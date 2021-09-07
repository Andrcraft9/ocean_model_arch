#include "macros/mpp_macros.fi"

program model
    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use mpp_module
    use mpp_sync_module, only: mpp_sync_init, mpp_sync_finalize
    use config_basinpar_module, only: load_config_basinpar, load_config_basinpar_from_file
    use config_parallel_module, only: load_config_parallel_from_file_and_cmd, dlb_balance_steps, dlb_model_steps
    use config_sw_module, only: load_config_sw
    use time_manager_module, only: num_step, num_step_max, tau,  &
                                   init_time_manager, time_manager_def, time_manager_print, is_local_print_step,  &
                                   year_loc, mon_loc, day_loc, hour_loc, min_loc, nrec_loc,  &
                                   loc_data_tstep, time_manager_update_nrec, yr_type
    use decomposition_module, only: domain_data
    use grid_module, only: grid_global_data, grid_data
    use ocean_module, only: ocean_data
    use io_module, only: read_global_mask
    use init_data_module, only: init_grid_data, init_ocean_data, init_device_data
    !use ocean_interface_module, only: envoke_div_velocity
    use output_module, only: local_output, output_init_buffers, output_clear_buffers
    use shallow_water_module, only: expl_shallow_water
    use tracer_control_module, only: expl_tracer
    use preprocess_module, only: dynamic_load_balance
    use errors_module
#ifdef _GPU_MODE_
    use shallow_water_gpu_module, only: expl_shallow_water_gpu
#endif

    implicit none

    real(wp8) :: t_local
    real(wp4) :: dev_local
    integer(wp8) :: local_num_step
#ifdef _GPU_MODE_
    integer :: istat
#endif

    call mpp_init()

    ! Load all configs
    call load_config_basinpar_from_file('basin.par')
    call load_config_sw()
    call load_config_parallel_from_file_and_cmd('parallel.par')

    ! Init Time Manager
    call init_time_manager('ocean_run.par')
    local_num_step = num_step

    ! Read global sea-land mask
    call grid_global_data%init()
    call read_global_mask(grid_global_data)

    ! Make decomposition and allocate data
    if (dlb_balance_steps <= 0 .or. dlb_model_steps <= 0) then
        call domain_data%init_from_config(grid_global_data%mask)

        ! Init sync buffers and patterns
        call mpp_sync_init(domain_data)
        
        ! Allocate data
        call ocean_data%init(domain_data)
        call grid_data%init(domain_data)
        
        ! Init data (read/set)
        call init_grid_data()
        call init_ocean_data()
    else
        ! Try to choose optimal decomposition
        if (mpp_is_master())  then
            print *, "MODEL: Start Dynamic Load Balance"
            call start_timer(t_local)
        endif
        call dynamic_load_balance(tau, dlb_balance_steps, dlb_model_steps)
        if (mpp_is_master()) then
            call end_timer(t_local)
            mpp_time_load_balance =  t_local
            print *, "MODEL: Dynamic Load Balance Time: ", mpp_time_load_balance
        endif
    endif
    ! Init other buffers
    call output_init_buffers()

    ! Device init data
    if (mpp_is_master_thread()) then
        call start_device_timer()
    endif
    call init_device_data
    if (mpp_is_master_thread()) then
        call end_device_timer(dev_local)
        mpp_device_init_data = mpp_device_init_data + dev_local
    endif

    call time_manager_def(local_num_step)
    if (mpp_is_master()) then
        print *,  '=================================================================='
        print *,  '------------ Eplicit shallow water scheme, version PSyKAl  -------'
        print *,  '=================================================================='
        print *,  '=================================================================='
        print *,  '----------- Starting shallow water model time integration --------'
        print *,  '=================================================================='
    endif

    call mpp_sync_output()

    if (is_local_print_step(local_num_step) > 0) then
        if (mpp_is_master()) print *, "MODEL: Output initial local data..."
        call local_output(1,  &
                          year_loc,  &
                          mon_loc,  &
                          day_loc,  &
                          hour_loc,  &
                          min_loc,  &
                          loc_data_tstep,  &
                          yr_type  )

        call time_manager_print(local_num_step)
    endif

    !debug
    !call abort_model("Stop")

    !$omp parallel default(shared) firstprivate(local_num_step, num_step_max, tau)

    ! Solver
    do while(local_num_step < num_step_max)
        if (mpp_is_master_thread()) then
            call start_timer(t_local)
            call start_device_timer()
        endif

        ! Computing one step of ocean dynamics
#ifdef _GPU_MODE_
        call expl_shallow_water_gpu(tau)
        istat = cudaDeviceSynchronize()
#else
        call expl_shallow_water(tau)
        call expl_tracer(tau)
#endif
        
        if (mpp_is_master_thread()) then
            call end_timer(t_local)
            call end_device_timer(dev_local)
            mpp_time_model_step = mpp_time_model_step + t_local
            mpp_device_time_model_step = mpp_device_time_model_step + dev_local
        endif

        ! Next time step
        local_num_step = local_num_step + 1
        if (mpp_is_master_thread()) then
            
            call time_manager_def(local_num_step)

            ! Output
            if (is_local_print_step(local_num_step) > 0) then
                call time_manager_update_nrec(local_num_step)
              
#ifdef _GPU_MODE_
                call ocean_data%ssh%sync_host_device(domain_data, .false.)
#endif
                call local_output(nrec_loc + 1,  &
                                  year_loc,  &
                                  mon_loc,  & 
                                  day_loc,  &
                                  hour_loc,  &
                                  min_loc,  &
                                  loc_data_tstep,  &
                                  yr_type  )

                call time_manager_print(local_num_step)
            endif
        endif

        !$omp barrier

    enddo

    !$omp end parallel

    ! Clear data
    call ocean_data%clear(domain_data)
    call grid_data%clear(domain_data)
    call grid_global_data%clear()
    call output_clear_buffers()

    ! Clear sync data
    call mpp_sync_finalize(domain_data)

    ! Clear decomposition
    call domain_data%clear()

    call mpp_finalize()

end program model