program model
    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use mpp_module, only: mpp_init, mpp_finalize, mpp_rank, start_timer, end_timer, mpp_time_model_step
    use mpp_sync_module, only: mpp_sync_init, mpp_sync_finalize
    use config_basinpar_module, only: load_config_basinpar
    use config_sw_module, only: load_config_sw
    use time_manager_module, only: num_step, num_step_max, tau,  &
                                   init_time_manager, time_manager_def, time_manager_print, is_local_print_step,  &
                                   year_loc, mon_loc, day_loc, hour_loc, min_loc, nrec_loc,  &
                                   loc_data_tstep, time_manager_update_nrec, yr_type
    use decomposition_module, only: domain_data
    use grid_module, only: grid_global_data, grid_data
    use ocean_module, only: ocean_data
    use io_module, only: read_global_mask
    use init_data_module, only: init_grid_data, init_ocean_data
    !use ocean_interface_module, only: envoke_div_velocity
    use output_module, only: local_output, output_init_buffers, output_clear_buffers
    use shallow_water_module, only: expl_shallow_water

    implicit none

    real(wp8) :: t_local

    call mpp_init()

    ! Load all configs
    call load_config_basinpar()
    call load_config_sw()

    ! Init Time Manager
    call init_time_manager('ocean_run.par')

    ! Read mask
    call grid_global_data%init()
    call read_global_mask(grid_global_data)

    ! Make decomposition
    call domain_data%init_from_config(grid_global_data%mask, 'parallel.par')

    ! Sync init
    call mpp_sync_init(domain_data)

    ! Allocate data
    call ocean_data%init(domain_data)
    call grid_data%init(domain_data)
    call output_init_buffers(domain_data)

    ! Init data (read/set)
    call init_grid_data(domain_data, grid_global_data, grid_data)
    call init_ocean_data(domain_data, grid_data, ocean_data)

    call time_manager_def()
    if (mpp_rank .eq. 0) then
        print *,  '=================================================================='
        print *,  '------------ Eplicit shallow water scheme, version PSyKAl  -------'
        print *,  '=================================================================='
        print *,  '=================================================================='
        print *,  '----------- Starting shallow water model time integration --------'
        print *,  '=================================================================='
    endif

    if (is_local_print_step() > 0) then
        if (mpp_rank == 0) print *, "Output initial local data..."
        call local_output(domain_data, &
                          grid_data,   &
                          ocean_data,  &
                          1,  &
                          year_loc,  &
                          mon_loc,  &
                          day_loc,  &
                          hour_loc,  &
                          min_loc,  &
                          loc_data_tstep,  &
                          yr_type  )

        call time_manager_print ()
    endif

    ! Solver
    do while(num_step<num_step_max)
        call start_timer(t_local)
        ! Computing one step of ocean dynamics
        call expl_shallow_water(tau, domain_data, grid_data, ocean_data)
        call end_timer(t_local)
        mpp_time_model_step = mpp_time_model_step + t_local
        
        ! Next time step
        num_step = num_step + 1
        call time_manager_def()

        ! Output
        if (is_local_print_step() > 0) then
            call time_manager_update_nrec ()
            
            call local_output(domain_data, &
                              grid_data,   &
                              ocean_data,  &
                             nrec_loc + 1,  &
                             year_loc,  &
                              mon_loc,  & 
                              day_loc,  &
                             hour_loc,  &
                              min_loc,  &
                       loc_data_tstep,  &
                              yr_type  )

            call time_manager_print ()
        endif
    enddo

    ! Clear data
    call ocean_data%clear(domain_data)
    call grid_data%clear(domain_data)
    call grid_global_data%clear()
    call output_clear_buffers(domain_data)

    ! Clear sync data
    call mpp_sync_finalize(domain_data)

    ! Clear decomposition
    call domain_data%clear()

    call mpp_finalize()

end program model