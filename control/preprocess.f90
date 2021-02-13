module preprocess_module

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use mpp_module
    use mpp_sync_module, only: mpp_sync_init, mpp_sync_finalize
    use decomposition_module, only: domain_data
    use grid_module, only: grid_global_data, grid_data
    use ocean_module, only: ocean_data
    use init_data_module, only: init_grid_data, init_ocean_data
    use shallow_water_module, only: shallow_water_kernels_compute_probe
    use errors_module
    
    implicit none
    save
    private

    public :: dynamic_load_balance

contains

    subroutine dynamic_load_balance(tau, balance_steps, model_steps)
        real(wp8), intent(in) :: tau
        integer, intent(in) :: balance_steps, model_steps
        
        integer :: balance_step, model_step
        real(wp8), allocatable :: compute_power(:)
        real(wp8) :: t_local
        integer :: ierr

        ! Make Init decomposition
        call domain_data%init_from_config(grid_global_data%mask)
        ! Init sync buffers and patterns
        call mpp_sync_init(domain_data)
        ! Allocate data
        call ocean_data%init(domain_data)
        call grid_data%init(domain_data)
        ! Init data (read/set)
        call init_grid_data()
        call init_ocean_data()

        allocate(compute_power(mpp_count))
        compute_power = 0.0
        ierr = 0

        do balance_step = 1, balance_steps

            ! Reset timers
            call mpp_reset()

            !$omp parallel default(shared) firstprivate(tau)

            do model_step = 1, model_steps
                if (mpp_is_master_thread())  call start_timer(t_local)
                ! Computing one step of ocean dynamics
                call shallow_water_kernels_compute_probe(tau)
                if (mpp_is_master_thread()) then
                    call end_timer(t_local)
                    mpp_time_model_step = mpp_time_model_step + t_local
                endif
                !$omp barrier
            enddo
    
            !$omp end parallel

            ! Fill Compute Powers of procs
            if (mpp_time_sync > 0.0 .and. mpp_time_model_step < 0.001) then
                ierr = 1
                print *, mpp_rank, 'Time sync and time model set are:', mpp_time_sync, mpp_time_model_step
            endif
            call check_error(ierr, 'PREP: ERROR: Cant compute Compute Powers of procs!')
            call mpi_allgather(domain_data%tot_weight / mpp_time_model_step, 1, mpi_real8, compute_power, 1, mpi_real8, mpp_cart_comm, ierr)
            compute_power = compute_power / MAXVAL(compute_power)
            if (mpp_is_master()) then
                print *, "PREP:  Dynamic Load Balance step is ", balance_step, " Compute Powers of procs (normalized): ", compute_power
            endif

            ! Clear data which are dependent on domain_data - this data must rebuild
            call ocean_data%clear(domain_data)
            call grid_data%clear(domain_data)
            ! Clear sync data
            call mpp_sync_finalize(domain_data)
            ! Clear decomposition
            call domain_data%clear()
            ! Reset timers
            call mpp_reset()

             ! Make new decomposition
            call domain_data%init_from_config_and_mpp_compute_power(grid_global_data%mask, compute_power)
            ! Init sync buffers and patterns
            call mpp_sync_init(domain_data)
            ! Allocate data
            call ocean_data%init(domain_data)
            call grid_data%init(domain_data)
            ! Init data (read/set)
            call init_grid_data()
            call init_ocean_data()
        enddo

        deallocate(compute_power)
    end subroutine

end module preprocess_module