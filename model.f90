program model
    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use mpp_module, only: mpp_init, mpp_finalize
    use mpp_sync_module, only: mpp_sync_init, mpp_sync_finalize
    use config_basinpar_module, only: load_config_basinpar
    use decomposition_module, only: domain_data
    use grid_module, only: grid_global_data, grid_data
    use ocean_module, only: ocean_data
    use io_module, only: read_global_mask
    use init_data_module, only: init_grid_data, init_ocean_data
    !use ocean_interface_module, only: envoke_div_velocity
    use output_module, only: local_output
    use shallow_water_module, only: expl_shallow_water

    implicit none

    integer, parameter :: num_step_max = 10
    real(wp8), parameter :: tau = 1
    integer :: num_step = 0

    call mpp_init()

    ! Load all configs
    call load_config_basinpar()

    ! Read mask
    call grid_global_data%init()
    call read_global_mask(grid_global_data)

    ! Make decomposition
    call domain_data%init(2, 2, grid_global_data%mask)

    ! Sync init
    call mpp_sync_init(domain_data)

    ! Allocate data
    call ocean_data%init(domain_data)
    call grid_data%init(domain_data)

    ! Init data (read/set)
    call init_ocean_data(domain_data, ocean_data)
    call init_grid_data(domain_data, grid_global_data, grid_data)

    ! Solver
    do while(num_step<num_step_max)
        ! Computing one step of ocean dynamics
        print *, "STEP: ", num_step
        call expl_shallow_water(tau, domain_data, grid_data, ocean_data)
        num_step = num_step + 1
    enddo
    

    ! Output
    call local_output(domain_data, grid_data, ocean_data, 1, 2019, 2, 10, 12, 30, 3600.0, 1)

    ! Clear data
    call ocean_data%clear(domain_data)
    call grid_data%clear(domain_data)
    call grid_global_data%clear()

    ! Clear sync data
    call mpp_sync_finalize(domain_data)

    ! Clear decomposition
    call domain_data%clear()

    call mpp_finalize()

end program model