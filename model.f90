program model
    use mpp_module, only: mpp_init, mpp_finalize
    use decomposition_module, only: domain_data
    use grid_module, only: grid_data
    use ocean_module, only: ocean_data
    use init_data_module, only: init_grid_data, init_ocean_data
    use ocean_model_module, only: envoke_div_velocity
    use output_module, only: local_output

    implicit none

    call mpp_init()

    ! Make decomposition
    call domain_data%init(2, 2)

    ! Allocate data
    call ocean_data%init(domain_data)
    call grid_data%init(domain_data)

    ! Init data (read/set)
    call init_ocean_data(domain_data, ocean_data)
    call init_grid_data(domain_data, grid_data)

    ! Solver
    call envoke_div_velocity(domain_data, grid_data, ocean_data)

    ! Output
    call local_output(domain_data, grid_data, ocean_data)

    ! Clear data
    call ocean_data%clear(domain_data)
    call grid_data%clear(domain_data)

    ! Clear decomposition
    call domain_data%clear()

    call mpp_finalize()

end program model