program model
    use parallel_module
    use decomposition_module
    use grid_module
    use ocean_module
    use init_data_module

    implicit none

    ! Make decomposition
    call procs_data%init()
    call domain_data%init(procs_data, 2, 2)

    ! Allocate data
    call ocean_data%init(domain_data)
    call grid_data%init(domain_data)

    ! Init data (read/set)
    call init_ocean_data(procs_data, domain_data, ocean_data)
    call init_grid_data(procs_data, domain_data, grid_data)

    ! Solver

    ! Clear data
    call ocean_data%clear(domain_data)
    call grid_data%clear(domain_data)

    ! Clear decomposition
    call domain_data%clear()
    call procs_data%finalize()

end program model