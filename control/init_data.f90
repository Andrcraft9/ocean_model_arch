module init_data_module
    ! Grid data

    use parallel_module, only: procs_type
    use basinpar_module, only: dxst, dyst
    use decomposition_module, only: domain_type
    use ocean_module, only: ocean_type
    use grid_module, only: grid_type
    use io_module

    implicit none
    save
    private

    public :: init_ocean_data
    public :: init_grid_data

contains

    subroutine init_ocean_data(procs, domain, ocean_data)
        ! Subroutine description
        type(procs_type), intent(in) :: procs
        type(domain_type), intent(in) :: domain
        type(ocean_type), intent(inout) :: ocean_data
        integer :: k

        call read_2d_real8(procs, domain, ocean_data%ssh)
        call read_2d_real8(procs, domain, ocean_data%ubrtr)
        call read_2d_real8(procs, domain, ocean_data%vbrtr)
        
        do k = 1, domain%bcount
            associate(ssh => ocean_data%ssh%block(k)%field,      &
                      ubrtr => ocean_data%ubrtr%block(k)%field,  &
                      vbrtr => ocean_data%vbrtr%block(k)%field)
            
                ! some computational on data
                ! ...
            end associate
        enddo

    end subroutine init_ocean_data

    subroutine init_grid_data(procs, domain, grid_data)
        ! Subroutine description
        type(procs_type), intent(in) :: procs
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        integer :: k

        call read_2d_real4(procs, domain, grid_data%hhq)
        call read_2d_real4(procs, domain, grid_data%hhu)
        call read_2d_real4(procs, domain, grid_data%hhv)

        do k = 1, domain%bcount
            associate(lu => grid_data%lu%block(k)%field,    &
                      hhq => grid_data%hhq%block(k)%field,  &
                      hhu => grid_data%hhu%block(k)%field,  &
                      hhv => grid_data%hhv%block(k)%field,  &
                      dxt => grid_data%dxt%block(k)%field,  &
                      dyt => grid_data%dyt%block(k)%field,  &
                      sqt => grid_data%sqt%block(k)%field)

                lu = 1.0
                dxt = dxst
                dyt = dyst
                sqt = dxst*dyst

            end associate
        enddo
    end subroutine init_grid_data

end module init_data_module