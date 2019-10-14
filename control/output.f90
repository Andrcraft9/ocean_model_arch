module output_module
    ! Module description

    use decomposition_module, only: domain_type
    use ocean_module, only: ocean_type
    use grid_module, only: grid_type
    use io_module

    implicit none
    save
    private

    public :: local_output

contains

    subroutine local_output(domain, grid_data, ocean_data)
        ! Subroutine description
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data

        call write_2d_real8(domain, ocean_data%div_btr)
    end subroutine local_output

end module output_module