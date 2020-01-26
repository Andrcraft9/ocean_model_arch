module output_module

    use decomposition_module, only: domain_type
    use ocean_module, only: ocean_type
    use grid_module, only: grid_type
    !use io_module, only:

    implicit none
    save
    private

    public :: local_output

contains

    subroutine local_output(domain, grid_data, ocean_data)
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(in) :: grid_data
        type(ocean_type), intent(in) :: ocean_data

        ! Nothing

    end subroutine local_output

end module output_module