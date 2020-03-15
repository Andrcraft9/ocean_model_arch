module ocean_interface_module
    ! Module description

    use decomposition_module, only: domain_type
    use ocean_module, only: ocean_type
    use grid_module, only: grid_type
    use kernel_interface_module, only: set_kernel_interface
    use velocity_module, only: div_velocity_kernel

    implicit none
    save
    private

    public :: envoke_div_velocity

contains

    subroutine envoke_div_velocity(domain, grid_data, ocean_data)
        ! Subroutine description
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(in) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data
        integer :: k

        do k = 1, domain%bcount
            call set_kernel_interface(domain, k)
            associate(  &
                ! Grid data
                lu => grid_data%lu%block(k)%field,    &
                !hhq => grid_data%hhq%block(k)%field,  &
                hu => grid_data%hhu%block(k)%field,  &
                hv => grid_data%hhv%block(k)%field,  &
                dxh => grid_data%dxt%block(k)%field,  &
                dyh => grid_data%dyt%block(k)%field,  &
                sqt => grid_data%sqt%block(k)%field,  &
                ! Ocean data
                ssh => ocean_data%ssh%block(k)%field,      &
                ub => ocean_data%ubrtr%block(k)%field,  &
                vb => ocean_data%vbrtr%block(k)%field,  &
                div_btr => ocean_data%div_btr%block(k)%field)

                ! Here kernel or kernels call
                call div_velocity_kernel(lu, dxh, dyh, sqt, hu, hv, ub, vb, div_btr)

            end associate
        enddo
        
    end subroutine envoke_div_velocity

end module ocean_interface_module