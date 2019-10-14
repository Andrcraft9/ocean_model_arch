module ocean_module
    ! Module description

    use decomposition_module, only: domain_type
    use ocean_module, only: ocean_type
    use grid_module, only: grid_type

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
            associate(
                ! Grid data
                lu => grid_data%lu%block(k)%field,    &
                hhq => grid_data%hhq%block(k)%field,  &
                hhu => grid_data%hhu%block(k)%field,  &
                hhv => grid_data%hhv%block(k)%field,  &
                dxt => grid_data%dxt%block(k)%field,  &
                dyt => grid_data%dyt%block(k)%field,  &
                sqt => grid_data%sqt%block(k)%field,  &
                ! Ocean data
                ssh => ocean_data%ssh%block(k)%field,      &
                ubrtr => ocean_data%ubrtr%block(k)%field,  &
                vbrtr => ocean_data%vbrtr%block(k)%field,  &
                div_btr = > ocean_data%div_btr%block(k)%field)

                call div_velocity_kernel()
            end associate
        enddo
        
    end subroutine envoke_div_velocity

end module ocean_module