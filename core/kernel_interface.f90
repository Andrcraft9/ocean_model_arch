module kernel_interface_module
    ! Module description

    use decomposition_module, only: domain_type

    implicit none
    save
    private

    integer, public :: bnd_x1, bnd_x2, bnd_y1, bnd_y2
    integer, public :: nx_start, nx_end, ny_start, ny_end

    public :: set_kernel_interface

contains

    subroutine set_kernel_interface(domain, k)
        ! Subroutine description
        type(domain_type), intent(in) :: domain
        integer, intent(in) :: k

        bnd_x1 = domain%bbnd_x1(k)
        bnd_x2 = domain%bbnd_x2(k)
        bnd_y1 = domain%bbnd_y1(k)
        bnd_y2 = domain%bbnd_y2(k)

        nx_start = domain%bnx_start(k)
        nx_end   = domain%bnx_end(k)
        ny_start = domain%bny_start(k)
        ny_end   = domain%bny_end(k)
    end subroutine 

end module kernel_interface_module