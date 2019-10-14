module kernel_interface_module
    ! Module description

    use decomposition_module, only: domain_type

    implicit none
    save
    private

    integer, public :: nxs, nxe
    integer, public :: nys, nye

    public :: set_kernel_interface

contains

    subroutine set_kernel_interface(domain, k)
        ! Subroutine description
        type(domain_type), intent(in) :: domain
        integer, intent(in) :: k

        nxs = domain%bnx_start(k)
        nxe = domain%bnx_end(k)
        nys = domain%bny_start(k)
        nye = domain%bny_end(k)
    end subroutine 

end module kernel_interface_module