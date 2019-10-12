module kernel_interface_module
    ! Module description

    implicit none
    save
    private

    integer, public :: nx_start, nx_end
    integer, public :: ny_start, ny_end

    public :: set_kernel_interface

contains

    subroutine set_kernel_interface(procs, domain, k)
        ! Subroutine description
        type(procs_type), intent(in) :: procs
        type(domain_type), intent(in) :: domain
        integer, intent(in) :: k

        nx_start = domain%bnx_start(k)
        nx_end = domain%bnx_end(k)
        ny_start = domain%bny_start(k)
        ny_end = domain%bny_end(k)
    end subroutine 

end module kernel_interface_module