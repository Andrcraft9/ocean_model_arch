module kernel_interface_module

    use decomposition_module, only: domain_type
    use kernel_runtime_module
    use mpp_module

    implicit none
    save
    public

    type, public :: block_bounds_type
        integer :: bnd_x1, bnd_x2, bnd_y1, bnd_y2
        integer :: nx_start, nx_end, ny_start, ny_end
    contains
        procedure, public :: set => set_block_bounds
    end type

    real(wp8) :: kernel_time_local

contains

    subroutine set_block_bounds(this, domain, k)
        class(block_bounds_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        integer, intent(in) :: k

        this%bnd_x1 = domain%bbnd_x1(k)
        this%bnd_x2 = domain%bbnd_x2(k)
        this%bnd_y1 = domain%bbnd_y1(k)
        this%bnd_y2 = domain%bbnd_y2(k)

        this%nx_start = domain%bnx_start(k)
        this%nx_end   = domain%bnx_end(k)
        this%ny_start = domain%bny_start(k)
        this%ny_end   = domain%bny_end(k)
    end subroutine

    subroutine end_kernel_timer(kernel_name)
        character(*), intent(in) :: kernel_name
        integer :: kernel_id

        call end_timer(kernel_time_local)
        call get_kernel_id(kernel_name, kernel_id)
        mpp_calls_kernels(kernel_id) = mpp_calls_kernels(kernel_id) + 1
        mpp_time_kernels(kernel_id) = mpp_time_kernels(kernel_id) + kernel_time_local
    end subroutine

    subroutine start_kernel_timer(kernel_name)
        character(*), intent(in) :: kernel_name
        integer :: kernel_id

        call start_timer(kernel_time_local)

    end subroutine

end module kernel_interface_module