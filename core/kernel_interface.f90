module kernel_interface_module

    use decomposition_module, only: domain_type
    use kernel_runtime_module
    use mpp_module

    implicit none
    save
    public

    real(wp8) :: kernel_time_local

contains

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