module kernel_interface_module

    use decomposition_module, only: domain_type
    use kernel_runtime_module
    use mpp_module

#include "macros/mpp_macros.fi"

    implicit none
    save
    public

    real(wp8) :: kernel_time_local
    real(wp8) :: kernel_time_local_threads(0 : _OMP_MAX_THREADS_ - 1)

contains

    subroutine end_kernel_timer(kernel_name)
        character(*), intent(in) :: kernel_name
        integer :: kernel_id

#ifdef _MPP_KERNEL_TIMER_ON_
        ! Master thread
        if (mpp_is_master_thread()) then
            call end_timer(kernel_time_local)
            call get_kernel_id(kernel_name, kernel_id)
            mpp_calls_kernels(kernel_id) = mpp_calls_kernels(kernel_id) + 1
            mpp_time_kernels(kernel_id) = mpp_time_kernels(kernel_id) + kernel_time_local
        endif

        ! Threads
        call end_timer(kernel_time_local_threads(mpp_get_thread()))
        call get_kernel_id(kernel_name, kernel_id)
        mpp_time_kernels_threads(mpp_get_thread(), kernel_id) = mpp_time_kernels_threads(mpp_get_thread(), kernel_id) + kernel_time_local_threads(mpp_get_thread())
#endif
    end subroutine

    subroutine start_kernel_timer(kernel_name)
        character(*), intent(in) :: kernel_name
        integer :: kernel_id

#ifdef _MPP_KERNEL_TIMER_ON_
        ! Master thread
        if (mpp_is_master_thread()) call start_timer(kernel_time_local)

        ! Threads
        call start_timer(kernel_time_local_threads(mpp_get_thread()))
#endif
    end subroutine

end module kernel_interface_module