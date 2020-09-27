module kernel_interface_module

    use kernel_runtime_module
    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use mpp_module
    use decomposition_module, only: domain_type
    use data_types_module, only: data2D_real8_type, data2D_real4_type
    use ocean_module, only: ocean_type
    use grid_module, only: grid_type
    use mpp_sync_module, only: hybrid_sync, sync, sync_parameters_type

#include "macros/mpp_macros.fi"

    implicit none
    save
    public

    ! Kernel timers for _MPP_KERNEL_TIMER_ON_
    real(wp8) :: kernel_time_local
    real(wp8) :: kernel_time_local_threads(0 : _OMP_MAX_THREADS_ - 1)

contains

!-----------------------------------------------------------------------------!
!-------------------------- Interface subroutines ----------------------------!
!-----------------------------------------------------------------------------!
    subroutine envoke_empty_kernel(k, domain, grid_data, ocean_data, param)
        integer, intent(in) :: k
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data
        real(wp8), intent(in) :: param
    end subroutine 

    subroutine envoke_empty_sync(sync_parameters, domain, grid_data, ocean_data)
        type(sync_parameters_type), intent(in) :: sync_parameters
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data
    end subroutine 

    subroutine envoke(domain, grid_data, ocean_data, sub_kernel, sub_sync, param)
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data
        procedure(envoke_empty_kernel), pointer :: sub_kernel
        procedure(envoke_empty_sync), pointer :: sub_sync
        real(wp8), intent(in) :: param
    
        integer :: k
        type(sync_parameters_type) :: sync_parameters_inner, sync_parameters_boundary, sync_parameters_intermediate, sync_parameters_all

        sync_parameters_inner%sync_mode = 0
        sync_parameters_boundary%sync_mode = 1
        sync_parameters_intermediate%sync_mode = 2
        sync_parameters_all%sync_mode = 3

#ifdef _MPP_NO_PARALLEL_MODE_

        do k = 1, domain%bcount
            call sub_kernel(k, domain, grid_data, ocean_data, param)
        enddo

        call sub_sync(sync_parameters_all, domain, grid_data, ocean_data)

#endif

#ifdef _MPP_BLOCK_MODE_

        !$omp do private(k) schedule(static, 1)
        do k = 1, domain%bcount
            call sub_kernel(k, domain, grid_data, ocean_data, param)
        enddo
        !$omp end do nowait

        call sub_sync(sync_parameters_all, domain, grid_data, ocean_data)

#endif

#ifdef _MPP_HYBRID_BLOCK_MODE_
    
        !$omp do private(k) schedule(static, 1)
        do k = 1, domain%bcount
            call sub_kernel(k, domain, grid_data, ocean_data, param)
        enddo
        !$omp end do nowait

        call sub_sync(sync_parameters_boundary, domain, grid_data, ocean_data)
        call sub_sync(sync_parameters_inner, domain, grid_data, ocean_data)
        call sub_sync(sync_parameters_intermediate, domain, grid_data, ocean_data)

#endif

    end subroutine

!-----------------------------------------------------------------------------!
!-------------------------- Timers subroutines -------------------------------!
!-----------------------------------------------------------------------------!
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