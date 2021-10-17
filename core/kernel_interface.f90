module kernel_interface_module

    use kernel_runtime_module
    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use mpp_module
    use mpp_sync_module, only: hybrid_sync, sync, sync_parameters_type
    use decomposition_module, only: domain_type, domain => domain_data

#include "macros/mpp_macros.fi"

    implicit none
    save
    public

    type, public :: kernel_parameters_type
        real(wp8) :: tau
        real(wp8) :: time_smooth
        integer :: data_id
    contains
        procedure, public  :: clear => clear_kernel_parameters_type
    end type kernel_parameters_type

    ! Kernel timers for _MPP_KERNEL_TIMER_ON_
    real(wp8) :: kernel_time_local
    real(wp8) :: kernel_time_local_threads(0 : _OMP_MAX_THREADS_ - 1)

contains

    subroutine clear_kernel_parameters_type(this)
        class(kernel_parameters_type), intent(inout) :: this
        this%tau = 0.0d0
        this%time_smooth = 0.0d0
        this%data_id = 0
    end subroutine
!-----------------------------------------------------------------------------!
!-------------------------- Interface subroutines ----------------------------!
!-----------------------------------------------------------------------------!
    subroutine envoke_empty_kernel(k, kernel_parameters)
        integer, intent(in) :: k
        type(kernel_parameters_type), intent(in) :: kernel_parameters
    end subroutine 

    subroutine envoke_empty_sync(k, sync_parameters)
        integer, intent(in) :: k
        type(sync_parameters_type), intent(in) :: sync_parameters
    end subroutine 

    subroutine envoke(sub_kernel, sub_sync, kernel_parameters)
        procedure(envoke_empty_kernel), pointer :: sub_kernel
        procedure(envoke_empty_sync), pointer :: sub_sync
        type(kernel_parameters_type), intent(in) :: kernel_parameters
    
        integer :: k, ierr
        type(sync_parameters_type) :: sync_parameters_inner, sync_parameters_boundary, sync_parameters_intermediate, sync_parameters_all
        real(wp8) :: t_local

        sync_parameters_inner%sync_mode = 0;         sync_parameters_inner%data_id = kernel_parameters%data_id
        sync_parameters_boundary%sync_mode = 1;      sync_parameters_boundary%data_id = kernel_parameters%data_id
        sync_parameters_intermediate%sync_mode = 2;  sync_parameters_intermediate%data_id = kernel_parameters%data_id
        sync_parameters_all%sync_mode = 3;           sync_parameters_all%data_id = kernel_parameters%data_id

#ifdef _MPP_NO_PARALLEL_MODE_

        do k = 1, domain%bcount
            call sub_kernel(k, kernel_parameters)
        enddo

        call sub_sync(-1, sync_parameters_all)

#endif

#ifdef _MPP_BLOCK_MODE_

#ifdef _DBG_TIME_PROFILE_
        !$omp master
        call mpi_barrier(mpp_cart_comm, ierr)
        !$omp end master
        !$omp barrier
        !$omp master
        call start_timer(t_local)
        !$omp end master
#endif

        !$omp do private(k) schedule(static, 1)
        do k = 1, domain%bcount
            call sub_kernel(k, kernel_parameters)
        enddo
        !$omp end do nowait

#ifdef _DBG_TIME_PROFILE_
        !$omp master
        call mpi_barrier(mpp_cart_comm, ierr)
        !$omp end master
        !$omp barrier
        !$omp master
        call end_timer(t_local)
        mpp_time_kernel_compute = mpp_time_kernel_compute + t_local
        !$omp end master
#endif

        call sub_sync(-1, sync_parameters_all)

#endif

#ifdef _MPP_HYBRID_BLOCK_MODE_
    
        !$omp do private(k) schedule(static, 1)
        do k = 1, domain%bcount
            call sub_kernel(k, kernel_parameters)
        enddo
        !$omp end do nowait

        call sub_sync(-1, sync_parameters_boundary)
        call sub_sync(-1, sync_parameters_inner)
        call sub_sync(-1, sync_parameters_intermediate)

#endif

    end subroutine

    subroutine envoke_device(sub_kernel, sub_sync, kernel_parameters)
        procedure(envoke_empty_kernel), pointer :: sub_kernel
        procedure(envoke_empty_sync), pointer :: sub_sync
        type(kernel_parameters_type), intent(in) :: kernel_parameters
        
        integer :: k, istat, ierr
        type(sync_parameters_type) :: sync_params_htod, sync_params_dtoh
        type(sync_parameters_type) :: sync_parameters_all
        real(wp8) :: t_local

#ifdef _GPU_MODE_

        sync_params_dtoh%sync_mode = 5; sync_params_dtoh%sync_device_host = 0; sync_params_dtoh%data_id = kernel_parameters%data_id
        sync_params_htod%sync_mode = 5; sync_params_htod%sync_device_host = 1; sync_params_htod%data_id = kernel_parameters%data_id
        sync_parameters_all%sync_mode = 3; sync_parameters_all%data_id = kernel_parameters%data_id

#ifdef _GPU_ASYNC_
        sync_params_dtoh%sync_device_host = 2
        sync_params_htod%sync_device_host = 3
#endif

#ifdef _GPU_FULL_
        do k = 1, domain%bcount
            call sub_kernel(k, kernel_parameters)
        enddo
#else
        ! If has GPU_ASYNC, GPU_MULTI or no modificator
        !

        ! Global debug sync before kernel compute
        !
#ifdef _DBG_TIME_PROFILE_
#ifdef _GPU_MULTI_
        !$omp do private(k) schedule(static, 1)
        do k = 1, domain%bcount
            istat = cudaSetDevice(k-1)
            istat = cudaDeviceSynchronize()
        enddo
        !$omp end do
#else
        !$omp master
        istat = cudaDeviceSynchronize()
        !$omp end master
#endif
        !$omp master
        call mpi_barrier(mpp_cart_comm, ierr)
        !$omp end master
        !$omp barrier
        !$omp master
        call start_timer(t_local)
        !$omp end master
#endif

        ! Kernel compute stage
        !
        !$omp do private(k) schedule(static, 1)
        do k = 1, domain%bcount
#ifdef _GPU_MULTI_
            istat = cudaSetDevice(k-1)
#endif
            call sub_kernel(k, kernel_parameters)
            call sub_sync(k, sync_params_dtoh)
        enddo
        !$omp end do nowait

        ! Global debug sync after kernel compute
        !
#ifdef _DBG_TIME_PROFILE_
#ifdef _GPU_MULTI_
        !$omp do private(k) schedule(static, 1)
        do k = 1, domain%bcount
            istat = cudaSetDevice(k-1)
            istat = cudaDeviceSynchronize()
        enddo
        !$omp end do
#else
        !$omp master
        istat = cudaDeviceSynchronize()
        !$omp end master
#endif
        !$omp master
        call mpi_barrier(mpp_cart_comm, ierr)
        !$omp end master
        !$omp barrier
        !$omp master
        call end_timer(t_local)
        mpp_time_kernel_compute = mpp_time_kernel_compute + t_local
        !$omp end master
#endif

        ! Sync devices before global sync call stage
        !
#ifdef _GPU_MULTI_
        !$omp do private(k) schedule(static, 1)
        do k = 1, domain%bcount
            istat = cudaSetDevice(k-1)
            istat = cudaDeviceSynchronize()
        enddo
        !$omp end do
#else
        !$omp master
        istat = cudaDeviceSynchronize()
        !$omp end master
        !$omp barrier
#endif

        ! Global sync call stage
        !
        call sub_sync(-1, sync_parameters_all)

        ! Global debug sync before host to device sync stage
        !
#ifdef _DBG_TIME_PROFILE_
#ifdef _GPU_MULTI_
        !$omp do private(k) schedule(static, 1)
        do k = 1, domain%bcount
            istat = cudaSetDevice(k-1)
            istat = cudaDeviceSynchronize()
        enddo
        !$omp end do
#else
        !$omp master
        istat = cudaDeviceSynchronize()
        !$omp end master
#endif
        !$omp master
        call mpi_barrier(mpp_cart_comm, ierr)
        !$omp end master
        !$omp barrier
        !$omp master
        call start_timer(t_local)
        !$omp end master
#endif

        ! Host to device sync stage
        !
        !$omp do private(k) schedule(static, 1)
        do k = 1, domain%bcount
#ifdef _GPU_MULTI_
            istat = cudaSetDevice(k-1)
#endif
            call sub_sync(k, sync_params_htod)
        enddo
        !$omp end do nowait

        ! Global debug sync after host to device sync stage
        !
#ifdef _DBG_TIME_PROFILE_
#ifdef _GPU_MULTI_
        !$omp do private(k) schedule(static, 1)
        do k = 1, domain%bcount
            istat = cudaSetDevice(k-1)
            istat = cudaDeviceSynchronize()
        enddo
        !$omp end do
#else
        !$omp master
        istat = cudaDeviceSynchronize()
        !$omp end master
#endif
        !$omp master
        call mpi_barrier(mpp_cart_comm, ierr)
        !$omp end master
        !$omp barrier
        !$omp master
        call end_timer(t_local)
        mpp_time_htod_sync = mpp_time_htod_sync + t_local
        !$omp end master
#endif

#endif

#endif
    end subroutine

    subroutine envoke_heterogeneous(sub_kernel, sub_sync, sub_kernel_gpu, sub_sync_gpu, kernel_parameters)
        procedure(envoke_empty_kernel), pointer :: sub_kernel
        procedure(envoke_empty_sync), pointer :: sub_sync
        procedure(envoke_empty_kernel), pointer :: sub_kernel_gpu
        procedure(envoke_empty_sync), pointer :: sub_sync_gpu
        type(kernel_parameters_type), intent(in) :: kernel_parameters
        
        integer :: k, istat, cpu_count, gpu_count
        type(sync_parameters_type) :: sync_params_htod, sync_params_dtoh
        type(sync_parameters_type) :: sync_parameters_all

#ifdef _GPU_MODE_

        sync_params_dtoh%sync_mode = 5; sync_params_dtoh%sync_device_host = 0; sync_params_dtoh%data_id = kernel_parameters%data_id
        sync_params_htod%sync_mode = 5; sync_params_htod%sync_device_host = 1; sync_params_htod%data_id = kernel_parameters%data_id
        sync_parameters_all%sync_mode = 3; sync_parameters_all%data_id = kernel_parameters%data_id

#ifdef _GPU_ASYNC_
        sync_params_dtoh%sync_device_host = 2
        sync_params_htod%sync_device_host = 3
#endif

        cpu_count = int(_CPU_GPU_RATIO_ * domain%bcount)
        gpu_count = domain%bcount - cpu_count

        do k = 1, gpu_count
            call sub_kernel_gpu(k, kernel_parameters)
            call sub_sync_gpu(k, sync_params_dtoh)
        enddo
        
        do k = gpu_count + 1, domain%bcount
            call sub_kernel(k, kernel_parameters)
        enddo

        istat = cudaDeviceSynchronize()
        call sub_sync(-1, sync_parameters_all)

        do k = 1, gpu_count
            call sub_sync(k, sync_params_htod)
        enddo

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