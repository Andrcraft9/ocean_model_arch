#include "macros/mpp_macros.fi"

module mpp_module
    ! MPI massively parallel processing library

    use mpi
    use kernel_runtime_module
    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use system_module, only: lrecl, lmpirecl, undef
#ifdef _GPU_MODE_
    use cudafor
#endif

    implicit none
    save
    public

!    include 'mpif.h'
    include "omp_lib.h"

    integer, public :: mpp_rank, mpp_count
    integer, public :: mpp_cart_comm
    integer, dimension(2), public :: mpp_size, mpp_coord
    logical, dimension(2), public :: mpp_period

    integer, public :: mpp_count_threads

    public :: mpp_sync_output
    public :: mpp_is_master, mpp_is_master_thread, mpp_is_master_process
    public :: start_timer, end_timer
    public :: start_device_timer, end_device_timer
    public :: mpp_init
    public :: mpp_reset
    public :: mpp_finalize

    ! Timers for master thread
    real(wp8) :: mpp_time_model_step, mpp_time_sync, mpp_time_sync_inner, mpp_time_sync_boundary, mpp_time_sync_intermediate
    real(wp8) :: mpp_time_sync_pack_mpi, mpp_time_sync_unpack_mpi, mpp_time_sync_isend_irecv, mpp_time_sync_wait
    real(wp8) :: mpp_time_load_balance
    real(wp8) :: mpp_time_sw, mpp_time_tracers
    ! Devcie timers
    real(wp4) :: mpp_device_time_model_step, mpp_device_init_data

    integer :: mpp_max_count_sync_send_recv, mpp_min_count_sync_send_recv
    real(wp8), allocatable :: mpp_time_kernels(:)
    integer, allocatable :: mpp_calls_kernels(:)

    ! Timers for threads
    real(wp8), allocatable :: mpp_time_kernels_threads(:, :)

    ! Memory profiling
    real(wp8) :: mpp_total_memory_allocations

#ifdef _GPU_MODE_
    type (cudaEvent), public :: device_start_event, device_stop_event
#endif
contains

    subroutine mpp_init()
        ! Subroutine description
        integer :: ierr, req, provided, tmp
        integer :: rank_cart
        character(len=mpi_max_library_version_string) :: version
#ifdef _GPU_MODE_
        type (cudaDeviceProp) :: prop_device
        integer :: nDevices = 0, idev
        integer :: istat
#endif

        mpp_period = (/.true., .true./)
        mpp_size = (/0,0/)
        !req = MPI_THREAD_MULTIPLE
        req = MPI_THREAD_FUNNELED
        ierr = 0

        !call mpi_init(ierr)
        call mpi_init_thread(req, provided, ierr)

        call mpi_get_library_version(version, tmp, ierr)

        call mpi_comm_rank(mpi_comm_world, mpp_rank, ierr)
        call mpi_comm_size(mpi_comm_world, mpp_count, ierr)
        call mpi_dims_create(mpp_count, 2, mpp_size, ierr)
        call mpi_cart_create(mpi_comm_world, 2, mpp_size, mpp_period, .false., mpp_cart_comm, ierr)
        call mpi_cart_coords(mpp_cart_comm, mpp_rank, 2, mpp_coord, ierr)
        
        call mpi_comm_rank(mpp_cart_comm, rank_cart, ierr)

        !$omp parallel
        mpp_count_threads = omp_get_num_threads()
        !$omp end parallel

        ! MPP INFO
        if (mpp_is_master()) then
            print *, "MPP INFO: Total Processes = ", mpp_count
            print *, "MPP INFO: Total Threads   = ", mpp_count_threads
            print *, "MPP INFO: Max possible threads per process = ", _OMP_MAX_THREADS_
            print *, "MPP INFO: required and provided thread level for mpi : ", req, provided
            print *, "MPP INFO: MPI byte, MPI real, MPI real8 : ", MPI_BYTE, MPI_REAL, MPI_REAL8
            print *, "MPP INFO: kind 4 and kind 8 : ", wp4, wp8
            print *, "MPP INFO: IO: lrecl, lmpirecl, undef : ", lrecl, lmpirecl, undef
            print *, "MPP INFO: MPI version: ", trim(version)

            print *, "MPP INFO: Check defines:"
#ifdef _MPP_KERNEL_TIMER_ON_
            print *, "MPP INFO: MPP_KERNEL_TIMER_ON: yes"
#endif

#ifdef _MPP_SORTED_BLOCKS_
            print *, "MPP INFO: MPP_SORTED_BLOCKS: yes"
#endif

#ifdef _MPP_NO_PARALLEL_MODE_
            print *, "MPP INFO: MPP_NO_PARALLEL_MODE: yes"
#endif

#ifdef _MPP_BLOCK_MODE_
            print *, "MPP INFO: MPP_BLOCK_MODE: yes"
#endif

#ifdef _MPP_HYBRID_BLOCK_MODE_
            print *, "MPP INFO: MPP_HYBRID_BLOCK_MODE: yes"
#endif

            print *, "MPP INFO: Version 3.0, only NO_PARALLEL, BLOCK mode:"
            print *, "MPP INFO: Single OMP PARALLEL, OMP for pack/unpack MPI buffers, nowait loops, one global sync at mpi_waitall"
            print *, "MPP INFO: dont use SIMUL syncs"
            print *, "MPP INFO: Tracers added"

#ifdef _GPU_MODE_
            ! CUDA INFO
            ! Number of CUDA-capable devices
            print *, "MPP INFO: GPU MODE: Block X, Block Y: ", _GPU_BLOCK_X_, _GPU_BLOCK_Y_
            print *, "MPP INFO: GPU INFO:"
            ierr = cudaGetDeviceCount(nDevices)
            if (nDevices == 0) then
                write(*,"(/,'No CUDA devices found',/)")
                stop
            else if (nDevices == 1) then
                write(*,"(/,'One CUDA device found',/)")
            else 
                write(*,"(/,i0,' CUDA devices found',/)") nDevices
            end if
            ! Loop over devices
            do idev = 0, nDevices-1
                write(*,"('Device Number: ',i0)") idev
                ierr = cudaGetDeviceProperties(prop_device, idev)
                if (ierr .eq. 0) then
                    write(*,"('  GetDeviceProperties for device ',i0,': Passed')") idev
                else
                    write(*,"('  GetDeviceProperties for device ',i0,': Failed')") idev
                endif
                ! General device info
                write(*,"('  Device Name: ',a)") trim(prop_device%name)
                write(*,"('  Compute Capability: ',i0,'.',i0)") prop_device%major, prop_device%minor
                write(*,"('  Number of Multiprocessors: ',i0)") prop_device%multiProcessorCount
                write(*,"('  Max Threads per Multiprocessor: ',i0)") prop_device%maxThreadsPerMultiprocessor
                write(*,"('  Global Memory (GB): ',f9.3,/)") prop_device%totalGlobalMem/1024.0**3
                ! Execution Configuration
                write(*,"('  Execution Configuration Limits')")
                write(*,"('    Max Grid Dims: ',2(i0,' x '),i0)") prop_device%maxGridSize
                write(*,"('    Max Block Dims: ',2(i0,' x '),i0)") prop_device%maxThreadsDim
                write(*,"('    Max Threads per Block: ',i0,/)") prop_device%maxThreadsPerBlock
                write(*,"('    Device Overlap: ',i0,/)") prop_device%deviceOverlap
            enddo
#endif

#ifdef _GPU_ASYNC_
            print *, "MPP INFO: GPU ASYNC: yes"
#endif
#ifdef _GPU_FULL_
            print *, "MPP INFO: GPU FULL: yes"
#endif
#ifdef _GPU_CPU_HETERO_
            print *, "MPP INFO: GPU_CPU_HETERO: yes. CPU_GPU_RATIO is ", _CPU_GPU_RATIO_
#endif
            print *, "------------------------------------------------------------"
        endif
        call mpi_barrier(mpp_cart_comm, ierr)

        !$omp parallel
        if (omp_get_thread_num() == 0) then 
            print *, 'rank_world, rank_cart, cord, provided, threads:', mpp_rank, rank_cart, mpp_coord, provided, omp_get_num_threads()
        endif
        !$omp end parallel
        call mpi_barrier(mpp_cart_comm, ierr)

        ! Timers
        call mpp_reset()

        ! Memory profile
        mpp_total_memory_allocations = 0.0d0

#ifdef _GPU_MODE_
        istat = cudaEventCreate(device_start_event)
        istat = cudaEventCreate(device_stop_event)
#endif

        call mpp_sync_output()
    end subroutine

    subroutine mpp_reset()
        ! Rest Timers
        mpp_time_model_step = 0.0d0
        mpp_time_sync = 0.0d0
        mpp_time_sync_inner = 0.0d0
        mpp_time_sync_boundary = 0.0d0
        mpp_time_sync_intermediate = 0.0d0

        mpp_time_sync_pack_mpi = 0.0d0
        mpp_time_sync_unpack_mpi = 0.0d0
        mpp_time_sync_isend_irecv = 0.0d0
        mpp_time_sync_wait = 0.0d0

        mpp_time_load_balance = 0.0d0

        mpp_time_sw = 0.0d0
        mpp_time_tracers = 0.0d0

        mpp_device_time_model_step = 0.0
        mpp_device_init_data = 0.0

        mpp_max_count_sync_send_recv = 0
        mpp_min_count_sync_send_recv = 0

#ifdef _MPP_KERNEL_TIMER_ON_
        mpp_time_kernels = 0.0d0
        mpp_calls_kernels = 0

        mpp_time_kernels_threads = 0
#endif
    end subroutine

    subroutine mpp_finalize()
        integer :: k, t, ierr
        real(wp8) :: maxtime_model_step, mintime_model_step
        real(wp8) :: maxtime_sync, mintime_sync
        real(wp8) :: maxtime_kernel, mintime_kernel
        real(wp8) :: maxtime_kernel_threads(0 : _OMP_MAX_THREADS_ - 1), mintime_kernel_threads(0 : _OMP_MAX_THREADS_ - 1)
        integer :: maxcount, mincount
        character(80) :: kernel_name
#ifdef _GPU_MODE_
        integer :: istat
#endif

        if (mpp_is_master()) then
            write(*,'(a50)') 'Times for master thread:'
        endif

        ! Timers for master thread
        call mpi_allreduce(mpp_time_model_step, maxtime_model_step, 1, mpi_real8, mpi_max, mpp_cart_comm, ierr)
        call mpi_allreduce(mpp_time_model_step, mintime_model_step, 1, mpi_real8, mpi_min, mpp_cart_comm, ierr)
        if (mpp_rank == 0) write(*,'(a50, F12.2, F12.2)') "Time full of model step (max and min): ", maxtime_model_step, mintime_model_step

        call mpi_allreduce(mpp_time_sw, maxtime_model_step, 1, mpi_real8, mpi_max, mpp_cart_comm, ierr)
        call mpi_allreduce(mpp_time_sw, mintime_model_step, 1, mpi_real8, mpi_min, mpp_cart_comm, ierr)
        if (mpp_rank == 0) write(*,'(a50, F12.2, F12.2)') "Time shallow water step (max and min): ", maxtime_model_step, mintime_model_step

        call mpi_allreduce(mpp_time_tracers, maxtime_model_step, 1, mpi_real8, mpi_max, mpp_cart_comm, ierr)
        call mpi_allreduce(mpp_time_tracers, mintime_model_step, 1, mpi_real8, mpi_min, mpp_cart_comm, ierr)
        if (mpp_rank == 0) write(*,'(a50, F12.2, F12.2)') "Time tracers step (max and min): ", maxtime_model_step, mintime_model_step
        
        call mpi_allreduce(mpp_time_sync, maxtime_sync, 1, mpi_real8, mpi_max, mpp_cart_comm, ierr)
        call mpi_allreduce(mpp_time_sync, mintime_sync, 1, mpi_real8, mpi_min, mpp_cart_comm, ierr)
        if (mpp_rank == 0) write(*,'(a50, F12.2, F12.2)') "Time sync (max and min): ", maxtime_sync, mintime_sync
        call mpi_allreduce(mpp_time_sync_inner, maxtime_sync, 1, mpi_real8, mpi_max, mpp_cart_comm, ierr)
        call mpi_allreduce(mpp_time_sync_inner, mintime_sync, 1, mpi_real8, mpi_min, mpp_cart_comm, ierr)
        if (mpp_rank == 0) write(*,'(a50, F12.2, F12.2)') "Time sync inner (max and min): ", maxtime_sync, mintime_sync
        call mpi_allreduce(mpp_time_sync_boundary, maxtime_sync, 1, mpi_real8, mpi_max, mpp_cart_comm, ierr)
        call mpi_allreduce(mpp_time_sync_boundary, mintime_sync, 1, mpi_real8, mpi_min, mpp_cart_comm, ierr)
        if (mpp_rank == 0) write(*,'(a50, F12.2, F12.2)') "Time sync boundary (max and min): ", maxtime_sync, mintime_sync
        call mpi_allreduce(mpp_time_sync_intermediate, maxtime_sync, 1, mpi_real8, mpi_max, mpp_cart_comm, ierr)
        call mpi_allreduce(mpp_time_sync_intermediate, mintime_sync, 1, mpi_real8, mpi_min, mpp_cart_comm, ierr)
        if (mpp_rank == 0) write(*,'(a50, F12.2, F12.2)') "Time sync intermediate (max and min): ", maxtime_sync, mintime_sync

        
        call mpi_allreduce(mpp_time_sync_pack_mpi, maxtime_sync, 1, mpi_real8, mpi_max, mpp_cart_comm, ierr)
        call mpi_allreduce(mpp_time_sync_pack_mpi, mintime_sync, 1, mpi_real8, mpi_min, mpp_cart_comm, ierr)
        if (mpp_rank == 0) write(*,'(a50, F12.2, F12.2)') "Time sync pack mpi (max and min): ", maxtime_sync, mintime_sync
        
        call mpi_allreduce(mpp_time_sync_unpack_mpi, maxtime_sync, 1, mpi_real8, mpi_max, mpp_cart_comm, ierr)
        call mpi_allreduce(mpp_time_sync_unpack_mpi, mintime_sync, 1, mpi_real8, mpi_min, mpp_cart_comm, ierr)
        if (mpp_rank == 0) write(*,'(a50, F12.2, F12.2)') "Time sync unpack mpi (max and min): ", maxtime_sync, mintime_sync
        
        call mpi_allreduce(mpp_time_sync_isend_irecv, maxtime_sync, 1, mpi_real8, mpi_max, mpp_cart_comm, ierr)
        call mpi_allreduce(mpp_time_sync_isend_irecv, mintime_sync, 1, mpi_real8, mpi_min, mpp_cart_comm, ierr)
        if (mpp_rank == 0) write(*,'(a50, F12.2, F12.2)') "Time sync isend and irecv (max and min): ", maxtime_sync, mintime_sync
        
        call mpi_allreduce(mpp_time_sync_wait, maxtime_sync, 1, mpi_real8, mpi_max, mpp_cart_comm, ierr)
        call mpi_allreduce(mpp_time_sync_wait, mintime_sync, 1, mpi_real8, mpi_min, mpp_cart_comm, ierr)
        if (mpp_rank == 0) write(*,'(a50, F12.2, F12.2)') "Time sync wait (max and min): ", maxtime_sync, mintime_sync

        if (mpp_rank == 0) write(*,'(a50, F12.2)') "Time load balance: ", mpp_time_load_balance

        if (mpp_rank == 0) write(*,'(a50, F12.2)') "Time device init data: ", mpp_device_init_data
        if (mpp_rank == 0) write(*,'(a50, F12.2)') "Time device model step: ", mpp_device_time_model_step

        call mpi_allreduce(mpp_max_count_sync_send_recv, maxcount, 1, mpi_integer, mpi_max, mpp_cart_comm, ierr)
        call mpi_allreduce(mpp_max_count_sync_send_recv, mincount, 1, mpi_integer, mpi_min, mpp_cart_comm, ierr)
        if (mpp_rank == 0) write(*,'(a50, I12, I12)') "Max Count of simul sync send, recv (max and min): ", maxcount, mincount

        call mpi_allreduce(mpp_min_count_sync_send_recv, maxcount, 1, mpi_integer, mpi_max, mpp_cart_comm, ierr)
        call mpi_allreduce(mpp_min_count_sync_send_recv, mincount, 1, mpi_integer, mpi_min, mpp_cart_comm, ierr)
        if (mpp_rank == 0) write(*,'(a50, I12, I12)') "Min Count of simul sync send, recv (max and min): ", maxcount, mincount


#ifdef _MPP_KERNEL_TIMER_ON_
        do k = 1, max_kernels
            call mpi_allreduce(mpp_time_kernels(k), maxtime_kernel, 1, mpi_real8, mpi_max, mpp_cart_comm, ierr)
            call mpi_allreduce(mpp_time_kernels(k), mintime_kernel, 1, mpi_real8, mpi_min, mpp_cart_comm, ierr)
            if (maxtime_kernel > 0.0d0) then
                call get_kernel_name(k, kernel_name)
                if (mpp_is_master()) then
                    write(*,'(a25, a25, I10, F12.2, F12.2)') kernel_name, "Time (calls, max and min): ", mpp_calls_kernels(k), maxtime_kernel, mintime_kernel
                endif
            endif
        enddo
        call mpi_barrier(mpp_cart_comm, ierr)

        if (mpp_is_master()) then
            write(*,'(a50)') 'Times for each thread:'
        endif

        ! Threads times
        do k = 1, max_kernels
            do t = 0, _OMP_MAX_THREADS_ - 1
                call mpi_allreduce(mpp_time_kernels_threads(t, k), maxtime_kernel_threads(t), 1, mpi_real8, mpi_max, mpp_cart_comm, ierr)
                call mpi_allreduce(mpp_time_kernels_threads(t, k), mintime_kernel_threads(t), 1, mpi_real8, mpi_min, mpp_cart_comm, ierr)
            enddo

            ! Max/min for thread
            maxtime_kernel = 0.0d0
            mintime_kernel = huge(0.0d0)
            do t = 0, _OMP_MAX_THREADS_ - 1
                if (maxtime_kernel_threads(t) > maxtime_kernel) maxtime_kernel = maxtime_kernel_threads(t)
                if (mintime_kernel_threads(t) > 0 .and. mintime_kernel_threads(t) < mintime_kernel) mintime_kernel = mintime_kernel_threads(t)
            enddo

            ! Print
            if (maxtime_kernel > 0.0d0) then
                call get_kernel_name(k, kernel_name)
                if (mpp_is_master()) then
                    write(*,'(a25, a50, F12.2, F12.2)') kernel_name, "(max and min for all threads): ", maxtime_kernel, mintime_kernel
                endif
            endif
        enddo
        call mpi_barrier(mpp_cart_comm, ierr)
        call deallocate(mpp_time_kernels_threads)
#endif

#ifdef _GPU_MODE_
        istat = cudaEventDestroy(device_start_event)
        istat = cudaEventDestroy(device_stop_event)
#endif

        call mpi_finalize(ierr)
    end subroutine

    subroutine start_timer(time)
        real(wp8), intent(inout) :: time
        time = mpi_wtime()
        return
    end subroutine

    subroutine end_timer(time)
        real(wp8), intent(inout) :: time
        time = mpi_wtime() - time
        return
    end subroutine

    subroutine start_device_timer()
        integer :: istat
#ifdef _GPU_MODE_
        istat = cudaEventRecord(device_start_event, 0)
#endif
    end subroutine

    subroutine end_device_timer(time)
        real(wp4), intent(inout) :: time
        integer :: istat
        time = 0.0
#ifdef _GPU_MODE_
        istat = cudaEventRecord(device_stop_event, 0)
        istat = cudaEventSynchronize(device_stop_event)
        istat = cudaEventElapsedTime(time, device_start_event, device_stop_event)
        time = time * 0.001 ! ms to sec
#endif
    end subroutine

    function mpp_get_thread() result(t)
        integer :: t
        t = omp_get_thread_num()
        return
    end function

    function mpp_is_master() result(is)
        logical :: is

        if (mpp_rank == 0 .and. omp_get_thread_num() == 0) then
            is = .true.
        else
            is = .false.
        endif

        return
    end function

    function mpp_is_master_thread() result(is)
        logical :: is

        if (omp_get_thread_num() == 0) then
            is = .true.
        else
            is = .false.
        endif

        return
    end function

    function mpp_is_master_process() result(is)
        logical :: is

        if (mpp_rank == 0) then
            is = .true.
        else
            is = .false.
        endif

        return
    end function

    subroutine mpp_sync_output()
        if (mpp_is_master()) then
            call flush(6)
            print *, " "
            call flush(6)
        endif
    end subroutine



end module mpp_module
