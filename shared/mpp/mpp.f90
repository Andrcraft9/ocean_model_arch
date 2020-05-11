module mpp_module
    ! MPI massively parallel processing library

    use mpi
    use kernel_runtime_module
    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4

    implicit none
    save
    public

#include "macros/mpp_macros.fi"

    !include 'mpif.h'
    include "omp_lib.h"

    integer, public :: mpp_rank, mpp_count
    integer, public :: mpp_cart_comm
    integer, dimension(2), public :: mpp_size, mpp_coord
    logical, dimension(2), public :: mpp_period

    integer, public :: mpp_count_threads

    public :: mpp_is_master, mpp_is_master_thread, mpp_is_master_process
    public :: start_timer, end_timer
    public :: mpp_init
    public :: mpp_finalize

    ! Timers for master thread
    real(wp8) :: mpp_time_model_step, mpp_time_sync
    real(wp8), allocatable :: mpp_time_kernels(:)
    integer, allocatable :: mpp_calls_kernels(:)

    ! Timers for threads
    real(wp8), allocatable :: mpp_time_kernels_threads(:, :)
contains

    subroutine mpp_init()
        ! Subroutine description
        integer :: ierr, req, provided, tmp
        integer :: rank_cart
        character(len=mpi_max_library_version_string) :: version

        mpp_period = (/.true., .true./)
        mpp_size = (/0,0/)
        req = MPI_THREAD_MULTIPLE
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
            print *, "MPP INFO: MPI version: ", trim(version)
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
        mpp_time_model_step = 0.0d0
        mpp_time_sync = 0.0d0

        allocate(mpp_time_kernels(max_kernels))
        allocate(mpp_calls_kernels(max_kernels))
        mpp_time_kernels = 0.0d0
        mpp_calls_kernels = 0

        allocate(mpp_time_kernels_threads(0 : _OMP_MAX_THREADS_ - 1, max_kernels))
        mpp_time_kernels_threads = 0
    end subroutine

    subroutine mpp_finalize()
        integer :: k, t, ierr
        real(wp8) :: maxtime_model_step, mintime_model_step
        real(wp8) :: maxtime_sync, mintime_sync
        real(wp8) :: maxtime_kernel, mintime_kernel
        real(wp8) :: maxtime_kernel_threads(0 : _OMP_MAX_THREADS_ - 1), mintime_kernel_threads(0 : _OMP_MAX_THREADS_ - 1)
        character(80) :: kernel_name

        if (mpp_is_master()) then
            write(*,'(a50)') 'Times for master thread:'
        endif

        ! Timers for master thread
        call mpi_allreduce(mpp_time_model_step, maxtime_model_step, 1, mpi_real8, mpi_max, mpp_cart_comm, ierr)
        call mpi_allreduce(mpp_time_model_step, mintime_model_step, 1, mpi_real8, mpi_min, mpp_cart_comm, ierr)
        if (mpp_rank == 0) write(*,'(a50, F12.2, F12.2)') "Time full of model step (max and min): ", maxtime_model_step, mintime_model_step
        call mpi_allreduce(mpp_time_sync, maxtime_sync, 1, mpi_real8, mpi_max, mpp_cart_comm, ierr)
        call mpi_allreduce(mpp_time_sync, mintime_sync, 1, mpi_real8, mpi_min, mpp_cart_comm, ierr)
        if (mpp_rank == 0) write(*,'(a50, F12.2, F12.2)') "Time sync (max and min): ", maxtime_sync, mintime_sync

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

end module mpp_module
