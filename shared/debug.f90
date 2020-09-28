module debug_module
    ! Configuration of debug in program, all utilits for debug

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use mpp_module
    use errors_module, only: abort_model, check_error

    implicit none
    save
    public

    integer, parameter :: debug_level = 6
    ! 0: no debug info
    ! 1: mean debug info per ranks, print only master rank
    ! 2: 
    ! 3: create debug files
    ! 4:
    ! 5: debug info for each rank
    ! 6: debug info for each thread
    ! 7: extension debug
    ! 8:
    ! 9: debug info for each block
    !
    ! 10: sync, 1 level
    ! 11: sync, 2 level
    ! 12: sync, 3 level

contains

    subroutine parallel_int_output(arr, x1, x2, y1, y2, msg)
        implicit none
        integer :: x1, x2, y1, y2
        integer :: arr(x1:x2, y1:y2)
        character*(*) msg
        integer :: k, m, n, ierr

        if (mpp_rank == 0) print *, msg
        call mpi_barrier(mpp_cart_comm, ierr)

        do k = 0, mpp_count-1
            if (mpp_rank == k) then
                print *, "rank: ", mpp_rank
                do m = x1, x2
                    do n = y1, y2
                        print *, mpp_rank, 'm, n, a(m, n)', m, n, arr(m, n)
                    enddo
                enddo
            endif
            call mpi_barrier(mpp_cart_comm, ierr)
        enddo
    end subroutine

end module debug_module
