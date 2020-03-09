module debug_module
    ! Configuration of debug in program, all utilits for debug

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use mpp_module, only: mpp_rank, mpp_count, mpp_cart_comm, mpp_size, mpp_coord, mpp_period
    use errors_module, only: abort_model, check_error

    implicit none
    save
    public

    integer, parameter :: debug_level = 5

contains

    subroutine parallel_int_output(arr, x1, x2, y1, y2, msg)
        implicit none
        integer :: x1, x2, y1, y2
        integer :: arr(x1:x2, y1:y2)
        character*(*) msg
        integer :: k, m, n, ierr

        if (mpp_rank .eq. 0) print *, msg
        call mpi_barrier(mpp_cart_comm, ierr)

        do k = 0, mpp_count-1
            if (mpp_rank .eq. k) then
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
