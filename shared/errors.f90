module errors_module
    ! Error processing module

    use mpp_module, only: mpp_rank, mpp_count, mpp_cart_comm

    implicit none
    save
    private

    public :: abort_model

contains

    subroutine abort_model(msg)
        ! Subroutine description
        character(len=*), intent(in) :: msg
        integer :: ierr
        
        print *, msg
        call mpi_abort(mpp_cart_comm, -1, ierr)
        stop
    end subroutine

end module errors_module