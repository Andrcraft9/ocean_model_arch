module errors_module
    ! Error processing module

    use mpi
    use mpp_module

    implicit none
    save
    private

    public :: abort_model
    public :: check_error

contains

    subroutine check_error(loc_err, msg)
        integer, intent(inout) :: loc_err
        character(len=*), intent(in) :: msg

        integer :: tot_err
        integer :: ierr

        call mpi_allreduce(loc_err, tot_err, 1, mpi_integer, mpi_sum, mpp_cart_comm, ierr)
        if (tot_err >= 1) then
            call abort_model(msg)
        endif
        loc_err = tot_err
    end subroutine

    subroutine abort_model(msg)
        character(len=*), intent(in) :: msg
        integer :: ierr
        
        print *, msg
        call mpi_abort(mpp_cart_comm, -1, ierr)
        stop
    end subroutine

end module errors_module