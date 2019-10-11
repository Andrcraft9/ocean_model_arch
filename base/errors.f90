module errors_module
    ! Module description
    !use kind_module, only: CL => SHR_KIND_CL
    use parallel_module

    implicit none
    save
    private

    public :: abort_model

contains

    subroutine abort_model(msg, procs)
        ! Subroutine description
        character(len=*), intent(in) :: msg
        type(procs_type), intent(in) :: procs
        integer :: ierr
        
        print *, msg
        call mpi_abort(procs%cart_comm, -1, ierr)
        stop
    end subroutine

end module errors_module