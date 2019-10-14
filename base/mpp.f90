module mpp_module
    ! MPI massively parallel processing library

    use mpi

    implicit none
    save
    private

    integer, public :: mpp_rank, mpp_count
    integer, public :: mpp_cart_comm
    integer, dimension(2), public :: mpp_size, mpp_coord
    logical, dimension(2), public :: mpp_period

    public :: mpp_init
    public :: mpp_finalize

contains

    subroutine mpp_init()
        ! Subroutine description
        integer :: ierr

        call mpi_init(ierr)

        mpp_period = (/.true., .true./)
        mpp_size = (/0,0/)
        ierr = 0

        call mpi_comm_rank(mpi_comm_world, mpp_rank, ierr)
        call mpi_comm_size(mpi_comm_world, mpp_count, ierr)
        call mpi_dims_create(mpp_count, 2, mpp_size, ierr)
        call mpi_cart_create(mpi_comm_world, 2, mpp_size, mpp_period, .false., mpp_cart_comm, ierr)
        call mpi_cart_coords(mpp_cart_comm, mpp_rank, 2, mpp_coord, ierr)
    end subroutine

    subroutine mpp_finalize()
        integer :: ierr

        call mpi_finalize(ierr)
    end subroutine

end module mpp_module
