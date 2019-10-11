module parallel_module
    ! Module description

    use mpi

    implicit none
    save
    private

    type, public :: procs_type
        integer, public :: rank, count
        integer, public :: cart_comm
        integer, dimension(2), public :: size, coord
        logical, dimension(2), public :: period
    contains
        procedure, public :: init
        procedure, public :: finalize
    end type procs_type

    type(procs_type), public, target :: procs_data

contains

    subroutine init(this)
        ! Subroutine description
        class(procs_type), intent(inout) :: this
        integer :: ierr

        call mpi_init(ierr)

        this%period = (/.true., .true./)
        this%size = (/0,0/)
        ierr = 0

        call mpi_comm_rank(mpi_comm_world, this%rank, ierr)
        call mpi_comm_size(mpi_comm_world, this%count, ierr)
        call mpi_dims_create(this%count, 2, this%size, ierr)
        call mpi_cart_create(mpi_comm_world, 2, this%size, this%period, .false., this%cart_comm, ierr)
        call mpi_cart_coords(this%cart_comm, this%rank, 2, this%coord, ierr)
    end subroutine

    subroutine finalize(this)
        class(procs_type), intent(inout) :: this
        integer :: ierr

        call mpi_finalize(ierr)
    end subroutine

end module parallel_module
