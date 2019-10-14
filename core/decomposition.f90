module decomposition_module
    ! Decomposition of computational area

    use basinpar_module, only: nx, ny
    use mpp_module, only: mpp_rank, mpp_count, mpp_cart_comm, mpp_size, mpp_coord
    use errors_module, only: abort_model
    
    implicit none
    save
    private

    type, public :: domain_type
        integer :: bcount
        integer, pointer :: bnx_start(:), bnx_end(:)
        integer, pointer :: bny_start(:), bny_end(:)
    contains
        procedure, public :: init
        procedure, public :: clear
    end type domain_type

    type(domain_type), public, target :: domain_data
    
contains

    subroutine init(this, bppnx, bppny)
        ! Initialization of each domain
        class(domain_type), intent(inout) :: this
        integer, intent(in) :: bppnx, bppny
        integer :: k, i, xloc, yloc, xs, ys, ierr

        this%bcount = bppnx*bppny
        allocate(this%bnx_start(this%bcount), this%bnx_end(this%bcount))
        allocate(this%bny_start(this%bcount), this%bny_end(this%bcount))

        if (mod(nx, (bppnx * mpp_size(1))) /= 0) then
            call abort_model('Can not decompose direct X axis')
        endif
        if (mod(ny, (bppny * mpp_size(2))) /= 0) then
            call abort_model('Can not decompose direct Y axis')
        endif

        xloc = nx / (bppnx * mpp_size(1))
        xs = 1
        yloc = ny / (bppny * mpp_size(2))
        ys = 1
        do k = 1, this%bcount
            this%bnx_start(k) = xs
            this%bnx_end(k) = xs + xloc - 1
            xs = xs + xloc
            
            this%bny_start(k) = ys
            this%bny_end(k) = ys + yloc - 1
            ys = ys + yloc
        enddo

        if (mpp_rank .eq. 0) print *, "MPI pocs: ", mpp_count, " Domain decomposition:"
        do i = 0, mpp_count - 1
            if (mpp_rank .eq. i) then
                print *, "rank, coord", mpp_rank, mpp_coord
                do k = 1, this%bcount
                    print *, "block, bnx bounds", k, this%bnx_start(k), this%bnx_end(k)
                    print *, "block, bny bounds", k, this%bny_start(k), this%bny_end(k)
                enddo
            endif
            call mpi_barrier(mpp_cart_comm, ierr)
        enddo
    end subroutine

    subroutine clear(this)
        ! Initialization of each domain
        class(domain_type), intent(inout) :: this
        
        deallocate(this%bnx_start, this%bnx_end)
        deallocate(this%bny_start, this%bny_end)
    end subroutine

end module decomposition_module