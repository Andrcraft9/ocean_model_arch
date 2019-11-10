module decomposition_module
    ! Decomposition of computational area

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use mpp_module, only: mpp_rank, mpp_count, mpp_cart_comm, mpp_size, mpp_coord
    use errors_module, only: abort_model, check_error
    
    implicit none
    save
    private

    ! Domain decomposition only in 2D area, s-levels are the same for each process, decomposition into rectangular blocks
    ! There are two types of block numeration:
    !   Local:
    !       for each proc k = 1, ..., bcount. In custom order!
    !   Global:
    !       Cart coords of blocks (1, 1), ..., (bnx, bny)
    type, public :: domain_type
        ! Local area
        integer :: bcount ! Number of data blocks at each process
        integer, pointer :: bnx_start(:), bnx_end(:) ! Significant point area in blocks (x direction)
        integer, pointer :: bny_start(:), bny_end(:) ! Significant point area in blocks (y direction)
        integer, pointer :: bbnd_x1(:), bbnd_x2(:) ! Array boundary in blocks (x direction)
        integer, pointer :: bbnd_y1(:), bbnd_y2(:) ! Array boundary in blocks (y direction)

        ! Global area
        integer :: bnx, bny ! Cart grid of blocks
        integer :: total_blocks ! Total active blocks
        integer :: bcount_max, bcount_min ! Min/Max of local blocks
        ! Map: block coords to proc
        ! If bglob_proc(m, n) == -1 => (m, n) block is land-block!
        integer, pointer :: bglob_proc(:, :)
        ! Map: local block number to block coords
        integer, pointer :: bindx(:, :)

        integer :: nx, ny ! Number of total points in x,y directions
        integer :: nz     ! Number of s-levels in vertical direction
        integer :: periodicity_x, periodicity_y ! Periodicity on x,y direction (0 - non-periodic, 1 - periodic)
        integer :: mmm, mm ! begin and end of significant area in x direction
        integer :: nnn, nn ! begin and end of significant area in y direction
        ! Note: full significant area for z direction (from 1 to nz) 
    contains
        procedure, public :: init
        procedure, public :: clear
        
        procedure, private :: read_config
        procedure, private :: block_uniform_decomposition
        procedure, private :: create_uniform_decomposition
    end type domain_type

    type(domain_type), public, target :: domain_data
    
contains

    subroutine read_config(this)
        
        ! Здесь нужно вместо того, чтобы брать значения из модуля, читать конфиг и брать значения оттуда.
        ! Это пока заглушка.
        use basinpar_module

        class(domain_type), intent(inout) :: this

        this%nx = nx
        this%ny = ny
        this%nz = nz
        
        this%periodicity_x = periodicity_x
        this%periodicity_y = periodicity_y
        
        this%mmm = mmm
        this%mm = mm
        this%nnn = nnn
        this%nn = nn
    end subroutine

    subroutine block_uniform_decomposition(this, lbasins,  &
                                           glob_bnx_start, glob_bnx_end, glob_bny_start, glob_bny_end,  & 
                                           glob_bbnd_x1, glob_bbnd_x2, glob_bbnd_y1, glob_bbnd_y2,  &
                                           bglob_weight, land_blocks)
        class(domain_type), intent(in) :: this
        integer, pointer, intent(in) :: lbasins(:,:) 
        integer, allocatable, intent(out) :: glob_bnx_start(:, :), glob_bnx_end(:, :),  &
                                             glob_bny_start(:, :), glob_bny_end(:, :)
        integer, allocatable, intent(out) :: glob_bbnd_x1(:, :), glob_bbnd_x2(:, :),  &
                                             glob_bbnd_y1(:, :), glob_bbnd_y2(:, :)
        real(wp8), allocatable, intent(out) :: bglob_weight(:, :)
        integer, intent(out) :: land_blocks
        integer :: m, n, locn
        
        associate(bnx => this%bnx,  &
                  bny => this%bny,  &
                  nx => this%nx,  &
                  ny => this%ny)
            
            land_blocks = 0
            bglob_weight = 0.0

            do m = 1, bnx
                do n = 1, bny
                    locn = floor(real(nx - 4)/real(bnx))
                    glob_bnx_start(m, n) = locn*(m-1) + 1 + 2
                    if (m .eq. bnx) then
                        locn = (nx - 2) - glob_bnx_start(m, n) + 1
                    endif
                    glob_bnx_end(m, n) = glob_bnx_start(m, n) + locn - 1
                    glob_bnx_start(m, n) = glob_bnx_start(m, n)
                    ! border area
                    glob_bbnd_x1(m, n) = glob_bnx_start(m, n) - 2
                    glob_bbnd_x2(m, n) = glob_bnx_end(m, n) + 2

                    locn = floor(real(ny - 4)/real(bny))
                    glob_bny_start(m, n) = locn*(n-1) + 1 + 2
                    if (n .eq. bny) then
                        locn = (ny - 2) - glob_bny_start(m, n) + 1
                    endif
                    glob_bny_end(m, n) = glob_bny_start(m, n) + locn - 1
                    glob_bny_start(m, n) = glob_bny_start(m, n)
                    ! border area
                    glob_bbnd_y1(m, n) = glob_bny_start(m, n) - 2
                    glob_bbnd_y2(m, n) = glob_bny_end(m, n) + 2

                    ! Compute load-balance
                    do i = glob_bnx_start(m, n), glob_bnx_end(m, n)
                        do j = glob_bny_start(m, n), glob_bny_end(m, n)
                            bglob_weight(m, n) = bglob_weight(m, n) + (1.0d0 - real(lbasins(i, j)))
                        enddo
                    enddo
                    ! Only-land blocks
                    if (bglob_weight(m, n) == 0.0d0) then
                        land_blocks = land_blocks + 1
                    endif
                enddo
            enddo
            if (mpp_rank == 0) print *, "Total land blocks:", land_blocks

            ierr = 0
            if (bnx*bny - land_blocks < procs) ierr = 1
            call check_error(ierr,  'procs > computational-blocks... Error!')
        end associate
    end subroutine

    subroutine create_uniform_decomposition(this, bglob_weight)
        class(domain_type), intent(inout) :: this
        real(wp8), allocatable, intent(in) :: bglob_weight(:, :)

        associate(bnx => this%bnx,  &
                  bny => this%bny)

        end associate
    end subroutine

    subroutine init(this, bppnx, bppny, lbasins)
        ! Initialization of each domain
        class(domain_type), intent(inout) :: this
        integer, intent(in) :: bppnx, bppny
        integer, pointer, intent(in) :: lbasins(:,:) 

        real(wp8), allocatable :: bglob_weight(:, :)
        real(wp8) :: bweight, max_bweight

        integer, allocatable :: glob_bnx_start(:, :), glob_bnx_end(:, :),  &
                                glob_bny_start(:, :), glob_bny_end(:, :)
        integer, allocatable :: glob_bbnd_x1(:, :), glob_bbnd_x2(:, :),  &
                                glob_bbnd_y1(:, :), glob_bbnd_y2(:, :)

        integer :: m, n, i, j, k, locn
        integer :: land_blocks
        integer :: bshared
        real(wp8) :: bcomm_metric, max_bcomm_metric
        integer, dimension(2) :: kblock
        integer :: ierr

        integer :: itot, ishared
        real(wp8) :: icomm_metric, max_icomm_metric

        associate(bcount => this%bcount,  &
                  bnx_start => this%bnx_start,  &
                  bnx_end => this%bnx_end,  &
                  bny_start => this%bny_start,  &
                  bny_end => this%bny_end,  &
                  bbnd_x1 => this%bbnd_x1,  &
                  bbnd_x2 => this%bbnd_x2,  &
                  bbnd_y1 => this%bbnd_y1,  &
                  bbnd_y2 => this%bbnd_y2,  &
                  bnx => this%bnx,  & 
                  bny => this%bny,  &
                  total_blocks => this%total_blocks,  &
                  bcount_max => this%bcount_max,  &
                  bcount_min => this%bcount_min,  &
                  bglob_proc => this%bglob_proc,  &
                  bindx => this%bindx,  &
                  nx => this%nx,  &
                  ny => this%ny)

            ! Set Cart grid of blocks
            bnx = bppnx * mpp_size(1)
            bny = bppny * mpp_size(2)

            if (mpp_rank == 0) print *, 'bnx, bny and Total blocks:', bnx, bny, bnx * bny
            if (mpp_rank == 0) print *, 'pnx, pny and procs:', mpp_size(1), mpp_size(2), mpp_count
            call mpi_barrier(mpp_cart_comm, ierr)

            ! Check bnx and bny
            ierr = 0
            if (mod(bnx, mpp_size(1)) /= 0) ierr = 1
            if (mod(bny, mpp_size(2)) /= 0) ierr = 1
            call check_error(ierr, 'Error in bny or bnx! mod(bnx, p_size(1)) or mod(bny, p_size(2)) not equal 0 !')

            allocate(bglob_weight(bnx, bny))
            allocate(glob_bnx_start(bnx, bny), glob_bnx_end(bnx, bny), glob_bny_start(bnx, bny), glob_bny_end(bnx, bny))
            allocate(glob_bbnd_x1(bnx, bny), glob_bbnd_x2(bnx, bny), glob_bbnd_y1(bnx, bny), glob_bbnd_y2(bnx, bny))
            bglob_weight = 0.0
            glob_bnx_start = 0; glob_bnx_end = 0; glob_bny_start = 0; glob_bny_end = 0
            glob_bbnd_x1 = 0; glob_bbnd_x2 = 0; glob_bbnd_y1 = 0; glob_bbnd_y2 = 0

            call this%block_uniform_decomposition(lbasins,  &
                                                  glob_bnx_start, glob_bnx_end, glob_bny_start, glob_bny_end,  & 
                                                  glob_bbnd_x1, glob_bbnd_x2, glob_bbnd_y1, glob_bbnd_y2,  &
                                                  bglob_weight, land_blocks)

            ! Compute bglob_proc
            allocate(bglob_proc(bnx, bny))
            if (mpp_rank == 0) print *, "Uniform blocks decomposition!..."
            call this%create_uniform_decomposition(bglob_weight)

            ! Compute blocks per proc
            bcount = 0; bweight = 0.0d0
            do m = 1, bnx
                do n = 1, bny
                    if (bglob_proc(m, n) == mpp_rank) then
                        bcount = bcount + 1
                        bweight = bweight + bglob_weight(m, n)
                    endif
                enddo
            enddo
            call mpi_allreduce(bcount, total_blocks, 1, mpi_integer, mpi_sum, mpp_cart_comm, ierr)
            call mpi_allreduce(bweight, max_bweight, 1, mpi_real8, mpi_max, mpp_cart_comm, ierr)
            call mpi_allreduce(bcount, bcount_max, 1, mpi_integer, mpi_max, mpp_cart_comm, ierr)
            call mpi_allreduce(bcount, bcount_min, 1, mpi_integer, mpi_min, mpp_cart_comm, ierr)

            ierr = 0
            if (bcount <= 0) ierr = 1
            call check_error(ierr, 'Proc with only land-blocks... Error!')

            ! Print information about blocks
            if (mpp_rank == 0) print *, 'Total blocks:', total_blocks, 'LB: ', max_bweight / (sum(bglob_weight) / real(mpp_count)),  &
                                        'max blocks per proc:', bcount_max, 'min blocks per proc:', bcount_min
            call mpi_barrier(mpp_cart_comm, ierr)
            
            ! Allocate blocks arrays per proc
            allocate(bindx(bcount, 2))
            allocate(bnx_start(bcount), bnx_end(bcount),  &
                     bny_start(bcount), bny_end(bcount))
            allocate(bbnd_x1(bcount), bbnd_x2(bcount),  &
                     bbnd_y1(bcount), bbnd_y2(bcount))
            bindx = 0
            bnx_start = 0; bnx_end = 0 
            bny_start = 0; bny_end = 0
            bbnd_x1 = 0; bbnd_x2 = 0 
            bbnd_y1 = 0; bbnd_y2 = 0

            k = 1
            do m = 1, bnx
                do n = 1, bny
                    if (bglob_proc(m, n) == mpp_rank) then
                        ! Map local block numeration to block coords
                        bindx(k, 1) = m
                        bindx(k, 2) = n

                        bnx_start(k) = glob_bnx_start(m, n)
                        bnx_end(k) = glob_bnx_end(m, n)
                        bny_start(k) = glob_bny_start(m, n)
                        bny_end(k) = glob_bny_end(m, n)

                        bbnd_x1(k) = glob_bbnd_x1(m, n)
                        bbnd_x2(k) = glob_bbnd_x2(m, n)
                        bbnd_y1(k) = glob_bbnd_y1(m, n)
                        bbnd_y2(k) = glob_bbnd_y2(m, n)

                        k = k + 1
                    endif
                enddo
            enddo

            deallocate(bglob_weight)
            deallocate(glob_bnx_start, glob_bnx_end, glob_bny_start, glob_bny_end)
            deallocate(glob_bbnd_x1, glob_bbnd_x2, glob_bbnd_y1, glob_bbnd_y2)

            ! DEBUG
            if (mpp_rank == 0) then
                do k = 1, bcount
                    write(*, '(i5, i5,i5,i5,i5)') k, bnx_start(k), bnx_end(k), bny_start(k), bny_end(k)
                enddo
            endif
        end associate
    end subroutine

    subroutine clear(this)
        ! Initialization of each domain
        class(domain_type), intent(inout) :: this
        
        deallocate(this%bnx_start, this%bnx_end)
        deallocate(this%bny_start, this%bny_end)
    end subroutine

end module decomposition_module