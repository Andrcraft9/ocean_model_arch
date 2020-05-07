module decomposition_module
    ! Decomposition of computational area

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use mpi
    use debug_module
    use mpp_module
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

        integer :: nz
    contains
        procedure, public :: init_from_config
        procedure, public :: init
        procedure, public :: clear
        
        procedure, private :: block_uniform_decomposition
        procedure, private :: create_uniform_decomposition
        procedure, private :: create_hilbert_curve_decomposition
    end type domain_type

!------------------------------------------------------------------------------

    type(domain_type), public, target :: domain_data
    
contains

    subroutine block_uniform_decomposition(this, lbasins,  &
                                           glob_bnx_start, glob_bnx_end, glob_bny_start, glob_bny_end,  & 
                                           glob_bbnd_x1, glob_bbnd_x2, glob_bbnd_y1, glob_bbnd_y2,  &
                                           bglob_weight, land_blocks)
        use config_basinpar_module, only: nx, ny
        
        class(domain_type), intent(in) :: this
        integer, pointer, intent(in) :: lbasins(:,:) 
        integer, allocatable, intent(inout) :: glob_bnx_start(:, :), glob_bnx_end(:, :),  &
                                             glob_bny_start(:, :), glob_bny_end(:, :)
        integer, allocatable, intent(inout) :: glob_bbnd_x1(:, :), glob_bbnd_x2(:, :),  &
                                             glob_bbnd_y1(:, :), glob_bbnd_y2(:, :)
        real(wp8), allocatable, intent(inout) :: bglob_weight(:, :)
        integer, intent(out) :: land_blocks
        integer :: m, n, i, j, locn, ierr
        
        associate(bnx => this%bnx,  &
                  bny => this%bny)
            
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
            if (mpp_is_master()) print *, "Total land blocks:", land_blocks

            ierr = 0
            if (bnx*bny - land_blocks < mpp_count) ierr = 1
            call check_error(ierr,  'procs > computational-blocks... Error!')
        end associate
    end subroutine

    subroutine create_hilbert_curve_decomposition(this, bglob_weight, land_blocks)
        
        use hilbert_curve_module

        class(domain_type), intent(in) :: this
        real(wp8), intent(in) :: bglob_weight(:, :)
        integer, intent(in) :: land_blocks

        integer :: k, i, ierr
        integer :: hilbert_index, hilbert_coord_x, hilbert_coord_y
        integer :: ks, sea_blocks
        real(wp8) :: weight, tot_weight, mean_weight, last_weight

        hilbert_index = int(log(real(this%bnx))/ log(2.0))
        ierr = 0
        if (this%bnx /= this%bny) then
            if (mpp_is_master()) print *, 'bnx not equal to bny! Can`t build Hilbert curve for this geometry!'
            ierr = 1
        endif
        if (2**hilbert_index /= this%bnx) then
            if (mpp_is_master()) print *, '2**M not eqal to bnx! Can`t build Hilbert curve for this geometry!'
            ierr = 1
        endif
        call check_error(ierr, 'Can`t build Hilbert curve for this geometry!')

        if (mpp_is_master()) print *, 'Hilbert curve index:', hilbert_index

        tot_weight = sum(bglob_weight)
        mean_weight = tot_weight / mpp_count
        sea_blocks = this%bnx * this%bny - land_blocks

        if (debug_level >= 1) then
            if (mpp_is_master()) print *, 'Total blocks weigth:', tot_weight, "Mean blocks weigth:", mean_weight
        endif
        call mpi_barrier(mpp_cart_comm, ierr)

        weight = 0.0d0
        last_weight = 0.0d0
        i = 0; ks = 0
        do k = 1, this%bnx * this%bny
            call hilbert_d2xy(hilbert_index, k-1, hilbert_coord_x, hilbert_coord_y)
            hilbert_coord_x = hilbert_coord_x + 1; hilbert_coord_y = hilbert_coord_y + 1

            ! Skip land blocks
            if (bglob_weight(hilbert_coord_x, hilbert_coord_y) == 0.0d0) then
                this%bglob_proc(hilbert_coord_x, hilbert_coord_y) = -1
                cycle
            else
                ks = ks + 1
            endif

            weight = weight + bglob_weight(hilbert_coord_x, hilbert_coord_y)

            !if (sea_blocks - ks >= mpp_count - i - 1) then
                if (weight + (weight - bglob_weight(hilbert_coord_x, hilbert_coord_y)) > 2.0*mean_weight) then
                    ! Recompute mean value
                    mean_weight = (tot_weight - last_weight) / (mpp_count - i - 1)
                    ! Go to next proc
                    i = i + 1
                    weight = bglob_weight(hilbert_coord_x, hilbert_coord_y)
                    if (i > mpp_count - 1) then
                        i = mpp_count - 1
                        !if (rank == 0 .and. parallel_dbg < 2) print *, 'Warning! Last procs ...'
                        if (mpp_is_master()) print *, k, 'Warning! Last procs ...'
                    endif
                endif
            !else
            !    i = i + 1
            !    weight = bglob_weight(hilbert_coord_x, hilbert_coord_y)
            !endif

            this%bglob_proc(hilbert_coord_x, hilbert_coord_y) = i
            last_weight = last_weight + bglob_weight(hilbert_coord_x, hilbert_coord_y)
        enddo

        if (debug_level >= 3) then
            call parallel_int_output(this%bglob_proc, 1, this%bnx, 1, this%bny, 'bglob_proc from load-balanced hilbert curve decomposition')
        endif
    end subroutine

    subroutine create_uniform_decomposition(this, bglob_weight)
        class(domain_type), intent(inout) :: this
        real(wp8), allocatable, intent(in) :: bglob_weight(:, :)
        integer :: m, n, ierr
        integer :: loc_bnx, loc_bny
        integer :: xblock_start, yblock_start
        integer, allocatable :: buf_int(:, :)

        associate(bnx => this%bnx,  &
                  bny => this%bny,  &
                  bgproc => this%bglob_proc)
    
            loc_bnx = bnx / mpp_size(1)
            loc_bny = bny / mpp_size(2)
            !bcount = loc_bnx*loc_bny
            !if (parallel_dbg >= 1) print *, rank, 'loc_bnx, loc_bny and Blocks per proc: ', loc_bnx, loc_bny, bcount
            !call mpi_barrier(cart_comm, ierr)
    
            xblock_start = 1 + mpp_coord(1)*loc_bnx
            yblock_start = 1 + mpp_coord(2)*loc_bny
            if (debug_level >= 2) then 
                print *, mpp_rank, 'xb_start, yb_start', xblock_start, yblock_start
                call mpi_barrier(mpp_cart_comm, ierr)
            endif
            bgproc = -1
            do m = xblock_start, xblock_start + loc_bnx - 1
                do n = yblock_start, yblock_start + loc_bny - 1
                    bgproc(m, n) = mpp_rank
                    if (bglob_weight(m, n) == 0.0) then
                        bgproc(m, n) = -1
                    endif
                enddo
            enddo
            ! Sync bglob_proc array
            allocate(buf_int(bnx, bny))
            buf_int = bgproc + 1
            call mpi_allreduce(buf_int, bgproc, bnx*bny, mpi_integer, mpi_sum, mpp_cart_comm, ierr)
            bgproc = bgproc - 1
    
            if (debug_level >= 3) then
                call parallel_int_output(bgproc, 1, bnx, 1, bny, 'bglob_proc from uniform decomposition')
            endif
    
            deallocate(buf_int)

        end associate
    end subroutine

    subroutine init(this, bppnx, bppny, mod_create, lbasins)
        use config_basinpar_module, only: nx, ny, nz

        ! Initialization of each domain
        class(domain_type), intent(inout) :: this
        integer, intent(in) :: bppnx, bppny
        integer, intent(in) :: mod_create
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
        integer :: max_points_per_block, min_points_per_block, reduced_max_points_per_block, reduced_min_points_per_block

        this%nz = nz

        associate(bcount => this%bcount,  &
                  bnx_start => this%bnx_start,  &
                  bnx_end => this%bnx_end,  &
                  bny_start => this%bny_start,  &
                  bny_end => this%bny_end,  &
                  bnx => this%bnx,  & 
                  bny => this%bny,  &
                  total_blocks => this%total_blocks,  &
                  bcount_max => this%bcount_max,  &
                  bcount_min => this%bcount_min)

            ! Set Cart grid of blocks
            bnx = bppnx * mpp_size(1)
            bny = bppny * mpp_size(2)

            if (mpp_is_master()) print *, 'bnx, bny and Total blocks:', bnx, bny, bnx * bny
            if (mpp_is_master()) print *, 'pnx, pny and procs:', mpp_size(1), mpp_size(2), mpp_count
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
            allocate(this%bglob_proc(bnx, bny))
            
            if (mod_create == 0) then
                if (mpp_is_master()) print *, "Uniform blocks decomposition!..."
                call this%create_uniform_decomposition(bglob_weight)
            elseif (mod_create == 1) then
                if (mpp_is_master()) print *, "Hilber Curve blocks decomposition!..."
                call this%create_hilbert_curve_decomposition(bglob_weight, land_blocks)
            else
                if (mpp_is_master()) print *, "Unknown mode!"
                call abort_model("Unknown decomposition mode!")
            endif

            ! Compute blocks per proc
            bcount = 0; bweight = 0.0d0
            do m = 1, bnx
                do n = 1, bny
                    if (this%bglob_proc(m, n) == mpp_rank) then
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
            if (mpp_is_master()) print *, 'Total blocks:', total_blocks, 'LB: ', max_bweight / (sum(bglob_weight) / real(mpp_count)),  &
                                        'max blocks per proc:', bcount_max, 'min blocks per proc:', bcount_min
            call mpi_barrier(mpp_cart_comm, ierr)
            
            ! Allocate blocks arrays per proc
            allocate(this%bindx(bcount, 2))
            allocate(this%bnx_start(bcount), this%bnx_end(bcount),  &
                     this%bny_start(bcount), this%bny_end(bcount))
            allocate(this%bbnd_x1(bcount), this%bbnd_x2(bcount),  &
                     this%bbnd_y1(bcount), this%bbnd_y2(bcount))
            this%bindx = 0
            this%bnx_start = 0; this%bnx_end = 0 
            this%bny_start = 0; this%bny_end = 0
            this%bbnd_x1 = 0; this%bbnd_x2 = 0 
            this%bbnd_y1 = 0; this%bbnd_y2 = 0

            k = 1
            do m = 1, bnx
                do n = 1, bny
                    if (this%bglob_proc(m, n) == mpp_rank) then
                        ! Map local block numeration to block coords
                        this%bindx(k, 1) = m
                        this%bindx(k, 2) = n

                        this%bnx_start(k) = glob_bnx_start(m, n)
                        this%bnx_end(k) = glob_bnx_end(m, n)
                        this%bny_start(k) = glob_bny_start(m, n)
                        this%bny_end(k) = glob_bny_end(m, n)

                        this%bbnd_x1(k) = glob_bbnd_x1(m, n)
                        this%bbnd_x2(k) = glob_bbnd_x2(m, n)
                        this%bbnd_y1(k) = glob_bbnd_y1(m, n)
                        this%bbnd_y2(k) = glob_bbnd_y2(m, n)

                        k = k + 1
                    endif
                enddo
            enddo

            deallocate(bglob_weight)
            deallocate(glob_bnx_start, glob_bnx_end, glob_bny_start, glob_bny_end)
            deallocate(glob_bbnd_x1, glob_bbnd_x2, glob_bbnd_y1, glob_bbnd_y2)

            ! DEBUG
            if (debug_level >= 1) then
                max_points_per_block = 0
                min_points_per_block = nx*ny
                do k = 1, bcount
                    if ( (this%bnx_end(k) - this%bnx_start(k) + 1) * (this%bny_end(k) - this%bny_start(k) + 1) > max_points_per_block) then
                        max_points_per_block = (this%bnx_end(k) - this%bnx_start(k) + 1) * (this%bny_end(k) - this%bny_start(k) + 1)
                    endif
                    if ( (this%bnx_end(k) - this%bnx_start(k) + 1) * (this%bny_end(k) - this%bny_start(k) + 1) < min_points_per_block) then
                        min_points_per_block = (this%bnx_end(k) - this%bnx_start(k) + 1) * (this%bny_end(k) - this%bny_start(k) + 1)
                    endif
                enddo
                print *, mpp_rank, 'Blocks', bcount, 'Max points per block:', max_points_per_block, 'Min points per blocks: ', min_points_per_block
                call mpi_barrier(mpp_cart_comm, ierr)
                call mpi_allreduce(min_points_per_block, reduced_min_points_per_block, 1, mpi_integer, mpi_min, mpp_cart_comm, ierr)
                call mpi_allreduce(max_points_per_block, reduced_max_points_per_block, 1, mpi_integer, mpi_max, mpp_cart_comm, ierr)
                if (mpp_is_master()) print *, 'Reduced Max points per block:', reduced_max_points_per_block, 'Reduced Min points per blocks: ', reduced_min_points_per_block
                call mpi_barrier(mpp_cart_comm, ierr)
            endif

            if (debug_level >= 2) then
                if (mpp_is_master()) then
                    do k = 1, bcount
                        write(*, '(i5, i5,i5,i5,i5)') k, this%bnx_start(k), this%bnx_end(k), this%bny_start(k), this%bny_end(k)
                        write(*, '(i5, i5,i5,i5,i5)') k, this%bbnd_x1(k), this%bbnd_x2(k), this%bbnd_y1(k), this%bbnd_y2(k)
                    enddo
                endif
            endif
        end associate
    end subroutine

    subroutine clear(this)
        ! Initialization of each domain
        class(domain_type), intent(inout) :: this
        
        deallocate(this%bnx_start, this%bnx_end)
        deallocate(this%bny_start, this%bny_end)

        deallocate(this%bbnd_x1, this%bbnd_x2)
        deallocate(this%bbnd_y1, this%bbnd_y2)

        deallocate(this%bglob_proc)
        deallocate(this%bindx)
    end subroutine

    subroutine init_from_config(this, lbasins, name)
        use rwpar_routes

        class(domain_type), intent(inout) :: this
        integer, pointer, intent(in) :: lbasins(:,:) 
        character(*) :: name

        integer :: nofcom
        character(128) :: comments(128)
        integer :: ierr

        integer :: mod_decomposition
        character(128) :: file_decomposition
        integer :: bppnx, bppny
        integer :: parallel_dbg
        integer :: parallel_mod
        character(128) :: file_output

        if (mpp_is_master()) then
            print *, 'Read decomposition config...'
            call readpar(name, comments, nofcom)
            read(comments(1),*) mod_decomposition
            call get_first_lexeme(comments(2), file_decomposition)
            read(comments(3),*) bppnx
            read(comments(4),*) bppny
            read(comments(5),*) parallel_dbg
            read(comments(6),*) parallel_mod
            call get_first_lexeme(comments(7), file_output)

            print *, 'mod_decomposition=', mod_decomposition
            print *, '(ignore this) decomposition file:', file_decomposition
            print *, 'bppnx=', bppnx
            print *, 'bppny=', bppny
            print *, '(ignore this) parallel_dbg=', parallel_dbg
            print *, '(ignore this) parallel_mod=', parallel_mod
            print *, '(ignore this) output file:', file_output
        endif

        call mpi_bcast(file_decomposition, 128, mpi_character, 0, mpp_cart_comm, ierr)
        call mpi_bcast(mod_decomposition, 1, mpi_integer, 0, mpp_cart_comm, ierr)
        call mpi_bcast(bppnx, 1, mpi_integer, 0, mpp_cart_comm, ierr)
        call mpi_bcast(bppny, 1, mpi_integer, 0, mpp_cart_comm, ierr)
        call mpi_bcast(parallel_dbg, 1, mpi_integer, 0, mpp_cart_comm, ierr)
        call mpi_bcast(parallel_mod, 1, mpi_integer, 0, mpp_cart_comm, ierr)
        call mpi_bcast(file_output, 128, mpi_character, 0, mpp_cart_comm, ierr)

        call this%init(bppnx, bppny, mod_decomposition, lbasins)

    end subroutine

end module decomposition_module