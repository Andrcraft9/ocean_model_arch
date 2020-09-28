module decomposition_module
    ! Decomposition of computational area

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    !use mpi
    use debug_module
    use mpp_module
    use errors_module, only: abort_model, check_error
    
    implicit none
    save
    private
    
#include "macros/mpp_macros.fi"
#include "macros/kernel_macros.fi"

    type, public :: block_info_type
        logical :: is_inner

        ! Local index of current block
        integer :: k
        ! Position in block grid
        integer :: bm, bn

        ! Info about boundary and edges
        !
        integer :: rank_nxp,     rank_nxm,     rank_nyp,     rank_nym
        integer :: rank_nxp_nyp, rank_nxp_nym, rank_nxm_nyp, rank_nxm_nym
        ! Local index of blocks near
        integer :: k_nxp,     k_nxm,     k_nyp,     k_nym
        integer :: k_nxp_nyp, k_nxp_nym, k_nxm_nyp, k_nxm_nym
    end type block_info_type

    ! Domain decomposition only in 2D area, s-levels are the same for each process, decomposition into rectangular blocks
    ! There are two types of block numeration:
    !   Local:
    !       for each proc k = 1, ..., bcount. In custom order!
    !   Global:
    !       Cart coords of blocks (1, 1), ..., (bnx, bny)
    type, public :: domain_type
        ! Local area
        integer :: bcount ! Number of data blocks at each process
        integer :: bcount_inner, bcount_boundary ! bcount_inner + bcount_boundary = bcount.

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
        ! Map: local block number to block coords: m=bindx(k,1), n=bindx(k,2)
        integer, pointer :: bindx(:, :)
        ! Map: block coords to local block number
        integer, pointer :: bglob_local_num(:, :)

        ! Info
        type(block_info_type), pointer :: blocks_info(:)
        integer, pointer :: ranks_near(:)
        integer :: amount_of_ranks_near
        integer, pointer :: boundary_blocks(:)

        ! Z dir
        integer :: nz
    contains
        procedure, public :: init_from_config
        procedure, public :: init
        procedure, public :: clear
        
        procedure, private :: block_uniform_decomposition
        procedure, private :: create_uniform_decomposition
        procedure, private :: create_hilbert_curve_decomposition
        procedure, private :: create_local_block_numearation
    end type domain_type

!------------------------------------------------------------------------------

    public :: get_boundary_points_of_block, get_halo_points_of_block
    public :: get_inverse_dir, get_rank_of_block, get_loc_num_of_block, fill_direction_array
    public :: is_block_on_current_proc, is_inner_block

    type(domain_type), public, target :: domain_data

contains

    subroutine get_boundary_points_of_block(domain, k, dir, nx_start, nx_end, ny_start, ny_end)
        type(domain_type), intent(in) :: domain
        integer, intent(in) :: k, dir
        integer, intent(out) :: nx_start, nx_end, ny_start, ny_end

        select case(dir)
            ! borders
            case(_NXP_)
                nx_start = domain%bnx_end(k)
                nx_end   = domain%bnx_end(k)
                ny_start = domain%bny_start(k)
                ny_end   = domain%bny_end(k)

            case(_NXM_)
                nx_start = domain%bnx_start(k)
                nx_end   = domain%bnx_start(k)
                ny_start = domain%bny_start(k)
                ny_end   = domain%bny_end(k)

            case(_NYP_)
                nx_start = domain%bnx_start(k)
                nx_end   = domain%bnx_end(k)
                ny_start = domain%bny_end(k)
                ny_end   = domain%bny_end(k)

            case(_NYM_)
                nx_start = domain%bnx_start(k)
                nx_end   = domain%bnx_end(k)
                ny_start = domain%bny_start(k)
                ny_end   = domain%bny_start(k)

            ! edges
            case (_NXP_NYP_)
                nx_start = domain%bnx_end(k)
                nx_end   = domain%bnx_end(k)
                ny_start = domain%bny_end(k)
                ny_end   = domain%bny_end(k)

            case (_NXP_NYM_)
                nx_start = domain%bnx_end(k)
                nx_end   = domain%bnx_end(k)
                ny_start = domain%bny_start(k)
                ny_end   = domain%bny_start(k)

            case (_NXM_NYP_)
                nx_start = domain%bnx_start(k)
                nx_end   = domain%bnx_start(k)
                ny_start = domain%bny_end(k)
                ny_end   = domain%bny_end(k)

            case (_NXM_NYM_)
                nx_start = domain%bnx_start(k)
                nx_end   = domain%bnx_start(k)
                ny_start = domain%bny_start(k)
                ny_end   = domain%bny_start(k)

            case default
                call abort_model("Unknown direction in get_boundary_points_of_block")
        end select
        
    end subroutine

    subroutine get_halo_points_of_block(domain, k, dir, nx_start, nx_end, ny_start, ny_end)
        type(domain_type), intent(in) :: domain
        integer, intent(in) :: k, dir
        integer, intent(out) :: nx_start, nx_end, ny_start, ny_end

        select case(dir)
            ! borders
            case(_NXP_)
                nx_start = domain%bnx_end(k) + 1
                nx_end   = domain%bnx_end(k) + 1
                ny_start = domain%bny_start(k)
                ny_end   = domain%bny_end(k)

            case(_NXM_)
                nx_start = domain%bnx_start(k) - 1
                nx_end   = domain%bnx_start(k) - 1
                ny_start = domain%bny_start(k)
                ny_end   = domain%bny_end(k)

            case(_NYP_)
                nx_start = domain%bnx_start(k)
                nx_end   = domain%bnx_end(k)
                ny_start = domain%bny_end(k) + 1
                ny_end   = domain%bny_end(k) + 1

            case(_NYM_)
                nx_start = domain%bnx_start(k)
                nx_end   = domain%bnx_end(k)
                ny_start = domain%bny_start(k) - 1
                ny_end   = domain%bny_start(k) - 1

            ! edges
            case (_NXP_NYP_)
                nx_start = domain%bnx_end(k) + 1
                nx_end   = domain%bnx_end(k) + 1
                ny_start = domain%bny_end(k) + 1
                ny_end   = domain%bny_end(k) + 1

            case (_NXP_NYM_)
                nx_start = domain%bnx_end(k) + 1
                nx_end   = domain%bnx_end(k) + 1
                ny_start = domain%bny_start(k) - 1
                ny_end   = domain%bny_start(k) - 1

            case (_NXM_NYP_)
                nx_start = domain%bnx_start(k) - 1
                nx_end   = domain%bnx_start(k) - 1
                ny_start = domain%bny_end(k) + 1
                ny_end   = domain%bny_end(k) + 1

            case (_NXM_NYM_)
                nx_start = domain%bnx_start(k) - 1
                nx_end   = domain%bnx_start(k) - 1
                ny_start = domain%bny_start(k) - 1
                ny_end   = domain%bny_start(k) - 1

            case default
                call abort_model("Unknown direction in get_halo_points_of_block")
        end select
        
    end subroutine

    function get_inverse_dir(dir) result(inverse_dir)
        integer, intent(in) :: dir
        integer :: inverse_dir

        select case(dir)
            ! borders
            case(_NXP_)
                inverse_dir = _NXM_
            case(_NXM_)
                inverse_dir = _NXP_
            case(_NYP_)
                inverse_dir = _NYM_
            case(_NYM_)
                inverse_dir = _NYP_
            ! edges
            case (_NXP_NYP_)
                inverse_dir = _NXM_NYM_
            case (_NXP_NYM_)
                inverse_dir = _NXM_NYP_
            case (_NXM_NYP_)
                inverse_dir = _NXP_NYM_
            case (_NXM_NYM_)
                inverse_dir = _NXP_NYP_
            case default
                call abort_model("Unknown direction in get_inverse_dir")
            end select

            return
    end function

    function get_rank_of_block(bm, bn, bglob_proc, bnx, bny) result(r)
        integer, intent(in) :: bm, bn
        integer, pointer, intent(in) :: bglob_proc(:, :)
        integer, intent(in) :: bnx, bny
        integer :: r

        if (bm > bnx .or. bm < 1 .or. bn > bny .or. bn < 1) then 
            r = -2
        else
            r = bglob_proc(bm, bn)
        endif
        return
    end function

    function get_loc_num_of_block(bm, bn, bglob_loc_num, bnx, bny) result(num)
        integer, intent(in) :: bm, bn
        integer, pointer, intent(in) :: bglob_loc_num(:, :)
        integer, intent(in) :: bnx, bny
        integer :: num

        if (bm > bnx .or. bm < 1 .or. bn > bny .or. bn < 1) then 
            num = -2
        else
            num = bglob_loc_num(bm, bn)
        endif
        return
    end function

    subroutine fill_direction_array(m, n, mdirs, ndirs)
        integer, intent(in) :: m, n
        integer, intent(inout) :: mdirs(8), ndirs(8)

        mdirs(1) = m + 1; ndirs(1) = n     ! nx+
        mdirs(2) = m - 1; ndirs(2) = n     ! nx-
        mdirs(3) = m;     ndirs(3) = n + 1 ! ny+
        mdirs(4) = m;     ndirs(4) = n - 1 ! ny-
        
        mdirs(5) = m + 1; ndirs(5) = n + 1 ! nx+ny+
        mdirs(6) = m + 1; ndirs(6) = n - 1 ! nx+ny-
        mdirs(7) = m - 1; ndirs(7) = n + 1 ! nx-ny+
        mdirs(8) = m - 1; ndirs(8) = n - 1 ! nx-ny-

        return
    end subroutine

    function is_block_on_current_proc(bm, bn, bglob_proc, bnx, bny) result(is)
        integer, intent(in) :: bm, bn
        integer, pointer, intent(in) :: bglob_proc(:, :)
        integer, intent(in) :: bnx, bny
        logical :: is

        is = .true.

        if (bm > bnx .or. bm < 1 .or. bn > bny .or. bn < 1) then 
            is = .false.
            return
        endif
        if (bglob_proc(bm, bn) /= mpp_rank) then
            is = .false.
            return
        endif

        return
    end function

    function is_inner_block(bm, bn, bglob_proc, bnx, bny) result(is)
        integer, intent(in) :: bm, bn
        integer, pointer, intent(in) :: bglob_proc(:, :)
        integer, intent(in) :: bnx, bny
        logical :: is

        integer :: bmdirs(8), bndirs(8)
        integer :: k

        is = .true.

        call fill_direction_array(bm, bn, bmdirs, bndirs)

        do k = 1, 8
            if (.not. is_block_on_current_proc(bmdirs(k), bndirs(k), bglob_proc, bnx, bny)) then
                is = .false.
                return
            endif
        enddo
        
        return
    end function

    subroutine check_domain_data(domain)
        type(domain_type), intent(in) :: domain

        integer :: k

        do k = 1, domain%bcount_boundary
            if (domain%blocks_info( domain%boundary_blocks(k) )%is_inner) then
                call abort_model("Incorrect boundaru blocks")
            endif
        enddo
        if (domain%bcount_inner + domain%bcount_boundary /= domain%bcount) then
            print *, 'DD INFO: bcount', domain%bcount
            call abort_model("Incorrect amount of inner and boundary blocks")
        endif
        
    end subroutine

    subroutine block_uniform_decomposition(this, lbasins,  &
                                           glob_bnx_start, glob_bnx_end, glob_bny_start, glob_bny_end,  & 
                                           glob_bbnd_x1, glob_bbnd_x2, glob_bbnd_y1, glob_bbnd_y2,  &
                                           bglob_weight, land_blocks)
        use config_basinpar_module, only: nx, ny
        
        class(domain_type), intent(in) :: this
        integer, allocatable, intent(in) :: lbasins(:,:) 
        integer, allocatable, intent(inout) :: glob_bnx_start(:, :), glob_bnx_end(:, :),  &
                                             glob_bny_start(:, :), glob_bny_end(:, :)
        integer, allocatable, intent(inout) :: glob_bbnd_x1(:, :), glob_bbnd_x2(:, :),  &
                                             glob_bbnd_y1(:, :), glob_bbnd_y2(:, :)
        real(wp8), allocatable, intent(inout) :: bglob_weight(:, :)
        integer, intent(out) :: land_blocks
        
        integer :: m, n, i, j, ierr, xsize, ysize, total_x_size, total_y_size
        integer, allocatable :: x_sizes(:), y_sizes(:), x_start(:), y_start(:)
        
        associate(bnx => this%bnx,  &
                  bny => this%bny)

            ! Compute optimal sizes x of blocks
            allocate(x_sizes(bnx), x_start(bnx))
            total_x_size = 0
            do m = 1, bnx
                if (m .eq. bnx) then
                    xsize = (nx - 4) - total_x_size
                else
                    xsize = floor(real((nx - 4) - total_x_size)/real(bnx - m + 1))
                endif

                if (xsize <= 0) then
                    call abort_model('Error in decomposition to uniform blocks: x size less or equal zero')
                endif

                x_sizes(m) = xsize
                x_start(m) = total_x_size
                total_x_size = total_x_size + xsize
            enddo
            if (SUM(x_sizes) /= (nx - 4)) call abort_model('Error in decomposition to uniform blocks: x total size not equal to nx - 4')

            ! Compute optimal sizes y of blocks
            allocate(y_sizes(bny), y_start(bny))
            total_y_size = 0
            do n = 1, bny
                if (n .eq. bny) then
                    ysize = (ny - 4) - total_y_size
                else
                    ysize = floor(real((ny - 4) - total_y_size)/real(bny - n + 1))
                endif

                if (ysize <= 0) then
                    call abort_model('Error in decomposition to uniform blocks: y size less or equal zero')
                endif

                y_sizes(n) = ysize
                y_start(n) = total_y_size
                total_y_size = total_y_size + ysize
            enddo
            if (SUM(y_sizes) /= (ny - 4)) call abort_model('Error in decomposition to uniform blocks: y total size not equal to ny - 4')
            
            ! Fill data with uniform decomposition and accumulate land blocks
            land_blocks = 0
            bglob_weight = 0.0
            do m = 1, bnx
                do n = 1, bny
                    glob_bnx_start(m, n) = 3 + x_start(m)
                    glob_bnx_end(m, n) = glob_bnx_start(m, n) + x_sizes(m) - 1
                    ! border area
                    glob_bbnd_x1(m, n) = glob_bnx_start(m, n) - 2
                    glob_bbnd_x2(m, n) = glob_bnx_end(m, n) + 2

                    glob_bny_start(m, n) = 3 + y_start(n)
                    glob_bny_end(m, n) = glob_bny_start(m, n) + y_sizes(n) - 1
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
            if (mpp_is_master()) print *, "DD INFO: Total land blocks:", land_blocks

            ierr = 0
            if (bnx*bny - land_blocks < mpp_count) ierr = 1
            call check_error(ierr,  'procs > computational-blocks... Error!')

            deallocate(x_sizes, y_sizes, x_start, y_start)
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
            if (mpp_is_master()) print *, 'DD INFO: bnx not equal to bny! Can`t build Hilbert curve for this geometry!'
            ierr = 1
        endif
        if (2**hilbert_index /= this%bnx) then
            if (mpp_is_master()) print *, 'DD INFO: 2**M not eqal to bnx! Can`t build Hilbert curve for this geometry!'
            ierr = 1
        endif
        call check_error(ierr, 'Can`t build Hilbert curve for this geometry!')

        if (mpp_is_master()) print *, 'DD INFO: Hilbert curve index:', hilbert_index

        tot_weight = sum(bglob_weight)
        mean_weight = tot_weight / mpp_count
        sea_blocks = this%bnx * this%bny - land_blocks

        if (debug_level >= 1) then
            if (mpp_is_master()) print *, 'DD INFO: Total blocks weigth:', tot_weight, "Mean blocks weigth:", mean_weight
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
                        if (mpp_is_master()) print *, 'DD INFO: for block', k, 'Warning! Last procs ...'
                    endif
                endif
            !else
            !    i = i + 1
            !    weight = bglob_weight(hilbert_coord_x, hilbert_coord_y)
            !endif

            this%bglob_proc(hilbert_coord_x, hilbert_coord_y) = i
            last_weight = last_weight + bglob_weight(hilbert_coord_x, hilbert_coord_y)
        enddo

        if (debug_level >= 7) then
            call parallel_int_output(this%bglob_proc, 1, this%bnx, 1, this%bny, 'bglob_proc from load-balanced hilbert curve decomposition')
        endif
    end subroutine

    subroutine create_uniform_decomposition(this, bglob_weight)
        class(domain_type), intent(inout) :: this
        real(wp8), allocatable, intent(in) :: bglob_weight(:, :)
        integer :: m, n, rank, ierr
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
            if (debug_level >= 1) then 
                do rank = 0, mpp_count - 1
                    if (mpp_rank == rank) then
                        print *, 'DD INFO: Uniform decomposition. rank, xb_start, yb_start', mpp_rank, xblock_start, yblock_start
                    endif
                    call mpi_barrier(mpp_cart_comm, ierr)
                enddo
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
    
            if (debug_level >= 7) then
                call parallel_int_output(bgproc, 1, bnx, 1, bny, 'DD INFO: bglob_proc from uniform decomposition')
            endif
    
            deallocate(buf_int)

        end associate
    end subroutine

    subroutine create_local_block_numearation(this,  &
                                              glob_bnx_start, glob_bnx_end, glob_bny_start, glob_bny_end,  & 
                                              glob_bbnd_x1, glob_bbnd_x2, glob_bbnd_y1, glob_bbnd_y2,  &
                                              bglob_weight)
        class(domain_type), intent(inout) :: this
        integer, allocatable, intent(inout) :: glob_bnx_start(:, :), glob_bnx_end(:, :),  &
                                               glob_bny_start(:, :), glob_bny_end(:, :)
        integer, allocatable, intent(inout) :: glob_bbnd_x1(:, :), glob_bbnd_x2(:, :),  &
                                               glob_bbnd_y1(:, :), glob_bbnd_y2(:, :)
        real(wp8), allocatable, intent(inout) :: bglob_weight(:, :)

        integer :: k, kk, k_inner, k_boundary, m, n, ierr
        integer, allocatable :: buf_int(:, :)
        
        logical, allocatable :: mask_blocks(:, :)
        integer, dimension(2) :: max_mn
        integer :: max_m, max_n

        ! Allocate blocks arrays per proc
        allocate(this%bindx(this%bcount, 2))
        allocate(this%bnx_start(this%bcount), this%bnx_end(this%bcount),  &
                 this%bny_start(this%bcount), this%bny_end(this%bcount))
        allocate(this%bbnd_x1(this%bcount), this%bbnd_x2(this%bcount),  &
                 this%bbnd_y1(this%bcount), this%bbnd_y2(this%bcount))
        this%bindx = 0
        this%bnx_start = 0; this%bnx_end = 0 
        this%bny_start = 0; this%bny_end = 0
        this%bbnd_x1 = 0; this%bbnd_x2 = 0 
        this%bbnd_y1 = 0; this%bbnd_y2 = 0

        ! Create local block numeration
        allocate(this%bglob_local_num(this%bnx, this%bny))
        this%bglob_local_num = -1

#ifdef _MPP_SORTED_BLOCKS_
            allocate(mask_blocks(this%bnx, this%bny))
            ! Set sorted by weight blocks
            k = 0; k_inner = 0; k_boundary = 0
            ! Mask only blocks on current proc
            mask_blocks = .false.
            do m = 1, this%bnx
                do n = 1, this%bny
                    if (this%bglob_proc(m, n) == mpp_rank) then
                        mask_blocks(m ,n) = .true.
                        ! Count inner and boundary blocks
                        if (is_inner_block(m, n, this%bglob_proc, this%bnx, this%bny)) then
                            k_inner = k_inner + 1
                        else
                            k_boundary = k_boundary + 1
                        endif
                    endif
                enddo
            enddo
                
            ! Sort blocks by weight
            do kk = 1, k_inner + k_boundary
                ! Count block
                k = k + 1

                max_mn = MAXLOC(bglob_weight, mask_blocks)
                max_m = max_mn(1)
                max_n = max_mn(2)

                mask_blocks(max_m, max_n) = .false.
                
                ! Map local block numeration to block coords
                this%bindx(k, 1) = max_m
                this%bindx(k, 2) = max_n

                this%bnx_start(k) = glob_bnx_start(max_m, max_n)
                this%bnx_end(k) = glob_bnx_end(max_m, max_n)
                this%bny_start(k) = glob_bny_start(max_m, max_n)
                this%bny_end(k) = glob_bny_end(max_m, max_n)

                this%bbnd_x1(k) = glob_bbnd_x1(max_m, max_n)
                this%bbnd_x2(k) = glob_bbnd_x2(max_m, max_n)
                this%bbnd_y1(k) = glob_bbnd_y1(max_m, max_n)
                this%bbnd_y2(k) = glob_bbnd_y2(max_m, max_n)
                
                if (debug_level >= 9) then
                    if (mpp_is_master()) print *, 'DD INFO: SORTED_BLOCKS: k, max_m, max_n, bw, is_inner?',  &
                                                    k, max_m, max_n, bglob_weight(max_m, max_n),  &
                                                    is_inner_block(max_m, max_n, this%bglob_proc, this%bnx, this%bny)
                endif

                this%bglob_local_num(max_m, max_n) = k
            enddo
            deallocate(mask_blocks)
#else
            k = 0; k_inner = 0; k_boundary = 0
            do m = 1, this%bnx
                do n = 1, this%bny
                    if (this%bglob_proc(m, n) == mpp_rank) then
                        ! Count block
                        k = k + 1
                        if (is_inner_block(m, n, this%bglob_proc, this%bnx, this%bny)) then
                            k_inner = k_inner + 1
                        else
                            k_boundary = k_boundary + 1
                        endif

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

                        this%bglob_local_num(m, n) = k
                    endif
                enddo
            enddo
#endif
        ! Sync bglob_local_num
        allocate(buf_int(this%bnx, this%bny))
        buf_int = this%bglob_local_num
        call mpi_allreduce(buf_int, this%bglob_local_num, this%bnx * this%bny, mpi_integer, mpi_max, mpp_cart_comm, ierr)
        if (debug_level >= 7) then
            call parallel_int_output(this%bglob_local_num, 1, this%bnx, 1, this%bny, 'DD INFO: bglob_local_num from uniform decomposition')
        endif
        deallocate(buf_int)

        if (k /= this%bcount .or. k_inner + k_boundary /= this%bcount) then
            print *, 'DD INFO:', mpp_rank, " k=", k, " bcount=", this%bcount, " k_inner=", k_inner, " k_boundary=", k_boundary
            call abort_model("Error in block numeration")
        endif

        ! Set values
        this%bcount_inner = k_inner
        this%bcount_boundary = k_boundary

    end subroutine

    subroutine init(this, bppnx, bppny, mod_create, lbasins)
        use config_basinpar_module, only: nx, ny, nz

        ! Initialization of each domain
        class(domain_type), intent(inout) :: this
        integer, intent(in) :: bppnx, bppny
        integer, intent(in) :: mod_create
        integer, allocatable, intent(in) :: lbasins(:,:) 

        real(wp8), allocatable :: bglob_weight(:, :)
        real(wp8) :: bweight, max_bweight

        integer, allocatable :: glob_bnx_start(:, :), glob_bnx_end(:, :),  &
                                glob_bny_start(:, :), glob_bny_end(:, :)
        integer, allocatable :: glob_bbnd_x1(:, :), glob_bbnd_x2(:, :),  &
                                glob_bbnd_y1(:, :), glob_bbnd_y2(:, :)

        integer :: m, n, k, k_boundary, kk, kkk, res, rank
        integer :: land_blocks
        integer :: ierr

        integer :: itot, ishared
        integer :: max_points_per_block, min_points_per_block, reduced_max_points_per_block, reduced_min_points_per_block
        integer :: buf_dirs(8)

        integer :: max_x, max_y, min_x, min_y

        real(wp8) :: max_w

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

            if (mpp_is_master()) print *, 'DD INFO: bnx, bny and Total blocks:', bnx, bny, bnx * bny
            if (mpp_is_master()) print *, 'DD INFO: pnx, pny and procs:', mpp_size(1), mpp_size(2), mpp_count
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
                if (mpp_is_master()) print *, "DD INFO: Uniform blocks decomposition!..."
                call this%create_uniform_decomposition(bglob_weight)
            elseif (mod_create == 1) then
                if (mpp_is_master()) print *, "DD INFO: Hilber Curve blocks decomposition!..."
                call this%create_hilbert_curve_decomposition(bglob_weight, land_blocks)
            else
                if (mpp_is_master()) print *, "DD INFO: Unknown mode!"
                call abort_model("Unknown decomposition mode!")
            endif

            ! Print decomposition in file for visualization
            if (debug_level >= 3) then
                if (mpp_rank == 0) then
                    print *, "DD INFO: Print decomposition in file decomposition.txt"
                    open(90, file = "decomposition.txt")
                    write(90, *)  bnx, bny, mpp_size(1), mpp_size(2)
                    do m = 1, bnx
                        do n = 1, bny
                            write(90, *) m, n, this%bglob_proc(m, n), bglob_weight(m, n)
                        enddo
                    enddo
                    close(90)
                    print *, "DD INFO: Print decomposition in file decomposition.txt is OK"
                endif
                call mpi_barrier(mpp_cart_comm, ierr)
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
            if (mpp_is_master()) print *, 'DD INFO: Total blocks:', total_blocks, 'LB: ', max_bweight / (sum(bglob_weight) / real(mpp_count)),  &
                                          'max blocks per proc:', bcount_max, 'min blocks per proc:', bcount_min
            call mpi_barrier(mpp_cart_comm, ierr)

            call this%create_local_block_numearation(glob_bnx_start, glob_bnx_end, glob_bny_start, glob_bny_end,  & 
                                                     glob_bbnd_x1, glob_bbnd_x2, glob_bbnd_y1, glob_bbnd_y2,  &
                                                     bglob_weight)

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
                call mpi_allreduce(min_points_per_block, reduced_min_points_per_block, 1, mpi_integer, mpi_min, mpp_cart_comm, ierr)
                call mpi_allreduce(max_points_per_block, reduced_max_points_per_block, 1, mpi_integer, mpi_max, mpp_cart_comm, ierr)
                if (mpp_is_master()) print *, 'DD INFO: Max points per block:', reduced_max_points_per_block, 'Min points per blocks: ', reduced_min_points_per_block
                call mpi_barrier(mpp_cart_comm, ierr)
            endif

            ! Init blocks_info: additional information about each block, especially for sync
            allocate(this%blocks_info(bcount))
            allocate(this%ranks_near(mpp_count))
            allocate(this%boundary_blocks(this%bcount_boundary))
            this%boundary_blocks = 0
            this%ranks_near = -1
            this%amount_of_ranks_near = 0
            k_boundary = 0
            do k = 1, bcount
                m = this%bindx(k, 1)
                n = this%bindx(k, 2)

                this%blocks_info(k)%is_inner = is_inner_block(m, n, this%bglob_proc, bnx, bny)
                this%blocks_info(k)%k = k
                this%blocks_info(k)%bm = m
                this%blocks_info(k)%bn = n
                
                ! Rank pf procs near
                this%blocks_info(k)%rank_nxp = get_rank_of_block(m + 1, n,     this%bglob_proc, bnx, bny)
                this%blocks_info(k)%rank_nxm = get_rank_of_block(m - 1, n,     this%bglob_proc, bnx, bny)
                this%blocks_info(k)%rank_nyp = get_rank_of_block(m,     n + 1, this%bglob_proc, bnx, bny)
                this%blocks_info(k)%rank_nym = get_rank_of_block(m,     n - 1, this%bglob_proc, bnx, bny)
                
                this%blocks_info(k)%rank_nxp_nyp = get_rank_of_block(m + 1, n + 1, this%bglob_proc, bnx, bny)
                this%blocks_info(k)%rank_nxp_nym = get_rank_of_block(m + 1, n - 1, this%bglob_proc, bnx, bny)
                this%blocks_info(k)%rank_nxm_nyp = get_rank_of_block(m - 1, n + 1, this%bglob_proc, bnx, bny)
                this%blocks_info(k)%rank_nxm_nym = get_rank_of_block(m - 1, n - 1, this%bglob_proc, bnx, bny)

                ! Local index of blocks near
                this%blocks_info(k)%k_nxp = get_loc_num_of_block(m + 1, n,     this%bglob_local_num, bnx, bny) 
                this%blocks_info(k)%k_nxm = get_loc_num_of_block(m - 1, n,     this%bglob_local_num, bnx, bny)
                this%blocks_info(k)%k_nyp = get_loc_num_of_block(m,     n + 1, this%bglob_local_num, bnx, bny)
                this%blocks_info(k)%k_nym = get_loc_num_of_block(m,     n - 1, this%bglob_local_num, bnx, bny)

                this%blocks_info(k)%k_nxp_nyp = get_loc_num_of_block(m + 1, n + 1, this%bglob_local_num, bnx, bny)
                this%blocks_info(k)%k_nxp_nym = get_loc_num_of_block(m + 1, n - 1, this%bglob_local_num, bnx, bny)
                this%blocks_info(k)%k_nxm_nyp = get_loc_num_of_block(m - 1, n + 1, this%bglob_local_num, bnx, bny)
                this%blocks_info(k)%k_nxm_nym = get_loc_num_of_block(m - 1, n - 1, this%bglob_local_num, bnx, bny)

                ! Gets ranks near
                buf_dirs(1) = this%blocks_info(k)%rank_nxp
                buf_dirs(2) = this%blocks_info(k)%rank_nxm
                buf_dirs(3) = this%blocks_info(k)%rank_nyp
                buf_dirs(4) = this%blocks_info(k)%rank_nym
                buf_dirs(5) = this%blocks_info(k)%rank_nxp_nyp
                buf_dirs(6) = this%blocks_info(k)%rank_nxp_nym
                buf_dirs(7) = this%blocks_info(k)%rank_nxm_nyp
                buf_dirs(8) = this%blocks_info(k)%rank_nxm_nym

                ! Boundary blocks
                if (.not. this%blocks_info(k)%is_inner) then
                    this%boundary_blocks(k_boundary + 1) = k
                    k_boundary = k_boundary + 1
                endif

                do kk = 1, 8
                    if (buf_dirs(kk) >= 0 .and. buf_dirs(kk) /= mpp_rank) then
                        ! Find this rank in ranks_near
                        res = 0
                        do kkk = 1, this%amount_of_ranks_near
                            if (this%ranks_near(kkk) == buf_dirs(kk)) then
                                res = 1
                            endif
                        enddo
                        ! If rank doesnt exist - count it
                        if (res == 0) then
                            this%amount_of_ranks_near = this%amount_of_ranks_near + 1
                            this%ranks_near(this%amount_of_ranks_near) = buf_dirs(kk)
                        endif
                    endif
                enddo

                if (debug_level >= 7 .and. mpp_is_master ()) then
                    print *, 'DD INFO: master: block k, is_inner', k, this%blocks_info(k)%is_inner
                    print *, 'DD INFO: master: block k, k', k, this%blocks_info(k)%k
                    print *, 'DD INFO: master: block k, bm', k, this%blocks_info(k)%bm
                    print *, 'DD INFO: master: block k, bn', k, this%blocks_info(k)%bn
                    print *, 'DD INFO: master: block k, rank_nxp', k, this%blocks_info(k)%rank_nxp
                    print *, 'DD INFO: master: block k, rank_nxm', k, this%blocks_info(k)%rank_nxm
                    print *, 'DD INFO: master: block k, rank_nyp', k, this%blocks_info(k)%rank_nyp
                    print *, 'DD INFO: master: block k, rank_nym', k, this%blocks_info(k)%rank_nym
                    print *, 'DD INFO: master: block k, rank_nxp_nyp', k, this%blocks_info(k)%rank_nxp_nyp
                    print *, 'DD INFO: master: block k, rank_nxp_nym', k, this%blocks_info(k)%rank_nxp_nym
                    print *, 'DD INFO: master: block k, rank_nxm_nyp', k, this%blocks_info(k)%rank_nxm_nyp
                    print *, 'DD INFO: master: block k, rank_nxm_nym', k, this%blocks_info(k)%rank_nxm_nym
                    print *, 'DD INFO: master: block k, k_nxp', k, this%blocks_info(k)%k_nxp
                    print *, 'DD INFO: master: block k, k_nxm', k, this%blocks_info(k)%k_nxm
                    print *, 'DD INFO: master: block k, k_nyp', k, this%blocks_info(k)%k_nyp
                    print *, 'DD INFO: master: block k, k_nym', k, this%blocks_info(k)%k_nym
                    print *, 'DD INFO: master: block k, k_nxp_nyp', k, this%blocks_info(k)%k_nxp_nyp
                    print *, 'DD INFO: master: block k, k_nxp_nym', k, this%blocks_info(k)%k_nxp_nym
                    print *, 'DD INFO: master: block k, k_nxm_nyp', k, this%blocks_info(k)%k_nxm_nyp
                    print *, 'DD INFO: master: block k, k_nxm_nym', k, this%blocks_info(k)%k_nxm_nym
                endif
            enddo

            ! Found max/min x and y points of compute area for each rank
            max_x = 0; max_y = 0
            min_x = nx; min_y = ny
            do k = 1, bcount
                if (this%bnx_start(k) < min_x) min_x = this%bnx_start(k)
                if (this%bnx_end(k) > max_x) max_x = this%bnx_end(k)
                if (this%bny_start(k) < min_y) min_y = this%bny_start(k)
                if (this%bny_end(k) > max_y) max_y = this%bny_end(k)
            enddo

            ! Info per rank
            if (debug_level >= 5) then
                do rank = 0, mpp_count - 1
                    if (mpp_rank == rank) then
                        print *, 'DD INFO: rank', mpp_rank, 'Blocks:', bcount, 'Min/Max points per block:', min_points_per_block, max_points_per_block,  &
                                 "Inner/Boundary blocks:", this%bcount_inner, this%bcount_boundary,  &
                                 "Amount of ranks near:", this%amount_of_ranks_near
                        print *, 'DD INFO: rank', mpp_rank, 'Min/Max nx', min_x, max_x, 'Min/Max ny', min_y, max_y
                    endif
                    call mpi_barrier(mpp_cart_comm, ierr)
                enddo
                call mpi_barrier(mpp_cart_comm, ierr)
            endif

            ! Info per thread
            if (debug_level >= 5) then
                max_w = 0.0d0
                res = 0
                !$omp parallel default(shared) private(kk, bweight)
                bweight = 0.0d0
                kk = 0
                !$omp do private(k) schedule(static, 1) reduction(max: max_w, res)
                do k = 1, bcount
                    bweight = bweight + bglob_weight(this%blocks_info(k)%bm, this%blocks_info(k)%bn)
                    kk = kk + 1

                    max_w = MAX(max_w, bweight)
                    res = MAX(res, kk)
                enddo
                !$omp end do

                if (debug_level >= 6) then
                    print *, 'DD INFO: rank, thread', mpp_rank, mpp_get_thread(), 'Blocks: ', kk, 'Total Weight: ', bweight
                endif
                
                !$omp end parallel

                print *, 'DD INFO: rank', mpp_rank, 'Max Blocks per threads: ', res, 'Max Total Weight per threads: ', max_w
            endif

            ! Info per block
            if (debug_level >= 9) then
                if (mpp_is_master()) then
                    do k = 1, bcount
                        print *, 'DD INFO: master: block:', k, ' is inner? ', this%blocks_info(k)%is_inner
                        write(*, '(A16, A42, i5, i5,i5,i5,i5)') 'DD INFO: master', 'k, nxs, nxe, nys, nye', k, this%bnx_start(k), this%bnx_end(k), this%bny_start(k), this%bny_end(k)
                        write(*, '(A16, A42, i5, i5,i5,i5,i5)') 'DD INFO: master', 'k, bndxs, bndxe, bndys, bndye', k, this%bbnd_x1(k), this%bbnd_x2(k), this%bbnd_y1(k), this%bbnd_y2(k)
                    enddo
                endif
            endif

            call mpp_sync_output()
            
            ! Check domain data
            call check_domain_data(this)

            deallocate(bglob_weight)

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

        deallocate(this%bglob_local_num)
        deallocate(this%blocks_info)

        deallocate(this%ranks_near)

        deallocate(this%boundary_blocks)
    end subroutine

    subroutine init_from_config(this, lbasins, name)
        use rwpar_routes

        class(domain_type), intent(inout) :: this
        integer, allocatable, intent(in) :: lbasins(:,:) 
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

        call mpp_sync_output()

        call this%init(bppnx, bppny, mod_decomposition, lbasins)

    end subroutine

end module decomposition_module