module mpp_sync_module
    ! Sync data of different types, depends on domain

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use mpi
    use debug_module
    use mpp_module, only: mpp_rank, mpp_count, mpp_cart_comm, mpp_size, mpp_coord, mpp_time_sync, start_timer, end_timer
    use errors_module, only: abort_model, check_error
    use data_types_module, only: block1D_real4_type, block1D_real8_type, block2D_real4_type, block2D_real8_type,  &
                                 data1D_real4_type, data1D_real8_type, data2D_real4_type, data2D_real8_type,      &
                                 data3D_real4_type, data3D_real8_type
    use decomposition_module, only: domain_type

#include "macros/mpp_macros.fi"

    implicit none
    save
    private

!------------------------------------------------------------------------------
    ! MPI buffers for 2D sync
    integer, allocatable :: reqsts(:), statuses(:, :)
    ! real8
    type(block1D_real8_type), dimension(:), pointer :: sync_buf8_send_nyp, sync_buf8_recv_nyp
    type(block1D_real8_type), dimension(:), pointer :: sync_buf8_send_nxp, sync_buf8_recv_nxp
    type(block1D_real8_type), dimension(:), pointer :: sync_buf8_send_nym, sync_buf8_recv_nym
    type(block1D_real8_type), dimension(:), pointer :: sync_buf8_send_nxm, sync_buf8_recv_nxm
    real(wp8), allocatable :: sync_edge_buf8_recv_nxp_nyp(:)
    real(wp8), allocatable :: sync_edge_buf8_recv_nxp_nym(:)
    real(wp8), allocatable :: sync_edge_buf8_recv_nxm_nyp(:)
    real(wp8), allocatable :: sync_edge_buf8_recv_nxm_nym(:)
    ! real4
    type(block1D_real4_type), dimension(:), pointer :: sync_buf4_send_nyp, sync_buf4_recv_nyp
    type(block1D_real4_type), dimension(:), pointer :: sync_buf4_send_nxp, sync_buf4_recv_nxp
    type(block1D_real4_type), dimension(:), pointer :: sync_buf4_send_nym, sync_buf4_recv_nym
    type(block1D_real4_type), dimension(:), pointer :: sync_buf4_send_nxm, sync_buf4_recv_nxm
    real(wp4), allocatable :: sync_edge_buf4_recv_nxp_nyp(:)
    real(wp4), allocatable :: sync_edge_buf4_recv_nxp_nym(:)
    real(wp4), allocatable :: sync_edge_buf4_recv_nxm_nyp(:)
    real(wp4), allocatable :: sync_edge_buf4_recv_nxm_nym(:)
!------------------------------------------------------------------------------
    ! MPI buffers for 3D sync
    ! real8
    type(block2D_real8_type), dimension(:), pointer :: sync_buf8_send_nyp_3D, sync_buf8_recv_nyp_3D
    type(block2D_real8_type), dimension(:), pointer :: sync_buf8_send_nxp_3D, sync_buf8_recv_nxp_3D
    type(block2D_real8_type), dimension(:), pointer :: sync_buf8_send_nym_3D, sync_buf8_recv_nym_3D
    type(block2D_real8_type), dimension(:), pointer :: sync_buf8_send_nxm_3D, sync_buf8_recv_nxm_3D
    type(block1D_real8_type), dimension(:), pointer :: sync_edge_buf8_recv_nxp_nyp_3D
    type(block1D_real8_type), dimension(:), pointer :: sync_edge_buf8_recv_nxp_nym_3D
    type(block1D_real8_type), dimension(:), pointer :: sync_edge_buf8_recv_nxm_nyp_3D
    type(block1D_real8_type), dimension(:), pointer :: sync_edge_buf8_recv_nxm_nym_3D
    ! real4
    type(block2D_real4_type), dimension(:), pointer :: sync_buf4_send_nyp_3D, sync_buf4_recv_nyp_3D
    type(block2D_real4_type), dimension(:), pointer :: sync_buf4_send_nxp_3D, sync_buf4_recv_nxp_3D
    type(block2D_real4_type), dimension(:), pointer :: sync_buf4_send_nym_3D, sync_buf4_recv_nym_3D
    type(block2D_real4_type), dimension(:), pointer :: sync_buf4_send_nxm_3D, sync_buf4_recv_nxm_3D
    type(block1D_real4_type), dimension(:), pointer :: sync_edge_buf4_recv_nxp_nyp_3D
    type(block1D_real4_type), dimension(:), pointer :: sync_edge_buf4_recv_nxp_nym_3D
    type(block1D_real4_type), dimension(:), pointer :: sync_edge_buf4_recv_nxm_nyp_3D
    type(block1D_real4_type), dimension(:), pointer :: sync_edge_buf4_recv_nxm_nym_3D
!------------------------------------------------------------------------------
    integer, private :: parallel_dbg
    !real(wp8), private :: time_sync
!------------------------------------------------------------------------------

    interface sync
        module procedure syncborder_data2D_real8
        module procedure syncborder_data2D_real4
        module procedure syncborder_data3D_real8
        module procedure syncborder_data3D_real4
    end interface

    public :: mpp_sync_init
    public :: mpp_sync_finalize
    public :: sync

    ! Private
    private :: allocate_mpp_sync_buffers
    private :: deallocate_mpp_sync_buffers

    private :: get_local_block_number
    private :: check_block_status
    private :: check_cart_coord

    private :: irecv_real8
    private :: irecv_real4
    private :: isend_real8
    private :: isend_real4
    private :: irecv3D_real8
    private :: irecv3D_real4
    private :: isend3D_real8
    private :: isend3D_real4

contains

    subroutine mpp_sync_init(domain)
        type(domain_type), intent(in) :: domain

        call allocate_mpp_sync_buffers(domain)
        parallel_dbg = debug_level
    end subroutine

    subroutine mpp_sync_finalize(domain)
        type(domain_type), intent(in) :: domain

        call deallocate_mpp_sync_buffers(domain)
    end subroutine

    subroutine allocate_mpp_sync_buffers(domain)
        type(domain_type), intent(in) :: domain
        integer :: k, nx_dir_size, ny_dir_size
        integer :: nzmax

        ! Set maximum nz levels
        nzmax = domain%nz + 2

        associate(bcount => domain%bcount,  &
                  bnx_start => domain%bnx_start,  &
                  bnx_end => domain%bnx_end,  &
                  bny_start => domain%bny_start,  &
                  bny_end => domain%bny_end)
            ! MPI buffers for 2D data
            ! real8
            allocate(sync_buf8_send_nyp(bcount), sync_buf8_recv_nyp(bcount))
            allocate(sync_buf8_send_nxp(bcount), sync_buf8_recv_nxp(bcount))
            allocate(sync_buf8_send_nym(bcount), sync_buf8_recv_nym(bcount))
            allocate(sync_buf8_send_nxm(bcount), sync_buf8_recv_nxm(bcount))
            do k = 1, bcount
                ny_dir_size = bnx_end(k) - bnx_start(k) + 1
                nx_dir_size = bny_end(k) - bny_start(k) + 1
                allocate(sync_buf8_send_nyp(k)%field(ny_dir_size), sync_buf8_recv_nyp(k)%field(ny_dir_size))
                allocate(sync_buf8_send_nxp(k)%field(nx_dir_size), sync_buf8_recv_nxp(k)%field(nx_dir_size))
                allocate(sync_buf8_send_nym(k)%field(ny_dir_size), sync_buf8_recv_nym(k)%field(ny_dir_size))
                allocate(sync_buf8_send_nxm(k)%field(nx_dir_size), sync_buf8_recv_nxm(k)%field(nx_dir_size))
            enddo
            allocate(sync_edge_buf8_recv_nxp_nyp(bcount))
            allocate(sync_edge_buf8_recv_nxp_nym(bcount))
            allocate(sync_edge_buf8_recv_nxm_nyp(bcount))
            allocate(sync_edge_buf8_recv_nxm_nym(bcount))
            ! real4
            allocate(sync_buf4_send_nyp(bcount), sync_buf4_recv_nyp(bcount))
            allocate(sync_buf4_send_nxp(bcount), sync_buf4_recv_nxp(bcount))
            allocate(sync_buf4_send_nym(bcount), sync_buf4_recv_nym(bcount))
            allocate(sync_buf4_send_nxm(bcount), sync_buf4_recv_nxm(bcount))
            do k = 1, bcount
                ny_dir_size = bnx_end(k) - bnx_start(k) + 1
                nx_dir_size = bny_end(k) - bny_start(k) + 1
                allocate(sync_buf4_send_nyp(k)%field(ny_dir_size), sync_buf4_recv_nyp(k)%field(ny_dir_size))
                allocate(sync_buf4_send_nxp(k)%field(nx_dir_size), sync_buf4_recv_nxp(k)%field(nx_dir_size))
                allocate(sync_buf4_send_nym(k)%field(ny_dir_size), sync_buf4_recv_nym(k)%field(ny_dir_size))
                allocate(sync_buf4_send_nxm(k)%field(nx_dir_size), sync_buf4_recv_nxm(k)%field(nx_dir_size))
            enddo
            allocate(sync_edge_buf4_recv_nxp_nyp(bcount))
            allocate(sync_edge_buf4_recv_nxp_nym(bcount))
            allocate(sync_edge_buf4_recv_nxm_nyp(bcount))
            allocate(sync_edge_buf4_recv_nxm_nym(bcount))

            ! MPI buffers for 3D data
            ! real8
            allocate(sync_buf8_send_nyp_3D(bcount), sync_buf8_recv_nyp_3D(bcount))
            allocate(sync_buf8_send_nxp_3D(bcount), sync_buf8_recv_nxp_3D(bcount))
            allocate(sync_buf8_send_nym_3D(bcount), sync_buf8_recv_nym_3D(bcount))
            allocate(sync_buf8_send_nxm_3D(bcount), sync_buf8_recv_nxm_3D(bcount))
            allocate(sync_edge_buf8_recv_nxp_nyp_3D(bcount))
            allocate(sync_edge_buf8_recv_nxp_nym_3D(bcount))
            allocate(sync_edge_buf8_recv_nxm_nyp_3D(bcount))
            allocate(sync_edge_buf8_recv_nxm_nym_3D(bcount))
            do k = 1, bcount
                ny_dir_size = bnx_end(k) - bnx_start(k) + 1
                nx_dir_size = bny_end(k) - bny_start(k) + 1
                allocate(sync_buf8_send_nyp_3D(k)%field(ny_dir_size, nzmax), sync_buf8_recv_nyp_3D(k)%field(ny_dir_size, nzmax))
                allocate(sync_buf8_send_nxp_3D(k)%field(nx_dir_size, nzmax), sync_buf8_recv_nxp_3D(k)%field(nx_dir_size, nzmax))
                allocate(sync_buf8_send_nym_3D(k)%field(ny_dir_size, nzmax), sync_buf8_recv_nym_3D(k)%field(ny_dir_size, nzmax))
                allocate(sync_buf8_send_nxm_3D(k)%field(nx_dir_size, nzmax), sync_buf8_recv_nxm_3D(k)%field(nx_dir_size, nzmax))
                allocate(sync_edge_buf8_recv_nxp_nyp_3D(k)%field(nzmax))
                allocate(sync_edge_buf8_recv_nxp_nym_3D(k)%field(nzmax))
                allocate(sync_edge_buf8_recv_nxm_nyp_3D(k)%field(nzmax))
                allocate(sync_edge_buf8_recv_nxm_nym_3D(k)%field(nzmax))
            enddo
            ! real4
            allocate(sync_buf4_send_nyp_3D(bcount), sync_buf4_recv_nyp_3D(bcount))
            allocate(sync_buf4_send_nxp_3D(bcount), sync_buf4_recv_nxp_3D(bcount))
            allocate(sync_buf4_send_nym_3D(bcount), sync_buf4_recv_nym_3D(bcount))
            allocate(sync_buf4_send_nxm_3D(bcount), sync_buf4_recv_nxm_3D(bcount))
            allocate(sync_edge_buf4_recv_nxp_nyp_3D(bcount))
            allocate(sync_edge_buf4_recv_nxp_nym_3D(bcount))
            allocate(sync_edge_buf4_recv_nxm_nyp_3D(bcount))
            allocate(sync_edge_buf4_recv_nxm_nym_3D(bcount))
            do k = 1, bcount
                ny_dir_size = bnx_end(k) - bnx_start(k) + 1
                nx_dir_size = bny_end(k) - bny_start(k) + 1
                allocate(sync_buf4_send_nyp_3D(k)%field(ny_dir_size, nzmax), sync_buf4_recv_nyp_3D(k)%field(ny_dir_size, nzmax))
                allocate(sync_buf4_send_nxp_3D(k)%field(nx_dir_size, nzmax), sync_buf4_recv_nxp_3D(k)%field(nx_dir_size, nzmax))
                allocate(sync_buf4_send_nym_3D(k)%field(ny_dir_size, nzmax), sync_buf4_recv_nym_3D(k)%field(ny_dir_size, nzmax))
                allocate(sync_buf4_send_nxm_3D(k)%field(nx_dir_size, nzmax), sync_buf4_recv_nxm_3D(k)%field(nx_dir_size, nzmax))
                allocate(sync_edge_buf4_recv_nxp_nyp_3D(k)%field(nzmax))
                allocate(sync_edge_buf4_recv_nxp_nym_3D(k)%field(nzmax))
                allocate(sync_edge_buf4_recv_nxm_nyp_3D(k)%field(nzmax))
                allocate(sync_edge_buf4_recv_nxm_nym_3D(k)%field(nzmax))
            enddo
            
            allocate(reqsts(2*8*bcount), statuses(MPI_STATUS_SIZE, 2*8*bcount))
        end associate
    end subroutine

    subroutine deallocate_mpp_sync_buffers(domain)
        implicit none
        type(domain_type), intent(in) :: domain
        integer :: k

        deallocate(reqsts)
        deallocate(statuses)

        associate(bcount => domain%bcount)
            ! MPI buffers for 2D data
            ! real8
            do k = 1, bcount
                deallocate(sync_buf8_send_nyp(k)%field, sync_buf8_recv_nyp(k)%field)
                deallocate(sync_buf8_send_nxp(k)%field, sync_buf8_recv_nxp(k)%field)
                deallocate(sync_buf8_send_nym(k)%field, sync_buf8_recv_nym(k)%field)
                deallocate(sync_buf8_send_nxm(k)%field, sync_buf8_recv_nxm(k)%field)
            enddo
            deallocate(sync_buf8_send_nyp, sync_buf8_recv_nyp)
            deallocate(sync_buf8_send_nxp, sync_buf8_recv_nxp)
            deallocate(sync_buf8_send_nym, sync_buf8_recv_nym)
            deallocate(sync_buf8_send_nxm, sync_buf8_recv_nxm)
            deallocate(sync_edge_buf8_recv_nxp_nyp)
            deallocate(sync_edge_buf8_recv_nxp_nym)
            deallocate(sync_edge_buf8_recv_nxm_nyp)
            deallocate(sync_edge_buf8_recv_nxm_nym)
            ! real4
            do k = 1, bcount
                deallocate(sync_buf4_send_nyp(k)%field, sync_buf4_recv_nyp(k)%field)
                deallocate(sync_buf4_send_nxp(k)%field, sync_buf4_recv_nxp(k)%field)
                deallocate(sync_buf4_send_nym(k)%field, sync_buf4_recv_nym(k)%field)
                deallocate(sync_buf4_send_nxm(k)%field, sync_buf4_recv_nxm(k)%field)
            enddo
            deallocate(sync_buf4_send_nyp, sync_buf4_recv_nyp)
            deallocate(sync_buf4_send_nxp, sync_buf4_recv_nxp)
            deallocate(sync_buf4_send_nym, sync_buf4_recv_nym)
            deallocate(sync_buf4_send_nxm, sync_buf4_recv_nxm)
            deallocate(sync_edge_buf4_recv_nxp_nyp)
            deallocate(sync_edge_buf4_recv_nxp_nym)
            deallocate(sync_edge_buf4_recv_nxm_nyp)
            deallocate(sync_edge_buf4_recv_nxm_nym)

            ! MPI buffers for 3D data
            ! real8
            do k = 1, bcount
                deallocate(sync_buf8_send_nyp_3D(k)%field, sync_buf8_recv_nyp_3D(k)%field)
                deallocate(sync_buf8_send_nxp_3D(k)%field, sync_buf8_recv_nxp_3D(k)%field)
                deallocate(sync_buf8_send_nym_3D(k)%field, sync_buf8_recv_nym_3D(k)%field)
                deallocate(sync_buf8_send_nxm_3D(k)%field, sync_buf8_recv_nxm_3D(k)%field)
                deallocate(sync_edge_buf8_recv_nxp_nyp_3D(k)%field)
                deallocate(sync_edge_buf8_recv_nxp_nym_3D(k)%field)
                deallocate(sync_edge_buf8_recv_nxm_nyp_3D(k)%field)
                deallocate(sync_edge_buf8_recv_nxm_nym_3D(k)%field)
            enddo
            deallocate(sync_buf8_send_nyp_3D, sync_buf8_recv_nyp_3D)
            deallocate(sync_buf8_send_nxp_3D, sync_buf8_recv_nxp_3D)
            deallocate(sync_buf8_send_nym_3D, sync_buf8_recv_nym_3D)
            deallocate(sync_buf8_send_nxm_3D, sync_buf8_recv_nxm_3D)
            deallocate(sync_edge_buf8_recv_nxp_nyp_3D)
            deallocate(sync_edge_buf8_recv_nxp_nym_3D)
            deallocate(sync_edge_buf8_recv_nxm_nyp_3D)
            deallocate(sync_edge_buf8_recv_nxm_nym_3D)
            ! real4
            do k = 1, bcount
                deallocate(sync_buf4_send_nyp_3D(k)%field, sync_buf4_recv_nyp_3D(k)%field)
                deallocate(sync_buf4_send_nxp_3D(k)%field, sync_buf4_recv_nxp_3D(k)%field)
                deallocate(sync_buf4_send_nym_3D(k)%field, sync_buf4_recv_nym_3D(k)%field)
                deallocate(sync_buf4_send_nxm_3D(k)%field, sync_buf4_recv_nxm_3D(k)%field)
                deallocate(sync_edge_buf4_recv_nxp_nyp_3D(k)%field)
                deallocate(sync_edge_buf4_recv_nxp_nym_3D(k)%field)
                deallocate(sync_edge_buf4_recv_nxm_nyp_3D(k)%field)
                deallocate(sync_edge_buf4_recv_nxm_nym_3D(k)%field)
            enddo
            deallocate(sync_buf4_send_nyp_3D, sync_buf4_recv_nyp_3D)
            deallocate(sync_buf4_send_nxp_3D, sync_buf4_recv_nxp_3D)
            deallocate(sync_buf4_send_nym_3D, sync_buf4_recv_nym_3D)
            deallocate(sync_buf4_send_nxm_3D, sync_buf4_recv_nxm_3D)
            deallocate(sync_edge_buf4_recv_nxp_nyp_3D)
            deallocate(sync_edge_buf4_recv_nxp_nym_3D)
            deallocate(sync_edge_buf4_recv_nxm_nyp_3D)
            deallocate(sync_edge_buf4_recv_nxm_nym_3D)
        end associate
    end subroutine

    integer function get_local_block_number(bindx, bcount, m, n)
        implicit none
        integer, pointer, intent(in) :: bindx(:, :)
        integer, intent(in) :: bcount
        integer, intent(in) :: m, n
        integer :: k
        get_local_block_number = -1
        do k = 1, bcount
            if (m == bindx(k, 1) .and. n == bindx(k, 2)) then
                get_local_block_number = k
                exit
            endif
        enddo
        return
    end function

    subroutine check_block_status(bglob_proc, bnx, bny, b_coords, p)
        implicit none
        integer, pointer, intent(in) :: bglob_proc(:, :)
        integer, intent(in) :: bnx, bny
        integer, dimension(2), intent(in) :: b_coords
        integer, intent(out) :: p
        integer, dimension(2) :: bgrid
        bgrid(1) = bnx; bgrid(2) = bny
        if (check_cart_coord(b_coords - 1, bgrid) /= 1) then
            ! Out of range, block does not exist
            p = -2
            return
        endif
        p = bglob_proc(b_coords(1), b_coords(2))
        return
    end subroutine

    integer function check_cart_coord(coord, grid_size)
        implicit none
        integer, dimension(2), intent(in) :: coord, grid_size
        check_cart_coord = 0
        if (all(coord .ge. 0) .and. all((grid_size - coord) .ge. 1)) then
            check_cart_coord = 1
        endif
        return
    end function

! ------------------------------------------------------------------------------- !
! ------------------------------------------------------------------------------- !
!                            Sync 2D data                                         !
! ------------------------------------------------------------------------------- !
! ------------------------------------------------------------------------------- !
    subroutine irecv_real8(k, src_block, src_proc, buff, buff_size, tag, reqst)
        implicit none
        integer, dimension(2) :: src_block
        integer :: k, src_proc, tag
        integer :: buff_size
        real(wp8) :: buff(:)
        integer :: reqst
        integer :: ierr

        !if (parallel_dbg >= 4) print *, rank, 'IRECV. block: ', bindx(k, 1), bindx(k, 2), 'src_block:', src_block(1), src_block(2), 'src_p:', src_proc, 'tag', tag
        if (parallel_dbg >= 4) print *, mpp_rank, 'IRECV. block: ', k, 'src_block:', src_block(1), src_block(2), 'src_p:', src_proc, 'tag', tag

        call mpi_irecv(buff, buff_size, mpi_real8, src_proc, tag, mpp_cart_comm, reqst, ierr)
    end subroutine

    subroutine isend_real8(k, dist_block, dist_proc, buff, buff_size, tag, reqst)
        implicit none
        integer, dimension(2) :: dist_block
        integer :: k, dist_proc, tag
        integer :: buff_size
        real(wp8) :: buff(:)
        integer :: reqst
        integer :: ierr

        !if (parallel_dbg >= 4) print *, rank, 'ISEND. block: ', bindx(k, 1), bindx(k, 2), 'dst_block:', dist_block(1), dist_block(2), 'dst_p:', dist_proc, 'tag', tag
        if (parallel_dbg >= 4) print *, mpp_rank, 'ISEND. block: ', k, 'dst_block:', dist_block(1), dist_block(2), 'dst_p:', dist_proc, 'tag', tag

        call mpi_isend(buff, buff_size, mpi_real8, dist_proc, tag, mpp_cart_comm, reqst, ierr)
    end subroutine

    subroutine syncborder_data2D_real8(domain, data2d)
#define _MPI_TYPE_ mpi_real8
#define _IRECV_ irecv_real8
#define _ISEND_ isend_real8
        
#define _SYNC_BUF_SEND_NYP_ sync_buf8_send_nyp
#define _SYNC_BUF_SEND_NXP_ sync_buf8_send_nxp
#define _SYNC_BUF_SEND_NYM_ sync_buf8_send_nym
#define _SYNC_BUF_SEND_NXM_ sync_buf8_send_nxm

#define _SYNC_BUF_RECV_NYP_ sync_buf8_recv_nyp
#define _SYNC_BUF_RECV_NXP_ sync_buf8_recv_nxp
#define _SYNC_BUF_RECV_NYM_ sync_buf8_recv_nym
#define _SYNC_BUF_RECV_NXM_ sync_buf8_recv_nxm

#define _SYNC_EDGE_BUF_RECV_NXP_NYP_ sync_edge_buf8_recv_nxp_nyp
#define _SYNC_EDGE_BUF_RECV_NXP_NYM_ sync_edge_buf8_recv_nxp_nym
#define _SYNC_EDGE_BUF_RECV_NXM_NYP_ sync_edge_buf8_recv_nxm_nyp
#define _SYNC_EDGE_BUF_RECV_NXM_NYM_ sync_edge_buf8_recv_nxm_nym
        
        implicit none
        type(data2D_real8_type), intent(inout) :: data2d
        type(domain_type), intent(in) :: domain

#include    "syncborder_block2D_gen.fi"

#undef _MPI_TYPE_
#undef _IRECV_
#undef _ISEND_
#undef _SYNC_BUF_SEND_NYP_
#undef _SYNC_BUF_SEND_NXP_
#undef _SYNC_BUF_SEND_NYM_
#undef _SYNC_BUF_SEND_NXM_
#undef _SYNC_BUF_RECV_NYP_
#undef _SYNC_BUF_RECV_NXP_
#undef _SYNC_BUF_RECV_NYM_
#undef _SYNC_BUF_RECV_NXM_
#undef _SYNC_EDGE_BUF_RECV_NXP_NYP_
#undef _SYNC_EDGE_BUF_RECV_NXP_NYM_
#undef _SYNC_EDGE_BUF_RECV_NXM_NYP_
#undef _SYNC_EDGE_BUF_RECV_NXM_NYM_
    end subroutine

    subroutine irecv_real4(k, src_block, src_proc, buff, buff_size, tag, reqst)
        implicit none
        integer, dimension(2) :: src_block
        integer :: k, src_proc, tag
        integer :: buff_size
        real(wp4) :: buff(:)
        integer :: reqst
        integer :: ierr

        !if (parallel_dbg >= 4) print *, rank, 'IRECV. block: ', bindx(k, 1), bindx(k, 2), 'src_block:', src_block(1), src_block(2), 'src_p:', src_proc, 'tag', tag
        if (parallel_dbg >= 4) print *, mpp_rank, 'IRECV. block: ', k, 'src_block:', src_block(1), src_block(2), 'src_p:', src_proc, 'tag', tag

        call mpi_irecv(buff, buff_size, mpi_real4, src_proc, tag, mpp_cart_comm, reqst, ierr)
    end subroutine

    subroutine isend_real4(k, dist_block, dist_proc, buff, buff_size, tag, reqst)
        implicit none
        integer, dimension(2) :: dist_block
        integer :: k, dist_proc, tag
        integer :: buff_size
        real(wp4) :: buff(:)
        integer :: reqst
        integer :: ierr

        !if (parallel_dbg >= 4) print *, rank, 'ISEND. block: ', bindx(k, 1), bindx(k, 2), 'dst_block:', dist_block(1), dist_block(2), 'dst_p:', dist_proc, 'tag', tag
        if (parallel_dbg >= 4) print *, mpp_rank, 'ISEND. block: ', k, 'dst_block:', dist_block(1), dist_block(2), 'dst_p:', dist_proc, 'tag', tag

        call mpi_isend(buff, buff_size, mpi_real4, dist_proc, tag, mpp_cart_comm, reqst, ierr)
    end subroutine

    subroutine syncborder_data2D_real4(domain, data2d)
#define _MPI_TYPE_ mpi_real4
#define _IRECV_ irecv_real4
#define _ISEND_ isend_real4
        
#define _SYNC_BUF_SEND_NYP_ sync_buf4_send_nyp
#define _SYNC_BUF_SEND_NXP_ sync_buf4_send_nxp
#define _SYNC_BUF_SEND_NYM_ sync_buf4_send_nym
#define _SYNC_BUF_SEND_NXM_ sync_buf4_send_nxm

#define _SYNC_BUF_RECV_NYP_ sync_buf4_recv_nyp
#define _SYNC_BUF_RECV_NXP_ sync_buf4_recv_nxp
#define _SYNC_BUF_RECV_NYM_ sync_buf4_recv_nym
#define _SYNC_BUF_RECV_NXM_ sync_buf4_recv_nxm

#define _SYNC_EDGE_BUF_RECV_NXP_NYP_ sync_edge_buf4_recv_nxp_nyp
#define _SYNC_EDGE_BUF_RECV_NXP_NYM_ sync_edge_buf4_recv_nxp_nym
#define _SYNC_EDGE_BUF_RECV_NXM_NYP_ sync_edge_buf4_recv_nxm_nyp
#define _SYNC_EDGE_BUF_RECV_NXM_NYM_ sync_edge_buf4_recv_nxm_nym
        
        implicit none
        type(data2D_real4_type), intent(inout) :: data2d
        type(domain_type), intent(in) :: domain

#include    "syncborder_block2D_gen.fi"

#undef _MPI_TYPE_
#undef _IRECV_
#undef _ISEND_
#undef _SYNC_BUF_SEND_NYP_
#undef _SYNC_BUF_SEND_NXP_
#undef _SYNC_BUF_SEND_NYM_
#undef _SYNC_BUF_SEND_NXM_
#undef _SYNC_BUF_RECV_NYP_
#undef _SYNC_BUF_RECV_NXP_
#undef _SYNC_BUF_RECV_NYM_
#undef _SYNC_BUF_RECV_NXM_
#undef _SYNC_EDGE_BUF_RECV_NXP_NYP_
#undef _SYNC_EDGE_BUF_RECV_NXP_NYM_
#undef _SYNC_EDGE_BUF_RECV_NXM_NYP_
#undef _SYNC_EDGE_BUF_RECV_NXM_NYM_
    end subroutine

! ------------------------------------------------------------------------------- !
! ------------------------------------------------------------------------------- !
!                            Sync 3D data                                         !
! ------------------------------------------------------------------------------- !
! ------------------------------------------------------------------------------- !
    subroutine irecv3D_real8(k, src_block, src_proc, buff, buff_size, tag, reqst)
        implicit none
        integer, dimension(2) :: src_block
        integer :: k, src_proc, tag
        integer :: buff_size
        real(wp8) :: buff(:, :)
        integer :: reqst
        integer :: ierr

        !if (parallel_dbg >= 4) print *, rank, 'IRECV. block: ', bindx(k, 1), bindx(k, 2), 'src_block:', src_block(1), src_block(2), 'src_p:', src_proc, 'tag', tag
        if (parallel_dbg >= 4) print *, mpp_rank, 'IRECV. block: ', k, 'src_block:', src_block(1), src_block(2), 'src_p:', src_proc, 'tag', tag

        call mpi_irecv(buff, buff_size, mpi_real8, src_proc, tag, mpp_cart_comm, reqst, ierr)
    end subroutine

    subroutine isend3D_real8(k, dist_block, dist_proc, buff, buff_size, tag, reqst)
        implicit none
        integer, dimension(2) :: dist_block
        integer :: k, dist_proc, tag
        integer :: buff_size
        real(wp8) :: buff(:, :)
        integer :: reqst
        integer :: ierr

        !if (parallel_dbg >= 4) print *, rank, 'ISEND. block: ', bindx(k, 1), bindx(k, 2), 'dst_block:', dist_block(1), dist_block(2), 'dst_p:', dist_proc, 'tag', tag
        if (parallel_dbg >= 4) print *, mpp_rank, 'ISEND. block: ', k, 'dst_block:', dist_block(1), dist_block(2), 'dst_p:', dist_proc, 'tag', tag

        call mpi_isend(buff, buff_size, mpi_real8, dist_proc, tag, mpp_cart_comm, reqst, ierr)
    end subroutine

    subroutine syncborder_data3D_real8(domain, data3d, nz)
#define _MPI_TYPE_ mpi_real8
#define _IRECV_ irecv3D_real8
#define _ISEND_ isend3D_real8
            
#define _SYNC_BUF_SEND_NYP_ sync_buf8_send_nyp_3D
#define _SYNC_BUF_SEND_NXP_ sync_buf8_send_nxp_3D
#define _SYNC_BUF_SEND_NYM_ sync_buf8_send_nym_3D
#define _SYNC_BUF_SEND_NXM_ sync_buf8_send_nxm_3D

#define _SYNC_BUF_RECV_NYP_ sync_buf8_recv_nyp_3D
#define _SYNC_BUF_RECV_NXP_ sync_buf8_recv_nxp_3D
#define _SYNC_BUF_RECV_NYM_ sync_buf8_recv_nym_3D
#define _SYNC_BUF_RECV_NXM_ sync_buf8_recv_nxm_3D

#define _SYNC_EDGE_BUF_RECV_NXP_NYP_ sync_edge_buf8_recv_nxp_nyp_3D
#define _SYNC_EDGE_BUF_RECV_NXP_NYM_ sync_edge_buf8_recv_nxp_nym_3D
#define _SYNC_EDGE_BUF_RECV_NXM_NYP_ sync_edge_buf8_recv_nxm_nyp_3D
#define _SYNC_EDGE_BUF_RECV_NXM_NYM_ sync_edge_buf8_recv_nxm_nym_3D
            
        implicit none
        type(data3D_real8_type), intent(inout) :: data3d
        type(domain_type), intent(in) :: domain
        integer, intent(in) :: nz

#include    "syncborder_block3D_gen.fi"

#undef _MPI_TYPE_
#undef _IRECV_
#undef _ISEND_
#undef _SYNC_BUF_SEND_NYP_
#undef _SYNC_BUF_SEND_NXP_
#undef _SYNC_BUF_SEND_NYM_
#undef _SYNC_BUF_SEND_NXM_
#undef _SYNC_BUF_RECV_NYP_
#undef _SYNC_BUF_RECV_NXP_
#undef _SYNC_BUF_RECV_NYM_
#undef _SYNC_BUF_RECV_NXM_
#undef _SYNC_EDGE_BUF_RECV_NXP_NYP_
#undef _SYNC_EDGE_BUF_RECV_NXP_NYM_
#undef _SYNC_EDGE_BUF_RECV_NXM_NYP_
#undef _SYNC_EDGE_BUF_RECV_NXM_NYM_
    end subroutine

    subroutine irecv3D_real4(k, src_block, src_proc, buff, buff_size, tag, reqst)
        implicit none
        integer, dimension(2) :: src_block
        integer :: k, src_proc, tag
        integer :: buff_size
        real(wp4) :: buff(:, :)
        integer :: reqst
        integer :: ierr

        !if (parallel_dbg >= 4) print *, rank, 'IRECV. block: ', bindx(k, 1), bindx(k, 2), 'src_block:', src_block(1), src_block(2), 'src_p:', src_proc, 'tag', tag
        if (parallel_dbg >= 4) print *, mpp_rank, 'IRECV. block: ', k, 'src_block:', src_block(1), src_block(2), 'src_p:', src_proc, 'tag', tag

        call mpi_irecv(buff, buff_size, mpi_real4, src_proc, tag, mpp_cart_comm, reqst, ierr)
    end subroutine

    subroutine isend3D_real4(k, dist_block, dist_proc, buff, buff_size, tag, reqst)
        implicit none
        integer, dimension(2) :: dist_block
        integer :: k, dist_proc, tag
        integer :: buff_size
        real(wp4) :: buff(:, :)
        integer :: reqst
        integer :: ierr

        !if (parallel_dbg >= 4) print *, rank, 'ISEND. block: ', bindx(k, 1), bindx(k, 2), 'dst_block:', dist_block(1), dist_block(2), 'dst_p:', dist_proc, 'tag', tag
        if (parallel_dbg >= 4) print *, mpp_rank, 'ISEND. block: ', k, 'dst_block:', dist_block(1), dist_block(2), 'dst_p:', dist_proc, 'tag', tag

        call mpi_isend(buff, buff_size, mpi_real4, dist_proc, tag, mpp_cart_comm, reqst, ierr)
    end subroutine

    subroutine syncborder_data3D_real4(domain, data3d, nz)
#define _MPI_TYPE_ mpi_real4
#define _IRECV_ irecv3D_real4
#define _ISEND_ isend3D_real4
        
#define _SYNC_BUF_SEND_NYP_ sync_buf4_send_nyp_3D
#define _SYNC_BUF_SEND_NXP_ sync_buf4_send_nxp_3D
#define _SYNC_BUF_SEND_NYM_ sync_buf4_send_nym_3D
#define _SYNC_BUF_SEND_NXM_ sync_buf4_send_nxm_3D

#define _SYNC_BUF_RECV_NYP_ sync_buf4_recv_nyp_3D
#define _SYNC_BUF_RECV_NXP_ sync_buf4_recv_nxp_3D
#define _SYNC_BUF_RECV_NYM_ sync_buf4_recv_nym_3D
#define _SYNC_BUF_RECV_NXM_ sync_buf4_recv_nxm_3D

#define _SYNC_EDGE_BUF_RECV_NXP_NYP_ sync_edge_buf4_recv_nxp_nyp_3D
#define _SYNC_EDGE_BUF_RECV_NXP_NYM_ sync_edge_buf4_recv_nxp_nym_3D
#define _SYNC_EDGE_BUF_RECV_NXM_NYP_ sync_edge_buf4_recv_nxm_nyp_3D
#define _SYNC_EDGE_BUF_RECV_NXM_NYM_ sync_edge_buf4_recv_nxm_nym_3D
        
        implicit none
        type(data3D_real4_type), intent(inout) :: data3d
        type(domain_type), intent(in) :: domain
        integer, intent(in) :: nz

#include    "syncborder_block3D_gen.fi"

#undef _MPI_TYPE_
#undef _IRECV_
#undef _ISEND_
#undef _SYNC_BUF_SEND_NYP_
#undef _SYNC_BUF_SEND_NXP_
#undef _SYNC_BUF_SEND_NYM_
#undef _SYNC_BUF_SEND_NXM_
#undef _SYNC_BUF_RECV_NYP_
#undef _SYNC_BUF_RECV_NXP_
#undef _SYNC_BUF_RECV_NYM_
#undef _SYNC_BUF_RECV_NXM_
#undef _SYNC_EDGE_BUF_RECV_NXP_NYP_
#undef _SYNC_EDGE_BUF_RECV_NXP_NYM_
#undef _SYNC_EDGE_BUF_RECV_NXM_NYP_
#undef _SYNC_EDGE_BUF_RECV_NXM_NYM_
    end subroutine

end module mpp_sync_module