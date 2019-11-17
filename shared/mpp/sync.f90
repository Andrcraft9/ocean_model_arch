module mpp_sync_module
    ! Module description

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use mpp_module, only: mpp_rank, mpp_count, mpp_cart_comm, mpp_size, mpp_coord
    use errors_module, only: abort_model, check_error
    use data_types_module
    use decomposition_module, only: domain_type

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

    public :: sync
    private :: allocate_mpp_sync_buffers
    private :: deallocate_mpp_sync_buffers

contains

    subroutine allocate_mpp_sync_buffers(domain)
        implicit none
        type(domain_type), intent(in) :: domain
        integer :: k, nx_dir_size, ny_dir_size

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
        type(domain), intent(in) :: domain
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

    integer function get_local_block_number(m, n)
        implicit none
        integer :: m, n
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

    subroutine check_block_status(b_coords, p)
        implicit none
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
        !write(*,*) coord,all(coord.ge.0),all((p_size-coord).ge.1)
        !print *, coord, p_size - coord, all((p_size-coord).ge.1)
        if (all(coord .ge. 0) .and. all((grid_size - coord) .ge. 1)) then
            check_cart_coord = 1
        endif
        return
    end function

    subroutine get_block_and_rank_by_point(m, n, out)
        implicit none
        integer :: m, n
        integer :: out(2)
        integer :: k, flag_r, r, flag_b, b, ierr

        flag_r = -1
        flag_b = -1
        do k = 1, bcount
            if (m >= bnx_start(k) .and. m <= bnx_end(k)) then
                if (n >= bny_start(k) .and. n <= bny_end(k)) then
                    flag_r = rank
                    flag_b = k
                endif
            endif
        enddo

        call mpi_allreduce(flag_r, r, 1, mpi_integer,      &
                            mpi_max, cart_comm, ierr)

        call mpi_allreduce(flag_b, b, 1, mpi_integer,      &
                            mpi_max, cart_comm, ierr)

        out(1) = r
        out(2) = b
    end subroutine

    integer function get_rank_by_point(m, n)
        implicit none
        integer :: m, n
        integer :: flag_r, r, ierr

        flag_r = -1
        if (m >= nx_start .and. m <= nx_end) then
            if (n >= ny_start .and. n <= ny_end) then
                flag_r = rank
            endif
        endif

        call mpi_allreduce(flag_r, r, 1, mpi_integer,      &
                            mpi_max, cart_comm, ierr)

        get_rank_by_point = r
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
        real*8 :: buff(:)
        integer :: reqst
        integer :: ierr

        if (parallel_dbg >= 4) print *, rank, 'IRECV. block: ', bindx(k, 1), bindx(k, 2), 'src_block:', src_block(1), src_block(2), 'src_p:', src_proc, 'tag', tag
        !call mpi_irecv(blks(k)%vals(dx1:dx2, dy1:dy2), (dx2 - dx1 + 1)*(dy2 - dy1 + 1), &
        !               mpi_real8, p_src, tag, cart_comm, reqst, ierr)
        call mpi_irecv(buff, buff_size, mpi_real8, src_proc, tag, cart_comm, reqst, ierr)
    end subroutine

    subroutine isend_real8(k, dist_block, dist_proc, buff, buff_size, tag, reqst)
        implicit none
        integer, dimension(2) :: dist_block
        integer :: k, dist_proc, tag
        integer :: buff_size
        real*8 :: buff(:)
        integer :: reqst
        integer :: ierr

        if (parallel_dbg >= 4) print *, rank, 'ISEND. block: ', bindx(k, 1), bindx(k, 2), 'dst_block:', dist_block(1), dist_block(2), 'dst_p:', dist_proc, 'tag', tag
        !call mpi_isend(blks(k)%vals(sx1:sx2, sy1:sy2), (sx2 - sx1 + 1)*(sy2 - sy1 + 1), &
        !               mpi_real8, p_dist, tag, cart_comm, reqst, ierr)
        call mpi_isend(buff, buff_size, mpi_real8, dist_proc, tag, cart_comm, reqst, ierr)
    end subroutine

    subroutine syncborder_block2D_real8(blks)
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
        type(block2D_real8), dimension(:), pointer :: blks

#include "syncborder_block2D_gen.fi"

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
        real*4 :: buff(:)
        integer :: reqst
        integer :: ierr

        if (parallel_dbg >= 4) print *, rank, 'IRECV. block: ', bindx(k, 1), bindx(k, 2), 'src_block:', src_block(1), src_block(2), 'src_p:', src_proc, 'tag', tag
        !call mpi_irecv(blks(k)%vals(dx1:dx2, dy1:dy2), (dx2 - dx1 + 1)*(dy2 - dy1 + 1), &
        !               mpi_real8, p_src, tag, cart_comm, reqst, ierr)
        call mpi_irecv(buff, buff_size, mpi_real4, src_proc, tag, cart_comm, reqst, ierr)
    end subroutine

    subroutine isend_real4(k, dist_block, dist_proc, buff, buff_size, tag, reqst)
        implicit none
        integer, dimension(2) :: dist_block
        integer :: k, dist_proc, tag
        integer :: buff_size
        real*4 :: buff(:)
        integer :: reqst
        integer :: ierr

        if (parallel_dbg >= 4) print *, rank, 'ISEND. block: ', bindx(k, 1), bindx(k, 2), 'dst_block:', dist_block(1), dist_block(2), 'dst_p:', dist_proc, 'tag', tag
        !call mpi_isend(blks(k)%vals(sx1:sx2, sy1:sy2), (sx2 - sx1 + 1)*(sy2 - sy1 + 1), &
        !               mpi_real8, p_dist, tag, cart_comm, reqst, ierr)
        call mpi_isend(buff, buff_size, mpi_real4, dist_proc, tag, cart_comm, reqst, ierr)
    end subroutine

    subroutine syncborder_block2D_real4(blks)
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
        type(block2D_real4), dimension(:), pointer :: blks
        
#include "syncborder_block2D_gen.fi"

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
        real*8 :: buff(:, :)
        integer :: reqst
        integer :: ierr

        if (parallel_dbg >= 4) print *, rank, 'IRECV. block: ', bindx(k, 1), bindx(k, 2), 'src_block:', src_block(1), src_block(2), 'src_p:', src_proc, 'tag', tag
        !call mpi_irecv(blks(k)%vals(dx1:dx2, dy1:dy2), (dx2 - dx1 + 1)*(dy2 - dy1 + 1), &
        !               mpi_real8, p_src, tag, cart_comm, reqst, ierr)
        call mpi_irecv(buff, buff_size, mpi_real8, src_proc, tag, cart_comm, reqst, ierr)
    end subroutine

    subroutine isend3D_real8(k, dist_block, dist_proc, buff, buff_size, tag, reqst)
        implicit none
        integer, dimension(2) :: dist_block
        integer :: k, dist_proc, tag
        integer :: buff_size
        real*8 :: buff(:, :)
        integer :: reqst
        integer :: ierr

        if (parallel_dbg >= 4) print *, rank, 'ISEND. block: ', bindx(k, 1), bindx(k, 2), 'dst_block:', dist_block(1), dist_block(2), 'dst_p:', dist_proc, 'tag', tag
        !call mpi_isend(blks(k)%vals(sx1:sx2, sy1:sy2), (sx2 - sx1 + 1)*(sy2 - sy1 + 1), &
        !               mpi_real8, p_dist, tag, cart_comm, reqst, ierr)
        call mpi_isend(buff, buff_size, mpi_real8, dist_proc, tag, cart_comm, reqst, ierr)
    end subroutine

    subroutine syncborder_block3D_real8(blks, nz)
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
        type(block3D_real8), dimension(:), pointer :: blks
        integer :: nz

#include "syncborder_block3D_gen.fi"

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
        real*4 :: buff(:, :)
        integer :: reqst
        integer :: ierr

        if (parallel_dbg >= 4) print *, rank, 'IRECV. block: ', bindx(k, 1), bindx(k, 2), 'src_block:', src_block(1), src_block(2), 'src_p:', src_proc, 'tag', tag
        !call mpi_irecv(blks(k)%vals(dx1:dx2, dy1:dy2), (dx2 - dx1 + 1)*(dy2 - dy1 + 1), &
        !               mpi_real8, p_src, tag, cart_comm, reqst, ierr)
        call mpi_irecv(buff, buff_size, mpi_real4, src_proc, tag, cart_comm, reqst, ierr)
    end subroutine

    subroutine isend3D_real4(k, dist_block, dist_proc, buff, buff_size, tag, reqst)
        implicit none
        integer, dimension(2) :: dist_block
        integer :: k, dist_proc, tag
        integer :: buff_size
        real*4 :: buff(:, :)
        integer :: reqst
        integer :: ierr

        if (parallel_dbg >= 4) print *, rank, 'ISEND. block: ', bindx(k, 1), bindx(k, 2), 'dst_block:', dist_block(1), dist_block(2), 'dst_p:', dist_proc, 'tag', tag
        !call mpi_isend(blks(k)%vals(sx1:sx2, sy1:sy2), (sx2 - sx1 + 1)*(sy2 - sy1 + 1), &
        !               mpi_real8, p_dist, tag, cart_comm, reqst, ierr)
        call mpi_isend(buff, buff_size, mpi_real4, dist_proc, tag, cart_comm, reqst, ierr)
    end subroutine

    subroutine syncborder_block3D_real4(blks, nz)
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
        type(block3D_real4), dimension(:), pointer :: blks
        integer :: nz
        
#include "syncborder_block3D_gen.fi"

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