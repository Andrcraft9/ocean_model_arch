module mpp_sync_module
    ! Sync data of different types, depends on domain

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    !use mpi
    use debug_module
    use mpp_module, only: mpp_rank, mpp_count, mpp_cart_comm, mpp_size, mpp_coord, mpp_time_sync, start_timer, end_timer
    use errors_module, only: abort_model, check_error
    use data_types_module, only: block1D_real4_type, block1D_real8_type, block2D_real4_type, block2D_real8_type,  &
                                 data1D_real4_type, data1D_real8_type, data2D_real4_type, data2D_real8_type,      &
                                 data3D_real4_type, data3D_real8_type
    use decomposition_module, only: domain_type, get_boundary_points_of_block, get_halo_points_of_block,             &
                                    get_inverse_dir, get_rank_of_block, get_loc_num_of_block, fill_direction_array,  &
                                    is_block_on_current_proc, is_inner_block

#include "macros/mpp_macros.fi"
#include "macros/kernel_macros.fi"

    implicit none
    save
    private

    interface sync
        module procedure syncborder_data2D_real8
        module procedure syncborder_data2D_real4
        !module procedure syncborder_data3D_real8
        !module procedure syncborder_data3D_real4
    end interface

    public :: mpp_sync_init
    public :: mpp_sync_finalize
    public :: sync

    ! Buffers for each near rank
    real(wp8), allocatable :: sync_send_buf_r8(:, :)
    real(wp4), allocatable :: sync_send_buf_r4(:, :)
    real(wp8), allocatable :: sync_recv_buf_r8(:, :)
    real(wp4), allocatable :: sync_recv_buf_r4(:, :)
    integer, allocatable :: sync_buf_pos(:)  ! Position of buffers, need for prepare buffers. Position of buffer for each near rank

!------------------------------------------------------------------------------

contains

    subroutine mpp_sync_init(domain)
        type(domain_type), intent(in) :: domain

        call allocate_mpp_sync_buffers(domain)
    end subroutine

    subroutine mpp_sync_finalize(domain)
        type(domain_type), intent(in) :: domain

        call deallocate_mpp_sync_buffers(domain)
    end subroutine

    subroutine allocate_mpp_sync_buffers(domain)
        type(domain_type), intent(in) :: domain

        integer :: k, kk
        integer :: sync_dir(8)
        integer :: k_dir(8)
        integer :: rank_dir(8)
        integer :: nxs, nxe, nys, nye     ! Boundary points
        integer, allocatable :: buf_sizes(:)
        integer :: max_buf_size = 1

        allocate(buf_sizes(mpp_count))
        buf_sizes = 0

        do k = domain%start_boundary, domain%start_boundary + domain%bcount_boundary - 1

            sync_dir(1) = _NXP_; k_dir(1) = domain%blocks_info(k)%k_nxp; rank_dir(1) = domain%blocks_info(k)%rank_nxp
            sync_dir(2) = _NXM_; k_dir(2) = domain%blocks_info(k)%k_nxm; rank_dir(2) = domain%blocks_info(k)%rank_nxm
            sync_dir(3) = _NYP_; k_dir(3) = domain%blocks_info(k)%k_nyp; rank_dir(3) = domain%blocks_info(k)%rank_nyp
            sync_dir(4) = _NYM_; k_dir(4) = domain%blocks_info(k)%k_nym; rank_dir(4) = domain%blocks_info(k)%rank_nym

            sync_dir(5) = _NXP_NYP_; k_dir(5) = domain%blocks_info(k)%k_nxp_nyp; rank_dir(5) = domain%blocks_info(k)%rank_nxp_nyp
            sync_dir(6) = _NXP_NYM_; k_dir(6) = domain%blocks_info(k)%k_nxp_nym; rank_dir(6) = domain%blocks_info(k)%rank_nxp_nym
            sync_dir(7) = _NXM_NYP_; k_dir(7) = domain%blocks_info(k)%k_nxm_nyp; rank_dir(7) = domain%blocks_info(k)%rank_nxm_nyp
            sync_dir(8) = _NXM_NYM_; k_dir(8) = domain%blocks_info(k)%k_nxm_nym; rank_dir(8) = domain%blocks_info(k)%rank_nxm_nym

            do kk = 1, 8
                if (rank_dir(kk) >= 0) then
                    if (rank_dir(kk) /= mpp_rank) then
                        call get_boundary_points_of_block(domain, k, sync_dir(kk), nxs, nxe, nys, nye)
                        buf_sizes(rank_dir(kk)) = buf_sizes(rank_dir(kk)) + 1 + (nxe - nxs + 1) + (nye - nys + 1)
                    endif
                endif
            enddo
        enddo

        do k = 1, mpp_count
            if (max_buf_size < buf_sizes(k)) max_buf_size = buf_sizes(k)
        enddo
        deallocate(buf_sizes)

        if (debug_level >= 5) then
            print *, "MPP SYNC INFO: rank, max buffer size = ", mpp_rank, max_buf_size
        endif

        allocate(sync_send_buf_r4(max_buf_size, domain%amount_of_ranks_near))
        allocate(sync_send_buf_r8(max_buf_size, domain%amount_of_ranks_near))
        allocate(sync_recv_buf_r4(max_buf_size, domain%amount_of_ranks_near))
        allocate(sync_recv_buf_r8(max_buf_size, domain%amount_of_ranks_near))
        allocate(sync_buf_pos(domain%amount_of_ranks_near))
        sync_send_buf_r4 = 0; sync_send_buf_r8 = 0
        sync_recv_buf_r4 = 0; sync_recv_buf_r8 = 0
        sync_buf_pos = 0

        call mpp_sync_output()
    end subroutine

    subroutine deallocate_mpp_sync_buffers(domain)
        type(domain_type), intent(in) :: domain

        deallocate(sync_send_buf_r4)
        deallocate(sync_send_buf_r8)
        deallocate(sync_recv_buf_r4)
        deallocate(sync_recv_buf_r8)
        deallocate(sync_buf_pos)
    end subroutine

!------------------------------------------------------------------------------
! Generic subroutines:

        subroutine syncborder_data2D_real8(domain, data2d)
#define _MPI_TYPE_ mpi_real8
#define _SYNC_SEND_BUF_ sync_send_buf_r8
#define _SYNC_RECV_BUF_ sync_recv_buf_r8
        
        implicit none
        type(data2D_real8_type), intent(inout) :: data2d
        type(domain_type), intent(in) :: domain

#include    "syncborder_block2D_gen.fi"

#undef _MPI_TYPE_
#undef _SYNC_SEND_BUF_
#undef _SYNC_RECV_BUF_
        end subroutine

        subroutine syncborder_data2D_real4(domain, data2d)
#define _MPI_TYPE_ mpi_real4
#define _SYNC_SEND_BUF_ sync_send_buf_r4
#define _SYNC_RECV_BUF_ sync_recv_buf_r4

        implicit none
        type(data2D_real4_type), intent(inout) :: data2d
        type(domain_type), intent(in) :: domain

#include    "syncborder_block2D_gen.fi"

#undef _MPI_TYPE_
#undef _SYNC_SEND_BUF_
#undef _SYNC_RECV_BUF_
        end subroutine

end module mpp_sync_module