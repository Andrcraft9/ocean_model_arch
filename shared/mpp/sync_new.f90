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
    use decomposition_module, only: domain_type

#include "macros/mpp_macros.fi"


    implicit none
    save
    private

    interface sync
        module procedure syncborder_data2D_real8
        module procedure syncborder_data2D_real4
        module procedure syncborder_data3D_real8
        module procedure syncborder_data3D_real4
    end interface

    public :: mpp_sync_init
    public :: mpp_sync_finalize
    public :: sync

!------------------------------------------------------------------------------

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

end module mpp_sync_module