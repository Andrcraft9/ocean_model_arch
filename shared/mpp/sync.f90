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
    real(wp8), allocatable :: sync_buf8(:, :)
    real(wp4), allocatable :: sync_buf4(:, :)

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

    end subroutine

    subroutine deallocate_mpp_sync_buffers(domain)
        type(domain_type), intent(in) :: domain
        
    end subroutine

        subroutine syncborder_data2D_real8(domain, data2d)
#define _MPI_TYPE_ mpi_real8
        
        implicit none
        type(data2D_real8_type), intent(inout) :: data2d
        type(domain_type), intent(in) :: domain

        integer :: k, kk
        integer :: sync_dir(8)
        integer :: k_dir(8)
        integer :: nxs, nxe, nys, nye     ! Boundary points
        integer :: hnxs, hnxe, hnys, hnye ! Halo points

        ! Update halo for inner blocks (only copy)
        !

        do k = domain%start_inner, domain%start_inner + domain%bcount_inner - 1
            
            sync_dir(1) = _NXP_; k_dir(1) = domain%blocks_info(k)%k_nxp
            sync_dir(2) = _NXM_; k_dir(2) = domain%blocks_info(k)%k_nxm
            sync_dir(3) = _NYP_; k_dir(3) = domain%blocks_info(k)%k_nyp
            sync_dir(4) = _NYM_; k_dir(4) = domain%blocks_info(k)%k_nym

            sync_dir(5) = _NXP_NYP_; k_dir(5) = domain%blocks_info(k)%k_nxp_nyp
            sync_dir(6) = _NXP_NYM_; k_dir(6) = domain%blocks_info(k)%k_nxp_nym
            sync_dir(7) = _NXM_NYP_; k_dir(7) = domain%blocks_info(k)%k_nxm_nyp
            sync_dir(8) = _NXM_NYM_; k_dir(8) = domain%blocks_info(k)%k_nxm_nym
            
            do kk = 1, 8
                call get_boundary_points_of_block(domain, k_dir(kk), get_inverse_dir(sync_dir(kk)), nxs, nxe, nys, nye)
                call get_halo_points_of_block(domain, k, sync_dir(kk), hnxs, hnxe, hnys, hnye)
                data2d%block(k)%field(hnxs : hnxe, hnys : hnye) = data2d%block(k_dir(kk))%field(nxs : nxe, nys : nye)
            enddo
        enddo

#undef _MPI_TYPE_
    end subroutine

    subroutine syncborder_data2D_real4(domain, data2d)
#define _MPI_TYPE_ mpi_real4
        
        implicit none
        type(data2D_real4_type), intent(inout) :: data2d
        type(domain_type), intent(in) :: domain

        integer :: k, kk
        integer :: sync_dir(8)
        integer :: k_dir(8)
        integer :: nxs, nxe, nys, nye     ! Boundary points
        integer :: hnxs, hnxe, hnys, hnye ! Halo points

        ! Update halo for inner blocks (only copy)
        !

        do k = domain%start_inner, domain%start_inner + domain%bcount_inner - 1
            
            sync_dir(1) = _NXP_; k_dir(1) = domain%blocks_info(k)%k_nxp
            sync_dir(2) = _NXM_; k_dir(2) = domain%blocks_info(k)%k_nxm
            sync_dir(3) = _NYP_; k_dir(3) = domain%blocks_info(k)%k_nyp
            sync_dir(4) = _NYM_; k_dir(4) = domain%blocks_info(k)%k_nym

            sync_dir(5) = _NXP_NYP_; k_dir(5) = domain%blocks_info(k)%k_nxp_nyp
            sync_dir(6) = _NXP_NYM_; k_dir(6) = domain%blocks_info(k)%k_nxp_nym
            sync_dir(7) = _NXM_NYP_; k_dir(7) = domain%blocks_info(k)%k_nxm_nyp
            sync_dir(8) = _NXM_NYM_; k_dir(8) = domain%blocks_info(k)%k_nxm_nym
            
            do kk = 1, 8
                call get_boundary_points_of_block(domain, k_dir(kk), get_inverse_dir(sync_dir(kk)), nxs, nxe, nys, nye)
                call get_halo_points_of_block(domain, k, sync_dir(kk), hnxs, hnxe, hnys, hnye)
                data2d%block(k)%field(hnxs : hnxe, hnys : hnye) = data2d%block(k_dir(kk))%field(nxs : nxe, nys : nye)
            enddo
        enddo

#undef _MPI_TYPE_
            end subroutine

end module mpp_sync_module