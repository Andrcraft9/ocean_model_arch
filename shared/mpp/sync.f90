#include "macros/mpp_macros.fi"
#include "macros/kernel_macros.fi"

module mpp_sync_module
    ! Sync data of different types, depends on domain

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    !use mpi
    use debug_module
    use mpp_module
    use errors_module, only: abort_model, check_error
    use data_types_module, only: block1D_real4_type, block1D_real8_type, block2D_real4_type, block2D_real8_type,  &
                                 data1D_real4_type, data1D_real8_type, data2D_real4_type, data2D_real8_type,      &
                                 data3D_real4_type, data3D_real8_type
    use decomposition_module, only: domain_type, get_boundary_points_of_block, get_halo_points_of_block,             &
                                    get_halo_and_boundary_points_of_block, get_expanded_halo_points_of_block,        &
                                    get_inverse_dir, get_rank_of_block, get_loc_num_of_block, fill_direction_array,  &
                                    is_block_on_current_proc, is_inner_block
#ifdef _GPU_MODE_
    use cudafor
#endif

    implicit none
    save
    private

    type, public :: sync_parameters_type
        integer :: sync_mode ! 0 - inner; 1 - boundary; 2 - intermediate; 3 - all, 5 - sync between device and host
        integer :: sync_device_host ! Check it only sync_mode is 5.
                                    ! 0 - dev to host (sync), 1 - host to dev (sync), 2 - dev to host (async), 3 - host to dev (async)
        integer :: data_id ! id of field to sync (temporary)
    end type sync_parameters_type

    interface sync_test
        module procedure syncborder_test_data2D_real8
        module procedure syncborder_test_data2D_real4
    end interface

    interface sync
        module procedure syncborder_data2D_real8
        module procedure syncborder_data2D_real4
    end interface

    interface hybrid_sync
        module procedure syncborder_data2D_hybrid_real8
        module procedure syncborder_data2D_hybrid_real4
    end interface

    interface sync_inner
        module procedure syncborder_data2D_inner_real8
        module procedure syncborder_data2D_inner_real4
    end interface

    interface sync_boundary
        module procedure syncborder_data2D_boundary_real8
        module procedure syncborder_data2D_boundary_real4
    end interface

    interface sync_intermediate
        module procedure syncborder_data2D_intermediate_real8
        module procedure syncborder_data2D_intermediate_real4
    end interface

    public :: mpp_sync_init
    public :: mpp_sync_finalize
    public :: sync, hybrid_sync

    public :: sync_test

    ! MPI info
    integer :: sync_count_send_recv
    integer, allocatable :: sync_requests(:), sync_statuses(:, :)
    ! Buffers for each near rank
    real(wp8), allocatable :: sync_send_buf_r8(:, :, :)
    real(wp4), allocatable :: sync_send_buf_r4(:, :, :)
    real(wp8), allocatable :: sync_recv_buf_r8(:, :, :)
    real(wp4), allocatable :: sync_recv_buf_r4(:, :, :)
    integer, allocatable :: sync_buf_size(:)      ! Size of buffer for each near rank
    integer, allocatable :: sync_buf_pos(:, :)    ! Position of buffer for each near rank and boundary block
    integer, allocatable :: sync_recv_buf_pos(:, :)


    integer, allocatable :: sync_map_rank(:)  ! Map: Gloabl rank to local rank number. 
                                              ! Global rank numbers - MPI ranks, from 0 to mpp_count - 1
                                              ! Local rank numbers  - indecies for all sync buffers, from 1 to domain%amount_of_ranks_near
    integer, allocatable :: sync_map_bb(:)    ! Map: Block number to boundary blocks number
                                              ! Block number: from 1 to bcount
                                              ! Boundary block number: from 1 to bcount_boundary

#ifdef _GPU_MODE_
    integer(kind=cuda_stream_kind), allocatable, public :: mpp_sync_cuda_streams(:)
#endif

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

        integer :: k, kk, rank, ierr, istat
        integer :: sync_dir(8)
        integer :: k_dir(8)
        integer :: rank_dir(8)
        integer :: nxs, nxe, nys, nye         ! Boundary points
        integer :: rk                         ! Local index of buffers for rank, see sync_map_rank
        integer :: bbk                        ! Local index of boundary block
        integer :: max_buf_size = 1, min_buf_size = huge(0)
        integer :: max_val, min_val
        integer :: leastpriority = 0, greatestpriority = 0, priority

        ! MPI info
        allocate(sync_requests(_MPP_MAX_SIMUL_SYNCS_ * 2 * domain%amount_of_ranks_near),  &
                 sync_statuses(MPI_STATUS_SIZE, _MPP_MAX_SIMUL_SYNCS_ * 2 * domain%amount_of_ranks_near))

        ! Create map: global rank numeration to local rank numeration used as indicies in MPI buffers
        allocate(sync_map_rank(0 : mpp_count - 1))
        sync_map_rank = -1
        do k = 1, domain%amount_of_ranks_near
            sync_map_rank(domain%ranks_near(k)) = k
        enddo

        ! Create map:
        allocate(sync_map_bb(domain%bcount))
        sync_map_bb = -1
        do k = 1, domain%bcount_boundary
            sync_map_bb(domain%boundary_blocks(k)) = k
        enddo

        allocate(sync_buf_size(domain%amount_of_ranks_near))
        allocate(sync_buf_pos(8, domain%bcount_boundary))
        allocate(sync_recv_buf_pos(8, domain%bcount_boundary))
        sync_buf_size = 0
        sync_buf_pos = 0
        sync_recv_buf_pos = 0

        do k = 1, domain%bcount
            if (.not. domain%blocks_info(k)%is_inner) then
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
                            rk = sync_map_rank(rank_dir(kk))
                            bbk = sync_map_bb(k)
                            call get_boundary_points_of_block(domain, k, sync_dir(kk), nxs, nxe, nys, nye)
                            sync_buf_pos(kk, bbk) = sync_buf_pos(kk, bbk) + sync_buf_size(rk) + 1
                            sync_buf_size(rk) = sync_buf_size(rk) + 1 + (nxe - nxs + 1)*(nye - nys + 1)
                            !sync_buf_size(rk) = sync_buf_size(rk) + 2 + (nxe - nxs + 1)*(nye - nys + 1)
                        endif
                    endif
                enddo
            endif
        enddo

        do rk = 1, domain%amount_of_ranks_near
            if (max_buf_size < sync_buf_size(rk)) max_buf_size = sync_buf_size(rk)
            if (min_buf_size > sync_buf_size(rk)) min_buf_size = sync_buf_size(rk)
        enddo

        if (mpp_count > 1) then
            call mpi_allreduce(max_buf_size, max_val, 1, mpi_integer, mpi_max, mpp_cart_comm, ierr)
            call mpi_allreduce(min_buf_size, min_val, 1, mpi_integer, mpi_min, mpp_cart_comm, ierr)
        else
            max_val = 0; max_buf_size = 0
            min_val = 0; min_buf_size = 0
        endif

        if (mpp_is_master()) then
            print *, "SYNC INFO: Allreduced Min/Max buffer size = ", min_val, max_val
        endif

        if (debug_level >= 5) then
            do rank = 0, mpp_count - 1
                if (mpp_rank == rank) then
                    print *, "SYNC INFO: rank, Min/Max buffer size = ", mpp_rank, min_buf_size, max_buf_size
                endif
                call mpi_barrier(mpp_cart_comm, ierr)
            enddo
        endif

        allocate(sync_send_buf_r4(max_buf_size, domain%amount_of_ranks_near, _MPP_MAX_SIMUL_SYNCS_))
        allocate(sync_send_buf_r8(max_buf_size, domain%amount_of_ranks_near, _MPP_MAX_SIMUL_SYNCS_))
        allocate(sync_recv_buf_r4(max_buf_size, domain%amount_of_ranks_near, _MPP_MAX_SIMUL_SYNCS_))
        allocate(sync_recv_buf_r8(max_buf_size, domain%amount_of_ranks_near, _MPP_MAX_SIMUL_SYNCS_))
        sync_send_buf_r4 = 0; sync_send_buf_r8 = 0
        sync_recv_buf_r4 = 0; sync_recv_buf_r8 = 0

        sync_count_send_recv = 0

#ifdef _GPU_MODE_
        allocate(mpp_sync_cuda_streams(domain%bcount))
        ! Set default CUDA stream
        mpp_sync_cuda_streams = 0
#endif

#ifdef _GPU_ASYNC_
        istat = cudaDeviceGetStreamPriorityRange(leastpriority, greatestpriority)
        priority = greatestpriority
        do k = 1, domain%bcount
            ! dont use priority streams
            !istat = cudaStreamCreateWithPriority(mpp_sync_cuda_streams(k), cudaStreamDefault, priority)
            !priority = priority + 1
            !if (priority > leastpriority) then
            !    priority = leastpriority
            !endif
            istat = cudaStreamCreate(mpp_sync_cuda_streams(k))
        enddo
#endif

        call mpp_sync_output()
    end subroutine

    subroutine deallocate_mpp_sync_buffers(domain)
        type(domain_type), intent(in) :: domain
        
        integer :: k, istat

        deallocate(sync_requests, sync_statuses)

        deallocate(sync_send_buf_r4)
        deallocate(sync_send_buf_r8)
        deallocate(sync_recv_buf_r4)
        deallocate(sync_recv_buf_r8)
        deallocate(sync_buf_size)
        deallocate(sync_buf_pos)
        deallocate(sync_recv_buf_pos)
        
        deallocate(sync_map_rank)
        deallocate(sync_map_bb)

#ifdef _GPU_ASYNC_
        do k = 1, domain%bcount
            istat = cudaStreamDestroy(mpp_sync_cuda_streams(k))
        enddo
#endif

#ifdef _GPU_MODE_
        deallocate(mpp_sync_cuda_streams)
#endif
    end subroutine

!------------------------------------------------------------------------------
! Test

    subroutine syncborder_test_data2D_real8(domain, data2d)
#define _MPI_TYPE_ mpi_real8
        
        implicit none
        type(data2D_real8_type), intent(inout) :: data2d
        type(domain_type), intent(in) :: domain

#include    "syncborder_block2D_gen_test.fi"

#undef _MPI_TYPE_

    end subroutine

    subroutine syncborder_test_data2D_real4(domain, data2d)
#define _MPI_TYPE_ mpi_real4
        
        implicit none
        type(data2D_real4_type), intent(inout) :: data2d
        type(domain_type), intent(in) :: domain

#include    "syncborder_block2D_gen_test.fi"

#undef _MPI_TYPE_

    end subroutine

!------------------------------------------------------------------------------
! Generic subroutines, hybrid:
        subroutine syncborder_data2D_hybrid_real8(k, sync_parameters, sync_tag, domain, data2d)
            implicit none
            integer, intent(in) :: k
            type(sync_parameters_type), intent(in) :: sync_parameters
            integer, intent(in) :: sync_tag
            type(data2D_real8_type), intent(inout) :: data2d
            type(domain_type), intent(in) :: domain

            if (sync_tag < 1) call abort_model("Sync Error: Invalid sync tag")

            if (sync_parameters%sync_mode == 5) then
                if (k < 1) call abort_model("Sync Error: Invalid block number for sync betweeb host and device")
                select case(sync_parameters%sync_device_host)
                    case(0)
                        call sync_device_to_host_data2D_real8(k, domain, data2d, .false.)
                    case(1)
                        call sync_host_to_device_data2D_real8(k, domain, data2d, .false.)
                    case(2)
                        call sync_device_to_host_data2D_real8(k, domain, data2d, .true.)
                    case(3)
                        call sync_host_to_device_data2D_real8(k, domain, data2d, .true.)
                    case default
                        call abort_model("Sync Error: Unknown sync between device and host")
                end select
                return
            endif

            select case(sync_parameters%sync_mode)
                case(0)
                    call syncborder_data2D_inner_real8(sync_tag, domain, data2d)
                case(1)
                    call syncborder_data2D_boundary_real8(sync_tag, domain, data2d)
                case(2)
                    call syncborder_data2D_intermediate_real8(sync_tag, domain, data2d)
                case(3)
                    call syncborder_data2D_real8(domain, data2d)
                case default
                    call abort_model("Sync Error: Unknown sync mode")
            end select
        end subroutine

        subroutine syncborder_data2D_hybrid_real4(k, sync_parameters, sync_tag, domain, data2d)
            implicit none
            integer, intent(in) :: k
            type(sync_parameters_type), intent(in) :: sync_parameters
            integer, intent(in) :: sync_tag
            type(data2D_real4_type), intent(inout) :: data2d
            type(domain_type), intent(in) :: domain

            if (sync_tag < 1) call abort_model("Sync Error: Invalid sync tag")

            if (sync_parameters%sync_mode == 5) then
                if (k < 1) call abort_model("Sync Error: Invalid block number for sync betweeb host and device")
                select case(sync_parameters%sync_device_host)
                    case(0)
                        call sync_device_to_host_data2D_real4(k, domain, data2d, .false.)
                    case(1)
                        call sync_host_to_device_data2D_real4(k, domain, data2d, .false.)
                    case(2)
                        call sync_device_to_host_data2D_real4(k, domain, data2d, .true.)
                    case(3)
                        call sync_host_to_device_data2D_real4(k, domain, data2d, .true.)
                    case default
                        call abort_model("Sync Error: Unknown sync between device and host")
                end select
                return
            endif

            select case(sync_parameters%sync_mode)
                case(0)
                    call syncborder_data2D_inner_real4(sync_tag, domain, data2d)
                case(1)
                    call syncborder_data2D_boundary_real4(sync_tag, domain, data2d)
                case(2)
                    call syncborder_data2D_intermediate_real4(sync_tag, domain, data2d)
                case(3)
                    call syncborder_data2D_real4(domain, data2d)
                case default
                    call abort_model("Sync Error: Unknown sync mode")
            end select
        end subroutine

!------------------------------------------------------------------------------
! Generic subroutines, sync between host and device:
        subroutine sync_device_to_host_data2D_real8(k, domain, data2d, is_async)
            implicit none
            integer, intent(in) :: k
            type(domain_type), intent(in) :: domain
            type(data2D_real8_type), intent(inout) :: data2d
            logical, intent(in) :: is_async

#ifdef _GPU_MODE_
            integer :: kk, istat
            integer :: nxs, nxe, nys, nye     ! Boundary points
            integer :: sync_dir(4)

            sync_dir(1) = _NXP_
            sync_dir(2) = _NXM_
            sync_dir(3) = _NYP_
            sync_dir(4) = _NYM_

            do kk = 1, 4
                call get_halo_and_boundary_points_of_block(domain, k, sync_dir(kk), nxs,  nxe,  nys,  nye)
                if (is_async) then
                    istat = cudaMemcpy2DAsync(data2d%block(k)%field(nxs, nys),            &
                                              domain%bbnd_x2(k) - domain%bbnd_x1(k) + 1,  &
                                              data2d%block(k)%field_device(nxs, nys),     &
                                              domain%bbnd_x2(k) - domain%bbnd_x1(k) + 1,  &
                                              nxe - nxs + 1,  &
                                              nye - nys + 1,  &
                                              stream=mpp_sync_cuda_streams(k))
                else
                    !data2d%block(k)%field(nxs : nxe, nys : nye) = data2d%block(k)%field_device(nxs : nxe, nys : nye)
                    istat = cudaMemcpy2D(data2d%block(k)%field(nxs, nys),            &
                                         domain%bbnd_x2(k) - domain%bbnd_x1(k) + 1,  &
                                         data2d%block(k)%field_device(nxs, nys),     &
                                         domain%bbnd_x2(k) - domain%bbnd_x1(k) + 1,  &
                                         nxe - nxs + 1,  &
                                         nye - nys + 1)
                endif
            enddo
#endif
        end subroutine

        subroutine sync_host_to_device_data2D_real8(k, domain, data2d, is_async)
            implicit none
            integer, intent(in) :: k
            type(domain_type), intent(in) :: domain
            type(data2D_real8_type), intent(inout) :: data2d
            logical, intent(in) :: is_async

#ifdef _GPU_MODE_
            integer :: kk, istat
            integer :: nxs, nxe, nys, nye     ! Boundary points
            integer :: sync_dir(4)

            sync_dir(1) = _NXP_
            sync_dir(2) = _NXM_
            sync_dir(3) = _NYP_
            sync_dir(4) = _NYM_

            do kk = 1, 4
                call get_expanded_halo_points_of_block(domain, k, sync_dir(kk), nxs,  nxe,  nys,  nye)
                if (is_async) then
                    istat = cudaMemcpy2DAsync(data2d%block(k)%field_device(nxs, nys),     &
                                              domain%bbnd_x2(k) - domain%bbnd_x1(k) + 1,  &
                                              data2d%block(k)%field(nxs, nys),            &
                                              domain%bbnd_x2(k) - domain%bbnd_x1(k) + 1,  &
                                              nxe - nxs + 1,  &
                                              nye - nys + 1,  &
                                              stream=mpp_sync_cuda_streams(k))
                else
                    !data2d%block(k)%field_device(nxs : nxe, nys : nye) = data2d%block(k)%field(nxs : nxe, nys : nye)
                    istat = cudaMemcpy2D(data2d%block(k)%field_device(nxs, nys),            &
                                                domain%bbnd_x2(k) - domain%bbnd_x1(k) + 1,  &
                                                data2d%block(k)%field(nxs, nys),            &
                                                domain%bbnd_x2(k) - domain%bbnd_x1(k) + 1,  &
                                                nxe - nxs + 1,  &
                                                nye - nys + 1)
                endif
            enddo
#endif
        end subroutine

        subroutine sync_device_to_host_data2D_real4(k, domain, data2d, is_async)
            implicit none
            integer, intent(in) :: k
            type(domain_type), intent(in) :: domain
            type(data2D_real4_type), intent(inout) :: data2d
            logical, intent(in) :: is_async

#ifdef _GPU_MODE_
            integer :: kk, istat
            integer :: nxs, nxe, nys, nye     ! Boundary points
            integer :: sync_dir(4)

            sync_dir(1) = _NXP_
            sync_dir(2) = _NXM_
            sync_dir(3) = _NYP_
            sync_dir(4) = _NYM_

            do kk = 1, 4
                call get_halo_and_boundary_points_of_block(domain, k, sync_dir(kk), nxs,  nxe,  nys,  nye)
                if (is_async) then
                    istat = cudaMemcpy2DAsync(data2d%block(k)%field(nxs, nys),            &
                                              domain%bbnd_x2(k) - domain%bbnd_x1(k) + 1,  &
                                              data2d%block(k)%field_device(nxs, nys),     &
                                              domain%bbnd_x2(k) - domain%bbnd_x1(k) + 1,  &
                                              nxe - nxs + 1,  &
                                              nye - nys + 1,  &
                                              stream=mpp_sync_cuda_streams(k))
                else
                    !data2d%block(k)%field(nxs : nxe, nys : nye) = data2d%block(k)%field_device(nxs : nxe, nys : nye)
                    istat = cudaMemcpy2D(data2d%block(k)%field(nxs, nys),            &
                                         domain%bbnd_x2(k) - domain%bbnd_x1(k) + 1,  &
                                         data2d%block(k)%field_device(nxs, nys),     &
                                         domain%bbnd_x2(k) - domain%bbnd_x1(k) + 1,  &
                                         nxe - nxs + 1,  &
                                         nye - nys + 1)
                endif
            enddo
#endif
        end subroutine

        subroutine sync_host_to_device_data2D_real4(k, domain, data2d, is_async)
            implicit none
            integer, intent(in) :: k
            type(domain_type), intent(in) :: domain
            type(data2D_real4_type), intent(inout) :: data2d
            logical, intent(in) :: is_async

#ifdef _GPU_MODE_
            integer :: kk, istat
            integer :: nxs, nxe, nys, nye     ! Boundary points
            integer :: sync_dir(4)

            sync_dir(1) = _NXP_
            sync_dir(2) = _NXM_
            sync_dir(3) = _NYP_
            sync_dir(4) = _NYM_

            do kk = 1, 4
                call get_expanded_halo_points_of_block(domain, k, sync_dir(kk), nxs,  nxe,  nys,  nye)
                if (is_async) then
                    istat = cudaMemcpy2DAsync(data2d%block(k)%field_device(nxs, nys),     &
                                              domain%bbnd_x2(k) - domain%bbnd_x1(k) + 1,  &
                                              data2d%block(k)%field(nxs, nys),            &
                                              domain%bbnd_x2(k) - domain%bbnd_x1(k) + 1,  &
                                              nxe - nxs + 1,  &
                                              nye - nys + 1,  &
                                              stream=mpp_sync_cuda_streams(k))
                else
                    !data2d%block(k)%field_device(nxs : nxe, nys : nye) = data2d%block(k)%field(nxs : nxe, nys : nye)
                    istat = cudaMemcpy2D(data2d%block(k)%field_device(nxs, nys),            &
                                                domain%bbnd_x2(k) - domain%bbnd_x1(k) + 1,  &
                                                data2d%block(k)%field(nxs, nys),            &
                                                domain%bbnd_x2(k) - domain%bbnd_x1(k) + 1,  &
                                                nxe - nxs + 1,  &
                                                nye - nys + 1)
                endif
            enddo
#endif
        end subroutine

!------------------------------------------------------------------------------
! Generic subroutines, all sync:

        subroutine syncborder_data2D_real8(domain, data2d)
#define _MPI_TYPE_ mpi_real8
#define _SYNC_SEND_BUF_ sync_send_buf_r8
#define _SYNC_RECV_BUF_ sync_recv_buf_r8
            
            implicit none
            type(data2D_real8_type), intent(inout) :: data2d
            type(domain_type), intent(in) :: domain
            integer, parameter :: sync_tag = 1

#include    "syncborder_block2D_gen_all.fi"

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
            integer, parameter :: sync_tag = 1

#include    "syncborder_block2D_gen_all.fi"

#undef _MPI_TYPE_
#undef _SYNC_SEND_BUF_
#undef _SYNC_RECV_BUF_
        end subroutine

!------------------------------------------------------------------------------
! Generic subroutines, inner:

        subroutine syncborder_data2D_inner_real8(sync_tag, domain, data2d)
#define _MPI_TYPE_ mpi_real8
#define _SYNC_SEND_BUF_ sync_send_buf_r8
#define _SYNC_RECV_BUF_ sync_recv_buf_r8
        
        implicit none
        integer, intent(in) :: sync_tag
        type(data2D_real8_type), intent(inout) :: data2d
        type(domain_type), intent(in) :: domain

#include    "syncborder_block2D_gen_inner.fi"

#undef _MPI_TYPE_
#undef _SYNC_SEND_BUF_
#undef _SYNC_RECV_BUF_
        end subroutine

        subroutine syncborder_data2D_inner_real4(sync_tag, domain, data2d)
#define _MPI_TYPE_ mpi_real4
#define _SYNC_SEND_BUF_ sync_send_buf_r4
#define _SYNC_RECV_BUF_ sync_recv_buf_r4

        implicit none
        integer, intent(in) :: sync_tag
        type(data2D_real4_type), intent(inout) :: data2d
        type(domain_type), intent(in) :: domain

#include    "syncborder_block2D_gen_inner.fi"

#undef _MPI_TYPE_
#undef _SYNC_SEND_BUF_
#undef _SYNC_RECV_BUF_
        end subroutine

!------------------------------------------------------------------------------
! Generic subroutines, boundary:

        subroutine syncborder_data2D_boundary_real8(sync_tag, domain, data2d)
#define _MPI_TYPE_ mpi_real8
#define _SYNC_SEND_BUF_ sync_send_buf_r8
#define _SYNC_RECV_BUF_ sync_recv_buf_r8
        
        implicit none
        integer, intent(in) :: sync_tag
        type(data2D_real8_type), intent(inout) :: data2d
        type(domain_type), intent(in) :: domain

#include    "syncborder_block2D_gen_boundary.fi"

#undef _MPI_TYPE_
#undef _SYNC_SEND_BUF_
#undef _SYNC_RECV_BUF_
        end subroutine

        subroutine syncborder_data2D_boundary_real4(sync_tag, domain, data2d)
#define _MPI_TYPE_ mpi_real4
#define _SYNC_SEND_BUF_ sync_send_buf_r4
#define _SYNC_RECV_BUF_ sync_recv_buf_r4

        implicit none
        integer, intent(in) :: sync_tag
        type(data2D_real4_type), intent(inout) :: data2d
        type(domain_type), intent(in) :: domain

#include    "syncborder_block2D_gen_boundary.fi"

#undef _MPI_TYPE_
#undef _SYNC_SEND_BUF_
#undef _SYNC_RECV_BUF_
        end subroutine

!------------------------------------------------------------------------------
! Generic subroutines, intermediate:

        subroutine syncborder_data2D_intermediate_real8(sync_tag, domain, data2d)
#define _MPI_TYPE_ mpi_real8
#define _SYNC_SEND_BUF_ sync_send_buf_r8
#define _SYNC_RECV_BUF_ sync_recv_buf_r8
        
        implicit none
        integer, intent(in) :: sync_tag
        type(data2D_real8_type), intent(inout) :: data2d
        type(domain_type), intent(in) :: domain

#include    "syncborder_block2D_gen_intermediate.fi"

#undef _MPI_TYPE_
#undef _SYNC_SEND_BUF_
#undef _SYNC_RECV_BUF_
        end subroutine

        subroutine syncborder_data2D_intermediate_real4(sync_tag, domain, data2d)
#define _MPI_TYPE_ mpi_real4
#define _SYNC_SEND_BUF_ sync_send_buf_r4
#define _SYNC_RECV_BUF_ sync_recv_buf_r4

        implicit none
        integer, intent(in) :: sync_tag
        type(data2D_real4_type), intent(inout) :: data2d
        type(domain_type), intent(in) :: domain

#include    "syncborder_block2D_gen_intermediate.fi"

#undef _MPI_TYPE_
#undef _SYNC_SEND_BUF_
#undef _SYNC_RECV_BUF_
        end subroutine

end module mpp_sync_module