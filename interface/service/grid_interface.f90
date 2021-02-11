module grid_interface_module

    use kernel_interface_module
    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use mpp_module
    use decomposition_module, only: domain_type, domain => domain_data
    use data_types_module, only: data2D_real8_type, data2D_real4_type
    use grid_module, only: grid_type, grid_data, grid_global_type, grid_global_data
    use mpp_sync_module, only: hybrid_sync, sync, sync_parameters_type
    use grid_kernels_module, only: lu_init_kernel, lu_lv_init_kernel, grid_base_init_kernel, grid_geo_init_kernel
    use errors_module, only: abort_model, check_error

#include "macros/kernel_macros.fi"
#include "macros/mpp_macros.fi"

    implicit none
    save
    private

    public :: envoke_lu_init_kernel
    public :: envoke_lu_lv_init_kernel
    public :: envoke_grid_base_init_kernel
    public :: envoke_grid_geo_init_kernel

contains

!-----------------------------------------------------------------------------!
!------------------------------- Kernels -------- ----------------------------!
!-----------------------------------------------------------------------------!

subroutine envoke_lu_init_kernel()
    integer :: k

    _OMP_BLOCKS_PARALLEL_BEGIN_
    do k = 1, domain%bcount

        call lu_init_kernel(domain%bbnd_x1(k),   domain%bbnd_x2(k), domain%bbnd_y1(k),   domain%bbnd_y2(k),  &
                            grid_global_data %mask,                  &
                            grid_data        %lu   %block(k)%field,  &
                            grid_data        %lu1  %block(k)%field)

    enddo
    _OMP_BLOCKS_PARALLEL_END_

    call sync(domain, grid_data%lu)
        
end subroutine

subroutine envoke_lu_lv_init_kernel()
    integer :: k

    _OMP_BLOCKS_PARALLEL_BEGIN_
    do k = 1, domain%bcount
        call lu_lv_init_kernel(domain%bbnd_x1(k),   domain%bbnd_x2(k), domain%bbnd_y1(k),   domain%bbnd_y2(k),  &
                               grid_data  %lu   %block(k)%field,  &
                               grid_data  %luh  %block(k)%field,  &
                               grid_data  %luu  %block(k)%field,  &
                               grid_data  %llu  %block(k)%field,  &
                               grid_data  %llv  %block(k)%field,  &
                               grid_data  %lcu  %block(k)%field,  &
                               grid_data  %lcv  %block(k)%field)
    enddo
    _OMP_BLOCKS_PARALLEL_END_

    call sync(domain, grid_data%luh)
    call sync(domain, grid_data%luu)
    call sync(domain, grid_data%lcu)
    call sync(domain, grid_data%llu)
    call sync(domain, grid_data%lcv)
    call sync(domain, grid_data%llv)

end subroutine

subroutine envoke_grid_base_init_kernel()
    integer :: k

    _OMP_BLOCKS_PARALLEL_BEGIN_
    do k = 1, domain%bcount
        call grid_base_init_kernel(domain%bnx_start(k), domain%bnx_end(k), domain%bny_start(k), domain%bny_end(k),  &
                                   domain%bbnd_x1(k),   domain%bbnd_x2(k), domain%bbnd_y1(k),   domain%bbnd_y2(k),  &
                                   grid_data  %xt     %block(k)%field,  &
                                   grid_data  %yt     %block(k)%field,  &
                                   grid_data  %xu     %block(k)%field,  &
                                   grid_data  %yv     %block(k)%field,  &
                                   grid_data  %dx     %block(k)%field,  &
                                   grid_data  %dy     %block(k)%field,  &
                                   grid_data  %dxt    %block(k)%field,  &
                                   grid_data  %dyt    %block(k)%field,  &
                                   grid_data  %dxh    %block(k)%field,  &
                                   grid_data  %dyh    %block(k)%field,  &
                                   grid_data  %dxb    %block(k)%field,  &
                                   grid_data  %dyb    %block(k)%field,  &
                                   grid_data  %rlh_s  %block(k)%field,  &
                                   grid_data  %rlh_c  %block(k)%field)
    enddo
    _OMP_BLOCKS_PARALLEL_END_

    call sync(domain, grid_data%dxt)
    call sync(domain, grid_data%dxb)
    call sync(domain, grid_data%dx)
    call sync(domain, grid_data%dxh)
    call sync(domain, grid_data%dyt)
    call sync(domain, grid_data%dyb)
    call sync(domain, grid_data%dy)
    call sync(domain, grid_data%dyh)

end subroutine

subroutine envoke_grid_geo_init_kernel()
    integer :: k

    _OMP_BLOCKS_PARALLEL_BEGIN_
    do k = 1, domain%bcount
        call grid_geo_init_kernel(domain%bnx_start(k), domain%bnx_end(k), domain%bny_start(k), domain%bny_end(k),  &
                                  domain%bbnd_x1(k),   domain%bbnd_x2(k), domain%bbnd_y1(k),   domain%bbnd_y2(k),  &
                                  grid_data  %xt            %block(k)%field,  &
                                  grid_data  %yt            %block(k)%field,  &
                                  grid_data  %xu            %block(k)%field,  &
                                  grid_data  %yv            %block(k)%field,  &
                                  grid_data  %dxt           %block(k)%field,  &
                                  grid_data  %dxb           %block(k)%field,  &
                                  grid_data  %dx            %block(k)%field,  &
                                  grid_data  %dxh           %block(k)%field,  &
                                  grid_data  %dyt           %block(k)%field,  &
                                  grid_data  %dyb           %block(k)%field,  &
                                  grid_data  %dy            %block(k)%field,  &
                                  grid_data  %dyh           %block(k)%field,  &
                                  grid_data  %geo_lon_t     %block(k)%field,  &
                                  grid_data  %geo_lat_t     %block(k)%field,  &
                                  grid_data  %geo_lon_u     %block(k)%field,  &
                                  grid_data  %geo_lat_u     %block(k)%field,  &
                                  grid_data  %geo_lon_v     %block(k)%field,  &
                                  grid_data  %geo_lat_v     %block(k)%field,  &
                                  grid_data  %geo_lon_h     %block(k)%field,  &
                                  grid_data  %geo_lat_h     %block(k)%field,  &
                                  grid_data  %rlh_s         %block(k)%field,  &
                                  grid_data  %rlh_c         %block(k)%field,  &
                                  grid_data  %rotvec_coeff  %block(k)%field,  &
                                  grid_data  %sqt           %block(k)%field,  &
                                  grid_data  %squ           %block(k)%field,  &
                                  grid_data  %sqv           %block(k)%field,  &
                                  grid_data  %sqh           %block(k)%field,  &
                                  grid_data  %rlh_sqh       %block(k)%field)
     enddo
    _OMP_BLOCKS_PARALLEL_END_

end subroutine

end module grid_interface_module