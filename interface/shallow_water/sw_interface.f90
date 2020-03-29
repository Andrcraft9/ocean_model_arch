module shallow_water_interface_module
#include "core/kernel_macros.fi"

    use kernel_interface_module, only: set_kernel_interface
    use mpp_module, only: mpp_rank, mpp_count, mpp_cart_comm, mpp_size, mpp_coord, mpp_period
    use decomposition_module, only: domain_type
    use data_types_module, only: data2D_real8_type, data3D_real8_type
    use ocean_module, only: ocean_type
    use grid_module, only: grid_type
    use mpp_sync_module, only: sync 
    use depth_module, only: hh_init_kernel, hh_update_kernel, hh_shift_kernel
    use velssh_sw_module, only: uv_trans_vort_kernel, uv_trans_kernel, uv_diff2_kernel


    implicit none
    save
    private

    public :: envoke_hh_init_kernel, envoke_hh_update_kernel, envoke_hh_shift_kernel
    public :: envoke_uv_trans_vort_kernel, envoke_uv_trans_kernel, envoke_uv_diff2_kernel

contains

    subroutine envoke_hh_init_kernel(domain, grid_data, ssh, sshp)

        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(data2D_real8_type), intent(in) :: ssh, sshp

        integer :: k

        do k = 1, domain%bcount
            
            call set_kernel_interface(domain, k)

            call hh_init_kernel(grid_data   %lu        %block(k)%field,  &
                                grid_data   %llu       %block(k)%field,  &
                                grid_data   %llv       %block(k)%field,  & 
                                grid_data   %luh       %block(k)%field,  &
                                grid_data   %dx        %block(k)%field,  &
                                grid_data   %dy        %block(k)%field,  &
                                grid_data   %dxt       %block(k)%field,  &
                                grid_data   %dyt       %block(k)%field,  &
                                grid_data   %dxh       %block(k)%field,  &
                                grid_data   %dyh       %block(k)%field,  &
                                grid_data   %dxb       %block(k)%field,  &
                                grid_data   %dyb       %block(k)%field,  &
                                grid_data   %hhq       %block(k)%field,  &
                                grid_data   %hhq_p     %block(k)%field,  & 
                                grid_data   %hhq_n     %block(k)%field,  &
                                grid_data   %hhu       %block(k)%field,  & 
                                grid_data   %hhu_p     %block(k)%field,  & 
                                grid_data   %hhu_n     %block(k)%field,  & 
                                grid_data   %hhv       %block(k)%field,  & 
                                grid_data   %hhv_p     %block(k)%field,  & 
                                grid_data   %hhv_n     %block(k)%field,  & 
                                grid_data   %hhh       %block(k)%field,  & 
                                grid_data   %hhh_p     %block(k)%field,  & 
                                grid_data   %hhh_n     %block(k)%field,  & 
                                             ssh       %block(k)%field,  &
                                             sshp      %block(k)%field,  &
                                grid_data   %hhq_rest  %block(k)%field)

        enddo

        call sync(domain, grid_data%hhu)
        call sync(domain, grid_data%hhu_p)
        call sync(domain, grid_data%hhu_n)
        call sync(domain, grid_data%hhv)
        call sync(domain, grid_data%hhv_p)
        call sync(domain, grid_data%hhv_n)
        call sync(domain, grid_data%hhh)
        call sync(domain, grid_data%hhh_p)
        call sync(domain, grid_data%hhh_n)
        
    end subroutine

    subroutine envoke_hh_update_kernel(domain, grid_data, ssh)

        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(data2D_real8_type), intent(in) :: ssh

        integer :: k

        do k = 1, domain%bcount

            call set_kernel_interface(domain, k)

            call hh_update_kernel(grid_data   %lu        %block(k)%field,  & 
                                  grid_data   %llu       %block(k)%field,  & 
                                  grid_data   %llv       %block(k)%field,  & 
                                  grid_data   %luh       %block(k)%field,  & 
                                  grid_data   %dx        %block(k)%field,  & 
                                  grid_data   %dy        %block(k)%field,  & 
                                  grid_data   %dxt       %block(k)%field,  & 
                                  grid_data   %dyt       %block(k)%field,  & 
                                  grid_data   %dxh       %block(k)%field,  & 
                                  grid_data   %dyh       %block(k)%field,  & 
                                  grid_data   %dxb       %block(k)%field,  & 
                                  grid_data   %dyb       %block(k)%field,  &
                                  grid_data   %hhq_n     %block(k)%field,  &
                                  grid_data   %hhu_n     %block(k)%field,  & 
                                  grid_data   %hhv_n     %block(k)%field,  & 
                                  grid_data   %hhh_n     %block(k)%field,  & 
                                               ssh       %block(k)%field,  & 
                                  grid_data   %hhq_rest  %block(k)%field)

        enddo

        call sync(domain, grid_data%hhu_n)
        call sync(domain, grid_data%hhv_n)
        call sync(domain, grid_data%hhh_n)
    
    end subroutine

    subroutine envoke_hh_shift_kernel(domain, grid_data)

        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data

        integer :: k

        do k = 1, domain%bcount

            call set_kernel_interface(domain, k)

            call hh_shift_kernel(grid_data  %lu     %block(k)%field,  &
                                 grid_data  %llu    %block(k)%field,  &
                                 grid_data  %llv    %block(k)%field,  &
                                 grid_data  %luh    %block(k)%field,  &
                                 grid_data  %hhq    %block(k)%field,  &
                                 grid_data  %hhq_p  %block(k)%field,  &
                                 grid_data  %hhq_n  %block(k)%field,  &
                                 grid_data  %hhu    %block(k)%field,  &
                                 grid_data  %hhu_p  %block(k)%field,  &
                                 grid_data  %hhu_n  %block(k)%field,  &
                                 grid_data  %hhv    %block(k)%field,  &
                                 grid_data  %hhv_p  %block(k)%field,  &
                                 grid_data  %hhv_n  %block(k)%field,  &
                                 grid_data  %hhh    %block(k)%field,  &
                                 grid_data  %hhh_p  %block(k)%field,  &
                                 grid_data  %hhh_n  %block(k)%field)

        enddo

    end subroutine

    subroutine envoke_uv_trans_vort_kernel(domain, grid_data, u, v, vort, nlev)

        type(domain_type), intent(in) :: domain
        type(grid_type), intent(in) :: grid_data
        type(data2D_real8_type), intent(in) :: u, v
        type(data3D_real8_type), intent(inout) :: vort
        integer, intent(in) :: nlev

        integer :: k

        do k = 1, domain%bcount
            
            call set_kernel_interface(domain, k)

            call uv_trans_vort_kernel(grid_data  %luu   %block(k)%field,  &
                                      grid_data  %dxt   %block(k)%field,  &
                                      grid_data  %dyt   %block(k)%field,  &
                                      grid_data  %dxb   %block(k)%field,  &
                                      grid_data  %dyb   %block(k)%field,  &
                                                  u     %block(k)%field,  &
                                                  v     %block(k)%field,  &
                                                  vort  %block(k)%field,  &
                                                  nlev)

        enddo

        call sync(domain, vort, nlev)

    end subroutine

    subroutine envoke_uv_trans_kernel(domain, grid_data, u, v, vort, RHSx, RHSy, nlev)

        type(domain_type), intent(in) :: domain
        type(grid_type), intent(in) :: grid_data
        type(data2D_real8_type), intent(in) :: u, v
        type(data3D_real8_type), intent(in) :: vort
        type(data2D_real8_type), intent(inout) :: RHSx, RHSy
        integer, intent(in) :: nlev

        integer :: k

        do k = 1, domain%bcount
            
            call set_kernel_interface(domain, k)

            call uv_trans_kernel(grid_data  %lcu   %block(k)%field,  &
                                 grid_data  %lcv   %block(k)%field,  &
                                 grid_data  %luu   %block(k)%field,  &
                                 grid_data  %dxh   %block(k)%field,  &
                                 grid_data  %dyh   %block(k)%field,  &
                                             u     %block(k)%field,  &
                                             v     %block(k)%field,  &
                                             vort  %block(k)%field,  &
                                 grid_data  %hhq   %block(k)%field,  &
                                 grid_data  %hhu   %block(k)%field,  &
                                 grid_data  %hhv   %block(k)%field,  &
                                 grid_data  %hhh   %block(k)%field,  &
                                             RHSx  %block(k)%field,  &
                                             RHSy  %block(k)%field,  &
                                 nlev)
        enddo
    
    end subroutine

    subroutine envoke_uv_diff2_kernel(domain, grid_data, mu, str_t, str_s, RHSx, RHSy, nlev)

        type(domain_type), intent(in) :: domain
        type(grid_type), intent(in) :: grid_data
        type(data2D_real8_type), intent(in) :: mu, str_t, str_s
        type(data2D_real8_type), intent(inout) :: RHSx, RHSy
        integer, intent(in) :: nlev

        integer :: k

        do k = 1, domain%bcount
            
            call set_kernel_interface(domain, k)
            
            call uv_diff2_kernel(grid_data  %lcu    %block(k)%field,  &
                                 grid_data  %lcv    %block(k)%field,  &
                                 grid_data  %dx     %block(k)%field,  &
                                 grid_data  %dy     %block(k)%field,  &
                                 grid_data  %dxt    %block(k)%field,  &
                                 grid_data  %dyt    %block(k)%field,  &
                                 grid_data  %dxh    %block(k)%field,  &
                                 grid_data  %dyh    %block(k)%field,  &
                                 grid_data  %dxb    %block(k)%field,  &
                                 grid_data  %dyb    %block(k)%field,  &
                                             mu     %block(k)%field,  &
                                             str_t  %block(k)%field,  &
                                             str_s  %block(k)%field,  &
                                 grid_data  %hhq    %block(k)%field,  &
                                 grid_data  %hhu    %block(k)%field,  &
                                 grid_data  %hhv    %block(k)%field,  &
                                 grid_data  %hhh    %block(k)%field,  &
                                             RHSx   %block(k)%field,  &
                                             RHSy   %block(k)%field,  &
                                 nlev)

        enddo
    
    end subroutine

end module shallow_water_interface_module