module shallow_water_interface_module
#include "macros/kernel_macros.fi"

    use kernel_interface_module
    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use mpp_module
    use decomposition_module, only: domain_type
    use data_types_module, only: data2D_real8_type, data2D_real4_type
    use ocean_module, only: ocean_type
    use grid_module, only: grid_type
    use mpp_sync_module, only: sync 
    use mixing_module, only: stress_components_kernel
    use depth_module, only: hh_init_kernel, hh_update_kernel, hh_shift_kernel
    use velssh_sw_module, only: check_ssh_err_kernel, sw_update_ssh_kernel, sw_update_uv, sw_next_step, uv_trans_vort_kernel, uv_trans_kernel, uv_diff2_kernel
    use errors_module, only: abort_model, check_error

    implicit none
    save
    private

    public :: envoke_check_ssh_err_kernel
    public :: envoke_stress_components_kernel
    public :: envoke_hh_init_kernel, envoke_hh_update_kernel, envoke_hh_shift_kernel
    public :: envoke_uv_trans_vort_kernel, envoke_uv_trans_kernel, envoke_uv_diff2_kernel
    public :: envoke_sw_update_ssh_kernel, envoke_sw_update_uv, envoke_sw_next_step

contains

    subroutine envoke_check_ssh_err_kernel(domain, grid_data, ssh, name)
    
    type(domain_type), intent(in) :: domain
    type(grid_type), intent(in) :: grid_data
    type(data2D_real8_type), intent(in) :: ssh
    character(*), intent(in) :: name

    integer :: k
    type(block_bounds_type) :: bb

    do k = 1, domain%bcount
            
        call bb%set(domain, k)

        call check_ssh_err_kernel(bb%nx_start, bb%nx_end, bb%ny_start, bb%ny_end,  &
                                  bb%bnd_x1,   bb%bnd_x2, bb%bnd_y1,   bb%bnd_y2,  &
                                  grid_data  %lu   %block(k)%field,  &
                                              ssh  %block(k)%field,  &
                                  name)

    enddo
    
    end subroutine

    subroutine envoke_stress_components_kernel(domain, grid_data, u, v, str_t, str_s)
    
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(data2D_real8_type), intent(in) :: u, v
        type(data2D_real8_type), intent(inout) :: str_t, str_s

        ! Interface only for 2D call
        integer, parameter :: nlev = 1

        integer :: k
        type(block_bounds_type) :: bb

        call start_kernel_timer('stress_components_kernel')
        do k = 1, domain%bcount
            
            call bb%set(domain, k)

            call stress_components_kernel(bb%nx_start, bb%nx_end, bb%ny_start, bb%ny_end,  &
                                          bb%bnd_x1,   bb%bnd_x2, bb%bnd_y1,   bb%bnd_y2,  &
                                          grid_data  %lu     %block(k)%field,  & 
                                          grid_data  %luu    %block(k)%field,  & 
                                          grid_data  %dx     %block(k)%field,  & 
                                          grid_data  %dy     %block(k)%field,  & 
                                          grid_data  %dxt    %block(k)%field,  & 
                                          grid_data  %dyt    %block(k)%field,  & 
                                          grid_data  %dxh    %block(k)%field,  & 
                                          grid_data  %dyh    %block(k)%field,  & 
                                          grid_data  %dxb    %block(k)%field,  & 
                                          grid_data  %dyb    %block(k)%field,  & 
                                                      u      %block(k)%field,  & 
                                                      v      %block(k)%field,  & 
                                                      str_t  %block(k)%field,  &
                                                      str_s  %block(k)%field,  &
                                          nlev)

        enddo
        call end_kernel_timer('stress_components_kernel')

        call sync(domain, str_t)
        call sync(domain, str_s)

    end subroutine

    subroutine envoke_hh_init_kernel(domain, grid_data, ssh, sshp)

        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(data2D_real8_type), intent(in) :: ssh, sshp

        integer :: k
        type(block_bounds_type) :: bb

        call start_kernel_timer('hh_init_kernel')
        do k = 1, domain%bcount
            
            call bb%set(domain, k)

            call hh_init_kernel(bb%nx_start, bb%nx_end, bb%ny_start, bb%ny_end,  &
                                bb%bnd_x1,   bb%bnd_x2, bb%bnd_y1,   bb%bnd_y2,  &
                                grid_data   %lu        %block(k)%field,  &
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
        call end_kernel_timer('hh_init_kernel')

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
        type(block_bounds_type) :: bb

        call start_kernel_timer('hh_update_kernel')
        do k = 1, domain%bcount

            call bb%set(domain, k)

            call hh_update_kernel(bb%nx_start, bb%nx_end, bb%ny_start, bb%ny_end,  &
                                  bb%bnd_x1,   bb%bnd_x2, bb%bnd_y1,   bb%bnd_y2,  &
                                  grid_data   %lu        %block(k)%field,  & 
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
        call end_kernel_timer('hh_update_kernel')

        call sync(domain, grid_data%hhu_n)
        call sync(domain, grid_data%hhv_n)
        call sync(domain, grid_data%hhh_n)
    
    end subroutine

    subroutine envoke_hh_shift_kernel(domain, grid_data)

        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data

        integer :: k
        type(block_bounds_type) :: bb

        call start_kernel_timer('hh_shift_kernel')
        do k = 1, domain%bcount

            call bb%set(domain, k)

            call hh_shift_kernel(bb%nx_start, bb%nx_end, bb%ny_start, bb%ny_end,  &
                                 bb%bnd_x1,   bb%bnd_x2, bb%bnd_y1,   bb%bnd_y2,  &
                                 grid_data  %lu     %block(k)%field,  &
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
        call end_kernel_timer('hh_shift_kernel')

    end subroutine

    subroutine envoke_uv_trans_vort_kernel(domain, grid_data, u, v, vort)

        type(domain_type), intent(in) :: domain
        type(grid_type), intent(in) :: grid_data
        type(data2D_real8_type), intent(in) :: u, v
        type(data2D_real8_type), intent(inout) :: vort

        ! Interface only for 2D call
        integer, parameter :: nlev = 1

        integer :: k
        type(block_bounds_type) :: bb

        call start_kernel_timer('uv_trans_vort_kernel')
        do k = 1, domain%bcount
            
            call bb%set(domain, k)

            call uv_trans_vort_kernel(bb%nx_start, bb%nx_end, bb%ny_start, bb%ny_end,  &
                                      bb%bnd_x1,   bb%bnd_x2, bb%bnd_y1,   bb%bnd_y2,  &
                                      grid_data  %luu   %block(k)%field,  &
                                      grid_data  %dxt   %block(k)%field,  &
                                      grid_data  %dyt   %block(k)%field,  &
                                      grid_data  %dxb   %block(k)%field,  &
                                      grid_data  %dyb   %block(k)%field,  &
                                                  u     %block(k)%field,  &
                                                  v     %block(k)%field,  &
                                                  vort  %block(k)%field,  &
                                      nlev)

        enddo
        call end_kernel_timer('uv_trans_vort_kernel')

        call sync(domain, vort)

    end subroutine

    subroutine envoke_uv_trans_kernel(domain, grid_data, u, v, vort, RHSx, RHSy)

        type(domain_type), intent(in) :: domain
        type(grid_type), intent(in) :: grid_data
        type(data2D_real8_type), intent(in) :: u, v
        type(data2D_real8_type), intent(in) :: vort
        type(data2D_real8_type), intent(inout) :: RHSx, RHSy
        
        ! Interface only for 2D call
        integer, parameter :: nlev = 1

        integer :: k
        type(block_bounds_type) :: bb

        call start_kernel_timer('uv_trans_kernel')
        do k = 1, domain%bcount
            
            call bb%set(domain, k)

            call uv_trans_kernel(bb%nx_start, bb%nx_end, bb%ny_start, bb%ny_end,  &
                                 bb%bnd_x1,   bb%bnd_x2, bb%bnd_y1,   bb%bnd_y2,  &
                                 grid_data  %lcu   %block(k)%field,  &
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
        call end_kernel_timer('uv_trans_kernel')
    
    end subroutine

    subroutine envoke_uv_diff2_kernel(domain, grid_data, mu, str_t, str_s, RHSx, RHSy)

        type(domain_type), intent(in) :: domain
        type(grid_type), intent(in) :: grid_data
        type(data2D_real8_type), intent(in) :: mu, str_t, str_s
        type(data2D_real8_type), intent(inout) :: RHSx, RHSy
        
        ! Interface only for 2D call
        integer, parameter :: nlev = 1

        integer :: k
        type(block_bounds_type) :: bb

        call start_kernel_timer('uv_diff2_kernel')
        do k = 1, domain%bcount
            
            call bb%set(domain, k)
            
            call uv_diff2_kernel(bb%nx_start, bb%nx_end, bb%ny_start, bb%ny_end,  &
                                 bb%bnd_x1,   bb%bnd_x2, bb%bnd_y1,   bb%bnd_y2,  &
                                 grid_data  %lcu    %block(k)%field,  &
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
        call end_kernel_timer('uv_diff2_kernel')
    
    end subroutine

    subroutine envoke_sw_update_ssh_kernel(domain, grid_data, tau, sshn, sshp, ubrtr, vbrtr)

        type(domain_type), intent(in) :: domain
        type(grid_type), intent(in) :: grid_data
        real(wp8), intent(in) :: tau
        type(data2D_real8_type), intent(in) :: ubrtr, vbrtr, sshp
        type(data2D_real8_type), intent(inout) :: sshn

        integer :: k
        type(block_bounds_type) :: bb

        call start_kernel_timer('sw_update_ssh_kernel')

        !$omp parallel do default(shared)    &
        !$omp private(k, bb)
        do k = 1, domain%bcount
            
            call bb%set(domain, k)
            
            call sw_update_ssh_kernel(bb%nx_start, bb%nx_end, bb%ny_start, bb%ny_end,  &
                                      bb%bnd_x1,   bb%bnd_x2, bb%bnd_y1,   bb%bnd_y2,  &
                                      tau,                                 &
                                      grid_data  %lu     %block(k)%field,  & 
                                      grid_data  %dx     %block(k)%field,  & 
                                      grid_data  %dy     %block(k)%field,  & 
                                      grid_data  %dxh    %block(k)%field,  & 
                                      grid_data  %dyh    %block(k)%field,  & 
                                      grid_data  %hhu    %block(k)%field,  & 
                                      grid_data  %hhv    %block(k)%field,  & 
                                                  sshn   %block(k)%field,  & 
                                                  sshp   %block(k)%field,  & 
                                                  ubrtr  %block(k)%field,  & 
                                                  vbrtr  %block(k)%field)

        enddo
        !$omp end parallel do
        call end_kernel_timer('sw_update_ssh_kernel')

        call sync(domain, sshn)
    
    end subroutine

    subroutine envoke_sw_update_uv(domain, grid_data, tau, ssh, ubrtr, ubrtrn, ubrtrp, vbrtr, vbrtrn, vbrtrp, rdis, &
                                   RHSx, RHSy, RHSx_adv, RHSy_adv, RHSx_dif, RHSy_dif)

        type(domain_type), intent(in) :: domain
        type(grid_type), intent(in) :: grid_data
        real(wp8), intent(in) :: tau
        type(data2D_real4_type), intent(in) :: rdis
        type(data2D_real8_type), intent(in) :: ssh, ubrtr, ubrtrp, vbrtr, vbrtrp,  &
                                               RHSx, RHSy, RHSx_adv, RHSy_adv, RHSx_dif, RHSy_dif
        type(data2D_real8_type), intent(inout) :: ubrtrn, vbrtrn

        integer :: k
        type(block_bounds_type) :: bb

        call start_kernel_timer('sw_update_uv')
        do k = 1, domain%bcount
            
            call bb%set(domain, k)
            
            call sw_update_uv(bb%nx_start, bb%nx_end, bb%ny_start, bb%ny_end,  &
                              bb%bnd_x1,   bb%bnd_x2, bb%bnd_y1,   bb%bnd_y2,  &
                              tau,                                    &
                              grid_data  %lcu       %block(k)%field,  &
                              grid_data  %lcv       %block(k)%field,  &
                              grid_data  %dxt       %block(k)%field,  &
                              grid_data  %dyt       %block(k)%field,  &
                              grid_data  %dxh       %block(k)%field,  &
                              grid_data  %dyh       %block(k)%field,  &
                              grid_data  %dxb       %block(k)%field,  &
                              grid_data  %dyb       %block(k)%field,  &
                              grid_data  %hhu       %block(k)%field,  &
                              grid_data  %hhu_n     %block(k)%field,  &
                              grid_data  %hhu_p     %block(k)%field,  &
                              grid_data  %hhv       %block(k)%field,  &
                              grid_data  %hhv_n     %block(k)%field,  &
                              grid_data  %hhv_p     %block(k)%field,  &
                              grid_data  %hhh       %block(k)%field,  &
                                          ssh       %block(k)%field,  &
                                          ubrtr     %block(k)%field,  &
                                          ubrtrn    %block(k)%field,  &
                                          ubrtrp    %block(k)%field,  &
                                          vbrtr     %block(k)%field,  &
                                          vbrtrn    %block(k)%field,  &
                                          vbrtrp    %block(k)%field,  &
                                          rdis      %block(k)%field,  &
                              grid_data  %rlh_s     %block(k)%field,  &
                                          RHSx      %block(k)%field,  &
                                          RHSy      %block(k)%field,  &
                                          RHSx_adv  %block(k)%field,  &
                                          RHSy_adv  %block(k)%field,  &
                                          RHSx_dif  %block(k)%field,  &
                                          RHSy_dif  %block(k)%field)

        enddo
        call end_kernel_timer('sw_update_uv')

        call sync(domain, ubrtrn)
        call sync(domain, vbrtrn)
    
    end subroutine

    subroutine envoke_sw_next_step(domain, grid_data, time_smooth, ssh, sshn, sshp, ubrtr, ubrtrn, ubrtrp, vbrtr, vbrtrn, vbrtrp)

        type(domain_type), intent(in) :: domain
        type(grid_type), intent(in) :: grid_data
        real(wp8), intent(in) :: time_smooth
        type(data2D_real8_type), intent(in) :: sshn, ubrtrn, vbrtrn
        type(data2D_real8_type), intent(inout) :: ssh, sshp, ubrtr, ubrtrp, vbrtr, vbrtrp

        integer :: k
        type(block_bounds_type) :: bb

        call start_kernel_timer('sw_next_step')
        do k = 1, domain%bcount
            
            call bb%set(domain, k)
            
            call sw_next_step(bb%nx_start, bb%nx_end, bb%ny_start, bb%ny_end,  &
                              bb%bnd_x1,   bb%bnd_x2, bb%bnd_y1,   bb%bnd_y2,  &
                              time_smooth,                          &
                              grid_data  %lu      %block(k)%field,  &
                              grid_data  %lcu     %block(k)%field,  &
                              grid_data  %lcv     %block(k)%field,  &
                                          ssh     %block(k)%field,  &
                                          sshn    %block(k)%field,  &
                                          sshp    %block(k)%field,  &
                                          ubrtr   %block(k)%field,  &
                                          ubrtrn  %block(k)%field,  &
                                          ubrtrp  %block(k)%field,  &
                                          vbrtr   %block(k)%field,  &
                                          vbrtrn  %block(k)%field,  &
                                          vbrtrp  %block(k)%field)
        enddo
        call end_kernel_timer('sw_next_step')

    end subroutine

end module shallow_water_interface_module