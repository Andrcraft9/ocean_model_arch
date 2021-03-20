#include "macros/mpp_macros.fi"
#ifdef _GPU_MODE_

module shallow_water_interface_gpu_module

    use kernel_interface_module
    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use decomposition_module, only: domain_type, domain => domain_data
    use data_types_module, only: data2D_real8_type, data2D_real4_type
    use ocean_module, only: ocean_type, ocean_data
    use grid_module, only: grid_type, grid_data
    use mpp_sync_module, only: hybrid_sync, sync, sync_parameters_type
    use mixing_gpu_module, only: stress_components_gpu_kernel
    use depth_gpu_module, only: hh_init_gpu_kernel, hh_update_gpu_kernel, hh_shift_gpu_kernel
    use velssh_sw_gpu_module, only: sw_update_ssh_gpu_kernel, sw_update_uv_gpu_kernel, sw_next_step_gpu_kernel, uv_trans_vort_gpu_kernel, uv_trans_gpu_kernel, uv_diff2_gpu_kernel

    implicit none
    save
    private

    public :: envoke_stress_components_gpu_kernel, envoke_stress_components_gpu_sync
    public :: envoke_hh_init_gpu_kernel,           envoke_hh_init_gpu_sync
    public :: envoke_hh_update_gpu_kernel,         envoke_hh_update_gpu_sync
    public :: envoke_hh_shift_gpu_kernel,          envoke_hh_shift_gpu_sync
    public :: envoke_uv_trans_vort_gpu_kernel,     envoke_uv_trans_vort_gpu_sync
    public :: envoke_uv_trans_gpu_kernel,          envoke_uv_trans_gpu_sync
    public :: envoke_uv_diff2_gpu_kernel,          envoke_uv_diff2_gpu_sync
    public :: envoke_sw_update_ssh_gpu_kernel,     envoke_sw_update_ssh_gpu_sync
    public :: envoke_sw_update_uv_gpu_kernel,      envoke_sw_update_uv_gpu_sync
    public :: envoke_sw_next_step_gpu_kernel,      envoke_sw_next_step_gpu_sync

contains

!-----------------------------------------------------------------------------!
!---------------------------- GPU Kernels ------------------------------------!
!-----------------------------------------------------------------------------!
    subroutine envoke_hh_init_kernel(k, param)
        integer, intent(in) :: k
        type(kernel_parameters_type), intent(in) :: param

        call hh_init_kernel(domain%bnx_start(k), domain%bnx_end(k), domain%bny_start(k), domain%bny_end(k),  &
                            domain%bbnd_x1(k),   domain%bbnd_x2(k), domain%bbnd_y1(k),   domain%bbnd_y2(k),  &
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
                            ocean_data  %ssh       %block(k)%field,  &
                            ocean_data  %sshp      %block(k)%field,  &
                            grid_data   %hhq_rest  %block(k)%field)
    end subroutine

    subroutine envoke_hh_init_sync(k, sync_parameters)
        integer, intent(in) :: k
        type(sync_parameters_type), intent(in) :: sync_parameters

        call hybrid_sync(sync_parameters, 1, domain, grid_data%hhu)
        call hybrid_sync(sync_parameters, 2, domain, grid_data%hhu_p)
        call hybrid_sync(sync_parameters, 3, domain, grid_data%hhu_n)
        call hybrid_sync(sync_parameters, 4, domain, grid_data%hhv)
        call hybrid_sync(sync_parameters, 5, domain, grid_data%hhv_p)
        call hybrid_sync(sync_parameters, 6, domain, grid_data%hhv_n)
        call hybrid_sync(sync_parameters, 7, domain, grid_data%hhh)
        call hybrid_sync(sync_parameters, 8, domain, grid_data%hhh_p)
        call hybrid_sync(sync_parameters, 9, domain, grid_data%hhh_n)
    end subroutine

!-----------------------------------------------------------------------------!
    subroutine envoke_stress_components_kernel(k, param)
        integer, intent(in) :: k
        type(kernel_parameters_type), intent(in) :: param

        ! Interface only for 2D call
        integer, parameter :: nlev = 1

        call stress_components_kernel(domain%bnx_start(k), domain%bnx_end(k), domain%bny_start(k), domain%bny_end(k),  &
                                      domain%bbnd_x1(k),   domain%bbnd_x2(k), domain%bbnd_y1(k),   domain%bbnd_y2(k),  &
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
                                      ocean_data %ubrtrp %block(k)%field,  & 
                                      ocean_data %vbrtrp %block(k)%field,  & 
                                      ocean_data %str_t  %block(k)%field,  &
                                      ocean_data %str_s  %block(k)%field,  &
                                      nlev)
    end subroutine

    subroutine envoke_stress_components_sync(k, sync_parameters)
        integer, intent(in) :: k
        type(sync_parameters_type), intent(in) :: sync_parameters

        call hybrid_sync(sync_parameters, 1, domain, ocean_data%str_t)
        call hybrid_sync(sync_parameters, 2, domain, ocean_data%str_s)
    end subroutine

!-----------------------------------------------------------------------------!
    subroutine envoke_hh_update_kernel(k, param)
        integer, intent(in) :: k
        type(kernel_parameters_type), intent(in) :: param

        call hh_update_kernel(domain%bnx_start(k), domain%bnx_end(k), domain%bny_start(k), domain%bny_end(k),  &
                              domain%bbnd_x1(k),   domain%bbnd_x2(k), domain%bbnd_y1(k),   domain%bbnd_y2(k),  &
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
                              ocean_data  %ssh       %block(k)%field,  & 
                              grid_data   %hhq_rest  %block(k)%field)
    end subroutine

    subroutine envoke_hh_update_sync(k, sync_parameters)
        integer, intent(in) :: k
        type(sync_parameters_type), intent(in) :: sync_parameters

        call hybrid_sync(sync_parameters, 1, domain, grid_data%hhu_n)
        call hybrid_sync(sync_parameters, 2, domain, grid_data%hhv_n)
        call hybrid_sync(sync_parameters, 3, domain, grid_data%hhh_n)
    end subroutine

!-----------------------------------------------------------------------------!
    subroutine envoke_hh_shift_kernel(k, param)
        integer, intent(in) :: k
        type(kernel_parameters_type), intent(in) :: param

        call hh_shift_kernel(domain%bnx_start(k), domain%bnx_end(k), domain%bny_start(k), domain%bny_end(k),  &
                             domain%bbnd_x1(k),   domain%bbnd_x2(k), domain%bbnd_y1(k),   domain%bbnd_y2(k),  &
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
    end subroutine

    subroutine envoke_hh_shift_sync(k, sync_parameters)
        integer, intent(in) :: k
        type(sync_parameters_type), intent(in) :: sync_parameters
    end subroutine

!-----------------------------------------------------------------------------!
    subroutine envoke_uv_trans_vort_kernel(k, param)
        integer, intent(in) :: k
        type(kernel_parameters_type), intent(in) :: param

        ! Interface only for 2D call
        integer, parameter :: nlev = 1

        call uv_trans_vort_kernel(domain%bnx_start(k), domain%bnx_end(k), domain%bny_start(k), domain%bny_end(k),  &
                                  domain%bbnd_x1(k),   domain%bbnd_x2(k), domain%bbnd_y1(k),   domain%bbnd_y2(k),  &
                                  grid_data  %luu   %block(k)%field,  &
                                  grid_data  %dxt   %block(k)%field,  &
                                  grid_data  %dyt   %block(k)%field,  &
                                  grid_data  %dxb   %block(k)%field,  &
                                  grid_data  %dyb   %block(k)%field,  &
                                  ocean_data %ubrtr %block(k)%field,  &
                                  ocean_data %vbrtr %block(k)%field,  &
                                  ocean_data %vort  %block(k)%field,  &
                                  nlev)
    end subroutine

    subroutine envoke_uv_trans_vort_sync(k, sync_parameters)
        integer, intent(in) :: k
        type(sync_parameters_type), intent(in) :: sync_parameters

        call hybrid_sync(sync_parameters, 1, domain, ocean_data%vort)
    end subroutine

!-----------------------------------------------------------------------------!
    subroutine envoke_uv_trans_kernel(k, param)
        integer, intent(in) :: k
        type(kernel_parameters_type), intent(in) :: param
        
        ! Interface only for 2D call
        integer, parameter :: nlev = 1

        call uv_trans_kernel(domain%bnx_start(k), domain%bnx_end(k), domain%bny_start(k), domain%bny_end(k),  &
                             domain%bbnd_x1(k),   domain%bbnd_x2(k), domain%bbnd_y1(k),   domain%bbnd_y2(k),  &
                             grid_data  %lcu   %block(k)%field,  &
                             grid_data  %lcv   %block(k)%field,  &
                             grid_data  %luu   %block(k)%field,  &
                             grid_data  %dxh   %block(k)%field,  &
                             grid_data  %dyh   %block(k)%field,  &
                             ocean_data %ubrtr %block(k)%field,  &
                             ocean_data %vbrtr %block(k)%field,  &
                             ocean_data %vort  %block(k)%field,  &
                             grid_data  %hhq   %block(k)%field,  &
                             grid_data  %hhu   %block(k)%field,  &
                             grid_data  %hhv   %block(k)%field,  &
                             grid_data  %hhh   %block(k)%field,  &
                             ocean_data %RHSx_adv  %block(k)%field,  &
                             ocean_data %RHSy_adv  %block(k)%field,  &
                             nlev)
    end subroutine

    subroutine envoke_uv_trans_sync(k, sync_parameters)
        integer, intent(in) :: k
        type(sync_parameters_type), intent(in) :: sync_parameters
    end subroutine

!-----------------------------------------------------------------------------!
    subroutine envoke_uv_diff2_kernel(k, param)
        integer, intent(in) :: k
        type(kernel_parameters_type), intent(in) :: param
        
        ! Interface only for 2D call
        integer, parameter :: nlev = 1
        
            call uv_diff2_kernel(domain%bnx_start(k), domain%bnx_end(k), domain%bny_start(k), domain%bny_end(k),  &
                                 domain%bbnd_x1(k),   domain%bbnd_x2(k), domain%bbnd_y1(k),   domain%bbnd_y2(k),  &
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
                                 ocean_data %mu     %block(k)%field,  &
                                 ocean_data %str_t  %block(k)%field,  &
                                 ocean_data %str_s  %block(k)%field,  &
                                 grid_data  %hhq    %block(k)%field,  &
                                 grid_data  %hhu    %block(k)%field,  &
                                 grid_data  %hhv    %block(k)%field,  &
                                 grid_data  %hhh    %block(k)%field,  &
                                 ocean_data %RHSx_dif %block(k)%field,  &
                                 ocean_data %RHSy_dif %block(k)%field,  &
                                 nlev)
    end subroutine

    subroutine envoke_uv_diff2_sync(k, sync_parameters)
        integer, intent(in) :: k
        type(sync_parameters_type), intent(in) :: sync_parameters
    end subroutine

!-----------------------------------------------------------------------------!
    subroutine envoke_sw_update_ssh_kernel(k, param)
        integer, intent(in) :: k
        type(kernel_parameters_type), intent(in) :: param

        real(wp8) :: tau
        tau = param%param_real8

        call sw_update_ssh_kernel(domain%bnx_start(k), domain%bnx_end(k), domain%bny_start(k), domain%bny_end(k),  &
                                  domain%bbnd_x1(k),   domain%bbnd_x2(k), domain%bbnd_y1(k),   domain%bbnd_y2(k),  &
                                  tau,                                 &
                                  grid_data  %lu     %block(k)%field,  & 
                                  grid_data  %dx     %block(k)%field,  & 
                                  grid_data  %dy     %block(k)%field,  & 
                                  grid_data  %dxh    %block(k)%field,  & 
                                  grid_data  %dyh    %block(k)%field,  & 
                                  grid_data  %hhu    %block(k)%field,  & 
                                  grid_data  %hhv    %block(k)%field,  & 
                                  ocean_data %sshn   %block(k)%field,  & 
                                  ocean_data %sshp   %block(k)%field,  & 
                                  ocean_data %ubrtr  %block(k)%field,  & 
                                  ocean_data %vbrtr  %block(k)%field)
    end subroutine

    subroutine envoke_sw_update_ssh_sync(k, sync_parameters)
        integer, intent(in) :: k
        type(sync_parameters_type), intent(in) :: sync_parameters

        call hybrid_sync(sync_parameters, 1, domain, ocean_data%sshn)
    end subroutine

!-----------------------------------------------------------------------------!
    subroutine envoke_sw_update_uv_kernel(k, param)
        integer, intent(in) :: k
        type(kernel_parameters_type), intent(in) :: param

        real(wp8) :: tau
        tau = param%param_real8

        call sw_update_uv(domain%bnx_start(k), domain%bnx_end(k), domain%bny_start(k), domain%bny_end(k),  &
                          domain%bbnd_x1(k),   domain%bbnd_x2(k), domain%bbnd_y1(k),   domain%bbnd_y2(k),  &
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
                          ocean_data %ssh       %block(k)%field,  &
                          ocean_data %ubrtr     %block(k)%field,  &
                          ocean_data %ubrtrn    %block(k)%field,  &
                          ocean_data %ubrtrp    %block(k)%field,  &
                          ocean_data %vbrtr     %block(k)%field,  &
                          ocean_data %vbrtrn    %block(k)%field,  &
                          ocean_data %vbrtrp    %block(k)%field,  &
                          ocean_data %r_diss    %block(k)%field,  &
                          grid_data  %rlh_s     %block(k)%field,  &
                          ocean_data %RHSx      %block(k)%field,  &
                          ocean_data %RHSy      %block(k)%field,  &
                          ocean_data %RHSx_adv  %block(k)%field,  &
                          ocean_data %RHSy_adv  %block(k)%field,  &
                          ocean_data %RHSx_dif  %block(k)%field,  &
                          ocean_data %RHSy_dif  %block(k)%field)
    end subroutine

    subroutine envoke_sw_update_uv_sync(k, sync_parameters)
        integer, intent(in) :: k
        type(sync_parameters_type), intent(in) :: sync_parameters

        call hybrid_sync(sync_parameters, 1, domain, ocean_data%vbrtrn)
        call hybrid_sync(sync_parameters, 2, domain, ocean_data%ubrtrn)
    end subroutine

!-----------------------------------------------------------------------------!
    subroutine envoke_sw_next_step_kernel(k, param)
        integer, intent(in) :: k
        type(kernel_parameters_type), intent(in) :: param

        real(wp8) :: time_smooth
        time_smooth = param%param_real8

        call sw_next_step(domain%bnx_start(k), domain%bnx_end(k), domain%bny_start(k), domain%bny_end(k),  &
                          domain%bbnd_x1(k),   domain%bbnd_x2(k), domain%bbnd_y1(k),   domain%bbnd_y2(k),  &
                          time_smooth,                          &
                          grid_data  %lu      %block(k)%field,  &
                          grid_data  %lcu     %block(k)%field,  &
                          grid_data  %lcv     %block(k)%field,  &
                          ocean_data %ssh     %block(k)%field,  &
                          ocean_data %sshn    %block(k)%field,  &
                          ocean_data %sshp    %block(k)%field,  &
                          ocean_data %ubrtr   %block(k)%field,  &
                          ocean_data %ubrtrn  %block(k)%field,  &
                          ocean_data %ubrtrp  %block(k)%field,  &
                          ocean_data %vbrtr   %block(k)%field,  &
                          ocean_data %vbrtrn  %block(k)%field,  &
                          ocean_data %vbrtrp  %block(k)%field)
    end subroutine

    subroutine envoke_sw_next_step_sync(k, sync_parameters)
        integer, intent(in) :: k
        type(sync_parameters_type), intent(in) :: sync_parameters
    end subroutine

end module shallow_water_interface_gpu_module

#endif