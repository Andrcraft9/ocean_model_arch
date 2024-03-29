module shallow_water_interface_module

    use kernel_interface_module
    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use mpp_module
    use decomposition_module, only: domain_type, domain => domain_data
    use data_types_module, only: data2D_real8_type, data2D_real4_type
    use ocean_module, only: ocean_type, ocean_data
    use grid_module, only: grid_type, grid_data
    use mpp_sync_module, only: hybrid_sync, sync, sync_parameters_type
    use mixing_module, only: stress_components_kernel
    use depth_module, only: hh_init_kernel, hh_update_kernel, hh_shift_kernel
    use velssh_sw_module, only: gaussian_elimination_kernel, check_ssh_err_kernel, sw_update_ssh_kernel, sw_update_uv, sw_next_step, uv_trans_vort_kernel, uv_trans_kernel, uv_diff2_kernel
    use errors_module, only: abort_model, check_error

#include "macros/kernel_macros.fi"
#include "macros/mpp_macros.fi"

    implicit none
    save
    private

    public :: envoke_gaussian_elimination

    public :: envoke_check_ssh_err_kernel,     envoke_check_ssh_err_sync
    public :: envoke_stress_components_kernel, envoke_stress_components_sync
    public :: envoke_hh_init_kernel,           envoke_hh_init_sync
    public :: envoke_hh_update_kernel,         envoke_hh_update_sync
    public :: envoke_hh_shift_kernel,          envoke_hh_shift_sync
    public :: envoke_uv_trans_vort_kernel,     envoke_uv_trans_vort_sync
    public :: envoke_uv_trans_kernel,          envoke_uv_trans_sync
    public :: envoke_uv_diff2_kernel,          envoke_uv_diff2_sync
    public :: envoke_sw_update_ssh_kernel,     envoke_sw_update_ssh_sync
    public :: envoke_sw_update_uv_kernel,      envoke_sw_update_uv_sync
    public :: envoke_sw_next_step_kernel,      envoke_sw_next_step_sync

contains

!-----------------------------------------------------------------------------!
!------------------------------- Kernels -------- ----------------------------!
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

        call hybrid_sync(k, sync_parameters, 1, domain, grid_data%hhu)
        call hybrid_sync(k, sync_parameters, 2, domain, grid_data%hhv)
        call hybrid_sync(k, sync_parameters, 3, domain, grid_data%hhh)
        !call hybrid_sync(k, sync_parameters, *, domain, grid_data%hhu_p) ! lazy update
        !call hybrid_sync(k, sync_parameters, *, domain, grid_data%hhu_n) ! not necessary update
        !call hybrid_sync(k, sync_parameters, *, domain, grid_data%hhv_p) ! lazy update
        !call hybrid_sync(k, sync_parameters, *, domain, grid_data%hhv_n) ! not necessary update
        !call hybrid_sync(k, sync_parameters, *, domain, grid_data%hhh_p) ! lazy update
        !call hybrid_sync(k, sync_parameters, *, domain, grid_data%hhh_n) ! not necessary update
    end subroutine

!-----------------------------------------------------------------------------!
    subroutine envoke_check_ssh_err_kernel(k, param)
        integer, intent(in) :: k
        type(kernel_parameters_type), intent(in) :: param

        call check_ssh_err_kernel(domain%bnx_start(k), domain%bnx_end(k), domain%bny_start(k), domain%bny_end(k),  &
                                  domain%bbnd_x1(k),   domain%bbnd_x2(k), domain%bbnd_y1(k),   domain%bbnd_y2(k),  &
                                  grid_data  %lu   %block(k)%field,  &
                                  ocean_data %ssh  %block(k)%field,  &
                                  'ssh')
    end subroutine

    subroutine envoke_check_ssh_err_sync(k, sync_parameters)
        integer, intent(in) :: k
        type(sync_parameters_type), intent(in) :: sync_parameters
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

        call hybrid_sync(k, sync_parameters, 1, domain, ocean_data%str_t)
        call hybrid_sync(k, sync_parameters, 2, domain, ocean_data%str_s)
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

        call hybrid_sync(k, sync_parameters, 1, domain, grid_data%hhu_n)
        call hybrid_sync(k, sync_parameters, 2, domain, grid_data%hhv_n)
        call hybrid_sync(k, sync_parameters, 3, domain, grid_data%hhh_n)
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
        call hybrid_sync(k, sync_parameters, 1, domain, ocean_data%vort)
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
        call hybrid_sync(k, sync_parameters, 1, domain, grid_data%hhu_p) ! Lazy update to hide with kernel computational
        call hybrid_sync(k, sync_parameters, 2, domain, grid_data%hhv_p) ! Lazy update to hide with kernel computational
        call hybrid_sync(k, sync_parameters, 3, domain, grid_data%hhh_p) ! Lazy update to hide with kernel computational
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

        call sw_update_ssh_kernel(domain%bnx_start(k), domain%bnx_end(k), domain%bny_start(k), domain%bny_end(k),  &
                                  domain%bbnd_x1(k),   domain%bbnd_x2(k), domain%bbnd_y1(k),   domain%bbnd_y2(k),  &
                                  param%tau,                           &
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
        call hybrid_sync(k, sync_parameters, 1, domain, ocean_data%sshn)
    end subroutine

!-----------------------------------------------------------------------------!
    subroutine envoke_sw_update_uv_kernel(k, param)
        integer, intent(in) :: k
        type(kernel_parameters_type), intent(in) :: param

        call sw_update_uv(domain%bnx_start(k), domain%bnx_end(k), domain%bny_start(k), domain%bny_end(k),  &
                          domain%bbnd_x1(k),   domain%bbnd_x2(k), domain%bbnd_y1(k),   domain%bbnd_y2(k),  &
                          param%tau,                              &
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
        call hybrid_sync(k, sync_parameters, 1, domain, ocean_data%vbrtrn)
        call hybrid_sync(k, sync_parameters, 2, domain, ocean_data%ubrtrn)
    end subroutine

!-----------------------------------------------------------------------------!
    subroutine envoke_sw_next_step_kernel(k, param)
        integer, intent(in) :: k
        type(kernel_parameters_type), intent(in) :: param

        call sw_next_step(domain%bnx_start(k), domain%bnx_end(k), domain%bny_start(k), domain%bny_end(k),  &
                          domain%bbnd_x1(k),   domain%bbnd_x2(k), domain%bbnd_y1(k),   domain%bbnd_y2(k),  &
                          param%time_smooth,                    &
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

!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!

    subroutine envoke_gaussian_elimination(ssh, sigma, nx0, ny0)
        
        type(data2D_real8_type), intent(inout) :: ssh
        real(wp8), intent(in) ::sigma
        integer, intent(in) :: nx0, ny0

        integer :: k

        _OMP_BLOCKS_PARALLEL_BEGIN_
        do k = 1, domain%bcount
    
            call gaussian_elimination_kernel(domain%bnx_start(k), domain%bnx_end(k), domain%bny_start(k), domain%bny_end(k),  &
                                             domain%bbnd_x1(k),   domain%bbnd_x2(k), domain%bbnd_y1(k),   domain%bbnd_y2(k),  &
                                             grid_data  %lu   %block(k)%field,  &
                                                         ssh  %block(k)%field,  &
                                             sigma, nx0, ny0)
    
        enddo
        _OMP_BLOCKS_PARALLEL_END_

        call sync(domain, ssh)
        
    end subroutine

end module shallow_water_interface_module