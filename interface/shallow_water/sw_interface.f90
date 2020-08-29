module shallow_water_interface_module

    use kernel_interface_module
    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use mpp_module
    use decomposition_module, only: domain_type
    use data_types_module, only: data2D_real8_type, data2D_real4_type
    use ocean_module, only: ocean_type
    use grid_module, only: grid_type
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

    public :: envoke, envoke_empty_kernel, envoke_empty_sync

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
!-------------------------- Interface subroutines ----------------------------!
    subroutine envoke_empty_kernel(k, domain, grid_data, ocean_data, param)
        integer, intent(in) :: k
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data
        real(wp8), intent(in) :: param
    end subroutine 

    subroutine envoke_empty_sync(sync_parameters, domain, grid_data, ocean_data)
        type(sync_parameters_type), intent(in) :: sync_parameters
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data
    end subroutine 

    subroutine envoke(domain, grid_data, ocean_data, sub_kernel, sub_sync, param)
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data
        procedure(envoke_empty_kernel), pointer :: sub_kernel
        procedure(envoke_empty_sync), pointer :: sub_sync
        real(wp8), intent(in) :: param
    
        integer :: k
        type(sync_parameters_type) :: sync_parameters_inner, sync_parameters_boundary

        sync_parameters_inner%is_inner_sync = .true.
        sync_parameters_boundary%is_inner_sync = .false.

#ifdef _MPP_HYBRID_BLOCK_MODE_

        !$omp parallel default(shared)
    
        !$omp do private(k) schedule(static, 1)
        do k = domain%start_boundary, domain%start_boundary + domain%bcount_boundary - 1
            call sub_kernel(k, domain, grid_data, ocean_data, param)
        enddo
        !$omp end do

        !$omp master
        call sub_sync(sync_parameters_boundary, domain, grid_data, ocean_data)
        !$omp end master

        !$omp do private(k) schedule(static, 1)
        do k = domain%start_inner, domain%start_inner + domain%bcount_inner - 1
            call sub_kernel(k, domain, grid_data, ocean_data, param)
        enddo
        !$omp end do

        call sub_sync(sync_parameters_inner, domain, grid_data, ocean_data)

        !$omp end parallel
#else

    _OMP_BLOCKS_PARALLEL_BEGIN_
    do k = 1, domain%bcount
        call sub_kernel(k, domain, grid_data, ocean_data, param)
    enddo
    _OMP_BLOCKS_PARALLEL_END_

    call sub_sync(sync_parameters_boundary, domain, grid_data, ocean_data)
    call sub_sync(sync_parameters_inner, domain, grid_data, ocean_data)

#endif
    end subroutine

!-----------------------------------------------------------------------------!
!------------------------------- Kernels -------- ----------------------------!

!-----------------------------------------------------------------------------!
    subroutine envoke_hh_init_kernel(k, domain, grid_data, ocean_data, param)
        integer, intent(in) :: k
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data
        real(wp8), intent(in) :: param

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

    subroutine envoke_hh_init_sync(sync_parameters, domain, grid_data, ocean_data)
        type(sync_parameters_type), intent(in) :: sync_parameters
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data

        call hybrid_sync(sync_parameters, domain, grid_data%hhu)
        call hybrid_sync(sync_parameters, domain, grid_data%hhu_p)
        call hybrid_sync(sync_parameters, domain, grid_data%hhu_n)
        call hybrid_sync(sync_parameters, domain, grid_data%hhv)
        call hybrid_sync(sync_parameters, domain, grid_data%hhv_p)
        call hybrid_sync(sync_parameters, domain, grid_data%hhv_n)
        call hybrid_sync(sync_parameters, domain, grid_data%hhh)
        call hybrid_sync(sync_parameters, domain, grid_data%hhh_p)
        call hybrid_sync(sync_parameters, domain, grid_data%hhh_n)
    end subroutine

!-----------------------------------------------------------------------------!
    subroutine envoke_check_ssh_err_kernel(k, domain, grid_data, ocean_data, param)
        integer, intent(in) :: k
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data
        real(wp8), intent(in) :: param

        call check_ssh_err_kernel(domain%bnx_start(k), domain%bnx_end(k), domain%bny_start(k), domain%bny_end(k),  &
                                  domain%bbnd_x1(k),   domain%bbnd_x2(k), domain%bbnd_y1(k),   domain%bbnd_y2(k),  &
                                  grid_data  %lu   %block(k)%field,  &
                                  ocean_data %ssh  %block(k)%field,  &
                                  'ssh')
    end subroutine

    subroutine envoke_check_ssh_err_sync(sync_parameters, domain, grid_data, ocean_data)
        type(sync_parameters_type), intent(in) :: sync_parameters
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data

    end subroutine

!-----------------------------------------------------------------------------!
    subroutine envoke_stress_components_kernel(k, domain, grid_data, ocean_data, param)
        integer, intent(in) :: k
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data
        real(wp8), intent(in) :: param

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

    subroutine envoke_stress_components_sync(sync_parameters, domain, grid_data, ocean_data)
        type(sync_parameters_type), intent(in) :: sync_parameters
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data

        call hybrid_sync(sync_parameters, domain, ocean_data%str_t)
        call hybrid_sync(sync_parameters, domain, ocean_data%str_s)
    end subroutine

!-----------------------------------------------------------------------------!
    subroutine envoke_hh_update_kernel(k, domain, grid_data, ocean_data, param)
        integer, intent(in) :: k
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data
        real(wp8), intent(in) :: param

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

    subroutine envoke_hh_update_sync(sync_parameters, domain, grid_data, ocean_data)
        type(sync_parameters_type), intent(in) :: sync_parameters
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data

        call hybrid_sync(sync_parameters, domain, grid_data%hhu_n)
        call hybrid_sync(sync_parameters, domain, grid_data%hhv_n)
        call hybrid_sync(sync_parameters, domain, grid_data%hhh_n)
    end subroutine

!-----------------------------------------------------------------------------!
    subroutine envoke_hh_shift_kernel(k, domain, grid_data, ocean_data, param)
        integer, intent(in) :: k
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data
        real(wp8), intent(in) :: param

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

    subroutine envoke_hh_shift_sync(sync_parameters, domain, grid_data, ocean_data)
        type(sync_parameters_type), intent(in) :: sync_parameters
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data

    end subroutine

!-----------------------------------------------------------------------------!
    subroutine envoke_uv_trans_vort_kernel(k, domain, grid_data, ocean_data, param)
        integer, intent(in) :: k
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data
        real(wp8), intent(in) :: param

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

    subroutine envoke_uv_trans_vort_sync(sync_parameters, domain, grid_data, ocean_data)
        type(sync_parameters_type), intent(in) :: sync_parameters
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data

        call hybrid_sync(sync_parameters, domain, ocean_data%vort)
    end subroutine

!-----------------------------------------------------------------------------!
    subroutine envoke_uv_trans_kernel(k, domain, grid_data, ocean_data, param)
        integer, intent(in) :: k
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data
        real(wp8), intent(in) :: param
        
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

    subroutine envoke_uv_trans_sync(sync_parameters, domain, grid_data, ocean_data)
        type(sync_parameters_type), intent(in) :: sync_parameters
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data

    end subroutine

!-----------------------------------------------------------------------------!
    subroutine envoke_uv_diff2_kernel(k, domain, grid_data, ocean_data, param)
        integer, intent(in) :: k
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data
        real(wp8), intent(in) :: param
        
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

    subroutine envoke_uv_diff2_sync(sync_parameters, domain, grid_data, ocean_data)
        type(sync_parameters_type), intent(in) :: sync_parameters
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data

    end subroutine

!-----------------------------------------------------------------------------!
    subroutine envoke_sw_update_ssh_kernel(k, domain, grid_data, ocean_data, tau)
        integer, intent(in) :: k
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data
        real(wp8), intent(in) :: tau

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

    subroutine envoke_sw_update_ssh_sync(sync_parameters, domain, grid_data, ocean_data)
        type(sync_parameters_type), intent(in) :: sync_parameters
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data

        call hybrid_sync(sync_parameters, domain, ocean_data%sshn)
    end subroutine

!-----------------------------------------------------------------------------!
    subroutine envoke_sw_update_uv_kernel(k, domain, grid_data, ocean_data, tau)
        integer, intent(in) :: k
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data
        real(wp8), intent(in) :: tau

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

    subroutine envoke_sw_update_uv_sync(sync_parameters, domain, grid_data, ocean_data)
        type(sync_parameters_type), intent(in) :: sync_parameters
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data

        call hybrid_sync(sync_parameters, domain, ocean_data%ubrtrn)
        call hybrid_sync(sync_parameters, domain, ocean_data%vbrtrn)
    end subroutine

!-----------------------------------------------------------------------------!
    subroutine envoke_sw_next_step_kernel(k, domain, grid_data, ocean_data, time_smooth)
        integer, intent(in) :: k
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data
        real(wp8), intent(in) :: time_smooth

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

    subroutine envoke_sw_next_step_sync(sync_parameters, domain, grid_data, ocean_data)
        type(sync_parameters_type), intent(in) :: sync_parameters
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data

    end subroutine

!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!

    subroutine envoke_gaussian_elimination(domain, grid_data, ssh, sigma, nx0, ny0)
        
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(in) :: grid_data
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