#include "macros/mpp_macros.fi"

module tracer_control_module
    
    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use decomposition_module, only: domain_type, domain => domain_data
    use ocean_module, only: ocean_type, ocean_data
    use grid_module, only: grid_type, grid_data
    use config_sw_module, only: use_tracers, tracer_num, time_smooth
    use mpp_module
#ifdef _GPU_MODE_
    use kernel_interface_module, only: envoke_empty_kernel, envoke_empty_sync, kernel_parameters_type, envoke => envoke_device
    use tracer_interface_gpu_module, only: envoke_tran_diff_fluxes_kernel => envoke_tran_diff_fluxes_kernel_gpu,  &
                                           envoke_tran_diff_fluxes_sync   => envoke_tran_diff_fluxes_sync_gpu,    &
                                           envoke_tran_diff_tracer_kernel => envoke_tran_diff_tracer_kernel_gpu,  &
                                           envoke_tran_diff_tracer_sync   => envoke_tran_diff_tracer_sync_gpu,    &
                                           envoke_tracer_next_step_kernel => envoke_tracer_next_step_kernel_gpu,  &
                                           envoke_tracer_next_step_sync   => envoke_tracer_next_step_sync_gpu
#else
    use kernel_interface_module, only: envoke_empty_kernel, envoke_empty_sync, kernel_parameters_type, envoke
    use tracer_interface_module, only: envoke_tran_diff_fluxes_kernel, envoke_tran_diff_fluxes_sync,  &
                                       envoke_tran_diff_tracer_kernel, envoke_tran_diff_tracer_sync,  &
                                       envoke_tracer_next_step_kernel, envoke_tracer_next_step_sync
#endif

    implicit none
    save
    public

contains

!===============================================================================
    subroutine expl_tracer(tau)
        real(wp8), intent(in) :: tau

        procedure(envoke_empty_kernel), pointer :: sub_kernel
        procedure(envoke_empty_sync), pointer :: sub_sync
        type(kernel_parameters_type) :: kernel_parameters
        integer :: k

        if (use_tracers > 0) then
            do k = 1, tracer_num
                ! Set kernel parameters for tracer kernels
                call kernel_parameters%clear()
                kernel_parameters%tau = tau
                kernel_parameters%time_smooth = time_smooth
                kernel_parameters%data_id = k

                sub_kernel => envoke_tran_diff_fluxes_kernel
                sub_sync   => envoke_tran_diff_fluxes_sync
                call envoke(sub_kernel, sub_sync, kernel_parameters)

                sub_kernel => envoke_tran_diff_tracer_kernel
                sub_sync   => envoke_tran_diff_tracer_sync
                call envoke(sub_kernel, sub_sync, kernel_parameters)

                sub_kernel => envoke_tracer_next_step_kernel
                sub_sync   => envoke_tracer_next_step_sync
                call envoke(sub_kernel, sub_sync, kernel_parameters)
            enddo
        endif
    endsubroutine expl_tracer

endmodule