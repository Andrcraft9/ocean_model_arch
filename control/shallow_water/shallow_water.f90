#include "macros/mpp_macros.fi"

module shallow_water_module
    
    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use decomposition_module, only: domain_type, domain => domain_data
    use ocean_module, only: ocean_type, ocean_data
    use grid_module, only: grid_type, grid_data
    use config_sw_module, only: full_free_surface, time_smooth, trans_terms, ksw_lat
    use mpp_module
    use kernel_interface_module
    use shallow_water_interface_module

    implicit none
    save
    public

contains

!===============================================================================
    ! explicit shallow water equation sloving
    subroutine expl_shallow_water(tau)

        real(wp8), intent(in) :: tau

        procedure(envoke_empty_kernel), pointer :: sub_kernel
        procedure(envoke_empty_sync), pointer :: sub_sync
        type(kernel_parameters_type) :: kernel_parameters

        ! Set kernel parameters for shallow water kernels
        call kernel_parameters%clear()
        kernel_parameters%tau = tau
        kernel_parameters%time_smooth = time_smooth

        !computing ssh        
        sub_kernel => envoke_sw_update_ssh_kernel
        sub_sync   => envoke_sw_update_ssh_sync
        call envoke(sub_kernel, sub_sync, kernel_parameters)

        if (full_free_surface>0) then
            sub_kernel => envoke_hh_update_kernel
            sub_sync   => envoke_hh_update_sync
            call envoke(sub_kernel, sub_sync, kernel_parameters)
        endif

        !computing advective and lateral-viscous terms for 2d-velocity
        if (trans_terms > 0) then
            sub_kernel => envoke_uv_trans_vort_kernel
            sub_sync   => envoke_uv_trans_vort_sync
            call envoke(sub_kernel, sub_sync, kernel_parameters)

            sub_kernel => envoke_uv_trans_kernel
            sub_sync   => envoke_uv_trans_sync
            call envoke(sub_kernel, sub_sync, kernel_parameters)
        endif

        if (ksw_lat > 0) then
            sub_kernel => envoke_stress_components_kernel
            sub_sync   => envoke_stress_components_sync
            call envoke(sub_kernel, sub_sync, kernel_parameters)

            sub_kernel => envoke_uv_diff2_kernel
            sub_sync   => envoke_uv_diff2_sync
            call envoke(sub_kernel, sub_sync, kernel_parameters)
        endif

        sub_kernel => envoke_sw_update_uv_kernel
        sub_sync   => envoke_sw_update_uv_sync
        call envoke(sub_kernel, sub_sync, kernel_parameters)

        !shifting time indices
        sub_kernel => envoke_sw_next_step_kernel
        sub_sync   => envoke_sw_next_step_sync
        call envoke(sub_kernel, sub_sync, kernel_parameters)

        if(full_free_surface>0) then
            sub_kernel => envoke_hh_shift_kernel
            sub_sync   => envoke_hh_shift_sync
            call envoke(sub_kernel, sub_sync,kernel_parameters)
        endif

        if(full_free_surface>0) then
            !initialize depth for external mode
            sub_kernel => envoke_hh_init_kernel
            sub_sync   => envoke_hh_init_sync
            call envoke(sub_kernel, sub_sync, kernel_parameters)
        endif

        ! Check error
        sub_kernel => envoke_check_ssh_err_kernel
        sub_sync   => envoke_check_ssh_err_sync
        call envoke(sub_kernel, sub_sync, kernel_parameters)

    endsubroutine expl_shallow_water

    subroutine shallow_water_kernels_compute_probe(tau)
        real(wp8), intent(in) :: tau

        procedure(envoke_empty_kernel), pointer :: sub_kernel
        procedure(envoke_empty_sync), pointer :: sub_sync
        type(kernel_parameters_type) :: kernel_parameters

        ! Set kernel parameters for shallow water kernels
        call kernel_parameters%clear()
        kernel_parameters%tau = tau
        kernel_parameters%time_smooth = time_smooth

        ! Compute probe without syncs
        sub_sync  => envoke_empty_sync
        
        ! Below - Body of expl_shallow_water
        ! But it can be just call kernels from sw
        ! We need just to measure time or compute some metrics for kernels

        !computing ssh
        sub_kernel => envoke_sw_update_ssh_kernel
        call envoke(sub_kernel, sub_sync, kernel_parameters)

        if (full_free_surface>0) then
            sub_kernel => envoke_hh_update_kernel
            call envoke(sub_kernel, sub_sync, kernel_parameters)
        endif

        !computing advective and lateral-viscous terms for 2d-velocity
        if (trans_terms > 0) then
            sub_kernel => envoke_uv_trans_vort_kernel
            call envoke(sub_kernel, sub_sync, kernel_parameters)

            sub_kernel => envoke_uv_trans_kernel
            call envoke(sub_kernel, sub_sync, kernel_parameters)
        endif

        if (ksw_lat > 0) then
            sub_kernel => envoke_stress_components_kernel
            call envoke(sub_kernel, sub_sync, kernel_parameters)

            sub_kernel => envoke_uv_diff2_kernel
            call envoke(sub_kernel, sub_sync, kernel_parameters)
        endif

        sub_kernel => envoke_sw_update_uv_kernel
        call envoke(sub_kernel, sub_sync, kernel_parameters)

        !shifting time indices
        sub_kernel => envoke_sw_next_step_kernel
        call envoke(sub_kernel, sub_sync, kernel_parameters)

        if(full_free_surface>0) then
            sub_kernel => envoke_hh_shift_kernel
            call envoke(sub_kernel, sub_sync, kernel_parameters)
        endif

        if(full_free_surface>0) then
            !initialize depth for external mode
            sub_kernel => envoke_hh_init_kernel
            call envoke(sub_kernel, sub_sync, kernel_parameters)
        endif

        ! Check error
        sub_kernel => envoke_check_ssh_err_kernel
        call envoke(sub_kernel, sub_sync, kernel_parameters)

    endsubroutine shallow_water_kernels_compute_probe
endmodule

#ifdef _GPU_MODE_

module shallow_water_gpu_module
    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use decomposition_module, only: domain_type, domain => domain_data
    use ocean_module, only: ocean_type, ocean_data
    use grid_module, only: grid_type, grid_data
    use config_sw_module, only: full_free_surface, time_smooth, trans_terms, ksw_lat
    use mpp_module
    use kernel_interface_module
    use shallow_water_interface_gpu_module

    implicit none
    save
    public

contains

    subroutine expl_shallow_water_gpu(tau)

        real(wp8), intent(in) :: tau

        procedure(envoke_empty_kernel), pointer :: sub_kernel
        procedure(envoke_empty_sync), pointer :: sub_sync
        type(kernel_parameters_type) :: kernel_parameters

        ! Set kernel parameters for shallow water kernels
        call kernel_parameters%clear()
        kernel_parameters%tau = tau
        kernel_parameters%time_smooth = time_smooth

        !computing ssh
        sub_kernel => envoke_sw_update_ssh_kernel_gpu
        sub_sync   => envoke_sw_update_ssh_sync_gpu
        call envoke_device(sub_kernel, sub_sync, kernel_parameters)

        !print *, "hh"
        if (full_free_surface>0) then
            sub_kernel => envoke_hh_update_kernel_gpu
            sub_sync   => envoke_hh_update_sync_gpu
            call envoke_device(sub_kernel, sub_sync, kernel_parameters)
        endif

        !computing advective and lateral-viscous terms for 2d-velocity
        !print *, "trans"
        if (trans_terms > 0) then
            sub_kernel => envoke_uv_trans_vort_kernel_gpu
            sub_sync   => envoke_uv_trans_vort_sync_gpu
            call envoke_device(sub_kernel, sub_sync, kernel_parameters)

            sub_kernel => envoke_uv_trans_kernel_gpu
            sub_sync   => envoke_uv_trans_sync_gpu
            call envoke_device(sub_kernel, sub_sync, kernel_parameters)
        endif

        !print *, "stress"
        if (ksw_lat > 0) then
            sub_kernel => envoke_stress_components_kernel_gpu
            sub_sync   => envoke_stress_components_sync_gpu
            call envoke_device(sub_kernel, sub_sync, kernel_parameters)

            sub_kernel => envoke_uv_diff2_kernel_gpu
            sub_sync   => envoke_uv_diff2_sync_gpu
            call envoke_device(sub_kernel, sub_sync, kernel_parameters)
        endif

        !print *, "uv"
        sub_kernel => envoke_sw_update_uv_kernel_gpu
        sub_sync   => envoke_sw_update_uv_sync_gpu
        call envoke_device(sub_kernel, sub_sync, kernel_parameters)

        !shifting time indices
        !print *, "ssh upd"
        sub_kernel => envoke_sw_next_step_kernel_gpu
        sub_sync   => envoke_sw_next_step_sync_gpu
        call envoke_device(sub_kernel, sub_sync, kernel_parameters)

        !print *, "hh upd"
        if(full_free_surface>0) then
            sub_kernel => envoke_hh_shift_kernel_gpu
            sub_sync   => envoke_hh_shift_sync_gpu
            call envoke_device(sub_kernel, sub_sync, kernel_parameters)
        endif

        !print *, "hh init"
        if(full_free_surface>0) then
            !initialize depth for external mode
            sub_kernel => envoke_hh_init_kernel_gpu
            sub_sync   => envoke_hh_init_sync_gpu
            call envoke_device(sub_kernel, sub_sync, kernel_parameters)
        endif

        ! Check error
        !sub_kernel => envoke_check_ssh_err_kernel
        !sub_sync   => envoke_check_ssh_err_sync
        !call envoke(sub_kernel, sub_sync, 0.0d0)

    endsubroutine expl_shallow_water_gpu

endmodule

#endif