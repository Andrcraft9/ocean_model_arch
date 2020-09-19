module shallow_water_module
    
    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use decomposition_module, only: domain_type
    use ocean_module, only: ocean_type
    use grid_module, only: grid_type
    use config_sw_module, only: full_free_surface, time_smooth, trans_terms, ksw_lat
    use mpp_module
    use kernel_interface_module
    use shallow_water_interface_module

    implicit none
    save
    private

    public :: expl_shallow_water

contains

!===============================================================================
    ! explicit shallow water equation sloving
    subroutine expl_shallow_water(tau, domain, grid_data, ocean_data)

        real(wp8), intent(in) :: tau
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data

        procedure(envoke_empty_kernel), pointer :: sub_kernel
        procedure(envoke_empty_sync), pointer :: sub_sync

        !computing ssh
        sub_kernel => envoke_sw_update_ssh_kernel
        sub_sync   => envoke_sw_update_ssh_sync
        call envoke(domain, grid_data, ocean_data, sub_kernel, sub_sync, tau)

        if (full_free_surface>0) then
            sub_kernel => envoke_hh_update_kernel
            sub_sync   => envoke_hh_update_sync
            call envoke(domain, grid_data, ocean_data, sub_kernel, sub_sync, 0.0d0)
        endif

        !computing advective and lateral-viscous terms for 2d-velocity
        if (trans_terms > 0) then
            sub_kernel => envoke_uv_trans_vort_kernel
            sub_sync   => envoke_uv_trans_vort_sync
            call envoke(domain, grid_data, ocean_data, sub_kernel, sub_sync, 0.0d0)

            sub_kernel => envoke_uv_trans_kernel
            sub_sync   => envoke_uv_trans_sync
            call envoke(domain, grid_data, ocean_data, sub_kernel, sub_sync, 0.0d0)
        endif

        if (ksw_lat > 0) then
            sub_kernel => envoke_stress_components_kernel
            sub_sync   => envoke_stress_components_sync
            call envoke(domain, grid_data, ocean_data, sub_kernel, sub_sync, 0.0d0)

            sub_kernel => envoke_uv_diff2_kernel
            sub_sync   => envoke_uv_diff2_sync
            call envoke(domain, grid_data, ocean_data, sub_kernel, sub_sync, 0.0d0)
        endif

        sub_kernel => envoke_sw_update_uv_kernel
        sub_sync   => envoke_sw_update_uv_sync
        call envoke(domain, grid_data, ocean_data, sub_kernel, sub_sync, tau)

        !shifting time indices
        sub_kernel => envoke_sw_next_step_kernel
        sub_sync   => envoke_sw_next_step_sync
        call envoke(domain, grid_data, ocean_data, sub_kernel, sub_sync, time_smooth)

        if(full_free_surface>0) then
            sub_kernel => envoke_hh_shift_kernel
            sub_sync   => envoke_hh_shift_sync
            call envoke(domain, grid_data, ocean_data, sub_kernel, sub_sync, 0.0d0)
        endif

        if(full_free_surface>0) then
            !initialize depth for external mode
            sub_kernel => envoke_hh_init_kernel
            sub_sync   => envoke_hh_init_sync
            call envoke(domain, grid_data, ocean_data, sub_kernel, sub_sync, 0.0d0)
        endif

        ! Check error
        sub_kernel => envoke_check_ssh_err_kernel
        sub_sync   => envoke_check_ssh_err_sync
        call envoke(domain, grid_data, ocean_data, sub_kernel, sub_sync, 0.0d0)

    endsubroutine expl_shallow_water

endmodule
