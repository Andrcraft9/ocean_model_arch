module shallow_water_module
    
    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use decomposition_module, only: domain_type
    use ocean_module, only: ocean_type
    use grid_module, only: grid_type
    use config_sw_module, only: full_free_surface, time_smooth, trans_terms, ksw_lat, lvisc_2
    use mpp_module, only: mpp_rank
    use shallow_water_interface_module

    implicit none
    save
    private

    public :: expl_shallow_water

contains

!===============================================================================
    ! explicit shallow water equation sloving
    subroutine expl_shallow_water(tau, domain, grid_data, ocean_data)

        implicit none

        real(wp8), intent(in) :: tau
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data

        ! init
        call ocean_data%mu%fill(domain, lvisc_2)

        !computing ssh
        call envoke_sw_update_ssh_kernel(domain, grid_data, tau, ocean_data%sshn, ocean_data%sshp, ocean_data%ubrtr, ocean_data%vbrtr)

        if (full_free_surface>0) then
            call envoke_hh_update_kernel(domain, grid_data, ocean_data%ssh)
        endif

        !computing advective and lateral-viscous terms for 2d-velocity
        if (trans_terms > 0) then
            call envoke_uv_trans_vort_kernel(domain, grid_data, ocean_data%ubrtr, ocean_data%vbrtr, ocean_data%vort)
            call envoke_uv_trans_kernel(domain, grid_data, ocean_data%ubrtr, ocean_data%vbrtr, ocean_data%vort, ocean_data%RHSx_adv, ocean_data%RHSy_adv)
        endif

        if (ksw_lat > 0) then
            call envoke_stress_components_kernel(domain, grid_data, ocean_data%ubrtrp, ocean_data%vbrtrp, ocean_data%str_t, ocean_data%str_s)
            call envoke_uv_diff2_kernel(domain, grid_data, ocean_data%mu, ocean_data%str_t, ocean_data%str_s, ocean_data%RHSx_dif, ocean_data%RHSy_dif)
        endif

        call envoke_sw_update_uv(domain, grid_data, tau, ocean_data%ssh, &
                                 ocean_data%ubrtr, ocean_data%ubrtrn, ocean_data%ubrtrp, &
                                 ocean_data%vbrtr, ocean_data%vbrtrn, ocean_data%vbrtrp, &
                                 ocean_data%r_diss, &
                                 ocean_data%RHSx, ocean_data%RHSy, ocean_data%RHSx_adv, ocean_data%RHSy_adv, ocean_data%RHSx_dif, ocean_data%RHSy_dif)

        !shifting time indices
        call envoke_sw_next_step(domain, grid_data, time_smooth, &
                                 ocean_data%ssh, ocean_data%sshn, ocean_data%sshp, &
                                 ocean_data%ubrtr, ocean_data%ubrtrn, ocean_data%ubrtrp, &
                                 ocean_data%vbrtr, ocean_data%vbrtrn, ocean_data%vbrtrp)

        if(full_free_surface>0) then
            call envoke_hh_shift_kernel(domain, grid_data)
        endif

        if(full_free_surface>0) then
            !initialize depth for external mode
            call envoke_hh_init_kernel(domain, grid_data, ocean_data%ssh, ocean_data%sshp)
        endif

        ! Check error
        call envoke_check_ssh_err_kernel(domain, grid_data, ocean_data%ssh, 'ssh')

    endsubroutine expl_shallow_water

endmodule
