module config_sw_module
    ! Module description

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4

    implicit none
    save
    public

    integer :: full_free_surface
    real(wp4) :: time_smooth

contains

    subroutine load_config_sw
        
        full_free_surface = 1
        time_smooth = 0.5
        
    end subroutine 

end module config_sw_module