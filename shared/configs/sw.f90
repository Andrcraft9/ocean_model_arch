module config_sw_module
    ! Module description

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4

    implicit none
    save
    public

    integer :: full_free_surface
    integer :: trans_terms
    integer :: ksw_lat
    real(wp8) :: time_smooth
    
    real(wp8) :: lvisc_2

    character(128) :: ssh_init_file_name

contains

    subroutine load_config_sw
        
        full_free_surface = 1
        time_smooth = 0.5d0
        trans_terms = 1
        ksw_lat = 1

        lvisc_2 = 1.0d+03

        ssh_init_file_name = 'slf.dat'
        
    end subroutine 

end module config_sw_module