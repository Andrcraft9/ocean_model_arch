module config_sw_module
    ! Module description

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use mpp_module
    use rwpar_routes

    implicit none
    save
    public

    integer :: full_free_surface
    integer :: trans_terms
    integer :: ksw_lat
    real(wp8) :: time_smooth
    real(wp8) :: lvisc_2
    integer :: use_tracers
    integer :: tracer_num
    character(128) :: ssh_init_file_name

contains

    subroutine load_config_sw_from_file(filepar)
        character(*), intent(in) :: filepar
        character(256) :: comments(256)
        integer :: nofcom, ierr

        ! reading parameters from file
        if (mpp_is_master()) then
            call readpar(filepar, comments, nofcom)
        endif
        call mpi_bcast(comments, 256*256, mpi_character, 0, mpp_cart_comm, ierr)

        read(comments( 1),*) full_free_surface
        read(comments( 2),*) trans_terms
        read(comments( 3),*) ksw_lat
        read(comments( 4),*) time_smooth
        read(comments( 5),*) lvisc_2
        read(comments( 6),*) use_tracers
        read(comments( 7),*) tracer_num
        call get_first_lexeme(comments(8), ssh_init_file_name)

        if (mpp_is_master()) then
            if (use_tracers > 0) then
                print *, "CONFIG: SW: Use tracers. Amount of tracers is ", tracer_num
            endif
        endif
    
        call mpp_sync_output()
    end subroutine

    subroutine load_config_sw_initial_ssh
        
        full_free_surface = 1
        time_smooth = 0.5d0
        trans_terms = 1
        ksw_lat = 1
        lvisc_2 = 1.0d+03
        use_tracers = 0
        tracer_num = 0
        ssh_init_file_name = 'slf.dat'
        
    end subroutine

    subroutine load_config_sw_test
        
        full_free_surface = 1
        time_smooth = 0.5d0
        trans_terms = 1
        ksw_lat = 1
        lvisc_2 = 1.0d+03
        use_tracers = 1
        tracer_num = 1
        ssh_init_file_name = 'none'
        
    end subroutine

end module config_sw_module