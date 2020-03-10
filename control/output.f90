module output_module

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use mpp_module, only: mpp_rank, mpp_count, mpp_cart_comm, mpp_size, mpp_coord
    use decomposition_module, only: domain_type
    use ocean_module, only: ocean_type
    use grid_module, only: grid_type
    use io_module, only: write_data

    implicit none
    save
    private

    public :: local_output

contains

    subroutine local_output(domain, grid_data, ocean_data,  &
                            nrec,    &
                            year,    &
                            month,   &
                            day,     &
                            hour,    &
                            minute,  &
                            tstep,   &
                            calendar)
        use config_basinpar_module, only: nx, ny, xgr_type, ygr_type, rlon, rlat, dxst, dyst
        use iodata_routes, only: fulfname, lrecl, lmpirecl, undef
        use rw_ctl_routes, only: ctl_file_write

        type(domain_type), intent(in) :: domain
        type(grid_type), intent(in) :: grid_data
        type(ocean_type), intent(in) :: ocean_data
        integer, intent(in) :: nrec, year, month, day, hour, minute, calendar
        real(wp4) :: tstep

        integer :: ierr
        real(wp8) :: z0(1), z1(1)
        real(wp8) :: xtm1loc(1), ytn1loc(1), xum1loc(1), yvn1loc(1)
        character(4096) :: fname

        xtm1loc = rlon
        ytn1loc = rlat
        !xum1loc = rlon - dxst/2.0
        !yvn1loc = rlat - dyst/2.0
        z0 = 0.0d0
        z1 = 1.0d0

        ierr = 0
        call write_data(domain, 'RESULTS/', 'hhq.dat', nrec, grid_data%hhq_rest, grid_data%lu, ierr)
        if(mpp_rank == 0) print *, 'hhq is written'
        
        call fulfname(fname, 'RESULTS/', 'hhq.dat', ierr)
        if (mpp_rank == 0) then
            call ctl_file_write(fname,     &     !file name
                                undef,     &     !value for undefined points
                                nx - 4,    &     !x-dimension
                                ny - 4,    &     !y-dimension
                                1,         &     !z-dimension
                                nrec,      &     !t-dimension
                                xgr_type,  &     !x-grid type (0 - linear, 1 - levels)
                                xtm1loc,   &     !first x-value (if linear) or x-array (if levels)
                                dxst,      &     !x-step (if linear)
                                ygr_type,  &     !y-grid type (0 - linear, 1 - levels)
                                ytn1loc,   &     !first y-value (if linear) or x-array (if levels)
                                dyst,      &     !y-step (if linear)
                                1,         &     !z-grid type (0 - linear, 1 - levels)
                                z0,        &     !first z-value (if linear) or x-array (if levels)
                                1.0d0,     &     !z-step (if linear)
                                calendar,  &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                                year,      &     !year   of the first field
                                month,     &     !month  of the first field
                                day,       &     !day    of the first field
                                hour,      &     !hour   of the first field
                                minute,    &     !minute of the first field
                                tstep,     &     !time step (in seconds)
                                'HHQ, m',  &     !title of dataset
                                'hhq')           !variable name
        endif
    end subroutine local_output

end module output_module