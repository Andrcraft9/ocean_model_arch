module time_manager_module

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use mpi
    use mpp_module
    use rwpar_routes
    use time_routes

    implicit none
    save
    private

    ! PUBLIC MEMBERS
    public :: init_time_manager, time_manager_def, time_manager_print, is_local_print_step, time_manager_update_nrec

    integer, parameter, public :: yr_type = 1      !0 - without leap-year, 1 - with leap-year

    integer, public :: year_loc,      &     !variables for writing local data
                       mon_loc,       &
                       day_loc,       &
                       hour_loc,      &
                       min_loc,       &
                       nrec_loc

    integer, public :: year_glob,     &     !variables for writing global data
                       mon_glob,      &
                       day_glob,      &
                       hour_glob,     &
                       min_glob,      &
                       nrec_glob

    real(wp4), public :: loc_data_tstep,       &  !Time step for writing  local data
                         glob_data_tstep          !Time step for writing global data

    integer(wp8), public :: num_step,        &  !Number of time step during the run
                            num_step_max        !The maximum number of time step for the run

    real(wp8), public :: tau

    character(256), public :: filepar,             &          !parameter file name
                              path2ocp,            &          !path to control points(results)
                              path2ocssdata,       &          !path to surface data on ocean grid
                              path2atmssdata,      &          !path to surface data on atm   grid
                              blank                
    character(128), public :: ss_ocfiles(8),   &          !files with sea surface data on oceanic grid
                        ss_ocfiles_fname(8),   &          !files with sea surface data on oceanic grid (fullnames)               
                            ss_atmfiles(14),   &          !files with sea surface data on atmospheric grid
                        ss_atmfiles_fname(14), &          !files with sea surface data on atmospheric grid (fullnames)  
                                atmask                    !file with atmospheric sea-land mask (1-land,0-ocean)
   
    ! PRIVATE MEMBERS
    real(wp4) :: time_step,            &  !Model time step (in seconds)
                 time_step_m,          &  !Model time step (in minutes)
                 time_step_h,          &  !Model time step (in hours)
                 time_step_d,          &  !Model time step (in days)
                 seconds_of_day           !Current seconds of the day

    integer :: start_type,        &  !Type of starting run (0 - from TS only, 1 - from the full checkpoint)
               nstep_icedyn,      &  !Number of internal time steps for ice dynamics
               nstep_barotrop        !Number of internal time steps for barotropic task

    integer :: init_year           !Initial year number for the run

    real(wp4) :: run_duration,         &  !Duration of the run in days
                 loc_data_wr_period,   &  !Period for writing local instantaneous data (in minutes; If >1440, then assumed 1440)
                                          !If <=0 then no local output is done
                 glob_data_wr_period      !Period for writing global time-mean data and checkpoints (in hours; If >1440, then assumed 1 month)
                                          !If <0 then only checkpoints are written

    integer :: loc_data_wr_period_step,       &   !period in steps to write to write local data
               glob_data_wr_period_step,      &   !period in steps to write to write global data
               loc_data_nstep,                &   !Number of record to write local data
               glob_data_nstep,               &   !Number of record to write global data
               monthly_output,                &   !Key for mohthly mean output (0 - no, 1 - yes)
               nofcom                             !number of lines in "ocean_run.par" (calculated)

    integer :: key_write_local,      &      !Key for write local  data(1 - yes, 0 - no)
               key_write_global,     &      !Key for write global data(1 - yes, 0 - no)
           time_write_global

    integer :: ndays_in_4yr(0:48),   &     !integer day distributions in 4-years (selected)
               ndays_noleap(0:48),   &     !integer day distributions in 4-years (without leap-year)
               ndays_leap  (0:48),   &     !integer day distributions in 4-years (with leap-year)
               m_sec_of_min,         &     !second counter in minute
               m_min_of_hour,        &     !minute counter in hour
               m_hour_of_day,        &     !hour counter in day     
               m_day_of_month,       &     !day counter in month
               m_day_of_year,        &     !day counter in year
               m_day_of_4yr,         &     !day counter in 4-years  
               m_month_of_year,      &     !mon counter in year
               m_month_of_4yr,       &     !mon counter in 4-years
               m_year_of_4yr,        &     !year counter in 4yrs
               m_day,                &     !model elapsed day counter starting from zero
               m_month,              &     !model elapsed month counter starting from zero
               m_year,               &     !year counter
               m_4yr,                &     !counter of 4-yr groups 
               m_time_changed(7),    &     !change indicator of time
               key_time_print,       &     !key of printing time:0-not,1-print
               nstep_per_day

    data ndays_noleap/ 0,  31,  59,  90, 120, 151, 181,   &    !1-st half year of 1-st year
                        212, 243, 273, 304, 334, 365,   &    !2-nd half year of 1-st year
                        396, 424, 455, 485, 516, 546,   &    !1-st half year of 2-nd year
                        577, 608, 638, 669, 699, 730,   &    !2-nd half year of 2-nd year
                        761, 789, 820, 850, 881, 911,   &    !1-st half year of 3-rd year
                        942, 973,1003,1034,1064,1095,   &    !2-nd half year of 3-rd year
                        1126,1154,1185,1215,1246,1276,   &    !1-st half year of 4-th year
                        1307,1338,1368,1399,1429,1460/        !2-nd half year of 4-th year

    data ndays_leap  / 0,  31,  59,  90, 120, 151, 181,   &    !1-st half year of 1-st year
                        212, 243, 273, 304, 334, 365,   &    !2-nd half year of 1-st year
                        396, 424, 455, 485, 516, 546,   &    !1-st half year of 2-nd year
                        577, 608, 638, 669, 699, 730,   &    !2-nd half year of 2-nd year
                        761, 789, 820, 850, 881, 911,   &    !1-st half year of 3-rd year
                        942, 973,1003,1034,1064,1095,   &    !2-nd half year of 3-rd year
                        1126,1155,1186,1216,1247,1277,   &    !1-st half year of 4-th year
                        1308,1339,1369,1400,1430,1461/        !2-nd half year of 4-th year

    ! input parameters of task and names of files with input data
    character(256) :: comments(256)

contains

subroutine load_config(name)
    
    character(*), intent(in) :: name

    integer :: ierr

    if (mpp_is_master_thread()) then

        blank = ' '
        filepar = name

        ! reading parameters from file
        if (mpp_is_master()) then
            call readpar(filepar, comments, nofcom)
        endif
        call mpi_bcast(comments, 256*256, mpi_character, 0, mpp_cart_comm, ierr)

        read(comments( 1),*) start_type          !Type of starting run (0 - from TS only, 1 - frim the full checkpoint)
        read(comments( 2),*) time_step           !Model time step (in seconds)
        read(comments( 3),*) run_duration        !Duration of the run in days
        read(comments( 4),*) num_step            !Number of time step during the run
        read(comments( 5),*) init_year           !Initial year number for the run
        read(comments( 6),*) loc_data_wr_period  !Period for writing local instantaneous data (minutes)
        read(comments( 7),*) glob_data_wr_period !Period for writing global time average data (minutes)
        read(comments( 8),*) nstep_icedyn        !Number of internal time steps for ice dynamics
        read(comments( 9),*) nstep_barotrop      !Number of internal time steps for barotropic task
        call get_first_lexeme(comments(10), path2ocp    )  !path to checkpoints(results)
        ! Files with data on oceanic grid:
        call get_first_lexeme(comments(11), path2ocssdata)   !path to ocean SS data
        call get_first_lexeme(comments(12), ss_ocfiles(1)  )  !file with SST 
        call get_first_lexeme(comments(13), ss_ocfiles(2)  )  !file with SSS 
        call get_first_lexeme(comments(14), ss_ocfiles(3)  )  !file with river runoff 
        call get_first_lexeme(comments(15), ss_ocfiles(4)  )  !file with TLBC
        call get_first_lexeme(comments(16), ss_ocfiles(5)  )  !file with SLBC
        call get_first_lexeme(comments(17), ss_ocfiles(6)  )  !file with ULBC
        call get_first_lexeme(comments(18), ss_ocfiles(7)  )  !file with VLBC
        call get_first_lexeme(comments(19), ss_ocfiles(8)  )  !file with SSHLBC
        ! Files with data on atmospheric grid:
        call get_first_lexeme(comments(20), path2atmssdata )  !path to atmospheric data
        call get_first_lexeme(comments(21), ss_atmfiles(1) )  !file with      zonal wind stress (1 and 2 condition)
        call get_first_lexeme(comments(22), ss_atmfiles(2) )  !file with meridional wind stress (1 and 2 condition)
        call get_first_lexeme(comments(23), ss_atmfiles(3) )  !file with          heat balance (2 condition) 
        call get_first_lexeme(comments(24), ss_atmfiles(4) )  !file with shortwave rad balance (2 condition)
        call get_first_lexeme(comments(25), ss_atmfiles(5) )  !file with    freshwater balance (2 condition)
        call get_first_lexeme(comments(26), ss_atmfiles(6) )  !file with       air temperature (3 condition)
        call get_first_lexeme(comments(27), ss_atmfiles(7) )  !file with          air humidity (3 condition)
        call get_first_lexeme(comments(28), ss_atmfiles(8) )  !file with wind      zonal speed (3 condition)
        call get_first_lexeme(comments(29), ss_atmfiles(9) )  !file with wind meridional speed (3 condition)
        call get_first_lexeme(comments(30), ss_atmfiles(10) )  !file with SLP (3 condition)
        call get_first_lexeme(comments(31), ss_atmfiles(11) )  !file with downwelling longwave radiation (3 condition)
        call get_first_lexeme(comments(32), ss_atmfiles(12) )  !file with downwelling shortwave radiation (3 condition)
        call get_first_lexeme(comments(33), ss_atmfiles(13) )  !file with wind liquid precipitation (rain)
        call get_first_lexeme(comments(34), ss_atmfiles(14) )  !file with wind solid precipitation (snow)
        call get_first_lexeme(comments(35), atmask)

    endif
    call mpp_barrier_threads()
    
end subroutine 


subroutine init_time_manager(name)
    
    character(*), intent(in) :: name

    call load_config(name)

    m_sec_of_min   = 0     !second counter in minute
    m_min_of_hour  = 0     !minute counter in hour
    m_hour_of_day  = 0     !hour counter in day     
    m_day_of_month = 0     !day counter in month
    m_day_of_year  = 0     !day counter in year
    m_day_of_4yr   = 0     !day counter in 4-years  
    m_month_of_year= 0     !mon counter in year
    m_month_of_4yr = 0     !mon counter in 4-years
    m_year_of_4yr  = 0     !year counter in 4yrs
    m_day          = 0     !model elapsed day counter starting from zero
    m_month        = 0     !model elapsed month counter starting from zero
    m_year         = 0     !year counter
    m_4yr          = 0     !counter of 4-yr groups 
    m_time_changed = 0     !change indicator of time
    key_time_print = 0     !key of printing time:0-not,1-print

    seconds_of_day = 0.0   !Current seconds of the day

    if(yr_type==0) then  !day distribution without leap-year
        ndays_in_4yr=ndays_noleap
    elseif(yr_type==1) then   !day distribution with leap-year
        ndays_in_4yr=ndays_leap
    endif

    if(loc_data_wr_period>0.0001) then
        key_write_local=1
    else
        key_write_local=0
    endif

    if(glob_data_wr_period>0.0001) then
        key_write_global=1
    else
        key_write_global=0
    endif

    time_write_global=0

    !computing some time parameters
    time_step_m = time_step/60.0
    time_step_h = time_step/3600.0
    time_step_d = time_step/86400.0
    nstep_per_day=nint(86400.0/time_step)       !NUMBER OF STEP PER DAY

    !for local output
    loc_data_wr_period = max(min(loc_data_wr_period,1440.0),time_step_m) !The maximum local output period is 1 day, minimum is time step
    loc_data_wr_period_step = nint(loc_data_wr_period/time_step_m)   !period in steps to write to write local data
    loc_data_tstep = loc_data_wr_period * 60.0                      !Time step in seconds for writing local data

    year_loc=init_year
    mon_loc=1
    day_loc=int(loc_data_wr_period/1440.0)+1
    hour_loc=mod(int(loc_data_wr_period/60.0),24)
    min_loc=mod(int(loc_data_wr_period),60)
    
    !for global output  
    if(abs(glob_data_wr_period)>1440.5) then
        monthly_output=1
        glob_data_tstep = 86400.0*30.0

        year_glob=init_year
        mon_glob=1
        day_glob=15
        hour_glob=0
        min_glob=0
    else
        monthly_output=0
        glob_data_wr_period = max(min(abs(glob_data_wr_period),1440.0),time_step_m)
        glob_data_wr_period_step = nint(glob_data_wr_period/time_step_m)   !period in steps to write to write global data
        glob_data_tstep = glob_data_wr_period * 60.0                        !Time step in seconds for writing global data

        year_glob=init_year
        mon_glob=1
        day_glob=1
        hour_glob=mod(int(glob_data_wr_period/2.0/60.0),24)
        min_glob=mod(int(glob_data_wr_period/2.0),60)
    endif
    
    num_step_max=int8(run_duration*nstep_per_day)

    key_time_print = 0

    tau = time_step
end subroutine

subroutine time_manager_def()

    call model_time_def(num_step,            &    !step counter,            input
                        time_step,           &    !time step in seconds,    input
                        ndays_in_4yr,        &    !integer day distribution in 4-years (49 months)
                        seconds_of_day,      &    !current seconds in day  ,output
                        m_sec_of_min,        &    !second counter in minute,output
                        m_min_of_hour,       &    !minute counter in hour  ,output
                        m_hour_of_day,       &    !hour counter in day     ,output
                        m_day_of_month,      &    !day counter in month    ,output
                        m_day_of_year,       &    !day counter in year     ,output
                        m_day_of_4yr,        &    !day counter in 4-years  ,output
                        m_month_of_year,     &    !mon counter in year     ,output
                        m_month_of_4yr,      &    !mon counter in 4-years  ,output
                        m_year_of_4yr,       &    !year counter in 4yrs    ,output
                        m_day,               &    !model elapsed day counter starting from zero
                        m_month,             &    !model elapsed month counter starting from zero
                        m_year,              &    !year counter            ,output
                        m_4yr,               &    !counter of 4-yr groups  ,output
                        m_time_changed,      &    !change indicator of time,output
                        key_time_print,      &    !key of printing time:0-not,1-print
                        init_year)                !initial real-time year
    
end subroutine

subroutine time_manager_print()
    
    call model_time_print(num_step,         &
                          m_sec_of_min,     &    !second counter in minute,output
                          m_min_of_hour,    &    !minute counter in hour  ,output
                          m_hour_of_day,    &    !hour counter in day     ,output
                          m_day_of_month,   &    !day counter in month    ,output
                          m_day_of_year,    &    !day counter in year     ,output
                          m_day_of_4yr,     &    !day counter in 4-years  ,output
                          m_month_of_year,  &    !mon counter in year     ,output
                          m_month,          &    !model elapsed month counter starting from zero
                          m_year )               !year counter            ,output
    
end subroutine

subroutine time_manager_update_nrec ()
    
    nrec_loc = num_step/loc_data_wr_period_step
end subroutine

function is_local_print_step() result(flag)
    integer :: flag

    flag = 0
    if( key_write_local>0) then
        if(mod(num_step, loc_data_wr_period_step)==0) then
            flag = 1
        endif
    endif

end function

end module time_manager_module