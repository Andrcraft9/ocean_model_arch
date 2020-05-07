module time_routes

use mpp_module

implicit none

contains
!======================================================================
subroutine model_time_def(   num_step,            &     !step counter,            input
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

!+++++++++++++++++time control parameters+++++++++++++++++++++++++++++++
real(4)  time_step,       &   !time step in seconds
         seconds_of_day       !current seconds in day

integer(8) num_step           !time step counter

integer ndays_in_4yr(0:48),     &   !integer day distributions in 4-years
        m_sec_of_min,           &   !second counter in minute
        m_min_of_hour,          &   !minute counter in hour
        m_hour_of_day,          &   !hour counter in day
        m_day_of_month,         &   !day  counter in month
        m_day_of_year,          &   !day  counter in year
        m_day_of_4yr,           &   !day  counter in 4-years
        m_month_of_year,        &   !mon  counter in year
        m_month_of_4yr,         &   !mon counter in 4-years
        m_year_of_4yr,          &   !year counter in 4yrs
        m_day,                  &   !model elapsed day counter starting from zero
        m_month,                &   !model elapsed month counter starting from zero
        m_year,                 &   !year counter 
        m_4yr,                  &   !counter of 4-yr groups
        m_time_changed(7),      &   !indicator of time changed (0-not,1-changed) for
                                    !1-sec,2-min,3-hour,4-day,5-month,6-year,7-4yrs
        key_time_print,         &   !key of printing time:0-not,1-print     
        init_year                    !initial real-time year

!------------------ internal variables: --------------------------------
! help variable preffix i denotes instant
integer i_day_of_year,    & !day counter in year
        i_day_of_4yr,     & !day counter in 4-year
        i_4yr,            & !counter of 4-yr groups
        i_sec_of_min,     & !second counter in minute
        i_min_of_hour,    & !minute counter in hour
        i_hour_of_day,    & !hour counter in day
        i_day_of_month,   & !day counter in month
        i_month_of_year,  & !mon counter in year
        i_year,           & !year counter 
        i_month_of_4yr,   & !mon counter in 4-years
        i_year_of_4yr       !year counter in 4yrs

character  month_name(12)*3
data month_name/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug',     &
                'Sep','Oct','Nov','Dec'/
integer i,k
integer nshift        ! number of years (0-3) for counter shift  
                      ! if the real-time initial year does not follow
                      ! the leap year
integer(8) ihelp8, khelp8
!-----------------------------------------------------------------------
  nshift=mod(init_year-1,4)
  
  m_time_changed=0
  
  ihelp8=int8(86400.0/time_step)

  m_day=int(num_step/ihelp8)+1
   
  khelp8=mod(num_step,ihelp8)
  seconds_of_day = time_step*dfloat(khelp8)
  
  i_4yr=(m_day-1+ndays_in_4yr(nshift*12))/ndays_in_4yr(48)+1
 
  i_day_of_4yr=mod(m_day-1+ndays_in_4yr(nshift*12), ndays_in_4yr(48))+1

          i_month_of_4yr=1
          i=49
1000      k=(i_month_of_4yr+i)/2
          if (i_day_of_4yr.le.ndays_in_4yr(k-1)) i=k
          if (i_day_of_4yr.gt.ndays_in_4yr(k-1)) i_month_of_4yr=k
          if (i.gt.i_month_of_4yr+1) go to 1000


      i_year_of_4yr = (i_month_of_4yr-1)/12 +1
      
      i_year        = (i_4yr-1)*4 + i_year_of_4yr +init_year-1-nshift

      i_month_of_year = mod(i_month_of_4yr-1,12)+1
      
      i_day_of_month= i_day_of_4yr-ndays_in_4yr(i_month_of_4yr-1)
      
      i_day_of_year = i_day_of_4yr-ndays_in_4yr((i_year_of_4yr-1)*12)
                
      i_hour_of_day = int(seconds_of_day/3600.0)

      i_min_of_hour = int(seconds_of_day/60.0)   ! min in day yet

      i_sec_of_min  = int(seconds_of_day) -  i_min_of_hour * 60

      i_min_of_hour = i_min_of_hour      -  i_hour_of_day * 60
                                        
      if(m_sec_of_min  .ne.i_sec_of_min  ) then
         m_sec_of_min    = i_sec_of_min
         m_time_changed(1)=1  
      end if  

      if(m_min_of_hour .ne.i_min_of_hour ) then
         m_min_of_hour   = i_min_of_hour
         m_time_changed(2)=1  
      end if 

      if(m_hour_of_day .ne.i_hour_of_day ) then
         m_hour_of_day   = i_hour_of_day
         m_time_changed(3)=1
      end if  
       
      if(m_day_of_month.ne.i_day_of_month) then
         m_day_of_4yr    = i_day_of_4yr
         m_day_of_year   = i_day_of_year
         m_day_of_month  = i_day_of_month
         m_time_changed(4)=1  
      end if
 
      if(m_month_of_year .ne.i_month_of_year ) then
         m_month_of_year   = i_month_of_year
       m_month_of_4yr    = i_month_of_4yr
         m_time_changed(5)=1  
      end if 
 
      if(m_year .ne.i_year ) then
         m_year        = i_year
         m_year_of_4yr = i_year_of_4yr
         m_time_changed(6)=1  
      end if 
      
      if(m_4yr .ne.i_4yr ) then
         m_4yr   = i_4yr
         m_time_changed(7)=1  
      end if 
         
         m_month=(m_year-init_year)*12+m_month_of_year
      
      if(key_time_print.ne.0) then
        if (mpp_is_master()) then
          write(*,'(a,i8,4(a,i2.2), a,i4.4, a,i3.3, a,i4.4,a,i5)')  &
             '   time step: ', num_step ,                      &
             '   model time: ',                           &
                            m_hour_of_day,':',            &
                            m_min_of_hour,':',            &
                            m_sec_of_min, '  ',           &
                            m_day_of_month,               &
                 month_name(m_month_of_year),m_year,        &
                 ';  day in year:',m_day_of_year,         &
                 ',  day in 4yrs:',m_day_of_4yr,          &
                 '   month of all:', m_month      
        endif
      endif 
end subroutine model_time_def
!======================================================================
subroutine model_time_print(num_step,         &
                            m_sec_of_min,     &    !second counter in minute,output
                            m_min_of_hour,    &    !minute counter in hour  ,output
                            m_hour_of_day,    &    !hour counter in day     ,output
                            m_day_of_month,   &    !day counter in month    ,output
                            m_day_of_year,    &    !day counter in year     ,output
                            m_day_of_4yr,     &    !day counter in 4-years  ,output
                            m_month_of_year,  &    !mon counter in year     ,output
                            m_month,          &    !model elapsed month counter starting from zero
                            m_year )               !year counter            ,output

integer(8) num_step
integer m_hour_of_day,      &
        m_min_of_hour,      &
        m_sec_of_min,       &
        m_day_of_month,     &
        m_month_of_year,    &
        m_month,            &
        m_year,             &
        m_day_of_year,      &
        m_day_of_4yr

character  month_name(12)*3
data month_name/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug',     &
                'Sep','Oct','Nov','Dec'/  
      
      if (mpp_is_master()) then   
          write(*,'(a,i8,4(a,i2.2), a,i4.4, a,i3.3, a,i4.4,a,i5)')  &
             '   time step: ', num_step ,                      &
             '   model time: ',                           &
                            m_hour_of_day,':',            &
                            m_min_of_hour,':',            &
                            m_sec_of_min, '  ',           &
                            m_day_of_month,               &
                 month_name(m_month_of_year),m_year,        &
                 ';  day in year:',m_day_of_year,         &
                 ',  day in 4yrs:',m_day_of_4yr,          &
                 '   month of all:', m_month 
      endif           
endsubroutine model_time_print

endmodule time_routes