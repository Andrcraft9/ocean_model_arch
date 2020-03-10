module rw_ctl_routes

  implicit none
  save
  public

contains

subroutine ctl_file_write(fname,    &     !file name
                          undef,    &     !value for undefined points
                          nx,       &     !x-dimension
                          ny,       &     !y-dimension
                          nz,       &     !z-dimension
                          nt,       &     !t-dimension
                          xtype,    &     !x-grid type (0 - linear, 1 - levels)
                          x0,       &     !first x-value (if linear) or x-array (if levels)
                          hx,       &     !x-step (if linear)
                          ytype,    &     !y-grid type (0 - linear, 1 - levels)
                          y0,       &     !first y-value (if linear) or x-array (if levels)
                          hy,       &     !y-step (if linear)
                          ztype,    &     !z-grid type (0 - linear, 1 - levels)
                          z0,       &     !first z-value (if linear) or x-array (if levels)
                          hz,       &     !z-step (if linear)
                          yr_type,  &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year0,    &     !year   of the first field
                          month0,   &     !month  of the first field
                          day0,     &     !day    of the first field
                          hour0,    &     !hour   of the first field
                          minute0,  &     !minute of the first field
                          ht,       &     !time step (in seconds)
                          title,    &     !title of dataset
                          varname   )     !variable name

character*(*) fname
character*(*) title
integer       nx, ny, nz, nt     !dimension of data

real(8)  x0(*),hx,      &   !x-coordinate parameters
         y0(*),hy,      &   !y-coordinate parameters 
         z0(*),hz           !z-coordinate parameters 

real     undef,         &   !undefinite value
         ht                 !time step

integer xtype,ytype,ztype,    &  !grid types
        yr_type,              &  !calendar type
        year0, month0, day0,  &  !Time parameters
        hour0, minute0

integer i
character*(*) varname
integer nyrht,nmoht,ndyht,nhrht,nmnht  !years,months,days,hours,minutes
!integer nyrt0,nmot0,ndyt0,nhrt0,nmnt0  !years,months,days,hours,minutes
character(128)  namectl,namedat,namedat2
character  tfstep*2,mon*3,month(12)*3
data month/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG',       &
           'SEP','OCT','NOV','DEC'/
integer k

!-------- create extensions ------------
namedat = fname
namectl = fname
call extname(namectl,'.ctl')
call extname(namedat,'.dat')

!--------  remove path from NAMEDAT ( .../.../xxxx.dat -> xxxx.dat )  ----
k  = len(namedat)
do while(namedat(k:k).NE.'/')
    k = k - 1
    if (k.le.0) exit
enddo
write(namedat2, '(a)') namedat( k+1: )

!-------- write to file     ------------
open (40,file=namectl,status='unknown')
write(40,'(a,a)') 'DSET    ^', namedat2
write(40,'(a,a)') 'TITLE    ', title

write(40,'(a,e12.5,a)'  ) 'UNDEF   ',undef,'  ! gap value'
 
if(xtype.eq.0) then
 write(40,'(a,i6,a,2(g15.8,5x))')         &
         'XDEF  ',nx,'  LINEAR   ', x0(1),hx
else
write(40,'(a,i6,a,7(g15.8,1x)/(22x,7(g15.8,1x)))')  &
         'XDEF  ',nx,'  LEVELS  ',(x0(i),i=1,nx)     
end if

if(ytype.eq.0) then
 write(40,'(a,i6,a,2(g15.8,5x))')         &
         'YDEF  ',ny,'  LINEAR   ', y0(1),hy
else
 write(40,'(a,i6,a,7(g15.8,1x)/(22x,7(g15.8,1x)))') &
         'YDEF  ',ny,'  LEVELS  ',(y0(i),i=1,ny)     
end if

if(ztype.eq.0) then
 write(40,'(a,i6,a,2(g15.8,5x))')         &
         'ZDEF  ',nz,'  LINEAR   ', z0(1),hz
else
 write(40,'(a,i6,a,7(g15.8,1x)/(22x,7(g15.8,1x)))') &
         'ZDEF  ',nz,'  LEVELS  ',(z0(i),i=1,nz)     
end if


call sec_to_yr_mo_hr_mn(ht,nyrht,nmoht,ndyht,nhrht,nmnht)

if(nyrht.ne.0) then
 tfstep='yr'
 k=nyrht
elseif(nmoht.ne.0) then
 tfstep='mo'
 k=nmoht
elseif(ndyht.ne.0) then
 tfstep='dy'
 k=ndyht
elseif(nhrht.ne.0) then
 tfstep='hr'
 k=nhrht
elseif(nmnht.ne.0) then
 tfstep='mn'
 k=nmnht
else
write(*,*)' Warning! time step error in creating ctl file ',namectl
end if


write(40,'(a,i6,a,i2.2,a,i2.2,a,i2.2,a,i4.4,5x,i4,a)')          &
          'TDEF  ',nt,'  LINEAR  ',hour0,':',minute0,'Z',       &
           day0,month(month0),year0,k,tfstep

if(yr_type.eq.0) then
write(40,'(a)') 'OPTIONS  365_DAY_CALENDAR'
endif


write(40,'(a)')  'VARS 1  ! Number of variables'
write(40,'(a,a,i5,a)')                    &
           varname, '        ', nz,  '  1 VARIABLES '
write(40,'(a)') 'ENDVARS'
write(40,*)
close (40,status='keep')

return
endsubroutine ctl_file_write

!====================================================================
subroutine extname(name,ext)
!--------------------------------------------------------------------
!     cONSTRUCTION OF EXTENSION  FOR THE  name
!      input : name
!              ext  (MUST CONSIST OF FOR SYMBOLS '.DAT' FOR EXAMPLE)
!         output: name = name + ext
!--------------------------------------------------------------------
   
      character*(*)  name
      character*(*)  ext
      integer        n
!     integer        nch, nchb
      do  n=len(name)-3,2,-1
      if (name(n:n).eq.'.') go to 5
      end do
      n=n+1
  5   name(n:n+4) = ext

endsubroutine extname

!======================================================================
subroutine sec_to_yr_mo_hr_mn(sec,nyr,nmo,ndy,nhr,nmn)

      integer nyr,nmo,ndy,nhr,nmn  !numbers of years,months,days,hours,minutes
      real sec, yr_sec,mo_sec,sec1 !time in seconds
     
      parameter (yr_sec=360.0*86400.0,mo_sec=30.0*86400.0)

      sec1=sec+0.10
      nyr  = int(sec1/yr_sec)
      sec1=(sec1-float(nyr)*yr_sec)

      nmo  = int(sec1/mo_sec)
      sec1=(sec1-float(nmo)*mo_sec)

      ndy  = int(sec1/86400.0)
      sec1=(sec1-float(ndy)*86400.0)

      nhr  = int(sec1/3600.0)
      sec1=(sec1-float(nhr)*3600.0)

      nmn  = int(sec1/60.0)

endsubroutine sec_to_yr_mo_hr_mn
!======================================================================
 subroutine ctl_file_read(ctlname,    &     !file name
                          datname,    &     !file name
                          undef,    &     !value for undefined points
                          nx,       &     !x-dimension
                          ny,       &     !y-dimension
                          nz,       &     !z-dimension
                          nt,       &     !t-dimension
                          xtype,    &     !x-grid type (0 - linear, 1 - levels)
                          x0,       &     !first x-value (if linear) or x-array (if levels)
                          hx,       &     !x-step (if linear)
                          ytype,    &     !y-grid type (0 - linear, 1 - levels)
                          y0,       &     !first y-value (if linear) or x-array (if levels)
                          hy,       &     !y-step (if linear)
                          ztype,    &     !z-grid type (0 - linear, 1 - levels)
                          z0,       &     !first z-value (if linear) or x-array (if levels)
                          hz,       &     !z-step (if linear)
                          nvar,     &     !number of variables
                          yr_type,  &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year0,    &     !year   of the first field
                          month0,   &     !month  of the first field
                          day0,     &     !day    of the first field
                          hour0,    &     !hour   of the first field
                          minute0,  &     !minute of the first field
                          ht,       &     !time step (in seconds)
                          title,    &     !title of dataset
                          varname,  &     !variable name
                          inform    )     

      character*(*) ctlname,datname
      character*(*) title
      character*(*) varname(1)
      integer   nvar              !number of variables
    character(len=32) vartype, gridtype, timeinit, timestep,option

      character(256) string
      integer i,l,l1,l2,nb,line,inform,lprint,k
!  external functions
!    integer nonblank
      integer       nx, ny, nz, nt     !dimension of data
character month1(12)*3,month2(12)*3,month3(12)*3
data month1/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug',       &
           'Sep','Oct','Nov','Dec'/
data month2/'jan','feb','mar','apr','may','jun','jul','aug',       &
           'sep','oct','nov','dec'/                 
data month3/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG',       &
           'SEP','OCT','NOV','DEC'/

      real(8)  x0(*),hx,      &   !x-coordinate parameters
               y0(*),hy,      &   !y-coordinate parameters 
               z0(*),hz           !z-coordinate parameters 
      
      real     undef,         &   !undefinite value
               ht                 !time step
      
      integer xtype,ytype,ztype,    &  !grid types
              yr_type,              &  !calendar type
              year0, month0, day0,  &  !Time parameters
              hour0, minute0

    if(inform.eq.0) then
          lprint=0
    else
        lprint=1
    end if

        inform=0

      open (40,file=ctlname,status='old',action='read')
      
    datname=''
      
      yr_type=1        
    do line=1,512
       
     read(40,'(a)',err=21) string
     nb=nonblank(string,1)
     if(nb.eq.0) go to 7
 
 if(string(nb:nb+6).eq.'ENDVARS'.or. string(nb:nb+6).eq.'endvars'.or.string(nb:nb+6).eq.'Endvars') then
        close (40)

! Printing information about CTL file if need
   if(lprint.ne.0) then

       write(*,'(a,a)') 'DSET    ^', datname
       write(*,'(a,a)') 'TITLE    ', title

       write(*,'(a,e12.5,a)'  ) 'UNDEF   ',undef,'  ! GAP VALUE'
      
     if(xtype==0) then
        write(*,'(a,i6,a,2(g15.8,5x))') 'XDEF  ',nx,'  LINEAR   ', x0(1),hx
     else
        write(*,'(a,i6,a,7(g15.8,1x)/(22x,7(g15.8,1x)))') 'XDEF  ',nx,'  LEVELS  ',(x0(i),i=1,nx)     
     end if

     if(ytype==0) then
        write(*,'(a,i6,a,2(g15.8,5x))')'YDEF  ',ny,'  LINEAR   ', y0(1),hy
     else
        write(*,'(a,i6,a,7(g14.7,1x)/(22x,7(g14.7,1x)))') 'YDEF  ',ny,'  LEVELS  ',(y0(i),i=1,ny)     
     end if
         
      write(*,*) 'TDEF  ',nt,'  LINEAR  ',timeinit,timestep

     if(ztype==0) then
        write(*,'(a,i6,a,2(g15.8,5x))') 'ZDEF  ',nz,'  LINEAR   ', z0(1),hz
       else
        write(*,'(a,i6,a,7(g15.8,1x)/(22x,7(g15.8,1x)))') 'ZDEF  ',nz,'  LEVELS  ',(z0(i),i=1,nz)
       endif
       
       if(yr_type.eq.0) then
        write(*,'(a)') 'OPTIONS  365_DAY_CALENDAR'
       endif

       write(*,'(a,i3,a)')  'VARS ', nvar ,  '! NUMBER OF VARIABLES'
       write(*,'(a,a,i5,a)') varname, '        ', nz,  '  1 VARIABLES '
       write(*,'(a)') 'ENDVARS'
       write(*,*)

    end if
      
    return
   end if
! define number of variables
     if(string(nb:nb+3).eq.'VARS'.or. string(nb:nb+3).eq.'Vars'.or. string(nb:nb+3).eq.'vars' ) then
         
        read(string,*,err=22) vartype,nvar 
        
        do l=1,nvar
        read(40,*,err=22) varname(l),i
        end do
        inform=inform+1

     end if 

! define name of file with data
     if(string(nb:nb+3).eq.'DSET'.or. string(nb:nb+3).eq.'dset'.or. string(nb:nb+3).eq.'Dset' ) then
          
         do  l1=nb+4,len(string)
          if (string(l1:l1).ne.' '.and.string(l1:l1).ne.'^') go to 5
         end do
5       continue
       do  l2=min(len(string),len(datname)+l1-1),nb+4,-1
          if (string(l2:l2).ne.' ') go to 6
         end do 
        l2=l2-1
6       continue
        
      if(l1.gt.l2) then               
         write(*,'(a,a)') ' Error in ',ctlname
         write(*,*) ' Error in definition name of file with data!'
         return
      else
           datname(1:l2-l1+1)=string(l1:l2)
         inform=inform+1
      end if

     end if

! define title of data
     if(string(nb:nb+4).eq.'TITLE'.or. string(nb:nb+4).eq.'Title'.or. string(nb:nb+4).eq.'title' ) then

        VARTYPE=string(nb:nb+4)
           read(string(NB+5:256),'(A)',err=14) title
          inform=inform+1

     end if

! define undefinit value
     if(string(nb:nb+4).eq.'UNDEF'.or. string(nb:nb+4).eq.'undef'.or. string(nb:nb+4).eq.'Undef' ) then
        read(string,*,err=8) vartype,undef 
        inform=inform+1
     end if   

       if(string(nb:nb+6).eq.'OPTIONS'.or. string(nb:nb+6).eq.'options'.or. string(nb:nb+6).eq.'Options' ) then
     
          read(string,*,err=9) vartype, option 
        inform=inform+1
          if(option(1:3)=='365') yr_type=0

     end if 
     
! define x grid
     if(string(nb:nb+3).eq.'xdef'.or. string(nb:nb+3).eq.'Xdef'.or.string(nb:nb+3).eq.'XDEF'  ) then

        read(string,*,err=10) vartype,nx,gridtype

          if(gridtype.eq.'linear'.or. gridtype.eq.'LINEAR'.or. gridtype.eq.'Linear'    ) then
        xtype=0
        read(string,*,err=10) vartype,nx,gridtype,x0(1),hx         
        
           do i=2,nx
           x0(i)=x0(1)+hx*DFLOAT(I-1)
           end do
           inform=inform+1

        else if(gridtype.eq.'levels'.or. gridtype.eq.'LEVELS'.or. gridtype.eq.'Levels'    ) then
                  xtype=1
              backspace (40)
                read(40,*,err=10) vartype,l,gridtype,(x0(i),i=1,nx)
                hx=0.0
                inform=inform+1
        end if

     end if
        

! define y grid
     if(string(nb:nb+3).eq.'ydef'.or. string(nb:nb+3).eq.'Ydef'.or. string(nb:nb+3).eq.'YDEF'    ) then

        read(string,*,err=11) vartype,ny,gridtype

          if(gridtype.eq.'linear'.or. gridtype.eq.'LINEAR'.or. gridtype.eq.'Linear'    ) then
        ytype=0
        read(string,*,err=11) vartype,ny,gridtype,y0(1),hy  
               
           do i=2,ny
           y0(i)=y0(1)+hy*DFLOAT(I-1)
           end do
           inform=inform+1

        else if(gridtype.eq.'levels'.or.  gridtype.eq.'LEVELS'.or.  gridtype.eq.'Levels'  ) then
            ytype=1
              backspace (40)
                read(40,*,err=11) vartype,l,gridtype,(y0(i),i=1,ny)
                  hy=0.0
                inform=inform+1
        end if

     end if

! define z grid
     if(string(nb:nb+3).eq.'zdef'.or. string(nb:nb+3).eq.'Zdef'.or.  string(nb:nb+3).eq.'ZDEF'  ) then

        read(string,*,err=12) vartype,nz,gridtype

          if(gridtype.eq.'linear'.or. gridtype.eq.'LINEAR'.or. gridtype.eq.'Linear'    ) then
            ztype=0
         read(string,*,err=12) vartype,nz,gridtype,z0(1),hz        

           do i=2,nz
           z0(i)=z0(1)+hz*DFLOAT(I)
           end do
           inform=inform+1

        else if(gridtype.eq.'levels'.or.  gridtype.eq.'LEVELS'.or. gridtype.eq.'Levels'    ) then
               ztype=1
              backspace (40)
                read(40,*,err=12) vartype,l,gridtype,(z0(i),i=1,nz)
                hz=0.0
                inform=inform+1
        end if

     end if
                
! define t grid
     if(string(nb:nb+3).eq.'tdef'.or. string(nb:nb+3).eq.'TDEF'.or. string(nb:nb+3).eq.'Tdef'    ) then
        read(string,*,err=13) vartype,nt,gridtype,timeinit,timestep

          if(gridtype.eq.'linear'.or. gridtype.eq.'LINEAR'.or. gridtype.eq.'Linear'    ) then
        
        do  l2=len(timestep),2,-1
          if (timestep(l2:l2).ne.' ') go to 26
          end do 

                              
26        continue            

           read(timeinit(1:2),*) hour0
           read(timeinit(4:5),*) minute0
           read(timeinit(7:8),*) day0
           read(timeinit(12:15),*) year0

           do k=1,12
            if(timeinit(9:11)==month1(k).or.timeinit(9:11)==month2(k).or.timeinit(9:11)==month3(k)) then
             month0=k
            endif
           enddo

        read(timestep(1:l2-2),*,err=13) ht

          if    (timestep(l2-1:l2).eq.'yr') then
                 ht=ht*365.0*24.0*3600.0

          elseif(timestep(l2-1:l2).eq.'mo') then
                 ht=ht*30.0*24.0*3600.0

          elseif(timestep(l2-1:l2).eq.'dy') then
                 ht=ht     *24.0*3600.0

          elseif(timestep(l2-1:l2).eq.'hr') then
                 ht=ht          *3600.0

          elseif(timestep(l2-1:l2).eq.'mn') then
                 ht=ht          *60.0


          else
           write(*,'(a,a)') ' error in ',ctlname
           write(*,*) ' in definition time step for t-grid!'
             close (40)
             return
          end if

          
!         read(string(nb+4:l1-1),*,err=11) nt
!        read(string(l1+7:len(string)),*,err=11) t0(1),ht
!             do i=2,nt
!             t0(i)=t0(i-1)+ht
!             end do
!        inform=inform+1

        else

        write(*,'(a,a)') ' warning in ',ctlname
        write(*,*) ' grid time is not linear!'

        end if

     end if   


7    continue
      end do


8       continue
      write(*,'(a,a)') ' error in ',ctlname
      write(*,*) ' in definition undefinit value!'
        close (40)
      return

9       continue
      write(*,'(a,a)') ' error in ',ctlname
      write(*,*) ' in definition Options!'
        close (40)
      return

10      continue
      write(*,'(a,a)') ' error in ',ctlname
      write(*,*) ' in definition x-grid!'
        close (40)
      return    

11      continue
      write(*,'(a,a)') ' error in ',ctlname
      write(*,*) ' in definition y-grid!'
        close (40)
      return
                        
12      continue
      write(*,'(a,a)') ' error in ',ctlname
      write(*,*) ' in definition z-grid!'
        close (40)
      return
                        
13      continue
      write(*,'(a,a)') ' error in ',ctlname
      write(*,*) ' in definition t-grid!'
        close (40)
      return

14      continue
      write(*,'(a,a)') ' error in ',ctlname
      write(*,*) ' in definition title of data!'
        close (40)
      return

21      continue
      write(*,*) ' warning in ',ctlname
      write(*,*) ' exit by eof, not by endvars!'
        close (40)          
    return

22      continue
      write(*,'(a,a)') ' error in ',ctlname
      write(*,*) ' in definition vars or variable name!'
        close (40)
      return

endsubroutine ctl_file_read

!====================================================================
integer function nonblank(string,nb1)
!--------------------------------------------------------------------
!     cALCULATE NON BLANK SYMBOL POSITION
!      input : string - STRING
!              nb1    - beginnig position
!      output: nonblank  - NON BLANK POSITION
!--------------------------------------------------------------------
      character*(*)  string
    integer nb,nb1

      do  nb=nb1,len(string)
      if (string(nb:nb).ne.' ') go to 5
      end do

    nonblank=0
    return

5     continue

      nonblank=nb

      if (nb.gt.len(string)-5) then 
     write(*,*) ' warning for non blank position in string:'
     write(*,'(a)') string
      nonblank=0
      end if

endfunction nonblank

!====================================================================
integer function charposition(string,nb1,char0)
!--------------------------------------------------------------------
!     Calculate non blank symbol position
!      INPUT : STRING - string
!              NB1    - BEGINNIG POSITION
!      OUTPUT: CHARPOSITION - position of character CHAR
!--------------------------------------------------------------------
      character*(*)  string
    character char0*1
    integer nb,nb1

      do  nb=nb1,len(string)
      if (string(nb:nb).eq.char0) go to 5
      end do

     write(*,*) ' warning in string:'
     write(*,'(a)') string
      write(*,'(a,a)') ' there is no character ', char0
     charposition=0
    return

5     continue

      charposition=nb

endfunction charposition

!======================================================================
subroutine compare_grids(maxgridlen,ierr,                                        &
                         undefa,nxa,nya,nza,nta,xa,hxa,ya,hya,ta,hta,za,hza,     &
                         undefb,nxb,nyb,nzb,ntb,xb,hxb,yb,hyb,tb,htb,zb,hzb)
    integer maxgridlen,ierr,i
! parameters for temperature
    real  undefa 
      real(8)      xa(maxgridlen),hxa,    &   
                   ya(maxgridlen),hya,    &
                   za(maxgridlen),hza
      real ta,hta 
    integer nxa,nya,nza,nta
 
! parameters for salinity
    real  undefb 
      real(8)      xb(maxgridlen),hxb,    &   
                   yb(maxgridlen),hyb,    &
                   zb(maxgridlen),hzb
      real tb,htb 
    integer nxb,nyb,nzb,ntb

    ierr=0

    if(nxa.gt.maxgridlen) then
    ierr=1
    write(*,'(a,i7,a,i7)') ' numbers of x-grid points ',nxa,' > ',' max be',maxgridlen
    return
    end if

    if(nya.gt.maxgridlen) then
    ierr=1
    write(*,'(a,i7,a,i7)') ' numbers of y-grid points ',nya,' > ',' max be',maxgridlen
    return
    end if
        
    if(nza.gt.maxgridlen) then
    ierr=1
    write(*,'(a,i7,a,i7)') ' numbers of z-grid points ',nza,' > ',' max be',maxgridlen
    return
    end if

    if(nxa.ne.nxb) then
    ierr=1
    write(*,'(a,i7,a,a,i7)') ' numbers of a x-grid points ',nxa,' is not equal to', ' numbers of b x-grid points ',nxb
    return
    end if
    
    if(nya.ne.nyb) then
    ierr=1
    write(*,'(a,i7,a,a,i7)') ' numbers of a y-grid points ',nya,' is not equal to', ' nmubers of b y-grid points ',nyb
    return
    end if

    if(nza.ne.nzb) then
    ierr=1
    write(*,'(a,i7,a,a,i7)') ' numbers of a z-grid points ',nza,' is not equal to', ' numbers of b z-grid points ',nzb
    return
    end if

    if(nta.ne.ntb) then
    ierr=1
    write(*,'(a,i7,a,a,i7)') ' numbers of a y-grid points ',nya,' is not equal to', ' numbers of b y-grid points ',nyb
    return
    end if

    do i=1,nxa
     if(xa(i).ne.xb(i)) then
     write(*,*)' x-grids (a) and (b) are different!'
     ierr=2
     end if
    end do

    do i=1,nya
     if(ya(i).ne.yb(i)) then
     write(*,*)' y-grids (a) and (b) are different!'
     ierr=2
     end if
    end do

    do i=1,nza
     if(za(i).ne.zb(i)) then
     write(*,*)' z-grids (a) and (b) are different!'
     ierr=2
     end if
    end do

endsubroutine compare_grids

endmodule rw_ctl_routes