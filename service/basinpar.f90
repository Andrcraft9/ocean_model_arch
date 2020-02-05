module basinpar_module
    ! Initializing basin grid parameters

    use mpp_module, only: mpp_rank, mpp_count, mpp_cart_comm, mpp_size, mpp_coord, mpp_period
    use constants_module, only: RadEarth, EarthAngVel, HeatCapWater, RefDen, FreeFallAcc
    use config_basinpar_module, only: nx, ny, mmm, nnn, dsxt, dyst, rlon, rlat,  &
                                      curve_grid, xgr_type, ygr_type, x_levels, y_levels,  &
                                      rotation_on_lon, rotation_on_lat, x_pole, y_pole, p_pole, q_pole
    use decomposition_module, only: domain_type
    use grid_module, only: grid_type
    use mpp_sync_module, only: sync

    implicit none
    save
    private

    public :: basinpar

contains

!===========================================================================================

    subroutine basinpar(domain, grid_data)
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        integer :: m, n, k, ierr

        ! temperature grid initialization

        ! x-coordinate (in degrees)
        ! in case of regular grid
        do k = 1, domain%bcount
            associate(bnd_x1 => domain%bbnd_x1(k), bnd_x2 => domain%bbnd_x2(k),  &
                      bnd_y1 => domain%bbnd_y1(k), bnd_y2 => domain%bbnd_y2(k),  &
                      xt => grid_data%xt%block(k)%field,  &
                      yt => grid_data%yt%block(k)%field)
  
            if (xgr_type==0) then
                do m = bnd_x1, bnd_x2
                    xt(m) = rlon + dfloat(m-mmm)*dxst
                end do
            else !in case of irregular grid
                do m = bnd_x1, bnd_x2
                    xt(m) = x_levels(m)
                end do
            endif

            ! y-coordinate (in degrees)
            ! in case of regular grid
            if (ygr_type==0) then
                do n = bnd_y1, bnd_y2
                    yt(n) = rlat + dfloat(n-nnn)*dyst
                end do
            else !in case of irregular grid
                do n = bnd_y1, bnd_y2
                    yt(n) = y_levels(n)
                end do
            endif

            end associate
        enddo

        ! parameters:
        if (mpp_rank == 0) then
            write(*,'(2x,a)')' Basin parameters from 1basinpar.inc:'

            if (curve_grid==0) then        ! Carthesian coordinates
                write(*,*) 'Coordinate system is carthesian'
            elseif (curve_grid==1) then
                write(*,*) 'Coordinate system is undistorted sphere'
                write(*,'(a,f10.3)') ' rotation angle on longitude is =',rotation_on_lon,    &
                                     ' rotation angle on  latitude is =',rotation_on_lat
            elseif (curve_grid==2) then
                write(*,*) 'Coordinate system is distorted sphere'
                write(*,'(a,f10.3)') ' geo longitude of new north pole is =',x_pole,    &
                                     ' geo  latitude of new north pole is =',y_pole,    &
                                     ' geo longitude of new south pole is =',p_pole,    &
                                     ' geo  latitude of new south pole is =',q_pole      
            endif

            if (xgr_type==0) then
                write(*,*) 'X-grid is uniform'
                write(*,'(2(a,f10.3),a)') ' initial x-coordinate (m=mmm) =',rlon,' step on x =',dxst,'[dgr] '
            else
                write(*,*) 'X-grid is non-uniform'      
                !write(*,'(a,f10.3)') ' minimal x-coordinate (m=mmm) =',xt(mmm),    &
                !                     ' maximal x-coordinate (m=mm ) =',xt(mm)
            endif

            if (ygr_type==0) then
                write(*,*) 'Y-grid is uniform'
                write(*,'(2(a,f10.3),a)') ' initial y-coordinate (n=nnn) =',rlat,' step on y =',dyst,'[dgr] '
            else
                write(*,*) 'Y-grid is non-uniform'      
                !write(*,'(a,f10.3)') ' minimal y-coordinate (n=nnn) =',yt(nnn),    &
                !                     ' maximal y-coordinate (n=nn ) =',yt(nn)
            endif

            !write(*,'(2(a,i2))') 'Periodicity on X =', periodicity_x,', Periodicity on Y =', periodicity_y
            write(*,'(4(a,i4))') '  nx=',nx, ';  ny=',ny,';  nz=',nz
            write(*,'(4(a,i4))') ' mmm=',mmm,';  mm=',mm,'; nnn=',nnn,';  nn=',nn

            write(*,'(2x,a,g14.7,a)')' Earth radius =',RadEarth,'(m)'
            write(*,'(2x,a,g14.7,a)') 'Earth angular velocity(omega) =',EarthAngVel,'[rad/sec]'
            write(*,'(2x,a,g14.7,a)') 'Heat capacity of water =',HeatCapWater,'[J/kg/�C] for 35%. sal'
            write(*,'(2x,a,g14.7,a)') 'reference density =',RefDen,'[kg/m**3]'
            write(*,'(2x,a,f10.3,a)') 'free fall acceleration(grv)=',FreeFallAcc,'[m/s**2]'
        endif     

        ! velocity grid initialization

        do k = 1, domain%bcount
            associate(bnd_x1 => domain%bbnd_x1(k), bnd_x2 => domain%bbnd_x2(k),  &
                      bnd_y1 => domain%bbnd_y1(k), bnd_y2 => domain%bbnd_y2(k),  &
                      nx_start => domain%bnx_start(k), nx_end => domain%bnx_end(k),  &
                      ny_start => domain%bny_start(k), ny_end => domain%bny_end(k),  &
                      xt => grid_data%xt%block(k)%field,  &
                      yt => grid_data%yt%block(k)%field,  &
                      xu => grid_data%xu%block(k)%field,  &
                      yv => grid_data%yv%block(k)%field,  &
                      ! dx
                      dxt => grid_data%dxt%block(k)%field,  &
                      dxb => grid_data%dxb%block(k)%field,  &
                      dx  => grid_data%dx%block(k)%field,   &
                      dxh => grid_data%dxh%block(k)%field,  &
                      ! dy
                      dyt => grid_data%dyt%block(k)%field,  &
                      dyb => grid_data%dyb%block(k)%field,  &
                      dy  => grid_data%dy%block(k)%field,   &
                      dyh => grid_data%dyh%block(k)%field,  &
                      ! rlh
                      rlh_s => grid_data%rlh_s%block(k)%field,  &
                      rlh_c => grid_data%rlh_c%block(k)%field)

            ! x-coordinate (in degrees)
            do m = bnd_x1, bnd_x2 - 1
                xu(m) = (xt(m)+xt(m+1))/2.0d0
            end do

            ! y-coordinate (in degrees)
            do n = bnd_y1, bnd_y2 - 1
                yv(n) = (yt(n)+yt(n+1))/2.0d0
            end do

            !Initialization of x-steps
            if (xgr_type>0) then
                do n = ny_start, ny_end
                    do m = nx_start, nx_end
                        !-----initialization of t- and v-grid x-steps in metres
                        dxt(m,n) = sngl(xt(m+1)-xt(m))*pip180*RadEarth
                        dxb(m,n) = sngl(xt(m+1)-xt(m))*pip180*RadEarth
                        !-----initialization of u- and h-grid x-steps in metres
                        dx(m,n) =  sngl(xu(m)-xu(m-1))*pip180*RadEarth
                        dxh(m,n) = sngl(xu(m)-xu(m-1))*pip180*RadEarth
                    end do
                end do
            else
                do n = ny_start, ny_end
                    do m = nx_start, nx_end
                        !-----initialization of t- and v-grid x-steps in centimeters
                        dxt(m,n) = sngl(dxst)*pip180*RadEarth
                        dxb(m,n) = sngl(dxst)*pip180*RadEarth
                        !-----initialization of u- and h-grid x-steps in centimeters
                        dx(m,n) =  sngl(dxst)*pip180*RadEarth
                        dxh(m,n) = sngl(dxst)*pip180*RadEarth
                    end do
                end do
            endif

            !Initialization of y-steps
            if (ygr_type>0) then      
                do n = ny_start, ny_end
                    do m = nx_start, nx_end
                        !-----initialization of t- and u-grid y-steps in centimeters
                        dyt(m,n) = sngl(yt(n+1)-yt(n))*pip180*RadEarth
                        dyb(m,n) = sngl(yt(n+1)-yt(n))*pip180*RadEarth
                        !-----initialization of v- and h-grid y-steps in centimeters
                        dy(m,n) =  sngl(yv(n)-yv(n-1))*pip180*RadEarth
                        dyh(m,n) = sngl(yv(n)-yv(n-1))*pip180*RadEarth
                    end do
                end do
            else
                do n = ny_start, ny_end
                    do m = nx_start, nx_end
                        !-----initialization of t- and u-grid y-steps in centimeters
                        dyt(m,n) = sngl(dyst)*pip180*RadEarth
                        dyb(m,n) = sngl(dyst)*pip180*RadEarth
                        !-----initialization of v- and h-grid y-steps in centimeters
                        dy(m,n) =  sngl(dyst)*pip180*RadEarth
                        dyh(m,n) = sngl(dyst)*pip180*RadEarth
                    end do
                end do
            endif

            !-----initialization of Coriolis terms--------------------------       
            rlh_s = 2.0*EarthAngVel
            rlh_c =-2.0*EarthAngVel

            end associate
        enddo

        call sync(domain, grid_data%dxt)
        call sync(domain, grid_data%dxb)
        call sync(domain, grid_data%dx)
        call sync(domain, grid_data%dxh)
        
        call sync(domain, grid_data%dyt)
        call sync(domain, grid_data%dyb)
        call sync(domain, grid_data%dy)
        call sync(domain, grid_data%dyh)


!-----metric initialization-------------------------------------------------------------- 
if(curve_grid==0) then   !in case of carthesian grid

!On T-grid 
call grid_parameters_carthesian(xt,   &   !model x-coordinate in degrees
                    yt,   &   !model y-coordinate in degrees
                bnd_x1,   &   !left   boundary of arrays
                bnd_x2,   &   !right  boundary of arrays
                bnd_y1,   &   !lower  boundary of arrays
                bnd_y2,   &   !upper  boundary of arrays
            geo_lon_t,   &   !geographical longitude in degrees
            geo_lat_t,   &   !geographical latitude  in degrees
                    dx,    &   !metrical coefficient on x
                    dy,    &   !metrical coefficient on x
            rotvec_coeff,  &   !rotation coefficients for vector transform
                rlh_s,   &   !coriolis main term (sin)
                rlh_c,   &   !coriolis second term (cos)
                    1,   &   !key to compute rotation coefficients (0/1)
                    0,   &   !key to compute coriolis terms (0/1)
                nx_start-1, &   !first significant point in x-direction (output)
                nx_end+1,   &   ! last significant point in x-direction (output)
                ny_start-1, &   !first significant point in y-direction (output)
                ny_end+1    )   ! last significant point in y-direction (output)   

!On U-grid    
call grid_parameters_carthesian(xu,   &   !model x-coordinate in degrees
                    yt,   &   !model y-coordinate in degrees
                bnd_x1,   &   !left   boundary of arrays
                bnd_x2,   &   !right  boundary of arrays
                bnd_y1,   &   !lower  boundary of arrays
                bnd_y2,   &   !upper  boundary of arrays
            geo_lon_u,   &   !geographical longitude in degrees
            geo_lat_u,   &   !geographical latitude  in degrees
                    dxt,    &   !metrical coefficient on x
                    dyh,    &   !metrical coefficient on x
            rotvec_coeff,  &   !rotation coefficients for vector transform
                rlh_s,   &   !coriolis main term (sin)
                rlh_c,   &   !coriolis second term (cos)
                    0,   &   !key to compute rotation coefficients (0/1)
                    0,   &   !key to compute coriolis terms (0/1)
                nx_start-1, &   !first significant point in x-direction (output)
                nx_end+1,   &   ! last significant point in x-direction (output)
                ny_start-1, &   !first significant point in y-direction (output)
                ny_end+1    )   ! last significant point in y-direction (output)   

!On V-grid    
call grid_parameters_carthesian(xt,   &   !model x-coordinate in degrees
                    yv,   &   !model y-coordinate in degrees
                bnd_x1,   &   !left   boundary of arrays
                bnd_x2,   &   !right  boundary of arrays
                bnd_y1,   &   !lower  boundary of arrays
                bnd_y2,   &   !upper  boundary of arrays
            geo_lon_v,   &   !geographical longitude in degrees
            geo_lat_v,   &   !geographical latitude  in degrees
                    dxh,    &   !metrical coefficient on x
                    dyt,    &   !metrical coefficient on x
            rotvec_coeff,  &   !rotation coefficients for vector transform
                rlh_s,   &   !coriolis main term (sin)
                rlh_c,   &   !coriolis second term (cos)
                    0,   &   !key to compute rotation coefficients (0/1)
                    0,   &   !key to compute coriolis terms (0/1)
                nx_start-1, &   !first significant point in x-direction (output)
                nx_end+1,   &   ! last significant point in x-direction (output)
                ny_start-1, &   !first significant point in y-direction (output)
                ny_end+1    )   ! last significant point in y-direction (output)         

!On H-grid
call grid_parameters_carthesian(xu,   &   !model x-coordinate in degrees
                    yv,   &   !model y-coordinate in degrees
                bnd_x1,   &   !left   boundary of arrays
                bnd_x2,   &   !right  boundary of arrays
                bnd_y1,   &   !lower  boundary of arrays
                bnd_y2,   &   !upper  boundary of arrays
            geo_lon_h,   &   !geographical longitude in degrees
            geo_lat_h,   &   !geographical latitude  in degrees
                    dxb,    &   !metrical coefficient on x
                    dyb,    &   !metrical coefficient on x
            rotvec_coeff,  &   !rotation coefficients for vector transform
                rlh_s,   &   !coriolis main term (sin)
                rlh_c,   &   !coriolis second term (cos)
                    0,   &   !key to compute rotation coefficients (0/1)
                    1,   &   !key to compute coriolis terms (0/1)
                nx_start-1, &   !first significant point in x-direction (output)
                nx_end+1,   &   ! last significant point in x-direction (output)
                ny_start-1, &   !first significant point in y-direction (output)
                ny_end+1    )   ! last significant point in y-direction (output)      

elseif(curve_grid==1) then !in case of spherical grid

!On T-grid 
call grid_parameters_spherical (xt,   &   !model x-coordinate in degrees
                    yt,   &   !model y-coordinate in degrees
                bnd_x1,   &   !left   boundary of arrays
                bnd_x2,   &   !right  boundary of arrays
                bnd_y1,   &   !lower  boundary of arrays
                bnd_y2,   &   !upper  boundary of arrays
        rotation_on_lon,   &   !euler angle on longitude
        rotation_on_lat,   &   !euler angle on latitude
            geo_lon_t,   &   !geographical longitude in degrees
            geo_lat_t,   &   !geographical latitude  in degrees
                    dx,    &   !metrical coefficient on x
                    dy,    &   !metrical coefficient on x
            rotvec_coeff,  &   !rotation coefficients for vector transform
                rlh_s,   &   !coriolis main term (sin)
                rlh_c,   &   !coriolis second term (cos)
                    1,   &   !key to compute rotation coefficients (0/1)
                    0,   &   !key to compute coriolis terms (0/1)
                nx_start-1, &   !first significant point in x-direction (output)
                nx_end+1,   &   ! last significant point in x-direction (output)
                ny_start-1, &   !first significant point in y-direction (output)
                ny_end+1    )   ! last significant point in y-direction (output)   

!On U-grid    
call grid_parameters_spherical (xu,   &   !model x-coordinate in degrees
                    yt,   &   !model y-coordinate in degrees
                bnd_x1,   &   !left   boundary of arrays
                bnd_x2,   &   !right  boundary of arrays
                bnd_y1,   &   !lower  boundary of arrays
                bnd_y2,   &   !upper  boundary of arrays
        rotation_on_lon,   &   !euler angle on longitude
        rotation_on_lat,   &   !euler angle on latitude
            geo_lon_u,   &   !geographical longitude in degrees
            geo_lat_u,   &   !geographical latitude  in degrees
                    dxt,    &   !metrical coefficient on x
                    dyh,    &   !metrical coefficient on x
            rotvec_coeff,  &   !rotation coefficients for vector transform
                rlh_s,   &   !coriolis main term (sin)
                rlh_c,   &   !coriolis second term (cos)
                    0,   &   !key to compute rotation coefficients (0/1)
                    0,   &   !key to compute coriolis terms (0/1)
                nx_start-1, &   !first significant point in x-direction (output)
                nx_end+1,   &   ! last significant point in x-direction (output)
                ny_start-1, &   !first significant point in y-direction (output)
                ny_end+1    )   ! last significant point in y-direction (output)   

!On V-grid    
call grid_parameters_spherical (xt,   &   !model x-coordinate in degrees
                    yv,   &   !model y-coordinate in degrees
                bnd_x1,   &   !left   boundary of arrays
                bnd_x2,   &   !right  boundary of arrays
                bnd_y1,   &   !lower  boundary of arrays
                bnd_y2,   &   !upper  boundary of arrays
        rotation_on_lon,   &   !euler angle on longitude
        rotation_on_lat,   &   !euler angle on latitude
            geo_lon_v,   &   !geographical longitude in degrees
            geo_lat_v,   &   !geographical latitude  in degrees
                    dxh,    &   !metrical coefficient on x
                    dyt,    &   !metrical coefficient on x
            rotvec_coeff,  &   !rotation coefficients for vector transform
                rlh_s,   &   !coriolis main term (sin)
                rlh_c,   &   !coriolis second term (cos)
                    0,   &   !key to compute rotation coefficients (0/1)
                    0,   &   !key to compute coriolis terms (0/1)
                nx_start-1, &   !first significant point in x-direction (output)
                nx_end+1,   &   ! last significant point in x-direction (output)
                ny_start-1, &   !first significant point in y-direction (output)
                ny_end+1    )   ! last significant point in y-direction (output)         

!On H-grid
call grid_parameters_spherical (xu,   &   !model x-coordinate in degrees
                    yv,   &   !model y-coordinate in degrees
                bnd_x1,   &   !left   boundary of arrays
                bnd_x2,   &   !right  boundary of arrays
                bnd_y1,   &   !lower  boundary of arrays
                bnd_y2,   &   !upper  boundary of arrays
        rotation_on_lon,   &   !euler angle on longitude
        rotation_on_lat,   &   !euler angle on latitude
            geo_lon_h,   &   !geographical longitude in degrees
            geo_lat_h,   &   !geographical latitude  in degrees
                    dxb,    &   !metrical coefficient on x
                    dyb,    &   !metrical coefficient on x
            rotvec_coeff,  &   !rotation coefficients for vector transform
                rlh_s,   &   !coriolis main term (sin)
                rlh_c,   &   !coriolis second term (cos)
                    0,   &   !key to compute rotation coefficients (0/1)
                    1,   &   !key to compute coriolis terms (0/1)
                nx_start-1, &   !first significant point in x-direction (output)
                nx_end+1,   &   ! last significant point in x-direction (output)
                ny_start-1, &   !first significant point in y-direction (output)
                ny_end+1    )   ! last significant point in y-direction (output)      

elseif(curve_grid==2) then   !in case of curve grid

!On T-grid 
call grid_parameters_curvilinear (xt,   &   !model x-coordinate in degrees
                        yt,   &   !model y-coordinate in degrees
                    bnd_x1,   &   !left   boundary of arrays
                    bnd_x2,   &   !right  boundary of arrays
                    bnd_y1,   &   !lower  boundary of arrays
                    bnd_y2,   &   !upper  boundary of arrays
                    x_pole,   &   !geo longitude of new north pole
                    y_pole,   &   !geo latitude  of new north pole
                    p_pole,   &   !geo longitude of new south pole
                    q_pole,   &   !geo latitude  of new south pole
                geo_lon_t,   &   !geographical longitude in degrees
                geo_lat_t,   &   !geographical latitude  in degrees
                    dx,    &   !metrical coefficient on x
                    dy,    &   !metrical coefficient on x
            rotvec_coeff,  &   !rotation coefficients for vector transform
                    rlh_s,   &   !coriolis main term (sin)
                    rlh_c,   &   !coriolis second term (cos)
                        1,   &   !key to compute rotation coefficients (0/1)
                        0,   &   !key to compute coriolis terms (0/1)
                nx_start-1, &   !first significant point in x-direction (output)
                nx_end+1,   &   ! last significant point in x-direction (output)
                ny_start-1, &   !first significant point in y-direction (output)
                ny_end+1    )   ! last significant point in y-direction (output)  

!On U-grid    
call grid_parameters_curvilinear(xu,   &   !model x-coordinate in degrees
                    yt,   &   !model y-coordinate in degrees
                bnd_x1,   &   !left   boundary of arrays
                bnd_x2,   &   !right  boundary of arrays
                bnd_y1,   &   !lower  boundary of arrays
                bnd_y2,   &   !upper  boundary of arrays
                x_pole,   &   !geo longitude of new north pole
                y_pole,   &   !geo latitude  of new north pole
                p_pole,   &   !geo longitude of new south pole
                q_pole,   &   !geo latitude  of new south pole
                geo_lon_u,   &   !geographical longitude in degrees
                geo_lat_u,   &   !geographical latitude  in degrees
                    dxt,    &   !metrical coefficient on x
                    dyh,    &   !metrical coefficient on x
            rotvec_coeff,  &   !rotation coefficients for vector transform
                    rlh_s,   &   !coriolis main term (sin)
                    rlh_c,   &   !coriolis second term (cos)
                        0,   &   !key to compute rotation coefficients (0/1)
                        0,   &   !key to compute coriolis terms (0/1)
                nx_start-1, &   !first significant point in x-direction (output)
                nx_end+1,   &   ! last significant point in x-direction (output)
                ny_start-1, &   !first significant point in y-direction (output)
                ny_end+1    )   ! last significant point in y-direction (output)
            
!On V-grid    
call grid_parameters_curvilinear(xt,   &   !model x-coordinate in degrees
                    yv,   &   !model y-coordinate in degrees
                bnd_x1,   &   !left   boundary of arrays
                bnd_x2,   &   !right  boundary of arrays
                bnd_y1,   &   !lower  boundary of arrays
                bnd_y2,   &   !upper  boundary of arrays
                x_pole,   &   !geo longitude of new north pole
                y_pole,   &   !geo latitude  of new north pole
                p_pole,   &   !geo longitude of new south pole
                q_pole,   &   !geo latitude  of new south pole
                geo_lon_v,   &   !geographical longitude in degrees
                geo_lat_v,   &   !geographical latitude  in degrees
                    dxh,    &   !metrical coefficient on x
                    dyt,    &   !metrical coefficient on x
            rotvec_coeff,  &   !rotation coefficients for vector transform
                    rlh_s,   &   !coriolis main term (sin)
                    rlh_c,   &   !coriolis second term (cos)
                        0,   &   !key to compute rotation coefficients (0/1)
                        0,   &   !key to compute coriolis terms (0/1)
                nx_start-1, &   !first significant point in x-direction (output)
                nx_end+1,   &   ! last significant point in x-direction (output)
                ny_start-1, &   !first significant point in y-direction (output)
                ny_end+1    )   ! last significant point in y-direction (output)         
        
!On H-grid
call grid_parameters_curvilinear(xu,   &   !model x-coordinate in degrees
                    yv,   &   !model y-coordinate in degrees
                bnd_x1,   &   !left   boundary of arrays
                bnd_x2,   &   !right  boundary of arrays
                bnd_y1,   &   !lower  boundary of arrays
                bnd_y2,   &   !upper  boundary of arrays
                x_pole,   &   !geo longitude of new north pole
                y_pole,   &   !geo latitude  of new north pole
                p_pole,   &   !geo longitude of new south pole
                q_pole,   &   !geo latitude  of new south pole
                geo_lon_h,   &   !geographical longitude in degrees
                geo_lat_h,   &   !geographical latitude  in degrees
                    dxb,    &   !metrical coefficient on x
                    dyb,    &   !metrical coefficient on x
            rotvec_coeff,  &   !rotation coefficients for vector transform
                    rlh_s,   &   !coriolis main term (sin)
                    rlh_c,   &   !coriolis second term (cos)
                        0,   &   !key to compute rotation coefficients (0/1)
                        1,   &   !key to compute coriolis terms (0/1)
                nx_start-1, &   !first significant point in x-direction (output)
                nx_end+1,   &   ! last significant point in x-direction (output)
                ny_start-1, &   !first significant point in y-direction (output)
                ny_end+1    )   ! last significant point in y-direction (output)      

end if

! Computing grid areas
do n=ny_start-1,ny_end+1
do m=nx_start-1,nx_end+1
sqt(m,n)=dx(m,n)*dy(m,n)
squ(m,n)=dxt(m,n)*dyh(m,n)
sqv(m,n)=dxh(m,n)*dyt(m,n)
sqh(m,n)=dxb(m,n)*dyb(m,n)
rlh_sqh(m,n)=rlh_s(m,n)*sqh(m,n)
end do
end do

!-----end of metric initialization------------------------------------------

endsubroutine basinpar

!====================================================================
subroutine grid_parameters_carthesian(x_mod,   &   !model x-coordinate in degrees
                    y_mod,   &   !model y-coordinate in degrees
                    bnd_x1,   &   !left   boundary of arrays
                    bnd_x2,   &   !right  boundary of arrays
                    bnd_y1,   &   !lower  boundary of arrays
                    bnd_y2,   &   !upper  boundary of arrays
                    geo_lon,   &   !geographical longitude in degrees
                    geo_lat,   &   !geographical latitude  in degrees
                    metr_x,    &   !metrical coefficient on x
                    metr_y,    &   !metrical coefficient on x
                    rot_coef,  &   !rotation coefficients for vector transform
                    cor_sin,   &   !coriolis main term (sin)
                    cor_cos,   &   !coriolis second term (cos)
                    key_rot,   &   !key to compute rotation coefficients (0/1)
                    key_cor,   &   !key to compute coriolis terms (0/1)
                    mmm_out,   &   !first significant point in x-direction (output)
                    mm_out,    &   ! last significant point in x-direction (output)
                    nnn_out,   &   !first significant point in y-direction (output)
                    nn_out)        !last significant point in y-direction (output)

integer bnd_x1,bnd_x2,bnd_y1,bnd_y2
integer mmm_out, mm_out, nnn_out, nn_out

real(8) x_mod(*),   &
y_mod(*),   &
geo_lon(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
geo_lat(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

real(4)   metr_x(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
metr_y(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
cor_sin(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
cor_cos(bnd_x1:bnd_x2,bnd_y1:bnd_y2)     

real(8) rot_coef(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4)

integer key_rot, key_cor
integer m,n

!$omp parallel do private (m,n)
do n=nnn_out,nn_out
do m=mmm_out,mm_out    

!         necessary latitude
geo_lat(m,n) = y_mod(n)

!         necessary longitude
geo_lon(m,n)= x_mod(m)

metr_x(m,n)=metr_x(m,n)*1.0
metr_y(m,n)=metr_y(m,n)*1.0

if(key_rot==1) then
!--------definition of angles between parallels-----------------------
rot_coef(m,n,1) = 1.0d0
rot_coef(m,n,2) = 0.0d0
rot_coef(m,n,3)=-rot_coef(m,n,2)
rot_coef(m,n,4)= rot_coef(m,n,1)
endif

if(key_cor==1) then
cor_sin(m,n)=cor_sin(m,n)/sqrt(2.0)
cor_cos(m,n)=cor_cos(m,n)/sqrt(2.0)
endif

enddo
enddo
!$omp end parallel do

endsubroutine grid_parameters_carthesian

!====================================================================
subroutine grid_parameters_spherical (x_mod,   &   !model x-coordinate in degrees
                    y_mod,   &   !model y-coordinate in degrees
                    bnd_x1,   &   !left   boundary of arrays
                    bnd_x2,   &   !right  boundary of arrays
                    bnd_y1,   &   !lower  boundary of arrays
                    bnd_y2,   &   !upper  boundary of arrays
            rotation_on_lon,   &   !euler angle on longitude
            rotation_on_lat,   &   !euler angle on latitude
                    geo_lon,   &   !geographical longitude in degrees
                    geo_lat,   &   !geographical latitude  in degrees
                    metr_x,    &   !metrical coefficient on x
                    metr_y,    &   !metrical coefficient on x
                    rot_coef,  &   !rotation coefficients for vector transform
                    cor_sin,   &   !coriolis main term (sin)
                    cor_cos,   &   !coriolis second term (cos)
                    key_rot,   &   !key to compute rotation coefficients (0/1)
                    key_cor,   &   !key to compute coriolis terms (0/1)
                    mmm_out,   &   !first significant point in x-direction (output)
                    mm_out,    &   ! last significant point in x-direction (output)
                    nnn_out,   &   !first significant point in y-direction (output)
                    nn_out)        !last significant point in y-direction (output)
use math_tools

include 'constants.fi'

integer bnd_x1,bnd_x2,bnd_y1,bnd_y2
integer mmm_out, mm_out, nnn_out, nn_out

real(8) x_mod(*),   &
y_mod(*),   &
geo_lon(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
geo_lat(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

real(4)   metr_x(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
metr_y(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
cor_sin(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
cor_cos(bnd_x1:bnd_x2,bnd_y1:bnd_y2)     

real(8) rot_coef(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4)

integer key_rot, key_cor

real(8) rotation_on_lon, rotation_on_lat

real(8) sin_lon, sin_lat, cos_lon, cos_lat, lat_mod       !auxilary variables
real(8) free_term_coslon, free_term_sinlon

real(8) sinlat_extr, coslat_extr, sum_rot_coef, sum_sincos

integer m,n

coslat_extr= dcosd(lat_extr)
sinlat_extr= dsind(lat_extr)

!$omp parallel do private (m,n,sin_lat,cos_lat,free_term_coslon,free_term_sinlon,cos_lon,sin_lon,  &
!$omp                      sum_rot_coef, sum_sincos, lat_mod)
do n=nnn_out,nn_out
lat_mod=max(min(y_mod(n),lat_extr),-lat_extr)
do m=mmm_out,mm_out    
sin_lat = dsind(y_mod(n)) * dcosd(rotation_on_lat)                           &
+ dcosd(x_mod(m)) * dcosd(y_mod(n)) * dsind(rotation_on_lat)

sin_lat=min(max(sin_lat, -sinlat_extr),sinlat_extr)
cos_lat = dsqrt(1d0-sin_lat**2)

!         necessary latitude
geo_lat(m,n) = dasind(sin_lat)

free_term_coslon =(  dcosd(x_mod(m)) * dcosd(y_mod(n)) * dcosd(rotation_on_lat)   &
        -  dsind(y_mod(n)) * dsind(rotation_on_lat)  )  / cos_lat

free_term_sinlon =(  dsind(x_mod(m)) * dcosd(y_mod(n))  ) / cos_lat

cos_lon=free_term_coslon*dcosd(rotation_on_lon)     & 
-free_term_sinlon*dsind(rotation_on_lon)

sin_lon=free_term_sinlon*dcosd(rotation_on_lon)     &
+free_term_coslon*dsind(rotation_on_lon)

sum_sincos=max(dsqrt(cos_lon**2+sin_lon**2),1d-10)
cos_lon=cos_lon/sum_sincos
sin_lon=sin_lon/sum_sincos

!         necessary longitude
geo_lon(m,n)=dsign(dacosd(cos_lon),sin_lon)

metr_x(m,n)=metr_x(m,n)*sngl(dcosd(lat_mod))
metr_y(m,n)=metr_y(m,n)*1.0

if(key_rot==1) then
!--------definition of angles between parallels-----------------------
rot_coef(m,n,1) = (  cos_lat*dcosd(rotation_on_lat) + sin_lat*dsind(rotation_on_lat)    &
        * ( cos_lon*dcosd(rotation_on_lon) + sin_lon*dsind(rotation_on_lon) )  &
                        )  / dcosd(lat_mod)

rot_coef(m,n,2) = (  -dsind(rotation_on_lat)                                                     &
        * ( sin_lon*dcosd(rotation_on_lon) - cos_lon*dsind(rotation_on_lon) )  )      &
                    / dcosd(lat_mod)

rot_coef(m,n,3)=-rot_coef(m,n,2)
rot_coef(m,n,4)= rot_coef(m,n,1)
sum_rot_coef=   max(dsqrt(rot_coef(m,n,1)*rot_coef(m,n,4)-rot_coef(m,n,2)*rot_coef(m,n,3)),1d-10)
rot_coef(m,n,:)=rot_coef(m,n,:)/sum_rot_coef
endif

if(key_cor==1) then
cor_sin(m,n)=cor_sin(m,n)*sngl(sin_lat)
cor_cos(m,n)=cor_cos(m,n)*sngl(cos_lat)
endif

enddo
enddo
!$omp end parallel do

endsubroutine grid_parameters_spherical

!====================================================================
subroutine grid_parameters_curvilinear (x_mod,   &   !model x-coordinate in degrees
                        y_mod,   &   !model y-coordinate in degrees
                    bnd_x1,   &   !left   boundary of arrays
                    bnd_x2,   &   !right  boundary of arrays
                    bnd_y1,   &   !lower  boundary of arrays
                    bnd_y2,   &   !upper  boundary of arrays
                    x_pole,   &   !geo longitude of new north pole
                    y_pole,   &   !geo latitude  of new north pole
                    p_pole,   &   !geo longitude of new south pole
                    q_pole,   &   !geo latitude  of new south pole
                    geo_lon,   &   !geographical longitude in degrees
                    geo_lat,   &   !geographical latitude  in degrees
                    metr_x,    &   !metrical coefficient on x
                    metr_y,    &   !metrical coefficient on x
                    rot_coef,  &   !rotation coefficients for vector transform
                    cor_sin,   &   !coriolis main term (sin)
                    cor_cos,   &   !coriolis second term (cos)
                    key_rot,   &   !key to compute rotation coefficients (0/1)
                    key_cor,   &   !key to compute coriolis terms (0/1)
                    mmm_out,  &   !first significant point in x-direction (output)
                    mm_out,   &   ! last significant point in x-direction (output)
                    nnn_out,  &   !first significant point in y-direction (output)
                    nn_out)       !last significant point in y-direction (output)
use math_tools
use mpi_parallel_tools, only: rank

include 'constants.fi'

integer bnd_x1,bnd_x2,bnd_y1,bnd_y2
integer mmm_out, mm_out, nnn_out, nn_out

real(8) x_mod(*),   &
y_mod(*),   &
geo_lon(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
geo_lat(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

real(4)   metr_x(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
metr_y(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
cor_sin(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
cor_cos(bnd_x1:bnd_x2,bnd_y1:bnd_y2)     

real(8) rot_coef(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4)

integer key_rot, key_cor

real(8) x_pole, y_pole, p_pole, q_pole
real(8) y_pole1, q_pole1, r3d,r2d

real(8) sin_lon, sin_lat, cos_lon, cos_lat, lat_mod    !auxilary variables
real(8) a,b,s,t,a0,b0,s0,t0                    !auxilary variables

real(8) num1,num2,numa,numb,denom1
real(8) numd1,numd2,numd3,numd4,numas,numat,numbs,numbt
real(8) alpha_scale
real(8) xn,yn,zn,xs,ys,zs,xm,ym,zm,lm,phm, sinphm,coslm,sinlm, phm1

real(8) dx_da, dx_db, dy_da, dy_db, da_ds, da_dt,   &
db_ds, db_dt, ds_dp, ds_dq, dt_dp, dt_dq,   &
da_dp, da_dq, db_dp, db_dq, dx_dp, dx_dq,   &
dy_dp, dy_dq, det, hp_divide_r, hq_divide_r
real(8) df(2,2), dfm1(2,2)

real(8) sinlat_extr, coslat_extr, sum_rot_coef, sum_sincos

integer m,n

y_pole1=min(max(y_pole, -lat_extr),lat_extr)
q_pole1=min(max(q_pole, -lat_extr),lat_extr)

coslat_extr= dcosd(lat_extr)
sinlat_extr= dsind(lat_extr)

xn=dcosd(x_pole)*dcosd(y_pole)
yn=dsind(x_pole)*dcosd(y_pole)
zn=dsind(y_pole)

xs=dcosd(p_pole)*dcosd(q_pole)
ys=dsind(p_pole)*dcosd(q_pole)
zs=dsind(q_pole)

xm=(xn+xs)/2.0d0
ym=(yn+ys)/2.0d0
zm=(zn+zs)/2.0d0

r3d=max(dsqrt(xm**2+ym**2+zm**2),1d-10)
r2d=max(dsqrt(xm**2+ym**2),1d-10)

sinphm=zm/r3d
sinlm= ym/r2d
coslm= xm/r2d

sinphm=min(max(sinphm, -sinlat_extr),sinlat_extr)
phm=dasind(sinphm)

sum_sincos=max(dsqrt(coslm**2+sinlm**2),1d-10)
coslm=coslm/sum_sincos
sinlm=sinlm/sum_sincos

lm=dsign(dacosd(coslm),sinlm)

s0 = 2.0d0*dtand(45.0d0 + y_pole1/2.0d0) *dcosd(x_pole)
t0 = 2.0d0*dtand(45.0d0 + y_pole1/2.0d0) *dsind(x_pole)
a0 = 2.0d0*dtand(45.0d0 + q_pole1/2.0d0) *dcosd(p_pole)
b0 = 2.0d0*dtand(45.0d0 + q_pole1/2.0d0) *dsind(p_pole)

alpha_scale=1.0d0

phm1=min(max(phm, -lat_extr),lat_extr)

S = 2.0d0*dtand(45.0d0 + phm1/2.0d0) *dcosd(lm)
T = 2.0d0*dtand(45.0d0 + phm1/2.0d0) *dsind(lm)

num1=(S-alpha_scale*A0)*(S-alpha_scale*S0) + (T-alpha_scale*B0)*(T-alpha_scale*T0)
num2=(T-alpha_scale*B0)*(S-alpha_scale*S0) - (S-alpha_scale*A0)*(T-alpha_scale*T0)

numa=s0*num1-t0*num2
numb=s0*num2+t0*num1

denom1=(S-alpha_scale*S0)**2+(T-alpha_scale*T0)**2

a=numa/denom1
b=numb/denom1

alpha_scale=2./dsqrt(a**2+b**2)

if(rank==0) then
write(*,'(a,f6.4)') 'alpha-scale is ', alpha_scale
endif

!$omp parallel do private(m,n,s,t,num1,num2,numa,numb,denom1,  &
!$omp      a,b,cos_lon,sin_lon,cos_lat,sin_lat, &
!$omp      dx_da,dx_db,dy_da,dy_db,numd1,numd2, &
!$omp      numd3,numd4,numas,numat,numbs,numbt, &
!$omp      da_ds,da_dt,db_ds,db_dt,ds_dp,ds_dq,dt_dp,dt_dq,  &
!$omp      da_dp,da_dq,db_dp,db_dq,dx_dp,dx_dq,dy_dp,dy_dq,  &
!$omp      dfm1,det,df,hp_divide_r,hq_divide_r, sum_rot_coef, sum_sincos, lat_mod)
do n=nnn_out,nn_out
lat_mod=max(min(y_mod(n),lat_extr),-lat_extr)
do m=mmm_out,mm_out     

!trasformation from new to old grid
s = 2.0d0*dtand(45.0d0 + lat_mod/2.0d0) *dcosd(x_mod(m))
t = 2.0d0*dtand(45.0d0 + lat_mod/2.0d0) *dsind(x_mod(m))

num1=(S-alpha_scale*A0)*(S-alpha_scale*S0) + (T-alpha_scale*B0)*(T-alpha_scale*T0)
num2=(T-alpha_scale*B0)*(S-alpha_scale*S0) - (S-alpha_scale*A0)*(T-alpha_scale*T0)

numa=s0*num1-t0*num2
numb=s0*num2+t0*num1

denom1=(S-alpha_scale*S0)**2+(T-alpha_scale*T0)**2

a=numa/denom1
b=numb/denom1

!     necessary latitude

sin_lat=(a**2+b**2-4.0d0)/(a**2+b**2+4.0d0)
sin_lat=min(max(sin_lat, -sinlat_extr),sinlat_extr)
cos_lat=dsqrt(1.0d0-sin_lat**2)

geo_lat(m,n)=dasind(sin_lat)

!     necessary longitude

cos_lon=a/dsqrt(a**2+b**2)
sin_lon=b/dsqrt(a**2+b**2)

sum_sincos=max(dsqrt(cos_lon**2+sin_lon**2),1d-10)
cos_lon=cos_lon/sum_sincos
sin_lon=sin_lon/sum_sincos

geo_lon(m,n)=dsign(dacosd(cos_lon),sin_lon)


!     differential of transformation

dx_da = -b / (a**2 + b**2)
dx_db =  a / (a**2 + b**2)

dy_da = a / ( dsqrt(a**2 + b**2) * (1.0d0+ (a**2 + b**2)/4.0d0))
dy_db = b / ( dsqrt(a**2 + b**2) * (1.0d0+ (a**2 + b**2)/4.0d0))


numd1=s-alpha_scale*s0+s-alpha_scale*a0
numd2=t-alpha_scale*t0+t-alpha_scale*b0
numd3=alpha_scale*(t0-b0)
numd4=alpha_scale*(a0-s0)

numas=s0*numd1-t0*numd3
numat=s0*numd2-t0*numd4
numbs=t0*numd1+s0*numd3
numbt=t0*numd2+s0*numd4

da_ds=numas/denom1-numa*2.0d0*(s-alpha_scale*s0)/(denom1*denom1)
da_dt=numat/denom1-numa*2.0d0*(t-alpha_scale*t0)/(denom1*denom1)
db_ds=numbs/denom1-numb*2.0d0*(s-alpha_scale*s0)/(denom1*denom1)
db_dt=numbt/denom1-numb*2.0d0*(t-alpha_scale*t0)/(denom1*denom1)


ds_dp = -2.0d0*dtand(45.0d0 + lat_mod/2.0d0) *dsind(x_mod(m))
ds_dq = dcosd(x_mod(m)) /(dcosd(45.0d0 + lat_mod/2.0d0))**2
dt_dp =  2.0d0*dtand(45.0d0 + lat_mod/2.0d0) *dcosd(x_mod(m))
dt_dq = dsind(x_mod(m)) /(dcosd(45.0d0 + lat_mod/2.0d0))**2

da_dp = da_ds*ds_dp + da_dt*dt_dp
da_dq = da_ds*ds_dq + da_dt*dt_dq
db_dp = db_ds*ds_dp + db_dt*dt_dp
db_dq = db_ds*ds_dq + db_dt*dt_dq

dx_dp = dx_da*da_Dp + dx_db*db_dp
dx_dq = dx_da*da_Dq + dx_db*db_dq
dy_dp = dy_da*da_Dp + dy_db*db_dp
dy_dq = dy_da*da_Dq + dy_db*db_dq

dfm1(1,1)=dx_dp   !*dcos(ret_lat*dpip180)
dfm1(1,2)=dx_dq   !*dcos(ret_lat*dpip180)
dfm1(2,1)=dy_dp
dfm1(2,2)=dy_dq


det=dfm1(2,2)*dfm1(1,1)-dfm1(1,2)*dfm1(2,1)

df(1,1)= dfm1(2,2)/det
df(1,2)=-dfm1(1,2)/det
df(2,1)=-dfm1(2,1)/det
df(2,2)= dfm1(1,1)/det

hp_divide_r = dsqrt((dx_dp*cos_lat)**2 + (dy_dp)**2)
hq_divide_r = dsqrt((dx_dq*cos_lat)**2 + (dy_dq)**2)

metr_x(m,n)=metr_x(m,n)*sngl(hp_divide_r)
metr_y(m,n)=metr_y(m,n)*sngl(hq_divide_r)

if(key_rot==1) then
!--------definition of angles between parallels-----------------------      
rot_coef(m,n,1)=df(1,1)*hp_divide_r/cos_lat
rot_coef(m,n,2)=df(1,2)*hp_divide_r
rot_coef(m,n,3)=df(2,1)*hq_divide_r/cos_lat
rot_coef(m,n,4)=df(2,2)*hq_divide_r      

sum_rot_coef=  max(dsqrt(rot_coef(m,n,1)*rot_coef(m,n,4)-rot_coef(m,n,2)*rot_coef(m,n,3)),1d-10)
rot_coef(m,n,:)=rot_coef(m,n,:)/sum_rot_coef
endif

if(key_cor==1) then
cor_sin(m,n)=cor_sin(m,n)*sngl(sin_lat)
cor_cos(m,n)=cor_cos(m,n)*sngl(cos_lat)
endif

enddo
enddo
!$omp end parallel do

endsubroutine grid_parameters_curvilinear

end module basinpar_module
