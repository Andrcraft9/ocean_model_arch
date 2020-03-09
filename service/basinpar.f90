module basinpar_module
#include "core/kernel_macros.fi"
    ! Initializing basin grid parameters

    use mpp_module, only: mpp_rank, mpp_count, mpp_cart_comm, mpp_size, mpp_coord, mpp_period
    use constants_module, only: RadEarth, EarthAngVel, HeatCapWater, RefDen, FreeFallAcc, pip180
    use config_basinpar_module, only: nx, ny, nz, mm, mmm, nn, nnn, dxst, dyst, rlon, rlat,    &
                                      curve_grid, xgr_type, ygr_type, x_levels, y_levels,      &
                                      rotation_on_lon, rotation_on_lat, x_pole, y_pole, p_pole, q_pole
    use decomposition_module, only: domain_type
    use grid_module, only: grid_type
    use mpp_sync_module, only: sync
    use grid_parameters_module, only: grid_parameters_carthesian, grid_parameters_spherical, grid_parameters_curvilinear

    implicit none
    save
    private

    public :: basinpar

contains

!===========================================================================================

    subroutine basinpar(domain, grid_data)
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        
        integer :: m, n, k

        ! temperature grid initialization

        ! x-coordinate (in degrees)
        ! in case of regular grid
        do k = 1, domain%bcount
            associate(_associate_domain_value_(bnd_x1, domain, bbnd_x1, k),  &
                      _associate_domain_value_(bnd_x2, domain, bbnd_x2, k),  &
                      _associate_domain_value_(bnd_y1, domain, bbnd_y1, k),  &
                      _associate_domain_value_(bnd_y2, domain, bbnd_y2, k),  &
                      _associate_data_(grid_data, xt, k),  &
                      _associate_data_(grid_data, yt, k))
  
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
                      ! geo
                      geo_lon_t => grid_data%geo_lon_t%block(k)%field, &
                      geo_lat_t => grid_data%geo_lat_t%block(k)%field, &
                      geo_lon_u => grid_data%geo_lon_u%block(k)%field, &
                      geo_lat_u => grid_data%geo_lat_u%block(k)%field, &
                      geo_lon_v => grid_data%geo_lon_v%block(k)%field, &
                      geo_lat_v => grid_data%geo_lat_v%block(k)%field, &
                      geo_lon_h => grid_data%geo_lon_h%block(k)%field, &
                      geo_lat_h => grid_data%geo_lat_h%block(k)%field, &
                      ! rlh
                      rlh_s => grid_data%rlh_s%block(k)%field,  &
                      rlh_c => grid_data%rlh_c%block(k)%field,  &
                      ! rotvec
                      rotvec_coeff => grid_data%rotvec_coeff%block(k)%field)

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
            end associate
        enddo

        do k = 1, domain%bcount
            associate(_associate_domain_value_(nx_start, domain, bnx_start, k),  &
                      _associate_domain_value_(nx_end,   domain, bnx_end,   k),  &
                      _associate_domain_value_(ny_start, domain, bny_start, k),  &
                      _associate_domain_value_(ny_end,   domain, bny_end,   k),  &
                      ! dx
                      _associate_data_(grid_data, dxt, k),  &
                      _associate_data_(grid_data, dxb, k),  &
                      _associate_data_(grid_data, dx,  k),   &
                      _associate_data_(grid_data, dxh, k),  &
                      ! dy
                      _associate_data_(grid_data, dyt, k),  &
                      _associate_data_(grid_data, dyb, k),  &
                      _associate_data_(grid_data, dy,  k),   &
                      _associate_data_(grid_data, dyh, k),  &
                      ! sq
                      _associate_data_(grid_data, sqt, k),  &
                      _associate_data_(grid_data, squ, k),  &
                      _associate_data_(grid_data, sqv, k),  &
                      _associate_data_(grid_data, sqh, k),  &
                      ! rlh
                      _associate_data_(grid_data, rlh_s,   k),  &
                      _associate_data_(grid_data, rlh_sqh, k))
                      
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
                
            end associate
        enddo
endsubroutine basinpar

end module basinpar_module
