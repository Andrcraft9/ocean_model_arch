module grid_kernels_module

  use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
  use constants_module, only: RadEarth, EarthAngVel, HeatCapWater, RefDen, FreeFallAcc, pip180
  use config_basinpar_module, only: nx, ny, nz, mm, mmm, nn, nnn, dxst, dyst, rlon, rlat,    &
                                    curve_grid, xgr_type, ygr_type, x_levels, y_levels,      &
                                    rotation_on_lon, rotation_on_lat, x_pole, y_pole, p_pole, q_pole
  use grid_parameters_module, only: grid_parameters_carthesian, grid_parameters_spherical, grid_parameters_curvilinear

  implicit none
  save
  public

#include "macros/mpp_macros.fi"

contains

subroutine lu_init_kernel(bnd_x1, bnd_x2, bnd_y1, bnd_y2,  &
                          mask, lu, lu1)

    integer, intent(in) :: bnd_x1, bnd_x2, bnd_y1, bnd_y2
    integer, intent(in) :: mask(nx, ny)
    real(wp4), intent(inout) :: lu(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                                lu1(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

    integer :: m, n

    do n = bnd_y1, bnd_y2
        do m = bnd_x1, bnd_x2
            if (mask(m,n) == 0) then
                lu(m,n) = 1.0
            endif
        end do
    end do

    lu1 = 1.0

end subroutine

subroutine lu_lv_init_kernel(bnd_x1, bnd_x2, bnd_y1, bnd_y2,  &
                             lu, luh, luu, llu, llv, lcu, lcv)

    integer, intent(in) :: bnd_x1, bnd_x2, bnd_y1, bnd_y2
    real(wp4), intent(in) :: lu(bnd_x1:bnd_x2, bnd_y1:bnd_y2)
    real(wp4), intent(inout) :: luh(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                                luu(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                                llu(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                                llv(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                                lcu(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                                lcv(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

    integer :: m, n

    !do n = ny_start-1,ny_end
    !do m = nx_start-1,nx_end
    do n = bnd_y1, bnd_y2-1
        do m = bnd_x1, bnd_x2-1
        if (lu(m,n)+lu(m+1,n)+lu(m,n+1)+lu(m+1,n+1)>0.5)  then
            luh(m,n)=1.0
        endif

        if(lu(m,n)*lu(m+1,n)*lu(m,n+1)*lu(m+1,n+1)>0.5)  then
            luu(m,n)=1.0
        endif
        enddo
    enddo

    !do n=ny_start-1,ny_end
    !do m=nx_start-1,nx_end
    do n = bnd_y1, bnd_y2-1
        do m = bnd_x1, bnd_x2-1    
        
        if (lu(m,n)+lu(m+1,n)>0.5) then
            llu(m,n)=1.0
        endif
    
        if (lu(m,n)+lu(m,n+1)>0.5) then
            llv(m,n)=1.0
        endif
    
        if (lu(m,n)*lu(m+1,n)>0.5) then
            lcu(m,n)=1.0
        endif
    
        if (lu(m,n)*lu(m,n+1)>0.5) then
            lcv(m,n)=1.0
        endif
                
        enddo
    enddo

end subroutine

subroutine grid_base_init_kernel(nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2,  &
                                 xt, yt, xu, yv,  &
                                 dx, dy, dxt, dyt, dxh, dyh, dxb, dyb,  &
                                 rlh_s, rlh_c)

    integer, intent(in) :: nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2
    real(wp8), intent(inout) :: xt(bnd_x1:bnd_x2), yt(bnd_y1:bnd_y2), xu(bnd_x1:bnd_x2), yv(bnd_y1:bnd_y2)
    real(wp4), intent(inout) :: dx(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                                dy(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                                dxt(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                                dyt(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                                dxh(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                                dyh(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                                dxb(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                                dyb(bnd_x1:bnd_x2, bnd_y1:bnd_y2)
    real(wp4), intent(inout) :: rlh_s(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  &
                                rlh_c(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

    integer :: m, n

    ! temperature grid initialization

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

    ! velocity grid initialization

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
        do n = ny_start - 1, ny_end + 1
            do m = nx_start - 1, nx_end + 1
                !-----initialization of t- and v-grid x-steps in metres
                dxt(m,n) = sngl(xt(m+1)-xt(m))*pip180*RadEarth
                dxb(m,n) = sngl(xt(m+1)-xt(m))*pip180*RadEarth
                !-----initialization of u- and h-grid x-steps in metres
                dx(m,n) =  sngl(xu(m)-xu(m-1))*pip180*RadEarth
                dxh(m,n) = sngl(xu(m)-xu(m-1))*pip180*RadEarth
            end do
        end do
    else
        do n = ny_start - 1, ny_end + 1
            do m = nx_start - 1, nx_end + 1
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
        do n = ny_start - 1, ny_end + 1
            do m = nx_start - 1, nx_end + 1
                !-----initialization of t- and u-grid y-steps in centimeters
                dyt(m,n) = sngl(yt(n+1)-yt(n))*pip180*RadEarth
                dyb(m,n) = sngl(yt(n+1)-yt(n))*pip180*RadEarth
                !-----initialization of v- and h-grid y-steps in centimeters
                dy(m,n) =  sngl(yv(n)-yv(n-1))*pip180*RadEarth
                dyh(m,n) = sngl(yv(n)-yv(n-1))*pip180*RadEarth
            end do
        end do
    else
        do n = ny_start - 1, ny_end + 1
            do m = nx_start - 1, nx_end + 1
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

end subroutine

subroutine grid_geo_init_kernel(nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2,  &
                                xt, yt, xu, yv,  &
                                dxt, dxb, dx, dxh, dyt, dyb, dy, dyh,  &
                                geo_lon_t, geo_lat_t, geo_lon_u, geo_lat_u, geo_lon_v, geo_lat_v, geo_lon_h, geo_lat_h,  &
                                rlh_s, rlh_c, rotvec_coeff,  &
                                sqt, squ, sqv, sqh, rlh_sqh)

    integer, intent(in) :: nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2
    real(wp8), intent(inout) :: xt(bnd_x1:bnd_x2), yt(bnd_y1:bnd_y2), xu(bnd_x1:bnd_x2), yv(bnd_y1:bnd_y2)
    real(wp4), intent(inout) :: dx(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                                dy(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                                dxt(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                                dyt(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                                dxh(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                                dyh(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                                dxb(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                                dyb(bnd_x1:bnd_x2, bnd_y1:bnd_y2)
    real(wp8), intent(inout) :: geo_lon_t(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                                geo_lat_t(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                                geo_lon_u(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                                geo_lat_u(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                                geo_lon_v(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                                geo_lat_v(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                                geo_lon_h(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                                geo_lat_h(bnd_x1:bnd_x2, bnd_y1:bnd_y2)
    real(wp4), intent(inout) :: rlh_s(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  &
                                rlh_c(bnd_x1:bnd_x2, bnd_y1:bnd_y2)
    real(wp8), intent(inout) :: rotvec_coeff(bnd_x1:bnd_x2, bnd_y1:bnd_y2, 4)
    real(wp4), intent(inout) :: sqt(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                                squ(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                                sqv(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                                sqh(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                                rlh_sqh(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

    integer :: m, n

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

end subroutine

end module grid_kernels_module