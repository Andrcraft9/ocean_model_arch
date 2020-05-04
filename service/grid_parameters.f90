module grid_parameters_module
    ! Module description

    use math_tools_module
    use mpp_module, only: mpp_rank
    use constants_module, only: RadEarth, EarthAngVel, HeatCapWater, RefDen, FreeFallAcc, lat_extr

    implicit none
    save
    private

    public :: grid_parameters_carthesian, grid_parameters_spherical, grid_parameters_curvilinear

contains

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

 integer :: bnd_x1,bnd_x2,bnd_y1,bnd_y2
 integer :: mmm_out, mm_out, nnn_out, nn_out
 real(wp8) :: x_mod(bnd_x1:bnd_x2),   &
              y_mod(bnd_y1:bnd_y2),   &
            geo_lon(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
            geo_lat(bnd_x1:bnd_x2,bnd_y1:bnd_y2)
 real(wp4) :: metr_x(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
              metr_y(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
             cor_sin(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
             cor_cos(bnd_x1:bnd_x2,bnd_y1:bnd_y2)     
 real(wp8) :: rot_coef(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4)
 integer :: key_rot, key_cor

 integer :: m,n

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
 endsubroutine grid_parameters_carthesian

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

 integer :: bnd_x1,bnd_x2,bnd_y1,bnd_y2
 integer :: mmm_out, mm_out, nnn_out, nn_out
 real(wp8) :: x_mod(bnd_x1:bnd_x2),   &
              y_mod(bnd_y1:bnd_y2),   &
  geo_lon(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
  geo_lat(bnd_x1:bnd_x2,bnd_y1:bnd_y2)
 real(wp4) ::  metr_x(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
               metr_y(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
              cor_sin(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
              cor_cos(bnd_x1:bnd_x2,bnd_y1:bnd_y2)     
 real(wp8) :: rot_coef(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4)
 integer :: key_rot, key_cor
 real(wp8) :: rotation_on_lon, rotation_on_lat

 ! auxilary variables
 real(wp8) :: sin_lon, sin_lat, cos_lon, cos_lat, lat_mod
 real(wp8) :: free_term_coslon, free_term_sinlon
 real(wp8) :: sinlat_extr, coslat_extr, sum_rot_coef, sum_sincos
 integer :: m,n

  coslat_extr= dcosd(lat_extr)
  sinlat_extr= dsind(lat_extr)

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
 endsubroutine grid_parameters_spherical
 
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

 integer :: bnd_x1,bnd_x2,bnd_y1,bnd_y2
 integer :: mmm_out, mm_out, nnn_out, nn_out
 real(wp8) :: x_mod(bnd_x1:bnd_x2),   &
              y_mod(bnd_y1:bnd_y2),   &
            geo_lon(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
            geo_lat(bnd_x1:bnd_x2,bnd_y1:bnd_y2)
 real(wp4) :: metr_x(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
              metr_y(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
             cor_sin(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
             cor_cos(bnd_x1:bnd_x2,bnd_y1:bnd_y2)     
 real(wp8) :: rot_coef(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4)
 integer key_rot, key_cor
 real(wp8) :: x_pole, y_pole, p_pole, q_pole
 real(wp8) :: y_pole1, q_pole1, r3d,r2d

 ! auxilary variables
 real(wp8) :: sin_lon, sin_lat, cos_lon, cos_lat, lat_mod
 real(wp8) :: a,b,s,t,a0,b0,s0,t0
 real(wp8) :: num1,num2,numa,numb,denom1
 real(wp8) :: numd1,numd2,numd3,numd4,numas,numat,numbs,numbt
 real(wp8) :: alpha_scale
 real(wp8) :: xn,yn,zn,xs,ys,zs,xm,ym,zm,lm,phm, sinphm,coslm,sinlm, phm1
 real(wp8) :: dx_da, dx_db, dy_da, dy_db, da_ds, da_dt,   &
              db_ds, db_dt, ds_dp, ds_dq, dt_dp, dt_dq,   &
              da_dp, da_dq, db_dp, db_dq, dx_dp, dx_dq,   &
              dy_dp, dy_dq, det, hp_divide_r, hq_divide_r
 real(wp8) :: df(2,2), dfm1(2,2)
 real(wp8) :: sinlat_extr, coslat_extr, sum_rot_coef, sum_sincos
 integer :: m,n

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
      
      if (mpp_rank==0) then
       write(*,'(a,f6.4)') 'alpha-scale is ', alpha_scale
      endif

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
 endsubroutine grid_parameters_curvilinear

end module grid_parameters_module