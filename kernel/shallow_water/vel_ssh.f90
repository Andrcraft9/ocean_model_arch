module velssh_sw_module

  use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
  use kernel_interface_module, only: nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2

  implicit none
  save
  private

  public :: uv_trans_vort_kernel, uv_trans_kernel, uv_diff2_kernel

contains

subroutine uv_trans_vort_kernel(luu,                 &
                                dxt, dyt, dxb, dyb,  &
                                u, v, vort,          &
                                nlev)

 real(wp4), intent(in) :: luu(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

 real(wp4), intent(in) :: dxt(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          dyt(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          dxb(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          dyb(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

 integer, intent(in) :: nlev
 real(wp8), intent(inout) :: vort(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev)
 real(wp8), intent(in) :: u(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev),        & !Transporting zonal velocity
                          v(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev)           !Transporting meridional velocity

 integer :: m, n, k

 do n=ny_start, ny_end
   do m=nx_start, nx_end
    if(luu(m,n)>0.5) then
     do k=1,nlev
      vort(m,n,k)= (v(m+1,n,k)*dyt(m+1,n)-v(m,n,k)*dyt(m,n))     &
                  -(u(m,n+1,k)*dxt(m,n+1)-u(m,n,k)*dxt(m,n))     &
                  -((v(m+1,n,k)-v(m,n,k))*dyb(m,n)-(u(m,n+1,k)-u(m,n,k))*dxb(m,n))
     enddo
    endif
   enddo
 enddo

end subroutine uv_trans_vort_kernel

subroutine uv_trans_kernel(lcu, lcv, luu,   &
                           dxh, dyh,        &
                           u, v, vort,      &
                           hq, hu, hv, hh,  &
                           RHSx, RHSy, nlev)

 integer, intent(in) :: nlev

 real(wp4), intent(in) :: lcu(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                          lcv(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          luu(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

 real(wp4), intent(in) :: dxh(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          dyh(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

 real(wp8), intent(in) :: u(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev),        & !Transporting zonal velocity
                          v(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev)           !Transporting meridional velocity

 real(wp8), intent(inout) :: RHSx(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev ),      & !Zonal source function
                             RHSy(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev )         !meridional source function

 real(wp8), intent(in) :: hq(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                          hu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                          hv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                          hh(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

 real(wp8), intent(in) :: vort(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev)

 real(wp8) :: fx_p, fx_m, fy_p, fy_m   !fluxes through cell edges

 integer :: m, n, k

  do n=ny_start,ny_end
    do m=nx_start,nx_end

!zonal velocity
      if(lcu(m,n)>0.5) then

        do k=1,nlev

         fx_p=(u(m  ,n  ,k)*dyh(m,n)*hu(m,n) + u(m+1,n  ,k)*dyh(m+1,n)*hu(m+1,n))/2.0d0   &
             *(u(m  ,n  ,k) + u(m+1,n  ,k))/2.0d0

         fx_m=(u(m  ,n  ,k)*dyh(m,n)*hu(m,n) + u(m-1,n  ,k)*dyh(m-1,n)*hu(m-1,n))/2.0d0   &
             *(u(m  ,n  ,k) + u(m-1,n  ,k))/2.0d0

         fy_p=(v(m  ,n  ,k)*dxh(m,n  )*hv(m,n  ) + v(m+1,n  ,k)*dxh(m+1,n  )*hv(m+1,n  ))/2.0d0   &
             *(u(m  ,n+1,k) + u(m  ,n  ,k))/2.0d0*dble(luu(m,n  ))

         fy_m=(v(m  ,n-1,k)*dxh(m,n-1)*hv(m,n-1) + v(m+1,n-1,k)*dxh(m+1,n-1)*hv(m+1,n-1))/2.0d0   &
             *(u(m  ,n-1,k) + u(m  ,n  ,k))/2.0d0*dble(luu(m,n-1))

         RHSx(m,n,k)= - (fx_p - fx_m + fy_p - fy_m)           &
            + ( vort(m,n  ,k)*hh(m,n  )*(v(m+1,n  ,k)+v(m,n  ,k))              &
            +   vort(m,n-1,k)*hh(m,n-1)*(v(m+1,n-1,k)+v(m,n-1,k))  )/4.0d0

        end do

      end if

!meridional velocity
      if(lcv(m,n)>0.5) then

        do k=1,nlev

         fy_p=(v(m  ,n  ,k)*dxh(m,n)*hv(m,n) + v(m  ,n+1,k)*dxh(m,n+1)*hv(m,n+1))/2.0d0    &
             *(v(m  ,n  ,k) + v(m  ,n+1,k))/2.0d0

         fy_m=(v(m  ,n  ,k)*dxh(m,n)*hv(m,n) + v(m  ,n-1,k)*dxh(m,n-1)*hv(m,n-1))/2.0d0    &
             *(v(m  ,n  ,k) + v(m  ,n-1,k))/2.0d0

         fx_p=(u(m  ,n  ,k)*dyh(m  ,n)*hu(m  ,n) + u(m  ,n+1,k)*dyh(m  ,n+1)*hu(m  ,n+1))/2.0d0    &
             *(v(m+1,n  ,k) + v(m  ,n  ,k))/2.0d0

         fx_m=(u(m-1,n  ,k)*dyh(m-1,n)*hu(m-1,n) + u(m-1,n+1,k)*dyh(m-1,n+1)*hu(m-1,n+1))/2.0d0    &
             *(v(m-1,n  ,k) + v(m  ,n  ,k))/2.0d0

         RHSy(m,n,k)= - (fx_p - fx_m + fy_p - fy_m)          &
             - ( vort(m  ,n,k)*hh(m  ,n)*(u(m  ,n+1,k)+u(m  ,n,k))               &
             +   vort(m-1,n,k)*hh(m-1,n)*(u(m-1,n+1,k)+u(m-1,n,k))  )/4.0d0
        end do

      end if

    end do
  end do

endsubroutine uv_trans_kernel

subroutine uv_diff2_kernel(lcu, lcv,                              &
                           dx, dy, dxt, dyt, dxh, dyh, dxb, dyb,  &
                           mu, str_t, str_s,                      &
                           hq, hu, hv, hh,                        &
                           RHSx, RHSy, nlev)

 integer, intent(in) :: nlev
 real(wp4), intent(in) :: lcu(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  &
                          lcv(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

 real(wp4), intent(in) :: dx(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                          dy(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                          dxt(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          dyt(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          dxh(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          dyh(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          dxb(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          dyb(bnd_x1:bnd_x2, bnd_y1:bnd_y2)


 real(wp8), intent(in) :: mu(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev ),      & !lateral viscosity coefficient
                       str_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev ),      & !Tension stress
                       str_s(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev )         !Shearing stress

 real(wp8), intent(inout) :: RHSx(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev),  & !Zonal source function
                             RHSy(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev)     !meridional source function

 real(wp8), intent(in) :: hq(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                          hu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                          hv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                          hh(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

 real(wp8) :: muh_p, muh_m
 integer :: m, n, k

  do n=ny_start,ny_end
    do m=nx_start,nx_end

!zonal velocity
      if(lcu(m,n)>0.5) then

        do k=1,nlev

         muh_p=(mu(m,n,k)+mu(m+1,n,k)+mu(m,n+1,k)+mu(m+1,n+1,k))/4.0d0
         muh_m=(mu(m,n,k)+mu(m+1,n,k)+mu(m,n-1,k)+mu(m+1,n-1,k))/4.0d0

         RHSx(m,n,k)=( dy(m+1,n)**2*mu(m+1,n,k)*hq(m+1,n)*str_t(m+1,n,k)             &
                      -dy(m  ,n)**2*mu(m  ,n,k)*hq(m  ,n)*str_t(m  ,n,k) )/dyh(m,n)  &
                   + (dxb(m,n  )**2*muh_p*hh(m,n  )*str_s(m,n  ,k)                   &
                     -dxb(m,n-1)**2*muh_m*hh(m,n-1)*str_s(m,n-1,k) )/dxt(m,n)
        end do

      end if

!meridional velocity
      if(lcv(m,n)>0.5) then

        do k=1,nlev

         muh_p=(mu(m,n,k)+mu(m+1,n,k)+mu(m,n+1,k)+mu(m+1,n+1,k))/4.0d0
         muh_m=(mu(m,n,k)+mu(m-1,n,k)+mu(m,n+1,k)+mu(m-1,n+1,k))/4.0d0

         RHSy(m,n,k)=-( dx(m,n+1)**2*mu(m,n+1,k)*hq(m,n+1)*str_t(m,n+1,k)              &
                       -dx(m,n  )**2*mu(m,n  ,k)*hq(m,n  )*str_t(m,n  ,k) ) /dxh(m,n)  &
                    + (dyb(m  ,n)**2*muh_p*hh(m  ,n)*str_s(m  ,n,k)                    &
                      -dyb(m-1,n)**2*muh_m*hh(m-1,n)*str_s(m-1,n,k) ) /dyt(m,n)
        end do

      end if

    end do
  end do

endsubroutine uv_diff2_kernel

endmodule velssh_sw_module
