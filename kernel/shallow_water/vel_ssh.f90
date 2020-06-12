module velssh_sw_module

  use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
  use constants_module, only: FreeFallAcc, dPi

  implicit none
  save
  private

#include "macros/mpp_macros.fi"

  public :: gaussian_elimination_kernel
  public :: check_ssh_err_kernel, sw_update_ssh_kernel, sw_update_uv, sw_next_step, uv_trans_vort_kernel, uv_trans_kernel, uv_diff2_kernel

contains


subroutine gaussian_elimination_kernel(nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2, lu, ssh, sigma, nx0, ny0)

  integer, intent(in) :: nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2

  real(wp4), intent(in) :: lu(bnd_x1:bnd_x2, bnd_y1:bnd_y2)
  real(wp8), intent(inout) :: ssh(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

  real(wp8), intent(in) ::sigma
  integer, intent(in) :: nx0, ny0

  integer :: m, n, ierr
  real(wp8) :: dx, dy

  _OMP_KERNEL_PARALLEL_BEGIN_
  do n = ny_start, ny_end
    do m = nx_start, nx_end
        if (lu(m,n)>0.5) then
          dx = real((m - nx0), kind=wp8) / (nx0*0.25d0)
          dy = real((n - ny0), kind=wp8) / (ny0*0.25d0)
          ssh(m, n) = (1.0 / (dsqrt(2*dPi) * sigma)) * dexp( -( (dx*dx + dy*dy) / (2*sigma*sigma )) )
        endif
    enddo
  enddo
  _OMP_KERNEL_PARALLEL_END_

end subroutine

subroutine check_ssh_err_kernel(nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2,  &
                                lu, ssh, name)
  use mpp_module
  use errors_module, only: abort_model
  
  integer, intent(in) :: nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2

  real(wp4), intent(in) :: lu(bnd_x1:bnd_x2, bnd_y1:bnd_y2)
  real(wp8), intent(in) :: ssh(bnd_x1:bnd_x2,bnd_y1:bnd_y2)
  character(*), intent(in) :: name
  integer :: m, n, ierr

  do n = ny_start, ny_end
      do m = nx_start, nx_end
          if (lu(m,n)>0.5) then
              if (ssh(m,n)<10000.0d0 .and. ssh(m,n)>-10000.0d0) then
                continue
              else
                  write(*,*) mpp_rank, 'ERROR!!! In the point m=', m, 'n=', n, name, '=', ssh(m,n)
                  !write(*,*) rank, 'ERR: Block k=', k, 'In the point m=', m, 'n=', n, 'ssh=', ssh(k)%vals(m,n),   &
                  !    'step: ', num_step, 'lon: ', geo_lon_t(k)%vals(m, n), 'lat: ', geo_lat_t(k)%vals(m, n)

                  call abort_model('SIGFPRE predict error')
              endif
          endif
      enddo
  enddo
end subroutine

subroutine sw_update_ssh_kernel(nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2,  &
                                tau, lu, dx, dy, dxh, dyh, hhu, hhv, sshn, sshp, ubrtr, vbrtr)

  integer, intent(in) :: nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2

  real(wp8), intent(in) :: tau

  real(wp4), intent(in) :: lu(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

  real(wp4), intent(in) :: dx(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                           dy(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                           dxh(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  &
                           dyh(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

 real(wp8), intent(in) :: hhu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),  &
                          hhv(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

  real(wp8), intent(inout) :: sshn(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

  real(wp8), intent(in) :: sshp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
                           ubrtr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),  &
                           vbrtr(bnd_x1:bnd_x2,bnd_y1:bnd_y2)
  
  integer :: m, n

  _OMP_KERNEL_PARALLEL_BEGIN_
  do n=ny_start,ny_end
    do m=nx_start,nx_end

        if(lu(m,n)>0.5) then
            sshn(m,n) = sshp(m,n) + 2.0d0*tau*(  &
            - ( ubrtr(m,n)*hhu(m,n)*dyh(m,n) - ubrtr(m-1,n)*hhu(m-1,n)*dyh(m-1,n)             &
              + vbrtr(m,n)*hhv(m,n)*dxh(m,n) - vbrtr(m,n-1)*hhv(m,n-1)*dxh(m,n-1) )/(dx(m,n)*dy(m,n))  )
        endif

    enddo
  enddo
  _OMP_KERNEL_PARALLEL_END_

end subroutine

subroutine sw_update_uv(nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2,  &
                        tau, lcu, lcv,  &
                        dxt, dyt, dxh, dyh, dxb, dyb,  &
                        hhu, hhun, hhup,  &
                        hhv, hhvn, hhvp,  &
                        hhh, ssh,  &
                        ubrtr, ubrtrn, ubrtrp, vbrtr, vbrtrn, vbrtrp,  &
                        rdis, rlh_s,  &
                        RHSx, RHSy, RHSx_adv, RHSy_adv, RHSx_dif, RHSy_dif)

  integer, intent(in) :: nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2

  real(wp8), intent(in) :: tau
  
  real(wp4), intent(in) :: lcu(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  &
                           lcv(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

  real(wp4), intent(in) :: dxt(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                           dyt(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                           dxh(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                           dyh(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                           dxb(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                           dyb(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

  real(wp8), intent(in) :: hhu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
                           hhun(bnd_x1:bnd_x2,bnd_y1:bnd_y2),  &
                           hhup(bnd_x1:bnd_x2,bnd_y1:bnd_y2),  &
                           hhv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
                           hhvn(bnd_x1:bnd_x2,bnd_y1:bnd_y2),  &
                           hhvp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),  &
                           hhh(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

  real(wp8), intent(in) :: ssh(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

  real(wp8), intent(inout) :: ubrtrn(bnd_x1:bnd_x2,bnd_y1:bnd_y2),  & 
                              vbrtrn(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

  real(wp8), intent(in) :: ubrtr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
                           vbrtr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
                           ubrtrp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),  &
                           vbrtrp(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

  real(wp4), intent(in) :: rdis(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                           rlh_s(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

 real(wp8), intent(in) :: RHSx(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                          RHSy(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                          RHSx_adv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),  &
                          RHSy_adv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),  &
                          RHSx_dif(bnd_x1:bnd_x2,bnd_y1:bnd_y2),  &
                          RHSy_dif(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

  integer :: m, n
  real(wp8) :: bp, bp0, slx, sly, grx, gry

  _OMP_KERNEL_PARALLEL_BEGIN_
  do n=ny_start,ny_end
    do m=nx_start,nx_end
        !zonal flux
        if(lcu(m,n)>0.5) then
            bp  = hhun(m,n)*dxt(m,n)*dyh(m,n)/2.0d0/tau
            bp0 = hhup(m,n)*dxt(m,n)*dyh(m,n)/2.0d0/tau

            slx = - FreeFallAcc * ( ssh(m+1,n) - ssh(m,n))*dyh(m,n)* hhu(m,n)
            grx = RHSx(m,n) + slx + RHSx_dif(m,n) + RHSx_adv(m,n)      &
                - (rdis(m,n)+rdis(m+1,n))/2.0d0 *ubrtrp(m,n)*dxt(m,n)*dyh(m,n)*hhu(m,n)        &
                + ( rlh_s(m,n  )*hhh(m,n  )*dxb(m,n  )*dyb(m,n  )*(vbrtr(m+1,n  ) + vbrtr(m,n  ))          &
                +   rlh_s(m,n-1)*hhh(m,n-1)*dxb(m,n-1)*dyb(m,n-1)*(vbrtr(m+1,n-1) + vbrtr(m,n-1)) )/4.0d0

            ubrtrn(m,n) = (ubrtrp(m,n)*bp0 + grx )/(bp)
        endif

        !meridional flux
        if(lcv(m,n)>0.5) then
            bp  = hhvn(m,n)*dyt(m,n)*dxh(m,n)/2.0d0/tau
            bp0 = hhvp(m,n)*dyt(m,n)*dxh(m,n)/2.0d0/tau

            sly = - FreeFallAcc * ( ssh(m,n+1)- ssh(m,n))*dxh(m,n)* hhv(m,n)
            gry = RHSy(m,n) + sly  + RHSy_dif(m,n) + RHSy_adv(m,n)      &
                - (rdis(m,n)+rdis(m,n+1))/2.0d0 *vbrtrp(m,n)*dxh(m,n)*dyt(m,n)*hhv(m,n)        &
                - ( rlh_s(m  ,n)*hhh(m  ,n)*dxb(m  ,n)*dyb(m  ,n)*(ubrtr(m  ,n+1)+ubrtr(m  ,n))             &
                +   rlh_s(m-1,n)*hhh(m-1,n)*dxb(m-1,n)*dyb(m-1,n)*(ubrtr(m-1,n+1)+ubrtr(m-1,n)) )/4.0d0

            vbrtrn(m,n) = (vbrtrp(m,n)*bp0 + gry )/(bp)
        endif
    enddo
  enddo
  _OMP_KERNEL_PARALLEL_END_

end subroutine

subroutine sw_next_step(nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2,  &
                        time_smooth,  &
                        lu, lcu, lcv,  &
                        ssh, sshn, sshp,  &
                        ubrtr, ubrtrn, ubrtrp,  &
                        vbrtr, vbrtrn, vbrtrp)

  integer, intent(in) :: nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2

  real(wp8), intent(in) :: time_smooth

  real(wp4), intent(in) :: lu(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                           lcu(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  &
                           lcv(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

  real(wp8), intent(inout) :: ssh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
                              sshn(bnd_x1:bnd_x2,bnd_y1:bnd_y2),  &
                              sshp(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

  real(wp8), intent(inout) :: ubrtr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
                              ubrtrn(bnd_x1:bnd_x2,bnd_y1:bnd_y2),  &
                              ubrtrp(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

  real(wp8), intent(inout) :: vbrtr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
                              vbrtrn(bnd_x1:bnd_x2,bnd_y1:bnd_y2),  &
                              vbrtrp(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

  integer :: m, n

  _OMP_KERNEL_PARALLEL_BEGIN_
  do n=ny_start-1,ny_end+1
    do m=nx_start-1,nx_end+1

        if(lu(m,n)>0.5) then
            sshp(m,n) = ssh(m,n)+time_smooth*(sshn(m,n)-2.0d0*ssh(m,n)+sshp(m,n))/2.0d0
            ssh(m,n) = sshn(m,n)
        endif
        if(lcu(m,n)>0.5) then
            ubrtrp(m,n) = ubrtr(m,n) + time_smooth*(ubrtrn(m,n)-2.0d0*ubrtr(m,n)+ubrtrp(m,n))/2.0d0
            ubrtr(m,n) = ubrtrn(m,n)
        endif
        if(lcv(m,n)>0.5) then
            vbrtrp(m,n) = vbrtr(m,n) + time_smooth*(vbrtrn(m,n)-2.0d0*vbrtr(m,n)+vbrtrp(m,n))/2.0d0
            vbrtr(m,n) = vbrtrn(m,n)
        endif

    enddo
  enddo
  _OMP_KERNEL_PARALLEL_END_

end subroutine

subroutine uv_trans_vort_kernel(nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2,  &
                                luu,                 &
                                dxt, dyt, dxb, dyb,  &
                                u, v, vort,          &
                                nlev)

 integer, intent(in) :: nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2

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

 _OMP_KERNEL_PARALLEL_BEGIN_
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
 _OMP_KERNEL_PARALLEL_END_

end subroutine uv_trans_vort_kernel

subroutine uv_trans_kernel(nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2,  &
                           lcu, lcv, luu,   &
                           dxh, dyh,        &
                           u, v, vort,      &
                           hq, hu, hv, hh,  &
                           RHSx, RHSy, nlev)

 integer, intent(in) :: nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2

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

 _OMP_KERNEL_PARALLEL_BEGIN_
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
  _OMP_KERNEL_PARALLEL_END_

endsubroutine uv_trans_kernel

subroutine uv_diff2_kernel(nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2,  &
                           lcu, lcv,                              &
                           dx, dy, dxt, dyt, dxh, dyh, dxb, dyb,  &
                           mu, str_t, str_s,                      &
                           hq, hu, hv, hh,                        &
                           RHSx, RHSy, nlev)

 integer, intent(in) :: nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2

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

  _OMP_KERNEL_PARALLEL_BEGIN_
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
  _OMP_KERNEL_PARALLEL_END_

endsubroutine uv_diff2_kernel

endmodule velssh_sw_module
