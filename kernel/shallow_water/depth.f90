module depth_module

  use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
  use config_sw_module, only: full_free_surface, time_smooth

  implicit none
  save
  private

#include "macros/mpp_macros.fi"

  public :: hh_init_kernel, hh_update_kernel, hh_shift_kernel

contains

subroutine hh_init_kernel(nx_start, nx_end, ny_start, ny_end,  &
                          lu, llu, llv, luh,  &
                          dx, dy, dxt, dyt, dxh, dyh, dxb, dyb,  &
                          hq, hqp, hqn,    &
                          hu, hup, hun,    &
                          hv, hvp, hvn,    &
                          hh, hhp, hhn,    &
                          sh, shp, h_r)

 integer, intent(in) :: nx_start, nx_end, ny_start, ny_end

 real(wp4), intent(in) :: lu (:, :),  &
                          llu(:, :),  & 
                          llv(:, :),  & 
                          luh(:, :)

 real(wp4), intent(in) :: dx (:, :),  &
                          dy (:, :),  &
                          dxt(:, :),  & 
                          dyt(:, :),  & 
                          dxh(:, :),  & 
                          dyh(:, :),  & 
                          dxb(:, :),  & 
                          dyb(:, :)

 real(wp8), intent(inout) :: hq(:, :), hqp(:, :), hqn(:, :),  &
                             hu(:, :), hup(:, :), hun(:, :),  &
                             hv(:, :), hvp(:, :), hvn(:, :),  &
                             hh(:, :), hhp(:, :), hhn(:, :),  &
                             sh(:, :), shp(:, :), h_r(:, :)

 real(wp8) :: slu
 integer :: m, n

       hq =h_r + sh *dfloat(full_free_surface)
       hqp=h_r + shp*dfloat(full_free_surface)
       hqn=h_r

      !do n=ny_start-2, ny_end+1
      !do m=nx_start-2, nx_end+1
      _OMP_KERNEL_PARALLEL_BEGIN_
      do n=ny_start-1,ny_end
       do m=nx_start-1,nx_end

        if(llu(m,n)>0.5) then
! interpolating hhq given on T-grid(lu) to hhu given on u-grid(lcu).
          slu=dble(lu(m,n)+lu(m+1,n))
          hu(m,n)=( hq(m  ,n)*dx(m  ,n)*dy(m  ,n)*dble(lu(m  ,n))   &
                  + hq(m+1,n)*dx(m+1,n)*dy(m+1,n)*dble(lu(m+1,n)) )/slu/dxt(m,n)/dyh(m,n)
         hup(m,n)=( hqp(m  ,n)*dx(m  ,n)*dy(m  ,n)*dble(lu(m  ,n))   &
                  + hqp(m+1,n)*dx(m+1,n)*dy(m+1,n)*dble(lu(m+1,n)) )/slu/dxt(m,n)/dyh(m,n)
         hun(m,n)=( hqn(m  ,n)*dx(m  ,n)*dy(m  ,n)*dble(lu(m  ,n))   &
                  + hqn(m+1,n)*dx(m+1,n)*dy(m+1,n)*dble(lu(m+1,n)) )/slu/dxt(m,n)/dyh(m,n)
        endif

        if(llv(m,n)>0.5) then
! interpolating hhq given on T-grid(lu) to hhv given on v-grid(lcv).
          slu=dble(lu(m,n)+lu(m,n+1))
          hv(m,n)=( hq(m,n  )*dx(m,n  )*dy(m,n  )*dble(lu(m,n  ))       &
                  + hq(m,n+1)*dx(m,n+1)*dy(m,n+1)*dble(lu(m,n+1)) )/slu/dxh(m,n)/dyt(m,n)
         hvp(m,n)=( hqp(m,n  )*dx(m,n  )*dy(m,n  )*dble(lu(m,n  ))       &
                  + hqp(m,n+1)*dx(m,n+1)*dy(m,n+1)*dble(lu(m,n+1)) )/slu/dxh(m,n)/dyt(m,n)
         hvn(m,n)=( hqn(m,n  )*dx(m,n  )*dy(m,n  )*dble(lu(m,n  ))       &
                  + hqn(m,n+1)*dx(m,n+1)*dy(m,n+1)*dble(lu(m,n+1)) )/slu/dxh(m,n)/dyt(m,n)
        endif

        if(luh(m,n)>0.5) then
! interpolating hhq given on T-grid(lu) to hhh given on h-grid(luu).
          slu=dble(lu(m,n)+lu(m+1,n)+lu(m,n+1)+lu(m+1,n+1))
          hh(m,n)=( hq(m  ,n  )*dx(m  ,n  )*dy(m  ,n  )*dble(lu(m  ,n  ))       &
                  + hq(m+1,n  )*dx(m+1,n  )*dy(m+1,n  )*dble(lu(m+1,n  ))       &
                   +hq(m  ,n+1)*dx(m  ,n+1)*dy(m  ,n+1)*dble(lu(m  ,n+1))       &
                  + hq(m+1,n+1)*dx(m+1,n+1)*dy(m+1,n+1)*dble(lu(m+1,n+1)) )/slu/dxb(m,n)/dyb(m,n)
         hhp(m,n)=( hqp(m  ,n  )*dx(m  ,n  )*dy(m  ,n  )*dble(lu(m  ,n  ))       &
                  + hqp(m+1,n  )*dx(m+1,n  )*dy(m+1,n  )*dble(lu(m+1,n  ))       &
                   +hqp(m  ,n+1)*dx(m  ,n+1)*dy(m  ,n+1)*dble(lu(m  ,n+1))       &
                  + hqp(m+1,n+1)*dx(m+1,n+1)*dy(m+1,n+1)*dble(lu(m+1,n+1)) )/slu/dxb(m,n)/dyb(m,n)
         hhn(m,n)=( hqn(m  ,n  )*dx(m  ,n  )*dy(m  ,n  )*dble(lu(m  ,n  ))       &
                  + hqn(m+1,n  )*dx(m+1,n  )*dy(m+1,n  )*dble(lu(m+1,n  ))       &
                   +hqn(m  ,n+1)*dx(m  ,n+1)*dy(m  ,n+1)*dble(lu(m  ,n+1))       &
                  + hqn(m+1,n+1)*dx(m+1,n+1)*dy(m+1,n+1)*dble(lu(m+1,n+1)) )/slu/dxb(m,n)/dyb(m,n)
        endif

       end do
    end do
    _OMP_KERNEL_PARALLEL_END_

endsubroutine hh_init_kernel

subroutine hh_update_kernel(nx_start, nx_end, ny_start, ny_end,  &
                            lu, llu, llv, luh,  &
                            dx, dy, dxt, dyt, dxh, dyh, dxb, dyb,  &
                            hqn, hun, hvn, hhn, sh, h_r)

 integer, intent(in) :: nx_start, nx_end, ny_start, ny_end

 real(wp4), intent(in) :: lu (:, :),  &
                          llu(:, :),  & 
                          llv(:, :),  & 
                          luh(:, :)

 real(wp4), intent(in) :: dx (:, :),  &
                          dy (:, :),  &
                          dxt(:, :),  & 
                          dyt(:, :),  & 
                          dxh(:, :),  & 
                          dyh(:, :),  & 
                          dxb(:, :),  & 
                          dyb(:, :)

 real(wp8), intent(inout) :: hqn(:, :), hun(:, :),  &
                             hvn(:, :), hhn(:, :),  &
                             sh (:, :), h_r(:, :)

 integer :: m,n
 real(wp8) :: slu

      hqn =h_r + sh

      !do n=ny_start-2, ny_end+1
      !do m=nx_start-2, nx_end+1
      _OMP_KERNEL_PARALLEL_BEGIN_
      do n=ny_start-1,ny_end
       do m=nx_start-1,nx_end

        if(llu(m,n)>0.5) then
! interpolating hhq given on T-grid(lu) to hhu given on u-grid(lcu).
          slu=dble(lu(m,n)+lu(m+1,n))
          hun(m,n)=( hqn(m  ,n)*dx(m  ,n)*dy(m  ,n)*dble(lu(m  ,n))   &
                   + hqn(m+1,n)*dx(m+1,n)*dy(m+1,n)*dble(lu(m+1,n)) )/slu/dxt(m,n)/dyh(m,n)
        endif

        if(llv(m,n)>0.5) then
! interpolating hhq given on T-grid(lu) to hhv given on v-grid(lcv).
          slu=dble(lu(m,n)+lu(m,n+1))
          hvn(m,n)=( hqn(m,n  )*dx(m,n  )*dy(m,n  )*dble(lu(m,n  ))       &
                   + hqn(m,n+1)*dx(m,n+1)*dy(m,n+1)*dble(lu(m,n+1)) )/slu/dxh(m,n)/dyt(m,n)
        endif

        if(luh(m,n)>0.5) then
! interpolating hhq given on T-grid(lu) to hhh given on h-grid(luu).
          slu=dble(lu(m,n)+lu(m+1,n)+lu(m,n+1)+lu(m+1,n+1))
          hhn(m,n)=( hqn(m  ,n  )*dx(m  ,n  )*dy(m  ,n  )*dble(lu(m  ,n  ))       &
                   + hqn(m+1,n  )*dx(m+1,n  )*dy(m+1,n  )*dble(lu(m+1,n  ))       &
                    +hqn(m  ,n+1)*dx(m  ,n+1)*dy(m  ,n+1)*dble(lu(m  ,n+1))       &
                   + hqn(m+1,n+1)*dx(m+1,n+1)*dy(m+1,n+1)*dble(lu(m+1,n+1)) )/slu/dxb(m,n)/dyb(m,n)
        endif

       end do
    end do
    _OMP_KERNEL_PARALLEL_END_

endsubroutine hh_update_kernel

subroutine hh_shift_kernel(nx_start, nx_end, ny_start, ny_end,  &
                           lu, llu, llv, luh,  &
                           hq, hqp, hqn,   &
                           hu, hup, hun,   &
                           hv, hvp, hvn,   &
                           hh, hhp, hhn)

 integer, intent(in) :: nx_start, nx_end, ny_start, ny_end

 real(wp4), intent(in) :: lu (:, :),  &
                          llu(:, :),  & 
                          llv(:, :),  & 
                          luh(:, :)

 real(wp8), intent(inout) :: hq(:, :), hqp(:, :), hqn(:, :),  &
                             hu(:, :), hup(:, :), hun(:, :),  &
                             hv(:, :), hvp(:, :), hvn(:, :),  &
                             hh(:, :), hhp(:, :), hhn(:, :)

 integer :: m, n

      _OMP_KERNEL_PARALLEL_BEGIN_
      do n=ny_start-1,ny_end+1
       do m=nx_start-1,nx_end+1

        if(llu(m,n)>0.5) then
          hup(m,n)= hu(m,n) + time_smooth*(hun(m,n)-2.0d0*hu(m,n)+hup(m,n))/2.0d0
           hu(m,n)= hun(m,n)
        endif

        if(llv(m,n)>0.5) then
          hvp(m,n)= hv(m,n) + time_smooth*(hvn(m,n)-2.0d0*hv(m,n)+hvp(m,n))/2.0d0
           hv(m,n)= hvn(m,n)
        endif

        if(lu(m,n)>0.5) then
          hqp(m,n)= hq(m,n) + time_smooth*(hqn(m,n)-2.0d0*hq(m,n)+hqp(m,n))/2.0d0
           hq(m,n)= hqn(m,n)
        endif

        if(luh(m,n)>0.5) then
          hhp(m,n)= hh(m,n) + time_smooth*(hhn(m,n)-2.0d0*hh(m,n)+hhp(m,n))/2.0d0
           hh(m,n)= hhn(m,n)
        endif

       end do
    end do
    _OMP_KERNEL_PARALLEL_END_

endsubroutine hh_shift_kernel

endmodule depth_module