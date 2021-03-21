#include "macros/mpp_macros.fi"
#ifdef _GPU_MODE_

module depth_gpu_module

  use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
  use cudafor

  implicit none
  save
  public

  integer, constant :: full_free_surface
  real(wp8), constant :: time_smooth

contains

subroutine depth_load_constant_mem()
  use config_sw_module, only: full_free_surface_host => full_free_surface, time_smooth_host => time_smooth

  full_free_surface = full_free_surface_host
  time_smooth = time_smooth_host
end subroutine

attributes(global) subroutine hh_init_kernel_gpu(nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2,  &
                                                 lu, llu, llv, luh,  &
                                                 dx, dy, dxt, dyt, dxh, dyh, dxb, dyb,  &
                                                 hq, hqp, hqn,    &
                                                 hu, hup, hun,    &
                                                 hv, hvp, hvn,    &
                                                 hh, hhp, hhn,    &
                                                 sh, shp, h_r)

 integer, intent(in), value :: nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2

 real(wp4), intent(in) :: lu(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                          llu(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          llv(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          luh(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

 real(wp4), intent(in) :: dx(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                          dy(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                          dxt(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          dyt(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          dxh(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          dyh(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          dxb(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          dyb(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

 real(wp8), intent(inout) :: hq(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hqp(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hqn(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  &
                             hu(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hup(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hun(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  &
                             hv(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hvp(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hvn(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  &
                             hh(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hhp(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hhn(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  &
                             sh(bnd_x1:bnd_x2, bnd_y1:bnd_y2), shp(bnd_x1:bnd_x2, bnd_y1:bnd_y2), h_r(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

 real(wp8), value :: slu
 integer, value :: m, n
 
 real(wp8), value ::  hq_mn, hq_mpn, hq_mnp, hq_mpnp
 real(wp8), value ::  hqp_mn, hqp_mpn, hqp_mnp, hqp_mpnp
 real(wp8), value ::  hqn_mn, hqn_mpn, hqn_mnp, hqn_mpnp

      ! do n=ny_start-1,ny_end
      ! do m=nx_start-1,nx_end
      m = (blockIdx%x-1)*blockDim%x + threadIdx%x + (nx_start - 1) - 1
      n = (blockIdx%y-1)*blockDim%y + threadIdx%y + (ny_start - 1) - 1
      if (m <= nx_end + 1 .and. n <= ny_end + 1) then
        ! Save values to avoid threadsync
        ! hqn=h_r
        hqn_mn   = h_r(m,  n  )
        hqn_mpn  = h_r(m+1,n  )
        hqn_mnp  = h_r(m,  n+1)
        hqn_mpnp = h_r(m+1,n+1)
        ! hqp=h_r + shp*dfloat(full_free_surface)
        hq_mn   = hqn_mn   + sh(m,  n  ) * dfloat(full_free_surface)
        hq_mpn  = hqn_mpn  + sh(m+1,n  ) * dfloat(full_free_surface)
        hq_mnp  = hqn_mnp  + sh(m,  n+1) * dfloat(full_free_surface)
        hq_mpnp = hqn_mpnp + sh(m+1,n+1) * dfloat(full_free_surface)
        ! hq =h_r + sh *dfloat(full_free_surface)
        hqp_mn   = hqn_mn   + shp(m,  n  ) * dfloat(full_free_surface)
        hqp_mpn  = hqn_mpn  + shp(m+1,n  ) * dfloat(full_free_surface)
        hqp_mnp  = hqn_mnp  + shp(m,  n+1) * dfloat(full_free_surface)
        hqp_mpnp = hqn_mpnp + shp(m+1,n+1) * dfloat(full_free_surface)

        ! Update for full area
        hq (m, n) = hq_mn
        hqp(m, n) = hqp_mn
        hqn(m, n) = hqn_mn

        if (m <= nx_end .and. n <= ny_end) then
          if(llu(m,n)>0.5) then
          ! interpolating hhq given on T-grid(lu) to hhu given on u-grid(lcu).
                slu=dble(lu(m,n)+lu(m+1,n))
                hu(m,n)=( hq_mn *dx(m  ,n)*dy(m  ,n)*dble(lu(m  ,n))   &
                        + hq_mpn*dx(m+1,n)*dy(m+1,n)*dble(lu(m+1,n)) )/slu/dxt(m,n)/dyh(m,n)
                hup(m,n)=(hqp_mn *dx(m  ,n)*dy(m  ,n)*dble(lu(m  ,n))   &
                        + hqp_mpn*dx(m+1,n)*dy(m+1,n)*dble(lu(m+1,n)) )/slu/dxt(m,n)/dyh(m,n)
                hun(m,n)=(hqn_mn *dx(m  ,n)*dy(m  ,n)*dble(lu(m  ,n))   &
                        + hqn_mpn*dx(m+1,n)*dy(m+1,n)*dble(lu(m+1,n)) )/slu/dxt(m,n)/dyh(m,n)
          endif

          if(llv(m,n)>0.5) then
          ! interpolating hhq given on T-grid(lu) to hhv given on v-grid(lcv).
                slu=dble(lu(m,n)+lu(m,n+1))
                hv(m,n)=( hq_mn *dx(m,n  )*dy(m,n  )*dble(lu(m,n  ))       &
                        + hq_mnp*dx(m,n+1)*dy(m,n+1)*dble(lu(m,n+1)) )/slu/dxh(m,n)/dyt(m,n)
                hvp(m,n)=( hqp_mn *dx(m,n  )*dy(m,n  )*dble(lu(m,n  ))       &
                         + hqp_mnp*dx(m,n+1)*dy(m,n+1)*dble(lu(m,n+1)) )/slu/dxh(m,n)/dyt(m,n)
                hvn(m,n)=( hqn_mn *dx(m,n  )*dy(m,n  )*dble(lu(m,n  ))       &
                         + hqn_mnp*dx(m,n+1)*dy(m,n+1)*dble(lu(m,n+1)) )/slu/dxh(m,n)/dyt(m,n)
          endif

          if(luh(m,n)>0.5) then
          ! interpolating hhq given on T-grid(lu) to hhh given on h-grid(luu).
                slu=dble(lu(m,n)+lu(m+1,n)+lu(m,n+1)+lu(m+1,n+1))
                hh(m,n)=( hq_mn *dx(m  ,n  )*dy(m  ,n  )*dble(lu(m  ,n  ))       &
                        + hq_mpn *dx(m+1,n  )*dy(m+1,n  )*dble(lu(m+1,n  ))       &
                        + hq_mnp *dx(m  ,n+1)*dy(m  ,n+1)*dble(lu(m  ,n+1))       &
                        + hq_mpnp*dx(m+1,n+1)*dy(m+1,n+1)*dble(lu(m+1,n+1)) )/slu/dxb(m,n)/dyb(m,n)
                hhp(m,n)=( hqp_mn  *dx(m  ,n  )*dy(m  ,n  )*dble(lu(m  ,n  ))       &
                        +  hqp_mpn *dx(m+1,n  )*dy(m+1,n  )*dble(lu(m+1,n  ))       &
                        +  hqp_mnp *dx(m  ,n+1)*dy(m  ,n+1)*dble(lu(m  ,n+1))       &
                        +  hqp_mpnp*dx(m+1,n+1)*dy(m+1,n+1)*dble(lu(m+1,n+1)) )/slu/dxb(m,n)/dyb(m,n)
                hhn(m,n)=( hqn_mn  *dx(m  ,n  )*dy(m  ,n  )*dble(lu(m  ,n  ))       &
                         + hqn_mpn *dx(m+1,n  )*dy(m+1,n  )*dble(lu(m+1,n  ))       &
                         + hqn_mnp *dx(m  ,n+1)*dy(m  ,n+1)*dble(lu(m  ,n+1))       &
                         + hqn_mpnp*dx(m+1,n+1)*dy(m+1,n+1)*dble(lu(m+1,n+1)) )/slu/dxb(m,n)/dyb(m,n)
          endif
        endif
      endif
endsubroutine hh_init_kernel_gpu

attributes(global) subroutine hh_update_kernel_gpu(nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2,  &
                                                   lu, llu, llv, luh,  &
                                                   dx, dy, dxt, dyt, dxh, dyh, dxb, dyb,  &
                                                   hqn, hun, hvn, hhn, sh, h_r)

 integer, intent(in), value :: nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2

 real(wp4), intent(in) :: lu(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                          llu(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          llv(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          luh(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

 real(wp4), intent(in) :: dx(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                          dy(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                          dxt(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          dyt(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          dxh(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          dyh(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          dxb(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          dyb(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

 real(wp8), intent(inout) :: hqn(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hun(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  &
                             hvn(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hhn(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  &
                             sh(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  h_r(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

 integer, value :: m,n
 real(wp8), value :: slu

 real(wp8), value :: hqn_mn, hqn_mpn, hqn_mnp, hqn_mpnp

      ! do n=ny_start-1,ny_end
      ! do m=nx_start-1,nx_end
      m = (blockIdx%x-1)*blockDim%x + threadIdx%x + (nx_start - 1) - 1
      n = (blockIdx%y-1)*blockDim%y + threadIdx%y + (ny_start - 1) - 1
      if (m <= nx_end + 1 .and. n <= ny_end + 1) then
        ! Save values to avoid threadsync
        ! hqn =h_r + sh
        hqn_mn   = h_r(m,  n  ) + sh(m,  n  )
        hqn_mpn  = h_r(m+1,n  ) + sh(m+1,n  )
        hqn_mnp  = h_r(m,  n+1) + sh(m,  n+1)
        hqn_mpnp = h_r(m+1,n+1) + sh(m+1,n+1)

        ! Update for full area
        hqn(m, n) = hqn_mn

        if (m <= nx_end .and. n <= ny_end) then
          if(llu(m,n)>0.5) then
            ! interpolating hhq given on T-grid(lu) to hhu given on u-grid(lcu).
            slu=dble(lu(m,n)+lu(m+1,n))
            hun(m,n)=( hqn_mn *dx(m  ,n)*dy(m  ,n)*dble(lu(m  ,n))   &
                     + hqn_mpn*dx(m+1,n)*dy(m+1,n)*dble(lu(m+1,n)) )/slu/dxt(m,n)/dyh(m,n)
          endif
  
          if(llv(m,n)>0.5) then
            ! interpolating hhq given on T-grid(lu) to hhv given on v-grid(lcv).
            slu=dble(lu(m,n)+lu(m,n+1))
            hvn(m,n)=( hqn_mn *dx(m,n  )*dy(m,n  )*dble(lu(m,n  ))       &
                     + hqn_mnp*dx(m,n+1)*dy(m,n+1)*dble(lu(m,n+1)) )/slu/dxh(m,n)/dyt(m,n)
          endif
  
          if(luh(m,n)>0.5) then
            ! interpolating hhq given on T-grid(lu) to hhh given on h-grid(luu).
            slu=dble(lu(m,n)+lu(m+1,n)+lu(m,n+1)+lu(m+1,n+1))
            hhn(m,n)=( hqn_mn  *dx(m  ,n  )*dy(m  ,n  )*dble(lu(m  ,n  ))       &
                     + hqn_mpn *dx(m+1,n  )*dy(m+1,n  )*dble(lu(m+1,n  ))       &
                     + hqn_mnp *dx(m  ,n+1)*dy(m  ,n+1)*dble(lu(m  ,n+1))       &
                     + hqn_mpnp*dx(m+1,n+1)*dy(m+1,n+1)*dble(lu(m+1,n+1)) )/slu/dxb(m,n)/dyb(m,n)
          endif
        endif
      endif

endsubroutine hh_update_kernel_gpu

attributes(global) subroutine hh_shift_kernel_gpu(nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2,  &
                                                  lu, llu, llv, luh,  &
                                                  hq, hqp, hqn,   &
                                                  hu, hup, hun,   &
                                                  hv, hvp, hvn,   &
                                                  hh, hhp, hhn)

 integer, intent(in), value :: nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2

 real(wp4), intent(in) :: lu(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                          llu(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          llv(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          luh(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

 real(wp8), intent(inout) :: hq(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hqp(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hqn(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  &
                             hu(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hup(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hun(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  &
                             hv(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hvp(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hvn(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  &
                             hh(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hhp(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hhn(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

 integer, value :: m, n

      ! do n=ny_start-1,ny_end+1
      ! do m=nx_start-1,nx_end+1
      m = (blockIdx%x-1)*blockDim%x + threadIdx%x + (nx_start - 1) - 1
      n = (blockIdx%y-1)*blockDim%y + threadIdx%y + (ny_start - 1) - 1
      if (m <= nx_end + 1 .and. n <= ny_end + 1) then
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
      endif

endsubroutine hh_shift_kernel_gpu

endmodule depth_gpu_module

#endif