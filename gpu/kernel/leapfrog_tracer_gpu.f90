#include "macros/mpp_macros.fi"
#ifdef _GPU_MODE_

module tracer_gpu_module

     use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
     use cudafor
   
     implicit none
     save
     public

     real(wp8), constant :: time_smooth
   
contains

subroutine tracer_load_constant_mem()
    use config_sw_module, only: time_smooth_host => time_smooth
    time_smooth = time_smooth_host
  end subroutine

attributes(global) subroutine tran_diff_fluxes_kernel_gpu(nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2,  &
                                                          lcu, lcv,            &
                                                          dxt, dyt, dxh, dyh,  &
                                                          hhu, hhv,            &
                                                          ff,                  &
                                                          ffp,                 &
                                                          uu,                  &
                                                          vv,                  &
                                                          mu,                  &
                                                          factor_mu,           &
                                                          flux_x,              &
                                                          flux_y)

 integer, intent(in), value :: nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2
 real(wp4), intent(in) :: lcu(bnd_x1:bnd_x2, bnd_y1:bnd_y2),   &
                          lcv(bnd_x1:bnd_x2, bnd_y1:bnd_y2)
 real(wp4), intent(in) :: dxt(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          dyt(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          dxh(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          dyh(bnd_x1:bnd_x2, bnd_y1:bnd_y2)
 real(wp8), intent(in) :: hhu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),  &
                          hhv(bnd_x1:bnd_x2,bnd_y1:bnd_y2)
 real(wp8), intent(in) :: ff(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &  !Transported tracer
                          ffp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),  &  !Transported tracer (previous step)
                          uu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
                          vv(bnd_x1:bnd_x2,bnd_y1:bnd_y2)       !Transporting velocities
 real(wp8), intent(in) :: mu(bnd_x1:bnd_x2,bnd_y1:bnd_y2)
 real(wp8), intent(in), value :: factor_mu
 real(wp8), intent(inout) :: flux_x(bnd_x1:bnd_x2,bnd_y1:bnd_y2),  &  !total flux on x-direction
                             flux_y(bnd_x1:bnd_x2,bnd_y1:bnd_y2)      !total flux on y-direction

 integer, value :: m, n
 real(wp8), value :: dfdx, dfdy
 real(wp8), value :: flux_1d, mu_1d
 real(wp8), value :: flux_adv, flux_diff, flux_gm

   !flux_x=0.0d0
   !flux_y=0.0d0
   !flux_f_x=0.0d0
   !flux_f_y=0.0d0

   ! Compute flux_x, flux_y, flux_f_x, flux_f_y
   !do n = ny_start, ny_end
   ! do m = nx_start, nx_end
    m = (blockIdx%x-1)*blockDim%x + threadIdx%x + (nx_start) - 1
    n = (blockIdx%y-1)*blockDim%y + threadIdx%y + (ny_start) - 1
    if (m <= nx_end .and. n <= ny_end) then
     if (lcu(m,n)>0.5) then

       flux_1d=0.0d0
       mu_1d=0.0d0

       dfdx = ff(m+1,n) - ff(m,n) ! Try ff instead of ffp
       mu_1d = (mu(m,n) + mu(m+1,n))/2.0d0*factor_mu*dyh(m,n)/dxt(m,n)
       flux_1d = mu_1d*hhu(m,n)*dfdx

       flux_diff =  flux_1d
       flux_adv  = -uu(m,n)*hhu(m,n)*dyh(m,n)*(ff(m,n) + ff(m+1,n))/2.0d0
       flux_gm = 0.0d0
       flux_x(m, n) = flux_adv + flux_diff + flux_gm
       !flux_f_x(m,n,1) = flux_f_x(m,n,1) + (flux_adv + flux_diff + flux_gm)
       !flux_f_x(m,n,2) = flux_f_x(m,n,2) + flux_adv
       !flux_f_x(m,n,3) = flux_f_x(m,n,3) + flux_diff
       !flux_f_x(m,n,4) = flux_f_x(m,n,4) + flux_gm
     endif

     if (lcv(m,n)>0.5) then

       flux_1d=0.0d0
       mu_1d=0.0d0

       dfdy = ff(m,n+1) - ff(m,n) ! Try ff instead of ffp
       mu_1d = (mu(m,n) + mu(m,n+1))/2.0d0*factor_mu*dxh(m,n)/dyt(m,n)
       flux_1d = mu_1d*hhv(m,n)*dfdy

       flux_diff = flux_1d
       flux_adv  = -vv(m,n)*hhv(m,n)*dxh(m,n)*(ff(m,n) + ff(m,n+1))/2.0d0
       flux_gm   = 0.0d0
       flux_y(m, n) = flux_adv + flux_diff + flux_gm
       !flux_f_y(m,n,1) = flux_f_y(m,n,1) + (flux_adv + flux_diff + flux_gm)
       !flux_f_y(m,n,2) = flux_f_y(m,n,2) + flux_adv
       !flux_f_y(m,n,3) = flux_f_y(m,n,3) + flux_diff
       !flux_f_y(m,n,4) = flux_f_y(m,n,4) + flux_gm
     endif
    endif
endsubroutine

attributes(global) subroutine tran_diff_tracer_kernel_gpu(nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2,  &
                                                          lu,                  &
                                                          dx, dy,              &
                                                          tau,                 &
                                                          hhqn, hhqp,          &
                                                          flux_x,              &
                                                          flux_y,              &
                                                          ffp,                 &
                                                          ffn)

 integer, intent(in), value :: nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2
 real(wp4), intent(in) :: lu(bnd_x1:bnd_x2, bnd_y1:bnd_y2)
 real(wp4), intent(in) :: dx(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                          dy(bnd_x1:bnd_x2, bnd_y1:bnd_y2)
 real(wp8), intent(in), value :: tau
 real(wp8), intent(in) :: hhqn(bnd_x1:bnd_x2,bnd_y1:bnd_y2),  &
                          hhqp(bnd_x1:bnd_x2,bnd_y1:bnd_y2)
 real(wp8), intent(in) :: flux_x(bnd_x1:bnd_x2,bnd_y1:bnd_y2),  &  !total flux on x-direction
                          flux_y(bnd_x1:bnd_x2,bnd_y1:bnd_y2)      !total flux on y-direction
 real(wp8), intent(in) :: ffp(bnd_x1:bnd_x2,bnd_y1:bnd_y2)  !Transported tracer (previous step)
 real(wp8), intent(inout) :: ffn(bnd_x1:bnd_x2,bnd_y1:bnd_y2)
 
integer, value :: m, n
real(wp8), value :: eta
real(wp8), value :: bp, bp0, rhs

      ! Compute fn
      !do n = ny_start, ny_end
      ! do m = nx_start, nx_end
      m = (blockIdx%x-1)*blockDim%x + threadIdx%x + (nx_start) - 1
      n = (blockIdx%y-1)*blockDim%y + threadIdx%y + (ny_start) - 1
      if (m <= nx_end .and. n <= ny_end) then
        if (lu(m,n)>0.5) then

         bp =  hhqn(m,n)*dx(m,n)*dy(m,n) / tau/2.0d0
         bp0=  hhqp(m,n)*dx(m,n)*dy(m,n) / tau/2.0d0

         rhs =  flux_x(m,n) - flux_x(m-1,n) + flux_y(m,n) - flux_y(m,n-1)
         eta = bp0*ffp(m,n) + rhs
         ffn(m, n) = eta / bp
        endif
      endif
endsubroutine

attributes(global) subroutine tracer_next_step_kernel_gpu(nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2,  &
                                                          time_smooth,  &
                                                          lu,   &
                                                          ffn,  &
                                                          ffp,  &
                                                          ff)

 integer, intent(in), value :: nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2
 real(wp8), intent(in), value :: time_smooth
 real(wp4), intent(in) :: lu(bnd_x1:bnd_x2, bnd_y1:bnd_y2)
 real(wp8), intent(in) :: ffn(bnd_x1:bnd_x2,bnd_y1:bnd_y2)
 real(wp8), intent(inout) :: ffp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),  &
                             ff(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

 integer, value :: m, n

      ! Update step: compute ffp, ff
      !do n = ny_start-1, ny_end+1
      ! do m = nx_start-1, nx_end+1
      m = (blockIdx%x-1)*blockDim%x + threadIdx%x + (nx_start - 1) - 1
      n = (blockIdx%y-1)*blockDim%y + threadIdx%y + (ny_start - 1) - 1
      if (m <= nx_end + 1 .and. n <= ny_end + 1) then
        if (lu(m,n)>0.5) then
          !ffp(m,n,k)=hhq(m,n)*ff(m,n,k)+time_smooth*(hhqn(m,n)*fn(m,n,k)-2.0d0*hhq(m,n)*ff(m,n,k)+hhqp(m,n)*ffp(m,n,k))/2.0d0
          ffp(m, n) = ff(m, n) + time_smooth*(ffn(m, n) - 2.0d0*ff(m, n) + ffp(m, n))/2.0d0
          ff(m, n) = ffn(m, n)
        endif
      endif
endsubroutine

end module tracer_gpu_module

#endif