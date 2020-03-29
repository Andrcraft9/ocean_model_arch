module tracer_module
    ! Module description

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use kernel_interface_module, only: nxs, nxe, nys, nye

    implicit none
    save
    private

    public :: tracer_transport

contains

subroutine tracer_predictor_fluxes_kernel(lu, lcu, lcv, dxh, dyh, ff, uh, vh, ww, flux_x_t, flux_y_t, flux_z_t)
    ! Predictor step, compute fluxes
    !
    ! INPUT/OUTPUT:
    ! Grid data:
    !
    ! Ocean data:
    real(4), intent(in) :: ff(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),      &  ! Transported tracer
                           uh(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),      &  ! Transporting velocities
                           vh(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),      &
                           ww(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1)
    real(4), intent(out) :: flux_x_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &
                            flux_y_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &
                            flux_z_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1)
    ! LOCAL:
    integer :: m, n, k 

    do n=ny_start,ny_end
        do m=nx_start,nx_end    

            if(lcu(m,n)>0.5) then
                do k=1,nz      
                    flux_x_t(m,n,k)= uh(m,n,k)*dyh(m,n)*(ff(m+1,n,k)-ff(m  ,n,k))/2.0
                enddo     
            endif

            if(lcv(m,n)>0.5) then
                do k=1,nz      
                    flux_y_t(m,n,k)= vh(m,n,k)*dxh(m,n)*(ff(m,n+1,k)-ff(m,n  ,k))/2.0
                enddo     
            endif

            if(lu(m,n)>0.5) then
                do k=2,nz
                    flux_z_t(m,n,k)= ww(m,n,k) * (ff(m,n,k)-ff(m,n,k-1))/2.0
                enddo
            endif

        enddo
    enddo

    call syncborder_real(flux_x_t, nz)
    call syncborder_real(flux_y_t, nz)
    if(periodicity_x/=0) then
      call cyclize_x(flux_x_t,nx,ny,nz,mmm,mm)
    end if
    if(periodicity_y/=0) then        
      call cyclize_y(flux_y_t,nx,ny,nz,nnn,nn)            
    end if
    
end subroutine

subroutine tracer_predictor_kernel(lu, sqt, dz, hhts_n, ff, flux_x_t, flux_y_t, flux_z_t, fm)
    ! Predictor, final step
    !
    ! INPUT/OUTPUT:
    ! Grid data:
    !
    ! Ocean data:
    real(4), intent(in) :: ff(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz)  ! Transported tracer
    real(4), intent(in) :: flux_x_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &
                           flux_y_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &
                           flux_z_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1)
    real(4), intent(out) :: fm(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz)
    ! LOCAL:
    integer :: m, n, k
    real(4) :: rhs

    do n=ny_start,ny_end
        do m=nx_start,nx_end    
            if(lu(m,n)>0.5) then
                do k=1,nz
                    rhs=(flux_x_t(m,n,k)+flux_x_t(m-1,n,k)              &
                        +flux_y_t(m,n,k)+flux_y_t(m,n-1,k)  )/sqt(m,n)  &
                        +(flux_z_t(m,n,k+1)+flux_z_t(m,n,k))/dz(k)

                    fm(m,n,k)=ff(m,n,k) - pred_corr_impl_ts*tau *rhs/hhts_n(m,n)
                enddo
            endif
        enddo
    enddo

    call syncborder_real(fm, nz)
    if(periodicity_x/=0) then
      call cyclize_x(fm,nx,ny,nz,mmm,mm)
    end if
    if(periodicity_y/=0) then
      call cyclize_y(fm,nx,ny,nz,nnn,nn)            
    end if
    
end subroutine

subroutine tracer_corrector_fluxes_kernel(lu, lcu, lcv, dxh, dyh, dz, uh, vh, ww, flux_top, flux_x_t, flux_y_t, flux_z_t, fm)
    ! Corrector step, compute fluxes
    !
    ! INPUT/OUTPUT:
    ! Grid data:
    !
    ! Ocean data:
    real(4), intent(in) :: uh(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),      &  ! Transporting velocities
                           vh(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),      &
                           ww(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1)
    real(4), intent(in) :: flux_top(bnd_x1:bnd_x2,bnd_y1:bnd_y2)
    real(4), intent(inout) :: flux_x_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &
                              flux_y_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &
                              flux_z_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1)
    real(4), intent(in) :: fm(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz)
    ! LOCAL:
    integer :: m, n, k

    do n=ny_start,ny_end
        do m=nx_start,nx_end    
            if(lcu(m,n)>0.5) then
                do k=1,nz
                    flux_x_t(m,n,k)= uh(m,n,k)*dyh(m,n)*(fm(m+1,n,k)+fm(m,n,k))/2.0
                    flux_f_x(m,n)=flux_f_x(m,n) + flux_x_t(m,n,k)*dz(k)
                enddo
            endif

            if(lcv(m,n)>0.5) then
                do k=1,nz
                    flux_y_t(m,n,k)= vh(m,n,k)*dxh(m,n)*(fm(m,n+1,k)+fm(m,n,k))/2.0
                    flux_f_y(m,n)=flux_f_y(m,n) + flux_y_t(m,n,k)*dz(k)       
                enddo
            endif  

            if(lu(m,n)>0.5) then
                flux_z_t(m,n,1)= flux_top(m,n)
                do k=2,nz
                    flux_z_t(m,n,k)= ww(m,n,k) * (fm(m,n,k)+fm(m,n,k-1))/2.0
                enddo
            endif
        enddo
    enddo

    call syncborder_real(flux_x_t, nz)
    call syncborder_real(flux_y_t, nz)
    if(periodicity_x/=0) then
      call cyclize_x(flux_x_t,nx,ny,nz  ,mmm,mm)
    end if
    if(periodicity_y/=0) then          
      call cyclize_y(flux_y_t,nx,ny,nz  ,nnn,nn)            
    end if
    call syncborder_real(flux_f_x, 1)
    call syncborder_real(flux_f_y, 1)
    if(periodicity_x/=0) then
      call cyclize_x(flux_f_x,nx,ny, 1,mmm,mm)    
    end if
    if(periodicity_y/=0) then
      call cyclize_y(flux_f_y,nx,ny, 1,nnn,nn)
    end if
end subroutine

subroutine tracer_corrector_kernel(lu, sqt, dz, hhts, hhts_n, ff, flux_x_t, flux_y_t, flux_z_t)
    ! Corrector, final step
    !
    ! INPUT/OUTPUT:
    ! Grid data:
    !
    ! Ocean data:
    real(4), intent(inout) :: ff(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz)  ! Transported tracer
    real(4), intent(in) :: flux_x_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &
                           flux_y_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &
                           flux_z_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1)
    ! LOCAL:
    integer :: m, n, k
    real(4) :: rhs

    do n=ny_start,ny_end
        do m=nx_start,nx_end    
            if(lu(m,n)>0.5) then
                do k=1,nz
                    rhs=(flux_x_t(m,n,k)-flux_x_t(m-1,n,k)                      &
                        +flux_y_t(m,n,k)-flux_y_t(m,n-1,k)  ) /sqt(m,n)         &
                        +(flux_z_t(m,n,k+1)-flux_z_t(m,n,k))/dz(k)

                    ff(m,n,k)=(ff(m,n,k)*hhts(m,n) - rhs*tau)/hhts_n(m,n)
                enddo
            endif
        enddo
    enddo

    call syncborder_real(ff, nz)
    if(periodicity_x/=0) then
      call cyclize_x(ff,nx,ny,nz,mmm,mm)
    end if
    if(periodicity_y/=0) then
      call cyclize_y(ff,nx,ny,nz,nnn,nn)            
    end if
end subroutine

!==========================================================================================================
subroutine tracer_transport(ff,     &
                           tau,     &
                            uh,     &
                            vh,     &
                            ww,     &
                        flux_top,   &
                        flux_f_x,   &
                        flux_f_y,   &
                           ksw_lbc, &
                           numlbc,  &
                           lqpx,    &
                           lqpy,    &
                           ind_lbc, &
                           flbc)

use basin_grid
use ocalg_routes

include 'transport.fi'
 
 real(4) tau, src_factor, factor_mu, factor_nu, rhs

 real(4) ff(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),      &   !Transported tracer
         uh(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),      &
         vh(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),      &
         ww(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1)         !Transporting velocities
 
 real(4) flux_top(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

 real(4) flux_f_x(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &   !total flux on x-direction(1-total, 2-advection, 3-diffusion, 4-GM)
         flux_f_y(bnd_x1:bnd_x2,bnd_y1:bnd_y2)           !total flux on y-direction(1-total, 2-advection, 3-diffusion, 4-GM)

 integer numlbc, ksw_lbc, lqp
 integer lqpx(numlbc), lqpy(numlbc)
 real(4) flbc(numlbc,nz)
 character ind_lbc(numlbc)

integer m,n,k

real(4), allocatable::  flux_x_t(:,:,:), flux_y_t(:,:,:),      &
                        flux_z_t(:,:,:), fm(:,:,:)

    allocate(flux_x_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &
             flux_y_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &
             flux_z_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),   &
                  fm(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz) ) 

   flux_x_t=0.0
   flux_y_t=0.0
   flux_z_t=0.0
         fm=0.0

   flux_f_x=0.0
   flux_f_y=0.0


   deallocate(fm,flux_z_t,flux_y_t,flux_x_t)

endsubroutine tracer_transport

end module tracer_module