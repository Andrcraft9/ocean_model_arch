module velocity_module

contains

!======================================================================
subroutine div_velocity(ub,vb,hu,hv,div_btr)
use basin_grid
use ocalg_routes

real(8)   ub(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
          vb(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
     div_btr(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

real(4) tau
integer m, n, k

!$omp parallel do private(m,n,k)
      do n=ny_start,ny_end
       do m=nx_start,nx_end
        if (lu(m,n)>0.5) then
          div_btr(m,n) = (  ub(m  ,n  )*hu(m  ,n  )*dyh(m  ,n  )                 &
                          - ub(m-1,n  )*hu(m-1,n  )*dyh(m-1,n  )                 &
                          + vb(m  ,n  )*hv(m  ,n  )*dxh(m  ,n  )                 &
                          - vb(m  ,n-1)*hv(m  ,n-1)*dxh(m  ,n-1)  )/sqt(m,n)     

        endif
       enddo
      enddo
!$omp end parallel do
  
  call syncborder_real(div_bcl, nz)
  call syncborder_real8(div_btr, 1)

  if(periodicity_x/=0) then
    call  cyclize_x(div_bcl,nx,ny,nz,mmm,mm)
    call cyclize8_x(div_btr,nx,ny,1,mmm,mm)
  endif

  if(periodicity_y/=0) then
    call  cyclize_y(div_bcl,nx,ny,nz,nnn,nn)
    call cyclize8_y(div_btr,nx,ny,1,nnn,nn)
  endif
endsubroutine div_velocity

end module velocity