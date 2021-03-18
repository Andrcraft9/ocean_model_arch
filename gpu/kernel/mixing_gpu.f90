#include "macros/mpp_macros.fi"
#ifdef _GPU_MODE_

module mixing_gpu_module

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
  
    implicit none
    save
    public

contains

!====================================================================================
attributes(global) subroutine stress_components_kernel(nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2,  &
                                                       lu, luu, dx, dy, dxt, dyt, dxh, dyh, dxb, dyb, u, v, str_t, str_s, nlev)

    integer, intent(in), value :: nx_start, nx_end, ny_start, ny_end, bnd_x1, bnd_x2, bnd_y1, bnd_y2
    integer, intent(in), value :: nlev
    
    real(wp4), intent(in) :: lu(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  &
                             luu(bnd_x1:bnd_x2, bnd_y1:bnd_y2)
                             
    real(wp4), intent(in) :: dx(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                             dy(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                             dxt(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                             dyt(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                             dxh(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                             dyh(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                             dxb(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  & 
                             dyb(bnd_x1:bnd_x2, bnd_y1:bnd_y2)
    real(wp8), intent(in) :: u(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev),    &
                              v(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev)
    real(wp8), intent(inout) :: str_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev),    &
                                str_s(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev)
     
    integer, value :: m, n, k
    
    ! do n=ny_start, ny_end
    ! do m=nx_start, nx_end
    m = (blockIdx%x-1)*blockDim%x + threadIdx%x + (nx_start) - 1
    n = (blockIdx%y-1)*blockDim%y + threadIdx%y + (ny_start) - 1
    if (m <= nx_end .and. n <= ny_end) then
        if(lu(m,n)>0.5) then
         do k=1,nlev 
          str_t(m,n,k)=dy(m,n)/dx(m,n)*(u(m,n,k)/dyh(m,n)-u(m-1,n,k)/dyh(m-1,n))     &
                      -dx(m,n)/dy(m,n)*(v(m,n,k)/dxh(m,n)-v(m,n-1,k)/dxh(m,n-1)) 
         enddo
        endif
    
        if(luu(m,n)>0.5) then
         do k=1,nlev
          str_s(m,n,k)=dxb(m,n)/dyb(m,n)*(u(m,n+1,k)/dxt(m,n+1)-u(m,n,k)/dxt(m,n))    &
                      +dyb(m,n)/dxb(m,n)*(v(m+1,n,k)/dyt(m+1,n)-v(m,n,k)/dyt(m,n))    
         enddo
        endif
    endif
endsubroutine

endmodule mixing_gpu_module

#endif