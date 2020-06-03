module mixing_module

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
  
    implicit none
    save
    private

#include "macros/mpp_macros.fi"

    public :: stress_components_kernel

contains

!====================================================================================
subroutine stress_components_kernel(nx_start, nx_end, ny_start, ny_end,  &
                                    lu, luu, dx, dy, dxt, dyt, dxh, dyh, dxb, dyb, u, v, str_t, str_s)

    integer, intent(in) :: nx_start, nx_end, ny_start, ny_end

    real(wp4), intent(in) :: lu (:, :),  &
                             luu(:, :)

    real(wp4), intent(in) :: dx (:, :),  & 
                             dy (:, :),  & 
                             dxt(:, :),  & 
                             dyt(:, :),  & 
                             dxh(:, :),  & 
                             dyh(:, :),  & 
                             dxb(:, :),  & 
                             dyb(:, :)
     real(wp8), intent(in) :: u(:, :),    &
                              v(:, :)
     real(wp8), intent(inout) :: str_t(:, :),    &
                                 str_s(:, :)
     
     integer :: m, n, k
    
     _OMP_KERNEL_PARALLEL_BEGIN_
     do n=ny_start, ny_end
       do m=nx_start, nx_end
    
        if(lu(m,n)>0.5) then
          str_t(m,n)  =dy(m,n)/dx(m,n)*(u(m,n)/dyh(m,n)-u(m-1,n)/dyh(m-1,n))     &
                      -dx(m,n)/dy(m,n)*(v(m,n)/dxh(m,n)-v(m,n-1)/dxh(m,n-1)) 
        endif
    
        if(luu(m,n)>0.5) then
          str_s(m,n)  =dxb(m,n)/dyb(m,n)*(u(m,n+1)/dxt(m,n+1)-u(m,n)/dxt(m,n))    &
                      +dyb(m,n)/dxb(m,n)*(v(m+1,n)/dyt(m+1,n)-v(m,n)/dyt(m,n))    
        endif
    
       enddo
     enddo
     _OMP_KERNEL_PARALLEL_END_
    
    endsubroutine

end