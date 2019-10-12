module velocity_module

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use kernel_interface_module, only: nx_start, nx_end, ny_start, ny_end

    implicit none
    save
    private

    public :: div_velocity_kernel

contains

    subroutine div_velocity_kernel(ub,vb,hu,hv,div_btr)
        real(wp4), intent(in) :: hu(x_start:nx_end, ny_start:ny_end),  &
                                 hv(x_start:nx_end, ny_start:ny_end)

        real(wp8), intent(in) :: ub(nx_start:nx_end, ny_start:ny_end),  &
                                 vb(nx_start:nx_end, ny_start:ny_end)
                
        real(wp8), intent(inout) :: div_btr(nx_start:nx_end, ny_start:ny_end)

        integer m, n

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
    endsubroutine

end module velocity_module