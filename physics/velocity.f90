module velocity_module

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use kernel_interface_module, only: nxs => bnd_x1, nxe => bnd_x2, nys => bnd_y1, nye => bnd_y2

    implicit none
    save
    private

    public :: div_velocity_kernel

contains

    subroutine div_velocity_kernel(lu, dxh, dyh, sqt, hu, hv, ub, vb, div_btr)
        ! Compute divergence of velocity vectors.
        !
        ! INPUT/OUTPUT:
        ! Grid data
        real(wp4), intent(in) :: lu(nxs:nxe, nys:nye),  &
                                 dxh(nxs:nxe, nys:nye), &
                                 dyh(nxs:nxe, nys:nye), &
                                 sqt(nxs:nxe, nys:nye)
        real(wp4), intent(in) :: hu(nxs:nxe, nys:nye),  &
                                 hv(nxs:nxe, nys:nye)
        ! Ocean data
        real(wp8), intent(in) :: ub(nxs:nxe, nys:nye),  &
                                 vb(nxs:nxe, nys:nye)
        real(wp8), intent(inout) :: div_btr(nxs:nxe, nys:nye)
        ! LOCAL:
        integer m, n

          !do n=nys,nye
          !  do m=nxs,nxe
          do n=nys+1,nye-1
              do m=nxs+1,nxe-1
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