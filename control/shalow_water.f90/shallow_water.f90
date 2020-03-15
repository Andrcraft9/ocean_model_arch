module shallow_water
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid

    use key_switches

    use ocalg_routes
    use depth_routes
    use mixing_routes
    use flux_routes
    use velssh_routes
    implicit none

    contains

!===============================================================================
    ! explicit shallow water equation sloving
    subroutine expl_shallow_water( tau,     &
                                   ubrtr,   &
                                   ubrtrp,  &
                                   ubrtrn,  &
                                   vbrtr,   &
                                   vbrtrp,  &
                                   vbrtrn,  &
                                   ssh,     &
                                   sshp,    &
                                   sshn,    &
                                     RHSx,  &
                                     RHSy,  &
                                    wflux,  &
                                       mu,  &
                                      mu4,  &
                                     vort,  &
                                  str_t2d,  &
                                  str_s2d,  &
                                       fx,  &
                                       fy,  &
                                     rdis,  &
                                 RHSx_adv,  &
                                 RHSy_adv,  &
                                 RHSx_dif,  &
                                 RHSy_dif,  &
                                 RHSx_bfc,  &
                                 RHSy_bfc)

        implicit none

        real(8) tau
        integer m, n

        real(8)  ubrtr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                ubrtrp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                ubrtrn(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                 vbrtr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                vbrtrp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                vbrtrn(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                   ssh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                  sshp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                  sshn(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                   wflux(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                    RHSx(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                    RHSy(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                      mu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                     mu4(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                    vort(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                 str_t2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                 str_s2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                      fx(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                      fy(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                    rdis(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                RHSx_adv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                RHSy_adv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                RHSx_dif(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                RHSy_dif(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                RHSx_bfc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                RHSy_bfc(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

        real(8) bp, bp0, grx, gry, slx, sly, slxn, slyn

        integer ierr

        !computing ssh
        !$omp parallel do
        do n=ny_start,ny_end
            do m=nx_start,nx_end

                if(lu(m,n)>0.5) then
                    sshn(m,n) = sshp(m,n) + 2.0d0*tau*( wflux(m,n)/RefDen*dfloat(full_free_surface)   &
                    - ( ubrtr(m,n)*hhu(m,n)*dyh(m,n) - ubrtr(m-1,n)*hhu(m-1,n)*dyh(m-1,n)             &
                      + vbrtr(m,n)*hhv(m,n)*dxh(m,n) - vbrtr(m,n-1)*hhv(m,n-1)*dxh(m,n-1) )/(dx(m,n)*dy(m,n))  )
                endif

            enddo
        enddo
        !$omp end parallel do

        call syncborder_real8(sshn, 1)
        if(periodicity_x/=0) then
            call cyclize8_x(sshn,nx,ny,1,mmm,mm)
        endif
        if(periodicity_y/=0) then
            call cyclize8_y(sshn,nx,ny,1,nnn,nn)
        endif

        if (full_free_surface>0) then
            call hh_update(hhqn, hhun, hhvn, hhhn, sshn, hhq_rest)
        endif

        !computing advective and lateral-viscous terms for 2d-velocity
        if (trans_terms > 0) then
            call uv_trans(ubrtr, vbrtr, vort,     &
                          hhq, hhu, hhv, hhh,     &
                          RHSx_adv, RHSy_adv, 1)
        endif

        if (ksw_lat > 0) then
            call stress_components(ubrtrp, vbrtrp, str_t2d, str_s2d, 1)
            call uv_diff2(mu, str_t2d, str_s2d,   &
                          hhq, hhu, hhv, hhh,     &
                          RHSx_dif, RHSy_dif, 1)
            
            if(ksw_lat4 > 0) then
                call uv_diff4( mu4, str_t2d, str_s2d,  &
                               fx, fy, hhq, hhu, hhv, hhh,    &
                               RHSx_dif, RHSy_dif, 1 )
            endif
        endif

        !compute Bottom Friction (bfc)
        if (ksw_bfc > 0) then
            call uv_bfc(ubrtrp, vbrtrp, hhq, hhu, hhv, hhh, RHSx_bfc, RHSy_bfc, nbfc)
        endif

        !$omp parallel do private(bp, bp0, grx, gry, slx, sly, slxn, slyn)
        do n=ny_start,ny_end
            do m=nx_start,nx_end
                !zonal flux
                if(lcu(m,n)>0.5) then
                    bp  = hhun(m,n)*dxt(m,n)*dyh(m,n)/2.0d0/tau
                    bp0 = hhup(m,n)*dxt(m,n)*dyh(m,n)/2.0d0/tau

                    slx = - FreeFallAcc * ( ssh(m+1,n) - ssh(m,n))*dyh(m,n)* hhu(m,n)
                    !slxn= - FreeFallAcc * (sshn(m+1,n) -sshn(m,n))*dyh(m,n)*hhun_e(m,n)
                    grx = RHSx(m,n) + slx + RHSx_dif(m,n) + RHSx_adv(m,n)      &
                        - (rdis(m,n)+rdis(m+1,n))/2.0d0 *ubrtrp(m,n)*dxt(m,n)*dyh(m,n)*hhu(m,n)        &
                        + ( rlh_s(m,n  )*hhh(m,n  )*dxb(m,n  )*dyb(m,n  )*(vbrtr(m+1,n  ) + vbrtr(m,n  ))          &
                        +   rlh_s(m,n-1)*hhh(m,n-1)*dxb(m,n-1)*dyb(m,n-1)*(vbrtr(m+1,n-1) + vbrtr(m,n-1)) )/4.0d0

                    ubrtrn(m,n) = (ubrtrp(m,n)*bp0 + grx )/(bp - RHSx_bfc(m, n))
                endif

                !meridional flux
                if(lcv(m,n)>0.5) then
                    bp  = hhvn(m,n)*dyt(m,n)*dxh(m,n)/2.0d0/tau
                    bp0 = hhvp(m,n)*dyt(m,n)*dxh(m,n)/2.0d0/tau

                    sly = - FreeFallAcc * ( ssh(m,n+1)- ssh(m,n))*dxh(m,n)* hhv(m,n)
                    !slyn= - FreeFallAcc * (sshn(m,n+1)-sshn(m,n))*dxh(m,n)*hhvn_e(m,n)
                    gry = RHSy(m,n) + sly  + RHSy_dif(m,n) + RHSy_adv(m,n)      &
                        - (rdis(m,n)+rdis(m,n+1))/2.0d0 *vbrtrp(m,n)*dxh(m,n)*dyt(m,n)*hhv(m,n)        &
                        - ( rlh_s(m  ,n)*hhh(m  ,n)*dxb(m  ,n)*dyb(m  ,n)*(ubrtr(m  ,n+1)+ubrtr(m  ,n))             &
                        +   rlh_s(m-1,n)*hhh(m-1,n)*dxb(m-1,n)*dyb(m-1,n)*(ubrtr(m-1,n+1)+ubrtr(m-1,n)) )/4.0d0

                    vbrtrn(m,n) = (vbrtrp(m,n)*bp0 + gry )/(bp - RHSy_bfc(m, n))
                endif
            enddo
        enddo
        !$omp end parallel do

        call syncborder_real8(ubrtrn, 1)
        call syncborder_real8(vbrtrn, 1)
        if(periodicity_x/=0) then
            call cyclize8_x(ubrtrn,nx,ny,1,mmm,mm)
            call cyclize8_x(vbrtrn,nx,ny,1,mmm,mm)
        endif
        if(periodicity_y/=0) then
            call cyclize8_y(ubrtrn,nx,ny,1,nnn,nn)
            call cyclize8_y(vbrtrn,nx,ny,1,nnn,nn)
        endif

        !shifting time indices
        !$omp parallel do private(m, n)
        do n=ny_start-1,ny_end+1
            do m=nx_start-1,nx_end+1

                if(lu(m,n)>0.5) then
                    sshp(m,n) = ssh(m,n)+time_smooth*(sshn(m,n)-2.0d0*ssh(m,n)+sshp(m,n))/2.0d0
                    ssh(m,n) = sshn(m,n)
                endif
                if(lcu(m,n)>0.5) then
                    !up(m,n) =  hhu_e(m,n)*u(m,n)+time_smooth*(hhun_e(m,n)*un(m,n)-2.0d0*hhu_e(m,n)*u(m,n)+hhup_e(m,n)*up(m,n))/2.0d0/dfloat(nstep)
                    ubrtrp(m,n) = ubrtr(m,n) + time_smooth*(ubrtrn(m,n)-2.0d0*ubrtr(m,n)+ubrtrp(m,n))/2.0d0
                    ubrtr(m,n) = ubrtrn(m,n)
                endif
                if(lcv(m,n)>0.5) then
                    !vp(m,n) =  hhv_e(m,n)*v(m,n)+time_smooth*(hhvn_e(m,n)*vn(m,n)-2.0d0*hhv_e(m,n)*v(m,n)+hhvp_e(m,n)*vp(m,n))/2.0d0/dfloat(nstep)
                    vbrtrp(m,n) = vbrtr(m,n) + time_smooth*(vbrtrn(m,n)-2.0d0*vbrtr(m,n)+vbrtrp(m,n))/2.0d0
                    vbrtr(m,n) = vbrtrn(m,n)
                endif

            enddo
        enddo
        !$omp end parallel do

        if(full_free_surface>0) then
            call hh_shift(hhq, hhqp, hhqn,   &
                          hhu, hhup, hhun,   &
                          hhv, hhvp, hhvn,   &
                          hhh, hhhp, hhhn)
        endif

        !call syncborder_real8(ubrtr_i, 1)
        !call syncborder_real8(vbrtr_i, 1)
        !if(periodicity_x/=0) then
        !    call cyclize8_x(ubrtr_i,nx,ny,1,mmm,mm)
        !    call cyclize8_x(vbrtr_i,nx,ny,1,mmm,mm)
        !endif
        !if(periodicity_y/=0) then
        !    call cyclize8_y(ubrtr_i,nx,ny,1,nnn,nn)
        !    call cyclize8_y(vbrtr_i,nx,ny,1,nnn,nn)
        !endif

        if(full_free_surface>0) then
            !initialize depth for external mode
            call hh_init(hhq, hhqp, hhqn,    &
                         hhu, hhup, hhun,    &
                         hhv, hhvp, hhvn,    &
                         hhh, hhhp, hhhn,    &
                         ssh, sshp, hhq_rest)
        endif

    endsubroutine expl_shallow_water

endmodule
