module basinpar_module
    ! Initializing basin grid parameters

    use mpp_module
    use constants_module, only: RadEarth, EarthAngVel, HeatCapWater, RefDen, FreeFallAcc, pip180
    use config_basinpar_module, only: nx, ny, nz, mm, mmm, nn, nnn, dxst, dyst, rlon, rlat,    &
                                      curve_grid, xgr_type, ygr_type, x_levels, y_levels,      &
                                      rotation_on_lon, rotation_on_lat, x_pole, y_pole, p_pole, q_pole
    use decomposition_module, only: domain_type
    use grid_module, only: grid_type
    use grid_interface_module, only: envoke_grid_base_init_kernel, envoke_grid_geo_init_kernel

    implicit none
    save
    private

    public :: basinpar

contains

!===========================================================================================

    subroutine basinpar(domain, grid_data)
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data

        ! parameters:
        if (mpp_is_master()) then
            write(*,'(2x,a)')' Basin parameters from 1basinpar.inc:'

            if (curve_grid==0) then        ! Carthesian coordinates
                write(*,*) 'Coordinate system is carthesian'
            elseif (curve_grid==1) then
                write(*,*) 'Coordinate system is undistorted sphere'
                write(*,'(a,f10.3)') ' rotation angle on longitude is =',rotation_on_lon,    &
                                     ' rotation angle on  latitude is =',rotation_on_lat
            elseif (curve_grid==2) then
                write(*,*) 'Coordinate system is distorted sphere'
                write(*,'(a,f10.3)') ' geo longitude of new north pole is =',x_pole,    &
                                     ' geo  latitude of new north pole is =',y_pole,    &
                                     ' geo longitude of new south pole is =',p_pole,    &
                                     ' geo  latitude of new south pole is =',q_pole      
            endif

            if (xgr_type==0) then
                write(*,*) 'X-grid is uniform'
                write(*,'(2(a,f10.3),a)') ' initial x-coordinate (m=mmm) =',rlon,' step on x =',dxst,'[dgr] '
            else
                write(*,*) 'X-grid is non-uniform'      
                !write(*,'(a,f10.3)') ' minimal x-coordinate (m=mmm) =',xt(mmm),    &
                !                     ' maximal x-coordinate (m=mm ) =',xt(mm)
            endif

            if (ygr_type==0) then
                write(*,*) 'Y-grid is uniform'
                write(*,'(2(a,f10.3),a)') ' initial y-coordinate (n=nnn) =',rlat,' step on y =',dyst,'[dgr] '
            else
                write(*,*) 'Y-grid is non-uniform'      
                !write(*,'(a,f10.3)') ' minimal y-coordinate (n=nnn) =',yt(nnn),    &
                !                     ' maximal y-coordinate (n=nn ) =',yt(nn)
            endif

            !write(*,'(2(a,i2))') 'Periodicity on X =', periodicity_x,', Periodicity on Y =', periodicity_y
            write(*,'(4(a,i4))') '  nx=',nx, ';  ny=',ny,';  nz=',nz
            write(*,'(4(a,i4))') ' mmm=',mmm,';  mm=',mm,'; nnn=',nnn,';  nn=',nn

            write(*,'(2x,a,g14.7,a)')' Earth radius =',RadEarth,'(m)'
            write(*,'(2x,a,g14.7,a)') 'Earth angular velocity(omega) =',EarthAngVel,'[rad/sec]'
            write(*,'(2x,a,g14.7,a)') 'Heat capacity of water =',HeatCapWater,'[J/kg/�C] for 35%. sal'
            write(*,'(2x,a,g14.7,a)') 'reference density =',RefDen,'[kg/m**3]'
            write(*,'(2x,a,f10.3,a)') 'free fall acceleration(grv)=',FreeFallAcc,'[m/s**2]'
        endif     

        call envoke_grid_base_init_kernel(domain, grid_data)
        call envoke_grid_geo_init_kernel(domain, grid_data)

endsubroutine basinpar

end module basinpar_module
