module config_basinpar_module
    ! Base parameters for model

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use mpp_module
    use rwpar_routes

    implicit none
    save
    public

    integer :: nx             ! number of total points in x-direction in base arrays
    integer :: ny             ! number of total points in y-direction in base arrays
    integer :: nz             ! number of s-levels in vertical direction for 3d arrays
    integer :: mmm            ! begin of significant area in x-direction
    integer :: nnn            ! begin of significant area in y-direction
    integer :: mm             ! end of significant area in x-direction
    integer :: nn             ! end of significant area in y-direction
    integer :: periodicity_x  ! Periodicity on x-direction (0 - non-periodic, 1 - periodic)
    integer :: periodicity_y  ! Periodicity on y-direction (0 - non-periodic, 1 - periodic) 
    character(len=256) :: mask_file_name ! name of file with temperature point sea-land mask
    character(len=256) :: bottom_topography_file_name ! name of file with bottom topography

    real(wp8) :: dxst     ! longitude step (in degrees) in case of regular grid
    real(wp8) :: dyst     ! latitude  step (in degrees) in case of regular grid
    real(wp8) :: rlon     ! first calcutable t-point (m=mmm) on latitude  (in degrees)
    real(wp8) :: rlat     ! first calcutable t-point (n=nnn) on longitude (in degrees)
    
    integer :: xgr_type
    integer :: ygr_type   ! grid types: 0 - regular, 1 - levels

    real(wp8), allocatable :: x_levels(:), y_levels(:)

    ! parameters of basin        ! =0 - Carthesian coordinates 
                                 !  in degrees, with metric coefficient rx=R, ry=R (R - Earth radius)
    integer :: curve_grid        ! =1 - undistorted spherical grid (poles are placed at the sphere axis, but not only Earth rotation axis)
                                 !  in degrees, with metric coefficients rx=R*cos(lat) and ry=R (R - Earth radius)
                                 ! =2 - distorted spherical grid (poles are placed beyond the sphere axis) 
                                 !  in degrees, with complicated metric coefficients

    ! Parameters of rotation (Euler angles in case of undistorted sphere)
    real(wp8) :: rotation_on_lon      ! geographic coordinates of the point  
    real(wp8) :: rotation_on_lat      ! that has coordinates (0,0) on rotated grid

    ! Parameters of curvature(in case of curvilinear cordinate system)
    real(wp8) :: x_pole,   &  ! lon of new north pole in geographical system
                 y_pole,   &  ! lat of new north pole in geographical system
                 p_pole,   &  ! lon of new south pole in geographical system
                 q_pole       ! lat of new south pole in geographical system
    
contains

    subroutine load_config_basinpar_from_file(filepar)
        character(*), intent(in) :: filepar
        character(256) :: comments(256)
        integer :: nofcom, ierr

        ! reading parameters from file
        if (mpp_is_master()) then
            call readpar(filepar, comments, nofcom)
        endif
        call mpi_bcast(comments, 256*256, mpi_character, 0, mpp_cart_comm, ierr)

        read(comments( 1),*) nx
        read(comments( 2),*) ny
        read(comments( 3),*) nz
        read(comments( 4),*) periodicity_x
        read(comments( 5),*) periodicity_y
        read(comments( 6),*) dxst
        read(comments( 7),*) dyst
        read(comments( 8),*) rlon
        read(comments( 9),*) rlat
        read(comments( 10),*) xgr_type
        read(comments( 11),*) ygr_type
        read(comments( 12),*) curve_grid
        read(comments( 13),*) rotation_on_lon
        read(comments( 14),*) rotation_on_lat
        read(comments( 15),*) x_pole
        read(comments( 16),*) y_pole
        read(comments( 17),*) p_pole
        read(comments( 18),*) q_pole
        call get_first_lexeme(comments(19), mask_file_name)
        call get_first_lexeme(comments(20), bottom_topography_file_name)

        ! Init config dependent parameters
        mmm = 3
        nnn = 3
        mm = nx-2
        nn = ny-2
        if (xgr_type > 0) allocate(x_levels(nx))
        if (ygr_type > 0) allocate(y_levels(ny))
    
        call mpp_sync_output()
    end subroutine

    subroutine load_config_basinpar()
        
        !call load_config_basinpar_as250m()
        call load_config_basinpar_as250m_test()
    end subroutine

    subroutine load_config_basinpar_bs4km()
        nx = 289
        ny = 163
        nz = 1
        mmm = 3
        nnn = 3
        mm = nx-2
        nn = ny-2
        periodicity_x = 0
        periodicity_y = 0

        rlon = 27.525d0
        rlat = 40.940d0
        dxst = 0.05d0
        dyst = 0.04d0

        xgr_type = 0
        ygr_type = 0

        if (xgr_type > 0) allocate(x_levels(nx))
        if (ygr_type > 0) allocate(y_levels(ny))

        curve_grid = 1

        rotation_on_lon = 0.0d0
        rotation_on_lat = 0.0d0

        x_pole = 90.0d0
        y_pole = 60.0d0
        p_pole = 90.0d0
        q_pole = -90.0d0

        mask_file_name = 'data/BS4km/mask.txt'
        bottom_topography_file_name = 'data/BS4km/topo.dat'
    end subroutine

    subroutine load_config_basinpar_as250m()
        nx = 1525
        ny = 1115
        nz = 1
        mmm = 3
        nnn = 3
        mm = nx-2
        nn = ny-2
        periodicity_x = 0
        periodicity_y = 0

        rlon = 34.751560d0
        rlat = 44.801125d0
        dxst = 0.00312d0
        dyst = 0.00225d0

        xgr_type = 0
        ygr_type = 0

        if (xgr_type > 0) allocate(x_levels(nx))
        if (ygr_type > 0) allocate(y_levels(ny))

        curve_grid = 1

        rotation_on_lon = 0.0d0
        rotation_on_lat = 0.0d0

        x_pole = 90.0d0
        y_pole = 60.0d0
        p_pole = 90.0d0
        q_pole = -90.0d0

        mask_file_name = 'data/AS/maskAzovCor.txt'
        bottom_topography_file_name = 'data/AS/TopoAzovSea_250x250m_3.5m.dat'
    end subroutine

    subroutine load_config_basinpar_as250m_test()
        nx = 1525
        ny = 1115
        nz = 1
        mmm = 3
        nnn = 3
        mm = nx-2
        nn = ny-2
        periodicity_x = 0
        periodicity_y = 0

        rlon = 34.751560d0
        rlat = 44.801125d0
        dxst = 0.00312d0
        dyst = 0.00225d0

        xgr_type = 0
        ygr_type = 0

        if (xgr_type > 0) allocate(x_levels(nx))
        if (ygr_type > 0) allocate(y_levels(ny))

        curve_grid = 1

        rotation_on_lon = 0.0d0
        rotation_on_lat = 0.0d0

        x_pole = 90.0d0
        y_pole = 60.0d0
        p_pole = 90.0d0
        q_pole = -90.0d0

        mask_file_name = 'none'
        bottom_topography_file_name = 'none'
    end subroutine

end module config_basinpar_module