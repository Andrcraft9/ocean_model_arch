module config_basinpar_module
    ! Base parameters for model

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4

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

    subroutine load_config_basinpar()
    ! Here must be reading from config file (parsing)

        nx = 289
        ny = 163
        nz = 20
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

        mask_file_name = '/home/andr/code/fortran/ocean_model_arch/data/BS4km/mask.txt'
        bottom_topography_file_name = 'none'
    end subroutine

end module config_basinpar_module