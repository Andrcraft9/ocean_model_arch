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
    real(wp8) :: dxst         ! longitude step (in degrees) in case of regular grid
    real(wp8) :: dyst         ! latitude  step (in degrees) in case of regular grid
    character(len=256) :: mask_file_name ! name of file with temperature point sea-land mask
    character(len=256) :: bottom_topography_file_name ! name of file with bottom topography
    
contains

    subroutine load_config_basinpar()
    ! Here must be reading from config file (parsing)

        nx = 20
        ny = 20
        nz = 10
        mmm = 3
        nnn = 3
        mm = nx-2
        nn = ny-2
        periodicity_x = 0
        periodicity_y = 0

        dxst = 0.5d0
        dyst = 0.5d0

        mask_file_name = './mask.txt'
        bottom_topography_file_name = 'none'
    end subroutine

end module config_basinpar_module