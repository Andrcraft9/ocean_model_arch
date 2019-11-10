module basinpar_module
    ! Base parameters for model
    !*** Все Inc файлы в модули, возможно брать все параметры из конфигов, оставить параметрами только констнаты, типо GravAcc и тд. ***!

    use kind_module, only: wp8 => SHR_KIND_R8

    implicit none
    save
    public
    
    integer, parameter :: nx = 20,            &  ! number of total points in x-direction in base arrays
                          ny = 20,            &  ! number of total points in y-direction in base arrays
                          nz = 10,            &  ! number of s-levels in vertical direction for 3d arrays
                          mmm = 3,            &  ! begin of significant area in x-direction
                          nnn = 3,            &  ! begin of significant area in y-direction
                          mm = nx-2,          &  ! end of significant area in x-direction
                          nn = ny-2,          &  ! end of significant area in y-direction
                          periodicity_x = 0,  &  ! Periodicity on x-direction (0 - non-periodic, 1 - periodic)
                          periodicity_y = 0      ! Periodicity on y-direction (0 - non-periodic, 1 - periodic) 

    real(wp8), parameter :: dxst = 0.5d0,  &  ! longitude step (in degrees) in case of regular grid
                            dyst = 0.5d0      ! latitude  step (in degrees) in case of regular grid

end module basinpar_module