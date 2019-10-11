module basinpar_module
    ! Base parameters for model
    !*** Все Inc файлы в модули, возможно брать все параметры из конфигов, оставить параметрами только констнаты, типо GravAcc и тд. ***!

    use kind_module, only: wp8 => SHR_KIND_R8

    implicit none
    save
    public
    
    integer, parameter :: nx = 64,  &  ! number of total points in x-direction in base arrays
                          ny = 64      ! number of total points in y-direction in base arrays

    real(wp8), parameter :: dxst = 0.5d0,  &  ! longitude step (in degrees) in case of regular grid
                            dyst = 0.5d0      ! latitude  step (in degrees) in case of regular grid

end module basinpar_module