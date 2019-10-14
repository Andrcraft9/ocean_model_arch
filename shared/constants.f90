module constants_module
    ! Module description

    implicit none
    save
    public

    real, parameter :: RADIUS = 6371.0e3   
    real, parameter :: GRAV   = 9.80    
    real, parameter :: RHO0    = 1.035e3

end module constants_module