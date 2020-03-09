module math_tools_module
    
    implicit none
    save
    public

#ifndef __INTEL_COMPILER

    real*8, parameter:: m_dpi = 3.14159265358979d0
    real, parameter :: m_pi = 3.141597

    contains

    function cosd(x)
        real :: x 
        real :: cosd

        cosd = cos((x/180.0)*m_pi)
    end function

    function sind(x)
        real :: x 
        real :: sind

        sind = sin((x/180.0)*m_pi)
    end function
    
    function dcosd(x)
        real*8 :: x 
        real*8 :: dcosd

        dcosd = dcos((x/180.0d0)*m_dpi)
    end function

    function dsind(x)
        real*8 :: x 
        real*8 :: dsind

        dsind = dsin((x/180.0d0)*m_dpi)
    end function

    function dasind(x)
        real*8 :: x 
        real*8 :: dasind

        dasind = dasin(x)/m_dpi*180.0d0
    end function

    function dacosd(x)
        real*8 :: x 
        real*8 :: dacosd

        dacosd = dacos(x)/m_dpi*180.0d0
    end function

    function dtand(x)
        real*8 :: x 
        real*8 :: dtand

        dtand = dtan((x/180.0d0)*m_dpi)
    end function

#endif

end module math_tools_module