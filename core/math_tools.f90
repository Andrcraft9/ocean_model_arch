module math_tools_module
    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use constants_module, only: m_dpi => dPi, m_pi => Pi

    implicit none
    save
    public

#ifndef __INTEL_COMPILER

    contains

    function cosd(x)
        real(wp4) :: x 
        real(wp4) :: cosd

        cosd = cos((x/180.0)*m_pi)
    end function

    function sind(x)
        real(wp4) :: x 
        real(wp4) :: sind

        sind = sin((x/180.0)*m_pi)
    end function
    
    function dcosd(x)
        real(wp8) :: x 
        real(wp8) :: dcosd

        dcosd = dcos((x/180.0d0)*m_dpi)
    end function

    function dsind(x)
        real(wp8) :: x 
        real(wp8) :: dsind

        dsind = dsin((x/180.0d0)*m_dpi)
    end function

    function dasind(x)
        real(wp8) :: x 
        real(wp8) :: dasind

        dasind = dasin(x)/m_dpi*180.0d0
    end function

    function dacosd(x)
        real(wp8) :: x 
        real(wp8) :: dacosd

        dacosd = dacos(x)/m_dpi*180.0d0
    end function

    function dtand(x)
        real(wp8) :: x 
        real(wp8) :: dtand

        dtand = dtan((x/180.0d0)*m_dpi)
    end function

#endif

end module math_tools_module