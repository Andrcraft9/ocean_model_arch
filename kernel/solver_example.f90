module solver
    ! Module description
    ! ...

    use data, only: wp => SHR_KIND_R8

    implicit none
    save
    private

    ! Parameters
    real(kind=wp), public :: alpha

    ! Public subroutines
    public :: solve

contains

    subroutine solve(x, y)
        ! Subroutine description
        !
        ! INPUT/OUTPUT:
        ! Grid data:
        !
        ! Ocean data:
        !
        ! LOCAL:
        
        real(kind=wp), intent(in) :: x
        real(kind=wp), intent(out) :: y

        y = y + alpha*x
    end subroutine

end module