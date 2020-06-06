module system_module

    implicit none
    save
    public

    ! Input/Output
    integer, parameter:: lrecl = 4  ! long of unique recl for PGI, gfortran, efc. For Intel use lrecl = 1 or flag -assume byterecl
    integer, parameter :: lmpirecl = 4  ! for mpi read/write
    real, parameter:: undef = -1.0e+32  ! missing data value

end module system_module