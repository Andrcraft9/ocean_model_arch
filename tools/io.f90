module io_module
    ! Input/Output module

    use parallel_module, only: procs_type
    use decomposition_module, only: domain_type
    use data_types_module, only: data2D_real4_type, data2D_real8_type

    implicit none
    save
    private

    public :: read_2D_real4, read_2D_real8
    public :: write_2D_real4, write_2D_real8

contains

    subroutine read_2D_real4(procs, domain, data_field)
        ! Read real4
        type(procs_type), intent(in) :: procs
        type(domain_type), intent(in) :: domain
        type(data2D_real4_type), intent(out) :: data_field
        integer :: k, m, n

        do k = 1, domain%bcount
            do m = domain%bnx_start(k), domain%bnx_end(k)
                do n = domain%bny_start(k), domain%bny_end(k)
                    associate(field => data_field%block(k)%field)
                        
                        field(m, n) = 0.0
                    
                    end associate
                enddo
            enddo
        enddo
    end subroutine

    subroutine read_2D_real8(procs, domain, data_field)
        ! Read real8
        type(procs_type), intent(in) :: procs
        type(domain_type), intent(in) :: domain
        type(data2D_real8_type), intent(out) :: data_field
        integer :: k, m, n

        do k = 1, domain%bcount
            do m = domain%bnx_start(k), domain%bnx_end(k)
                do n = domain%bny_start(k), domain%bny_end(k)
                    associate(field => data_field%block(k)%field)
                        
                        field(m, n) = 0.0
                    
                    end associate
                enddo
            enddo
        enddo
    end subroutine

    subroutine write_2D_real4(procs, domain, data_field)
        ! Write real4
        type(procs_type), intent(in) :: procs
        type(domain_type), intent(in) :: domain
        type(data2D_real4_type), intent(in) :: data_field
        integer :: k, m, n

        do k = 1, domain%bcount
            do m = domain%bnx_start(k), domain%bnx_end(k)
                do n = domain%bny_start(k), domain%bny_end(k)
                    associate(field => data_field%block(k)%field)
                        
                        print *, m, n, field(m, n)
                    
                    end associate
                enddo
            enddo
        enddo
    end subroutine

    subroutine write_2D_real8(procs, domain, data_field)
        ! Write real4
        type(procs_type), intent(in) :: procs
        type(domain_type), intent(in) :: domain
        type(data2D_real8_type), intent(in) :: data_field
        integer :: k, m, n

        do k = 1, domain%bcount
            do m = domain%bnx_start(k), domain%bnx_end(k)
                do n = domain%bny_start(k), domain%bny_end(k)
                    associate(field => data_field%block(k)%field)
                        
                        print *, m, n, field(m, n)
                    
                    end associate
                enddo
            enddo
        enddo
    end subroutine

end module io_module