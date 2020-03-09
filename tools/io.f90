module io_module
    ! Input/Output module

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use mpi
    use mpp_module, only: mpp_rank, mpp_count, mpp_cart_comm, mpp_size, mpp_coord
    use errors_module, only: abort_model, check_error
    use data_types_module, only: data1D_real4_type, data1D_real8_type, data2D_real4_type, data2D_real8_type,      &
                                 data3D_real4_type, data3D_real8_type
    use decomposition_module, only: domain_type
    use grid_module, only: grid_global_type
    use, intrinsic :: iso_fortran_env, only: file_storage_size

    implicit none
    save
    private

    interface read_data
        module procedure read_data2D_real4
        module procedure read_data2D_real8
    end interface

    public :: read_data
    public :: read_global_mask

contains

    subroutine read_global_mask(grid_global_data)
        use config_basinpar_module, only: nx, ny, mask_file_name

        type(grid_global_type), intent(inout) :: grid_global_data

        character :: frmt*16, comment*80
        integer :: m, n, ierr
        
        write(frmt, 1000) nx
        1000  format('(',i9,'i1)')
        
        ! reading mask from:
        if (mpp_rank == 0) then
          print *, "Reading mask of computational area..."
          print *, mask_file_name, nx, ny, file_storage_size
          open (11, file=mask_file_name, status='old', recl=nx*file_storage_size, err=98)
            read (11,  '(a)') comment(1:min(80,nx))
            write(*,'(1x,a)') comment
            do n = ny, 1, -1
              read(11, frmt, end=99, err=99) (grid_global_data%mask(m,n), m=1,nx)
            enddo
          close(11)
          print *, "...OK"
        endif
        call mpi_bcast(grid_global_data%mask, nx*ny, mpi_integer, 0, mpp_cart_comm, ierr)
        
        return
        98   call abort_model('Error in opening mask file')
        99   call abort_model('Error in reading mask file')
    end subroutine

    subroutine read_data2D_real4(domain, path, fname, nfild, data2d, lu, ierr)
        type(domain_type), intent(in) :: domain
        type(data2D_real4_type), intent(inout) :: data2d
        type(data2D_real4_type), intent(in) :: lu
        character*(*), intent(in) :: path, fname
        integer, intent(in) :: nfild
        integer, intent(out) :: ierr

        integer :: m, n, k

        do k = 1, domain%bcount
            associate(nx_start => domain%bnx_start(k), nx_end => domain%bnx_end(k),  &
                      ny_start => domain%bny_start(k), ny_end => domain%bny_end(k),  &
                      field => data2d%block(k)%field)
            
            do n = ny_start, ny_end
                do m = nx_start, nx_end
                    field(m, n) = 0.0
                enddo
            enddo

            end associate
        enddo

    end subroutine

    subroutine read_data2D_real8(domain, path, fname, nfild, data2d, lu, ierr)
        type(domain_type), intent(in) :: domain
        type(data2D_real8_type), intent(inout) :: data2d
        type(data2D_real4_type), intent(in) :: lu
        character*(*), intent(in) :: path, fname
        integer, intent(in) :: nfild
        integer, intent(out) :: ierr

        integer :: m, n, k

        do k = 1, domain%bcount
            associate(nx_start => domain%bnx_start(k), nx_end => domain%bnx_end(k),  &
                      ny_start => domain%bny_start(k), ny_end => domain%bny_end(k),  &
                      field => data2d%block(k)%field)
            
            do n = ny_start, ny_end
                do m = nx_start, nx_end
                    field(m, n) = 0.0d0
                enddo
            enddo

            end associate
        enddo

    end subroutine

end module io_module