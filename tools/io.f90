module io_module
    ! Input/Output module

#include "macros/kernel_macros.fi"

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    !use mpi
    use mpp_module
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

    interface write_data
        module procedure write_data2D_real4
        module procedure write_data2D_real8
    end interface

    public :: read_data
    public :: write_data
    public :: read_global_mask

contains

    subroutine read_global_mask(grid_global_data)
        use config_basinpar_module, only: nx, ny, mask_file_name

        type(grid_global_type), intent(inout) :: grid_global_data

        character :: frmt*16, comment*80
        integer :: m, n, ierr
        
        if (mpp_is_master_thread()) then

            write(frmt, 1000) nx
            1000  format('(',i9,'i1)')
            
            ! reading mask from:
            if (mpp_is_master()) then
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
        
        endif

        return
        98    call abort_model('Error in opening mask file')
        99    call abort_model('Error in reading mask file')
    end subroutine

    subroutine read_data2D_real4(domain, path, fname, nfild, data2d, lu2d, ierr)
        use config_basinpar_module, only: nx, ny, nxb => mmm, nxe => mm, nyb => nnn, nye => nn
        use iodata_routes, only: fulfname, lmpirecl

        type(domain_type), intent(in) :: domain
        type(data2D_real4_type), intent(inout) :: data2d
        type(data2D_real4_type), intent(in) :: lu2d
        character*(*), intent(in) :: path, fname
        integer, intent(in) :: nfild
        integer, intent(out) :: ierr
        
        integer :: m, n, k, kfull, i, j
        integer(kind=mpi_offset_kind) :: disp
        character(4096) :: namofile
        integer :: hfile
        integer :: tsubarr, totsize
        integer :: sizes2(2), locsizes2(2), offset2(2)

        if (mpp_is_master_thread()) then

            ! Definition full file name
            call fulfname(namofile, path, fname, ierr)
            if (ierr .ne. 0) then
                call abort_model("Error in full name of file for reading")
            endif

            disp = (nxe - nxb + 1)*(nye - nyb + 1) * int(lmpirecl, mpi_offset_kind) * (nfild - 1)
            !locsizes = (/nx_end - nx_start + 1, ny_end - ny_start + 1, nze - nzb + 1/)
            !totsize = locsizes(1)*locsizes(2)*locsizes(3)

            call mpi_file_open(mpp_cart_comm, namofile, mpi_mode_rdonly, mpi_info_null, hfile, ierr)
            if (ierr .ne. mpi_success) then
                call abort_model("Error in open file for reading: "//namofile(1 : len_trim(namofile)))
            endif

            do kfull = 1, domain%bcount_max
                if (kfull <= domain%bcount) then
                    k = kfull
                else 
                    k = 1
                endif

                associate(_associate_domain_value_(nx_start, domain, bnx_start, k),  &
                        _associate_domain_value_(nx_end,   domain, bnx_end,   k),  &
                        _associate_domain_value_(ny_start, domain, bny_start, k),  &
                        _associate_domain_value_(ny_end,   domain, bny_end,   k),  &
                        field => data2d%block(k)%field, lu => lu2d%block(k)%field)

                    offset2 = (/nx_start - nxb, ny_start - nyb/)
                    locsizes2 = (/nx_end - nx_start + 1, ny_end - ny_start + 1/)
                    sizes2 = (/nxe - nxb + 1, nye - nyb + 1/)

                    totsize = locsizes2(1)*locsizes2(2)
                    
                    !print *, sizes2, locsizes2, offset2
                    call mpi_type_create_subarray(2, sizes2, locsizes2, offset2,        &
                                                mpi_order_fortran, mpi_real, tsubarr, ierr)
                    
                    call mpi_type_commit(tsubarr, ierr)

                    call mpi_file_set_view(hfile, disp, mpi_real, tsubarr, "native", mpi_info_null, ierr)

                    call mpi_file_read_all(hfile, field(nx_start:nx_end, ny_start:ny_end),  &
                                        totsize, mpi_real, mpi_status_ignore, ierr)
                    if (ierr .ne. mpi_success) then
                        call abort_model("Error in reading from file: "//namofile(1 : len_trim(namofile)))
                    endif

                    call mpi_type_free(tsubarr, ierr)
                end associate
            enddo 
            call mpi_file_close(hfile, ierr)

            ! Filling undefinite points by zero instead undef
            do k = 1, domain%bcount
                associate(_associate_domain_value_(bnd_x1, domain, bbnd_x1, k),  &
                        _associate_domain_value_(bnd_x2, domain, bbnd_x2, k),  &
                        _associate_domain_value_(bnd_y1, domain, bbnd_y1, k),  &
                        _associate_domain_value_(bnd_y2, domain, bbnd_y2, k),  &
                        field => data2d%block(k)%field, lu => lu2d%block(k)%field)

                    do j = bnd_y1, bnd_y2
                        do i = bnd_x1, bnd_x2
                            if (abs(lu(i,j)) < 0.5) then
                                field(i, j) = 0.0
                            end if
                        enddo
                    enddo
                end associate
            enddo

        endif

        return
    end subroutine

    subroutine read_data2D_real8(domain, path, fname, nfild, data2d, lu2d, ierr)
        use config_basinpar_module, only: nx, ny, nxb => mmm, nxe => mm, nyb => nnn, nye => nn
        use iodata_routes, only: fulfname, lmpirecl

        type(domain_type), intent(in) :: domain
        type(data2D_real8_type), intent(inout) :: data2d
        type(data2D_real4_type), intent(in) :: lu2d
        character*(*), intent(in) :: path, fname
        integer, intent(in) :: nfild
        integer, intent(out) :: ierr
        
        integer :: m, n, k, kfull, i, j
        integer(kind=mpi_offset_kind) :: disp
        character(4096) :: namofile
        integer :: hfile
        integer :: tsubarr, totsize
        integer :: sizes2(2), locsizes2(2), offset2(2)

        if (mpp_is_master_thread()) then

            ! Definition full file name
            call fulfname(namofile, path, fname, ierr)
            if (ierr .ne. 0) then
                call abort_model("Error in full name of file for reading")
            endif

            disp = (nxe - nxb + 1)*(nye - nyb + 1) * 2*int(lmpirecl, mpi_offset_kind) * (nfild - 1)
            !locsizes = (/nx_end - nx_start + 1, ny_end - ny_start + 1, nze - nzb + 1/)
            !totsize = locsizes(1)*locsizes(2)*locsizes(3)

            call mpi_file_open(mpp_cart_comm, namofile, mpi_mode_rdonly, mpi_info_null, hfile, ierr)
            if (ierr .ne. mpi_success) then
                call abort_model("Error in open file for reading: "//namofile(1 : len_trim(namofile)))
            endif

            do kfull = 1, domain%bcount_max
                if (kfull <= domain%bcount) then
                    k = kfull
                else 
                    k = 1
                endif

                associate(_associate_domain_value_(nx_start, domain, bnx_start, k),  &
                        _associate_domain_value_(nx_end,   domain, bnx_end,   k),  &
                        _associate_domain_value_(ny_start, domain, bny_start, k),  &
                        _associate_domain_value_(ny_end,   domain, bny_end,   k),  &
                        field => data2d%block(k)%field, lu => lu2d%block(k)%field)

                    offset2 = (/nx_start - nxb, ny_start - nyb/)
                    locsizes2 = (/nx_end - nx_start + 1, ny_end - ny_start + 1/)
                    sizes2 = (/nxe - nxb + 1, nye - nyb + 1/)

                    totsize = locsizes2(1)*locsizes2(2)
                    
                    !print *, sizes2, locsizes2, offset2
                    call mpi_type_create_subarray(2, sizes2, locsizes2, offset2,        &
                                                mpi_order_fortran, mpi_real8, tsubarr, ierr)
                    
                    call mpi_type_commit(tsubarr, ierr)

                    call mpi_file_set_view(hfile, disp, mpi_real8, tsubarr, "native", mpi_info_null, ierr)

                    call mpi_file_read_all(hfile, field(nx_start:nx_end, ny_start:ny_end),  &
                                        totsize, mpi_real8, mpi_status_ignore, ierr)
                    if (ierr .ne. mpi_success) then
                        call abort_model("Error in reading from file: "//namofile(1 : len_trim(namofile)))
                    endif

                    call mpi_type_free(tsubarr, ierr)
                end associate
            enddo 
            call mpi_file_close(hfile, ierr)

            ! Filling undefinite points by zero instead undef
            do k = 1, domain%bcount
                associate(_associate_domain_value_(bnd_x1, domain, bbnd_x1, k),  &
                        _associate_domain_value_(bnd_x2, domain, bbnd_x2, k),  &
                        _associate_domain_value_(bnd_y1, domain, bbnd_y1, k),  &
                        _associate_domain_value_(bnd_y2, domain, bbnd_y2, k),  &
                        field => data2d%block(k)%field, lu => lu2d%block(k)%field)

                    do j = bnd_y1, bnd_y2
                        do i = bnd_x1, bnd_x2
                            if (abs(lu(i,j)) < 0.5) then
                                field(i, j) = 0.0d0
                            end if
                        enddo
                    enddo
                end associate
            enddo

        endif

        return
    end subroutine

    subroutine write_data2D_real4(domain, path, fname, nfild, data2d, lu2d, ierr)
        use config_basinpar_module, only: nx, ny, nxb => mmm, nxe => mm, nyb => nnn, nye => nn
        use iodata_routes, only: fulfname, lrecl, lmpirecl, undef

        type(domain_type), intent(in) :: domain
        type(data2D_real4_type), intent(in) :: data2d
        type(data2D_real4_type), intent(in) :: lu2d
        character*(*), intent(in) :: path, fname
        integer, intent(in) :: nfild
        integer, intent(out) :: ierr
        
        integer :: m, n, k, kfull, i, j
        integer(kind=mpi_offset_kind) :: disp
        character(4096) :: namofile
        integer :: hfile
        integer :: tsubarr, totsize
        integer :: sizes2(2), locsizes2(2), offset2(2)

        if (mpp_is_master_thread()) then

            ! Definition full file name
            call fulfname(namofile, path, fname, ierr)
            if (ierr .ne. 0) then
                call abort_model("Error in full name of file for reading")
            endif

            disp = (nxe - nxb + 1)*(nye - nyb + 1) * int(lmpirecl, mpi_offset_kind) * (nfild - 1)
            !locsizes = (/nx_end - nx_start + 1, ny_end - ny_start + 1, nze - nzb + 1/)
            !totsize = locsizes(1)*locsizes(2)*locsizes(3)

            if (mpp_is_master()) then
                open(40, file=namofile, status='unknown', access='direct', form='unformatted',recl=(nxe - nxb + 1)*(nye - nyb + 1)*lrecl, err=101)
                ! writing on the file
                write(40, rec=nfild, err=102) ((undef, i=nxb,nxe), j=nyb,nye)
                close(40)
            endif

            call mpi_file_open(mpp_cart_comm, namofile, ior(mpi_mode_wronly, mpi_mode_create), mpi_info_null, hfile, ierr)
            if (ierr .ne. mpi_success) then
                call abort_model("Error in open file for writing: "//namofile(1 : len_trim(namofile)))
            endif

            do kfull = 1, domain%bcount_max
                if (kfull <= domain%bcount) then
                    k = kfull
                else 
                    k = 1
                endif

                associate(_associate_domain_value_(nx_start, domain, bnx_start, k),  &
                        _associate_domain_value_(nx_end,   domain, bnx_end,   k),  &
                        _associate_domain_value_(ny_start, domain, bny_start, k),  &
                        _associate_domain_value_(ny_end,   domain, bny_end,   k),  &
                        field => data2d%block(k)%field, lu => lu2d%block(k)%field)
                
                    offset2 = (/nx_start - nxb, ny_start - nyb/)
                    locsizes2 = (/nx_end - nx_start + 1, ny_end - ny_start + 1/)
                    sizes2 = (/nxe - nxb + 1, nye - nyb + 1/)
                    totsize = locsizes2(1)*locsizes2(2)

                    call mpi_type_create_subarray(2, sizes2, locsizes2, offset2,               &
                                                mpi_order_fortran, mpi_real, tsubarr, ierr)
                
                    call mpi_type_commit(tsubarr, ierr)

                    call mpi_file_set_view(hfile, disp, mpi_real, tsubarr, "native", mpi_info_null, ierr)
                    
                    do j = ny_start, ny_end
                        do i = nx_start, nx_end
                            if (abs(lu(i,j)) .lt. 0.5) then
                                field(i, j) = undef
                            end if
                        enddo
                    enddo
                    
                    call mpi_file_write_all(hfile, field(nx_start:nx_end, ny_start:ny_end),  &
                                            totsize, mpi_real, mpi_status_ignore, ierr)
                    if (ierr .ne. mpi_success) then
                        call abort_model("Error in writig in file: "//namofile(1 : len_trim(namofile)))
                    endif
                    
                    call mpi_type_free(tsubarr, ierr)
                end associate
            enddo

            call mpi_file_close(hfile, ierr)

            do k = 1, domain%bcount
                associate(_associate_domain_value_(nx_start, domain, bnx_start, k),  &
                    _associate_domain_value_(nx_end,   domain, bnx_end,   k),  &
                    _associate_domain_value_(ny_start, domain, bny_start, k),  &
                    _associate_domain_value_(ny_end,   domain, bny_end,   k),  &
                    field => data2d%block(k)%field, lu => lu2d%block(k)%field)

                    do j = ny_start, ny_end
                        do i = nx_start, nx_end
                            if (abs(lu(i,j)) .lt. 0.5) then
                                field(i, j) = 0.0
                            end if
                        enddo
                    enddo

                end associate
            enddo

        endif

        return
        101    call abort_model("Error in pre-opening file for writing undefs: "//namofile(1 : len_trim(namofile)))
        102    call abort_model("Error in pre-writing undefs in file: "//namofile(1 : len_trim(namofile)))
    end subroutine

    subroutine write_data2D_real8(domain, path, fname, nfild, data2d, lu2d, ierr)
        use config_basinpar_module, only: nx, ny, nxb => mmm, nxe => mm, nyb => nnn, nye => nn
        use iodata_routes, only: fulfname, lrecl, lmpirecl, undef

        type(domain_type), intent(in) :: domain
        type(data2D_real8_type), intent(in) :: data2d
        type(data2D_real4_type), intent(in) :: lu2d
        character*(*), intent(in) :: path, fname
        integer, intent(in) :: nfild
        integer, intent(out) :: ierr
        
        integer :: m, n, k, kfull, i, j
        integer(kind=mpi_offset_kind) :: disp
        character(4096) :: namofile
        integer :: hfile
        integer :: tsubarr, totsize
        integer :: sizes2(2), locsizes2(2), offset2(2)

        if (mpp_is_master_thread()) then

            ! Definition full file name
            call fulfname(namofile, path, fname, ierr)
            if (ierr .ne. 0) then
                call abort_model("Error in full name of file for reading")
            endif

            disp = (nxe - nxb + 1)*(nye - nyb + 1) * 2*int(lmpirecl, mpi_offset_kind) * (nfild - 1)
            !locsizes = (/nx_end - nx_start + 1, ny_end - ny_start + 1, nze - nzb + 1/)
            !totsize = locsizes(1)*locsizes(2)*locsizes(3)

            if (mpp_is_master()) then
                open(40, file=namofile, status='unknown', access='direct', form='unformatted',recl=(nxe - nxb + 1)*(nye - nyb + 1)*2*lrecl, err=101)
                ! writing on the file
                write(40, rec=nfild, err=102) ((undef, i=nxb,nxe), j=nyb,nye)
                close(40)
            endif

            call mpi_file_open(mpp_cart_comm, namofile, ior(mpi_mode_wronly, mpi_mode_create), mpi_info_null, hfile, ierr)
            if (ierr .ne. mpi_success) then
                call abort_model("Error in open file for writing: "//namofile(1 : len_trim(namofile)))
            endif

            do kfull = 1, domain%bcount_max
                if (kfull <= domain%bcount) then
                    k = kfull
                else 
                    k = 1
                endif

                associate(_associate_domain_value_(nx_start, domain, bnx_start, k),  &
                        _associate_domain_value_(nx_end,   domain, bnx_end,   k),  &
                        _associate_domain_value_(ny_start, domain, bny_start, k),  &
                        _associate_domain_value_(ny_end,   domain, bny_end,   k),  &
                        field => data2d%block(k)%field, lu => lu2d%block(k)%field)
                
                    offset2 = (/nx_start - nxb, ny_start - nyb/)
                    locsizes2 = (/nx_end - nx_start + 1, ny_end - ny_start + 1/)
                    sizes2 = (/nxe - nxb + 1, nye - nyb + 1/)
                    totsize = locsizes2(1)*locsizes2(2)

                    call mpi_type_create_subarray(2, sizes2, locsizes2, offset2,               &
                                                mpi_order_fortran, mpi_real8, tsubarr, ierr)
                
                    call mpi_type_commit(tsubarr, ierr)

                    call mpi_file_set_view(hfile, disp, mpi_real8, tsubarr, "native", mpi_info_null, ierr)
                    
                    do j = ny_start, ny_end
                        do i = nx_start, nx_end
                            if (abs(lu(i,j)) .lt. 0.5) then
                                field(i, j) = undef
                            end if
                        enddo
                    enddo
                    
                    call mpi_file_write_all(hfile, field(nx_start:nx_end, ny_start:ny_end),  &
                                            totsize, mpi_real8, mpi_status_ignore, ierr)
                    if (ierr .ne. mpi_success) then
                        call abort_model("Error in writig in file: "//namofile(1 : len_trim(namofile)))
                    endif
                    
                    call mpi_type_free(tsubarr, ierr)
                end associate
            enddo

            call mpi_file_close(hfile, ierr)

            do k = 1, domain%bcount
                associate(_associate_domain_value_(nx_start, domain, bnx_start, k),  &
                    _associate_domain_value_(nx_end,   domain, bnx_end,   k),  &
                    _associate_domain_value_(ny_start, domain, bny_start, k),  &
                    _associate_domain_value_(ny_end,   domain, bny_end,   k),  &
                    field => data2d%block(k)%field, lu => lu2d%block(k)%field)

                    do j = ny_start, ny_end
                        do i = nx_start, nx_end
                            if (abs(lu(i,j)) .lt. 0.5) then
                                field(i, j) = 0.0d0
                            end if
                        enddo
                    enddo

                end associate
            enddo

        endif
        
        return
        101    call abort_model("Error in pre-opening file for writing undefs: "//namofile(1 : len_trim(namofile)))
        102    call abort_model("Error in pre-writing undefs in file: "//namofile(1 : len_trim(namofile)))
    end subroutine


end module io_module