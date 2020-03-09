module iodata_routes
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

    !integer, parameter:: lrecl=1         !long of unique recl for Compaq, Intel fortran
    integer, parameter:: lrecl=4         !long of unique recl for PGI, gfortran, efc

    integer, parameter :: lmpirecl=4  !for mpi read/write
    real, parameter:: undef=-1.0e+32  !missing data value
    integer, parameter:: optimized_io = 0  !Optimized input/output, see pdrstd and pwdstd

contains

!======================================================================
!  version#1 of input-output subroutines
!======================================================================
subroutine fulfname(name,path,filen,ierr)

!  definit sum: name = path + filen removing blanks from path
      character*(*) path, filen, name
      character*1 fndevider
      integer n, lp, ierr, np1, np2, lf, nf1, nf2
!  fndevider - devider in full file name with path
!  change it according to operation sistem
!  for dos:
!     fndevider='\'
!  for unix(dos also):
      fndevider='/'

      ierr=0
!     name=filen
!     return
      do n=1,len(name)
         name(n:n)=' '
      end do
      lp=len(path)
      if(lp<1) then
         ierr=1
         write(*,'(2x,a)')'error in subroutine fulfname:'
         write(*,'(2x,a,a)')' error in path to file ',   name(1:len_trim(name))
         return
      end if
!  fined real initial and last definit position in path without blank
          np1=0
      if(lp>1) then
          n=1
          do while(n<lp)
             if(path(n:n)==' '.and.path(n+1:n+1)/=' ') np1=n
             if(path(n:n)/=' '.and.path(n+1:n+1)==' ') exit
             n=n+1
          end do
          if(n==lp.and.path(lp:lp)==' ') then
             np2 = 0
          else
             np2 = n
          end if
      else
          if(path(1:1)==' ') then
              np2=0
          else
              np2=1
          end if
      end if
      lp=np2-np1                !long of real path witout blanks

      lf=len(filen)
      if(lf<1) then
         ierr=1
         write(*,'(2x,a)') 'error in subroutine fulfname:'
         write(*,'(2x,a)') 'error in file name: '
         write(*,*)'path: ', path(1:len_trim(path)), ';  name:',  name(1:len_trim(path))
         return
      end if
!  fined real initial and last definit position in filename without blank
          nf1=0
      if(lf>1) then
          n=1
          do while(n<lf)
             if(filen(n:n)==' '.and.filen(n+1:n+1)/=' ') nf1=n
             if(filen(n:n)/=' '.and.filen(n+1:n+1)==' ') exit
             n=n+1
          end do
          if(n==lf.and.filen(lf:lf)==' ') then
             nf2=0
          else
             nf2 = n
          end if
      else
          if(filen(1:1)==' ') then
              nf2=0
          else
              nf2=1
          end if
      end if
      lf=nf2-nf1                !long of real filename witout blanks

      if(lf<=0) then
      ierr=1
      write(*,'(2x,a)')'error in subroutine fulfname:'
      write(*,'(2x,a)')'there is no file name!'
      return
      end if

      if(lp+lf>len(name)) then
      ierr=1
      write(*,'(2x,a)')'error in subroutine fulfname:'
      write(*,'(2x,a)') 'error in file name: '
      write(*,*)'path: ', path(1:len_trim(path)), ';  name:',  name(1:len_trim(name))
      write(*,'(2x,a)')'len of fulname < path+filename:'
      return
      end if

      if(lp>0) then
          name(1:lp) = path(np1+1:np2)
          if(name(lp:lp)==fndevider.and.filen(nf1+1:nf1+1)==fndevider) then
            lp=lp-1
          end if
          if(name(lp:lp)/=fndevider.and.filen(nf1+1:nf1+1)/=fndevider) then
            lp=lp+1
            name(lp:lp)=fndevider
          end if
          name(lp+1:lp+lf) = filen(nf1+1:nf2)
      else
          name(1:lf) = filen(nf1+1:nf2)
      end if

endsubroutine fulfname

!======================================================================
subroutine rdstd(path,fname,nfield,field,lu,nx,ny,nz,nxb,nxe,nyb,nye,nzb,nze,ierr)
!  nx,ny,nz - general dimesion of fild
!  nxb,nxe,nyb,nye,nzb,nze - grid coordinates of treat array subdomain
!      where index b denotes begin, and e - end

!     This subroutine fills (read) array FILD from unformatted deirect
!     file of DIOGIN standard
!     path           - path to file (i.g. 'F:\ARAB')
!     fname          - name of file (i.g.: 'taux.std')
!     nfild          - number of field in file (on t)
!     fild(nx,ny,nz) - field array
!     lu(nx,ny)    - ocean mask
!---------------------------------------------------------------------  
      character*(*) path, fname
      integer       nfield,nrecf,nx,ny,nz,i,j,k,ierr
      real          field(nx,ny,nz), lu(nx,ny)
      character(4096) namofile
      integer  nxe, nxb, nye, nyb, nzb, nze, l, kr

!  definition full file name
      call fulfname(namofile,path,fname,ierr)
      if (ierr/=0) go to 100

!  check of correctness of grid coordinates of treat array part
      if(nxe>nx.or.nxb<1.or.nxb>nxe) then
          ierr=1
          goto 103
      end if
      if(nye>ny.or.nyb<1.or.nyb>nye) then
          ierr=2
          goto 103
      end if
      if(nze>nz.or.nzb<1.or.nzb>nze) then
          ierr=3
          goto 103
      end if
!      write(*,'(2x,a)')'open direct file for reading:'
!      write(*,'(2x,a)')  namofile
!      write(*,'(2x,a,i4,a,i4)')'recl lenth in words:', nx,' x',ny
      open(40,file=namofile,status='old',access='direct',form='unformatted',recl=(nxe-nxb+1)*(nye-nyb+1)*lrecl,err=101)

      nrecf=(nfield-1)*(nze-nzb+1)     !initial number of record

      l=0
      kr=0
      do k=nzb,nze
      kr=kr+1
      read(40,rec=nrecf+kr,err=102) ((field(i,j,k),i=nxb,nxe),j=nyb,nye)

!  filling undefinite points by zero instead undef
          do j=nyb,nye
            do i=nxb,nxe

            if(abs(lu(i,j))<0.5) then
               l=l+1
               field(i,j,k)=0.0
            end if

            enddo
         enddo

      end do

      close(40)

      l=(nxe-nxb+1)*(nye-nyb+1)-l/kr
      write(*,'(1x,a,a)') 'input data from ',namofile(1:len_trim (namofile))
      write(*,'(7x,a,i7,a,i7,a,i8,a)') 'dimension of field =',nxe-nxb+1,' *',nye-nyb+1, ' (',l,'-ocean points)'
!     write(*,'(2x,a)')'close direct file:'
!     write(*,'(2x,a)')  namofile
      return
100   write(*,'(2x,a)')'error in full name of file for reading: '
      write(*,'(2x,a)') namofile(1:len_trim(namofile))  
      stop
101   write(*,'(2x,a)')'error in open file for reading: '
      write(*,'(2x,a)') namofile(1:len_trim(namofile)) 
      stop
102   write(*,'(2x,a)')'error in reading from file: '
      write(*,'(2x,a)') namofile(1:len_trim(namofile)) 
      write(*,'(18h error in reading ,i3,6h level)') k
      stop
103   write(*,'(2x,a)')'error in reading from file: '
      write(*,'(2x,a)') namofile(1:len_trim(namofile)) 
      write(*,'(2x,a,i3,a)')'error in grid diapason of ', ierr,' - coordinate'
      stop
endsubroutine rdstd

!======================================================================
subroutine wdstd(path,fname,nfield,field,lu,nx,ny,nz,nxb,nxe,nyb,nye,nzb,nze,ierr)
!  nx,ny,nz - general dimesion of fild
!  nxb,nxe,nyb,nye,nzb,nze - grid coordinates of treat array subdomain
!      where index b denotes begin, and e - end

!     This subroutine fills (WRITE) array FILD to unformatted deirect
!     file of DIOGIN standard
!     path           - path to file (i.g. 'F:\ARAB')
!     fname          - name of file (i.g.: 'taux.std')
!     nfild          - number of field in file (on t)
!     fild(nx,ny,nz) - field array
!     lu(nx,ny)    - ocean mask
!     ierr - error information. if ierr ne 0 in input, then no printing
!---------------------------------------------------------------------
      character*(*) path, fname
      integer       nfield,nrecf,nx,ny,nz,i,j,k,ierr
      real          field(nx,ny,nz), lu(nx,ny)
      character(4096) namofile
      integer      nxe, nxb, nye, nyb, nzb, nze, l, kr, lprint

    if(ierr==0) then
      lprint=1
    else
      lprint=0    
    end if

    ierr=0

!  definition full file name
      call fulfname(namofile,path,fname,ierr)
      if (ierr/=0) go to 100

!  check of correctness of grid coordinates of treat array part
      if(nxe>nx.or.nxb<1.or.nxb>nxe) then
          ierr=1
          goto 103
      end if
      if(nye>ny.or.nyb<1.or.nyb>nye) then
          ierr=2
          goto 103
      end if
      if(nze>nz.or.nzb<1.or.nzb>nze) then
          ierr=3
          goto 103
      end if
!      write(*,'(2x,a)')'open direct file for writing:'
!      write(*,'(2x,a)')  namofile
!      write(*,'(2x,a,i4,a,i4)')'recl lenth in words:', nx,' x',ny
      open(40,file=namofile,status='unknown',access='direct', form='unformatted',recl=(nxe-nxb+1)*(nye-nyb+1)*lrecl,err=101)

      nrecf=(nfield-1)*(nze-nzb+1)     !initial number of record

      l   =0
      ierr=0
      kr=0
      do k=nzb,nze
         kr=kr+1
!  fulling undefinite points by 0ver insted zero
          do j=nyb,nye
            do i=nxb,nxe

            if(abs(lu(i,j))<0.5) then
               l=l+1
               field(i,j,k)=undef
            end if

            enddo
         enddo

!  writing on the file
       write(40,rec=nrecf+kr,err=102)((field(i,j,k),i=nxb,nxe),j=nyb,nye)

!  fulling undefinite points by zero insted undef
          do j=nyb,nye
            do i=nxb,nxe

            if(abs(lu(i,j))<0.5) then
               ierr=ierr+1
               field(i,j,k)=0.0
            end if

            enddo
         enddo
      end do


      l=(nxe-nxb+1)*(nye-nyb+1)-l/kr
    
    if(lprint==1) then
! print information on terminal     
      write(*,'(1x,a,a)')  'output data to ',namofile(1:len_trim (namofile))
      write(*,'(8x,a,i7,a,i7,a,i8,a)') 'dimension of field =',nxe-nxb+1,' *',nye-nyb+1, ' (',l,'-ocean points)'
      end if
     
      close(40)
!      write(*,'(2x,a)')'close direct file:'
!      write(*,'(2x,a)') namofile

      ierr=(nxe-nxb+1)*(nye-nyb+1)-ierr/kr-l
      if (ierr/=0) then
            write(*,'(2x,a)')  namofile
            write(*,'(i7,a)') ierr, 'errors in number of ocean horizontal grid points.'
      endif

      return
100   write(*,'(2x,a)')'error in full name of file for writing: '
      write(*,'(2x,a)') namofile(1:len_trim(namofile)) 
      stop
101   write(*,'(2x,a)')'error in open file for writing: '
      write(*,'(2x,a)') namofile(1:len_trim(namofile)) 
      stop
102   write(*,'(2x,a)')'error in writing on file: '
      write(*,'(2x,a)') namofile(1:len_trim(namofile)) 
      write(*,'(18h error in writing ,i3,6h level)') k
      stop
103   write(*,'(2x,a)')'error in writing to file: '
      write(*,'(2x,a)') namofile(1:len_trim(namofile)) 
      write(*,'(2x,a,i3,a)')'error in grid diapason of ',  ierr,' - coordinate'
      stop
endsubroutine wdstd

!======================================================================
subroutine rdstd8(path,fname,nfield,field,lu,nx,ny,nz,nxb,nxe,nyb,nye,nzb,nze,ierr)
!  nx,ny,nz - general dimesion of fild
!  nxb,nxe,nyb,nye,nzb,nze - grid coordinates of treat array subdomain
!      where index b denotes begin, and e - end

!     This subroutine fills (read) array FILD from unformatted deirect
!     file of DIOGIN standard
!     path           - path to file (i.g. 'F:\ARAB')
!     fname          - name of file (i.g.: 'taux.std')
!     nfild          - number of field in file (on t)
!     fild(nx,ny,nz) - field array
!     lu(nx,ny)    - ocean mask
!---------------------------------------------------------------------
      character*(*) path, fname
      integer       nfield,nrecf,nx,ny,nz,i,j,k,ierr
      real*8        field(nx,ny,nz)
      real            lu(nx,ny)
      character(4096) namofile
      integer  nxe, nxb, nye, nyb, nzb, nze, l, kr

!  definition full file name
      call fulfname(namofile,path,fname,ierr)
      if (ierr/=0) go to 100

!  check of correctness of grid coordinates of treat array part
      if(nxe>nx.or.nxb<1.or.nxb>nxe) then
          ierr=1
          goto 103
      end if
      if(nye>ny.or.nyb<1.or.nyb>nye) then
          ierr=2
          goto 103
      end if
      if(nze>nz.or.nzb<1.or.nzb>nze) then
          ierr=3
          goto 103
      end if
!      write(*,'(2x,a)')'open direct file for reading:'
!      write(*,'(2x,a)')  namofile
!      write(*,'(2x,a,i4,a,i4)')'recl lenth in words:', nx,' x',ny
      open(40,file=namofile,status='old',access='direct', form='unformatted',recl=(nxe-nxb+1)*(nye-nyb+1)*lrecl*2,err=101)

      nrecf=(nfield-1)*(nze-nzb+1)     !initial number of record

      l=0
      kr=0
      do k=nzb,nze
      kr=kr+1
      read(40,rec=nrecf+kr,err=102) ((field(i,j,k),i=nxb,nxe),j=nyb,nye)

!  fulling undefinite points by zero insted undef
          do j=nyb,nye
            do i=nxb,nxe

            if(abs(lu(i,j))<0.5) then
               l=l+1
               field(i,j,k)=0.0
            end if

            enddo
         enddo

      end do

      close(40)

      l=(nxe-nxb+1)*(nye-nyb+1)-l/kr
      write(*,'(1x,a,a)') 'input data from ',namofile(1:len_trim (namofile))
      write(*,'(7x,a,i7,a,i7,a,i8,a)') 'dimension of field =',nxe-nxb+1,' *',nye-nyb+1, ' (',l,'-ocean points)'
!     write(*,'(2x,a)')'close direct file:'
!     write(*,'(2x,a)')  namofile
      return
100   write(*,'(2x,a)')'error in full name of file for reading: '
      write(*,'(2x,a)') namofile(1:len_trim(namofile)) 
      stop
101   write(*,'(2x,a)')'error in open file for reading: '
      write(*,'(2x,a)') namofile(1:len_trim(namofile)) 
      stop
102   write(*,'(2x,a)')'error in reading from file: '
      write(*,'(2x,a)') namofile(1:len_trim(namofile)) 
      write(*,'(18h error in reading ,i3,6h level)') k
      stop
103   write(*,'(2x,a)')'error in reading from file: '
      write(*,'(2x,a)') namofile(1:len_trim(namofile)) 
      write(*,'(2x,a,i3,a)')'error in grid diapason of ',  ierr,' - coordinate'
      stop
endsubroutine rdstd8

!======================================================================
subroutine wdstd8(path,fname,nfield,field,lu,nx,ny,nz, nxb,nxe,nyb,nye,nzb,nze,ierr)
!  nx,ny,nz - general dimesion of fild
!  nxb,nxe,nyb,nye,nzb,nze - grid coordinates of treat array subdomain
!      where index b denotes begin, and e - end

!     This subroutine fills (WRITE) array FILD to unformatted deirect
!     file of DIOGIN standard
!     path           - path to file (i.g. 'F:\ARAB')
!     fname          - name of file (i.g.: 'taux.std')
!     nfield          - number of field in file (on t)
!     field(nx,ny,nz) - field array
!     lu(nx,ny)    - ocean mask
!---------------------------------------------------------------------
      character*(*) path, fname
      integer       nfield,nrecf,nx,ny,nz,i,j,k,ierr
      real*8        field(nx,ny,nz)
      real          lu(nx,ny)
      character(4096) namofile
      integer      nxe, nxb, nye, nyb, nzb, nze, l, kr

!  definition full file name
      call fulfname(namofile,path,fname,ierr)
      if (ierr/=0) go to 100

!  check of correctness of grid coordinates of treat array part
      if(nxe>nx.or.nxb<1.or.nxb>nxe) then
          ierr=1
          goto 103
      end if
      if(nye>ny.or.nyb<1.or.nyb>nye) then
          ierr=2
          goto 103
      end if
      if(nze>nz.or.nzb<1.or.nzb>nze) then
          ierr=3
          goto 103
      end if
!      write(*,'(2x,a)')'open direct file for writing:'
!      write(*,'(2x,a)')  namofile
!      write(*,'(2x,a,i4,a,i4)')'recl lenth in words:', nx,' x',ny
      open(40,file=namofile,status='unknown',access='direct', form='unformatted',recl=(nxe-nxb+1)*(nye-nyb+1)*lrecl*2,err=101)

      nrecf=(nfield-1)*(nze-nzb+1)     !initial number of record

      l   =0
      ierr=0
      kr=0
      do k=nzb,nze
         kr=kr+1
!  fulling undefinite points by 0ver insted zero
          do j=nyb,nye
            do i=nxb,nxe

            if(abs(lu(i,j))<0.5) then
               l=l+1
               field(i,j,k)=dble(undef)
            end if

            enddo
         enddo

!  writing on the file
       write(40,rec=nrecf+kr,err=102)((field(i,j,k),i=nxb,nxe),j=nyb,nye)

!  fulling undefinite points by zero insted undef
          do j=nyb,nye
            do i=nxb,nxe

            if(abs(lu(i,j))<0.5) then
               ierr=ierr+1
               field(i,j,k)=0.0
            end if

            enddo
         enddo
      end do


      l=(nxe-nxb+1)*(nye-nyb+1)-l/kr
      write(*,'(1x,a,a)')'output data to ',namofile(1:len_trim (namofile))
      write(*,'(8x,a,i7,a,i7,a,i8,a)')'dimension of field =',nxe-nxb+1,' *',nye-nyb+1, ' (',l,'-ocean points)'
      close(40)
!      write(*,'(2x,a)')'close direct file:'
!      write(*,'(2x,a)') namofile

      ierr=(nxe-nxb+1)*(nye-nyb+1)-ierr/kr-l
      if (ierr/=0) then
            write(*,'(2x,a)')  namofile
            write(*,'(i7,a)') ierr,'errors in number of ocean horizontal grid points.'
      endif

      return
100   write(*,'(2x,a)')'error in full name of file for writing: '
      write(*,'(2x,a)') namofile(1:len_trim(namofile)) 
      stop
101   write(*,'(2x,a)')'error in open file for writing: '
      write(*,'(2x,a)') namofile(1:len_trim(namofile)) 
      stop
102   write(*,'(2x,a)')'error in writing on file: '
      write(*,'(2x,a)') namofile(1:len_trim(namofile)) 
      write(*,'(18h error in writing ,i3,6h level)') k
      stop
103   write(*,'(2x,a)')'error in writing to file: '
      write(*,'(2x,a)') namofile(1:len_trim(namofile)) 
      write(*,'(2x,a,i3,a)')'error in grid diapason of ', ierr,' - coordinate'
      stop
endsubroutine wdstd8

!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
!---------------------- Simple MPI subroutines ---------------------------------!
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
subroutine prdstd2D(path, fname, nfild, fild, lu, nx, ny, nxb, nxe, nyb, nye, ierr)

    character*(*) :: path, fname
    integer :: nx, ny
    integer :: nfild, i, j, k, kfull, ierr
    type(block2D_real4), dimension(:), pointer :: fild
    type(block2D_real4), dimension(:), pointer :: lu
    character(4096) :: namofile
    integer :: nxe, nxb, nye, nyb
    integer :: hfile
    integer(kind=mpi_offset_kind) :: disp
    integer :: tsubarr, totsize
    integer :: sizes2(2), locsizes2(2), offset2(2)

    !  definition full file name
    call fulfname(namofile,path,fname,ierr)
    if (ierr.ne.0) go to 100
    !  check of correctness of grid coordinates of treat array part
    if (nxe.gt.nx .or. nxb.lt.1 .or. nxb.gt.nxe) then
        ierr=1
        goto 103
    end if
    if (nye.gt.ny .or. nyb.lt.1 .or. nyb.gt.nye) then
        ierr=2
        goto 103
    end if

    disp = (nxe-nxb+1)*(nye-nyb+1)*int(lmpirecl, mpi_offset_kind)*(nfild-1)
    !locsizes = (/nx_end - nx_start + 1, ny_end - ny_start + 1, nze - nzb + 1/)
    !totsize = locsizes(1)*locsizes(2)*locsizes(3)

    call mpi_file_open(cart_comm, namofile, mpi_mode_rdonly, mpi_info_null, hfile, ierr)
    if (ierr .ne. mpi_success) goto 101    

    do kfull = 1, bcount_max
        if (kfull <= bcount) then
            k = kfull
            call set_block(k)
        else 
            k = 1
            call set_block(1)
        endif

        offset2 = (/nx_start - nxb, ny_start - nyb/)
        locsizes2 = (/nx_end - nx_start + 1, ny_end - ny_start + 1/)
        sizes2 = (/nxe - nxb + 1, nye - nyb + 1/)

        totsize = locsizes2(1)*locsizes2(2)
        
        !print *, sizes2, locsizes2, offset2
        call mpi_type_create_subarray(2, sizes2, locsizes2, offset2,        &
                                    mpi_order_fortran, mpi_real, tsubarr, ierr)
        
        call mpi_type_commit(tsubarr, ierr)

        call mpi_file_set_view(hfile, disp, mpi_real, tsubarr, "native", mpi_info_null, ierr)

        call mpi_file_read_all(hfile, fild(k)%vals(nx_start:nx_end, ny_start:ny_end),  &
                               totsize, mpi_real, mpi_status_ignore, ierr)
        if (ierr .ne. mpi_success) goto 102

        call mpi_type_free(tsubarr, ierr)
    enddo 
    call mpi_file_close(hfile, ierr)

    !  filling undefinite points by zero instead undef
    do k = 1, bcount
        call set_block(k)
        do j = ny_start-2, ny_end+2
            do i = nx_start-2, nx_end+2
                if (abs(lu(k)%vals(i,j)) < 0.5) then
                    fild(k)%vals(i, j) = 0.0
                end if
            enddo
        enddo
    enddo

    return

100   write(*,'(a,i5)') 'on rank ',rank
    write(*,'(2x,a)')'error in full name of file for reading: '
    write(*,'(2x,a)') namofile(1:len_trim(namofile))
    call mpi_abort(cart_comm,5,ierr)
    stop
101   write(*,'(a,i5)') 'on rank ',rank
    write(*,'(2x,a)')'error in open file for reading: '
    write(*,'(2x,a)') namofile(1:len_trim(namofile))
    call mpi_abort(cart_comm,5,ierr)
    stop
102   write(*,'(a,i5)') 'on rank ',rank
    write(*,'(2x,a)')'error in reading from file: '
    write(*,'(2x,a)') namofile(1:len_trim(namofile))
    write(*,'(18h error in reading ,i3,6h level)') k
    call mpi_abort(cart_comm,5,ierr)
    stop
103   write(*,'(a,i5)') 'on rank ',rank
    write(*,'(2x,a)')'error in reading from file: '
    write(*,'(2x,a)') namofile(1:len_trim(namofile))
    write(*,'(2x,a,i3,a)')'error in grid diapason of ', ierr,' - coordinate'
    call mpi_abort(cart_comm,5,ierr)
    stop
end subroutine prdstd2D

subroutine pwdstd2D(path, fname, nfild, fild, lu, nx, ny, nxb, nxe, nyb, nye, ierr)
    !  nx,ny,nz - general dimesion of fild
    !  nxb,nxe,nyb,nye,nzb,nze - grid coordinates of treat array subdomain
    !      where index b denotes begin, and e - end
    !     this subroutine fills (write) array fild to unformatted deirect
    !     file of diogin standard
    !     path           - path to file (i.g. 'f:\arab')
    !     fname          - name of file (i.g.: 'taux.std')
    !     nfild          - number of field in file (on t)
    !     fild(nx,ny,nz) - field array
    !     lu(nx,ny)    - ocean mask
    !---------------------------------------------------------------------
    character*(*) :: path, fname
    integer :: hfile
    integer(kind=mpi_offset_kind) :: disp
    integer :: nfild, nrecf, nx, ny, i, j, k, kfull, ierr
    type(block2D_real4), dimension(:), pointer :: fild
    type(block2D_real4), dimension(:), pointer :: lu
    character*4096 :: namofile
    integer  :: nxe, nxb, nye, nyb
    integer :: tsubarr, totsize
    integer :: sizes2(2), locsizes2(2), offset2(2)
    
    !  definition full file name
    call fulfname(namofile,path,fname,ierr)
    !      write(*,*) "writing to ", trim(namofile)
    if (ierr.ne.0) go to 100
    ! check of correctness of grid coordinates of treat array part
    if (nxe.gt.nx .or. nxb.lt.1 .or. nxb.gt.nxe) then
          ierr=1
          goto 103
    end if
    if (nye.gt.ny .or. nyb.lt.1 .or. nyb.gt.nye) then
          ierr=2
          goto 103
    end if
    
    disp = (nxe-nxb+1)*(nye-nyb+1)*int(lmpirecl, mpi_offset_kind)*(nfild-1)
    !locsizes = (/nx_end - nx_start + 1, ny_end - ny_start + 1, nze - nzb + 1/)
    !totsize = locsizes(1)*locsizes(2)*locsizes(3)

    if (rank == 0) then
        open(40,file=namofile,status='unknown',access='direct', form='unformatted',recl=(nxe-nxb+1)*(nye-nyb+1)*lrecl,err=101)
        ! initial number of record
        nrecf = nfild  
        ! writing on the file
        write(40,rec=nrecf,err=102)(( undef, i=nxb,nxe), j=nyb,nye )
        close(40)
    endif

    call mpi_file_open(cart_comm, namofile, ior(mpi_mode_wronly,mpi_mode_create), mpi_info_null, hfile, ierr)
    if (ierr .ne. mpi_success) goto 101

    do kfull = 1, bcount_max
        if (kfull <= bcount) then
            k = kfull
            call set_block(k)
        else 
            k = 1
            call set_block(1)
        endif
        
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
                    if (abs(lu(k)%vals(i,j)) .lt. 0.5) then
                        fild(k)%vals(i, j) = undef
                    end if
            enddo
        enddo
        
        call mpi_file_write_all(hfile, fild(k)%vals(nx_start:nx_end, ny_start:ny_end),  &
                                totsize, mpi_real, mpi_status_ignore, ierr)
        if (ierr.ne.mpi_success) goto 102
        
        call mpi_type_free(tsubarr, ierr)
    enddo

    call mpi_file_close(hfile, ierr)

    do k = 1, bcount
        call set_block(k)
        do j = ny_start, ny_end
            do i = nx_start, nx_end
                    if (abs(lu(k)%vals(i,j)) .lt. 0.5) then
                        fild(k)%vals(i, j) = 0.0
                    end if
            enddo
        enddo
    enddo

    return
    
100   write(*,'(2x,a)')'error in full name of file for writing: '
    write(*,'(2x,a)') namofile(1:len_trim(namofile))
    call mpi_abort(cart_comm,-1,ierr)
    stop
101   write(*,'(2x,a)')'error in open file for writing: '
    write(*,'(2x,a)') namofile(1:len_trim(namofile))
    call mpi_abort(cart_comm,-1,ierr)
    stop
102   write(*,'(2x,a)')'error in writing on file: '
    write(*,'(2x,a)') namofile(1:len_trim(namofile))
    write(*,'(18h error in writing ,i3,6h level)') k
    call mpi_abort(cart_comm,-1,ierr)
    stop
103   write(*,'(2x,a)')'error in writing to file: '
    write(*,'(2x,a)') namofile(1:len_trim(namofile))
    write(*,'(2x,a,i3,a)')'error in grid diapason of ', ierr,' - coordinate'
    call mpi_abort(cart_comm,-1,ierr)
    stop
end subroutine pwdstd2D

subroutine prdstd3D(path, fname, nfild, fild, lu, nx, ny, nz, nxb, nxe, nyb, nye, nzb, nze, ierr)

    character*(*) :: path, fname
    integer :: nx, ny, nz
    integer :: nfild, i, j, k, kfull, ierr
    type(block3D_real4), dimension(:), pointer :: fild
    type(block2D_real4), dimension(:), pointer :: lu
    character(4096) :: namofile
    integer :: nxe, nxb, nye, nyb, nzb, nze
    integer :: hfile
    integer(kind=mpi_offset_kind) :: disp
    integer :: tsubarr, totsize
    integer :: sizes3(3), locsizes3(3), offset3(3)

    !  definition full file name
    call fulfname(namofile,path,fname,ierr)
    if (ierr.ne.0) go to 100
    !  check of correctness of grid coordinates of treat array part
    if (nxe.gt.nx .or. nxb.lt.1 .or. nxb.gt.nxe) then
        ierr=1
        goto 103
    end if
    if (nye.gt.ny .or. nyb.lt.1 .or. nyb.gt.nye) then
        ierr=2
        goto 103
    end if
    if (nze.gt.nz .or. nzb.lt.1 .or. nzb.gt.nze) then
        ierr=3
        goto 103
    end if

    disp = (nxe-nxb+1)*(nye-nyb+1)*int(lmpirecl, mpi_offset_kind)*(nfild-1)*(nze-nzb+1)
    !locsizes = (/nx_end - nx_start + 1, ny_end - ny_start + 1, nze - nzb + 1/)
    !totsize = locsizes(1)*locsizes(2)*locsizes(3)

    call mpi_file_open(cart_comm, namofile, mpi_mode_rdonly, mpi_info_null, hfile, ierr)
    if (ierr .ne. mpi_success) goto 101    

    do kfull = 1, bcount_max
        if (kfull <= bcount) then
            k = kfull
            call set_block(k)
        else 
            k = 1
            call set_block(1)
        endif
        offset3 = (/nx_start - nxb, ny_start - nyb, 0/)
        locsizes3 = (/nx_end - nx_start + 1, ny_end - ny_start + 1, nze - nzb + 1/)
        sizes3 = (/nxe - nxb + 1, nye - nyb + 1, nze - nzb + 1/)
        totsize = locsizes3(1)*locsizes3(2)*locsizes3(3)
        
        call mpi_type_create_subarray(3, sizes3, locsizes3, offset3, mpi_order_fortran, mpi_real, tsubarr, ierr)
        call mpi_type_commit(tsubarr, ierr)

        call mpi_file_set_view(hfile, disp, mpi_real, tsubarr, "native", mpi_info_null, ierr)

        call mpi_file_read_all(hfile, fild(k)%vals(nx_start:nx_end, ny_start:ny_end, nzb:nze),  &
                               totsize, mpi_real, mpi_status_ignore, ierr)
        if (ierr .ne. mpi_success) goto 102

        call mpi_type_free(tsubarr, ierr)
    enddo 
    call mpi_file_close(hfile, ierr)

    !  filling undefinite points by zero instead undef
    do k = 1, bcount
        call set_block(k)
        do j = ny_start-2, ny_end+2
            do i = nx_start-2, nx_end+2
                if (abs(lu(k)%vals(i,j)) < 0.5) then
                    fild(k)%vals(i, j, nzb:nze) = 0.0
                end if
            enddo
        enddo
    enddo

    return

100   write(*,'(a,i5)') 'on rank ',rank
    write(*,'(2x,a)')'error in full name of file for reading: '
    write(*,'(2x,a)') namofile(1:len_trim(namofile))
    call mpi_abort(cart_comm,5,ierr)
    stop
101   write(*,'(a,i5)') 'on rank ',rank
    write(*,'(2x,a)')'error in open file for reading: '
    write(*,'(2x,a)') namofile(1:len_trim(namofile))
    call mpi_abort(cart_comm,5,ierr)
    stop
102   write(*,'(a,i5)') 'on rank ',rank
    write(*,'(2x,a)')'error in reading from file: '
    write(*,'(2x,a)') namofile(1:len_trim(namofile))
    write(*,'(18h error in reading ,i3,6h level)') k
    call mpi_abort(cart_comm,5,ierr)
    stop
103   write(*,'(a,i5)') 'on rank ',rank
    write(*,'(2x,a)')'error in reading from file: '
    write(*,'(2x,a)') namofile(1:len_trim(namofile))
    write(*,'(2x,a,i3,a)')'error in grid diapason of ', ierr,' - coordinate'
    call mpi_abort(cart_comm,5,ierr)
    stop
end subroutine prdstd3D

subroutine pwdstd3D(path, fname, nfild, fild, lu, nx, ny, nz, nxb, nxe, nyb, nye, nzb, nze, ierr)
    !  nx,ny,nz - general dimesion of fild
    !  nxb,nxe,nyb,nye,nzb,nze - grid coordinates of treat array subdomain
    !      where index b denotes begin, and e - end
    !     this subroutine fills (write) array fild to unformatted deirect
    !     file of diogin standard
    !     path           - path to file (i.g. 'f:\arab')
    !     fname          - name of file (i.g.: 'taux.std')
    !     nfild          - number of field in file (on t)
    !     fild(nx,ny,nz) - field array
    !     lu(nx,ny)    - ocean mask
    !---------------------------------------------------------------------
    character*(*) :: path, fname
    integer :: hfile
    integer(kind=mpi_offset_kind) :: disp
    integer :: nfild, nrecf, nx, ny, nz, i, j, k, kfull, ierr, kr
    type(block3D_real4), dimension(:), pointer :: fild
    type(block2D_real4), dimension(:), pointer :: lu
    character*4096 :: namofile
    integer  :: nxe, nxb, nye, nyb, nzb, nze
    integer :: tsubarr, totsize
    integer :: sizes3(3), locsizes3(3), offset3(3)
    
    !  definition full file name
    call fulfname(namofile,path,fname,ierr)
    !      write(*,*) "writing to ", trim(namofile)
    if (ierr.ne.0) go to 100
    ! check of correctness of grid coordinates of treat array part
    if (nxe.gt.nx .or. nxb.lt.1 .or. nxb.gt.nxe) then
          ierr=1
          goto 103
    end if
    if (nye.gt.ny .or. nyb.lt.1 .or. nyb.gt.nye) then
          ierr=2
          goto 103
    end if
    if (nze.gt.nz .or. nzb.lt.1 .or. nzb.gt.nze) then
        ierr=3
        goto 103
    end if
    
    disp = (nxe-nxb+1)*(nye-nyb+1)*int(lmpirecl, mpi_offset_kind)*(nfild-1)*(nze-nzb+1)
    !locsizes = (/nx_end - nx_start + 1, ny_end - ny_start + 1, nze - nzb + 1/)
    !totsize = locsizes(1)*locsizes(2)*locsizes(3)

    if (rank == 0) then
        open(40,file=namofile,status='unknown',access='direct', form='unformatted',recl=(nxe-nxb+1)*(nye-nyb+1)*lrecl,err=101)
        ! initial number of record
        nrecf = (nfild - 1)*(nze - nzb + 1)
        ! writing on the file 
        kr = 0
        do k = nzb, nze
            kr = kr + 1
            write(40, rec=nrecf+kr, err=102)(( undef, i=nxb,nxe), j=nyb,nye)
        enddo
        close(40)
    endif

    call mpi_file_open(cart_comm, namofile, ior(mpi_mode_wronly,mpi_mode_create), mpi_info_null, hfile, ierr)
    if (ierr .ne. mpi_success) goto 101

    do kfull = 1, bcount_max
        if (kfull <= bcount) then
            k = kfull
            call set_block(k)
        else 
            k = 1
            call set_block(1)
        endif
        
        offset3 = (/nx_start - nxb, ny_start - nyb, 0/)
        locsizes3 = (/nx_end - nx_start + 1, ny_end - ny_start + 1, nze - nzb + 1/)
        sizes3 = (/nxe - nxb + 1, nye - nyb + 1, nze - nzb + 1/)
        totsize = locsizes3(1)*locsizes3(2)*locsizes3(3)
        call mpi_type_create_subarray(3, sizes3, locsizes3, offset3, mpi_order_fortran, mpi_real, tsubarr, ierr)
        call mpi_type_commit(tsubarr, ierr)

        call mpi_file_set_view(hfile, disp, mpi_real, tsubarr, "native", mpi_info_null, ierr)
        
        do j = ny_start, ny_end
            do i = nx_start, nx_end
                    if (abs(lu(k)%vals(i,j)) .lt. 0.5) then
                        fild(k)%vals(i, j, nzb:nze) = undef
                    end if
            enddo
        enddo
        
        call mpi_file_write_all(hfile, fild(k)%vals(nx_start:nx_end, ny_start:ny_end, nzb:nze),  &
                                totsize, mpi_real, mpi_status_ignore, ierr)
        if (ierr.ne.mpi_success) goto 102
        
        call mpi_type_free(tsubarr, ierr)
    enddo

    call mpi_file_close(hfile, ierr)

    do k = 1, bcount
        call set_block(k)
        do j = ny_start, ny_end
            do i = nx_start, nx_end
                    if (abs(lu(k)%vals(i,j)) .lt. 0.5) then
                        fild(k)%vals(i, j, nzb:nze) = 0.0
                    end if
            enddo
        enddo
    enddo

    return
    
100   write(*,'(2x,a)')'error in full name of file for writing: '
    write(*,'(2x,a)') namofile(1:len_trim(namofile))
    call mpi_abort(cart_comm,-1,ierr)
    stop
101   write(*,'(2x,a)')'error in open file for writing: '
    write(*,'(2x,a)') namofile(1:len_trim(namofile))
    call mpi_abort(cart_comm,-1,ierr)
    stop
102   write(*,'(2x,a)')'error in writing on file: '
    write(*,'(2x,a)') namofile(1:len_trim(namofile))
    write(*,'(18h error in writing ,i3,6h level)') k
    call mpi_abort(cart_comm,-1,ierr)
    stop
103   write(*,'(2x,a)')'error in writing to file: '
    write(*,'(2x,a)') namofile(1:len_trim(namofile))
    write(*,'(2x,a,i3,a)')'error in grid diapason of ', ierr,' - coordinate'
    call mpi_abort(cart_comm,-1,ierr)
    stop
end subroutine pwdstd3D

endmodule iodata_routes
