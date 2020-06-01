module iodata_routes
    
    use, intrinsic :: iso_fortran_env, only: file_storage_size

    implicit none
    save
    public


    !integer, parameter:: lrecl=1         !long of unique recl for Compaq, Intel fortran
    integer, parameter:: lrecl=4         !long of unique recl for PGI, gfortran, efc

    integer, parameter :: lmpirecl=4  !for mpi read/write
    real, parameter:: undef=-1.0e+32  !missing data value

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

endmodule iodata_routes
