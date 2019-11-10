module init_data_module
    ! Grid data

    use basinpar_module, only: dxst, dyst
    use decomposition_module, only: domain_type
    use ocean_module, only: ocean_type
    use grid_module, only: grid_type
    use io_module

    implicit none
    save
    private

    public :: init_ocean_data
    public :: init_grid_data

contains

    subroutine init_ocean_data(domain, ocean_data)
        ! Subroutine description
        type(domain_type), intent(in) :: domain
        type(ocean_type), intent(inout) :: ocean_data
        integer :: k

        call read_2d_real8(domain, ocean_data%ssh)
        call read_2d_real8(domain, ocean_data%ubrtr)
        call read_2d_real8(domain, ocean_data%vbrtr)
        
        do k = 1, domain%bcount
            associate(ssh => ocean_data%ssh%block(k)%field,      &
                      ubrtr => ocean_data%ubrtr%block(k)%field,  &
                      vbrtr => ocean_data%vbrtr%block(k)%field)
            
                ! some computational on data
                ! ...
            end associate
        enddo

    end subroutine init_ocean_data

    subroutine init_grid_data(domain, grid_data)
        ! Subroutine description
        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        integer :: k

        call read_2d_real4(domain, grid_data%hhq)
        call read_2d_real4(domain, grid_data%hhu)
        call read_2d_real4(domain, grid_data%hhv)

        do k = 1, domain%bcount
            associate(lu => grid_data%lu%block(k)%field,    &
                      hhq => grid_data%hhq%block(k)%field,  &
                      hhu => grid_data%hhu%block(k)%field,  &
                      hhv => grid_data%hhv%block(k)%field,  &
                      dxt => grid_data%dxt%block(k)%field,  &
                      dyt => grid_data%dyt%block(k)%field,  &
                      sqt => grid_data%sqt%block(k)%field)

                lu = 1.0
                dxt = dxst
                dyt = dyst
                sqt = dxst*dyst

            end associate
        enddo
    end subroutine init_grid_data

    subroutine init_mask
        ! Subroutine description
        implicit none
        
    end subroutine init_mask

    ! grid consruction module by temperature mask.
subroutine gridcon(ftemask)
use basin_grid
use ocalg_routes
include 'reclen.fi'
! subroutin for construction pass boundary, velosity and bottom masks
! using temperature mask in diogin standart
!--------------------------------------------------------------------
character*(*) ftemask
character frmt*16,comment*80
! temporary integer indexes
integer m, n, ierr

write(frmt,1000) nx
1000  format('(',i9,'i1)')

! reading mask from:
if (rank == 0) then
  open (11,file=ftemask,status='old',recl=nx*lrecl)
    read (11,  '(a)') comment(1:min(80,nx))
    write(*,'(1x,a)') comment
    do n=ny,1,-1
      read(11,frmt,end=99) (lbasins(m,n),m=1,nx)
    enddo
  close(11)
endif
call mpi_bcast(lbasins, nx*ny, mpi_integer, 0, cart_comm, ierr)

! conversion integer diogin mask to real model mask
do n = bnd_y1, bnd_y2
  do m = bnd_x1, bnd_x2
    if (lbasins(m,n)==0) then
      lu(m,n)=1.0
    endif
  end do
end do

lu1=1.0

!  forming mask for depth grid points
!  forming luh from lu, which have land neibours in luh.
! constructing array luh for relief hh.
 
if(periodicity_x/=0) then
  call cyclize_x(lu,nx,ny,1,mmm,mm)
endif

if(periodicity_y/=0) then
  call cyclize_y(lu,nx,ny,1,nnn,nn)
endif
       
if (rank == 0) then
  write(*,*) 'Construction of H-grid masks: '
  write(*,*) 'LUH (includes boundary) and LUU (does not include boundary)'
endif 

!do n = ny_start-1,ny_end
!do m = nx_start-1,nx_end
do n = bnd_y1, bnd_y2-1
  do m = bnd_x1, bnd_x2-1
    if (lu(m,n)+lu(m+1,n)+lu(m,n+1)+lu(m+1,n+1)>0.5)  then
      luh(m,n)=1.0
    endif
  enddo
enddo

!do n=ny_start-1,ny_end
!do m=nx_start-1,nx_end
do n = bnd_y1, bnd_y2-1
  do m = bnd_x1, bnd_x2-1
    if(lu(m,n)*lu(m+1,n)*lu(m,n+1)*lu(m+1,n+1)>0.5)  then
      luu(m,n)=1.0
    endif
  enddo
enddo

if(periodicity_x/=0) then
  call cyclize_x(luh,nx,ny,1,mmm,mm)
  call cyclize_x(luu,nx,ny,1,mmm,mm)
endif

if(periodicity_y/=0) then
  call cyclize_y(luh,nx,ny,1,nnn,nn)
  call cyclize_y(luu,nx,ny,1,nnn,nn)
endif

if (rank == 0) then
  write(*,*) 'Construction of U- and V-grid masks: '
  write(*,*) 'LCU and LCV (do not include boundary) and LLU and LLV (include boundary)'
endif

!do n=ny_start-1,ny_end
!do m=nx_start-1,nx_end
do n = bnd_y1, bnd_y2-1
  do m = bnd_x1, bnd_x2-1    
    
    if (lu(m,n)+lu(m+1,n)>0.5) then
      llu(m,n)=1.0
    endif

    if (lu(m,n)+lu(m,n+1)>0.5) then
      llv(m,n)=1.0
    endif

    if (lu(m,n)*lu(m+1,n)>0.5) then
      lcu(m,n)=1.0
    endif

    if (lu(m,n)*lu(m,n+1)>0.5) then
      lcv(m,n)=1.0
    endif
          
  enddo
enddo

if (periodicity_x/=0) then
  if (rank == 0) write(*,*)'  set periodicity to u-grid mask(lcu,llu).'
  call cyclize_x(lcu,nx,ny,1,mmm,mm)
  call cyclize_x(llu,nx,ny,1,mmm,mm)
  if (rank == 0) write(*,*)'  set periodicity to v-grid mask(lcv,llv).'
  call cyclize_x(lcv,nx,ny,1,mmm,mm)
  call cyclize_x(llv,nx,ny,1,mmm,mm)
endif

if (periodicity_y/=0) then
  call cyclize_y(lcu,nx,ny,1,nnn,nn)
  call cyclize_y(llu,nx,ny,1,nnn,nn)
  call cyclize_y(lcv,nx,ny,1,nnn,nn)
  call cyclize_y(llv,nx,ny,1,nnn,nn)
endif

return
99    write(*,*)'  error in reading file ',ftemask(1:len_trim(ftemask))
call mpi_abort(cart_comm, 1, ierr)
stop 1
endsubroutine gridcon

end module init_data_module