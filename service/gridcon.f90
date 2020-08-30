module gridcon_module
    ! Grid consruction module

    use mpp_module
    use decomposition_module, only: domain_type
    use grid_module, only: grid_type, grid_global_type
    use mpp_sync_module, only: sync
    use errors_module, only: abort_model

#include "macros/mpp_macros.fi"

    implicit none
    save
    private

    public :: gridcon

contains

    subroutine gridcon(domain, grid_global_data, grid_data)
        ! Grid consruction module by temperature mask.
        ! subroutin for construction pass boundary, velosity and bottom masks
        ! using temperature mask in diogin standart
        use config_basinpar_module, only: nx, ny, mask_file_name

        type(domain_type), intent(in) :: domain
        type(grid_global_type), intent(in) :: grid_global_data
        type(grid_type), intent(inout) :: grid_data
        
        integer :: m, n, k, ierr

        ! conversion integer diogin mask to real model mask
        do k = 1, domain%bcount
          associate(bnd_x1 => domain%bbnd_x1(k), bnd_x2 => domain%bbnd_x2(k),  &
                    bnd_y1 => domain%bbnd_y1(k), bnd_y2 => domain%bbnd_y2(k),  &
                    mask => grid_global_data%mask,       &
                    lu => grid_data%lu%block(k)%field,   &
                    lu1 =>grid_data%lu1%block(k)%field)

          do n = bnd_y1, bnd_y2
            do m = bnd_x1, bnd_x2
              if (mask(m,n) == 0) then
                lu(m,n) = 1.0
              endif
            end do
          end do
          
          lu1 = 1.0
          
          end associate
        enddo

        call sync(domain, grid_data%lu)

        ! forming mask for depth grid points
        ! forming luh from lu, which have land neibours in luh.
        ! constructing array luh for relief hh.
               
        if (mpp_is_master()) then
          write(*,*) 'Construction of H-grid masks: '
          write(*,*) 'LUH (includes boundary) and LUU (does not include boundary)'
        endif 
        
        do k = 1, domain%bcount
          associate(bnd_x1 => domain%bbnd_x1(k), bnd_x2 => domain%bbnd_x2(k),  &
                    bnd_y1 => domain%bbnd_y1(k), bnd_y2 => domain%bbnd_y2(k),  &
                    lu => grid_data%lu%block(k)%field,   &
                    luh =>grid_data%luh%block(k)%field,  &
                    luu =>grid_data%luu%block(k)%field)

          !do n = ny_start-1,ny_end
          !do m = nx_start-1,nx_end
          do n = bnd_y1, bnd_y2-1
            do m = bnd_x1, bnd_x2-1
              if (lu(m,n)+lu(m+1,n)+lu(m,n+1)+lu(m+1,n+1)>0.5)  then
                luh(m,n)=1.0
              endif

              if(lu(m,n)*lu(m+1,n)*lu(m,n+1)*lu(m+1,n+1)>0.5)  then
                luu(m,n)=1.0
              endif
            enddo
          enddo

          end associate
        enddo

        call sync(domain, grid_data%luh)
        call sync(domain, grid_data%luu)
        
        if (mpp_is_master()) then
          write(*,*) 'Construction of U- and V-grid masks: '
          write(*,*) 'LCU and LCV (do not include boundary) and LLU and LLV (include boundary)'
        endif
        
        do k = 1, domain%bcount
          associate(bnd_x1 => domain%bbnd_x1(k), bnd_x2 => domain%bbnd_x2(k),  &
                    bnd_y1 => domain%bbnd_y1(k), bnd_y2 => domain%bbnd_y2(k),  &
                    lu => grid_data%lu%block(k)%field,   &
                    llu =>grid_data%llu%block(k)%field,  &
                    llv =>grid_data%llv%block(k)%field,  &
                    lcu =>grid_data%lcu%block(k)%field,  &
                    lcv =>grid_data%lcv%block(k)%field)

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

          end associate
        enddo

        call sync(domain, grid_data%lcu)
        call sync(domain, grid_data%llu)
        call sync(domain, grid_data%lcv)
        call sync(domain, grid_data%llv)
        
    end subroutine gridcon

end module gridcon_module