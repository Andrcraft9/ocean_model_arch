#include "macros/mpp_macros.fi"

module ocean_module
    ! Ocean data

    use decomposition_module, only: domain_type
    use data_types_module, only: data2D_real8_type, data2D_real4_type

    implicit none
    save
    private

    type, public :: ocean_type
        ! Barotropic dynamics arrays
        type(data2D_real8_type) :: ssh,      &  !sea surface height (SSH) at current  time step [m]
                                   pgrx,     &  !pressure gradient x-component for RHS
                                   pgry,     &  !pressure gradient y-component for RHS
                                   ubrtr,    &  !barotropic velocity      zonal[m/s] at current time step
                                   vbrtr,    &  !barotropic velocity meridional[m/s] at current time step
                                   RHSx2d,   &  !x-component of external force(barotropic)
                                   RHSy2d       !y-component of external force(barotropic)
        
        ! Data from previous time step
        type(data2D_real8_type) :: sshn,     &
                                   sshp,     &  !sea surface height (SSH) at previous time step [m]
                                   ubrtrn,   &
                                   ubrtrp,   &  !barotropic velocity      zonal[m/s] at previous time step
                                   vbrtrn,   &
                                   vbrtrp       !barotropic velocity meridional[m/s] at previous time step

        type(data2D_real4_type) :: r_diss       !Rayleigh friction scale (1/s)

        type(data2D_real8_type) :: RHSx, RHSy, RHSx_dif, RHSy_dif, RHSx_adv, RHSy_adv
        type(data2D_real8_type) :: mu, str_t, str_s, vort

        ! Tracers
        type(data2D_real8_type) :: ff1, ff1p, ff1n

        ! Fluxes for computations
        type(data2D_real8_type) :: flux_x, flux_y
    contains
        procedure, public  :: init
        procedure, public  :: clear
#ifdef _GPU_MODE_
        procedure, public  :: sync_host_device
#endif
    end type ocean_type

!------------------------------------------------------------------------------

    type(ocean_type), public, target :: ocean_data

contains

    subroutine init(this, domain)
        ! Initialization of grid data
        class(ocean_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain

        call this%ssh%init(domain)
        call this%pgrx%init(domain)
        call this%pgry%init(domain)
        call this%ubrtr%init(domain)
        call this%vbrtr%init(domain)
        call this%RHSx2d%init(domain)
        call this%RHSy2d%init(domain)

        call this%sshn%init(domain)
        call this%sshp%init(domain)
        call this%ubrtrn%init(domain)
        call this%ubrtrp%init(domain)
        call this%vbrtrn%init(domain)
        call this%vbrtrp%init(domain)

        call this%r_diss%init(domain)

        call this%RHSx%init(domain) 
        call this%RHSy%init(domain) 
        call this%RHSx_dif%init(domain) 
        call this%RHSy_dif%init(domain) 
        call this%RHSx_adv%init(domain) 
        call this%RHSy_adv%init(domain)
        call this%mu%init(domain) 
        call this%str_t%init(domain) 
        call this%str_s%init(domain) 
        call this%vort%init(domain)

        call this%ff1%init(domain)
        call this%ff1p%init(domain)
        call this%ff1n%init(domain)

        call this%flux_x%init(domain)
        call this%flux_y%init(domain)
    end subroutine

    subroutine clear(this, domain)
        ! Clear grid data
        class(ocean_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain

        call this%ssh%clear(domain)
        call this%pgrx%clear(domain)
        call this%pgry%clear(domain)
        call this%ubrtr%clear(domain)
        call this%vbrtr%clear(domain)
        call this%RHSx2d%clear(domain)
        call this%RHSy2d%clear(domain)

        call this%sshn%clear(domain)
        call this%sshp%clear(domain)
        call this%ubrtrn%clear(domain)
        call this%ubrtrp%clear(domain)
        call this%vbrtrn%clear(domain)
        call this%vbrtrp%clear(domain)

        call this%r_diss%clear(domain)

        call this%RHSx%clear(domain) 
        call this%RHSy%clear(domain) 
        call this%RHSx_dif%clear(domain) 
        call this%RHSy_dif%clear(domain) 
        call this%RHSx_adv%clear(domain) 
        call this%RHSy_adv%clear(domain)
        call this%mu%clear(domain) 
        call this%str_t%clear(domain) 
        call this%str_s%clear(domain) 
        call this%vort%clear(domain)

        call this%ff1%clear(domain)
        call this%ff1p%clear(domain)
        call this%ff1n%clear(domain)

        call this%flux_x%clear(domain)
        call this%flux_y%clear(domain)
    end subroutine

#ifdef _GPU_MODE_
    subroutine sync_host_device(this, domain, is_htod)
        ! Initialization of grid data
        class(ocean_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        logical, intent(in) :: is_htod

        call this%ssh%sync_host_device(domain, is_htod)
        call this%pgrx%sync_host_device(domain, is_htod)
        call this%pgry%sync_host_device(domain, is_htod)
        call this%ubrtr%sync_host_device(domain, is_htod)
        call this%vbrtr%sync_host_device(domain, is_htod)
        call this%RHSx2d%sync_host_device(domain, is_htod)
        call this%RHSy2d%sync_host_device(domain, is_htod)

        call this%sshn%sync_host_device(domain, is_htod)
        call this%sshp%sync_host_device(domain, is_htod)
        call this%ubrtrn%sync_host_device(domain, is_htod)
        call this%ubrtrp%sync_host_device(domain, is_htod)
        call this%vbrtrn%sync_host_device(domain, is_htod)
        call this%vbrtrp%sync_host_device(domain, is_htod)

        call this%r_diss%sync_host_device(domain, is_htod)
        
        call this%RHSx%sync_host_device(domain, is_htod) 
        call this%RHSy%sync_host_device(domain, is_htod) 
        call this%RHSx_dif%sync_host_device(domain, is_htod) 
        call this%RHSy_dif%sync_host_device(domain, is_htod) 
        call this%RHSx_adv%sync_host_device(domain, is_htod) 
        call this%RHSy_adv%sync_host_device(domain, is_htod)
        call this%mu%sync_host_device(domain, is_htod) 
        call this%str_t%sync_host_device(domain, is_htod) 
        call this%str_s%sync_host_device(domain, is_htod) 
        call this%vort%sync_host_device(domain, is_htod)

        call this%ff1%sync_host_device(domain, is_htod)
        call this%ff1p%sync_host_device(domain, is_htod)
        call this%ff1n%sync_host_device(domain, is_htod)

        call this%flux_x%sync_host_device(domain, is_htod)
        call this%flux_y%sync_host_device(domain, is_htod)
    end subroutine
#endif

endmodule ocean_module
