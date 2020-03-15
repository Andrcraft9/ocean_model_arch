module ocean_module
    ! Ocean data

    use decomposition_module, only: domain_type
    use data_types_module, only: data2D_real8_type

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
    contains
        procedure, public  :: init
        procedure, public  :: clear
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
    end subroutine

endmodule ocean_module
