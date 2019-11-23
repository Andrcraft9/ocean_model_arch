module ocean_module
    ! Ocean data

    use decomposition_module, only: domain_type
    use data_types_module, only: data2D_real8_type

    implicit none
    save
    private

    type, public :: ocean_type
        type(data2D_real8_type) :: ssh,    &  ! sea surface height (SSH) at current  time step [m] (internal mode)
                                   ubrtr,  &  ! barotropic velocity      zonal[m/s] at current time step (internal mode)
                                   vbrtr      ! barotropic velocity meridional[m/s] at current time step (internal mode)
        type(data2D_real8_type) :: div_btr
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
        call this%ubrtr%init(domain)
        call this%vbrtr%init(domain)

        call this%div_btr%init(domain)
    end subroutine

    subroutine clear(this, domain)
        ! Clear grid data
        class(ocean_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain

        call this%ssh%clear(domain)
        call this%ubrtr%clear(domain)
        call this%vbrtr%clear(domain)

        call this%div_btr%clear(domain)
    end subroutine

endmodule ocean_module
