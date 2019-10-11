module grid_module
    ! Grid data
    !** Обернуть все в типы, зачем нужен target? см. E3SM

    use decomposition_module, only: domain_type
    use data_types_module, only: data2D_real4_type

    implicit none
    save
    private
    
    type, public :: grid_type
        type(data2D_real4_type) :: lu,   &   ! mask of t-grid
                                   hhq,  &   ! ocean depth on lu  (t-points)
                                   hhu,  &   ! ocean depth on lcu (u-points)
                                   hhv       ! ocean depth on lcv (v-points)
        type(data2D_real4_type) :: dxt, dyt  ! horizontal grid steps between t-points (in meters)
        type(data2D_real4_type) :: sqt       ! grid area in T-points
    contains
        procedure, public :: init
        procedure, public :: clear
    end type grid_type

    type(grid_type), public, target :: grid_data
    
contains

    subroutine init(this, domain)
        ! Initialization of grid data
        class(grid_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        
        call this%lu%init(domain)
        call this%hhq%init(domain)
        call this%hhu%init(domain)
        call this%hhv%init(domain)

        call this%dxt%init(domain)
        call this%dyt%init(domain)
        call this%sqt%init(domain)
    end subroutine

    subroutine clear(this, domain)
        ! Clear grid data
        class(grid_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain

        call this%lu%clear(domain)
        call this%hhq%clear(domain)
        call this%hhu%clear(domain)
        call this%hhv%clear(domain)

        call this%dxt%clear(domain)
        call this%dyt%clear(domain)
        call this%sqt%clear(domain)
    end subroutine

end module grid_module