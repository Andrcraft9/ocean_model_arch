module init_data_module
#include "core/kernel_macros.fi"
    ! Initialize data in model

    use mpp_module, only: mpp_rank, mpp_count, mpp_cart_comm, mpp_size, mpp_coord, mpp_period
    use data_types_module, only: data2D_real4_type
    use decomposition_module, only: domain_type
    use ocean_module, only: ocean_type
    use grid_module, only: grid_type, grid_global_type
    use io_module, only: read_data
    use gridcon_module, only: gridcon
    use basinpar_module, only: basinpar
    use mpp_sync_module, only: sync 

    implicit none
    save
    private

    public :: init_ocean_data
    public :: init_grid_data

contains

    subroutine init_ocean_data(domain, ocean_data)
        type(domain_type), intent(in) :: domain
        type(ocean_type), intent(inout) :: ocean_data
        integer :: k
        
        do k = 1, domain%bcount
            associate(_associate_data_(ocean_data, ssh,   k),  &
                      _associate_data_(ocean_data, ubrtr, k),  &
                      _associate_data_(ocean_data, vbrtr, k))
            
              ssh = 0.0
              ubrtr = 0.0
              vbrtr = 0.0

            end associate
        enddo

    end subroutine init_ocean_data

    subroutine init_grid_data(domain, grid_global_data, grid_data)
        use config_basinpar_module, only: bottom_topography_file_name

        type(domain_type), intent(in) :: domain
        type(grid_global_type), intent(inout) :: grid_global_data
        type(grid_type), intent(inout) :: grid_data

        type(data2D_real4_type) :: tmp_data
        integer :: ierr

        ! area mask initialization
        call gridcon(domain, grid_global_data, grid_data)

        ! setting vertical t-,w- grid levels
        !call vgrid
        
        ! define grid geographical coordinates, steps and coriolis parameters
        call basinpar(domain, grid_data)

        ! Read bottom topograhy (real4 binary file)
        call tmp_data%init(domain)
        call read_data(domain, ' ', bottom_topography_file_name, 1, tmp_data, grid_data%lu, ierr)
        call grid_data%hhq_rest%copy_real4(domain, tmp_data)
        call sync(domain, grid_data%hhq_rest)
        call tmp_data%clear(domain)

    end subroutine init_grid_data

end module init_data_module