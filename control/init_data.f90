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
    use errors_module, only: abort_model

    implicit none
    save
    private

    public :: init_ocean_data
    public :: init_grid_data

contains

    subroutine init_ocean_data(domain, grid_data, ocean_data)
        use config_sw_module, only: ssh_init_file_name
        use shallow_water_interface_module, only: envoke_hh_init_kernel, envoke_check_ssh_err_kernel

        type(domain_type), intent(in) :: domain
        type(grid_type), intent(inout) :: grid_data
        type(ocean_type), intent(inout) :: ocean_data
        
        type(data2D_real4_type) :: tmp_data
        integer :: k, ierr

        call tmp_data%init(domain)
        call read_data(domain, 'INIT/', ssh_init_file_name, 1, tmp_data, grid_data%lu, ierr)
        call ocean_data%ssh%copy_from_real4(domain, tmp_data)
        call sync(domain, ocean_data%ssh)
        call tmp_data%clear(domain)

        call ocean_data%sshn%copy_from(domain, ocean_data%ssh)
        call ocean_data%sshp%copy_from(domain, ocean_data%ssh)

        call envoke_hh_init_kernel(domain, grid_data, ocean_data%ssh, ocean_data%sshp)

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
        call grid_data%hhq_rest%copy_from_real4(domain, tmp_data)
        call sync(domain, grid_data%hhq_rest)
        call tmp_data%clear(domain)

    end subroutine init_grid_data

end module init_data_module