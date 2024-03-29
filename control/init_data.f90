#include "macros/kernel_macros.fi"
#include "macros/mpp_macros.fi"

module init_data_module
    ! Initialize data in model

    use mpp_module
    use data_types_module, only: data2D_real4_type, data2D_real8_type
    use decomposition_module, only: domain_type, domain => domain_data
    use ocean_module, only: ocean_type, ocean_data
    use grid_module, only: grid_type, grid_data, grid_global_type, grid_global_data
    use io_module, only: read_data
    use gridcon_module, only: gridcon
    use basinpar_module, only: basinpar
    use config_sw_module, only: use_tracers, tracer_num
    use mpp_sync_module, only: sync, sync_test
    use errors_module, only: abort_model

    implicit none
    save
    private

    public :: init_ocean_data
    public :: init_grid_data
    public :: init_device_data

contains

    subroutine init_ocean_data()
        use config_sw_module, only: ssh_init_file_name, lvisc_2
        use kernel_interface_module
        use shallow_water_interface_module, only: envoke_gaussian_elimination, envoke_hh_init_kernel, envoke_hh_init_sync
        use config_basinpar_module, only: nx, ny
        
        type(data2D_real4_type) :: tmp_data
        procedure(envoke_empty_kernel), pointer :: sub_kernel
        procedure(envoke_empty_sync), pointer :: sub_sync
        type(kernel_parameters_type) :: kernel_parameters
        integer :: k, ierr

        !debug
        !call sync_test(domain, ocean_data%ssh)
        !call sync_test(domain, grid_data%lu)
        !return

        if (ssh_init_file_name .eq. 'none') then
            if (mpp_is_master()) print *, "INIT: WARNING: SSH init file is none. Set gaussian elimination as init ssh level."
            call envoke_gaussian_elimination(ocean_data%ssh, 1.0d0, nx / 2, ny / 2)
        else
            call tmp_data%init(domain)
            call read_data(domain, 'INIT/', ssh_init_file_name, 1, tmp_data, grid_data%lu, ierr)
            call ocean_data%ssh%copy_from_real4(domain, tmp_data)
            call sync(domain, ocean_data%ssh)
            call tmp_data%clear(domain)
        endif

        call ocean_data%sshn%copy_from(domain, ocean_data%ssh)
        call ocean_data%sshp%copy_from(domain, ocean_data%ssh)

        call kernel_parameters%clear()
        sub_kernel => envoke_hh_init_kernel
        sub_sync   => envoke_hh_init_sync
        call envoke(sub_kernel, sub_sync, kernel_parameters)

        ! Zero velocity
        ! U
        call ocean_data%ubrtr%fill(domain, 0.0d0)
        call ocean_data%ubrtrn%fill(domain, 0.0d0)
        call ocean_data%ubrtrp%fill(domain, 0.0d0)
        ! V
        call ocean_data%vbrtr%fill(domain, 0.0d0)
        call ocean_data%vbrtrn%fill(domain, 0.0d0)
        call ocean_data%vbrtrp%fill(domain, 0.0d0)

        ! init
        call ocean_data%mu%fill(domain, lvisc_2)
        call ocean_data%mu%fill(domain, 0.0d0)

        ! Tracers
        if (use_tracers > 0) then
            if (mpp_is_master()) print *, "INIT: Tracers init as gaussian elimination. Tracers num ", tracer_num
            do k = 1, tracer_num        
                call envoke_gaussian_elimination(ocean_data%ff1(k), 0.5d0, nx / 2, ny / 2)
                call ocean_data%ff1n(k)%copy_from(domain, ocean_data%ff1(k))
                call ocean_data%ff1p(k)%copy_from(domain, ocean_data%ff1(k))
            enddo
            ! Zero fluxes
            call ocean_data%flux_x%fill(domain, 0.0d0)
            call ocean_data%flux_y%fill(domain, 0.0d0)
        endif

        call mpp_sync_output()

    end subroutine init_ocean_data

    subroutine init_grid_data()
        use config_basinpar_module, only: bottom_topography_file_name

        type(data2D_real4_type) :: tmp_data
        integer :: ierr

        ! area mask initialization
        call gridcon

        ! setting vertical t-,w- grid levels
        !call vgrid
        
        ! define grid geographical coordinates, steps and coriolis parameters
        call basinpar

        ! Read bottom topograhy (real4 binary file)
        if (bottom_topography_file_name .eq. 'none') then
            if (mpp_is_master()) print *, 'INIT: WARNING: Bottom topography is none. Set 100m as ocean depth.'
            call grid_data%hhq_rest%fill(domain, 100.0d0)
        else
            call tmp_data%init(domain)
            call read_data(domain, ' ', bottom_topography_file_name, 1, tmp_data, grid_data%lu, ierr)
            call grid_data%hhq_rest%copy_from_real4(domain, tmp_data)
            call sync(domain, grid_data%hhq_rest)
            call tmp_data%clear(domain)
        endif

        call mpp_sync_output()

    end subroutine init_grid_data

    subroutine init_device_data
#ifdef _GPU_MODE_
        use depth_gpu_module, only: depth_load_constant_mem
        use velssh_sw_gpu_module, only: velssh_sw_load_constant_mem
        use tracer_gpu_module, only: tracer_load_constant_mem

        ! Sync init data between host and device
        call grid_data%sync_host_device(domain, .true.)
        call ocean_data%sync_host_device(domain, .true.)

        ! Prepare constant id device memory
        call depth_load_constant_mem
        call velssh_sw_load_constant_mem
        if (use_tracers) then
            call tracer_load_constant_mem
        endif
#endif
        return
    end subroutine init_device_data

end module init_data_module