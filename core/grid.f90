module grid_module
    ! Grid data

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use mpp_module, only: mpp_rank, mpp_count, mpp_cart_comm, mpp_size, mpp_coord, mpp_period
    use decomposition_module, only: domain_type
    use data_types_module, only: data1D_real4_type, data1D_real8_type, data2D_real4_type, data2D_real8_type, data3D_real8_type
    use errors_module, only: abort_model, check_error

    implicit none
    save
    private

    type, public :: grid_global_type
        integer, allocatable :: mask(:, :)
    contains
        procedure, public :: init => init_global_grid_type
        procedure, public :: clear => clear_global_grid_type
    end type grid_global_type
    
    type, public :: grid_type
        type(data2D_real4_type) :: lu,        &  ! Mask of t-grid
                                   lu1,       &  ! Mask of t-grid (1 everywhere)
                                   luu,       &  ! Mask of h-grid (0 on boundary)
                                   luh,       &  ! Mask of h-grid (1 on boundary)
                                   lcu,       &  ! Mask of u-grid (0 on boundary)
                                   lcv,       &  ! Mask of v-grid (0 on boundary)
                                   llu,       &  ! Mask of u-grid (1 on boundary)
                                   llv           ! Mask of v-grid (1 on boundary)

        ! Temporairly using of real8 instead of real4 for depths
        !type(data2D_real4_type) :: hhh,       &  ! Ocean depth on luh (h-points)
        type(data2D_real8_type) :: hhh,       &  ! Ocean depth on luh (h-points)
                                   hhh_n,     &  ! Ocean depth on luh (h-points) at previous step
                                   hhh_p,     &  ! Ocean depth on luh (h-points) at pre-previous step
                                   hhq_rest,  &  ! Ocean depth on lu  (t-points) at rest state
                                   hhu_rest,  &  ! Ocean depth on lcu (u-points) at rest state
                                   hhv_rest,  &  ! Ocean depth on lcv (v-points) at rest state
                                   hhq,       &  ! Ocean depth on lu  (t-points) 
                                   hhq_n,     &  ! Ocean depth on lu  (t-points) at previous step
                                   hhq_p,     &  ! Ocean depth on lu  (t-points) at pre-previous step
                                   hhts,      &
                                   hhts_n,    &
                                   hhu,       &  ! Ocean depth on lcu (u-points)
                                   hhu_n,     &  ! Ocean depth on lcu (u-points) at previous step
                                   hhu_p,     &  ! Ocean depth on lcu (u-points) at pre-previous step
                                   hhv,       &  ! Ocean depth on lcv (v-points)
                                   hhv_n,     &  ! Ocean depth on lcv (v-points) at previous step
                                   hhv_p         ! Ocean depth on lcv (v-points) at pre-previous step

        type(data2D_real4_type) :: rlh_s,     &  ! Main (sin) coriolis parameter on h-points
                                   rlh_sqh,   &  ! Main (sin) coriolis parameter on h-points by area
                                   rlh_c         ! 2-nd (cos) coriolis parameter on h-points

        type(data1D_real4_type) :: z, zw,  &   ! Vertical sigma-levels (t-points and w-points)
                                   hzt, dz     ! Steps between t-levels and w-levels

        type(data2D_real4_type) :: dxt, dyt,  &  ! Horizontal grid steps between   t-points (in meters)
                                   dx , dy ,  &  ! Horizontal grid steps between u,v-points (in meters)
                                   dxh, dyh,  &  ! Horizontal grid steps between   h-points (in meters)
                                   dxb, dyb      ! Horizontal grid steps between v,u-points (in meters)

        type(data2D_real4_type) :: sqt,  &  ! Grid area in T-points
                                   squ,  &  ! Grid area in U-points
                                   sqv,  &  ! Grid area in V-points
                                   sqh      ! Grid area in H-points

        type(data1D_real8_type) :: xt, yt,  &  ! Horizontal t-grid            x- and y-coordinates (in degrees)
                                   xu, yv      ! Horizontal u-grid and v-grid x- and y-coordinates (in degrees)


        type(data2D_real8_type) :: geo_lon_t,  &  ! Geographical longitudes of T-points
                                   geo_lat_t,  &  ! Geographical latitudes  of T-points
                                   geo_lon_u,  &  ! Geographical longitudes of U-points
                                   geo_lat_u,  &  ! Geographical latitudes  of U-points
                                   geo_lon_v,  &  ! Geographical longitudes of V-points
                                   geo_lat_v,  &  ! Geographical latitudes  of V-points
                                   geo_lon_h,  &  ! Geographical longitudes of H-points
                                   geo_lat_h      ! Geographical latitudes  of H-points

        type(data3D_real8_type) :: rotvec_coeff  ! Cos and sin of angles between coordinate lines
    contains
        procedure, public :: init => init_grid_type
        procedure, public :: clear => clear_grid_type
    end type grid_type

!------------------------------------------------------------------------------

    type(grid_global_type), public, target :: grid_global_data
    type(grid_type), public, target :: grid_data
    
contains

    subroutine init_global_grid_type(this)
        use config_basinpar_module, only: nx, ny

        class(grid_global_type), intent(inout) :: this

        allocate(this%mask(nx, ny))
    end subroutine

    subroutine clear_global_grid_type(this)
        class(grid_global_type), intent(inout) :: this

        deallocate(this%mask)
    end subroutine

    subroutine init_grid_type(this, domain)
        use config_basinpar_module, only: nz

        ! Initialization of grid data
        class(grid_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        
        call this%lu %init(domain) 
        call this%lu1%init(domain)
        call this%luu%init(domain)
        call this%luh%init(domain)
        call this%lcu%init(domain)
        call this%lcv%init(domain)
        call this%llu%init(domain)
        call this%llv%init(domain)

        call this%hhh     %init(domain)
        call this%hhh_n   %init(domain)
        call this%hhq_rest%init(domain)
        call this%hhu_rest%init(domain)
        call this%hhv_rest%init(domain)
        call this%hhq     %init(domain)
        call this%hhq_n   %init(domain)
        call this%hhts    %init(domain)
        call this%hhts_n  %init(domain)
        call this%hhu     %init(domain)
        call this%hhu_n   %init(domain)
        call this%hhv     %init(domain)
        call this%hhv_n   %init(domain)
        call this%rlh_s   %init(domain)
        call this%rlh_sqh %init(domain)
        call this%rlh_c   %init(domain)

        call this%z  %init(domain, 1, domain%nz)
        call this%zw %init(domain, 1, domain%nz + 1)
        call this%hzt%init(domain, 1, domain%nz + 1)
        call this%dz %init(domain, 1, domain%nz)

        call this%dxt%init(domain)
        call this%dyt%init(domain)
        call this%dx %init(domain) 
        call this%dy %init(domain) 
        call this%dxh%init(domain)
        call this%dyh%init(domain)
        call this%dxb%init(domain)
        call this%dyb%init(domain)

        call this%sqt%init(domain)
        call this%squ%init(domain)
        call this%sqv%init(domain)
        call this%sqh%init(domain)

        call this%xt%init_from_domain(domain, .true.)  ! nx
        call this%yt%init_from_domain(domain, .false.) ! ny
        call this%xu%init_from_domain(domain, .true.)  ! nx
        call this%yv%init_from_domain(domain, .false.) ! ny

        call this%geo_lon_t%init(domain)
        call this%geo_lat_t%init(domain)
        call this%geo_lon_u%init(domain)
        call this%geo_lat_u%init(domain)
        call this%geo_lon_v%init(domain)
        call this%geo_lat_v%init(domain)
        call this%geo_lon_h%init(domain)
        call this%geo_lat_h%init(domain)

        call this%rotvec_coeff%init(domain, 1, 4)

        ! Data from pre-previous time step
        call this%hhh_p%init(domain)  ! Ocean depth on luh (h-points) at pre-previous step
        call this%hhq_p%init(domain)  ! Ocean depth on lu  (t-points) at pre-previous step
        call this%hhu_p%init(domain)  ! Ocean depth on lcu (u-points) at pre-previous step
        call this%hhv_p%init(domain)  ! Ocean depth on lcv (v-points) at pre-previous step
        

    end subroutine

    subroutine clear_grid_type(this, domain)
        ! Clear grid data
        class(grid_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        
        call this%lu %clear(domain) 
        call this%lu1%clear(domain)
        call this%luu%clear(domain)
        call this%luh%clear(domain)
        call this%lcu%clear(domain)
        call this%lcv%clear(domain)
        call this%llu%clear(domain)
        call this%llv%clear(domain)

        call this%hhh     %clear(domain)
        call this%hhh_n   %clear(domain)
        call this%hhq_rest%clear(domain)
        call this%hhu_rest%clear(domain)
        call this%hhv_rest%clear(domain)
        call this%hhq     %clear(domain)
        call this%hhq_n   %clear(domain)
        call this%hhts    %clear(domain)
        call this%hhts_n  %clear(domain)
        call this%hhu     %clear(domain)
        call this%hhu_n   %clear(domain)
        call this%hhv     %clear(domain)
        call this%hhv_n   %clear(domain)
        call this%rlh_s   %clear(domain)
        call this%rlh_sqh %clear(domain)
        call this%rlh_c   %clear(domain)

        call this%z  %clear(domain)
        call this%zw %clear(domain)
        call this%hzt%clear(domain)
        call this%dz %clear(domain)

        call this%dxt%clear(domain)
        call this%dyt%clear(domain)
        call this%dx %clear(domain) 
        call this%dy %clear(domain) 
        call this%dxh%clear(domain)
        call this%dyh%clear(domain)
        call this%dxb%clear(domain)
        call this%dyb%clear(domain)

        call this%sqt%clear(domain)
        call this%squ%clear(domain)
        call this%sqv%clear(domain)
        call this%sqh%clear(domain)

        call this%xt%clear(domain)
        call this%yt%clear(domain)
        call this%xu%clear(domain)
        call this%yv%clear(domain)

        call this%geo_lon_t%clear(domain)
        call this%geo_lat_t%clear(domain)
        call this%geo_lon_u%clear(domain)
        call this%geo_lat_u%clear(domain)
        call this%geo_lon_v%clear(domain)
        call this%geo_lat_v%clear(domain)
        call this%geo_lon_h%clear(domain)
        call this%geo_lat_h%clear(domain)

        call this%rotvec_coeff%clear(domain)

        ! Data from pre-previous time step
        call this%hhh_p%clear(domain)  ! Ocean depth on luh (h-points) at pre-previous step
        call this%hhq_p%clear(domain)  ! Ocean depth on lu  (t-points) at pre-previous step
        call this%hhu_p%clear(domain)  ! Ocean depth on lcu (u-points) at pre-previous step
        call this%hhv_p%clear(domain)  ! Ocean depth on lcv (v-points) at pre-previous step

    end subroutine

end module grid_module