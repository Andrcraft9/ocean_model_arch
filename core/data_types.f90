#include "macros/mpp_macros.fi"

module data_types_module
    ! Module description

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use decomposition_module, only: domain_type
#ifdef _GPU_MODE_
    use cudafor
#endif

    implicit none
    save
    private

!------------------------------------------------------------------------------
    ! Data types for Parallel System layer
    ! 1D blocks
#ifdef _GPU_MODE_
    type, public :: block1D_real4_type
        real(wp4), pinned, allocatable :: field(:)
        real(wp4), allocatable, device :: field_device(:)
    end type block1D_real4_type
#else
    type, public :: block1D_real4_type
        real(wp4), allocatable :: field(:)
    end type block1D_real4_type
#endif

#ifdef _GPU_MODE_
    type, public :: block1D_real8_type
        real(wp8), pinned, allocatable :: field(:)
        real(wp8), allocatable, device :: field_device(:)
    end type block1D_real8_type
#else
    type, public :: block1D_real8_type
        real(wp8), allocatable :: field(:)
    end type block1D_real8_type
#endif

    ! 2D blocks
#ifdef _GPU_MODE_
    type, public :: block2D_real4_type
        real(wp4), pinned, allocatable :: field(:, :)
        real(wp4), allocatable, device :: field_device(:, :)
    end type block2D_real4_type
#else
    type, public :: block2D_real4_type
        real(wp4), allocatable :: field(:, :)
    end type block2D_real4_type
#endif

#ifdef _GPU_MODE_
    type, public :: block2D_real8_type
        real(wp8), pinned, allocatable :: field(:, :)
        real(wp8), allocatable, device :: field_device(:, :)
    end type block2D_real8_type
#else
    type, public :: block2D_real8_type
        real(wp8), allocatable :: field(:, :)
    end type block2D_real8_type
#endif

    ! 3D blocks
#ifdef _GPU_MODE_
    type, public :: block3D_real4_type
        real(wp4), pinned, allocatable :: field(:, :, :)
        real(wp4), allocatable, device :: field_device(:, :, :)
    end type block3D_real4_type
#else
    type, public :: block3D_real4_type
        real(wp4), allocatable :: field(:, :, :)
    end type block3D_real4_type
#endif

#ifdef _GPU_MODE_
    type, public :: block3D_real8_type
        real(wp8), pinned, allocatable :: field(:, :, :)
        real(wp8), allocatable, device :: field_device(:, :, :)
    end type block3D_real8_type
#else
    type, public :: block3D_real8_type
        real(wp8), allocatable :: field(:, :, :)
    end type block3D_real8_type
#endif

!------------------------------------------------------------------------------
    ! Data types for Algorithm layer
    ! 1D data
    type, public :: data1D_real4_type
        type(block1D_real4_type), pointer :: block(:)
    contains
        procedure, public :: init => init_data1D_real4
        procedure, public :: init_from_domain => init_from_domain_data1D_real4
        procedure, public :: clear => clear_data1D_real4
    end type data1D_real4_type

    type, public :: data1D_real8_type
        type(block1D_real8_type), pointer :: block(:)
    contains
        procedure, public :: init => init_data1D_real8
        procedure, public :: init_from_domain => init_from_domain_data1D_real8
        procedure, public :: clear => clear_data1D_real8
    end type data1D_real8_type

    ! 2D data
    type, public :: data2D_real4_type
        type(block2D_real4_type), pointer :: block(:)
    contains
        procedure, public :: init => init_data2D_real4
        procedure, public :: clear => clear_data2D_real4
        procedure, public :: copy_from => copy_data2D_real4_from_real4
        procedure, public :: copy_from_real8 => copy_data2D_real4_from_real8
        procedure, public :: fill => fill_data2D_real4
        ! debug
        procedure, public :: fill_debug => fill_debug_data2D_real4
    end type data2D_real4_type

    type, public :: data2D_real8_type
        type(block2D_real8_type), pointer :: block(:)
    contains
        procedure, public :: init => init_data2D_real8
        procedure, public :: clear => clear_data2D_real8
        procedure, public :: copy_from => copy_data2D_real8_from_real8
        procedure, public :: copy_from_real4 => copy_data2D_real8_from_real4
        procedure, public :: fill => fill_data2D_real8
        ! debug
        procedure, public :: fill_debug => fill_debug_data2D_real8
        procedure, public :: init_nans => init_nans_data2D_real8
    end type data2D_real8_type

    ! 3D data
    type, public :: data3D_real4_type
        type(block3D_real4_type), pointer :: block(:)
    contains
        procedure, public :: init => init_data3D_real4
        procedure, public :: clear => clear_data3D_real4
    end type data3D_real4_type

    type, public :: data3D_real8_type
        type(block3D_real8_type), pointer :: block(:)
    contains
        procedure, public :: init => init_data3D_real8
        procedure, public :: clear => clear_data3D_real8
    end type data3D_real8_type

    !interface copy_data2D_real4
    !    module procedure copy_data2D_real4_from_real4
    !    module procedure copy_data2D_real4_from_real8
    !end interface copy_data2D_real4
    !interface copy_data2D_real8
    !    module procedure copy_data2D_real8_from_real4
    !    module procedure copy_data2D_real8_from_real8
    !end interface copy_data2D_real8

contains

!------------------------------------------------------------------------------
    ! Init and clear for 1D data type
    ! real4 base type
    subroutine init_data1D_real4(this, domain, n_start, n_end)
        class(data1D_real4_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        integer, intent(in) :: n_start, n_end
        integer :: k

        allocate(this%block(domain%bcount))

        _OMP_BLOCKS_PARALLEL_BEGIN_
        do k = 1, domain%bcount
            allocate(this%block(k)%field(n_start : n_end))
            this%block(k)%field = 0.0
        enddo
        _OMP_BLOCKS_PARALLEL_END_

#ifdef _GPU_MODE_
        do k = 1, domain%bcount
            allocate(this%block(k)%field_device(n_start : n_end))
            this%block(k)%field_device = 0.0
        enddo
#endif
    end subroutine
    
    subroutine init_from_domain_data1D_real4(this, domain, is_x_direction)
        class(data1D_real4_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        logical, intent(in) :: is_x_direction
        integer :: k

        allocate(this%block(domain%bcount))

        _OMP_BLOCKS_PARALLEL_BEGIN_
        do k = 1, domain%bcount
            if (is_x_direction) then
                allocate(this%block(k)%field(domain%bbnd_x1(k) : domain%bbnd_x2(k)))
            else
                allocate(this%block(k)%field(domain%bbnd_y1(k) : domain%bbnd_y2(k)))
            endif
            this%block(k)%field = 0.0
        enddo
        _OMP_BLOCKS_PARALLEL_END_

#ifdef _GPU_MODE_
        do k = 1, domain%bcount
            if (is_x_direction) then
                allocate(this%block(k)%field_device(domain%bbnd_x1(k) : domain%bbnd_x2(k)))
            else
                allocate(this%block(k)%field_device(domain%bbnd_y1(k) : domain%bbnd_y2(k)))
            endif
            this%block(k)%field_device = 0.0
        enddo
#endif
    end subroutine

    subroutine clear_data1D_real4(this, domain)
        class(data1D_real4_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        integer :: k

#ifdef _GPU_MODE_
        do k = 1, domain%bcount
            deallocate(this%block(k)%field_device)
        enddo
#endif

        do k = 1, domain%bcount
            deallocate(this%block(k)%field)
        enddo

        deallocate(this%block)
    end subroutine

    ! real8 base type
    subroutine init_data1D_real8(this, domain, n_start, n_end)
        class(data1D_real8_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        integer, intent(in) :: n_start, n_end
        integer :: k

        allocate(this%block(domain%bcount))

        _OMP_BLOCKS_PARALLEL_BEGIN_
        do k = 1, domain%bcount
            allocate(this%block(k)%field(n_start : n_end))
            this%block(k)%field = 0.0d0
        enddo
        _OMP_BLOCKS_PARALLEL_END_

#ifdef _GPU_MODE_
        do k = 1, domain%bcount
            allocate(this%block(k)%field_device(n_start : n_end))
            this%block(k)%field_device = 0.0d0
        enddo
#endif
    end subroutine
    
    subroutine init_from_domain_data1D_real8(this, domain, is_x_direction)
        class(data1D_real8_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        logical, intent(in) :: is_x_direction
        integer :: k

        allocate(this%block(domain%bcount))

        _OMP_BLOCKS_PARALLEL_BEGIN_
        do k = 1, domain%bcount
            if (is_x_direction) then
                allocate(this%block(k)%field(domain%bbnd_x1(k) : domain%bbnd_x2(k)))
            else
                allocate(this%block(k)%field(domain%bbnd_y1(k) : domain%bbnd_y2(k)))
            endif
            this%block(k)%field = 0.0d0
        enddo
        _OMP_BLOCKS_PARALLEL_END_

#ifdef _GPU_MODE_
        do k = 1, domain%bcount
            if (is_x_direction) then
                allocate(this%block(k)%field_device(domain%bbnd_x1(k) : domain%bbnd_x2(k)))
            else
                allocate(this%block(k)%field_device(domain%bbnd_y1(k) : domain%bbnd_y2(k)))
            endif
            this%block(k)%field = 0.0d0
        enddo
#endif
    end subroutine

    subroutine clear_data1D_real8(this, domain)
        class(data1D_real8_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        integer :: k

#ifdef _GPU_MODE_
        do k = 1, domain%bcount
            deallocate(this%block(k)%field_device)
        enddo
#endif

        do k = 1, domain%bcount
            deallocate(this%block(k)%field)
        enddo

        deallocate(this%block)
    end subroutine

!------------------------------------------------------------------------------
    ! Init and clear for 2D data type
    ! real4 base type
    subroutine init_data2D_real4(this, domain)
        class(data2D_real4_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        integer :: k, m, n

        allocate(this%block(domain%bcount))

        _OMP_BLOCKS_PARALLEL_BEGIN_
        do k = 1, domain%bcount
            allocate(this%block(k)%field(domain%bbnd_x1(k) : domain%bbnd_x2(k), domain%bbnd_y1(k) : domain%bbnd_y2(k)))
            do n = domain%bbnd_y1(k), domain%bbnd_y2(k)
                do m = domain%bbnd_x1(k), domain%bbnd_x2(k)
                    this%block(k)%field(m, n) = 0.0
                enddo
            enddo
        enddo
        _OMP_BLOCKS_PARALLEL_END_

#ifdef _GPU_MODE_
        do k = 1, domain%bcount
            allocate(this%block(k)%field_device(domain%bbnd_x1(k) : domain%bbnd_x2(k), domain%bbnd_y1(k) : domain%bbnd_y2(k)))
            this%block(k)%field_device = 0.0
        enddo
#endif
    end subroutine

    subroutine clear_data2D_real4(this, domain)
        class(data2D_real4_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        integer :: k

#ifdef _GPU_MODE_
        do k = 1, domain%bcount
            deallocate(this%block(k)%field_device)
        enddo
#endif

        do k = 1, domain%bcount
            deallocate(this%block(k)%field)
        enddo

        deallocate(this%block)
    end subroutine

    subroutine copy_data2D_real4_from_real4(this, domain, copy_data)
        class(data2D_real4_type), intent(inout) :: this
        type(data2D_real4_type), intent(in) :: copy_data
        type(domain_type), intent(in) :: domain
        integer :: k, m, n

        _OMP_BLOCKS_PARALLEL_BEGIN_
        do k = 1, domain%bcount
            do n = domain%bbnd_y1(k), domain%bbnd_y2(k)
                do m = domain%bbnd_x1(k), domain%bbnd_x2(k)
                    this%block(k)%field(m, n) = copy_data%block(k)%field(m, n)
                enddo
            enddo
        enddo
        _OMP_BLOCKS_PARALLEL_END_

#ifdef _GPU_MODE_
        ! Copy from host to device, because device to device may be unsupported
        do k = 1, domain%bcount
            this%block(k)%field_device = this%block(k)%field
        enddo
#endif
    end subroutine

    subroutine copy_data2D_real4_from_real8(this, domain, copy_data)
        class(data2D_real4_type), intent(inout) :: this
        type(data2D_real8_type), intent(in) :: copy_data
        type(domain_type), intent(in) :: domain
        integer :: k, m, n

        _OMP_BLOCKS_PARALLEL_BEGIN_
        do k = 1, domain%bcount
            do n = domain%bbnd_y1(k), domain%bbnd_y2(k)
                do m = domain%bbnd_x1(k), domain%bbnd_x2(k)
                    this%block(k)%field(m, n) = real(copy_data%block(k)%field(m, n), wp4)
                enddo
            enddo
        enddo
        _OMP_BLOCKS_PARALLEL_END_

#ifdef _GPU_MODE_
        ! Copy from host to device, because device to device may be unsupported
        do k = 1, domain%bcount
            this%block(k)%field_device = this%block(k)%field
        enddo
#endif
    end subroutine

    subroutine fill_data2D_real4(this, domain, val)
        class(data2D_real4_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        real(wp4), intent(in) :: val
        integer :: k, m, n

        _OMP_BLOCKS_PARALLEL_BEGIN_
        do k = 1, domain%bcount
            do n = domain%bbnd_y1(k), domain%bbnd_y2(k)
                do m = domain%bbnd_x1(k), domain%bbnd_x2(k)
                    this%block(k)%field(m, n) = val
                enddo
            enddo
        enddo
        _OMP_BLOCKS_PARALLEL_END_

#ifdef _GPU_MODE_
        do k = 1, domain%bcount
            this%block(k)%field_device = val
        enddo
#endif
    end subroutine

    subroutine fill_debug_data2D_real4(this, domain)
        class(data2D_real4_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        integer :: k, m, n

        _OMP_BLOCKS_PARALLEL_BEGIN_
        do k = 1, domain%bcount
            do n = domain%bbnd_y1(k), domain%bbnd_y2(k)
                do m = domain%bbnd_x1(k), domain%bbnd_x2(k)
                    this%block(k)%field(m, n) = real(k, wp4)
                enddo
            enddo
        enddo
        _OMP_BLOCKS_PARALLEL_END_

#ifdef _GPU_MODE_
        do k = 1, domain%bcount
            this%block(k)%field_device = real(k, wp4)
        enddo
#endif
    end subroutine

    ! real8 base type
    subroutine init_data2D_real8(this, domain)
        class(data2D_real8_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        integer :: k, m, n

        allocate(this%block(domain%bcount))

        _OMP_BLOCKS_PARALLEL_BEGIN_
        do k = 1, domain%bcount
            allocate(this%block(k)%field(domain%bbnd_x1(k) : domain%bbnd_x2(k), domain%bbnd_y1(k) : domain%bbnd_y2(k)))
            do n = domain%bbnd_y1(k), domain%bbnd_y2(k)
                do m = domain%bbnd_x1(k), domain%bbnd_x2(k)
                    this%block(k)%field(m, n) = 0.0d0
                enddo
            enddo
        enddo
        _OMP_BLOCKS_PARALLEL_END_

#ifdef _GPU_MODE_
        do k = 1, domain%bcount
            allocate(this%block(k)%field_device(domain%bbnd_x1(k) : domain%bbnd_x2(k), domain%bbnd_y1(k) : domain%bbnd_y2(k)))
            this%block(k)%field_device = 0.0d0
        enddo
#endif
    end subroutine

    subroutine init_nans_data2D_real8(this, domain)
        use mpp_module
        use, intrinsic :: iso_fortran_env
        use, intrinsic :: ieee_arithmetic

        class(data2D_real8_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        integer :: k

        integer :: m, n
        real(wp8) :: val_nan

        val_nan = ieee_value(val_nan, ieee_quiet_nan)

        allocate(this%block(domain%bcount))
        do k = 1, domain%bcount
            allocate(this%block(k)%field(domain%bbnd_x1(k) : domain%bbnd_x2(k), domain%bbnd_y1(k) : domain%bbnd_y2(k)))
            this%block(k)%field = val_nan

! May be useful for MPI debug:            
!            this%block(k)%field = 0.0
!            if (mpp_rank == 5) then
!                    associate(bnd_x1 => domain%bbnd_x1(k), bnd_x2 => domain%bbnd_x2(k),  &
!                            bnd_y1 => domain%bbnd_y1(k), bnd_y2 => domain%bbnd_y2(k),  &
!                            nx_start => domain%bnx_start(k), nx_end => domain%bnx_end(k),  &
!                            ny_start => domain%bny_start(k), ny_end => domain%bny_end(k))
!            
!                        do n = ny_start-1, ny_end+1
!                            do m = nx_start-1, nx_end+1
!                                if (m == nx_start - 1) then
!                                    this%block(k)%field(m, n) = val_nan
!                                endif
!                                if (m == nx_end + 1) then
!                                    this%block(k)%field(m, n) = val_nan
!                                endif
!                                if (n == ny_start - 1) then
!                                    this%block(k)%field(m, n) = val_nan
!                                endif
!                                if (n == ny_end + 1) then
!                                    this%block(k)%field(m, n) = val_nan
!                                endif
!                            enddo
!                        enddo
!                    end associate
!                endif
        enddo

#ifdef _GPU_MODE_
        do k = 1, domain%bcount
            this%block(k)%field_device = val_nan
        enddo
#endif
    end subroutine

    subroutine clear_data2D_real8(this, domain)
        class(data2D_real8_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        integer :: k
        
#ifdef _GPU_MODE_
        do k = 1, domain%bcount
            deallocate(this%block(k)%field_device)
        enddo
#endif

        do k = 1, domain%bcount
            deallocate(this%block(k)%field)
        enddo
        
        deallocate(this%block)
    end subroutine

    subroutine copy_data2D_real8_from_real4(this, domain, copy_data)
        class(data2D_real8_type), intent(inout) :: this
        type(data2D_real4_type), intent(in) :: copy_data
        type(domain_type), intent(in) :: domain
        integer :: k, m, n

        _OMP_BLOCKS_PARALLEL_BEGIN_
        do k = 1, domain%bcount
            do n = domain%bbnd_y1(k), domain%bbnd_y2(k)
                do m = domain%bbnd_x1(k), domain%bbnd_x2(k)
                    this%block(k)%field(m, n) = real(copy_data%block(k)%field(m, n), wp8)
                enddo
            enddo
        enddo
        _OMP_BLOCKS_PARALLEL_END_

#ifdef _GPU_MODE_
        ! Copy from host to device, because device to device may be unsupported
        do k = 1, domain%bcount
            this%block(k)%field_device = this%block(k)%field
        enddo
#endif
    end subroutine

    subroutine copy_data2D_real8_from_real8(this, domain, copy_data)
        class(data2D_real8_type), intent(inout) :: this
        type(data2D_real8_type), intent(in) :: copy_data
        type(domain_type), intent(in) :: domain
        integer :: k, m, n

        _OMP_BLOCKS_PARALLEL_BEGIN_
        do k = 1, domain%bcount
            do n = domain%bbnd_y1(k), domain%bbnd_y2(k)
                do m = domain%bbnd_x1(k), domain%bbnd_x2(k)
                    this%block(k)%field(m, n) = copy_data%block(k)%field(m, n)
                enddo
            enddo
        enddo
        _OMP_BLOCKS_PARALLEL_END_

#ifdef _GPU_MODE_
        ! Copy from host to device, because device to device may be unsupported
        do k = 1, domain%bcount
            this%block(k)%field_device = this%block(k)%field
        enddo
#endif
    end subroutine

    subroutine fill_data2D_real8(this, domain, val)
        class(data2D_real8_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        real(wp8), intent(in) :: val
        integer :: k, m, n

        _OMP_BLOCKS_PARALLEL_BEGIN_
        do k = 1, domain%bcount
            do n = domain%bbnd_y1(k), domain%bbnd_y2(k)
                do m = domain%bbnd_x1(k), domain%bbnd_x2(k)
                    this%block(k)%field(m, n) = val
                enddo
            enddo
        enddo
        _OMP_BLOCKS_PARALLEL_END_

#ifdef _GPU_MODE_
        do k = 1, domain%bcount
            this%block(k)%field_device = val
        enddo
#endif
    end subroutine

    subroutine fill_debug_data2D_real8(this, domain)
        class(data2D_real8_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        integer :: k, m, n

        _OMP_BLOCKS_PARALLEL_BEGIN_
        do k = 1, domain%bcount
            do n = domain%bbnd_y1(k), domain%bbnd_y2(k)
                do m = domain%bbnd_x1(k), domain%bbnd_x2(k)
                    this%block(k)%field(m, n) = real(k, wp8)
                enddo
            enddo
        enddo
        _OMP_BLOCKS_PARALLEL_END_

#ifdef _GPU_MODE_
        do k = 1, domain%bcount
            this%block(k)%field_device = real(k, wp8)
        enddo
#endif
    end subroutine

!------------------------------------------------------------------------------
    ! Init and clear for 3D data type
    ! real4 base type
    subroutine init_data3D_real4(this, domain, nz_start, nz_end)
        class(data3D_real4_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        integer, intent(in) :: nz_start, nz_end
        integer :: k, m, n

        allocate(this%block(domain%bcount))

        _OMP_BLOCKS_PARALLEL_BEGIN_
        do k = 1, domain%bcount
            allocate(this%block(k)%field(domain%bbnd_x1(k) : domain%bbnd_x2(k), domain%bbnd_y1(k) : domain%bbnd_y2(k), nz_start : nz_end))
            do n = domain%bbnd_y1(k), domain%bbnd_y2(k)
                do m = domain%bbnd_x1(k), domain%bbnd_x2(k)
                    this%block(k)%field(m, n, :) = 0.0
                enddo
            enddo
        enddo
        _OMP_BLOCKS_PARALLEL_END_

#ifdef _GPU_MODE_
        do k = 1, domain%bcount
            allocate(this%block(k)%field_device(domain%bbnd_x1(k) : domain%bbnd_x2(k), domain%bbnd_y1(k) : domain%bbnd_y2(k), nz_start : nz_end))
            this%block(k)%field_device = 0.0
        enddo
#endif
    end subroutine

    subroutine clear_data3D_real4(this, domain)
        class(data3D_real4_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        integer :: k

#ifdef _GPU_MODE_
        do k = 1, domain%bcount
            deallocate(this%block(k)%field_device)
        enddo
#endif

        do k = 1, domain%bcount
            deallocate(this%block(k)%field)
        enddo

        deallocate(this%block)
    end subroutine

    ! real8 base type
    subroutine init_data3D_real8(this, domain, nz_start, nz_end)
        class(data3D_real8_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        integer, intent(in) :: nz_start, nz_end
        integer :: k, m, n

        allocate(this%block(domain%bcount))
        
        _OMP_BLOCKS_PARALLEL_BEGIN_
        do k = 1, domain%bcount
            allocate(this%block(k)%field(domain%bbnd_x1(k) : domain%bbnd_x2(k), domain%bbnd_y1(k) : domain%bbnd_y2(k), nz_start : nz_end))
            do n = domain%bbnd_y1(k), domain%bbnd_y2(k)
                do m = domain%bbnd_x1(k), domain%bbnd_x2(k)
                    this%block(k)%field(m, n, :) = 0.0d0
                enddo
            enddo
        enddo
        _OMP_BLOCKS_PARALLEL_END_

#ifdef _GPU_MODE_
        do k = 1, domain%bcount
            allocate(this%block(k)%field_device(domain%bbnd_x1(k) : domain%bbnd_x2(k), domain%bbnd_y1(k) : domain%bbnd_y2(k), nz_start : nz_end))
            this%block(k)%field_device = 0.0d0
        enddo
#endif
    end subroutine

    subroutine clear_data3D_real8(this, domain)
        class(data3D_real8_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        integer :: k

#ifdef _GPU_MODE_
        do k = 1, domain%bcount
            deallocate(this%block(k)%field_device)
        enddo
#endif

        do k = 1, domain%bcount
            deallocate(this%block(k)%field)
        enddo
        
        deallocate(this%block)
    end subroutine

end module data_types_module