module data_types_module
    ! Module description

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use decomposition_module

    implicit none
    save
    private

!------------------------------------------------------------------------------
    ! Data types for Parallel System layer
    ! 1D blocks
    type, public :: block1D_real4_type
        real(wp4), pointer :: field(:)
    end type block1D_real4_type

    type, public :: block1D_real8_type
        real(wp8), pointer :: field(:)
    end type block1D_real8_type

    ! 2D blocks
    type, public :: block2D_real4_type
        real(wp4), pointer :: field(:, :)
    end type block2D_real4_type

    type, public :: block2D_real8_type
        real(wp8), pointer :: field(:, :)
    end type block2D_real8_type

    ! 3D blocks
    type, public :: block3D_real4_type
        real(wp4), pointer :: field(:, :, :)
    end type block3D_real4_type

    type, public :: block3D_real8_type
        real(wp8), pointer :: field(:, :, :)
    end type block3D_real8_type

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
        procedure, public :: init_from_domian => init_from_domain_data1D_real8
        procedure, public :: clear => clear_data1D_real8
    end type data1D_real8_type

    ! 2D data
    type, public :: data2D_real4_type
        type(block2D_real4_type), pointer :: block(:)
    contains
        procedure, public :: init => init_data2D_real4
        procedure, public :: clear => clear_data2D_real4
    end type data2D_real4_type

    type, public :: data2D_real8_type
        type(block2D_real8_type), pointer :: block(:)
    contains
        procedure, public :: init => init_data2D_real8
        procedure, public :: clear => clear_data2D_real8
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
        do k = 1, domain%bcount
            allocate(this%block(k)%field(n_start : n_end))
            this%block(k)%field = 0.0
        enddo
    end subroutine
    
    subroutine init_from_domain_data1D_real4(this, domain, is_x_direction)
        class(data1D_real4_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        logical, intent(in) :: is_x_direction
        integer :: k
        allocate(this%block(domain%bcount))
        do k = 1, domain%bcount
            if (is_x_direction) then
                allocate(this%block(k)%field(domain%bnx_start(k) : domain%bnx_end(k)))
            else
                allocate(this%block(k)%field(domain%bny_start(k) : domain%bny_end(k)))
            endif
            this%block(k)%field = 0.0
        enddo
    end subroutine

    subroutine clear_data1D_real4(this, domain)
        class(data1D_real4_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        integer :: k
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
        do k = 1, domain%bcount
            allocate(this%block(k)%field(n_start : n_end))
            this%block(k)%field = 0.0
        enddo
    end subroutine
    
    subroutine init_from_domain_data1D_real8(this, domain, is_x_direction)
        class(data1D_real8_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        logical, intent(in) :: is_x_direction
        integer :: k
        allocate(this%block(domain%bcount))
        do k = 1, domain%bcount
            if (is_x_direction) then
                allocate(this%block(k)%field(domain%bnx_start(k) : domain%bnx_end(k)))
            else
                allocate(this%block(k)%field(domain%bny_start(k) : domain%bny_end(k)))
            endif
            this%block(k)%field = 0.0
        enddo
    end subroutine

    subroutine clear_data1D_real8(this, domain)
        class(data1D_real8_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        integer :: k
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
        integer :: k
        allocate(this%block(domain%bcount))
        do k = 1, domain%bcount
            allocate(this%block(k)%field(domain%bnx_start(k) : domain%bnx_end(k), domain%bny_start(k) : domain%bny_end(k)))
            this%block(k)%field = 0.0
        enddo
    end subroutine

    subroutine clear_data2D_real4(this, domain)
        class(data2D_real4_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        integer :: k
        do k = 1, domain%bcount
            deallocate(this%block(k)%field)
        enddo
        deallocate(this%block)
    end subroutine

    ! real8 base type
    subroutine init_data2D_real8(this, domain)
        class(data2D_real8_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        integer :: k
        allocate(this%block(domain%bcount))
        do k = 1, domain%bcount
            allocate(this%block(k)%field(domain%bnx_start(k) : domain%bnx_end(k), domain%bny_start(k) : domain%bny_end(k)))
            this%block(k)%field = 0.0
        enddo
    end subroutine

    subroutine clear_data2D_real8(this, domain)
        class(data2D_real8_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        integer :: k
        do k = 1, domain%bcount
            deallocate(this%block(k)%field)
        enddo
        deallocate(this%block)
    end subroutine

!------------------------------------------------------------------------------
    ! Init and clear for 3D data type
    ! real4 base type
    subroutine init_data3D_real4(this, domain, nz_start, nz_end)
        class(data3D_real4_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        integer, intent(in) :: nz_start, nz_end
        integer :: k
        allocate(this%block(domain%bcount))
        do k = 1, domain%bcount
            allocate(this%block(k)%field(domain%bnx_start(k) : domain%bnx_end(k), domain%bny_start(k) : domain%bny_end(k), nz_start : nz_end))
            this%block(k)%field = 0.0
        enddo
    end subroutine

    subroutine clear_data3D_real4(this, domain)
        class(data3D_real4_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        integer :: k
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
        integer :: k
        allocate(this%block(domain%bcount))
        do k = 1, domain%bcount
            allocate(this%block(k)%field(domain%bnx_start(k) : domain%bnx_end(k), domain%bny_start(k) : domain%bny_end(k), nz_start : nz_end))
            this%block(k)%field = 0.0
        enddo
    end subroutine

    subroutine clear_data3D_real8(this, domain)
        class(data3D_real8_type), intent(inout) :: this
        type(domain_type), intent(in) :: domain
        integer :: k
        do k = 1, domain%bcount
            deallocate(this%block(k)%field)
        enddo
        deallocate(this%block)
    end subroutine

end module data_types_module