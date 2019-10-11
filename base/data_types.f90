module data_types_module
    ! Module description

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use decomposition_module

    implicit none
    save
    private

    ! Data types for Parallel System layer
    type, public :: block2D_real4_type
        real(wp4), pointer :: field(:, :)
    end type block2D_real4_type

    type, public :: block2D_real8_type
        real(wp8), pointer :: field(:, :)
    end type block2D_real8_type

    ! Data types for Algorithm layer
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

contains

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

end module data_types_module