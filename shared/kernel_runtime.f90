module kernel_runtime_module

    implicit none
    save
    public

    integer, parameter :: max_kernels = 15

contains

subroutine get_kernel_id(name, kernel_id)
    character(*), intent(in) :: name
    integer, intent(inout) :: kernel_id

    select case (name)
        case ('stress_components_kernel')
            kernel_id = 1
        case ('hh_init_kernel')
            kernel_id = 2
        case ('hh_update_kernel')
            kernel_id = 3 
        case ('hh_shift_kernel')
            kernel_id = 4
        case ('uv_trans_vort_kernel')
            kernel_id = 5
        case ('uv_trans_kernel')
            kernel_id = 6
        case ('uv_diff2_kernel')
            kernel_id = 7
        case ('sw_update_ssh_kernel')
            kernel_id = 8
        case ('sw_update_uv')
            kernel_id = 9
        case ('sw_next_step')
            kernel_id = 10
        case default
            kernel_id = 15
    end select
end subroutine 

subroutine get_kernel_name(kernel_id, name)
    integer, intent(in) :: kernel_id
    character(*), intent(inout) :: name

    select case (kernel_id)
        case (1)
            name = 'stress_components_kernel'
        case (2)
            name = 'hh_init_kernel'
        case (3)
            name = 'hh_update_kernel'
        case (4)
            name = 'hh_shift_kernel'
        case (5)
            name = 'uv_trans_vort_kernel'
        case (6)
            name = 'uv_trans_kernel'
        case (7)
            name = 'uv_diff2_kernel'
        case (8)
            name = 'sw_update_ssh_kernel'
        case (9)
            name = 'sw_update_uv'
        case (10)
            name = 'sw_next_step'
        case default
            name = 'unknown_kernel'
    end select
end subroutine

end module kernel_runtime_module