module tracer_interface_module

    use kernel_interface_module
    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use mpp_module
    use decomposition_module, only: domain_type, domain => domain_data
    use data_types_module, only: data2D_real8_type, data2D_real4_type
    use ocean_module, only: ocean_type, ocean_data
    use grid_module, only: grid_type, grid_data
    use mpp_sync_module, only: hybrid_sync, sync, sync_parameters_type
    use tracer_module, only: tran_diff_fluxes_kernel, tran_diff_tracer_kernel, tracer_next_step_kernel

#include "macros/mpp_macros.fi"

    implicit none
    save
    private

    public :: envoke_tran_diff_fluxes_kernel, envoke_tran_diff_fluxes_sync
    public :: envoke_tran_diff_tracer_kernel, envoke_tran_diff_tracer_sync
    public :: envoke_tracer_next_step_kernel, envoke_tracer_next_step_sync

contains

!-----------------------------------------------------------------------------!
!------------------------------- Kernels -------- ----------------------------!
!-----------------------------------------------------------------------------!
subroutine envoke_tran_diff_fluxes_kernel(k, param)
    integer, intent(in) :: k
    type(kernel_parameters_type), intent(in) :: param

    call tran_diff_fluxes_kernel(domain%bnx_start(k), domain%bnx_end(k), domain%bny_start(k), domain%bny_end(k),  &
                                 domain%bbnd_x1(k),   domain%bbnd_x2(k), domain%bbnd_y1(k),   domain%bbnd_y2(k),  &
                                 grid_data %lcu %block(k)%field,  &
                                 grid_data %lcv %block(k)%field,  &
                                 grid_data %dxt %block(k)%field,  &
                                 grid_data %dyt %block(k)%field,  &
                                 grid_data %dxh %block(k)%field,  &
                                 grid_data %dyh %block(k)%field,  &
                                 grid_data %hhu %block(k)%field,  &
                                 grid_data %hhv %block(k)%field,  &
                                 ocean_data%ff1(param%data_id)  %block(k)%field,  &
                                 ocean_data%ff1p(param%data_id) %block(k)%field,  &
                                 ocean_data%ubrtr %block(k)%field,  &
                                 ocean_data%vbrtr %block(k)%field,  &
                                 ocean_data%mu    %block(k)%field,  &
                                 1.0d0,                             &
                                 ocean_data%flux_x%block(k)%field,  &
                                 ocean_data%flux_y%block(k)%field)
end subroutine

subroutine envoke_tran_diff_fluxes_sync(k, sync_parameters)
    integer, intent(in) :: k
    type(sync_parameters_type), intent(in) :: sync_parameters

    call hybrid_sync(k, sync_parameters, 1, domain, ocean_data%flux_x)
    call hybrid_sync(k, sync_parameters, 2, domain, ocean_data%flux_y)
end subroutine

subroutine envoke_tran_diff_tracer_kernel(k, param)
    integer, intent(in) :: k
    type(kernel_parameters_type), intent(in) :: param

    call tran_diff_tracer_kernel(domain%bnx_start(k), domain%bnx_end(k), domain%bny_start(k), domain%bny_end(k),  &
                                 domain%bbnd_x1(k),   domain%bbnd_x2(k), domain%bbnd_y1(k),   domain%bbnd_y2(k),  &
                                 grid_data %lu %block(k)%field,   &
                                 grid_data %dx %block(k)%field,   &
                                 grid_data %dy %block(k)%field,   &
                                 param%tau,                       &
                                 grid_data %hhq_n  %block(k)%field,  &
                                 grid_data %hhq_p  %block(k)%field,  &
                                 ocean_data%flux_x %block(k)%field,  &
                                 ocean_data%flux_y %block(k)%field,  &
                                 ocean_data%ff1p(param%data_id) %block(k)%field,  &
                                 ocean_data%ff1n(param%data_id) %block(k)%field)
end subroutine

subroutine envoke_tran_diff_tracer_sync(k, sync_parameters)
    integer, intent(in) :: k
    type(sync_parameters_type), intent(in) :: sync_parameters

    call hybrid_sync(k, sync_parameters, 1, domain, ocean_data%ff1n(sync_parameters%data_id))
end subroutine

subroutine envoke_tracer_next_step_kernel(k, param)
    integer, intent(in) :: k
    type(kernel_parameters_type), intent(in) :: param

    call tracer_next_step_kernel(domain%bnx_start(k), domain%bnx_end(k), domain%bny_start(k), domain%bny_end(k),  &
                                 domain%bbnd_x1(k),   domain%bbnd_x2(k), domain%bbnd_y1(k),   domain%bbnd_y2(k),  &
                                 param%time_smooth,              &
                                 grid_data %lu %block(k)%field,  &
                                 ocean_data%ff1n(param%data_id)%block(k)%field,  &
                                 ocean_data%ff1p(param%data_id)%block(k)%field,  &
                                 ocean_data%ff1 (param%data_id)%block(k)%field)
end subroutine

subroutine envoke_tracer_next_step_sync(k, sync_parameters)
    integer, intent(in) :: k
    type(sync_parameters_type), intent(in) :: sync_parameters

end subroutine

end module tracer_interface_module