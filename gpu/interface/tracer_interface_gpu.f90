#include "macros/mpp_macros.fi"
#ifdef _GPU_MODE_

module tracer_interface_gpu_module

    use kernel_interface_module
    use cudafor
    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use mpp_module
    use decomposition_module, only: domain_type, domain => domain_data
    use data_types_module, only: data2D_real8_type, data2D_real4_type
    use ocean_module, only: ocean_type, ocean_data
    use grid_module, only: grid_type, grid_data
    use mpp_sync_module, only: hybrid_sync, sync, sync_parameters_type, mpp_sync_cuda_streams
    use tracer_gpu_module, only: tran_diff_fluxes_kernel_gpu, tran_diff_tracer_kernel_gpu, tracer_next_step_kernel_gpu

    implicit none
    save
    private

    public :: envoke_tran_diff_fluxes_kernel_gpu, envoke_tran_diff_fluxes_sync_gpu
    public :: envoke_tran_diff_tracer_kernel_gpu, envoke_tran_diff_tracer_sync_gpu
    public :: envoke_tracer_next_step_kernel_gpu, envoke_tracer_next_step_sync_gpu

contains

!-----------------------------------------------------------------------------!
!------------------------------- Kernels -------- ----------------------------!
!-----------------------------------------------------------------------------!
subroutine envoke_tran_diff_fluxes_kernel_gpu(k, param)
    integer, intent(in) :: k
    type(kernel_parameters_type), intent(in) :: param

    type(dim3) :: tBlock, grid

    tBlock = dim3(_GPU_BLOCK_X_, _GPU_BLOCK_Y_, 1)
    grid = dim3(ceiling(real( (domain%bnx_end(k) + 1) - (domain%bnx_start(k) - 1) + 1 ) / tBlock%x),  &
                ceiling(real( (domain%bny_end(k) + 1) - (domain%bny_start(k) - 1) + 1 ) / tBlock%y),  &
                1)

    call tran_diff_fluxes_kernel_gpu<<<grid, tBlock, 0, mpp_sync_cuda_streams(k)>>>(    &
                                 domain%bnx_start(k), domain%bnx_end(k), domain%bny_start(k), domain%bny_end(k),  &
                                 domain%bbnd_x1(k),   domain%bbnd_x2(k), domain%bbnd_y1(k),   domain%bbnd_y2(k),  &
                                 grid_data %lcu %block(k)%field_device,  &
                                 grid_data %lcv %block(k)%field_device,  &
                                 grid_data %dxt %block(k)%field_device,  &
                                 grid_data %dyt %block(k)%field_device,  &
                                 grid_data %dxh %block(k)%field_device,  &
                                 grid_data %dyh %block(k)%field_device,  &
                                 grid_data %hhu %block(k)%field_device,  &
                                 grid_data %hhv %block(k)%field_device,  &
                                 ocean_data%ff1(param%data_id)  %block(k)%field_device,  &
                                 ocean_data%ff1p(param%data_id) %block(k)%field_device,  &
                                 ocean_data%ubrtr %block(k)%field_device,  &
                                 ocean_data%vbrtr %block(k)%field_device,  &
                                 ocean_data%mu    %block(k)%field_device,  &
                                 1.0d0,                             &
                                 ocean_data%flux_x%block(k)%field_device,  &
                                 ocean_data%flux_y%block(k)%field_device)
end subroutine

subroutine envoke_tran_diff_fluxes_sync_gpu(k, sync_parameters)
    integer, intent(in) :: k
    type(sync_parameters_type), intent(in) :: sync_parameters

    call hybrid_sync(k, sync_parameters, 1, domain, ocean_data%flux_x)
    call hybrid_sync(k, sync_parameters, 2, domain, ocean_data%flux_y)
end subroutine

subroutine envoke_tran_diff_tracer_kernel_gpu(k, param)
    integer, intent(in) :: k
    type(kernel_parameters_type), intent(in) :: param

    type(dim3) :: tBlock, grid

    tBlock = dim3(_GPU_BLOCK_X_, _GPU_BLOCK_Y_, 1)
    grid = dim3(ceiling(real( (domain%bnx_end(k) + 1) - (domain%bnx_start(k) - 1) + 1 ) / tBlock%x),  &
                ceiling(real( (domain%bny_end(k) + 1) - (domain%bny_start(k) - 1) + 1 ) / tBlock%y),  &
                1)

    call tran_diff_tracer_kernel_gpu<<<grid, tBlock, 0, mpp_sync_cuda_streams(k)>>>(    &
                                 domain%bnx_start(k), domain%bnx_end(k), domain%bny_start(k), domain%bny_end(k),  &
                                 domain%bbnd_x1(k),   domain%bbnd_x2(k), domain%bbnd_y1(k),   domain%bbnd_y2(k),  &
                                 grid_data %lu %block(k)%field_device,   &
                                 grid_data %dx %block(k)%field_device,   &
                                 grid_data %dy %block(k)%field_device,   &
                                 param%tau,                              &
                                 grid_data %hhq_n  %block(k)%field_device,  &
                                 grid_data %hhq_p  %block(k)%field_device,  &
                                 ocean_data%flux_x %block(k)%field_device,  &
                                 ocean_data%flux_y %block(k)%field_device,  &
                                 ocean_data%ff1p(param%data_id) %block(k)%field_device,  &
                                 ocean_data%ff1n(param%data_id) %block(k)%field_device)
end subroutine

subroutine envoke_tran_diff_tracer_sync_gpu(k, sync_parameters)
    integer, intent(in) :: k
    type(sync_parameters_type), intent(in) :: sync_parameters

    call hybrid_sync(k, sync_parameters, 1, domain, ocean_data%ff1n(sync_parameters%data_id))
end subroutine

subroutine envoke_tracer_next_step_kernel_gpu(k, param)
    integer, intent(in) :: k
    type(kernel_parameters_type), intent(in) :: param

    type(dim3) :: tBlock, grid

    tBlock = dim3(_GPU_BLOCK_X_, _GPU_BLOCK_Y_, 1)
    grid = dim3(ceiling(real( (domain%bnx_end(k) + 1) - (domain%bnx_start(k) - 1) + 1 ) / tBlock%x),  &
                ceiling(real( (domain%bny_end(k) + 1) - (domain%bny_start(k) - 1) + 1 ) / tBlock%y),  &
                1)

    call tracer_next_step_kernel_gpu<<<grid, tBlock, 0, mpp_sync_cuda_streams(k)>>>(    &
                                 domain%bnx_start(k), domain%bnx_end(k), domain%bny_start(k), domain%bny_end(k),  &
                                 domain%bbnd_x1(k),   domain%bbnd_x2(k), domain%bbnd_y1(k),   domain%bbnd_y2(k),  &
                                 param%time_smooth,                     &
                                 grid_data %lu %block(k)%field_device,  &
                                 ocean_data%ff1n(param%data_id)%block(k)%field_device,  &
                                 ocean_data%ff1p(param%data_id)%block(k)%field_device,  &
                                 ocean_data%ff1 (param%data_id)%block(k)%field_device)
end subroutine

subroutine envoke_tracer_next_step_sync_gpu(k, sync_parameters)
    integer, intent(in) :: k
    type(sync_parameters_type), intent(in) :: sync_parameters

end subroutine

end module tracer_interface_gpu_module

#endif