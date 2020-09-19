module gridcon_module
    ! Grid consruction module

    use mpp_module
    use decomposition_module, only: domain_type
    use grid_module, only: grid_type, grid_global_type
    use mpp_sync_module, only: sync
    use grid_interface_module, only: envoke_lu_init_kernel, envoke_lu_lv_init_kernel
    use errors_module, only: abort_model

#include "macros/mpp_macros.fi"

    implicit none
    save
    private

    public :: gridcon

contains

    subroutine gridcon(domain, grid_global_data, grid_data)
        ! Grid consruction module by temperature mask.
        ! subroutin for construction pass boundary, velosity and bottom masks
        ! using temperature mask in diogin standart

        type(domain_type), intent(in) :: domain
        type(grid_global_type), intent(in) :: grid_global_data
        type(grid_type), intent(inout) :: grid_data

        ! conversion integer diogin mask to real model mask
        call envoke_lu_init_kernel(domain, grid_data, grid_global_data)

        ! forming mask for depth grid points
        ! forming luh from lu, which have land neibours in luh.
        ! constructing array luh for relief hh.

        if (mpp_is_master()) then
          write(*,*) 'Construction of H-grid masks: '
          write(*,*) 'LUH (includes boundary) and LUU (does not include boundary)'
          write(*,*) 'Construction of U- and V-grid masks: '
          write(*,*) 'LCU and LCV (do not include boundary) and LLU and LLV (include boundary)'
        endif

        call envoke_lu_lv_init_kernel(domain, grid_data)
        
    end subroutine gridcon

end module gridcon_module