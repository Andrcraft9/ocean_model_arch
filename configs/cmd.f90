module config_cmd_module

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use mpp_module
    use config_basinpar_module
    use config_parallel_module
    use config_sw_module

    implicit none
    save
    public

contains

    subroutine overwrite_configs_by_cmd()
        integer :: ierr
        character(len=32) :: arg
        
        if (mpp_is_master()) then
            if (command_argument_count() > 2) then
                call get_command_argument(1, arg)
                read(arg, *) mod_decomposition
                call get_command_argument(2, arg)
                read(arg, *) bppnx
                call get_command_argument(3, arg)
                read(arg, *) bppny
                print *, 'CMD: Overwtite by command line: mod_decomp(1), bppnx(2), bppny(3):', mod_decomposition, bppnx, bppny
            else
                print *, 'CMD: Ignore arguments from command line'
                print *, 'CMD: Possible command line arguments: mod_decomp(1), bppnx(2), bppny(3)'
            endif
        endif
        call mpi_bcast(mod_decomposition, 1, mpi_integer, 0, mpp_cart_comm, ierr)
        call mpi_bcast(bppnx, 1, mpi_integer, 0, mpp_cart_comm, ierr)
        call mpi_bcast(bppny, 1, mpi_integer, 0, mpp_cart_comm, ierr)

        call mpp_sync_output()
    end subroutine

end module
