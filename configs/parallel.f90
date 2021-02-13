module config_parallel_module
    ! Parallel config for model

    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4
    use mpp_module
    use rwpar_routes

    implicit none
    save
    public

    integer :: mod_decomposition
    character(128) :: file_decomposition
    integer :: bppnx, bppny
    integer :: parallel_dbg
    integer :: parallel_mod
    character(128) :: file_output
    integer :: dlb_balance_steps
    integer :: dlb_model_steps
    
contains

    subroutine load_config_parallel_from_file_and_cmd(filepar)
        character(*), intent(in) :: filepar
        character(256) :: comments(256)
        integer :: nofcom, ierr
        character(len=32) :: arg

        ! reading parameters from file
        if (mpp_is_master()) then
            call readpar(filepar, comments, nofcom)
        endif
        call mpi_bcast(comments, 256*256, mpi_character, 0, mpp_cart_comm, ierr)

        read(comments(1),*) mod_decomposition
        call get_first_lexeme(comments(2), file_decomposition)
        read(comments(3),*) bppnx
        read(comments(4),*) bppny
        read(comments(5),*) parallel_dbg
        read(comments(6),*) parallel_mod
        call get_first_lexeme(comments(7), file_output)
        read(comments(8),*) dlb_balance_steps
        read(comments(9),*) dlb_model_steps

        if (command_argument_count() > 2) then
            call get_command_argument(1, arg)
            read(arg, *) mod_decomposition
            call get_command_argument(2, arg)
            read(arg, *) bppnx
            call get_command_argument(3, arg)
            read(arg, *) bppny
            print *, 'Warning: Overwtite mode, bppnx and bppny:', mod_decomposition, bppnx, bppny
        endif
        
        if (mpp_is_master()) then
            print *, 'Parallel decomposition config:'
            print *, 'mod_decomposition=', mod_decomposition
            print *, '(ignore this) decomposition file: ', file_decomposition
            print *, 'bppnx=', bppnx
            print *, 'bppny=', bppny
            print *, '(ignore this) parallel_dbg=', parallel_dbg
            print *, '(ignore this) parallel_mod=', parallel_mod
            print *, '(ignore this) output file: ', file_output
            print *, 'Balance steps for DLB: ', dlb_balance_steps
            print *, 'Model steps for DLB  : ', dlb_model_steps
        endif
    
        call mpp_sync_output()
    end subroutine

end module config_parallel_module