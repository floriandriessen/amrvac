!> This module defines the handling of the eos of amrvac

module mod_eos
    use mod_global_parameters

    implicit none
    public

    type eos_fluid

        integer :: eos_type
    
    end type eos_fluid

    !> eos fluid object
    type(eos_fluid), public, allocatable     :: eos_fl


    contains
        !> Read this module"s parameters from a file
        subroutine eos_read_params(files)
            character(len=*), intent(in) :: files(:)
            integer                      :: n

            namelist /eos_list/ eos_type

            do n = 1, size(files)
                open(unitpar, file=trim(files(n)), status="old")
                read(unitpar, eos_list, end=111)
        111     close(unitpar)
            end do
        
        end subroutine eos_read_params

        subroutine eos_init(eos_fl)
            
            allocate(eos_fl)

            call eos_read_params(par_files)

        end subroutine eos_init

end module mod_eos