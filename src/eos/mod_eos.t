!> This module defines the handling of the eos of amrvac

module mod_eos
    use mod_global_parameters

    implicit none
    private

    type eos_container

        integer :: eos_type
        double precision :: He_abundance
    
    end type eos_container

    !> eos object
    type(eos_container), public, allocatable     :: eos


    contains
        !> Read this module"s parameters from a file
        subroutine eos_read_params(files)
            character(len=*), intent(in) :: files(:)
            integer                      :: n

            !> Default values
            character(len=std_len)  :: eos_type
            double precision, parameter :: He_abundance = 0.1d0

            namelist /eos_list/ eos_type, He_abundance

            do n = 1, size(files)
                open(unitpar, file=trim(files(n)), status="old")
                read(unitpar, eos_list, end=111)
        111     close(unitpar)
            end do

            eos%eos_type = eos_type
            eos%He_abundance = He_abundance
        
        end subroutine eos_read_params

        subroutine eos_init()
            
            allocate(eos)
            call eos_read_params(par_files)

        end subroutine eos_init

        !> I think that the eos stuff is initialised like this but it has to be tested.
        !> Next steps are to add the methods where the module sets the relevant params inside the eos

end module mod_eos