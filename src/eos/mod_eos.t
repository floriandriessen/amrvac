!> This module defines the handling of the eos of amrvac

module mod_eos
    use mod_global_parameters
    use mod_eos_container
    implicit none
    private

    !> eos object
    type(eos_container), public, allocatable     :: eos

    public :: eos_init!, link_eos

    contains
        !> Read this module"s parameters from a file
        subroutine eos_read_params(files)
            character(len=*), intent(in) :: files(:)
            integer                      :: n

            !> Default values
            character(len=std_len)  :: eos_type
            double precision :: He_abundance
            double precision :: gamma
            double precision :: gamma_minus_1
            double precision :: inv_gamma
            double precision :: inv_gamma_minus_1

            eos_type = 'FI'
            He_abundance = 0.1d0
            gamma = 5.0d0/3.0d0
            gamma_minus_1 = gamma - 1.0d0
            inv_gamma = 1.0d0 / gamma
            inv_gamma_minus_1 = 1.0d0 / gamma_minus_1

            namelist /eos_list/ eos_type, He_abundance, gamma

            do n = 1, size(files)
                open(unitpar, file=trim(files(n)), status="old")
                read(unitpar, eos_list, end=111)
        111     close(unitpar)
            end do

            eos%eos_type = eos_type
            eos%He_abundance = He_abundance
            eos%gamma = gamma
            eos%gamma_minus_1 = eos%gamma - 1.0d0
            eos%inv_gamma = 1.0d0 / eos%gamma
            eos%inv_gamma_minus_1 = 1.0d0 / eos%gamma_minus_1
        
        end subroutine eos_read_params

        subroutine eos_init()
            
            allocate(eos)
            call eos_read_params(par_files)

        end subroutine eos_init

end module mod_eos
!> Needs a line after to pass the preprocesor