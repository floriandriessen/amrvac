module mod_eos_container
    use mod_global_parameters
    use mod_phys_dict
    implicit none
    private

    abstract interface
        subroutine convert_condition(ixI^L, ixO^L, w, x)
            use mod_global_parameters
            integer, intent(in)             :: ixI^L, ixO^L
            double precision, intent(inout) :: w(ixI^S, nw)
            double precision, intent(in)    :: x(ixI^S, 1:^ND)
        end subroutine convert_condition
    end interface
    
    type eos_container

        character(len=std_len) :: eos_type
        double precision :: He_abundance
        double precision :: gamma
        double precision :: gamma_minus_1
        double precision :: inv_gamma
        double precision :: inv_gamma_minus_1
        double precision :: inv_squared_c0
        double precision :: inv_squared_c

        type(phys_dict) :: phys_dict

        procedure (convert_condition), pointer, nopass :: to_conserved => null()
        procedure (convert_condition), pointer, nopass :: to_primitive => null()
    
    end type eos_container

    public :: eos_container

end module mod_eos_container
!> Needs a line after to pass the preprocesor