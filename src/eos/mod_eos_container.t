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
    
    type eos_table_container
    
        character(len=std_len) :: filename
        double precision, allocatable :: table(:,:)
        integer :: dim1, dim2
        double precision :: var1_min, var1_max, var2_min, var2_max

    end type eos_table_container

    type eos_container

        character(len=std_len) :: eos_type
        logical :: ionE
        logical :: table_check
        double precision :: He_abundance
        double precision :: gamma
        double precision :: gamma_minus_1
        double precision :: inv_gamma
        double precision :: inv_gamma_minus_1
        double precision :: inv_squared_c0
        double precision :: inv_squared_c
        
        double precision, allocatable :: T_space(:^D&)
        double precision, allocatable :: y_space(:^D&)
        !> EoS cast as Eint = f * pressure
        !> f = 1/(gamma-1) + ionE/pressure
        double precision, allocatable :: p_to_eint_factor(:^D&)

        type(phys_dict) :: phys_dict
        type(eos_table_container) :: neOnH

        procedure (convert_condition), pointer, nopass :: to_conserved => null()
        procedure (convert_condition), pointer, nopass :: to_primitive => null()
    
    end type eos_container

    public :: eos_container

end module mod_eos_container
!> Needs a line after to pass the preprocesor