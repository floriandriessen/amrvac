module mod_eos_container
    use mod_global_parameters
    implicit none
    private

    public :: eos_container, eos_table_container

    abstract interface
        subroutine convert_condition(ixI^L, ixO^L, w, x)
            use mod_global_parameters
            integer, intent(in)             :: ixI^L, ixO^L
            double precision, intent(inout) :: w(ixI^S, nw)
            double precision, intent(in)    :: x(ixI^S, 1:ndim)
        end subroutine convert_condition

        subroutine update_eos_spaces(ixI^L, ixO^L, w, x)
            use mod_global_parameters
            integer, intent(in)             :: ixI^L, ixO^L
            double precision, intent(inout) :: w(ixI^S, nw)
            double precision, intent(in)    :: x(ixI^S, 1:ndim)
        end subroutine update_eos_spaces

        subroutine update_eos_spaces_alt(ixI^L, ixO^L, wCT, w, x)
            use mod_global_parameters
            integer, intent(in)             :: ixI^L, ixO^L
            double precision, intent(inout) :: w(ixI^S, nw)
            double precision, intent(in)    :: wCT(ixI^S, nw), x(ixI^S, 1:ndim)
        end subroutine update_eos_spaces_alt

        subroutine get_eos_variable(ixI^L, ixO^L, w, x, res)
            use mod_global_parameters
            integer, intent(in)             :: ixI^L, ixO^L
            double precision, intent(in) :: w(ixI^S, nw)
            double precision, intent(in)    :: x(ixI^S, 1:ndim)
            double precision, intent(out)   :: res(ixI^S)
        end subroutine get_eos_variable

        subroutine get_eos_variable_alt(w, x, ixI^L, ixO^L,  res)
            use mod_global_parameters
            integer, intent(in)             :: ixI^L, ixO^L
            double precision, intent(in) :: w(ixI^S, nw)
            double precision, intent(in)    :: x(ixI^S, 1:ndim)
            double precision, intent(out)   :: res(ixI^S)
        end subroutine get_eos_variable_alt
    end interface
    
    type eos_table_container
    
        character(len=std_len) :: filename
        double precision, allocatable :: table(:,:)
        integer :: dim1, dim2
        double precision :: var1_min, var1_max, var2_min, var2_max
        double precision :: step_inv_1, step_inv_2  !> precomputed (dim-1)/(max-min) for each axis

    end type eos_table_container

    type eos_container

        character(len=std_len) :: eos_type
        character(len=20)     :: method = 'tables'        !> 'tables' or 'analytic'
        character(len=20)     :: gamma1_method = 'exact'   !> 'exact' or 'effective'
        character(len=20)     :: inversion = 'bisect'      !> 'bisect' or 'newton' (analytic only)
        logical :: ionE
        character(len=std_len) :: table_location
        logical :: table_check
        double precision :: He_abundance
        double precision :: gamma
        double precision :: gamma_minus_1
        double precision :: inv_gamma
        double precision :: inv_gamma_minus_1
        double precision :: inv_squared_c0
        double precision :: inv_squared_c
        double precision :: nH2rhoFactor

        !> Fully-ionised regime bypass constants (precomputed in eos_finalise)
        double precision :: eion_per_nH          !> Total ionisation energy per nH [code units]
        double precision :: eint_rho_FI_threshold !> eint/rho above which gas is fully ionised [code units]
        double precision :: p_rho_FI_threshold    !> p/rho above which gas is fully ionised [code units]
        double precision :: n_per_nH_FI          !> Total particles per nH when FI = 2 + 3*A_He
        double precision :: neOnH_FI             !> ne/nH when fully ionised = 1 + 2*A_He
        logical :: disable_FI_bypass = .false.  !> Set .true. to force all cells through table path

        !>Leaving a comment here to remind me that it's a bad idea to have anything outside
        !> of the w() array if one wants to /ensure/ the fields are consistent 
        !> given OMP directives and load balancing

        ! Expected components of the EoS
        type(eos_table_container) :: neOnH
        type(eos_table_container) :: T
        type(eos_table_container) :: p2eint
        type(eos_table_container) :: gamma1    !> First adiabatic index Gamma_1(nH, eint/nH)
        type(eos_table_container) :: gamma1_p  !> First adiabatic index Gamma_1(nH, p/nH) for fast csound2
        type(eos_table_container) :: eint_from_T  !> Inverse T table: eint/nH(nH, T)
        type(eos_table_container) :: log_p        !> Merged log10(p/nH)(nH, eint/nH) for WB bisection

        procedure (convert_condition), pointer, nopass :: to_conserved => null()
        procedure (convert_condition), pointer, nopass :: to_primitive => null()
        procedure (convert_condition), pointer, nopass :: p_to_e => null()
        procedure (update_eos_spaces), pointer, nopass :: update_eos => null()
        procedure (update_eos_spaces_alt), pointer, nopass :: update_PI_temperature => null()
        procedure (get_eos_variable_alt), pointer, nopass :: get_thermal_pressure => null()
        procedure (get_eos_variable_alt), pointer, nopass :: get_temperature_from_eint => null()
        procedure (get_eos_variable_alt), pointer, nopass :: get_temperature_from_etot => null()
        procedure (get_eos_variable_alt), pointer, nopass :: get_Rfactor => null()
        procedure (get_eos_variable_alt), pointer, nopass :: get_rho => null()
        procedure (get_eos_variable_alt), pointer, nopass :: get_nH => null()
        procedure (get_eos_variable_alt), pointer, nopass :: get_Te => null()
        procedure (get_eos_variable_alt), pointer, nopass :: get_csound2 => null()
        procedure (get_eos_variable_alt), pointer, nopass :: get_gamma1 => null()
        procedure (convert_condition), pointer, nopass :: to_prolong => null()
        procedure (convert_condition), pointer, nopass :: from_prolong => null()

    end type eos_container

end module mod_eos_container
!> Needs a line after to pass the preprocesor