module mod_hd_eos
    use mod_global_parameters
    use mod_physics
    use mod_eos
    use mod_eos_container
    use mod_hd_phys
    use mod_timing

    use mod_comm_lib, only: mpistop

    implicit none
    private

    public :: hd_link_eos

    contains

        !> Link the appropriate EOS conversion routines based on the selected EoS type
        subroutine hd_link_eos()
            use mod_usr_methods

            if (eos%eos_type == 'FI') then
                eos%to_conserved => hd_to_conserved
                eos%to_primitive => hd_to_primitive
            else if (eos%eos_type == 'LTE') then
                eos%to_conserved => hd_to_conserved_LTE
                eos%to_primitive => hd_to_primitive_LTE
            else
                call mpistop('Error: Unknown HD EOS type: ' // trim(eos%eos_type))
            end if

            !> Going to keep this approach for now until fully tested 
            !> - then refactor to remove phys_to_prim/con references
            phys_to_primitive => eos%to_primitive
            phys_to_conserved => eos%to_conserved
            phys_get_rho      => eos%get_rho
            phys_get_pthermal => eos%get_thermal_pressure
            phys_bind_eos_to_source  => bind_eos_to_source
            
            eos%p_to_e => p_to_e !> suitable for both FI and LTE

            ! choose Rfactor in ideal gas law
            if(hd_partial_ionization) then
                eos%get_Rfactor=>Rfactor_from_PI_temperature !> defined in eos file
                phys_update_temperature => update_PI_temperature !> defined in eos file
            else if(associated(usr_Rfactor)) then
                eos%get_Rfactor=>usr_Rfactor
            else
                eos%get_Rfactor=>Rfactor_from_constant_ionization
            end if

        end subroutine hd_link_eos

        subroutine bind_eos_to_source() !> this is called in eos_finalise through mod_physics procedure linking
            if (allocated(tc_fl)) then
                tc_fl%get_temperature_from_conserved => eos%get_temperature_from_etot
                tc_fl%get_temperature_from_eint => eos%get_temperature_from_eint
                tc_fl%get_rho => eos%get_rho
            end if

            if (allocated(rc_fl)) then
                rc_fl%get_rho => eos%get_rho
                rc_fl%get_pthermal => eos%get_thermal_pressure
                rc_fl%get_var_Rfactor => eos%get_Rfactor
                rc_fl%get_Te => eos%get_Te
            end if
        end subroutine bind_eos_to_source

        !> Transform primitive variables into conservative ones
        subroutine hd_to_conserved(ixI^L, ixO^L, w, x)
            use mod_global_parameters
            use mod_dust, only: dust_to_conserved
            integer, intent(in)             :: ixI^L, ixO^L
            double precision, intent(inout) :: w(ixI^S, nw)
            double precision, intent(in)    :: x(ixI^S, 1:ndim)

            integer :: ix^D

            timeeos0 = MPI_WTIME() !> For monitoring cost of eos module

            {do ix^DB=ixOmin^DB,ixOmax^DB\}
                if (hd_energy) then
                    ! Calculate total energy from pressure and kinetic energy
                    w(ix^D,e_)=w(ix^D,p_)*eos%inv_gamma_minus_1+&
                        half*(^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)
                end if
                ! Convert velocity to momentum
                ^C&w(ix^D,m^C_)=w(ix^D,rho_)*w(ix^D,m^C_)\
            {end do\}

            if (hd_dust) then
                call dust_to_conserved(ixI^L, ixO^L, w, x)
            end if

            timeeos_tot=timeeos_tot+(MPI_WTIME()-timeeos0) !> For monitoring cost of eos module

        end subroutine hd_to_conserved

        !> Transform conservative variables into primitive ones
        subroutine hd_to_primitive(ixI^L, ixO^L, w, x)
            use mod_global_parameters
            use mod_dust, only: dust_to_primitive
            integer, intent(in)             :: ixI^L, ixO^L
            double precision, intent(inout) :: w(ixI^S, nw)
            double precision, intent(in)    :: x(ixI^S, 1:ndim)

            double precision                :: inv_rho
            integer :: ix^D

            timeeos0 = MPI_WTIME() !> For monitoring cost of eos module

            {do ix^DB=ixOmin^DB,ixOmax^DB\}
                inv_rho = 1.d0/w(ix^D,rho_)
                ! Convert momentum to velocity
                ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\
                ! Calculate pressure = (gamma-1) * (e-ek)
                if(hd_energy) then
                    ! Compute pressure
                    w(ix^D,p_)=(eos%gamma_minus_1)*(w(ix^D,e_)&
                        -half*w(ix^D,rho_)*(^C&w(ix^D,m^C_)**2+))
                end if
            {end do\}

            if (fix_small_values) then
                call hd_handle_small_values(.false., w, x, ixI^L, ixO^L, 'hd_to_primitive')
            end if

            ! Convert dust momentum to dust velocity
            if (hd_dust) then
                call dust_to_primitive(ixI^L, ixO^L, w, x)
            end if

            timeeos_tot=timeeos_tot+(MPI_WTIME()-timeeos0) !> For monitoring cost of eos module

        end subroutine hd_to_primitive

        subroutine hd_to_conserved_LTE(ixI^L, ixO^L, w, x)
            use mod_global_parameters
            use mod_dust, only: dust_to_conserved
            integer, intent(in)             :: ixI^L, ixO^L
            double precision, intent(inout) :: w(ixI^S, nw)
            double precision, intent(in)    :: x(ixI^S, 1:ndim)

            timeeos0 = MPI_WTIME() !> For monitoring cost of eos module

            ! Calculate total energy from pressure and kinetic energy
            call p_to_e(ixI^L, ixO^L, w, x)

            ! Convert velocity to momentum
            ^C&w(ixO^S,m^C_)=w(ixO^S,rho_)*w(ixO^S,m^C_)\

            if (hd_dust) then
                call dust_to_conserved(ixI^L, ixO^L, w, x)
            end if

            timeeos_tot=timeeos_tot+(MPI_WTIME()-timeeos0) !> For monitoring cost of eos module
            
        end subroutine hd_to_conserved_LTE

        subroutine p_to_e(ixI^L, ixO^L, w, x)
            !> Needed to separate this for compatibility with hd_small_values_check
            !> i.e., to avoid circular calls
            use mod_global_parameters
            integer, intent(in)             :: ixI^L, ixO^L
            double precision, intent(inout) :: w(ixI^S, nw)
            double precision, intent(in)    :: x(ixI^S, 1:ndim)

            integer :: ix^D
            double precision :: p_to_eint
            double precision :: nH(ixI^S), nH_in(ixI^S), p_in(ixI^S)
            
            if (eos%ionE) then
                call eos%get_nH(w, x, ixI^L, ixO^L, nH)
                nH_in(ixO^S) = dlog10(nH(ixO^S))
                p_in(ixO^S) = dlog10(w(ixO^S,p_)) - nH_in(ixO^S)
            endif

            p_to_eint = eos%inv_gamma_minus_1 !> default assumed conversion value
            {do ix^DB=ixOmin^DB,ixOmax^DB\}
                if (hd_energy) then
                    !> calculate scaling according to tables
                    if (eos%ionE) p_to_eint = p2eint_from_nH_p(nH_in(ix^D), p_in(ix^D))
                    w(ix^D,e_)=w(ix^D,p_)*p_to_eint+&
                        half*(^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)
                end if
            {end do\}

        end subroutine p_to_e

        subroutine hd_to_primitive_LTE(ixI^L, ixO^L, w, x)
            use mod_global_parameters
            use mod_dust, only: dust_to_primitive
            integer, intent(in)             :: ixI^L, ixO^L
            double precision, intent(inout) :: w(ixI^S, nw)
            double precision, intent(in)    :: x(ixI^S, 1:ndim)

            double precision                :: inv_rho
            double precision                :: nH(ixI^S)
            double precision                :: log_nH(ixI^S)
            double precision                :: eint_val, p_guess
            integer :: ix^D

            timeeos0 = MPI_WTIME() !> For monitoring cost of eos module

            call eos%get_nH(w, x, ixI^L, ixO^L, nH)

            ! Cache log10(nH) for all cells (used by root-finder for IonE)
            if (eos%ionE) then
                log_nH(ixO^S) = dlog10(nH(ixO^S))
            end if

            {do ix^DB=ixOmin^DB,ixOmax^DB\}
                inv_rho = 1.d0/w(ix^D,rho_)
                ! Convert momentum to velocity
                ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\
                ! Calculate pressure
                if(hd_energy) then
                    if (eos%ionE) then
                        ! OLD: use stored Ne/Te from previous update_eos (closure gap with p2eint table)
                        w(ix^D,p_)=nH(ix^D) * (1.0d0 + eos%He_abundance + (w(ix^D,Ne_) / nH(ix^D))) * w(ix^D,Te_)
                        ! ! NEW: Root-find p by inverting the same p2eint table used in p_to_e (exact closure)
                        ! eint_val = w(ix^D,e_) - half*w(ix^D,rho_)*(^C&w(ix^D,m^C_)**2+)
                        ! ! Initial guess from stored Ne and Te (close to true p)
                        ! p_guess = nH(ix^D) * (1.0d0 + eos%He_abundance + w(ix^D,Ne_)/nH(ix^D)) * w(ix^D,Te_)
                        ! w(ix^D,p_) = p_from_eint_IonE(eint_val, nH(ix^D), log_nH(ix^D), p_guess)
                    else
                        w(ix^D,p_)=(eos%gamma_minus_1)*(w(ix^D,e_)&
                            -half*w(ix^D,rho_)*(^C&w(ix^D,m^C_)**2+))
                    end if
                end if
            {end do\}

            if (fix_small_values) then
                call hd_handle_small_values(.false., w, x, ixI^L, ixO^L, 'hd_to_primitive_LTE')
            end if

            ! Convert dust momentum to dust velocity
            if (hd_dust) then
                call dust_to_primitive(ixI^L, ixO^L, w, x)
            end if

            timeeos_tot=timeeos_tot+(MPI_WTIME()-timeeos0) !> For monitoring cost of eos module

        end subroutine hd_to_primitive_LTE

        subroutine Rfactor_from_constant_ionization(w,x,ixI^L,ixO^L,Rfactor)
            use mod_global_parameters
            integer, intent(in) :: ixI^L, ixO^L
            double precision, intent(in) :: w(ixI^S,1:nw)
            double precision, intent(in) :: x(ixI^S,1:ndim)
            double precision, intent(out):: Rfactor(ixI^S)

            Rfactor(ixO^S)=RR

        end subroutine Rfactor_from_constant_ionization

    end module mod_hd_eos
!> Needs a line after to pass the preprocesor