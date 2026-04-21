module mod_hd_eos
    use mod_global_parameters
    use mod_physics
    use mod_eos
    use mod_eos_container
    use mod_hd_phys
    use mod_timing
    use mod_radiative_cooling, only: build_Y_mod_table

    use mod_comm_lib, only: mpistop

    implicit none
    private

    public :: hd_link_eos, hd_get_pthermal

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
            phys_get_pthermal => hd_get_pthermal
            eos%get_thermal_pressure => hd_get_pthermal
            phys_bind_eos_to_source  => bind_eos_to_source
            
            eos%p_to_e => p_to_e !> suitable for both FI and LTE

            ! Link sound speed and gamma1 computation
            if (eos%eos_type == 'LTE' .and. eos%ionE) then
                eos%get_csound2 => hd_get_csound2_LTE
                phys_get_gamma1 => hd_get_gamma1_LTE
                ! EoS-aware prolongation: interpolate in (rho, v, T) space
                phys_to_prolong   => hd_to_prolong_LTE
                phys_from_prolong => hd_from_prolong_LTE
            else
                eos%get_csound2 => hd_get_csound2_FI
                phys_get_gamma1 => hd_get_gamma1_FI
            end if

            ! choose Rfactor in ideal gas law (usr_Rfactor takes priority)
            if(associated(usr_Rfactor)) then
                eos%get_Rfactor=>usr_Rfactor
            else if (eos%eos_type == 'LTE') then
                eos%get_Rfactor => Rfactor_from_LTE
            else if(hd_partial_ionization) then
                eos%get_Rfactor=>Rfactor_from_PI_temperature !> defined in eos file
                phys_update_temperature => update_PI_temperature !> defined in eos file
            else
                eos%get_Rfactor=>Rfactor_from_constant_ionization
            end if

        end subroutine hd_link_eos

        subroutine bind_eos_to_source() !> this is called in eos_finalise through mod_physics procedure linking
            if (allocated(tc_fl)) then
                tc_fl%get_temperature_from_conserved => eos%get_temperature_from_etot
                !> Use fast bilinear T lookup for TC STS substeps (density fixed,
                !> TC flux dominated by corona where T(eint) is smooth)
                if (eos%eos_type == 'LTE' .and. eos%ionE) then
                    tc_fl%get_temperature_from_eint => get_temperature_from_eint_LTE_fast
                else
                    tc_fl%get_temperature_from_eint => eos%get_temperature_from_eint
                end if
                tc_fl%get_rho => eos%get_rho
                tc_fl%get_ne_nH => eos%get_ne_nH
                tc_fl%get_var_Rfactor => eos%get_Rfactor
            end if

            if (allocated(rc_fl)) then
                rc_fl%get_rho => eos%get_rho
                rc_fl%get_pthermal => eos%get_thermal_pressure
                rc_fl%get_var_Rfactor => eos%get_Rfactor
                rc_fl%get_Te => eos%get_Te
                rc_fl%get_ne_nH => eos%get_ne_nH
                !> Build the variable-c_V Townsend Y_mod table now that all
                !> EoS tables (eint_from_T, T, neOnH) are in code units.
                !> build_Y_mod_table checks coolmethod=='exact' and .not.isPPL
                !> internally and early-returns otherwise.
                if (eos%ionE) call build_Y_mod_table(rc_fl)
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

            timeeos_conv=timeeos_conv+(MPI_WTIME()-timeeos0)

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

            timeeos_conv=timeeos_conv+(MPI_WTIME()-timeeos0)

        end subroutine hd_to_primitive

        !> LTE primitive -> conserved conversion.
        !>
        !> On entry: rho_ = density, m_ = velocity, p_ = pressure.
        !> On exit:  rho_ = density (unchanged), m_ = momentum, e_ = total energy.
        subroutine hd_to_conserved_LTE(ixI^L, ixO^L, w, x)
            use mod_global_parameters
            use mod_dust, only: dust_to_conserved
            integer, intent(in)             :: ixI^L, ixO^L
            double precision, intent(inout) :: w(ixI^S, nw)
            double precision, intent(in)    :: x(ixI^S, 1:ndim)

            timeeos0 = MPI_WTIME()

            call p_to_e(ixI^L, ixO^L, w, x)

            ! Convert velocity to momentum
            ^C&w(ixO^S,m^C_)=w(ixO^S,rho_)*w(ixO^S,m^C_)\

            if (hd_dust) then
                call dust_to_conserved(ixI^L, ixO^L, w, x)
            end if

            timeeos_conv=timeeos_conv+(MPI_WTIME()-timeeos0)

        end subroutine hd_to_conserved_LTE

        subroutine p_to_e(ixI^L, ixO^L, w, x)
            !> Convert pressure to total energy: E = eint(rho, p) + KE.
            !>
            !> On entry: w(rho_) = density, w(m_) = velocity, w(p_) = pressure.
            !> On exit:  w(rho_) unchanged, w(m_) unchanged, w(e_) = total energy.
            !>
            !> Four paths for ionE (LTE with ionisation):
            !>   1. FI bypass (p/rho > threshold): analytic eint = p/(gamma-1) + eion*nH
            !>   2. Analytical Saha (eos_method='analytic'): direct Saha solve for T,y from p
            !>   3. WB mode (hd_well_balanced): cached bisection on forward tables
            !>   4. Standard: p2eint inverse table lookup (fast, ~0.01% round-trip error)
            use mod_global_parameters
            use mod_eos, only: saha_p_to_T
            integer, intent(in)             :: ixI^L, ixO^L
            double precision, intent(inout) :: w(ixI^S, nw)
            double precision, intent(in)    :: x(ixI^S, 1:ndim)

            integer :: ix^D
            double precision :: p_to_eint, p_over_rho
            double precision :: nH(ixI^S), nH_in(ixI^S), p_in(ixI^S)
            double precision :: log_eint_mid, eint_total

            if (eos%ionE) then
                call eos%get_nH(w, x, ixI^L, ixO^L, nH)
                nH_in(ixO^S) = dlog10(nH(ixO^S))
                if (eos%p2eint_method /= 'bisect') then
                    p_in(ixO^S) = dlog10(w(ixO^S,p_)) - nH_in(ixO^S)
                end if
            endif

            p_to_eint = eos%inv_gamma_minus_1
            {do ix^DB=ixOmin^DB,ixOmax^DB\}
                if (hd_energy) then
                    if (eos%ionE) then
                        p_over_rho = w(ix^D,p_) / w(ix^D,rho_)
                        if (p_over_rho > eos%p_rho_FI_threshold) then
                            !> FI bypass: exact inverse of to_primitive
                            p_to_eint = eos%inv_gamma_minus_1 &
                                + eos%eion_per_nH * nH(ix^D) / w(ix^D,p_)
                            w(ix^D,e_) = w(ix^D,p_)*p_to_eint + &
                                half*(^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)
                        else if (eos%method == 'analytic') then
                            !> Analytical Saha: bisect for T from p, return eint directly
                            block
                                double precision :: T_solve, y_solve, eint_nH_solve
                                call saha_p_to_T(nH(ix^D), w(ix^D,p_), &
                                    T_solve, y_solve, eint_nH_solve)
                                w(ix^D,e_) = eint_nH_solve * nH(ix^D) + &
                                    half*(^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)
                            end block
                        else if (eos%p2eint_method == 'bisect') then
                            !> Bisection on forward p table: finds eint such that
                            !> p(nH, eint) = p_target exactly. More accurate than
                            !> the p2eint table in the ionisation zone.
                            call eint_from_p_bisect(nH_in(ix^D), &
                                dlog10(w(ix^D,p_)), log_eint_mid)
                            eint_total = nH(ix^D) * 10.0d0**log_eint_mid
                            eint_total = max(eint_total, &
                                nH(ix^D) * 10.0d0**eos%T%var2_min)
                            w(ix^D,e_) = eint_total + &
                                half*(^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)
                        else
                            !> Standard table lookup (fast, default)
                            p_to_eint = p2eint_from_nH_p(nH_in(ix^D), p_in(ix^D))
                            w(ix^D,e_) = w(ix^D,p_)*p_to_eint + &
                                half*(^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)
                        end if
                    else
                        w(ix^D,e_) = w(ix^D,p_)*p_to_eint + &
                            half*(^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)
                    end if
                end if
            {end do\}

        end subroutine p_to_e

        !> LTE conserved -> primitive conversion.
        !>
        !> On entry: rho_ = density, m_ = momentum, e_ = total energy.
        !> On exit:  rho_ = density (unchanged), m_ = velocity, p_ = pressure.
        !>
        !> Pressure is computed energy-consistently from actual eint via EoS
        !> table lookups (T and ne/nH). Cannot use stored Ne_/Te_ because
        !> they may be stale after AMR prolongation/coarsening.
        subroutine hd_to_primitive_LTE(ixI^L, ixO^L, w, x)
            use mod_global_parameters
            use mod_dust, only: dust_to_primitive
            use mod_eos, only: T_from_nH_eint, y_from_nH_eint, &
                saha_T_from_nH_eint, saha_y_from_nH_T, p_nH_from_eint
            integer, intent(in)             :: ixI^L, ixO^L
            double precision, intent(inout) :: w(ixI^S, nw)
            double precision, intent(in)    :: x(ixI^S, 1:ndim)

            double precision                :: inv_rho
            double precision                :: nH(ixI^S)
            double precision                :: log_nH(ixI^S)
            double precision                :: eint_val, eint_in, T_loc, y_loc
            integer :: ix^D

            timeeos0 = MPI_WTIME()

            call eos%get_nH(w, x, ixI^L, ixO^L, nH)

            ! Cache log10(nH) for all cells (used by table lookups for IonE)
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
                        ! Energy-consistent pressure from actual eint via EoS tables.
                        ! Cannot use stored Ne_/Te_ because they may be stale after
                        ! AMR prolongation/coarsening (nonlinear EoS breaks averaging).
                        eint_val = w(ix^D,e_) - half*w(ix^D,rho_)*(^C&w(ix^D,m^C_)**2+)
                        ! Floor eint to prevent unphysical values
                        if (eos%method /= 'analytic') then
                            eint_val = max(eint_val, nH(ix^D) * 10.0d0**eos%T%var2_min)
                        end if
                        eint_val = max(eint_val, smalldouble)
                        if (eint_val * inv_rho > eos%eint_rho_FI_threshold) then
                            ! FI bypass: p = (gamma-1)*(eint - eion*nH)
                            w(ix^D,p_) = eos%gamma_minus_1 &
                                * (eint_val - eos%eion_per_nH * nH(ix^D))
                        else
                            ! Ionisation zone: single p/nH lookup
                            eint_in = dlog10(eint_val) - log_nH(ix^D)
                            w(ix^D,p_) = nH(ix^D) &
                                * p_nH_from_eint(log_nH(ix^D), eint_in)
                        end if
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

            timeeos_conv=timeeos_conv+(MPI_WTIME()-timeeos0)

        end subroutine hd_to_primitive_LTE

        !> Calculate thermal pressure within ixO^L.
        !> For energy runs delegates to eos%get_thermal_pressure; for no-energy
        !> (isothermal) uses the adiabatic relation p = adiab * rho^gamma.
        subroutine hd_get_pthermal(w, x, ixI^L, ixO^L, pth)
        use mod_global_parameters
        use mod_physics, only: phys_get_ei
        use mod_usr_methods, only: usr_set_pthermal
        use mod_small_values, only: trace_small_values

        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S, 1:nw)
        double precision, intent(in) :: x(ixI^S, 1:ndim)
        double precision, intent(out):: pth(ixI^S)
        integer                      :: iw, ix^D
        double precision :: nH(ixI^S)

        if (hd_energy) then
            if (eos%eos_type == 'LTE') then
                ! LTE: p = nH * (1 + He + ne/nH) * T from stored state
                call eos%get_nH(w, x, ixI^L, ixO^L, nH)
                pth(ixO^S) = nH(ixO^S) * (1.0d0 + eos%He_abundance &
                    + (w(ixO^S,iw_ne) / nH(ixO^S))) * w(ixO^S,iw_te)
            else
                ! FI: p = (gamma-1) * (e - KE)
                pth(ixO^S) = eos%gamma_minus_1 * phys_get_ei(w, ixI^L, ixO^L)
            end if
        else
            if (.not. associated(usr_set_pthermal)) then
                pth(ixO^S) = hd_adiab * w(ixO^S, rho_)**eos%gamma
            else
                call usr_set_pthermal(w,x,ixI^L,ixO^L,pth)
            end if
        end if

        if (fix_small_values) then
            {do ix^DB= ixO^LIM^DB\}
            if(pth(ix^D)<small_pressure) then
                pth(ix^D)=small_pressure
            endif
            {enddo^D&\}
        else if (check_small_values) then
            {do ix^DB= ixO^LIM^DB\}
            if(pth(ix^D)<small_pressure) then
                write(*,*) "Error: small value of gas pressure",pth(ix^D),&
                    " encountered when call hd_get_pthermal"
                write(*,*) "Iteration: ", it, " Time: ", global_time
                write(*,*) "Location: ", x(ix^D,:)
                write(*,*) "Cell number: ", ix^D
                do iw=1,nw
                write(*,*) trim(cons_wnames(iw)),": ",w(ix^D,iw)
                end do
                if(trace_small_values) write(*,*) dsqrt(pth(ix^D)-bigdouble)
                write(*,*) "Saving status at the previous time step"
                crash=.true.
            end if
            {enddo^D&\}
        end if

        end subroutine hd_get_pthermal

        !> Sound speed squared for FI (fully ionized / constant gamma) EoS.
        !> Expects w in primitive form: w(p_) = pressure, w(rho_) = density.
        subroutine hd_get_csound2_FI(w, x, ixI^L, ixO^L, cs2)
            use mod_global_parameters
            integer, intent(in)             :: ixI^L, ixO^L
            double precision, intent(in)    :: w(ixI^S, nw)
            double precision, intent(in)    :: x(ixI^S, 1:ndim)
            double precision, intent(out)   :: cs2(ixI^S)

            timeeos0 = MPI_WTIME()

            cs2(ixO^S) = eos%gamma * w(ixO^S, p_) / w(ixO^S, rho_)

            timeeos_csound = timeeos_csound + (MPI_WTIME()-timeeos0)

        end subroutine hd_get_csound2_FI

        !> Sound speed squared for LTE+IonE EoS.
        !> Delegates Gamma_1 computation to hd_get_gamma1_LTE, then cs2 = Gamma_1 * p/rho.
        subroutine hd_get_csound2_LTE(w, x, ixI^L, ixO^L, cs2)
            use mod_global_parameters
            integer, intent(in)             :: ixI^L, ixO^L
            double precision, intent(in)    :: w(ixI^S, nw)
            double precision, intent(in)    :: x(ixI^S, 1:ndim)
            double precision, intent(out)   :: cs2(ixI^S)
            integer :: ix^D

            timeeos0 = MPI_WTIME()

            call hd_get_gamma1_LTE(w, x, ixI^L, ixO^L, cs2)
            {do ix^DB=ixOmin^DB,ixOmax^DB\}
                cs2(ix^D) = cs2(ix^D) * w(ix^D, p_) / w(ix^D, rho_)
            {end do\}

            timeeos_csound = timeeos_csound + (MPI_WTIME()-timeeos0)

        end subroutine hd_get_csound2_LTE

        !> Return constant gamma for FI EoS
        subroutine hd_get_gamma1_FI(w, x, ixI^L, ixO^L, gamma1)
            use mod_global_parameters
            integer, intent(in)             :: ixI^L, ixO^L
            double precision, intent(in)    :: w(ixI^S, nw)
            double precision, intent(in)    :: x(ixI^S, 1:ndim)
            double precision, intent(out)   :: gamma1(ixI^S)

            gamma1(ixO^S) = eos%gamma

        end subroutine hd_get_gamma1_FI

        !> Return effective adiabatic index for LTE+IonE EoS.
        !> Dispatches on gamma1_method and eos%method.
        subroutine hd_get_gamma1_LTE(w, x, ixI^L, ixO^L, gamma1)
            use mod_global_parameters
            use mod_eos, only: gamma1_from_nH_p, gamma1_from_nH_T_analytic
            integer, intent(in)             :: ixI^L, ixO^L
            double precision, intent(in)    :: w(ixI^S, nw)
            double precision, intent(in)    :: x(ixI^S, 1:ndim)
            double precision, intent(out)   :: gamma1(ixI^S)

            double precision :: nH_val, p_over_rho
            integer :: ix^D

            if (eos%gamma1_method == 'constant') then
                gamma1(ixO^S) = eos%gamma
                return
            end if

            {do ix^DB=ixOmin^DB,ixOmax^DB\}
                p_over_rho = w(ix^D, p_) / w(ix^D, rho_)
                if (p_over_rho > eos%p_rho_FI_threshold) then
                    gamma1(ix^D) = eos%gamma
                else
                    nH_val = w(ix^D, rho_) / eos%nH2rhoFactor
                    if (eos%gamma1_method == 'effective') then
                        if (iw_te > 0 .and. w(ix^D,iw_te) > 0.0d0 .and. &
                            iw_ne > 0) then
                            block
                                double precision :: y_g, eint_g
                                y_g = w(ix^D,iw_ne) / nH_val
                                eint_g = eos%inv_gamma_minus_1 * (1.0d0 + y_g) * nH_val &
                                    * w(ix^D,iw_te)
                                if (eos%ionE) eint_g = eint_g &
                                    + y_g * eos%eion_per_nH * nH_val
                                if (eint_g > 0.0d0) then
                                    gamma1(ix^D) = 1.0d0 + w(ix^D,p_) / eint_g
                                else
                                    gamma1(ix^D) = eos%gamma
                                end if
                            end block
                        else
                            gamma1(ix^D) = eos%gamma
                        end if
                    else if (eos%method == 'analytic') then
                        if (iw_te > 0 .and. w(ix^D,iw_te) > 0.0d0) then
                            gamma1(ix^D) = gamma1_from_nH_T_analytic(nH_val, w(ix^D,iw_te))
                        else
                            gamma1(ix^D) = eos%gamma
                        end if
                    else
                        gamma1(ix^D) = gamma1_from_nH_p(dlog10(nH_val), &
                            dlog10(w(ix^D, p_) / nH_val))
                    end if
                end if
            {end do\}

        end subroutine hd_get_gamma1_LTE

        subroutine Rfactor_from_constant_ionization(w,x,ixI^L,ixO^L,Rfactor)
            use mod_global_parameters
            integer, intent(in) :: ixI^L, ixO^L
            double precision, intent(in) :: w(ixI^S,1:nw)
            double precision, intent(in) :: x(ixI^S,1:ndim)
            double precision, intent(out):: Rfactor(ixI^S)

            Rfactor(ixO^S)=RR

        end subroutine Rfactor_from_constant_ionization

        !> Rfactor from LTE EOS: R = (1 + He + ne/nH) / nH2rhoFactor
        !> Uses stored Ne_ wextra field from update_eos_LTE.
        subroutine Rfactor_from_LTE(w,x,ixI^L,ixO^L,Rfactor)
            use mod_global_parameters
            integer, intent(in) :: ixI^L, ixO^L
            double precision, intent(in) :: w(ixI^S,1:nw)
            double precision, intent(in) :: x(ixI^S,1:ndim)
            double precision, intent(out):: Rfactor(ixI^S)

            double precision :: y_nH
            integer :: ix^D

            {do ix^DB=ixOmin^DB,ixOmax^DB\}
              y_nH = w(ix^D, iw_ne) * eos%nH2rhoFactor / w(ix^D, iw_rho)
              Rfactor(ix^D) = (1.0d0 + eos%He_abundance + y_nH) / eos%nH2rhoFactor
            {end do\}

        end subroutine Rfactor_from_LTE

        !> Convert conserved (rho, rho*v, E) to prolong form (rho, v, T).
        !> T is stored in the e_ slot. Interpolation in this space avoids
        !> Jensen's inequality across the ionisation plateau.
        subroutine hd_to_prolong_LTE(ixI^L, ixO^L, w, x)
        use mod_eos, only: T_from_nH_eint
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        double precision :: inv_rho, eint_val, nH_val, log_nH
        integer :: ix^D

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            inv_rho = 1.d0 / w(ix^D, rho_)
            nH_val = w(ix^D, rho_) / eos%nH2rhoFactor
            log_nH = dlog10(nH_val)

            ! Convert momentum to velocity
            w(ix^D, mom(1)) = w(ix^D, mom(1)) * inv_rho

            ! Compute eint = E - 0.5*rho*v^2
            eint_val = w(ix^D, e_) - half * w(ix^D, rho_) * w(ix^D, mom(1))**2
            eint_val = max(eint_val, smalldouble)

            ! Convert eint to T via EoS table
            if (eint_val * inv_rho > eos%eint_rho_FI_threshold) then
            ! FI bypass: T = (gamma-1)*(eint - eion*nH) / (nH * n_per_nH_FI)
            w(ix^D, e_) = eos%gamma_minus_1 &
                * (eint_val - eos%eion_per_nH * nH_val) &
                / (nH_val * eos%n_per_nH_FI)
            else
            ! Ionisation zone: T from table
            w(ix^D, e_) = T_from_nH_eint(log_nH, &
                dlog10(eint_val) - log_nH)
            end if
        {end do\}

        end subroutine hd_to_prolong_LTE

        !> Convert prolong form (rho, v, T) back to conserved (rho, rho*v, E).
        !> T is read from the e_ slot. Uses eint_nH_from_T for back-conversion.
        subroutine hd_from_prolong_LTE(ixI^L, ixO^L, w, x)
        use mod_eos, only: eint_nH_from_T
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        double precision :: T_val, eint_val, T_FI, nH_val, log_nH
        integer :: ix^D

        T_FI = (eos%eint_rho_FI_threshold &
            * eos%nH2rhoFactor - eos%eion_per_nH) &
            * eos%gamma_minus_1 / eos%n_per_nH_FI

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            T_val = w(ix^D, e_)  ! T stored in e_ slot
            nH_val = w(ix^D, rho_) / eos%nH2rhoFactor
            log_nH = dlog10(nH_val)

            if (T_val > T_FI) then
            ! FI: eint = nH*(n_per_nH*T/(gamma-1) + eion)
            eint_val = nH_val &
                * (eos%n_per_nH_FI * T_val * eos%inv_gamma_minus_1 &
                + eos%eion_per_nH)
            else
            ! Ionisation zone: eint/nH from T table
            eint_val = eint_nH_from_T(log_nH, &
                dlog10(max(T_val, 10.0d0**eos%eint_from_T%var2_min))) &
                * nH_val
            end if

            ! E = eint + 0.5*rho*v^2
            w(ix^D, e_) = eint_val + half * w(ix^D, rho_) * w(ix^D, mom(1))**2
            ! Convert velocity to momentum
            w(ix^D, mom(1)) = w(ix^D, rho_) * w(ix^D, mom(1))
        {end do\}

        end subroutine hd_from_prolong_LTE

    end module mod_hd_eos
!> Needs a line after to pass the preprocesor