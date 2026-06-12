!=============================================================================
!> Shared EoS accessors (EoS-type-agnostic) for the eos% family.
!>
!> Carved out of mod_eos.t (the per-type split). These operate on `eos` and the
!> conserved/primitive state without committing to FI / LTE / PI, so they are
!> the common dependency that mod_eos_FI / mod_eos_LTE / mod_eos_PI and the
!> orchestrating mod_eos all draw on -- keeping the type-specific modules free
!> of any back-dependency on mod_eos (no circular use). mod_eos re-exports the
!> public names so existing `use mod_eos` callers are unaffected.
!=============================================================================
module mod_eos_shared_functions
    use mod_global_parameters
    use mod_eos_container

    implicit none
    private

    public :: get_rho, get_nH, get_ne_nH
    public :: get_temperature_from_etot
    public :: eos_get_log_T_floor

contains

    !> Mass density (code units).
    subroutine get_rho(w,x,ixI^L,ixO^L,rho)
        integer, intent(in)           :: ixI^L, ixO^L
        double precision, intent(in)  :: w(ixI^S,1:nw)
        double precision, intent(in)  :: x(ixI^S,1:ndim)
        double precision, intent(out) :: rho(ixI^S)

        rho(ixO^S) = w(ixO^S,iw_rho)

    end subroutine get_rho

    !> Hydrogen number density: nH = rho / nH2rhoFactor.
    subroutine get_nH(w,x,ixI^L,ixO^L,nH)
        integer, intent(in)           :: ixI^L, ixO^L
        double precision, intent(in)  :: w(ixI^S,1:nw)
        double precision, intent(in)  :: x(ixI^S,1:ndim)
        double precision, intent(out) :: nH(ixI^S)

        nH(ixO^S) = w(ixO^S,iw_rho) / eos%nH2rhoFactor

    end subroutine get_nH

    !> Return electron and hydrogen number densities in code units.
    !> For LTE (iw_ne allocated): ne from Saha EoS, nH from rho/nH2rhoFactor.
    !> For FI  (iw_ne not allocated): ne = nH * neOnH_FI (full ionisation).
    subroutine get_ne_nH(ixI^L, ixO^L, w, ne, nH)
        integer, intent(in)           :: ixI^L, ixO^L
        double precision, intent(in)  :: w(ixI^S, 1:nw)
        double precision, intent(out) :: ne(ixI^S), nH(ixI^S)

        nH(ixO^S) = w(ixO^S, iw_rho) / eos%nH2rhoFactor
        if (iw_ne > 0) then
            ne(ixO^S) = w(ixO^S, iw_ne)
        else
            ne(ixO^S) = nH(ixO^S) * eos%neOnH_FI
        end if
    end subroutine get_ne_nH

    !> Temperature from the TOTAL energy: strip kinetic (+magnetic) energy via
    !> phys_e_to_ei, then delegate to the type-specific eos%get_temperature_from_eint.
    subroutine get_temperature_from_etot(w, x, ixI^L, ixO^L, res)
        use mod_physics
        integer, intent(in)             :: ixI^L,ixO^L
        double precision, intent(in)    :: x(ixI^S,1:ndim)
        double precision, intent(in)    :: w(ixI^S,1:nw)
        double precision, intent(out)   :: res(ixI^S)
        double precision :: wlocal(ixI^S,1:nw)

        !> No timing here: get_temperature_from_eint is already timed.
        !> The array copy and e_to_ei subtraction are trivial.
        wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
        call phys_e_to_ei(ixI^L, ixO^L, wlocal, x)
        call eos%get_temperature_from_eint(wlocal, x, ixI^L, ixO^L, res)

    end subroutine get_temperature_from_etot

    !> Lower log10(T) bound for the (rho,T) inverse table; FI (no tables) floors
    !> at log10(smalldouble). Used by the TC fluid-port log_T_floor.
    double precision function eos_get_log_T_floor() result(log_T_min)
        if (eos%method_id == EOS_ENTROPY .and. allocated(eos%eintT%table)) then
            log_T_min = eos%eintT%var2_min
        else if (allocated(eos%eint_from_T%table)) then
            log_T_min = eos%eint_from_T%var2_min
        else
            log_T_min = dlog10(smalldouble)
        end if
    end function eos_get_log_T_floor

end module mod_eos_shared_functions
