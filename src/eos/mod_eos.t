!> Equation of state for AMRVAC, handled through a single eos_container object.
!>
!> Two gas models (eos_type):
!>   FI  -- fully ionised ideal gas (constant gamma, mu).
!>   LTE -- local thermodynamic equilibrium with partial ionisation of H (+He),
!>          via one of three interchangeable methods (eos_method):
!>            tables   -- bilinear/PCHIP lookup of pre-built tables.
!>            entropy  -- bicubic-Hermite reconstruction (see mod_eos_entropy).
!>            analytic -- on-the-fly H-only Saha solve (see saha_* routines).
!>          Tables are generated offline and read as binary; their axes are
!>          (log10 nH, log10 eint/nH) in code units after the unit shift applied
!>          in eos_finalise. ionE optionally folds ionisation energy into eint.
!>
!> The code-unit normalisation is kept so an FI atmosphere recovers R = 1;
!> update_eos_LTE relies on this. update_eos refreshes the cached EoS fields
!> (Te, ne) once per RK substep, and for boundary cells via update_eos_4_bc;
!> call eos%update_eos manually before any extra output.
!>
!> Jack Jenkins (27.01.2026); H(+He) LTE tables by Chris Osborne, in
!> consultation with Damien Przybylski (MURaM).

module mod_eos
    use mod_global_parameters
    use mod_eos_container
    use mod_eos_interp
    use mod_eos_saha
    use mod_eos_tables
    use mod_timing
    use mod_comm_lib, only: mpistop

    implicit none
    private

    character(len=std_len) :: AMRVAC_DIR

    !> The EoS state object (eos) and the cached-log10(nH) wextra index
    !> (iw_log_nH) now live in mod_eos_container so every EoS sub-module can
    !> reach them; re-export here so existing `use mod_eos` callers still see them.
    public :: eos, iw_log_nH

    !> Lifecycle (called from amrvac.t)
    public :: eos_init, eos_finalise, prepare_eos_w_fields
    !> Scalar EoS kernels used by the hd/mhd EoS layer and sub-physics ports
    public :: y_from_nH_eint, T_from_nH_eint, p2eint_from_nH_p
    public :: gamma1_from_nH_p, p_nH_from_eint
    public :: eint_from_p_bisect, eint_nH_from_T
    public :: eos_get_log_T_floor, eos_get_eintT_grid
    !> Analytic Saha entry points (analytic method)
    public :: saha_eint_from_nH_T, saha_T_from_nH_eint, saha_p_to_T
    public :: gamma1_from_nH_T_analytic
    !> Partial-ionization temperature (Leenaarts-style; hd only)
    public :: Rfactor_from_PI_temperature, update_PI_temperature
    !> Fast TC temperature kernel + electron/H density getter
    public :: get_temperature_from_eint_LTE_fast
    public :: get_ne_nH

contains
    !> Read this module"s parameters from a file
    subroutine eos_read_params(files)
        use mod_eos_entropy, only: entropy_set_bilinear
        character(len=*), intent(in) :: files(:)
        integer                      :: n

        !> Default values
        character(len=std_len)  :: eos_type
        logical :: ionE
        character(len=std_len) :: table_location
        logical :: table_check
        double precision :: He_abundance
        double precision :: gamma
        double precision :: gamma_minus_1
        double precision :: inv_gamma
        double precision :: inv_gamma_minus_1
        double precision :: nH2rhoFactor

        logical :: disable_FI_bypass
        character(len=20) :: eos_method, gamma1_method, lte_h_inversion, p2eint_method
        character(len=20) :: entropy_interp

        namelist /eos_list/ eos_type, table_location, ionE
        namelist /eos_list/ table_check, He_abundance, gamma
        namelist /eos_list/ disable_FI_bypass
        namelist /eos_list/ eos_method, gamma1_method, lte_h_inversion, p2eint_method
        namelist /eos_list/ entropy_interp

        eos_type = 'FI'
        eos_method = 'tables'
        gamma1_method = 'exact'
        lte_h_inversion = 'bisect'
        p2eint_method = 'table'
        entropy_interp = 'bicubic'
        disable_FI_bypass = .false.
        call get_environment_variable("AMRVAC_DIR", AMRVAC_DIR)
        table_location = trim(AMRVAC_DIR)//"/src/tables/eos_tables/uniform256/"
        ionE = .false.
        table_check = .false.
        He_abundance = 0.1d0
        gamma = 5.0d0/3.0d0

        do n = 1, size(files)
            open(unitpar, file=trim(files(n)), status="old")
            read(unitpar, eos_list, end=111)
    111     close(unitpar)
        end do

        eos%eos_type = eos_type
        if(mype==0) write(*,*) "EoS type: "//trim(eos%eos_type)
        eos%ionE = ionE
        if ((eos%eos_type == 'FI') .and. eos%ionE) then
            eos%ionE = .false.
            if(mype==0) write(*,*) "WARNING: eos_type = 'FI' requires ionE = .false."
        endif
        eos%table_location = table_location
        if((mype==0) .and. (eos%eos_type == 'LTE')) write(*,*) "EoS table location: "//trim(eos%table_location)
        eos%table_check = table_check
        eos%He_abundance = He_abundance
        eos%gamma = gamma
        eos%gamma_minus_1 = eos%gamma - 1.0d0
        eos%inv_gamma = 1.0d0 / eos%gamma
        if (abs(eos%gamma_minus_1) > 1.0d-12) then
            eos%inv_gamma_minus_1 = 1.0d0 / eos%gamma_minus_1
        else
            eos%inv_gamma_minus_1 = 0.0d0
        end if
        eos%nH2rhoFactor = 1.0d0 !> Assume the default is 1.0 i.e., FI
        eos%disable_FI_bypass = disable_FI_bypass
        eos%method = eos_method
        select case (eos%method)
        case ('analytic'); eos%method_id = EOS_ANALYTIC
        case ('entropy');  eos%method_id = EOS_ENTROPY
        case default;      eos%method_id = EOS_TABLES
        end select
        eos%gamma1_method = gamma1_method
        eos%p2eint_method = p2eint_method
        eos%inversion = lte_h_inversion

        !> Entropy-method interpolation order. 'bicubic' (default) uses the
        !> 4-derivative Hermite kernel. 'bilinear' ignores the derivative
        !> tables and uses plain corner-value bilinear -- faster per call,
        !> less accurate. Bench-only knob; the bicubic path remains the
        !> production default.
        if (entropy_interp == 'bilinear') then
            call entropy_set_bilinear(.true.)
            if (mype == 0) write(*,*) &
                'NOTE: entropy_interp=bilinear -- bicubic Hermite disabled, '// &
                'derivative tables ignored.'
        else if (entropy_interp == 'bicubic') then
            call entropy_set_bilinear(.false.)
        else
            call mpistop("Unknown entropy_interp: "//trim(entropy_interp)// &
                " (expected 'bicubic' or 'bilinear')")
        end if

        ! Validate method combinations
        if (eos%eos_type == 'FI' .and. eos%method_id == EOS_ANALYTIC) then
            call mpistop("eos_method='analytic' requires eos_type='LTE'")
        end if
        if (eos%method_id == EOS_ANALYTIC .and. eos%He_abundance > 0.0d0) then
            if (mype == 0) write(*,*) &
                "WARNING: eos_method='analytic' is H-only. He_abundance=", &
                eos%He_abundance, " used in pressure law but not in Saha ionization."
        end if
        if (eos%method /= 'analytic' .and. eos%method /= 'tables' &
            .and. eos%method_id /= EOS_ENTROPY) then
            call mpistop("Unknown eos_method: "//trim(eos%method)// &
                " (expected 'tables', 'analytic', or 'entropy')")
        end if
        if (eos%gamma1_method /= 'exact' .and. eos%gamma1_method /= 'effective' &
            .and. eos%gamma1_method /= 'constant') then
            call mpistop("Unknown gamma1_method: "//trim(eos%gamma1_method)// &
                " (expected 'exact', 'effective', or 'constant')")
        end if
        if (eos%inversion /= 'bisect' .and. eos%inversion /= 'newton') then
            call mpistop("Unknown lte_h_inversion: "//trim(eos%inversion)// &
                " (expected 'bisect' or 'newton')")
        end if
        ! effective gamma1 requires table-based EoS (needs T, neOnH tables to build gamma_eff table)
        if (eos%method_id == EOS_ANALYTIC .and. eos%gamma1_method == 'effective') then
            if (mype == 0) write(*,*) &
                "WARNING: gamma1_method='effective' not applicable with analytic EoS."
            if (mype == 0) write(*,*) &
                "         Falling back to gamma1_method='exact' (analytic 2D table)."
            eos%gamma1_method = 'exact'
        end if
        ! analytic with ionE=false is just ideal gas -- warn
        if (eos%method_id == EOS_ANALYTIC .and. .not. eos%ionE) then
            if (mype == 0) then
                write(*,*) "WARNING: eos_method='analytic' without ionE"
                write(*,*) "  reduces to ideal gas. Consider eos_type='FI'."
            end if
        end if
        ! effective gamma1: table built with gamma_eff = 1 + p/eint approximation
        if (eos%eos_type == 'LTE' .and. eos%ionE .and. &
            eos%gamma1_method == 'effective') then
            if (mype == 0) write(*,*) &
                "NOTE: gamma1_method='effective' builds gamma_eff = 1+p/eint table."
        end if

    end subroutine eos_read_params

    !> Phase 'create' (before units are known): allocate the EoS object, read
    !> its &eos_list parameters, and load the raw tables for the chosen method.
    !> The getters wired here must be live before (m)hd_link_eos runs.
    subroutine eos_init()
        allocate(eos)
        call eos_read_params(par_files)
        if (eos%eos_type == 'LTE') call eos_load_method_tables()

        eos%get_rho   => get_rho
        eos%get_nH    => get_nH
        eos%get_ne_nH => get_ne_nH
        !> FI temperature-from-pressure: wired here (not in eos_finalise) so it is
        !> non-null before (m)hd_link_eos captures it into the radiation fluid
        !> object. The LTE variant is a later pass.
        if (eos%eos_type == 'FI') &
            eos%get_temperature_from_pressure => get_temperature_from_pressure_FI
    end subroutine eos_init

    !> Load the binary LTE tables for the chosen method (no files for analytic).
    subroutine eos_load_method_tables()
        integer :: iq, id
        character(len=8), parameter :: quantity(6) = &
            [character(8):: 'neOnH', 'Tfwd', 'pfwd', 'eintP', 'g1p', 'eintT']
        character(len=3), parameter :: deriv(4) = &
            [character(3):: '', '_x', '_y', '_xy']

        select case (eos%method_id)
        case (EOS_ANALYTIC)
            if (mype == 0) then
                write(*,*) 'EoS method: analytic (H-only Saha)'
                write(*,*) '  inversion:     ', trim(eos%inversion)
                write(*,*) '  gamma1_method: ', trim(eos%gamma1_method)
            end if
        case (EOS_ENTROPY)
            !> Six quantities, each a value table plus its three derivative tables,
            !> on the forward (nH, eint/nH), inverse-p (nH, p/nH) and inverse-T
            !> (nH, T) grids. Every runtime query is one bicubic-Hermite evaluation;
            !> see mod_eos_entropy and entropy/generate_all_tables.py.
            if (mype == 0) write(*,*) &
                'EoS method: entropy (bicubic Hermite, 4-derivative corners, no closure)'
            do iq = 1, size(quantity)
                do id = 1, size(deriv)
                    call load_lte_tables(trim(quantity(iq))//trim(deriv(id)))
                end do
            end do
        case default   ! EOS_TABLES
            call load_lte_tables("T")
            call load_lte_tables("neOnH")
            if (eos%ionE) then
                call load_lte_tables("p2eint")
                call try_load_lte_tables("gamma1")
                call try_load_lte_tables("eint_from_T")
            end if
        end select
    end subroutine eos_load_method_tables

    !> Phase 'commit' (after units are known): finalise the dispatch for the
    !> loaded physics -- wire the method's runtime pointers and finalise its
    !> tables. The cross-module wiring that follows (into the ghost-cell update
    !> and the physics source terms) is done by the caller in amrvac.t, so the
    !> EoS never reaches back into those modules.
    subroutine eos_finalise()
        if (allocated(eos)) then
            if (eos%eos_type == 'LTE') then
                eos%update_eos => update_eos_LTE
                ! eos%get_thermal_pressure is set by (m)hd_link_eos
                eos%get_temperature_from_eint => get_temperature_from_eint_LTE
                eos%get_Te => get_Te_LTE
                call eos_finalise_lte_tables()
            else if (eos%eos_type == 'FI') then
                eos%update_eos => update_eos_FI   !> nothing cached for ideal gas
                eos%get_temperature_from_eint => get_temperature_from_eint_FI
                eos%get_Te => get_Te_FI
                !> Fully-ionised particle counts (2 + 3*A_He per H; ne/nH = 1 + 2*A_He)
                eos%n_per_nH_FI = 2.0d0 + 3.0d0 * eos%He_abundance
                eos%neOnH_FI    = 1.0d0 + 2.0d0 * eos%He_abundance
            else
                call mpistop('eos_finalise: eos_type '//trim(eos%eos_type)//' not recognised')
            endif

            eos%get_temperature_from_etot => get_temperature_from_etot
        else
            call mpistop('eos_finalise: eos not allocated')
        endif
    end subroutine eos_finalise

    !> Finalise the loaded LTE tables for the selected method: shift axes to
    !> code units, build any derived tables, validate. The table mechanics live
    !> in mod_eos_tables; this routine only orders them for the chosen method.
    subroutine eos_finalise_lte_tables()
        if (eos%method_id == EOS_ANALYTIC) then
            !> Analytical H-only Saha: no table unit conversion needed.
            !> Build Gamma1 table if 'exact' mode requested.
            if (eos%gamma1_method == 'exact') then
                call build_gamma1_analytic_table()
            end if

            !> Precompute FI bypass constants (same as table path)
            call precompute_FI_bypass_constants()

        else if (eos%method_id == EOS_ENTROPY) then
            !> Entropy finalise: shift each loaded table's axes from CGS
            !> to code units, then build its lookup-support arrays. No aux
            !> builds -- every runtime quantity is a single lookup. Forward
            !> (axis2 = log eint/nH) and inverse-p (axis2 = log p/nH) share
            !> the erg/particle unit (shift_axis_to_code); inverse-T uses
            !> log T (shift_axis_to_code_T).
            call shift_axis_to_code(eos%neOnH)
            call shift_axis_to_code(eos%neOnH_x)
            call shift_axis_to_code(eos%neOnH_y)
            call shift_axis_to_code(eos%neOnH_xy)
            call shift_axis_to_code(eos%Tfwd)
            call shift_axis_to_code(eos%Tfwd_x)
            call shift_axis_to_code(eos%Tfwd_y)
            call shift_axis_to_code(eos%Tfwd_xy)
            call shift_axis_to_code(eos%pfwd)
            call shift_axis_to_code(eos%pfwd_x)
            call shift_axis_to_code(eos%pfwd_y)
            call shift_axis_to_code(eos%pfwd_xy)
            call shift_axis_to_code(eos%eintP)
            call shift_axis_to_code(eos%eintP_x)
            call shift_axis_to_code(eos%eintP_y)
            call shift_axis_to_code(eos%eintP_xy)
            call shift_axis_to_code(eos%g1p)
            call shift_axis_to_code(eos%g1p_x)
            call shift_axis_to_code(eos%g1p_y)
            call shift_axis_to_code(eos%g1p_xy)
            call shift_axis_to_code_T(eos%eintT)
            call shift_axis_to_code_T(eos%eintT_x)
            call shift_axis_to_code_T(eos%eintT_y)
            call shift_axis_to_code_T(eos%eintT_xy)
            call entropy_table_prepare(eos%neOnH);    call entropy_table_prepare(eos%neOnH_x)
            call entropy_table_prepare(eos%neOnH_y);  call entropy_table_prepare(eos%neOnH_xy)
            call entropy_table_prepare(eos%Tfwd);     call entropy_table_prepare(eos%Tfwd_x)
            call entropy_table_prepare(eos%Tfwd_y);   call entropy_table_prepare(eos%Tfwd_xy)
            call entropy_table_prepare(eos%pfwd);     call entropy_table_prepare(eos%pfwd_x)
            call entropy_table_prepare(eos%pfwd_y);   call entropy_table_prepare(eos%pfwd_xy)
            call entropy_table_prepare(eos%eintP);    call entropy_table_prepare(eos%eintP_x)
            call entropy_table_prepare(eos%eintP_y);  call entropy_table_prepare(eos%eintP_xy)
            call entropy_table_prepare(eos%g1p);      call entropy_table_prepare(eos%g1p_x)
            call entropy_table_prepare(eos%g1p_y);    call entropy_table_prepare(eos%g1p_xy)
            call entropy_table_prepare(eos%eintT);    call entropy_table_prepare(eos%eintT_x)
            call entropy_table_prepare(eos%eintT_y);  call entropy_table_prepare(eos%eintT_xy)
            call precompute_FI_bypass_constants()
            !> The entropy method never loads eos%T, but several callers
            !> read eos%T%var{1,2}_{min,max} (the cold-cell eint floor,
            !> FI-fallback checks). Source them from eos%Tfwd, which is on
            !> the same forward (log nH, log eint/nH) grid by construction.
            eos%T%var2_min = eos%Tfwd%var2_min
            eos%T%var2_max = eos%Tfwd%var2_max
            eos%T%var1_min = eos%Tfwd%var1_min
            eos%T%var1_max = eos%Tfwd%var1_max
            if (mype == 0) write(*,'(A)') &
                ' EoS entropy: forward + inverse tables loaded; no aux builds.'

        else
            !> Table-based path: unit conversion and table building
            eos%T%table = eos%T%table - dlog10(unit_temperature)

            eos%T%var1_min = eos%T%var1_min - dlog10(unit_numberdensity)
            eos%T%var1_max = eos%T%var1_max - dlog10(unit_numberdensity)
            eos%T%var2_min = eos%T%var2_min - dlog10(unit_pressure/unit_numberdensity)
            eos%T%var2_max = eos%T%var2_max - dlog10(unit_pressure/unit_numberdensity)
            if (.not. eos%T%is_uniform) then
                eos%T%var1_nodes = eos%T%var1_nodes - dlog10(unit_numberdensity)
                eos%T%var2_nodes = eos%T%var2_nodes - dlog10(unit_pressure/unit_numberdensity)
            end if

            eos%neOnH%var1_min = eos%neOnH%var1_min - dlog10(unit_numberdensity)
            eos%neOnH%var1_max = eos%neOnH%var1_max - dlog10(unit_numberdensity)
            eos%neOnH%var2_min = eos%neOnH%var2_min - dlog10(unit_pressure/unit_numberdensity)
            eos%neOnH%var2_max = eos%neOnH%var2_max - dlog10(unit_pressure/unit_numberdensity)
            if (.not. eos%neOnH%is_uniform) then
                eos%neOnH%var1_nodes = eos%neOnH%var1_nodes - dlog10(unit_numberdensity)
                eos%neOnH%var2_nodes = eos%neOnH%var2_nodes - dlog10(unit_pressure/unit_numberdensity)
            end if

            !> Ensure axis-node arrays exist for T and neOnH so that
            !> the build_*_table routines below can use a single
            !> code path (var1_nodes(i) / var2_nodes(j)) regardless
            !> of whether the source tables are uniform or adaptive.
            call ensure_axis_nodes(eos%T)
            call ensure_axis_nodes(eos%neOnH)

            if (eos%ionE) then
                if (allocated(eos%p2eint%table)) then
                    eos%p2eint%var1_min = eos%p2eint%var1_min - dlog10(unit_numberdensity)
                    eos%p2eint%var1_max = eos%p2eint%var1_max - dlog10(unit_numberdensity)
                    eos%p2eint%var2_min = eos%p2eint%var2_min - dlog10(unit_pressure/unit_numberdensity)
                    eos%p2eint%var2_max = eos%p2eint%var2_max - dlog10(unit_pressure/unit_numberdensity)
                    if (.not. eos%p2eint%is_uniform) then
                        eos%p2eint%var1_nodes = eos%p2eint%var1_nodes - dlog10(unit_numberdensity)
                        eos%p2eint%var2_nodes = eos%p2eint%var2_nodes - dlog10(unit_pressure/unit_numberdensity)
                    end if
                    call ensure_axis_nodes(eos%p2eint)
                end if

                !> Build derived tables from loaded T and neOnH.
                !> For tables that have a loaded adaptive version
                !> on disk (trailer present, is_uniform == .false.),
                !> *keep the loaded one* -- this enables the hybrid
                !> mode where T+neOnH ship uniform (interleaved
                !> fast path preserved) while gamma1/p2eint/
                !> eint_from_T ship adaptive (gamma_1 peak resolved,
                !> tighter closure). The loaded versions were
                !> built by the same Saha solver as T/neOnH so
                !> they are physically consistent.
                ! Skip rebuild for any derived table that was successfully
                ! loaded from disk (whether uniform or adaptive). This
                ! enables shipping pre-built tables and avoids the
                ! one-time runtime rebuild cost.
                if (.not. allocated(eos%gamma1%table)) then
                    call build_gamma1_table()
                else
                    if (mype == 0) write(*,'(A)') " Using loaded gamma1 table (skip rebuild)"
                    !> gamma1 stored on disk on (log_nH, log_eint/nH) in CGS log10;
                    !> shift axes to code units to match runtime queries.
                    !> Values are dimensionless, no value shift needed.
                    eos%gamma1%var1_min = eos%gamma1%var1_min - dlog10(unit_numberdensity)
                    eos%gamma1%var1_max = eos%gamma1%var1_max - dlog10(unit_numberdensity)
                    eos%gamma1%var2_min = eos%gamma1%var2_min - dlog10(unit_pressure/unit_numberdensity)
                    eos%gamma1%var2_max = eos%gamma1%var2_max - dlog10(unit_pressure/unit_numberdensity)
                    if (.not. eos%gamma1%is_uniform) then
                        eos%gamma1%var1_nodes = eos%gamma1%var1_nodes - dlog10(unit_numberdensity)
                        eos%gamma1%var2_nodes = eos%gamma1%var2_nodes - dlog10(unit_pressure/unit_numberdensity)
                    end if
                    call ensure_axis_nodes(eos%gamma1)
                end if

                !> Order matters: build_gamma1_p_table reads eos%p2eint%var*_nodes,
                !> so p2eint must be allocated first.  For 'tables' it's loaded
                !> from disk earlier; for 'entropy' it has to be built right here.
                if (.not. allocated(eos%p2eint%table)) then
                    call build_p2eint_table()
                else
                    if (mype == 0) write(*,'(A)') " Using loaded p2eint table (skip rebuild)"
                    call ensure_axis_nodes(eos%p2eint)
                end if

                call build_gamma1_p_table()
                call build_log_p_table()
                call build_p_over_nH_table()

                if (.not. allocated(eos%eint_from_T%table)) then
                    call build_eint_from_T_table()
                else
                    if (mype == 0) write(*,'(A)') " Using loaded eint_from_T table (skip rebuild)"
                    !> eint_from_T stored on (log_nH, log_T) in CGS log10
                    !> with values log10(eint/nH) in CGS. Shift axes AND values
                    !> to code units.
                    eos%eint_from_T%var1_min = eos%eint_from_T%var1_min - dlog10(unit_numberdensity)
                    eos%eint_from_T%var1_max = eos%eint_from_T%var1_max - dlog10(unit_numberdensity)
                    eos%eint_from_T%var2_min = eos%eint_from_T%var2_min - dlog10(unit_temperature)
                    eos%eint_from_T%var2_max = eos%eint_from_T%var2_max - dlog10(unit_temperature)
                    eos%eint_from_T%table = eos%eint_from_T%table - dlog10(unit_pressure/unit_numberdensity)
                    if (.not. eos%eint_from_T%is_uniform) then
                        eos%eint_from_T%var1_nodes = eos%eint_from_T%var1_nodes - dlog10(unit_numberdensity)
                        eos%eint_from_T%var2_nodes = eos%eint_from_T%var2_nodes - dlog10(unit_temperature)
                    end if
                    call ensure_axis_nodes(eos%eint_from_T)
                end if
            endif

            call precompute_step_inv(eos%T)
            call precompute_step_inv(eos%neOnH)
            call eos_build_guards(eos%T)
            call eos_build_guards(eos%neOnH)
            if (eos%ionE) then
                call precompute_step_inv(eos%p2eint)
                call precompute_step_inv(eos%gamma1)
                call precompute_step_inv(eos%gamma1_p)
                call precompute_step_inv(eos%eint_from_T)
                call precompute_step_inv(eos%log_p)
                call precompute_step_inv(eos%p_over_nH)
                call eos_build_guards(eos%p2eint)
                call eos_build_guards(eos%gamma1)
                call eos_build_guards(eos%gamma1_p)
                call eos_build_guards(eos%eint_from_T)
                call eos_build_guards(eos%log_p)
                call eos_build_guards(eos%p_over_nH)

                ! Build interleaved Group A table: T, neOnH, p_over_nH
                call build_interleaved_eint_table()
            end if

            ! Defensive check: invariants on every populated table.
            ! Catches silent-failure modes (forgotten guard, mismatched
            ! node array length, non-monotonic nodes) at startup
            ! rather than at first lookup.
            call eos_validate_table(eos%T,           "eos%T")
            call eos_validate_table(eos%neOnH,       "eos%neOnH")
            if (eos%ionE) then
                call eos_validate_table(eos%p2eint,      "eos%p2eint")
                call eos_validate_table(eos%gamma1,      "eos%gamma1")
                call eos_validate_table(eos%gamma1_p,    "eos%gamma1_p")
                call eos_validate_table(eos%eint_from_T, "eos%eint_from_T")
                call eos_validate_table(eos%log_p,       "eos%log_p")
                call eos_validate_table(eos%p_over_nH,   "eos%p_over_nH")
            end if

            call precompute_FI_bypass_constants()

            if (eos%ionE) call verify_eos_round_trips()
        end if  ! method == 'analytic' or 'tables' or 'entropy'
    end subroutine eos_finalise_lte_tables

    subroutine prepare_eos_w_fields()
        integer :: iigrid, igrid

        ! fill eos space of all grids, inc boundaries
        do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
            call eos%update_eos(ixG^LL,ixG^LL,ps(igrid)%w,ps(igrid)%x)
        end do

    end subroutine prepare_eos_w_fields

    subroutine update_eos_LTE(ixI^L, ixO^L, w, x)
        use mod_physics
        !> This routine is called before each RK substep
        integer, intent(in)             :: ixI^L,ixO^L
        double precision, intent(in)    :: x(ixI^S,1:ndim)
        double precision, intent(inout) :: w(ixI^S,1:nw)

        double precision :: wlocal(ixI^S,1:nw)
        double precision :: deltaEfactor(ixI^S)
        double precision :: prev_y(ixI^S), new_y(ixI^S)
        double precision :: pth(ixI^S)
        double precision :: nH_in(ixI^S), nH(ixI^S), eint_in(ixI^S)
        double precision :: Rfactor_FI
        double precision :: yy, eint_nH_floor
        double precision :: time0
        integer :: ix^D
        
        timeeos0 = MPI_WTIME() !> For monitoring cost of eos module

        wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
        call phys_e_to_ei(ixI^L,ixO^L,wlocal,x) !> wlocal now contains internal energy, NOT total energy

        pth(ixO^S) = eos%gamma_minus_1 * wlocal(ixO^S,iw_e) !> pressure from internal energy only - SHOULD ONLY BE USED WHEN IonE IS UNIMPORTANT

        call eos%get_nH(w, x, ixI^L, ixO^L, nH)
        ! Enforce internal energy floor for EoS lookup only.
        ! Prevents NaN from dlog10(eint<0) after strong rarefactions where
        ! kinetic energy can numerically exceed total energy.
        ! The conserved energy w(iw_e) is NOT modified: the physical fluxes
        ! based on the small positive pressure will naturally restore the cell.
        if (eos%method_id == EOS_ANALYTIC) then
            ! Analytic: floor to ~100 K equivalent (no tables loaded)
            ! In code units: eint/nH = T_code / (gamma-1) for neutral gas
            {do ix^DB=ixOmin^DB,ixOmax^DB\}
                eint_nH_floor = nH(ix^D) * eos%inv_gamma_minus_1 * 100.0d0 / unit_temperature
                wlocal(ix^D,iw_e) = max(wlocal(ix^D,iw_e), eint_nH_floor)
            {end do\}
        else
            nH_in(ixO^S) = dlog10(nH(ixO^S))
            {do ix^DB=ixOmin^DB,ixOmax^DB\}
                !> eos%T%var2_min holds the forward log(eint/nH) lower bound for
                !> all table methods (entropy sources it from Tfwd in eos_finalise).
                eint_nH_floor = nH(ix^D) * 10.0d0**eos%T%var2_min
                if (wlocal(ix^D,iw_e) < eint_nH_floor) then
                    wlocal(ix^D,iw_e) = eint_nH_floor
                end if
            {end do\}
            eint_in(ixO^S) = dlog10(wlocal(ixO^S,iw_e)) - nH_in(ixO^S)
        end if

        !> Constant FI Rfactor for the fully-ionised fast path.
        !> Must NOT use eos%get_Rfactor which reads iw_ne (uninitialised at IC).
        Rfactor_FI = eos%n_per_nH_FI / eos%nH2rhoFactor

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            if (wlocal(ix^D,iw_e) / w(ix^D,iw_rho) > eos%eint_rho_FI_threshold) then
                !> Fully ionised: analytical formulae, no table lookups
                new_y(ix^D) = eos%neOnH_FI
                if (eos%ionE) then
                    !> T = (eint - eion*nH) * (gamma-1) / (Rfactor_FI * rho)
                    w(ix^D,iw_te) = eos%gamma_minus_1 &
                        * (wlocal(ix^D,iw_e) - eos%eion_per_nH * nH(ix^D)) &
                        / (Rfactor_FI * w(ix^D,iw_rho))
                else
                    w(ix^D,iw_te) = pth(ix^D) / (nH(ix^D) * (1.0d0 + eos%He_abundance + new_y(ix^D)))
                end if
            else
                !> Ionisation zone
                if (eos%method_id == EOS_ANALYTIC) then
                    !> Analytical Saha: solve quadratic for T and y
                    call saha_T_from_nH_eint(nH(ix^D), &
                        wlocal(ix^D,iw_e) / nH(ix^D), &
                        w(ix^D,iw_te), yy)
                    new_y(ix^D) = yy
                else
                    !> Table lookups: fused T+y for shared index computation
                    if (.not. eos%ionE) then
                        new_y(ix^D) = y_from_nH_eint(nH_in(ix^D),eint_in(ix^D))
                        w(ix^D,iw_te) = pth(ix^D) / (nH(ix^D) * (1.0d0 + eos%He_abundance + new_y(ix^D)))
                    else if (eos%method_id == EOS_ENTROPY) then
                        !> Entropy method dispatches T_and_y_from_nH_eint to the
                        !> biquintic-Hermite forward at any eint; no FI fallback
                        !> needed (and eos%T isn't loaded so var2_max is 0 here).
                        call T_and_y_from_nH_eint(nH_in(ix^D), eint_in(ix^D), &
                            w(ix^D,iw_te), yy)
                        new_y(ix^D) = yy
                    else
                        if (eint_in(ix^D) < eos%T%var2_max) then
                            call T_and_y_from_nH_eint(nH_in(ix^D), eint_in(ix^D), &
                                w(ix^D,iw_te), yy)
                            new_y(ix^D) = yy
                        else
                            !> Above-table fallback with ionisation energy correction
                            new_y(ix^D) = eos%neOnH_FI
                            w(ix^D,iw_te) = eos%gamma_minus_1 &
                                * (wlocal(ix^D,iw_e) - eos%eion_per_nH * nH(ix^D)) &
                                / (Rfactor_FI * w(ix^D,iw_rho))
                        end if
                    end if
                end if
            end if
        {end do\}

        w(ixO^S,iw_ne) = new_y(ixO^S) * nH(ixO^S)

        timeeos_update=timeeos_update+(MPI_WTIME()-timeeos0)

    end subroutine update_eos_LTE

    subroutine update_eos_FI(ixI^L, ixO^L, w, x)
        !> This routine is called before each RK substep
        integer, intent(in)             :: ixI^L,ixO^L
        double precision, intent(in)    :: x(ixI^S,1:ndim)
        double precision, intent(inout) :: w(ixI^S,1:nw)

        !> Nothing needed for the FI variant

    end subroutine update_eos_FI

    subroutine get_rho(w,x,ixI^L,ixO^L,rho)
        use mod_global_parameters
        integer, intent(in)           :: ixI^L, ixO^L
        double precision, intent(in)  :: w(ixI^S,1:nw)
        double precision, intent(in)  :: x(ixI^S,1:ndim)
        double precision, intent(out) :: rho(ixI^S)

        rho(ixO^S) = w(ixO^S,iw_rho)

    end subroutine get_rho

    subroutine get_nH(w,x,ixI^L,ixO^L,nH)
        use mod_global_parameters
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
        use mod_global_parameters
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

    subroutine get_Te_LTE(w,x,ixI^L,ixO^L,T)
        use mod_global_parameters
        integer, intent(in)           :: ixI^L, ixO^L
        double precision, intent(in)  :: w(ixI^S,1:nw)
        double precision, intent(in)  :: x(ixI^S,1:ndim)
        double precision, intent(out) :: T(ixI^S)

        T(ixO^S) = w(ixO^S,iw_te)

    end subroutine get_Te_LTE

    subroutine get_Te_FI(w,x,ixI^L,ixO^L,T)
        use mod_global_parameters
        integer, intent(in)           :: ixI^L, ixO^L
        double precision, intent(in)  :: w(ixI^S,1:nw)
        double precision, intent(in)  :: x(ixI^S,1:ndim)
        double precision, intent(out) :: T(ixI^S)
        double precision :: Rfactor(ixI^S), pth(ixI^S)

        !> No timing here: get_thermal_pressure is already timed.
        !> The division and Rfactor call are trivial.
        call eos%get_thermal_pressure(w, x, ixI^L, ixO^L, pth)
        call eos%get_Rfactor(w,x,ixI^L,ixO^L,Rfactor)
        T(ixO^S) = pth(ixO^S) / (w(ixO^S,iw_rho) * Rfactor(ixO^S))

    end subroutine get_Te_FI

    !> The next four subroutines get the temperature from the energy variable depending on the eos type
    !> These distinctions are necessary primarily for the STS methods (dt -> conserved, set_source -> eint)
    !> See mod_thermal_conduction for more details
    subroutine get_temperature_from_eint_FI(w, x, ixI^L, ixO^L, res)
        !> Assumes input energy is internal energy
        use mod_physics
        integer, intent(in)             :: ixI^L,ixO^L
        double precision, intent(in)    :: x(ixI^S,1:ndim)
        double precision, intent(in) :: w(ixI^S,1:nw)
        double precision, intent(out)   :: res(ixI^S)

        double precision :: Rfactor(ixI^S)

        timeeos0 = MPI_WTIME()

        call eos%get_Rfactor(w,x,ixI^L,ixI^L,Rfactor)
        res(ixO^S) = (eos%gamma_minus_1 * w(ixO^S,iw_e) / (Rfactor(ixO^S) * w(ixO^S,iw_rho))) !> pth/rho

        timeeos_Tfromei=timeeos_Tfromei+(MPI_WTIME()-timeeos0)
    end subroutine get_temperature_from_eint_FI

    !> FI temperature from primitive pressure: T = p / (R * rho).
    !> w is PRIMITIVE here, so iw_e holds the thermal pressure. Lives in the EoS
    !> so hd and mhd share one routine (replaces the orphaned, never-called
    !> hd_get_temperature_from_prim that used to sit in mod_hd_phys).
    subroutine get_temperature_from_pressure_FI(w, x, ixI^L, ixO^L, res)
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(in)    :: x(ixI^S,1:ndim)
        double precision, intent(in)    :: w(ixI^S,1:nw)
        double precision, intent(out)   :: res(ixI^S)

        double precision :: Rfactor(ixI^S)

        call eos%get_Rfactor(w, x, ixI^L, ixO^L, Rfactor)
        res(ixO^S) = w(ixO^S,iw_e) / (Rfactor(ixO^S) * w(ixO^S,iw_rho))
    end subroutine get_temperature_from_pressure_FI

    !> Floor on log10(T) for the (rho,T) inverse tables. Encapsulates the
    !> method->container choice so callers need not know EoS internals: the
    !> 'tables' method populates eint_from_T, the entropy method eintT; picking
    !> the wrong one leaves var2_min=0 and clobbers cold cells. FI has no such
    !> table, so the floor is log10(smalldouble). Snapshot this into the sub-
    !> physics fluid objects at bind time (it is fixed once tables are loaded).
    double precision function eos_get_log_T_floor() result(log_T_min)
        if (eos%method_id == EOS_ENTROPY .and. allocated(eos%eintT%table)) then
            log_T_min = eos%eintT%var2_min
        else if (allocated(eos%eint_from_T%table)) then
            log_T_min = eos%eint_from_T%var2_min
        else
            log_T_min = dlog10(smalldouble)
        end if
    end function eos_get_log_T_floor

    !> log_nH grid metadata of the (log_nH, log_T) inverse table (eint from T),
    !> choosing the container by method ('tables'->eint_from_T, 'entropy'->eintT).
    !> n_nH=0 if no such table (analytic/FI). Lets cooling's build_Y_mod_table
    !> get its grid without reaching into EoS table internals.
    subroutine eos_get_eintT_grid(n_nH, lg_nH_min, lg_nH_max)
        integer, intent(out)          :: n_nH
        double precision, intent(out) :: lg_nH_min, lg_nH_max

        if (allocated(eos%eint_from_T%table)) then
            n_nH = eos%eint_from_T%dim1
            lg_nH_min = eos%eint_from_T%var1_min
            lg_nH_max = eos%eint_from_T%var1_max
        else if (allocated(eos%eintT%table)) then
            n_nH = eos%eintT%dim1
            lg_nH_min = eos%eintT%var1_min
            lg_nH_max = eos%eintT%var1_max
        else
            n_nH = 0
            lg_nH_min = 0.0d0
            lg_nH_max = 0.0d0
        end if
    end subroutine eos_get_eintT_grid

    subroutine get_temperature_from_eint_LTE(w, x, ixI^L, ixO^L, res)
        !> Assumes input energy is internal energy.
        !> Includes FI bypass: cells with eint/rho above threshold skip table lookups.
        use mod_physics
        integer, intent(in)             :: ixI^L,ixO^L
        double precision, intent(in)    :: x(ixI^S,1:ndim)
        double precision, intent(in) :: w(ixI^S,1:nw)
        double precision, intent(out)   :: res(ixI^S)

        double precision :: nH(ixI^S),nH_in(ixI^S), eint_in(ixI^S)
        double precision :: Rfactor(ixI^S), Rfactor_FI
        integer :: ix^D

        timeeos0 = MPI_WTIME()

        Rfactor_FI = eos%n_per_nH_FI / (1.0d0 + 4.0d0*eos%He_abundance)

        call eos%get_nH(w, x, ixI^L, ixO^L, nH)
        nH_in(ixO^S) = dlog10(nH(ixO^S))
        eint_in(ixO^S) = dlog10(w(ixO^S,iw_e)) - nH_in(ixO^S)

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            if (w(ix^D,iw_e) / w(ix^D,iw_rho) > eos%eint_rho_FI_threshold) then
                !> Fully ionised: T = (eint - eion*nH) * (gamma-1) / (Rfactor_FI * rho)
                res(ix^D) = eos%gamma_minus_1 &
                    * (w(ix^D,iw_e) - eos%eion_per_nH * nH(ix^D)) &
                    / (Rfactor_FI * w(ix^D,iw_rho))
            else if (eos%method_id == EOS_ANALYTIC) then
                block0: block
                    double precision :: T_loc, y_loc
                    call saha_T_from_nH_eint(nH(ix^D), &
                        w(ix^D,iw_e) / nH(ix^D), T_loc, y_loc)
                    res(ix^D) = T_loc
                end block block0
            else
                res(ix^D) = T_from_nH_eint(nH_in(ix^D),eint_in(ix^D))
            endif
        {end do\}

        timeeos_Tfromei=timeeos_Tfromei+(MPI_WTIME()-timeeos0)
    end subroutine get_temperature_from_eint_LTE

    subroutine get_temperature_from_eint_LTE_fast(w, x, ixI^L, ixO^L, res)
        !> Fast TC variant: two-pass regime-aware bypass.
        !>
        !> Pass 1 (vectorised): compute FI formula for ALL cells as array ops.
        !>   For the ~97% fully-ionised cells (T > 50 kK), this is exact.
        !>   For the ~3% ionisation-zone cells, the result is immediately overwritten.
        !>   The compiler vectorises this with SVML (no branches, pure arithmetic).
        !>
        !> Pass 2 (scalar, ~3% of cells): overwrite ionisation-zone cells with
        !>   bilinear table lookup + dexp.  Threshold check uses multiply
        !>   (eint <= threshold * rho) to avoid a 14-cycle scalar division.
        use mod_physics
        integer, intent(in)             :: ixI^L,ixO^L
        double precision, intent(in)    :: x(ixI^S,1:ndim)
        double precision, intent(in)    :: w(ixI^S,1:nw)
        double precision, intent(out)   :: res(ixI^S)

        double precision :: inv_Rfactor_FI, eion_rho_inv
        double precision :: log_nH_val, log_eint_nH_val
        double precision :: fx, fy, rx, ry
        integer :: jx, jy, jx1, jy1, ix^D
        double precision, parameter :: ln10 = 2.302585092994046d0

        timeeos0 = MPI_WTIME()

        !> Precompute scalar constants (avoid per-cell divisions)
        inv_Rfactor_FI = (1.0d0 + 4.0d0*eos%He_abundance) / eos%n_per_nH_FI
        eion_rho_inv = eos%eion_per_nH / eos%nH2rhoFactor

        !> Pass 1: FI formula for ALL cells (vectorisable array operations).
        !> T = (gamma-1) * (eint - eion_per_rho * rho) / (Rfactor_FI * rho)
        res(ixO^S) = eos%gamma_minus_1 * inv_Rfactor_FI &
            * (w(ixO^S,iw_e) - eion_rho_inv * w(ixO^S,iw_rho)) &
            / w(ixO^S,iw_rho)

        !> Pass 2: overwrite ionisation-zone cells.
        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            if (w(ix^D,iw_e) <= eos%eint_rho_FI_threshold * w(ix^D,iw_rho)) then
                if (eos%method_id == EOS_ANALYTIC) then
                    block1: block
                        double precision :: nH_loc, T_loc, y_loc
                        nH_loc = w(ix^D, iw_rho) / eos%nH2rhoFactor
                        call saha_T_from_nH_eint(nH_loc, &
                            w(ix^D,iw_e) / nH_loc, T_loc, y_loc)
                        res(ix^D) = T_loc
                    end block block1
                else if (eos%method_id == EOS_ENTROPY) then
                    !> Entropy method dispatches T_from_nH_eint to the
                    !> bicubic-Hermite forward (Tfwd) -- eos%T is NOT loaded
                    !> for entropy, so the legacy bilinear path below would
                    !> read an unallocated table -> SIGSEGV. Use the entropy
                    !> bicubic kernel instead. Reuses block%wextra log_nH
                    !> cache when present (same fast path as legacy).
                    block2: block
                        use mod_eos_entropy, only: entropy_T_from_log_nH_eint
                        if (iw_log_nH > 0) then
                            log_nH_val = block%wextra(ix^D, iw_log_nH)
                        else
                            log_nH_val = dlog10(w(ix^D, iw_rho) / eos%nH2rhoFactor)
                        end if
                        log_eint_nH_val = dlog10(w(ix^D, iw_e)) - log_nH_val
                        res(ix^D) = entropy_T_from_log_nH_eint(eos%Tfwd, eos%Tfwd_x, &
                            eos%Tfwd_y, eos%Tfwd_xy, log_nH_val, log_eint_nH_val)
                    end block block2
                else
                    if (iw_log_nH > 0) then
                        log_nH_val = block%wextra(ix^D, iw_log_nH)
                    else
                        log_nH_val = dlog10(w(ix^D, iw_rho) / eos%nH2rhoFactor)
                    end if
                    log_eint_nH_val = dlog10(w(ix^D, iw_e)) - log_nH_val

                    if (eos%T%is_uniform) then
                        !> Inlined uniform bilinear (hot path: no function call).
                        ry = max(0.0d0, min((log_nH_val - eos%T%var1_min) * eos%T%step_inv_1, &
                                             dble(eos%T%dim1-1)))
                        rx = max(0.0d0, min((log_eint_nH_val - eos%T%var2_min) * eos%T%step_inv_2, &
                                             dble(eos%T%dim2-1)))
                        jy = int(ry); jx = int(rx)
                        jy1 = min(jy+1, eos%T%dim1-1)
                        jx1 = min(jx+1, eos%T%dim2-1)
                        fy = ry - dble(jy); fx = rx - dble(jx)
                        res(ix^D) = dexp(ln10 * ( &
                            (1.0d0-fy)*((1.0d0-fx)*eos%T%table(jy+1,jx+1)  &
                                              + fx *eos%T%table(jy+1,jx1+1)) &
                           +      fy *((1.0d0-fx)*eos%T%table(jy1+1,jx+1)   &
                                              + fx *eos%T%table(jy1+1,jx1+1))))
                    else
                        !> Non-uniform: dispatch to the binary-search bilinear.
                        !> Slightly slower per cell (binary search adds ~8 cmps per axis)
                        !> but unavoidable when the grid is non-uniform.
                        res(ix^D) = dexp(ln10 * bilinear_lookup( &
                            log_nH_val, log_eint_nH_val, eos%T))
                    end if
                end if
            end if
        {end do\}

        timeeos_Tfromei = timeeos_Tfromei + (MPI_WTIME()-timeeos0)
    end subroutine get_temperature_from_eint_LTE_fast

    subroutine get_temperature_from_etot(w, x, ixI^L, ixO^L,  res)
        !> Assumes input energy is total energy
        use mod_physics
        integer, intent(in)             :: ixI^L,ixO^L
        double precision, intent(in)    :: x(ixI^S,1:ndim)
        double precision, intent(in) :: w(ixI^S,1:nw)
        double precision, intent(out)   :: res(ixI^S)
        double precision :: wlocal(ixI^S,1:nw)

        !> No timing here: get_temperature_from_eint is already timed.
        !> The array copy and e_to_ei subtraction are trivial.
        wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
        call phys_e_to_ei(ixI^L, ixO^L, wlocal, x)
        call eos%get_temperature_from_eint(wlocal, x, ixI^L, ixO^L, res)

    end subroutine get_temperature_from_etot

    !> Ionization fraction from (log10 nH, log10 eint/nH) in code units.
    !> Dispatches: analytic -> Saha quadratic, tables -> PCHIP interpolation.
    double precision function y_from_nH_eint(nH, eint_nh) result(result_val)
        use mod_eos_entropy, only: entropy_y_from_log_nH_eint
        double precision, intent(in) :: nH, eint_nh
        double precision, parameter :: ln10 = 2.302585092994046d0
        double precision :: T_loc, y_loc, eint_rho

        if (eos%method_id == EOS_ANALYTIC) then
            ! FI bypass: skip Saha solve for fully ionized cells
            eint_rho = 10.0d0**eint_nh / eos%nH2rhoFactor
            if (eint_rho > eos%eint_rho_FI_threshold) then
                result_val = eos%neOnH_FI
                return
            end if
            call saha_T_from_nH_eint(10.0d0**nH, 10.0d0**eint_nh, T_loc, y_loc)
            result_val = y_loc
        else if (eos%method_id == EOS_ENTROPY) then
            ! neOnH stored linearly (not log10); single bicubic Hermite lookup.
            result_val = entropy_y_from_log_nH_eint(eos%neOnH, eos%neOnH_x, &
                eos%neOnH_y, eos%neOnH_xy, nH, eint_nh)
        else
            result_val = dexp(ln10 * bicubic_lookup(nH, eint_nh, eos%neOnH))
        end if
    end function y_from_nH_eint

    !> Temperature from (log10 nH, log10 eint/nH) in code units.
    !> Dispatches: analytic -> Saha bisection/Newton, tables -> PCHIP interpolation.
    double precision function T_from_nH_eint(nH, eint_nh) result(result_val)
        use mod_eos_entropy, only: entropy_T_from_log_nH_eint
        double precision, intent(in) :: nH, eint_nh
        double precision, parameter :: ln10 = 2.302585092994046d0
        double precision :: T_loc, y_loc, eint_rho, Rfactor_FI

        if (eos%method_id == EOS_ANALYTIC) then
            ! FI bypass: skip Saha solve for fully ionized cells
            eint_rho = 10.0d0**eint_nh / eos%nH2rhoFactor
            if (eint_rho > eos%eint_rho_FI_threshold) then
                Rfactor_FI = eos%n_per_nH_FI / (1.0d0 + 4.0d0*eos%He_abundance)
                result_val = eos%gamma_minus_1 &
                    * (10.0d0**eint_nh - eos%eion_per_nH) / Rfactor_FI
                return
            end if
            call saha_T_from_nH_eint(10.0d0**nH, 10.0d0**eint_nh, T_loc, y_loc)
            result_val = T_loc
        else if (eos%method_id == EOS_ENTROPY) then
            result_val = entropy_T_from_log_nH_eint(eos%Tfwd, eos%Tfwd_x, &
                eos%Tfwd_y, eos%Tfwd_xy, nH, eint_nh)
        else
            result_val = dexp(ln10 * bicubic_lookup(nH, eint_nh, eos%T))
        end if
    end function T_from_nH_eint

    !> Fused T+y lookup from (log10 nH, log10 eint/nH) in code units.
    !> Computes grid indices once, evaluates both T and y tables.
    !> Saves one index computation + better cache utilisation vs separate calls.
    subroutine T_and_y_from_nH_eint(log_nH, log_eint_nH, T_out, y_out)
        use mod_eos_entropy, only: entropy_T_and_y_from_log_nH_eint
        double precision, intent(in) :: log_nH, log_eint_nH
        double precision, intent(out) :: T_out, y_out
        double precision, parameter :: ln10 = 2.302585092994046d0
        double precision :: results(3)

        if (eos%method_id == EOS_ANALYTIC) then
            call saha_T_from_nH_eint(10.0d0**log_nH, 10.0d0**log_eint_nH, T_out, y_out)
        else if (eos%method_id == EOS_ENTROPY) then
            call entropy_T_and_y_from_log_nH_eint(eos%Tfwd, eos%Tfwd_x, &
                eos%Tfwd_y, eos%Tfwd_xy, &
                eos%neOnH, eos%neOnH_x, eos%neOnH_y, eos%neOnH_xy, &
                log_nH, log_eint_nH, T_out, y_out)
        else if (allocated(eos%table_eint_il)) then
            ! Interleaved PCHIP fast-path. For uniform tables, the affine
            ! index calculation is inlined; for adaptive tables, the same
            ! kernel runs with binary-search index calculation. Both keep
            ! the stride-1 access pattern over the 3 packed slots
            ! (T, n_e/n_H, p/n_H), so cache behaviour is preserved.
            if (eos%T%is_uniform) then
                call interp_pchip_interleaved(log_nH, log_eint_nH, &
                    eos%table_eint_il, 3, eos%T%dim2, eos%T%dim1, &
                    eos%T%var1_min, eos%T%var1_max, &
                    eos%T%var2_min, eos%T%var2_max, results)
            else
                call interp_pchip_interleaved_nu(log_nH, log_eint_nH, &
                    eos%table_eint_il, 3, eos%T%dim2, eos%T%dim1, &
                    eos%T%var1_nodes, eos%T%var2_nodes, &
                    eos%T%guard_1, eos%T%guard_M_1, eos%T%guard_scale_1, &
                    eos%T%guard_2, eos%T%guard_M_2, eos%T%guard_scale_2, &
                    results)
            end if
            T_out = dexp(ln10 * results(1))
            y_out = dexp(ln10 * results(2))
        else
            ! No interleaved table built -- separate dispatcher calls
            T_out = dexp(ln10 * bicubic_lookup(log_nH, log_eint_nH, eos%T))
            y_out = dexp(ln10 * bicubic_lookup(log_nH, log_eint_nH, eos%neOnH))
        end if
    end subroutine T_and_y_from_nH_eint

    !> Pressure-to-eint ratio from (log10 nH, log10 p/nH) in code units.
    !> Dispatches: analytic -> Saha solve for eint/p, tables -> PCHIP interpolation.
    double precision function p2eint_from_nH_p(nH, ponH) result(result_val)
        use mod_eos_entropy, only: entropy_eint_from_nH_p
        double precision, intent(in) :: nH, ponH
        double precision :: nH_code, p_code, T_loc, y_loc, eint_nH_loc, p_rho

        if (eos%method_id == EOS_ANALYTIC) then
            nH_code = 10.0d0**nH
            p_code = nH_code * 10.0d0**ponH
            ! FI bypass: skip Saha solve for fully ionized cells
            p_rho = p_code / (nH_code * eos%nH2rhoFactor)
            if (p_rho > eos%p_rho_FI_threshold) then
                result_val = eos%inv_gamma_minus_1 &
                    + eos%eion_per_nH * nH_code / p_code
                return
            end if
            call saha_p_to_T(nH_code, p_code, T_loc, y_loc, eint_nH_loc)
            result_val = eint_nH_loc * nH_code / p_code
        else if (eos%method_id == EOS_ENTROPY) then
            ! One bicubic-Hermite lookup on the eintP table (built at
            ! generation time by inverting the pfwd polynomial, exact at nodes).
            result_val = entropy_eint_from_nH_p(eos%eintP, eos%eintP_x, &
                eos%eintP_y, eos%eintP_xy, nH, ponH)
        else
            result_val = bicubic_lookup(nH, ponH, eos%p2eint)
        end if
    end function p2eint_from_nH_p

    !> Gamma_1 from pressure-indexed table: (log10 nH, log10 p/nH) -> Gamma_1.
    !> For 'entropy' the conversion p -> eint -> Gamma_1 via formula keeps Maxwell
    !> consistency with the forward at runtime. The p->eint inverse is one
    !> bisection per cell -- non-trivial cost; this function is in the hot path
    !> via hd_get_csound2_LTE.
    double precision function gamma1_from_nH_p(log_nH, log_p_nH) result(g1)
        use mod_eos_entropy, only: entropy_gamma1_from_nH_p
        double precision, intent(in) :: log_nH, log_p_nH
        if (eos%method_id == EOS_ENTROPY) then
            ! Single bicubic Hermite lookup on the pre-built Gamma_1(rho, p) table.
            ! ZERO runtime iterations; no p->eint intermediate inversion.
            g1 = entropy_gamma1_from_nH_p(eos%g1p, eos%g1p_x, eos%g1p_y, &
                eos%g1p_xy, log_nH, log_p_nH)
        else
            g1 = bicubic_lookup(log_nH, log_p_nH, eos%gamma1_p)
        end if
    end function gamma1_from_nH_p

    !> Merged log10(p/nH) lookup: (log10 nH, log10 eint/nH) -> log10(p/nH)
    !> Single PCHIP evaluation replacing separate T + neOnH lookups.
    double precision function log_p_from_nH_eint(log_nH, log_eint_nH) result(lp)
        double precision, intent(in) :: log_nH, log_eint_nH
        lp = bicubic_lookup(log_nH, log_eint_nH, eos%log_p)
    end function log_p_from_nH_eint

    !> Cached bisection on the log_p table for WB p->eint inversion.
    !> Precomputes nH-direction indices and table values once, then
    !> bisects using only cheap PCHIP evaluations with varying tx.
    !> Expects a narrow initial bracket [lo, hi] (e.g. from p2eint guess).
    !>
    !> Adaptive-grid support: when eos%log_p%is_uniform is .false., the
    !> y- and x-direction index calculations switch to binary search on
    !> the explicit node arrays. The cached 4x4 stencil mechanism is
    !> unchanged; only the index-to-cell mapping changes.
    subroutine log_p_bisect_cached(log_nH, log_p_target, &
        log_eint_lo, log_eint_hi, max_iter, log_eint_result)
        use mod_lookup_table, only: find_index_bsearch
        double precision, intent(in)    :: log_nH, log_p_target
        double precision, intent(inout) :: log_eint_lo, log_eint_hi
        integer, intent(in)             :: max_iter
        double precision, intent(out)   :: log_eint_result

        integer :: nx, ny, ix, iy, iter, ii
        integer :: i0, i1, i2, i3, j0, j1, j2, j3
        double precision :: tx, ty, rx, ry
        double precision :: xstep_inv, ystep_inv
        double precision :: vmin_x, vmax_x, vmin_y, vmax_y
        double precision :: log_eint_mid, log_p_eval
        logical :: is_unif

        !> Table values: tv(y_row, x_col) for 4 y-rows x 4 x-cols
        double precision :: tv(4,4)
        !> Cached ix for detecting grid cell change
        integer :: ix_cached

        nx = eos%log_p%dim2   ! eint axis (varx)
        ny = eos%log_p%dim1   ! nH axis (vary)
        is_unif = eos%log_p%is_uniform

        if (is_unif) then
            vmin_x = eos%log_p%var2_min
            vmax_x = eos%log_p%var2_max
            vmin_y = eos%log_p%var1_min
            vmax_y = eos%log_p%var1_max
            xstep_inv = dble(nx-1) / (vmax_x - vmin_x)
            ystep_inv = dble(ny-1) / (vmax_y - vmin_y)
        end if

        !> Precompute nH indices (constant across all iterations)
        if (is_unif) then
            ry = (log_nH - vmin_y) * ystep_inv
            ry = max(0.0d0, min(ry, dble(ny-1)))
            iy = int(ry)
            ty = ry - dble(iy)
        else
            !> Non-uniform: binary search on var1_nodes (length ny = dim1)
            if (log_nH <= eos%log_p%var1_nodes(1)) then
                iy = 0; ty = 0.0d0
            else if (log_nH >= eos%log_p%var1_nodes(ny)) then
                iy = ny - 2; ty = 1.0d0
            else
                ii = find_index_guard(eos%log_p%var1_nodes, ny, log_nH, &
                    eos%log_p%guard_1, eos%log_p%guard_M_1, eos%log_p%guard_scale_1)
                iy = max(0, min(ii - 2, ny - 2))
                ty = (log_nH - eos%log_p%var1_nodes(iy+1)) &
                   / (eos%log_p%var1_nodes(iy+2) - eos%log_p%var1_nodes(iy+1))
                ty = max(0.0d0, min(ty, 1.0d0))
            end if
        end if

        !> Clamped y-row indices
        j0 = max(0, min(ny-1, iy-1))
        j1 = max(0, min(ny-1, iy  ))
        j2 = max(0, min(ny-1, iy+1))
        j3 = max(0, min(ny-1, iy+2))

        ix_cached = -1

        do iter = 1, max_iter
            log_eint_mid = 0.5d0 * (log_eint_lo + log_eint_hi)

            !> Compute eint-direction index
            if (is_unif) then
                rx = (log_eint_mid - vmin_x) * xstep_inv
                rx = max(0.0d0, min(rx, dble(nx-1)))
                ix = int(rx)
                tx = rx - dble(ix)
            else
                !> Non-uniform: binary search on var2_nodes (length nx = dim2)
                if (log_eint_mid <= eos%log_p%var2_nodes(1)) then
                    ix = 0; tx = 0.0d0
                else if (log_eint_mid >= eos%log_p%var2_nodes(nx)) then
                    ix = nx - 2; tx = 1.0d0
                else
                    ii = find_index_guard(eos%log_p%var2_nodes, nx, log_eint_mid, &
                        eos%log_p%guard_2, eos%log_p%guard_M_2, eos%log_p%guard_scale_2)
                    ix = max(0, min(ii - 2, nx - 2))
                    tx = (log_eint_mid - eos%log_p%var2_nodes(ix+1)) &
                       / (eos%log_p%var2_nodes(ix+2) - eos%log_p%var2_nodes(ix+1))
                    tx = max(0.0d0, min(tx, 1.0d0))
                end if
            end if

            !> Load 16 table values only when ix changes
            if (ix /= ix_cached) then
                i0 = max(0, min(nx-1, ix-1))
                i1 = max(0, min(nx-1, ix  ))
                i2 = max(0, min(nx-1, ix+1))
                i3 = max(0, min(nx-1, ix+2))

                tv(1,1) = eos%log_p%table(j0+1, i0+1)
                tv(1,2) = eos%log_p%table(j0+1, i1+1)
                tv(1,3) = eos%log_p%table(j0+1, i2+1)
                tv(1,4) = eos%log_p%table(j0+1, i3+1)
                tv(2,1) = eos%log_p%table(j1+1, i0+1)
                tv(2,2) = eos%log_p%table(j1+1, i1+1)
                tv(2,3) = eos%log_p%table(j1+1, i2+1)
                tv(2,4) = eos%log_p%table(j1+1, i3+1)
                tv(3,1) = eos%log_p%table(j2+1, i0+1)
                tv(3,2) = eos%log_p%table(j2+1, i1+1)
                tv(3,3) = eos%log_p%table(j2+1, i2+1)
                tv(3,4) = eos%log_p%table(j2+1, i3+1)
                tv(4,1) = eos%log_p%table(j3+1, i0+1)
                tv(4,2) = eos%log_p%table(j3+1, i1+1)
                tv(4,3) = eos%log_p%table(j3+1, i2+1)
                tv(4,4) = eos%log_p%table(j3+1, i3+1)

                ix_cached = ix
            end if

            !> Evaluate PCHIP: 4 x-rows then 1 y-interp
            log_p_eval = pchip_2d_from_cache(tv, tx, ty)

            if (log_p_eval < log_p_target) then
                log_eint_lo = log_eint_mid
            else
                log_eint_hi = log_eint_mid
            end if

            if (dabs(log_eint_hi - log_eint_lo) < 1.0d-14) exit
        end do

        log_eint_result = 0.5d0 * (log_eint_lo + log_eint_hi)

    contains

        pure double precision function pchip_1d(p0, p1, p2, p3, t) result(v)
            double precision, intent(in) :: p0, p1, p2, p3, t
            double precision :: d0, d1, d2, m1, m2, s, a1, a2, lim
            double precision :: tt, ttt, h00, h10, h01, h11

            d0 = p1 - p0;  d1 = p2 - p1;  d2 = p3 - p2

            if (d1 == 0.0d0) then
                m1 = 0.0d0;  m2 = 0.0d0
            else
                if (d0*d1 <= 0.0d0) then
                    m1 = 0.0d0
                else
                    m1 = 2.0d0*d0*d1 / (d0 + d1)
                end if
                if (d1*d2 <= 0.0d0) then
                    m2 = 0.0d0
                else
                    m2 = 2.0d0*d1*d2 / (d1 + d2)
                end if
                s = sign(1.0d0, d1)
                a1 = s*m1;  a2 = s*m2
                if (a1 < 0.0d0) a1 = 0.0d0
                if (a2 < 0.0d0) a2 = 0.0d0
                lim = 3.0d0*abs(d1)
                if (a1 > lim) a1 = lim
                if (a2 > lim) a2 = lim
                m1 = s*a1;  m2 = s*a2
            end if

            tt = t*t;  ttt = tt*t
            h00 = 2.0d0*ttt - 3.0d0*tt + 1.0d0
            h10 = ttt - 2.0d0*tt + t
            h01 = -2.0d0*ttt + 3.0d0*tt
            h11 = ttt - tt
            v = h00*p1 + h10*m1 + h01*p2 + h11*m2
        end function pchip_1d

        pure double precision function pchip_2d_from_cache(c, tx, ty) result(z)
            double precision, intent(in) :: c(4,4), tx, ty
            double precision :: g0, g1, g2, g3
            !> 4 x-direction interpolations (one per y-row)
            g0 = pchip_1d(c(1,1), c(1,2), c(1,3), c(1,4), tx)
            g1 = pchip_1d(c(2,1), c(2,2), c(2,3), c(2,4), tx)
            g2 = pchip_1d(c(3,1), c(3,2), c(3,3), c(3,4), tx)
            g3 = pchip_1d(c(4,1), c(4,2), c(4,3), c(4,4), tx)
            !> 1 y-direction interpolation
            z = pchip_1d(g0, g1, g2, g3, ty)
        end function pchip_2d_from_cache

    end subroutine log_p_bisect_cached

    !> Given log10(nH) and log10(p), find log10(eint/nH) by table-guessed
    !> bisection on the forward pressure table.  Wraps the bracketing logic
    !> and log_p_bisect_cached into a single call for use by both HD and MHD
    !> from_conserved routines.
    subroutine eint_from_p_bisect(log_nH_val, log_p_val, log_eint_nH_out)
        double precision, intent(in)  :: log_nH_val, log_p_val
        double precision, intent(out) :: log_eint_nH_out

        double precision :: log_p_target, p2eint_ratio, log_eint_guess
        double precision :: log_p_at_guess, log_eint_lo, log_eint_hi
        double precision :: f_bracket, margin
        integer :: max_iter

        log_p_target = log_p_val - log_nH_val

        ! Initial guess from p2eint table
        p2eint_ratio = p2eint_from_nH_p(log_nH_val, log_p_target)
        log_eint_guess = dlog10(p2eint_ratio) + log_p_target

        log_eint_guess = max(log_eint_guess, eos%log_p%var2_min)
        log_eint_guess = min(log_eint_guess, eos%log_p%var2_max)

        log_p_at_guess = log_p_from_nH_eint(log_nH_val, log_eint_guess)

        ! Establish bracket around the guess
        margin = 5.0d-4
        if (log_p_at_guess < log_p_target) then
            log_eint_lo = log_eint_guess
            log_eint_hi = min(log_eint_guess + margin, eos%log_p%var2_max)
            f_bracket = log_p_from_nH_eint(log_nH_val, log_eint_hi) - log_p_target
        else
            log_eint_lo = max(log_eint_guess - margin, eos%log_p%var2_min)
            log_eint_hi = log_eint_guess
            f_bracket = -(log_p_from_nH_eint(log_nH_val, log_eint_lo) - log_p_target)
        end if

        if (f_bracket >= 0.0d0) then
            max_iter = 8
        else
            log_eint_lo = eos%log_p%var2_min
            log_eint_hi = eos%log_p%var2_max
            max_iter = 20
        end if

        call log_p_bisect_cached(log_nH_val, log_p_target, &
            log_eint_lo, log_eint_hi, max_iter, log_eint_nH_out)

    end subroutine eint_from_p_bisect

    !> Internal energy per nH from (log10 nH, log10 T) in code units.
    !> Uses the bisection-built inverse table (H+He, machine precision).
    !> Fallback: H-only Saha if table not built.
    double precision function eint_nH_from_T(log_nH, log_T) result(eint_nH)
        use mod_eos_entropy, only: entropy_eint_from_nH_T
        double precision, intent(in) :: log_nH, log_T
        double precision, parameter :: ln10 = 2.302585092994046d0
        double precision :: log_e_nh

        if (eos%method_id == EOS_ENTROPY) then
            ! Single bicubic Hermite lookup on the pre-built eint(rho, T)
            ! inverse table. ZERO runtime iterations.
            log_e_nh = entropy_eint_from_nH_T(eos%eintT, eos%eintT_x, &
                eos%eintT_y, eos%eintT_xy, log_nH, log_T)
            eint_nH = 10.0d0**log_e_nh
        else if (allocated(eos%eint_from_T%table)) then
            eint_nH = dexp(ln10 * bicubic_lookup(log_nH, log_T, eos%eint_from_T))
        else
            eint_nH = saha_eint_from_nH_T(10.0d0**log_nH, 10.0d0**log_T)
        end if
    end function eint_nH_from_T

    !> p/nH from (log10 nH, log10 eint/nH) in code units.
    !> Returns (1+He+y)*T directly -- single lookup replaces T + y lookups.
    double precision function p_nH_from_eint(log_nH, log_eint_nH) result(p_nH)
        use mod_eos_entropy, only: entropy_p_nH_from_eint
        double precision, intent(in) :: log_nH, log_eint_nH
        double precision, parameter :: ln10 = 2.302585092994046d0

        if (eos%method_id == EOS_ANALYTIC) then
            block
                double precision :: T_loc, y_loc
                call saha_T_from_nH_eint(10.0d0**log_nH, 10.0d0**log_eint_nH, T_loc, y_loc)
                p_nH = (1.0d0 + eos%He_abundance + y_loc) * T_loc
            end block
        else if (eos%method_id == EOS_ENTROPY) then
            p_nH = entropy_p_nH_from_eint(eos%pfwd, eos%pfwd_x, eos%pfwd_y, &
                                           eos%pfwd_xy, log_nH, log_eint_nH)
        else
            p_nH = dexp(ln10 * bicubic_lookup(log_nH, log_eint_nH, eos%p_over_nH))
        end if
    end function p_nH_from_eint

    !> routines needed for Chun's partially ionised module (From Leenaarts et al. 2012?)
    subroutine Rfactor_from_PI_temperature(w,x,ixI^L,ixO^L,Rfactor)
        use mod_global_parameters
        use mod_ionization_degree
        integer, intent(in) :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S,1:nw)
        double precision, intent(in) :: x(ixI^S,1:ndim)
        double precision, intent(out):: Rfactor(ixI^S)

        double precision :: iz_H(ixO^S),iz_He(ixO^S)

        call ionization_degree_from_temperature(ixI^L,ixO^L,w(ixI^S,iw_te),iz_H,iz_He)
        ! assume the first and second ionization of Helium have the same degree
        Rfactor(ixO^S)=(1.d0+iz_H(ixO^S)+0.1d0*(1.d0+iz_He(ixO^S)*(1.d0+iz_He(ixO^S))))/(1.d0+4.d0*eos%He_abundance)

    end subroutine Rfactor_from_PI_temperature
    
    subroutine update_PI_temperature(ixI^L,ixO^L,wCT,w,x)
        use mod_global_parameters
        use mod_ionization_degree

        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
        double precision, intent(inout) :: w(ixI^S,1:nw)

        double precision :: iz_H(ixO^S),iz_He(ixO^S), pth(ixI^S)

        call ionization_degree_from_temperature(ixI^L,ixO^L,wCT(ixI^S,iw_te),iz_H,iz_He)

        call eos%get_thermal_pressure(w,x,ixI^L,ixO^L,pth)

        w(ixO^S,iw_te)=(2.d0+3.d0*eos%He_abundance)*pth(ixO^S)/(w(ixO^S,iw_rho)*(1.d0+iz_H(ixO^S)+&
        eos%He_abundance*(iz_He(ixO^S)*(iz_He(ixO^S)+1.d0)+1.d0)))

    end subroutine update_PI_temperature

end module mod_eos
!> Needs a line after to pass the preprocesor