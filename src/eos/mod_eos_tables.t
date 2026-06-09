!=============================================================================
!> EoS table lifecycle (build/IO; not on the simulation hot path).
!>
!> Loads the binary LTE tables from disk (read_eos_from_file / load_lte_tables),
!> prepares each container for lookup (axis nodes, code-unit shift, O(1) guard
!> arrays, validation), builds the derived tables that the 'tables' method needs
!> (gamma1, log_p, p_over_nH, interleaved, pressure-indexed gamma1, the bisected
!> p2eint and eint-from-T inverses), and precomputes the fully-ionised bypass
!> constants. Everything here runs once, during eos_init/eos_finalise.
!=============================================================================
module mod_eos_tables
    use mod_global_parameters
    use mod_eos_container
    use mod_eos_interp
    use mod_comm_lib, only: mpistop
    implicit none
    private

    !> Table loading
    public :: load_lte_tables, try_load_lte_tables
    !> Per-container preparation
    public :: ensure_axis_nodes, shift_axis_to_code, shift_axis_to_code_T
    public :: entropy_table_prepare, eos_build_guards, eos_validate_table
    public :: precompute_FI_bypass_constants
    !> Derived-table builders ('tables' method) + round-trip check
    public :: build_gamma1_table, build_log_p_table, build_p_over_nH_table
    public :: build_interleaved_eint_table, build_gamma1_p_table
    public :: build_p2eint_table, build_eint_from_T_table, verify_eos_round_trips

    character(len=std_len), parameter :: eos_table_prefix = "LTEeos_"

contains

    !> Ensure tc%var1_nodes / tc%var2_nodes are allocated and populated so
    !> that build_*_table routines and other code can reference axis positions
    !> uniformly (no is_uniform branching). For uniform tables we generate
    !> evenly-spaced node arrays from var{1,2}_min/max.
    !>
    !> Lookup performance is unaffected: bicubic_lookup / bilinear_lookup
    !> still branch on tc%is_uniform -- uniform tables go through the affine
    !> fast path; only adaptive tables consult var{1,2}_nodes.
    subroutine ensure_axis_nodes(tc)
        type(eos_table_container), intent(inout) :: tc
        integer :: k
        double precision :: dx
        if (.not. allocated(tc%table)) return
        if (.not. allocated(tc%var1_nodes)) then
            allocate(tc%var1_nodes(tc%dim1))
            dx = (tc%var1_max - tc%var1_min) / dble(tc%dim1 - 1)
            do k = 1, tc%dim1
                tc%var1_nodes(k) = tc%var1_min + (k - 1) * dx
            end do
        end if
        if (.not. allocated(tc%var2_nodes)) then
            allocate(tc%var2_nodes(tc%dim2))
            dx = (tc%var2_max - tc%var2_min) / dble(tc%dim2 - 1)
            do k = 1, tc%dim2
                tc%var2_nodes(k) = tc%var2_min + (k - 1) * dx
            end do
        end if
    end subroutine ensure_axis_nodes

    !> Shift a table's (log_nH, log_eint/nH) axes -- and adaptive node arrays
    !> if present -- from CGS storage to code units.  Used by the entropy
    !> method so eos_finalise can apply the shift AFTER unit_* are populated.
    subroutine shift_axis_to_code(tc)
        !> CGS -> code shift assuming axis-1 = log(nH) and axis-2 = log(eint/nH)
        !> (or log(p/nH); both share the same erg/particle unit conversion).
        type(eos_table_container), intent(inout) :: tc
        tc%var1_min = tc%var1_min - dlog10(unit_numberdensity)
        tc%var1_max = tc%var1_max - dlog10(unit_numberdensity)
        tc%var2_min = tc%var2_min - dlog10(unit_pressure/unit_numberdensity)
        tc%var2_max = tc%var2_max - dlog10(unit_pressure/unit_numberdensity)
        if (.not. tc%is_uniform) then
            tc%var1_nodes = tc%var1_nodes - dlog10(unit_numberdensity)
            tc%var2_nodes = tc%var2_nodes - dlog10(unit_pressure/unit_numberdensity)
        end if
    end subroutine shift_axis_to_code

    subroutine entropy_table_prepare(tc)
        !> Bundle ensure_axis_nodes + precompute_step_inv + build_guards + validate,
        !> the identical four-step setup each loaded entropy table container needs.
        type(eos_table_container), intent(inout) :: tc
        call ensure_axis_nodes(tc)
        call precompute_step_inv(tc)
        call eos_build_guards(tc)
        call eos_validate_table(tc, trim(tc%filename))
    end subroutine entropy_table_prepare

    subroutine shift_axis_to_code_T(tc)
        !> CGS -> code shift for tables whose axis-2 = log(T) (K).
        !> Axis-1 = log(nH) shift as in the standard case.
        type(eos_table_container), intent(inout) :: tc
        tc%var1_min = tc%var1_min - dlog10(unit_numberdensity)
        tc%var1_max = tc%var1_max - dlog10(unit_numberdensity)
        tc%var2_min = tc%var2_min - dlog10(unit_temperature)
        tc%var2_max = tc%var2_max - dlog10(unit_temperature)
        if (.not. tc%is_uniform) then
            tc%var1_nodes = tc%var1_nodes - dlog10(unit_numberdensity)
            tc%var2_nodes = tc%var2_nodes - dlog10(unit_temperature)
        end if
    end subroutine shift_axis_to_code_T

    !> Build guard (bucket) arrays so that adaptive index lookup is O(1).
    !> No-op for uniform tables. Safe to call repeatedly.
    subroutine eos_build_guards(tc)
        type(eos_table_container), intent(inout) :: tc
        if (tc%is_uniform) return
        if (.not. allocated(tc%var1_nodes)) return
        if (.not. allocated(tc%var2_nodes)) return
        call build_guard_for_axis(tc%var1_nodes, tc%dim1, &
            tc%guard_1, tc%guard_M_1, tc%guard_scale_1)
        call build_guard_for_axis(tc%var2_nodes, tc%dim2, &
            tc%guard_2, tc%guard_M_2, tc%guard_scale_2)
    end subroutine eos_build_guards

    subroutine build_guard_for_axis(nodes, n, guard, M_out, scale_out)
        double precision, intent(in) :: nodes(:)
        integer, intent(in) :: n
        integer, allocatable, intent(inout) :: guard(:)
        integer, intent(out) :: M_out
        double precision, intent(out) :: scale_out

        integer, parameter :: SAFETY = 4
        integer :: k, ii, M
        double precision :: xmin, xmax, span, min_gap, bk

        xmin = nodes(1); xmax = nodes(n)
        span = xmax - xmin
        if (span <= 0.0d0 .or. n < 2) then
            M_out = 0; scale_out = 0.0d0
            if (allocated(guard)) deallocate(guard)
            return
        end if
        min_gap = huge(1.0d0)
        do ii = 1, n-1
            if (nodes(ii+1) - nodes(ii) < min_gap) min_gap = nodes(ii+1) - nodes(ii)
        end do
        if (min_gap <= 0.0d0) then
            !> Degenerate node array -- fall back to bsearch by leaving M=0
            M_out = 0; scale_out = 0.0d0
            if (allocated(guard)) deallocate(guard)
            return
        end if
        M = max(2*n, SAFETY * int(ceiling(span / min_gap)))
        if (allocated(guard)) deallocate(guard)
        allocate(guard(M))
        scale_out = dble(M) / span
        M_out = M

        ii = 1
        do k = 1, M
            bk = xmin + dble(k-1) * span / dble(M)
            do while (ii < n .and. nodes(ii) < bk)
                ii = ii + 1
            end do
            guard(k) = ii
        end do
    end subroutine build_guard_for_axis

    !> Defensive consistency check for one table after init. Aborts with a
    !> diagnostic message if any silent-failure invariant is violated.
    subroutine eos_validate_table(tc, name)
        type(eos_table_container), intent(in) :: tc
        character(len=*), intent(in) :: name
        integer :: i

        if (.not. allocated(tc%table)) then
            if (mype == 0) write(*,'(A,A)') " EOS validate FAIL: table not allocated for ", trim(name)
            call mpistop("eos_validate_table: missing table data")
        end if
        if (size(tc%table, 1) /= tc%dim1 .or. size(tc%table, 2) /= tc%dim2) then
            if (mype == 0) write(*,'(A,A)') " EOS validate FAIL: dim mismatch for ", trim(name)
            call mpistop("eos_validate_table: dim mismatch")
        end if
        if (.not. tc%is_uniform) then
            if (.not. allocated(tc%var1_nodes) .or. .not. allocated(tc%var2_nodes)) then
                if (mype == 0) write(*,'(A,A)') &
                    " EOS validate FAIL: adaptive table missing node arrays for ", trim(name)
                call mpistop("eos_validate_table: adaptive table without nodes")
            end if
            if (size(tc%var1_nodes) /= tc%dim1 .or. size(tc%var2_nodes) /= tc%dim2) then
                if (mype == 0) write(*,'(A,A)') &
                    " EOS validate FAIL: node array length mismatch for ", trim(name)
                call mpistop("eos_validate_table: node length mismatch")
            end if
            do i = 1, tc%dim1 - 1
                if (tc%var1_nodes(i+1) <= tc%var1_nodes(i)) then
                    if (mype == 0) write(*,'(A,A,A,I0)') &
                        " EOS validate FAIL: var1_nodes not strictly monotonic for ", &
                        trim(name), " at index ", i
                    call mpistop("eos_validate_table: non-monotonic nodes")
                end if
            end do
            do i = 1, tc%dim2 - 1
                if (tc%var2_nodes(i+1) <= tc%var2_nodes(i)) then
                    if (mype == 0) write(*,'(A,A,A,I0)') &
                        " EOS validate FAIL: var2_nodes not strictly monotonic for ", &
                        trim(name), " at index ", i
                    call mpistop("eos_validate_table: non-monotonic nodes")
                end if
            end do
            if (.not. allocated(tc%guard_1) .or. .not. allocated(tc%guard_2) .or. &
                tc%guard_M_1 <= 0 .or. tc%guard_M_2 <= 0 .or. &
                tc%guard_scale_1 <= 0.0d0 .or. tc%guard_scale_2 <= 0.0d0) then
                if (mype == 0) write(*,'(A,A,A)') &
                    " EOS validate WARN: adaptive table without guard for ", trim(name), &
                    " - lookup will fall back to bsearch"
            end if
        end if
        if (mype == 0) write(*,'(A,A)') " EOS validate OK: ", trim(name)
    end subroutine eos_validate_table

    subroutine precompute_FI_bypass_constants()
        !> Precompute constants for bypassing table lookups when gas is fully ionised.
        !> Physics: above T ~ 50 kK, H is >99.99% ionised and He is >90% ionised.
        !> The EoS reduces to ideal gas with constant mu, Gamma_1 = 5/3, ne/nH = 1+2*A_He.
        !> We check eint/rho > threshold (1 division + 1 comparison) to skip all table lookups.
        double precision :: chi_H, chi_HeI, chi_HeII
        double precision :: T_thresh, eint_nH_thresh, p_nH_thresh

        !> Ionisation energies in CGS (erg)
        chi_H   = 13.6d0  * 1.602d-12
        chi_HeI = 24.6d0  * 1.602d-12
        chi_HeII = 54.4d0 * 1.602d-12

        !> Fully-ionised particle counts and electron fraction
        eos%n_per_nH_FI = 2.0d0 + 3.0d0 * eos%He_abundance
        eos%neOnH_FI    = 1.0d0 + 2.0d0 * eos%He_abundance

        !> Total ionisation energy per hydrogen atom (CGS), converted to code units
        !> eion_per_nH [code] = eion_per_nH [erg] / (unit_pressure / unit_numberdensity)
        eos%eion_per_nH = (chi_H + eos%He_abundance * (chi_HeI + chi_HeII)) &
                        * unit_numberdensity / unit_pressure

        !> Bypass threshold: eint/rho value corresponding to T = 200,000 K fully ionised.
        !> Raised from 50 kK to 200 kK so that blended LTE tables (which converge to
        !> FI by 200 kK) match the bypass formula exactly at the threshold -- no
        !> discontinuity in ne/nH, T, Gamma_1, or p when cells cross the boundary.
        !> eint/nH = 1.5 * n_per_nH * kB * T + eion_per_nH  (all in code units)
        !> eint/rho = eint/nH / nH2rhoFactor
        T_thresh = 200000.0d0 / unit_temperature
        eint_nH_thresh = 1.5d0 * eos%n_per_nH_FI * T_thresh + eos%eion_per_nH
        eos%eint_rho_FI_threshold = eint_nH_thresh / eos%nH2rhoFactor

        !> Pressure-based threshold for primitive-variable checks:
        !> For FI gas, T = p / (rho * Rfactor_FI), so p/rho > Rfactor_FI * T_thresh
        !> where Rfactor_FI = n_per_nH_FI / (1 + 4*A_He) in code units
        eos%p_rho_FI_threshold = eos%n_per_nH_FI &
            / (1.0d0 + 4.0d0*eos%He_abundance) * T_thresh

        !> Allow disabling the bypass for verification (forces all cells through table path)
        if (eos%disable_FI_bypass) then
            eos%eint_rho_FI_threshold = HUGE(1.0d0)
            eos%p_rho_FI_threshold    = HUGE(1.0d0)
        end if

        if (mype == 0) then
            if (eos%disable_FI_bypass) then
                write(*,'(a)') ' FI bypass: DISABLED (all cells use table path)'
            else
                write(*,'(a,es12.4)') ' FI bypass: eion_per_nH (code) = ', eos%eion_per_nH
                write(*,'(a,es12.4)') ' FI bypass: eint/rho threshold = ', eos%eint_rho_FI_threshold
                write(*,'(a,f8.4)')   ' FI bypass: n_per_nH_FI        = ', eos%n_per_nH_FI
                write(*,'(a,f8.4)')   ' FI bypass: neOnH_FI           = ', eos%neOnH_FI
            end if
        end if

    end subroutine precompute_FI_bypass_constants

    !> Read one named table file into its eos% container. The filename encodes
    !> the composition (H or HHe) and whether ionisation energy is included.
    subroutine load_lte_tables(fieldname)
        character(len=*), intent(in) :: fieldname
        character(len=std_len) :: subname, tablename, filename, filepath
        subname = "H" !> Always consider at least hydrogen

        if (eos%He_abundance > 0.0d0) then
            write(subname,"(A,A)") trim(subname), "He"
        endif

        if (eos%ionE) then !> Consider the ionisation energy?
            write(subname,"(A,A)") trim(subname), "_IonE"
        else
            write(subname,"(A,A)") trim(subname), "_NoIonE"
        endif

        write(tablename, "(A,A,A)") trim(eos_table_prefix), fieldname, "_"
        write(filename,"(A,A,A)") trim(tablename), trim(subname), '.bin'
        
        if (mype==0) then
            print*, "Reading EoS tables from: ", trim(filename)
        endif

        write(filepath,"(A,A)") trim(eos%table_location), trim(filename)

        select case (fieldname)
        case("T")
            eos%T%filename = filename
            call read_eos_from_file(trim(filepath), eos%T, .true.)
        case("neOnH")
            eos%neOnH%filename = filename
            ! Legacy tables (method='tables'/'analytic') store log10(neOnH); the
            ! entropy method's tables store raw (linear) neOnH because they go
            ! through bicubic Hermite (which interpolates the value directly,
            ! no log10 transform needed and would break stored derivatives).
            call read_eos_from_file(trim(filepath), eos%neOnH, &
                                     eos%method_id /= EOS_ENTROPY)
        case("p2eint")
            eos%p2eint%filename = filename
            call read_eos_from_file(trim(filepath), eos%p2eint, .false.)
        case("gamma1")
            eos%gamma1%filename = filename
            call read_eos_from_file(trim(filepath), eos%gamma1, .false.)
        case("eint_from_T")
            eos%eint_from_T%filename = filename
            call read_eos_from_file(trim(filepath), eos%eint_from_T, .false.)
        !> Entropy-method tables. Unit shifts MUST happen in eos_finalise (where
        !> unit_* are set), not here -- units aren't initialised at load time.
        case("neOnH_x")
            eos%neOnH_x%filename = filename
            call read_eos_from_file(trim(filepath), eos%neOnH_x, .false.)
        case("neOnH_y")
            eos%neOnH_y%filename = filename
            call read_eos_from_file(trim(filepath), eos%neOnH_y, .false.)
        case("neOnH_xy")
            eos%neOnH_xy%filename = filename
            call read_eos_from_file(trim(filepath), eos%neOnH_xy, .false.)
        case("eintP")
            eos%eintP%filename = filename
            call read_eos_from_file(trim(filepath), eos%eintP, .false.)
        case("eintP_x")
            eos%eintP_x%filename = filename
            call read_eos_from_file(trim(filepath), eos%eintP_x, .false.)
        case("eintP_y")
            eos%eintP_y%filename = filename
            call read_eos_from_file(trim(filepath), eos%eintP_y, .false.)
        case("eintP_xy")
            eos%eintP_xy%filename = filename
            call read_eos_from_file(trim(filepath), eos%eintP_xy, .false.)
        case("g1p")
            eos%g1p%filename = filename
            call read_eos_from_file(trim(filepath), eos%g1p, .false.)
        case("g1p_x")
            eos%g1p_x%filename = filename
            call read_eos_from_file(trim(filepath), eos%g1p_x, .false.)
        case("g1p_y")
            eos%g1p_y%filename = filename
            call read_eos_from_file(trim(filepath), eos%g1p_y, .false.)
        case("g1p_xy")
            eos%g1p_xy%filename = filename
            call read_eos_from_file(trim(filepath), eos%g1p_xy, .false.)
        case("eintT")
            eos%eintT%filename = filename
            call read_eos_from_file(trim(filepath), eos%eintT, .false.)
        case("eintT_x")
            eos%eintT_x%filename = filename
            call read_eos_from_file(trim(filepath), eos%eintT_x, .false.)
        case("eintT_y")
            eos%eintT_y%filename = filename
            call read_eos_from_file(trim(filepath), eos%eintT_y, .false.)
        case("eintT_xy")
            eos%eintT_xy%filename = filename
            call read_eos_from_file(trim(filepath), eos%eintT_xy, .false.)
        ! 'yT' tables dropped -- zero runtime references.
        case("Tfwd")
            eos%Tfwd%filename = filename
            call read_eos_from_file(trim(filepath), eos%Tfwd, .false.)
        case("Tfwd_x")
            eos%Tfwd_x%filename = filename
            call read_eos_from_file(trim(filepath), eos%Tfwd_x, .false.)
        case("Tfwd_y")
            eos%Tfwd_y%filename = filename
            call read_eos_from_file(trim(filepath), eos%Tfwd_y, .false.)
        case("Tfwd_xy")
            eos%Tfwd_xy%filename = filename
            call read_eos_from_file(trim(filepath), eos%Tfwd_xy, .false.)
        case("pfwd")
            eos%pfwd%filename = filename
            call read_eos_from_file(trim(filepath), eos%pfwd, .false.)
        case("pfwd_x")
            eos%pfwd_x%filename = filename
            call read_eos_from_file(trim(filepath), eos%pfwd_x, .false.)
        case("pfwd_y")
            eos%pfwd_y%filename = filename
            call read_eos_from_file(trim(filepath), eos%pfwd_y, .false.)
        case("pfwd_xy")
            eos%pfwd_xy%filename = filename
            call read_eos_from_file(trim(filepath), eos%pfwd_xy, .false.)
        case default
            call mpistop("eos table name "//trim(fieldname)//" not recognised in load_lte_tables")
        end select

    end subroutine load_lte_tables

    !> Attempt to load pre-computed table from binary file.
    !> If the file does not exist, silently skip -- the table will be built at runtime.
    subroutine try_load_lte_tables(fieldname)
        character(len=*), intent(in) :: fieldname
        character(len=std_len) :: subname, tablename, filename, filepath
        logical :: file_exists

        subname = "H"
        if (eos%He_abundance > 0.0d0) then
            write(subname,"(A,A)") trim(subname), "He"
        endif
        if (eos%ionE) then
            write(subname,"(A,A)") trim(subname), "_IonE"
        else
            write(subname,"(A,A)") trim(subname), "_NoIonE"
        endif
        write(tablename, "(A,A,A)") trim(eos_table_prefix), fieldname, "_"
        write(filename,"(A,A,A)") trim(tablename), trim(subname), '.bin'
        write(filepath,"(A,A)") trim(eos%table_location), trim(filename)

        inquire(file=trim(filepath), exist=file_exists)
        if (file_exists) then
            call load_lte_tables(fieldname)
        else
            if (mype == 0) write(*,*) &
                trim(fieldname)//' table not found, will build at runtime: '//trim(filename)
        endif
    end subroutine try_load_lte_tables

    !> Read an EoS table from a binary file written by generate_lte_tables.py.
    !>
    !> File format (uniform -- legacy):
    !>   [2 x int32]   : (dim1, dim2) = (dimy, dimx)   in Fortran column-major
    !>   [4 x float64] : (var1_min, var1_max, var2_min, var2_max)
    !>   [dim1*dim2 x float64] : table values
    !>
    !> File format (adaptive -- extension): same prefix, then optional trailer
    !>   [1 x int32]      : trailer_flag (= 1)
    !>   [dim1 x float64] : x_nodes (axis-1 node positions, log10 space)
    !>   [dim2 x float64] : y_nodes (axis-2 node positions, log10 space)
    !>
    !> The trailer is detected via end-of-file on the trailer_flag read:
    !> uniform tables (no trailer) leave is_uniform = .true. and produce
    !> ios /= 0 on the read; adaptive tables set is_uniform = .false. and
    !> populate x_nodes, y_nodes.
    subroutine read_eos_from_file(filename, tc, logtable)
        character(len=*), intent(in) :: filename
        type(eos_table_container), intent(inout) :: tc
        logical, intent(in) :: logtable

        integer :: ios, trailer_flag

        open(unit=10, file=filename, access='stream', form='unformatted', action='read')

        !> First read the header information, size from python file, and allocate
        read(10) tc%dim1, tc%dim2

        !> Then read the header information, bounds from python file
        read(10) tc%var1_min, tc%var1_max, tc%var2_min, tc%var2_max

        if (allocated(tc%table)) deallocate(tc%table)
        allocate(tc%table(tc%dim1, tc%dim2))

        !> Then read the 2D cube from the python file
        read(10) tc%table

        if (logtable) then
            tc%table = dlog10(tc%table) !> Store tables in log10 space for interpolation
        endif

        !> Optional adaptive-grid trailer. Uniform (legacy) tables hit EOF here
        !> and the iostat-guarded read keeps is_uniform = .true.
        tc%is_uniform = .true.
        if (allocated(tc%var1_nodes)) deallocate(tc%var1_nodes)
        if (allocated(tc%var2_nodes)) deallocate(tc%var2_nodes)

        read(10, iostat=ios) trailer_flag
        if (ios == 0) then
            if (trailer_flag == 1) then
                tc%is_uniform = .false.
                allocate(tc%var1_nodes(tc%dim1))
                allocate(tc%var2_nodes(tc%dim2))
                read(10) tc%var1_nodes
                read(10) tc%var2_nodes
                if (mype == 0) then
                    print*, trim(filename), " (adaptive grid trailer detected)"
                endif
            else
                if (mype == 0) then
                    print*, "WARNING: unrecognised trailer flag ", trailer_flag, &
                        " in ", trim(filename), " - treating as uniform"
                endif
            endif
        end if

        if (mype==0) then
            if (eos%table_check) then
                print*, filename, " table check"
                print*, tc%var1_min, tc%var1_max, tc%var2_min, tc%var2_max
                print*, 'corners', ' interpolated'
                print*, tc%table(1,1), interp_clamped_bilinear_table(tc%var1_min, &
                    tc%var2_min, tc%table, tc%dim1, tc%dim2, tc%var1_min, tc%var1_max, tc%var2_min, tc%var2_max)
                print*, tc%table(1,tc%dim2), interp_clamped_bilinear_table(tc%var1_min, &
                    tc%var2_max, tc%table, tc%dim1, tc%dim2, tc%var1_min, tc%var1_max, tc%var2_min, tc%var2_max)
                print*, tc%table(tc%dim1,1), interp_clamped_bilinear_table(tc%var1_max, &
                    tc%var2_min, tc%table, tc%dim1, tc%dim2, tc%var1_min, tc%var1_max, tc%var2_min, tc%var2_max)
                print*, tc%table(tc%dim1,tc%dim2), interp_clamped_bilinear_table(tc%var1_max, &
                    tc%var2_max, tc%table, tc%dim1, tc%dim2, tc%var1_min, tc%var1_max, tc%var2_min, tc%var2_max)
                print*, "These should be identical, if not then the table structure is inconsistent"
                print*, "###########"
            endif
        endif

        close(10)

    end subroutine read_eos_from_file

    !> Gamma_1 from (log nH, log eint/nH). Used only by build_gamma1_p_table,
    !> i.e. the 'tables' method; the entropy method reads its own g1p table.
    double precision function gamma1_from_nH_eint(log_nH, log_eint_nH) result(g1)
        double precision, intent(in) :: log_nH, log_eint_nH
        g1 = bicubic_lookup(log_nH, log_eint_nH, eos%gamma1)
    end function gamma1_from_nH_eint

    !> Build the first adiabatic index Gamma_1 table from T and neOnH tables.
    !> Called during eos_finalise after unit conversion, when tables are in code units.
    !>
    !> Physics: Gamma_1 = (rho/p) * (dp/drho)_s = a + b * (p / eint_vol)
    !>   where a = d(log10 p)/d(log10 nH) at constant eint/nH
    !>         b = d(log10 p)/d(log10 eint/nH) at constant nH
    !>   and p = nH * T * (1 + He + y), eint_vol = nH * (eint/nH)
    !>
    !> Derivation: from the thermodynamic identity
    !>   a^2 = (dp/drho)_{eint} + ((eint+p)/rho) * (dp/deint)_rho
    !> expressed in table coordinates (log10 nH, log10 eint/nH).
    !> For ideal gas: a=1, b=1, p/eint=gamma-1, so Gamma_1 = 1+(gamma-1) = gamma.
    subroutine build_gamma1_table()
        integer :: n1, n2, i, j, im, ip, jm, jp
        double precision :: x1, x2
        double precision, allocatable :: log_p(:,:)
        double precision :: T_val, y_val, p_val, eint_vol
        double precision :: a_val, b_val, g1_val
        double precision :: g1_min, g1_max

        n1 = eos%T%dim1
        n2 = eos%T%dim2

        !> Gamma1 table shares the same grid as T and neOnH (inherits
        !> uniform-or-adaptive node layout from the source).
        eos%gamma1%dim1 = n1
        eos%gamma1%dim2 = n2
        eos%gamma1%var1_min = eos%T%var1_min
        eos%gamma1%var1_max = eos%T%var1_max
        eos%gamma1%var2_min = eos%T%var2_min
        eos%gamma1%var2_max = eos%T%var2_max
        eos%gamma1%filename = 'computed_gamma1'
        eos%gamma1%is_uniform = eos%T%is_uniform
        if (allocated(eos%gamma1%var1_nodes)) deallocate(eos%gamma1%var1_nodes)
        if (allocated(eos%gamma1%var2_nodes)) deallocate(eos%gamma1%var2_nodes)
        allocate(eos%gamma1%var1_nodes(n1)); eos%gamma1%var1_nodes = eos%T%var1_nodes
        allocate(eos%gamma1%var2_nodes(n2)); eos%gamma1%var2_nodes = eos%T%var2_nodes

        if (allocated(eos%gamma1%table)) deallocate(eos%gamma1%table)
        allocate(eos%gamma1%table(n1, n2))
        allocate(log_p(n1, n2))

        !> Step 1: Compute log10(p) at each grid point (using actual node positions)
        !> In code units: p = nH * T * (1 + He + y), kB absorbed
        do j = 1, n2
            x2 = eos%T%var2_nodes(j)
            do i = 1, n1
                x1 = eos%T%var1_nodes(i)
                T_val = 10.0d0**eos%T%table(i, j)         !> T in code units
                y_val = 10.0d0**eos%neOnH%table(i, j)     !> Ne/nH (dimensionless)
                !> log10(p) = log10(nH) + log10(T) + log10(1 + He + y)
                log_p(i, j) = x1 + eos%T%table(i, j) + dlog10(1.0d0 + eos%He_abundance + y_val)
            end do
        end do

        !> Step 2: Compute Gamma_1 from finite differences of log_p using
        !> actual node spacings (correct for both uniform and adaptive grids).
        g1_min = 1.0d30
        g1_max = -1.0d30

        do j = 1, n2
            x2 = eos%T%var2_nodes(j)
            do i = 1, n1
                x1 = eos%T%var1_nodes(i)

                !> Centered finite differences with one-sided at boundaries
                im = max(i - 1, 1)
                ip = min(i + 1, n1)
                jm = max(j - 1, 1)
                jp = min(j + 1, n2)

                !> p/eint_vol = T * (1+He+y) / (eint/nH)
                T_val = 10.0d0**eos%T%table(i, j)
                y_val = 10.0d0**eos%neOnH%table(i, j)
                eint_vol = 10.0d0**x2  !> eint/nH in code units
                p_val = T_val * (1.0d0 + eos%He_abundance + y_val)

                if (eos%gamma1_method == 'effective') then
                    !> Approximate: gamma_eff = 1 + p/eint_vol
                    g1_val = 1.0d0 + p_val / eint_vol
                else
                    !> Full formula: Gamma_1 = a + b * (p / eint_vol)
                    !> a = d(log10 p) / d(log10 nH) at constant eint/nH
                    a_val = (log_p(ip, j) - log_p(im, j)) &
                          / (eos%T%var1_nodes(ip) - eos%T%var1_nodes(im))
                    !> b = d(log10 p) / d(log10 eint/nH) at constant nH
                    b_val = (log_p(i, jp) - log_p(i, jm)) &
                          / (eos%T%var2_nodes(jp) - eos%T%var2_nodes(jm))
                    g1_val = a_val + b_val * (p_val / eint_vol)
                endif

                !> Clamp to physical range [1, 5/3] (monatomic H+He gas)
                g1_val = max(1.001d0, min(g1_val, 5.0d0/3.0d0))

                eos%gamma1%table(i, j) = g1_val

                g1_min = min(g1_min, g1_val)
                g1_max = max(g1_max, g1_val)
            end do
        end do

        deallocate(log_p)

        if (mype == 0) then
            write(*, '(A,A,A,F8.4,A,F8.4)') &
                ' Gamma1 table built (', trim(eos%gamma1_method), '): min = ', &
                g1_min, ', max = ', g1_max
        end if

    end subroutine build_gamma1_table

    !> Build merged log10(p/nH) table: log10(p/nH)(nH, eint/nH).
    !> Same grid as T and neOnH tables. Each entry is:
    !>   log10(p/nH) = log10( 10^T * (1 + He + 10^y) )
    !> where T and y are from the forward tables at the same grid point.
    !> Used by the WB bisection to replace two separate PCHIP lookups
    !> (T + neOnH) with a single lookup per iteration.
    subroutine build_log_p_table()
        integer :: n1, n2, i, j
        double precision :: log_nH, log_eint_nH, T_val, y_val

        n1 = eos%T%dim1
        n2 = eos%T%dim2

        eos%log_p%dim1 = n1
        eos%log_p%dim2 = n2
        eos%log_p%var1_min = eos%T%var1_min
        eos%log_p%var1_max = eos%T%var1_max
        eos%log_p%var2_min = eos%T%var2_min
        eos%log_p%var2_max = eos%T%var2_max
        eos%log_p%is_uniform = eos%T%is_uniform
        if (allocated(eos%log_p%var1_nodes)) deallocate(eos%log_p%var1_nodes)
        if (allocated(eos%log_p%var2_nodes)) deallocate(eos%log_p%var2_nodes)
        allocate(eos%log_p%var1_nodes(n1)); eos%log_p%var1_nodes = eos%T%var1_nodes
        allocate(eos%log_p%var2_nodes(n2)); eos%log_p%var2_nodes = eos%T%var2_nodes

        if (allocated(eos%log_p%table)) deallocate(eos%log_p%table)
        allocate(eos%log_p%table(n1, n2))

        do j = 1, n2
            do i = 1, n1
                T_val = eos%T%table(i, j)
                y_val = eos%neOnH%table(i, j)
                eos%log_p%table(i, j) = dlog10(10.0d0**T_val &
                    * (1.0d0 + eos%He_abundance + 10.0d0**y_val))
            end do
        end do

        if (mype == 0) then
            write(*, '(A,I4,A,I4)') &
                ' Merged log_p table built: ', n1, ' x ', n2
        end if

    end subroutine build_log_p_table

    !> Build p/nH table: (1+He+y)*T at each (nH, eint/nH) grid point.
    !> Shares axes with the T and neOnH tables (Group A).
    !> Used by to_primitive_LTE for single-lookup pressure computation.
    subroutine build_p_over_nH_table()
        integer :: n1, n2, i, j
        double precision :: T_val, y_val, p_nH_val
        double precision :: p_min, p_max

        n1 = eos%T%dim1
        n2 = eos%T%dim2

        eos%p_over_nH%dim1 = n1
        eos%p_over_nH%dim2 = n2
        eos%p_over_nH%var1_min = eos%T%var1_min
        eos%p_over_nH%var1_max = eos%T%var1_max
        eos%p_over_nH%var2_min = eos%T%var2_min
        eos%p_over_nH%var2_max = eos%T%var2_max
        eos%p_over_nH%filename = 'computed_p_over_nH'
        eos%p_over_nH%is_uniform = eos%T%is_uniform
        if (allocated(eos%p_over_nH%var1_nodes)) deallocate(eos%p_over_nH%var1_nodes)
        if (allocated(eos%p_over_nH%var2_nodes)) deallocate(eos%p_over_nH%var2_nodes)
        allocate(eos%p_over_nH%var1_nodes(n1)); eos%p_over_nH%var1_nodes = eos%T%var1_nodes
        allocate(eos%p_over_nH%var2_nodes(n2)); eos%p_over_nH%var2_nodes = eos%T%var2_nodes

        if (allocated(eos%p_over_nH%table)) deallocate(eos%p_over_nH%table)
        allocate(eos%p_over_nH%table(n1, n2))

        p_min = 1.0d30
        p_max = -1.0d30

        do j = 1, n2
            do i = 1, n1
                !> T table stores log10(T_code), neOnH stores log10(y)
                T_val = 10.0d0**eos%T%table(i, j)
                y_val = 10.0d0**eos%neOnH%table(i, j)
                !> p/nH = (1 + He + y) * T in code units
                p_nH_val = (1.0d0 + eos%He_abundance + y_val) * T_val
                !> Store as log10 for PCHIP interpolation consistency
                eos%p_over_nH%table(i, j) = dlog10(p_nH_val)
                p_min = min(p_min, p_nH_val)
                p_max = max(p_max, p_nH_val)
            end do
        end do

        call precompute_step_inv(eos%p_over_nH)

        if (mype == 0) then
            write(*, '(A,ES10.3,A,ES10.3)') &
                ' p/nH table built: min = ', p_min, ', max = ', p_max
        end if

    end subroutine build_p_over_nH_table

    !> Build interleaved Group A table from existing T, neOnH, p_over_nH tables.
    !> Layout: table_eint_il(3, ny, nx) where:
    !>   slot 1 = T (log10), slot 2 = neOnH (log10), slot 3 = p_over_nH (log10)
    !> All three are on the same (nH, eint/nH) axes.
    subroutine build_interleaved_eint_table()
        integer :: n1, n2, i, j

        n1 = eos%T%dim1
        n2 = eos%T%dim2

        if (allocated(eos%table_eint_il)) deallocate(eos%table_eint_il)
        allocate(eos%table_eint_il(3, n1, n2))

        do j = 1, n2
            do i = 1, n1
                eos%table_eint_il(1, i, j) = eos%T%table(i, j)
                eos%table_eint_il(2, i, j) = eos%neOnH%table(i, j)
                eos%table_eint_il(3, i, j) = eos%p_over_nH%table(i, j)
            end do
        end do

        if (mype == 0) then
            write(*, '(A,I0,A,I0,A,I0,A,F6.1,A)') &
                ' Interleaved Group A table built: ', 3, ' x ', n1, ' x ', n2, &
                ' (', 3.0d0*n1*n2*8.0d0/1024.0d0, ' KB)'
        end if

    end subroutine build_interleaved_eint_table

    !> Build pressure-indexed Gamma_1 table: Gamma_1(nH, p/nH).
    !> Re-indexes the eint-based gamma1 table into pressure space,
    !> sharing the same grid axes as the p2eint table.
    !> This eliminates the intermediate p2eint lookup from hd_get_csound2_LTE.
    subroutine build_gamma1_p_table()
        integer :: n1, n2, i, j
        double precision :: log_nH, log_p_nH
        double precision :: p2eint_ratio, log_eint_nH
        double precision :: g1_val, g1_min, g1_max

        n1 = eos%p2eint%dim1
        n2 = eos%p2eint%dim2

        ! gamma1_p shares the same axis layout as p2eint (inherits adaptive
        ! v1 nodes from T via p2eint; v2 is uniform over the log_p/nH range).
        eos%gamma1_p%dim1 = n1
        eos%gamma1_p%dim2 = n2
        eos%gamma1_p%var1_min = eos%p2eint%var1_min
        eos%gamma1_p%var1_max = eos%p2eint%var1_max
        eos%gamma1_p%var2_min = eos%p2eint%var2_min
        eos%gamma1_p%var2_max = eos%p2eint%var2_max
        eos%gamma1_p%filename = 'computed_gamma1_p'
        eos%gamma1_p%is_uniform = eos%p2eint%is_uniform
        if (allocated(eos%gamma1_p%var1_nodes)) deallocate(eos%gamma1_p%var1_nodes)
        if (allocated(eos%gamma1_p%var2_nodes)) deallocate(eos%gamma1_p%var2_nodes)
        allocate(eos%gamma1_p%var1_nodes(n1)); eos%gamma1_p%var1_nodes = eos%p2eint%var1_nodes
        allocate(eos%gamma1_p%var2_nodes(n2)); eos%gamma1_p%var2_nodes = eos%p2eint%var2_nodes

        if (allocated(eos%gamma1_p%table)) deallocate(eos%gamma1_p%table)
        allocate(eos%gamma1_p%table(n1, n2))

        g1_min = 1.0d30
        g1_max = -1.0d30

        do j = 1, n2
            log_p_nH = eos%p2eint%var2_nodes(j)
            do i = 1, n1
                log_nH = eos%p2eint%var1_nodes(i)

                ! Convert from (nH, p/nH) to (nH, eint/nH) via the p2eint table
                ! p2eint table stores the ratio eint/(p), so:
                ! log10(eint/nH) = log10(p/nH) + log10(p2eint_ratio)
                p2eint_ratio = eos%p2eint%table(i, j)
                log_eint_nH = log_p_nH + dlog10(p2eint_ratio)

                ! Look up Gamma_1 from the eint-indexed table
                g1_val = gamma1_from_nH_eint(log_nH, log_eint_nH)

                eos%gamma1_p%table(i, j) = g1_val

                g1_min = min(g1_min, g1_val)
                g1_max = max(g1_max, g1_val)
            end do
        end do

        if (mype == 0) then
            write(*, '(A,F8.4,A,F8.4)') &
                ' Gamma1_p table built (pressure-indexed): min = ', g1_min, ', max = ', g1_max
        end if

    end subroutine build_gamma1_p_table

    !> Build p2eint table at runtime by inverting the forward (T, neOnH) tables.
    !> For each (nH, p/nH) grid point, find eint/nH by bisection such that
    !> the PCHIP-interpolated forward tables give p(nH, eint/nH) = p_target.
    !> This guarantees exact round-trip: to_conserved(p) -> eint -> to_primitive -> p.
    !>
    !> The p/nH grid bounds are computed from the forward tables (not hardcoded),
    !> ensuring coverage of all physically realisable pressures.
    subroutine build_p2eint_table()
        integer :: n1, n2_fwd, n2_p, i, j, iter
        double precision :: log_nH_i, log_eint_lo, log_eint_hi, log_eint_mid
        double precision :: log_p_target, log_p_eval
        double precision :: T_val, y_val
        double precision :: p_global_min, p_global_max
        double precision :: dx_p, log_p_nH_ij
        double precision :: max_err, err
        double precision :: log_p_lo_i, log_p_hi_i
        integer :: i_worst, j_worst

        if (mype == 0) write(*,*) 'Building p2eint table from forward tables...'

        n1 = eos%T%dim1
        n2_fwd = eos%T%dim2

        !> Compute p/nH at every forward grid node to find the pressure range
        p_global_min =  1.0d30
        p_global_max = -1.0d30
        do j = 1, n2_fwd
            do i = 1, n1
                T_val = 10.0d0**eos%T%table(i, j)
                y_val = 10.0d0**eos%neOnH%table(i, j)
                log_p_nH_ij = dlog10(T_val * (1.0d0 + eos%He_abundance + y_val))
                p_global_min = min(p_global_min, log_p_nH_ij)
                p_global_max = max(p_global_max, log_p_nH_ij)
            end do
        end do

        !> Pad pressure range to avoid edge clamping at off-grid PCHIP evaluations.
        !> The PCHIP polynomial can overshoot grid-node extrema by ~h^2 * f'',
        !> which for steep ionisation-zone gradients needs generous padding.
        p_global_min = p_global_min - 0.2d0
        p_global_max = p_global_max + 0.2d0

        n2_p = n2_fwd  !> Same resolution as forward table
        if (allocated(eos%p2eint%table)) deallocate(eos%p2eint%table)
        eos%p2eint%dim1 = n1
        eos%p2eint%dim2 = n2_p
        eos%p2eint%var1_min = eos%T%var1_min
        eos%p2eint%var1_max = eos%T%var1_max
        eos%p2eint%var2_min = p_global_min
        eos%p2eint%var2_max = p_global_max
        eos%p2eint%filename = 'computed_p2eint'
        !> v1 axis inherited from T (adaptive if T is); v2 axis is uniform
        !> over the (different) log_p/nH range. When T is adaptive, we still
        !> store explicit var2_nodes so lookups dispatch via _nu uniformly.
        eos%p2eint%is_uniform = eos%T%is_uniform
        if (allocated(eos%p2eint%var1_nodes)) deallocate(eos%p2eint%var1_nodes)
        if (allocated(eos%p2eint%var2_nodes)) deallocate(eos%p2eint%var2_nodes)
        allocate(eos%p2eint%var1_nodes(n1)); eos%p2eint%var1_nodes = eos%T%var1_nodes
        allocate(eos%p2eint%var2_nodes(n2_p))
        do j = 1, n2_p
            eos%p2eint%var2_nodes(j) = p_global_min &
                + (j - 1) * (p_global_max - p_global_min) / dble(n2_p - 1)
        end do
        allocate(eos%p2eint%table(n1, n2_p))

        dx_p = (p_global_max - p_global_min) / dble(n2_p - 1)
        max_err = 0.0d0

        !> For each (nH_i, p_j/nH), bisect on log10(eint/nH) to find the
        !> value where the PCHIP-interpolated forward tables give this pressure.
        i_worst = 1
        j_worst = 1
        do i = 1, n1
            log_nH_i = eos%T%var1_nodes(i)

            !> Compute achievable p range for this nH from forward table endpoints
            T_val = bicubic_lookup(log_nH_i, eos%T%var2_min, eos%T)
            y_val = bicubic_lookup(log_nH_i, eos%T%var2_min, eos%neOnH)
            log_p_lo_i = dlog10(10.0d0**T_val &
                * (1.0d0 + eos%He_abundance + 10.0d0**y_val))
            T_val = bicubic_lookup(log_nH_i, eos%T%var2_max, eos%T)
            y_val = bicubic_lookup(log_nH_i, eos%T%var2_max, eos%neOnH)
            log_p_hi_i = dlog10(10.0d0**T_val &
                * (1.0d0 + eos%He_abundance + 10.0d0**y_val))

            do j = 1, n2_p
                log_p_target = p_global_min + (j - 1) * dx_p

                !> Clamp target to achievable range -- outside this range,
                !> map to the table boundary (best available)
                if (log_p_target <= log_p_lo_i) then
                    eos%p2eint%table(i, j) = 10.0d0**(eos%T%var2_min - log_p_target)
                    cycle
                else if (log_p_target >= log_p_hi_i) then
                    eos%p2eint%table(i, j) = 10.0d0**(eos%T%var2_max - log_p_target)
                    cycle
                end if

                !> Bisection on log10(eint/nH)
                log_eint_lo = eos%T%var2_min
                log_eint_hi = eos%T%var2_max

                do iter = 1, 52  !> 52 bisections -> 2^-52 ~ 10^-15.7 precision
                    log_eint_mid = 0.5d0 * (log_eint_lo + log_eint_hi)

                    !> Evaluate p from forward tables at (nH_i, eint_mid)
                    T_val = bicubic_lookup(log_nH_i, log_eint_mid, eos%T)
                    y_val = bicubic_lookup(log_nH_i, log_eint_mid, eos%neOnH)
                    log_p_eval = dlog10(10.0d0**T_val &
                        * (1.0d0 + eos%He_abundance + 10.0d0**y_val))

                    if (log_p_eval < log_p_target) then
                        log_eint_lo = log_eint_mid
                    else
                        log_eint_hi = log_eint_mid
                    end if

                    if (dabs(log_eint_hi - log_eint_lo) < 1.0d-14) exit
                end do

                eos%p2eint%table(i, j) = 10.0d0**(log_eint_mid - log_p_target)

                !> Track round-trip error (only for in-range points)
                err = dabs(log_p_eval - log_p_target)
                if (err > max_err) then
                    max_err = err
                    i_worst = i
                    j_worst = j
                end if
            end do
        end do

        if (mype == 0) then
            write(*, '(A,ES10.3,A,I4,A,I4)') &
                ' p2eint table built: max bisection err = ', max_err, &
                ' at (i,j)=(', i_worst, ',', j_worst, ')'
            write(*, '(A,F8.4,A,F8.4)') &
                '   log10(p/nH) range = ', p_global_min, ' to ', p_global_max
        end if

    end subroutine build_p2eint_table

    !> Build inverse T table: given (nH, T), return log10(eint/nH).
    !> Inverts the forward T(nH, eint/nH) table by bisection using the SAME
    !> PCHIP interpolation kernel as T_from_nH_eint. This guarantees that
    !> T_from_nH_eint(nH, eint_from_T(nH, T)) = T to machine precision.
    !> Called during eos_finalise after unit conversion.
    subroutine build_eint_from_T_table()
        integer :: n1, n2_fwd, n2_inv, i, j, iter
        double precision :: dx2_inv
        double precision :: T_global_min, T_global_max, T_FI_threshold
        double precision :: log_nH_i, log_eint_lo, log_eint_hi, log_eint_mid
        double precision :: T_target, T_eval, T_max_at_nH, eint_nH_FI
        double precision :: eint_nH_min, eint_nH_max, max_err, err
        integer :: n_FI_filled

        n1 = eos%T%dim1
        n2_fwd = eos%T%dim2

        !> Find global T range from the forward T table
        T_global_min =  1.0d30
        T_global_max = -1.0d30
        do j = 1, n2_fwd
            do i = 1, n1
                T_global_min = min(T_global_min, eos%T%table(i, j))
                T_global_max = max(T_global_max, eos%T%table(i, j))
            end do
        end do

        !> FI threshold temperature in code units (200,000 K)
        T_FI_threshold = eos%p_rho_FI_threshold &
            * (1.0d0 + 4.0d0*eos%He_abundance) / eos%n_per_nH_FI

        !> Extend T axis to cover the FI regime up to 10 MK (or T_global_max if larger)
        T_global_max = max(T_global_max, dlog10(10.0d0))  !> 10 MK in code units (unit_T=1MK)

        !> Set up inverse table with same nH axis, extended T axis
        n2_inv = n2_fwd  !> Same resolution as forward table
        if (allocated(eos%eint_from_T%table)) deallocate(eos%eint_from_T%table)
        eos%eint_from_T%dim1 = n1
        eos%eint_from_T%dim2 = n2_inv
        eos%eint_from_T%var1_min = eos%T%var1_min
        eos%eint_from_T%var1_max = eos%T%var1_max
        eos%eint_from_T%var2_min = T_global_min
        eos%eint_from_T%var2_max = T_global_max
        eos%eint_from_T%filename = 'computed_eint_from_T'
        !> Inherit v1 axis from T (adaptive if T is); uniform v2 axis covers
        !> [T_global_min, T_global_max] with explicit nodes.
        eos%eint_from_T%is_uniform = eos%T%is_uniform
        if (allocated(eos%eint_from_T%var1_nodes)) deallocate(eos%eint_from_T%var1_nodes)
        if (allocated(eos%eint_from_T%var2_nodes)) deallocate(eos%eint_from_T%var2_nodes)
        allocate(eos%eint_from_T%var1_nodes(n1)); eos%eint_from_T%var1_nodes = eos%T%var1_nodes
        allocate(eos%eint_from_T%var2_nodes(n2_inv))
        do j = 1, n2_inv
            eos%eint_from_T%var2_nodes(j) = T_global_min &
                + (j - 1) * (T_global_max - T_global_min) / dble(n2_inv - 1)
        end do

        allocate(eos%eint_from_T%table(n1, n2_inv))

        dx2_inv = (T_global_max - T_global_min) / dble(n2_inv - 1)

        eint_nH_min =  1.0d30
        eint_nH_max = -1.0d30
        max_err = 0.0d0
        n_FI_filled = 0

        !> For each (nH_i, T_j): use bisection on the forward table where it
        !> has coverage, and the analytic FI formula above the table range.
        do i = 1, n1
            log_nH_i = eos%T%var1_nodes(i)

            !> Find the maximum T in the forward table at this nH
            T_max_at_nH = eos%T%table(i, n2_fwd)

            do j = 1, n2_inv
                T_target = T_global_min + (j - 1) * dx2_inv

                if (10.0d0**T_target > T_FI_threshold .or. &
                    T_target > T_max_at_nH) then
                    !> Above table coverage or FI regime: use analytic formula
                    !> eint/nH = n_per_nH_FI / (gamma-1) * T + neOnH_FI * eion_per_nH
                    eint_nH_FI = eos%n_per_nH_FI * eos%inv_gamma_minus_1 &
                        * 10.0d0**T_target
                    if (eos%ionE) eint_nH_FI = eint_nH_FI &
                        + eos%neOnH_FI * eos%eion_per_nH
                    eos%eint_from_T%table(i, j) = dlog10(eint_nH_FI)
                    n_FI_filled = n_FI_filled + 1
                else
                    !> Bisection on log10(eint/nH) using the forward T table
                    log_eint_lo = eos%T%var2_min
                    log_eint_hi = eos%T%var2_max

                    do iter = 1, 52
                        log_eint_mid = 0.5d0 * (log_eint_lo + log_eint_hi)

                        T_eval = bicubic_lookup(log_nH_i, log_eint_mid, eos%T)

                        if (T_eval < T_target) then
                            log_eint_lo = log_eint_mid
                        else
                            log_eint_hi = log_eint_mid
                        end if

                        if (dabs(log_eint_hi - log_eint_lo) < 1.0d-14) exit
                    end do

                    eos%eint_from_T%table(i, j) = log_eint_mid

                    err = dabs(T_eval - T_target)
                    max_err = max(max_err, err)
                end if

                eint_nH_min = min(eint_nH_min, eos%eint_from_T%table(i, j))
                eint_nH_max = max(eint_nH_max, eos%eint_from_T%table(i, j))
            end do
        end do

        if (mype == 0) then
            write(*, '(A,ES10.3,A,I0,A,I0)') &
                ' Inverse T table built: max bisection err = ', max_err, &
                ', FI-filled = ', n_FI_filled, ' / ', n1*n2_inv
            write(*, '(A,F8.4,A,F8.4)') &
                '   log10(T) range = ', T_global_min, ' to ', T_global_max
            write(*, '(A,F8.4,A,F8.4)') &
                '   log10(eint/nH) range = ', eint_nH_min, ' to ', eint_nH_max
        end if

    end subroutine build_eint_from_T_table

    !> Verify round-trip consistency of ALL EoS table pathways.
    !> Tests at off-grid random points to measure real interpolation error.
    !> Called after all tables are built.
    subroutine verify_eos_round_trips()
        integer :: n1, n2, i, j, n_tested
        double precision :: log_nH, log_eint_nH
        double precision :: dx1, dx2
        double precision :: T_val, y_val, log_p_nH, p2eint_val, log_eint_recovered
        double precision :: T_recovered, p_from_Ty, log_p_recovered
        double precision :: eint_from_T_val, log_eint_from_T_recovered
        double precision :: err_p, err_eint_T
        double precision :: max_err_p, max_err_eint_T, mean_err_p, mean_err_eint_T
        double precision :: log_nH_worst_p, log_eint_worst_p

        n1 = eos%T%dim1
        n2 = eos%T%dim2
        dx1 = (eos%T%var1_max - eos%T%var1_min) / dble(n1 - 1)
        dx2 = (eos%T%var2_max - eos%T%var2_min) / dble(n2 - 1)

        max_err_p = 0.0d0
        max_err_eint_T = 0.0d0
        mean_err_p = 0.0d0
        mean_err_eint_T = 0.0d0
        n_tested = 0

        !> Test at cell-centre-like points (offset by 0.5*dx from grid nodes)
        !> to exercise the PCHIP interpolation between nodes.
        do i = 2, n1 - 1
            log_nH = eos%T%var1_min + (dble(i) - 0.5d0) * dx1
            do j = 2, n2 - 1
                log_eint_nH = eos%T%var2_min + (dble(j) - 0.5d0) * dx2
                n_tested = n_tested + 1

                !> === Round-trip 1: eint -> T,y -> p -> p2eint -> eint' ===
                !> Forward: get T and y from eint
                T_val = bicubic_lookup(log_nH, log_eint_nH, eos%T)
                y_val = bicubic_lookup(log_nH, log_eint_nH, eos%neOnH)
                !> Compute p/nH
                log_p_nH = dlog10(10.0d0**T_val &
                    * (1.0d0 + eos%He_abundance + 10.0d0**y_val))
                !> Inverse: get eint from p via p2eint table
                p2eint_val = bicubic_lookup(log_nH, log_p_nH, eos%p2eint)
                log_eint_recovered = log_p_nH + dlog10(p2eint_val)

                err_p = dabs(log_eint_recovered - log_eint_nH)
                mean_err_p = mean_err_p + err_p
                if (err_p > max_err_p) then
                    max_err_p = err_p
                    log_nH_worst_p = log_nH
                    log_eint_worst_p = log_eint_nH
                end if

                !> === Round-trip 2: eint -> T -> eint_from_T -> eint' ===
                eint_from_T_val = bicubic_lookup(log_nH, T_val, eos%eint_from_T)

                err_eint_T = dabs(eint_from_T_val - log_eint_nH)
                mean_err_eint_T = mean_err_eint_T + err_eint_T
            end do
        end do

        mean_err_p = mean_err_p / dble(max(n_tested, 1))
        mean_err_eint_T = mean_err_eint_T / dble(max(n_tested, 1))

        if (mype == 0) then
            write(*, '(A)') ' === EoS round-trip verification (off-grid points) ==='
            write(*, '(A,I6,A)') '   Tested ', n_tested, ' points (half-cell offsets)'
            write(*, '(A,ES10.3,A,ES10.3)') &
                '   p round-trip:      max err = ', max_err_p, &
                '  mean err = ', mean_err_p
            write(*, '(A,F8.4,A,F8.4)') &
                '     worst at log_nH = ', log_nH_worst_p, &
                ', log_eint = ', log_eint_worst_p
            write(*, '(A,ES10.3,A,ES10.3)') &
                '   T round-trip:      max err = ', max_err_eint_T, &
                '  mean err = ', mean_err_eint_T
            write(*, '(A)') ' ====================================================='
        end if

    end subroutine verify_eos_round_trips

end module mod_eos_tables
!> Needs a line after to pass the preprocessor
