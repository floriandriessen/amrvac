!> Module for escape probability radiative cooling modification.
!>
!> Computes column mass along the 1D domain and provides an escape probability
!> factor beta(tau) = (1 - exp(-tau))/tau to multiply the radiative cooling rate.
!> This replaces the ad-hoc density taper with a physically-motivated opacity proxy
!> based on the formalism of Carlsson & Leenaarts (2012).
!>
!> Column mass is measured from the domain centre (apex of a symmetric loop)
!> TOWARD the footpoints. This gives zero column mass at the apex (full coronal
!> cooling) and large column mass at the footpoints (suppressed chromospheric
!> cooling).
!>
!> Currently supports 1D slab_uniform geometry with AMR (Phase 1).
module mod_escape_probability
  use mod_global_parameters, only: std_len
  implicit none
  private

  ! Public routines
  public :: escape_prob_init
  public :: escape_prob_compute_colmass

  ! Module state
  integer :: n_aux = 0                !< Total cells at finest AMR level
  double precision :: dx_aux = 0.0d0  !< Cell width at finest level
  integer :: iw_colmass_ = -1         !< Index into wextra for column mass
  logical :: escape_sym_ = .false.    !< Symmetric loop (apex at domain center)
  logical :: arrays_allocated = .false. !< Deferred allocation flag

  ! Work arrays (allocated on first use, after mesh params are set)
  double precision, allocatable :: rho_dx_local(:)
  double precision, allocatable :: rho_dx_global(:)
  double precision, allocatable :: colmass_left(:)
  double precision, allocatable :: colmass_right(:)

contains

  !> Register escape probability parameters.
  !> Called during hd_phys_init (before mesh parameters are available).
  !> Grid-dependent allocation is deferred to the first compute call.
  subroutine escape_prob_init(iw_colmass, escape_sym)
    use mod_global_parameters
    integer, intent(in) :: iw_colmass
    logical, intent(in) :: escape_sym

    iw_colmass_ = iw_colmass
    escape_sym_ = escape_sym

    {^NOONED
    call mpistop("escape probability currently only supports 1D")
    }

  end subroutine escape_prob_init

  !> Allocate work arrays once mesh parameters are available.
  subroutine escape_prob_setup_arrays()
    use mod_global_parameters
    integer :: level_factor

    {^IFONED
    level_factor = 2**(refine_max_level - 1)
    n_aux = domain_nx1 * level_factor
    dx_aux = (xprobmax1 - xprobmin1) / dble(n_aux)
    }

    allocate(rho_dx_local(n_aux))
    allocate(rho_dx_global(n_aux))
    allocate(colmass_left(n_aux))
    allocate(colmass_right(n_aux))
    arrays_allocated = .true.

    if (mype == 0) then
      write(*, '(A,I8,A,ES10.3,A,L1)') &
        ' Escape probability: n_aux=', n_aux, &
        ' dx_aux=', dx_aux, ' symmetric=', escape_sym_
    end if

  end subroutine escape_prob_setup_arrays

  !> Compute column mass from the corona toward the footpoints and scatter
  !> to wextra.  Called once per timestep before advance().
  !>
  !> For symmetric loops (escape_sym_=.true.): uses the domain centre as
  !> the apex and integrates outward toward both footpoints. This avoids
  !> the apex jumping when condensations form (density minimum is not stable).
  !> For non-symmetric atmospheres: integrates from the top of the domain
  !> (xprobmax) downward.
  subroutine escape_prob_compute_colmass()
    use mod_global_parameters
    use mod_variables, only: iw_rho
    use mod_connectivity, only: igrids_active, igridstail_active

    integer :: iigrid, igrid, level, ratio
    integer :: ix1, k, klo, khi, kk, k_apex
    double precision :: x_cell, rho_cell, running_sum
    double precision :: cm_value

    ! Deferred allocation on first call
    if (.not. arrays_allocated) call escape_prob_setup_arrays()
    if (n_aux == 0) return

    {^IFONED
    ! 1. Zero local accumulator
    rho_dx_local = 0.0d0

    ! 2. Gather: loop over active blocks, deposit density onto aux array
    do iigrid = 1, igridstail_active
      igrid = igrids_active(iigrid)
      level = node(plevel_, igrid)
      ratio = 2**(refine_max_level - level)

      do ix1 = ixMlo1, ixMhi1
        x_cell = ps(igrid)%x(ix1, 1)
        rho_cell = ps(igrid)%w(ix1, iw_rho)

        ! Map cell center to aux index (1-based)
        k = floor((x_cell - xprobmin1) / dx_aux) + 1
        k = max(1, min(k, n_aux))

        if (ratio == 1) then
          ! Fine level: 1:1 mapping
          rho_dx_local(k) = rho_dx_local(k) + rho_cell * dx_aux
        else
          ! Coarse level: spread density across 'ratio' fine cells.
          ! The coarse cell spans 'ratio' fine cells centered on k.
          klo = k - ratio/2 + 1
          khi = klo + ratio - 1
          klo = max(1, klo)
          khi = min(n_aux, khi)
          do kk = klo, khi
            rho_dx_local(kk) = rho_dx_local(kk) + rho_cell * dx_aux
          end do
        end if
      end do
    end do

    ! 3. MPI global reduction
    call MPI_ALLREDUCE(rho_dx_local, rho_dx_global, n_aux, &
         MPI_DOUBLE_PRECISION, MPI_SUM, icomm, ierrmpi)

    ! 3b. Fill gaps left by AMR level transitions.
    !     At coarse/fine boundaries, the interior-cell gather can miss
    !     a few aux cells (ghost cell region between blocks at different
    !     levels).  Forward-fill from the nearest populated neighbour.
    do k = 2, n_aux
      if (rho_dx_global(k) == 0.0d0) rho_dx_global(k) = rho_dx_global(k-1)
    end do
    ! Backward pass for any leading zeros (cells before the first block)
    do k = n_aux - 1, 1, -1
      if (rho_dx_global(k) == 0.0d0) rho_dx_global(k) = rho_dx_global(k+1)
    end do

    ! 4. Prefix sums for column mass from both boundaries
    ! Left: colmass_left(i) = integral of rho*dx from x=xprobmin to cell i
    running_sum = 0.0d0
    do k = 1, n_aux
      running_sum = running_sum + rho_dx_global(k)
      colmass_left(k) = running_sum
    end do

    ! Right: colmass_right(i) = integral of rho*dx from cell i to x=xprobmax
    running_sum = 0.0d0
    do k = n_aux, 1, -1
      running_sum = running_sum + rho_dx_global(k)
      colmass_right(k) = running_sum
    end do

    ! 5. Set apex index for column mass reference point
    if (escape_sym_) then
      ! Symmetric loop: apex is at the domain centre (fixed)
      k_apex = n_aux / 2
    else
      ! Non-symmetric: "corona" is at the top of the domain
      k_apex = n_aux
    end if

    ! 6. Scatter: write column mass FROM APEX TOWARD FOOTPOINTS to wextra.
    !    For symmetric: cm(x) = integral from apex to x (outward from corona)
    !    This gives 0 at apex, increasing toward both footpoints.
    do iigrid = 1, igridstail_active
      igrid = igrids_active(iigrid)

      do ix1 = ixGlo1, ixGhi1
        x_cell = ps(igrid)%x(ix1, 1)
        k = floor((x_cell - xprobmin1) / dx_aux) + 1
        k = max(1, min(k, n_aux))

        if (escape_sym_) then
          ! Symmetric loop: column mass from apex toward each footpoint
          if (k <= k_apex) then
            ! Left half: cm = colmass from k to k_apex
            cm_value = colmass_left(k_apex) - colmass_left(k)
          else
            ! Right half: cm = colmass from k_apex to k
            cm_value = colmass_right(k_apex) - colmass_right(k)
          end if
        else
          ! Non-symmetric: cm from top of domain downward = colmass_right
          cm_value = colmass_right(k)
        end if

        ps(igrid)%wextra(ix1, iw_colmass_) = max(cm_value, 0.0d0)
      end do
    end do
    }

  end subroutine escape_prob_compute_colmass

end module mod_escape_probability
