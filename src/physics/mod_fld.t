!> Module for flux limited diffusion (FLD)-approximation in Radiation-(Magneto)hydrodynamics simulations
!>
!> Full description of RHD-FLD in
!> Moens N., Sundqvist J.O., El Mellah I., Poniatowski L., Teunissen J. & Keppens R. 2022, A&A 657, A81
!> Radiation-hydrodynamics with MPI-AMRVAC . Flux-limited diffusion
!> doi:10.1051/0004-6361/202141023
!>
!> Full description for RMHD-FLD in
!> N. Narechania, R. Keppens, A. ud-Doula, N. Moens & J. Sundqvist 2025, A&A 696, A131
!> doi:10.1051/0004-6361/202452208
!> Radiation-magnetohydrodynamics with MPI-AMRVAC using flux-limited diffusion

module mod_fld
    use mod_comm_lib, only: mpistop
    use mod_geometry
    implicit none
    !> source split for energy interact and radforce:
    logical :: fld_Radforce_split = .false.
    !> Opacity per unit of unit_density
    double precision, public :: fld_kappa0 = 0.34d0
    !> Tolerance for bisection method for Energy sourceterms
    !> This is a percentage of the minimum of gas- and radiation energy
    double precision, public :: fld_bisect_tol = 1.d-4
    !> Tolerance for adi method for radiative Energy diffusion
    double precision, public :: fld_diff_tol = 1.d-4
    !> switches for opacity
    character(len=8)  :: fld_opacity_law = 'const'
    character(len=50) :: fld_opal_table = 'Y09800' 
    !> flux limiter choice
    character(len=16) :: fld_fluxlimiter = 'Pomraning'
    !> diffusion coefficient for multigrid method
    integer :: i_diff_mg
    !> Which method to find the root for the energy interaction polynomial
    character(len=8) :: fld_interaction_method = 'Halley'
    !> Resume run when multigrid returns error
    logical :: diffcrash_resume = .true.
    !> A copy of (m)hd_gamma
    double precision, private, protected :: fld_gamma
    !> running timestep for diffusion solver, initialised as zero
    double precision :: dt_diff = 0.d0
    !> public methods
    !> these are called in mod_hd_phys or mod_mhd_phys
    public :: fld_init
    public :: fld_get_radpress
    public :: fld_get_diffcoef_central
    public :: add_fld_rad_force
    public :: fld_radforce_get_dt
    ! these are made public for mod_usr purposes and diagnostics
    public :: fld_get_radflux
    public :: fld_get_fluxlimiter
    public :: fld_get_fluxlimiter_prim
    public :: fld_get_opacity_prim
  contains

  !> Reading in fld-list parameters from .par file
  subroutine fld_params_read(files)
    use mod_global_parameters, only: unitpar
    use mod_constants
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /fld_list/ fld_kappa0, fld_Radforce_split, &
    fld_bisect_tol, fld_diff_tol, fld_opacity_law, fld_fluxlimiter, &
    fld_interaction_method, diffcrash_resume, fld_opal_table

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, fld_list, end=111)
       111    close(unitpar)
    end do
  end subroutine fld_params_read

  !> Initialising FLD-module
  !> Read opacities
  !> Initialise Multigrid and adimensionalise kappa
  subroutine fld_init(He_abundance, r_gamma)
    use mod_global_parameters
    use mod_variables
    use mod_physics
    use mod_opal_opacity, only: init_opal_table
    use mod_multigrid_coupling

    double precision, intent(in) :: He_abundance, r_gamma
    double precision :: sigma_thomson

    call fld_params_read(par_files)
    phys_implicit_update => fld_implicit_update
    phys_evaluate_implicit => fld_evaluate_implicit
    use_multigrid = .true.
    mg%n_extra_vars = 1
    mg%operator_type = mg_vhelmholtz
    i_diff_mg = var_set_extravar("D", "D")
    !> set gamma
    fld_gamma = r_gamma
    !> Read in opacity table if necesary
    if(fld_opacity_law .eq. 'opal') call init_opal_table(fld_opal_table)
    if(fld_opacity_law .eq. 'thomson')  then
      ! special fixed value for thomson opacity, assuming cgs units
      sigma_thomson = 6.6524585d-25
      fld_kappa0 = sigma_thomson/const_mp * (1.+2.*He_abundance)/(1.+4.*He_abundance)
    endif
  end subroutine fld_init

  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  !> This subroutine handles the radiation force and its work added explicitly
  !> and the energy interaction term combined with photon tiring using an implicit update
  subroutine add_fld_rad_force(qdt,ixI^L,ixO^L,wCT,wCTprim,w,x,qsourcesplit,active)
    use mod_constants
    use mod_global_parameters
    use mod_geometry
    use mod_physics, only: phys_get_Rfactor
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw),wCTprim(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in) :: qsourcesplit
    logical, intent(inout) :: active

    integer :: idir,jdir,nth_for_fld,ix^D
    double precision :: sigma_b,cc
    double precision :: grad_E(ixI^S)
    double precision, dimension(ixO^S) :: a1,a2,a3,c0,c1,kappa
    double precision, dimension(ixI^S) :: e_gas,E_rad,tmp
    double precision, dimension(ixI^S,1:ndim,1:ndim) :: div_v,edd

    !> Calculate and add sourceterms
    if(qsourcesplit .eqv. fld_Radforce_split) then
      active = .true.
      nth_for_fld=2
      ! store here lambda in a1 and fld_R in a2
      call fld_get_eddington(wCTprim,x,ixI^L,ixO^L,edd,a1,a2,nth_for_fld)
      do idir = 1,ndim
        call gradient(wCTprim(ixI^S,iw_r_e),ixI^L,ixO^L,idir,grad_E,nth_for_fld)
        ! Radiation force = kappa*rho/c *Flux = lambda gradE
        ! recycle grad E to store -lambda (grad E)_i
        grad_E(ixO^S) = -a1(ixO^S)*grad_E(ixO^S)
        !> Momentum equation source term
        w(ixO^S,iw_mom(idir)) = w(ixO^S,iw_mom(idir))+ qdt*grad_E(ixO^S)
        !> Energy equation source term 
        w(ixO^S,iw_e) = w(ixO^S,iw_e) + qdt*wCTprim(ixO^S,iw_mom(idir))*grad_E(ixO^S)
      enddo

      !> Photon tiring : calculate tensor grad v (named div_v here)
      do idir = 1,ndim
        do jdir = 1,ndim
          call gradient(wCTprim(ixI^S,iw_mom(jdir)),ixI^L,ixO^L,idir,tmp)
          div_v(ixO^S,idir,jdir) = tmp(ixO^S)
        enddo
      enddo
      ! perform contraction fe : grad(v) with fe eddington tensor
      {^IFONED
      a3(ixO^S) = div_v(ixO^S,1,1)*edd(ixO^S,1,1)
      }
      {^IFTWOD
      a3(ixO^S) = div_v(ixO^S,1,1)*edd(ixO^S,1,1) &
                      + div_v(ixO^S,1,2)*edd(ixO^S,1,2) &
                      + div_v(ixO^S,2,1)*edd(ixO^S,2,1) &
                      + div_v(ixO^S,2,2)*edd(ixO^S,2,2)
      }
      {^IFTHREED
      a3(ixO^S) = div_v(ixO^S,1,1)*edd(ixO^S,1,1) &
                      + div_v(ixO^S,1,2)*edd(ixO^S,1,2) &
                      + div_v(ixO^S,1,3)*edd(ixO^S,1,3) &
                      + div_v(ixO^S,2,1)*edd(ixO^S,2,1) &
                      + div_v(ixO^S,2,2)*edd(ixO^S,2,2) &
                      + div_v(ixO^S,2,3)*edd(ixO^S,2,3) &
                      + div_v(ixO^S,3,1)*edd(ixO^S,3,1) &
                      + div_v(ixO^S,3,2)*edd(ixO^S,3,2) &
                      + div_v(ixO^S,3,3)*edd(ixO^S,3,3)
      }
      call fld_get_opacity_prim(wCTprim,x,ixI^L,ixO^L,kappa)

      cc = const_c/unit_velocity
      sigma_b = const_rad_a*cc/4.d0*(unit_temperature**4.d0)/(unit_pressure)
      !> e_gas is the INTERNAL ENERGY without KINETIC ENERGY
      e_gas(ixO^S) = w(ixO^S,iw_e)-half*sum(w(ixO^S,iw_mom(:))**2,dim=ndim+1)/w(ixO^S,iw_rho)
      if(allocated(iw_mag)) then
        e_gas(ixO^S) = e_gas(ixO^S)-half*sum(w(ixO^S,iw_mag(:))**2,dim=ndim+1)
      endif
      {do ix^D = ixOmin^D,ixOmax^D\ }
        e_gas(ix^D) = max(e_gas(ix^D),small_pressure/(fld_gamma-1.d0))
      {enddo\}
      E_rad(ixO^S) = w(ixO^S,iw_r_e)
      call phys_get_Rfactor(wCT,x,ixI^L,ixO^L,tmp)
      !> Coefficients for the polynomial in Moens et al. 2022, eq 37. but with photon tiring (a3)
      a1(ixO^S) = 4.d0*sigma_b*qdt*kappa(ixO^S)*(fld_gamma-one)**4/(wCT(ixO^S,iw_rho)**3*tmp(ixO^S)**4)
      a2(ixO^S) = cc*kappa(ixO^S)*wCT(ixO^S,iw_rho)*qdt
      a3(ixO^S) = a3(ixO^S)*qdt

      c0(ixO^S) = ((one+a2(ixO^S)+a3(ixO^S))*e_gas(ixO^S)+a2(ixO^S)*E_rad(ixO^S))/a1(ixO^S)/(one+a3(ixO^S))
      c1(ixO^S) = (one+a2(ixO^S)+a3(ixO^S))/a1(ixO^S)/(one+a3(ixO^S))

      !> Loop over every cell for rootfinding method
      {do ix^D = ixOmin^D,ixOmax^D\}
        select case(fld_interaction_method)
        case('Bisect')
        call Bisection_method(e_gas(ix^D),E_rad(ix^D),c0(ix^D),c1(ix^D))
        case('Newton')
        call Newton_method(e_gas(ix^D),E_rad(ix^D),c0(ix^D),c1(ix^D))
        case('Halley')
        call Halley_method(e_gas(ix^D),E_rad(ix^D),c0(ix^D),c1(ix^D))
        case default
        call mpistop('root-method not known')
        end select 
      {enddo\}

      E_rad(ixO^S) = (a1(ixO^S)*e_gas(ixO^S)**4.d0+E_rad(ixO^S))/(one+a2(ixO^S)+a3(ixO^S))
      !> Update gas-energy in w, internal + kinetic
      w(ixO^S,iw_e) = e_gas(ixO^S)
      w(ixO^S,iw_e) = w(ixO^S,iw_e)+half*sum(w(ixO^S,iw_mom(:))**2,dim=ndim+1)/w(ixO^S,iw_rho)
      if(allocated(iw_mag)) then
        w(ixO^S,iw_e) = w(ixO^S,iw_e)+half*sum(w(ixO^S,iw_mag(:))**2,dim=ndim+1)
      endif
      !> Update rad-energy in w
      w(ixO^S,iw_r_e) = E_rad(ixO^S)
    end if
  end subroutine add_fld_rad_force

  !> get dt limit for radiation force and FLD explicit source additions
  !> NOTE: w is primitive on entry
  subroutine fld_radforce_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry
    use mod_physics, only: phys_get_csrad2

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S,1:ndim), w(ixI^S,1:nw)
    double precision, intent(inout) :: dtnew

    integer          :: idim,idims,nth_for_fld
    double precision :: dxinv(1:ndim), max_grav
    double precision :: lambda(ixO^S),fld_R(ixO^S)
    double precision :: tmp(ixI^S)
    double precision :: cmax(ixI^S),cmaxtot(ixI^S),courantmaxtots

    nth_for_fld=2
    call fld_get_fluxlimiter_prim(w,x,ixI^L,ixO^L,lambda,fld_R,nth_for_fld)
    if(slab_uniform) then
      ^D&dxinv(^D)=one/dx^D;
      do idim = 1, ndim
        call gradient(w(ixI^S,iw_r_e),ixI^L,ixO^L,idim,tmp,nth_for_fld)
        max_grav = maxval(dabs(-lambda(ixO^S)*tmp(ixO^S)/w(ixO^S,iw_rho)))
        max_grav = max(max_grav, epsilon(1.0d0))
        dtnew = min(dtnew, 1.0d0 / dsqrt(max_grav * dxinv(idim)))
      end do
    else
      do idim = 1, ndim
        call gradient(w(ixI^S,iw_r_e),ixI^L,ixO^L,idim,tmp,nth_for_fld)
        max_grav = maxval(dabs(-lambda(ixO^S)*tmp(ixO^S)/w(ixO^S,iw_rho))/block%ds(ixO^S,idim))
        max_grav = max(max_grav, epsilon(1.0d0))
        dtnew = min(dtnew, 1.0d0 / dsqrt(max_grav))
      end do
    endif

    ! here we interface back to fld_get_radpress
    call phys_get_csrad2(w,x,ixI^L,ixO^L,tmp)
    if(slab_uniform) then
       ^D&dxinv(^D)=one/dx^D;
       do idims=1,ndim
          cmax(ixO^S)=dabs(w(ixO^S,iw_mom(idims)))+dsqrt(tmp(ixO^S))
          if(idims==1) then
            cmaxtot(ixO^S)=cmax(ixO^S)*dxinv(idims)
          else
            cmaxtot(ixO^S)=cmaxtot(ixO^S)+cmax(ixO^S)*dxinv(idims)
          end if
       end do
    else
       do idims=1,ndim
          cmax(ixO^S)=dabs(w(ixO^S,iw_mom(idims)))+dsqrt(tmp(ixO^S))
          if(idims==1) then
            cmaxtot(ixO^S)=cmax(ixO^S)/block%ds(ixO^S,idims)
          else
            cmaxtot(ixO^S)=cmaxtot(ixO^S)+cmax(ixO^S)/block%ds(ixO^S,idims)
          end if
       end do
    end if
    ! courantmaxtots='max(summed c/dx)'
    courantmaxtots=maxval(cmaxtot(ixO^S))
    if(courantmaxtots>smalldouble) dtnew=min(dtnew,courantpar/courantmaxtots)

  end subroutine fld_radforce_get_dt

  !> Sets the opacity in the w-array
  !> by calling mod_opal_opacity
  !> NOTE: assumes primitives in w
  !> NOTE: assuming opacity is local, not ok with cak line force
  subroutine fld_get_opacity_prim(w, x, ixI^L, ixO^L, fld_kappa)
    use mod_global_parameters
    use mod_physics, only: phys_get_tgas
    use mod_usr_methods
    use mod_opal_opacity
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out) :: fld_kappa(ixO^S)

    integer :: ix^D
    double precision :: rho0,Temp0,n
    double precision :: Temp(ixI^S)

    select case (fld_opacity_law)
      case('const','thomson')
        fld_kappa(ixO^S) = fld_kappa0/unit_opacity
      case('opal')
        call phys_get_tgas(w,x,ixI^L,ixO^L,Temp)
        {do ix^D=ixOmin^D,ixOmax^D\ }
          rho0 = w(ix^D,iw_rho)*unit_density
          Temp0 = Temp(ix^D)*unit_temperature
          call set_opal_opacity(rho0,Temp0,n)
          fld_kappa(ix^D) = n/unit_opacity
        {enddo\ }
      case('special')
        if (.not. associated(usr_special_opacity)) then
          call mpistop("special opacity not defined")
        endif
        call usr_special_opacity(ixI^L, ixO^L, w, x, fld_kappa)
      case default
        call mpistop("Doesn't know opacity law")
      end select
  end subroutine fld_get_opacity_prim

  !> Returns Radiation Pressure as tensor
  !> NOTE: w is primitive on entry
  subroutine fld_get_radpress(w, x, ixI^L, ixO^L, rad_pressure, nth)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L, nth
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: rad_pressure(ixO^S,1:ndim,1:ndim)
    integer :: i,j
    double precision             :: eddington_tensor(ixI^S,1:ndim,1:ndim)
    double precision             :: lambda(ixO^S),fld_R(ixO^S)

    call fld_get_eddington(w, x, ixI^L, ixO^L, eddington_tensor, lambda, fld_R, nth)
    do i=1,ndim
      do j=1,ndim
        rad_pressure(ixO^S,i,j) = eddington_tensor(ixO^S,i,j)*w(ixO^S,iw_r_e)
      enddo
    enddo
  end subroutine fld_get_radpress

  !> This subroutine calculates flux limiter lambda according to fld_fluxlimiter
  !> It also calculates fld_R which is ratio of radiation scaleheight and mean free path
  !> NOTE: nth and ixI and ixO not free to choose here: TODO
  subroutine fld_get_fluxlimiter(w,x,ixI^L,ixO^L,fld_lambda,fld_R,nth)
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods
    use mod_physics, only: phys_to_primitive
    integer, intent(in)           :: ixI^L,ixO^L,nth
    double precision, intent(in)  :: w(ixI^S,1:nw)
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(out) :: fld_R(ixO^S),fld_lambda(ixO^S)

    double precision :: wprim(ixI^S,1:nw)

    wprim=w
    call phys_to_primitive(ixI^L,ixI^L,wprim,x)
    call fld_get_fluxlimiter_prim(wprim,x,ixI^L,ixO^L,fld_lambda,fld_R,nth)

  end subroutine fld_get_fluxlimiter

  !> This subroutine calculates flux limiter lambda according to fld_fluxlimiter
  !> It also calculates fld_R which is ratio of radiation scaleheight and mean free path
  !> NOTE: this one operates on primitives
  !> NOTE: nth and ixI and ixO not free to choose 
  subroutine fld_get_fluxlimiter_prim(w,x,ixI^L,ixO^L,fld_lambda,fld_R,nth)
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods
    integer, intent(in)           :: ixI^L,ixO^L,nth
    double precision, intent(in)  :: w(ixI^S,1:nw)
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(out) :: fld_R(ixO^S),fld_lambda(ixO^S)

    integer :: idir, ix^D
    double precision :: kappa(ixO^S),normgrad2(ixO^S)
    double precision :: grad_r_e(ixI^S)

    select case(fld_fluxlimiter)
    case('Diffusion')
      ! optically thick limit
      fld_lambda(ixO^S) = 1.d0/3.d0
      fld_R(ixO^S) = zero
    case('FreeStream')
      ! optically thin limit
      normgrad2(ixO^S) = zero
      do idir=1,ndim
        call gradient(w(ixI^S,iw_r_e),ixI^L,ixO^L,idir,grad_r_e,nth)
        normgrad2(ixO^S) = normgrad2(ixO^S)+grad_r_e(ixO^S)**2
      end do
      call fld_get_opacity_prim(w,x,ixI^L,ixO^L,kappa)
      ! Calculate R everywhere
      ! |grad E|/(rho kappa E)
      fld_R(ixO^S) = dsqrt(normgrad2(ixO^S))/(kappa(ixO^S)*w(ixO^S,iw_rho)*w(ixO^S,iw_r_e))
      where(normgrad2(ixO^S)<smalldouble**2)
        ! treat uniform case as diffusion limit
        fld_R(ixO^S)=zero
        fld_lambda(ixO^S) = 1.0d0/3.0d0
      elsewhere
        fld_lambda(ixO^S) = one/fld_R(ixO^S)
      endwhere
    case('Pomraning')
      ! Calculate R everywhere
      ! |grad E|/(rho kappa E)
      normgrad2(ixO^S) = zero
      do idir = 1,ndim
        call gradient(w(ixI^S,iw_r_e),ixI^L,ixO^L,idir,grad_r_e,nth)
        normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_e(ixO^S)**2
      end do
      call fld_get_opacity_prim(w,x,ixI^L,ixO^L,kappa)
      fld_R(ixO^S) = dsqrt(normgrad2(ixO^S))/(kappa(ixO^S)*w(ixO^S,iw_rho)*w(ixO^S,iw_r_e))
      ! Calculate the flux limiter, lambda
      ! Levermore and Pomraning: lambda = (2 + R)/(6 + 3R + R^2)
      fld_lambda(ixO^S) = (2.d0+fld_R(ixO^S))/(6.d0+3*fld_R(ixO^S)+fld_R(ixO^S)**2)
    case('Minerbo')
      ! Calculate R everywhere
      ! |grad E|/(rho kappa E)
      normgrad2(ixO^S) = zero
      do idir = 1,ndim
        call gradient(w(ixI^S,iw_r_e),ixI^L,ixO^L,idir,grad_r_e,nth)
        normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_e(ixO^S)**2
      end do
      call fld_get_opacity_prim(w, x, ixI^L, ixO^L, kappa)
      fld_R(ixO^S) = dsqrt(normgrad2(ixO^S))/(kappa(ixO^S)*w(ixO^S,iw_rho)*w(ixO^S,iw_r_e))
      ! Calculate the flux limiter, lambda
      ! Minerbo:
      {do ix^D = ixOmin^D,ixOmax^D\ }
        if(fld_R(ix^D) .lt. 3.d0/2.d0) then
          fld_lambda(ix^D) = 2.d0/(3.d0+dsqrt(9.d0+12.d0*fld_R(ix^D)**2))
        else
          fld_lambda(ix^D) = 1.d0/(1.d0+fld_R(ix^D)+dsqrt(1.d0+2.d0*fld_R(ix^D)))
        endif
      {enddo\}
    case('special')
      if (.not. associated(usr_special_fluxlimiter)) then
        call mpistop("special fluxlimiter not defined")
      endif
      call usr_special_fluxlimiter(ixI^L,ixO^L,w,x,fld_lambda,fld_R)
    case default
      call mpistop('Fluxlimiter unknown')
    end select

  end subroutine fld_get_fluxlimiter_prim

  !> Calculate Radiation Flux (use of cgs units)
  !> NOTE: only for diagnostics purposes (w conservative on entry)
  !> This returns cell centered values for radiation flux 
  subroutine fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux)
    use mod_global_parameters
    use mod_geometry
    use mod_physics, only: phys_to_primitive
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out) :: rad_flux(ixO^S, 1:ndim)

    integer :: idir,nth_for_fld
    double precision :: cc
    double precision :: wprim(ixI^S,1:nw)
    double precision :: grad_r_e(ixI^S)
    double precision :: kappa(ixO^S), lambda(ixO^S), fld_R(ixO^S)

    wprim=w
    call phys_to_primitive(ixI^L,ixI^L,wprim,x)

    cc = const_c/unit_velocity
    ! NOTE: assuming opacity is local, not ok with cak line force
    call fld_get_opacity_prim(wprim, x, ixI^L, ixO^L, kappa)
    ! CHECK ranges
    nth_for_fld=2
    call fld_get_fluxlimiter_prim(wprim, x, ixI^L, ixO^L, lambda, fld_R, nth_for_fld)
    !> Calculate the Flux using the fld closure relation
    !> F = -c*lambda/(kappa*rho) *grad E
    do idir = 1,ndim
      call gradient(wprim(ixI^S,iw_r_e),ixI^L,ixO^L,idir,grad_r_e,nth_for_fld)
      rad_flux(ixO^S, idir) = -cc*lambda(ixO^S)/(kappa(ixO^S)*wprim(ixO^S,iw_rho))*grad_r_e(ixO^S)
    end do
  end subroutine fld_get_radflux

  !> Calculate Eddington-tensor (where w is primitive)
  !> also feeds back the flux limiter lambda and R
  subroutine fld_get_eddington(w, x, ixI^L, ixO^L, eddington_tensor, lambda, fld_R, nth)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)          :: ixI^L, ixO^L, nth
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out) :: eddington_tensor(ixI^S,1:ndim,1:ndim)
    double precision, intent(out) :: lambda(ixO^S),fld_R(ixO^S)

    integer :: idir,jdir
    double precision :: normgrad2(ixO^S), f(ixO^S)
    double precision :: grad_r_e(ixI^S,1:ndim)

    normgrad2(ixO^S) = zero
    do idir = 1,ndim
      call gradient(w(ixI^S, iw_r_e),ixI^L,ixO^L,idir,grad_r_e(ixI^S,idir),nth)
      normgrad2(ixO^S)=normgrad2(ixO^S)+grad_r_e(ixO^S,idir)**2
    end do
    ! get lambda and R  
    call fld_get_fluxlimiter_prim(w,x,ixI^L,ixO^L,lambda,fld_R,nth)
    ! store f_e= lambda + lambda^2 R^2
    f(ixO^S) = lambda(ixO^S)+(lambda(ixO^S)*fld_R(ixO^S))**2
    do idir = 1,ndim
      ! first compute the isotropic (diagonal) part
      eddington_tensor(ixO^S,idir,idir) = half*(one-f(ixO^S))
    enddo
    do idir = 1,ndim
      do jdir = 1,ndim
        ! initialize off-diagonal part here
        if(idir .ne. jdir) eddington_tensor(ixO^S,idir,jdir) = zero
        ! add part depending on unit vectors along gradient E
        where(normgrad2(ixO^S) .gt. smalldouble**2)
          eddington_tensor(ixO^S,idir,jdir) = eddington_tensor(ixO^S,idir,jdir)+&
            half*(3.d0*f(ixO^S)-one)*grad_r_e(ixO^S,idir)*grad_r_e(ixO^S,jdir)/normgrad2(ixO^S)
        endwhere
      enddo
    enddo
  end subroutine fld_get_eddington

  !> Calling all subroutines to perform the multigrid method
  !> Communicates rad_e and diff_coeff to multigrid library
  !> Advance psa=psb+dtfactor*qdt*F_im(psa)
  subroutine fld_implicit_update(dtfactor,qdt,qtC,psa,psb)
    use mod_global_parameters
    use mod_forest
    use mod_multigrid_coupling
    type(state), target :: psa(max_blocks)
    type(state), target :: psb(max_blocks)
    double precision, intent(in) :: qdt
    double precision, intent(in) :: qtC
    double precision, intent(in) :: dtfactor

    integer                      :: n
    double precision             :: res, max_residual, lambda
    integer, parameter           :: max_its = 100
    integer                      :: iw_to,iw_from
    integer                      :: iigrid, igrid, id
    integer                      :: nc, lvl, i
    type(tree_node), pointer     :: pnode
    double precision             :: fac, facD

    !> Set diffusion timestep, add previous timestep if mg did not converge:
    if(it == 0) dt_diff = 0
    dt_diff = dt_diff + qdt
    ! Avoid setting a very restrictive limit to the residual when the time step
    ! is small (as the operator is ~ 1/(D * qdt))
    if(qdt < dtmin) then
      if(mype==0)then
          print *,'skipping implicit solve: dt too small!'
          print *,'Currently at time=',global_time,' time step=',qdt,' dtmin=',dtmin
      endif
      return
    endif
    max_residual = fld_diff_tol
    mg%operator_type = mg_vhelmholtz
    mg%smoother_type = mg_smoother_gs
    call mg_set_methods(mg)
    if(.not. mg%is_allocated) call mpistop("multigrid tree not allocated yet")
    lambda = 1.d0/(dtfactor * dt_diff)
    call vhelmholtz_set_lambda(lambda)
    call update_diffcoeff(psb)
    fac = 1.d0
    facD = 1.d0
    call mg_copy_to_tree(i_diff_mg, mg_iveps, factor=facD, state_from=psb)
    call mg_copy_to_tree(iw_r_e, mg_iphi, factor=fac, state_from=psb)
    call mg_copy_to_tree(iw_r_e, mg_irhs, factor=-1/(dtfactor*dt_diff), state_from=psb)
    if(time_advance)then
      call mg_restrict(mg, mg_iveps)
      call mg_fill_ghost_cells(mg, mg_iveps)
    endif
    call mg_fas_fmg(mg, .true., max_res=res)
    do n = 1, max_its
      if(res < max_residual) exit
      call mg_fas_vcycle(mg, max_res=res)
    end do
    if(res .le. 0.d0) then
      if (diffcrash_resume) then
        if (mg%my_rank == 0) &
        write(*,*) it, ' residual zero ', res
        return
      endif
      if(mg%my_rank == 0) then
        print*, res
        error stop "Diffusion residual to zero"
      endif
    endif
!    if(n==max_its+1) then
!       if(diffcrash_resume) then
!         if(mg%my_rank == 0) &
!         write(*,*) it, ' residual high ', res
!         return
!       endif
!       if(mg%my_rank == 0) then
!          print *, "Did you specify boundary conditions correctly?"
!          print *, "Or is the variation in diffusion too large?"
!          print*, n, res
!          print *, mg%bc(1, mg_iphi)%bc_value, mg%bc(2, mg_iphi)%bc_value
!       end if
!       error stop "diffusion_solve: no convergence"
!    end if
    !> Reset dt_diff when diffusion worked out
    dt_diff = 0.d0
    call mg_copy_from_tree_gc(mg_iphi, iw_r_e, state_to=psa)
  end subroutine fld_implicit_update

  !> inplace update of psa==>F_im(psa)
  subroutine fld_evaluate_implicit(qtC,psa)
    use mod_global_parameters
    use mod_ghostcells_update
    type(state), target :: psa(max_blocks)
    double precision, intent(in) :: qtC
    integer :: iigrid, igrid, level
    integer :: ixO^L

    call update_diffcoeff(psa)
    ixO^L=ixM^LL^LADD1;
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
       call put_diffterm_onegrid(ixG^LL,ixO^L,psa(igrid)%w)
    end do
    !$OMP END PARALLEL DO
    ! enforce boundary conditions for psa
    call getbc(qtC,0.d0,psa,1,nwflux+nwaux)
  end subroutine fld_evaluate_implicit

  !> inplace update of psa==>F_im(psa)
  subroutine put_diffterm_onegrid(ixI^L,ixO^L,w)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision :: gradE(ixO^S),divF(ixO^S)
    double precision :: divF_h(ixO^S),divF_j(ixO^S)
    double precision :: diff_term(ixO^S)
    integer :: idir, jxO^L, hxO^L

    divF(ixO^S) = 0.d0
    do idir = 1,ndim
      hxO^L=ixO^L-kr(idir,^D);
      jxO^L=ixO^L+kr(idir,^D);
      divF_h(ixO^S) = w(ixO^S,i_diff_mg)*w(hxO^S,i_diff_mg)/(w(ixO^S,i_diff_mg) + w(hxO^S,i_diff_mg))*(w(hxO^S,iw_r_e) - w(ixO^S,iw_r_e))
      divF_j(ixO^S) = w(ixO^S,i_diff_mg)*w(jxO^S,i_diff_mg)/(w(ixO^S,i_diff_mg) + w(jxO^S,i_diff_mg))*(w(jxO^S,iw_r_e) - w(ixO^S,iw_r_e))
      divF(ixO^S) = divF(ixO^S) + 2.d0*(divF_h(ixO^S) + divF_j(ixO^S))/dxlevel(idir)**2
    enddo
    w(ixO^S,iw_r_e) = divF(ixO^S)
  end subroutine put_diffterm_onegrid

  !> Calculates cell-centered diffusion coefficient to be used in multigrid
  subroutine fld_get_diffcoef_central(w, x, ixI^L, ixO^L)
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods
    use mod_physics, only: phys_to_primitive
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)

    integer :: idir,i,j,ix^D
    double precision :: wprim(ixI^S,1:nw)
    double precision :: cc
    double precision :: kappa(ixO^S),lambda(ixO^S),fld_R(ixO^S)
    double precision :: max_D(ixI^S),grad_r_e(ixI^S),rad_e(ixI^S)

    cc = const_c/unit_velocity
    wprim=w
    call phys_to_primitive(ixI^L,ixI^L,wprim,x)
    call fld_get_opacity_prim(wprim, x, ixI^L, ixO^L, kappa)
    call fld_get_fluxlimiter_prim(wprim, x, ixI^L, ixO^L, lambda, fld_R, nghostcells-1)
    w(ixO^S,i_diff_mg) = cc*lambda(ixO^S)/(kappa(ixO^S)*wprim(ixO^S,iw_rho))
    where(w(ixO^S,i_diff_mg) .lt. 0.d0) w(ixO^S,i_diff_mg) = smalldouble
    if(associated(usr_special_diffcoef)) call usr_special_diffcoef(w, wprim, x, ixI^L, ixO^L)

  end subroutine fld_get_diffcoef_central

  subroutine update_diffcoeff(psa)
    use mod_global_parameters
    type(state), target :: psa(max_blocks)
    integer :: iigrid, igrid, level
    integer :: ixO^L

    ixO^L=ixM^LL^LADD1;
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
        call fld_get_diffcoef_central(psa(igrid)%w, psa(igrid)%x, ixG^LL, ixO^L)
    end do
    !$OMP END PARALLEL DO
  end subroutine update_diffcoeff

  !> Find the root of the 4th degree polynomial using the bisection method
  subroutine Bisection_method(e_gas, E_rad, c0, c1)
    use mod_global_parameters
    double precision, intent(in)    :: c0, c1
    double precision, intent(in)    :: E_rad
    double precision, intent(inout) :: e_gas
    double precision :: bisect_a, bisect_b, bisect_c
    integer :: n, max_its

    n = 0
    max_its = 100
    bisect_a = zero
    bisect_b = min(abs(c0/c1),abs(c0)**(1.d0/4.d0))+smalldouble
    do while(dabs(bisect_b-bisect_a) .ge. fld_bisect_tol*min(e_gas,E_rad))
      bisect_c = (bisect_a + bisect_b)/two
      n = n +1
      if(n .gt. max_its) then
        call mpistop('No convergece in bisection scheme')
      endif
      if(Polynomial_Bisection(bisect_a, c0, c1)*&
      Polynomial_Bisection(bisect_b, c0, c1) .lt. zero) then
        if(Polynomial_Bisection(bisect_a, c0, c1)*&
        Polynomial_Bisection(bisect_c, c0, c1) .lt. zero) then
          bisect_b = bisect_c
        elseif(Polynomial_Bisection(bisect_b, c0, c1)*&
        Polynomial_Bisection(bisect_c, c0, c1) .lt. zero) then
          bisect_a = bisect_c
        elseif(Polynomial_Bisection(bisect_a, c0, c1) .eq. zero) then
          bisect_b = bisect_a
          bisect_c = bisect_a
          goto 2435
        elseif(Polynomial_Bisection(bisect_b, c0, c1) .eq. zero) then
          bisect_a = bisect_b
          bisect_c = bisect_b
          goto 2435
        elseif(Polynomial_Bisection(bisect_c, c0, c1) .eq. zero) then
          bisect_a = bisect_c
          bisect_b = bisect_c
          goto 2435
        else
          call mpistop("Problem with fld bisection method")
        endif
      elseif(Polynomial_Bisection(bisect_a, c0, c1) &
        - Polynomial_Bisection(bisect_b, c0, c1) .lt. fld_bisect_tol*Polynomial_Bisection(bisect_a, c0, c1)) then
        goto 2435
      else
        bisect_a = e_gas
        bisect_b = e_gas
        print*, "IGNORING GAS-RAD ENERGY EXCHANGE ", c0, c1
        print*, Polynomial_Bisection(bisect_a, c0, c1), Polynomial_Bisection(bisect_b, c0, c1)
        call mpistop('issues in bisection scheme')
        if(Polynomial_Bisection(bisect_a, c0, c1) .le. smalldouble) then
          bisect_b = bisect_a
        elseif(Polynomial_Bisection(bisect_a, c0, c1) .le. smalldouble) then
          bisect_a = bisect_b
        endif
        goto 2435
      endif
    enddo
      2435 e_gas = (bisect_a + bisect_b)/two
  end subroutine Bisection_method

  !> Find the root of the 4th degree polynomial using the Newton method
  subroutine Newton_method(e_gas, E_rad, c0, c1)
    use mod_global_parameters
    double precision, intent(in)    :: c0, c1
    double precision, intent(in)    :: E_rad
    double precision, intent(inout) :: e_gas
    double precision :: xval, yval, der, deltax
    integer :: ii

    yval = bigdouble
    xval = e_gas
    der = one
    deltax = one
    ii = 0
    !> Compare error with dx = dx/dy dy
    do while(dabs(deltax) .gt. fld_bisect_tol)
      yval = Polynomial_Bisection(xval, c0, c1)
      der  = dPolynomial_Bisection(xval, c0, c1)
      deltax = -yval/der
      xval = xval + deltax
      ii = ii + 1
      if(ii .gt. 1d3) then
        print*, 'skip to bisection algorithm'
        call Bisection_method(e_gas, E_rad, c0, c1)
        return
      endif
    enddo
    e_gas = xval
  end subroutine Newton_method

  !> Find the root of the 4th degree polynomial using the Halley method
  subroutine Halley_method(e_gas, E_rad, c0, c1)
    use mod_global_parameters
    double precision, intent(in)    :: c0, c1
    double precision, intent(in)    :: E_rad
    double precision, intent(inout) :: e_gas
    double precision :: xval, yval, der, dder, deltax
    integer :: ii

    yval = bigdouble
    xval = e_gas
    der = one
    dder = one
    deltax = one
    ii = 0
    !> Compare error with dx = dx/dy dy
    do while (dabs(deltax) .gt. fld_bisect_tol)
      yval = Polynomial_Bisection(xval, c0, c1)
      der  = dPolynomial_Bisection(xval, c0, c1)
      dder = ddPolynomial_Bisection(xval, c0, c1)
      deltax = -two*yval*der/(two*der**2 - yval*dder)
      xval = xval + deltax
      ii = ii + 1
      if(ii .gt. 1d3) then
        print*, 'skip to Newton algorithm'
        !call mpistop("iterates to 1000, entering Newton method")
        call Newton_method(e_gas, E_rad, c0, c1)
        return
      endif
    enddo
    e_gas = xval
  end subroutine Halley_method

  !> Evaluate polynomial at argument e_gas
  function Polynomial_Bisection(e_gas, c0, c1) result(val)
    use mod_global_parameters
    double precision, intent(in) :: e_gas
    double precision, intent(in) :: c0, c1
    double precision :: val

    val = e_gas**4.d0 + c1*e_gas - c0
  end function Polynomial_Bisection

  !> Evaluate first derivative of polynomial at argument e_gas
  function dPolynomial_Bisection(e_gas, c0, c1) result(der)
    use mod_global_parameters
    double precision, intent(in) :: e_gas
    double precision, intent(in) :: c0, c1
    double precision :: der

    der = 4.d0*e_gas**3.d0 + c1
  end function dPolynomial_Bisection

  !> Evaluate second derivative of polynomial at argument e_gas
  function ddPolynomial_Bisection(e_gas, c0, c1) result(dder)
    use mod_global_parameters
    double precision, intent(in) :: e_gas
    double precision, intent(in) :: c0, c1
    double precision :: dder

    dder = 4.d0*3.d0*e_gas**2.d0
  end function ddPolynomial_Bisection
end module mod_fld
