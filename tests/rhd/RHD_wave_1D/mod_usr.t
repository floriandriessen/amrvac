module mod_usr

  ! Include a physics module
  use mod_hd
  use mod_fld

  implicit none

  ! input variables
  double precision :: rho0,eg0,Er0,wvlength,A_rho

  ! derived vars
  double precision :: p0,Trad0,wvl,omega,wavenumber,csound
  double precision :: A_v, A_p, A_e

  double precision :: rho0_norm,eg0_norm,Er0_norm,p0_norm,omega_norm

  ! Storing additional var in the dat file
  integer :: Tgas_,Trad_,pres_,vel_,amr_


contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Note how we here must set three values that in turn define M-L-T
    unit_density       = 3.216d-9    ! g cm^-3
    unit_length        = 1.24d11  ! cm
    unit_velocity      = 2.32d6   ! cm/s

    csound=dsqrt(26020.d0*2.0/3.0/3.216d-9)
    wvl=1000.0/3.216d-9/0.4
    print *,'wavelength then=',wvl
    omega = 2.d0*dpi*csound/wvl
    print *,'omega then=',omega

    call usr_params_read(par_files)

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

   ! to add selected variables to the .dat file
    usr_modify_output => set_output_vars

    ! Drive the wave using an internal boundary
    usr_internal_bc => Initialize_Wave

    ! Choose coordinate system as n-D Cartesian
    {^IFONED call set_coordinate_system("Cartesian_1D")}
    {^IFTWOD call set_coordinate_system("Cartesian_2D")}
    {^IFTHREED call set_coordinate_system("Cartesian_3D")}

    ! Active the physics module
    call hd_activate()

    ! to add selected variables to the .dat file
    Tgas_ = var_set_extravar("Tgas", "Tgas")
    Trad_ = var_set_extravar("Trad", "Trad")
    pres_ = var_set_extravar("pres", "pres")
    vel_ = var_set_extravar("vel", "vel")
    amr_ = var_set_extravar("level", "level")

  end subroutine usr_init


  !> Read parameters from a file
  subroutine usr_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /wave_list/ rho0, eg0, Er0, wvlength, A_rho

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       rewind(unitpar)
       read(unitpar, wave_list, end=113)
113    close(unitpar)
    end do

    p0 = eg0*(hd_gamma - one)
    wavenumber = 2.d0*dpi/wvlength
  if(mype==0)then
    print *,'============================================'
    write(*,*) 'INPUT GIVEN IN cgs units is'
    write(*,*) 'input density, gas energy density, radiation energy density=',rho0,eg0,Er0
    write(*,*) 'derived gas pressure is =',p0
    write(*,*) 'derived dimensionless rad/gas energy ratio =',Er0/eg0
    write(*,*) 'derived ratio r=E/(4gamma e) =',Er0/(4.d0*hd_gamma*eg0)
    csound=dsqrt(hd_gamma*p0/rho0)
    ! note: assuming cgs in line below
    write(*,*) 'derived Boltzmann ratio Bo=4gamma ca e/(cE) =',4.0d0*hd_gamma*csound*eg0/(const_c*Er0)
    write(*,*) 'input dimensionless amplitude for density is   =',A_rho
    write(*,*) 'input dimensionless wavelength is              =',wvlength
    write(*,*) 'derived dimensionless wavenumber is            =',wavenumber
    print *,'============================================'
  endif

  end subroutine usr_params_read

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    use mod_fld

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: temp(ixI^S)
    double precision :: a,b,Xfrac,Yfrac
    logical, save:: first=.true.

    rho0_norm = rho0/unit_density
    p0_norm   = p0/unit_pressure
    eg0_norm  = eg0/unit_pressure
    Er0_norm  = Er0/unit_pressure

    omega_norm=omega*unit_time
    A_v   = omega_norm/(wavenumber*rho0_norm)*A_rho
    A_p   = omega_norm**2/wavenumber**2*A_rho
    A_e   = 1.d0/(hd_gamma-one)*p0_norm/rho0_norm*A_rho

    ! Set initial values for
    w(ixI^S, rho_)  = rho0_norm
    w(ixI^S, mom(:))= 0.d0
    w(ixI^S, e_)    = eg0_norm
    w(ixI^S, r_e)   = Er0_norm

  if(mype==0.and.first)then
    print *,'===IN INITIAL CONDITIONS========================================='
    write(*,*) 'converted to normalized values'
    write(*,*) 'normalized density, pressure, energy density gas-rad =',rho0_norm,p0_norm,eg0_norm,Er0_norm
    write(*,*) 'and amplitudes become A_rho A_v A_p A_e=',A_rho,A_v,A_p,A_e
    write(*,*) 'omega and wavenumber=',omega,wavenumber
    print *,'===========================================-====================='
  print *,'========GLOBAL values==========='
  print *,'rho0=',rho0
  print *,'p0=',p0
  print *,'eg0=',eg0
  print *,'Er0=',Er0
  print *,'arad_norm=',arad_norm
  write(*,*) 'derived gas temperature is       =',unit_temperature*(p0_norm/(rho0_norm*RR)),' or dimensionless=',(p0_norm/(rho0_norm*RR))
  Trad0=(Er0/const_rad_a)**(1.0d0/4.0d0)
  write(*,*) 'derived radiation temperature is =',Trad0,' or dimensionless=',Trad0/unit_temperature
  print *,'c_norm=',c_norm
  print *,'const_kappae  =',const_kappae
  if(trim(fld_opacity_law).eq.'const_norm')then
      print *,'normalized fld_kappa0    =',fld_kappa0
      print *,'physical value           =',fld_kappa0*unit_opacity
  endif
  if(trim(fld_opacity_law).eq.'const')then
      print *,'physical fld_kappa in cgs =',fld_kappa0
      print *,'normalized value          =',fld_kappa0/unit_opacity
  endif
  print *,'gamma=',hd_gamma
  print *,'========UNITS==========='
  print *,'SI_unit       =',SI_unit
  print *,'const_rad_a   =',const_rad_a
  print *,'eq_state_units=',eq_state_units
  print *,'He_abundance  =',He_abundance
  print *,'RR            =',RR
  print *,'unit_time          =',unit_time
  print *,'unit_length        =',unit_length
  print *,'unit_velocity      =',unit_velocity
  print *,'unit_pressure      =',unit_pressure
  print *,'unit_Erad          =',unit_Erad
  print *,'unit_numberdensity =',unit_numberdensity
  print *,'unit_density       =',unit_density
  print *,'unit_mass          =',unit_mass
  print *,'unit_temperature   =',unit_temperature
  print *,'unit_radflux       =',unit_radflux
  print *, 'CHECK that ',unit_pressure,' equals ',unit_density*unit_velocity**2
  print *, 'CHECK that ',unit_length,' equals ',unit_velocity*unit_time
  print *, 'CHECK that ',unit_mass,' equals ',unit_density*unit_length**3
  print *, 'density to numberdensity has factor   ',unit_density/unit_numberdensity
  print *, '                     compare  this to ',mp_cgs*(1.d0+4.d0*He_abundance)
  print *, 'pressure to n T has factor            ',unit_pressure/(unit_numberdensity*unit_temperature)
  print *, '                     compare  this to ',kB_cgs*(2.d0+3.d0*He_abundance)
  a=unit_density/unit_numberdensity/mp_cgs
  b=unit_pressure/(unit_numberdensity*unit_temperature*kB_cgs)
  print *, 'mean molecular weight mu adopted is =',a/b,' and this equals ', (1.d0+4.d0*He_abundance)/(2.d0+3.d0*He_abundance)
  Xfrac=1.d0/a
  Yfrac=4.d0*He_abundance/(1.d0+4.d0*He_abundance)
  print *, 'mass fraction hydrogen X is =',1/a,' and this equals ', 1.d0/(1.d0+4.d0*He_abundance)
  print *, 'mass fraction helium   Y is =',Yfrac
  print *, ' check that 1/mu', b/a,' is equal to 2X+3Y/4=',2.d0*Xfrac+3.d0*Yfrac/4.d0
  print *, ' ratio n_e/n_p=',1.d0+2.0d0*He_abundance
  print *,'========UNITS==========='
    first=.false.
  endif

  end subroutine initial_conditions

  subroutine Initialize_Wave(level,qt,ixI^L,ixO^L,w,x)
    use mod_global_parameters
    use mod_fld
    integer, intent(in)             :: ixI^L,ixO^L,level
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    double precision :: temp(ixI^S),vel(ixI^S),pres(ixI^S)

    where (x(ixI^S,1) .lt. one)
      vel(ixI^S)       = A_v*dsin(wavenumber*x(ixI^S,1)-omega_norm*qt)
      w(ixI^S, rho_)   = rho0_norm + A_rho*dsin(wavenumber*x(ixI^S,1)-omega_norm*qt)
      w(ixI^S, mom(1)) = w(ixI^S, rho_)*vel(ixI^S)
      w(ixI^S, e_)     = (p0_norm + A_p*dsin(wavenumber*x(ixI^S,1)-omega_norm*qt))/(hd_gamma-1.0d0)+half*w(ixI^S,rho_)*vel(ixI^S)**2
      w(ixI^S, r_e)    = Er0_norm
    endwhere

  end subroutine Initialize_Wave

  subroutine set_output_vars(ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: Trad(ixI^S),Tgas(ixI^S),pth(ixI^S)

    call hd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    call hd_get_trad(w,x,ixI^L,ixO^L,Trad)
    call hd_get_temperature_from_etot(w,x,ixI^L,ixO^L,Tgas)
    w(ixO^S,Tgas_)=Tgas(ixO^S)
    w(ixO^S,Trad_)=Trad(ixO^S)
    w(ixO^S,pres_)=pth(ixO^S)
    w(ixO^S,vel_)=w(ixO^S,mom(1))/w(ixO^S,rho_)
    ! output the AMR level (assuming uniform grid blocks)
    w(ixO^S,amr_)=dlog(((xprobmax1-xprobmin1)/domain_nx1)/dxlevel(1))/dlog(2.0d0)+1.0d0

  end subroutine set_output_vars

end module mod_usr
