module mod_usr

  ! Include a physics module
  use mod_hd
  use mod_fld

  implicit none

  ! input in cgs
  double precision :: rho1,rho2,v1,T1

  ! derived in cgs from steady shock condition
  double precision :: v2,T2,momc,momc_n
  double precision :: v2_orig,T2_orig

  ! normalized and derived values
  double precision :: rho1_norm,rho2_norm
  double precision :: v1_norm,v2_norm
  double precision :: T1_norm,T2_norm
  double precision :: p1_norm,p2_norm
  double precision :: eg1_norm,eg2_norm,Er1_norm,Er2_norm

  ! Storing additional var in the dat file
  integer :: Tgas_,Trad_,pres_,vel_,amr_


contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Note how we here must set three values that in turn define M-L-T
    unit_density       =0.01d0       ! 0.01 g cm^-3
    unit_velocity      =1.d8         ! 10^8  cm/s
    unit_length        =1.d5         ! 10^5 cm 

    call usr_params_read(par_files)

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

   ! to add selected variables to the .dat file
    usr_modify_output => set_output_vars

    ! Manually refine grid near shock
    usr_refine_grid => refine_shock

    ! Boundary conditions needed for MG solver and for HD part
    usr_special_bc => boundary_conditions
    usr_special_mg_bc => mg_boundary_conditions

    ! Choose coordinate system as 2D Cartesian with three components for vectors
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

    namelist /shock_list/ rho1, rho2, v1, T1

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       rewind(unitpar)
       read(unitpar, shock_list, end=113)
113    close(unitpar)
    end do

  v2=rho1*v1/rho2
  momc=rho1*v1
  v2_orig=1.458051139d8
  T2_orig=4.239d7
  if(mype==0)then
    print *,'============================================'
    write(*,*) 'INPUT GIVEN IN cgs units is'
    write(*,*) 'input density 1-2     =',rho1,rho2
    write(*,*) 'input Temperature 1   =',T1
    write(*,*) 'input velocity 1    =',v1
    write(*,*) 'derived velocity 2  =',v2
    write(*,*) 'original velocity 2  =',v2_orig
    write(*,*) 'original temperature 2  =',T2_orig
    print *,'============================================'
  endif

  end subroutine usr_params_read

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: a,b,Xfrac,Yfrac,Tval
    double precision :: c1A_norm,c0A_norm
    double precision :: c1B_norm,c0B_norm
    double precision :: T2A_norm,T2B_norm,T2C_norm
    logical, save:: first=.true.

    rho1_norm = rho1/unit_density
    rho2_norm = rho2/unit_density
    T1_norm = T1/unit_temperature
    p1_norm=T1_norm*rho1_norm*RR
    Er1_norm = arad_norm*T1_norm**4.d0
    v1_norm = v1/unit_velocity
    momc_n=momc/(unit_density*unit_velocity)
    eg1_norm=p1_norm/(hd_gamma-1.0d0)+half*rho1_norm*v1_norm**2

    !!v2=v2_orig
    v2_norm = v2/unit_velocity
    ! use the Halley root finder for the 4th order polynomial in T2 (normalized)
    c1A_norm=3.0d0*rho2_norm/arad_norm
    c0A_norm=((rho2_norm-rho1_norm)*momc_n**2/(rho2_norm*rho1_norm) &
             +rho1_norm*T1_norm+Er1_norm/3.0d0)*3.0d0/arad_norm
    T2A_norm=T1_norm ! this is a bad guess
    call Halley_method(T2A_norm,T2A_norm,c0A_norm,c1A_norm)

    if(mype==0.and.first)then
       print *,'CHECK ON ZERO of',T2A_norm**4+c1A_norm*T2A_norm-c0A_norm
    endif
    if(mype==0.and.first)print *,'HERE WE GET',T2A_norm

    Tval=0.5d0*T2A_norm

    ! use the Halley root finder for the 4th order polynomial in T2 (normalized)
    c1B_norm=3.0d0*rho2_norm*hd_gamma/(arad_norm*4.0d0*(hd_gamma-1.d0))
    c0B_norm=(3.0d0/(4.0d0*arad_norm))* &
       ((hd_gamma*rho1_norm*T1_norm/(hd_gamma-1.d0) &
         +4.0d0*Er1_norm/3.0d0+half*momc_n**2/rho1_norm)*rho2_norm/rho1_norm &
       -half*momc_n**2/rho2_norm)
    T2B_norm=T1_norm ! this is a bad guess
    call Halley_method(T2B_norm,T2B_norm,c0B_norm,c1B_norm)
    if(mype==0.and.first)then
       print *,'CHECK ON ZERO of',T2B_norm**4+c1B_norm*T2B_norm-c0B_norm
    endif

    if(mype==0.and.first)print *,'HERE WE use INSTEAD',T2B_norm
    Tval=Tval+0.5d0*T2B_norm
    if(mype==0.and.first)print *,'average is',Tval
    T2=T2_orig
    T2_norm=T2_orig/unit_temperature
    if(mype==0.and.first)print *,'WHILE we had original value',T2_norm

    if(mype==0.and.first)then
       print *,'difference in T is',T2A_norm-T2B_norm
       print *,'makes the factor ',(c1A_norm*T2A_norm-c0A_norm-c1B_norm*T2B_norm+c0B_norm)
       print *,'equal to',T2A_norm**4-T2B_norm**4
       print *,'c0A_norm, c0B_norm, difference=',c0A_norm,c0B_norm,c0A_norm-c0B_norm
       print *,'c1A_norm, c1B_norm, difference=',c1A_norm,c1B_norm,c1A_norm-c1B_norm
       print *,'versus ratio=',(c0A_norm-c0B_norm)/(c1A_norm-c1B_norm)
    endif
    T2_norm=T2B_norm
    if(mype==0.and.first)print *,'USING',T2_norm
     
    T2=T2_norm*unit_temperature
    p2_norm=T2_norm*rho2_norm*RR
    eg2_norm=p2_norm/(hd_gamma-1.0d0)+half*rho2_norm*v2_norm**2
    Er2_norm = arad_norm*T2_norm**4.d0

    w(ixI^S,rho_)   = rho1_norm
    w(ixI^S,mom(:)) = 0.d0
    w(ixI^S,mom(1)) = rho1_norm*v1_norm
    w(ixI^S,e_)     = eg1_norm
    w(ixI^S,r_e)    = Er1_norm

    where (x(ixI^S,1) .gt. 0.d0)
      w(ixI^S,rho_)   = rho2_norm
      w(ixI^S,mom(1)) = rho2_norm*v2_norm
      w(ixI^S,e_)     = eg2_norm
      w(ixI^S,r_e)    = Er2_norm
    end where

  if(mype==0.and.first)then
    print *,'===IN INITIAL CONDITIONS========================================='
    write(*,*) 'converted to normalized values'
    write(*,*) 'normalized densities          =',rho1_norm,rho2_norm
    write(*,*) 'normalized velocities         =',v1_norm,v2_norm
    write(*,*) 'normalized temperatures       =',T1_norm,T2_norm
    write(*,*) 'normalized gas energies       =',eg1_norm,eg2_norm
    write(*,*) 'normalized radiation energies =',Er1_norm,Er2_norm
    print *,'================================================================='
  print *,'========GLOBAL values==========='
  print *,'rho 1-2=',rho1,rho2
  print *,'T   1-2=',T1,T2
  print *,'v   1-2=',v1,v2
    print *,'================================================================='
    print*, 'RHD-fluxes for stationary shock: ', 'Left', ' | ', 'Right'
    print*, 'density flux must be exactly equal  :', rho1_norm*v1_norm, ' | ', rho2_norm*v2_norm
    print*, 'momentum flux must be exact   equal :',rho1_norm*v1_norm**2+p1_norm+Er1_norm/3,' | ',rho2_norm*v2_norm**2+p2_norm+Er2_norm/3
    print*, 'energy flux must be equal           :', (p1_norm+eg1_norm+4.0d0*Er1_norm/3.0d0)*v1_norm, ' | ', (p2_norm+eg2_norm+4.0d0*Er2_norm/3.0d0)*v2_norm
    print*, 'with radiation and gas temperatures equal on each side: LEFT  is', T1_norm,(Er1_norm/arad_norm)**0.25d0
    print*, 'with radiation and gas temperatures equal on each side: RIGHT is', T2_norm,(Er2_norm/arad_norm)**0.25d0
    print *,'================================================================='
  print *,'arad_norm=',arad_norm
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
  print *, 'CHECK that p_u ',unit_pressure,' equals ',unit_density*unit_velocity**2
  print *, 'CHECK that L_u ',unit_length,' equals ',unit_velocity*unit_time
  print *, 'CHECK that M_u',unit_mass,' equals ',unit_density*unit_length**3
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


  subroutine boundary_conditions(qt,ixI^L,ixB^L,iB,w,x)
    use mod_global_parameters
    use mod_fld

    integer, intent(in)             :: ixI^L, ixB^L, iB
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    integer :: ii

    select case (iB)

    case(1)
      w(ixB^S,rho_) = rho1_norm
      w(ixB^S,mom(:)) = 0.d0
      w(ixB^S,mom(1)) = rho1_norm*v1_norm
      w(ixB^S,e_) = eg1_norm
      !! w(ixB^S,r_e) = Er1_norm
    case(2)  
      w(ixB^S,rho_) = rho2_norm
      w(ixB^S,mom(:)) = 0.d0
      w(ixB^S,mom(1)) = rho2_norm*v2_norm
      w(ixB^S,e_) = eg2_norm
      !!w(ixB^S,r_e) = Er2_norm
    case default
      call mpistop('boundary not known')
    end select
  end subroutine boundary_conditions


  subroutine mg_boundary_conditions(iB)

    use mod_global_parameters
    use mod_multigrid_coupling

    integer, intent(in)             :: iB

    select case (iB)
    case (1)
        mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
        mg%bc(iB, mg_iphi)%bc_value = Er1_norm
    case (2)
        mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
        mg%bc(iB, mg_iphi)%bc_value = Er2_norm
    case default
        call mpistop("no such boundary for MG in mod_usr")
    end select
  end subroutine mg_boundary_conditions

  subroutine refine_shock(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)
    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    ! you must set consistent values for integers refine/coarsen:
    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement
    ! coarsen = -1 enforce to not coarsen
    ! coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixG^L, ix^L
    double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    if (any(dabs(x(ixG^S,1)) < 0.1d0)) refine=1

  end subroutine refine_shock

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
