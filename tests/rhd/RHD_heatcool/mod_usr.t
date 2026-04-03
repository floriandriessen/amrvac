module mod_usr

  ! Include a physics module
  use mod_hd
  use mod_fld

  implicit none

  double precision :: e_eq
  double precision :: rho0
  double precision :: t0
  double precision :: e0
  double precision :: E_r0
  double precision :: fld_mu
  ! Storing additional var in the dat file
  integer :: delta_e_, amr_


contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods
    use mod_constants

    call set_coordinate_system("Cartesian_2D")

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

   ! to add selected variables to the .dat file
    usr_modify_output => set_output_vars

    ! Active the physics module
    call hd_activate()

    ! to add selected variables to the .dat file
    delta_e_ = var_set_extravar("delta_e", "delta_e")
    !amr_ = var_set_extravar("amr", "amr")

  end subroutine usr_init

subroutine initglobaldata_usr
  use mod_global_parameters

  !Units
  unit_numberdensity = 1.d0/((1.d0+4.d0*He_abundance)*const_mp)
  unit_pressure = 1.d0
  unit_time = 1.d0

  unit_density=(1.d0+4.d0*He_abundance)*const_mp*unit_numberdensity
  unit_velocity = dsqrt(unit_pressure/unit_density)
  unit_length = unit_time*unit_velocity
  unit_temperature=unit_pressure/((2.d0+3.d0*He_abundance)*unit_numberdensity*const_kB)

  unit_radflux = unit_velocity*unit_pressure
  unit_opacity = one/(unit_density*unit_length)

  call usr_params_read(par_files)
  fld_mu=(1.0d0+4.0d0*He_abundance)/(2.0d0+3.0d0*He_abundance)

  e_eq = (E_r0/(const_rad_a))**(1.d0/4.d0) &
        *one/(hd_gamma-one)*const_kB*rho0 &
        /(fld_mu*const_mp)

end subroutine initglobaldata_usr


subroutine usr_params_read(files)
  use mod_global_parameters, only: unitpar
  character(len=*), intent(in) :: files(:)
  integer                      :: n

  namelist /test_list/ rho0, E_r0, t0

  do n = 1, size(files)
     open(unitpar, file=trim(files(n)), status="old")
     read(unitpar, test_list, end=111)
     111    close(unitpar)
  end do


end subroutine usr_params_read

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    use mod_constants

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S, ndim)
    double precision, intent(inout) :: w(ixI^S, nw)
    logical, save:: first=.true.

    ! Set initial values for w
    w(ixI^S, rho_) = rho0
    w(ixI^S, mom(:)) = zero
    w(ixI^S,r_e) = E_r0
    w(ixI^S, e_) = t0*e_eq

  if(mype==0.and.first)then
    write(*,*)'fld_mu mean molecular weight=',fld_mu
    write(*,*)'initial density=',rho0
    write(*,*)'initial E0=',E_r0
    write(*,*)'initial e=',t0*e_eq
    write(*,*)'equilibrium e=',e_eq
    write(*,*)'forced deviation of factor t0=',t0
    first=.false.
  endif
  end subroutine initial_conditions

  subroutine set_output_vars(ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    w(ixO^S,delta_e_)=w(ixO^S,e_)/e_eq
    ! output the AMR level (assuming uniform grid blocks)
    !w(ixO^S,amr_)=dlog(((xprobmax1-xprobmin1)/domain_nx1)/dxlevel(1))/dlog(2.0d0)+1.0d0

  end subroutine set_output_vars

end module mod_usr
