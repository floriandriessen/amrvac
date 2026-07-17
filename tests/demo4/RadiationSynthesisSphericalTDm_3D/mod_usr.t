module mod_usr
  use mod_tdfluxrope
  use mod_eos, only: eos
  use mod_mhd
  implicit none

  double precision :: tcorona, rhocorona, prom_density_factor
  double precision :: atm_rho_floor, atm_rho0, atm_scale_height
  double precision :: atm_tchrom, atm_htra, atm_wtra
  double precision :: patch_rmin, patch_r0, patch_theta0, patch_phi0
  double precision :: prom_x0, prom_y0, prom_z0
  double precision :: prom_ax, prom_ay, prom_az

contains

  subroutine usr_init()
    call set_coordinate_system("spherical_3D")

    unit_length        = 1.d9
    unit_temperature   = 1.d6
    unit_numberdensity = 1.d10

    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid  => initonegrid_usr
    usr_refine_grid    => specialrefine_grid
    usr_aux_output     => specialvar_output
    usr_add_aux_names  => specialvarnames_output

    call mhd_activate()
  end subroutine usr_init

  subroutine initglobaldata_usr()
    double precision :: Itube, Nt_TD99

    Li_TD99 = 5.0d-1
    tcorona = 10.d0**5.8d0/unit_temperature
    rhocorona = 1.0d10/unit_numberdensity
    p_Bt_ratio = 0.94d0

    ! Local TD99 rope embedded in the small spherical patch.  The
    ! spherical cell centers are projected onto a tangent plane before
    ! TD99 is evaluated, so these parameters are in code units.
    d_TD99 = 5.0d0
    L_TD99 = 3.5d0
    R_TD99 = 8.05219014004735d0
    a_TD99 = 2.25d0
    q_TD99 = 1.0d0
    Izero_TD99 = -1.0d0

    prom_density_factor = 100.0d0
    atm_rho_floor = 0.08d0
    atm_rho0 = 1.07d0
    atm_scale_height = 6.5d0
    atm_tchrom = 2.0d4/unit_temperature
    atm_htra = 0.25d0
    atm_wtra = 0.08d0
    patch_rmin = 69.61d0
    patch_r0 = 0.5d0*(69.61d0+79.61d0)
    patch_theta0 = 0.5d0*dpi
    patch_phi0 = 0.d0

    prom_x0 = 0.0d0
    prom_y0 = 0.0d0
    prom_z0 = R_TD99-d_TD99
    prom_ax = 1.8d0
    prom_ay = 6.0d0
    prom_az = 1.2d0

    Itube = 2.d0*q_TD99*L_TD99*R_TD99 &
         /(L_TD99**2+R_TD99**2) &
         /sqrt(L_TD99**2+R_TD99**2) &
         /(log(8.d0*R_TD99/a_TD99)-1.5d0+Li_TD99/2.d0)
    Nt_TD99 = abs(Itube/Izero_TD99)*R_TD99**2/a_TD99**2

    if (mype == 0) then
      write(*,*) 'Spherical radiation-synthesis TD99 prominence test'
      write(*,*) 'stratified rho(h=50 Mm), T_corona:', &
           stratified_density(5.0d0), tcorona
      write(*,*) 'prominence density factor:', prom_density_factor
      write(*,*) 'local TD99 turns:', Nt_TD99
    end if
  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    logical :: prom(ixI^S)
    double precision :: xTD(ixI^S,1:ndim), BTD(ixI^S,1:ndir)

    call spherical_to_local_td99(ixI^L,ixO^L,x,xTD)
    call TD99(ixI^L,ixO^L,xTD,BTD)

    call fill_stratified_atmosphere(ixI^L,ixO^L,x,w)
    w(ixO^S,mom(:)) = zero
    call local_td99_field_to_spherical(ixI^L,ixO^L,BTD,w(ixI^S,mag(1):mag(ndir)))

    prom(ixI^S) = .false.
    call get_prominence_mask(ixI^L,ixO^L,x,prom)
    where (prom(ixO^S))
      w(ixO^S,rho_) = prom_density_factor*w(ixO^S,rho_)
    end where

    call eos%to_conserved(ixI^L,ixO^L,w,x)
  end subroutine initonegrid_usr

  subroutine fill_stratified_atmosphere(ixI^L,ixO^L,x,w)
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: height(ixI^S), temperature(ixI^S)

    height(ixO^S) = max(zero,x(ixO^S,1)-patch_rmin)
    w(ixO^S,rho_) = stratified_density(height(ixO^S))
    temperature(ixO^S) = stratified_temperature(height(ixO^S))
    if (mhd_energy) w(ixO^S,p_) = w(ixO^S,rho_)*temperature(ixO^S)
  end subroutine fill_stratified_atmosphere

  elemental double precision function stratified_density(height)
    double precision, intent(in) :: height

    stratified_density = atm_rho_floor + atm_rho0*exp(-max(zero,height)/atm_scale_height)
  end function stratified_density

  elemental double precision function stratified_temperature(height)
    double precision, intent(in) :: height

    stratified_temperature = atm_tchrom + half*(tcorona-atm_tchrom) * &
         (tanh((max(zero,height)-atm_htra)/atm_wtra)+one)
  end function stratified_temperature

  subroutine specialrefine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    logical :: prom(ixI^S)

    prom(ixI^S) = .false.
    call get_prominence_mask(ixI^L,ixO^L,x,prom)

    if (any(prom(ixO^S))) then
      if (level < refine_max_level) then
        refine = 1
        coarsen = 0
      else
        refine = -1
        coarsen = -1
      end if
    end if
  end subroutine specialrefine_grid

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision :: w(ixI^S,nw+nwauxio)
    double precision :: normconv(0:nw+nwauxio)

    double precision :: pth(ixI^S)
    logical :: prom(ixI^S)

    call eos%get_thermal_pressure(w,x,ixI^L,ixO^L,pth)
    prom(ixI^S) = .false.
    call get_prominence_mask(ixI^L,ixO^L,x,prom)

    w(ixO^S,nw+1) = pth(ixO^S)/w(ixO^S,rho_)
    w(ixO^S,nw+2) = zero
    where (prom(ixO^S)) w(ixO^S,nw+2) = one
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    character(len=*) :: varnames

    varnames = 'T prominence_mask'
  end subroutine specialvarnames_output

  subroutine get_prominence_mask(ixI^L,ixO^L,x,prom)
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    logical, intent(inout) :: prom(ixI^S)

    double precision :: ell(ixI^S), rvertical(ixI^S), rhoadius(ixI^S)
    double precision :: xTD(ixI^S,1:ndim)

    call spherical_to_local_td99(ixI^L,ixO^L,x,xTD)
    rvertical(ixO^S) = dsqrt(xTD(ixO^S,2)**2 + (xTD(ixO^S,3)+d_TD99)**2)
    rhoadius(ixO^S) = dsqrt(xTD(ixO^S,1)**2 + (rvertical(ixO^S)-R_TD99)**2)
    ell(ixO^S) = ((xTD(ixO^S,1)-prom_x0)/prom_ax)**2 &
         + ((xTD(ixO^S,2)-prom_y0)/prom_ay)**2 &
         + ((xTD(ixO^S,3)-prom_z0)/prom_az)**2

    prom(ixO^S) = ell(ixO^S) <= one .and. rhoadius(ixO^S) <= 1.35d0*a_TD99
  end subroutine get_prominence_mask

  subroutine spherical_to_local_td99(ixI^L,ixO^L,x,xTD)
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out) :: xTD(ixI^S,1:ndim)

    xTD(ixO^S,1) = patch_r0*dsin(patch_theta0)*(x(ixO^S,3)-patch_phi0)
    xTD(ixO^S,2) = patch_r0*(x(ixO^S,2)-patch_theta0)
    xTD(ixO^S,3) = x(ixO^S,1)-patch_rmin
  end subroutine spherical_to_local_td99

  subroutine local_td99_field_to_spherical(ixI^L,ixO^L,BTD,Bsph)
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: BTD(ixI^S,1:ndir)
    double precision, intent(out) :: Bsph(ixI^S,1:ndir)

    Bsph(ixO^S,1) = BTD(ixO^S,3)
    Bsph(ixO^S,2) = BTD(ixO^S,2)
    Bsph(ixO^S,3) = BTD(ixO^S,1)
  end subroutine local_td99_field_to_spherical

end module mod_usr
