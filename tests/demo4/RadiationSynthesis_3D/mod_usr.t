module mod_usr
  use mod_tdfluxrope
  use mod_eos, only: eos
  use mod_mhd
  implicit none

  double precision :: tcorona, rhocorona, prom_density_factor
  double precision :: prom_x0, prom_y0, prom_z0
  double precision :: prom_ax, prom_ay, prom_az

contains

  subroutine usr_init()
    call set_coordinate_system("Cartesian_3D")

    unit_length        = 5.d9
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

    d_TD99 = 3.0d9/unit_length
    L_TD99 = 3.5d9/unit_length
    R_TD99 = 11.0d9/unit_length
    a_TD99 = 2.39168d9/unit_length
    q_TD99 = 4.0d21/unit_magneticfield/unit_length**2
    Izero_TD99 = 1.0d0/(unit_magneticfield*unit_length)/const_c * &
         (-1.0d6)*2.99792456d9

    prom_density_factor = 100.0d0
    prom_x0 = 0.0d0
    prom_y0 = 0.0d0
    prom_z0 = 1.45d0
    prom_ax = 0.55d0
    prom_ay = 2.1d0
    prom_az = 0.45d0

    Itube = 2.d0*q_TD99*L_TD99*R_TD99 &
         /(L_TD99**2+R_TD99**2) &
         /sqrt(L_TD99**2+R_TD99**2) &
         /(log(8.d0*R_TD99/a_TD99)-1.5d0+Li_TD99/2.d0)
    Nt_TD99 = abs(Itube/Izero_TD99)*R_TD99**2/a_TD99**2

    if (mype == 0) then
      write(*,*) 'Radiation-synthesis TDm prominence test'
      write(*,*) 'rho_corona, T_corona:', rhocorona, tcorona
      write(*,*) 'prominence density factor:', prom_density_factor
      write(*,*) 'TD99 turns:', Nt_TD99
    end if
  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: Bf(ixI^S,1:ndir)
    logical :: prom(ixI^S)

    call TD99(ixI^L,ixO^L,x,Bf)

    w(ixO^S,rho_) = rhocorona
    w(ixO^S,mom(:)) = zero
    w(ixO^S,mag(1:ndir)) = Bf(ixO^S,1:ndir)
    if (mhd_energy) w(ixO^S,p_) = rhocorona*tcorona

    prom(ixI^S) = .false.
    call get_prominence_mask(ixI^L,ixO^L,x,prom)
    where (prom(ixO^S))
      w(ixO^S,rho_) = prom_density_factor*w(ixO^S,rho_)
    end where

    call eos%to_conserved(ixI^L,ixO^L,w,x)
  end subroutine initonegrid_usr

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
    integer, intent(in) :: ixI^L, ixO^L
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

    rvertical(ixO^S) = dsqrt(x(ixO^S,2)**2 + (x(ixO^S,3)+d_TD99)**2)
    rhoadius(ixO^S) = dsqrt(x(ixO^S,1)**2 + (rvertical(ixO^S)-R_TD99)**2)
    ell(ixO^S) = ((x(ixO^S,1)-prom_x0)/prom_ax)**2 &
         + ((x(ixO^S,2)-prom_y0)/prom_ay)**2 &
         + ((x(ixO^S,3)-prom_z0)/prom_az)**2

    prom(ixO^S) = ell(ixO^S) <= one .and. rhoadius(ixO^S) <= 1.35d0*a_TD99
  end subroutine get_prominence_mask

end module mod_usr
