module mod_mhd_eos
    use mod_global_parameters
    use mod_physics
    use mod_eos
    use mod_eos_container
    use mod_mhd_phys
    use mod_timing

    use mod_comm_lib, only: mpistop

    implicit none
    private

    public :: mhd_link_eos
contains

    !> Link the appropriate EOS conversion routines based on the selected EoS type
    subroutine mhd_link_eos()
      use mod_usr_methods

      ! === Primitive <-> Conserved routing ===
      if (eos%eos_type == 'FI') then
        if(mhd_hydrodynamic_e) then
          eos%to_primitive        => mhd_to_primitive_hde
          eos%to_conserved        => mhd_to_conserved_hde
        else if(mhd_semirelativistic) then
          if(mhd_energy) then
            eos%to_primitive        => mhd_to_primitive_semirelati
            eos%to_conserved        => mhd_to_conserved_semirelati
          else
            eos%to_primitive        => mhd_to_primitive_semirelati_noe
            eos%to_conserved        => mhd_to_conserved_semirelati_noe
          end if
        else
          if(has_equi_rho0) then
            eos%to_primitive        => mhd_to_primitive_split_rho
            eos%to_conserved        => mhd_to_conserved_split_rho
          else if(mhd_internal_e) then
            eos%to_primitive        => mhd_to_primitive_inte
            eos%to_conserved        => mhd_to_conserved_inte
          else if(mhd_energy) then
            eos%to_primitive         => mhd_to_primitive_origin
            eos%to_conserved         => mhd_to_conserved_origin
          else
            eos%to_primitive         => mhd_to_primitive_origin_noe
            eos%to_conserved         => mhd_to_conserved_origin_noe
          end if
        end if
      else if (eos%eos_type == 'LTE') then
        ! Guard unsupported combinations
        if (mhd_hydrodynamic_e) &
          call mpistop('LTE EoS not supported with mhd_hydrodynamic_e')
        if (mhd_semirelativistic) &
          call mpistop('LTE EoS not supported with mhd_semirelativistic')
        if (has_equi_rho0 .or. has_equi_pe0) &
          call mpistop('LTE EoS not supported with equilibrium splitting')
        if (.not. mhd_energy) &
          call mpistop('LTE EoS requires mhd_energy=.true.')
        ! Route LTE conversion based on energy formulation
        if (mhd_internal_e) then
          eos%to_primitive  => mhd_to_primitive_inte_LTE
          eos%to_conserved  => mhd_to_conserved_inte_LTE
        else
          eos%to_primitive  => mhd_to_primitive_origin_LTE
          eos%to_conserved  => mhd_to_conserved_origin_LTE
        end if
      else
        call mpistop('Error: Unknown MHD EOS type: ' // trim(eos%eos_type))
      end if

      phys_to_primitive       => eos%to_primitive
      phys_to_conserved       => eos%to_conserved
      phys_get_rho            => eos%get_rho
      phys_bind_eos_to_source => bind_eos_to_source

      ! === p_to_e (pressure -> total energy for origin, or eint for inte) ===
      if (mhd_internal_e) then
        eos%p_to_e => mhd_p_to_eint
      else
        eos%p_to_e => mhd_p_to_e
      end if

      ! === Sound speed and Gamma1 ===
      if (eos%eos_type == 'LTE' .and. eos%ionE) then
        eos%get_csound2 => mhd_get_csound2_LTE
        phys_get_gamma1 => mhd_get_gamma1_LTE
      else
        eos%get_csound2 => mhd_get_csound2_FI
        phys_get_gamma1 => mhd_get_gamma1_FI
      end if

      ! === Prolongation ===
      if (eos%eos_type == 'LTE' .and. eos%ionE) then
        eos%to_prolong    => mhd_to_prolong_LTE
        eos%from_prolong  => mhd_from_prolong_LTE
        phys_to_prolong   => mhd_to_prolong_LTE
        phys_from_prolong => mhd_from_prolong_LTE
      end if

      ! === Rfactor (matching HD: routing in link_eos, not phys_init) ===
      if (eos%eos_type == 'LTE') then
        eos%get_Rfactor => Rfactor_from_LTE
      else if(mhd_partial_ionization) then
        eos%get_Rfactor => Rfactor_from_temperature_ionization
        phys_update_temperature => mhd_update_temperature
      else if(associated(usr_Rfactor)) then
        eos%get_Rfactor => usr_Rfactor
      else
        eos%get_Rfactor => Rfactor_from_constant_ionization
      end if

      ! === Thermal pressure routing ===
      if(mhd_internal_e) then
        phys_get_pthermal        => mhd_get_pthermal_inte
        mhd_get_pthermal         => mhd_get_pthermal_inte
      else if(mhd_hydrodynamic_e) then
        phys_get_pthermal        => mhd_get_pthermal_hde
        mhd_get_pthermal         => mhd_get_pthermal_hde
      else if(mhd_semirelativistic) then
        phys_get_pthermal        => mhd_get_pthermal_semirelati
        mhd_get_pthermal         => mhd_get_pthermal_semirelati
      else if(mhd_energy) then
        phys_get_pthermal        => mhd_get_pthermal_origin
        mhd_get_pthermal         => mhd_get_pthermal_origin
      else
        phys_get_pthermal        => mhd_get_pthermal_noe
        mhd_get_pthermal         => mhd_get_pthermal_noe
      end if

      ! Override pthermal for LTE (delegates to eos%get_thermal_pressure at runtime)
      if (eos%eos_type == 'LTE') then
        phys_get_pthermal => mhd_get_pthermal_LTE
        mhd_get_pthermal  => mhd_get_pthermal_LTE
      end if

      ! === Temperature routing ===
      if(mhd_partial_ionization .or. eos%eos_type == 'LTE') then
        mhd_get_temperature => mhd_get_temperature_from_Te
      else
        if(mhd_internal_e) then
          if(has_equi_pe0 .and. has_equi_rho0) then
            mhd_get_temperature => mhd_get_temperature_from_eint_with_equi
          else
            mhd_get_temperature => mhd_get_temperature_from_eint
          end if
        else
          if(has_equi_pe0 .and. has_equi_rho0) then
            mhd_get_temperature => mhd_get_temperature_from_etot_with_equi
          else
            mhd_get_temperature => mhd_get_temperature_from_etot
          end if
        end if
      end if

      ! === Internal energy extraction ===
      if(mhd_internal_e) then
        phys_get_ei => mhd_get_ei_inte
      else if(mhd_energy) then
        phys_get_ei => mhd_get_ei_origin
      end if

      ! phys_e_to_ei / phys_ei_to_e assigned in mhd_phys_init (mod_mhd_phys.t)

    end subroutine mhd_link_eos

    !> Called by eos_finalise via phys_bind_eos_to_source to link
    !> TC and RC modules to EoS-aware function pointers.
    subroutine bind_eos_to_source()

      if (allocated(tc_fl)) then
        tc_fl%get_temperature_from_conserved => eos%get_temperature_from_etot
        if (eos%eos_type == 'LTE' .and. eos%ionE) then
          tc_fl%get_temperature_from_eint => get_temperature_from_eint_LTE_fast
        else
          tc_fl%get_temperature_from_eint => eos%get_temperature_from_eint
        end if
        tc_fl%get_rho => eos%get_rho
        ! Equilibrium-specific pointers
        if(has_equi_pe0 .and. has_equi_rho0 .and. mhd_equi_thermal) then
          tc_fl%has_equi = .true.
          tc_fl%get_temperature_equi => mhd_get_temperature_equi
          tc_fl%get_rho_equi => mhd_get_rho_equi
        else
          tc_fl%has_equi = .false.
        end if
      end if

      if (allocated(rc_fl)) then
        rc_fl%get_rho          => eos%get_rho
        rc_fl%get_pthermal     => eos%get_thermal_pressure
        rc_fl%get_var_Rfactor  => eos%get_Rfactor
        rc_fl%get_Te           => eos%get_Te
        ! Equilibrium-specific pointers
        if(has_equi_pe0 .and. has_equi_rho0 .and. mhd_equi_thermal) then
          rc_fl%has_equi = .true.
          rc_fl%get_rho_equi => mhd_get_rho_equi
          rc_fl%get_pthermal_equi => mhd_get_pe_equi
        else
          rc_fl%has_equi = .false.
        end if
      end if

    end subroutine bind_eos_to_source

    ! ========================================================================
    ! FI conversion routines (signatures updated for convert_condition interface)
    ! ========================================================================

    !> Transform primitive variables into conservative ones (origin, total energy)
    subroutine mhd_to_conserved_origin(ixI^L,ixO^L,w,x)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      double precision, intent(in)    :: x(ixI^S, 1:ndim)

      integer :: ix^D

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! Calculate total energy from pressure, kinetic and magnetic energy
        w(ix^D,e_)=w(ix^D,p_)*eos%inv_gamma_minus_1&
                  +half*((^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)&
                  +(^C&w(ix^D,b^C_)**2+))
        ! Convert velocity to momentum
        ^C&w(ix^D,m^C_)=w(ix^D,rho_)*w(ix^D,m^C_)\
      {end do\}

    end subroutine mhd_to_conserved_origin

    !> Transform primitive variables into conservative ones (no energy)
    subroutine mhd_to_conserved_origin_noe(ixI^L,ixO^L,w,x)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      double precision, intent(in)    :: x(ixI^S, 1:ndim)

      integer :: ix^D

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! Convert velocity to momentum
        ^C&w(ix^D,m^C_)=w(ix^D,rho_)*w(ix^D,m^C_)\
      {end do\}

    end subroutine mhd_to_conserved_origin_noe

    !> Transform primitive variables into conservative ones (hydrodynamic energy)
    subroutine mhd_to_conserved_hde(ixI^L,ixO^L,w,x)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      double precision, intent(in)    :: x(ixI^S, 1:ndim)

      integer :: ix^D

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! Calculate total energy from pressure, kinetic and magnetic energy
        w(ix^D,e_)=w(ix^D,p_)*eos%inv_gamma_minus_1&
                  +half*(^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)
        ! Convert velocity to momentum
        ^C&w(ix^D,m^C_)=w(ix^D,rho_)*w(ix^D,m^C_)\
      {end do\}

    end subroutine mhd_to_conserved_hde

    !> Transform primitive variables into conservative ones (internal energy)
    subroutine mhd_to_conserved_inte(ixI^L,ixO^L,w,x)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      double precision, intent(in)    :: x(ixI^S, 1:ndim)

      integer :: ix^D

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! Calculate internal energy from pressure
        w(ix^D,e_)=w(ix^D,p_)*eos%inv_gamma_minus_1
        ! Convert velocity to momentum
        ^C&w(ix^D,m^C_)=w(ix^D,rho_)*w(ix^D,m^C_)\
      {end do\}

    end subroutine mhd_to_conserved_inte

    !> Transform primitive variables into conservative ones (split rho)
    subroutine mhd_to_conserved_split_rho(ixI^L,ixO^L,w,x)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      double precision, intent(in)    :: x(ixI^S, 1:ndim)

      double precision :: rho
      integer :: ix^D

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        rho=w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,b0i)
        ! Calculate total energy from pressure, kinetic and magnetic energy
        w(ix^D,e_)=w(ix^D,p_)*eos%inv_gamma_minus_1&
                  +half*((^C&w(ix^D,m^C_)**2+)*rho&
                        +(^C&w(ix^D,b^C_)**2+))
        ! Convert velocity to momentum
        ^C&w(ix^D,m^C_)=rho*w(ix^D,m^C_)\
      {end do\}

    end subroutine mhd_to_conserved_split_rho

    !> Transform primitive variables into conservative ones (semirelativistic)
    subroutine mhd_to_conserved_semirelati(ixI^L,ixO^L,w,x)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      double precision, intent(in)    :: x(ixI^S, 1:ndim)

      double precision :: E(ixO^S,1:ndir), S(ixO^S,1:ndir)
      integer :: ix^D

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        {^IFTHREEC
        E(ix^D,1)=w(ix^D,b2_)*w(ix^D,m3_)-w(ix^D,b3_)*w(ix^D,m2_)
        E(ix^D,2)=w(ix^D,b3_)*w(ix^D,m1_)-w(ix^D,b1_)*w(ix^D,m3_)
        E(ix^D,3)=w(ix^D,b1_)*w(ix^D,m2_)-w(ix^D,b2_)*w(ix^D,m1_)
        S(ix^D,1)=E(ix^D,2)*w(ix^D,b3_)-E(ix^D,3)*w(ix^D,b2_)
        S(ix^D,2)=E(ix^D,3)*w(ix^D,b1_)-E(ix^D,1)*w(ix^D,b3_)
        S(ix^D,3)=E(ix^D,1)*w(ix^D,b2_)-E(ix^D,2)*w(ix^D,b1_)
        }
        {^IFTWOC
        E(ix^D,1)=zero
        ! switch 3 with 2 to add 3 when ^C from 1 to 2
        E(ix^D,2)=w(ix^D,b1_)*w(ix^D,m2_)-w(ix^D,b2_)*w(ix^D,m1_)
        S(ix^D,1)=-E(ix^D,2)*w(ix^D,b2_)
        S(ix^D,2)=E(ix^D,2)*w(ix^D,b1_)
        }
        {^IFONEC
        E(ix^D,1)=zero
        S(ix^D,1)=zero
        }
        if(mhd_internal_e) then
          ! internal energy
          w(ix^D,e_)=w(ix^D,p_)*eos%inv_gamma_minus_1
        else
          ! equation (9)
          ! Calculate total energy from internal, kinetic and magnetic energy
          w(ix^D,e_)=w(ix^D,p_)*eos%inv_gamma_minus_1&
                    +half*((^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)&
                    +(^C&w(ix^D,b^C_)**2+)&
                    +(^C&e(ix^D,^C)**2+)*eos%inv_squared_c)
        end if

        ! Convert velocity to momentum, equation (9)
        ^C&w(ix^D,m^C_)=w(ix^D,rho_)*w(ix^D,m^C_)+S(ix^D,^C)*eos%inv_squared_c\

      {end do\}

    end subroutine mhd_to_conserved_semirelati

    subroutine mhd_to_conserved_semirelati_noe(ixI^L,ixO^L,w,x)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      double precision, intent(in)    :: x(ixI^S, 1:ndim)

      double precision :: E(ixO^S,1:ndir), S(ixO^S,1:ndir)
      integer :: ix^D

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        {^IFTHREEC
        E(ix^D,1)=w(ix^D,b2_)*w(ix^D,m3_)-w(ix^D,b3_)*w(ix^D,m2_)
        E(ix^D,2)=w(ix^D,b3_)*w(ix^D,m1_)-w(ix^D,b1_)*w(ix^D,m3_)
        E(ix^D,3)=w(ix^D,b1_)*w(ix^D,m2_)-w(ix^D,b2_)*w(ix^D,m1_)
        S(ix^D,1)=E(ix^D,2)*w(ix^D,b3_)-E(ix^D,3)*w(ix^D,b2_)
        S(ix^D,2)=E(ix^D,3)*w(ix^D,b1_)-E(ix^D,1)*w(ix^D,b3_)
        S(ix^D,3)=E(ix^D,1)*w(ix^D,b2_)-E(ix^D,2)*w(ix^D,b1_)
        }
        {^IFTWOC
        E(ix^D,1)=zero
        ! switch 3 with 2 to add 3 when ^C from 1 to 2
        E(ix^D,2)=w(ix^D,b1_)*w(ix^D,m2_)-w(ix^D,b2_)*w(ix^D,m1_)
        S(ix^D,1)=-E(ix^D,2)*w(ix^D,b2_)
        S(ix^D,2)=E(ix^D,2)*w(ix^D,b1_)
        }
        {^IFONEC
        S(ix^D,1)=zero
        }
        ! Convert velocity to momentum, equation (9)
        ^C&w(ix^D,m^C_)=w(ix^D,rho_)*w(ix^D,m^C_)+S(ix^D,^C)*eos%inv_squared_c\

      {end do\}

    end subroutine mhd_to_conserved_semirelati_noe

    !> Transform conservative variables into primitive ones (origin)
    subroutine mhd_to_primitive_origin(ixI^L,ixO^L,w,x)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      double precision, intent(in)    :: x(ixI^S, 1:ndim)

      double precision                :: inv_rho
      integer :: ix^D

      if (fix_small_values) then
        call phys_handle_small_values(.false., w, x, ixI^L, ixO^L, 'mhd_to_primitive_origin')
      end if

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        inv_rho = 1.d0/w(ix^D,rho_)
        ! Convert momentum to velocity
        ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\
        ! Calculate pressure = (gamma-1) * (e-ek-eb)
        w(ix^D,p_)=eos%gamma_minus_1*(w(ix^D,e_)&
                  -half*(w(ix^D,rho_)*(^C&w(ix^D,m^C_)**2+)&
                    +(^C&w(ix^D,b^C_)**2+)))
      {end do\}

    end subroutine mhd_to_primitive_origin

    !> Transform conservative variables into primitive ones (no energy)
    subroutine mhd_to_primitive_origin_noe(ixI^L,ixO^L,w,x)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      double precision, intent(in)    :: x(ixI^S, 1:ndim)

      double precision                :: inv_rho
      integer :: ix^D

      if (fix_small_values) then
        call phys_handle_small_values(.false., w, x, ixI^L, ixO^L, 'mhd_to_primitive_origin_noe')
      end if

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        inv_rho = 1.d0/w(ix^D,rho_)
        ! Convert momentum to velocity
        ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\
      {end do\}

    end subroutine mhd_to_primitive_origin_noe

    !> Transform conservative variables into primitive ones (hde)
    subroutine mhd_to_primitive_hde(ixI^L,ixO^L,w,x)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      double precision, intent(in)    :: x(ixI^S, 1:ndim)

      double precision                :: inv_rho
      integer                         :: ix^D

      if (fix_small_values) then
        call phys_handle_small_values(.false., w, x, ixI^L, ixO^L, 'mhd_to_primitive_hde')
      end if

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        inv_rho = 1.d0/w(ix^D,rho_)
        ! Convert momentum to velocity
        ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\
        ! Calculate pressure = (gamma-1) * (e-ek)
        w(ix^D,p_)=eos%gamma_minus_1*(w(ix^D,e_)-half*w(ix^D,rho_)*(^C&w(ix^D,m^C_)**2+))
      {end do\}

    end subroutine mhd_to_primitive_hde

    !> Transform conservative variables into primitive ones (internal energy)
    subroutine mhd_to_primitive_inte(ixI^L,ixO^L,w,x)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      double precision, intent(in)    :: x(ixI^S, 1:ndim)

      double precision                :: inv_rho
      integer                         :: ix^D

      if (fix_small_values) then
        call phys_handle_small_values(.false., w, x, ixI^L, ixO^L, 'mhd_to_primitive_inte')
      end if

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! Calculate pressure = (gamma-1) * e_internal
        w(ix^D,p_)=w(ix^D,e_)*eos%gamma_minus_1
        ! Convert momentum to velocity
        inv_rho = 1.d0/w(ix^D,rho_)
        ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\
      {end do\}

    end subroutine mhd_to_primitive_inte

    !> Transform conservative variables into primitive ones (split rho)
    subroutine mhd_to_primitive_split_rho(ixI^L,ixO^L,w,x)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      double precision, intent(in)    :: x(ixI^S, 1:ndim)

      double precision :: inv_rho
      integer :: ix^D

      if (fix_small_values) then
        call phys_handle_small_values(.false., w, x, ixI^L, ixO^L, 'mhd_to_primitive_split_rho')
      end if

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        inv_rho=1.d0/(w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,b0i))
        ! Convert momentum to velocity
        ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\
        ! Calculate pressure = (gamma-1) * (e-ek-eb)
        w(ix^D,p_)=eos%gamma_minus_1*(w(ix^D,e_)&
                    -half*((w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,b0i))*&
                    (^C&w(ix^D,m^C_)**2+)+(^C&w(ix^D,b^C_)**2+)))
      {end do\}

    end subroutine mhd_to_primitive_split_rho

    !> Transform conservative variables into primitive ones (semirelativistic)
    subroutine mhd_to_primitive_semirelati(ixI^L,ixO^L,w,x)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      double precision, intent(in)    :: x(ixI^S, 1:ndim)

      double precision :: b(ixO^S,1:ndir), tmp, b2, gamma2, inv_rho
      integer :: ix^D

      if (fix_small_values) then
        call phys_handle_small_values(.false., w, x, ixI^L, ixO^L, 'mhd_to_primitive_semirelati')
      end if

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        b2=(^C&w(ix^D,b^C_)**2+)
        if(b2>smalldouble) then
          tmp=1.d0/sqrt(b2)
        else
          tmp=0.d0
        end if
        ^C&b(ix^D,^C)=w(ix^D,b^C_)*tmp\
        tmp=(^C&b(ix^D,^C)*w(ix^D,m^C_)+)

        inv_rho=1.d0/w(ix^D,rho_)
        ! Va^2/c^2
        b2=b2*inv_rho*eos%inv_squared_c
        ! equation (15)
        gamma2=1.d0/(1.d0+b2)
        ! Convert momentum to velocity
        ^C&w(ix^D,m^C_)=gamma2*(w(ix^D,m^C_)+b2*b(ix^D,^C)*tmp)*inv_rho\

        if(mhd_internal_e) then
          ! internal energy to pressure
          w(ix^D,p_)=eos%gamma_minus_1*w(ix^D,e_)
        else
          ! E=Bxv
          {^IFTHREEC
          b(ix^D,1)=w(ix^D,b2_)*w(ix^D,m3_)-w(ix^D,b3_)*w(ix^D,m2_)
          b(ix^D,2)=w(ix^D,b3_)*w(ix^D,m1_)-w(ix^D,b1_)*w(ix^D,m3_)
          b(ix^D,3)=w(ix^D,b1_)*w(ix^D,m2_)-w(ix^D,b2_)*w(ix^D,m1_)
          }
          {^IFTWOC
          b(ix^D,1)=zero
          b(ix^D,2)=w(ix^D,b1_)*w(ix^D,m2_)-w(ix^D,b2_)*w(ix^D,m1_)
          }
          {^IFONEC
          b(ix^D,1)=zero
          }
          ! Calculate pressure = (gamma-1) * (e-eK-eB-eE)
          w(ix^D,p_)=eos%gamma_minus_1*(w(ix^D,e_)&
                    -half*((^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)&
                    +(^C&w(ix^D,b^C_)**2+)&
                    +(^C&b(ix^D,^C)**2+)*eos%inv_squared_c))
        end if
      {end do\}

    end subroutine mhd_to_primitive_semirelati

    !> Transform conservative variables into primitive ones (semirelativistic noe)
    subroutine mhd_to_primitive_semirelati_noe(ixI^L,ixO^L,w,x)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      double precision, intent(in)    :: x(ixI^S, 1:ndim)

      double precision :: b(ixO^S,1:ndir),tmp,b2,gamma2,inv_rho
      integer :: ix^D, idir

      if (fix_small_values) then
        call phys_handle_small_values(.false., w, x, ixI^L, ixO^L, 'mhd_to_primitive_semirelati_noe')
      end if

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        b2=(^C&w(ix^D,b^C_)**2+)
        if(b2>smalldouble) then
          tmp=1.d0/sqrt(b2)
        else
          tmp=0.d0
        end if
        ^C&b(ix^D,^C)=w(ix^D,b^C_)*tmp\
        tmp=(^C&b(ix^D,^C)*w(ix^D,m^C_)+)

        inv_rho=1.d0/w(ix^D,rho_)
        ! Va^2/c^2
        b2=b2*inv_rho*eos%inv_squared_c
        ! equation (15)
        gamma2=1.d0/(1.d0+b2)
        ! Convert momentum to velocity
        ^C&w(ix^D,m^C_)=gamma2*(w(ix^D,m^C_)+b2*b(ix^D,^C)*tmp)*inv_rho\
      {end do\}

    end subroutine mhd_to_primitive_semirelati_noe

    ! ========================================================================
    ! LTE conversion routines
    ! ========================================================================

    !> LTE: Primitive (rho, v, p, B) -> Conserved (rho, rho*v, E_total, B)
    !> E_total = eint + 0.5*rho*v^2 + 0.5*B^2
    subroutine mhd_to_conserved_origin_LTE(ixI^L,ixO^L,w,x)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      double precision, intent(in)    :: x(ixI^S, 1:ndim)

      timeeos0 = MPI_WTIME()

      ! Convert p -> E_total (eint + eK + eB) via EoS table
      call mhd_p_to_e(ixI^L, ixO^L, w, x)
      ! Convert velocity to momentum
      ^C&w(ixO^S,m^C_)=w(ixO^S,rho_)*w(ixO^S,m^C_)\

      timeeos_conv=timeeos_conv+(MPI_WTIME()-timeeos0)

    end subroutine mhd_to_conserved_origin_LTE

    !> LTE: Primitive (rho, v, p, B) -> Conserved (rho, rho*v, eint, B)
    !> Internal energy formulation: e_ stores eint only
    subroutine mhd_to_conserved_inte_LTE(ixI^L,ixO^L,w,x)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      double precision, intent(in)    :: x(ixI^S, 1:ndim)

      timeeos0 = MPI_WTIME()

      ! Convert p -> eint via EoS table (no mechanical energy)
      call mhd_p_to_eint(ixI^L, ixO^L, w, x)
      ! Convert velocity to momentum
      ^C&w(ixO^S,m^C_)=w(ixO^S,rho_)*w(ixO^S,m^C_)\

      timeeos_conv=timeeos_conv+(MPI_WTIME()-timeeos0)

    end subroutine mhd_to_conserved_inte_LTE

    !> Convert pressure to total energy (eint + eK + eB) for origin formulation.
    !> Uses p2eint table for LTE ionE, with FI bypass for hot cells.
    subroutine mhd_p_to_e(ixI^L,ixO^L,w,x)
      use mod_global_parameters
      use mod_eos, only: p2eint_from_nH_p, saha_p_to_T
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      double precision, intent(in)    :: x(ixI^S, 1:ndim)

      integer :: ix^D
      double precision :: p_to_eint, p_over_rho
      double precision :: nH(ixI^S), nH_in(ixI^S), p_in(ixI^S)
      double precision :: T_solve, y_solve, eint_nH_solve

      if (eos%ionE) then
        call eos%get_nH(w, x, ixI^L, ixO^L, nH)
        if (eos%method /= 'analytic') then
          nH_in(ixO^S) = dlog10(nH(ixO^S))
          p_in(ixO^S) = dlog10(w(ixO^S,p_)) - nH_in(ixO^S)
        end if
      endif

      p_to_eint = eos%inv_gamma_minus_1
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        if (eos%ionE) then
          p_over_rho = w(ix^D,p_) / w(ix^D,rho_)
          if (p_over_rho > eos%p_rho_FI_threshold) then
            p_to_eint = eos%inv_gamma_minus_1 &
                + eos%eion_per_nH * nH(ix^D) / w(ix^D,p_)
          else if (eos%method == 'analytic') then
            call saha_p_to_T(nH(ix^D), w(ix^D,p_), &
                T_solve, y_solve, eint_nH_solve)
            w(ix^D,e_) = eint_nH_solve * nH(ix^D) &
                +half*((^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)&
                +(^C&w(ix^D,b^C_)**2+))
            cycle
          else
            p_to_eint = p2eint_from_nH_p(nH_in(ix^D), p_in(ix^D))
          end if
        end if
        w(ix^D,e_)=w(ix^D,p_)*p_to_eint&
                  +half*((^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)&
                  +(^C&w(ix^D,b^C_)**2+))
      {end do\}

    end subroutine mhd_p_to_e

    !> Convert pressure to internal energy only (for internal_e formulation).
    !> Uses p2eint table for LTE ionE, with FI bypass for hot cells.
    subroutine mhd_p_to_eint(ixI^L,ixO^L,w,x)
      use mod_global_parameters
      use mod_eos, only: p2eint_from_nH_p, saha_p_to_T
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      double precision, intent(in)    :: x(ixI^S, 1:ndim)

      integer :: ix^D
      double precision :: p_to_eint, p_over_rho
      double precision :: nH(ixI^S), nH_in(ixI^S), p_in(ixI^S)
      double precision :: T_solve, y_solve, eint_nH_solve

      if (eos%ionE) then
        call eos%get_nH(w, x, ixI^L, ixO^L, nH)
        if (eos%method /= 'analytic') then
          nH_in(ixO^S) = dlog10(nH(ixO^S))
          p_in(ixO^S) = dlog10(w(ixO^S,p_)) - nH_in(ixO^S)
        end if
      endif

      p_to_eint = eos%inv_gamma_minus_1
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        if (eos%ionE) then
          p_over_rho = w(ix^D,p_) / w(ix^D,rho_)
          if (p_over_rho > eos%p_rho_FI_threshold) then
            p_to_eint = eos%inv_gamma_minus_1 &
                + eos%eion_per_nH * nH(ix^D) / w(ix^D,p_)
          else if (eos%method == 'analytic') then
            call saha_p_to_T(nH(ix^D), w(ix^D,p_), &
                T_solve, y_solve, eint_nH_solve)
            w(ix^D,e_) = eint_nH_solve * nH(ix^D)
            cycle
          else
            p_to_eint = p2eint_from_nH_p(nH_in(ix^D), p_in(ix^D))
          end if
        end if
        w(ix^D,e_)=w(ix^D,p_)*p_to_eint
      {end do\}

    end subroutine mhd_p_to_eint

    !> LTE: Conserved (rho, rho*v, E_total, B) -> Primitive (rho, v, p, B)
    !> Uses T and neOnH tables with FI bypass for hot cells.
    subroutine mhd_to_primitive_origin_LTE(ixI^L,ixO^L,w,x)
      use mod_global_parameters
      use mod_eos, only: T_from_nH_eint, y_from_nH_eint, &
          saha_T_from_nH_eint
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      double precision, intent(in)    :: x(ixI^S, 1:ndim)

      double precision :: inv_rho, eint_val, eint_in, T_loc, y_loc
      double precision :: nH(ixI^S), log_nH(ixI^S)
      integer :: ix^D

      timeeos0 = MPI_WTIME()

      call eos%get_nH(w, x, ixI^L, ixO^L, nH)
      if (eos%ionE) then
        log_nH(ixO^S) = dlog10(nH(ixO^S))
      end if

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        inv_rho = 1.d0/w(ix^D,rho_)
        ! Convert momentum to velocity
        ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\
        ! Extract internal energy: eint = E - 0.5*rho*v^2 - 0.5*B^2
        eint_val = w(ix^D,e_) &
                  - half*(w(ix^D,rho_)*(^C&w(ix^D,m^C_)**2+) &
                  + (^C&w(ix^D,b^C_)**2+))
        ! Floor eint to prevent unphysical values
        if (eos%method /= 'analytic') then
          eint_val = max(eint_val, nH(ix^D) * 10.0d0**eos%T%var2_min)
        end if
        eint_val = max(eint_val, smalldouble)

        if (eos%ionE) then
          if (eint_val * inv_rho > eos%eint_rho_FI_threshold) then
            ! FI bypass: p = (gamma-1) * (eint - eion*nH)
            w(ix^D,p_) = eos%gamma_minus_1 &
                * (eint_val - eos%eion_per_nH * nH(ix^D))
          else if (eos%method == 'analytic') then
            ! Analytical Saha: solve for T and y from eint
            call saha_T_from_nH_eint(nH(ix^D), &
                eint_val / nH(ix^D), T_loc, y_loc)
            w(ix^D,p_) = nH(ix^D) &
                * (1.0d0 + eos%He_abundance + y_loc) * T_loc
          else
            ! Ionisation zone: table lookup
            eint_in = dlog10(eint_val) - log_nH(ix^D)
            T_loc = T_from_nH_eint(log_nH(ix^D), eint_in)
            y_loc = y_from_nH_eint(log_nH(ix^D), eint_in)
            w(ix^D,p_) = nH(ix^D) &
                * (1.0d0 + eos%He_abundance + y_loc) * T_loc
          end if
        else
          ! FI without ionE tables
          w(ix^D,p_) = eos%gamma_minus_1 * eint_val
        end if
      {end do\}

      if (fix_small_values) then
        call phys_handle_small_values(.false., w, x, ixI^L, ixO^L, &
            'mhd_to_primitive_origin_LTE')
      end if

      timeeos_conv=timeeos_conv+(MPI_WTIME()-timeeos0)

    end subroutine mhd_to_primitive_origin_LTE

    !> LTE: Conserved (rho, rho*v, eint, B) -> Primitive (rho, v, p, B)
    !> Internal energy formulation: eint = w(e_) directly.
    subroutine mhd_to_primitive_inte_LTE(ixI^L,ixO^L,w,x)
      use mod_global_parameters
      use mod_eos, only: T_from_nH_eint, y_from_nH_eint, &
          saha_T_from_nH_eint
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      double precision, intent(in)    :: x(ixI^S, 1:ndim)

      double precision :: inv_rho, eint_val, eint_in, T_loc, y_loc
      double precision :: nH(ixI^S), log_nH(ixI^S)
      integer :: ix^D

      timeeos0 = MPI_WTIME()

      call eos%get_nH(w, x, ixI^L, ixO^L, nH)
      if (eos%ionE) then
        log_nH(ixO^S) = dlog10(nH(ixO^S))
      end if

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! eint = w(e_) directly for internal energy formulation
        eint_val = w(ix^D,e_)
        eint_val = max(eint_val, nH(ix^D) * 10.0d0**eos%T%var2_min)

        if (eos%ionE) then
          if (eint_val / w(ix^D,rho_) > eos%eint_rho_FI_threshold) then
            w(ix^D,p_) = eos%gamma_minus_1 &
                * (eint_val - eos%eion_per_nH * nH(ix^D))
          else if (eos%method == 'analytic') then
            call saha_T_from_nH_eint(nH(ix^D), &
                eint_val / nH(ix^D), T_loc, y_loc)
            w(ix^D,p_) = nH(ix^D) &
                * (1.0d0 + eos%He_abundance + y_loc) * T_loc
          else
            eint_in = dlog10(eint_val) - log_nH(ix^D)
            T_loc = T_from_nH_eint(log_nH(ix^D), eint_in)
            y_loc = y_from_nH_eint(log_nH(ix^D), eint_in)
            w(ix^D,p_) = nH(ix^D) &
                * (1.0d0 + eos%He_abundance + y_loc) * T_loc
          end if
        else
          w(ix^D,p_) = eos%gamma_minus_1 * eint_val
        end if

        ! Convert momentum to velocity
        inv_rho = 1.d0/w(ix^D,rho_)
        ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\
      {end do\}

      if (fix_small_values) then
        call phys_handle_small_values(.false., w, x, ixI^L, ixO^L, &
            'mhd_to_primitive_inte_LTE')
      end if

      timeeos_conv=timeeos_conv+(MPI_WTIME()-timeeos0)

    end subroutine mhd_to_primitive_inte_LTE

    ! ========================================================================
    ! LTE prolongation routines
    ! ========================================================================

    !> Convert conserved (rho, rho*v, E_total, B) to prolongation form (rho, v, T, B).
    !> Temperature is stored in the p_ slot for AMR interpolation.
    !> T is smoothest through the ionisation zone (conduction-dominated, linear gradient).
    subroutine mhd_to_prolong_LTE(ixI^L,ixO^L,w,x)
      use mod_global_parameters
      use mod_eos, only: T_from_nH_eint, saha_T_from_nH_eint
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      double precision, intent(in)    :: x(ixI^S, 1:ndim)

      double precision :: inv_rho, eint_val
      double precision :: nH(ixI^S), log_nH(ixI^S)
      integer :: ix^D

      call eos%get_nH(w, x, ixI^L, ixO^L, nH)
      log_nH(ixO^S) = dlog10(nH(ixO^S))

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        inv_rho = 1.d0/w(ix^D,rho_)
        ! Convert momentum to velocity
        ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\
        if (mhd_internal_e) then
          eint_val = w(ix^D,e_)
        else
          ! Compute eint = E - 0.5*rho*v^2 - 0.5*B^2
          eint_val = w(ix^D,e_) &
              - half*(w(ix^D,rho_)*(^C&w(ix^D,m^C_)**2+) &
              + (^C&w(ix^D,b^C_)**2+))
        end if
        ! Floor eint to prevent unphysical values
        if (eos%method /= 'analytic') then
          eint_val = max(eint_val, nH(ix^D) * 10.0d0**eos%T%var2_min)
        end if
        eint_val = max(eint_val, smalldouble)

        if (eint_val * inv_rho > eos%eint_rho_FI_threshold) then
          ! FI: T = (gamma-1)*(eint - eion*nH) / (nH * n_per_nH_FI)
          w(ix^D,p_) = eos%gamma_minus_1 &
              * (eint_val - eos%eion_per_nH * nH(ix^D)) &
              / (nH(ix^D) * eos%n_per_nH_FI)
        else if (eos%method == 'analytic') then
          ! Analytical Saha: solve for T from eint
          block
            double precision :: T_loc, y_loc
            call saha_T_from_nH_eint(nH(ix^D), &
                eint_val / nH(ix^D), T_loc, y_loc)
            w(ix^D,p_) = T_loc
          end block
        else
          ! Ionisation zone: T from table
          w(ix^D,p_) = T_from_nH_eint( &
              log_nH(ix^D), &
              dlog10(eint_val) - log_nH(ix^D))
        end if
      {end do\}

    end subroutine mhd_to_prolong_LTE

    !> Convert prolongation form (rho, v, T, B) to conserved (rho, rho*v, E_total, B).
    !> T is read from the p_ slot. Uses eint_nH_from_T table for back-conversion.
    subroutine mhd_from_prolong_LTE(ixI^L,ixO^L,w,x)
      use mod_global_parameters
      use mod_eos, only: eint_nH_from_T, saha_eint_from_nH_T
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      double precision, intent(in)    :: x(ixI^S, 1:ndim)

      double precision :: T_val, eint_val, T_FI
      double precision :: nH(ixI^S), log_nH(ixI^S)
      integer :: ix^D

      ! Temperature above which gas is fully ionised
      T_FI = (eos%eint_rho_FI_threshold &
          * eos%nH2rhoFactor - eos%eion_per_nH) &
          * eos%gamma_minus_1 / eos%n_per_nH_FI

      call eos%get_nH(w, x, ixI^L, ixO^L, nH)
      log_nH(ixO^S) = dlog10(nH(ixO^S))

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        T_val = w(ix^D,p_)  ! T stored in p_ slot
        if (T_val > T_FI) then
          ! FI: eint = nH*(n_per_nH*T/(gamma-1) + eion)
          eint_val = nH(ix^D) &
              * (eos%n_per_nH_FI * T_val * eos%inv_gamma_minus_1 &
              + eos%eion_per_nH)
        else if (eos%method == 'analytic') then
          ! Analytical Saha: eint from T directly
          eint_val = saha_eint_from_nH_T(nH(ix^D), T_val) * nH(ix^D)
        else
          ! Ionisation zone: eint/nH from T table
          eint_val = eint_nH_from_T( &
              log_nH(ix^D), &
              dlog10(max(T_val, 10.0d0**eos%eint_from_T%var2_min))) &
              * nH(ix^D)
        end if
        if (mhd_internal_e) then
          w(ix^D,e_) = eint_val
        else
          ! E = eint + 0.5*rho*v^2 + 0.5*B^2
          w(ix^D,e_) = eint_val &
              + half*(w(ix^D,rho_)*(^C&w(ix^D,m^C_)**2+) &
              + (^C&w(ix^D,b^C_)**2+))
        end if
        ! Convert velocity to momentum
        ^C&w(ix^D,m^C_)=w(ix^D,rho_)*w(ix^D,m^C_)\
      {end do\}

    end subroutine mhd_from_prolong_LTE

    ! ========================================================================
    ! Sound speed routines
    ! ========================================================================

    !> Acoustic sound speed squared for FI (constant gamma) EoS.
    !> Expects w in primitive form: w(p_) = pressure, w(rho_) = density.
    subroutine mhd_get_csound2_FI(w, x, ixI^L, ixO^L, cs2)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(in)    :: w(ixI^S, nw)
      double precision, intent(in)    :: x(ixI^S, 1:ndim)
      double precision, intent(out)   :: cs2(ixI^S)

      timeeos0 = MPI_WTIME()

      cs2(ixO^S) = eos%gamma * w(ixO^S, p_) / w(ixO^S, rho_)

      timeeos_csound = timeeos_csound + (MPI_WTIME()-timeeos0)

    end subroutine mhd_get_csound2_FI

    !> Acoustic sound speed squared for LTE+IonE EoS using Gamma_1 table.
    !> Expects w in primitive form. Uses FI bypass for hot cells.
    subroutine mhd_get_csound2_LTE(w, x, ixI^L, ixO^L, cs2)
      use mod_global_parameters
      use mod_eos, only: gamma1_from_nH_p, gamma1_from_nH_T_analytic
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(in)    :: w(ixI^S, nw)
      double precision, intent(in)    :: x(ixI^S, 1:ndim)
      double precision, intent(out)   :: cs2(ixI^S)

      double precision :: nH_val, log_nH, log_p_nH, g1, p_over_rho
      integer :: ix^D

      timeeos0 = MPI_WTIME()

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        p_over_rho = w(ix^D, p_) / w(ix^D, rho_)
        if (p_over_rho > eos%p_rho_FI_threshold) then
          cs2(ix^D) = eos%gamma * p_over_rho
        else
          nH_val = w(ix^D, rho_) / eos%nH2rhoFactor
          if (eos%gamma1_method == 'effective') then
            if (iw_te > 0 .and. w(ix^D,iw_te) > 0.0d0 .and. iw_ne > 0) then
              block
                double precision :: y_loc, eint_loc
                y_loc = w(ix^D,iw_ne) / nH_val
                eint_loc = eos%inv_gamma_minus_1 * (1.0d0 + y_loc) * nH_val * w(ix^D,iw_te)
                if (eos%ionE) eint_loc = eint_loc &
                    + y_loc * eos%eion_per_nH * nH_val
                if (eint_loc > 0.0d0) then
                  g1 = 1.0d0 + w(ix^D,p_) / eint_loc
                else
                  g1 = eos%gamma
                end if
              end block
            else
              g1 = eos%gamma
            end if
          else if (eos%method == 'analytic') then
            if (iw_te > 0 .and. w(ix^D,iw_te) > 0.0d0) then
              g1 = gamma1_from_nH_T_analytic(nH_val, w(ix^D,iw_te))
            else
              g1 = eos%gamma
            end if
          else
            log_nH = dlog10(nH_val)
            log_p_nH = dlog10(w(ix^D, p_) / nH_val)
            g1 = gamma1_from_nH_p(log_nH, log_p_nH)
          end if
          cs2(ix^D) = g1 * p_over_rho
        end if
      {end do\}

      timeeos_csound = timeeos_csound + (MPI_WTIME()-timeeos0)

    end subroutine mhd_get_csound2_LTE

    subroutine mhd_get_gamma1_FI(w, x, ixI^L, ixO^L, gamma1)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(in)    :: w(ixI^S, nw)
      double precision, intent(in)    :: x(ixI^S, 1:ndim)
      double precision, intent(out)   :: gamma1(ixI^S)

      gamma1(ixO^S) = eos%gamma

    end subroutine mhd_get_gamma1_FI

    subroutine mhd_get_gamma1_LTE(w, x, ixI^L, ixO^L, gamma1)
      use mod_global_parameters
      use mod_eos, only: gamma1_from_nH_p, gamma1_from_nH_T_analytic
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(in)    :: w(ixI^S, nw)
      double precision, intent(in)    :: x(ixI^S, 1:ndim)
      double precision, intent(out)   :: gamma1(ixI^S)

      double precision :: nH_val, p_over_rho
      integer :: ix^D

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        p_over_rho = w(ix^D, p_) / w(ix^D, rho_)
        if (p_over_rho > eos%p_rho_FI_threshold) then
          gamma1(ix^D) = eos%gamma
        else
          nH_val = w(ix^D, rho_) / eos%nH2rhoFactor
          if (eos%gamma1_method == 'effective') then
            if (iw_te > 0 .and. w(ix^D,iw_te) > 0.0d0 .and. iw_ne > 0) then
              block
                double precision :: y_g, eint_g
                y_g = w(ix^D,iw_ne) / nH_val
                eint_g = eos%inv_gamma_minus_1 * (1.0d0 + y_g) * nH_val * w(ix^D,iw_te)
                if (eos%ionE) eint_g = eint_g &
                    + y_g * eos%eion_per_nH * nH_val
                if (eint_g > 0.0d0) then
                  gamma1(ix^D) = 1.0d0 + w(ix^D,p_) / eint_g
                else
                  gamma1(ix^D) = eos%gamma
                end if
              end block
            else
              gamma1(ix^D) = eos%gamma
            end if
          else if (eos%method == 'analytic') then
            if (iw_te > 0 .and. w(ix^D,iw_te) > 0.0d0) then
              gamma1(ix^D) = gamma1_from_nH_T_analytic(nH_val, w(ix^D,iw_te))
            else
              gamma1(ix^D) = eos%gamma
            end if
          else
            gamma1(ix^D) = gamma1_from_nH_p(dlog10(nH_val), &
                dlog10(w(ix^D, p_) / nH_val))
          end if
        end if
      {end do\}

    end subroutine mhd_get_gamma1_LTE

    ! ========================================================================
    ! Rfactor routines (matching HD: EoS functions live in EoS module)
    ! ========================================================================

    !> Rfactor = p/(rho*T) for constant ionisation degree
    subroutine Rfactor_from_constant_ionization(w,x,ixI^L,ixO^L,Rfactor)
      use mod_global_parameters
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S,1:nw)
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision, intent(out):: Rfactor(ixI^S)

      Rfactor(ixO^S)=RR

    end subroutine Rfactor_from_constant_ionization

    !> Rfactor from LTE EOS: R = (1 + He + ne/nH) / nH2rhoFactor
    !> Uses stored Ne_ wextra field from update_eos_LTE.
    subroutine Rfactor_from_LTE(w,x,ixI^L,ixO^L,Rfactor)
        use mod_global_parameters
        integer, intent(in) :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S,1:nw)
        double precision, intent(in) :: x(ixI^S,1:ndim)
        double precision, intent(out):: Rfactor(ixI^S)

        double precision :: y_nH
        integer :: ix^D

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
          y_nH = w(ix^D, iw_ne) * eos%nH2rhoFactor / w(ix^D, iw_rho)
          Rfactor(ix^D) = (1.0d0 + eos%He_abundance + y_nH) / eos%nH2rhoFactor
        {end do\}

    end subroutine Rfactor_from_LTE

    !> Rfactor from partial ionisation temperature lookup (Leenaarts et al. 2012)
    subroutine Rfactor_from_temperature_ionization(w,x,ixI^L,ixO^L,Rfactor)
      use mod_global_parameters
      use mod_ionization_degree
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S,1:nw)
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision, intent(out):: Rfactor(ixI^S)

      double precision :: iz_H(ixO^S),iz_He(ixO^S)

      call ionization_degree_from_temperature(ixI^L,ixO^L,w(ixI^S,Te_),iz_H,iz_He)
      ! assume the first and second ionization of Helium have the same degree
      Rfactor(ixO^S)=(1.d0+iz_H(ixO^S)+0.1d0*(1.d0+iz_He(ixO^S)*(1.d0+iz_He(ixO^S))))/(2.d0+3.d0*He_abundance)

    end subroutine Rfactor_from_temperature_ionization

    ! ========================================================================
    ! Thermal pressure routines (moved from mod_mhd_phys.t)
    ! ========================================================================

    !> Calculate isothermal thermal pressure
    subroutine mhd_get_pthermal_noe(w,x,ixI^L,ixO^L,pth)
      use mod_global_parameters

      integer, intent(in)          :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S,nw)
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision, intent(out):: pth(ixI^S)

      if(has_equi_rho0) then
        pth(ixO^S)=mhd_adiab*(w(ixO^S,rho_)+block%equi_vars(ixO^S,equi_rho0_,0))**mhd_gamma
      else
        pth(ixO^S)=mhd_adiab*w(ixO^S,rho_)**mhd_gamma
      end if

    end subroutine mhd_get_pthermal_noe

    !> Calculate thermal pressure from internal energy
    subroutine mhd_get_pthermal_inte(w,x,ixI^L,ixO^L,pth)
      use mod_global_parameters
      use mod_small_values, only: trace_small_values

      integer, intent(in)          :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S,nw)
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision, intent(out):: pth(ixI^S)

      integer :: iw, ix^D

     {do ix^DB= ixOmin^DB,ixOmax^DB\}
        if(has_equi_pe0) then
          pth(ix^D)=eos%gamma_minus_1*w(ix^D,e_)+block%equi_vars(ix^D,equi_pe0_,0)
        else
          pth(ix^D)=eos%gamma_minus_1*w(ix^D,e_)
        end if
        if(fix_small_values.and.pth(ix^D)<small_pressure) pth(ix^D)=small_pressure
     {end do\}

      if(check_small_values.and..not.fix_small_values) then
       {do ix^DB= ixOmin^DB,ixOmax^DB\}
          if(pth(ix^D)<small_pressure) then
            write(*,*) "Error: small value of gas pressure",pth(ix^D),&
                 " encountered when call mhd_get_pthermal_inte"
            write(*,*) "Iteration: ", it, " Time: ", global_time
            write(*,*) "Location: ", x(ix^D,:)
            write(*,*) "Cell number: ", ix^D
            do iw=1,nw
              write(*,*) trim(cons_wnames(iw)),": ",w(ix^D,iw)
            end do
            if(trace_small_values) write(*,*) sqrt(pth(ix^D)-bigdouble)
            write(*,*) "Saving status at the previous time step"
            crash=.true.
          end if
       {end do\}
      end if

    end subroutine mhd_get_pthermal_inte

    !> Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho-b**2/2) within ixO^L
    subroutine mhd_get_pthermal_origin(w,x,ixI^L,ixO^L,pth)
      use mod_global_parameters
      use mod_small_values, only: trace_small_values

      integer, intent(in)          :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S,nw)
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision, intent(out):: pth(ixI^S)

      integer :: iw, ix^D

     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        if(has_equi_rho0) then
          pth(ix^D)=eos%gamma_minus_1*(w(ix^D,e_)-half*((^C&w(ix^D,m^C_)**2+)/(w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,0))&
               +(^C&w(ix^D,b^C_)**2+)))+block%equi_vars(ix^D,equi_pe0_,0)
        else
          pth(ix^D)=eos%gamma_minus_1*(w(ix^D,e_)-half*((^C&w(ix^D,m^C_)**2+)/w(ix^D,rho_)&
               +(^C&w(ix^D,b^C_)**2+)))
        end if
        if(fix_small_values.and.pth(ix^D)<small_pressure) pth(ix^D)=small_pressure
     {end do\}

      if(check_small_values.and..not.fix_small_values) then
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
          if(pth(ix^D)<small_pressure) then
            write(*,*) "Error: small value of gas pressure",pth(ix^D),&
                 " encountered when call mhd_get_pthermal"
            write(*,*) "Iteration: ", it, " Time: ", global_time
            write(*,*) "Location: ", x(ix^D,:)
            write(*,*) "Cell number: ", ix^D
            do iw=1,nw
              write(*,*) trim(cons_wnames(iw)),": ",w(ix^D,iw)
            end do
            if(trace_small_values) write(*,*) sqrt(pth(ix^D)-bigdouble)
            write(*,*) "Saving status at the previous time step"
            crash=.true.
          end if
       {end do\}
      end if

    end subroutine mhd_get_pthermal_origin

    !> Calculate thermal pressure for LTE EoS (delegates to eos%get_thermal_pressure)
    subroutine mhd_get_pthermal_LTE(w,x,ixI^L,ixO^L,pth)
      use mod_global_parameters
      use mod_small_values, only: trace_small_values

      integer, intent(in)          :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S,nw)
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision, intent(out):: pth(ixI^S)

      integer :: iw, ix^D

      call eos%get_thermal_pressure(w, x, ixI^L, ixO^L, pth)

      if(fix_small_values) then
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
          if(pth(ix^D)<small_pressure) pth(ix^D)=small_pressure
       {end do\}
      else if(check_small_values) then
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
          if(pth(ix^D)<small_pressure) then
            write(*,*) "Error: small value of gas pressure",pth(ix^D),&
                 " encountered when call mhd_get_pthermal_LTE"
            write(*,*) "Iteration: ", it, " Time: ", global_time
            write(*,*) "Location: ", x(ix^D,:)
            write(*,*) "Cell number: ", ix^D
            do iw=1,nw
              write(*,*) trim(cons_wnames(iw)),": ",w(ix^D,iw)
            end do
            if(trace_small_values) write(*,*) sqrt(pth(ix^D)-bigdouble)
            write(*,*) "Saving status at the previous time step"
            crash=.true.
          end if
       {end do\}
      end if

    end subroutine mhd_get_pthermal_LTE

    !> Calculate thermal pressure for semirelativistic MHD
    subroutine mhd_get_pthermal_semirelati(w,x,ixI^L,ixO^L,pth)
      use mod_global_parameters
      use mod_small_values, only: trace_small_values

      integer, intent(in)          :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S,nw)
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision, intent(out):: pth(ixI^S)

      double precision :: b(ixO^S,1:ndir), v(ixO^S,1:ndir), tmp, b2, gamma2, inv_rho
      integer :: iw, ix^D

     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        b2=(^C&w(ix^D,b^C_)**2+)
        if(b2>smalldouble) then
          tmp=1.d0/sqrt(b2)
        else
          tmp=0.d0
        end if
        ^C&b(ix^D,^C)=w(ix^D,b^C_)*tmp\
        tmp=(^C&b(ix^D,^C)*w(ix^D,m^C_)+)

        inv_rho=1.d0/w(ix^D,rho_)
        ! Va^2/c^2
        b2=b2*inv_rho*eos%inv_squared_c
        ! equation (15)
        gamma2=1.d0/(1.d0+b2)
        ! Convert momentum to velocity
        ^C&v(ix^D,^C)=gamma2*(w(ix^D,m^C_)+b2*b(ix^D,^C)*tmp)*inv_rho\

        ! E=Bxv
        {^IFTHREEC
        b(ix^D,1)=w(ix^D,b2_)*v(ix^D,3)-w(ix^D,b3_)*v(ix^D,2)
        b(ix^D,2)=w(ix^D,b3_)*v(ix^D,1)-w(ix^D,b1_)*v(ix^D,3)
        b(ix^D,3)=w(ix^D,b1_)*v(ix^D,2)-w(ix^D,b2_)*v(ix^D,1)
        }
        {^IFTWOC
        b(ix^D,1)=zero
        b(ix^D,2)=w(ix^D,b1_)*v(ix^D,2)-w(ix^D,b2_)*v(ix^D,1)
        }
        {^IFONEC
        b(ix^D,1)=zero
        }
        ! Calculate pressure = (gamma-1) * (e-eK-eB-eE)
        pth(ix^D)=eos%gamma_minus_1*(w(ix^D,e_)&
                   -half*((^C&v(ix^D,^C)**2+)*w(ix^D,rho_)&
                   +(^C&w(ix^D,b^C_)**2+)&
                   +(^C&b(ix^D,^C)**2+)*eos%inv_squared_c))
        if(fix_small_values.and.pth(ix^D)<small_pressure) pth(ix^D)=small_pressure
     {end do\}

      if(check_small_values.and..not.fix_small_values) then
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
          if(pth(ix^D)<small_pressure) then
            write(*,*) "Error: small value of gas pressure",pth(ix^D),&
                 " encountered when call mhd_get_pthermal_semirelati"
            write(*,*) "Iteration: ", it, " Time: ", global_time
            write(*,*) "Location: ", x(ix^D,:)
            write(*,*) "Cell number: ", ix^D
            do iw=1,nw
              write(*,*) trim(cons_wnames(iw)),": ",w(ix^D,iw)
            end do
            if(trace_small_values) write(*,*) sqrt(pth(ix^D)-bigdouble)
            write(*,*) "Saving status at the previous time step"
            crash=.true.
          end if
       {end do\}
      end if

    end subroutine mhd_get_pthermal_semirelati

    !> Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho) within ixO^L
    subroutine mhd_get_pthermal_hde(w,x,ixI^L,ixO^L,pth)
      use mod_global_parameters
      use mod_small_values, only: trace_small_values

      integer, intent(in)          :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S,nw)
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision, intent(out):: pth(ixI^S)

      integer :: iw, ix^D

     {do ix^DB= ixOmin^DB,ixOmax^DB\}
        pth(ix^D)=eos%gamma_minus_1*(w(ix^D,e_)-half*((^C&w(ix^D,m^C_)**2+)/w(ix^D,rho_)))
        if(fix_small_values.and.pth(ix^D)<small_pressure) pth(ix^D)=small_pressure
     {end do\}
      if(check_small_values.and..not.fix_small_values) then
        {do ix^DB= ixOmin^DB,ixOmax^DB\}
          if(pth(ix^D)<small_pressure) then
            write(*,*) "Error: small value of gas pressure",pth(ix^D),&
                 " encountered when call mhd_get_pthermal_hde"
            write(*,*) "Iteration: ", it, " Time: ", global_time
            write(*,*) "Location: ", x(ix^D,:)
            write(*,*) "Cell number: ", ix^D
            do iw=1,nw
              write(*,*) trim(cons_wnames(iw)),": ",w(ix^D,iw)
            end do
            if(trace_small_values) write(*,*) sqrt(pth(ix^D)-bigdouble)
            write(*,*) "Saving status at the previous time step"
            crash=.true.
          end if
       {end do\}
      end if

    end subroutine mhd_get_pthermal_hde

    ! ========================================================================
    ! Temperature routines (moved from mod_mhd_phys.t)
    ! ========================================================================

    !> Copy temperature from stored Te variable
    subroutine mhd_get_temperature_from_Te(w, x, ixI^L, ixO^L, res)
      use mod_global_parameters
      integer, intent(in)          :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S, 1:nw)
      double precision, intent(in) :: x(ixI^S, 1:ndim)
      double precision, intent(out):: res(ixI^S)
      res(ixO^S) = w(ixO^S, Te_)
    end subroutine mhd_get_temperature_from_Te

    !> Calculate temperature=p/rho when in e_ the internal energy is stored
    subroutine mhd_get_temperature_from_eint(w, x, ixI^L, ixO^L, res)
      use mod_global_parameters
      integer, intent(in)          :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S, 1:nw)
      double precision, intent(in) :: x(ixI^S, 1:ndim)
      double precision, intent(out):: res(ixI^S)

      double precision :: R(ixI^S)

      call eos%get_Rfactor(w,x,ixI^L,ixO^L,R)
      res(ixO^S) = eos%gamma_minus_1 * w(ixO^S, e_)/(w(ixO^S,rho_)*R(ixO^S))
    end subroutine mhd_get_temperature_from_eint

    !> Calculate temperature=p/rho when in e_ the total energy is stored
    subroutine mhd_get_temperature_from_etot(w, x, ixI^L, ixO^L, res)
      use mod_global_parameters
      integer, intent(in)          :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S, 1:nw)
      double precision, intent(in) :: x(ixI^S, 1:ndim)
      double precision, intent(out):: res(ixI^S)

      double precision :: R(ixI^S)

      call eos%get_Rfactor(w,x,ixI^L,ixO^L,R)
      call mhd_get_pthermal(w,x,ixI^L,ixO^L,res)
      res(ixO^S)=res(ixO^S)/(R(ixO^S)*w(ixO^S,rho_))

    end subroutine mhd_get_temperature_from_etot

    subroutine mhd_get_temperature_from_etot_with_equi(w, x, ixI^L, ixO^L, res)
      use mod_global_parameters
      integer, intent(in)          :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S, 1:nw)
      double precision, intent(in) :: x(ixI^S, 1:ndim)
      double precision, intent(out):: res(ixI^S)

      double precision :: R(ixI^S)

      call eos%get_Rfactor(w,x,ixI^L,ixO^L,R)
      call mhd_get_pthermal(w,x,ixI^L,ixO^L,res)
      res(ixO^S)=res(ixO^S)/(R(ixO^S)*(w(ixO^S,rho_)+block%equi_vars(ixO^S,equi_rho0_,b0i)))

    end subroutine mhd_get_temperature_from_etot_with_equi

    subroutine mhd_get_temperature_from_eint_with_equi(w, x, ixI^L, ixO^L, res)
      use mod_global_parameters
      integer, intent(in)          :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S, 1:nw)
      double precision, intent(in) :: x(ixI^S, 1:ndim)
      double precision, intent(out):: res(ixI^S)

      double precision :: R(ixI^S)

      call eos%get_Rfactor(w,x,ixI^L,ixO^L,R)
      res(ixO^S) = (eos%gamma_minus_1 * w(ixO^S, e_) + block%equi_vars(ixO^S,equi_pe0_,b0i)) /&
                  ((w(ixO^S,rho_) +block%equi_vars(ixO^S,equi_rho0_,b0i))*R(ixO^S))

    end subroutine mhd_get_temperature_from_eint_with_equi

    subroutine mhd_get_temperature_equi(w,x, ixI^L, ixO^L, res)
      use mod_global_parameters
      integer, intent(in)          :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S, 1:nw)
      double precision, intent(in) :: x(ixI^S, 1:ndim)
      double precision, intent(out):: res(ixI^S)

      double precision :: R(ixI^S)

      call eos%get_Rfactor(w,x,ixI^L,ixO^L,R)
      res(ixO^S)= block%equi_vars(ixO^S,equi_pe0_,b0i)/(block%equi_vars(ixO^S,equi_rho0_,b0i)*R(ixO^S))

    end subroutine mhd_get_temperature_equi

    ! ========================================================================
    ! Equilibrium extraction routines (moved from mod_mhd_phys.t)
    ! ========================================================================

    subroutine mhd_get_rho_equi(w, x, ixI^L, ixO^L, res)
      use mod_global_parameters
      integer, intent(in)          :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S, 1:nw)
      double precision, intent(in) :: x(ixI^S, 1:ndim)
      double precision, intent(out):: res(ixI^S)
      res(ixO^S) = block%equi_vars(ixO^S,equi_rho0_,b0i)
    end subroutine mhd_get_rho_equi

    subroutine mhd_get_pe_equi(w,x, ixI^L, ixO^L, res)
      use mod_global_parameters
      integer, intent(in)          :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S, 1:nw)
      double precision, intent(in) :: x(ixI^S, 1:ndim)
      double precision, intent(out):: res(ixI^S)
      res(ixO^S) = block%equi_vars(ixO^S,equi_pe0_,b0i)
    end subroutine mhd_get_pe_equi

    ! ========================================================================
    ! Internal energy extraction + small value handling (moved from mod_mhd_phys.t)
    ! ========================================================================

    !> Get internal energy from conserved state (origin formulation)
    function mhd_get_ei_origin(w, ixI^L, ixO^L) result(ei)
      use mod_global_parameters
      integer, intent(in)          :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S, nw)
      double precision             :: ei(ixO^S)

      ! ei = E_total - e_kinetic - e_magnetic
      ei(ixO^S) = w(ixO^S,e_) - half*((^C&w(ixO^S,m^C_)**2+)/w(ixO^S,rho_) &
           + (^C&w(ixO^S,b^C_)**2+))
    end function mhd_get_ei_origin

    !> Get internal energy from conserved state (internal energy formulation)
    function mhd_get_ei_inte(w, ixI^L, ixO^L) result(ei)
      use mod_global_parameters
      integer, intent(in)          :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S, nw)
      double precision             :: ei(ixO^S)

      ! e_ already stores internal energy
      ei(ixO^S) = w(ixO^S,e_)
    end function mhd_get_ei_inte

    ! mhd_handle_small_ei stays in mod_mhd_phys.t (callers are there)

    !> Update temperature Te_ from pressure using partial ionization
    subroutine mhd_update_temperature(ixI^L,ixO^L,wCT,w,x)
      use mod_global_parameters
      use mod_ionization_degree

      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
      double precision, intent(inout) :: w(ixI^S,1:nw)

      double precision :: iz_H(ixO^S),iz_He(ixO^S), pth(ixI^S)

      call ionization_degree_from_temperature(ixI^L,ixO^L,wCT(ixI^S,Te_),iz_H,iz_He)

      call mhd_get_pthermal(w,x,ixI^L,ixO^L,pth)

      w(ixO^S,Te_)=(2.d0+3.d0*He_abundance)*pth(ixO^S)/(w(ixO^S,rho_)*(1.d0+iz_H(ixO^S)+&
       He_abundance*(iz_He(ixO^S)*(iz_He(ixO^S)+1.d0)+1.d0)))

    end subroutine mhd_update_temperature

end module mod_mhd_eos
!> Needs a line after to pass the preprocesor
