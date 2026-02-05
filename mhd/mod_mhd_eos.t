module mod_mhd_eos
    use mod_global_parameters
    use mod_physics
    use mod_eos
    use mod_eos_container
    use mod_mhd_phys

    use mod_comm_lib, only: mpistop

    implicit none
    private

    public :: mhd_link_eos
contains

    !> Link the appropriate EOS conversion routines based on the selected EoS type
    subroutine mhd_link_eos()

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
        call mpistop('Error: MHD EOS type LTE not yet implemented.')
      else
        call mpistop('Error: Unknown MHD EOS type: ' // trim(eos%eos_type))
      end if

      !> Going to keep this approach for now until fully tested 
      !> - then refactor to remove phys_to_prim/con references
      phys_to_primitive => eos%to_primitive
      phys_to_conserved => eos%to_conserved

    end subroutine mhd_link_eos

    !> Transform primitive variables into conservative ones
    subroutine mhd_to_conserved_origin(ixI^L,ixO^L,w)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      ! double precision, intent(in)    :: x(ixI^S, 1:ndim)

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

    !> Transform primitive variables into conservative ones
    subroutine mhd_to_conserved_origin_noe(ixI^L,ixO^L,w)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      ! double precision, intent(in)    :: x(ixI^S, 1:ndim)

      integer :: ix^D

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! Convert velocity to momentum
        ^C&w(ix^D,m^C_)=w(ix^D,rho_)*w(ix^D,m^C_)\
      {end do\}

    end subroutine mhd_to_conserved_origin_noe

    !> Transform primitive variables into conservative ones
    subroutine mhd_to_conserved_hde(ixI^L,ixO^L,w)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      ! double precision, intent(in)    :: x(ixI^S, 1:ndim)

      integer :: ix^D

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! Calculate total energy from pressure, kinetic and magnetic energy
        w(ix^D,e_)=w(ix^D,p_)*eos%inv_gamma_minus_1&
                  +half*(^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)
        ! Convert velocity to momentum
        ^C&w(ix^D,m^C_)=w(ix^D,rho_)*w(ix^D,m^C_)\
      {end do\}

    end subroutine mhd_to_conserved_hde

    !> Transform primitive variables into conservative ones
    subroutine mhd_to_conserved_inte(ixI^L,ixO^L,w)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      ! double precision, intent(in)    :: x(ixI^S, 1:ndim)

      integer :: ix^D

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! Calculate total energy from pressure, kinetic and magnetic energy
        w(ix^D,e_)=w(ix^D,p_)*eos%inv_gamma_minus_1
        ! Convert velocity to momentum
        ^C&w(ix^D,m^C_)=w(ix^D,rho_)*w(ix^D,m^C_)\
      {end do\}

    end subroutine mhd_to_conserved_inte

    !> Transform primitive variables into conservative ones
    subroutine mhd_to_conserved_split_rho(ixI^L,ixO^L,w)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      ! double precision, intent(in)    :: x(ixI^S, 1:ndim)

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

    !> Transform primitive variables into conservative ones
    subroutine mhd_to_conserved_semirelati(ixI^L,ixO^L,w)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      ! double precision, intent(in)    :: x(ixI^S, 1:ndim)

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

    subroutine mhd_to_conserved_semirelati_noe(ixI^L,ixO^L,w)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      ! double precision, intent(in)    :: x(ixI^S, 1:ndim)

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

    !> Transform conservative variables into primitive ones
    subroutine mhd_to_primitive_origin(ixI^L,ixO^L,w)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      ! double precision, intent(in)    :: x(ixI^S, 1:ndim)

      double precision                :: inv_rho
      integer :: ix^D

      ! if (fix_small_values) then
      !   ! fix small values preventing NaN numbers in the following converting
      !   call mhd_handle_small_values(.false., w, x, ixI^L, ixO^L, 'mhd_to_primitive_origin')
      ! end if

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

    !> Transform conservative variables into primitive ones
    subroutine mhd_to_primitive_origin_noe(ixI^L,ixO^L,w)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      ! double precision, intent(in)    :: x(ixI^S, 1:ndim)

      double precision                :: inv_rho
      integer :: ix^D

      ! if (fix_small_values) then
      !   ! fix small values preventing NaN numbers in the following converting
      !   call mhd_handle_small_values(.false., w, x, ixI^L, ixO^L, 'mhd_to_primitive_origin_noe')
      ! end if

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        inv_rho = 1.d0/w(ix^D,rho_)
        ! Convert momentum to velocity
        ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\
      {end do\}

    end subroutine mhd_to_primitive_origin_noe

    !> Transform conservative variables into primitive ones
    subroutine mhd_to_primitive_hde(ixI^L,ixO^L,w)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      ! double precision, intent(in)    :: x(ixI^S, 1:ndim)

      double precision                :: inv_rho
      integer                         :: ix^D

      ! if (fix_small_values) then
      !   ! fix small values preventing NaN numbers in the following converting
      !   call mhd_handle_small_values(.false., w, x, ixI^L, ixO^L, 'mhd_to_primitive_hde')
      ! end if

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        inv_rho = 1.d0/w(ix^D,rho_)
        ! Convert momentum to velocity
        ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\
        ! Calculate pressure = (gamma-1) * (e-ek)
        w(ix^D,p_)=eos%gamma_minus_1*(w(ix^D,e_)-half*w(ix^D,rho_)*(^C&w(ix^D,m^C_)**2+))
      {end do\}

    end subroutine mhd_to_primitive_hde

    !> Transform conservative variables into primitive ones
    subroutine mhd_to_primitive_inte(ixI^L,ixO^L,w)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      ! double precision, intent(in)    :: x(ixI^S, 1:ndim)

      double precision                :: inv_rho
      integer                         :: ix^D

      ! if (fix_small_values) then
      !   ! fix small values preventing NaN numbers in the following converting
      !   call mhd_handle_small_values(.false., w, x, ixI^L, ixO^L, 'mhd_to_primitive_inte')
      ! end if

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! Calculate pressure = (gamma-1) * e_internal
        w(ix^D,p_)=w(ix^D,e_)*eos%gamma_minus_1
        ! Convert momentum to velocity
        inv_rho = 1.d0/w(ix^D,rho_)
        ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\
      {end do\}

    end subroutine mhd_to_primitive_inte

    !> Transform conservative variables into primitive ones
    subroutine mhd_to_primitive_split_rho(ixI^L,ixO^L,w)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      ! double precision, intent(in)    :: x(ixI^S, 1:ndim)

      double precision :: inv_rho
      integer :: ix^D

      ! if (fix_small_values) then
      !   ! fix small values preventing NaN numbers in the following converting
      !   call mhd_handle_small_values(.false., w, x, ixI^L, ixO^L, 'mhd_to_primitive_split_rho')
      ! end if

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

    !> Transform conservative variables into primitive ones
    subroutine mhd_to_primitive_semirelati(ixI^L,ixO^L,w)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      ! double precision, intent(in)    :: x(ixI^S, 1:ndim)

      double precision :: b(ixO^S,1:ndir), tmp, b2, gamma2, inv_rho
      integer :: ix^D

      ! if (fix_small_values) then
      !   ! fix small values preventing NaN numbers in the following converting
      !   call mhd_handle_small_values(.false., w, x, ixI^L, ixO^L, 'mhd_to_primitive_semirelati')
      ! end if

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

    !> Transform conservative variables into primitive ones
    subroutine mhd_to_primitive_semirelati_noe(ixI^L,ixO^L,w)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(inout) :: w(ixI^S, nw)
      ! double precision, intent(in)    :: x(ixI^S, 1:ndim)

      double precision :: b(ixO^S,1:ndir),tmp,b2,gamma2,inv_rho
      integer :: ix^D, idir

      ! if (fix_small_values) then
      !   ! fix small values preventing NaN numbers in the following converting
      !   call mhd_handle_small_values(.false., w, x, ixI^L, ixO^L, 'mhd_to_primitive_semirelati_noe')
      ! end if

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


end module mod_mhd_eos
!> Needs a line after to pass the preprocesor