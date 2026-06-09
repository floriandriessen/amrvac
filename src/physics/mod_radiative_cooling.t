!> module radiative cooling -- add optically thin radiative cooling
!>
!> only uses the (Townsend) exact integration method, can be used in HD, ffhd, MHD, twofl
!>
!> Cooling curves assume an H/He plasma and ionization equilibrium.
!> The gas temperature is supplied by the active physics EOS; the cooling
!> table composition and abundance assumptions remain curve-dependent.
!> Formula: Q=-n_H*n_e*f(T), where f(T) is tabulated or piecewise power law.
!>
module mod_radiative_cooling

  use mod_global_parameters, only: std_len
  use mod_physics
  use mod_radloss_tables
  use mod_comm_lib, only: mpistop
  implicit none

  !> Helium abundance over Hydrogen
  double precision, private    :: He_abundance

  !> The adiabatic index
  double precision, private :: rc_gamma

  !> The adiabatic index minus 1
  double precision, private :: rc_gamma_1

  !> inverse of the adiabatic index minus 1
  double precision, private :: invgam

  abstract interface
    subroutine get_subr1(w,x,ixI^L,ixO^L,res)
      use mod_global_parameters
      integer, intent(in)          :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S,nw)
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision, intent(out):: res(ixI^S)
    end subroutine get_subr1

    subroutine get_Rfactor_T(T, Rfactor)
      double precision, intent(in)  :: T
      double precision, intent(out) :: Rfactor
    end subroutine get_Rfactor_T
  end interface

  type rc_fluid

    double precision :: rad_damp_height
    double precision :: rad_damp_scale

    ! these are set in init method
    double precision, allocatable :: tcool(:), Lcool(:), dLdtcool(:)
    double precision, allocatable :: Yc(:)
    double precision  :: tref, lref, tcoolmin,tcoolmax
    double precision  :: lgtcoolmin, lgtcoolmax, lgstep

    ! The piecewise powerlaw (PPL) tabels and variabels
    ! x_* en t_* are given as log_10
    double precision, allocatable :: y_PPL(:), t_PPL(:), l_PPL(:), a_PPL(:)

    !> Lower limit of temperature
    double precision   :: tlow

    !> Index of the energy density
    integer              :: e_
    !> Index of cut off temperature for TRAC
    integer              :: Tcoff_

    ! these are set as parameters
    !> Resolution of temperature in interpolated tables
    integer :: ncool

    integer :: n_PPL

    !> Fixed temperature not lower than tlow
    logical   :: Tfix

    !> Add cooling source in a split way (.true.) or un-split way (.false.)
    logical    :: rc_split

    logical :: isPPL = .false.

    !> Apply radiative damping near both x1 boundaries for 1D loop models
    logical :: rc_is_1d_loop = .false.
    !> cutoff radiative cooling below rad_damp_height
    logical :: rad_damp
    !> whether background equilibrium contribution is split off
    logical :: has_equi = .false.
    !> whether background equilibrium is compensated in thermal balance
    logical :: subtract_equi = .false.

    double precision, allocatable :: frac_lowFIP(:)
    !> Index of primitive FIP abundance variable, -1 if disabled
    integer :: fip_ = -1
    !> Enable local Newton cooling/heating approximation for optically thick losses
    logical :: rad_newton = .false.
    double precision :: rad_newton_pthick = 25.d0
    double precision :: rad_newton_trad = 0.006d0
    double precision :: rad_newton_rhosurf = 1.d4

    !> Name of cooling curve
    character(len=std_len)  :: coolcurve

    procedure (get_subr1), pointer, nopass :: get_rho => null()
    procedure (get_subr1), pointer, nopass :: get_rho_equi => null()
    procedure (get_subr1), pointer, nopass :: get_pthermal => null()
    procedure (get_subr1), pointer, nopass :: get_temperature => null()
    procedure (get_subr1), pointer, nopass :: get_pthermal_equi => null()
    procedure (get_subr1), pointer, nopass :: get_var_Rfactor => null()
    procedure (get_Rfactor_T), pointer, nopass :: get_Rfactor_from_temperature => null()
    procedure (get_subr1), pointer, nopass :: get_temperature_equi => null()

  end type rc_fluid

  contains

    !> Radiative cooling initialization
    subroutine radiative_cooling_init_params(phys_gamma,He_abund)
      use mod_global_parameters
      double precision, intent(in) :: phys_gamma,He_abund

      rc_gamma=phys_gamma
      He_abundance=He_abund
    end subroutine radiative_cooling_init_params

    subroutine radiative_cooling_init(fl,read_params)
      use mod_global_parameters
      interface
        subroutine read_params(fl)
          use mod_global_parameters, only: unitpar,par_files
          import rc_fluid
          type(rc_fluid), intent(inout) :: fl

        end subroutine read_params
      end interface

      type(rc_fluid), intent(inout) :: fl

      double precision, dimension(:), allocatable :: t_table
      double precision, dimension(:), allocatable :: L_table
      double precision, dimension(:), allocatable :: f_table
      double precision :: ratt, fact1, fact2, fact3, dL1, dL2
      double precision :: tstep, Lstep
      integer :: ntable, i, j
      logical :: jump
      Character(len=65) :: PPL_curves(1:6)

      fl%ncool=4000
      fl%coolcurve='JCcorona'
      fl%tlow=bigdouble
      fl%Tfix=.false.
      fl%rc_split=.false.
      fl%rc_is_1d_loop = .false.
      fl%rad_damp=.false.
      fl%rad_damp_height=0.5d0
      fl%rad_damp_scale=0.15d0
      call read_params(fl)

      if (fl%fip_ > 0) then
        select case (trim(fl%coolcurve))
        case ('Dere_photo', 'Dere_photo_DM')
        case default
          call mpistop("FIP cooling requires coolcurve='Dere_photo' or 'Dere_photo_DM'")
        end select
      end if

      if(fl%rc_split) any_source_split=.true.

      ! Checks if coolcurve is a piecewise power law (PPL)
      PPL_curves = [Character(len=65) :: 'Hildner','FM', 'Rosner', 'Klimchuk','SPEX_DM_rough','SPEX_DM_fine']
      do i=1,size(PPL_curves)
         if (PPL_curves(i)==fl%coolcurve) then
            fl%isPPL = .true.
         end if
      end do

      ! Init for PPL
      if (fl%isPPL) then
         ! Read in tables and create t_PPL, l_PPL, a_PPL
         select case(fl%coolcurve)

         case('Hildner')
            if(mype ==0) &
            print *,'Use Hildner (1974) piecewise power law'
            fl%n_PPL = n_Hildner
            allocate(fl%t_PPL(1:fl%n_PPL+1), fl%l_PPL(1:fl%n_PPL+1))
            allocate(fl%a_PPL(1:fl%n_PPL))
            fl%t_PPL(1:fl%n_PPL+1) = t_Hildner(1:n_Hildner+1)
            fl%a_PPL(1:fl%n_PPL) = a_Hildner(1:n_Hildner)
            fl%l_PPL(1:fl%n_PPL) = 10.d0**x_Hildner(1:n_Hildner) * (10.d0**fl%t_PPL(1:fl%n_PPL))**fl%a_PPL(1:fl%n_PPL)

         case('FM')
            if(mype==0) &
            print *,'Use Forbes and Malherbe (1991)-like piecewise power law'
            fl%n_PPL = n_FM
            allocate(fl%t_PPL(1:fl%n_PPL+1), fl%l_PPL(1:fl%n_PPL+1))
            allocate(fl%a_PPL(1:fl%n_PPL))
            fl%t_PPL(1:fl%n_PPL+1) = t_FM(1:n_FM+1)
            fl%a_PPL(1:fl%n_PPL) = a_FM(1:n_FM)
            fl%l_PPL(1:fl%n_PPL) = 10.d0**x_FM(1:n_FM) * (10.d0**fl%t_PPL(1:fl%n_PPL))**fl%a_PPL(1:fl%n_PPL)

         case('Rosner')
            if(mype==0) &
            print *,'Use piecewise power law according to Rosner (1978)'
            if(mype ==0) &
            print *,'and extended by Priest (1982) from Van Der Linden (1991)'
            fl%n_PPL = n_Rosner
            allocate(fl%t_PPL(1:fl%n_PPL+1), fl%l_PPL(1:fl%n_PPL+1))
            allocate(fl%a_PPL(1:fl%n_PPL))
            fl%t_PPL(1:fl%n_PPL+1) = t_Rosner(1:n_Rosner+1)
            fl%a_PPL(1:fl%n_PPL) = a_Rosner(1:n_Rosner)
            fl%l_PPL(1:fl%n_PPL) = 10.d0**x_Rosner(1:n_Rosner) * (10.d0**fl%t_PPL(1:fl%n_PPL))**fl%a_PPL(1:fl%n_PPL)

         case('Klimchuk')
            if(mype==0) &
            print *,'Use Klimchuk (2008) piecewise power law'
            fl%n_PPL = n_Klimchuk
            allocate(fl%t_PPL(1:fl%n_PPL+1), fl%l_PPL(1:fl%n_PPL+1))
            allocate(fl%a_PPL(1:fl%n_PPL))
            fl%t_PPL(1:fl%n_PPL+1) = t_Klimchuk(1:n_Klimchuk+1)
            fl%a_PPL(1:fl%n_PPL) = a_Klimchuk(1:n_Klimchuk)
            fl%l_PPL(1:fl%n_PPL) = 10.d0**x_Klimchuk(1:n_Klimchuk) * (10.d0**fl%t_PPL(1:fl%n_PPL))**fl%a_PPL(1:fl%n_PPL)

         case('SPEX_DM_rough')
            if(mype==0) &
            print *,'Use the rough piece wise power law fit to the SPEX_DM curve (2009)'
            fl%n_PPL = n_SPEX_DM_rough
            allocate(fl%t_PPL(1:fl%n_PPL+1), fl%l_PPL(1:fl%n_PPL+1))
            allocate(fl%a_PPL(1:fl%n_PPL))
            fl%t_PPL(1:fl%n_PPL+1) = t_SPEX_DM_rough(1:n_SPEX_DM_rough+1)
            fl%a_PPL(1:fl%n_PPL) = a_SPEX_DM_rough(1:n_SPEX_DM_rough)
            fl%l_PPL(1:fl%n_PPL) = 10.d0**x_SPEX_DM_rough(1:n_SPEX_DM_rough) * (10.d0**fl%t_PPL(1:fl%n_PPL))**fl%a_PPL(1:fl%n_PPL)

         case('SPEX_DM_fine')
            if(mype==0) &
            print *,'Use the fine, detailed piece wise power law fit to the SPEX_DM curve (2009)'
            fl%n_PPL = n_SPEX_DM_fine
            allocate(fl%t_PPL(1:fl%n_PPL+1), fl%l_PPL(1:fl%n_PPL+1))
            allocate(fl%a_PPL(1:fl%n_PPL))
            fl%t_PPL(1:fl%n_PPL+1) = t_SPEX_DM_fine(1:n_SPEX_DM_fine+1)
            fl%a_PPL(1:fl%n_PPL) = a_SPEX_DM_fine(1:n_SPEX_DM_fine)
            fl%l_PPL(1:fl%n_PPL) = 10.d0**x_SPEX_DM_fine(1:n_SPEX_DM_fine) * (10.d0**fl%t_PPL(1:fl%n_PPL))**fl%a_PPL(1:fl%n_PPL)

         case default
            call mpistop("This piecewise power law is unknown")
         end select

         ! Go from logarithmic to actual values.
         fl%t_PPL(1:fl%n_PPL+1) = 10.d0**fl%t_PPL(1:fl%n_PPL+1)
         ! Change unit of table if SI is used instead of cgs
         if (SI_unit) fl%l_PPL(1:fl%n_PPL) = fl%l_PPL(1:fl%n_PPL) * 10.0d0**(-13)

         ! Make dimensionless
         fl%t_PPL(1:fl%n_PPL+1) = fl%t_PPL(1:fl%n_PPL+1) / unit_temperature
         fl%l_PPL(1:fl%n_PPL) = fl%l_PPL(1:fl%n_PPL) * unit_numberdensity**2 * unit_time / unit_pressure * (1.d0+2.d0*He_abundance)

         ! Set tref en lref
         fl%l_PPL(fl%n_PPL+1) = fl%l_PPL(fl%n_PPL) * ( fl%t_PPL(fl%n_PPL+1) / fl%t_PPL(fl%n_PPL) )**fl%a_PPL(fl%n_PPL)
         fl%lref = fl%l_PPL(fl%n_PPL+1)
         fl%tref = fl%t_PPL(fl%n_PPL+1)

         ! Set tcoolmin and tcoolmax
         fl%tcoolmin = fl%t_PPL(1)
         fl%tcoolmax = fl%t_PPL(fl%n_PPL+1)
         ! smaller value for lowest temperatures from cooling table and user's choice
         if (fl%tlow==bigdouble) fl%tlow=fl%tcoolmin
         !create y_PPL
         call create_y_PPL(fl)

      else

         ! Init for interpolatable tables
         allocate(fl%tcool(1:fl%ncool), fl%Lcool(1:fl%ncool), fl%dLdtcool(1:fl%ncool))
         allocate(fl%Yc(1:fl%ncool))
         if(fl%fip_ > 0) allocate(fl%frac_lowFIP(1:fl%ncool))

         fl%tcool(1:fl%ncool)    = zero
         fl%Lcool(1:fl%ncool)    = zero
         fl%dLdtcool(1:fl%ncool) = zero

         ! Read in the selected cooling curve
         select case(fl%coolcurve)

         case('JCcorona')
            if(mype ==0) &
            print *,'Use Colgan & Feldman (2008) cooling curve'
            if(mype ==0) &
            print *,'This version only till 10000 K, beware for floor T treatment'
            ntable = n_JCcorona
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            t_table(1:ntable) = t_JCcorona(1:n_JCcorona)
            L_table(1:ntable) = l_JCcorona(1:n_JCcorona)

         case('DM')
            if(mype ==0) &
            print *,'Use Dalgarno & McCray (1972) cooling curve'
            ntable = n_DM
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            t_table(1:ntable) = t_DM(1:n_DM)
            L_table(1:ntable) = l_DM(1:n_DM)

         case('MB')
            if(mype ==0) &
            write(*,'(3a)') 'Use MacDonald & Bailey (1981) cooling curve '&
                 ,'as implemented in ZEUS-3D, with the values '&
                 ,'from Dalgarno & McCRay (1972) for low temperatures.'
            ntable = n_MB + 20
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            t_table(1:ntable) = t_DM(1:21)
            L_table(1:ntable) = l_DM(1:21)
            t_table(22:ntable) = t_MB(2:n_MB)
            L_table(22:ntable) = l_MB(2:n_MB)

         case('MLcosmol')
            if(mype ==0) &
            print *,'Use Mellema & Lundqvist (2002) cooling curve '&
                 ,'for zero metallicity '
            ntable = n_MLcosmol
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            t_table(1:ntable) = t_MLcosmol(1:n_MLcosmol)
            L_table(1:ntable) = l_MLcosmol(1:n_MLcosmol)

         case('MLwc')
            if(mype ==0) &
            print *,'Use Mellema & Lundqvist (2002) cooling curve '&
                 ,'for WC-star metallicity '
            ntable = n_MLwc
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            t_table(1:ntable) = t_MLwc(1:n_MLwc)
            L_table(1:ntable) = l_MLwc(1:n_MLwc)

         case('MLsolar1')
            if(mype ==0) &
            print *,'Use Mellema & Lundqvist (2002) cooling curve '&
                 ,'for solar metallicity '
            ntable = n_MLsolar1
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            t_table(1:ntable) = t_MLsolar1(1:n_MLsolar1)
            L_table(1:ntable) = l_MLsolar1(1:n_MLsolar1)

         case('cloudy_ism')
            if(mype ==0) &
            print *,'Use Cloudy based cooling curve '&
                 ,'for ism metallicity '
            ntable = n_cl_ism
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            t_table(1:ntable) = t_cl_ism(1:n_cl_ism)
            L_table(1:ntable) = l_cl_ism(1:n_cl_ism)

         case('cloudy_solar')
            if(mype ==0) &
            print *,'Use Cloudy based cooling curve '&
                 ,'for solar metallicity '
            ntable = n_cl_solar
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            t_table(1:ntable) = t_cl_solar(1:n_cl_solar)
            L_table(1:ntable) = l_cl_solar(1:n_cl_solar)

         case('SPEX')
            if(mype ==0) &
            print *,'Use SPEX cooling curve (Schure et al. 2009) '&
                 ,'for solar metallicity '
            ntable = n_SPEX
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            t_table(1:ntable) = t_SPEX(1:n_SPEX)
            L_table(1:ntable) = l_SPEX(1:n_SPEX) + log10(nenh_SPEX(1:n_SPEX))

         case('SPEX_DM')
            if(mype ==0) then
               print *, 'Use SPEX cooling curve for solar metallicity above 10^4 K. '
               print *, 'At lower temperatures,use Dalgarno & McCray (1972), '
               print *, 'with a pre-set ionization fraction of 10^-3. '
               print *, 'as described by Schure et al. (2009). '
            endif
            ntable = n_SPEX + n_DM_2 - 6
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            t_table(1:n_DM_2-1) = t_DM_2(1:n_DM_2-1)
            L_table(1:n_DM_2-1) = L_DM_2(1:n_DM_2-1)
            t_table(n_DM_2:ntable) = t_SPEX(6:n_SPEX)
            L_table(n_DM_2:ntable) = l_SPEX(6:n_SPEX) + log10(nenh_SPEX(6:n_SPEX))

         case('Dere_corona')
            if(mype ==0) &
            print *,'Use Dere (2009) cooling curve for solar corona'
            ntable = n_Dere
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            t_table(1:ntable) = t_Dere(1:n_Dere)
            L_table(1:ntable) = l_Dere_corona(1:n_Dere)

         case('Dere_corona_DM')
            if(mype==0)&
            print *, 'Combination of Dere_corona (2009) for high temperatures and'
            if(mype==0)&
            print *, 'Dalgarno & McCray (1972), DM2, for low temperatures'
            ntable = n_Dere + n_DM_2 - 1
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            t_table(1:n_DM_2-1) = t_DM_2(1:n_DM_2-1)
            L_table(1:n_DM_2-1) = L_DM_2(1:n_DM_2-1)
            t_table(n_DM_2:ntable) = t_Dere(1:n_Dere)
            L_table(n_DM_2:ntable) = l_Dere_corona(1:n_Dere)

         case('Dere_photo')
            if(mype ==0) &
            print *,'Use Dere (2009) cooling curve for solar photophere'
            ntable = n_Dere
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            if (fl%fip_ > 0) allocate(f_table(1:ntable))
            t_table(1:ntable) = t_Dere(1:n_Dere)
            L_table(1:ntable) = l_Dere_photo(1:n_Dere)
            if (fl%fip_ > 0) f_table(1:ntable) = lowFIP_frac(1:n_Dere)

         case('Dere_photo_DM')
            if(mype==0)&
            print *, 'Combination of Dere_photo (2009) for high temperatures and'
            if(mype==0)&
            print *, 'Dalgarno & McCray (1972), DM2, for low temperatures'
            ntable = n_Dere + n_DM_2 - 1
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            if (fl%fip_ > 0) allocate(f_table(1:ntable))
            t_table(1:n_DM_2-1) = t_DM_2(1:n_DM_2-1)
            L_table(1:n_DM_2-1) = L_DM_2(1:n_DM_2-1)
            t_table(n_DM_2:ntable) = t_Dere(1:n_Dere)
            L_table(n_DM_2:ntable) = l_Dere_photo(1:n_Dere)
            if (fl%fip_ > 0) then
              f_table(1:n_DM_2-1) = zero
              f_table(n_DM_2:ntable) = lowFIP_frac(1:n_Dere)
            end if

         case('Colgan')
            if(mype==0) &
            print *, 'Use Colgan (2008) cooling curve'
            ntable = n_Colgan
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            t_table(1:ntable) = t_Colgan(1:n_Colgan)
            L_table(1:ntable) = l_Colgan(1:n_Colgan)

         case('Colgan_DM')
            if(mype==0)&
            print *, 'Combination of Colgan (2008) for high temperatures and'
            if(mype==0)&
            print *, 'Dalgarno & McCray (1972), DM2, for low temperatures'
            ntable = n_Colgan + n_DM_2
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            t_table(1:n_DM_2) = t_DM_2(1:n_DM_2)
            L_table(1:n_DM_2) = L_DM_2(1:n_DM_2)
            t_table(n_DM_2+1:ntable) = t_Colgan(1:n_Colgan)
            L_table(n_DM_2+1:ntable) = l_Colgan(1:n_Colgan)

         case default
            call mpistop("This coolingcurve is unknown")
         end select

         ! create cooling table(s) for use in amrvac
         fl%tcoolmax = t_table(ntable)
         fl%tcoolmin = t_table(1)
         ratt = (fl%tcoolmax-fl%tcoolmin)/( dble(fl%ncool-1) + smalldouble)

         fl%tcool(1) = fl%tcoolmin
         fl%Lcool(1) = L_table(1)

         fl%tcool(fl%ncool) = fl%tcoolmax
         fl%Lcool(fl%ncool) = L_table(ntable)

         if (fl%fip_ > 0) then
           fl%frac_lowFIP(1) = f_table(1)
           fl%frac_lowFIP(fl%ncool) = f_table(ntable)
         end if

         do i=2,fl%ncool        ! loop to create one table
           fl%tcool(i) = fl%tcool(i-1)+ratt
           do j=1,ntable-1   ! loop to create one spot on a table
           ! Second order polynomial interpolation, except at the outer edge,
           ! or in case of a large jump.
             if(fl%tcool(i) < t_table(j+1)) then
                if(j.eq. ntable-1 )then
                  fact1 = (fl%tcool(i)-t_table(j+1))     &
                        /(t_table(j)-t_table(j+1))
                  fact2 = (fl%tcool(i)-t_table(j))       &
                        /(t_table(j+1)-t_table(j))
                  fl%Lcool(i) = L_table(j)*fact1 + L_table(j+1)*fact2
                  if (fl%fip_ > 0) then
                    fl%frac_lowFIP(i) = f_table(j)*fact1 + f_table(j+1)*fact2
                  end if
                  exit
                else
                  dL1 = L_table(j+1)-L_table(j)
                  dL2 = L_table(j+2)-L_table(j+1)
                  jump =(max(dabs(dL1),dabs(dL2)) > 2*min(dabs(dL1),dabs(dL2)))
                end if
                if( jump ) then
                  fact1 = (fl%tcool(i)-t_table(j+1))     &
                        /(t_table(j)-t_table(j+1))
                  fact2 = (fl%tcool(i)-t_table(j))       &
                        /(t_table(j+1)-t_table(j))
                  fl%Lcool(i) = L_table(j)*fact1 + L_table(j+1)*fact2
                  if (fl%fip_ > 0) then
                    fl%frac_lowFIP(i) = f_table(j)*fact1 + f_table(j+1)*fact2
                  end if
                  exit
                else
                  fact1 = ((fl%tcool(i)-t_table(j+1))     &
                        * (fl%tcool(i)-t_table(j+2)))   &
                        / ((t_table(j)-t_table(j+1)) &
                        * (t_table(j)-t_table(j+2)))
                  fact2 = ((fl%tcool(i)-t_table(j))       &
                        * (fl%tcool(i)-t_table(j+2)))   &
                        / ((t_table(j+1)-t_table(j)) &
                        * (t_table(j+1)-t_table(j+2)))
                  fact3 = ((fl%tcool(i)-t_table(j))       &
                        * (fl%tcool(i)-t_table(j+1)))   &
                        / ((t_table(j+2)-t_table(j)) &
                        * (t_table(j+2)-t_table(j+1)))
                  fl%Lcool(i) = L_table(j)*fact1 + L_table(j+1)*fact2 &
                           + L_table(j+2)*fact3
                  if (fl%fip_ > 0) then
                    fl%frac_lowFIP(i) = f_table(j)*fact1 + f_table(j+1)*fact2 &
                                      + f_table(j+2)*fact3
                  end if
                  exit
                end if
             end if
           end do  ! end loop to find create one spot on a table
         end do    ! end loop to create one table

         ! Go from logarithmic to actual values.
         fl%tcool(1:fl%ncool) = 10.0D0**fl%tcool(1:fl%ncool)
         fl%Lcool(1:fl%ncool) = 10.0D0**fl%Lcool(1:fl%ncool)

         ! Change unit of table if SI is used instead of cgs
         if (SI_unit) fl%Lcool(1:fl%ncool) = fl%Lcool(1:fl%ncool) * 10.0d0**(-13)

         ! Scale both T and Lambda
         fl%tcool(1:fl%ncool) = fl%tcool(1:fl%ncool) / unit_temperature
         fl%Lcool(1:fl%ncool) = fl%Lcool(1:fl%ncool) * unit_numberdensity**2 * unit_time / unit_pressure * (1.d0+2.d0*He_abundance)

         fl%tcoolmin       = fl%tcool(1)+smalldouble  ! avoid pointless interpolation
         ! smaller value for lowest temperatures from cooling table and user's choice
         if (fl%tlow==bigdouble) fl%tlow=fl%tcoolmin
         fl%tcoolmax       = fl%tcool(fl%ncool)
         fl%lgtcoolmin = dlog10(fl%tcoolmin)
         fl%lgtcoolmax = dlog10(fl%tcoolmax)
         fl%lgstep = (fl%lgtcoolmax-fl%lgtcoolmin) * 1.d0 / (fl%ncool-1)
         fl%dLdtcool(1)     = (fl%Lcool(2)-fl%Lcool(1))/(fl%tcool(2)-fl%tcool(1))
         fl%dLdtcool(fl%ncool) = (fl%Lcool(fl%ncool)-fl%Lcool(fl%ncool-1))/(fl%tcool(fl%ncool)-fl%tcool(fl%ncool-1))

         do i=2,fl%ncool-1
           fl%dLdtcool(i) = (fl%Lcool(i+1)-fl%Lcool(i-1))/(fl%tcool(i+1)-fl%tcool(i-1))
         end do

         deallocate(t_table)
         deallocate(L_table)
         if (allocated(f_table)) deallocate(f_table)

         fl%tref = fl%tcoolmax
         fl%lref = fl%Lcool(fl%ncool)
         fl%Yc(fl%ncool) = zero
         do i=fl%ncool-1, 1, -1
            fl%Yc(i) = fl%Yc(i+1)
            do j=1,100
               tstep = 1.0d-2*(fl%tcool(i+1)-fl%tcool(i))
               call findL(fl%tcool(i+1)-j*tstep, Lstep, fl)
               fl%Yc(i) = fl%Yc(i) + fl%lref/fl%tref*tstep/Lstep
            end do
         end do
      end if

      rc_gamma_1=rc_gamma-1.d0
      invgam = 1.d0/rc_gamma_1

    end subroutine radiative_cooling_init

    subroutine create_y_PPL(fl)
    !  creates the constants of integration needed for solving
    !  the cooling law exact for a piecewise power law
    !  In correspondence with eq. A6 of Townsend (2009)
      use mod_global_parameters
      type(rc_fluid) :: fl
      double precision :: y_extra, factor
      integer :: i

      allocate(fl%y_PPL(1:fl%n_PPL+1))

      fl%y_PPL(1:fl%n_PPL+1) = zero

      do i=fl%n_PPL, 1, -1
          factor = fl%l_PPL(fl%n_PPL+1) * fl%t_PPL(i) / (fl%l_PPL(i) * fl%t_PPL(fl%n_PPL+1))
          if (fl%a_PPL(i) == 1.d0) then
             y_extra =  log( fl%t_PPL(i) / fl%t_PPL(i+1) )
          else
             y_extra = 1 / (1 - fl%a_PPL(i)) * (1 - ( fl%t_PPL(i) / fl%t_PPL(i+1) )**(fl%a_PPL(i)-1) )
          end if
          fl%y_PPL(i) = fl%y_PPL(i+1) - factor*y_extra
      end do
    end subroutine create_y_PPL

    subroutine getvar_cooling(ixI^L,ixO^L,w,x,coolrate,fl)
    ! Create extra variable to show cooling rate in the output.
    ! This diagnostic returns the effective instantaneous optically-thin
    ! radiative cooling rate.
    ! Newton cooling is intentionally excluded here because it is not an
    ! optically-thin radiative loss term.
      use mod_global_parameters

      integer, intent(in)          :: ixI^L,ixO^L
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision             :: w(ixI^S,1:nw)
      double precision, intent(out):: coolrate(ixI^S)
      type(rc_fluid), intent(in)   :: fl

      double precision :: rho(ixI^S)
      double precision :: L1, Te(ixI^S)
      double precision :: rho_safe
      double precision :: fip_prim, frac_lowFIP, fip_factor
      double precision :: rad_damp_factor
      integer :: ix^D

      call fl%get_rho(w,x,ixI^L,ixO^L,rho)
      call fl%get_temperature(w,x,ixI^L,ixO^L,Te)

      {do ix^DB = ixO^LIM^DB\}
        ! Ordinary optically-thin cooling curve contribution: rho^2 Lambda(T)
        if(Te(ix^D) <= fl%tcoolmin) then
          L1 = zero
        else if(Te(ix^D) >= fl%tcoolmax) then
          call calc_l_extended(Te(ix^D),L1,fl)
          L1 = L1*rho(ix^D)**2
        else
          call findL(Te(ix^D),L1,fl)
          L1 = L1*rho(ix^D)**2
        end if

        ! FIP-dependent correction.
        ! In this routine w is conserved, so the tracer is rho*fip rather than fip.
        fip_factor = one
        if (fl%fip_ > 0) then
          rho_safe = max(rho(ix^D), small_density)
          fip_prim = min(maxfip, max(minfip, w(ix^D,fl%fip_) / rho_safe))
          frac_lowFIP = lowFIP_fraction(Te(ix^D), fl)
          fip_factor = one - frac_lowFIP + fip_prim * frac_lowFIP
        end if

        ! Geometric damping of the optically-thin cooling, matching cool_exact.
        rad_damp_factor = one
        {^IFONED
        if (slab_uniform .and. fl%rad_damp .and. x(ix^D,ndim) <= xprobmin1 + fl%rad_damp_height) then
          rad_damp_factor = exp(-(x(ix^D,ndim)-xprobmin1-fl%rad_damp_height)**2/fl%rad_damp_scale**2)
        end if
        if (fl%rc_is_1d_loop .and. slab_uniform .and. fl%rad_damp &
               .and. x(ix^D,ndim) >= xprobmax1 - fl%rad_damp_height) then
          rad_damp_factor = exp(-(x(ix^D,ndim)-xprobmax1+fl%rad_damp_height)**2/fl%rad_damp_scale**2)
        end if
        }
        {^IFTWOD
        if (slab_uniform .and. fl%rad_damp .and. x(ix^D,ndim) <= xprobmin2 + fl%rad_damp_height) then
          rad_damp_factor = exp(-(x(ix^D,ndim)-xprobmin2-fl%rad_damp_height)**2/fl%rad_damp_scale**2)
        end if
        }
        {^IFTHREED
        if (slab_uniform .and. fl%rad_damp .and. x(ix^D,ndim) <= xprobmin3 + fl%rad_damp_height) then
          rad_damp_factor = exp(-(x(ix^D,ndim)-xprobmin3-fl%rad_damp_height)**2/fl%rad_damp_scale**2)
        end if
        }

        coolrate(ix^D) = L1 * fip_factor * rad_damp_factor
      {end do\}
    end subroutine getvar_cooling

    subroutine getvar_cooling_exact(qdt, ixI^L, ixO^L, wCT, w, x, coolrate, fl)
    ! TODO: getvar_cooling_exact is a diagnostic routine and is not currently
    ! guaranteed to reproduce cool_exact after the FIP-dependent cooling,
    ! rad_damp, and Newton-cooling extensions.
      use mod_global_parameters

      integer, intent(in)           :: ixI^L, ixO^L
      double precision, intent(in)  :: qdt, x(ixI^S, 1:ndim), wCT(ixI^S, 1:nw)
      double precision              :: w(ixI^S, 1:nw)
      double precision, intent(out) :: coolrate(ixI^S)
      type(rc_fluid), intent(in)    :: fl
      double precision              :: y1, y2, l1, tlocal2
      double precision              :: Te(ixI^S), pnew(ixI^S), rho(ixI^S), rhonew(ixI^S)
      double precision              :: emin, Lmax, fact, Rfactor(ixI^S)
      double precision              :: Rfloor, Rnew
      integer                       :: ix^D
      logical :: use_Rfactor_T
      double precision :: Rfloor_const

      call fl%get_rho(wCT, x, ixI^L, ixO^L, rho)
      call fl%get_temperature(wCT, x, ixI^L, ixO^L, Te)
      call fl%get_pthermal(w, x, ixI^L, ixO^L, pnew)
      call fl%get_rho(w, x, ixI^L, ixO^L, rhonew)

      fact=fl%lref*qdt/fl%tref
      use_Rfactor_T = associated(fl%get_Rfactor_from_temperature)
      if (use_Rfactor_T) then
        {do ix^DB = ixO^LIM^DB\}
          call fl%get_Rfactor_from_temperature(Te(ix^D), Rfactor(ix^D))
        {end do\}
      else
        call fl%get_var_Rfactor(wCT, x, ixI^L, ixO^L, Rfactor)
      end if
      if (use_Rfactor_T) then
        call cooling_get_Rfactor_T(fl, fl%tlow, 1.d0, Rfloor_const)
      end if
      {do ix^DB = ixO^LIM^DB\}
         if (use_Rfactor_T) then
           Rfloor = Rfloor_const
         else
           Rfloor = Rfactor(ix^D)
         end if
         emin = rhonew(ix^D) * fl%tlow * Rfloor * invgam
         lmax = max(zero, ( pnew(ix^D)*invgam - emin ) / qdt)

         ! No cooling if temperature is below floor level.
         ! Assuming Bremsstrahlung if temperature is higher than maximum.
         if( Te(ix^D)<= fl%tcoolmin) then
           l1 = zero
         else if( Te(ix^D)>= fl%tcoolmax ) then
           call calc_l_extended(Te(ix^D), l1, fl)
           l1 = l1 * rho(ix^D)**2
           l1 = min(l1, lmax)
         else
           call findY(Te(ix^D), y1, fl)
           y2   = y1 +  fact * rho(ix^D)*rc_gamma_1
           call findT(tlocal2, y2, fl)
           if( tlocal2 <= fl%tcoolmin ) then
             l1 = lmax
           else
             call cooling_get_Rfactor_T( &
                fl, tlocal2, Rfactor(ix^D), Rnew)
             l1 = rho(ix^D) * (Rfactor(ix^D)*Te(ix^D) - Rnew*tlocal2) * invgam/qdt
             l1 = max(zero, l1)
           end if
           l1 = min(l1, lmax)
         end if
         if(slab_uniform .and. fl%rad_damp .and. x(ix^D,ndim) .le. fl%rad_damp_height) then
           l1 = l1*exp(-(x(ix^D,ndim)-fl%rad_damp_height)**2/fl%rad_damp_scale**2)
         end if
        coolrate(ix^D) = l1
      {end do\}
    end subroutine getvar_cooling_exact

    subroutine radiative_cooling_add_source(qdt,ixI^L,ixO^L,wCT,wCTprim,w,x,&
         qsourcesplit,active,fl)
    ! w[iw]=w[iw]+qdt*S[wCT,x] where S is the source based on wCT within ixO
      use mod_global_parameters
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: qdt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw), wCTprim(ixI^S,1:nw)
      double precision, intent(inout) :: w(ixI^S,1:nw)
      logical, intent(in) :: qsourcesplit
      logical, intent(inout) :: active
      type(rc_fluid), intent(in) :: fl
      double precision, allocatable, dimension(:^D&) :: Lequi

      if(qsourcesplit .eqv.fl%rc_split) then
       active = .true.
       call cool_exact(qdt,ixI^L,ixO^L,wCT,wCTprim,w,x,fl)
       if(fl%subtract_equi) then
          allocate(Lequi(ixI^S))
          call get_cool_equi(qdt,ixI^L,ixO^L,wCT,w,x,fl,Lequi)
          w(ixO^S,fl%e_) = w(ixO^S,fl%e_)+Lequi(ixO^S)
          deallocate(Lequi)
       endif
       if( fl%Tfix ) call floortemperature(qdt,ixI^L,ixO^L,wCT,w,x,fl)
      end if
    end subroutine radiative_cooling_add_source

    subroutine floortemperature(qdt,ixI^L,ixO^L,wCT,w,x,fl)
    !  Force minimum temperature to a fixed temperature
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(in)    :: qdt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
      double precision, intent(inout) :: w(ixI^S,1:nw)
      type(rc_fluid), intent(in) :: fl
      double precision :: etherm(ixI^S), rho(ixI^S)
      double precision :: Rfactor(ixI^S), emin, Rfloor
      integer :: ix^D
      logical :: use_Rfactor_T
      double precision :: Rfloor_const

      call fl%get_pthermal(w,x,ixI^L,ixO^L,etherm)
      call fl%get_rho(w,x,ixI^L,ixO^L,rho)
      call fl%get_var_Rfactor(wCT,x,ixI^L,ixO^L,Rfactor)
      use_Rfactor_T = associated(fl%get_Rfactor_from_temperature)
      if (use_Rfactor_T) then
        call fl%get_Rfactor_from_temperature(fl%tlow, Rfloor_const)
      end if
      {do ix^DB = ixO^LIM^DB\}
         if (use_Rfactor_T) then
           Rfloor = Rfloor_const
         else
           Rfloor = Rfactor(ix^D)
         end if
         emin = rho(ix^D)*fl%tlow*Rfloor
         if (etherm(ix^D) < emin) then
           w(ix^D,fl%e_) = w(ix^D,fl%e_) + &
                (emin-etherm(ix^D))*invgam
         end if
      {end do\}
    end subroutine floortemperature

    subroutine get_cool_equi(qdt,ixI^L,ixO^L,wCT,w,x,fl,res)
      use mod_global_parameters

      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(in)    :: qdt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
      double precision, intent(inout) :: w(ixI^S,1:nw)
      type(rc_fluid), intent(in) :: fl
      double precision, intent(out) :: res(ixI^S)

      double precision :: pth(ixI^S),rho(ixI^S),Rfactor(ixI^S),L1,Tlocal2
      double precision :: Te(ixI^S)
      double precision :: emin, Lmax
      double precision :: Y1, Y2
      double precision :: de, emax,fact
      double precision :: Rfloor, Rnew
      integer :: ix^D
      logical :: use_Rfactor_T
      double precision :: Rfloor_const

      call fl%get_pthermal_equi(wCT,x,ixI^L,ixO^L,pth)
      call fl%get_rho_equi(wCT,x,ixI^L,ixO^L,rho)
      call fl%get_temperature_equi(wCT,x,ixI^L,ixO^L,Te)
      Rfactor(ixO^S)=pth(ixO^S)/(rho(ixO^S)*Te(ixO^S))
      res=0d0
      fact = fl%lref*qdt/fl%tref
      use_Rfactor_T = associated(fl%get_Rfactor_from_temperature)
      if (use_Rfactor_T) then
        call cooling_get_Rfactor_T(fl, fl%tlow, 1.d0, Rfloor_const)
      end if
      {do ix^DB = ixO^LIM^DB\}
           if (use_Rfactor_T) then
             Rfloor = Rfloor_const
           else
             Rfloor = Rfactor(ix^D)
           end if
           emin = rho(ix^D) * fl%tlow * Rfloor * invgam
           Lmax = max(zero,(pth(ix^D)*invgam-emin)/qdt)
           emax = max(zero, pth(ix^D)*invgam-emin)
           if( Te(ix^D)<=fl%tcoolmin ) then
             ! temperature is below floor level, no cooling.
             L1 = zero
           else if( Te(ix^D)>=fl%tcoolmax )then
             call calc_l_extended(Te(ix^D), L1,fl)
             L1 = L1*rho(ix^D)**2
             if(phys_trac) then
               if(Te(ix^D)<block%wextra(ix^D,fl%Tcoff_)) then
                 L1=L1*sqrt((Te(ix^D)/block%wextra(ix^D,fl%Tcoff_))**5)
               end if
             end if
             L1 = min(L1,Lmax)
             res(ix^D) =  L1*qdt
           else
             call findY(Te(ix^D),Y1,fl)
             Y2 = Y1 + fact * rho(ix^D)*rc_gamma_1
             call findT(Tlocal2,Y2,fl)
             if(Tlocal2<=fl%tcoolmin) then
               de = emax
             else
               call cooling_get_Rfactor_T( &
                    fl, Tlocal2, Rfactor(ix^D), Rnew)
               de = rho(ix^D) * &
                    (Rfactor(ix^D)*Te(ix^D) - Rnew*Tlocal2) * &
                    invgam
             end if
             if(phys_trac) then
               if(Te(ix^D)<block%wextra(ix^D,fl%Tcoff_)) then
                 de=de*sqrt((Te(ix^D)/block%wextra(ix^D,fl%Tcoff_))**5)
               end if
             end if
             de = min(max(zero, de), emax)
             res(ix^D) = de
           end if
      {end do\}
    end subroutine get_cool_equi

    subroutine cooling_get_Rfactor_T(fl, T, Rold, Rnew)
      type(rc_fluid), intent(in) :: fl
      double precision, intent(in) :: T, Rold
      double precision, intent(out) :: Rnew

      if (associated(fl%get_Rfactor_from_temperature)) then
        call fl%get_Rfactor_from_temperature(T, Rnew)
      else
        Rnew = Rold
      end if
    end subroutine cooling_get_Rfactor_T

    subroutine cool_exact(qdt,ixI^L,ixO^L,wCT,wCTprim,w,x,fl)
    !  Cooling routine using exact integration method from Townsend 2009
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(in)    :: qdt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw), wCTprim(ixI^S,1:nw)
      double precision, intent(inout) :: w(ixI^S,1:nw)
      type(rc_fluid), intent(in) :: fl
      double precision :: Y1, Y2
      double precision :: L1, Tlocal2, pnew(ixI^S)
      double precision :: Rfloor, Rnew, R1, R2
      double precision :: rho(ixI^S), Te(ixI^S), rhonew(ixI^S), Rfactor(ixI^S)
      double precision :: emin, Lmax, fact
      double precision :: de_thin, de_thick, emax
      double precision :: T1, T2, Tnew(ixI^S), tau, xi
      double precision :: xi_arr(ixI^S), emax_rem_arr(ixI^S)
      double precision :: cool_fac, fip_prim, frac_lowFIP, fip_factor
      integer :: ix^D
      logical :: use_Rfactor_T
      double precision :: Rfloor_const

      call fl%get_rho(wCT,x,ixI^L,ixO^L,rho)
      call fl%get_temperature(wCT,x,ixI^L,ixO^L,Te)
      call fl%get_pthermal(w,x,ixI^L,ixO^L,pnew)
      call fl%get_rho(w,x,ixI^L,ixO^L,rhonew)
      use_Rfactor_T = associated(fl%get_Rfactor_from_temperature)
      if (use_Rfactor_T) then
        {do ix^DB = ixO^LIM^DB\}
          call fl%get_Rfactor_from_temperature(Te(ix^D), Rfactor(ix^D))
        {end do\}
      else
        call fl%get_var_Rfactor(wCT, x, ixI^L, ixO^L, Rfactor)
      end if
      fact = fl%lref*qdt/fl%tref
      if (use_Rfactor_T) then
        call cooling_get_Rfactor_T(fl, fl%tlow, 1.d0, Rfloor_const)
      end if
      xi_arr = one
      emax_rem_arr = zero
      {do ix^DB = ixO^LIM^DB\}
        if (use_Rfactor_T) then
          Rfloor = Rfloor_const
        else
          Rfloor = Rfactor(ix^D)
        end if
        emin = rhonew(ix^D)*fl%tlow*Rfloor*invgam
        Lmax = max(zero,pnew(ix^D)*invgam-emin)/qdt
        emax = max(zero,pnew(ix^D)*invgam-emin)
        if (Te(ix^D)<=fl%tcoolmin) cycle
        if (fl%rad_newton) then
          xi = exp(-pnew(ix^D) / fl%rad_newton_pthick)
          xi = min(max(xi, zero), one)
        else
          xi = one
        end if
        cool_fac = xi

        ! ------- (A) OPTICALLY-THIN PART --------
        ! --- FIP factor with T-dependent r(T) ---
        if (fl%fip_ > 0) then
          fip_prim = min(maxfip, max(minfip, wCTprim(ix^D,fl%fip_)))
          ! frac_lowFIP(T) in [0,1]: low-FIP contribution to total Lambda at T
          frac_lowFIP = lowFIP_fraction(Te(ix^D), fl)
          ! Effective multiplicative factor on Lambda(T) for this time step:
          !   Lambda_eff(T) = g_fip * Lambda(T), with g_fip frozen at T^n
          fip_factor = one - frac_lowFIP + fip_prim * frac_lowFIP
          cool_fac = cool_fac * fip_factor
        end if
        ! --- geometric damping of optically thin radiative losses ---
        ! Suppress optically thin cooling near the lower boundary (chromosphere/photosphere)
        ! using a phenomenological Gaussian height-dependent damping.
        {^IFONED
        if (slab_uniform .and. fl%rad_damp .and. x(ix^D,ndim) <= xprobmin1 + fl%rad_damp_height) then
          cool_fac = cool_fac * exp(-(x(ix^D,ndim)-xprobmin1-fl%rad_damp_height)**2/fl%rad_damp_scale**2)
        end if
        if (fl%rc_is_1d_loop .and. slab_uniform .and. fl%rad_damp &
               .and. x(ix^D,ndim) >= xprobmax1 - fl%rad_damp_height) then
          cool_fac = cool_fac * exp(-(x(ix^D,ndim)-xprobmax1+fl%rad_damp_height)**2/fl%rad_damp_scale**2)
        end if
        }
        {^IFTWOD
        if (slab_uniform .and. fl%rad_damp .and. x(ix^D,ndim) <= xprobmin2 + fl%rad_damp_height) then
          cool_fac = cool_fac * exp(-(x(ix^D,ndim)-xprobmin2-fl%rad_damp_height)**2/fl%rad_damp_scale**2)
        end if
        }
        {^IFTHREED
        if (slab_uniform .and. fl%rad_damp .and. x(ix^D,ndim) <= xprobmin3 + fl%rad_damp_height) then
          cool_fac = cool_fac * exp(-(x(ix^D,ndim)-xprobmin3-fl%rad_damp_height)**2/fl%rad_damp_scale**2)
        end if
        }
        !  If the temperature is higher than the maximum table value,
        !  assume Bremsstrahlung and cool explicitly.
        if (Te(ix^D) >= fl%tcoolmax) then
          call calc_l_extended(Te(ix^D), L1, fl)
          L1 = L1 * rho(ix^D)**2
          L1 = min(cool_fac * L1, Lmax)
          de_thin = L1*qdt
          w(ix^D,fl%e_) = w(ix^D,fl%e_) - de_thin
        else
          ! Map T^n to Y-space
          call findY(Te(ix^D), Y1, fl)
          ! FIP and cutoff are included here as multiplicative factors in Lambda_eff.
          Y2 = Y1 + cool_fac * fact * rho(ix^D) * rc_gamma_1
          ! Invert Y(T) to get T^{n+1}
          call findT(Tlocal2, Y2, fl)
          ! Convert delta_T to an energy loss delta_e, respecting internal-energy floor.
          if (Tlocal2 <= fl%tcoolmin) then
            de_thin = emax
          else
            call cooling_get_Rfactor_T(fl, Tlocal2, Rfactor(ix^D), Rnew)
            de_thin = rho(ix^D) * &
                 (Rfactor(ix^D)*Te(ix^D) - Rnew*Tlocal2) * invgam
          end if
          ! --- TRAC modification: approximate, NOT strictly EI ---
          ! This rescaling is applied *after* the EI step and depends on T^n,
          ! so it does not correspond to integrating dT/dt with a modified
          ! Lambda(T). Kept here for practical reason, but should be
          ! understood as an approximate correction on top of EI.
          if (phys_trac) then
            if (Te(ix^D) < block%wextra(ix^D,fl%Tcoff_)) then
              de_thin = de_thin * sqrt( (Te(ix^D)/block%wextra(ix^D,fl%Tcoff_))**5 )
            end if
          end if
          de_thin = min(max(zero, de_thin), emax)
          w(ix^D,fl%e_) = w(ix^D,fl%e_) - de_thin
          if (fl%rad_newton) then
            xi_arr(ix^D) = xi
            emax_rem_arr(ix^D) = max(zero, emax - de_thin)
          end if
        end if
      {end do\}

      ! ------- (B) OPTICALLY-THICK (NEWTON) PART --------
      if (fl%rad_newton) then
        call fl%get_temperature(w, x, ixI^L, ixO^L, Tnew)
        {do ix^DB = ixO^LIM^DB\}
          T1 = Tnew(ix^D)
          tau = max(0.1d0 * sqrt(fl%rad_newton_rhosurf/rho(ix^D)), 4.d0*qdt)
          T2 = fl%rad_newton_trad + (T1-fl%rad_newton_trad)*exp(-qdt/tau)
          call cooling_get_Rfactor_T(fl, T1, Rfactor(ix^D), R1)
          call cooling_get_Rfactor_T(fl, T2, R1, R2)
          de_thick = (one-xi_arr(ix^D)) * rho(ix^D) * (R1*T1-R2*T2)*invgam
          de_thick = min(de_thick, emax_rem_arr(ix^D))
          w(ix^D,fl%e_) = w(ix^D,fl%e_) - de_thick
        {end do\}
      end if
    end subroutine cool_exact

    subroutine calc_l_extended (tpoint, lpoint,fl)
    !  Calculate l for t beyond tcoolmax
    !  Assumes Bremsstrahlung for the interpolated tables
    !  Uses the power law for piecewise power laws
      double precision, intent(IN)  :: tpoint
      double precision, intent(OUT) :: lpoint
      type(rc_fluid), intent(in) :: fl

      if(fl%isPPL) then
        lpoint =fl%l_PPL(fl%n_PPL) * ( tpoint / fl%t_PPL(fl%n_PPL) )**fl%a_PPL(fl%n_PPL)
      else
        lpoint = fl%Lcool(fl%ncool) * sqrt( tpoint / fl%tcoolmax)
      end if
    end subroutine calc_l_extended

    double precision function lowFIP_fraction(tpoint, fl)
      use mod_global_parameters

      double precision, intent(in) :: tpoint
      type(rc_fluid), intent(in)   :: fl

      double precision :: lgtp
      integer :: jl

      if (tpoint <= fl%tcool(1)) then
        lowFIP_fraction = fl%frac_lowFIP(1)
        return
      else if (tpoint >= fl%tcool(fl%ncool)) then
        lowFIP_fraction = fl%frac_lowFIP(fl%ncool)
        return
      end if

      lgtp = dlog10(tpoint)
      jl   = int((lgtp - fl%lgtcoolmin) / fl%lgstep) + 1
      jl   = max(1, min(fl%ncool-1, jl))

      lowFIP_fraction = fl%frac_lowFIP(jl) &
           + (tpoint - fl%tcool(jl)) &
           * (fl%frac_lowFIP(jl+1) - fl%frac_lowFIP(jl)) &
           / (fl%tcool(jl+1) - fl%tcool(jl))
    end function lowFIP_fraction

    subroutine findL (tpoint,Lpoint,fl)
    !  Fast search option to find correct point
    !  in cooling curve
      use mod_global_parameters

      double precision,intent(IN)   :: tpoint
      double precision, intent(OUT) :: Lpoint
      type(rc_fluid), intent(in) :: fl

      double precision :: lgtp
      integer :: jl,i

      if(fl%isPPL) then
        i = maxloc(fl%t_PPL, dim=1, mask=fl%t_PPL<tpoint)
        Lpoint = fl%l_PPL(i) * (tpoint / fl%t_PPL(i))**fl%a_PPL(i)
      else
        lgtp = dlog10(tpoint)
        jl = int((lgtp - fl%lgtcoolmin) /fl%lgstep) + 1
        Lpoint = fl%Lcool(jl)+ (tpoint-fl%tcool(jl)) &
                  * (fl%Lcool(jl+1)-fl%Lcool(jl)) &
                  / (fl%tcool(jl+1)-fl%tcool(jl))
      end if

    end subroutine findL

    subroutine findY (tpoint,Ypoint,fl)
    !  Fast search option to find correct point in cooling time
      use mod_global_parameters

      double precision,intent(IN)   :: tpoint
      double precision, intent(OUT) :: Ypoint
      type(rc_fluid), intent(in) :: fl

      double precision :: lgtp
      double precision :: y_extra,factor
      integer :: jl,i

      if(fl%isPPL) then
        i = maxloc(fl%t_PPL, dim=1, mask=fl%t_PPL<tpoint)
        factor = fl%l_PPL(fl%n_PPL+1) * fl%t_PPL(i) / (fl%l_PPL(i) * fl%t_PPL(fl%n_PPL+1))
        if(fl%a_PPL(i)==1.d0) then
          y_extra = log( fl%t_PPL(i) / tpoint )
        else
          y_extra = 1 / (1 - fl%a_PPL(i)) * (1 - ( fl%t_PPL(i) / tpoint )**(fl%a_PPL(i)-1) )
        end if
        Ypoint = fl%y_PPL(i) + factor*y_extra
      else
        lgtp = dlog10(tpoint)
        jl = int((lgtp - fl%lgtcoolmin) / fl%lgstep) + 1
        Ypoint = fl%Yc(jl)+ (tpoint-fl%tcool(jl)) &
                  * (fl%Yc(jl+1)-fl%Yc(jl)) &
                  / (fl%tcool(jl+1)-fl%tcool(jl))
      end if

    end subroutine findY

    subroutine findT (tpoint,Ypoint,fl)
    !  Fast search option to find correct temperature
    !  from temporal evolution function. Only possible this way because T is a monotonously
    !  decreasing function for the interpolated tables
    !  Uses eq. A7 from Townsend 2009 for piecewise power laws
      use mod_global_parameters

      double precision,intent(OUT)   :: tpoint
      double precision, intent(IN) :: Ypoint
      type(rc_fluid), intent(in) :: fl

      double precision :: factor
      integer :: jl,jc,jh,i

      if(fl%isPPL) then
        i = minloc(fl%y_PPL, dim=1, mask=fl%y_PPL>Ypoint)
        factor =  fl%l_PPL(i) * fl%t_PPL(fl%n_PPL+1) / (fl%l_PPL(fl%n_PPL+1) * fl%t_PPL(i))
        if(fl%a_PPL(i)==1.d0) then
          tpoint = fl%t_PPL(i) * exp( -1.d0 * factor * ( Ypoint - fl%y_PPL(i)))
        else
          tpoint = fl%t_PPL(i) * (1 - (1 - fl%a_PPL(i)) * factor * (Ypoint - fl%y_PPL(i)))**(1 / (1 - fl%a_PPL(i)))
        end if
      else
        if(Ypoint >= fl%Yc(1)) then
          tpoint = fl%tcoolmin
        else if (Ypoint == fl%Yc(fl%ncool)) then
          tpoint = fl%tcoolmax
        else
          jl=0
          jh=fl%ncool+1
          do
            if(jh-jl <= 1) exit
            jc=(jh+jl)/2
            if(Ypoint <= fl%Yc(jc)) then
              jl=jc
            else
              jh=jc
            end if
          end do
          ! Linear interpolation to obtain correct temperature
          tpoint = fl%tcool(jl)+ (Ypoint-fl%Yc(jl)) &
                 * (fl%tcool(jl+1)-fl%tcool(jl)) &
                 / (fl%Yc(jl+1)-fl%Yc(jl))
        end if
      end if
    end subroutine findT

end module mod_radiative_cooling
