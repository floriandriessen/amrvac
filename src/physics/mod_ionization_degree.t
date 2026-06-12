!> module ionization degree - get ionization degree for given temperature
module mod_ionization_degree
  implicit none
  private
  double precision :: iz_H_CL(1:100)
  double precision :: Te_H_CL(1:100)
  double precision :: ionization_He_abundance = 0.1d0
  double precision :: ionization_Rfactor_norm = 2.3d0
  double precision, allocatable :: prom_T_table(:)
  double precision, allocatable :: prom_p_table(:)
  double precision, allocatable :: prom_T_eos_table(:)
  double precision, allocatable :: prom_Y_eos_table(:,:)
  double precision :: Te_table_min, Te_table_max, inv_Te_table_step
  double precision :: Te_low_iz_He=5413.d0
  double precision :: Y_table_min, Y_table_max
  double precision :: ionization_H_energy_unit = 0.d0
  double precision, dimension(:), allocatable :: Rfactor_table
  double precision, dimension(:), allocatable :: Y_table
  double precision, dimension(:), allocatable :: eps_ion_H_table
  double precision, dimension(:), allocatable :: Te_H_table
  double precision, dimension(:), allocatable :: iz_H_table
  integer, parameter :: ionization_mode_carlsson = 1
  integer, parameter :: ionization_mode_gunar2025 = 2
  integer, parameter :: n_prom_T = 10
  integer, parameter :: n_prom_p = 7
  integer, parameter :: n_prom_eos = 4000
  integer :: ionization_table_mode = ionization_mode_carlsson
  ! Initialize
  public :: ionization_degree_init
  ! Current rho,p -> complete EOS state
  public :: ionization_get_state
  public :: ionization_get_state_scalar
  ! Prescribed T,p -> Rfactor and optional ionization degrees
  public :: ionization_state_Tp
  ! Prescribed rho,T -> consistent pthermal and Rfactor
  public :: ionization_solve_p_Rfactor
  ! T-only cooling callback and mode query
  public :: ionization_get_Rfactor_from_temperature
  public :: ionization_is_temperature_only
  public :: ionization_includes_energy
  public :: ionization_get_state_from_eint
  public :: ionization_get_p_eint_from_rho_T
  public :: ionization_get_eps_derivative_T
  public :: ionization_get_csound2_T
  public :: ionization_check_eint_table
  ! Compatibility interfaces
  public :: ionization_state_from_temperature_scalar
  ! Legacy compatibility API
  public :: ionization_degree_from_temperature
  logical :: ionization_initialized = .false.
  logical :: ionization_include_energy = .false.

  double precision, parameter :: prom_T_K(n_prom_T) = (/ &
       6000.d0, 7000.d0, 8000.d0, 9000.d0, 10000.d0, &
       12000.d0, 15000.d0, 20000.d0, 25000.d0, 30000.d0 /)

  double precision, parameter :: prom_p_cgs(n_prom_p) = (/ &
        0.03d0, 0.05d0, 0.10d0, 0.30d0, 0.50d0, 0.70d0, 1.00d0 /)

  ! First index: temperature. Second index: pressure.
  ! Data are ordered one complete temperature column per pressure.
  double precision, parameter :: &
      prom_iz_H_data(n_prom_T,n_prom_p) = reshape((/ &
      ! p = 0.03 dyn cm^-2
      0.53d0, 0.59d0, 0.65d0, 0.73d0, 0.80d0, &
      0.90d0, 0.95d0, 0.98d0, 0.99d0, 1.00d0, &
      ! p = 0.05 dyn cm^-2
      0.40d0, 0.47d0, 0.56d0, 0.69d0, 0.79d0, &
      0.91d0, 0.96d0, 0.98d0, 0.99d0, 1.00d0, &
      ! p = 0.10 dyn cm^-2
      0.25d0, 0.32d0, 0.47d0, 0.67d0, 0.81d0, &
      0.93d0, 0.97d0, 0.99d0, 0.99d0, 1.00d0, &
      ! p = 0.30 dyn cm^-2
      0.10d0, 0.17d0, 0.41d0, 0.69d0, 0.85d0, &
      0.96d0, 0.98d0, 0.99d0, 1.00d0, 1.00d0, &
      ! p = 0.50 dyn cm^-2
      0.06d0, 0.14d0, 0.38d0, 0.70d0, 0.87d0, &
      0.97d0, 0.99d0, 1.00d0, 1.00d0, 1.00d0, &
      ! p = 0.70 dyn cm^-2
      0.04d0, 0.12d0, 0.36d0, 0.70d0, 0.88d0, &
      0.97d0, 0.99d0, 1.00d0, 1.00d0, 1.00d0, &
      ! p = 1.00 dyn cm^-2
      0.03d0, 0.10d0, 0.34d0, 0.71d0, 0.89d0, &
      0.98d0, 0.99d0, 1.00d0, 1.00d0, 1.00d0 /), &
      (/ n_prom_T, n_prom_p /))

  ! Carlsson, M., & Leenaarts, J. 2012, A&A, 539, A39
  data Te_H_CL /2.1000005d+03, 2.2349495d+03, 2.3785708d+03, 2.5314199d+03, 2.6940928d+03, &
                2.8672190d+03, 3.0514692d+03, 3.2475613d+03, 3.4562544d+03, 3.6783584d+03, &
                3.9147351d+03, 4.1662998d+03, 4.4340322d+03, 4.7189697d+03, 5.0222153d+03, &
                5.3449502d+03, 5.6884248d+03, 6.0539712d+03, 6.4430083d+03, 6.8570420d+03, &
                7.2976860d+03, 7.7666460d+03, 8.2657383d+03, 8.7969062d+03, 9.3622090d+03, &
                9.9638330d+03, 1.0604130d+04, 1.1285561d+04, 1.2010781d+04, 1.2782619d+04, &
                1.3604041d+04, 1.4478250d+04, 1.5408652d+04, 1.6398826d+04, 1.7452648d+04, &
                1.8574172d+04, 1.9767766d+04, 2.1038082d+04, 2.2390010d+04, 2.3828838d+04, &
                2.5360100d+04, 2.6989764d+04, 2.8724182d+04, 3.0570023d+04, 3.2534518d+04, &
                3.4625215d+04, 3.6850262d+04, 3.9218293d+04, 4.1738543d+04, 4.4420699d+04, &
                4.7275266d+04, 5.0313219d+04, 5.3546391d+04, 5.6987395d+04, 6.0649453d+04, &
                6.4546914d+04, 6.8694758d+04, 7.3109148d+04, 7.7807289d+04, 8.2807258d+04, &
                8.8128625d+04, 9.3791852d+04, 9.9819000d+04, 1.0623346d+05, 1.1306024d+05, &
                1.2032560d+05, 1.2805797d+05, 1.3628709d+05, 1.4504503d+05, 1.5436592d+05, &
                1.6428561d+05, 1.7484295d+05, 1.8607852d+05, 1.9803609d+05, 2.1076208d+05, &
                2.2430608d+05, 2.3872045d+05, 2.5406084d+05, 2.7038703d+05, 2.8776234d+05, &
                3.0625456d+05, 3.2593475d+05, 3.4688000d+05, 3.6917081d+05, 3.9289406d+05, &
                4.1814181d+05, 4.4501247d+05, 4.7360988d+05, 5.0404450d+05, 5.3643488d+05, &
                5.7090725d+05, 6.0759431d+05, 6.4663956d+05, 6.8819319d+05, 7.3241712d+05, &
                7.7948294d+05, 8.2957412d+05, 8.8288331d+05, 9.3961919d+05, 1.0000000d+06  /

  ! Carlsson, M., & Leenaarts, J. 2012, A&A, 539, A39
  data iz_H_CL /1.0000000d+00, 1.0000000d+00, 1.0000000d+00, 1.0000000d+00, 1.0000000d+00, &
                1.0000000d+00, 1.0000000d+00, 9.9999928d-01, 9.9999774d-01, 9.9999624d-01, &
                9.9999428d-01, 9.9999219d-01, 9.9998862d-01, 9.9998313d-01, 9.9997944d-01, &
                9.9995774d-01, 1.0000324d+00, 9.9822283d-01, 9.6391600d-01, 8.8995951d-01, &
                8.0025935d-01, 7.4219799d-01, 7.1378660d-01, 6.9565988d-01, 6.7326778d-01, &
                6.5032905d-01, 6.2487972d-01, 5.9925246d-01, 5.7028764d-01, 5.3673893d-01, &
                4.9144468d-01, 4.4332290d-01, 3.9158964d-01, 3.4621361d-01, 2.7826238d-01, &
                1.9633104d-01, 1.2593637d-01, 8.1383914d-02, 5.0944358d-02, 2.8770180d-02, &
                1.4449068d-02, 7.6071853d-03, 4.3617180d-03, 2.7694483d-03, 1.8341533d-03, &
                1.2467797d-03, 8.6159195d-04, 6.0232455d-04, 4.2862241d-04, 3.0666622d-04, &
                2.2256430d-04, 1.6328457d-04, 1.2179965d-04, 9.2065704d-05, 7.0357783d-05, &
                5.4168402d-05, 4.2026968d-05, 3.2807184d-05, 2.5801779d-05, 2.0414085d-05, &
                1.6271379d-05, 1.3041737d-05, 1.0538992d-05, 8.5329248d-06, 6.9645471d-06, &
                5.7069005d-06, 4.7155136d-06, 3.9203474d-06, 3.2734958d-06, 2.7415674d-06, &
                2.2911979d-06, 1.9045801d-06, 1.5692771d-06, 1.2782705d-06, 1.0281080d-06, &
                8.1567373d-07, 6.4200310d-07, 5.0668103d-07, 4.0577518d-07, 3.3245340d-07, &
                2.8006056d-07, 2.4119515d-07, 2.1033340d-07, 1.8429903d-07, 1.6189058d-07, &
                1.4264889d-07, 1.2616209d-07, 1.1193421d-07, 9.9408616d-08, 8.8335305d-08, &
                7.8593857d-08, 7.0042269d-08, 6.2532941d-08, 5.5898361d-08, 5.0049657d-08, &
                4.4888207d-08, 4.0317172d-08, 3.6192560d-08, 3.2330732d-08, 2.9387850d-08  /

  contains

    subroutine ionization_degree_init(He_abundance, Rfactor_norm, table_name, include_energy)
      use mod_global_parameters, only: unit_temperature, unit_numberdensity, &
           unit_pressure, SI_unit
      use mod_constants, only: const_eV
      use mod_comm_lib, only: mpistop
      double precision, intent(in) :: He_abundance, Rfactor_norm
      character(len=*), intent(in), optional :: table_name
      character(len=32) :: selected_table
      logical, intent(in), optional :: include_energy
      double precision, parameter :: erg_to_joule = 1.0d-7

      if (ionization_initialized) then
        call mpistop("ionization_degree_init called more than once")
      end if

      if (He_abundance < 0.d0) call mpistop("negative He abundance")
      if (Rfactor_norm <= 0.d0) call mpistop("invalid Rfactor normalization")

      ionization_He_abundance = He_abundance
      ionization_Rfactor_norm = Rfactor_norm
      Te_low_iz_He = 5413.d0/unit_temperature
      ionization_include_energy = .false.
      if (present(include_energy)) ionization_include_energy = include_energy

      selected_table = "carlsson2012"
      if (present(table_name)) selected_table = trim(adjustl(table_name))

      select case (trim(selected_table))
      case ("carlsson2012")
        ionization_table_mode = ionization_mode_carlsson
      case ("gunar2025_a_l500")
        ionization_table_mode = ionization_mode_gunar2025
      case default
        call mpistop("unknown partial-ionization table")
      end select

      if (ionization_include_energy .and. &
          ionization_table_mode /= ionization_mode_carlsson) then
        call mpistop("ionization energy requires a T-only ionization table")
      end if

      if (ionization_include_energy) then
        if (unit_numberdensity <= 0.d0 .or. unit_pressure <= 0.d0) then
          call mpistop("invalid units for hydrogen ionization energy")
        end if
        if (SI_unit) then
          ! const_eV is in erg/eV; convert to joule for SI pressure units.
          ionization_H_energy_unit = &
               13.6d0*const_eV*erg_to_joule * &
               unit_numberdensity/unit_pressure
        else
          ionization_H_energy_unit = &
               13.6d0*const_eV*unit_numberdensity/unit_pressure
        end if
      else
        ionization_H_energy_unit = 0.d0
      end if

      select case (ionization_table_mode)
      case (ionization_mode_carlsson)
        call ionization_build_carlsson_table()
      case (ionization_mode_gunar2025)
        call ionization_build_prominence_table()
      end select
      ionization_initialized = .true.
    end subroutine ionization_degree_init

    subroutine ionization_build_carlsson_table()
      use mod_global_parameters, only: unit_temperature, zero

      double precision :: fact1, fact2, fact3
      double precision :: dL1, dL2, Te_table_step
      integer :: i, j, ntable, n_interpolated_table
      logical :: jump

      ntable=100
      n_interpolated_table=4000

      ! change neutral fraction to ionization degree
      iz_H_CL=1.d0-iz_H_CL

      allocate(Te_H_table(1:n_interpolated_table))
      allocate(iz_H_table(1:n_interpolated_table))
      Te_H_table=zero
      iz_H_table=zero

      ! create cooling table(s) for use in amrvac
      Te_table_step = (Te_H_CL(ntable)-Te_H_CL(1))/dble(n_interpolated_table-1)

      Te_H_table(1) = Te_H_CL(1)
      iz_H_table(1) = iz_H_CL(1)

      Te_H_table(n_interpolated_table) = Te_H_CL(ntable)
      iz_H_table(n_interpolated_table) = iz_H_CL(ntable)

      do i=2,n_interpolated_table        ! loop to create one table
        Te_H_table(i) = Te_H_table(i-1)+Te_table_step
        do j=1,ntable-1   ! loop to create one spot on a table
        ! Second order polynomial interpolation, except at the outer edge,
        ! or in case of a large jump.
          if(Te_H_table(i) < Te_H_CL(j+1)) then
             if(j.eq. ntable-1 )then
               fact1 = (Te_H_table(i)-Te_H_CL(j+1))     &
                     /(Te_H_CL(j)-Te_H_CL(j+1))

               fact2 = (Te_H_table(i)-Te_H_CL(j))       &
                     /(Te_H_CL(j+1)-Te_H_CL(j))

               iz_H_table(i) = iz_H_CL(j)*fact1 + iz_H_CL(j+1)*fact2
               exit
             else
               dL1 = iz_H_CL(j+1)-iz_H_CL(j)
               dL2 = iz_H_CL(j+2)-iz_H_CL(j+1)
               jump =(max(dabs(dL1),dabs(dL2)) > 2.d0*min(dabs(dL1),dabs(dL2)))
             end if

             if( jump ) then
               fact1 = (Te_H_table(i)-Te_H_CL(j+1))     &
                     /(Te_H_CL(j)-Te_H_CL(j+1))

               fact2 = (Te_H_table(i)-Te_H_CL(j))       &
                     /(Te_H_CL(j+1)-Te_H_CL(j))

               iz_H_table(i) = iz_H_CL(j)*fact1 + iz_H_CL(j+1)*fact2
               exit
             else
               fact1 = ((Te_H_table(i)-Te_H_CL(j+1))     &
                     * (Te_H_table(i)-Te_H_CL(j+2)))   &
                     / ((Te_H_CL(j)-Te_H_CL(j+1)) &
                     * (Te_H_CL(j)-Te_H_CL(j+2)))

               fact2 = ((Te_H_table(i)-Te_H_CL(j))       &
                     * (Te_H_table(i)-Te_H_CL(j+2)))   &
                     / ((Te_H_CL(j+1)-Te_H_CL(j)) &
                     * (Te_H_CL(j+1)-Te_H_CL(j+2)))

               fact3 = ((Te_H_table(i)-Te_H_CL(j))       &
                     * (Te_H_table(i)-Te_H_CL(j+1)))   &
                     / ((Te_H_CL(j+2)-Te_H_CL(j)) &
                     * (Te_H_CL(j+2)-Te_H_CL(j+1)))

               iz_H_table(i) = iz_H_CL(j)*fact1 + iz_H_CL(j+1)*fact2 &
                        + iz_H_CL(j+2)*fact3
               exit
             endif
          endif
        enddo  ! end loop to find create one spot on a table
      enddo    ! end loop to create one table

      ! nondimensionalize temperature
      Te_H_table=Te_H_table/unit_temperature
      Te_table_min=Te_H_table(1)
      Te_table_max=Te_H_table(n_interpolated_table)
      inv_Te_table_step=dble(n_interpolated_table-1)/(Te_table_max-Te_table_min)

      call ionization_build_eos_table()
    end subroutine ionization_build_carlsson_table

    subroutine ionization_build_prominence_table()
      use mod_global_parameters, only: unit_temperature, unit_pressure, SI_unit
      use mod_comm_lib, only: mpistop

      integer :: iT, ip
      double precision :: pressure_unit_cgs
      double precision :: dT, Rfactor

      if (unit_temperature <= 0.d0 .or. unit_pressure <= 0.d0) then
        call mpistop("invalid units for prominence ionization table")
      end if

      if (SI_unit) then
        ! unit_pressure is Pa; 1 Pa = 10 dyn cm^-2.
        pressure_unit_cgs = 10.d0*unit_pressure
      else
        ! In cgs mode unit_pressure is already dyn cm^-2.
        pressure_unit_cgs = unit_pressure
      end if

      allocate(prom_T_table(n_prom_T))
      allocate(prom_p_table(n_prom_p))
      allocate(prom_T_eos_table(n_prom_eos))
      allocate(prom_Y_eos_table(n_prom_eos,n_prom_p))

      prom_T_table = prom_T_K/unit_temperature
      prom_p_table = prom_p_cgs/pressure_unit_cgs

      dT = (1.d6 - 2.1000005d3) / dble(n_prom_eos-1)
      do iT = 1, n_prom_eos
        prom_T_eos_table(iT) = &
              (2.1000005d3 + dble(iT-1)*dT)/unit_temperature
      end do

      do ip = 1, n_prom_p
        do iT = 1, n_prom_eos
          call ionization_state_Tp(prom_T_eos_table(iT), prom_p_table(ip), Rfactor)
          prom_Y_eos_table(iT,ip) = &
                Rfactor*prom_T_eos_table(iT)
        end do
      end do

      call ionization_check_prominence_y_table()
    end subroutine ionization_build_prominence_table

    subroutine ionization_prominence_pressure_bracket(p, ip, weight)
      double precision, intent(in) :: p
      integer, intent(out) :: ip
      double precision, intent(out) :: weight

      integer :: j
      double precision :: p_use

      p_use = min(prom_p_table(n_prom_p), &
            max(prom_p_table(1), p))

      if (p_use <= prom_p_table(1)) then
        ip = 1
        weight = 0.d0
        return
      else if (p_use >= prom_p_table(n_prom_p)) then
        ip = n_prom_p-1
        weight = 1.d0
        return
      end if

      do j = 1, n_prom_p-1
        if (p_use <= prom_p_table(j+1)) then
          ip = j
          weight = (p_use-prom_p_table(j)) / &
                (prom_p_table(j+1)-prom_p_table(j))
          return
        end if
      end do
    end subroutine ionization_prominence_pressure_bracket

    !> Legacy array interface for T -> ionization degrees.
    !> Not used anymore; kept for compatibility and may be removed later.
    subroutine ionization_degree_from_temperature(ixI^L,ixO^L,Te,iz_H,iz_He)
      use mod_global_parameters
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: Te(ixI^S)
      double precision, intent(out) :: iz_H(ixO^S), iz_He(ixO^S)
      integer :: ix^D
      double precision :: Rtmp

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        call ionization_state_from_temperature_scalar( &
             Te(ix^D), Rtmp, iz_H(ix^D), iz_He(ix^D))
      {end do\}
    end subroutine ionization_degree_from_temperature

    subroutine ionization_H_degree_carlsson_scalar(T, iz_H)
      use mod_comm_lib, only: mpistop

      double precision, intent(in) :: T
      double precision, intent(out) :: iz_H

      integer :: i

      if (.not. allocated(Te_H_table) .or. &
          .not. allocated(iz_H_table)) then
        call mpistop("Carlsson ionization table is not initialized")
      end if

      if (T < Te_table_min) then
        iz_H = 0.d0
      else if (T >= Te_table_max) then
        iz_H = 1.d0
      else
        i = int((T-Te_table_min)*inv_Te_table_step) + 1
        i = max(1, min(size(Te_H_table)-1, i))

        iz_H = iz_H_table(i) + &
              (T-Te_H_table(i)) * &
              (iz_H_table(i+1)-iz_H_table(i)) * &
              inv_Te_table_step
      end if

      iz_H = min(1.d0, max(0.d0, iz_H))
    end subroutine ionization_H_degree_carlsson_scalar

    subroutine ionization_H_degree_prominence_scalar(T, p, iz_H)
      use mod_comm_lib, only: mpistop

      double precision, intent(in) :: T, p
      double precision, intent(out) :: iz_H

      integer :: iT, ip, j
      double precision :: weight_T, weight_p
      double precision :: iz_low, iz_high

      if (.not. allocated(prom_T_table) .or. &
          .not. allocated(prom_p_table)) then
        call mpistop("prominence ionization table is not initialized")
      end if

      call ionization_prominence_pressure_bracket(p, ip, weight_p)

      if (T <= prom_T_table(1)) then
        iT = 1
        weight_T = 0.d0
      else if (T >= prom_T_table(n_prom_T)) then
        iz_H = 1.d0
        return
      else
        iT = 1
        do j = 1, n_prom_T-1
          if (T <= prom_T_table(j+1)) then
            iT = j
            exit
          end if
        end do

        weight_T = (T-prom_T_table(iT)) / &
              (prom_T_table(iT+1)-prom_T_table(iT))
      end if

      if (iT == 1 .and. weight_T == 0.d0) then
        iz_H = (1.d0-weight_p)*prom_iz_H_data(1,ip) + &
              weight_p*prom_iz_H_data(1,ip+1)
      else
        iz_low = (1.d0-weight_p)*prom_iz_H_data(iT,ip) + &
              weight_p*prom_iz_H_data(iT,ip+1)

        iz_high = (1.d0-weight_p)*prom_iz_H_data(iT+1,ip) + &
              weight_p*prom_iz_H_data(iT+1,ip+1)

        iz_H = (1.d0-weight_T)*iz_low + weight_T*iz_high
      end if

      iz_H = min(1.d0, max(0.d0, iz_H))
    end subroutine ionization_H_degree_prominence_scalar

    subroutine ionization_get_state_scalar(rho, p, T, Rfactor, iz_H, iz_He)
      use mod_comm_lib, only: mpistop
      double precision, intent(in) :: rho, p
      double precision, intent(out) :: T, Rfactor
      double precision, intent(out), optional :: iz_H, iz_He

      double precision :: y
      double precision :: iz_H_loc, iz_He_loc

      if (.not. ionization_initialized) then
        call mpistop("ionization tables are not initialized")
      end if

      if (rho <= 0.d0 .or. p <= 0.d0) then
        select case (ionization_table_mode)
        case (ionization_mode_carlsson)
          T = Te_H_table(1)
        case (ionization_mode_gunar2025)
          T = prom_T_eos_table(1)
        end select
      else
        y = p/rho

        select case (ionization_table_mode)
        case (ionization_mode_carlsson)
          call ionization_invert_y_scalar(y, T)

        case (ionization_mode_gunar2025)
          call ionization_invert_y_prominence_scalar(y, p, T)

        case default
          call mpistop("invalid ionization table mode")
        end select
      end if

      call ionization_state_Tp( T, p, Rfactor, iz_H_loc, iz_He_loc)

      if (present(iz_H)) iz_H = iz_H_loc
      if (present(iz_He)) iz_He = iz_He_loc
    end subroutine ionization_get_state_scalar

    subroutine ionization_He_degree_scalar(T, iz_He)
      use mod_global_parameters, only: unit_temperature

      double precision, intent(in) :: T
      double precision, intent(out) :: iz_He

      if (T < Te_low_iz_He) then
        iz_He = 1.0084814d-4
      else
        iz_He = 1.d0 - &
              10.d0**(0.322571d0-5.96d-5*T*unit_temperature)
      end if

      iz_He = min(1.d0, max(0.d0, iz_He))
    end subroutine ionization_He_degree_scalar

    subroutine ionization_state_Tp(T, p, Rfactor, iz_H, iz_He)
      use mod_comm_lib, only: mpistop

      double precision, intent(in) :: T, p
      double precision, intent(out) :: Rfactor
      double precision, intent(out), optional :: iz_H, iz_He

      double precision :: iz_H_loc, iz_He_loc

      select case (ionization_table_mode)
      case (ionization_mode_carlsson)
        call ionization_H_degree_carlsson_scalar(T, iz_H_loc)

      case (ionization_mode_gunar2025)
        call ionization_H_degree_prominence_scalar( &
              T, p, iz_H_loc)

      case default
        call mpistop("invalid ionization table mode")
      end select

      call ionization_He_degree_scalar(T, iz_He_loc)

      iz_H_loc = min(1.d0, max(0.d0, iz_H_loc))
      iz_He_loc = min(1.d0, max(0.d0, iz_He_loc))

      Rfactor = ionization_Rfactor_from_degrees( &
            iz_H_loc, iz_He_loc)

      if (present(iz_H)) iz_H = iz_H_loc
      if (present(iz_He)) iz_He = iz_He_loc
    end subroutine ionization_state_Tp

    subroutine ionization_solve_p_Rfactor(rho, T, p, Rfactor)
      use mod_comm_lib, only: mpistop

      integer, parameter :: max_expand = 32, max_iter = 80
      double precision, intent(in) :: rho, T
      double precision, intent(out) :: p, Rfactor

      integer :: iter
      double precision :: plo, phi, pmid
      double precision :: Rlo, Rhi, Rmid
      double precision :: fhi, fmid, scale

      if (.not. ionization_initialized) then
        call mpistop("ionization pressure solver called before initialization")
      end if
      if (rho <= 0.d0 .or. T <= 0.d0) then
        call mpistop("invalid rho or T in ionization pressure solver")
      end if

      if (ionization_is_temperature_only()) then
        call ionization_state_from_temperature_scalar(T, Rfactor)
        p = rho*Rfactor*T
        return
      end if

      ! Low-pressure clamp branch
      plo = prom_p_table(1)
      call ionization_state_Tp(T, plo, Rlo)
      p = rho*Rlo*T
      if (p <= plo) then
        Rfactor = Rlo
        return
      end if

      ! High-pressure clamp branch
      phi = prom_p_table(n_prom_p)
      call ionization_state_Tp(T, phi, Rhi)
      p = rho*Rhi*T
      if (p >= phi) then
        Rfactor = Rhi
        return
      end if

      ! Now the root must lie inside [plo, phi]
      do iter = 1, max_iter
        pmid = 0.5d0*(plo+phi)
        call ionization_state_Tp(T, pmid, Rmid)
        fmid = pmid-rho*Rmid*T
        scale = max(abs(pmid), abs(rho*Rmid*T), tiny(1.d0))

        if (abs(fmid) <= 1.d-12*scale) then
          p = pmid
          Rfactor = Rmid
          return
        end if

        if (fmid > 0.d0) then
          phi = pmid
        else
          plo = pmid
        end if
      end do

      call mpistop("pressure-dependent EOS solve did not converge")
    end subroutine ionization_solve_p_Rfactor

    subroutine ionization_get_state(ixI^L, ixO^L, rho, p, T, Rfactor, iz_H, iz_He)
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: rho(ixI^S), p(ixI^S)
      double precision, intent(out) :: T(ixI^S), Rfactor(ixI^S)
      double precision, intent(out), optional :: iz_H(ixI^S), iz_He(ixI^S)

      double precision :: iz_H_loc, iz_He_loc
      integer :: ix^D

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        call ionization_get_state_scalar(rho(ix^D), p(ix^D), T(ix^D), &
             Rfactor(ix^D), iz_H_loc, iz_He_loc)

        if (present(iz_H)) iz_H(ix^D) = iz_H_loc
        if (present(iz_He)) iz_He(ix^D) = iz_He_loc
      {end do\}
    end subroutine ionization_get_state

    subroutine ionization_invert_y_scalar(y, T)
      use mod_comm_lib, only: mpistop
      double precision, intent(in) :: y
      double precision, intent(out) :: T

      integer :: ilo, ihi, imid, n
      double precision :: denom

      if (.not. allocated(Y_table) .or. .not. allocated(Rfactor_table)) then
        call mpistop("ionization_invert_y_scalar: EOS tables are not initialized")
      end if

      n = size(Y_table)

      if (n < 2) then
        call mpistop("ionization_invert_y_scalar: EOS table has fewer than 2 points")
      end if

      if (y <= 0.d0) then
        T = Te_H_table(1)
        return
      end if

      if (y < Y_table_min) then
        if (Rfactor_table(1) > 0.d0) then
          T = max(tiny(1.d0), y/Rfactor_table(1))
        else
          T = Te_H_table(1)
        end if
        return
      end if

      if (y >= Y_table_max) then
        if (Rfactor_table(n) > 0.d0) then
          T = y/Rfactor_table(n)
        else
          T = Te_H_table(n)
        end if
        return
      end if

      ilo = 1
      ihi = n

      do while (ihi - ilo > 1)
        imid = (ilo + ihi)/2
        if (y >= Y_table(imid)) then
          ilo = imid
        else
          ihi = imid
        end if
      end do

      denom = Y_table(ilo+1) - Y_table(ilo)

      if (denom <= max(tiny(denom), 1.d-12*dabs(Y_table(ilo)))) then
        call mpistop("ionization_invert_y_scalar: non-increasing Y_table interval")
      end if

      T = Te_H_table(ilo) + (y - Y_table(ilo)) * &
          (Te_H_table(ilo+1) - Te_H_table(ilo))/denom
    end subroutine ionization_invert_y_scalar

    subroutine ionization_state_from_temperature_scalar(T, Rfactor, iz_H, iz_He)
      use mod_comm_lib, only: mpistop

      double precision, intent(in) :: T
      double precision, intent(out) :: Rfactor
      double precision, intent(out), optional :: iz_H, iz_He
      double precision :: iz_H_loc, iz_He_loc

      if (ionization_table_mode /= ionization_mode_carlsson) then
        call mpistop("temperature-only ionization state requested "// &
              "for a pressure-dependent table")
      end if

      call ionization_H_degree_carlsson_scalar(T, iz_H_loc)
      call ionization_He_degree_scalar(T, iz_He_loc)

      Rfactor = ionization_Rfactor_from_degrees( &
            iz_H_loc, iz_He_loc)

      if (present(iz_H)) iz_H = iz_H_loc
      if (present(iz_He)) iz_He = iz_He_loc
    end subroutine ionization_state_from_temperature_scalar

    subroutine ionization_check_prominence_y_table()
      use mod_comm_lib, only: mpistop
      integer :: iT, ip

      do ip = 1, n_prom_p
        do iT = 1, n_prom_eos-1
          if (prom_Y_eos_table(iT+1,ip) <= &
              prom_Y_eos_table(iT,ip)) then
            call mpistop("prominence EOS Y(T,p) is not monotonic")
          end if
        end do
      end do
    end subroutine ionization_check_prominence_y_table

    subroutine ionization_invert_y_prominence_scalar(y, p, T)
      use mod_comm_lib, only: mpistop

      double precision, intent(in) :: y, p
      double precision, intent(out) :: T

      integer :: ip, ilo, ihi, imid
      double precision :: wp, Ylo, Yhi, Ymid
      double precision :: Ymin, Ymax, Rend, denom

      if (.not. allocated(prom_Y_eos_table)) then
        call mpistop("prominence EOS inversion table is missing")
      end if

      call ionization_prominence_pressure_bracket(p, ip, wp)

      Ymin = (1.d0-wp)*prom_Y_eos_table(1,ip) + &
            wp*prom_Y_eos_table(1,ip+1)

      Ymax = (1.d0-wp)*prom_Y_eos_table(n_prom_eos,ip) + &
            wp*prom_Y_eos_table(n_prom_eos,ip+1)

      if (y <= 0.d0) then
        T = prom_T_eos_table(1)
        return
      end if

      if (y < Ymin) then
        Rend = Ymin/prom_T_eos_table(1)
        T = max(tiny(1.d0), y/Rend)
        return
      end if

      if (y >= Ymax) then
        Rend = Ymax/prom_T_eos_table(n_prom_eos)
        T = y/Rend
        return
      end if

      ilo = 1
      ihi = n_prom_eos

      do while (ihi-ilo > 1)
        imid = (ilo+ihi)/2
        Ymid = (1.d0-wp)*prom_Y_eos_table(imid,ip) + &
              wp*prom_Y_eos_table(imid,ip+1)

        if (y >= Ymid) then
          ilo = imid
        else
          ihi = imid
        end if
      end do

      Ylo = (1.d0-wp)*prom_Y_eos_table(ilo,ip) + &
            wp*prom_Y_eos_table(ilo,ip+1)

      Yhi = (1.d0-wp)*prom_Y_eos_table(ilo+1,ip) + &
            wp*prom_Y_eos_table(ilo+1,ip+1)

      denom = Yhi-Ylo
      if (denom <= max(tiny(denom),1.d-12*abs(Ylo))) then
        call mpistop("invalid prominence EOS inversion interval")
      end if

      T = prom_T_eos_table(ilo) + (y-Ylo) * &
            (prom_T_eos_table(ilo+1)-prom_T_eos_table(ilo)) / &
            denom
    end subroutine ionization_invert_y_prominence_scalar

    subroutine ionization_build_eos_table()
      use mod_comm_lib, only: mpistop
      integer :: i, n
      double precision :: iz_H

      if (.not. allocated(Te_H_table) .or. .not. allocated(iz_H_table)) then
        call mpistop("ionization_build_eos_table: ionization tables are not initialized")
      end if

      if (allocated(Rfactor_table)) deallocate(Rfactor_table)
      if (allocated(Y_table)) deallocate(Y_table)
      if (allocated(eps_ion_H_table)) deallocate(eps_ion_H_table)

      n = size(Te_H_table)

      if (n < 2) then
        call mpistop("ionization_build_eos_table: Te_H_table has fewer than 2 points")
      end if

      allocate(Rfactor_table(1:n))
      allocate(Y_table(1:n))
      if (ionization_include_energy) allocate(eps_ion_H_table(1:n))

      do i = 1, n
        call ionization_state_from_temperature_scalar( &
             Te_H_table(i), Rfactor_table(i), iz_H)
        Y_table(i) = Rfactor_table(i)*Te_H_table(i)
        if (ionization_include_energy) then
          eps_ion_H_table(i) = ionization_H_energy_unit*iz_H
        end if
      end do

      Y_table_min = Y_table(1)
      Y_table_max = Y_table(n)

      call ionization_check_y_table()
    end subroutine ionization_build_eos_table

    subroutine ionization_check_y_table()
      use mod_comm_lib, only: mpistop
      integer :: i

      if (.not. allocated(Y_table)) then
        call mpistop("ionization_check_y_table: Y_table is not allocated")
      end if

      if (size(Y_table) < 2) then
        call mpistop("ionization_check_y_table: Y_table has fewer than 2 points")
      end if

      do i = 1, size(Y_table)-1
        if (Y_table(i+1) <= Y_table(i)) then
          call mpistop("ionization_check_y_table: Y(T)=Rfactor(T)*T is not strictly increasing")
        end if
      end do
    end subroutine ionization_check_y_table

    double precision function ionization_Rfactor_from_degrees(iz_H, iz_He)
      double precision, intent(in) :: iz_H, iz_He

      ionization_Rfactor_from_degrees = &
           (1.d0 + iz_H + ionization_He_abundance * &
           (1.d0 + iz_He*(1.d0 + iz_He))) / ionization_Rfactor_norm
    end function ionization_Rfactor_from_degrees

    subroutine ionization_get_Rfactor_from_temperature(T, Rfactor)
      double precision, intent(in)  :: T
      double precision, intent(out) :: Rfactor

      call ionization_state_from_temperature_scalar(T, Rfactor)
    end subroutine ionization_get_Rfactor_from_temperature

    logical function ionization_is_temperature_only()
      ionization_is_temperature_only = &
            ionization_table_mode == ionization_mode_carlsson
    end function ionization_is_temperature_only

    logical function ionization_includes_energy()
      ionization_includes_energy = ionization_include_energy
    end function ionization_includes_energy

    double precision function ionization_eint_table_value(i, invgam)
      integer, intent(in) :: i
      double precision, intent(in) :: invgam

      ionization_eint_table_value = &
           Y_table(i)*invgam + eps_ion_H_table(i)
    end function ionization_eint_table_value

    subroutine ionization_check_eint_table(invgam)
      use mod_comm_lib, only: mpistop
      double precision, intent(in) :: invgam
      integer :: i
      double precision :: eps1, eps2

      if (.not. ionization_include_energy) return
      if (.not. allocated(eps_ion_H_table)) then
        call mpistop("ionization-energy table is not initialized")
      end if
      if (invgam <= 0.d0) then
        call mpistop("invalid inverse gamma minus one for ionization EOS")
      end if

      do i = 1, size(Y_table)-1
        eps1 = ionization_eint_table_value(i, invgam)
        eps2 = ionization_eint_table_value(i+1, invgam)
        if (eps2 <= eps1) then
          call mpistop("ionization specific-energy table is not monotonic")
        end if
      end do
    end subroutine ionization_check_eint_table

    subroutine ionization_invert_eint_scalar(eps, invgam, T)
      use mod_comm_lib, only: mpistop
      double precision, intent(in) :: eps, invgam
      double precision, intent(out) :: T

      integer :: ilo, ihi, imid, n, iter
      double precision :: eps_min, eps_max, eps_lo, eps_hi, denom
      double precision :: eps_eval, slope, Rfactor, iz_H, Tnew

      if (.not. ionization_include_energy .or. &
          .not. allocated(eps_ion_H_table)) then
        call mpistop("ionization-energy EOS is not initialized")
      end if
      if (invgam <= 0.d0) then
        call mpistop("invalid inverse gamma minus one in energy inversion")
      end if

      n = size(Y_table)
      eps_min = ionization_eint_table_value(1, invgam)
      eps_max = ionization_eint_table_value(n, invgam)

      if (eps < eps_min) then
        T = max(tiny(1.d0), &
             (eps-eps_ion_H_table(1))/(Rfactor_table(1)*invgam))
        return
      end if

      if (eps >= eps_max) then
        T = max(tiny(1.d0), &
             (eps-eps_ion_H_table(n))/(Rfactor_table(n)*invgam))
        return
      end if

      ilo = 1
      ihi = n
      do while (ihi-ilo > 1)
        imid = (ilo+ihi)/2
        if (eps >= ionization_eint_table_value(imid, invgam)) then
          ilo = imid
        else
          ihi = imid
        end if
      end do

      eps_lo = ionization_eint_table_value(ilo, invgam)
      eps_hi = ionization_eint_table_value(ilo+1, invgam)
      denom = eps_hi-eps_lo
      if (denom <= max(tiny(denom), 1.d-12*abs(eps_lo))) then
        call mpistop("invalid ionization-energy inversion interval")
      end if

      T = Te_H_table(ilo) + (eps-eps_lo) * &
           (Te_H_table(ilo+1)-Te_H_table(ilo))/denom

      ! The tabulated value supplies a close initial guess. Correct it against
      ! the actual interpolated ionization state so rho,eint -> T and
      ! rho,T -> eint remain mutually consistent inside each table interval.
      slope = denom/(Te_H_table(ilo+1)-Te_H_table(ilo))
      do iter = 1, 3
        call ionization_state_from_temperature_scalar( &
             T, Rfactor, iz_H)
        eps_eval = Rfactor*T*invgam + ionization_H_energy_unit*iz_H
        if (abs(eps_eval-eps) <= &
            1.d-12*max(abs(eps), tiny(1.d0))) exit
        Tnew = T-(eps_eval-eps)/slope
        T = max(Te_H_table(ilo), min(Te_H_table(ilo+1), Tnew))
      end do
    end subroutine ionization_invert_eint_scalar

    subroutine ionization_get_state_from_eint( &
         rho, eint, invgam, T, p, Rfactor, iz_H, iz_He)
      use mod_comm_lib, only: mpistop
      double precision, intent(in) :: rho, eint, invgam
      double precision, intent(out) :: T, p, Rfactor
      double precision, intent(out), optional :: iz_H, iz_He
      double precision :: iz_H_loc, iz_He_loc

      if (.not. ionization_include_energy) &
           call mpistop("ionization-energy EOS is not enabled")
      if (.not. ionization_is_temperature_only()) &
           call mpistop("ionization-energy EOS requires a T-only table")
      if (rho <= 0.d0 .or. eint <= 0.d0) &
           call mpistop("invalid rho or eint in ionization EOS")
      call ionization_invert_eint_scalar(eint/rho, invgam, T)
      call ionization_state_from_temperature_scalar( &
           T, Rfactor, iz_H_loc, iz_He_loc)
      p = rho*Rfactor*T
      if (present(iz_H)) iz_H = iz_H_loc
      if (present(iz_He)) iz_He = iz_He_loc
    end subroutine ionization_get_state_from_eint

    subroutine ionization_get_p_eint_from_rho_T( &
         rho, T, invgam, p, eint, Rfactor, iz_H, iz_He)
      use mod_comm_lib, only: mpistop
      double precision, intent(in) :: rho, T, invgam
      double precision, intent(out) :: p, eint, Rfactor
      double precision, intent(out), optional :: iz_H, iz_He
      double precision :: iz_H_loc, iz_He_loc

      if (rho <= 0.d0 .or. T <= 0.d0 .or. invgam <= 0.d0) then
        call mpistop("invalid state in rho,T ionization EOS mapping")
      end if

      if (ionization_include_energy) then
        if (.not. ionization_is_temperature_only()) then
          call mpistop("ionization energy requires a T-only ionization table")
        end if
        call ionization_state_from_temperature_scalar( &
             T, Rfactor, iz_H_loc, iz_He_loc)
        p = rho*Rfactor*T
        ! Only hydrogen ionization energy is included in this first version;
        ! helium ionization energy is intentionally deferred.
        eint = rho*(Rfactor*T*invgam + &
             ionization_H_energy_unit*iz_H_loc)
      else
        call ionization_solve_p_Rfactor(rho, T, p, Rfactor)
        call ionization_state_Tp(T, p, Rfactor, iz_H_loc, iz_He_loc)
        eint = p*invgam
      end if

      if (present(iz_H)) iz_H = iz_H_loc
      if (present(iz_He)) iz_He = iz_He_loc
    end subroutine ionization_get_p_eint_from_rho_T

    subroutine ionization_get_eps_derivative_T( &
         T, invgam, eps, deps_dT, dq_dT_out)
      use mod_comm_lib, only: mpistop
      double precision, intent(in) :: T, invgam
      double precision, intent(out) :: eps, deps_dT
      double precision, intent(out), optional :: dq_dT_out

      integer :: i, n
      double precision :: Rfactor, iz_H, q, dq_dT, depsion_dT

      if (.not. ionization_include_energy) then
        call mpistop("ionization-energy derivative requested while disabled")
      end if
      if (.not. ionization_is_temperature_only()) then
        call mpistop("ionization-energy derivative requires a T-only table")
      end if
      if (T <= 0.d0 .or. invgam <= 0.d0) then
        call mpistop("invalid state in ionization-energy derivative")
      end if

      n = size(Te_H_table)
      call ionization_state_from_temperature_scalar(T, Rfactor, iz_H)
      q = Rfactor*T
      eps = q*invgam + ionization_H_energy_unit*iz_H

      if (T < Te_table_min) then
        dq_dT = Rfactor_table(1)
        depsion_dT = 0.d0
      else if (T >= Te_table_max) then
        dq_dT = Rfactor_table(n)
        depsion_dT = 0.d0
      else
        i = int((T-Te_table_min)*inv_Te_table_step) + 1
        i = max(1, min(n-1, i))
        dq_dT = (Y_table(i+1)-Y_table(i)) / &
             (Te_H_table(i+1)-Te_H_table(i))
        depsion_dT = (eps_ion_H_table(i+1)-eps_ion_H_table(i)) / &
             (Te_H_table(i+1)-Te_H_table(i))
      end if

      deps_dT = dq_dT*invgam + depsion_dT
      if (deps_dT <= 0.d0) then
        call mpistop("non-positive heat capacity in ionization EOS")
      end if
      if (present(dq_dT_out)) dq_dT_out = dq_dT
    end subroutine ionization_get_eps_derivative_T

    subroutine ionization_get_csound2_T(T, invgam, csound2)
      double precision, intent(in) :: T, invgam
      double precision, intent(out) :: csound2

      double precision :: Rfactor, q, eps, deps_dT, dq_dT

      call ionization_get_eps_derivative_T( &
           T, invgam, eps, deps_dT, dq_dT)
      call ionization_state_from_temperature_scalar(T, Rfactor)
      q = Rfactor*T
      csound2 = q*(1.d0+dq_dT/deps_dT)
    end subroutine ionization_get_csound2_T

end module mod_ionization_degree
