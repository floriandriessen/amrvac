module mod_magnetic_topology
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan, &
       ieee_is_finite
  use mod_global_parameters, only: ndim,npe,slab_uniform, &
       xprobmin1,xprobmax1,xprobmin2,xprobmax2,xprobmin3,xprobmax3
  use mod_comm_lib, only: mpistop
  use mod_trace_field, only: trace_length_result,trace_twist_result, &
       trace_mapping_result,trace_qperp_result,trace_topology_result, &
       trace_field_length_single,trace_field_length_multi, &
       trace_field_twist_single,trace_field_twist_multi, &
       trace_field_mapping_single,trace_field_mapping_multi, &
       trace_field_topology_multi, &
       trace_field_qperp_single,trace_field_qperp_multi, &
       trace_status_boundary, &
       trace_face_none,trace_face_xmin,trace_face_xmax, &
       trace_face_ymin,trace_face_ymax,trace_face_zmin, &
       trace_face_zmax,trace_face_ambiguous
  implicit none
  private

  public :: mt_length_single,mt_length_seeds,mt_length_plane_xy
  public :: mt_length_plane_xz,mt_length_plane_yz
  public :: mt_fieldline_products_seeds
  public :: mt_twist_single,mt_twist_seeds,mt_twist_plane_xy
  public :: mt_twist_plane_xz,mt_twist_plane_yz
  public :: mt_mapping_plane_xy
  public :: mt_qperp_plane_xy,mt_qperp_plane_xz,mt_qperp_plane_yz
  public :: mt_qperp_plane_arbitrary
  public :: mt_fieldline_products_plane_arbitrary
  public :: mt_topology_plane_xy,mt_topology_plane_xz,mt_topology_plane_yz
  public :: mt_qsl_plane_vtu_xy,mt_qsl_plane_vtu_xz,mt_qsl_plane_vtu_yz
  public :: mt_write_cartesian_vti_pointdata
  public :: mt_fieldline_products_volume_vti
  public :: mt_params_read,mt_run_topology_task

  integer, parameter :: mt_vti_kind_float64 = 1
  integer, parameter :: mt_vti_kind_int32   = 2
  integer, parameter :: mt_vti_name_len     = 64
  integer, parameter :: mt_task_name_len    = 128
  double precision, parameter :: mt_unset_real = huge(1.d0)

  logical :: mt_enable = .false.
  character(len=mt_task_name_len) :: mt_mode = ''
  character(len=mt_task_name_len) :: mt_output_file = ''
  character(len=mt_task_name_len) :: mt_output_prefix = ''
  character(len=mt_task_name_len) :: mt_seed_file = ''
  character(len=mt_task_name_len) :: mt_vtk_detail = 'minimal'
  double precision :: mt_dL = -1.d0
  integer :: mt_max_steps = -1
  double precision :: mt_b_min = -1.d0
  logical :: mt_compute_twist = .false.
  logical :: mt_compute_qperp = .false.
  logical :: mt_write_csv = .false.
  character(len=mt_task_name_len) :: mt_plane = ''
  double precision :: mt_origin(3) = mt_unset_real
  double precision :: mt_e1(3) = mt_unset_real
  double precision :: mt_e2(3) = mt_unset_real
  double precision :: mt_xmin = mt_unset_real
  double precision :: mt_xmax = mt_unset_real
  double precision :: mt_ymin = mt_unset_real
  double precision :: mt_ymax = mt_unset_real
  double precision :: mt_zmin = mt_unset_real
  double precision :: mt_zmax = mt_unset_real
  double precision :: mt_x0 = mt_unset_real
  double precision :: mt_y0 = mt_unset_real
  double precision :: mt_z0 = mt_unset_real
  integer :: mt_nx = -1
  integer :: mt_ny = -1
  integer :: mt_nz = -1
  double precision :: mt_s1min = mt_unset_real
  double precision :: mt_s1max = mt_unset_real
  double precision :: mt_s2min = mt_unset_real
  double precision :: mt_s2max = mt_unset_real
  integer :: mt_n1 = -1
  integer :: mt_n2 = -1
  integer :: mt_chunk_nz = -1

  type, private :: mt_vti_array_desc
    character(len=mt_vti_name_len) :: name
    integer :: kind
    integer(kind=8) :: nbytes
    integer(kind=8) :: offset
  end type mt_vti_array_desc

  type, private :: mt_volume_products
    double precision, allocatable :: length_total(:)
    double precision, allocatable :: length_backward(:)
    double precision, allocatable :: length_forward(:)
    integer, allocatable :: nstep_backward_length(:)
    integer, allocatable :: nstep_forward_length(:)
    integer, allocatable :: status_backward_length(:)
    integer, allocatable :: status_forward_length(:)
    double precision, allocatable :: twist_total(:)
    double precision, allocatable :: twist_backward(:)
    double precision, allocatable :: twist_forward(:)
    integer, allocatable :: nstep_backward_twist(:)
    integer, allocatable :: nstep_forward_twist(:)
    integer, allocatable :: status_backward_twist(:)
    integer, allocatable :: status_forward_twist(:)
    double precision, allocatable :: qperp(:)
    double precision, allocatable :: logqperp(:)
    double precision, allocatable :: N2(:)
    double precision, allocatable :: bfactor(:)
    double precision, allocatable :: length_forward_qperp(:)
    double precision, allocatable :: length_backward_qperp(:)
    double precision, allocatable :: Bseed_norm(:)
    double precision, allocatable :: Bf_norm(:)
    double precision, allocatable :: Bb_norm(:)
    integer, allocatable :: valid_qperp(:)
    integer, allocatable :: status_qperp(:)
    integer, allocatable :: face_forward_qperp(:)
    integer, allocatable :: face_backward_qperp(:)
    integer, allocatable :: status_forward_qperp(:)
    integer, allocatable :: status_backward_qperp(:)
  end type mt_volume_products

contains

  subroutine mt_params_read(files)
    ! Read the optional one-task magnetic-topology namelist.
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)

    integer :: n

    namelist /magnetic_topology_list/ mt_enable,mt_mode,mt_output_file, &
         mt_output_prefix,mt_seed_file,mt_vtk_detail, &
         mt_dL,mt_max_steps,mt_b_min, &
         mt_compute_twist,mt_compute_qperp, &
         mt_write_csv,mt_plane, &
         mt_origin,mt_e1,mt_e2,mt_xmin,mt_xmax,mt_ymin,mt_ymax, &
         mt_zmin,mt_zmax,mt_nx,mt_ny,mt_nz,mt_x0,mt_y0,mt_z0, &
         mt_s1min,mt_s1max,mt_n1,mt_s2min,mt_s2max,mt_n2, &
         mt_chunk_nz

    call mt_set_default_params()
    do n=1,size(files)
      open(unitpar,file=trim(files(n)),status='old')
      read(unitpar,magnetic_topology_list,end=111)
111   close(unitpar)
    enddo
  end subroutine mt_params_read

  subroutine mt_run_topology_task()
    ! Dispatch the namelist-selected topology/QSL postprocessing task.
    character(len=mt_task_name_len) :: mode

    if (.not.mt_enable) return

    mode=mt_lowercase(trim(mt_mode))
    call mt_validate_common_params(mode)

    select case (trim(mode))
    case ('axis_plane_full_vtu')
      call mt_run_axis_plane_full_vtu_task()
    case ('volume_vti')
      call mt_run_volume_vti_task()
    case ('arbitrary_plane_products')
      call mt_run_arbitrary_plane_products_task()
    case ('seed_products')
      call mt_run_seed_products_task()
    case ('axis_plane_csv')
      call mt_run_axis_plane_csv_task()
    case default
      call mpistop('mt_run_topology_task: unknown mt_mode='//trim(mt_mode))
    end select
  end subroutine mt_run_topology_task

  subroutine mt_set_default_params()
    mt_enable=.false.
    mt_mode=''
    mt_output_file=''
    mt_output_prefix=''
    mt_seed_file=''
    mt_vtk_detail='minimal'
    mt_dL=-1.d0
    mt_max_steps=-1
    mt_b_min=-1.d0
    mt_compute_twist=.false.
    mt_compute_qperp=.false.
    mt_write_csv=.false.
    mt_plane=''
    mt_origin=mt_unset_real
    mt_e1=mt_unset_real
    mt_e2=mt_unset_real
    mt_xmin=mt_unset_real
    mt_xmax=mt_unset_real
    mt_ymin=mt_unset_real
    mt_ymax=mt_unset_real
    mt_zmin=mt_unset_real
    mt_zmax=mt_unset_real
    mt_x0=mt_unset_real
    mt_y0=mt_unset_real
    mt_z0=mt_unset_real
    mt_nx=-1
    mt_ny=-1
    mt_nz=-1
    mt_s1min=mt_unset_real
    mt_s1max=mt_unset_real
    mt_s2min=mt_unset_real
    mt_s2max=mt_unset_real
    mt_n1=-1
    mt_n2=-1
    mt_chunk_nz=-1
  end subroutine mt_set_default_params

  subroutine mt_validate_common_params(mode)
    character(len=*), intent(in) :: mode

    if (len_trim(mode)==0) then
      call mpistop('mt_run_topology_task requires mt_mode')
    endif
    if (len_trim(mt_output_file)==0 .and. &
        len_trim(mt_output_prefix)==0) then
      call mpistop('mt_run_topology_task requires mt_output_file or mt_output_prefix')
    endif
    if (mt_dL<=0.d0) then
      call mpistop('mt_run_topology_task requires mt_dL > 0')
    endif
    if (mt_max_steps<=0) then
      call mpistop('mt_run_topology_task requires mt_max_steps > 0')
    endif
    select case (mt_lowercase(trim(mt_vtk_detail)))
    case ('minimal','full')
    case default
      call mpistop('mt_run_topology_task requires mt_vtk_detail=minimal or full')
    end select
  end subroutine mt_validate_common_params

  subroutine mt_run_axis_plane_full_vtu_task()
    character(len=mt_task_name_len) :: plane
    character(len=mt_task_name_len) :: output_file

    plane=mt_lowercase(trim(mt_plane))
    call mt_resolve_output_file('axis_plane_full_vtu','.vtu', &
         '_'//trim(plane)//'_full.vtu',output_file)
    if (mt_compute_twist .or. mt_compute_qperp) then
      write(*,'(a)') 'mt_run_topology_task: axis_plane_full_vtu ignores compute flags'
      write(*,'(a)') 'mt_run_topology_task: it always writes full plane products'
    endif

    select case (trim(plane))
    case ('xy')
      call mt_require_real('mt_xmin',mt_xmin)
      call mt_require_real('mt_xmax',mt_xmax)
      call mt_require_real('mt_ymin',mt_ymin)
      call mt_require_real('mt_ymax',mt_ymax)
      call mt_require_real('mt_z0',mt_z0)
      call mt_require_positive_int('mt_nx',mt_nx)
      call mt_require_positive_int('mt_ny',mt_ny)
      call mt_require_ordered('mt_xmin','mt_xmax',mt_xmin,mt_xmax)
      call mt_require_ordered('mt_ymin','mt_ymax',mt_ymin,mt_ymax)
      write(*,'(a)') 'mt_run_topology_task: writing axis-plane full VTU '//trim(output_file)
      if (mt_b_min>0.d0) then
        call mt_qsl_plane_vtu_xy(mt_xmin,mt_xmax,mt_nx,mt_ymin,mt_ymax, &
             mt_ny,mt_z0,mt_dL,mt_max_steps,trim(output_file), &
             b_min=mt_b_min)
      else
        call mt_qsl_plane_vtu_xy(mt_xmin,mt_xmax,mt_nx,mt_ymin,mt_ymax, &
             mt_ny,mt_z0,mt_dL,mt_max_steps,trim(output_file))
      endif
    case ('xz')
      call mt_require_real('mt_xmin',mt_xmin)
      call mt_require_real('mt_xmax',mt_xmax)
      call mt_require_real('mt_zmin',mt_zmin)
      call mt_require_real('mt_zmax',mt_zmax)
      call mt_require_real('mt_y0',mt_y0)
      call mt_require_positive_int('mt_nx',mt_nx)
      call mt_require_positive_int('mt_nz',mt_nz)
      call mt_require_ordered('mt_xmin','mt_xmax',mt_xmin,mt_xmax)
      call mt_require_ordered('mt_zmin','mt_zmax',mt_zmin,mt_zmax)
      write(*,'(a)') 'mt_run_topology_task: writing axis-plane full VTU '//trim(output_file)
      if (mt_b_min>0.d0) then
        call mt_qsl_plane_vtu_xz(mt_xmin,mt_xmax,mt_nx,mt_zmin,mt_zmax, &
             mt_nz,mt_y0,mt_dL,mt_max_steps,trim(output_file), &
             b_min=mt_b_min)
      else
        call mt_qsl_plane_vtu_xz(mt_xmin,mt_xmax,mt_nx,mt_zmin,mt_zmax, &
             mt_nz,mt_y0,mt_dL,mt_max_steps,trim(output_file))
      endif
    case ('yz')
      call mt_require_real('mt_ymin',mt_ymin)
      call mt_require_real('mt_ymax',mt_ymax)
      call mt_require_real('mt_zmin',mt_zmin)
      call mt_require_real('mt_zmax',mt_zmax)
      call mt_require_real('mt_x0',mt_x0)
      call mt_require_positive_int('mt_ny',mt_ny)
      call mt_require_positive_int('mt_nz',mt_nz)
      call mt_require_ordered('mt_ymin','mt_ymax',mt_ymin,mt_ymax)
      call mt_require_ordered('mt_zmin','mt_zmax',mt_zmin,mt_zmax)
      write(*,'(a)') 'mt_run_topology_task: writing axis-plane full VTU '//trim(output_file)
      if (mt_b_min>0.d0) then
        call mt_qsl_plane_vtu_yz(mt_ymin,mt_ymax,mt_ny,mt_zmin,mt_zmax, &
             mt_nz,mt_x0,mt_dL,mt_max_steps,trim(output_file), &
             b_min=mt_b_min)
      else
        call mt_qsl_plane_vtu_yz(mt_ymin,mt_ymax,mt_ny,mt_zmin,mt_zmax, &
             mt_nz,mt_x0,mt_dL,mt_max_steps,trim(output_file))
      endif
    case default
      call mpistop('mt_run_topology_task: axis_plane_full_vtu requires mt_plane=xy, xz, or yz')
    end select
  end subroutine mt_run_axis_plane_full_vtu_task

  subroutine mt_run_axis_plane_csv_task()
    character(len=mt_task_name_len) :: plane
    character(len=mt_task_name_len) :: length_csv,twist_csv,mapping_csv
    character(len=mt_task_name_len) :: qperp_csv

    if (len_trim(mt_output_prefix)==0) then
      call mpistop('axis_plane_csv requires mt_output_prefix')
    endif
    if (len_trim(mt_output_file)>0) then
      write(*,'(a)') 'mt_run_topology_task: axis_plane_csv ignores mt_output_file'
      write(*,'(a)') 'mt_run_topology_task: axis_plane_csv uses mt_output_prefix for CSV outputs'
    endif

    plane=mt_lowercase(trim(mt_plane))
    call mt_resolve_prefix_file('_'//trim(plane)//'_length.csv', &
         length_csv,'axis_plane_csv length requires mt_output_prefix')
    twist_csv=''
    mapping_csv=''
    qperp_csv=''
    if (mt_compute_twist) then
      call mt_resolve_prefix_file('_'//trim(plane)//'_twist.csv', &
           twist_csv,'axis_plane_csv twist requires mt_output_prefix')
    endif
    if (mt_compute_qperp) then
      call mt_resolve_prefix_file('_'//trim(plane)//'_qperp_method2.csv', &
           qperp_csv,'axis_plane_csv Qperp requires mt_output_prefix')
    endif

    select case (trim(plane))
    case ('xy')
      call mt_require_real('mt_xmin',mt_xmin)
      call mt_require_real('mt_xmax',mt_xmax)
      call mt_require_real('mt_ymin',mt_ymin)
      call mt_require_real('mt_ymax',mt_ymax)
      call mt_require_real('mt_z0',mt_z0)
      call mt_require_positive_int('mt_nx',mt_nx)
      call mt_require_positive_int('mt_ny',mt_ny)
      call mt_require_ordered('mt_xmin','mt_xmax',mt_xmin,mt_xmax)
      call mt_require_ordered('mt_ymin','mt_ymax',mt_ymin,mt_ymax)
      if (mt_b_min>0.d0) then
        call mt_topology_plane_xy(mt_xmin,mt_xmax,mt_nx,mt_ymin, &
             mt_ymax,mt_ny,mt_z0,mt_dL,mt_max_steps,trim(length_csv), &
             trim(twist_csv),trim(mapping_csv),b_min=mt_b_min)
        if (mt_compute_qperp) call mt_qperp_plane_xy(mt_xmin,mt_xmax, &
             mt_nx,mt_ymin,mt_ymax,mt_ny,mt_z0,mt_dL,mt_max_steps, &
             trim(qperp_csv),b_min=mt_b_min)
      else
        call mt_topology_plane_xy(mt_xmin,mt_xmax,mt_nx,mt_ymin, &
             mt_ymax,mt_ny,mt_z0,mt_dL,mt_max_steps,trim(length_csv), &
             trim(twist_csv),trim(mapping_csv))
        if (mt_compute_qperp) call mt_qperp_plane_xy(mt_xmin,mt_xmax, &
             mt_nx,mt_ymin,mt_ymax,mt_ny,mt_z0,mt_dL,mt_max_steps, &
             trim(qperp_csv))
      endif
    case ('xz')
      call mt_require_real('mt_xmin',mt_xmin)
      call mt_require_real('mt_xmax',mt_xmax)
      call mt_require_real('mt_zmin',mt_zmin)
      call mt_require_real('mt_zmax',mt_zmax)
      call mt_require_real('mt_y0',mt_y0)
      call mt_require_positive_int('mt_nx',mt_nx)
      call mt_require_positive_int('mt_nz',mt_nz)
      call mt_require_ordered('mt_xmin','mt_xmax',mt_xmin,mt_xmax)
      call mt_require_ordered('mt_zmin','mt_zmax',mt_zmin,mt_zmax)
      if (mt_b_min>0.d0) then
        call mt_topology_plane_xz(mt_xmin,mt_xmax,mt_nx,mt_zmin, &
             mt_zmax,mt_nz,mt_y0,mt_dL,mt_max_steps,trim(length_csv), &
             trim(twist_csv),trim(mapping_csv),b_min=mt_b_min)
        if (mt_compute_qperp) call mt_qperp_plane_xz(mt_xmin,mt_xmax, &
             mt_nx,mt_zmin,mt_zmax,mt_nz,mt_y0,mt_dL,mt_max_steps, &
             trim(qperp_csv),b_min=mt_b_min)
      else
        call mt_topology_plane_xz(mt_xmin,mt_xmax,mt_nx,mt_zmin, &
             mt_zmax,mt_nz,mt_y0,mt_dL,mt_max_steps,trim(length_csv), &
             trim(twist_csv),trim(mapping_csv))
        if (mt_compute_qperp) call mt_qperp_plane_xz(mt_xmin,mt_xmax, &
             mt_nx,mt_zmin,mt_zmax,mt_nz,mt_y0,mt_dL,mt_max_steps, &
             trim(qperp_csv))
      endif
    case ('yz')
      call mt_require_real('mt_ymin',mt_ymin)
      call mt_require_real('mt_ymax',mt_ymax)
      call mt_require_real('mt_zmin',mt_zmin)
      call mt_require_real('mt_zmax',mt_zmax)
      call mt_require_real('mt_x0',mt_x0)
      call mt_require_positive_int('mt_ny',mt_ny)
      call mt_require_positive_int('mt_nz',mt_nz)
      call mt_require_ordered('mt_ymin','mt_ymax',mt_ymin,mt_ymax)
      call mt_require_ordered('mt_zmin','mt_zmax',mt_zmin,mt_zmax)
      if (mt_b_min>0.d0) then
        call mt_topology_plane_yz(mt_ymin,mt_ymax,mt_ny,mt_zmin, &
             mt_zmax,mt_nz,mt_x0,mt_dL,mt_max_steps,trim(length_csv), &
             trim(twist_csv),trim(mapping_csv),b_min=mt_b_min)
        if (mt_compute_qperp) call mt_qperp_plane_yz(mt_ymin,mt_ymax, &
             mt_ny,mt_zmin,mt_zmax,mt_nz,mt_x0,mt_dL,mt_max_steps, &
             trim(qperp_csv),b_min=mt_b_min)
      else
        call mt_topology_plane_yz(mt_ymin,mt_ymax,mt_ny,mt_zmin, &
             mt_zmax,mt_nz,mt_x0,mt_dL,mt_max_steps,trim(length_csv), &
             trim(twist_csv),trim(mapping_csv))
        if (mt_compute_qperp) call mt_qperp_plane_yz(mt_ymin,mt_ymax, &
             mt_ny,mt_zmin,mt_zmax,mt_nz,mt_x0,mt_dL,mt_max_steps, &
             trim(qperp_csv))
      endif
    case default
      call mpistop('mt_run_topology_task: axis_plane_csv requires mt_plane=xy, xz, or yz')
    end select

    write(*,'(a)') 'mt_run_topology_task: wrote axis-plane CSV length '// &
         trim(length_csv)
    if (mt_compute_twist) write(*,'(a)') &
         'mt_run_topology_task: wrote axis-plane CSV twist '//trim(twist_csv)
    if (mt_compute_qperp) write(*,'(a)') &
         'mt_run_topology_task: wrote axis-plane CSV Qperp '//trim(qperp_csv)
  end subroutine mt_run_axis_plane_csv_task

  subroutine mt_run_volume_vti_task()
    character(len=mt_task_name_len) :: output_file

    call mt_resolve_output_file('volume_vti','.vti','_volume.vti', &
         output_file)
    call mt_require_real('mt_xmin',mt_xmin)
    call mt_require_real('mt_xmax',mt_xmax)
    call mt_require_real('mt_ymin',mt_ymin)
    call mt_require_real('mt_ymax',mt_ymax)
    call mt_require_real('mt_zmin',mt_zmin)
    call mt_require_real('mt_zmax',mt_zmax)
    call mt_require_positive_int('mt_nx',mt_nx)
    call mt_require_positive_int('mt_ny',mt_ny)
    call mt_require_positive_int('mt_nz',mt_nz)
    call mt_require_ordered('mt_xmin','mt_xmax',mt_xmin,mt_xmax)
    call mt_require_ordered('mt_ymin','mt_ymax',mt_ymin,mt_ymax)
    call mt_require_ordered('mt_zmin','mt_zmax',mt_zmin,mt_zmax)

    write(*,'(a)') 'mt_run_topology_task: writing volume VTI '//trim(output_file)
    if (mt_b_min>0.d0) then
      if (mt_chunk_nz>0) then
        call mt_fieldline_products_volume_vti(mt_xmin,mt_xmax,mt_nx, &
             mt_ymin,mt_ymax,mt_ny,mt_zmin,mt_zmax,mt_nz, &
             mt_dL,mt_max_steps,trim(output_file),b_min=mt_b_min, &
             compute_twist=mt_compute_twist, &
             compute_qperp=mt_compute_qperp,chunk_nz=mt_chunk_nz)
      else
        call mt_fieldline_products_volume_vti(mt_xmin,mt_xmax,mt_nx, &
             mt_ymin,mt_ymax,mt_ny,mt_zmin,mt_zmax,mt_nz, &
             mt_dL,mt_max_steps,trim(output_file),b_min=mt_b_min, &
             compute_twist=mt_compute_twist, &
             compute_qperp=mt_compute_qperp)
      endif
    else
      if (mt_chunk_nz>0) then
        call mt_fieldline_products_volume_vti(mt_xmin,mt_xmax,mt_nx, &
             mt_ymin,mt_ymax,mt_ny,mt_zmin,mt_zmax,mt_nz, &
             mt_dL,mt_max_steps,trim(output_file), &
             compute_twist=mt_compute_twist, &
             compute_qperp=mt_compute_qperp,chunk_nz=mt_chunk_nz)
      else
        call mt_fieldline_products_volume_vti(mt_xmin,mt_xmax,mt_nx, &
             mt_ymin,mt_ymax,mt_ny,mt_zmin,mt_zmax,mt_nz, &
             mt_dL,mt_max_steps,trim(output_file), &
             compute_twist=mt_compute_twist, &
             compute_qperp=mt_compute_qperp)
      endif
    endif
  end subroutine mt_run_volume_vti_task

  subroutine mt_run_arbitrary_plane_products_task()
    character(len=mt_task_name_len) :: vtu_file,csv_file

    call mt_resolve_output_file('arbitrary_plane_products','.vtu', &
         '_arbitrary_plane_products.vtu',vtu_file)
    csv_file=''
    if (mt_write_csv) then
      call mt_resolve_prefix_file('_arbitrary_plane_products.csv',csv_file, &
           'arbitrary_plane_products CSV requires mt_output_prefix')
    endif

    call mt_require_real('mt_origin(1)',mt_origin(1))
    call mt_require_real('mt_origin(2)',mt_origin(2))
    call mt_require_real('mt_origin(3)',mt_origin(3))
    call mt_require_real('mt_e1(1)',mt_e1(1))
    call mt_require_real('mt_e1(2)',mt_e1(2))
    call mt_require_real('mt_e1(3)',mt_e1(3))
    call mt_require_real('mt_e2(1)',mt_e2(1))
    call mt_require_real('mt_e2(2)',mt_e2(2))
    call mt_require_real('mt_e2(3)',mt_e2(3))
    call mt_require_real('mt_s1min',mt_s1min)
    call mt_require_real('mt_s1max',mt_s1max)
    call mt_require_real('mt_s2min',mt_s2min)
    call mt_require_real('mt_s2max',mt_s2max)
    call mt_require_positive_int('mt_n1',mt_n1)
    call mt_require_positive_int('mt_n2',mt_n2)
    call mt_require_ordered('mt_s1min','mt_s1max',mt_s1min,mt_s1max)
    call mt_require_ordered('mt_s2min','mt_s2max',mt_s2min,mt_s2max)

    write(*,'(a)') 'mt_run_topology_task: writing arbitrary-plane VTU '// &
         trim(vtu_file)
    if (mt_write_csv) then
      write(*,'(a)') 'mt_run_topology_task: writing arbitrary-plane CSV '// &
           trim(csv_file)
    endif

    if (mt_b_min>0.d0) then
      call mt_fieldline_products_plane_arbitrary(mt_origin(1:ndim), &
           mt_e1(1:ndim),mt_e2(1:ndim),mt_s1min,mt_s1max,mt_n1, &
           mt_s2min,mt_s2max,mt_n2,mt_dL,mt_max_steps,trim(csv_file), &
           b_min=mt_b_min,compute_twist=mt_compute_twist, &
           compute_qperp=mt_compute_qperp,vtu_file=trim(vtu_file), &
           write_csv=mt_write_csv)
    else
      call mt_fieldline_products_plane_arbitrary(mt_origin(1:ndim), &
           mt_e1(1:ndim),mt_e2(1:ndim),mt_s1min,mt_s1max,mt_n1, &
           mt_s2min,mt_s2max,mt_n2,mt_dL,mt_max_steps,trim(csv_file), &
           compute_twist=mt_compute_twist, &
           compute_qperp=mt_compute_qperp,vtu_file=trim(vtu_file), &
           write_csv=mt_write_csv)
    endif
  end subroutine mt_run_arbitrary_plane_products_task

  subroutine mt_run_seed_products_task()
    character(len=mt_task_name_len) :: csv_file,vtu_file
    double precision, allocatable :: seeds(:,:)
    integer :: nseed

    if (len_trim(mt_output_prefix)==0) then
      call mpistop('seed_products requires mt_output_prefix')
    endif
    if (len_trim(mt_output_file)>0) then
      write(*,'(a)') 'mt_run_topology_task: seed_products ignores mt_output_file'
      write(*,'(a)') 'mt_run_topology_task: seed_products uses mt_output_prefix for CSV and VTU'
    endif
    if (len_trim(mt_seed_file)==0) then
      call mpistop('seed_products requires mt_seed_file')
    endif

    call mt_resolve_prefix_file('_seed_products.csv',csv_file, &
         'seed_products CSV requires mt_output_prefix')
    call mt_resolve_prefix_file('_seed_products.vtu',vtu_file, &
         'seed_products VTU requires mt_output_prefix')
    call mt_read_seed_file(trim(mt_seed_file),seeds,nseed)

    write(*,'(a)') 'mt_run_topology_task: writing seed-products CSV '// &
         trim(csv_file)
    write(*,'(a)') 'mt_run_topology_task: writing seed-products VTU '// &
         trim(vtu_file)

    if (mt_write_csv) then
      write(*,'(a)') 'mt_run_topology_task: seed_products always writes CSV'
    endif

    if (mt_b_min>0.d0) then
      call mt_fieldline_products_seeds(seeds,nseed,mt_dL,mt_max_steps, &
           trim(csv_file),b_min=mt_b_min,compute_twist=mt_compute_twist, &
           compute_qperp=mt_compute_qperp,vtu_file=trim(vtu_file))
    else
      call mt_fieldline_products_seeds(seeds,nseed,mt_dL,mt_max_steps, &
           trim(csv_file),compute_twist=mt_compute_twist, &
           compute_qperp=mt_compute_qperp,vtu_file=trim(vtu_file))
    endif

    deallocate(seeds)
  end subroutine mt_run_seed_products_task

  subroutine mt_resolve_output_file(mode,extension,suffix,output_file)
    character(len=*), intent(in) :: mode,extension,suffix
    character(len=*), intent(out) :: output_file
    integer :: required_len

    output_file=''
    if (len_trim(mt_output_file)>0) then
      output_file=trim(mt_output_file)
      if (len_trim(mt_output_prefix)>0) then
        write(*,'(a)') 'mt_run_topology_task: mt_output_file takes precedence over mt_output_prefix'
      endif
      call mt_warn_extension(output_file,extension,mode)
    else
      required_len=len_trim(mt_output_prefix)+len_trim(suffix)
      if (required_len>len(output_file)) then
        call mpistop('mt_run_topology_task output filename exceeds internal length')
      endif
      output_file=trim(mt_output_prefix)//trim(suffix)
    endif
  end subroutine mt_resolve_output_file

  subroutine mt_resolve_prefix_file(suffix,output_file,missing_message)
    character(len=*), intent(in) :: suffix,missing_message
    character(len=*), intent(out) :: output_file
    integer :: required_len

    if (len_trim(mt_output_prefix)==0) then
      call mpistop(trim(missing_message))
    endif
    required_len=len_trim(mt_output_prefix)+len_trim(suffix)
    if (required_len>len(output_file)) then
      call mpistop('mt_run_topology_task output filename exceeds internal length')
    endif
    output_file=trim(mt_output_prefix)//trim(suffix)
  end subroutine mt_resolve_prefix_file

  logical function mt_is_set_real(value) result(is_set)
    double precision, intent(in) :: value

    is_set=(value/=mt_unset_real)
  end function mt_is_set_real

  subroutine mt_require_real(name,value)
    character(len=*), intent(in) :: name
    double precision, intent(in) :: value

    if (.not.mt_is_set_real(value)) then
      call mpistop('mt_run_topology_task requires '//trim(name))
    endif
  end subroutine mt_require_real

  subroutine mt_require_positive_int(name,value)
    character(len=*), intent(in) :: name
    integer, intent(in) :: value

    if (value<=0) then
      call mpistop('mt_run_topology_task requires '//trim(name)//' > 0')
    endif
  end subroutine mt_require_positive_int

  subroutine mt_require_ordered(name_min,name_max,value_min,value_max)
    character(len=*), intent(in) :: name_min,name_max
    double precision, intent(in) :: value_min,value_max

    if (value_max<value_min) then
      call mpistop('mt_run_topology_task requires ordered '// &
           trim(name_min)//'/'//trim(name_max))
    endif
  end subroutine mt_require_ordered

  subroutine mt_read_seed_file(filename,seeds,nseed)
    character(len=*), intent(in) :: filename
    double precision, allocatable, intent(out) :: seeds(:,:)
    integer, intent(out) :: nseed

    character(len=1024) :: line
    double precision :: xyz(3)
    integer :: io,ios,line_no,iseed

    nseed=0
    open(newunit=io,file=trim(filename),status='old',action='read', &
         form='formatted',iostat=ios)
    if (ios/=0) then
      call mpistop('mt_read_seed_file could not open '//trim(filename))
    endif
    line_no=0
    do
      read(io,'(a)',iostat=ios) line
      if (ios/=0) exit
      line_no=line_no+1
      if (.not.mt_seed_line_is_data(line)) cycle
      call mt_validate_seed_line(filename,line,line_no)
      nseed=nseed+1
    enddo
    close(io)
    if (nseed<=0) then
      call mpistop('mt_read_seed_file found no seeds in '//trim(filename))
    endif

    allocate(seeds(nseed,ndim))
    seeds=0.d0
    open(newunit=io,file=trim(filename),status='old',action='read', &
         form='formatted',iostat=ios)
    if (ios/=0) then
      call mpistop('mt_read_seed_file could not reopen '//trim(filename))
    endif
    line_no=0
    iseed=0
    do
      read(io,'(a)',iostat=ios) line
      if (ios/=0) exit
      line_no=line_no+1
      if (.not.mt_seed_line_is_data(line)) cycle
      iseed=iseed+1
      read(line,*,iostat=ios) xyz
      if (ios/=0) then
        call mt_fail_seed_line(filename,line_no)
      endif
      seeds(iseed,1:ndim)=xyz(1:ndim)
    enddo
    close(io)
  end subroutine mt_read_seed_file

  logical function mt_seed_line_is_data(line) result(is_data)
    character(len=*), intent(in) :: line
    character(len=len(line)) :: trimmed

    trimmed=adjustl(line)
    is_data=len_trim(trimmed)>0
    if (is_data) is_data=trimmed(1:1)/='#'
  end function mt_seed_line_is_data

  subroutine mt_validate_seed_line(filename,line,line_no)
    character(len=*), intent(in) :: filename,line
    integer, intent(in) :: line_no

    double precision :: x,y,z
    integer :: ios

    if (index(line,',')>0) then
      call mt_fail_seed_line(filename,line_no)
    endif
    if (mt_count_tokens(line)/=3) then
      call mt_fail_seed_line(filename,line_no)
    endif
    read(line,*,iostat=ios) x,y,z
    if (ios/=0) then
      call mt_fail_seed_line(filename,line_no)
    endif
  end subroutine mt_validate_seed_line

  integer function mt_count_tokens(line) result(ntoken)
    character(len=*), intent(in) :: line
    integer :: i
    logical :: in_token

    ntoken=0
    in_token=.false.
    do i=1,len_trim(line)
      if (line(i:i)==' ' .or. line(i:i)==achar(9)) then
        in_token=.false.
      else
        if (.not.in_token) then
          ntoken=ntoken+1
          in_token=.true.
        endif
      endif
    enddo
  end function mt_count_tokens

  subroutine mt_fail_seed_line(filename,line_no)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: line_no

    call mpistop('mt_read_seed_file malformed line '// &
         trim(mt_int_to_string(line_no))//' in '//trim(filename))
  end subroutine mt_fail_seed_line

  function mt_int_to_string(value) result(text)
    integer, intent(in) :: value
    character(len=32) :: text

    write(text,'(i0)') value
  end function mt_int_to_string

  subroutine mt_warn_extension(filename,extension,mode)
    character(len=*), intent(in) :: filename,extension,mode

    if (.not.mt_has_extension(filename,extension)) then
      write(*,'(a)') 'mt_run_topology_task warning: '//trim(mode)// &
           ' output file does not end in '//trim(extension)
    endif
  end subroutine mt_warn_extension

  logical function mt_has_extension(filename,extension) result(has_ext)
    character(len=*), intent(in) :: filename,extension
    character(len=mt_task_name_len) :: fname,ext
    integer :: lf,le

    fname=mt_lowercase(trim(filename))
    ext=mt_lowercase(trim(extension))
    lf=len_trim(fname)
    le=len_trim(ext)
    has_ext=.false.
    if (lf>=le) has_ext=(fname(lf-le+1:lf)==ext(1:le))
  end function mt_has_extension

  function mt_lowercase(input) result(output)
    character(len=*), intent(in) :: input
    character(len=len(input)) :: output
    integer :: i,code

    output=input
    do i=1,len(input)
      code=iachar(output(i:i))
      if (code>=iachar('A') .and. code<=iachar('Z')) then
        output(i:i)=achar(code+iachar('a')-iachar('A'))
      endif
    enddo
  end function mt_lowercase

  logical function mt_vtk_detail_is_full() result(is_full)
    is_full=(mt_lowercase(trim(mt_vtk_detail))=='full')
  end function mt_vtk_detail_is_full

  double precision function mt_visual_float(value) result(out_value)
    double precision, intent(in) :: value

    if (ieee_is_finite(value)) then
      out_value=value
    else
      out_value=0.d0
    endif
  end function mt_visual_float

  double precision function mt_visual_valid_float(value,is_valid) &
       result(out_value)
    double precision, intent(in) :: value
    logical, intent(in) :: is_valid

    if (is_valid .and. ieee_is_finite(value)) then
      out_value=value
    else
      out_value=0.d0
    endif
  end function mt_visual_valid_float

  subroutine mt_length_single(seed,dL,max_steps,csv_file,b_min)
    ! Trace one magnetic field line and write its summary to a CSV file.
    double precision, intent(in) :: seed(ndim),dL
    integer, intent(in) :: max_steps
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    type(trace_length_result) :: result
    double precision :: seed_xyz(3)
    integer :: csv_unit,io_status

    if (npe/=1) then
      call mpistop('mt_length_single currently requires npe=1')
    endif

    if (present(b_min)) then
      call trace_field_length_single(seed,dL,max_steps,result,b_min)
    else
      call trace_field_length_single(seed,dL,max_steps,result)
    endif

    seed_xyz=0.d0
    seed_xyz(1:ndim)=result%seed

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop('mt_length_single could not open CSV file')
    endif

    write(csv_unit,'(a)',iostat=io_status) &
         'seed_x,seed_y,seed_z,'// &
         'length_total,length_backward,length_forward,'// &
         'nstep_backward,nstep_forward,'// &
         'status_backward,status_forward'
    if (io_status/=0) then
      close(csv_unit)
      call mpistop('mt_length_single could not write CSV header')
    endif

    write(csv_unit,'(es24.16,5(",",es24.16),4(",",i0))',iostat=io_status) &
         seed_xyz,result%total_length, &
         result%backward_length,result%forward_length, &
         result%backward_nstep,result%forward_nstep, &
         result%backward_status,result%forward_status
    if (io_status/=0) then
      close(csv_unit)
      call mpistop('mt_length_single could not write CSV data')
    endif

    close(csv_unit)
  end subroutine mt_length_single

  subroutine mt_length_seeds(seeds,nseed,dL,max_steps,csv_file,b_min)
    ! Trace multiple magnetic field lines and write one summary row per seed.
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    type(trace_length_result), allocatable :: results(:)
    double precision :: seed_xyz(3)
    integer :: csv_unit,io_status,iseed

    if (npe/=1) then
      call mpistop('mt_length_seeds currently requires npe=1')
    endif
    if (nseed<0) then
      call mpistop('mt_length_seeds requires nseed>=0')
    endif

    allocate(results(nseed))
    if (present(b_min)) then
      call trace_field_length_multi(seeds,nseed,dL,max_steps,results,b_min)
    else
      call trace_field_length_multi(seeds,nseed,dL,max_steps,results)
    endif

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      deallocate(results)
      call mpistop('mt_length_seeds could not open CSV file')
    endif

    write(csv_unit,'(a)',iostat=io_status) &
         'seed_id,seed_x,seed_y,seed_z,'// &
         'length_total,length_backward,length_forward,'// &
         'nstep_backward,nstep_forward,'// &
         'status_backward,status_forward'
    if (io_status/=0) then
      close(csv_unit)
      deallocate(results)
      call mpistop('mt_length_seeds could not write CSV header')
    endif

    do iseed=1,nseed
      seed_xyz=0.d0
      seed_xyz(1:ndim)=results(iseed)%seed
      write(csv_unit,'(i0,6(",",es24.16),4(",",i0))',iostat=io_status) &
           iseed,seed_xyz,results(iseed)%total_length, &
           results(iseed)%backward_length,results(iseed)%forward_length, &
           results(iseed)%backward_nstep,results(iseed)%forward_nstep, &
           results(iseed)%backward_status,results(iseed)%forward_status
      if (io_status/=0) then
        close(csv_unit)
        deallocate(results)
        call mpistop('mt_length_seeds could not write CSV data')
      endif
    enddo

    close(csv_unit)
    deallocate(results)
  end subroutine mt_length_seeds

  subroutine mt_fieldline_products_seeds(seeds,nseed,dL,max_steps,csv_file, &
       b_min,compute_twist,compute_qperp,vtu_file)
    ! Write selected per-field-line diagnostics for an arbitrary seed set.
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min
    logical, intent(in), optional :: compute_twist,compute_qperp
    character(len=*), intent(in), optional :: vtu_file

    type(trace_length_result), allocatable :: length_results(:)
    type(trace_twist_result), allocatable :: twist_results(:)
    type(trace_qperp_result), allocatable :: qperp_results(:)
    logical :: do_twist,do_qperp

    if (npe/=1) then
      call mpistop('mt_fieldline_products_seeds currently requires npe=1')
    endif
    if (ndim/=3 .or. .not.slab_uniform) then
      call mpistop('mt_fieldline_products_seeds requires 3D uniform Cartesian geometry')
    endif
    if (nseed<0) then
      call mpistop('mt_fieldline_products_seeds requires nseed>=0')
    endif

    do_twist=.false.
    do_qperp=.false.
    if (present(compute_twist)) do_twist=compute_twist
    if (present(compute_qperp)) do_qperp=compute_qperp

    allocate(length_results(nseed))
    if (do_twist) then
      allocate(twist_results(nseed))
    else
      allocate(twist_results(0))
    endif
    if (do_qperp) then
      allocate(qperp_results(nseed))
    else
      allocate(qperp_results(0))
    endif

    if (present(b_min)) then
      call mt_trace_fieldline_products_seedset(seeds,nseed,dL,max_steps, &
           length_results,twist_results,qperp_results,do_twist,do_qperp, &
           b_min)
    else
      call mt_trace_fieldline_products_seedset(seeds,nseed,dL,max_steps, &
           length_results,twist_results,qperp_results,do_twist,do_qperp)
    endif

    call mt_write_fieldline_products_seeds_csv(length_results, &
         twist_results,qperp_results,nseed,csv_file,do_twist,do_qperp)
    if (present(vtu_file)) then
      if (len_trim(vtu_file)>0) then
        call mt_write_fieldline_products_vtu_vertices(vtu_file, &
             length_results,twist_results,qperp_results,nseed,do_twist, &
             do_qperp,'mt_fieldline_products_seeds')
      endif
    endif

    if (allocated(qperp_results)) deallocate(qperp_results)
    if (allocated(twist_results)) deallocate(twist_results)
    deallocate(length_results)
  end subroutine mt_fieldline_products_seeds

  subroutine mt_length_plane_xy(xmin,xmax,nx,ymin,ymax,ny,z0,dL, &
       max_steps,csv_file,b_min)
    ! Trace a uniform seed grid on a constant-z Cartesian plane.
    integer, intent(in) :: nx,ny,max_steps
    double precision, intent(in) :: xmin,xmax,ymin,ymax,z0,dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_length_plane_axis(xmin,xmax,nx,ymin,ymax,ny,z0,1,2,3, &
           dL,max_steps,csv_file,'mt_length_plane_xy','ix,iy',b_min)
    else
      call mt_length_plane_axis(xmin,xmax,nx,ymin,ymax,ny,z0,1,2,3, &
           dL,max_steps,csv_file,'mt_length_plane_xy','ix,iy')
    endif
  end subroutine mt_length_plane_xy

  subroutine mt_mapping_plane_xy(xmin,xmax,nx,ymin,ymax,ny,z0,dL, &
       max_steps,csv_file,b_min)
    ! Trace a constant-z seed grid and write endpoint metadata for Q mapping.
    integer, intent(in) :: nx,ny,max_steps
    double precision, intent(in) :: xmin,xmax,ymin,ymax,z0,dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    type(trace_mapping_result), allocatable :: results(:)
    double precision, allocatable :: seeds(:,:)
    double precision :: source_normal(3)
    integer :: nseed

    call mt_validate_axis_plane(xmin,xmax,nx,ymin,ymax,ny, &
         'mt_mapping_plane_xy')
    call mt_build_axis_plane_seeds(xmin,xmax,nx,ymin,ymax,ny,z0, &
         1,2,3,seeds)

    source_normal=0.d0
    source_normal(3)=1.d0
    nseed=nx*ny
    allocate(results(nseed))
    if (present(b_min)) then
      call trace_field_mapping_multi(seeds,nseed,dL,max_steps,results, &
           b_min=b_min,source_normal=source_normal)
    else
      call trace_field_mapping_multi(seeds,nseed,dL,max_steps,results, &
           source_normal=source_normal)
    endif
    call mt_write_mapping_plane_xy_csv(results,nx,ny,csv_file)

    deallocate(seeds,results)
  end subroutine mt_mapping_plane_xy

  subroutine mt_qperp_plane_xy(xmin,xmax,nx,ymin,ymax,ny,z0,dL, &
       max_steps,csv_file,b_min)
    ! Compute Method-II Qperp diagnostics on a constant-z source plane.
    integer, intent(in) :: nx,ny,max_steps
    double precision, intent(in) :: xmin,xmax,ymin,ymax,z0,dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_qperp_plane_axis(xmin,xmax,nx,ymin,ymax,ny,z0,1,2,3, &
           dL,max_steps,csv_file,'mt_qperp_plane_xy','ix,iy',b_min)
    else
      call mt_qperp_plane_axis(xmin,xmax,nx,ymin,ymax,ny,z0,1,2,3, &
           dL,max_steps,csv_file,'mt_qperp_plane_xy','ix,iy')
    endif
  end subroutine mt_qperp_plane_xy

  subroutine mt_qperp_plane_xz(xmin,xmax,nx,zmin,zmax,nz,y0,dL, &
       max_steps,csv_file,b_min)
    ! Compute Method-II Qperp diagnostics on a constant-y source plane.
    integer, intent(in) :: nx,nz,max_steps
    double precision, intent(in) :: xmin,xmax,zmin,zmax,y0,dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_qperp_plane_axis(xmin,xmax,nx,zmin,zmax,nz,y0,1,3,2, &
           dL,max_steps,csv_file,'mt_qperp_plane_xz','ix,iz',b_min)
    else
      call mt_qperp_plane_axis(xmin,xmax,nx,zmin,zmax,nz,y0,1,3,2, &
           dL,max_steps,csv_file,'mt_qperp_plane_xz','ix,iz')
    endif
  end subroutine mt_qperp_plane_xz

  subroutine mt_qperp_plane_yz(ymin,ymax,ny,zmin,zmax,nz,x0,dL, &
       max_steps,csv_file,b_min)
    ! Compute Method-II Qperp diagnostics on a constant-x source plane.
    integer, intent(in) :: ny,nz,max_steps
    double precision, intent(in) :: ymin,ymax,zmin,zmax,x0,dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_qperp_plane_axis(ymin,ymax,ny,zmin,zmax,nz,x0,2,3,1, &
           dL,max_steps,csv_file,'mt_qperp_plane_yz','iy,iz',b_min)
    else
      call mt_qperp_plane_axis(ymin,ymax,ny,zmin,zmax,nz,x0,2,3,1, &
           dL,max_steps,csv_file,'mt_qperp_plane_yz','iy,iz')
    endif
  end subroutine mt_qperp_plane_yz

  subroutine mt_qperp_plane_arbitrary(origin,e1,e2,s1min,s1max,n1, &
       s2min,s2max,n2,dL,max_steps,csv_file,b_min)
    ! Compute Method-II Qperp diagnostics on an orthonormal arbitrary plane.
    integer, intent(in) :: n1,n2,max_steps
    double precision, intent(in) :: origin(ndim),e1(ndim),e2(ndim)
    double precision, intent(in) :: s1min,s1max,s2min,s2max,dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    type(trace_qperp_result), allocatable :: results(:)
    double precision, allocatable :: seeds(:,:),s1(:),s2(:)
    integer :: nseed

    if (.not.mt_validate_arbitrary_plane_basis(e1,e2,s1min,s1max,n1, &
         s2min,s2max,n2,'mt_qperp_plane_arbitrary')) then
      call mt_write_qperp_arbitrary_header(csv_file, &
           'mt_qperp_plane_arbitrary')
      return
    endif

    call mt_build_arbitrary_plane_seeds(origin,e1,e2,s1min,s1max,n1, &
         s2min,s2max,n2,seeds,s1,s2)

    nseed=n1*n2
    allocate(results(nseed))
    if (present(b_min)) then
      call trace_field_qperp_multi(seeds,nseed,dL,max_steps,results,b_min)
    else
      call trace_field_qperp_multi(seeds,nseed,dL,max_steps,results)
    endif
    call mt_write_qperp_arbitrary_csv(results,s1,s2,n1,n2,csv_file, &
         'mt_qperp_plane_arbitrary')

    deallocate(seeds,s1,s2,results)
  end subroutine mt_qperp_plane_arbitrary

  subroutine mt_fieldline_products_plane_arbitrary(origin,e1,e2,s1min, &
       s1max,n1,s2min,s2max,n2,dL,max_steps,csv_file,b_min, &
       compute_twist,compute_qperp,vtu_file,write_csv)
    ! Write selected per-field-line diagnostics on an arbitrary seed plane.
    integer, intent(in) :: n1,n2,max_steps
    double precision, intent(in) :: origin(ndim),e1(ndim),e2(ndim)
    double precision, intent(in) :: s1min,s1max,s2min,s2max,dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min
    logical, intent(in), optional :: compute_twist,compute_qperp
    character(len=*), intent(in), optional :: vtu_file
    logical, intent(in), optional :: write_csv

    type(trace_length_result), allocatable :: length_results(:)
    type(trace_twist_result), allocatable :: twist_results(:)
    type(trace_qperp_result), allocatable :: qperp_results(:)
    double precision, allocatable :: seeds(:,:),s1(:),s2(:)
    integer :: nseed
    logical :: do_twist,do_qperp,do_csv

    do_twist=.false.
    do_qperp=.false.
    do_csv=.true.
    if (present(compute_twist)) do_twist=compute_twist
    if (present(compute_qperp)) do_qperp=compute_qperp
    if (present(write_csv)) do_csv=write_csv

    if (.not.mt_validate_arbitrary_plane_basis(e1,e2,s1min,s1max,n1, &
         s2min,s2max,n2,'mt_fieldline_products_plane_arbitrary')) then
      if (do_csv) then
        call mt_write_fieldline_products_plane_arbitrary_header(csv_file, &
             'mt_fieldline_products_plane_arbitrary',do_twist,do_qperp)
      endif
      if (present(vtu_file)) then
        if (len_trim(vtu_file)>0) then
          call mt_write_fieldline_products_vtu_empty(vtu_file,do_twist, &
               do_qperp,'mt_fieldline_products_plane_arbitrary')
        endif
      endif
      return
    endif

    call mt_build_arbitrary_plane_seeds(origin,e1,e2,s1min,s1max,n1, &
         s2min,s2max,n2,seeds,s1,s2)

    nseed=n1*n2
    allocate(length_results(nseed))
    if (do_twist) then
      allocate(twist_results(nseed))
    else
      allocate(twist_results(0))
    endif
    if (do_qperp) then
      allocate(qperp_results(nseed))
    else
      allocate(qperp_results(0))
    endif

    if (present(b_min)) then
      call mt_trace_fieldline_products_seedset(seeds,nseed,dL,max_steps, &
           length_results,twist_results,qperp_results,do_twist,do_qperp, &
           b_min)
    else
      call mt_trace_fieldline_products_seedset(seeds,nseed,dL,max_steps, &
           length_results,twist_results,qperp_results,do_twist,do_qperp)
    endif

    if (do_csv) then
      call mt_write_fieldline_products_plane_arbitrary_csv(length_results, &
           twist_results,qperp_results,s1,s2,n1,n2,csv_file,do_twist, &
           do_qperp,'mt_fieldline_products_plane_arbitrary')
    endif
    if (present(vtu_file)) then
      if (len_trim(vtu_file)>0) then
        call mt_write_fieldline_products_vtu_plane(vtu_file, &
             length_results,twist_results,qperp_results,n1,n2,do_twist, &
             do_qperp,'mt_fieldline_products_plane_arbitrary')
      endif
    endif

    deallocate(qperp_results,twist_results,length_results,seeds,s1,s2)
  end subroutine mt_fieldline_products_plane_arbitrary

  subroutine mt_topology_plane_xy(xmin,xmax,nx,ymin,ymax,ny,z0,dL, &
       max_steps,length_csv,twist_csv,mapping_csv,b_min)
    ! Trace a constant-z seed plane once and write selected topology products.
    integer, intent(in) :: nx,ny,max_steps
    double precision, intent(in) :: xmin,xmax,ymin,ymax,z0,dL
    character(len=*), intent(in) :: length_csv,twist_csv,mapping_csv
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_topology_plane_axis(xmin,xmax,nx,ymin,ymax,ny,z0,1,2,3, &
           dL,max_steps,length_csv,twist_csv,mapping_csv, &
           'mt_topology_plane_xy','ix,iy','ix,iy',b_min)
    else
      call mt_topology_plane_axis(xmin,xmax,nx,ymin,ymax,ny,z0,1,2,3, &
           dL,max_steps,length_csv,twist_csv,mapping_csv, &
           'mt_topology_plane_xy','ix,iy','ix,iy')
    endif
  end subroutine mt_topology_plane_xy

  subroutine mt_topology_plane_xz(xmin,xmax,nx,zmin,zmax,nz,y0,dL, &
       max_steps,length_csv,twist_csv,mapping_csv,b_min)
    ! Trace a constant-y seed plane once and write selected topology products.
    integer, intent(in) :: nx,nz,max_steps
    double precision, intent(in) :: xmin,xmax,zmin,zmax,y0,dL
    character(len=*), intent(in) :: length_csv,twist_csv,mapping_csv
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_topology_plane_axis(xmin,xmax,nx,zmin,zmax,nz,y0,1,3,2, &
           dL,max_steps,length_csv,twist_csv,mapping_csv, &
           'mt_topology_plane_xz','ix,iz','i,j',b_min)
    else
      call mt_topology_plane_axis(xmin,xmax,nx,zmin,zmax,nz,y0,1,3,2, &
           dL,max_steps,length_csv,twist_csv,mapping_csv, &
           'mt_topology_plane_xz','ix,iz','i,j')
    endif
  end subroutine mt_topology_plane_xz

  subroutine mt_topology_plane_yz(ymin,ymax,ny,zmin,zmax,nz,x0,dL, &
       max_steps,length_csv,twist_csv,mapping_csv,b_min)
    ! Trace a constant-x seed plane once and write selected topology products.
    integer, intent(in) :: ny,nz,max_steps
    double precision, intent(in) :: ymin,ymax,zmin,zmax,x0,dL
    character(len=*), intent(in) :: length_csv,twist_csv,mapping_csv
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_topology_plane_axis(ymin,ymax,ny,zmin,zmax,nz,x0,2,3,1, &
           dL,max_steps,length_csv,twist_csv,mapping_csv, &
           'mt_topology_plane_yz','iy,iz','i,j',b_min)
    else
      call mt_topology_plane_axis(ymin,ymax,ny,zmin,zmax,nz,x0,2,3,1, &
           dL,max_steps,length_csv,twist_csv,mapping_csv, &
           'mt_topology_plane_yz','iy,iz','i,j')
    endif
  end subroutine mt_topology_plane_yz

  subroutine mt_qsl_plane_vtu_xy(xmin,xmax,nx,ymin,ymax,ny,z0,dL, &
       max_steps,vtu_file,b_min)
    ! Write a full topology/QSL visualization file for a constant-z plane.
    integer, intent(in) :: nx,ny,max_steps
    double precision, intent(in) :: xmin,xmax,ymin,ymax,z0,dL
    character(len=*), intent(in) :: vtu_file
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_qsl_plane_vtu_axis(xmin,xmax,nx,ymin,ymax,ny,z0,1,2,3, &
           dL,max_steps,vtu_file,'mt_qsl_plane_vtu_xy',b_min)
    else
      call mt_qsl_plane_vtu_axis(xmin,xmax,nx,ymin,ymax,ny,z0,1,2,3, &
           dL,max_steps,vtu_file,'mt_qsl_plane_vtu_xy')
    endif
  end subroutine mt_qsl_plane_vtu_xy

  subroutine mt_qsl_plane_vtu_xz(xmin,xmax,nx,zmin,zmax,nz,y0,dL, &
       max_steps,vtu_file,b_min)
    ! Write a full topology/QSL visualization file for a constant-y plane.
    integer, intent(in) :: nx,nz,max_steps
    double precision, intent(in) :: xmin,xmax,zmin,zmax,y0,dL
    character(len=*), intent(in) :: vtu_file
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_qsl_plane_vtu_axis(xmin,xmax,nx,zmin,zmax,nz,y0,1,3,2, &
           dL,max_steps,vtu_file,'mt_qsl_plane_vtu_xz',b_min)
    else
      call mt_qsl_plane_vtu_axis(xmin,xmax,nx,zmin,zmax,nz,y0,1,3,2, &
           dL,max_steps,vtu_file,'mt_qsl_plane_vtu_xz')
    endif
  end subroutine mt_qsl_plane_vtu_xz

  subroutine mt_qsl_plane_vtu_yz(ymin,ymax,ny,zmin,zmax,nz,x0,dL, &
       max_steps,vtu_file,b_min)
    ! Write a full topology/QSL visualization file for a constant-x plane.
    integer, intent(in) :: ny,nz,max_steps
    double precision, intent(in) :: ymin,ymax,zmin,zmax,x0,dL
    character(len=*), intent(in) :: vtu_file
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_qsl_plane_vtu_axis(ymin,ymax,ny,zmin,zmax,nz,x0,2,3,1, &
           dL,max_steps,vtu_file,'mt_qsl_plane_vtu_yz',b_min)
    else
      call mt_qsl_plane_vtu_axis(ymin,ymax,ny,zmin,zmax,nz,x0,2,3,1, &
           dL,max_steps,vtu_file,'mt_qsl_plane_vtu_yz')
    endif
  end subroutine mt_qsl_plane_vtu_yz

  subroutine mt_write_cartesian_vti_pointdata(vti_file,xmin,xmax,nx, &
       ymin,ymax,ny,zmin,zmax,nz,length_total,qperp,status)
    ! Low-level appended-binary VTI writer for future uniform seed-volume
    ! products. PointData values use i-fastest ordering:
    ! ipoint = (k-1)*nx*ny + (j-1)*nx + i.
    character(len=*), intent(in) :: vti_file
    integer, intent(in) :: nx,ny,nz
    double precision, intent(in) :: xmin,xmax,ymin,ymax,zmin,zmax
    double precision, intent(in) :: length_total(:),qperp(:)
    integer, intent(in) :: status(:)

    integer :: npoint
    double precision :: origin(3),spacing(3)

    if (len_trim(vti_file)==0) then
      call mpistop('mt_write_cartesian_vti_pointdata requires a VTI file')
    endif
    if (nx<=0 .or. ny<=0 .or. nz<=0) then
      call mpistop('mt_write_cartesian_vti_pointdata requires positive dimensions')
    endif

    npoint=nx*ny*nz
    if (size(length_total)/=npoint .or. size(qperp)/=npoint .or. &
         size(status)/=npoint) then
      call mpistop('mt_write_cartesian_vti_pointdata array size mismatch')
    endif

    origin=(/ xmin,ymin,zmin /)
    spacing(1)=mt_vti_axis_spacing(xmin,xmax,nx)
    spacing(2)=mt_vti_axis_spacing(ymin,ymax,ny)
    spacing(3)=mt_vti_axis_spacing(zmin,zmax,nz)
    call mt_write_vti_pointdata_fixed(vti_file,origin,spacing,nx,ny,nz, &
         length_total,qperp,status,'mt_write_cartesian_vti_pointdata')
  end subroutine mt_write_cartesian_vti_pointdata

  subroutine mt_fieldline_products_volume_vti(xmin,xmax,nx, &
       ymin,ymax,ny,zmin,zmax,nz,dL,max_steps,vti_file,b_min, &
       compute_twist,compute_qperp,chunk_nz)
    ! Compute per-seed field-line products on a user-defined uniform
    ! Cartesian sampling volume and write PointData to appended-binary VTI.
    integer, intent(in) :: nx,ny,nz,max_steps
    double precision, intent(in) :: xmin,xmax,ymin,ymax,zmin,zmax,dL
    character(len=*), intent(in) :: vti_file
    double precision, intent(in), optional :: b_min
    logical, intent(in), optional :: compute_twist,compute_qperp
    integer, intent(in), optional :: chunk_nz

    type(trace_length_result), allocatable :: length_slab(:)
    type(trace_twist_result), allocatable :: twist_slab(:)
    type(trace_qperp_result), allocatable :: qperp_slab(:)
    type(mt_volume_products) :: products
    double precision, allocatable :: seeds_slab(:,:)
    double precision :: origin(3),spacing(3)
    integer :: k_start,k_end,slab_nz,slab_nseed,nseed
    integer :: chunk_nz_eff
    integer(kind=8) :: nseed64
    logical :: do_twist,do_qperp

    if (len_trim(vti_file)==0) then
      call mpistop('mt_fieldline_products_volume_vti requires a VTI file')
    endif
    if (npe/=1) then
      call mpistop('mt_fieldline_products_volume_vti currently requires npe=1')
    endif
    if (ndim/=3 .or. .not.slab_uniform) then
      call mpistop('mt_fieldline_products_volume_vti requires 3D uniform Cartesian geometry')
    endif
    if (nx<=0 .or. ny<=0 .or. nz<=0) then
      call mpistop('mt_fieldline_products_volume_vti requires positive dimensions')
    endif
    if (xmax<xmin .or. ymax<ymin .or. zmax<zmin) then
      call mpistop('mt_fieldline_products_volume_vti requires ordered bounds')
    endif

    nseed64=int(nx,kind=8)*int(ny,kind=8)*int(nz,kind=8)
    if (nseed64>int(huge(nseed),kind=8)) then
      call mpistop('mt_fieldline_products_volume_vti seed count overflows integer')
    endif
    if (nseed64>1000000_8) then
      write(*,'(a,i0,a)') &
           'mt_fieldline_products_volume_vti warning: volume has ', &
           nseed64,' seeds; final VTI arrays still scale with nseed'
    endif
    nseed=int(nseed64)

    do_twist=.false.
    do_qperp=.false.
    if (present(compute_twist)) do_twist=compute_twist
    if (present(compute_qperp)) do_qperp=compute_qperp

    chunk_nz_eff=nz
    if (present(chunk_nz)) chunk_nz_eff=max(1,min(nz,chunk_nz))

    origin=(/ xmin,ymin,zmin /)
    spacing(1)=mt_vti_axis_spacing(xmin,xmax,nx)
    spacing(2)=mt_vti_axis_spacing(ymin,ymax,ny)
    spacing(3)=mt_vti_axis_spacing(zmin,zmax,nz)

    call mt_allocate_volume_products(products,nseed,do_twist,do_qperp)

    do k_start=1,nz,chunk_nz_eff
      k_end=min(nz,k_start+chunk_nz_eff-1)
      slab_nz=k_end-k_start+1
      slab_nseed=nx*ny*slab_nz

      allocate(seeds_slab(slab_nseed,ndim))
      call mt_build_volume_slab_seeds(seeds_slab,xmin,ymin,zmin, &
           spacing,nx,ny,k_start,slab_nz)

      allocate(length_slab(slab_nseed))
      if (do_twist) then
        allocate(twist_slab(slab_nseed))
      else
        allocate(twist_slab(0))
      endif

      if (do_qperp) then
        allocate(qperp_slab(slab_nseed))
      else
        allocate(qperp_slab(0))
      endif

      if (present(b_min)) then
        call mt_trace_fieldline_products_seedset(seeds_slab,slab_nseed, &
             dL,max_steps,length_slab,twist_slab,qperp_slab,do_twist, &
             do_qperp,b_min)
      else
        call mt_trace_fieldline_products_seedset(seeds_slab,slab_nseed, &
             dL,max_steps,length_slab,twist_slab,qperp_slab,do_twist, &
             do_qperp)
      endif

      call mt_copy_volume_length_slab(products,length_slab,nx,ny, &
           k_start,slab_nz)
      deallocate(length_slab)

      if (do_twist) then
        call mt_copy_volume_twist_slab(products,twist_slab,nx,ny, &
             k_start,slab_nz)
      endif
      deallocate(twist_slab)

      if (do_qperp) then
        call mt_copy_volume_qperp_slab(products,qperp_slab,nx,ny, &
             k_start,slab_nz)
      endif
      deallocate(qperp_slab)

      deallocate(seeds_slab)
    enddo

    call mt_write_fieldline_products_volume_vti(vti_file,origin,spacing, &
         nx,ny,nz,products, &
         do_twist,do_qperp,'mt_fieldline_products_volume_vti')
    call mt_deallocate_volume_products(products)
  end subroutine mt_fieldline_products_volume_vti

  subroutine mt_trace_fieldline_products_seedset(seeds,nseed,dL,max_steps, &
       length_results,twist_results,qperp_results,do_twist,do_qperp,b_min)
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    type(trace_length_result), intent(out) :: length_results(nseed)
    type(trace_twist_result), intent(out) :: twist_results(:)
    type(trace_qperp_result), intent(out) :: qperp_results(:)
    logical, intent(in) :: do_twist,do_qperp
    double precision, intent(in), optional :: b_min

    double precision :: seed_local(ndim)
    integer :: iseed

    if (present(b_min)) then
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iseed,seed_local) SCHEDULE(DYNAMIC,16)
      do iseed=1,nseed
        seed_local=seeds(iseed,:)
        call trace_field_length_single(seed_local,dL,max_steps, &
             length_results(iseed),b_min)
        if (do_twist) then
          call trace_field_twist_single(seed_local,dL,max_steps, &
               twist_results(iseed),b_min)
        endif
        if (do_qperp) then
          call trace_field_qperp_single(seed_local,dL,max_steps, &
               qperp_results(iseed),b_min)
        endif
      enddo
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iseed,seed_local) SCHEDULE(DYNAMIC,16)
      do iseed=1,nseed
        seed_local=seeds(iseed,:)
        call trace_field_length_single(seed_local,dL,max_steps, &
             length_results(iseed))
        if (do_twist) then
          call trace_field_twist_single(seed_local,dL,max_steps, &
               twist_results(iseed))
        endif
        if (do_qperp) then
          call trace_field_qperp_single(seed_local,dL,max_steps, &
               qperp_results(iseed))
        endif
      enddo
      !$OMP END PARALLEL DO
    endif
  end subroutine mt_trace_fieldline_products_seedset

  subroutine mt_twist_single(seed,dL,max_steps,csv_file,b_min)
    ! Trace one magnetic field line and write its length and twist summary.
    double precision, intent(in) :: seed(ndim),dL
    integer, intent(in) :: max_steps
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    type(trace_twist_result) :: result
    double precision :: seed_xyz(3)
    integer :: csv_unit,io_status

    if (npe/=1) then
      call mpistop('mt_twist_single currently requires npe=1')
    endif
    if (ndim/=3 .or. .not.slab_uniform) then
      call mpistop('mt_twist_single requires 3D uniform Cartesian geometry')
    endif

    if (present(b_min)) then
      call trace_field_twist_single(seed,dL,max_steps,result,b_min)
    else
      call trace_field_twist_single(seed,dL,max_steps,result)
    endif

    seed_xyz=0.d0
    seed_xyz(1:ndim)=result%line%seed

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop('mt_twist_single could not open CSV file')
    endif

    write(csv_unit,'(a)',iostat=io_status) &
         'seed_id,seed_x,seed_y,seed_z,'// &
         'length_total,length_backward,length_forward,'// &
         'twist_total,twist_backward,twist_forward,'// &
         'nstep_backward,nstep_forward,'// &
         'status_backward,status_forward'
    if (io_status/=0) then
      close(csv_unit)
      call mpistop('mt_twist_single could not write CSV header')
    endif

    write(csv_unit,'(i0,9(",",es24.16),4(",",i0))',iostat=io_status) &
         1,seed_xyz,result%line%total_length, &
         result%line%backward_length,result%line%forward_length, &
         result%total_twist,result%backward_twist,result%forward_twist, &
         result%line%backward_nstep,result%line%forward_nstep, &
         result%line%backward_status,result%line%forward_status
    if (io_status/=0) then
      close(csv_unit)
      call mpistop('mt_twist_single could not write CSV data')
    endif

    close(csv_unit)
  end subroutine mt_twist_single

  subroutine mt_twist_seeds(seeds,nseed,dL,max_steps,csv_file,b_min)
    ! Trace multiple magnetic field lines and write length and twist summaries.
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    type(trace_twist_result), allocatable :: results(:)
    double precision :: seed_xyz(3)
    integer :: csv_unit,io_status,iseed

    if (npe/=1) then
      call mpistop('mt_twist_seeds currently requires npe=1')
    endif
    if (ndim/=3 .or. .not.slab_uniform) then
      call mpistop('mt_twist_seeds requires 3D uniform Cartesian geometry')
    endif
    if (nseed<0) then
      call mpistop('mt_twist_seeds requires nseed>=0')
    endif

    allocate(results(nseed))
    if (present(b_min)) then
      call trace_field_twist_multi(seeds,nseed,dL,max_steps,results,b_min)
    else
      call trace_field_twist_multi(seeds,nseed,dL,max_steps,results)
    endif

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      deallocate(results)
      call mpistop('mt_twist_seeds could not open CSV file')
    endif

    write(csv_unit,'(a)',iostat=io_status) &
         'seed_id,seed_x,seed_y,seed_z,'// &
         'length_total,length_backward,length_forward,'// &
         'twist_total,twist_backward,twist_forward,'// &
         'nstep_backward,nstep_forward,'// &
         'status_backward,status_forward'
    if (io_status/=0) then
      close(csv_unit)
      deallocate(results)
      call mpistop('mt_twist_seeds could not write CSV header')
    endif

    do iseed=1,nseed
      seed_xyz=0.d0
      seed_xyz(1:ndim)=results(iseed)%line%seed
      write(csv_unit,'(i0,9(",",es24.16),4(",",i0))',iostat=io_status) &
           iseed,seed_xyz,results(iseed)%line%total_length, &
           results(iseed)%line%backward_length, &
           results(iseed)%line%forward_length, &
           results(iseed)%total_twist,results(iseed)%backward_twist, &
           results(iseed)%forward_twist, &
           results(iseed)%line%backward_nstep, &
           results(iseed)%line%forward_nstep, &
           results(iseed)%line%backward_status, &
           results(iseed)%line%forward_status
      if (io_status/=0) then
        close(csv_unit)
        deallocate(results)
        call mpistop('mt_twist_seeds could not write CSV data')
      endif
    enddo

    close(csv_unit)
    deallocate(results)
  end subroutine mt_twist_seeds

  subroutine mt_twist_plane_xy(xmin,xmax,nx,ymin,ymax,ny,z0,dL, &
       max_steps,csv_file,b_min)
    ! Trace a uniform seed grid on a constant-z plane and write twist summaries.
    integer, intent(in) :: nx,ny,max_steps
    double precision, intent(in) :: xmin,xmax,ymin,ymax,z0,dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_twist_plane_axis(xmin,xmax,nx,ymin,ymax,ny,z0,1,2,3, &
           dL,max_steps,csv_file,'mt_twist_plane_xy','ix,iy',b_min)
    else
      call mt_twist_plane_axis(xmin,xmax,nx,ymin,ymax,ny,z0,1,2,3, &
           dL,max_steps,csv_file,'mt_twist_plane_xy','ix,iy')
    endif
  end subroutine mt_twist_plane_xy

  subroutine mt_length_plane_xz(xmin,xmax,nx,zmin,zmax,nz,y0,dL, &
       max_steps,csv_file,b_min)
    integer, intent(in) :: nx,nz,max_steps
    double precision, intent(in) :: xmin,xmax,zmin,zmax,y0,dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_length_plane_axis(xmin,xmax,nx,zmin,zmax,nz,y0,1,3,2, &
           dL,max_steps,csv_file,'mt_length_plane_xz','i,j',b_min)
    else
      call mt_length_plane_axis(xmin,xmax,nx,zmin,zmax,nz,y0,1,3,2, &
           dL,max_steps,csv_file,'mt_length_plane_xz','i,j')
    endif
  end subroutine mt_length_plane_xz

  subroutine mt_length_plane_yz(ymin,ymax,ny,zmin,zmax,nz,x0,dL, &
       max_steps,csv_file,b_min)
    integer, intent(in) :: ny,nz,max_steps
    double precision, intent(in) :: ymin,ymax,zmin,zmax,x0,dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_length_plane_axis(ymin,ymax,ny,zmin,zmax,nz,x0,2,3,1, &
           dL,max_steps,csv_file,'mt_length_plane_yz','i,j',b_min)
    else
      call mt_length_plane_axis(ymin,ymax,ny,zmin,zmax,nz,x0,2,3,1, &
           dL,max_steps,csv_file,'mt_length_plane_yz','i,j')
    endif
  end subroutine mt_length_plane_yz

  subroutine mt_twist_plane_xz(xmin,xmax,nx,zmin,zmax,nz,y0,dL, &
       max_steps,csv_file,b_min)
    integer, intent(in) :: nx,nz,max_steps
    double precision, intent(in) :: xmin,xmax,zmin,zmax,y0,dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_twist_plane_axis(xmin,xmax,nx,zmin,zmax,nz,y0,1,3,2, &
           dL,max_steps,csv_file,'mt_twist_plane_xz','i,j',b_min)
    else
      call mt_twist_plane_axis(xmin,xmax,nx,zmin,zmax,nz,y0,1,3,2, &
           dL,max_steps,csv_file,'mt_twist_plane_xz','i,j')
    endif
  end subroutine mt_twist_plane_xz

  subroutine mt_twist_plane_yz(ymin,ymax,ny,zmin,zmax,nz,x0,dL, &
       max_steps,csv_file,b_min)
    integer, intent(in) :: ny,nz,max_steps
    double precision, intent(in) :: ymin,ymax,zmin,zmax,x0,dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_twist_plane_axis(ymin,ymax,ny,zmin,zmax,nz,x0,2,3,1, &
           dL,max_steps,csv_file,'mt_twist_plane_yz','i,j',b_min)
    else
      call mt_twist_plane_axis(ymin,ymax,ny,zmin,zmax,nz,x0,2,3,1, &
           dL,max_steps,csv_file,'mt_twist_plane_yz','i,j')
    endif
  end subroutine mt_twist_plane_yz

  subroutine mt_length_plane_axis(c1min,c1max,n1,c2min,c2max,n2, &
       fixed_value,axis1,axis2,fixed_axis,dL,max_steps,csv_file, &
       caller,index_header,b_min)
    integer, intent(in) :: n1,n2,axis1,axis2,fixed_axis,max_steps
    double precision, intent(in) :: c1min,c1max,c2min,c2max
    double precision, intent(in) :: fixed_value,dL
    character(len=*), intent(in) :: csv_file,caller,index_header
    double precision, intent(in), optional :: b_min

    type(trace_length_result), allocatable :: results(:)
    double precision, allocatable :: seeds(:,:)
    integer :: nseed

    call mt_validate_axis_plane(c1min,c1max,n1,c2min,c2max,n2,caller)
    call mt_build_axis_plane_seeds(c1min,c1max,n1,c2min,c2max,n2, &
         fixed_value,axis1,axis2,fixed_axis,seeds)

    nseed=n1*n2
    allocate(results(nseed))
    if (present(b_min)) then
      call trace_field_length_multi(seeds,nseed,dL,max_steps,results,b_min)
    else
      call trace_field_length_multi(seeds,nseed,dL,max_steps,results)
    endif
    call mt_write_length_plane_csv(results,n1,n2,csv_file,caller, &
         index_header)

    deallocate(seeds,results)
  end subroutine mt_length_plane_axis

  subroutine mt_twist_plane_axis(c1min,c1max,n1,c2min,c2max,n2, &
       fixed_value,axis1,axis2,fixed_axis,dL,max_steps,csv_file, &
       caller,index_header,b_min)
    integer, intent(in) :: n1,n2,axis1,axis2,fixed_axis,max_steps
    double precision, intent(in) :: c1min,c1max,c2min,c2max
    double precision, intent(in) :: fixed_value,dL
    character(len=*), intent(in) :: csv_file,caller,index_header
    double precision, intent(in), optional :: b_min

    type(trace_twist_result), allocatable :: results(:)
    double precision, allocatable :: seeds(:,:)
    integer :: nseed

    call mt_validate_axis_plane(c1min,c1max,n1,c2min,c2max,n2,caller)
    call mt_build_axis_plane_seeds(c1min,c1max,n1,c2min,c2max,n2, &
         fixed_value,axis1,axis2,fixed_axis,seeds)

    nseed=n1*n2
    allocate(results(nseed))
    if (present(b_min)) then
      call trace_field_twist_multi(seeds,nseed,dL,max_steps,results,b_min)
    else
      call trace_field_twist_multi(seeds,nseed,dL,max_steps,results)
    endif
    call mt_write_twist_plane_csv(results,n1,n2,csv_file,caller, &
         index_header)

    deallocate(seeds,results)
  end subroutine mt_twist_plane_axis

  subroutine mt_qperp_plane_axis(c1min,c1max,n1,c2min,c2max,n2, &
       fixed_value,axis1,axis2,fixed_axis,dL,max_steps,csv_file, &
       caller,index_header,b_min)
    integer, intent(in) :: n1,n2,axis1,axis2,fixed_axis,max_steps
    double precision, intent(in) :: c1min,c1max,c2min,c2max
    double precision, intent(in) :: fixed_value,dL
    character(len=*), intent(in) :: csv_file,caller,index_header
    double precision, intent(in), optional :: b_min

    type(trace_qperp_result), allocatable :: results(:)
    double precision, allocatable :: seeds(:,:)
    integer :: nseed

    call mt_validate_axis_plane(c1min,c1max,n1,c2min,c2max,n2,caller)
    call mt_build_axis_plane_seeds(c1min,c1max,n1,c2min,c2max,n2, &
         fixed_value,axis1,axis2,fixed_axis,seeds)

    nseed=n1*n2
    allocate(results(nseed))
    if (present(b_min)) then
      call trace_field_qperp_multi(seeds,nseed,dL,max_steps,results,b_min)
    else
      call trace_field_qperp_multi(seeds,nseed,dL,max_steps,results)
    endif
    call mt_write_qperp_plane_csv(results,n1,n2,csv_file,caller, &
         index_header)

    deallocate(seeds,results)
  end subroutine mt_qperp_plane_axis

  subroutine mt_topology_plane_axis(c1min,c1max,n1,c2min,c2max,n2, &
       fixed_value,axis1,axis2,fixed_axis,dL,max_steps,length_csv, &
       twist_csv,mapping_csv,caller,index_header,length_index_header, &
       b_min)
    integer, intent(in) :: n1,n2,axis1,axis2,fixed_axis,max_steps
    double precision, intent(in) :: c1min,c1max,c2min,c2max
    double precision, intent(in) :: fixed_value,dL
    character(len=*), intent(in) :: length_csv,twist_csv,mapping_csv
    character(len=*), intent(in) :: caller,index_header,length_index_header
    double precision, intent(in), optional :: b_min

    type(trace_topology_result), allocatable :: topology(:)
    type(trace_length_result), allocatable :: length_results(:)
    type(trace_twist_result), allocatable :: twist_results(:)
    type(trace_mapping_result), allocatable :: mapping_results(:)
    double precision, allocatable :: seeds(:,:)
    double precision :: source_normal(3)
    integer :: nseed
    logical :: write_twist,write_mapping,need_mapping

    call mt_validate_axis_plane(c1min,c1max,n1,c2min,c2max,n2,caller)
    if (len_trim(length_csv)==0) then
      call mpistop(trim(caller)//' requires a length CSV file')
    endif

    write_twist=len_trim(twist_csv)>0
    write_mapping=len_trim(mapping_csv)>0
    need_mapping=write_mapping

    call mt_build_axis_plane_seeds(c1min,c1max,n1,c2min,c2max,n2, &
         fixed_value,axis1,axis2,fixed_axis,seeds)
    nseed=n1*n2
    allocate(topology(nseed),length_results(nseed))
    source_normal=0.d0
    source_normal(fixed_axis)=1.d0
    if (present(b_min)) then
      call trace_field_topology_multi(seeds,nseed,dL,max_steps,topology, &
           need_twist=write_twist,need_mapping=need_mapping,b_min=b_min, &
           source_normal=source_normal)
    else
      call trace_field_topology_multi(seeds,nseed,dL,max_steps,topology, &
           need_twist=write_twist,need_mapping=need_mapping, &
           source_normal=source_normal)
    endif

    call mt_topology_to_length(topology,nseed,length_results)
    call mt_write_length_plane_csv(length_results,n1,n2,length_csv, &
         caller,length_index_header)

    if (write_twist) then
      allocate(twist_results(nseed))
      call mt_topology_to_twist(topology,nseed,twist_results)
      call mt_write_twist_plane_csv(twist_results,n1,n2,twist_csv, &
           caller,length_index_header)
      deallocate(twist_results)
    endif

    if (need_mapping) then
      allocate(mapping_results(nseed))
      call mt_topology_to_mapping(topology,nseed,mapping_results)
      if (write_mapping) then
        call mt_write_mapping_plane_csv(mapping_results,n1,n2,mapping_csv, &
             caller,index_header)
      endif
      deallocate(mapping_results)
    endif

    deallocate(length_results,topology,seeds)
  end subroutine mt_topology_plane_axis

  subroutine mt_qsl_plane_vtu_axis(c1min,c1max,n1,c2min,c2max,n2, &
       fixed_value,axis1,axis2,fixed_axis,dL,max_steps,vtu_file,caller, &
       b_min)
    integer, intent(in) :: n1,n2,axis1,axis2,fixed_axis,max_steps
    double precision, intent(in) :: c1min,c1max,c2min,c2max
    double precision, intent(in) :: fixed_value,dL
    character(len=*), intent(in) :: vtu_file,caller
    double precision, intent(in), optional :: b_min

    type(trace_length_result), allocatable :: length_results(:)
    type(trace_twist_result), allocatable :: twist_results(:)
    type(trace_mapping_result), allocatable :: mapping_results(:)
    type(trace_qperp_result), allocatable :: qperp_results(:)
    double precision, allocatable :: seeds(:,:)
    double precision :: source_normal(3)
    integer :: nseed
    logical :: need_mapping

    call mt_validate_axis_plane(c1min,c1max,n1,c2min,c2max,n2,caller)
    if (len_trim(vtu_file)==0) then
      call mpistop(trim(caller)//' requires a VTU file')
    endif

    call mt_build_axis_plane_seeds(c1min,c1max,n1,c2min,c2max,n2, &
         fixed_value,axis1,axis2,fixed_axis,seeds)
    nseed=n1*n2
    allocate(length_results(nseed),twist_results(nseed), &
         mapping_results(nseed),qperp_results(nseed))

    source_normal=0.d0
    source_normal(fixed_axis)=1.d0
    need_mapping=mt_vtk_detail_is_full()
    if (need_mapping) then
      if (present(b_min)) then
        call mt_trace_qsl_plane_full(seeds,nseed,dL,max_steps,source_normal, &
             qperp_results,twist_results,mapping_results,b_min)
      else
        call mt_trace_qsl_plane_full(seeds,nseed,dL,max_steps,source_normal, &
             qperp_results,twist_results,mapping_results)
      endif
      call mt_qperp_to_length(qperp_results,nseed,length_results)
    else
      if (present(b_min)) then
        call mt_trace_qsl_plane_minimal(seeds,nseed,dL,max_steps, &
             qperp_results,twist_results,b_min)
      else
        call mt_trace_qsl_plane_minimal(seeds,nseed,dL,max_steps, &
             qperp_results,twist_results)
      endif
      call mt_qperp_to_length(qperp_results,nseed,length_results)
    endif

    call mt_write_qsl_plane_vtu(vtu_file,length_results,twist_results, &
         mapping_results,qperp_results,n1,n2,caller)

    deallocate(qperp_results,mapping_results,twist_results, &
         length_results,seeds)
  end subroutine mt_qsl_plane_vtu_axis

  subroutine mt_trace_qsl_plane_full(seeds,nseed,dL,max_steps,source_normal, &
       qperp_results,twist_results,mapping_results,b_min)
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL,source_normal(3)
    type(trace_qperp_result), intent(out) :: qperp_results(nseed)
    type(trace_twist_result), intent(out) :: twist_results(nseed)
    type(trace_mapping_result), intent(out) :: mapping_results(nseed)
    double precision, intent(in), optional :: b_min

    double precision :: seed_local(ndim)
    integer :: iseed

    if (present(b_min)) then
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iseed,seed_local) SCHEDULE(DYNAMIC,16)
      do iseed=1,nseed
        seed_local=seeds(iseed,:)
        call trace_field_qperp_single(seed_local,dL,max_steps, &
             qperp_results(iseed),b_min)
        call trace_field_twist_single(seed_local,dL,max_steps, &
             twist_results(iseed),b_min)
        call trace_field_mapping_single(seed_local,dL,max_steps, &
             mapping_results(iseed),b_min,source_normal)
      enddo
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iseed,seed_local) SCHEDULE(DYNAMIC,16)
      do iseed=1,nseed
        seed_local=seeds(iseed,:)
        call trace_field_qperp_single(seed_local,dL,max_steps, &
             qperp_results(iseed))
        call trace_field_twist_single(seed_local,dL,max_steps, &
             twist_results(iseed))
        call trace_field_mapping_single(seed_local,dL,max_steps, &
             mapping_results(iseed),source_normal=source_normal)
      enddo
      !$OMP END PARALLEL DO
    endif
  end subroutine mt_trace_qsl_plane_full

  subroutine mt_trace_qsl_plane_minimal(seeds,nseed,dL,max_steps, &
       qperp_results,twist_results,b_min)
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    type(trace_qperp_result), intent(out) :: qperp_results(nseed)
    type(trace_twist_result), intent(out) :: twist_results(nseed)
    double precision, intent(in), optional :: b_min

    double precision :: seed_local(ndim)
    integer :: iseed

    if (present(b_min)) then
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iseed,seed_local) SCHEDULE(DYNAMIC,16)
      do iseed=1,nseed
        seed_local=seeds(iseed,:)
        call trace_field_qperp_single(seed_local,dL,max_steps, &
             qperp_results(iseed),b_min)
        call trace_field_twist_single(seed_local,dL,max_steps, &
             twist_results(iseed),b_min)
      enddo
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iseed,seed_local) SCHEDULE(DYNAMIC,16)
      do iseed=1,nseed
        seed_local=seeds(iseed,:)
        call trace_field_qperp_single(seed_local,dL,max_steps, &
             qperp_results(iseed))
        call trace_field_twist_single(seed_local,dL,max_steps, &
             twist_results(iseed))
      enddo
      !$OMP END PARALLEL DO
    endif
  end subroutine mt_trace_qsl_plane_minimal

  subroutine mt_qperp_to_length(qperp_results,nseed,results)
    integer, intent(in) :: nseed
    type(trace_qperp_result), intent(in) :: qperp_results(nseed)
    type(trace_length_result), intent(out) :: results(nseed)

    integer :: iseed

    do iseed=1,nseed
      results(iseed)%seed=qperp_results(iseed)%seed
      results(iseed)%forward_footpoint=qperp_results(iseed)%forward_endpoint
      results(iseed)%backward_footpoint=qperp_results(iseed)%backward_endpoint
      results(iseed)%forward_status=qperp_results(iseed)%forward_status
      results(iseed)%backward_status=qperp_results(iseed)%backward_status
      results(iseed)%forward_nstep=0
      results(iseed)%backward_nstep=0
      if (qperp_results(iseed)%valid_q0) then
        results(iseed)%forward_length=qperp_results(iseed)%forward_length
        results(iseed)%backward_length=qperp_results(iseed)%backward_length
        results(iseed)%total_length=results(iseed)%forward_length &
             +results(iseed)%backward_length
      else
        results(iseed)%forward_length=0.d0
        results(iseed)%backward_length=0.d0
        results(iseed)%total_length=0.d0
      endif
    enddo
  end subroutine mt_qperp_to_length

  subroutine mt_topology_to_length(topology,nseed,results)
    integer, intent(in) :: nseed
    type(trace_topology_result), intent(in) :: topology(nseed)
    type(trace_length_result), intent(out) :: results(nseed)

    integer :: iseed

    do iseed=1,nseed
      results(iseed)%seed=topology(iseed)%seed
      results(iseed)%forward_footpoint=topology(iseed)%forward_endpoint
      results(iseed)%backward_footpoint=topology(iseed)%backward_endpoint
      results(iseed)%forward_length=topology(iseed)%length_forward
      results(iseed)%backward_length=topology(iseed)%length_backward
      results(iseed)%total_length=topology(iseed)%length_total
      results(iseed)%forward_nstep=topology(iseed)%forward_nstep
      results(iseed)%backward_nstep=topology(iseed)%backward_nstep
      results(iseed)%forward_status=topology(iseed)%forward_status
      results(iseed)%backward_status=topology(iseed)%backward_status
    enddo
  end subroutine mt_topology_to_length

  subroutine mt_topology_to_twist(topology,nseed,results)
    integer, intent(in) :: nseed
    type(trace_topology_result), intent(in) :: topology(nseed)
    type(trace_twist_result), intent(out) :: results(nseed)

    integer :: iseed

    do iseed=1,nseed
      results(iseed)%line%seed=topology(iseed)%seed
      results(iseed)%line%forward_footpoint=topology(iseed)%forward_endpoint
      results(iseed)%line%backward_footpoint=topology(iseed)%backward_endpoint
      results(iseed)%line%forward_length=topology(iseed)%length_forward
      results(iseed)%line%backward_length=topology(iseed)%length_backward
      results(iseed)%line%total_length=topology(iseed)%length_total
      results(iseed)%line%forward_nstep=topology(iseed)%forward_nstep
      results(iseed)%line%backward_nstep=topology(iseed)%backward_nstep
      results(iseed)%line%forward_status=topology(iseed)%forward_status
      results(iseed)%line%backward_status=topology(iseed)%backward_status
      results(iseed)%forward_twist=topology(iseed)%twist_forward
      results(iseed)%backward_twist=topology(iseed)%twist_backward
      results(iseed)%total_twist=topology(iseed)%twist_total
    enddo
  end subroutine mt_topology_to_twist

  subroutine mt_topology_to_mapping(topology,nseed,results)
    integer, intent(in) :: nseed
    type(trace_topology_result), intent(in) :: topology(nseed)
    type(trace_mapping_result), intent(out) :: results(nseed)

    integer :: iseed

    do iseed=1,nseed
      results(iseed)%seed=topology(iseed)%seed
      results(iseed)%source_B=topology(iseed)%source_B
      results(iseed)%forward_footpoint=topology(iseed)%map_forward_endpoint
      results(iseed)%backward_footpoint=topology(iseed)%map_backward_endpoint
      results(iseed)%forward_B=topology(iseed)%forward_B
      results(iseed)%backward_B=topology(iseed)%backward_B
      results(iseed)%forward_length=topology(iseed)%map_forward_length
      results(iseed)%backward_length=topology(iseed)%map_backward_length
      results(iseed)%source_Bn=topology(iseed)%source_Bn
      results(iseed)%forward_Bn=topology(iseed)%forward_Bn
      results(iseed)%backward_Bn=topology(iseed)%backward_Bn
      results(iseed)%forward_face=topology(iseed)%map_forward_face
      results(iseed)%backward_face=topology(iseed)%map_backward_face
      results(iseed)%forward_status=topology(iseed)%map_forward_status
      results(iseed)%backward_status=topology(iseed)%map_backward_status
      results(iseed)%valid=topology(iseed)%has_mapping .and. &
           topology(iseed)%valid
    enddo
  end subroutine mt_topology_to_mapping

  subroutine mt_validate_axis_plane(c1min,c1max,n1,c2min,c2max,n2,caller)
    integer, intent(in) :: n1,n2
    double precision, intent(in) :: c1min,c1max,c2min,c2max
    character(len=*), intent(in) :: caller

    if (npe/=1) then
      call mpistop(trim(caller)//' currently requires npe=1')
    endif
    if (ndim/=3) then
      call mpistop(trim(caller)//' requires ndim=3')
      return
    endif
    if (.not.slab_uniform) then
      call mpistop(trim(caller)//' requires uniform Cartesian geometry')
    endif
    if (n1<1 .or. n2<1) then
      call mpistop(trim(caller)//' requires both sample counts >=1')
    endif
    if (c1max<c1min .or. c2max<c2min) then
      call mpistop(trim(caller)//' requires ordered plane bounds')
    endif
  end subroutine mt_validate_axis_plane

  subroutine mt_build_axis_plane_seeds(c1min,c1max,n1,c2min,c2max,n2, &
       fixed_value,axis1,axis2,fixed_axis,seeds)
    integer, intent(in) :: n1,n2,axis1,axis2,fixed_axis
    double precision, intent(in) :: c1min,c1max,c2min,c2max,fixed_value
    double precision, allocatable, intent(out) :: seeds(:,:)

    double precision :: dc1,dc2
    integer :: i,j,iseed

    dc1=0.d0
    dc2=0.d0
    if (n1>1) dc1=(c1max-c1min)/dble(n1-1)
    if (n2>1) dc2=(c2max-c2min)/dble(n2-1)

    allocate(seeds(n1*n2,ndim))
    seeds=0.d0
    do j=1,n2
      do i=1,n1
        iseed=(j-1)*n1+i
        seeds(iseed,axis1)=c1min+dble(i-1)*dc1
        seeds(iseed,axis2)=c2min+dble(j-1)*dc2
        seeds(iseed,fixed_axis)=fixed_value
      enddo
    enddo
  end subroutine mt_build_axis_plane_seeds

  logical function mt_validate_arbitrary_plane_basis(e1,e2,s1min,s1max, &
       n1,s2min,s2max,n2,caller) result(valid)
    integer, intent(in) :: n1,n2
    double precision, intent(in) :: e1(ndim),e2(ndim)
    double precision, intent(in) :: s1min,s1max,s2min,s2max
    character(len=*), intent(in) :: caller

    double precision, parameter :: basis_tol=1.d-10
    double precision :: e1_norm,e2_norm,e12_dot

    valid=.false.
    if (npe/=1) then
      write(*,'(a)') trim(caller)//' currently requires npe=1'
      return
    endif
    if (ndim/=3) then
      write(*,'(a)') trim(caller)//' requires ndim=3'
      return
    endif
    if (.not.slab_uniform) then
      write(*,'(a)') trim(caller)//' requires uniform Cartesian geometry'
      return
    endif
    if (n1<1 .or. n2<1) then
      write(*,'(a)') trim(caller)//' requires both sample counts >=1'
      return
    endif
    if (s1max<s1min .or. s2max<s2min) then
      write(*,'(a)') trim(caller)//' requires ordered plane bounds'
      return
    endif
    e1_norm=dsqrt(sum(e1**2))
    e2_norm=dsqrt(sum(e2**2))
    e12_dot=sum(e1*e2)
    if (abs(e1_norm-1.d0)>basis_tol .or. &
         abs(e2_norm-1.d0)>basis_tol .or. &
         abs(e12_dot)>basis_tol) then
      write(*,'(a)') trim(caller)//' requires orthonormal e1/e2'
      return
    endif

    valid=.true.
  end function mt_validate_arbitrary_plane_basis

  subroutine mt_build_arbitrary_plane_seeds(origin,e1,e2,s1min,s1max,n1, &
       s2min,s2max,n2,seeds,s1,s2)
    integer, intent(in) :: n1,n2
    double precision, intent(in) :: origin(ndim),e1(ndim),e2(ndim)
    double precision, intent(in) :: s1min,s1max,s2min,s2max
    double precision, allocatable, intent(out) :: seeds(:,:),s1(:),s2(:)

    double precision :: ds1,ds2
    integer :: i,j,iseed

    ds1=0.d0
    ds2=0.d0
    if (n1>1) ds1=(s1max-s1min)/dble(n1-1)
    if (n2>1) ds2=(s2max-s2min)/dble(n2-1)

    allocate(seeds(n1*n2,ndim),s1(n1*n2),s2(n1*n2))
    do j=1,n2
      do i=1,n1
        iseed=(j-1)*n1+i
        s1(iseed)=s1min+dble(i-1)*ds1
        s2(iseed)=s2min+dble(j-1)*ds2
        seeds(iseed,:)=origin+s1(iseed)*e1+s2(iseed)*e2
      enddo
    enddo
  end subroutine mt_build_arbitrary_plane_seeds

  subroutine mt_write_length_plane_csv(results,n1,n2,csv_file,caller, &
       index_header)
    integer, intent(in) :: n1,n2
    type(trace_length_result), intent(in) :: results(n1*n2)
    character(len=*), intent(in) :: csv_file,caller,index_header

    double precision :: seed_xyz(3)
    integer :: csv_unit,io_status,i,j,iseed

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open CSV file')
    endif

    write(csv_unit,'(a)',iostat=io_status) &
         trim(index_header)//',seed_x,seed_y,seed_z,'// &
         'length_total,length_backward,length_forward,'// &
         'nstep_backward,nstep_forward,'// &
         'status_backward,status_forward'
    if (io_status/=0) then
      close(csv_unit)
      call mpistop(trim(caller)//' could not write CSV header')
    endif

    do j=1,n2
      do i=1,n1
        iseed=(j-1)*n1+i
        seed_xyz=0.d0
        seed_xyz(1:ndim)=results(iseed)%seed
        write(csv_unit,'(i0,",",i0,6(",",es24.16),4(",",i0))', &
             iostat=io_status) i,j,seed_xyz,results(iseed)%total_length, &
             results(iseed)%backward_length, &
             results(iseed)%forward_length, &
             results(iseed)%backward_nstep, &
             results(iseed)%forward_nstep, &
             results(iseed)%backward_status, &
             results(iseed)%forward_status
        if (io_status/=0) then
          close(csv_unit)
          call mpistop(trim(caller)//' could not write CSV data')
        endif
      enddo
    enddo

    close(csv_unit)
  end subroutine mt_write_length_plane_csv

  subroutine mt_write_fieldline_products_seeds_csv(length_results, &
       twist_results,qperp_results,nseed,csv_file,do_twist,do_qperp)
    integer, intent(in) :: nseed
    type(trace_length_result), intent(in) :: length_results(nseed)
    type(trace_twist_result), intent(in) :: twist_results(:)
    type(trace_qperp_result), intent(in) :: qperp_results(:)
    character(len=*), intent(in) :: csv_file
    logical, intent(in) :: do_twist,do_qperp

    double precision :: seed_xyz(3)
    double precision :: Bseed_norm,Bf_norm,Bb_norm
    character(len=2048) :: header
    integer :: csv_unit,io_status,iseed

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop('mt_fieldline_products_seeds could not open CSV file')
    endif

    header='seed_id,seed_x,seed_y,seed_z,'// &
         'length_total,length_backward,length_forward,'// &
         'nstep_backward_length,nstep_forward_length,'// &
         'status_backward_length,status_forward_length'
    if (do_twist) then
      header=trim(header)//','// &
           'twist_total,twist_backward,twist_forward,'// &
           'nstep_backward_twist,nstep_forward_twist,'// &
           'status_backward_twist,status_forward_twist'
    endif
    if (do_qperp) then
      header=trim(header)//','// &
           'qperp,logqperp,valid_qperp,status_qperp,'// &
           'N2,bfactor,'// &
           'face_forward_qperp,face_backward_qperp,'// &
           'status_forward_qperp,status_backward_qperp,'// &
           'length_forward_qperp,length_backward_qperp,'// &
           'Bseed_norm,Bf_norm,Bb_norm,'// &
           'x_f_qperp,y_f_qperp,z_f_qperp,'// &
           'x_b_qperp,y_b_qperp,z_b_qperp'
    endif
    write(csv_unit,'(a)',iostat=io_status) trim(header)
    if (io_status/=0) then
      close(csv_unit)
      call mpistop('mt_fieldline_products_seeds could not write CSV header')
    endif

    do iseed=1,nseed
      seed_xyz=0.d0
      seed_xyz(1:ndim)=length_results(iseed)%seed
      write(csv_unit,'(i0,6(",",es24.16),4(",",i0))', &
           advance='no',iostat=io_status) &
           iseed,seed_xyz,length_results(iseed)%total_length, &
           length_results(iseed)%backward_length, &
           length_results(iseed)%forward_length, &
           length_results(iseed)%backward_nstep, &
           length_results(iseed)%forward_nstep, &
           length_results(iseed)%backward_status, &
           length_results(iseed)%forward_status
      if (io_status/=0) exit

      if (do_twist) then
        write(csv_unit,'(3(",",es24.16),4(",",i0))', &
             advance='no',iostat=io_status) &
             twist_results(iseed)%total_twist, &
             twist_results(iseed)%backward_twist, &
             twist_results(iseed)%forward_twist, &
             twist_results(iseed)%line%backward_nstep, &
             twist_results(iseed)%line%forward_nstep, &
             twist_results(iseed)%line%backward_status, &
             twist_results(iseed)%line%forward_status
        if (io_status/=0) exit
      endif

      if (do_qperp) then
        Bseed_norm=dsqrt(sum(qperp_results(iseed)%B_seed**2))
        Bf_norm=dsqrt(sum(qperp_results(iseed)%forward_B**2))
        Bb_norm=dsqrt(sum(qperp_results(iseed)%backward_B**2))
        write(csv_unit, &
             '(2(",",es24.16),",",l1,",",i0,2(",",es24.16),'// &
             '4(",",i0),11(",",es24.16))', &
             advance='no',iostat=io_status) &
             qperp_results(iseed)%qperp, &
             qperp_results(iseed)%logqperp, &
             qperp_results(iseed)%valid,qperp_results(iseed)%status, &
             qperp_results(iseed)%N2,qperp_results(iseed)%bfactor, &
             qperp_results(iseed)%forward_face, &
             qperp_results(iseed)%backward_face, &
             qperp_results(iseed)%forward_status, &
             qperp_results(iseed)%backward_status, &
             qperp_results(iseed)%forward_length, &
             qperp_results(iseed)%backward_length, &
             Bseed_norm,Bf_norm,Bb_norm, &
             qperp_results(iseed)%forward_endpoint, &
             qperp_results(iseed)%backward_endpoint
        if (io_status/=0) exit
      endif

      write(csv_unit,'()',iostat=io_status)
      if (io_status/=0) exit
    enddo

    if (io_status/=0) then
      close(csv_unit)
      call mpistop('mt_fieldline_products_seeds could not write CSV data')
    endif
    close(csv_unit)
  end subroutine mt_write_fieldline_products_seeds_csv

  subroutine mt_fieldline_products_append_header(header,do_twist,do_qperp)
    character(len=*), intent(inout) :: header
    logical, intent(in) :: do_twist,do_qperp

    header=trim(header)//','// &
         'length_total,length_backward,length_forward,'// &
         'nstep_backward_length,nstep_forward_length,'// &
         'status_backward_length,status_forward_length'
    if (do_twist) then
      header=trim(header)//','// &
           'twist_total,twist_backward,twist_forward,'// &
           'nstep_backward_twist,nstep_forward_twist,'// &
           'status_backward_twist,status_forward_twist'
    endif
    if (do_qperp) then
      header=trim(header)//','// &
           'qperp,logqperp,valid_qperp,status_qperp,'// &
           'N2,bfactor,'// &
           'face_forward_qperp,face_backward_qperp,'// &
           'status_forward_qperp,status_backward_qperp,'// &
           'length_forward_qperp,length_backward_qperp,'// &
           'Bseed_norm,Bf_norm,Bb_norm,'// &
           'x_f_qperp,y_f_qperp,z_f_qperp,'// &
           'x_b_qperp,y_b_qperp,z_b_qperp'
    endif
  end subroutine mt_fieldline_products_append_header

  subroutine mt_write_fieldline_products_plane_arbitrary_header(csv_file, &
       caller,do_twist,do_qperp)
    character(len=*), intent(in) :: csv_file,caller
    logical, intent(in) :: do_twist,do_qperp

    character(len=2048) :: header
    integer :: csv_unit,io_status

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open CSV file')
    endif

    header='i,j,s1,s2,seed_x,seed_y,seed_z'
    call mt_fieldline_products_append_header(header,do_twist,do_qperp)
    write(csv_unit,'(a)',iostat=io_status) trim(header)
    if (io_status/=0) then
      close(csv_unit)
      call mpistop(trim(caller)//' could not write CSV header')
    endif

    close(csv_unit)
  end subroutine mt_write_fieldline_products_plane_arbitrary_header

  subroutine mt_write_fieldline_products_plane_arbitrary_csv( &
       length_results,twist_results,qperp_results,s1,s2,n1,n2,csv_file, &
       do_twist,do_qperp,caller)
    integer, intent(in) :: n1,n2
    type(trace_length_result), intent(in) :: length_results(n1*n2)
    type(trace_twist_result), intent(in) :: twist_results(:)
    type(trace_qperp_result), intent(in) :: qperp_results(:)
    double precision, intent(in) :: s1(n1*n2),s2(n1*n2)
    character(len=*), intent(in) :: csv_file,caller
    logical, intent(in) :: do_twist,do_qperp

    double precision :: seed_xyz(3)
    double precision :: Bseed_norm,Bf_norm,Bb_norm
    character(len=2048) :: header
    integer :: csv_unit,io_status,i,j,iseed

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open CSV file')
    endif

    header='i,j,s1,s2,seed_x,seed_y,seed_z'
    call mt_fieldline_products_append_header(header,do_twist,do_qperp)
    write(csv_unit,'(a)',iostat=io_status) trim(header)
    if (io_status/=0) then
      close(csv_unit)
      call mpistop(trim(caller)//' could not write CSV header')
    endif

    do j=1,n2
      do i=1,n1
        iseed=(j-1)*n1+i
        seed_xyz=0.d0
        seed_xyz(1:ndim)=length_results(iseed)%seed
        write(csv_unit,'(i0,",",i0,8(",",es24.16),4(",",i0))', &
             advance='no',iostat=io_status) &
             i,j,s1(iseed),s2(iseed),seed_xyz, &
             length_results(iseed)%total_length, &
             length_results(iseed)%backward_length, &
             length_results(iseed)%forward_length, &
             length_results(iseed)%backward_nstep, &
             length_results(iseed)%forward_nstep, &
             length_results(iseed)%backward_status, &
             length_results(iseed)%forward_status
        if (io_status/=0) exit

        if (do_twist) then
          write(csv_unit,'(3(",",es24.16),4(",",i0))', &
               advance='no',iostat=io_status) &
               twist_results(iseed)%total_twist, &
               twist_results(iseed)%backward_twist, &
               twist_results(iseed)%forward_twist, &
               twist_results(iseed)%line%backward_nstep, &
               twist_results(iseed)%line%forward_nstep, &
               twist_results(iseed)%line%backward_status, &
               twist_results(iseed)%line%forward_status
          if (io_status/=0) exit
        endif

        if (do_qperp) then
          Bseed_norm=dsqrt(sum(qperp_results(iseed)%B_seed**2))
          Bf_norm=dsqrt(sum(qperp_results(iseed)%forward_B**2))
          Bb_norm=dsqrt(sum(qperp_results(iseed)%backward_B**2))
          write(csv_unit, &
               '(2(",",es24.16),",",l1,",",i0,2(",",es24.16),'// &
               '4(",",i0),11(",",es24.16))', &
               advance='no',iostat=io_status) &
               qperp_results(iseed)%qperp, &
               qperp_results(iseed)%logqperp, &
               qperp_results(iseed)%valid,qperp_results(iseed)%status, &
               qperp_results(iseed)%N2,qperp_results(iseed)%bfactor, &
               qperp_results(iseed)%forward_face, &
               qperp_results(iseed)%backward_face, &
               qperp_results(iseed)%forward_status, &
               qperp_results(iseed)%backward_status, &
               qperp_results(iseed)%forward_length, &
               qperp_results(iseed)%backward_length, &
               Bseed_norm,Bf_norm,Bb_norm, &
               qperp_results(iseed)%forward_endpoint, &
               qperp_results(iseed)%backward_endpoint
          if (io_status/=0) exit
        endif

        write(csv_unit,'()',iostat=io_status)
        if (io_status/=0) exit
      enddo
      if (io_status/=0) exit
    enddo

    if (io_status/=0) then
      close(csv_unit)
      call mpistop(trim(caller)//' could not write CSV data')
    endif
    close(csv_unit)
  end subroutine mt_write_fieldline_products_plane_arbitrary_csv

  subroutine mt_write_fieldline_products_vtu_empty(vtu_file,do_twist, &
       do_qperp,caller)
    character(len=*), intent(in) :: vtu_file,caller
    logical, intent(in) :: do_twist,do_qperp

    integer :: vtu_unit,io_status

    open(newunit=vtu_unit,file=trim(vtu_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open VTU file')
    endif

    call mt_write_vtu_file_header(vtu_unit,0,0)
    call mt_write_vtu_empty_pointdata(vtu_unit,do_twist,do_qperp)
    write(vtu_unit,'(a)',iostat=io_status) '<CellData>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '</CellData>'
    if (io_status==0) then
      write(vtu_unit,'(a)',iostat=io_status) '<Points>'
    endif
    if (io_status==0) then
      write(vtu_unit,'(a)',iostat=io_status) &
           '<DataArray type="Float64" NumberOfComponents="3" format="ascii">'
    endif
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '</DataArray>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '</Points>'
    if (io_status==0) call mt_write_vtu_cells_empty(vtu_unit,io_status)
    if (io_status==0) call mt_write_vtu_file_footer(vtu_unit,io_status)
    if (io_status/=0) then
      close(vtu_unit)
      call mpistop(trim(caller)//' could not write VTU file')
    endif

    close(vtu_unit)
  end subroutine mt_write_fieldline_products_vtu_empty

  subroutine mt_write_fieldline_products_vtu_vertices(vtu_file, &
       length_results,twist_results,qperp_results,npoint,do_twist, &
       do_qperp,caller)
    integer, intent(in) :: npoint
    type(trace_length_result), intent(in) :: length_results(npoint)
    type(trace_twist_result), intent(in) :: twist_results(:)
    type(trace_qperp_result), intent(in) :: qperp_results(:)
    character(len=*), intent(in) :: vtu_file,caller
    logical, intent(in) :: do_twist,do_qperp

    integer :: vtu_unit,io_status

    open(newunit=vtu_unit,file=trim(vtu_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open VTU file')
    endif

    call mt_write_vtu_file_header(vtu_unit,npoint,npoint)
    call mt_write_vtu_product_pointdata(vtu_unit,length_results, &
         twist_results,qperp_results,npoint,do_twist,do_qperp,io_status)
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '<CellData>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '</CellData>'
    if (io_status==0) call mt_write_vtu_points(vtu_unit,length_results, &
         npoint,io_status)
    if (io_status==0) call mt_write_vtu_vertex_cells(vtu_unit,npoint, &
         io_status)
    if (io_status==0) call mt_write_vtu_file_footer(vtu_unit,io_status)
    if (io_status/=0) then
      close(vtu_unit)
      call mpistop(trim(caller)//' could not write VTU file')
    endif

    close(vtu_unit)
  end subroutine mt_write_fieldline_products_vtu_vertices

  subroutine mt_write_fieldline_products_vtu_plane(vtu_file, &
       length_results,twist_results,qperp_results,n1,n2,do_twist, &
       do_qperp,caller)
    integer, intent(in) :: n1,n2
    type(trace_length_result), intent(in) :: length_results(n1*n2)
    type(trace_twist_result), intent(in) :: twist_results(:)
    type(trace_qperp_result), intent(in) :: qperp_results(:)
    character(len=*), intent(in) :: vtu_file,caller
    logical, intent(in) :: do_twist,do_qperp

    integer :: npoint,ncell,vtu_unit,io_status

    npoint=n1*n2
    if (n1>1 .and. n2>1) then
      ncell=(n1-1)*(n2-1)
    else
      ncell=npoint
    endif

    open(newunit=vtu_unit,file=trim(vtu_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open VTU file')
    endif

    call mt_write_vtu_file_header(vtu_unit,npoint,ncell)
    call mt_write_vtu_product_pointdata(vtu_unit,length_results, &
         twist_results,qperp_results,npoint,do_twist,do_qperp,io_status)
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '<CellData>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '</CellData>'
    if (io_status==0) call mt_write_vtu_points(vtu_unit,length_results, &
         npoint,io_status)
    if (io_status==0) then
      if (n1>1 .and. n2>1) then
        call mt_write_vtu_quad_cells(vtu_unit,n1,n2,io_status)
      else
        call mt_write_vtu_vertex_cells(vtu_unit,npoint,io_status)
      endif
    endif
    if (io_status==0) call mt_write_vtu_file_footer(vtu_unit,io_status)
    if (io_status/=0) then
      close(vtu_unit)
      call mpistop(trim(caller)//' could not write VTU file')
    endif

    close(vtu_unit)
  end subroutine mt_write_fieldline_products_vtu_plane

  subroutine mt_write_qsl_plane_vtu(vtu_file,length_results,twist_results, &
       mapping_results,qperp_results,n1,n2,caller)
    integer, intent(in) :: n1,n2
    type(trace_length_result), intent(in) :: length_results(n1*n2)
    type(trace_twist_result), intent(in) :: twist_results(n1*n2)
    type(trace_mapping_result), intent(in) :: mapping_results(n1*n2)
    type(trace_qperp_result), intent(in) :: qperp_results(n1*n2)
    character(len=*), intent(in) :: vtu_file,caller

    integer :: npoint,ncell,vtu_unit,io_status

    npoint=n1*n2
    if (n1>1 .and. n2>1) then
      ncell=(n1-1)*(n2-1)
    else
      ncell=npoint
    endif

    open(newunit=vtu_unit,file=trim(vtu_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open VTU file')
    endif

    call mt_write_vtu_file_header(vtu_unit,npoint,ncell)
    call mt_write_qsl_plane_vtu_pointdata(vtu_unit,length_results, &
         twist_results,mapping_results,qperp_results,npoint,io_status)
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '<CellData>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '</CellData>'
    if (io_status==0) call mt_write_vtu_points(vtu_unit,length_results, &
         npoint,io_status)
    if (io_status==0) then
      if (n1>1 .and. n2>1) then
        call mt_write_vtu_quad_cells(vtu_unit,n1,n2,io_status)
      else
        call mt_write_vtu_vertex_cells(vtu_unit,npoint,io_status)
      endif
    endif
    if (io_status==0) call mt_write_vtu_file_footer(vtu_unit,io_status)
    if (io_status/=0) then
      close(vtu_unit)
      call mpistop(trim(caller)//' could not write VTU file')
    endif

    close(vtu_unit)
  end subroutine mt_write_qsl_plane_vtu

  subroutine mt_write_qsl_plane_vtu_pointdata(vtu_unit,length_results, &
       twist_results,mapping_results,qperp_results,npoint,io_status)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_length_result), intent(in) :: length_results(npoint)
    type(trace_twist_result), intent(in) :: twist_results(npoint)
    type(trace_mapping_result), intent(in) :: mapping_results(npoint)
    type(trace_qperp_result), intent(in) :: qperp_results(npoint)
    integer, intent(inout) :: io_status

    write(vtu_unit,'(a)',iostat=io_status) '<PointData>'
    if (io_status/=0) return

    if (mt_vtk_detail_is_full()) then
      call mt_write_vtu_length_pointdata(vtu_unit,length_results,npoint, &
           io_status)
      call mt_write_vtu_twist_pointdata(vtu_unit,twist_results,npoint, &
           io_status)
      call mt_write_vtu_mapping_pointdata(vtu_unit,mapping_results,npoint, &
           io_status)
      call mt_write_vtu_q_product_pointdata(vtu_unit,qperp_results, &
           npoint,.true.,io_status)
      call mt_write_vtu_qperp_public_pointdata(vtu_unit,qperp_results, &
           npoint,io_status)
      call mt_write_vtu_qperp_method2_pointdata(vtu_unit,qperp_results, &
           npoint,io_status)
    else
      call mt_write_vtu_qsl_minimal_pointdata(vtu_unit,length_results, &
           twist_results,qperp_results,npoint,io_status)
    endif

    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '</PointData>'
  end subroutine mt_write_qsl_plane_vtu_pointdata

  subroutine mt_write_vtu_qsl_minimal_pointdata(vtu_unit,length_results, &
       twist_results,qperp_results,npoint,io_status)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_length_result), intent(in) :: length_results(npoint)
    type(trace_twist_result), intent(in) :: twist_results(npoint)
    type(trace_qperp_result), intent(in) :: qperp_results(npoint)
    integer, intent(inout) :: io_status

    integer :: ipoint

    call mt_write_vtu_float_array_start(vtu_unit,'length_total',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_visual_float(length_results(ipoint)%total_length), &
         ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    call mt_write_vtu_float_array_start(vtu_unit,'twist_total',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_visual_float(twist_results(ipoint)%total_twist), &
         ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    call mt_write_vtu_float_array_start(vtu_unit,'logQ',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_visual_valid_float(qperp_results(ipoint)%logq0, &
         qperp_results(ipoint)%valid_q0),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    call mt_write_vtu_float_array_start(vtu_unit,'logQperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_visual_valid_float(qperp_results(ipoint)%logqperp, &
         qperp_results(ipoint)%valid),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    call mt_write_vtu_int_array_start(vtu_unit,'connection_type_Q', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (mt_q_result_connection_type(qperp_results(ipoint)),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

  end subroutine mt_write_vtu_qsl_minimal_pointdata

  double precision function mt_vti_axis_spacing(xmin,xmax,n)
    integer, intent(in) :: n
    double precision, intent(in) :: xmin,xmax

    if (n>1) then
      mt_vti_axis_spacing=(xmax-xmin)/dble(n-1)
    else
      mt_vti_axis_spacing=0.d0
    endif
  end function mt_vti_axis_spacing

  subroutine mt_allocate_volume_products(products,nseed,do_twist,do_qperp)
    type(mt_volume_products), intent(inout) :: products
    integer, intent(in) :: nseed
    logical, intent(in) :: do_twist,do_qperp

    allocate(products%length_total(nseed))
    allocate(products%length_backward(nseed))
    allocate(products%length_forward(nseed))
    allocate(products%nstep_backward_length(nseed))
    allocate(products%nstep_forward_length(nseed))
    allocate(products%status_backward_length(nseed))
    allocate(products%status_forward_length(nseed))

    if (do_twist) then
      allocate(products%twist_total(nseed))
      allocate(products%twist_backward(nseed))
      allocate(products%twist_forward(nseed))
      allocate(products%nstep_backward_twist(nseed))
      allocate(products%nstep_forward_twist(nseed))
      allocate(products%status_backward_twist(nseed))
      allocate(products%status_forward_twist(nseed))
    endif

    if (do_qperp) then
      allocate(products%qperp(nseed))
      allocate(products%logqperp(nseed))
      allocate(products%N2(nseed))
      allocate(products%bfactor(nseed))
      allocate(products%length_forward_qperp(nseed))
      allocate(products%length_backward_qperp(nseed))
      allocate(products%Bseed_norm(nseed))
      allocate(products%Bf_norm(nseed))
      allocate(products%Bb_norm(nseed))
      allocate(products%valid_qperp(nseed))
      allocate(products%status_qperp(nseed))
      allocate(products%face_forward_qperp(nseed))
      allocate(products%face_backward_qperp(nseed))
      allocate(products%status_forward_qperp(nseed))
      allocate(products%status_backward_qperp(nseed))
    endif
  end subroutine mt_allocate_volume_products

  subroutine mt_deallocate_volume_products(products)
    type(mt_volume_products), intent(inout) :: products

    if (allocated(products%length_total)) deallocate(products%length_total)
    if (allocated(products%length_backward)) deallocate(products%length_backward)
    if (allocated(products%length_forward)) deallocate(products%length_forward)
    if (allocated(products%nstep_backward_length)) &
         deallocate(products%nstep_backward_length)
    if (allocated(products%nstep_forward_length)) &
         deallocate(products%nstep_forward_length)
    if (allocated(products%status_backward_length)) &
         deallocate(products%status_backward_length)
    if (allocated(products%status_forward_length)) &
         deallocate(products%status_forward_length)

    if (allocated(products%twist_total)) deallocate(products%twist_total)
    if (allocated(products%twist_backward)) deallocate(products%twist_backward)
    if (allocated(products%twist_forward)) deallocate(products%twist_forward)
    if (allocated(products%nstep_backward_twist)) &
         deallocate(products%nstep_backward_twist)
    if (allocated(products%nstep_forward_twist)) &
         deallocate(products%nstep_forward_twist)
    if (allocated(products%status_backward_twist)) &
         deallocate(products%status_backward_twist)
    if (allocated(products%status_forward_twist)) &
         deallocate(products%status_forward_twist)

    if (allocated(products%qperp)) deallocate(products%qperp)
    if (allocated(products%logqperp)) deallocate(products%logqperp)
    if (allocated(products%N2)) deallocate(products%N2)
    if (allocated(products%bfactor)) deallocate(products%bfactor)
    if (allocated(products%length_forward_qperp)) &
         deallocate(products%length_forward_qperp)
    if (allocated(products%length_backward_qperp)) &
         deallocate(products%length_backward_qperp)
    if (allocated(products%Bseed_norm)) deallocate(products%Bseed_norm)
    if (allocated(products%Bf_norm)) deallocate(products%Bf_norm)
    if (allocated(products%Bb_norm)) deallocate(products%Bb_norm)
    if (allocated(products%valid_qperp)) deallocate(products%valid_qperp)
    if (allocated(products%status_qperp)) deallocate(products%status_qperp)
    if (allocated(products%face_forward_qperp)) &
         deallocate(products%face_forward_qperp)
    if (allocated(products%face_backward_qperp)) &
         deallocate(products%face_backward_qperp)
    if (allocated(products%status_forward_qperp)) &
         deallocate(products%status_forward_qperp)
    if (allocated(products%status_backward_qperp)) &
         deallocate(products%status_backward_qperp)
  end subroutine mt_deallocate_volume_products

  subroutine mt_build_volume_slab_seeds(seeds,xmin,ymin,zmin,spacing, &
       nx,ny,k_start,slab_nz)
    double precision, intent(out) :: seeds(:,:)
    double precision, intent(in) :: xmin,ymin,zmin,spacing(3)
    integer, intent(in) :: nx,ny,k_start,slab_nz

    integer :: i,j,k,kk,iseed

    do kk=1,slab_nz
      k=k_start+kk-1
      do j=1,ny
        do i=1,nx
          iseed=(kk-1)*nx*ny+(j-1)*nx+i
          seeds(iseed,1)=xmin+dble(i-1)*spacing(1)
          seeds(iseed,2)=ymin+dble(j-1)*spacing(2)
          seeds(iseed,3)=zmin+dble(k-1)*spacing(3)
        enddo
      enddo
    enddo
  end subroutine mt_build_volume_slab_seeds

  subroutine mt_copy_volume_length_slab(products,results,nx,ny,k_start, &
       slab_nz)
    type(mt_volume_products), intent(inout) :: products
    type(trace_length_result), intent(in) :: results(:)
    integer, intent(in) :: nx,ny,k_start,slab_nz

    integer :: i,j,kk,k,ilocal,iglobal

    do kk=1,slab_nz
      k=k_start+kk-1
      do j=1,ny
        do i=1,nx
          ilocal=(kk-1)*nx*ny+(j-1)*nx+i
          iglobal=(k-1)*nx*ny+(j-1)*nx+i
          products%length_total(iglobal)=results(ilocal)%total_length
          products%length_backward(iglobal)=results(ilocal)%backward_length
          products%length_forward(iglobal)=results(ilocal)%forward_length
          products%nstep_backward_length(iglobal)= &
               results(ilocal)%backward_nstep
          products%nstep_forward_length(iglobal)= &
               results(ilocal)%forward_nstep
          products%status_backward_length(iglobal)= &
               results(ilocal)%backward_status
          products%status_forward_length(iglobal)= &
               results(ilocal)%forward_status
        enddo
      enddo
    enddo
  end subroutine mt_copy_volume_length_slab

  subroutine mt_copy_volume_twist_slab(products,results,nx,ny,k_start, &
       slab_nz)
    type(mt_volume_products), intent(inout) :: products
    type(trace_twist_result), intent(in) :: results(:)
    integer, intent(in) :: nx,ny,k_start,slab_nz

    integer :: i,j,kk,k,ilocal,iglobal

    do kk=1,slab_nz
      k=k_start+kk-1
      do j=1,ny
        do i=1,nx
          ilocal=(kk-1)*nx*ny+(j-1)*nx+i
          iglobal=(k-1)*nx*ny+(j-1)*nx+i
          products%twist_total(iglobal)=results(ilocal)%total_twist
          products%twist_backward(iglobal)=results(ilocal)%backward_twist
          products%twist_forward(iglobal)=results(ilocal)%forward_twist
          products%nstep_backward_twist(iglobal)= &
               results(ilocal)%line%backward_nstep
          products%nstep_forward_twist(iglobal)= &
               results(ilocal)%line%forward_nstep
          products%status_backward_twist(iglobal)= &
               results(ilocal)%line%backward_status
          products%status_forward_twist(iglobal)= &
               results(ilocal)%line%forward_status
        enddo
      enddo
    enddo
  end subroutine mt_copy_volume_twist_slab

  subroutine mt_copy_volume_qperp_slab(products,results,nx,ny,k_start, &
       slab_nz)
    type(mt_volume_products), intent(inout) :: products
    type(trace_qperp_result), intent(in) :: results(:)
    integer, intent(in) :: nx,ny,k_start,slab_nz

    integer :: i,j,kk,k,ilocal,iglobal

    do kk=1,slab_nz
      k=k_start+kk-1
      do j=1,ny
        do i=1,nx
          ilocal=(kk-1)*nx*ny+(j-1)*nx+i
          iglobal=(k-1)*nx*ny+(j-1)*nx+i
          products%qperp(iglobal)=results(ilocal)%qperp
          products%logqperp(iglobal)=results(ilocal)%logqperp
          products%N2(iglobal)=results(ilocal)%N2
          products%bfactor(iglobal)=results(ilocal)%bfactor
          products%length_forward_qperp(iglobal)= &
               results(ilocal)%forward_length
          products%length_backward_qperp(iglobal)= &
               results(ilocal)%backward_length
          products%Bseed_norm(iglobal)=dsqrt(sum(results(ilocal)%B_seed**2))
          products%Bf_norm(iglobal)=dsqrt(sum(results(ilocal)%forward_B**2))
          products%Bb_norm(iglobal)=dsqrt(sum(results(ilocal)%backward_B**2))
          products%valid_qperp(iglobal)=merge(1,0,results(ilocal)%valid)
          products%status_qperp(iglobal)=results(ilocal)%status
          products%face_forward_qperp(iglobal)=results(ilocal)%forward_face
          products%face_backward_qperp(iglobal)=results(ilocal)%backward_face
          products%status_forward_qperp(iglobal)= &
               results(ilocal)%forward_status
          products%status_backward_qperp(iglobal)= &
               results(ilocal)%backward_status
        enddo
      enddo
    enddo
  end subroutine mt_copy_volume_qperp_slab

  subroutine mt_build_vti_desc(desc,name,array_kind,npoint)
    type(mt_vti_array_desc), intent(out) :: desc
    character(len=*), intent(in) :: name
    integer, intent(in) :: array_kind,npoint

    desc%name=''
    desc%name=trim(name)
    desc%kind=array_kind
    select case (array_kind)
    case (mt_vti_kind_float64)
      desc%nbytes=int(npoint,kind=8)*8_8
    case (mt_vti_kind_int32)
      desc%nbytes=int(npoint,kind=8)*4_8
    case default
      desc%nbytes=0_8
    end select
    desc%offset=0_8
  end subroutine mt_build_vti_desc

  subroutine mt_finalize_vti_desc_offsets(descs,ndesc)
    integer, intent(in) :: ndesc
    type(mt_vti_array_desc), intent(inout) :: descs(ndesc)

    integer :: idesc
    integer(kind=8) :: offset

    offset=0_8
    do idesc=1,ndesc
      descs(idesc)%offset=offset
      offset=offset+descs(idesc)%nbytes+4_8
    enddo
  end subroutine mt_finalize_vti_desc_offsets

  subroutine mt_append_vti_desc(descs,idesc,name,array_kind,npoint)
    type(mt_vti_array_desc), intent(inout) :: descs(:)
    integer, intent(inout) :: idesc
    character(len=*), intent(in) :: name
    integer, intent(in) :: array_kind,npoint

    idesc=idesc+1
    call mt_build_vti_desc(descs(idesc),name,array_kind,npoint)
  end subroutine mt_append_vti_desc

  subroutine mt_build_cartesian_vti_pointdata_descs(npoint,descs,ndesc)
    integer, intent(in) :: npoint
    type(mt_vti_array_desc), allocatable, intent(out) :: descs(:)
    integer, intent(out) :: ndesc

    ndesc=3
    allocate(descs(ndesc))
    call mt_build_vti_desc(descs(1),'length_total',mt_vti_kind_float64, &
         npoint)
    call mt_build_vti_desc(descs(2),'qperp',mt_vti_kind_float64,npoint)
    call mt_build_vti_desc(descs(3),'status',mt_vti_kind_int32,npoint)
    call mt_finalize_vti_desc_offsets(descs,ndesc)
  end subroutine mt_build_cartesian_vti_pointdata_descs

  subroutine mt_build_volume_vti_descs(npoint,do_twist,do_qperp,descs, &
       ndesc)
    integer, intent(in) :: npoint
    logical, intent(in) :: do_twist,do_qperp
    type(mt_vti_array_desc), allocatable, intent(out) :: descs(:)
    integer, intent(out) :: ndesc

    integer :: idesc,max_desc

    if (.not.mt_vtk_detail_is_full()) then
      max_desc=1
      if (do_twist) max_desc=max_desc+1
      if (do_qperp) max_desc=max_desc+3
    else
      max_desc=7
      if (do_twist) max_desc=max_desc+7
      if (do_qperp) max_desc=max_desc+15
    endif
    allocate(descs(max_desc))

    idesc=0
    call mt_append_vti_desc(descs,idesc,'length_total', &
         mt_vti_kind_float64,npoint)
    if (.not.mt_vtk_detail_is_full()) then
      if (do_twist) then
        call mt_append_vti_desc(descs,idesc,'twist_total', &
             mt_vti_kind_float64,npoint)
      endif
      if (do_qperp) then
        call mt_append_vti_desc(descs,idesc,'logQperp', &
             mt_vti_kind_float64,npoint)
        call mt_append_vti_desc(descs,idesc,'valid_Qperp', &
             mt_vti_kind_int32,npoint)
        call mt_append_vti_desc(descs,idesc,'status_Qperp', &
             mt_vti_kind_int32,npoint)
      endif
      ndesc=idesc
      call mt_finalize_vti_desc_offsets(descs,ndesc)
      return
    endif

    call mt_append_vti_desc(descs,idesc,'length_backward', &
         mt_vti_kind_float64,npoint)
    call mt_append_vti_desc(descs,idesc,'length_forward', &
         mt_vti_kind_float64,npoint)
    call mt_append_vti_desc(descs,idesc,'nstep_backward_length', &
         mt_vti_kind_int32,npoint)
    call mt_append_vti_desc(descs,idesc,'nstep_forward_length', &
         mt_vti_kind_int32,npoint)
    call mt_append_vti_desc(descs,idesc,'status_backward_length', &
         mt_vti_kind_int32,npoint)
    call mt_append_vti_desc(descs,idesc,'status_forward_length', &
         mt_vti_kind_int32,npoint)

    if (do_twist) then
      call mt_append_vti_desc(descs,idesc,'twist_total', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'twist_backward', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'twist_forward', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'nstep_backward_twist', &
           mt_vti_kind_int32,npoint)
      call mt_append_vti_desc(descs,idesc,'nstep_forward_twist', &
           mt_vti_kind_int32,npoint)
      call mt_append_vti_desc(descs,idesc,'status_backward_twist', &
           mt_vti_kind_int32,npoint)
      call mt_append_vti_desc(descs,idesc,'status_forward_twist', &
           mt_vti_kind_int32,npoint)
    endif

    if (do_qperp) then
      call mt_append_vti_desc(descs,idesc,'qperp', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'logqperp', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'N2',mt_vti_kind_float64, &
           npoint)
      call mt_append_vti_desc(descs,idesc,'bfactor', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'length_forward_qperp', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'length_backward_qperp', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'Bseed_norm', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'Bf_norm', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'Bb_norm', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'valid_qperp', &
           mt_vti_kind_int32,npoint)
      call mt_append_vti_desc(descs,idesc,'status_qperp', &
           mt_vti_kind_int32,npoint)
      call mt_append_vti_desc(descs,idesc,'face_forward_qperp', &
           mt_vti_kind_int32,npoint)
      call mt_append_vti_desc(descs,idesc,'face_backward_qperp', &
           mt_vti_kind_int32,npoint)
      call mt_append_vti_desc(descs,idesc,'status_forward_qperp', &
           mt_vti_kind_int32,npoint)
      call mt_append_vti_desc(descs,idesc,'status_backward_qperp', &
           mt_vti_kind_int32,npoint)
    endif

    ndesc=idesc
    call mt_finalize_vti_desc_offsets(descs,ndesc)
  end subroutine mt_build_volume_vti_descs

  character(len=8) function mt_vti_type_name(array_kind)
    integer, intent(in) :: array_kind

    select case (array_kind)
    case (mt_vti_kind_float64)
      mt_vti_type_name='Float64'
    case (mt_vti_kind_int32)
      mt_vti_type_name='Int32'
    case default
      mt_vti_type_name='Unknown'
    end select
  end function mt_vti_type_name

  subroutine mt_write_vti_image_header(vti_unit,origin,spacing,nx,ny,nz, &
       descs,ndesc,io_status)
    integer, intent(in) :: vti_unit,nx,ny,nz,ndesc
    double precision, intent(in) :: origin(3),spacing(3)
    type(mt_vti_array_desc), intent(in) :: descs(ndesc)
    integer, intent(out) :: io_status

    integer :: extent(6),idesc

    extent=(/ 0,nx-1,0,ny-1,0,nz-1 /)
    write(vti_unit,'(a)',iostat=io_status) '<?xml version="1.0"?>'
    if (io_status/=0) return
    write(vti_unit,'(a)',iostat=io_status) &
         '<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">'
    if (io_status/=0) return
    write(vti_unit,'(a,3(1pe24.16),a,6(i0,1x),a,3(1pe24.16),a)', &
         iostat=io_status) '  <ImageData Origin="',origin, &
         '" WholeExtent="',extent,'" Spacing="',spacing,'">'
    if (io_status/=0) return
    write(vti_unit,'(a,6(i0,1x),a)',iostat=io_status) &
         '    <Piece Extent="',extent,'">'
    if (io_status/=0) return
    write(vti_unit,'(a)',iostat=io_status) '      <PointData>'
    if (io_status/=0) return

    do idesc=1,ndesc
      call mt_write_vti_appended_array(vti_unit, &
           mt_vti_type_name(descs(idesc)%kind),trim(descs(idesc)%name), &
           descs(idesc)%offset,io_status)
      if (io_status/=0) return
    enddo

    write(vti_unit,'(a)',iostat=io_status) '      </PointData>'
    if (io_status/=0) return
    write(vti_unit,'(a)',iostat=io_status) '    </Piece>'
    if (io_status/=0) return
    write(vti_unit,'(a)',iostat=io_status) '  </ImageData>'
    if (io_status/=0) return
    write(vti_unit,'(a)',iostat=io_status) '<AppendedData encoding="raw">'
  end subroutine mt_write_vti_image_header

  subroutine mt_write_fieldline_products_volume_vti(vti_file,origin, &
       spacing,nx,ny,nz,products,do_twist,do_qperp,caller)
    integer, intent(in) :: nx,ny,nz
    double precision, intent(in) :: origin(3),spacing(3)
    type(mt_volume_products), intent(in) :: products
    logical, intent(in) :: do_twist,do_qperp
    character(len=*), intent(in) :: vti_file,caller

    integer :: vti_unit,io_status,npoint
    character(len=1) :: marker
    type(mt_vti_array_desc), allocatable :: descs(:)
    integer :: ndesc

    npoint=nx*ny*nz
    call mt_check_vti_byte_count(npoint,caller)
    call mt_build_volume_vti_descs(npoint,do_twist,do_qperp,descs,ndesc)

    open(newunit=vti_unit,file=trim(vti_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open VTI file')
    endif

    call mt_write_vti_image_header(vti_unit,origin,spacing,nx,ny,nz, &
         descs,ndesc,io_status)
    if (io_status/=0) then
      close(vti_unit)
      call mpistop(trim(caller)//' could not write VTI header')
    endif
    close(vti_unit)

    open(newunit=vti_unit,file=trim(vti_file),access='stream', &
         form='unformatted',status='old',position='append', &
         action='write',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not append VTI payload')
    endif

    marker='_'
    write(vti_unit,iostat=io_status) marker
    if (io_status==0) then
      call mt_write_volume_vti_payload(vti_unit,products,npoint,descs, &
           ndesc,io_status)
    endif
    close(vti_unit)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not write VTI payload')
    endif

    open(newunit=vti_unit,file=trim(vti_file),status='old', &
         action='write',form='formatted',position='append',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not append VTI footer')
    endif
    write(vti_unit,'(a)',iostat=io_status) '</AppendedData>'
    if (io_status==0) write(vti_unit,'(a)',iostat=io_status) '</VTKFile>'
    close(vti_unit)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not write VTI footer')
    endif
    deallocate(descs)
  end subroutine mt_write_fieldline_products_volume_vti

  subroutine mt_write_volume_vti_payload(vti_unit,products,npoint,descs, &
       ndesc,io_status)
    integer, intent(in) :: vti_unit,npoint,ndesc
    type(mt_volume_products), intent(in) :: products
    type(mt_vti_array_desc), intent(in) :: descs(ndesc)
    integer, intent(inout) :: io_status

    integer :: idesc

    do idesc=1,ndesc
      select case (trim(descs(idesc)%name))
      case ('length_total')
        if (mt_vtk_detail_is_full()) then
          call mt_write_vti_payload_float64(vti_unit,products%length_total, &
               npoint,io_status)
        else
          call mt_write_vti_payload_float64_visual(vti_unit, &
               products%length_total,npoint,io_status)
        endif
      case ('length_backward')
        call mt_write_vti_payload_float64(vti_unit,products%length_backward, &
             npoint,io_status)
      case ('length_forward')
        call mt_write_vti_payload_float64(vti_unit,products%length_forward, &
             npoint,io_status)
      case ('nstep_backward_length')
        call mt_write_vti_payload_int32(vti_unit, &
             products%nstep_backward_length,npoint,io_status)
      case ('nstep_forward_length')
        call mt_write_vti_payload_int32(vti_unit, &
             products%nstep_forward_length,npoint,io_status)
      case ('status_backward_length')
        call mt_write_vti_payload_int32(vti_unit, &
             products%status_backward_length,npoint,io_status)
      case ('status_forward_length')
        call mt_write_vti_payload_int32(vti_unit, &
             products%status_forward_length,npoint,io_status)
      case ('twist_total')
        if (mt_vtk_detail_is_full()) then
          call mt_write_vti_payload_float64(vti_unit,products%twist_total, &
               npoint,io_status)
        else
          call mt_write_vti_payload_float64_visual(vti_unit, &
               products%twist_total,npoint,io_status)
        endif
      case ('twist_backward')
        call mt_write_vti_payload_float64(vti_unit,products%twist_backward, &
             npoint,io_status)
      case ('twist_forward')
        call mt_write_vti_payload_float64(vti_unit,products%twist_forward, &
             npoint,io_status)
      case ('nstep_backward_twist')
        call mt_write_vti_payload_int32(vti_unit, &
             products%nstep_backward_twist,npoint,io_status)
      case ('nstep_forward_twist')
        call mt_write_vti_payload_int32(vti_unit, &
             products%nstep_forward_twist,npoint,io_status)
      case ('status_backward_twist')
        call mt_write_vti_payload_int32(vti_unit, &
             products%status_backward_twist,npoint,io_status)
      case ('status_forward_twist')
        call mt_write_vti_payload_int32(vti_unit, &
             products%status_forward_twist,npoint,io_status)
      case ('qperp')
        call mt_write_vti_payload_float64(vti_unit,products%qperp,npoint, &
             io_status)
      case ('logqperp')
        call mt_write_vti_payload_float64(vti_unit,products%logqperp, &
             npoint,io_status)
      case ('logQperp')
        call mt_write_vti_payload_logqperp_visual(vti_unit,products,npoint, &
             io_status)
      case ('N2')
        call mt_write_vti_payload_float64(vti_unit,products%N2,npoint, &
             io_status)
      case ('bfactor')
        call mt_write_vti_payload_float64(vti_unit,products%bfactor,npoint, &
             io_status)
      case ('length_forward_qperp')
        call mt_write_vti_payload_float64(vti_unit, &
             products%length_forward_qperp,npoint,io_status)
      case ('length_backward_qperp')
        call mt_write_vti_payload_float64(vti_unit, &
             products%length_backward_qperp,npoint,io_status)
      case ('Bseed_norm')
        call mt_write_vti_payload_float64(vti_unit,products%Bseed_norm, &
             npoint,io_status)
      case ('Bf_norm')
        call mt_write_vti_payload_float64(vti_unit,products%Bf_norm,npoint, &
             io_status)
      case ('Bb_norm')
        call mt_write_vti_payload_float64(vti_unit,products%Bb_norm,npoint, &
             io_status)
      case ('valid_qperp')
        call mt_write_vti_payload_int32(vti_unit,products%valid_qperp, &
             npoint,io_status)
      case ('valid_Qperp')
        call mt_write_vti_payload_int32(vti_unit,products%valid_qperp, &
             npoint,io_status)
      case ('status_qperp')
        call mt_write_vti_payload_int32(vti_unit,products%status_qperp, &
             npoint,io_status)
      case ('status_Qperp')
        call mt_write_vti_payload_int32(vti_unit,products%status_qperp, &
             npoint,io_status)
      case ('face_forward_qperp')
        call mt_write_vti_payload_int32(vti_unit, &
             products%face_forward_qperp,npoint,io_status)
      case ('face_backward_qperp')
        call mt_write_vti_payload_int32(vti_unit, &
             products%face_backward_qperp,npoint,io_status)
      case ('status_forward_qperp')
        call mt_write_vti_payload_int32(vti_unit, &
             products%status_forward_qperp,npoint,io_status)
      case ('status_backward_qperp')
        call mt_write_vti_payload_int32(vti_unit, &
             products%status_backward_qperp,npoint,io_status)
      case default
        io_status=1
      end select
      if (io_status/=0) exit
    enddo
  end subroutine mt_write_volume_vti_payload

  subroutine mt_write_vti_payload_float64_visual(vti_unit,values,npoint, &
       io_status)
    integer, intent(in) :: vti_unit,npoint
    double precision, intent(in) :: values(npoint)
    integer, intent(inout) :: io_status

    double precision, allocatable :: visual_values(:)
    integer :: ipoint

    if (io_status/=0) return
    allocate(visual_values(npoint))
    do ipoint=1,npoint
      visual_values(ipoint)=mt_visual_float(values(ipoint))
    enddo
    call mt_write_vti_payload_float64(vti_unit,visual_values,npoint, &
         io_status)
    deallocate(visual_values)
  end subroutine mt_write_vti_payload_float64_visual

  subroutine mt_write_vti_payload_logqperp_visual(vti_unit,products,npoint, &
       io_status)
    integer, intent(in) :: vti_unit,npoint
    type(mt_volume_products), intent(in) :: products
    integer, intent(inout) :: io_status

    double precision, allocatable :: visual_values(:)
    integer :: ipoint

    if (io_status/=0) return
    allocate(visual_values(npoint))
    do ipoint=1,npoint
      visual_values(ipoint)=mt_visual_valid_float(products%logqperp(ipoint), &
           products%valid_qperp(ipoint)==1)
    enddo
    call mt_write_vti_payload_float64(vti_unit,visual_values,npoint, &
         io_status)
    deallocate(visual_values)
  end subroutine mt_write_vti_payload_logqperp_visual

  subroutine mt_check_vti_byte_count(npoint,caller)
    integer, intent(in) :: npoint
    character(len=*), intent(in) :: caller

    if (int(npoint,kind=8)*8_8>int(huge(0),kind=8)) then
      call mpistop(trim(caller)//' exceeds 32-bit VTI byte count')
    endif
  end subroutine mt_check_vti_byte_count

  subroutine mt_write_vti_payload_float64(vti_unit,values,npoint,io_status)
    integer, intent(in) :: vti_unit,npoint
    double precision, intent(in) :: values(npoint)
    integer, intent(inout) :: io_status

    integer :: byte_count

    if (io_status/=0) return
    byte_count=npoint*8
    write(vti_unit,iostat=io_status) byte_count
    if (io_status==0) write(vti_unit,iostat=io_status) values
  end subroutine mt_write_vti_payload_float64

  subroutine mt_write_vti_payload_int32(vti_unit,values,npoint,io_status)
    integer, intent(in) :: vti_unit,npoint
    integer, intent(in) :: values(npoint)
    integer, intent(inout) :: io_status

    integer :: byte_count

    if (io_status/=0) return
    byte_count=npoint*4
    write(vti_unit,iostat=io_status) byte_count
    if (io_status==0) write(vti_unit,iostat=io_status) values
  end subroutine mt_write_vti_payload_int32

  subroutine mt_write_vti_pointdata_fixed(vti_file,origin,spacing,nx,ny,nz, &
       length_total,qperp,status,caller)
    character(len=*), intent(in) :: vti_file,caller
    integer, intent(in) :: nx,ny,nz
    double precision, intent(in) :: origin(3),spacing(3)
    double precision, intent(in) :: length_total(nx*ny*nz),qperp(nx*ny*nz)
    integer, intent(in) :: status(nx*ny*nz)

    integer :: vti_unit,io_status,npoint
    integer :: idesc,ndesc
    character(len=1) :: marker
    type(mt_vti_array_desc), allocatable :: descs(:)

    npoint=nx*ny*nz
    call mt_check_vti_byte_count(npoint,caller)
    call mt_build_cartesian_vti_pointdata_descs(npoint,descs,ndesc)

    open(newunit=vti_unit,file=trim(vti_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open VTI file')
    endif

    call mt_write_vti_image_header(vti_unit,origin,spacing,nx,ny,nz, &
         descs,ndesc,io_status)
    if (io_status/=0) then
      close(vti_unit)
      call mpistop(trim(caller)//' could not write VTI header')
    endif
    close(vti_unit)

    open(newunit=vti_unit,file=trim(vti_file),access='stream', &
         form='unformatted',status='old',position='append', &
         action='write',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not append VTI payload')
    endif

    marker='_'
    write(vti_unit,iostat=io_status) marker
    do idesc=1,ndesc
      if (io_status/=0) exit
      select case (trim(descs(idesc)%name))
      case ('length_total')
        call mt_write_vti_payload_float64(vti_unit,length_total,npoint, &
             io_status)
      case ('qperp')
        call mt_write_vti_payload_float64(vti_unit,qperp,npoint,io_status)
      case ('status')
        call mt_write_vti_payload_int32(vti_unit,status,npoint,io_status)
      case default
        io_status=1
      end select
    enddo
    close(vti_unit)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not write VTI payload')
    endif

    open(newunit=vti_unit,file=trim(vti_file),status='old', &
         action='write',form='formatted',position='append',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not append VTI footer')
    endif
    write(vti_unit,'(a)',iostat=io_status) '</AppendedData>'
    if (io_status==0) write(vti_unit,'(a)',iostat=io_status) '</VTKFile>'
    close(vti_unit)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not write VTI footer')
    endif
    deallocate(descs)
  end subroutine mt_write_vti_pointdata_fixed

  subroutine mt_write_vti_appended_array(vti_unit,vtk_type,name,offset, &
       io_status)
    integer, intent(in) :: vti_unit
    character(len=*), intent(in) :: vtk_type,name
    integer(kind=8), intent(in) :: offset
    integer, intent(out) :: io_status

    write(vti_unit,'(a)',advance='no',iostat=io_status) &
         '        <DataArray type="'
    if (io_status==0) write(vti_unit,'(a)',advance='no', &
         iostat=io_status) trim(vtk_type)
    if (io_status==0) write(vti_unit,'(a)',advance='no', &
         iostat=io_status) '" Name="'
    if (io_status==0) write(vti_unit,'(a)',advance='no', &
         iostat=io_status) trim(name)
    if (io_status==0) write(vti_unit,'(a)',advance='no', &
         iostat=io_status) '" format="appended" offset="'
    if (io_status==0) write(vti_unit,'(i0)',advance='no', &
         iostat=io_status) offset
    if (io_status==0) write(vti_unit,'(a)',iostat=io_status) '"/>'
  end subroutine mt_write_vti_appended_array

  subroutine mt_write_vtu_length_pointdata(vtu_unit,length_results,npoint, &
       io_status)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_length_result), intent(in) :: length_results(npoint)
    integer, intent(inout) :: io_status

    integer :: ipoint

    call mt_write_vtu_float_array_start(vtu_unit,'length_total',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (length_results(ipoint)%total_length,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'length_backward', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (length_results(ipoint)%backward_length,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'length_forward', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (length_results(ipoint)%forward_length,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    call mt_write_vtu_int_array_start(vtu_unit,'nstep_backward_length', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (length_results(ipoint)%backward_nstep,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'nstep_forward_length', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (length_results(ipoint)%forward_nstep,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_backward_length', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (length_results(ipoint)%backward_status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_forward_length', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (length_results(ipoint)%forward_status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
  end subroutine mt_write_vtu_length_pointdata

  subroutine mt_write_vtu_mapping_pointdata(vtu_unit,mapping_results, &
       npoint,io_status)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_mapping_result), intent(in) :: mapping_results(npoint)
    integer, intent(inout) :: io_status

    integer :: ipoint

    call mt_write_vtu_float_array_start(vtu_unit,'x_f_mapping',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mapping_results(ipoint)%forward_footpoint(1),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'y_f_mapping',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mapping_results(ipoint)%forward_footpoint(2),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'z_f_mapping',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mapping_results(ipoint)%forward_footpoint(3),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'x_b_mapping',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mapping_results(ipoint)%backward_footpoint(1),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'y_b_mapping',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mapping_results(ipoint)%backward_footpoint(2),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'z_b_mapping',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mapping_results(ipoint)%backward_footpoint(3),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    call mt_write_vtu_float_array_start(vtu_unit,'source_Bn_mapping', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mapping_results(ipoint)%source_Bn,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'forward_Bn_mapping', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mapping_results(ipoint)%forward_Bn,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'backward_Bn_mapping', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mapping_results(ipoint)%backward_Bn,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    call mt_write_vtu_int_array_start(vtu_unit,'face_forward_mapping', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (mapping_results(ipoint)%forward_face,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'face_backward_mapping', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (mapping_results(ipoint)%backward_face,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_forward_mapping', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (mapping_results(ipoint)%forward_status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_backward_mapping', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (mapping_results(ipoint)%backward_status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'valid_mapping',io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (merge(1,0,mapping_results(ipoint)%valid),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
  end subroutine mt_write_vtu_mapping_pointdata

  subroutine mt_write_vtu_qperp_public_pointdata(vtu_unit,qperp_results, &
       npoint,io_status)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_qperp_result), intent(in) :: qperp_results(npoint)
    integer, intent(inout) :: io_status

    integer :: ipoint

    call mt_write_vtu_float_array_start(vtu_unit,'logQperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_visual_valid_float(qperp_results(ipoint)%logqperp, &
         qperp_results(ipoint)%valid),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'valid_Qperp',io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (merge(1,0,qperp_results(ipoint)%valid),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_Qperp',io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (qperp_results(ipoint)%status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
  end subroutine mt_write_vtu_qperp_public_pointdata

  subroutine mt_write_vtu_qperp_method2_pointdata(vtu_unit,qperp_results, &
       npoint,io_status)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_qperp_result), intent(in) :: qperp_results(npoint)
    integer, intent(inout) :: io_status

    integer :: ipoint

    call mt_write_vtu_float_array_start(vtu_unit,'qperp_method2',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%qperp,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'logqperp_method2', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%logqperp,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'N2_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%N2,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'bfactor_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%bfactor,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'length_forward_qperp', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%forward_length,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'length_backward_qperp', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%backward_length,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    call mt_write_vtu_float_array_start(vtu_unit,'Bseed_norm_qperp', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (dsqrt(sum(qperp_results(ipoint)%B_seed**2)),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'Bf_norm_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (dsqrt(sum(qperp_results(ipoint)%forward_B**2)),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'Bb_norm_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (dsqrt(sum(qperp_results(ipoint)%backward_B**2)),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    call mt_write_vtu_float_array_start(vtu_unit,'x_f_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%forward_endpoint(1),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'y_f_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%forward_endpoint(2),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'z_f_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%forward_endpoint(3),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'x_b_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%backward_endpoint(1),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'y_b_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%backward_endpoint(2),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'z_b_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%backward_endpoint(3),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    call mt_write_vtu_int_array_start(vtu_unit,'valid_qperp_method2', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (merge(1,0,qperp_results(ipoint)%valid),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_qperp_method2', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (qperp_results(ipoint)%status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'face_forward_qperp', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (qperp_results(ipoint)%forward_face,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'face_backward_qperp', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (qperp_results(ipoint)%backward_face,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_forward_qperp', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (qperp_results(ipoint)%forward_status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_backward_qperp', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (qperp_results(ipoint)%backward_status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
  end subroutine mt_write_vtu_qperp_method2_pointdata

  subroutine mt_write_vtu_q_product_pointdata(vtu_unit,qperp_results, &
       npoint,write_raw_q,io_status)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_qperp_result), intent(in) :: qperp_results(npoint)
    logical, intent(in) :: write_raw_q
    integer, intent(inout) :: io_status

    integer :: ipoint

    if (write_raw_q) then
      call mt_write_vtu_float_array_start(vtu_unit,'q',io_status)
      if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
           (qperp_results(ipoint)%q0,ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
    endif

    call mt_write_vtu_float_array_start(vtu_unit,'logQ',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_visual_valid_float(qperp_results(ipoint)%logq0, &
         qperp_results(ipoint)%valid_q0),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'valid_Q',io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (merge(1,0,qperp_results(ipoint)%valid_q0),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_Q',io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (qperp_results(ipoint)%status_q0,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'face_pair_Q',io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (mt_q_result_face_pair(qperp_results(ipoint)),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'connection_type_Q', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (mt_q_result_connection_type(qperp_results(ipoint)),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
  end subroutine mt_write_vtu_q_product_pointdata

  subroutine mt_write_vtu_file_header(vtu_unit,npoint,ncell)
    integer, intent(in) :: vtu_unit,npoint,ncell

    write(vtu_unit,'(a)') '<?xml version="1.0"?>'
    write(vtu_unit,'(a)') &
         '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
    write(vtu_unit,'(a)') '<UnstructuredGrid>'
    write(vtu_unit,'(a,i0,a,i0,a)') '<Piece NumberOfPoints="',npoint, &
         '" NumberOfCells="',ncell,'">'
  end subroutine mt_write_vtu_file_header

  subroutine mt_write_vtu_file_footer(vtu_unit,io_status)
    integer, intent(in) :: vtu_unit
    integer, intent(inout) :: io_status

    write(vtu_unit,'(a)',iostat=io_status) '</Piece>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '</UnstructuredGrid>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '</VTKFile>'
  end subroutine mt_write_vtu_file_footer

  subroutine mt_write_vtu_product_pointdata(vtu_unit,length_results, &
       twist_results,qperp_results,npoint,do_twist,do_qperp,io_status)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_length_result), intent(in) :: length_results(npoint)
    type(trace_twist_result), intent(in) :: twist_results(:)
    type(trace_qperp_result), intent(in) :: qperp_results(:)
    logical, intent(in) :: do_twist,do_qperp
    integer, intent(inout) :: io_status

    integer :: ipoint

    write(vtu_unit,'(a)',iostat=io_status) '<PointData>'
    if (io_status/=0) return

    if (.not.mt_vtk_detail_is_full()) then
      call mt_write_vtu_product_minimal_pointdata(vtu_unit,length_results, &
           twist_results,qperp_results,npoint,do_twist,do_qperp,io_status)
      if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
           '</PointData>'
      return
    endif

    call mt_write_vtu_float_array_start(vtu_unit,'length_total',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (length_results(ipoint)%total_length,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'length_backward', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (length_results(ipoint)%backward_length,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'length_forward', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (length_results(ipoint)%forward_length,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    call mt_write_vtu_int_array_start(vtu_unit,'nstep_backward_length', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (length_results(ipoint)%backward_nstep,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'nstep_forward_length', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (length_results(ipoint)%forward_nstep,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_backward_length', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (length_results(ipoint)%backward_status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_forward_length', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (length_results(ipoint)%forward_status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    if (do_twist) then
      call mt_write_vtu_twist_pointdata(vtu_unit,twist_results,npoint, &
           io_status)
    endif
    if (do_qperp) then
      call mt_write_vtu_qperp_pointdata(vtu_unit,qperp_results,npoint, &
           io_status)
    endif

    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '</PointData>'
  end subroutine mt_write_vtu_product_pointdata

  subroutine mt_write_vtu_product_minimal_pointdata(vtu_unit, &
       length_results,twist_results,qperp_results,npoint,do_twist, &
       do_qperp,io_status)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_length_result), intent(in) :: length_results(npoint)
    type(trace_twist_result), intent(in) :: twist_results(:)
    type(trace_qperp_result), intent(in) :: qperp_results(:)
    logical, intent(in) :: do_twist,do_qperp
    integer, intent(inout) :: io_status

    integer :: ipoint

    call mt_write_vtu_float_array_start(vtu_unit,'length_total',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_visual_float(length_results(ipoint)%total_length), &
         ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    if (do_twist) then
      call mt_write_vtu_float_array_start(vtu_unit,'twist_total',io_status)
      if (io_status==0) write(vtu_unit,'(4(1pe24.16))', &
           iostat=io_status) &
           (mt_visual_float(twist_results(ipoint)%total_twist), &
           ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
    endif

    if (do_qperp) then
      call mt_write_vtu_float_array_start(vtu_unit,'logQperp',io_status)
      if (io_status==0) write(vtu_unit,'(4(1pe24.16))', &
           iostat=io_status) &
           (mt_visual_valid_float(qperp_results(ipoint)%logqperp, &
           qperp_results(ipoint)%valid),ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
      call mt_write_vtu_int_array_start(vtu_unit,'valid_Qperp', &
           io_status)
      if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
           (merge(1,0,qperp_results(ipoint)%valid),ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
      call mt_write_vtu_int_array_start(vtu_unit,'status_Qperp', &
           io_status)
      if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
           (qperp_results(ipoint)%status,ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
    endif
  end subroutine mt_write_vtu_product_minimal_pointdata

  subroutine mt_write_vtu_twist_pointdata(vtu_unit,twist_results,npoint, &
       io_status)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_twist_result), intent(in) :: twist_results(npoint)
    integer, intent(inout) :: io_status

    integer :: ipoint

    call mt_write_vtu_float_array_start(vtu_unit,'twist_total',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (twist_results(ipoint)%total_twist,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'twist_backward', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (twist_results(ipoint)%backward_twist,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'twist_forward', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (twist_results(ipoint)%forward_twist,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    call mt_write_vtu_int_array_start(vtu_unit,'nstep_backward_twist', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (twist_results(ipoint)%line%backward_nstep,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'nstep_forward_twist', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (twist_results(ipoint)%line%forward_nstep,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_backward_twist', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (twist_results(ipoint)%line%backward_status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_forward_twist', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (twist_results(ipoint)%line%forward_status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
  end subroutine mt_write_vtu_twist_pointdata

  subroutine mt_write_vtu_qperp_pointdata(vtu_unit,qperp_results,npoint, &
       io_status)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_qperp_result), intent(in) :: qperp_results(npoint)
    integer, intent(inout) :: io_status

    integer :: ipoint

    call mt_write_vtu_float_array_start(vtu_unit,'qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%qperp,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'logqperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%logqperp,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'N2',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%N2,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'bfactor',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%bfactor,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'length_forward_qperp', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%forward_length,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'length_backward_qperp', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%backward_length,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    call mt_write_vtu_float_array_start(vtu_unit,'Bseed_norm',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (dsqrt(sum(qperp_results(ipoint)%B_seed**2)),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'Bf_norm',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (dsqrt(sum(qperp_results(ipoint)%forward_B**2)),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'Bb_norm',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (dsqrt(sum(qperp_results(ipoint)%backward_B**2)),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    call mt_write_vtu_qperp_endpoint_pointdata(vtu_unit,qperp_results, &
         npoint,io_status)
    call mt_write_vtu_qperp_int_pointdata(vtu_unit,qperp_results,npoint, &
         io_status)
  end subroutine mt_write_vtu_qperp_pointdata

  subroutine mt_write_vtu_qperp_endpoint_pointdata(vtu_unit,qperp_results, &
       npoint,io_status)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_qperp_result), intent(in) :: qperp_results(npoint)
    integer, intent(inout) :: io_status

    integer :: ipoint

    call mt_write_vtu_float_array_start(vtu_unit,'x_f_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%forward_endpoint(1),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'y_f_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%forward_endpoint(2),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'z_f_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%forward_endpoint(3),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'x_b_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%backward_endpoint(1),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'y_b_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%backward_endpoint(2),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'z_b_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%backward_endpoint(3),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
  end subroutine mt_write_vtu_qperp_endpoint_pointdata

  subroutine mt_write_vtu_qperp_int_pointdata(vtu_unit,qperp_results, &
       npoint,io_status)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_qperp_result), intent(in) :: qperp_results(npoint)
    integer, intent(inout) :: io_status

    integer :: ipoint

    call mt_write_vtu_int_array_start(vtu_unit,'valid_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (merge(1,0,qperp_results(ipoint)%valid),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (qperp_results(ipoint)%status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'face_forward_qperp', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (qperp_results(ipoint)%forward_face,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'face_backward_qperp', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (qperp_results(ipoint)%backward_face,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_forward_qperp', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (qperp_results(ipoint)%forward_status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_backward_qperp', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (qperp_results(ipoint)%backward_status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
  end subroutine mt_write_vtu_qperp_int_pointdata

  subroutine mt_write_vtu_empty_pointdata(vtu_unit,do_twist,do_qperp)
    integer, intent(in) :: vtu_unit
    logical, intent(in) :: do_twist,do_qperp

    write(vtu_unit,'(a)') '<PointData>'
    if (.not.mt_vtk_detail_is_full()) then
      call mt_write_vtu_empty_float_array(vtu_unit,'length_total')
      if (do_twist) then
        call mt_write_vtu_empty_float_array(vtu_unit,'twist_total')
      endif
      if (do_qperp) then
        call mt_write_vtu_empty_float_array(vtu_unit,'logQperp')
        call mt_write_vtu_empty_int_array(vtu_unit,'valid_Qperp')
        call mt_write_vtu_empty_int_array(vtu_unit,'status_Qperp')
      endif
      write(vtu_unit,'(a)') '</PointData>'
      return
    endif
    call mt_write_vtu_empty_float_array(vtu_unit,'length_total')
    call mt_write_vtu_empty_float_array(vtu_unit,'length_backward')
    call mt_write_vtu_empty_float_array(vtu_unit,'length_forward')
    call mt_write_vtu_empty_int_array(vtu_unit,'nstep_backward_length')
    call mt_write_vtu_empty_int_array(vtu_unit,'nstep_forward_length')
    call mt_write_vtu_empty_int_array(vtu_unit,'status_backward_length')
    call mt_write_vtu_empty_int_array(vtu_unit,'status_forward_length')
    if (do_twist) then
      call mt_write_vtu_empty_float_array(vtu_unit,'twist_total')
      call mt_write_vtu_empty_float_array(vtu_unit,'twist_backward')
      call mt_write_vtu_empty_float_array(vtu_unit,'twist_forward')
      call mt_write_vtu_empty_int_array(vtu_unit,'nstep_backward_twist')
      call mt_write_vtu_empty_int_array(vtu_unit,'nstep_forward_twist')
      call mt_write_vtu_empty_int_array(vtu_unit,'status_backward_twist')
      call mt_write_vtu_empty_int_array(vtu_unit,'status_forward_twist')
    endif
    if (do_qperp) then
      call mt_write_vtu_empty_float_array(vtu_unit,'qperp')
      call mt_write_vtu_empty_float_array(vtu_unit,'logqperp')
      call mt_write_vtu_empty_float_array(vtu_unit,'N2')
      call mt_write_vtu_empty_float_array(vtu_unit,'bfactor')
      call mt_write_vtu_empty_float_array(vtu_unit,'length_forward_qperp')
      call mt_write_vtu_empty_float_array(vtu_unit,'length_backward_qperp')
      call mt_write_vtu_empty_float_array(vtu_unit,'Bseed_norm')
      call mt_write_vtu_empty_float_array(vtu_unit,'Bf_norm')
      call mt_write_vtu_empty_float_array(vtu_unit,'Bb_norm')
      call mt_write_vtu_empty_float_array(vtu_unit,'x_f_qperp')
      call mt_write_vtu_empty_float_array(vtu_unit,'y_f_qperp')
      call mt_write_vtu_empty_float_array(vtu_unit,'z_f_qperp')
      call mt_write_vtu_empty_float_array(vtu_unit,'x_b_qperp')
      call mt_write_vtu_empty_float_array(vtu_unit,'y_b_qperp')
      call mt_write_vtu_empty_float_array(vtu_unit,'z_b_qperp')
      call mt_write_vtu_empty_int_array(vtu_unit,'valid_qperp')
      call mt_write_vtu_empty_int_array(vtu_unit,'status_qperp')
      call mt_write_vtu_empty_int_array(vtu_unit,'face_forward_qperp')
      call mt_write_vtu_empty_int_array(vtu_unit,'face_backward_qperp')
      call mt_write_vtu_empty_int_array(vtu_unit,'status_forward_qperp')
      call mt_write_vtu_empty_int_array(vtu_unit,'status_backward_qperp')
    endif
    write(vtu_unit,'(a)') '</PointData>'
  end subroutine mt_write_vtu_empty_pointdata

  subroutine mt_write_vtu_empty_float_array(vtu_unit,name)
    integer, intent(in) :: vtu_unit
    character(len=*), intent(in) :: name

    write(vtu_unit,'(a,a,a)') '<DataArray type="Float64" Name="', &
         trim(name),'" format="ascii">'
    write(vtu_unit,'(a)') '</DataArray>'
  end subroutine mt_write_vtu_empty_float_array

  subroutine mt_write_vtu_empty_int_array(vtu_unit,name)
    integer, intent(in) :: vtu_unit
    character(len=*), intent(in) :: name

    write(vtu_unit,'(a,a,a)') '<DataArray type="Int32" Name="', &
         trim(name),'" format="ascii">'
    write(vtu_unit,'(a)') '</DataArray>'
  end subroutine mt_write_vtu_empty_int_array

  subroutine mt_write_vtu_float_array_start(vtu_unit,name,io_status)
    integer, intent(in) :: vtu_unit
    character(len=*), intent(in) :: name
    integer, intent(inout) :: io_status

    if (io_status/=0) return
    write(vtu_unit,'(a,a,a)',iostat=io_status) &
         '<DataArray type="Float64" Name="',trim(name),'" format="ascii">'
  end subroutine mt_write_vtu_float_array_start

  subroutine mt_write_vtu_int_array_start(vtu_unit,name,io_status)
    integer, intent(in) :: vtu_unit
    character(len=*), intent(in) :: name
    integer, intent(inout) :: io_status

    if (io_status/=0) return
    write(vtu_unit,'(a,a,a)',iostat=io_status) &
         '<DataArray type="Int32" Name="',trim(name),'" format="ascii">'
  end subroutine mt_write_vtu_int_array_start

  subroutine mt_write_vtu_data_array_end(vtu_unit,io_status)
    integer, intent(in) :: vtu_unit
    integer, intent(inout) :: io_status

    if (io_status/=0) return
    write(vtu_unit,'(a)',iostat=io_status) '</DataArray>'
  end subroutine mt_write_vtu_data_array_end

  subroutine mt_write_vtu_points(vtu_unit,length_results,npoint,io_status)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_length_result), intent(in) :: length_results(npoint)
    integer, intent(inout) :: io_status

    double precision :: seed_xyz(3)
    integer :: ipoint

    write(vtu_unit,'(a)',iostat=io_status) '<Points>'
    if (io_status==0) then
      write(vtu_unit,'(a)',iostat=io_status) &
           '<DataArray type="Float64" NumberOfComponents="3" format="ascii">'
    endif
    do ipoint=1,npoint
      seed_xyz=0.d0
      seed_xyz(1:ndim)=length_results(ipoint)%seed
      if (io_status==0) write(vtu_unit,'(3(1pe24.16))', &
           iostat=io_status) seed_xyz
    enddo
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '</DataArray>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '</Points>'
  end subroutine mt_write_vtu_points

  subroutine mt_write_vtu_vertex_cells(vtu_unit,npoint,io_status)
    integer, intent(in) :: vtu_unit,npoint
    integer, intent(inout) :: io_status

    integer :: ipoint

    write(vtu_unit,'(a)',iostat=io_status) '<Cells>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '<DataArray type="Int32" Name="connectivity" format="ascii">'
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (ipoint-1,ipoint=1,npoint)
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '</DataArray>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '<DataArray type="Int32" Name="offsets" format="ascii">'
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (ipoint,ipoint=1,npoint)
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '</DataArray>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '<DataArray type="UInt8" Name="types" format="ascii">'
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (1,ipoint=1,npoint)
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '</DataArray>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '</Cells>'
  end subroutine mt_write_vtu_vertex_cells

  subroutine mt_write_vtu_quad_cells(vtu_unit,n1,n2,io_status)
    integer, intent(in) :: vtu_unit,n1,n2
    integer, intent(inout) :: io_status

    integer :: i,j,icell,p,ncell

    ncell=(n1-1)*(n2-1)
    write(vtu_unit,'(a)',iostat=io_status) '<Cells>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '<DataArray type="Int32" Name="connectivity" format="ascii">'
    do j=1,n2-1
      do i=1,n1-1
        p=(j-1)*n1+(i-1)
        if (io_status==0) write(vtu_unit,'(4(i0,1x))',iostat=io_status) &
             p,p+1,p+1+n1,p+n1
      enddo
    enddo
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '</DataArray>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '<DataArray type="Int32" Name="offsets" format="ascii">'
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (4*icell,icell=1,ncell)
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '</DataArray>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '<DataArray type="UInt8" Name="types" format="ascii">'
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (9,icell=1,ncell)
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '</DataArray>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '</Cells>'
  end subroutine mt_write_vtu_quad_cells

  subroutine mt_write_vtu_cells_empty(vtu_unit,io_status)
    integer, intent(in) :: vtu_unit
    integer, intent(inout) :: io_status

    write(vtu_unit,'(a)',iostat=io_status) '<Cells>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '<DataArray type="Int32" Name="connectivity" format="ascii">'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '</DataArray>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '<DataArray type="Int32" Name="offsets" format="ascii">'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '</DataArray>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '<DataArray type="UInt8" Name="types" format="ascii">'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '</DataArray>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '</Cells>'
  end subroutine mt_write_vtu_cells_empty

  subroutine mt_write_twist_plane_csv(results,n1,n2,csv_file,caller, &
       index_header)
    integer, intent(in) :: n1,n2
    type(trace_twist_result), intent(in) :: results(n1*n2)
    character(len=*), intent(in) :: csv_file,caller,index_header

    double precision :: seed_xyz(3)
    integer :: csv_unit,io_status,i,j,iseed

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open CSV file')
    endif

    write(csv_unit,'(a)',iostat=io_status) &
         trim(index_header)//',seed_x,seed_y,seed_z,'// &
         'length_total,length_backward,length_forward,'// &
         'twist_total,twist_backward,twist_forward,'// &
         'nstep_backward,nstep_forward,'// &
         'status_backward,status_forward'
    if (io_status/=0) then
      close(csv_unit)
      call mpistop(trim(caller)//' could not write CSV header')
    endif

    do j=1,n2
      do i=1,n1
        iseed=(j-1)*n1+i
        seed_xyz=0.d0
        seed_xyz(1:ndim)=results(iseed)%line%seed
        write(csv_unit,'(i0,",",i0,9(",",es24.16),4(",",i0))', &
             iostat=io_status) i,j,seed_xyz, &
             results(iseed)%line%total_length, &
             results(iseed)%line%backward_length, &
             results(iseed)%line%forward_length, &
             results(iseed)%total_twist, &
             results(iseed)%backward_twist, &
             results(iseed)%forward_twist, &
             results(iseed)%line%backward_nstep, &
             results(iseed)%line%forward_nstep, &
             results(iseed)%line%backward_status, &
             results(iseed)%line%forward_status
        if (io_status/=0) then
          close(csv_unit)
          call mpistop(trim(caller)//' could not write CSV data')
        endif
      enddo
    enddo

    close(csv_unit)
  end subroutine mt_write_twist_plane_csv

  subroutine mt_write_qperp_plane_csv(results,n1,n2,csv_file,caller, &
       index_header)
    integer, intent(in) :: n1,n2
    type(trace_qperp_result), intent(in) :: results(n1*n2)
    character(len=*), intent(in) :: csv_file,caller,index_header

    double precision :: Bseed_norm,Bf_norm,Bb_norm
    integer :: csv_unit,io_status,i,j,iseed

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open CSV file')
    endif

    write(csv_unit,'(a)',iostat=io_status) &
         trim(index_header)//',seed_x,seed_y,seed_z,'// &
         'qperp,logqperp,valid,status,'// &
         'N2,bfactor,'// &
         'face_forward,face_backward,'// &
         'status_forward,status_backward,'// &
         'length_forward,length_backward,'// &
         'Bseed_norm,Bf_norm,Bb_norm,'// &
         'x_f,y_f,z_f,x_b,y_b,z_b'
    if (io_status/=0) then
      close(csv_unit)
      call mpistop(trim(caller)//' could not write CSV header')
    endif

    do j=1,n2
      do i=1,n1
        iseed=(j-1)*n1+i
        Bseed_norm=dsqrt(sum(results(iseed)%B_seed**2))
        Bf_norm=dsqrt(sum(results(iseed)%forward_B**2))
        Bb_norm=dsqrt(sum(results(iseed)%backward_B**2))
        write(csv_unit, &
             '(i0,",",i0,5(",",es24.16),",",l1,",",i0,'// &
             '2(",",es24.16),4(",",i0),11(",",es24.16))', &
             iostat=io_status) &
             i,j,results(iseed)%seed, &
             results(iseed)%qperp,results(iseed)%logqperp, &
             results(iseed)%valid,results(iseed)%status, &
             results(iseed)%N2,results(iseed)%bfactor, &
             results(iseed)%forward_face,results(iseed)%backward_face, &
             results(iseed)%forward_status,results(iseed)%backward_status, &
             results(iseed)%forward_length,results(iseed)%backward_length, &
             Bseed_norm,Bf_norm,Bb_norm, &
             results(iseed)%forward_endpoint, &
             results(iseed)%backward_endpoint
        if (io_status/=0) then
          close(csv_unit)
          call mpistop(trim(caller)//' could not write CSV data')
        endif
      enddo
    enddo

    close(csv_unit)
  end subroutine mt_write_qperp_plane_csv

  subroutine mt_write_qperp_arbitrary_header(csv_file,caller)
    character(len=*), intent(in) :: csv_file,caller

    integer :: csv_unit,io_status

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open CSV file')
    endif

    write(csv_unit,'(a)',iostat=io_status) &
         'i,j,s1,s2,seed_x,seed_y,seed_z,'// &
         'qperp,logqperp,valid,status,'// &
         'N2,bfactor,'// &
         'face_forward,face_backward,'// &
         'status_forward,status_backward,'// &
         'length_forward,length_backward,'// &
         'Bseed_norm,Bf_norm,Bb_norm,'// &
         'x_f,y_f,z_f,x_b,y_b,z_b'
    if (io_status/=0) then
      close(csv_unit)
      call mpistop(trim(caller)//' could not write CSV header')
    endif

    close(csv_unit)
  end subroutine mt_write_qperp_arbitrary_header

  subroutine mt_write_qperp_arbitrary_csv(results,s1,s2,n1,n2,csv_file, &
       caller)
    integer, intent(in) :: n1,n2
    type(trace_qperp_result), intent(in) :: results(n1*n2)
    double precision, intent(in) :: s1(n1*n2),s2(n1*n2)
    character(len=*), intent(in) :: csv_file,caller

    double precision :: Bseed_norm,Bf_norm,Bb_norm
    integer :: csv_unit,io_status,i,j,iseed

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open CSV file')
    endif

    write(csv_unit,'(a)',iostat=io_status) &
         'i,j,s1,s2,seed_x,seed_y,seed_z,'// &
         'qperp,logqperp,valid,status,'// &
         'N2,bfactor,'// &
         'face_forward,face_backward,'// &
         'status_forward,status_backward,'// &
         'length_forward,length_backward,'// &
         'Bseed_norm,Bf_norm,Bb_norm,'// &
         'x_f,y_f,z_f,x_b,y_b,z_b'
    if (io_status/=0) then
      close(csv_unit)
      call mpistop(trim(caller)//' could not write CSV header')
    endif

    do j=1,n2
      do i=1,n1
        iseed=(j-1)*n1+i
        Bseed_norm=dsqrt(sum(results(iseed)%B_seed**2))
        Bf_norm=dsqrt(sum(results(iseed)%forward_B**2))
        Bb_norm=dsqrt(sum(results(iseed)%backward_B**2))
        write(csv_unit, &
             '(i0,",",i0,7(",",es24.16),",",l1,",",i0,'// &
             '2(",",es24.16),4(",",i0),11(",",es24.16))', &
             iostat=io_status) &
             i,j,s1(iseed),s2(iseed),results(iseed)%seed, &
             results(iseed)%qperp,results(iseed)%logqperp, &
             results(iseed)%valid,results(iseed)%status, &
             results(iseed)%N2,results(iseed)%bfactor, &
             results(iseed)%forward_face,results(iseed)%backward_face, &
             results(iseed)%forward_status,results(iseed)%backward_status, &
             results(iseed)%forward_length,results(iseed)%backward_length, &
             Bseed_norm,Bf_norm,Bb_norm, &
             results(iseed)%forward_endpoint, &
             results(iseed)%backward_endpoint
        if (io_status/=0) then
          close(csv_unit)
          call mpistop(trim(caller)//' could not write CSV data')
        endif
      enddo
    enddo

    close(csv_unit)
  end subroutine mt_write_qperp_arbitrary_csv

  subroutine mt_write_mapping_plane_xy_csv(results,nx,ny,csv_file)
    integer, intent(in) :: nx,ny
    type(trace_mapping_result), intent(in) :: results(nx*ny)
    character(len=*), intent(in) :: csv_file

    call mt_write_mapping_plane_csv(results,nx,ny,csv_file, &
         'mt_mapping_plane_xy','ix,iy')
  end subroutine mt_write_mapping_plane_xy_csv

  subroutine mt_write_mapping_plane_csv(results,n1,n2,csv_file,caller, &
       index_header)
    integer, intent(in) :: n1,n2
    type(trace_mapping_result), intent(in) :: results(n1*n2)
    character(len=*), intent(in) :: csv_file,caller,index_header

    integer :: csv_unit,io_status,i,j,iseed

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open CSV file')
    endif

    write(csv_unit,'(a)',iostat=io_status) &
         trim(index_header)//','// &
         'seed_x,seed_y,seed_z,'// &
         'source_Bx,source_By,source_Bz,source_Bn,'// &
         'backward_x,backward_y,backward_z,'// &
         'backward_Bx,backward_By,backward_Bz,backward_Bn,'// &
         'backward_face,backward_length,backward_status,'// &
         'forward_x,forward_y,forward_z,'// &
         'forward_Bx,forward_By,forward_Bz,forward_Bn,'// &
         'forward_face,forward_length,forward_status,valid'
    if (io_status/=0) then
      close(csv_unit)
      call mpistop(trim(caller)//' could not write CSV header')
    endif

    do j=1,n2
      do i=1,n1
        iseed=(j-1)*n1+i
        write(csv_unit, &
             '(i0,",",i0,14(",",es24.16),",",i0,",",es24.16,'// &
             '",",i0,7(",",es24.16),",",i0,",",es24.16,'// &
             '",",i0,",",l1)',iostat=io_status) &
             i,j, &
             results(iseed)%seed, &
             results(iseed)%source_B,results(iseed)%source_Bn, &
             results(iseed)%backward_footpoint, &
             results(iseed)%backward_B,results(iseed)%backward_Bn, &
             results(iseed)%backward_face, &
             results(iseed)%backward_length, &
             results(iseed)%backward_status, &
             results(iseed)%forward_footpoint, &
             results(iseed)%forward_B,results(iseed)%forward_Bn, &
             results(iseed)%forward_face, &
             results(iseed)%forward_length, &
             results(iseed)%forward_status, &
             results(iseed)%valid
        if (io_status/=0) then
          close(csv_unit)
          call mpistop(trim(caller)//' could not write CSV data')
        endif
      enddo
    enddo

    close(csv_unit)
  end subroutine mt_write_mapping_plane_csv

  integer function mt_q_result_face_pair(qperp_result) result(face_pair)
    type(trace_qperp_result), intent(in) :: qperp_result

    face_pair=0
    if (.not.qperp_result%valid_q0) return
    face_pair=mt_q_face_pair_from_faces(qperp_result%backward_face, &
         qperp_result%forward_face)
  end function mt_q_result_face_pair

  integer function mt_q_result_connection_type(qperp_result) &
       result(connection_type)
    type(trace_qperp_result), intent(in) :: qperp_result

    connection_type=0
    if (.not.qperp_result%valid_q0) return
    connection_type=mt_q_connection_type_from_faces( &
         qperp_result%backward_face,qperp_result%forward_face)
  end function mt_q_result_connection_type

  integer function mt_q_face_pair_from_faces(face_b,face_f) result(face_pair)
    integer, intent(in) :: face_b,face_f

    face_pair=0
    if (.not.mt_q_face_valid(face_b)) return
    if (.not.mt_q_face_valid(face_f)) return
    face_pair=10*face_b+face_f
  end function mt_q_face_pair_from_faces

  integer function mt_q_connection_type_from_faces(face_b,face_f) &
       result(connection_type)
    integer, intent(in) :: face_b,face_f

    connection_type=0
    if (.not.mt_q_face_valid(face_b)) return
    if (.not.mt_q_face_valid(face_f)) return
    if (face_b==trace_face_zmin .and. face_f==trace_face_zmin) then
      connection_type=1
    else if (face_b==trace_face_zmax .and. face_f==trace_face_zmax) then
      connection_type=2
    else if ((face_b==trace_face_zmin .and. face_f==trace_face_zmax) .or. &
             (face_b==trace_face_zmax .and. face_f==trace_face_zmin)) then
      connection_type=3
    else if (face_b==face_f) then
      connection_type=4
    else
      connection_type=5
    endif
  end function mt_q_connection_type_from_faces

  logical function mt_q_face_valid(face)
    integer, intent(in) :: face

    mt_q_face_valid=face>=trace_face_xmin .and. face<=trace_face_zmax
  end function mt_q_face_valid

end module mod_magnetic_topology
