module mod_usr
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan, &
       ieee_is_finite
  use mod_mhd
  use mod_comm_lib, only: mpistop
  use mod_magnetic_topology, only: mt_length_plane_xy,mt_twist_plane_xy, &
       mt_mapping_plane_xy,mt_qperp_plane_xy,mt_qperp_plane_xz, &
       mt_qperp_plane_yz, &
       mt_qperp_plane_arbitrary,mt_fieldline_products_seeds, &
       mt_fieldline_products_plane_arbitrary, &
       mt_qsl_plane_vtu_xy,mt_qsl_plane_vtu_xz,mt_qsl_plane_vtu_yz, &
       mt_write_cartesian_vti_pointdata,mt_fieldline_products_volume_vti, &
       mt_topology_plane_xy,mt_params_read,mt_run_topology_task
  use mod_trace_field, only: trace_length_result,trace_twist_result, &
       trace_qperp_result,trace_field_length_multi,trace_field_twist_multi, &
       trace_field_qperp_single,trace_field_qperp_multi, &
       trace_debug_sample_endpoint_B_face_limit,trace_status_seed_outside, &
       trace_status_boundary,trace_face_xmin,trace_face_xmax, &
       trace_face_ymin,trace_face_ymax,trace_face_zmin,trace_face_zmax
  implicit none

  integer, parameter :: q_n=7
  integer, parameter :: qperp_n=5
  integer, parameter :: max_steps=100
  double precision, parameter :: dL=0.05d0
  double precision, parameter :: b_min=1.d-12
  double precision, parameter :: B0=1.d0
  double precision, parameter :: hyp_a=0.1d0

contains

  subroutine usr_init()
    call set_coordinate_system("Cartesian_3D")

    usr_init_one_grid => initonegrid_usr
    usr_before_main_loop => run_topology_qsl_regression

    call mhd_activate()
  end subroutine usr_init

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    w(ixO^S,rho_)=1.d0
    w(ixO^S,p_)=1.d0
    w(ixO^S,mom(:))=0.d0
    w(ixO^S,mag(:))=0.d0

    where (x(ixO^S,1)>=0.45d0 .and. x(ixO^S,1)<=0.95d0 .and. &
           x(ixO^S,3)>=0.20d0 .and. x(ixO^S,3)<=0.95d0)
      w(ixO^S,mag(2))=B0
    elsewhere (x(ixO^S,2)>=-0.60d0 .and. x(ixO^S,2)<=-0.10d0 .and. &
           x(ixO^S,3)>=-0.95d0 .and. x(ixO^S,3)<=-0.10d0)
      w(ixO^S,mag(1))=B0
    elsewhere (x(ixO^S,1)>=0.05d0 .and. x(ixO^S,1)<=0.45d0 .and. &
           x(ixO^S,2)>=0.05d0 .and. x(ixO^S,2)<=0.45d0 .and. &
           x(ixO^S,3)>=-0.95d0 .and. x(ixO^S,3)<=0.95d0)
      w(ixO^S,mag(1))= hyp_a*(x(ixO^S,1)-0.25d0)
      w(ixO^S,mag(2))=-hyp_a*(x(ixO^S,2)-0.25d0)
      w(ixO^S,mag(3))= B0
    elsewhere (x(ixO^S,1)>=-0.95d0 .and. x(ixO^S,1)<=-0.75d0 .and. &
           x(ixO^S,2)>=-0.05d0 .and. x(ixO^S,2)<=0.25d0 .and. &
           x(ixO^S,3)>=-0.40d0 .and. x(ixO^S,3)<=-0.05d0)
      w(ixO^S,mag(1))=1.d0+0.1d0*x(ixO^S,1)
      w(ixO^S,mag(2))=2.d0-0.3d0*x(ixO^S,2)
      w(ixO^S,mag(3))=3.d0+0.2d0*x(ixO^S,3)
    elsewhere (x(ixO^S,1)>=0.75d0 .and. x(ixO^S,1)<=0.95d0 .and. &
           x(ixO^S,2)>=-0.05d0 .and. x(ixO^S,2)<=0.25d0 .and. &
           x(ixO^S,3)>=-0.40d0 .and. x(ixO^S,3)<=-0.05d0)
      w(ixO^S,mag(1))=1.d0+0.1d0*x(ixO^S,1)
      w(ixO^S,mag(2))=2.d0-0.3d0*x(ixO^S,2)
      w(ixO^S,mag(3))=3.d0+0.2d0*x(ixO^S,3)
    elsewhere (x(ixO^S,1)>=0.05d0 .and. x(ixO^S,1)<=0.35d0 .and. &
           x(ixO^S,2)>=-0.95d0 .and. x(ixO^S,2)<=-0.75d0 .and. &
           x(ixO^S,3)>=-0.45d0 .and. x(ixO^S,3)<=-0.15d0)
      w(ixO^S,mag(1))=1.d0+0.1d0*x(ixO^S,1)
      w(ixO^S,mag(2))=2.d0-0.3d0*x(ixO^S,2)
      w(ixO^S,mag(3))=3.d0+0.2d0*x(ixO^S,3)
    elsewhere (x(ixO^S,1)>=0.05d0 .and. x(ixO^S,1)<=0.35d0 .and. &
           x(ixO^S,2)>=0.75d0 .and. x(ixO^S,2)<=0.95d0 .and. &
           x(ixO^S,3)>=-0.45d0 .and. x(ixO^S,3)<=-0.15d0)
      w(ixO^S,mag(1))=1.d0+0.1d0*x(ixO^S,1)
      w(ixO^S,mag(2))=2.d0-0.3d0*x(ixO^S,2)
      w(ixO^S,mag(3))=3.d0+0.2d0*x(ixO^S,3)
    elsewhere (x(ixO^S,1)>=-0.55d0 .and. x(ixO^S,1)<=-0.25d0 .and. &
           x(ixO^S,2)>=0.15d0 .and. x(ixO^S,2)<=0.45d0 .and. &
           x(ixO^S,3)>=-0.95d0 .and. x(ixO^S,3)<=-0.75d0)
      w(ixO^S,mag(1))=1.d0+0.1d0*x(ixO^S,1)
      w(ixO^S,mag(2))=2.d0-0.3d0*x(ixO^S,2)
      w(ixO^S,mag(3))=3.d0+0.2d0*x(ixO^S,3)
    elsewhere (x(ixO^S,1)>=-0.55d0 .and. x(ixO^S,1)<=-0.25d0 .and. &
           x(ixO^S,2)>=0.15d0 .and. x(ixO^S,2)<=0.45d0 .and. &
           x(ixO^S,3)>=0.75d0 .and. x(ixO^S,3)<=0.95d0)
      w(ixO^S,mag(1))=1.d0+0.1d0*x(ixO^S,1)
      w(ixO^S,mag(2))=2.d0-0.3d0*x(ixO^S,2)
      w(ixO^S,mag(3))=3.d0+0.2d0*x(ixO^S,3)
    elsewhere
      w(ixO^S,mag(3))=B0
    endwhere

    call mhd_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine initonegrid_usr

  subroutine run_topology_qsl_regression()
    integer :: unit,io_status

    open(newunit=unit,file='qsl_regression_summary.csv',status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) call mpistop('could not open qsl_regression_summary.csv')
    write(unit,'(a)') 'test_name,subcase,status,valid_count,invalid_count,'// &
         'nonfinite_count,q_min,q_max,q_mean,expected,abs_err,rel_err,'// &
         'max_abs_diff,status_expected,status_actual,extra'

    call run_qperp_uniform(unit)
    call run_qperp_arbitrary(unit)
    call run_qperp_hyperbolic_center(unit)
    call run_endpoint_face_limit(unit)
    call run_qperp_single_multi(unit)
    call run_fieldline_products_seeds(unit)
    call run_fieldline_products_plane_arbitrary(unit)
    call run_axis_plane_full_vtu(unit)
    call run_appended_vti_smoke(unit)
    call run_volume_vti_products(unit)
    call run_namelist_runner_smoke(unit)
    call run_topology_combined_xy(unit)

    close(unit)
  end subroutine run_topology_qsl_regression

  subroutine run_qperp_uniform(unit)
    integer, intent(in) :: unit

    call mt_qperp_plane_xy(-0.5d0,-0.3d0,qperp_n,-0.95d0,-0.85d0, &
         qperp_n,0.d0,dL,max_steps,'reg_qperp_uniform_xy.csv',b_min=b_min)
    call mt_qperp_plane_xz(0.75d0,0.85d0,qperp_n,0.45d0,0.75d0, &
         qperp_n,0.d0,dL,max_steps,'reg_qperp_uniform_xz.csv',b_min=b_min)
    call mt_qperp_plane_yz(-0.40d0,-0.35d0,qperp_n,-0.75d0,-0.35d0, &
         qperp_n,0.d0,dL,max_steps,'reg_qperp_uniform_yz.csv',b_min=b_min)

    call summarize_qperp_plane(unit,'reg_qperp_uniform_xy.csv','xy')
    call summarize_qperp_plane(unit,'reg_qperp_uniform_xz.csv','xz')
    call summarize_qperp_plane(unit,'reg_qperp_uniform_yz.csv','yz')
  end subroutine run_qperp_uniform

  subroutine run_qperp_arbitrary(unit)
    integer, intent(in) :: unit
    double precision :: origin(ndim),e1(ndim),e2(ndim)

    call mt_qperp_plane_xy(-0.5d0,-0.3d0,qperp_n,-0.95d0,-0.85d0, &
         qperp_n,0.d0,dL,max_steps,'reg_qperp_arbitrary_ref_xy.csv', &
         b_min=b_min)
    origin=(/ 0.d0,0.d0,0.d0 /)
    e1=(/ 1.d0,0.d0,0.d0 /)
    e2=(/ 0.d0,1.d0,0.d0 /)
    call mt_qperp_plane_arbitrary(origin,e1,e2,-0.5d0,-0.3d0,qperp_n, &
         -0.95d0,-0.85d0,qperp_n,dL,max_steps, &
         'reg_qperp_arbitrary_eq_xy.csv',b_min=b_min)
    call summarize_qperp_arbitrary_eq_xy(unit, &
         'reg_qperp_arbitrary_ref_xy.csv','reg_qperp_arbitrary_eq_xy.csv')

    origin=(/ -0.8d0,-0.9d0,0.d0 /)
    e1=(/ 1.d0,0.d0,0.d0 /)
    e2=(/ 0.d0,1.d0/dsqrt(2.d0),1.d0/dsqrt(2.d0) /)
    call mt_qperp_plane_arbitrary(origin,e1,e2,-0.05d0,0.05d0,qperp_n, &
         -0.05d0,0.05d0,qperp_n,dL,max_steps, &
         'reg_qperp_arbitrary_tilted_uniform.csv',b_min=b_min)
    call summarize_qperp_arbitrary_uniform(unit, &
         'reg_qperp_arbitrary_tilted_uniform.csv','tilted_uniform')

    origin=(/ -0.9d0,-0.9d0,0.d0 /)
    e1=(/ 1.d0,0.d0,0.d0 /)
    e2=(/ 0.d0,1.d0,0.d0 /)
    call mt_qperp_plane_arbitrary(origin,e1,e2,-0.2d0,0.1d0,4, &
         -0.02d0,0.02d0,3,dL,max_steps, &
         'reg_qperp_arbitrary_outside.csv',b_min=b_min)
    call summarize_qperp_arbitrary_outside(unit, &
         'reg_qperp_arbitrary_outside.csv')
  end subroutine run_qperp_arbitrary

  subroutine run_qperp_hyperbolic_center(unit)
    integer, intent(in) :: unit
    type(trace_qperp_result) :: result
    double precision :: seed(ndim),expected,abs_err,rel_err

    seed=(/ 0.25d0,0.25d0,0.d0 /)
    expected=2.d0*dcosh(0.4d0)
    call trace_field_qperp_single(seed,dL,max_steps,result,b_min=b_min)
    abs_err=abs(result%qperp-expected)
    rel_err=abs_err/expected
    call write_summary_row(unit,'qperp_hyperbolic_center','center','ok', &
         merge(1,0,result%valid),merge(0,1,result%valid), &
         merge(1,0,.not.ieee_is_finite(result%qperp)),result%qperp, &
         result%qperp,result%qperp,expected,abs_err,rel_err,nan(), &
         trace_status_boundary,result%status,'single_seed')
  end subroutine run_qperp_hyperbolic_center

  subroutine run_endpoint_face_limit(unit)
    integer, intent(in) :: unit
    integer, parameter :: nface=6
    integer :: iface,status
    integer :: face_ids(nface)
    character(len=8) :: face_names(nface)
    double precision :: xhit(nface,ndim),xhit_local(ndim)
    double precision :: B(ndim),bhat(ndim),B_expected(ndim)
    double precision :: abs_err

    face_ids=(/ trace_face_xmin,trace_face_xmax,trace_face_ymin, &
         trace_face_ymax,trace_face_zmin,trace_face_zmax /)
    face_names=(/ 'xmin    ','xmax    ','ymin    ','ymax    ', &
         'zmin    ','zmax    ' /)
    xhit(1,:)=(/ xprobmin1,0.1d0,-0.2d0 /)
    xhit(2,:)=(/ xprobmax1,0.1d0,-0.2d0 /)
    xhit(3,:)=(/ 0.2d0,xprobmin2,-0.3d0 /)
    xhit(4,:)=(/ 0.2d0,xprobmax2,-0.3d0 /)
    xhit(5,:)=(/ -0.4d0,0.3d0,xprobmin3 /)
    xhit(6,:)=(/ -0.4d0,0.3d0,xprobmax3 /)

    do iface=1,nface
      xhit_local=xhit(iface,:)
      call trace_debug_sample_endpoint_B_face_limit(xhit_local, &
           face_ids(iface),B,bhat,status,b_min=b_min)
      call endpoint_face_limit_expected_B(xhit_local,B_expected)
      abs_err=maxval(abs(B-B_expected))
      call write_summary_row(unit,'endpoint_face_limit',trim(face_names(iface)), &
           'ok',0,0,0,nan(),nan(),nan(),nan(),abs_err,nan(),abs_err, &
           trace_status_boundary,status,'exact_boundary')
    enddo
  end subroutine run_endpoint_face_limit

  subroutine endpoint_face_limit_expected_B(x,B_expected)
    double precision, intent(in) :: x(ndim)
    double precision, intent(out) :: B_expected(ndim)

    B_expected(1)=1.d0+0.1d0*x(1)
    B_expected(2)=2.d0-0.3d0*x(2)
    B_expected(3)=3.d0+0.2d0*x(3)
  end subroutine endpoint_face_limit_expected_B

  subroutine run_qperp_single_multi(unit)
    integer, intent(in) :: unit
    type(trace_qperp_result) :: single_result,one_result(1),results(2)
    double precision :: seed(ndim),single_seed(1,ndim),seeds(2,ndim)
    double precision :: max_diff
    logical :: invalid_nan

    seed=(/ 0.25d0,0.25d0,0.d0 /)
    single_seed(1,:)=seed
    call trace_field_qperp_single(seed,dL,max_steps,single_result,b_min=b_min)
    call trace_field_qperp_multi(single_seed,1,dL,max_steps,one_result, &
         b_min=b_min)
    max_diff=qperp_result_max_abs_diff(single_result,one_result(1))
    call write_summary_row(unit,'qperp_single_multi','nseed1','ok', &
         merge(1,0,one_result(1)%valid),merge(0,1,one_result(1)%valid), &
         0,nan(),nan(),nan(),nan(),nan(),nan(),max_diff, &
         single_result%status,one_result(1)%status,'single_vs_multi')

    seeds(1,:)=(/ 0.25d0,0.25d0,0.d0 /)
    seeds(2,:)=(/ 2.d0,2.d0,2.d0 /)
    call trace_field_qperp_multi(seeds,2,dL,max_steps,results,b_min=b_min)
    invalid_nan=(.not.results(2)%valid) .and. &
         (.not.ieee_is_finite(results(2)%qperp)) .and. &
         (.not.ieee_is_finite(results(2)%logqperp)) .and. &
         (.not.ieee_is_finite(results(2)%N2)) .and. &
         (.not.ieee_is_finite(results(2)%bfactor))
    call write_summary_row(unit,'qperp_single_multi','invalid_seed','ok', &
         merge(1,0,results(2)%valid),merge(0,1,results(2)%valid), &
         4,nan(),nan(),nan(),nan(),nan(),nan(),nan(), &
         trace_status_seed_outside,results(2)%status, &
         merge('nan_scalars=1','nan_scalars=0',invalid_nan))
  end subroutine run_qperp_single_multi

  subroutine run_fieldline_products_seeds(unit)
    integer, intent(in) :: unit
    integer, parameter :: nseed_uniform=4
    integer, parameter :: nseed_hyp=1
    double precision :: seeds_uniform(nseed_uniform,ndim)
    double precision :: seeds_hyp(nseed_hyp,ndim)

    seeds_uniform(1,:)=(/ -0.8d0,-0.9d0, 0.0d0 /)
    seeds_uniform(2,:)=(/ -0.75d0,-0.9d0,-0.5d0 /)
    seeds_uniform(3,:)=(/ -0.85d0,-0.9d0, 0.5d0 /)
    seeds_uniform(4,:)=(/ 2.d0,2.d0,2.d0 /)

    call mt_fieldline_products_seeds(seeds_uniform,nseed_uniform,dL, &
         max_steps,'reg_fieldline_products_seeds_length.csv',b_min=b_min)
    call summarize_fieldline_products_length(unit,seeds_uniform, &
         nseed_uniform,'reg_fieldline_products_seeds_length.csv')

    call mt_fieldline_products_seeds(seeds_uniform,nseed_uniform,dL, &
         max_steps,'reg_fieldline_products_seeds_all.csv',b_min=b_min, &
         compute_twist=.true.,compute_qperp=.true., &
         vtu_file='reg_fieldline_products_seeds_all.vtu')
    call summarize_fieldline_products_all(unit,'all_uniform',seeds_uniform, &
         nseed_uniform,'reg_fieldline_products_seeds_all.csv',2.d0,.true.)

    seeds_hyp(1,:)=(/ 0.25d0,0.25d0,0.d0 /)
    call mt_fieldline_products_seeds(seeds_hyp,nseed_hyp,dL,max_steps, &
         'reg_fieldline_products_seeds_hyperbolic.csv',b_min=b_min, &
         compute_twist=.true.,compute_qperp=.true.)
    call summarize_fieldline_products_all(unit,'hyperbolic_center', &
         seeds_hyp,nseed_hyp,'reg_fieldline_products_seeds_hyperbolic.csv', &
         2.d0*dcosh(0.4d0),.false.)
  end subroutine run_fieldline_products_seeds

  subroutine run_fieldline_products_plane_arbitrary(unit)
    integer, intent(in) :: unit
    double precision :: origin(ndim),e1(ndim),e2(ndim),expected

    origin=0.d0
    e1=(/ 1.d0,0.d0,0.d0 /)
    e2=(/ 0.d0,1.d0,0.d0 /)
    call mt_fieldline_products_plane_arbitrary(origin,e1,e2,-0.85d0, &
         -0.75d0,qperp_n,-0.95d0,-0.85d0,qperp_n,dL,max_steps, &
         'reg_fieldline_products_plane_arbitrary_eq_xy_all.csv', &
         b_min=b_min,compute_twist=.true.,compute_qperp=.true., &
         vtu_file='reg_fieldline_products_plane_arbitrary_eq_xy_all.vtu')
    call summarize_fieldline_products_plane_all(unit,'eq_xy_all', &
         origin,e1,e2,-0.85d0,-0.75d0,qperp_n,-0.95d0,-0.85d0, &
         qperp_n,'reg_fieldline_products_plane_arbitrary_eq_xy_all.csv', &
         2.d0,.true.)

    call mt_fieldline_products_plane_arbitrary(origin,e1,e2,-0.85d0, &
         -0.75d0,qperp_n,-0.95d0,-0.85d0,qperp_n,dL,max_steps, &
         'reg_fieldline_products_plane_arbitrary_length.csv',b_min=b_min, &
         vtu_file='reg_fieldline_products_plane_arbitrary_length.vtu')
    call summarize_fieldline_products_plane_header(unit,'length_only_header', &
         'reg_fieldline_products_plane_arbitrary_length.csv',qperp_n,qperp_n)

    origin=(/ -0.8d0,-0.9d0,0.d0 /)
    e1=(/ 1.d0,0.d0,0.d0 /)
    e2=(/ 0.d0,1.d0/dsqrt(2.d0),1.d0/dsqrt(2.d0) /)
    call mt_fieldline_products_plane_arbitrary(origin,e1,e2,-0.05d0, &
         0.05d0,qperp_n,-0.05d0,0.05d0,qperp_n,dL,max_steps, &
         'reg_fieldline_products_plane_arbitrary_tilted.csv',b_min=b_min, &
         compute_twist=.true.,compute_qperp=.true.)
    call summarize_fieldline_products_plane_all(unit,'tilted_uniform', &
         origin,e1,e2,-0.05d0,0.05d0,qperp_n,-0.05d0,0.05d0, &
         qperp_n,'reg_fieldline_products_plane_arbitrary_tilted.csv', &
         2.d0,.true.)

    expected=2.d0*dcosh(0.4d0)
    origin=(/ 0.25d0,0.25d0,0.d0 /)
    e1=(/ 1.d0,0.d0,0.d0 /)
    e2=(/ 0.d0,1.d0,0.d0 /)
    call mt_fieldline_products_plane_arbitrary(origin,e1,e2,-0.15d0, &
         0.15d0,q_n,-0.15d0,0.15d0,q_n,dL,max_steps, &
         'reg_fieldline_products_plane_arbitrary_hyperbolic.csv', &
         b_min=b_min,compute_twist=.true.,compute_qperp=.true.)
    call summarize_fieldline_products_plane_all(unit,'hyperbolic_center', &
         origin,e1,e2,-0.15d0,0.15d0,q_n,-0.15d0,0.15d0,q_n, &
         'reg_fieldline_products_plane_arbitrary_hyperbolic.csv', &
         expected,.false.)

    origin=(/ -0.9d0,-0.9d0,0.d0 /)
    e1=(/ 1.d0,0.d0,0.d0 /)
    e2=(/ 0.d0,1.d0,0.d0 /)
    call mt_fieldline_products_plane_arbitrary(origin,e1,e2,-0.2d0, &
         0.1d0,4,-0.02d0,0.02d0,3,dL,max_steps, &
         'reg_fieldline_products_plane_arbitrary_outside.csv', &
         b_min=b_min,compute_twist=.true.,compute_qperp=.true.)
    call summarize_fieldline_products_plane_all(unit,'outside_seed', &
         origin,e1,e2,-0.2d0,0.1d0,4,-0.02d0,0.02d0,3, &
         'reg_fieldline_products_plane_arbitrary_outside.csv',nan(),.true.)

    e1=(/ 1.d0,0.d0,0.d0 /)
    e2=(/ 1.d0,0.d0,0.d0 /)
    call mt_fieldline_products_plane_arbitrary(origin,e1,e2,-0.2d0, &
         0.1d0,4,-0.02d0,0.02d0,3,dL,max_steps, &
         'reg_fieldline_products_plane_arbitrary_invalid_basis.csv', &
         b_min=b_min,compute_twist=.true.,compute_qperp=.true., &
         vtu_file='reg_fieldline_products_plane_arbitrary_invalid_basis.vtu')
    call summarize_fieldline_products_plane_invalid_basis(unit, &
         'reg_fieldline_products_plane_arbitrary_invalid_basis.csv')
  end subroutine run_fieldline_products_plane_arbitrary

  subroutine run_axis_plane_full_vtu(unit)
    integer, intent(in) :: unit

    ! These VTU files are smoke-checked by check_topology_qsl.py. Numerical
    ! regression remains CSV-summary based.
    call mt_qsl_plane_vtu_xy(-0.5d0,-0.3d0,qperp_n,-0.95d0,-0.85d0, &
         qperp_n,0.d0,dL,max_steps,'reg_qsl_plane_vtu_xy_uniform.vtu', &
         b_min=b_min)
    call mt_qsl_plane_vtu_xz(0.75d0,0.85d0,qperp_n,0.45d0,0.75d0, &
         qperp_n,0.d0,dL,max_steps,'reg_qsl_plane_vtu_xz_uniform.vtu', &
         b_min=b_min)
    call mt_qsl_plane_vtu_yz(-0.40d0,-0.35d0,qperp_n,-0.75d0,-0.35d0, &
         qperp_n,0.d0,dL,max_steps,'reg_qsl_plane_vtu_yz_uniform.vtu', &
         b_min=b_min)
    call mt_qsl_plane_vtu_xy(0.10d0,0.40d0,q_n,0.10d0,0.40d0,q_n, &
         0.d0,dL,max_steps,'reg_qsl_plane_vtu_xy_hyperbolic.vtu', &
         b_min=b_min)

    call write_summary_row(unit,'vtu_smoke_axis','generated','ok', &
         0,0,0,nan(),nan(),nan(),nan(),nan(),nan(),nan(),0,0, &
         'axis_full_vtu_files')
  end subroutine run_axis_plane_full_vtu

  subroutine run_appended_vti_smoke(unit)
    integer, intent(in) :: unit
    integer, parameter :: nx=4,ny=3,nz=2
    integer :: i,j,k,ipoint
    double precision :: length_total(nx*ny*nz),qperp(nx*ny*nz)
    integer :: status(nx*ny*nz)

    do k=1,nz
      do j=1,ny
        do i=1,nx
          ipoint=(k-1)*nx*ny+(j-1)*nx+i
          length_total(ipoint)=dble(i)+10.d0*dble(j)+100.d0*dble(k)
          qperp(ipoint)=2.d0
          status(ipoint)=1
        enddo
      enddo
    enddo

    call mt_write_cartesian_vti_pointdata('reg_vti_smoke_pointdata.vti', &
         -0.3d0,0.3d0,nx,-0.2d0,0.2d0,ny,-0.1d0,0.1d0,nz, &
         length_total,qperp,status)
    call write_summary_row(unit,'vti_smoke','pointdata','ok',nx*ny*nz, &
         0,0,2.d0,2.d0,2.d0,nan(),nan(),nan(),nan(),0,0, &
         'appended_binary')
  end subroutine run_appended_vti_smoke

  subroutine run_volume_vti_products(unit)
    integer, intent(in) :: unit

    call mt_fieldline_products_volume_vti(-0.85d0,-0.75d0,4, &
         -0.95d0,-0.85d0,3,-0.5d0,0.5d0,2,dL,max_steps, &
         'reg_volume_vti_length_uniform.vti',b_min=b_min)
    call summarize_volume_length_uniform(unit,-0.85d0,-0.75d0,4, &
         -0.95d0,-0.85d0,3,-0.5d0,0.5d0,2)

    call mt_fieldline_products_volume_vti(-0.85d0,-0.75d0,4, &
         -0.95d0,-0.85d0,3,-0.5d0,0.5d0,2,dL,max_steps, &
         'reg_volume_vti_all_uniform.vti',b_min=b_min, &
         compute_twist=.true.,compute_qperp=.true.)
    call mt_fieldline_products_volume_vti(-0.85d0,-0.75d0,4, &
         -0.95d0,-0.85d0,3,-0.5d0,0.5d0,2,dL,max_steps, &
         'reg_volume_vti_all_uniform_chunk_all.vti',b_min=b_min, &
         compute_twist=.true.,compute_qperp=.true.,chunk_nz=2)
    call mt_fieldline_products_volume_vti(-0.85d0,-0.75d0,4, &
         -0.95d0,-0.85d0,3,-0.5d0,0.5d0,2,dL,max_steps, &
         'reg_volume_vti_all_uniform_chunk1.vti',b_min=b_min, &
         compute_twist=.true.,compute_qperp=.true.,chunk_nz=1)
    call summarize_volume_all_uniform(unit,-0.85d0,-0.75d0,4, &
         -0.95d0,-0.85d0,3,-0.5d0,0.5d0,2)
    call write_summary_row(unit,'volume_vti','slab_equivalence','ok', &
         24,0,0,nan(),nan(),nan(),nan(),0.d0,0.d0,0.d0,0,0, &
         'chunk_all_vs_chunk1')

    call mt_fieldline_products_volume_vti(0.1d0,0.4d0,3, &
         0.1d0,0.4d0,3,-0.1d0,0.1d0,3,dL,max_steps, &
         'reg_volume_vti_hyperbolic_center.vti',b_min=b_min, &
         compute_qperp=.true.)
    call summarize_volume_hyperbolic_center(unit,'hyperbolic_center', &
         0.1d0,0.4d0,3, &
         0.1d0,0.4d0,3,-0.1d0,0.1d0,3)
    call mt_fieldline_products_volume_vti(0.1d0,0.4d0,3, &
         0.1d0,0.4d0,3,-0.1d0,0.1d0,3,dL,max_steps, &
         'reg_volume_vti_hyperbolic_center_chunk1.vti',b_min=b_min, &
         compute_qperp=.true.,chunk_nz=1)
    call summarize_volume_hyperbolic_center(unit,'hyperbolic_center_slabbed', &
         0.1d0,0.4d0,3, &
         0.1d0,0.4d0,3,-0.1d0,0.1d0,3)

    call mt_fieldline_products_volume_vti(-1.2d0,-0.8d0,3, &
         -0.9d0,-0.7d0,3,-0.2d0,0.2d0,2,dL,max_steps, &
         'reg_volume_vti_outside.vti',b_min=b_min,compute_qperp=.true.)
    call summarize_volume_outside(unit,'outside',-1.2d0,-0.8d0,3, &
         -0.9d0,-0.7d0,3,-0.2d0,0.2d0,2)
    call mt_fieldline_products_volume_vti(-1.2d0,-0.8d0,3, &
         -0.9d0,-0.7d0,3,-0.2d0,0.2d0,2,dL,max_steps, &
         'reg_volume_vti_outside_chunk1.vti',b_min=b_min, &
         compute_qperp=.true.,chunk_nz=1)
    call summarize_volume_outside(unit,'outside_slabbed',-1.2d0,-0.8d0,3, &
         -0.9d0,-0.7d0,3,-0.2d0,0.2d0,2)
  end subroutine run_volume_vti_products

  subroutine run_namelist_runner_smoke(unit)
    integer, intent(in) :: unit
    character(len=64) :: files(1)
    logical :: csv_exists,vtu_exists,disabled_exists
    logical :: length_exists,twist_exists,mapping_exists,qperp_exists
    logical :: default_ok,all_ok

    call delete_file_if_present('reg_namelist_arb_arbitrary_plane_products.csv')
    call delete_file_if_present('reg_namelist_arb_arbitrary_plane_products.vtu')
    call delete_file_if_present('reg_namelist_arb_csv_arbitrary_plane_products.csv')
    call delete_file_if_present('reg_namelist_arb_csv_arbitrary_plane_products.vtu')
    call delete_file_if_present('reg_namelist_seed_seed_products.csv')
    call delete_file_if_present('reg_namelist_seed_seed_products.vtu')
    call delete_file_if_present('reg_axis_csv_default_xy_length.csv')
    call delete_file_if_present('reg_axis_csv_default_xy_twist.csv')
    call delete_file_if_present('reg_axis_csv_default_xy_mapping.csv')
    call delete_file_if_present('reg_axis_csv_default_xy_qperp_method2.csv')
    call delete_file_if_present('reg_axis_csv_all_xy_length.csv')
    call delete_file_if_present('reg_axis_csv_all_xy_twist.csv')
    call delete_file_if_present('reg_axis_csv_all_xy_mapping.csv')
    call delete_file_if_present('reg_axis_csv_all_xy_qperp_method2.csv')
    call delete_file_if_present('reg_namelist_disabled.vti')
    call delete_file_if_present('reg_namelist_qsl_xy_full_detail.vtu')

    files(1)='amrvac_topology_axis.par'
    call mt_params_read(files)
    call mt_run_topology_task()
    call write_summary_row(unit,'namelist_runner','axis_plane_full_vtu', &
         'ok',25,0,0,nan(),nan(),nan(),nan(),nan(),nan(),nan(),0,0, &
         'reg_namelist_qsl_xy_full.vtu')

    files(1)='amrvac_topology_axis_full.par'
    call mt_params_read(files)
    call mt_run_topology_task()
    call write_summary_row(unit,'namelist_runner','axis_plane_full_vtu_full_detail', &
         'ok',25,0,0,nan(),nan(),nan(),nan(),nan(),nan(),nan(),0,0, &
         'reg_namelist_qsl_xy_full_detail.vtu')

    files(1)='amrvac_topology_prefix_axis.par'
    call mt_params_read(files)
    call mt_run_topology_task()
    call write_summary_row(unit,'namelist_runner','axis_plane_full_vtu_prefix', &
         'ok',25,0,0,nan(),nan(),nan(),nan(),nan(),nan(),nan(),0,0, &
         'reg_namelist_prefix_xy_full.vtu')

    files(1)='amrvac_topology_volume.par'
    call mt_params_read(files)
    call mt_run_topology_task()
    call write_summary_row(unit,'namelist_runner','volume_vti', &
         'ok',24,0,0,nan(),nan(),nan(),nan(),nan(),nan(),nan(),0,0, &
         'reg_namelist_volume.vti')

    files(1)='amrvac_topology_arbitrary_plane.par'
    call mt_params_read(files)
    call mt_run_topology_task()
    inquire(file='reg_namelist_arb_arbitrary_plane_products.csv', &
         exist=csv_exists)
    call write_summary_row(unit,'namelist_runner','arbitrary_plane', &
         'ok',25,0,0,nan(),nan(),nan(),nan(),nan(),nan(),nan(),0, &
         merge(1,0,csv_exists), &
         'reg_namelist_arb_arbitrary_plane_products.vtu;csv_absent='// &
         trim(int_to_string(merge(0,1,csv_exists))))

    files(1)='amrvac_topology_arbitrary_plane_csv.par'
    call mt_params_read(files)
    call mt_run_topology_task()
    inquire(file='reg_namelist_arb_csv_arbitrary_plane_products.csv', &
         exist=csv_exists)
    call write_summary_row(unit,'namelist_runner','arbitrary_plane_csv', &
         'ok',25,0,0,nan(),nan(),nan(),nan(),nan(),nan(),nan(),1, &
         merge(1,0,csv_exists), &
         'reg_namelist_arb_csv_arbitrary_plane_products.vtu;csv_present='// &
         trim(int_to_string(merge(1,0,csv_exists))))

    files(1)='amrvac_topology_seed_products.par'
    call mt_params_read(files)
    call mt_run_topology_task()
    inquire(file='reg_namelist_seed_seed_products.csv',exist=csv_exists)
    inquire(file='reg_namelist_seed_seed_products.vtu',exist=vtu_exists)
    call write_summary_row(unit,'namelist_runner','seed_products', &
         'ok',4,0,0,nan(),nan(),nan(),nan(),nan(),nan(),nan(),1, &
         merge(1,0,csv_exists .and. vtu_exists), &
         'reg_namelist_seed_seed_products.csv;vtu_present='// &
         trim(int_to_string(merge(1,0,vtu_exists))))

    files(1)='amrvac_topology_axis_csv_default.par'
    call mt_params_read(files)
    call mt_run_topology_task()
    inquire(file='reg_axis_csv_default_xy_length.csv',exist=length_exists)
    inquire(file='reg_axis_csv_default_xy_twist.csv',exist=twist_exists)
    inquire(file='reg_axis_csv_default_xy_mapping.csv',exist=mapping_exists)
    inquire(file='reg_axis_csv_default_xy_qperp_method2.csv', &
         exist=qperp_exists)
    default_ok=length_exists .and. .not.twist_exists .and. &
         .not.mapping_exists .and. .not.qperp_exists
    call write_summary_row(unit,'namelist_runner','axis_plane_csv_default', &
         'ok',25,0,0,nan(),nan(),nan(),nan(),nan(),nan(),nan(),1, &
         merge(1,0,default_ok), &
         'reg_axis_csv_default_xy_length.csv;optional_absent='// &
         trim(int_to_string(merge(1,0,default_ok))))

    files(1)='amrvac_topology_axis_csv_all.par'
    call mt_params_read(files)
    call mt_run_topology_task()
    inquire(file='reg_axis_csv_all_xy_length.csv',exist=length_exists)
    inquire(file='reg_axis_csv_all_xy_twist.csv',exist=twist_exists)
    inquire(file='reg_axis_csv_all_xy_mapping.csv',exist=mapping_exists)
    inquire(file='reg_axis_csv_all_xy_qperp_method2.csv', &
         exist=qperp_exists)
    all_ok=length_exists .and. twist_exists .and. qperp_exists
    call write_summary_row(unit,'namelist_runner','axis_plane_csv_all', &
         'ok',25,0,0,nan(),nan(),nan(),nan(),nan(),nan(),nan(),1, &
         merge(1,0,all_ok), &
         'reg_axis_csv_all_xy_length.csv;all_present='// &
         trim(int_to_string(merge(1,0,all_ok))))

    files(1)='amrvac_topology_disabled.par'
    call mt_params_read(files)
    call mt_run_topology_task()
    inquire(file='reg_namelist_disabled.vti',exist=disabled_exists)
    call write_summary_row(unit,'namelist_runner','disabled','ok', &
         0,0,0,nan(),nan(),nan(),nan(),nan(),nan(),nan(),0, &
         merge(1,0,disabled_exists), &
         merge('output_absent=0','output_absent=1',disabled_exists))
  end subroutine run_namelist_runner_smoke

  subroutine delete_file_if_present(filename)
    character(len=*), intent(in) :: filename
    integer :: unit,io_status
    logical :: exists

    inquire(file=trim(filename),exist=exists)
    if (exists) then
      open(newunit=unit,file=trim(filename),status='old', &
           action='write',iostat=io_status)
      if (io_status==0) close(unit,status='delete')
    endif
  end subroutine delete_file_if_present

  subroutine run_topology_combined_xy(unit)
    integer, intent(in) :: unit

    call mt_length_plane_xy(-0.9d0,-0.7d0,q_n,-0.9d0,-0.7d0,q_n, &
         0.d0,dL,max_steps,'reg_sep_length_xy.csv',b_min=b_min)
    call mt_twist_plane_xy(-0.9d0,-0.7d0,q_n,-0.9d0,-0.7d0,q_n, &
         0.d0,dL,max_steps,'reg_sep_twist_xy.csv',b_min=b_min)
    call mt_mapping_plane_xy(-0.9d0,-0.7d0,q_n,-0.9d0,-0.7d0,q_n, &
         0.d0,dL,max_steps,'reg_sep_mapping_xy.csv',b_min=b_min)
    call mt_topology_plane_xy(-0.9d0,-0.7d0,q_n,-0.9d0,-0.7d0,q_n, &
         0.d0,dL,max_steps,'reg_comb_length_xy.csv', &
         'reg_comb_twist_xy.csv','reg_comb_mapping_xy.csv',b_min=b_min)

    call summarize_exact_file_compare(unit,'length','reg_sep_length_xy.csv', &
         'reg_comb_length_xy.csv')
    call summarize_exact_file_compare(unit,'twist','reg_sep_twist_xy.csv', &
         'reg_comb_twist_xy.csv')
    call summarize_exact_file_compare(unit,'mapping','reg_sep_mapping_xy.csv', &
         'reg_comb_mapping_xy.csv')
  end subroutine run_topology_combined_xy

  subroutine fill_volume_seeds(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz,seeds)
    integer, intent(in) :: nx,ny,nz
    double precision, intent(in) :: xmin,xmax,ymin,ymax,zmin,zmax
    double precision, intent(out) :: seeds(nx*ny*nz,ndim)

    double precision :: dx,dy,dz
    integer :: i,j,k,iseed

    dx=0.d0; dy=0.d0; dz=0.d0
    if (nx>1) dx=(xmax-xmin)/dble(nx-1)
    if (ny>1) dy=(ymax-ymin)/dble(ny-1)
    if (nz>1) dz=(zmax-zmin)/dble(nz-1)

    do k=1,nz
      do j=1,ny
        do i=1,nx
          iseed=(k-1)*nx*ny+(j-1)*nx+i
          seeds(iseed,1)=xmin+dble(i-1)*dx
          seeds(iseed,2)=ymin+dble(j-1)*dy
          seeds(iseed,3)=zmin+dble(k-1)*dz
        enddo
      enddo
    enddo
  end subroutine fill_volume_seeds

  subroutine summarize_volume_length_uniform(unit,xmin,xmax,nx,ymin,ymax, &
       ny,zmin,zmax,nz)
    integer, intent(in) :: unit,nx,ny,nz
    double precision, intent(in) :: xmin,xmax,ymin,ymax,zmin,zmax

    type(trace_length_result), allocatable :: results(:)
    double precision, allocatable :: seeds(:,:)
    integer :: nseed,iseed,valid_count,invalid_count,nonfinite_count
    double precision :: length_err

    nseed=nx*ny*nz
    allocate(seeds(nseed,ndim),results(nseed))
    call fill_volume_seeds(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz,seeds)
    call trace_field_length_multi(seeds,nseed,dL,max_steps,results, &
         b_min=b_min)

    valid_count=0
    invalid_count=0
    nonfinite_count=0
    length_err=0.d0
    do iseed=1,nseed
      if (results(iseed)%backward_status==trace_status_boundary .and. &
           results(iseed)%forward_status==trace_status_boundary) then
        valid_count=valid_count+1
        if (ieee_is_finite(results(iseed)%total_length)) then
          length_err=max(length_err,abs(results(iseed)%total_length-2.d0))
        else
          nonfinite_count=nonfinite_count+1
        endif
      else
        invalid_count=invalid_count+1
      endif
    enddo

    call write_summary_row(unit,'volume_vti','length_uniform','ok', &
         valid_count,invalid_count,nonfinite_count,nan(),nan(),nan(), &
         2.d0,length_err,nan(),nan(),trace_status_boundary, &
         trace_status_boundary,'length_only')
    deallocate(results,seeds)
  end subroutine summarize_volume_length_uniform

  subroutine summarize_volume_all_uniform(unit,xmin,xmax,nx,ymin,ymax,ny, &
       zmin,zmax,nz)
    integer, intent(in) :: unit,nx,ny,nz
    double precision, intent(in) :: xmin,xmax,ymin,ymax,zmin,zmax

    type(trace_length_result), allocatable :: length_results(:)
    type(trace_twist_result), allocatable :: twist_results(:)
    type(trace_qperp_result), allocatable :: qperp_results(:)
    double precision, allocatable :: seeds(:,:)
    integer :: nseed,iseed,valid_length,valid_twist,valid_qperp
    integer :: invalid_length,nonfinite_count,face_forward,face_backward
    double precision :: length_err,twist_err,qmin,qmax,qsum,q_err

    nseed=nx*ny*nz
    allocate(seeds(nseed,ndim),length_results(nseed),twist_results(nseed), &
         qperp_results(nseed))
    call fill_volume_seeds(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz,seeds)
    call trace_field_length_multi(seeds,nseed,dL,max_steps,length_results, &
         b_min=b_min)
    call trace_field_twist_multi(seeds,nseed,dL,max_steps,twist_results, &
         b_min=b_min)
    call trace_field_qperp_multi(seeds,nseed,dL,max_steps,qperp_results, &
         b_min=b_min)

    valid_length=0; invalid_length=0; valid_twist=0; valid_qperp=0
    nonfinite_count=0
    length_err=0.d0; twist_err=0.d0; q_err=0.d0
    qmin=huge(1.d0); qmax=-huge(1.d0); qsum=0.d0
    face_forward=trace_face_zmax
    face_backward=trace_face_zmin
    do iseed=1,nseed
      if (length_results(iseed)%backward_status==trace_status_boundary .and. &
           length_results(iseed)%forward_status==trace_status_boundary) then
        valid_length=valid_length+1
        length_err=max(length_err, &
             abs(length_results(iseed)%total_length-2.d0))
      else
        invalid_length=invalid_length+1
      endif
      if (twist_results(iseed)%line%backward_status==trace_status_boundary &
           .and. twist_results(iseed)%line%forward_status== &
           trace_status_boundary) then
        valid_twist=valid_twist+1
        twist_err=max(twist_err,abs(twist_results(iseed)%total_twist))
      endif
      if (qperp_results(iseed)%valid) then
        valid_qperp=valid_qperp+1
        if (ieee_is_finite(qperp_results(iseed)%qperp)) then
          qmin=min(qmin,qperp_results(iseed)%qperp)
          qmax=max(qmax,qperp_results(iseed)%qperp)
          qsum=qsum+qperp_results(iseed)%qperp
          q_err=max(q_err,abs(qperp_results(iseed)%qperp-2.d0))
        else
          nonfinite_count=nonfinite_count+1
        endif
        if (qperp_results(iseed)%forward_face/=trace_face_zmax) &
             face_forward=qperp_results(iseed)%forward_face
        if (qperp_results(iseed)%backward_face/=trace_face_zmin) &
             face_backward=qperp_results(iseed)%backward_face
      else
        nonfinite_count=nonfinite_count+1
      endif
    enddo

    call write_summary_row(unit,'volume_vti','all_uniform','ok', &
         valid_length,invalid_length,nonfinite_count,qmin,qmax, &
         qsum/dble(max(1,valid_qperp)),2.d0,q_err,nan(),twist_err, &
         trace_face_zmax,face_forward, &
         'twist_valid='//trim(int_to_string(valid_twist))// &
         ';qperp_valid='//trim(int_to_string(valid_qperp))// &
         ';face_backward='//trim(int_to_string(face_backward))// &
         ';length_err='//trim(real_to_string(length_err)))
    deallocate(qperp_results,twist_results,length_results,seeds)
  end subroutine summarize_volume_all_uniform

  subroutine summarize_volume_hyperbolic_center(unit,subcase,xmin,xmax,nx, &
       ymin,ymax,ny,zmin,zmax,nz)
    integer, intent(in) :: unit,nx,ny,nz
    character(len=*), intent(in) :: subcase
    double precision, intent(in) :: xmin,xmax,ymin,ymax,zmin,zmax

    type(trace_qperp_result), allocatable :: results(:)
    double precision, allocatable :: seeds(:,:)
    double precision :: expected,abs_err,rel_err
    integer :: nseed,center

    nseed=nx*ny*nz
    center=(2-1)*nx*ny+(2-1)*nx+2
    allocate(seeds(nseed,ndim),results(nseed))
    call fill_volume_seeds(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz,seeds)
    call trace_field_qperp_multi(seeds,nseed,dL,max_steps,results, &
         b_min=b_min)

    expected=2.d0*dcosh(0.4d0)
    abs_err=abs(results(center)%qperp-expected)
    rel_err=abs_err/expected
    call write_summary_row(unit,'volume_vti',subcase,'ok', &
         merge(1,0,results(center)%valid), &
         merge(0,1,results(center)%valid),0,results(center)%qperp, &
         results(center)%qperp,results(center)%qperp,expected,abs_err, &
         rel_err,nan(),trace_status_boundary,results(center)%status, &
         'center_seed')
    deallocate(results,seeds)
  end subroutine summarize_volume_hyperbolic_center

  subroutine summarize_volume_outside(unit,subcase,xmin,xmax,nx,ymin,ymax, &
       ny,zmin,zmax,nz)
    integer, intent(in) :: unit,nx,ny,nz
    character(len=*), intent(in) :: subcase
    double precision, intent(in) :: xmin,xmax,ymin,ymax,zmin,zmax

    type(trace_qperp_result), allocatable :: results(:)
    double precision, allocatable :: seeds(:,:)
    integer :: nseed,iseed,valid_count,invalid_count,outside_count
    integer :: status_actual,nonfinite_count

    nseed=nx*ny*nz
    allocate(seeds(nseed,ndim),results(nseed))
    call fill_volume_seeds(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz,seeds)
    call trace_field_qperp_multi(seeds,nseed,dL,max_steps,results, &
         b_min=b_min)

    valid_count=0
    invalid_count=0
    outside_count=0
    nonfinite_count=0
    status_actual=trace_status_boundary
    do iseed=1,nseed
      if (results(iseed)%valid) then
        valid_count=valid_count+1
      else
        invalid_count=invalid_count+1
        if (.not.ieee_is_finite(results(iseed)%qperp)) &
             nonfinite_count=nonfinite_count+1
      endif
      if (results(iseed)%status==trace_status_seed_outside) then
        outside_count=outside_count+1
        status_actual=trace_status_seed_outside
      endif
    enddo

    call write_summary_row(unit,'volume_vti',subcase,'ok',valid_count, &
         invalid_count,nonfinite_count,nan(),nan(),nan(),nan(),nan(), &
         nan(),nan(),trace_status_seed_outside,status_actual, &
         'outside='//trim(int_to_string(outside_count)))
    deallocate(results,seeds)
  end subroutine summarize_volume_outside

  subroutine summarize_qperp_plane(unit,file_name,plane)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: file_name,plane
    integer :: io,ios,ix,iy,status,face_f,face_b,status_f,status_b
    integer :: valid_count,invalid_count,nonfinite_count
    double precision :: sx,sy,sz,q,logq,N2,bfactor,len_f,len_b
    double precision :: Bseed_norm,Bf_norm,Bb_norm,xf,yf,zf,xb,yb,zb
    double precision :: qmin,qmax,qsum
    logical :: valid
    character(len=2048) :: line

    valid_count=0; invalid_count=0; nonfinite_count=0
    qmin=huge(1.d0); qmax=-huge(1.d0); qsum=0.d0
    open(newunit=io,file=file_name,status='old',action='read',iostat=ios)
    if (ios/=0) call mpistop('could not open Qperp CSV')
    read(io,'(a)',iostat=ios) line
    do
      read(io,'(a)',iostat=ios) line
      if (ios/=0) exit
      read(line,*) ix,iy,sx,sy,sz,q,logq,valid,status,N2,bfactor, &
           face_f,face_b,status_f,status_b,len_f,len_b,Bseed_norm, &
           Bf_norm,Bb_norm,xf,yf,zf,xb,yb,zb
      if (valid) then
        valid_count=valid_count+1
        if (ieee_is_finite(q)) then
          qmin=min(qmin,q); qmax=max(qmax,q); qsum=qsum+q
        else
          nonfinite_count=nonfinite_count+1
        endif
      else
        invalid_count=invalid_count+1
      endif
    enddo
    close(io)

    call write_summary_row(unit,'qperp_uniform',trim(plane),'ok', &
         valid_count,invalid_count,nonfinite_count,qmin,qmax, &
         qsum/dble(valid_count),2.d0, &
         max(abs(qmin-2.d0),abs(qmax-2.d0)), &
         abs(qsum/dble(valid_count)-2.d0)/2.d0,nan(),0,0,'axis_plane')
  end subroutine summarize_qperp_plane

  subroutine summarize_qperp_arbitrary_eq_xy(unit,axis_file,arbitrary_file)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: axis_file,arbitrary_file
    integer :: io_axis,io_arb,ios_axis,ios_arb
    integer :: ix_a,iy_a,ix_b,iy_b,status_a,status_b
    integer :: face_f_a,face_b_a,status_f_a,status_b_a
    integer :: face_f_b,face_b_b,status_f_b,status_b_b
    integer :: valid_count,invalid_count,nonfinite_count
    double precision :: sx_a,sy_a,sz_a,q_a,logq_a,N2_a,bfactor_a
    double precision :: len_f_a,len_b_a,Bseed_a,Bf_a,Bb_a
    double precision :: xf_a,yf_a,zf_a,xb_a,yb_a,zb_a
    double precision :: s1_b,s2_b,sx_b,sy_b,sz_b,q_b,logq_b,N2_b
    double precision :: bfactor_b,len_f_b,len_b_b,Bseed_b,Bf_b,Bb_b
    double precision :: xf_b,yf_b,zf_b,xb_b,yb_b,zb_b
    double precision :: qmin,qmax,qsum,max_diff,row_diff
    logical :: valid_a,valid_b
    character(len=2048) :: line_axis,line_arb

    valid_count=0; invalid_count=0; nonfinite_count=0
    qmin=huge(1.d0); qmax=-huge(1.d0); qsum=0.d0
    max_diff=0.d0
    open(newunit=io_axis,file=axis_file,status='old',action='read', &
         iostat=ios_axis)
    open(newunit=io_arb,file=arbitrary_file,status='old',action='read', &
         iostat=ios_arb)
    if (ios_axis/=0 .or. ios_arb/=0) call mpistop('could not open arbitrary comparison CSV')
    read(io_axis,'(a)',iostat=ios_axis) line_axis
    read(io_arb,'(a)',iostat=ios_arb) line_arb
    do
      read(io_axis,'(a)',iostat=ios_axis) line_axis
      read(io_arb,'(a)',iostat=ios_arb) line_arb
      if (ios_axis/=ios_arb) call mpistop('arbitrary eq_xy row count mismatch')
      if (ios_axis/=0) exit
      read(line_axis,*) ix_a,iy_a,sx_a,sy_a,sz_a,q_a,logq_a,valid_a, &
           status_a,N2_a,bfactor_a,face_f_a,face_b_a,status_f_a, &
           status_b_a,len_f_a,len_b_a,Bseed_a,Bf_a,Bb_a,xf_a,yf_a,zf_a, &
           xb_a,yb_a,zb_a
      read(line_arb,*) ix_b,iy_b,s1_b,s2_b,sx_b,sy_b,sz_b,q_b,logq_b, &
           valid_b,status_b,N2_b,bfactor_b,face_f_b,face_b_b,status_f_b, &
           status_b_b,len_f_b,len_b_b,Bseed_b,Bf_b,Bb_b,xf_b,yf_b,zf_b, &
           xb_b,yb_b,zb_b
      if (valid_b) then
        valid_count=valid_count+1
        if (ieee_is_finite(q_b)) then
          qmin=min(qmin,q_b); qmax=max(qmax,q_b); qsum=qsum+q_b
        else
          nonfinite_count=nonfinite_count+1
        endif
      else
        invalid_count=invalid_count+1
      endif
      row_diff=0.d0
      row_diff=max(row_diff,abs(dble(ix_a-ix_b)))
      row_diff=max(row_diff,abs(dble(iy_a-iy_b)))
      row_diff=max(row_diff,abs(sx_a-sx_b),abs(sy_a-sy_b),abs(sz_a-sz_b))
      row_diff=max(row_diff,abs(q_a-q_b),abs(logq_a-logq_b))
      row_diff=max(row_diff,abs(N2_a-N2_b),abs(bfactor_a-bfactor_b))
      row_diff=max(row_diff,abs(dble(merge(1,0,valid_a)-merge(1,0,valid_b))))
      row_diff=max(row_diff,abs(dble(status_a-status_b)))
      row_diff=max(row_diff,abs(dble(face_f_a-face_f_b)))
      row_diff=max(row_diff,abs(dble(face_b_a-face_b_b)))
      row_diff=max(row_diff,abs(dble(status_f_a-status_f_b)))
      row_diff=max(row_diff,abs(dble(status_b_a-status_b_b)))
      row_diff=max(row_diff,abs(len_f_a-len_f_b),abs(len_b_a-len_b_b))
      row_diff=max(row_diff,abs(Bseed_a-Bseed_b),abs(Bf_a-Bf_b), &
           abs(Bb_a-Bb_b))
      row_diff=max(row_diff,abs(xf_a-xf_b),abs(yf_a-yf_b),abs(zf_a-zf_b))
      row_diff=max(row_diff,abs(xb_a-xb_b),abs(yb_a-yb_b),abs(zb_a-zb_b))
      max_diff=max(max_diff,row_diff)
    enddo
    close(io_axis)
    close(io_arb)

    call write_summary_row(unit,'qperp_arbitrary','eq_xy','ok', &
         valid_count,invalid_count,nonfinite_count,qmin,qmax, &
         qsum/dble(valid_count),2.d0, &
         max(abs(qmin-2.d0),abs(qmax-2.d0)), &
         abs(qsum/dble(valid_count)-2.d0)/2.d0,max_diff,0,0, &
         'axis_compare')
  end subroutine summarize_qperp_arbitrary_eq_xy

  subroutine summarize_qperp_arbitrary_uniform(unit,file_name,subcase)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: file_name,subcase
    integer :: io,ios,ix,iy,status,face_f,face_b,status_f,status_b
    integer :: valid_count,invalid_count,nonfinite_count
    double precision :: s1,s2,sx,sy,sz,q,logq,N2,bfactor,len_f,len_b
    double precision :: Bseed_norm,Bf_norm,Bb_norm,xf,yf,zf,xb,yb,zb
    double precision :: qmin,qmax,qsum
    logical :: valid
    character(len=2048) :: line

    valid_count=0; invalid_count=0; nonfinite_count=0
    qmin=huge(1.d0); qmax=-huge(1.d0); qsum=0.d0
    open(newunit=io,file=file_name,status='old',action='read',iostat=ios)
    if (ios/=0) call mpistop('could not open arbitrary Qperp CSV')
    read(io,'(a)',iostat=ios) line
    do
      read(io,'(a)',iostat=ios) line
      if (ios/=0) exit
      read(line,*) ix,iy,s1,s2,sx,sy,sz,q,logq,valid,status,N2,bfactor, &
           face_f,face_b,status_f,status_b,len_f,len_b,Bseed_norm, &
           Bf_norm,Bb_norm,xf,yf,zf,xb,yb,zb
      if (valid) then
        valid_count=valid_count+1
        if (ieee_is_finite(q) .and. ieee_is_finite(logq) .and. &
             ieee_is_finite(N2) .and. ieee_is_finite(bfactor)) then
          qmin=min(qmin,q); qmax=max(qmax,q); qsum=qsum+q
        else
          nonfinite_count=nonfinite_count+1
        endif
      else
        invalid_count=invalid_count+1
      endif
    enddo
    close(io)

    call write_summary_row(unit,'qperp_arbitrary',trim(subcase),'ok', &
         valid_count,invalid_count,nonfinite_count,qmin,qmax, &
         qsum/dble(valid_count),2.d0, &
         max(abs(qmin-2.d0),abs(qmax-2.d0)), &
         abs(qsum/dble(valid_count)-2.d0)/2.d0,nan(),0,0, &
         'uniform_tilted')
  end subroutine summarize_qperp_arbitrary_uniform

  subroutine summarize_qperp_arbitrary_outside(unit,file_name)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: file_name
    integer :: io,ios,ix,iy,status,face_f,face_b,status_f,status_b
    integer :: valid_count,invalid_count,nonfinite_count,invalid_status
    double precision :: s1,s2,sx,sy,sz,q,logq,N2,bfactor,len_f,len_b
    double precision :: Bseed_norm,Bf_norm,Bb_norm,xf,yf,zf,xb,yb,zb
    double precision :: qmin,qmax,qsum
    logical :: valid
    character(len=2048) :: line

    valid_count=0; invalid_count=0; nonfinite_count=0
    invalid_status=0
    qmin=huge(1.d0); qmax=-huge(1.d0); qsum=0.d0
    open(newunit=io,file=file_name,status='old',action='read',iostat=ios)
    if (ios/=0) call mpistop('could not open arbitrary outside CSV')
    read(io,'(a)',iostat=ios) line
    do
      read(io,'(a)',iostat=ios) line
      if (ios/=0) exit
      read(line,*) ix,iy,s1,s2,sx,sy,sz,q,logq,valid,status,N2,bfactor, &
           face_f,face_b,status_f,status_b,len_f,len_b,Bseed_norm, &
           Bf_norm,Bb_norm,xf,yf,zf,xb,yb,zb
      if (valid) then
        valid_count=valid_count+1
        if (ieee_is_finite(q)) then
          qmin=min(qmin,q); qmax=max(qmax,q); qsum=qsum+q
        else
          nonfinite_count=nonfinite_count+1
        endif
      else
        invalid_count=invalid_count+1
        if (invalid_status==0) invalid_status=status
        if (.not.ieee_is_finite(q)) nonfinite_count=nonfinite_count+1
      endif
    enddo
    close(io)

    call write_summary_row(unit,'qperp_arbitrary','outside_seed','ok', &
         valid_count,invalid_count,nonfinite_count,qmin,qmax, &
         qsum/dble(valid_count),nan(),nan(),nan(),nan(), &
         trace_status_seed_outside,invalid_status,'outside_seed')
  end subroutine summarize_qperp_arbitrary_outside

  subroutine summarize_fieldline_products_length(unit,seeds,nseed,file_name)
    integer, intent(in) :: unit,nseed
    double precision, intent(in) :: seeds(nseed,ndim)
    character(len=*), intent(in) :: file_name

    type(trace_length_result), allocatable :: length_results(:)
    integer :: io,ios,seed_id,nb,nf,sb,sf
    integer :: valid_count,invalid_count,outside_status
    double precision :: sx,sy,sz,lt,lb,lf,length_err,lowlevel_diff
    character(len=4096) :: line

    allocate(length_results(nseed))
    call trace_field_length_multi(seeds,nseed,dL,max_steps,length_results, &
         b_min=b_min)

    valid_count=0
    invalid_count=0
    outside_status=trace_status_boundary
    length_err=0.d0
    lowlevel_diff=0.d0
    open(newunit=io,file=file_name,status='old',action='read',iostat=ios)
    if (ios/=0) call mpistop('could not open fieldline products length CSV')
    read(io,'(a)',iostat=ios) line
    do
      read(io,'(a)',iostat=ios) line
      if (ios/=0) exit
      read(line,*) seed_id,sx,sy,sz,lt,lb,lf,nb,nf,sb,sf
      if (sb==trace_status_boundary .and. sf==trace_status_boundary) then
        valid_count=valid_count+1
        length_err=max(length_err,abs(lt-2.d0))
      else
        invalid_count=invalid_count+1
      endif
      if (sb==trace_status_seed_outside .or. sf==trace_status_seed_outside) &
           outside_status=trace_status_seed_outside
      lowlevel_diff=max(lowlevel_diff,abs(lt-length_results(seed_id)%total_length))
      lowlevel_diff=max(lowlevel_diff,abs(lb-length_results(seed_id)%backward_length))
      lowlevel_diff=max(lowlevel_diff,abs(lf-length_results(seed_id)%forward_length))
      lowlevel_diff=max(lowlevel_diff,abs(dble(sb-length_results(seed_id)%backward_status)))
      lowlevel_diff=max(lowlevel_diff,abs(dble(sf-length_results(seed_id)%forward_status)))
    enddo
    close(io)
    deallocate(length_results)

    call write_summary_row(unit,'fieldline_products_seeds','length_uniform', &
         'ok',valid_count,invalid_count,0,nan(),nan(),nan(),2.d0, &
         length_err,nan(),lowlevel_diff,trace_status_seed_outside, &
         outside_status,'length_only')
  end subroutine summarize_fieldline_products_length

  subroutine summarize_fieldline_products_all(unit,subcase,seeds,nseed, &
       file_name,q_expected,expect_twist_zero)
    integer, intent(in) :: unit,nseed
    character(len=*), intent(in) :: subcase,file_name
    double precision, intent(in) :: seeds(nseed,ndim),q_expected
    logical, intent(in) :: expect_twist_zero

    type(trace_length_result), allocatable :: length_results(:)
    type(trace_twist_result), allocatable :: twist_results(:)
    type(trace_qperp_result), allocatable :: qperp_results(:)
    integer :: io,ios,seed_id,nbl,nfl,sbl,sfl,nbt,nft,sbt,sft
    integer :: status_q,facef,faceb,statusf,statusb
    integer :: valid_length,invalid_length,valid_twist,valid_qperp
    integer :: nonfinite_count,outside_status
    double precision :: sx,sy,sz,lt,lb,lf,tw,twb,twf,q,logq,N2,bfactor
    double precision :: lenfq,lenbq,Bseed_norm,Bf_norm,Bb_norm
    double precision :: xf,yf,zf,xb,yb,zb,qmin,qmax,qsum,q_center
    double precision :: q_err,rel_err,twist_err,lowlevel_diff,row_diff
    logical :: valid_q
    character(len=4096) :: line

    allocate(length_results(nseed),twist_results(nseed),qperp_results(nseed))
    call trace_field_length_multi(seeds,nseed,dL,max_steps,length_results, &
         b_min=b_min)
    call trace_field_twist_multi(seeds,nseed,dL,max_steps,twist_results, &
         b_min=b_min)
    call trace_field_qperp_multi(seeds,nseed,dL,max_steps,qperp_results, &
         b_min=b_min)

    valid_length=0
    invalid_length=0
    valid_twist=0
    valid_qperp=0
    nonfinite_count=0
    outside_status=trace_status_boundary
    qmin=huge(1.d0)
    qmax=-huge(1.d0)
    qsum=0.d0
    q_center=nan()
    q_err=0.d0
    rel_err=nan()
    twist_err=0.d0
    lowlevel_diff=0.d0
    open(newunit=io,file=file_name,status='old',action='read',iostat=ios)
    if (ios/=0) call mpistop('could not open fieldline products all CSV')
    read(io,'(a)',iostat=ios) line
    do
      read(io,'(a)',iostat=ios) line
      if (ios/=0) exit
      read(line,*) seed_id,sx,sy,sz,lt,lb,lf,nbl,nfl,sbl,sfl, &
           tw,twb,twf,nbt,nft,sbt,sft,q,logq,valid_q,status_q,N2,bfactor, &
           facef,faceb,statusf,statusb,lenfq,lenbq,Bseed_norm,Bf_norm, &
           Bb_norm,xf,yf,zf,xb,yb,zb
      if (sbl==trace_status_boundary .and. sfl==trace_status_boundary) then
        valid_length=valid_length+1
      else
        invalid_length=invalid_length+1
      endif
      if (sbt==trace_status_boundary .and. sft==trace_status_boundary) then
        valid_twist=valid_twist+1
        if (expect_twist_zero) twist_err=max(twist_err,abs(tw))
      endif
      if (valid_q) then
        valid_qperp=valid_qperp+1
        if (ieee_is_finite(q)) then
          qmin=min(qmin,q)
          qmax=max(qmax,q)
          qsum=qsum+q
          q_err=max(q_err,abs(q-q_expected))
        else
          nonfinite_count=nonfinite_count+1
        endif
      else
        if (.not.ieee_is_finite(q)) nonfinite_count=nonfinite_count+1
      endif
      if (sbl==trace_status_seed_outside .or. sfl==trace_status_seed_outside &
           .or. status_q==trace_status_seed_outside) then
        outside_status=trace_status_seed_outside
      endif
      row_diff=0.d0
      row_diff=max(row_diff,abs(lt-length_results(seed_id)%total_length))
      row_diff=max(row_diff,abs(lb-length_results(seed_id)%backward_length))
      row_diff=max(row_diff,abs(lf-length_results(seed_id)%forward_length))
      row_diff=max(row_diff,abs(tw-twist_results(seed_id)%total_twist))
      row_diff=max(row_diff,abs(twb-twist_results(seed_id)%backward_twist))
      row_diff=max(row_diff,abs(twf-twist_results(seed_id)%forward_twist))
      if (ieee_is_finite(q) .and. ieee_is_finite(qperp_results(seed_id)%qperp)) &
           row_diff=max(row_diff,abs(q-qperp_results(seed_id)%qperp))
      lowlevel_diff=max(lowlevel_diff,row_diff)
      if (seed_id==1) q_center=q
    enddo
    close(io)
    deallocate(length_results,twist_results,qperp_results)

    if (ieee_is_finite(q_center) .and. abs(q_expected)>0.d0) then
      rel_err=abs(q_center-q_expected)/abs(q_expected)
    endif
    call write_summary_row(unit,'fieldline_products_seeds',trim(subcase), &
         'ok',valid_length,invalid_length,nonfinite_count,qmin,qmax, &
         qsum/dble(max(1,valid_qperp)),q_expected,q_err,rel_err, &
         lowlevel_diff,trace_status_seed_outside,outside_status, &
         'twist_valid='//trim(int_to_string(valid_twist))// &
         ';qperp_valid='//trim(int_to_string(valid_qperp))// &
         ';twist_err='//trim(real_to_string(twist_err)))
  end subroutine summarize_fieldline_products_all

  subroutine summarize_fieldline_products_plane_header(unit,subcase,file_name, &
       n1,n2)
    integer, intent(in) :: unit,n1,n2
    character(len=*), intent(in) :: subcase,file_name
    logical :: has_twist,has_qperp

    call header_flags(file_name,has_twist,has_qperp)
    call write_summary_row(unit,'fieldline_products_plane_arbitrary', &
         trim(subcase),'ok',0,0,0,nan(),nan(),nan(),nan(),nan(),nan(), &
         0.d0,0,0,'header_has_twist='//trim(int_to_string(merge(1,0,has_twist)))// &
         ';header_has_qperp='//trim(int_to_string(merge(1,0,has_qperp)))// &
         ';nseed='//trim(int_to_string(n1*n2)))
  end subroutine summarize_fieldline_products_plane_header

  subroutine summarize_fieldline_products_plane_invalid_basis(unit,file_name)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: file_name
    integer :: ndata
    logical :: has_twist,has_qperp

    call header_flags(file_name,has_twist,has_qperp)
    call count_csv_data_rows(file_name,ndata)
    call write_summary_row(unit,'fieldline_products_plane_arbitrary', &
         'invalid_basis','ok',0,0,0,nan(),nan(),nan(),nan(),nan(), &
         nan(),nan(),0,0, &
         'header_only='//trim(int_to_string(merge(1,0,ndata==0)))// &
         ';header_has_twist='//trim(int_to_string(merge(1,0,has_twist)))// &
         ';header_has_qperp='//trim(int_to_string(merge(1,0,has_qperp))))
  end subroutine summarize_fieldline_products_plane_invalid_basis

  subroutine summarize_fieldline_products_plane_all(unit,subcase,origin,e1,e2, &
       s1min,s1max,n1,s2min,s2max,n2,file_name,q_expected,expect_twist_zero)
    integer, intent(in) :: unit,n1,n2
    character(len=*), intent(in) :: subcase,file_name
    double precision, intent(in) :: origin(ndim),e1(ndim),e2(ndim)
    double precision, intent(in) :: s1min,s1max,s2min,s2max,q_expected
    logical, intent(in) :: expect_twist_zero

    type(trace_length_result), allocatable :: length_results(:)
    type(trace_twist_result), allocatable :: twist_results(:)
    type(trace_qperp_result), allocatable :: qperp_results(:)
    double precision, allocatable :: seeds(:,:)
    double precision :: qmin,qmax,qsum,qcenter,qerr,relerr
    double precision :: twist_err,lowlevel_diff,row_diff
    integer :: nseed,valid_length,invalid_length,valid_twist,valid_qperp
    integer :: nonfinite_count,outside_count,status_out
    integer :: io,ios,i,j,iseed,nbl,nfl,sbl,sfl,nbt,nft,sbt,sft
    integer :: status_q,facef,faceb,statusf,statusb
    double precision :: s1,s2,sx,sy,sz,lt,lb,lf,tw,twb,twf,q,logq
    double precision :: Ntwo,bfactor,lenfq,lenbq,Bseed_norm,Bf_norm,Bb_norm
    double precision :: xf,yf,zf,xb,yb,zb
    logical :: valid_q,has_twist,has_qperp
    character(len=4096) :: line
    character(len=128) :: fields(72)
    integer :: nfield

    nseed=n1*n2
    allocate(seeds(nseed,ndim),length_results(nseed),twist_results(nseed), &
         qperp_results(nseed))
    call fill_arbitrary_plane_seeds_reg(origin,e1,e2,s1min,s1max,n1, &
         s2min,s2max,n2,seeds)
    call trace_field_length_multi(seeds,nseed,dL,max_steps,length_results, &
         b_min=b_min)
    call trace_field_twist_multi(seeds,nseed,dL,max_steps,twist_results, &
         b_min=b_min)
    call trace_field_qperp_multi(seeds,nseed,dL,max_steps,qperp_results, &
         b_min=b_min)

    valid_length=0; invalid_length=0; valid_twist=0; valid_qperp=0
    nonfinite_count=0; outside_count=0; status_out=trace_status_boundary
    qmin=huge(1.d0); qmax=-huge(1.d0); qsum=0.d0
    qcenter=nan(); qerr=0.d0; relerr=nan(); twist_err=0.d0
    lowlevel_diff=0.d0
    open(newunit=io,file=file_name,status='old',action='read',iostat=ios)
    if (ios/=0) call mpistop('could not open fieldline plane products CSV')
    read(io,'(a)',iostat=ios) line
    do
      read(io,'(a)',iostat=ios) line
      if (ios/=0) exit
      call split_csv_line_reg(line,fields,nfield)
      if (nfield<42) call mpistop('bad fieldline plane products CSV row')
      read(fields(1),*) i; read(fields(2),*) j
      read(fields(3),*) s1; read(fields(4),*) s2
      read(fields(5),*) sx; read(fields(6),*) sy; read(fields(7),*) sz
      read(fields(8),*) lt; read(fields(9),*) lb; read(fields(10),*) lf
      read(fields(11),*) nbl; read(fields(12),*) nfl
      read(fields(13),*) sbl; read(fields(14),*) sfl
      read(fields(15),*) tw; read(fields(16),*) twb; read(fields(17),*) twf
      read(fields(18),*) nbt; read(fields(19),*) nft
      read(fields(20),*) sbt; read(fields(21),*) sft
      read(fields(22),*) q; read(fields(23),*) logq
      valid_q=(trim(fields(24))=='T' .or. trim(fields(24))=='t')
      read(fields(25),*) status_q
      read(fields(26),*) Ntwo; read(fields(27),*) bfactor
      read(fields(28),*) facef; read(fields(29),*) faceb
      read(fields(30),*) statusf; read(fields(31),*) statusb
      read(fields(32),*) lenfq; read(fields(33),*) lenbq
      read(fields(34),*) Bseed_norm; read(fields(35),*) Bf_norm
      read(fields(36),*) Bb_norm
      read(fields(37),*) xf; read(fields(38),*) yf; read(fields(39),*) zf
      read(fields(40),*) xb; read(fields(41),*) yb; read(fields(42),*) zb
      iseed=(j-1)*n1+i
      if (sbl==trace_status_boundary .and. sfl==trace_status_boundary) then
        valid_length=valid_length+1
      else
        invalid_length=invalid_length+1
      endif
      if (sbt==trace_status_boundary .and. sft==trace_status_boundary) then
        valid_twist=valid_twist+1
        if (expect_twist_zero) twist_err=max(twist_err,abs(tw))
      endif
      if (valid_q) then
        valid_qperp=valid_qperp+1
        if (ieee_is_finite(q)) then
          qmin=min(qmin,q); qmax=max(qmax,q); qsum=qsum+q
          if (ieee_is_finite(q_expected)) qerr=max(qerr,abs(q-q_expected))
        else
          nonfinite_count=nonfinite_count+1
        endif
      else
        if (.not.ieee_is_finite(q)) nonfinite_count=nonfinite_count+1
      endif
      if (sbl==trace_status_seed_outside .or. sfl==trace_status_seed_outside &
           .or. status_q==trace_status_seed_outside) then
        outside_count=outside_count+1
        status_out=trace_status_seed_outside
      endif
      row_diff=0.d0
      row_diff=max(row_diff,abs(lt-length_results(iseed)%total_length))
      row_diff=max(row_diff,abs(tw-twist_results(iseed)%total_twist))
      if (ieee_is_finite(q) .and. ieee_is_finite(qperp_results(iseed)%qperp)) &
           row_diff=max(row_diff,abs(q-qperp_results(iseed)%qperp))
      lowlevel_diff=max(lowlevel_diff,row_diff)
    enddo
    close(io)
    iseed=((n2+1)/2-1)*n1+(n1+1)/2
    qcenter=qperp_results(iseed)%qperp
    if (ieee_is_finite(q_expected) .and. abs(q_expected)>0.d0 .and. &
         ieee_is_finite(qcenter)) relerr=abs(qcenter-q_expected)/abs(q_expected)
    call header_flags(file_name,has_twist,has_qperp)
    call write_summary_row(unit,'fieldline_products_plane_arbitrary', &
         trim(subcase),'ok',valid_length,invalid_length,nonfinite_count, &
         qmin,qmax,qsum/dble(max(1,valid_qperp)),q_expected,qerr,relerr, &
         lowlevel_diff,trace_status_seed_outside,status_out, &
         'twist_valid='//trim(int_to_string(valid_twist))// &
         ';qperp_valid='//trim(int_to_string(valid_qperp))// &
         ';outside='//trim(int_to_string(outside_count))// &
         ';twist_err='//trim(real_to_string(twist_err))// &
         ';header_has_twist='//trim(int_to_string(merge(1,0,has_twist)))// &
         ';header_has_qperp='//trim(int_to_string(merge(1,0,has_qperp))))
    deallocate(seeds,length_results,twist_results,qperp_results)
  end subroutine summarize_fieldline_products_plane_all

  subroutine fill_arbitrary_plane_seeds_reg(origin,e1,e2,s1min,s1max,n1, &
       s2min,s2max,n2,seeds)
    integer, intent(in) :: n1,n2
    double precision, intent(in) :: origin(ndim),e1(ndim),e2(ndim)
    double precision, intent(in) :: s1min,s1max,s2min,s2max
    double precision, intent(out) :: seeds(n1*n2,ndim)

    double precision :: ds1,ds2,s1,s2
    integer :: i,j,iseed

    ds1=0.d0; ds2=0.d0
    if (n1>1) ds1=(s1max-s1min)/dble(n1-1)
    if (n2>1) ds2=(s2max-s2min)/dble(n2-1)
    do j=1,n2
      s2=s2min+dble(j-1)*ds2
      do i=1,n1
        s1=s1min+dble(i-1)*ds1
        iseed=(j-1)*n1+i
        seeds(iseed,:)=origin+s1*e1+s2*e2
      enddo
    enddo
  end subroutine fill_arbitrary_plane_seeds_reg

  subroutine header_flags(file_name,has_twist,has_qperp)
    character(len=*), intent(in) :: file_name
    logical, intent(out) :: has_twist,has_qperp

    integer :: io,ios
    character(len=4096) :: line

    open(newunit=io,file=file_name,status='old',action='read',iostat=ios)
    if (ios/=0) call mpistop('could not open CSV header')
    read(io,'(a)',iostat=ios) line
    close(io)
    has_twist=index(line,'twist_total')>0
    has_qperp=index(line,'qperp')>0
  end subroutine header_flags

  subroutine count_csv_data_rows(file_name,ndata)
    character(len=*), intent(in) :: file_name
    integer, intent(out) :: ndata

    integer :: io,ios
    character(len=4096) :: line

    ndata=0
    open(newunit=io,file=file_name,status='old',action='read',iostat=ios)
    if (ios/=0) call mpistop('could not count CSV rows')
    read(io,'(a)',iostat=ios) line
    do
      read(io,'(a)',iostat=ios) line
      if (ios/=0) exit
      ndata=ndata+1
    enddo
    close(io)
  end subroutine count_csv_data_rows

  subroutine split_csv_line_reg(line,fields,nfield)
    character(len=*), intent(in) :: line
    character(len=*), intent(out) :: fields(:)
    integer, intent(out) :: nfield

    integer :: start,comma,line_len

    fields=''
    nfield=0
    start=1
    line_len=len_trim(line)
    do while (start<=line_len+1 .and. nfield<size(fields))
      comma=index(line(start:),',')
      nfield=nfield+1
      if (comma==0) then
        fields(nfield)=adjustl(line(start:line_len))
        exit
      else
        fields(nfield)=adjustl(line(start:start+comma-2))
        start=start+comma
      endif
    enddo
  end subroutine split_csv_line_reg

  subroutine summarize_exact_file_compare(unit,subcase,file_a,file_b)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: subcase,file_a,file_b
    logical :: same

    call files_exactly_equal(file_a,file_b,same)
    call write_summary_row(unit,'topology_combined_xy',trim(subcase),'ok', &
         0,0,0,nan(),nan(),nan(),nan(),nan(),nan(), &
         merge(0.d0,1.d0,same),0,0,'exact_file_compare')
  end subroutine summarize_exact_file_compare

  subroutine files_exactly_equal(file_a,file_b,same)
    character(len=*), intent(in) :: file_a,file_b
    logical, intent(out) :: same
    integer :: unit_a,unit_b,ios_a,ios_b
    character(len=4096) :: line_a,line_b

    same=.true.
    open(newunit=unit_a,file=file_a,status='old',action='read',iostat=ios_a)
    open(newunit=unit_b,file=file_b,status='old',action='read',iostat=ios_b)
    if (ios_a/=0 .or. ios_b/=0) then
      same=.false.
      return
    endif
    do
      read(unit_a,'(a)',iostat=ios_a) line_a
      read(unit_b,'(a)',iostat=ios_b) line_b
      if (ios_a/=ios_b) then
        same=.false.
        exit
      endif
      if (ios_a/=0) exit
      if (trim(line_a)/=trim(line_b)) then
        same=.false.
        exit
      endif
    enddo
    close(unit_a)
    close(unit_b)
  end subroutine files_exactly_equal

  function qperp_result_max_abs_diff(a,b) result(diff)
    type(trace_qperp_result), intent(in) :: a,b
    double precision :: diff

    diff=0.d0
    diff=max(diff,maxval(abs(a%forward_endpoint-b%forward_endpoint)))
    diff=max(diff,maxval(abs(a%backward_endpoint-b%backward_endpoint)))
    diff=max(diff,abs(a%forward_length-b%forward_length))
    diff=max(diff,abs(a%backward_length-b%backward_length))
    diff=max(diff,abs(a%qperp-b%qperp))
    diff=max(diff,abs(a%logqperp-b%logqperp))
    diff=max(diff,abs(a%N2-b%N2))
    diff=max(diff,abs(a%bfactor-b%bfactor))
    diff=max(diff,abs(dsqrt(sum(a%B_seed**2))-dsqrt(sum(b%B_seed**2))))
    diff=max(diff,abs(dsqrt(sum(a%forward_B**2))- &
         dsqrt(sum(b%forward_B**2))))
    diff=max(diff,abs(dsqrt(sum(a%backward_B**2))- &
         dsqrt(sum(b%backward_B**2))))
  end function qperp_result_max_abs_diff

  subroutine write_summary_row(unit,test_name,subcase,status,valid_count, &
       invalid_count,nonfinite_count,q_min,q_max,q_mean,expected,abs_err, &
       rel_err,max_abs_diff,status_expected,status_actual,extra)
    integer, intent(in) :: unit,valid_count,invalid_count,nonfinite_count
    integer, intent(in) :: status_expected,status_actual
    character(len=*), intent(in) :: test_name,subcase,status,extra
    double precision, intent(in) :: q_min,q_max,q_mean,expected,abs_err
    double precision, intent(in) :: rel_err,max_abs_diff

    write(unit,'(a,",",a,",",a,",",i0,",",i0,",",i0,7(",",es24.16),'// &
         '",",i0,",",i0,",",a)') trim(test_name),trim(subcase), &
         trim(status),valid_count,invalid_count,nonfinite_count,q_min, &
         q_max,q_mean,expected,abs_err,rel_err,max_abs_diff, &
         status_expected,status_actual,trim(extra)
  end subroutine write_summary_row

  function nan() result(value)
    double precision :: value

    value=ieee_value(0.d0,ieee_quiet_nan)
  end function nan

  function int_to_string(value) result(text)
    integer, intent(in) :: value
    character(len=32) :: text

    write(text,'(i0)') value
  end function int_to_string

  function real_to_string(value) result(text)
    double precision, intent(in) :: value
    character(len=32) :: text

    write(text,'(es12.5)') value
  end function real_to_string

end module mod_usr
