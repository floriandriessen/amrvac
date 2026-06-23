module mod_trace_field
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan, &
       ieee_is_finite
  use mod_global_parameters
  use mod_physics
  implicit none

  integer, parameter, public :: trace_status_active               = 0
  integer, parameter, public :: trace_status_boundary             = 1
  integer, parameter, public :: trace_status_weak_field           = 2
  integer, parameter, public :: trace_status_max_steps            = 3
  integer, parameter, public :: trace_status_seed_outside         = 4
  integer, parameter, public :: trace_status_out_of_domain        = 5
  integer, parameter, public :: trace_status_invalid_input        = 6
  integer, parameter, public :: trace_status_unsupported_geometry = 7
  integer, parameter, public :: trace_status_trac_stop            = 8
  integer, parameter, public :: trace_status_mpi_unsupported      = 9
  integer, parameter, public :: trace_status_bad_curl_stencil     = 10
  integer, parameter, public :: trace_status_bad_grad_stencil     = 11
  integer, parameter, public :: trace_status_bad_face_limit_sample= 12

  integer, parameter, public :: trace_face_none      = 0
  integer, parameter, public :: trace_face_xmin      = 1
  integer, parameter, public :: trace_face_xmax      = 2
  integer, parameter, public :: trace_face_ymin      = 3
  integer, parameter, public :: trace_face_ymax      = 4
  integer, parameter, public :: trace_face_zmin      = 5
  integer, parameter, public :: trace_face_zmax      = 6
  integer, parameter, public :: trace_face_ambiguous = 7

  public :: trace_field_length_single,trace_field_length_multi
  public :: trace_field_twist_single,trace_field_twist_multi
  public :: trace_field_mapping_single,trace_field_mapping_multi
  public :: trace_field_topology_multi
  public :: trace_field_qperp_single,trace_field_qperp_multi
  public :: trace_debug_sample_bhat_gradbhat
  public :: trace_debug_compare_gradbhat
  public :: trace_debug_compare_gradbhat_methods
  public :: trace_debug_transport_tangents
  public :: trace_debug_transport_tangents_to_boundary
  public :: trace_debug_sample_endpoint_B_face_limit
  public :: trace_debug_qperp_single
  private :: trace_summary_seed,trace_summary_multi
  private :: trace_summary_twist_seed,trace_summary_twist_multi
  private :: trace_summary_mapping_seed,trace_summary_mapping_multi
  private :: trace_summary_topology_multi
  private :: trace_summary_validate,trace_summary_init_result
  private :: trace_summary_init_twist_result
  private :: trace_summary_init_mapping_result
  private :: trace_summary_init_topology_result
  private :: trace_summary_locate_seed,trace_summary_trace_seed
  private :: trace_summary_trace_twist_seed
  private :: trace_summary_init_state,trace_summary_advance_state
  private :: trace_summary_trace_state_to_end,trace_summary_trace_states_grouped
  private :: trace_summary_advance_state_in_grid
  private :: trace_summary_rk2_step,trace_intersect_domain
  private :: trace_summary_accumulate_twist,sample_B_curlB_at_point
  private :: sample_B_at_point,trace_summary_fill_mapping
  private :: trace_total_B_at_cell,trace_bhat_at_cell
  private :: sample_bhat_gradbhat_at_point
  private :: sample_bhat_gradbhat_cellfd_at_point
  private :: sample_bhat_gradbhat_xeps_at_point
  private :: sample_bhat_gradbhat_interpderiv_at_point
  private :: trace_debug_locate_point,trace_locate_point_with_hint
  private :: trace_tangent_rhs
  private :: trace_advance_tangent_state_rk2,trace_project_to_perp_bhat
  private :: trace_tangent_rk2_trial,trace_endpoint_B_bhat
  private :: trace_tangent_rk2_trial_from_rhs
  private :: trace_sample_B_on_domain_face_limit
  private :: trace_locate_face_limit_grid
  private :: trace_face_normal,trace_project_to_boundary_face
  private :: trace_make_perp_basis,trace_debug_nan
  private :: trace_face_from_dim_side,trace_face_Bn,trace_face_is_boundary
  private :: trace_boundary_face_at_point
  private :: trace_init_qperp_result,trace_qperp_single_core
  private :: trace_qperp_compute_q0_scalars
  private :: trace_qperp_project_to_boundary_face
  private :: trace_tangent_state_init
  private :: trace_tangent_trace_states_grouped
  private :: trace_tangent_state_advance_in_grid
  private :: trace_tangent_state_finalize_boundary
  private :: trace_qperp_prepare_seed_result
  private :: trace_qperp_compute_scalars
  private :: trace_qperp_finalize_from_states

  type, private :: trace_summary_state
    double precision :: x(ndim)
    double precision :: footpoint(ndim)
    double precision :: length
    double precision :: twist
    integer :: nstep
    integer :: status
    integer :: igrid
    integer :: seed_id
    integer :: face
    logical :: forward
    logical :: active
    logical :: accumulate_twist
  end type trace_summary_state

  type, private :: trace_tangent_state
    integer :: seed_id
    integer :: direction
    logical :: active
    logical :: complete
    double precision :: x(ndim)
    double precision :: u(ndim)
    double precision :: v(ndim)
    double precision :: length
    integer :: nstep
    integer :: igrid
    integer :: status
    integer :: face_id
    double precision :: endpoint(ndim)
    double precision :: endpoint_B(ndim)
    double precision :: endpoint_bhat(ndim)
    double precision :: u_perp(ndim)
    double precision :: v_perp(ndim)
  end type trace_tangent_state

  type, public :: trace_length_result
    double precision :: seed(ndim)
    double precision :: forward_footpoint(ndim)
    double precision :: backward_footpoint(ndim)
    double precision :: forward_length
    double precision :: backward_length
    double precision :: total_length
    integer :: forward_nstep
    integer :: backward_nstep
    integer :: forward_status
    integer :: backward_status
  end type trace_length_result

  type, public :: trace_twist_result
    type(trace_length_result) :: line
    double precision :: forward_twist
    double precision :: backward_twist
    double precision :: total_twist
  end type trace_twist_result

  type, public :: trace_mapping_result
    double precision :: seed(ndim)
    double precision :: source_B(3)
    double precision :: forward_footpoint(ndim)
    double precision :: backward_footpoint(ndim)
    double precision :: forward_B(3)
    double precision :: backward_B(3)
    double precision :: forward_length
    double precision :: backward_length
    double precision :: source_Bn
    double precision :: forward_Bn
    double precision :: backward_Bn
    integer :: forward_face
    integer :: backward_face
    integer :: forward_status
    integer :: backward_status
    logical :: valid
  end type trace_mapping_result

  type, public :: trace_qperp_result
    double precision :: seed(ndim)
    double precision :: qperp
    double precision :: logqperp
    double precision :: N2
    double precision :: bfactor
    double precision :: q0
    double precision :: logq0
    double precision :: qperp0
    double precision :: logqperp0
    double precision :: N2_qperp0
    double precision :: bfactor_qperp0
    double precision :: forward_Bn_q0
    double precision :: backward_Bn_q0
    double precision :: B_seed(ndim)
    double precision :: bhat_seed(ndim)
    double precision :: forward_endpoint(ndim)
    double precision :: backward_endpoint(ndim)
    double precision :: forward_B(ndim)
    double precision :: backward_B(ndim)
    double precision :: forward_bhat(ndim)
    double precision :: backward_bhat(ndim)
    double precision :: forward_length
    double precision :: backward_length
    double precision :: u0(ndim)
    double precision :: v0(ndim)
    double precision :: u_forward_perp(ndim)
    double precision :: v_forward_perp(ndim)
    double precision :: u_backward_perp(ndim)
    double precision :: v_backward_perp(ndim)
    integer :: forward_face
    integer :: backward_face
    integer :: forward_status
    integer :: backward_status
    integer :: status_q0
    integer :: status_qperp0
    integer :: status
    logical :: valid
    logical :: valid_q0
    logical :: valid_qperp0
  end type trace_qperp_result

  type, public :: trace_topology_result
    double precision :: seed(ndim)
    double precision :: length_forward
    double precision :: length_backward
    double precision :: length_total
    double precision :: twist_forward
    double precision :: twist_backward
    double precision :: twist_total
    double precision :: forward_endpoint(ndim)
    double precision :: backward_endpoint(ndim)
    integer :: forward_nstep
    integer :: backward_nstep
    integer :: forward_face
    integer :: backward_face
    integer :: forward_status
    integer :: backward_status
    double precision :: map_forward_endpoint(ndim)
    double precision :: map_backward_endpoint(ndim)
    double precision :: map_forward_length
    double precision :: map_backward_length
    integer :: map_forward_face
    integer :: map_backward_face
    integer :: map_forward_status
    integer :: map_backward_status
    double precision :: source_B(3)
    double precision :: forward_B(3)
    double precision :: backward_B(3)
    double precision :: source_Bn
    double precision :: forward_Bn
    double precision :: backward_Bn
    logical :: has_twist
    logical :: has_mapping
    logical :: valid
    integer :: status
  end type trace_topology_result

contains

  subroutine trace_field_length_single(seed,dL,max_steps,result,b_min)
    ! Trace one magnetic field line in both directions without storing its path.
    ! This first implementation supports uniform Cartesian grids and npe=1.
    double precision, intent(in) :: seed(ndim),dL
    integer, intent(in) :: max_steps
    type(trace_length_result), intent(out) :: result
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call trace_summary_seed(seed,dL,max_steps,result,b_min)
    else
      call trace_summary_seed(seed,dL,max_steps,result)
    endif
  end subroutine trace_field_length_single

  subroutine trace_field_twist_single(seed,dL,max_steps,result,b_min)
    ! Trace one magnetic field line and integrate twist without storing its path.
    double precision, intent(in) :: seed(ndim),dL
    integer, intent(in) :: max_steps
    type(trace_twist_result), intent(out) :: result
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call trace_summary_twist_seed(seed,dL,max_steps,result,b_min)
    else
      call trace_summary_twist_seed(seed,dL,max_steps,result)
    endif
  end subroutine trace_field_twist_single

  subroutine trace_field_mapping_single(seed,dL,max_steps,result,b_min, &
       source_normal)
    ! Trace one seed and return endpoint metadata for future mapping/Q work.
    double precision, intent(in) :: seed(ndim),dL
    integer, intent(in) :: max_steps
    type(trace_mapping_result), intent(out) :: result
    double precision, intent(in), optional :: b_min
    double precision, intent(in), optional :: source_normal(3)

    if (present(b_min)) then
      if (present(source_normal)) then
        call trace_summary_mapping_seed(seed,dL,max_steps,result,b_min, &
             source_normal)
      else
        call trace_summary_mapping_seed(seed,dL,max_steps,result,b_min)
      endif
    else
      if (present(source_normal)) then
        call trace_summary_mapping_seed(seed,dL,max_steps,result, &
             source_normal=source_normal)
      else
        call trace_summary_mapping_seed(seed,dL,max_steps,result)
      endif
    endif
  end subroutine trace_field_mapping_single

  subroutine trace_field_qperp_single(seed,dL,max_steps,result,b_min)
    ! Trace one seed with Method-II tangent transport and return Q_perp.
    double precision, intent(in) :: seed(ndim),dL
    integer, intent(in) :: max_steps
    type(trace_qperp_result), intent(out) :: result
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call trace_qperp_single_core(seed,dL,max_steps,result,b_min)
    else
      call trace_qperp_single_core(seed,dL,max_steps,result)
    endif
  end subroutine trace_field_qperp_single

  subroutine trace_debug_sample_bhat_gradbhat(seed,bhat,grad_bhat,status,b_min)
    ! Debug/validation hook for the future Method-II sampler.
    double precision, intent(in) :: seed(ndim)
    double precision, intent(out) :: bhat(3),grad_bhat(3,3)
    integer, intent(out) :: status
    double precision, intent(in), optional :: b_min

    type(trace_length_result) :: located
    double precision :: field_min
    integer :: igrid

    bhat=zero
    grad_bhat=zero
    status=trace_status_active
    if (npe/=1) then
      status=trace_status_mpi_unsupported
      return
    else if (ndim/=3 .or. .not.slab_uniform) then
      status=trace_status_unsupported_geometry
      return
    endif

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)

    call trace_summary_locate_seed(seed,located,igrid)
    status=located%forward_status
    if (status/=trace_status_active) return

    call sample_bhat_gradbhat_at_point(seed,igrid,field_min,bhat, &
         grad_bhat,status)
  end subroutine trace_debug_sample_bhat_gradbhat

  subroutine trace_debug_sample_endpoint_B_face_limit(xhit,face_id,B,bhat, &
       status,b_min)
    ! Debug/validation hook for the Qperp endpoint face-limit B sampler.
    double precision, intent(in) :: xhit(ndim)
    integer, intent(in) :: face_id
    double precision, intent(out) :: B(ndim),bhat(ndim)
    integer, intent(out) :: status
    double precision, intent(in), optional :: b_min

    double precision :: field_min,xprobe(ndim),B3(3),bhat3(3)
    double precision :: xface,span
    integer :: normal_dim,side,igrid

    B=zero
    bhat=zero
    status=trace_status_active
    if (npe/=1) then
      status=trace_status_mpi_unsupported
      return
    else if (ndim/=3 .or. .not.slab_uniform) then
      status=trace_status_unsupported_geometry
      return
    endif
    if (.not.trace_face_is_boundary(face_id)) then
      status=trace_status_bad_face_limit_sample
      return
    endif

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)

    {^IFTHREED
    select case(face_id)
    case(trace_face_xmin)
      normal_dim=1
      side=-1
      xface=xprobmin1
    case(trace_face_xmax)
      normal_dim=1
      side=1
      xface=xprobmax1
    case(trace_face_ymin)
      normal_dim=2
      side=-1
      xface=xprobmin2
    case(trace_face_ymax)
      normal_dim=2
      side=1
      xface=xprobmax2
    case(trace_face_zmin)
      normal_dim=3
      side=-1
      xface=xprobmin3
    case(trace_face_zmax)
      normal_dim=3
      side=1
      xface=xprobmax3
    case default
      status=trace_status_bad_face_limit_sample
      return
    end select

    span=max(abs(xprobmax1-xprobmin1), &
         max(abs(xprobmax2-xprobmin2),abs(xprobmax3-xprobmin3)))
    xprobe=xhit
    xprobe(normal_dim)=xface-side*1.d-6*max(one,span)
    call trace_debug_locate_point(xprobe,igrid,status)
    if (status/=trace_status_active) return

    call trace_sample_B_on_domain_face_limit(xhit,face_id,igrid,field_min, &
         B3,bhat3,status)
    if (status/=trace_status_active) return
    B=B3(1:ndim)
    bhat=bhat3(1:ndim)
    status=trace_status_boundary
    }
  end subroutine trace_debug_sample_endpoint_B_face_limit

  subroutine trace_debug_compare_gradbhat(seed,eps,bhat_current, &
       grad_current,status_current,bhat_xeps,grad_xeps,status_xeps,b_min)
    ! Debug/validation hook comparing production and x+-eps grad(bhat).
    double precision, intent(in) :: seed(ndim),eps
    double precision, intent(out) :: bhat_current(3),grad_current(3,3)
    double precision, intent(out) :: bhat_xeps(3),grad_xeps(3,3)
    integer, intent(out) :: status_current,status_xeps
    double precision, intent(in), optional :: b_min

    type(trace_length_result) :: located
    double precision :: field_min
    integer :: igrid

    bhat_current=zero
    grad_current=zero
    bhat_xeps=zero
    grad_xeps=zero
    status_current=trace_status_active
    status_xeps=trace_status_active
    if (npe/=1) then
      status_current=trace_status_mpi_unsupported
      status_xeps=trace_status_mpi_unsupported
      return
    else if (ndim/=3 .or. .not.slab_uniform) then
      status_current=trace_status_unsupported_geometry
      status_xeps=trace_status_unsupported_geometry
      return
    endif

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)

    call trace_summary_locate_seed(seed,located,igrid)
    status_current=located%forward_status
    if (status_current==trace_status_active) then
      call sample_bhat_gradbhat_at_point(seed,igrid,field_min, &
           bhat_current,grad_current,status_current)
    endif

    call sample_bhat_gradbhat_xeps_at_point(seed,eps,field_min,bhat_xeps, &
         grad_xeps,status_xeps)
  end subroutine trace_debug_compare_gradbhat

  subroutine trace_debug_compare_gradbhat_methods(seed,eps,bhat_current, &
       grad_current,status_current,bhat_xeps,grad_xeps,status_xeps, &
       bhat_interp,grad_interp,status_interp,b_min)
    ! Debug/validation hook comparing all grad(bhat) sampler candidates.
    double precision, intent(in) :: seed(ndim),eps
    double precision, intent(out) :: bhat_current(3),grad_current(3,3)
    double precision, intent(out) :: bhat_xeps(3),grad_xeps(3,3)
    double precision, intent(out) :: bhat_interp(3),grad_interp(3,3)
    integer, intent(out) :: status_current,status_xeps,status_interp
    double precision, intent(in), optional :: b_min

    type(trace_length_result) :: located
    double precision :: field_min
    integer :: igrid

    bhat_current=zero
    grad_current=zero
    bhat_xeps=zero
    grad_xeps=zero
    bhat_interp=zero
    grad_interp=zero
    status_current=trace_status_active
    status_xeps=trace_status_active
    status_interp=trace_status_active
    if (npe/=1) then
      status_current=trace_status_mpi_unsupported
      status_xeps=trace_status_mpi_unsupported
      status_interp=trace_status_mpi_unsupported
      return
    else if (ndim/=3 .or. .not.slab_uniform) then
      status_current=trace_status_unsupported_geometry
      status_xeps=trace_status_unsupported_geometry
      status_interp=trace_status_unsupported_geometry
      return
    endif

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)

    call trace_summary_locate_seed(seed,located,igrid)
    status_current=located%forward_status
    status_interp=located%forward_status
    if (status_current==trace_status_active) then
      call sample_bhat_gradbhat_cellfd_at_point(seed,igrid,field_min, &
           bhat_current,grad_current,status_current)
    endif
    if (status_interp==trace_status_active) then
      call sample_bhat_gradbhat_interpderiv_at_point(seed,igrid,field_min, &
           bhat_interp,grad_interp,status_interp)
    endif

    call sample_bhat_gradbhat_xeps_at_point(seed,eps,field_min,bhat_xeps, &
         grad_xeps,status_xeps)
  end subroutine trace_debug_compare_gradbhat_methods

  subroutine trace_debug_transport_tangents(seed,u0,v0,dL,nstep,forward, &
       x_end,u_end,v_end,bhat_end,u_perp_end,v_perp_end,status,b_min)
    ! Debug/validation hook for future Method-II tangent transport.
    double precision, intent(in) :: seed(ndim),u0(ndim),v0(ndim),dL
    integer, intent(in) :: nstep
    logical, intent(in) :: forward
    double precision, intent(out) :: x_end(ndim),u_end(ndim),v_end(ndim)
    double precision, intent(out) :: bhat_end(ndim)
    double precision, intent(out) :: u_perp_end(ndim),v_perp_end(ndim)
    integer, intent(out) :: status
    double precision, intent(in), optional :: b_min

    double precision :: field_min,h,bhat3(3),grad_end(3,3)
    integer :: igrid,istep

    x_end=seed
    u_end=u0
    v_end=v0
    bhat_end=zero
    u_perp_end=zero
    v_perp_end=zero

    status=trace_status_active
    if (npe/=1) then
      status=trace_status_mpi_unsupported
      return
    else if (ndim/=3 .or. .not.slab_uniform) then
      status=trace_status_unsupported_geometry
      return
    else if (dL<=zero .or. nstep<0) then
      status=trace_status_invalid_input
      return
    endif

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)
    h=dL
    if (.not.forward) h=-dL

    call trace_debug_locate_point(x_end,igrid,status)
    if (status/=trace_status_active) return

    do istep=1,nstep
      call trace_advance_tangent_state_rk2(x_end,u_end,v_end,igrid,h, &
           field_min,status)
      if (status/=trace_status_active) exit
    enddo

    if (status==trace_status_active) then
      {^IFTHREED
      call sample_bhat_gradbhat_at_point(x_end,igrid,field_min,bhat3, &
           grad_end,status)
      if (status==trace_status_active) bhat_end=bhat3(1:ndim)
      }
    endif
    if (status/=trace_status_active) return

    call trace_project_to_perp_bhat(u_end,bhat_end,u_perp_end)
    call trace_project_to_perp_bhat(v_end,bhat_end,v_perp_end)
  end subroutine trace_debug_transport_tangents

  subroutine trace_debug_transport_tangents_to_boundary(seed,u0,v0,dL, &
       max_steps,forward,x_end,u_end,v_end,bhat_end,b_end,u_perp_end, &
       v_perp_end,u_face_end,v_face_end,length,face_id,status,b_min)
    ! Debug/validation hook for Method-II style tangent transport to boundary.
    double precision, intent(in) :: seed(ndim),u0(ndim),v0(ndim),dL
    integer, intent(in) :: max_steps
    logical, intent(in) :: forward
    double precision, intent(out) :: x_end(ndim),u_end(ndim),v_end(ndim)
    double precision, intent(out) :: bhat_end(ndim),b_end(ndim)
    double precision, intent(out) :: u_perp_end(ndim),v_perp_end(ndim)
    double precision, intent(out) :: u_face_end(ndim),v_face_end(ndim)
    double precision, intent(out) :: length
    integer, intent(out) :: face_id,status
    double precision, intent(in), optional :: b_min

    double precision :: field_min,h,h_partial,partial_length
    double precision :: xprobe(ndim),xtrial(ndim),utrial(ndim),vtrial(ndim)
    double precision :: xhit(ndim),xold(ndim),uold(ndim),vold(ndim)
    double precision :: kx1(ndim),ku1(ndim),kv1(ndim)
    integer :: igrid,igrid_trial,istep,point_domain
    logical :: hit_ok

    x_end=seed
    u_end=u0
    v_end=v0
    bhat_end=zero
    b_end=zero
    u_perp_end=zero
    v_perp_end=zero
    u_face_end=zero
    v_face_end=zero
    length=zero
    face_id=trace_face_none

    status=trace_status_active
    if (npe/=1) then
      status=trace_status_mpi_unsupported
      return
    else if (ndim/=3 .or. .not.slab_uniform) then
      status=trace_status_unsupported_geometry
      return
    else if (dL<=zero .or. max_steps<0) then
      status=trace_status_invalid_input
      return
    endif

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)
    h=dL
    if (.not.forward) h=-dL

    call trace_debug_locate_point(x_end,igrid,status)
    if (status/=trace_status_active) return

    do istep=1,max_steps
      xold=x_end
      uold=u_end
      vold=v_end

      call trace_tangent_rhs(xold,uold,vold,igrid,field_min,kx1,ku1,kv1, &
           status)
      if (status/=trace_status_active) return
      xprobe=xold+h*kx1

      point_domain=0
      {if (xprobe(^DB)>=xprobmin^DB .and. xprobe(^DB)<xprobmax^DB) point_domain=point_domain+1\}
      if (point_domain/=ndim) then
        call trace_intersect_domain(xold,xprobe,xhit,hit_ok,face_id)
        if (.not.hit_ok) then
          call trace_boundary_face_at_point(xold,face_id,hit_ok)
          if (.not.hit_ok) then
            status=trace_status_boundary
            x_end=xold
            return
          endif
          x_end=xold
          u_end=uold
          v_end=vold
          status=trace_status_boundary
          exit
        endif
        partial_length=dsqrt(sum((xhit-xold)**2))
        if (partial_length<=100.d0*epsilon(one)*max(one,abs(h))) then
          x_end=xhit
          u_end=uold
          v_end=vold
          status=trace_status_boundary
          exit
        endif
        h_partial=sign(partial_length,h)
        call trace_tangent_rk2_trial_from_rhs(xold,uold,vold,igrid, &
             h_partial,field_min,kx1,ku1,kv1,xtrial,utrial,vtrial,status)
        if (status/=trace_status_active) return
        x_end=xhit
        u_end=utrial
        v_end=vtrial
        length=length+partial_length
        status=trace_status_boundary
        exit
      endif

      call trace_tangent_rk2_trial_from_rhs(xold,uold,vold,igrid,h, &
           field_min,kx1,ku1,kv1,xtrial,utrial,vtrial,status)
      if (status/=trace_status_active) return

      point_domain=0
      {if (xtrial(^DB)>=xprobmin^DB .and. xtrial(^DB)<xprobmax^DB) point_domain=point_domain+1\}
      if (point_domain==ndim) then
        call trace_locate_point_with_hint(xtrial,igrid,igrid_trial,status)
        if (status/=trace_status_active) return
        x_end=xtrial
        u_end=utrial
        v_end=vtrial
        igrid=igrid_trial
        length=length+abs(h)
      else
        call trace_intersect_domain(xold,xtrial,xhit,hit_ok,face_id)
        if (.not.hit_ok) then
          call trace_boundary_face_at_point(xold,face_id,hit_ok)
          if (.not.hit_ok) then
            status=trace_status_boundary
            x_end=xold
            return
          endif
          x_end=xold
          u_end=uold
          v_end=vold
          status=trace_status_boundary
          exit
        endif
        partial_length=dsqrt(sum((xhit-xold)**2))
        if (partial_length<=100.d0*epsilon(one)*max(one,abs(h))) then
          x_end=xhit
          u_end=uold
          v_end=vold
          status=trace_status_boundary
          exit
        endif
        h_partial=sign(partial_length,h)
        call trace_tangent_rk2_trial_from_rhs(xold,uold,vold,igrid, &
             h_partial,field_min,kx1,ku1,kv1,xtrial,utrial,vtrial,status)
        if (status/=trace_status_active) return
        x_end=xhit
        u_end=utrial
        v_end=vtrial
        length=length+partial_length
        status=trace_status_boundary
        exit
      endif
    enddo

    if (status==trace_status_active) then
      status=trace_status_max_steps
      return
    endif
    if (status/=trace_status_boundary) return

    call trace_endpoint_B_bhat(x_end,face_id,igrid,field_min,b_end, &
         bhat_end,status)
    if (status/=trace_status_boundary) return

    call trace_project_to_perp_bhat(u_end,bhat_end,u_perp_end)
    call trace_project_to_perp_bhat(v_end,bhat_end,v_perp_end)
    call trace_project_to_boundary_face(u_end,bhat_end,face_id,u_face_end)
    call trace_project_to_boundary_face(v_end,bhat_end,face_id,v_face_end)
  end subroutine trace_debug_transport_tangents_to_boundary

  subroutine trace_debug_qperp_single(seed,dL,max_steps,qperp,logqperp, &
       x_f,x_b,u_f_perp,v_f_perp,u_b_perp,v_b_perp,b_seed,b_f,b_b, &
       length_f,length_b,face_f,face_b,status,N2,bfactor,b_min)
    ! Debug/validation hook for single-seed Method-II Q_perp diagnostics.
    double precision, intent(in) :: seed(ndim),dL
    integer, intent(in) :: max_steps
    double precision, intent(out) :: qperp,logqperp,N2,bfactor
    double precision, intent(out) :: x_f(ndim),x_b(ndim)
    double precision, intent(out) :: u_f_perp(ndim),v_f_perp(ndim)
    double precision, intent(out) :: u_b_perp(ndim),v_b_perp(ndim)
    double precision, intent(out) :: b_seed(ndim),b_f(ndim),b_b(ndim)
    double precision, intent(out) :: length_f,length_b
    integer, intent(out) :: face_f,face_b,status
    double precision, intent(in), optional :: b_min

    type(trace_qperp_result) :: result

    if (present(b_min)) then
      call trace_field_qperp_single(seed,dL,max_steps,result,b_min)
    else
      call trace_field_qperp_single(seed,dL,max_steps,result)
    endif

    qperp=result%qperp
    logqperp=result%logqperp
    N2=result%N2
    bfactor=result%bfactor
    x_f=result%forward_endpoint
    x_b=result%backward_endpoint
    u_f_perp=result%u_forward_perp
    v_f_perp=result%v_forward_perp
    u_b_perp=result%u_backward_perp
    v_b_perp=result%v_backward_perp
    b_seed=result%B_seed
    b_f=result%forward_B
    b_b=result%backward_B
    length_f=result%forward_length
    length_b=result%backward_length
    face_f=result%forward_face
    face_b=result%backward_face
    status=result%status
  end subroutine trace_debug_qperp_single

  subroutine trace_init_qperp_result(seed,result)
    double precision, intent(in) :: seed(ndim)
    type(trace_qperp_result), intent(out) :: result

    result%seed=seed
    result%qperp=trace_debug_nan()
    result%logqperp=trace_debug_nan()
    result%N2=trace_debug_nan()
    result%bfactor=trace_debug_nan()
    result%q0=trace_debug_nan()
    result%logq0=trace_debug_nan()
    result%qperp0=trace_debug_nan()
    result%logqperp0=trace_debug_nan()
    result%N2_qperp0=trace_debug_nan()
    result%bfactor_qperp0=trace_debug_nan()
    result%forward_Bn_q0=trace_debug_nan()
    result%backward_Bn_q0=trace_debug_nan()
    result%B_seed=zero
    result%bhat_seed=zero
    result%forward_endpoint=seed
    result%backward_endpoint=seed
    result%forward_B=zero
    result%backward_B=zero
    result%forward_bhat=zero
    result%backward_bhat=zero
    result%forward_length=zero
    result%backward_length=zero
    result%u0=zero
    result%v0=zero
    result%u_forward_perp=zero
    result%v_forward_perp=zero
    result%u_backward_perp=zero
    result%v_backward_perp=zero
    result%forward_face=trace_face_none
    result%backward_face=trace_face_none
    result%forward_status=trace_status_active
    result%backward_status=trace_status_active
    result%status_q0=trace_status_active
    result%status_qperp0=trace_status_active
    result%status=trace_status_active
    result%valid=.false.
    result%valid_q0=.false.
    result%valid_qperp0=.false.
  end subroutine trace_init_qperp_result

  subroutine trace_qperp_single_core(seed,dL,max_steps,result,b_min)
    double precision, intent(in) :: seed(ndim),dL
    integer, intent(in) :: max_steps
    type(trace_qperp_result), intent(out) :: result
    double precision, intent(in), optional :: b_min

    double precision :: field_min,B3(3),Bseed_norm
    double precision :: u_end(ndim),v_end(ndim),u_face(ndim),v_face(ndim)
    integer :: igrid,status,indomain

    call trace_init_qperp_result(seed,result)

    if (npe/=1) then
      result%status=trace_status_mpi_unsupported
      return
    else if (ndim/=3 .or. .not.slab_uniform) then
      result%status=trace_status_unsupported_geometry
      return
    else if (dL<=zero .or. max_steps<=0) then
      result%status=trace_status_invalid_input
      return
    endif

    indomain=0
    {if (seed(^DB)>=xprobmin^DB .and. seed(^DB)<xprobmax^DB) indomain=indomain+1\}
    if (indomain/=ndim) then
      result%status=trace_status_seed_outside
      return
    endif

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)

    call trace_debug_locate_point(seed,igrid,status)
    if (status/=trace_status_active) then
      result%status=status
      return
    endif
    call sample_B_at_point(seed,igrid,field_min,B3,status)
    if (status/=trace_status_active) then
      result%status=status
      return
    endif
    result%B_seed=B3(1:ndim)
    Bseed_norm=dsqrt(sum(result%B_seed**2))
    if (Bseed_norm<=zero .or. Bseed_norm<field_min) then
      result%status=trace_status_weak_field
      return
    endif
    result%bhat_seed=result%B_seed/Bseed_norm
    call trace_make_perp_basis(B3/Bseed_norm,result%u0,result%v0,status)
    if (status/=trace_status_active) then
      result%status=status
      return
    endif

    call trace_debug_transport_tangents_to_boundary(seed,result%u0, &
         result%v0,dL,max_steps,.true.,result%forward_endpoint,u_end, &
         v_end,result%forward_bhat,result%forward_B,result%u_forward_perp, &
         result%v_forward_perp,u_face,v_face,result%forward_length, &
         result%forward_face,status,b_min=field_min)
    result%forward_status=status
    if (status/=trace_status_boundary) then
      result%status=status
      return
    endif

    call trace_debug_transport_tangents_to_boundary(seed,result%u0, &
         result%v0,dL,max_steps,.false.,result%backward_endpoint,u_end, &
         v_end,result%backward_bhat,result%backward_B,result%u_backward_perp, &
         result%v_backward_perp,u_face,v_face,result%backward_length, &
         result%backward_face,status,b_min=field_min)
    result%backward_status=status
    if (status/=trace_status_boundary) then
      result%status=status
      return
    endif

    call trace_qperp_compute_scalars(result,field_min)
  end subroutine trace_qperp_single_core

  subroutine trace_tangent_state_init(seed,u0,v0,igrid,seed_id,direction,state)
    double precision, intent(in) :: seed(ndim),u0(ndim),v0(ndim)
    integer, intent(in) :: igrid,seed_id,direction
    type(trace_tangent_state), intent(out) :: state

    state%seed_id=seed_id
    state%direction=direction
    state%active=.true.
    state%complete=.false.
    state%x=seed
    state%u=u0
    state%v=v0
    state%length=zero
    state%nstep=0
    state%igrid=igrid
    state%status=trace_status_active
    state%face_id=trace_face_none
    state%endpoint=seed
    state%endpoint_B=zero
    state%endpoint_bhat=zero
    state%u_perp=zero
    state%v_perp=zero
  end subroutine trace_tangent_state_init

  subroutine trace_tangent_trace_states_grouped(states,nstate,dL,max_steps, &
       threshold)
    integer, intent(in) :: nstate,max_steps
    type(trace_tangent_state), intent(inout) :: states(nstate)
    double precision, intent(in) :: dL,threshold

    logical, allocatable :: processed(:)
    integer :: istate,jstate,target_grid
    logical :: any_active

    allocate(processed(nstate))
    do
      any_active=.false.
      do istate=1,nstate
        if (states(istate)%active) then
          any_active=.true.
          exit
        endif
      enddo
      if (.not.any_active) exit

      processed=.false.
      do istate=1,nstate
        if (.not.states(istate)%active .or. processed(istate)) cycle
        target_grid=states(istate)%igrid
        do jstate=1,nstate
          if (.not.states(jstate)%active .or. processed(jstate)) cycle
          if (states(jstate)%igrid/=target_grid) cycle
          call trace_tangent_state_advance_in_grid(states(jstate), &
               target_grid,dL,max_steps,threshold)
          processed(jstate)=.true.
        enddo
      enddo
    enddo
    deallocate(processed)
  end subroutine trace_tangent_trace_states_grouped

  subroutine trace_tangent_state_advance_in_grid(state,igrid,dL,max_steps, &
       threshold)
    type(trace_tangent_state), intent(inout) :: state
    integer, intent(in) :: igrid,max_steps
    double precision, intent(in) :: dL,threshold

    double precision :: h,h_partial,partial_length
    double precision :: xprobe(ndim),xtrial(ndim),utrial(ndim),vtrial(ndim)
    double precision :: xhit(ndim),xold(ndim),uold(ndim),vold(ndim)
    double precision :: kx1(ndim),ku1(ndim),kv1(ndim)
    integer :: igrid_trial,point_domain
    logical :: hit_ok

    if (.not.state%active) return

    h=dL
    if (state%direction<0) h=-dL

    do while(state%active .and. state%igrid==igrid .and. &
         state%nstep<max_steps)
      xold=state%x
      uold=state%u
      vold=state%v

      call trace_tangent_rhs(xold,uold,vold,state%igrid,threshold,kx1,ku1, &
           kv1,state%status)
      if (state%status/=trace_status_active) then
        state%active=.false.
        state%complete=.true.
        return
      endif
      xprobe=xold+h*kx1

      point_domain=0
      {if (xprobe(^DB)>=xprobmin^DB .and. xprobe(^DB)<xprobmax^DB) point_domain=point_domain+1\}
      if (point_domain/=ndim) then
        call trace_intersect_domain(xold,xprobe,xhit,hit_ok,state%face_id)
        if (.not.hit_ok) then
          call trace_boundary_face_at_point(xold,state%face_id,hit_ok)
          state%status=trace_status_boundary
          state%active=.false.
          state%x=xold
          state%u=uold
          state%v=vold
          state%endpoint=xold
          if (hit_ok) call trace_tangent_state_finalize_boundary(state, &
               threshold)
          return
        endif
        partial_length=dsqrt(sum((xhit-xold)**2))
        if (partial_length<=100.d0*epsilon(one)*max(one,abs(h))) then
          state%x=xhit
          state%u=uold
          state%v=vold
          state%endpoint=xhit
          state%status=trace_status_boundary
          state%active=.false.
          call trace_tangent_state_finalize_boundary(state,threshold)
          return
        endif
        h_partial=sign(partial_length,h)
        call trace_tangent_rk2_trial_from_rhs(xold,uold,vold,state%igrid, &
             h_partial,threshold,kx1,ku1,kv1,xtrial,utrial,vtrial, &
             state%status)
        if (state%status/=trace_status_active) then
          state%active=.false.
          state%complete=.true.
          return
        endif
        state%x=xhit
        state%u=utrial
        state%v=vtrial
        state%length=state%length+partial_length
        state%endpoint=xhit
        state%status=trace_status_boundary
        state%active=.false.
        call trace_tangent_state_finalize_boundary(state,threshold)
        return
      endif

      call trace_tangent_rk2_trial_from_rhs(xold,uold,vold,state%igrid,h, &
           threshold,kx1,ku1,kv1,xtrial,utrial,vtrial,state%status)
      if (state%status/=trace_status_active) then
        state%active=.false.
        state%complete=.true.
        return
      endif

      point_domain=0
      {if (xtrial(^DB)>=xprobmin^DB .and. xtrial(^DB)<xprobmax^DB) point_domain=point_domain+1\}
      if (point_domain==ndim) then
        call trace_locate_point_with_hint(xtrial,state%igrid,igrid_trial, &
             state%status)
        if (state%status/=trace_status_active) then
          state%active=.false.
          state%complete=.true.
          return
        endif
        state%x=xtrial
        state%u=utrial
        state%v=vtrial
        state%igrid=igrid_trial
        state%length=state%length+abs(h)
        state%nstep=state%nstep+1
      else
        call trace_intersect_domain(xold,xtrial,xhit,hit_ok,state%face_id)
        if (.not.hit_ok) then
          call trace_boundary_face_at_point(xold,state%face_id,hit_ok)
          state%status=trace_status_boundary
          state%active=.false.
          state%x=xold
          state%u=uold
          state%v=vold
          state%endpoint=xold
          if (hit_ok) call trace_tangent_state_finalize_boundary(state, &
               threshold)
          return
        endif
        partial_length=dsqrt(sum((xhit-xold)**2))
        if (partial_length<=100.d0*epsilon(one)*max(one,abs(h))) then
          state%x=xhit
          state%u=uold
          state%v=vold
          state%endpoint=xhit
          state%status=trace_status_boundary
          state%active=.false.
          call trace_tangent_state_finalize_boundary(state,threshold)
          return
        endif
        h_partial=sign(partial_length,h)
        call trace_tangent_rk2_trial_from_rhs(xold,uold,vold,state%igrid, &
             h_partial,threshold,kx1,ku1,kv1,xtrial,utrial,vtrial, &
             state%status)
        if (state%status/=trace_status_active) then
          state%active=.false.
          state%complete=.true.
          return
        endif
        state%x=xhit
        state%u=utrial
        state%v=vtrial
        state%length=state%length+partial_length
        state%endpoint=xhit
        state%status=trace_status_boundary
        state%active=.false.
        call trace_tangent_state_finalize_boundary(state,threshold)
        return
      endif
    enddo

    if (state%active .and. state%nstep>=max_steps) then
      state%status=trace_status_max_steps
      state%active=.false.
      state%complete=.true.
    endif
  end subroutine trace_tangent_state_advance_in_grid

  subroutine trace_tangent_state_finalize_boundary(state,threshold)
    type(trace_tangent_state), intent(inout) :: state
    double precision, intent(in) :: threshold

    if (state%status/=trace_status_boundary) return

    call trace_endpoint_B_bhat(state%endpoint,state%face_id,state%igrid, &
         threshold,state%endpoint_B,state%endpoint_bhat,state%status)
    if (state%status/=trace_status_boundary) return

    call trace_project_to_perp_bhat(state%u,state%endpoint_bhat, &
         state%u_perp)
    call trace_project_to_perp_bhat(state%v,state%endpoint_bhat, &
         state%v_perp)
    state%complete=.true.
  end subroutine trace_tangent_state_finalize_boundary

  subroutine trace_qperp_prepare_seed_result(seed,field_min,result,igrid, &
       status)
    double precision, intent(in) :: seed(ndim),field_min
    type(trace_qperp_result), intent(out) :: result
    integer, intent(out) :: igrid,status

    double precision :: B3(3),Bseed_norm
    integer :: indomain

    call trace_init_qperp_result(seed,result)
    igrid=-1
    status=trace_status_active

    indomain=0
    {if (seed(^DB)>=xprobmin^DB .and. seed(^DB)<xprobmax^DB) indomain=indomain+1\}
    if (indomain/=ndim) then
      status=trace_status_seed_outside
      result%status=status
      return
    endif

    call trace_debug_locate_point(seed,igrid,status)
    if (status/=trace_status_active) then
      result%status=status
      return
    endif

    call sample_B_at_point(seed,igrid,field_min,B3,status)
    if (status/=trace_status_active) then
      result%status=status
      return
    endif
    result%B_seed=B3(1:ndim)
    Bseed_norm=dsqrt(sum(result%B_seed**2))
    if (Bseed_norm<=zero .or. Bseed_norm<field_min) then
      status=trace_status_weak_field
      result%status=status
      return
    endif

    result%bhat_seed=result%B_seed/Bseed_norm
    call trace_make_perp_basis(B3/Bseed_norm,result%u0,result%v0,status)
    if (status/=trace_status_active) result%status=status
  end subroutine trace_qperp_prepare_seed_result

  subroutine trace_qperp_compute_scalars(result,field_min)
    type(trace_qperp_result), intent(inout) :: result
    double precision, intent(in) :: field_min

    double precision :: Bseed_norm,Bf_norm,Bb_norm,dot_f,dot_b,qtmp

    Bseed_norm=dsqrt(sum(result%B_seed**2))
    Bf_norm=dsqrt(sum(result%forward_B**2))
    Bb_norm=dsqrt(sum(result%backward_B**2))
    if (Bseed_norm<=zero .or. Bf_norm<=zero .or. Bb_norm<=zero .or. &
         Bseed_norm<field_min .or. Bf_norm<field_min .or. &
         Bb_norm<field_min) then
      result%status=trace_status_weak_field
      return
    endif

    dot_f=sum(result%u_forward_perp*result%v_forward_perp)
    dot_b=sum(result%u_backward_perp*result%v_backward_perp)
    result%N2=sum(result%u_forward_perp**2)* &
         sum(result%v_backward_perp**2)+ &
         sum(result%u_backward_perp**2)* &
         sum(result%v_forward_perp**2)-2.d0*dot_f*dot_b
    result%bfactor=abs(Bf_norm*Bb_norm)/(Bseed_norm**2)
    qtmp=result%N2*result%bfactor
    if (qtmp>zero .and. ieee_is_finite(qtmp)) then
      result%qperp=qtmp
      result%logqperp=dlog10(qtmp)
      result%valid=.true.
      result%status=trace_status_boundary
      result%qperp0=qtmp
      result%logqperp0=result%logqperp
      result%valid_qperp0=.true.
      result%status_qperp0=trace_status_boundary
    else
      result%N2=trace_debug_nan()
      result%bfactor=trace_debug_nan()
      result%qperp=trace_debug_nan()
      result%logqperp=trace_debug_nan()
      result%status=trace_status_invalid_input
      result%valid=.false.
      result%qperp0=trace_debug_nan()
      result%logqperp0=trace_debug_nan()
      result%valid_qperp0=.false.
      result%status_qperp0=trace_status_invalid_input
    endif
    call trace_qperp_compute_q0_scalars(result,field_min)
  end subroutine trace_qperp_compute_scalars

  subroutine trace_qperp_compute_q0_scalars(result,field_min)
    type(trace_qperp_result), intent(inout) :: result
    double precision, intent(in) :: field_min

    double precision :: Bseed_norm,Bnf,Bnb,N2bnd,qtmp
    double precision :: uf_face(ndim),vf_face(ndim)
    double precision :: ub_face(ndim),vb_face(ndim)
    integer :: status_f,status_b

    result%q0=trace_debug_nan()
    result%logq0=trace_debug_nan()
    result%N2_qperp0=trace_debug_nan()
    result%bfactor_qperp0=trace_debug_nan()
    result%forward_Bn_q0=trace_debug_nan()
    result%backward_Bn_q0=trace_debug_nan()
    result%valid_q0=.false.

    if (result%forward_status/=trace_status_boundary) then
      result%status_q0=result%forward_status
      return
    endif
    if (result%backward_status/=trace_status_boundary) then
      result%status_q0=result%backward_status
      return
    endif

    Bseed_norm=dsqrt(sum(result%B_seed**2))
    if (Bseed_norm<=zero .or. Bseed_norm<field_min) then
      result%status_q0=trace_status_weak_field
      return
    endif

    call trace_qperp_project_to_boundary_face(result%u_forward_perp, &
         result%forward_B,result%forward_face,uf_face,Bnf,status_f)
    call trace_qperp_project_to_boundary_face(result%v_forward_perp, &
         result%forward_B,result%forward_face,vf_face,Bnf,status_f)
    call trace_qperp_project_to_boundary_face(result%u_backward_perp, &
         result%backward_B,result%backward_face,ub_face,Bnb,status_b)
    call trace_qperp_project_to_boundary_face(result%v_backward_perp, &
         result%backward_B,result%backward_face,vb_face,Bnb,status_b)
    if (status_f/=trace_status_active) then
      result%status_q0=status_f
      return
    endif
    if (status_b/=trace_status_active) then
      result%status_q0=status_b
      return
    endif
    if (abs(Bnf)<field_min .or. abs(Bnb)<field_min) then
      result%status_q0=trace_status_weak_field
      return
    endif

    N2bnd=sum(uf_face**2)*sum(vb_face**2)+ &
         sum(ub_face**2)*sum(vf_face**2)- &
         2.d0*sum(uf_face*vf_face)*sum(ub_face*vb_face)
    result%N2_qperp0=abs(N2bnd)
    result%bfactor_qperp0=abs(Bnf*Bnb)/(Bseed_norm**2)
    result%forward_Bn_q0=Bnf
    result%backward_Bn_q0=Bnb
    qtmp=result%N2_qperp0*result%bfactor_qperp0
    if (qtmp>zero .and. ieee_is_finite(qtmp)) then
      result%q0=qtmp
      result%logq0=dlog10(qtmp)
      result%valid_q0=.true.
      result%status_q0=trace_status_boundary
    else
      result%q0=trace_debug_nan()
      result%logq0=trace_debug_nan()
      result%valid_q0=.false.
      result%status_q0=trace_status_invalid_input
    endif
  end subroutine trace_qperp_compute_q0_scalars

  subroutine trace_qperp_project_to_boundary_face(vec,B,face_id,vec_face, &
       Bn,status)
    double precision, intent(in) :: vec(ndim),B(ndim)
    integer, intent(in) :: face_id
    double precision, intent(out) :: vec_face(ndim),Bn
    integer, intent(out) :: status

    double precision :: normal(ndim),vecn
    logical :: ok

    vec_face=zero
    Bn=zero
    status=trace_status_invalid_input
    call trace_face_normal(face_id,normal,ok)
    if (.not.ok) return
    Bn=sum(B*normal)
    if (Bn==zero) then
      status=trace_status_weak_field
      return
    endif
    vecn=sum(vec*normal)
    vec_face=vec-vecn/Bn*B
    status=trace_status_active
  end subroutine trace_qperp_project_to_boundary_face

  subroutine trace_qperp_finalize_from_states(forward_state,backward_state, &
       result,field_min)
    type(trace_tangent_state), intent(in) :: forward_state,backward_state
    type(trace_qperp_result), intent(inout) :: result
    double precision, intent(in) :: field_min

    result%forward_endpoint=forward_state%endpoint
    result%forward_B=forward_state%endpoint_B
    result%forward_bhat=forward_state%endpoint_bhat
    result%forward_length=forward_state%length
    result%forward_face=forward_state%face_id
    result%forward_status=forward_state%status
    result%u_forward_perp=forward_state%u_perp
    result%v_forward_perp=forward_state%v_perp

    result%backward_endpoint=backward_state%endpoint
    result%backward_B=backward_state%endpoint_B
    result%backward_bhat=backward_state%endpoint_bhat
    result%backward_length=backward_state%length
    result%backward_face=backward_state%face_id
    result%backward_status=backward_state%status
    result%u_backward_perp=backward_state%u_perp
    result%v_backward_perp=backward_state%v_perp

    if (forward_state%status/=trace_status_boundary) then
      result%status=forward_state%status
      return
    endif
    if (backward_state%status/=trace_status_boundary) then
      result%status=backward_state%status
      return
    endif

    call trace_qperp_compute_scalars(result,field_min)
  end subroutine trace_qperp_finalize_from_states

  subroutine trace_debug_locate_point(x,igrid,status)
    use mod_particle_base, only: find_particle_ipe

    double precision, intent(in) :: x(ndim)
    integer, intent(out) :: igrid,status

    double precision :: x3d(3)
    integer :: indomain,ipe,j

    igrid=-1
    status=trace_status_active
    indomain=0
    {if (x(^DB)>=xprobmin^DB .and. x(^DB)<xprobmax^DB) indomain=indomain+1\}
    if (indomain/=ndim) then
      status=trace_status_boundary
      return
    endif

    x3d=zero
    do j=1,ndim
      x3d(j)=x(j)
    enddo
    call find_particle_ipe(x3d,igrid,ipe)
    if (igrid<0 .or. ipe/=mype) then
      status=trace_status_out_of_domain
    endif
  end subroutine trace_debug_locate_point

  subroutine trace_locate_point_with_hint(x,igrid_hint,igrid,status)
    double precision, intent(in) :: x(ndim)
    integer, intent(in) :: igrid_hint
    integer, intent(out) :: igrid,status

    double precision :: xbmin^D,xbmax^D
    integer :: inblock

    igrid=-1
    status=trace_status_active
    if (igrid_hint>=0) then
      ^D&xbmin^D=rnode(rpxmin^D_,igrid_hint);
      ^D&xbmax^D=rnode(rpxmax^D_,igrid_hint);
      inblock=0
      {if (x(^DB)>xbmin^DB .and. x(^DB)<xbmax^DB) inblock=inblock+1\}
      if (inblock==ndim) then
        igrid=igrid_hint
        return
      endif
    endif

    call trace_debug_locate_point(x,igrid,status)
  end subroutine trace_locate_point_with_hint

  subroutine trace_tangent_rhs(x,u,v,igrid,threshold,kx,ku,kv,status)
    double precision, intent(in) :: x(ndim),u(ndim),v(ndim),threshold
    integer, intent(in) :: igrid
    double precision, intent(out) :: kx(ndim),ku(ndim),kv(ndim)
    integer, intent(out) :: status

    double precision :: bhat(3),grad_bhat(3,3)
    double precision :: u3(3),v3(3),ku3(3),kv3(3)

    kx=zero
    ku=zero
    kv=zero
    status=trace_status_unsupported_geometry
    if (ndim/=3) return

    {^IFTHREED
    call sample_bhat_gradbhat_at_point(x,igrid,threshold,bhat,grad_bhat,status)
    if (status/=trace_status_active) return

    u3=zero
    v3=zero
    u3(1:ndim)=u
    v3(1:ndim)=v
    ku3=matmul(grad_bhat,u3)
    kv3=matmul(grad_bhat,v3)
    kx=bhat(1:ndim)
    ku=ku3(1:ndim)
    kv=kv3(1:ndim)
    }
  end subroutine trace_tangent_rhs

  subroutine trace_advance_tangent_state_rk2(x,u,v,igrid,h,threshold,status)
    double precision, intent(inout) :: x(ndim),u(ndim),v(ndim)
    integer, intent(inout) :: igrid
    double precision, intent(in) :: h,threshold
    integer, intent(out) :: status

    double precision :: xnew(ndim),unew(ndim),vnew(ndim)
    integer :: igrid_mid,igrid_new

    call trace_tangent_rk2_trial(x,u,v,igrid,h,threshold,xnew,unew,vnew, &
         status)
    if (status/=trace_status_active) return

    call trace_debug_locate_point(xnew,igrid_new,status)
    if (status/=trace_status_active) return

    x=xnew
    u=unew
    v=vnew
    igrid=igrid_new
  end subroutine trace_advance_tangent_state_rk2

  subroutine trace_tangent_rk2_trial(x,u,v,igrid,h,threshold,xnew,unew, &
       vnew,status)
    double precision, intent(in) :: x(ndim),u(ndim),v(ndim),h,threshold
    integer, intent(in) :: igrid
    double precision, intent(out) :: xnew(ndim),unew(ndim),vnew(ndim)
    integer, intent(out) :: status

    double precision :: kx1(ndim),ku1(ndim),kv1(ndim)

    call trace_tangent_rhs(x,u,v,igrid,threshold,kx1,ku1,kv1,status)
    if (status/=trace_status_active) return

    call trace_tangent_rk2_trial_from_rhs(x,u,v,igrid,h,threshold, &
         kx1,ku1,kv1,xnew,unew,vnew,status)
  end subroutine trace_tangent_rk2_trial

  subroutine trace_tangent_rk2_trial_from_rhs(x,u,v,igrid,h,threshold, &
       kx1,ku1,kv1,xnew,unew,vnew,status)
    double precision, intent(in) :: x(ndim),u(ndim),v(ndim),h,threshold
    double precision, intent(in) :: kx1(ndim),ku1(ndim),kv1(ndim)
    integer, intent(in) :: igrid
    double precision, intent(out) :: xnew(ndim),unew(ndim),vnew(ndim)
    integer, intent(out) :: status

    double precision :: kx2(ndim),ku2(ndim),kv2(ndim)
    double precision :: xmid(ndim),umid(ndim),vmid(ndim)
    integer :: igrid_mid

    xmid=x+half*h*kx1
    umid=u+half*h*ku1
    vmid=v+half*h*kv1
    call trace_locate_point_with_hint(xmid,igrid,igrid_mid,status)
    if (status/=trace_status_active) return

    call trace_tangent_rhs(xmid,umid,vmid,igrid_mid,threshold,kx2,ku2,kv2, &
         status)
    if (status/=trace_status_active) return

    xnew=x+h*kx2
    unew=u+h*ku2
    vnew=v+h*kv2
  end subroutine trace_tangent_rk2_trial_from_rhs

  subroutine trace_project_to_perp_bhat(vec,bhat,vec_perp)
    double precision, intent(in) :: vec(ndim),bhat(ndim)
    double precision, intent(out) :: vec_perp(ndim)

    vec_perp=vec-sum(vec*bhat)*bhat
  end subroutine trace_project_to_perp_bhat

  subroutine trace_endpoint_B_bhat(x,face_id,igrid,threshold,B,bhat,status)
    double precision, intent(in) :: x(ndim),threshold
    integer, intent(in) :: face_id,igrid
    double precision, intent(out) :: B(ndim),bhat(ndim)
    integer, intent(out) :: status

    double precision :: B3(3),bhat3(3)
    integer :: sample_igrid

    B=zero
    bhat=zero
    call trace_locate_face_limit_grid(x,face_id,igrid,sample_igrid,status)
    if (status/=trace_status_active) return
    call trace_sample_B_on_domain_face_limit(x,face_id,sample_igrid,threshold, &
         B3,bhat3,status)
    if (status/=trace_status_active) return
    B=B3(1:ndim)
    bhat=bhat3(1:ndim)
    status=trace_status_boundary
  end subroutine trace_endpoint_B_bhat

  subroutine trace_locate_face_limit_grid(x,face_id,igrid_in,igrid_out, &
       status)
    double precision, intent(in) :: x(ndim)
    integer, intent(in) :: face_id,igrid_in
    integer, intent(out) :: igrid_out,status

    double precision :: xprobe(ndim),dxn,xface
    integer :: normal_dim,side

    igrid_out=-1
    status=trace_status_bad_face_limit_sample
    if (ndim/=3 .or. .not.slab_uniform .or. igrid_in<0) then
      status=trace_status_unsupported_geometry
      return
    endif
    if (.not.trace_face_is_boundary(face_id)) return

    {^IFTHREED
    select case(face_id)
    case(trace_face_xmin)
      normal_dim=1
      side=-1
      xface=xprobmin1
      dxn=rnode(rpdx1_,igrid_in)
    case(trace_face_xmax)
      normal_dim=1
      side=1
      xface=xprobmax1
      dxn=rnode(rpdx1_,igrid_in)
    case(trace_face_ymin)
      normal_dim=2
      side=-1
      xface=xprobmin2
      dxn=rnode(rpdx2_,igrid_in)
    case(trace_face_ymax)
      normal_dim=2
      side=1
      xface=xprobmax2
      dxn=rnode(rpdx2_,igrid_in)
    case(trace_face_zmin)
      normal_dim=3
      side=-1
      xface=xprobmin3
      dxn=rnode(rpdx3_,igrid_in)
    case(trace_face_zmax)
      normal_dim=3
      side=1
      xface=xprobmax3
      dxn=rnode(rpdx3_,igrid_in)
    case default
      return
    end select

    xprobe=x
    xprobe(normal_dim)=xface-half*side*dxn
    call trace_debug_locate_point(xprobe,igrid_out,status)
    }
  end subroutine trace_locate_face_limit_grid

  subroutine trace_sample_B_on_domain_face_limit(xhit,face_id,igrid, &
       threshold,B,bhat,status)
    ! Sample total B at a rectangular domain face from the inside-side limit.
    ! The normal direction is linearly extrapolated from the two nearest
    ! interior cell centers; tangential directions use bilinear interpolation.
    double precision, intent(in) :: xhit(ndim),threshold
    integer, intent(in) :: face_id,igrid
    double precision, intent(out) :: B(3),bhat(3)
    integer, intent(out) :: status

    double precision :: Bnear(3),Binner(3)
    double precision :: dxb(3),xface,xnear,xinner
    double precision :: xnear_pt(ndim),xinner_pt(ndim)
    double precision :: Bnorm
    integer :: ixO^L
    integer :: normal_dim,side
    integer :: igrid_near,igrid_inner

    B=zero
    bhat=zero
    status=trace_status_bad_face_limit_sample
    if (ndim/=3 .or. .not.slab_uniform .or. igrid<0) then
      status=trace_status_unsupported_geometry
      return
    endif
    if (.not.trace_face_is_boundary(face_id)) return

    {^IFTHREED
    select case(face_id)
    case(trace_face_xmin)
      normal_dim=1
      side=-1
      xface=xprobmin1
    case(trace_face_xmax)
      normal_dim=1
      side=1
      xface=xprobmax1
    case(trace_face_ymin)
      normal_dim=2
      side=-1
      xface=xprobmin2
    case(trace_face_ymax)
      normal_dim=2
      side=1
      xface=xprobmax2
    case(trace_face_zmin)
      normal_dim=3
      side=-1
      xface=xprobmin3
    case(trace_face_zmax)
      normal_dim=3
      side=1
      xface=xprobmax3
    case default
      return
    end select

    ixO^L=ixM^LL;
    dxb=(/ rnode(rpdx1_,igrid),rnode(rpdx2_,igrid), &
         rnode(rpdx3_,igrid) /)

    xnear=xface-half*side*dxb(normal_dim)
    xinner=xface-1.5d0*side*dxb(normal_dim)
    xnear_pt=xhit
    xinner_pt=xhit
    xnear_pt(normal_dim)=xnear
    xinner_pt(normal_dim)=xinner

    call trace_debug_locate_point(xnear_pt,igrid_near,status)
    if (status/=trace_status_active) return
    call sample_B_at_point(xnear_pt,igrid_near,threshold,Bnear,status)
    if (status/=trace_status_active) return
    call trace_debug_locate_point(xinner_pt,igrid_inner,status)
    if (status/=trace_status_active) return
    call sample_B_at_point(xinner_pt,igrid_inner,threshold,Binner,status)
    if (status/=trace_status_active) return

    B=Bnear+(Bnear-Binner)/(xnear-xinner)*(xface-xnear)
    Bnorm=dsqrt(sum(B**2))
    if (Bnorm<=zero .or. Bnorm<threshold) then
      status=trace_status_weak_field
      return
    endif
    bhat=B/Bnorm
    status=trace_status_active
    }
  end subroutine trace_sample_B_on_domain_face_limit

  subroutine trace_face_normal(face_id,normal,ok)
    integer, intent(in) :: face_id
    double precision, intent(out) :: normal(ndim)
    logical, intent(out) :: ok

    normal=zero
    ok=.true.
    select case(face_id)
    case(trace_face_xmin)
      normal(1)=-one
    case(trace_face_xmax)
      normal(1)=one
    {^IFONED
    case default
      ok=.false.
    }
    {^IFTWOD
    case(trace_face_ymin)
      normal(2)=-one
    case(trace_face_ymax)
      normal(2)=one
    case default
      ok=.false.
    }
    {^IFTHREED
    case(trace_face_ymin)
      normal(2)=-one
    case(trace_face_ymax)
      normal(2)=one
    case(trace_face_zmin)
      normal(3)=-one
    case(trace_face_zmax)
      normal(3)=one
    case default
      ok=.false.
    }
    end select
  end subroutine trace_face_normal

  subroutine trace_project_to_boundary_face(vec,bhat,face_id,vec_face)
    double precision, intent(in) :: vec(ndim),bhat(ndim)
    integer, intent(in) :: face_id
    double precision, intent(out) :: vec_face(ndim)

    double precision :: normal(ndim),bn
    logical :: ok

    vec_face=zero
    call trace_face_normal(face_id,normal,ok)
    if (.not.ok) return
    bn=sum(bhat*normal)
    if (abs(bn)<=smalldouble) return
    vec_face=vec-sum(vec*normal)/bn*bhat
  end subroutine trace_project_to_boundary_face

  subroutine trace_make_perp_basis(bhat,u0,v0,status)
    double precision, intent(in) :: bhat(3)
    double precision, intent(out) :: u0(ndim),v0(ndim)
    integer, intent(out) :: status

    double precision :: ref(3),u3(3),v3(3),b3(3),unorm,bnorm
    integer :: iref

    u0=zero
    v0=zero
    status=trace_status_active
    if (ndim/=3) then
      status=trace_status_unsupported_geometry
      return
    endif

    b3=bhat
    bnorm=dsqrt(sum(b3**2))
    if (bnorm<=zero) then
      status=trace_status_weak_field
      return
    endif
    b3=b3/bnorm

    iref=1
    if (abs(b3(2))<abs(b3(iref))) iref=2
    if (abs(b3(3))<abs(b3(iref))) iref=3
    ref=zero
    ref(iref)=one

    u3=ref-sum(ref*b3)*b3
    unorm=dsqrt(sum(u3**2))
    if (unorm<=smalldouble) then
      status=trace_status_weak_field
      return
    endif
    u3=u3/unorm

    v3(1)=b3(2)*u3(3)-b3(3)*u3(2)
    v3(2)=b3(3)*u3(1)-b3(1)*u3(3)
    v3(3)=b3(1)*u3(2)-b3(2)*u3(1)
    v3=v3/dsqrt(sum(v3**2))

    u0=u3(1:ndim)
    v0=v3(1:ndim)
  end subroutine trace_make_perp_basis

  double precision function trace_debug_nan()
    trace_debug_nan=ieee_value(0.d0,ieee_quiet_nan)
  end function trace_debug_nan

  subroutine trace_summary_seed(seed,dL,max_steps,result,b_min)
    double precision, intent(in) :: seed(ndim),dL
    integer, intent(in) :: max_steps
    type(trace_length_result), intent(out) :: result
    double precision, intent(in), optional :: b_min

    double precision :: field_min
    integer :: common_status,igrid

    call trace_summary_validate(dL,max_steps,common_status)
    if (common_status/=trace_status_active) then
      call trace_summary_init_result(seed,result,common_status)
      return
    endif

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)

    call trace_summary_locate_seed(seed,result,igrid)
    if (result%forward_status/=trace_status_active) return
    call trace_summary_trace_seed(igrid,dL,max_steps,field_min,result)
  end subroutine trace_summary_seed

  subroutine trace_summary_twist_seed(seed,dL,max_steps,result,b_min)
    double precision, intent(in) :: seed(ndim),dL
    integer, intent(in) :: max_steps
    type(trace_twist_result), intent(out) :: result
    double precision, intent(in), optional :: b_min

    double precision :: field_min
    integer :: common_status,igrid

    call trace_summary_validate(dL,max_steps,common_status)
    if (common_status==trace_status_active .and. ndim/=3) then
      common_status=trace_status_unsupported_geometry
    endif
    if (common_status/=trace_status_active) then
      call trace_summary_init_twist_result(seed,result,common_status)
      return
    endif

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)

    call trace_summary_locate_seed(seed,result%line,igrid)
    result%forward_twist=zero
    result%backward_twist=zero
    result%total_twist=zero
    if (result%line%forward_status/=trace_status_active) return
    call trace_summary_trace_twist_seed(igrid,dL,max_steps,field_min,result)
  end subroutine trace_summary_twist_seed

  subroutine trace_summary_mapping_seed(seed,dL,max_steps,result,b_min, &
       source_normal)
    double precision, intent(in) :: seed(ndim),dL
    integer, intent(in) :: max_steps
    type(trace_mapping_result), intent(out) :: result
    double precision, intent(in), optional :: b_min
    double precision, intent(in), optional :: source_normal(3)

    type(trace_summary_state) :: forward_state,backward_state
    type(trace_length_result) :: located
    double precision :: field_min,normal(3)
    integer :: common_status,igrid
    logical :: have_normal

    call trace_summary_validate(dL,max_steps,common_status)
    if (common_status==trace_status_active .and. ndim/=3) then
      common_status=trace_status_unsupported_geometry
    endif
    if (common_status/=trace_status_active) then
      call trace_summary_init_mapping_result(seed,result,common_status)
      return
    endif

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)
    normal=zero
    have_normal=present(source_normal)
    if (have_normal) normal=source_normal

    call trace_summary_locate_seed(seed,located,igrid)
    call trace_summary_init_mapping_result(seed,result,located%forward_status)
    if (located%forward_status/=trace_status_active) return

    call trace_summary_init_state(seed,igrid,.true.,1,.false.,forward_state)
    call trace_summary_trace_state_to_end(forward_state,dL,max_steps,field_min)
    call trace_summary_init_state(seed,igrid,.false.,1,.false.,backward_state)
    call trace_summary_trace_state_to_end(backward_state,dL,max_steps,field_min)
    call trace_summary_fill_mapping(seed,igrid,forward_state,backward_state, &
         field_min,have_normal,normal,result)
  end subroutine trace_summary_mapping_seed

  subroutine trace_summary_multi(seeds,nseed,dL,max_steps,results,field_min)
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL,field_min
    type(trace_length_result), intent(out) :: results(nseed)

    type(trace_summary_state), allocatable :: states(:)
    double precision :: seed_local(ndim)
    integer :: common_status,iseed,igrid,iforward,ibackward

    if (nseed<=0) return

    call trace_summary_validate(dL,max_steps,common_status)
    if (common_status/=trace_status_active) then
      do iseed=1,nseed
        seed_local=seeds(iseed,:)
        call trace_summary_init_result(seed_local,results(iseed),common_status)
      enddo
      return
    endif

    allocate(states(2*nseed))
    do iseed=1,nseed
      seed_local=seeds(iseed,:)
      call trace_summary_locate_seed(seed_local,results(iseed),igrid)
      iforward=2*iseed-1
      ibackward=2*iseed
      call trace_summary_init_state(seed_local,igrid,.true.,iseed,.false., &
           states(iforward))
      call trace_summary_init_state(seed_local,igrid,.false.,iseed,.false., &
           states(ibackward))
      if (results(iseed)%forward_status/=trace_status_active) then
        states(iforward)%status=results(iseed)%forward_status
        states(iforward)%active=.false.
        states(ibackward)%status=results(iseed)%backward_status
        states(ibackward)%active=.false.
      endif
    enddo

    call trace_summary_trace_states_grouped(states,2*nseed,dL,max_steps, &
         field_min)

    do iseed=1,nseed
      iforward=2*iseed-1
      ibackward=2*iseed
      results(iseed)%forward_footpoint=states(iforward)%footpoint
      results(iseed)%forward_length=states(iforward)%length
      results(iseed)%forward_nstep=states(iforward)%nstep
      results(iseed)%forward_status=states(iforward)%status
      results(iseed)%backward_footpoint=states(ibackward)%footpoint
      results(iseed)%backward_length=states(ibackward)%length
      results(iseed)%backward_nstep=states(ibackward)%nstep
      results(iseed)%backward_status=states(ibackward)%status
      results(iseed)%total_length=results(iseed)%forward_length &
           +results(iseed)%backward_length
    enddo
    deallocate(states)
  end subroutine trace_summary_multi

  subroutine trace_summary_twist_multi(seeds,nseed,dL,max_steps,results, &
       field_min)
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL,field_min
    type(trace_twist_result), intent(out) :: results(nseed)

    type(trace_summary_state), allocatable :: states(:)
    double precision :: seed_local(ndim)
    integer :: common_status,iseed,igrid,iforward,ibackward

    if (nseed<=0) return

    call trace_summary_validate(dL,max_steps,common_status)
    if (common_status==trace_status_active .and. ndim/=3) then
      common_status=trace_status_unsupported_geometry
    endif
    if (common_status/=trace_status_active) then
      do iseed=1,nseed
        seed_local=seeds(iseed,:)
        call trace_summary_init_twist_result(seed_local,results(iseed), &
             common_status)
      enddo
      return
    endif

    allocate(states(2*nseed))
    do iseed=1,nseed
      seed_local=seeds(iseed,:)
      call trace_summary_locate_seed(seed_local,results(iseed)%line,igrid)
      results(iseed)%forward_twist=zero
      results(iseed)%backward_twist=zero
      results(iseed)%total_twist=zero
      iforward=2*iseed-1
      ibackward=2*iseed
      call trace_summary_init_state(seed_local,igrid,.true.,iseed,.true., &
           states(iforward))
      call trace_summary_init_state(seed_local,igrid,.false.,iseed,.true., &
           states(ibackward))
      if (results(iseed)%line%forward_status/=trace_status_active) then
        states(iforward)%status=results(iseed)%line%forward_status
        states(iforward)%active=.false.
        states(ibackward)%status=results(iseed)%line%backward_status
        states(ibackward)%active=.false.
      endif
    enddo

    call trace_summary_trace_states_grouped(states,2*nseed,dL,max_steps, &
         field_min)

    do iseed=1,nseed
      iforward=2*iseed-1
      ibackward=2*iseed
      results(iseed)%line%forward_footpoint=states(iforward)%footpoint
      results(iseed)%line%forward_length=states(iforward)%length
      results(iseed)%line%forward_nstep=states(iforward)%nstep
      results(iseed)%line%forward_status=states(iforward)%status
      results(iseed)%line%backward_footpoint=states(ibackward)%footpoint
      results(iseed)%line%backward_length=states(ibackward)%length
      results(iseed)%line%backward_nstep=states(ibackward)%nstep
      results(iseed)%line%backward_status=states(ibackward)%status
      results(iseed)%line%total_length=results(iseed)%line%forward_length &
           +results(iseed)%line%backward_length
      results(iseed)%forward_twist=states(iforward)%twist
      results(iseed)%backward_twist=states(ibackward)%twist
      results(iseed)%total_twist=results(iseed)%forward_twist &
           +results(iseed)%backward_twist
    enddo
    deallocate(states)
  end subroutine trace_summary_twist_multi

  subroutine trace_summary_mapping_multi(seeds,nseed,dL,max_steps,results, &
       field_min,have_source_normal,source_normal)
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL,field_min
    logical, intent(in) :: have_source_normal
    double precision, intent(in) :: source_normal(3)
    type(trace_mapping_result), intent(out) :: results(nseed)

    type(trace_summary_state), allocatable :: states(:)
    integer, allocatable :: source_igrids(:)
    type(trace_length_result) :: located
    double precision :: seed_local(ndim)
    integer :: common_status,iseed,igrid,iforward,ibackward

    if (nseed<=0) return

    call trace_summary_validate(dL,max_steps,common_status)
    if (common_status==trace_status_active .and. ndim/=3) then
      common_status=trace_status_unsupported_geometry
    endif
    if (common_status/=trace_status_active) then
      do iseed=1,nseed
        seed_local=seeds(iseed,:)
        call trace_summary_init_mapping_result(seed_local,results(iseed), &
             common_status)
      enddo
      return
    endif

    allocate(states(2*nseed))
    allocate(source_igrids(nseed))
    source_igrids=-1
    do iseed=1,nseed
      seed_local=seeds(iseed,:)
      call trace_summary_locate_seed(seed_local,located,igrid)
      source_igrids(iseed)=igrid
      call trace_summary_init_mapping_result(seed_local,results(iseed), &
           located%forward_status)
      iforward=2*iseed-1
      ibackward=2*iseed
      call trace_summary_init_state(seed_local,igrid,.true.,iseed,.false., &
           states(iforward))
      call trace_summary_init_state(seed_local,igrid,.false.,iseed,.false., &
           states(ibackward))
      if (located%forward_status/=trace_status_active) then
        states(iforward)%status=located%forward_status
        states(iforward)%active=.false.
        states(ibackward)%status=located%backward_status
        states(ibackward)%active=.false.
      endif
    enddo

    call trace_summary_trace_states_grouped(states,2*nseed,dL,max_steps, &
         field_min)

    do iseed=1,nseed
      iforward=2*iseed-1
      ibackward=2*iseed
      if (source_igrids(iseed)>=0) then
        call trace_summary_fill_mapping(results(iseed)%seed, &
             source_igrids(iseed),states(iforward),states(ibackward), &
             field_min,have_source_normal,source_normal,results(iseed))
      endif
    enddo
    deallocate(source_igrids,states)
  end subroutine trace_summary_mapping_multi

  subroutine trace_summary_topology_multi(seeds,nseed,dL,max_steps,results, &
       field_min,need_twist,need_mapping,have_source_normal,source_normal)
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL,field_min
    logical, intent(in) :: need_twist,need_mapping,have_source_normal
    double precision, intent(in) :: source_normal(3)
    type(trace_topology_result), intent(out) :: results(nseed)

    type(trace_summary_state), allocatable :: states(:)
    type(trace_summary_state), allocatable :: map_states(:)
    integer, allocatable :: source_igrids(:)
    type(trace_length_result) :: located
    type(trace_mapping_result) :: mapping
    double precision :: seed_local(ndim)
    integer :: common_status,iseed,igrid,iforward,ibackward,idim

    if (nseed<=0) return

    call trace_summary_validate(dL,max_steps,common_status)
    if (common_status==trace_status_active .and. &
         (need_twist .or. need_mapping) .and. ndim/=3) then
      common_status=trace_status_unsupported_geometry
    endif
    if (common_status/=trace_status_active) then
      do iseed=1,nseed
        do idim=1,ndim
          seed_local(idim)=seeds(iseed,idim)
        enddo
        call trace_summary_init_topology_result(seed_local,results(iseed), &
             common_status,need_twist,need_mapping)
      enddo
      return
    endif

    allocate(states(2*nseed))
    ! Mapping summaries need only B/Bn; keep them independent of
    ! curl-stencil failures that can stop twist accumulation near boundaries.
    if (need_mapping .and. need_twist) allocate(map_states(2*nseed))
    allocate(source_igrids(nseed))
    source_igrids=-1
    do iseed=1,nseed
      do idim=1,ndim
        seed_local(idim)=seeds(iseed,idim)
      enddo
      call trace_summary_locate_seed(seed_local,located,igrid)
      source_igrids(iseed)=igrid
      call trace_summary_init_topology_result(seed_local,results(iseed), &
           located%forward_status,need_twist,need_mapping)
      iforward=2*iseed-1
      ibackward=2*iseed
      call trace_summary_init_state(seed_local,igrid,.true.,iseed, &
           need_twist,states(iforward))
      call trace_summary_init_state(seed_local,igrid,.false.,iseed, &
           need_twist,states(ibackward))
      if (allocated(map_states)) then
        call trace_summary_init_state(seed_local,igrid,.true.,iseed, &
             .false.,map_states(iforward))
        call trace_summary_init_state(seed_local,igrid,.false.,iseed, &
             .false.,map_states(ibackward))
      endif
      if (located%forward_status/=trace_status_active) then
        states(iforward)%status=located%forward_status
        states(iforward)%active=.false.
        states(ibackward)%status=located%backward_status
        states(ibackward)%active=.false.
        if (allocated(map_states)) then
          map_states(iforward)%status=located%forward_status
          map_states(iforward)%active=.false.
          map_states(ibackward)%status=located%backward_status
          map_states(ibackward)%active=.false.
        endif
      endif
    enddo

    call trace_summary_trace_states_grouped(states,2*nseed,dL,max_steps, &
         field_min)
    if (allocated(map_states)) then
      call trace_summary_trace_states_grouped(map_states,2*nseed,dL, &
           max_steps,field_min)
    endif

    do iseed=1,nseed
      iforward=2*iseed-1
      ibackward=2*iseed
      results(iseed)%forward_endpoint=states(iforward)%footpoint
      results(iseed)%backward_endpoint=states(ibackward)%footpoint
      results(iseed)%length_forward=states(iforward)%length
      results(iseed)%length_backward=states(ibackward)%length
      results(iseed)%length_total=results(iseed)%length_forward &
           +results(iseed)%length_backward
      results(iseed)%forward_nstep=states(iforward)%nstep
      results(iseed)%backward_nstep=states(ibackward)%nstep
      results(iseed)%forward_face=states(iforward)%face
      results(iseed)%backward_face=states(ibackward)%face
      results(iseed)%forward_status=states(iforward)%status
      results(iseed)%backward_status=states(ibackward)%status
      if (need_twist) then
        results(iseed)%twist_forward=states(iforward)%twist
        results(iseed)%twist_backward=states(ibackward)%twist
        results(iseed)%twist_total=results(iseed)%twist_forward &
             +results(iseed)%twist_backward
      endif

      if (states(iforward)%status==trace_status_boundary .and. &
           states(ibackward)%status==trace_status_boundary) then
        results(iseed)%status=trace_status_boundary
        results(iseed)%valid=.true.
      else
        results(iseed)%valid=.false.
        if (states(iforward)%status/=trace_status_boundary) then
          results(iseed)%status=states(iforward)%status
        else
          results(iseed)%status=states(ibackward)%status
        endif
      endif

      if (need_mapping .and. source_igrids(iseed)>=0) then
        if (allocated(map_states)) then
          call trace_summary_fill_mapping(results(iseed)%seed, &
               source_igrids(iseed),map_states(iforward), &
               map_states(ibackward),field_min,have_source_normal, &
               source_normal,mapping)
        else
          call trace_summary_fill_mapping(results(iseed)%seed, &
               source_igrids(iseed),states(iforward),states(ibackward), &
               field_min,have_source_normal,source_normal,mapping)
        endif
        results(iseed)%source_B=mapping%source_B
        results(iseed)%forward_B=mapping%forward_B
        results(iseed)%backward_B=mapping%backward_B
        results(iseed)%source_Bn=mapping%source_Bn
        results(iseed)%forward_Bn=mapping%forward_Bn
        results(iseed)%backward_Bn=mapping%backward_Bn
        results(iseed)%map_forward_endpoint=mapping%forward_footpoint
        results(iseed)%map_backward_endpoint=mapping%backward_footpoint
        results(iseed)%map_forward_length=mapping%forward_length
        results(iseed)%map_backward_length=mapping%backward_length
        results(iseed)%map_forward_face=mapping%forward_face
        results(iseed)%map_backward_face=mapping%backward_face
        results(iseed)%map_forward_status=mapping%forward_status
        results(iseed)%map_backward_status=mapping%backward_status
        results(iseed)%forward_status=mapping%forward_status
        results(iseed)%backward_status=mapping%backward_status
        results(iseed)%valid=mapping%valid
        if (mapping%valid) then
          results(iseed)%status=trace_status_boundary
        else if (mapping%forward_status/=trace_status_boundary) then
          results(iseed)%status=mapping%forward_status
        else
          results(iseed)%status=mapping%backward_status
        endif
      endif
    enddo
    if (allocated(map_states)) deallocate(map_states)
    deallocate(source_igrids,states)
  end subroutine trace_summary_topology_multi

  subroutine trace_summary_validate(dL,max_steps,status)
    double precision, intent(in) :: dL
    integer, intent(in) :: max_steps
    integer, intent(out) :: status

    status=trace_status_active
    if (npe/=1) then
      status=trace_status_mpi_unsupported
    else if (dL<=zero .or. max_steps<0) then
      status=trace_status_invalid_input
    else if (.not.slab_uniform) then
      status=trace_status_unsupported_geometry
    endif
  end subroutine trace_summary_validate

  subroutine trace_summary_init_result(seed,result,status)
    double precision, intent(in) :: seed(ndim)
    type(trace_length_result), intent(out) :: result
    integer, intent(in) :: status

    result%seed=seed
    result%forward_footpoint=seed
    result%backward_footpoint=seed
    result%forward_length=zero
    result%backward_length=zero
    result%total_length=zero
    result%forward_nstep=0
    result%backward_nstep=0
    result%forward_status=status
    result%backward_status=status
  end subroutine trace_summary_init_result

  subroutine trace_summary_init_twist_result(seed,result,status)
    double precision, intent(in) :: seed(ndim)
    type(trace_twist_result), intent(out) :: result
    integer, intent(in) :: status

    call trace_summary_init_result(seed,result%line,status)
    result%forward_twist=zero
    result%backward_twist=zero
    result%total_twist=zero
  end subroutine trace_summary_init_twist_result

  subroutine trace_summary_init_mapping_result(seed,result,status)
    double precision, intent(in) :: seed(ndim)
    type(trace_mapping_result), intent(out) :: result
    integer, intent(in) :: status

    result%seed=seed
    result%source_B=zero
    result%forward_footpoint=seed
    result%backward_footpoint=seed
    result%forward_B=zero
    result%backward_B=zero
    result%forward_length=zero
    result%backward_length=zero
    result%source_Bn=zero
    result%forward_Bn=zero
    result%backward_Bn=zero
    result%forward_face=trace_face_none
    result%backward_face=trace_face_none
    result%forward_status=status
    result%backward_status=status
    result%valid=.false.
  end subroutine trace_summary_init_mapping_result

  subroutine trace_summary_init_topology_result(seed,result,status, &
       need_twist,need_mapping)
    double precision, intent(in) :: seed(ndim)
    type(trace_topology_result), intent(out) :: result
    integer, intent(in) :: status
    logical, intent(in) :: need_twist,need_mapping

    result%seed=seed
    result%length_forward=zero
    result%length_backward=zero
    result%length_total=zero
    result%twist_forward=zero
    result%twist_backward=zero
    result%twist_total=zero
    result%forward_endpoint=seed
    result%backward_endpoint=seed
    result%forward_nstep=0
    result%backward_nstep=0
    result%forward_face=trace_face_none
    result%backward_face=trace_face_none
    result%forward_status=status
    result%backward_status=status
    result%map_forward_endpoint=seed
    result%map_backward_endpoint=seed
    result%map_forward_length=zero
    result%map_backward_length=zero
    result%map_forward_face=trace_face_none
    result%map_backward_face=trace_face_none
    result%map_forward_status=status
    result%map_backward_status=status
    result%source_B=zero
    result%forward_B=zero
    result%backward_B=zero
    result%source_Bn=zero
    result%forward_Bn=zero
    result%backward_Bn=zero
    result%has_twist=need_twist
    result%has_mapping=need_mapping
    result%valid=.false.
    result%status=status
  end subroutine trace_summary_init_topology_result

  subroutine trace_summary_locate_seed(seed,result,igrid)
    use mod_particle_base, only: find_particle_ipe

    double precision, intent(in) :: seed(ndim)
    type(trace_length_result), intent(out) :: result
    integer, intent(out) :: igrid

    double precision :: x3d(3)
    integer :: indomain,ipe,j

    call trace_summary_init_result(seed,result,trace_status_active)
    igrid=-1
    indomain=0
    {if (seed(^DB)>=xprobmin^DB .and. seed(^DB)<xprobmax^DB) indomain=indomain+1\}
    if (indomain/=ndim) then
      result%forward_status=trace_status_seed_outside
      result%backward_status=trace_status_seed_outside
      return
    endif

    x3d=zero
    do j=1,ndim
      x3d(j)=seed(j)
    enddo
    call find_particle_ipe(x3d,igrid,ipe)
    if (igrid<0 .or. ipe/=mype) then
      result%forward_status=trace_status_out_of_domain
      result%backward_status=trace_status_out_of_domain
      return
    endif
  end subroutine trace_summary_locate_seed

  subroutine trace_summary_trace_seed(igrid,dL,max_steps,field_min,result)
    integer, intent(in) :: igrid,max_steps
    double precision, intent(in) :: dL,field_min
    type(trace_length_result), intent(inout) :: result

    type(trace_summary_state) :: forward_state,backward_state

    call trace_summary_init_state(result%seed,igrid,.true.,1,.false.,forward_state)
    call trace_summary_trace_state_to_end(forward_state,dL,max_steps,field_min)
    call trace_summary_init_state(result%seed,igrid,.false.,1,.false.,backward_state)
    call trace_summary_trace_state_to_end(backward_state,dL,max_steps,field_min)

    result%forward_footpoint=forward_state%footpoint
    result%forward_length=forward_state%length
    result%forward_nstep=forward_state%nstep
    result%forward_status=forward_state%status
    result%backward_footpoint=backward_state%footpoint
    result%backward_length=backward_state%length
    result%backward_nstep=backward_state%nstep
    result%backward_status=backward_state%status
    result%total_length=result%forward_length+result%backward_length
  end subroutine trace_summary_trace_seed

  subroutine trace_summary_trace_twist_seed(igrid,dL,max_steps,field_min, &
       result)
    integer, intent(in) :: igrid,max_steps
    double precision, intent(in) :: dL,field_min
    type(trace_twist_result), intent(inout) :: result

    type(trace_summary_state) :: forward_state,backward_state

    call trace_summary_init_state(result%line%seed,igrid,.true.,1,.true., &
         forward_state)
    call trace_summary_trace_state_to_end(forward_state,dL,max_steps,field_min)
    call trace_summary_init_state(result%line%seed,igrid,.false.,1,.true., &
         backward_state)
    call trace_summary_trace_state_to_end(backward_state,dL,max_steps,field_min)

    result%line%forward_footpoint=forward_state%footpoint
    result%line%forward_length=forward_state%length
    result%line%forward_nstep=forward_state%nstep
    result%line%forward_status=forward_state%status
    result%line%backward_footpoint=backward_state%footpoint
    result%line%backward_length=backward_state%length
    result%line%backward_nstep=backward_state%nstep
    result%line%backward_status=backward_state%status
    result%line%total_length=result%line%forward_length &
         +result%line%backward_length
    result%forward_twist=forward_state%twist
    result%backward_twist=backward_state%twist
    result%total_twist=result%forward_twist+result%backward_twist
  end subroutine trace_summary_trace_twist_seed

  subroutine trace_summary_init_state(xseed,igrid,forward,seed_id, &
       accumulate_twist,state)
    double precision, intent(in) :: xseed(ndim)
    integer, intent(in) :: igrid,seed_id
    logical, intent(in) :: forward,accumulate_twist
    type(trace_summary_state), intent(out) :: state

    state%x=xseed
    state%footpoint=xseed
    state%length=zero
    state%twist=zero
    state%nstep=0
    state%status=trace_status_active
    state%igrid=igrid
    state%seed_id=seed_id
    state%face=trace_face_none
    state%forward=forward
    state%active=.true.
    state%accumulate_twist=accumulate_twist
  end subroutine trace_summary_init_state

  subroutine trace_summary_trace_state_to_end(state,ds,max_nstep,threshold)
    type(trace_summary_state), intent(inout) :: state
    double precision, intent(in) :: ds,threshold
    integer, intent(in) :: max_nstep

    do while(state%active .and. state%nstep<max_nstep)
      call trace_summary_advance_state(state,ds,threshold)
    enddo

    if (state%active) then
      state%status=trace_status_max_steps
      state%active=.false.
    endif
    state%footpoint=state%x
  end subroutine trace_summary_trace_state_to_end

  subroutine trace_summary_trace_states_grouped(states,nstate,ds,max_nstep,threshold)
    integer, intent(in) :: nstate,max_nstep
    type(trace_summary_state), intent(inout) :: states(nstate)
    double precision, intent(in) :: ds,threshold

    logical, allocatable :: processed(:)
    integer :: istate,jstate,target_grid
    logical :: any_active

    allocate(processed(nstate))
    do
      any_active=.false.
      do istate=1,nstate
        if (states(istate)%active) then
          any_active=.true.
          exit
        endif
      enddo
      if (.not.any_active) exit

      processed=.false.
      do istate=1,nstate
        if (.not.states(istate)%active .or. processed(istate)) cycle
        target_grid=states(istate)%igrid
        do jstate=1,nstate
          if (.not.states(jstate)%active .or. processed(jstate)) cycle
          if (states(jstate)%igrid/=target_grid) cycle
          call trace_summary_advance_state_in_grid(states(jstate),target_grid, &
               ds,max_nstep,threshold)
          processed(jstate)=.true.
        enddo
      enddo
    enddo
    deallocate(processed)
  end subroutine trace_summary_trace_states_grouped

  subroutine trace_summary_advance_state_in_grid(state,igrid,ds,max_nstep,threshold)
    type(trace_summary_state), intent(inout) :: state
    integer, intent(in) :: igrid,max_nstep
    double precision, intent(in) :: ds,threshold

    do while(state%active .and. state%igrid==igrid .and. &
         state%nstep<max_nstep)
      call trace_summary_advance_state(state,ds,threshold)
    enddo

    if (state%active .and. state%nstep>=max_nstep) then
      state%status=trace_status_max_steps
      state%active=.false.
      state%footpoint=state%x
    endif
  end subroutine trace_summary_advance_state_in_grid

  subroutine trace_summary_advance_state(state,ds,threshold)
    type(trace_summary_state), intent(inout) :: state
    double precision, intent(in) :: ds,threshold

    double precision :: xnext(ndim),xhit(ndim),ds_actual,segment_length
    double precision :: xbmin^D,xbmax^D
    integer :: igrid_next,ipe_next,inblock,step_status,face_id
    integer :: point_domain
    logical :: newpe,stop_trace,hit_ok

    if (.not.state%active) return

    call trace_summary_rk2_step(state%x,state%igrid,ds,state%forward,threshold, &
         xnext,ds_actual,step_status)
    if (step_status/=trace_status_active) then
      state%status=step_status
      state%active=.false.
      return
    endif

    point_domain=0
    {if (xnext(^DB)>=xprobmin^DB .and. xnext(^DB)<xprobmax^DB) point_domain=point_domain+1\}
    if (point_domain/=ndim) then
      call trace_intersect_domain(state%x,xnext,xhit,hit_ok,face_id)
      if (hit_ok) then
        segment_length=dsqrt(sum((xhit-state%x)**2))
        if (state%accumulate_twist) then
          call trace_summary_accumulate_twist(state,xhit,segment_length, &
               threshold,step_status)
        endif
        state%length=state%length+segment_length
        ! Count the nonzero terminal boundary segment as one completed step.
        state%nstep=state%nstep+1
        state%x=xhit
        state%footpoint=state%x
        state%face=face_id
      endif
      state%status=trace_status_boundary
      state%active=.false.
      return
    endif

    if (state%accumulate_twist) then
      call trace_summary_accumulate_twist(state,xnext,ds_actual,threshold, &
           step_status)
    endif
    state%length=state%length+ds_actual
    state%nstep=state%nstep+1
    state%x=xnext
    state%footpoint=state%x

    ^D&xbmin^D=rnode(rpxmin^D_,state%igrid)\
    ^D&xbmax^D=rnode(rpxmax^D_,state%igrid)\
    inblock=0
    {if (state%x(^DB)>=xbmin^DB .and. state%x(^DB)<xbmax^DB) inblock=inblock+1\}
    if (inblock==ndim) return

    igrid_next=state%igrid
    ipe_next=mype
    newpe=.false.
    stop_trace=.false.
    call find_next_grid(state%igrid,igrid_next,ipe_next,state%x,newpe,stop_trace)
    if (stop_trace) then
      state%status=trace_status_boundary
      state%active=.false.
      return
    endif
    if (newpe .or. ipe_next/=mype) then
      state%status=trace_status_out_of_domain
      state%active=.false.
      return
    endif
    state%igrid=igrid_next
  end subroutine trace_summary_advance_state

  subroutine trace_summary_accumulate_twist(state,xend,segment_length, &
       threshold,status)
    type(trace_summary_state), intent(inout) :: state
    double precision, intent(in) :: xend(ndim),segment_length,threshold
    integer, intent(out) :: status

    double precision :: xmid(ndim),B(3),curlB(3),B2,twist_density

    status=trace_status_active
    if (segment_length<=zero) return

    xmid=half*(state%x+xend)
    call sample_B_curlB_at_point(xmid,state%igrid,threshold,B,curlB,status)
    if (status/=trace_status_active) return

    B2=sum(B**2)
    twist_density=sum(curlB*B)/(4.d0*dpi*B2)
    ! Both directions contribute with positive arc length to the full-line Tw.
    state%twist=state%twist+twist_density*segment_length
  end subroutine trace_summary_accumulate_twist

  subroutine trace_total_B_at_cell(igrid,ix1,ix2,ix3,B,status)
    integer, intent(in) :: igrid,ix1,ix2,ix3
    double precision, intent(out) :: B(3)
    integer, intent(out) :: status

    B=zero
    status=trace_status_out_of_domain
    if (ndim/=3 .or. igrid<0) return

    {^IFTHREED
    if (.not.B0field .and. .not.allocated(iw_mag)) return
    if (allocated(iw_mag)) then
      B(1:3)=ps(igrid)%w(ix1,ix2,ix3,iw_mag(1:3))
    endif
    if (B0field) then
      B(1:3)=B(1:3)+ps(igrid)%B0(ix1,ix2,ix3,1:3,0)
    endif
    status=trace_status_active
    }
  end subroutine trace_total_B_at_cell

  subroutine trace_bhat_at_cell(igrid,ix1,ix2,ix3,threshold,bhat,status)
    integer, intent(in) :: igrid,ix1,ix2,ix3
    double precision, intent(in) :: threshold
    double precision, intent(out) :: bhat(3)
    integer, intent(out) :: status

    double precision :: B(3),Bnorm

    bhat=zero
    call trace_total_B_at_cell(igrid,ix1,ix2,ix3,B,status)
    if (status/=trace_status_active) return

    Bnorm=dsqrt(sum(B**2))
    if (Bnorm<=zero .or. Bnorm<threshold) then
      status=trace_status_weak_field
      return
    endif
    bhat=B/Bnorm
  end subroutine trace_bhat_at_cell

  subroutine sample_bhat_gradbhat_at_point(x,igrid,threshold,bhat, &
       grad_bhat,status)
    ! Production Method-II sampler: trilinear total-B derivative + chain rule.
    double precision, intent(in) :: x(ndim),threshold
    integer, intent(in) :: igrid
    double precision, intent(out) :: bhat(3),grad_bhat(3,3)
    integer, intent(out) :: status

    call sample_bhat_gradbhat_interpderiv_at_point(x,igrid,threshold,bhat, &
         grad_bhat,status)
  end subroutine sample_bhat_gradbhat_at_point

  subroutine sample_bhat_gradbhat_cellfd_at_point(x,igrid,threshold,bhat, &
       grad_bhat,status)
    ! Legacy diagnostic sampler: interpolate cell-centered finite-difference grad(bhat).
    double precision, intent(in) :: x(ndim),threshold
    integer, intent(in) :: igrid
    double precision, intent(out) :: bhat(3),grad_bhat(3,3)
    integer, intent(out) :: status

    double precision :: dxb^D,xd^D
    double precision :: bhat_cell(0:1^D&,3)
    double precision :: grad_cell(0:1^D&,3,3)
    double precision :: factor(0:1^D&)
    double precision :: bhat_center(3),bhat_plus(3),bhat_minus(3)
    integer :: ixI^L,ixbl^D,ix^D,j,k

    bhat=zero
    grad_bhat=zero
    status=trace_status_bad_grad_stencil
    if (ndim/=3 .or. igrid<0) return

    {^IFTHREED
    ixI^L=ixG^LL;
    ^D&dxb^D=rnode(rpdx^D_,igrid);
    ^D&ixbl^D=floor((x(^D)-ps(igrid)%x(ixImin^DD,^D))/dxb^D)+ixImin^D;
    ^D&xd^D=(x(^D)-ps(igrid)%x(ixbl^DD,^D))/dxb^D;

    if (ixbl1-1<ixImin1 .or. ixbl1+2>ixImax1 .or. &
         ixbl2-1<ixImin2 .or. ixbl2+2>ixImax2 .or. &
         ixbl3-1<ixImin3 .or. ixbl3+2>ixImax3) return

    bhat_cell=zero
    grad_cell=zero
    {do ix^D=0,1\}
      call trace_bhat_at_cell(igrid,ixbl1+ix1,ixbl2+ix2, &
           ixbl3+ix3,threshold,bhat_center,status)
      if (status/=trace_status_active) return
      bhat_cell(ix^D,1:3)=bhat_center

      call trace_bhat_at_cell(igrid,ixbl1+ix1+1,ixbl2+ix2, &
           ixbl3+ix3,threshold,bhat_plus,status)
      if (status/=trace_status_active) return
      call trace_bhat_at_cell(igrid,ixbl1+ix1-1,ixbl2+ix2, &
           ixbl3+ix3,threshold,bhat_minus,status)
      if (status/=trace_status_active) return
      grad_cell(ix^D,1:3,1)=(bhat_plus-bhat_minus)/(2.d0*dxb1)

      call trace_bhat_at_cell(igrid,ixbl1+ix1,ixbl2+ix2+1, &
           ixbl3+ix3,threshold,bhat_plus,status)
      if (status/=trace_status_active) return
      call trace_bhat_at_cell(igrid,ixbl1+ix1,ixbl2+ix2-1, &
           ixbl3+ix3,threshold,bhat_minus,status)
      if (status/=trace_status_active) return
      grad_cell(ix^D,1:3,2)=(bhat_plus-bhat_minus)/(2.d0*dxb2)

      call trace_bhat_at_cell(igrid,ixbl1+ix1,ixbl2+ix2, &
           ixbl3+ix3+1,threshold,bhat_plus,status)
      if (status/=trace_status_active) return
      call trace_bhat_at_cell(igrid,ixbl1+ix1,ixbl2+ix2, &
           ixbl3+ix3-1,threshold,bhat_minus,status)
      if (status/=trace_status_active) return
      grad_cell(ix^D,1:3,3)=(bhat_plus-bhat_minus)/(2.d0*dxb3)
    {enddo\}

    {do ix^D=0,1\}
      factor(ix^D)={abs(1-ix^D-xd^D)*}
    {enddo\}
    {do ix^D=0,1\}
      do j=1,3
        bhat(j)=bhat(j)+bhat_cell(ix^D,j)*factor(ix^D)
        do k=1,3
          grad_bhat(j,k)=grad_bhat(j,k)+grad_cell(ix^D,j,k)*factor(ix^D)
        enddo
      enddo
    {enddo\}

    status=trace_status_active
    }
  end subroutine sample_bhat_gradbhat_cellfd_at_point

  subroutine sample_bhat_gradbhat_xeps_at_point(x,eps,threshold,bhat, &
       grad_bhat,status)
    ! Reference-only grad(bhat) by central differences of interpolated bhat.
    double precision, intent(in) :: x(ndim),eps,threshold
    double precision, intent(out) :: bhat(3),grad_bhat(3,3)
    integer, intent(out) :: status

    double precision :: xp(ndim),xm(ndim)
    double precision :: B(3),Bp(3),Bm(3),bhat_p(3),bhat_m(3)
    double precision :: Bnorm,Bpnorm,Bmnorm
    integer :: igrid,idir,point_domain

    bhat=zero
    grad_bhat=zero
    status=trace_status_bad_grad_stencil
    if (ndim/=3 .or. .not.slab_uniform) then
      status=trace_status_unsupported_geometry
      return
    endif
    if (eps<=zero) then
      status=trace_status_invalid_input
      return
    endif

    call trace_debug_locate_point(x,igrid,status)
    if (status/=trace_status_active) return
    call sample_B_at_point(x,igrid,threshold,B,status)
    if (status/=trace_status_active) return
    Bnorm=dsqrt(sum(B**2))
    if (Bnorm<=zero .or. Bnorm<threshold) then
      status=trace_status_weak_field
      return
    endif
    bhat=B/Bnorm

    do idir=1,ndim
      xp=x
      xm=x
      xp(idir)=xp(idir)+eps
      xm(idir)=xm(idir)-eps

      point_domain=0
      {if (xp(^DB)>=xprobmin^DB .and. xp(^DB)<xprobmax^DB) point_domain=point_domain+1\}
      if (point_domain/=ndim) then
        status=trace_status_bad_grad_stencil
        return
      endif
      call trace_debug_locate_point(xp,igrid,status)
      if (status/=trace_status_active) return
      call sample_B_at_point(xp,igrid,threshold,Bp,status)
      if (status/=trace_status_active) return
      Bpnorm=dsqrt(sum(Bp**2))
      if (Bpnorm<=zero .or. Bpnorm<threshold) then
        status=trace_status_weak_field
        return
      endif
      bhat_p=Bp/Bpnorm

      point_domain=0
      {if (xm(^DB)>=xprobmin^DB .and. xm(^DB)<xprobmax^DB) point_domain=point_domain+1\}
      if (point_domain/=ndim) then
        status=trace_status_bad_grad_stencil
        return
      endif
      call trace_debug_locate_point(xm,igrid,status)
      if (status/=trace_status_active) return
      call sample_B_at_point(xm,igrid,threshold,Bm,status)
      if (status/=trace_status_active) return
      Bmnorm=dsqrt(sum(Bm**2))
      if (Bmnorm<=zero .or. Bmnorm<threshold) then
        status=trace_status_weak_field
        return
      endif
      bhat_m=Bm/Bmnorm

      grad_bhat(1:3,idir)=(bhat_p-bhat_m)/(2.d0*eps)
    enddo

    status=trace_status_active
  end subroutine sample_bhat_gradbhat_xeps_at_point

  subroutine sample_bhat_gradbhat_interpderiv_at_point(x,igrid,threshold, &
       bhat,grad_bhat,status)
    ! Candidate sampler: differentiate trilinear total-B interpolation.
    double precision, intent(in) :: x(ndim),threshold
    integer, intent(in) :: igrid
    double precision, intent(out) :: bhat(3),grad_bhat(3,3)
    integer, intent(out) :: status

    double precision :: dxb^D,xd^D
    double precision :: field(0:1^D&,3)
    double precision :: B(3),dBdx(3,3)
    double precision :: wx(0:1),wy(0:1),wz(0:1)
    double precision :: dwx(0:1),dwy(0:1),dwz(0:1)
    double precision :: weight,Bnorm,projected
    integer :: ixI^L,ixbl^D,ix^D,j

    bhat=zero
    grad_bhat=zero
    status=trace_status_bad_grad_stencil
    if (ndim/=3 .or. .not.slab_uniform .or. igrid<0) then
      status=trace_status_unsupported_geometry
      return
    endif

    {^IFTHREED
    ixI^L=ixG^LL;
    ^D&dxb^D=rnode(rpdx^D_,igrid);
    ^D&ixbl^D=floor((x(^D)-ps(igrid)%x(ixImin^DD,^D))/dxb^D)+ixImin^D;
    ^D&xd^D=(x(^D)-ps(igrid)%x(ixbl^DD,^D))/dxb^D;

    if (ixbl1<ixImin1 .or. ixbl1+1>ixImax1 .or. &
         ixbl2<ixImin2 .or. ixbl2+1>ixImax2 .or. &
         ixbl3<ixImin3 .or. ixbl3+1>ixImax3) return

    status=trace_status_out_of_domain
    if (.not.B0field .and. .not.allocated(iw_mag)) return

    field=zero
    {do ix^D=0,1\}
      if (allocated(iw_mag)) then
        field(ix^D,1:3)=ps(igrid)%w(ixbl1+ix1,ixbl2+ix2, &
             ixbl3+ix3,iw_mag(1:3))
      endif
      if (B0field) then
        field(ix^D,1:3)=field(ix^D,1:3) &
             +ps(igrid)%B0(ixbl1+ix1,ixbl2+ix2,ixbl3+ix3,1:3,0)
      endif
    {enddo\}

    wx(0)=one-xd1
    wx(1)=xd1
    wy(0)=one-xd2
    wy(1)=xd2
    wz(0)=one-xd3
    wz(1)=xd3
    dwx(0)=-one/dxb1
    dwx(1)= one/dxb1
    dwy(0)=-one/dxb2
    dwy(1)= one/dxb2
    dwz(0)=-one/dxb3
    dwz(1)= one/dxb3

    B=zero
    dBdx=zero
    {do ix^D=0,1\}
      weight=wx(ix1)*wy(ix2)*wz(ix3)
      B(1:3)=B(1:3)+field(ix^D,1:3)*weight
      dBdx(1:3,1)=dBdx(1:3,1)+field(ix^D,1:3)* &
           dwx(ix1)*wy(ix2)*wz(ix3)
      dBdx(1:3,2)=dBdx(1:3,2)+field(ix^D,1:3)* &
           wx(ix1)*dwy(ix2)*wz(ix3)
      dBdx(1:3,3)=dBdx(1:3,3)+field(ix^D,1:3)* &
           wx(ix1)*wy(ix2)*dwz(ix3)
    {enddo\}

    Bnorm=dsqrt(sum(B**2))
    if (Bnorm<=zero .or. Bnorm<threshold) then
      status=trace_status_weak_field
      return
    endif

    bhat=B/Bnorm
    do j=1,3
      projected=sum(bhat*dBdx(1:3,j))
      grad_bhat(1:3,j)=(dBdx(1:3,j)-bhat(1:3)*projected)/Bnorm
    enddo

    status=trace_status_active
    }
  end subroutine sample_bhat_gradbhat_interpderiv_at_point

  subroutine sample_B_curlB_at_point(x,igrid,threshold,B,curlB,status)
    ! Interpolate total B and its cell-centered second-order curl locally.
    double precision, intent(in) :: x(ndim),threshold
    integer, intent(in) :: igrid
    double precision, intent(out) :: B(3),curlB(3)
    integer, intent(out) :: status

    double precision :: dxb^D,xd^D
    double precision :: field(0:1^D&,3),current(0:1^D&,3)
    double precision :: factor(0:1^D&),B2
    integer :: ixI^L,ixJ^L,ixbl^D,ix^D,j

    B=zero
    curlB=zero
    status=trace_status_bad_curl_stencil
    if (ndim/=3) return

    {^IFTHREED
    ixI^L=ixG^LL;
    ixJ^L=ixM^LL^LADD1;
    ^D&dxb^D=rnode(rpdx^D_,igrid);
    ^D&ixbl^D=floor((x(^D)-ps(igrid)%x(ixImin^DD,^D))/dxb^D)+ixImin^D;
    ^D&xd^D=(x(^D)-ps(igrid)%x(ixbl^DD,^D))/dxb^D;

    if (ixbl1-1<ixImin1 .or. ixbl1+2>ixImax1 .or. &
         ixbl2-1<ixImin2 .or. ixbl2+2>ixImax2 .or. &
         ixbl3-1<ixImin3 .or. ixbl3+2>ixImax3) return
    if (B0field) then
      if (ixbl1<ixJmin1 .or. ixbl1+1>ixJmax1 .or. &
           ixbl2<ixJmin2 .or. ixbl2+1>ixJmax2 .or. &
           ixbl3<ixJmin3 .or. ixbl3+1>ixJmax3) return
    endif
    status=trace_status_out_of_domain
    if (.not.B0field .and. .not.allocated(iw_mag)) return

    field=zero
    current=zero
    {do ix^D=0,1\}
      if (allocated(iw_mag)) then
        field(ix^D,1:3)=ps(igrid)%w(ixbl1+ix1,ixbl2+ix2, &
             ixbl3+ix3,iw_mag(1:3))
      endif
      if (B0field) then
        field(ix^D,1:3)=field(ix^D,1:3) &
             +ps(igrid)%B0(ixbl1+ix1,ixbl2+ix2,ixbl3+ix3,1:3,0)
      endif
      if (allocated(iw_mag)) then
        current(ix^D,1)=half*( &
             (ps(igrid)%w(ixbl1+ix1,ixbl2+ix2+1,ixbl3+ix3,iw_mag(3)) &
             -ps(igrid)%w(ixbl1+ix1,ixbl2+ix2-1,ixbl3+ix3,iw_mag(3)))/dxb2 &
             -(ps(igrid)%w(ixbl1+ix1,ixbl2+ix2,ixbl3+ix3+1,iw_mag(2)) &
             -ps(igrid)%w(ixbl1+ix1,ixbl2+ix2,ixbl3+ix3-1,iw_mag(2)))/dxb3)
        current(ix^D,2)=half*( &
             (ps(igrid)%w(ixbl1+ix1,ixbl2+ix2,ixbl3+ix3+1,iw_mag(1)) &
             -ps(igrid)%w(ixbl1+ix1,ixbl2+ix2,ixbl3+ix3-1,iw_mag(1)))/dxb3 &
             -(ps(igrid)%w(ixbl1+ix1+1,ixbl2+ix2,ixbl3+ix3,iw_mag(3)) &
             -ps(igrid)%w(ixbl1+ix1-1,ixbl2+ix2,ixbl3+ix3,iw_mag(3)))/dxb1)
        current(ix^D,3)=half*( &
             (ps(igrid)%w(ixbl1+ix1+1,ixbl2+ix2,ixbl3+ix3,iw_mag(2)) &
             -ps(igrid)%w(ixbl1+ix1-1,ixbl2+ix2,ixbl3+ix3,iw_mag(2)))/dxb1 &
             -(ps(igrid)%w(ixbl1+ix1,ixbl2+ix2+1,ixbl3+ix3,iw_mag(1)) &
             -ps(igrid)%w(ixbl1+ix1,ixbl2+ix2-1,ixbl3+ix3,iw_mag(1)))/dxb2)
      endif
      if (B0field) then
        current(ix^D,1:3)=current(ix^D,1:3) &
             +ps(igrid)%J0(ixbl^D+ix^D,1:3)
      endif
    {enddo\}

    {do ix^D=0,1\}
      factor(ix^D)={abs(1-ix^D-xd^D)*}
    {enddo\}
    {do ix^D=0,1\}
      do j=1,3
        B(j)=B(j)+field(ix^D,j)*factor(ix^D)
        curlB(j)=curlB(j)+current(ix^D,j)*factor(ix^D)
      enddo
    {enddo\}

    B2=sum(B**2)
    if (B2<=zero .or. dsqrt(B2)<threshold) then
      status=trace_status_weak_field
      return
    endif
    status=trace_status_active
    }
  end subroutine sample_B_curlB_at_point

  subroutine sample_B_at_point(x,igrid,threshold,B,status)
    ! Interpolate total magnetic field at x using the local block stencil.
    double precision, intent(in) :: x(ndim),threshold
    integer, intent(in) :: igrid
    double precision, intent(out) :: B(3)
    integer, intent(out) :: status

    double precision :: dxb^D,xd^D
    double precision :: field(0:1^D&,3),factor(0:1^D&),B2
    double precision :: Bcell(3)
    integer :: ixI^L,ixbl^D,ix^D,j

    B=zero
    status=trace_status_out_of_domain
    if (ndim/=3 .or. igrid<0) return

    {^IFTHREED
    ixI^L=ixG^LL;
    ^D&dxb^D=rnode(rpdx^D_,igrid);
    ^D&ixbl^D=floor((x(^D)-ps(igrid)%x(ixImin^DD,^D))/dxb^D)+ixImin^D;
    ^D&xd^D=(x(^D)-ps(igrid)%x(ixbl^DD,^D))/dxb^D;

    if (ixbl1<ixImin1 .or. ixbl1+1>ixImax1 .or. &
         ixbl2<ixImin2 .or. ixbl2+1>ixImax2 .or. &
         ixbl3<ixImin3 .or. ixbl3+1>ixImax3) return
    field=zero
    {do ix^D=0,1\}
      call trace_total_B_at_cell(igrid,ixbl1+ix1,ixbl2+ix2, &
           ixbl3+ix3,Bcell,status)
      if (status/=trace_status_active) return
      field(ix^D,1:3)=Bcell
    {enddo\}

    {do ix^D=0,1\}
      factor(ix^D)={abs(1-ix^D-xd^D)*}
    {enddo\}
    {do ix^D=0,1\}
      do j=1,3
        B(j)=B(j)+field(ix^D,j)*factor(ix^D)
      enddo
    {enddo\}

    B2=sum(B**2)
    if (B2<=zero .or. dsqrt(B2)<threshold) then
      status=trace_status_weak_field
      return
    endif
    status=trace_status_active
    }
  end subroutine sample_B_at_point

  subroutine trace_summary_fill_mapping(seed,source_igrid,forward_state, &
       backward_state,threshold,have_source_normal,source_normal,result)
    double precision, intent(in) :: seed(ndim),threshold,source_normal(3)
    integer, intent(in) :: source_igrid
    logical, intent(in) :: have_source_normal
    type(trace_summary_state), intent(in) :: forward_state,backward_state
    type(trace_mapping_result), intent(inout) :: result

    integer :: status

    call trace_summary_init_mapping_result(seed,result,trace_status_active)

    result%forward_footpoint=forward_state%footpoint
    result%backward_footpoint=backward_state%footpoint
    result%forward_length=forward_state%length
    result%backward_length=backward_state%length
    result%forward_face=forward_state%face
    result%backward_face=backward_state%face
    result%forward_status=forward_state%status
    result%backward_status=backward_state%status

    call sample_B_at_point(seed,source_igrid,threshold,result%source_B,status)
    if (status/=trace_status_active) then
      result%forward_status=status
      result%backward_status=status
      result%valid=.false.
      return
    endif
    if (have_source_normal) then
      result%source_Bn=sum(result%source_B*source_normal)
    endif

    if (forward_state%status==trace_status_boundary) then
      call sample_B_at_point(forward_state%footpoint,forward_state%igrid, &
           threshold,result%forward_B,status)
      if (status/=trace_status_active) result%forward_status=status
    endif
    if (backward_state%status==trace_status_boundary) then
      call sample_B_at_point(backward_state%footpoint,backward_state%igrid, &
           threshold,result%backward_B,status)
      if (status/=trace_status_active) result%backward_status=status
    endif

    result%forward_Bn=trace_face_Bn(result%forward_face,result%forward_B)
    result%backward_Bn=trace_face_Bn(result%backward_face,result%backward_B)
    result%valid=result%forward_status==trace_status_boundary .and. &
         result%backward_status==trace_status_boundary .and. &
         trace_face_is_boundary(result%forward_face) .and. &
         trace_face_is_boundary(result%backward_face)
  end subroutine trace_summary_fill_mapping

  double precision function trace_face_Bn(face_id,B) result(Bn)
    integer, intent(in) :: face_id
    double precision, intent(in) :: B(3)

    Bn=zero
    select case(face_id)
    case(trace_face_xmin)
      Bn=-B(1)
    case(trace_face_xmax)
      Bn=B(1)
    case(trace_face_ymin)
      Bn=-B(2)
    case(trace_face_ymax)
      Bn=B(2)
    case(trace_face_zmin)
      Bn=-B(3)
    case(trace_face_zmax)
      Bn=B(3)
    end select
  end function trace_face_Bn

  logical function trace_face_is_boundary(face_id) result(is_boundary)
    integer, intent(in) :: face_id

    is_boundary=face_id>=trace_face_xmin .and. face_id<=trace_face_zmax
  end function trace_face_is_boundary

  subroutine trace_boundary_face_at_point(x,face_id,on_boundary)
    double precision, intent(in) :: x(ndim)
    integer, intent(out) :: face_id
    logical, intent(out) :: on_boundary

    double precision :: domain_min(ndim),domain_max(ndim),tol
    integer :: idim,nface,candidate_face

    ^D&domain_min(^D)=xprobmin^D;
    ^D&domain_max(^D)=xprobmax^D;
    tol=100.d0*epsilon(one)*max(one,maxval(abs(domain_max-domain_min)))
    face_id=trace_face_none
    nface=0
    do idim=1,ndim
      if (abs(x(idim)-domain_min(idim))<=tol) then
        candidate_face=trace_face_from_dim_side(idim,-1)
        nface=nface+1
        if (nface==1) then
          face_id=candidate_face
        else if (candidate_face/=face_id) then
          face_id=trace_face_ambiguous
        endif
      endif
      if (abs(x(idim)-domain_max(idim))<=tol) then
        candidate_face=trace_face_from_dim_side(idim,1)
        nface=nface+1
        if (nface==1) then
          face_id=candidate_face
        else if (candidate_face/=face_id) then
          face_id=trace_face_ambiguous
        endif
      endif
    enddo
    on_boundary=nface>0
  end subroutine trace_boundary_face_at_point

  subroutine trace_intersect_domain(xinside,xoutside,xhit,hit_ok,face_id)
    double precision, intent(in) :: xinside(ndim),xoutside(ndim)
    double precision, intent(out) :: xhit(ndim)
    logical, intent(out) :: hit_ok
    integer, intent(out) :: face_id

    double precision :: alpha,alpha_hit,alpha_tol,delta
    double precision :: domain_min(ndim),domain_max(ndim)
    integer :: idim,hit_dim,hit_side,candidate_face

    ^D&domain_min(^D)=xprobmin^D;
    ^D&domain_max(^D)=xprobmax^D;
    alpha_hit=huge(one)
    alpha_tol=100.d0*epsilon(one)
    hit_dim=0
    hit_side=0
    face_id=trace_face_none
    do idim=1,ndim
      delta=xoutside(idim)-xinside(idim)
      if (xoutside(idim)<domain_min(idim) .and. abs(delta)>smalldouble) then
        alpha=(domain_min(idim)-xinside(idim))/delta
        if (alpha>zero .and. alpha<=one) then
          candidate_face=trace_face_from_dim_side(idim,-1)
          if (alpha<alpha_hit-alpha_tol) then
            alpha_hit=alpha
            hit_dim=idim
            hit_side=-1
            face_id=candidate_face
          else if (abs(alpha-alpha_hit)<=alpha_tol .and. &
               candidate_face/=face_id) then
            face_id=trace_face_ambiguous
          endif
        endif
      else if (xoutside(idim)>=domain_max(idim) .and. &
           abs(delta)>smalldouble) then
        alpha=(domain_max(idim)-xinside(idim))/delta
        if (alpha>zero .and. alpha<=one) then
          candidate_face=trace_face_from_dim_side(idim,1)
          if (alpha<alpha_hit-alpha_tol) then
            alpha_hit=alpha
            hit_dim=idim
            hit_side=1
            face_id=candidate_face
          else if (abs(alpha-alpha_hit)<=alpha_tol .and. &
               candidate_face/=face_id) then
            face_id=trace_face_ambiguous
          endif
        endif
      endif
    enddo

    hit_ok=hit_dim>0
    xhit=xinside
    if (.not.hit_ok) return

    xhit=xinside+alpha_hit*(xoutside-xinside)
    do idim=1,ndim
      xhit(idim)=min(max(xhit(idim),domain_min(idim)),domain_max(idim))
    enddo
    if (hit_side<0) then
      xhit(hit_dim)=domain_min(hit_dim)
    else
      xhit(hit_dim)=domain_max(hit_dim)
    endif
  end subroutine trace_intersect_domain

  integer function trace_face_from_dim_side(idim,side) result(face_id)
    integer, intent(in) :: idim,side

    face_id=trace_face_none
    select case(idim)
    case(1)
      if (side<0) then
        face_id=trace_face_xmin
      else
        face_id=trace_face_xmax
      endif
    case(2)
      if (side<0) then
        face_id=trace_face_ymin
      else
        face_id=trace_face_ymax
      endif
    case(3)
      if (side<0) then
        face_id=trace_face_zmin
      else
        face_id=trace_face_zmax
      endif
    end select
  end function trace_face_from_dim_side

  subroutine trace_summary_rk2_step(xnow,igrid,ds,forward,threshold, &
       xnext,ds_actual,status)
    double precision, intent(in) :: xnow(ndim),ds,threshold
    integer, intent(in) :: igrid
    logical, intent(in) :: forward
    double precision, intent(out) :: xnext(ndim),ds_actual
    integer, intent(out) :: status

    double precision :: xstage(ndim),xstage_sample(ndim),K1(ndim),K2(ndim)
    double precision :: dxb^D
    integer :: ixI^L
    integer :: igrid_stage
    logical :: field_ok
    character(len=std_len) :: field_type

    ixI^L=ixG^LL;
    ^D&dxb^D=rnode(rpdx^D_,igrid);
    field_type='Bfield'
    xnext=xnow
    ds_actual=zero
    status=trace_status_active

    call get_K(xnow,igrid,K1,ixI^L,dxb^D,field_type,threshold,field_ok)
    if (.not.field_ok) then
      status=trace_status_weak_field
      return
    endif

    if (forward) then
      xstage=xnow+ds*K1
    else
      xstage=xnow-ds*K1
    endif
    xstage_sample=xstage
    ^D&xstage_sample(^D)=max(xprobmin^D,min(xstage_sample(^D), &
         xprobmax^D-100.d0*epsilon(one)*max(one,abs(xprobmax^D))));
    call trace_locate_point_with_hint(xstage_sample,igrid,igrid_stage,status)
    if (status/=trace_status_active) return
    ^D&dxb^D=rnode(rpdx^D_,igrid_stage);
    call get_K(xstage_sample,igrid_stage,K2,ixI^L,dxb^D,field_type, &
         threshold,field_ok)
    if (.not.field_ok) then
      status=trace_status_weak_field
      return
    endif

    if (forward) then
      xnext=xnow+ds*(half*K1+half*K2)
    else
      xnext=xnow-ds*(half*K1+half*K2)
    endif
    ds_actual=dsqrt(sum((xnext-xnow)**2))
  end subroutine trace_summary_rk2_step

  subroutine trace_field_length_multi(seeds,nseed,dL,max_steps,results,b_min)
    ! Trace multiple seeds without storing field-line paths.
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    type(trace_length_result), intent(out) :: results(nseed)
    double precision, intent(in), optional :: b_min

    double precision :: field_min

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)
    call trace_summary_multi(seeds,nseed,dL,max_steps,results,field_min)
  end subroutine trace_field_length_multi

  subroutine trace_field_twist_multi(seeds,nseed,dL,max_steps,results,b_min)
    ! Trace multiple seeds and integrate twist without storing field-line paths.
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    type(trace_twist_result), intent(out) :: results(nseed)
    double precision, intent(in), optional :: b_min

    double precision :: field_min

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)
    call trace_summary_twist_multi(seeds,nseed,dL,max_steps,results,field_min)
  end subroutine trace_field_twist_multi

  subroutine trace_field_mapping_multi(seeds,nseed,dL,max_steps,results, &
       b_min,source_normal)
    ! Trace multiple seeds and return endpoint metadata without storing paths.
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    type(trace_mapping_result), intent(out) :: results(nseed)
    double precision, intent(in), optional :: b_min
    double precision, intent(in), optional :: source_normal(3)

    double precision :: field_min,normal(3)
    logical :: have_normal

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)
    normal=zero
    have_normal=present(source_normal)
    if (have_normal) normal=source_normal
    call trace_summary_mapping_multi(seeds,nseed,dL,max_steps,results, &
         field_min,have_normal,normal)
  end subroutine trace_field_mapping_multi

  subroutine trace_field_topology_multi(seeds,nseed,dL,max_steps,results, &
       need_twist,need_mapping,b_min,source_normal)
    ! Trace multiple seeds once and return length plus optional Tw/mapping.
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    type(trace_topology_result), intent(out) :: results(nseed)
    logical, intent(in), optional :: need_twist,need_mapping
    double precision, intent(in), optional :: b_min
    double precision, intent(in), optional :: source_normal(3)

    double precision :: field_min,normal(3)
    logical :: do_twist,do_mapping,have_normal

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)
    do_twist=.false.
    do_mapping=.false.
    if (present(need_twist)) do_twist=need_twist
    if (present(need_mapping)) do_mapping=need_mapping
    normal=zero
    have_normal=present(source_normal)
    if (have_normal) normal=source_normal

    call trace_summary_topology_multi(seeds,nseed,dL,max_steps,results, &
         field_min,do_twist,do_mapping,have_normal,normal)
  end subroutine trace_field_topology_multi

  subroutine trace_field_qperp_multi(seeds,nseed,dL,max_steps,results,b_min)
    ! Block-grouped tangent-state tracing for multiple Method-II/Q_perp seeds.
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    type(trace_qperp_result), intent(out) :: results(nseed)
    double precision, intent(in), optional :: b_min

    type(trace_tangent_state), allocatable :: states(:)
    double precision :: field_min
    double precision :: seed_local(ndim),zero_vec(ndim)
    integer :: common_status,iseed,idim,igrid,iforward,ibackward
    integer :: seed_status

    if (nseed<=0) return

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)
    zero_vec=zero

    common_status=trace_status_active
    if (npe/=1) then
      common_status=trace_status_mpi_unsupported
    else if (ndim/=3 .or. .not.slab_uniform) then
      common_status=trace_status_unsupported_geometry
    else if (dL<=zero .or. max_steps<=0) then
      common_status=trace_status_invalid_input
    endif
    if (common_status/=trace_status_active) then
      do iseed=1,nseed
        do idim=1,ndim
          seed_local(idim)=seeds(iseed,idim)
        enddo
        call trace_init_qperp_result(seed_local,results(iseed))
        results(iseed)%status=common_status
      enddo
      return
    endif

    allocate(states(2*nseed))
    do iseed=1,nseed
      do idim=1,ndim
        seed_local(idim)=seeds(iseed,idim)
      enddo
      call trace_qperp_prepare_seed_result(seed_local,field_min, &
           results(iseed),igrid,seed_status)
      iforward=2*iseed-1
      ibackward=2*iseed
      if (seed_status==trace_status_active) then
        call trace_tangent_state_init(seed_local,results(iseed)%u0, &
             results(iseed)%v0,igrid,iseed,1,states(iforward))
        call trace_tangent_state_init(seed_local,results(iseed)%u0, &
             results(iseed)%v0,igrid,iseed,-1,states(ibackward))
      else
        call trace_tangent_state_init(seed_local,zero_vec,zero_vec,igrid, &
             iseed,1,states(iforward))
        call trace_tangent_state_init(seed_local,zero_vec,zero_vec,igrid, &
             iseed,-1,states(ibackward))
        states(iforward)%active=.false.
        states(ibackward)%active=.false.
        states(iforward)%status=seed_status
        states(ibackward)%status=seed_status
      endif
    enddo

    call trace_tangent_trace_states_grouped(states,2*nseed,dL,max_steps, &
         field_min)

    do iseed=1,nseed
      iforward=2*iseed-1
      ibackward=2*iseed
      if (results(iseed)%status==trace_status_active) then
        call trace_qperp_finalize_from_states(states(iforward), &
             states(ibackward),results(iseed),field_min)
      endif
    enddo
    deallocate(states)
  end subroutine trace_field_qperp_multi

  subroutine trace_field_multi(xfm,wPm,wLm,dL,numL,numP,nwP,nwL,forwardm,ftype,tcondi)
    ! trace multiple field lines
    ! xfm: locations of points at the field lines. User should provide xfm(1:numL,1,1:ndim)
    !   as seed points, then field line wills be traced from the seed points. xfm(1:numL,1,:) 
    !   are user given and other points are given by the subroutine as feedback.
    ! numL: number of field lines user wants to trace. user given
    ! numP: maximum number of points at the field line. User defined. Note that not every
    !   point of numP is valid, as the tracing will stop when the field line leave the
    !   simulation box. The number of valid point is given in wLm(1)
    ! wPm: point variables, which have values at each point of xfm. The way to get wP is
    !   user defined, with help of IO given by subroutine set_field_w in mod_usr_methods.
    !   User can calculate the density/temperature at the point and then store the values
    !   into wPm, or do something else.
    ! nwP: number of point variables. user given
    ! wLm: line variables, variables for the lines rather than for each point of the lines.
    !   For example, wLm(1) stores the number of valid points at the field lines. The way to 
    !   get wLm is also user defined with set_field_w, the same as wPm. User can calculate
    !   the maximum density/temperature at the field lines and stores it in wLm
    ! nwL: number of line variables. user given
    ! dL: length step for field line tracing. user given
    ! forwardm: true--trace field forward; false--trace field line backward. user given
    ! ftype: type of field user wants to trace. User can trace velocity field by setting
    !   ftype='Vfield' or trace magnetic field by setting ftype='Bfield'. It is possible 
    !   to trace other fields, e.g. electric field, where user can define the field with
    !   IO given by subroutine set_field in mod_usr_methods. user given
    ! tcondi: user given
    use mod_particle_base

    integer, intent(in) :: numL,numP,nwP,nwL
    double precision, intent(inout) :: xfm(numL,numP,ndim),wPm(numL,numP,nwP),wLm(numL,1+nwL)
    double precision, intent(in) :: dL
    logical, intent(in) :: forwardm(numL)
    character(len=std_len), intent(in) :: ftype,tcondi

    double precision :: x3d(3),statusF(4+ndim),statusL(numL,4+ndim+nwL),statusS(numL,4+ndim+nwL)
    double precision :: xf(numP,ndim),wP(numP,nwP),wL(1+nwL)
    double precision, allocatable :: data_send(:,:,:),data_recv(:,:,:)
    integer :: indomain,ipe_now,igrid_now,igrid,j,iL
    integer :: ipoint_in,ipoint_out,numSend,nRT,nRTmax
    integer :: ipointm(numL),igridm(numL)
    logical :: continueL(numL),myL(numL)
    logical :: stopT,forward

    if (tcondi/='TRAC') then
      wPm=zero
    else
      wPm=-1
    endif
    wLm=zero
    xfm(1:numL,2:numP,:)=zero
    stopT=.TRUE.
    myL=.FALSE.
    xf=zero
    wP=zero

    ! find the pe and igrid for the first point
    do iL=1,numL
      indomain=0
      wLm(iL,1)=0
      {if (xfm(iL,1,^DB)>=xprobmin^DB .and. xfm(iL,1,^DB)<xprobmax^DB) indomain=indomain+1\}
      if (indomain==ndim) then
        if (tcondi/='TRAC') wLm(iL,1)=1
        continueL(iL)=.TRUE.
        ! find pe and igrid
        x3d=0.d0
        do j=1,ndim
          x3d(j)=xfm(iL,1,j)
        enddo
        call find_particle_ipe(x3d,igrid_now,ipe_now)
        stopT=.FALSE.
        ipointm(iL)=1
        igridm(iL)=igrid_now
        if (mype==ipe_now) then 
          myL(iL)=.TRUE.
        else
          xfm(iL,1,:)=zero
        endif
      else
        continueL(iL)=.FALSE.
        wLm(iL,1)=zero
      endif
    enddo

    do while(stopT .eqv. .FALSE.)
      ! tracing multiple field lines inside pe
      statusS=zero
      do iL=1,numL
        if (myL(iL) .and. continueL(iL)) then
          igrid=igridm(iL)
          ipoint_in=ipointm(iL)
          xf(ipoint_in,:)=xfm(iL,ipoint_in,:)
          wL(:)=wLm(iL,:)
          forward=forwardm(iL)
          statusF=zero
          call find_points_in_pe(igrid,ipoint_in,xf,wP,wL,dL,numP,nwP,nwL,forward,ftype,tcondi,statusF)
          ipoint_out=int(statusF(1))
          xfm(iL,ipoint_in:ipoint_out-1,:)=xf(ipoint_in:ipoint_out-1,:)
          wPm(iL,ipoint_in:ipoint_out-1,:)=wP(ipoint_in:ipoint_out-1,:)
          ! status for each field line
          ! 1: index of next point
          ! 2: ipe of next point
          ! 3: igrid of next point
          ! 4: trace_status_active -> continue; otherwise stop tracing
          ! 5:4+ndim: coordinate of next point
          ! 4+ndim+1:4+ndim+nwL: wL(2:1+nwL)
          ! for TRAC nwL=2 -> wL(2): current Tcoff; wL(3): Tmax
          ! for TRAC nwP=2 -> wP(:,1):ipe; wP(:,2):igrid
          statusS(iL,1:4+ndim)=statusF(1:4+ndim)
          statusS(iL,4+ndim+1:4+ndim+nwL)=wL(2:1+nwL)
          if (tcondi=='TRAC') wLm(iL,1)=ipoint_out-1
        endif
      enddo

      ! comunicating tracing results
      numSend=numL*(4+ndim+nwL)
      call MPI_ALLREDUCE(statusS,statusL,numSend,MPI_DOUBLE_PRECISION,&
                       MPI_SUM,icomm,ierrmpi)

      ! for next step
      stopT=.TRUE.
      myL=.FALSE.
      do iL=1,numL
        if (continueL(iL)) then
          ipointm(iL)=int(statusL(iL,1))
          if (mype==int(statusL(iL,2))) myL(iL)=.TRUE.
          igridm(iL)=int(statusL(iL,3))
          if (int(statusL(iL,4))==trace_status_active) then
            stopT=.FALSE.
          else
            continueL(iL)=.FALSE.
          endif
          if (myL(iL)) xfm(iL,ipointm(iL),1:ndim)=statusL(iL,4+1:4+ndim)
          if (tcondi/='TRAC') then
            if (int(statusL(iL,4))==trace_status_weak_field) then
              wLm(iL,1)=ipointm(iL)
            else
              wLm(iL,1)=ipointm(iL)-1
            endif
          endif
          wLm(iL,2:1+nwL)=statusL(iL,4+ndim+1:4+ndim+nwL)
        endif
      enddo
    enddo

    ! communication after tracing
    if (tcondi/='TRAC') then
      nRTmax=0
      do iL=1,numL
        if (nRTmax<int(wLm(iL,1))) nRTmax=int(wLm(iL,1))
      enddo
      numSend=numL*nRTmax*(ndim+nwP)

      allocate(data_send(numL,nRTmax,ndim+nwP),data_recv(numL,nRTmax,ndim+nwP))
      data_send(:,:,:)=zero
      do iL=1,numL
        nRT=int(wLm(iL,1))
        data_send(iL,1:nRT,1:ndim)=xfm(iL,1:nRT,1:ndim)
        if (nwP>0) data_send(iL,1:nRT,1+ndim:ndim+nwP)=wPm(iL,1:nRT,1:nwP)
      enddo
      call MPI_ALLREDUCE(data_send,data_recv,numSend,MPI_DOUBLE_PRECISION,&
                         MPI_SUM,icomm,ierrmpi)
      do iL=1,numL
        nRT=int(wLm(iL,1))
        xfm(iL,1:nRT,1:ndim)=data_recv(iL,1:nRT,1:ndim)
        if (nwP>0) wPm(iL,1:nRT,1:nwP)=data_recv(iL,1:nRT,1+ndim:ndim+nwP)
      enddo
      deallocate(data_send,data_recv)
    endif

  end subroutine trace_field_multi

  subroutine trace_field_single(xf,wP,wL,dL,numP,nwP,nwL,forward,ftype,tcondi)
    ! trace a field line
    ! xf: locations of points at the field line. User should provide xf(1,1:ndim)
    !   as seed point, then field line will be traced from the seed point. xf(1,:) is 
    !   user given and xf(2:wL(1),:) are given by the subroutine as feedback.
    ! numP: maximum number of points at the field line. User defined. Note that not every
    !   point of numP is valid, as the tracing will stop when the field line leave the
    !   simulation box. The number of valid point is given in wL(1)
    ! wP: point variables, which have values at each point of xf. The way to get wP is
    !   user defined, with help of IO given by subroutine set_field_w in mod_usr_methods.
    !   User can calculate the density/temperature at the point and then store the values
    !   into wP, or do something else.
    ! nwP: number of point variables. user given
    ! wL: line variables, variables for the line rather than for each point of the line.
    !   For example, wL(1) stores the number of valid points at the field line. The way to 
    !   get wL is also user defined with set_field_w, the same as wP. User can calculate
    !   the maximum density/temperature at the field line and stores it in wL
    ! nwL: number of line variables. user given
    ! dL: length step for field line tracing. user given
    ! forward: true--trace field forward; false--trace field line backward. user given
    ! ftype: type of field user wants to trace. User can trace velocity field by setting
    !   ftype='Vfield' or trace magnetic field by setting ftype='Bfield'. It is possible 
    !   to trace other fields, e.g. electric field, where user can define the field with
    !   IO given by subroutine set_field in mod_usr_methods. user given
    ! tcondi: user given
    use mod_usr_methods
    use mod_particle_base

    integer, intent(in) :: numP,nwP,nwL
    double precision, intent(inout) :: xf(numP,ndim),wP(numP,nwP),wL(1+nwL)
    double precision, intent(in) :: dL
    logical, intent(in) :: forward
    character(len=std_len), intent(in) :: ftype,tcondi

    double precision :: x3d(3),statusF(4+ndim),status_bcast(4+ndim+nwL)
    double precision, allocatable :: data_send(:,:),data_recv(:,:)
    integer :: indomain,ipoint_in,ipe_now,igrid_now,igrid,j
    integer :: ipoint_out,ipe_next,igrid_next,numRT
    logical :: stopT

    wP=zero
    wL=zero
    xf(2:numP,:)=zero

    ! check whether or the first point is inside simulation box. if yes, find
    ! the pe and igrid for the point
    indomain=0
    wL(1)=0
    {if (xf(1,^DB)>=xprobmin^DB .and. xf(1,^DB)<xprobmax^DB) indomain=indomain+1\}
    if (indomain==ndim) then
      wL(1)=1

       ! find pe and igrid
       x3d=0.d0
       do j=1,ndim
         x3d(j)=xf(1,j)
       enddo
      call find_particle_ipe(x3d,igrid_now,ipe_now)
      stopT=.FALSE.
      ipoint_in=1
      if (mype/=ipe_now) xf(1,:)=zero
    else
      if (mype==0) then
        call MPISTOP('Field tracing error: given point is not in simulation box!')
      endif
    endif


    ! other points in field line
    do while(stopT .eqv. .FALSE.)

      if (mype==ipe_now) then
        igrid=igrid_now
        ! looking for points in one pe
        call find_points_in_pe(igrid,ipoint_in,xf,wP,wL,dL,numP,nwP,nwL,forward,ftype,tcondi,statusF)
        status_bcast(1:4+ndim)=statusF(1:4+ndim)
        status_bcast(4+ndim+1:4+ndim+nwL)=wL(2:1+nwL)
      endif
      ! comunication
      call MPI_BCAST(status_bcast,4+ndim+nwL,MPI_DOUBLE_PRECISION,ipe_now,icomm,ierrmpi)
      statusF(1:4+ndim)=status_bcast(1:4+ndim)
      wL(2:1+nwL)=status_bcast(4+ndim+1:4+ndim+nwL)

      ! prepare for next step
      ipoint_out=int(statusF(1))
      ipe_next=int(statusF(2))
      igrid_next=int(statusF(3))
      if (int(statusF(4))/=trace_status_active) then
        stopT=.TRUE.
        if (int(statusF(4))==trace_status_weak_field) then
          wL(1)=ipoint_out
        else
          wL(1)=ipoint_out-1
        endif
      endif
      if (mype==ipe_next) then
        do j=1,ndim
          xf(ipoint_out,j)=statusF(4+j)
        enddo
      else
        xf(ipoint_out,:)=zero
      endif

      ! pe and grid of next point
      ipe_now=ipe_next
      igrid_now=igrid_next
      ipoint_in=ipoint_out
    enddo

    if (tcondi/='TRAC') then
      numRT=int(wL(1))
      allocate(data_send(numRT,ndim+nwP),data_recv(numRT,ndim+nwP))
      data_send(:,:)=zero
      data_recv(:,:)=zero
      data_send(1:numRT,1:ndim)=xf(1:numRT,1:ndim)
      if (nwP>0) data_send(1:numRT,1+ndim:ndim+nwP)=wP(1:numRT,1:nwP)
      call MPI_ALLREDUCE(data_send,data_recv,numRT*(ndim+nwP),MPI_DOUBLE_PRECISION,&
                         MPI_SUM,icomm,ierrmpi)
      xf(1:numRT,1:ndim)=data_recv(1:numRT,1:ndim)
      if (nwP>0) wP(1:numRT,1:nwP)=data_recv(1:numRT,1+ndim:ndim+nwP)
      deallocate(data_send,data_recv)
    endif

  end subroutine trace_field_single

  subroutine find_points_in_pe(igrid,ipoint_in,xf,wP,wL,dL,numP,nwP,nwL,forward,ftype,tcondi,statusF)

    integer, intent(inout) :: igrid
    integer, intent(in) :: ipoint_in,numP,nwP,nwL
    double precision, intent(inout) :: xf(numP,ndim),wP(numP,nwP),wL(1+nwL)
    double precision, intent(in) :: dL
    logical, intent(in) :: forward
    character(len=std_len), intent(in) :: ftype,tcondi
    double precision, intent(inout) :: statusF(4+ndim)

    double precision :: xfout(ndim)
    integer :: ipe_next,igrid_next,ip_in,ip_out,j,indomain,trace_status
    logical :: newpe,stopT

    ip_in=ipoint_in
    newpe=.FALSE.
    ipe_next=mype
    igrid_next=igrid
    statusF=zero

    do while(newpe .eqv. .FALSE.)
      ! looking for points in given grid    
      call find_points_interp(igrid,ip_in,ip_out,xf,wP,wL,numP,nwP,nwL, &
           dL,forward,ftype,tcondi,trace_status)
      ip_in=ip_out

      ! when next point is out of given grid, find next grid  
      if (trace_status/=trace_status_active) then
        newpe=.TRUE.
        stopT=.TRUE.
      else
        indomain=0
        {if (xf(ip_out,^DB)>=xprobmin^DB .and. xf(ip_out,^DB)<=xprobmax^DB) indomain=indomain+1\}
        if (ip_out<numP .and. indomain==ndim) then
          if (tcondi/='TRAC') then
            stopT=.FALSE.
            xfout=xf(ip_out,:)
            call find_next_grid(igrid,igrid_next,ipe_next,xfout,newpe,stopT)
            if (stopT) trace_status=trace_status_boundary
          else
            if (xf(ip_out,ndim)>phys_trac_mask) then
              newpe=.TRUE.
              stopT=.TRUE.
              trace_status=trace_status_trac_stop
            else
              stopT=.FALSE.
              xfout=xf(ip_out,:)
              call find_next_grid(igrid,igrid_next,ipe_next,xfout,newpe,stopT)
              if (stopT) trace_status=trace_status_boundary
            endif
          endif
        else
          newpe=.TRUE.
          stopT=.TRUE.
          if (ip_out>=numP) then
            trace_status=trace_status_max_steps
          else
            trace_status=trace_status_boundary
          endif
        endif
      endif

      if (newpe) then
        statusF(1)=ip_out
        statusF(2)=ipe_next
        statusF(3)=igrid_next
        statusF(4)=trace_status_active
        if (stopT) statusF(4)=trace_status
        do j=1,ndim
          statusF(4+j)=xf(ip_out,j)
        enddo
      endif

      if (newpe .eqv. .FALSE.) igrid=igrid_next
    enddo

  end subroutine find_points_in_pe

  subroutine find_next_grid(igrid,igrid_next,ipe_next,xf1,newpe,stopT)
    ! check the grid and pe of next point
    use mod_usr_methods
    use mod_global_parameters
    use mod_forest

    integer, intent(inout) :: igrid,igrid_next,ipe_next
    double precision, intent(in) :: xf1(ndim)
    logical, intent(inout) :: newpe,stopT

    double precision :: dxb^D,xb^L,xbmid^D
    double precision :: xbn^L
    integer :: idn^D,my_neighbor_type,inblock
    integer :: ic^D,inc^D,ipe_neighbor,igrid_neighbor

    igrid_next=igrid
    ipe_next=mype

    ^D&xbmin^D=rnode(rpxmin^D_,igrid)\
    ^D&xbmax^D=rnode(rpxmax^D_,igrid)\
    inblock=0

    ! direction of next grid
    idn^D=0\
    {if (xf1(^D)<=xbmin^D) idn^D=-1\}
    {if (xf1(^D)>=xbmax^D) idn^D=1\}
    my_neighbor_type=neighbor_type(idn^D,igrid)
    igrid_neighbor=neighbor(1,idn^D,igrid)
    ipe_neighbor=neighbor(2,idn^D,igrid)

    ! ipe and igrid of next grid
    select case(my_neighbor_type)
    case (neighbor_boundary)
      ! next point is not in simulation box
      newpe=.TRUE.
      stopT=.TRUE.

    case(neighbor_coarse)
      ! neighbor grid has lower refinement level      
      igrid_next=igrid_neighbor
      ipe_next=ipe_neighbor
      if (mype==ipe_neighbor) then
        newpe=.FALSE.
      else
        newpe=.TRUE.
      endif

    case(neighbor_sibling)
      ! neighbor grid has lower refinement level 
      igrid_next=igrid_neighbor
      ipe_next=ipe_neighbor
      if (mype==ipe_neighbor) then
        newpe=.FALSE.
      else
        newpe=.TRUE.
      endif

    case(neighbor_fine)
      ! neighbor grid has higher refinement level 
      {xbmid^D=(xbmin^D+xbmax^D)/2.d0\}
      ^D&inc^D=1\
      {if (xf1(^D)<=xbmin^D) inc^D=0\}
      {if (xf1(^D)>xbmin^D .and. xf1(^D)<=xbmid^D) inc^D=1\}
      {if (xf1(^D)>xbmid^D .and. xf1(^D)<xbmax^D) inc^D=2\}
      {if (xf1(^D)>=xbmax^D) inc^D=3\}
      ipe_next=neighbor_child(2,inc^D,igrid)
      igrid_next=neighbor_child(1,inc^D,igrid)
      if (mype==ipe_next) then
        newpe=.FALSE.
      else
        newpe=.TRUE.
      endif
    end select

  end subroutine find_next_grid

  subroutine find_points_interp(igrid,ip_in,ip_out,xf,wP,wL,numP,nwP,nwL, &
       dL,forward,ftype,tcondi,trace_status)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in) :: igrid,ip_in,numP,nwP,nwL
    integer, intent(inout) :: ip_out
    double precision, intent(inout) :: xf(numP,ndim),wP(numP,nwP),wL(1+nwL)
    double precision, intent(in) :: dL
    logical, intent(in) :: forward
    character(len=std_len), intent(in) :: ftype,tcondi
    integer, intent(out) :: trace_status

    double precision :: dxb^D,xb^L
    double precision :: field(ixg^T,ndir)
    double precision :: xs1(ndim),xs2(ndim),K1(ndim),K2(ndim)
    double precision :: xfpre(ndim),xfnow(ndim),xfnext(ndim)
    double precision :: Tpre,Tnow,Tnext,dTds,Lt,Lr,ds,T_bott,trac_delta
    integer          :: ip,inblock,ixI^L,ixO^L,j
    logical :: field_ok

    ixI^L=ixG^LL;
    ixO^L=ixM^LL;
    ^D&dxb^D=rnode(rpdx^D_,igrid);
    ^D&xbmin^D=rnode(rpxmin^D_,igrid);
    ^D&xbmax^D=rnode(rpxmax^D_,igrid);
    ip_out=ip_in
    trace_status=trace_status_active

    if (tcondi/='TRAC') then
      ds=dL
    else
      ds=dxb^ND
      T_bott=2.d4/unit_temperature
      trac_delta=0.25d0
    endif

    ! main loop
    MAINLOOP: do ip=ip_in,numP-1

      ! integrate magnetic field with Runge-Kutta method
      xs1(:)=xf(ip,:)
      call get_K(xs1,igrid,K1,ixI^L,dxb^D,ftype,smalldouble,field_ok)
      if (.not.field_ok) then
        ip_out=ip
        trace_status=trace_status_weak_field
        return
      endif
      if (forward) then
        xs2(:)=xf(ip,:)+ds*K1(:)
      else
        xs2(:)=xf(ip,:)-ds*K1(:)
      endif
      call get_K(xs2,igrid,K2,ixI^L,dxb^D,ftype,smalldouble,field_ok)
      if (.not.field_ok) then
        ip_out=ip
        trace_status=trace_status_weak_field
        return
      endif
      if (forward) then
        xf(ip+1,:)=xf(ip,:)+ds*(0.5*K1(:)+0.5*K2(:))
      else
        xf(ip+1,:)=xf(ip,:)-ds*(0.5*K1(:)+0.5*K2(:))
      endif
      ip_out=ip+1

      ! get local values for variable via interpolation
      if (tcondi/='TRAC') then
        if (associated(usr_set_field_w)) then 
          call usr_set_field_w(igrid,ip,xf,wP,wL,numP,nwP,nwL,dL,forward,ftype,tcondi)
        endif
      else  ! get TRAC Tcoff
        wP(ip,1)=mype
        wP(ip,2)=igrid
        if (ip==ip_in) then
          if (forward) then
            xfpre(:)=xf(ip,:)-ds*K1(:)
          else
            xfpre(:)=xf(ip,:)+ds*K1(:)
          endif
          xfnow(:)=xf(ip,:)
          xfnext(:)=xf(ip+1,:)
          call get_T_loc_TRAC(igrid,xfpre,Tpre,ixI^L,dxb^D)
          call get_T_loc_TRAC(igrid,xfnow,Tnow,ixI^L,dxb^D)
          call get_T_loc_TRAC(igrid,xfnext,Tnext,ixI^L,dxb^D)
        else
          xfpre=xf(ip-1,:)
          xfnow(:)=xf(ip,:)
          xfnext(:)=xf(ip+1,:)
          Tpre=Tnow
          Tnow=Tnext
          call get_T_loc_TRAC(igrid,xfnext,Tnext,ixI^L,dxb^D)
        endif
        dTds=abs(Tnext-Tpre)/(2*ds)
        if (ip==1) then
          Lt=0.d0
          wL(2)=T_bott   ! current Tcofl
          wL(3)=Tnow     ! current Tlmax
        else
          Lt=0.d0
          if (dTds>0.d0) then
            Lt=Tnow/dTds
            Lr=ds
            ! renew cutoff temperature
            if(Lr>trac_delta*Lt) then
              if (Tnow>wL(2)) wL(2)=Tnow
            endif
          endif
          if (Tnow>wL(3)) wL(3)=Tnow
        endif
      endif

      ! exit loop if next point is not in this block
      inblock=0
      {if (xf(ip+1,^DB)>=xbmin^DB .and. xf(ip+1,^DB)<xbmax^DB) inblock=inblock+1\}
      if (tcondi=='TRAC' .and. xf(ip+1,ndim)>phys_trac_mask) inblock=0
      if (inblock/=ndim) exit MAINLOOP

    enddo MAINLOOP

  end subroutine find_points_interp

  subroutine get_K(xfn,igrid,K,ixI^L,dxb^D,ftype,b_min,field_ok)
    use mod_usr_methods

    integer :: ixI^L,igrid
    double precision :: dxb^D
    double precision :: xfn(ndim),K(ndim)
    character(len=std_len) :: ftype
    double precision, intent(in), optional :: b_min
    logical, intent(out), optional :: field_ok

    double precision :: xd^D
    double precision :: field(0:1^D&,ndir),Fx(ndim),factor(0:1^D&)
    double precision :: Ftotal,field_min
    logical :: valid_field
    integer          :: ixb^D,ix^D,ixbl^D,j

    ^D&ixbl^D=floor((xfn(^D)-ps(igrid)%x(ixImin^DD,^D))/dxb^D)+ixImin^D;
    ^D&xd^D=(xfn(^D)-ps(igrid)%x(ixbl^DD,^D))/dxb^D;

    field=zero
    if(ftype=='Bfield') then
      if(B0field) then
        if(allocated(iw_mag)) then
          {do ix^D=0,1\}
            field(ix^D,1:ndir)=ps(igrid)%w(ixbl^D+ix^D,iw_mag(1:ndir)) &
                 +ps(igrid)%B0(ixbl^D+ix^D,1:ndir,0)
          {enddo\}
        else
          {do ix^D=0,1\}
            field(ix^D,1:ndir)=ps(igrid)%B0(ixbl^D+ix^D,1:ndir,0)
          {enddo\}
        endif
      else
        {do ix^D=0,1\}
          field(ix^D,1:ndir)=ps(igrid)%w(ixbl^D+ix^D,iw_mag(1:ndir))
        {enddo\}
      endif
    else if (ftype=='Vfield') then
      block
        double precision :: w_stencil(ixbl^D:ixbl^D+1^D&,1:nw)
        double precision :: x_stencil(ixbl^D:ixbl^D+1^D&,1:ndim)
        double precision :: vector_stencil(ixbl^D:ixbl^D+1^D&,1:ndir)
        integer :: ixS^L

        ^D&ixSmin^D=ixbl^D;
        ^D&ixSmax^D=ixbl^D+1;
        w_stencil(ixS^S,1:nw)=ps(igrid)%w(ixS^S,1:nw)
        x_stencil(ixS^S,1:ndim)=ps(igrid)%x(ixS^S,1:ndim)
        call phys_get_v(w_stencil,x_stencil,ixS^L,ixS^L,vector_stencil)
        {do ix^D=0,1\}
          field(ix^D,1:ndir)=vector_stencil(ixbl^D+ix^D,1:ndir)
        {enddo\}
      end block
    endif
    {do ix^D=0,1\}
      factor(ix^D)={abs(1-ix^D-xd^D)*}
    {enddo\}

    if (ftype=='Bfield' .or. ftype=='Vfield') then
      Fx=0.d0
      {do ix^DB=0,1\}
        do j=1,ndim
          Fx(j)=Fx(j)+field(ix^D,j)*factor(ix^D)
        enddo
      {enddo\}
    else if (associated(usr_set_field)) then
      call usr_set_field(xfn,igrid,Fx,ftype)
    else
      call MPISTOP('Field tracing error: wrong field type!')
    endif

    Ftotal=zero
    do j=1,ndim
      Ftotal=Ftotal+(Fx(j))**2
    enddo
    Ftotal=dsqrt(Ftotal)

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)
    valid_field=Ftotal>zero .and. Ftotal>=field_min
    if (present(field_ok)) field_ok=valid_field

    K=zero
    if (valid_field) then
      K(1:ndim)=Fx(1:ndim)/Ftotal
    endif

  end subroutine get_K

  subroutine get_T_loc_TRAC(igrid,xloc,Tloc,ixI^L,dxb^D)
    ! grid T has been calculated and stored in wextra(ixI^S,Tcoff_)
    ! see mod_trac
    integer, intent(in) :: igrid,ixI^L
    double precision, intent(inout) :: xloc(ndim)
    double precision, intent(inout) :: Tloc
    double precision, intent(in) :: dxb^D

    double precision :: xd^D
    double precision :: factor(0:1^D&),Tnear(0:1^D&)
    integer          :: ixb^D,ix^D,ixbl^D,j,ixO^L

    ^D&ixbl^D=floor((xloc(^D)-ps(igrid)%x(ixImin^DD,^D))/dxb^D)+ixImin^D;
    ^D&xd^D=(xloc(^D)-ps(igrid)%x(ixbl^DD,^D))/dxb^D;
    ^D&ixOmin^D=ixbl^D;
    ^D&ixOmax^D=ixOmin^D+1;

    {do ix^D=0,1\}
      factor(ix^D)={abs(1-ix^D-xd^D)*}
      Tnear(ix^D)=ps(igrid)%wextra(ixbl^D+ix^D,iw_tcoff)
    {enddo\}

    Tloc=0.d0
    ! interpolation
    {do ix^DB=0,1\}
      Tloc=Tloc+Tnear(ix^D)*factor(ix^D)
    {enddo\}

  end subroutine get_T_loc_TRAC

end module mod_trace_field
