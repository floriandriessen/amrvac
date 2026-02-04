!>   The module add viscous source terms and check time step
!>
!>   Viscous forces in the momentum equations:
!>   d m_i/dt +=  - div (vc_mu * PI)
!>   !! Viscous work in the energy equation:
!>   !! de/dt    += - div (v . vc_mu * PI)
!>   where the PI stress tensor is
!>   PI_i,j = - (dv_j/dx_i + dv_i/dx_j) + Kronecker delta_i,j*(2/3)*Sum_k dv_k/dx_k
!>   where vc_mu is the dynamic viscosity coefficient (g cm^-1 s^-1).
module mod_viscosity
  use mod_comm_lib, only: mpistop
  implicit none

  !> Viscosity coefficient
  double precision, public :: vc_mu = 1.d0

  !> Index of the density (in the w array)
  integer, private, parameter              :: rho_ = 1

  !> Indices of the momentum density
  integer, allocatable, private, protected :: mom(:)

  !> Index of the energy density (-1 if not present)
  integer, private, protected              :: e_
  !> Indices of the velocity for the form of better vectorization
  integer, public, protected              :: v1_,v2_,v3_

  !> fourth order
  logical :: vc_4th_order = .false.

  !> source split or not
  logical :: vc_split= .false.

  !> whether to compute the viscous terms as
  !> fluxes (ie in the div on the LHS), or not (by default)
  logical :: viscInDiv= .false.

  procedure(sub_add_source), pointer :: viscosity_add_source => null()
  ! Public methods
  public :: viscosity_add_source
  public :: visc_get_flux_prim

contains
  !> Read this module"s parameters from a file
  subroutine vc_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /vc_list/ vc_mu, vc_4th_order, vc_split, viscInDiv

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, vc_list, end=111)
111    close(unitpar)
    end do


  end subroutine vc_params_read

  !> Initialize the module
  subroutine viscosity_init(phys_wider_stencil)
    use mod_global_parameters
    use mod_geometry
    integer, intent(inout) :: phys_wider_stencil
    integer :: nwx,idir

    call vc_params_read(par_files)

    if(vc_split) any_source_split=.true.

    ! Determine flux variables
    nwx = 1                  ! rho (density)

    allocate(mom(ndir))
    do idir = 1, ndir
       nwx    = nwx + 1
       mom(idir) = nwx       ! momentum density
    end do
    v^C_=mom(^C);

    nwx = nwx + 1
    e_     = nwx          ! energy density

    if (viscInDiv) then
      ! to compute the derivatives from left and right upwinded values
      phys_wider_stencil = 1
    end if
    select case (coordinate)
      case (Cartesian,Cartesian_stretched,Cartesian_expansion)
        viscosity_add_source => viscosity_add_source_Cartesian
      case (cylindrical)
        viscosity_add_source => viscosity_add_source_cylinder
      case (spherical)
        viscosity_add_source => viscosity_add_source_sphere
    end select

  end subroutine viscosity_init

  subroutine viscosity_add_source_Cartesian(qdt,ixI^L,ixO^L,wCT,wp,w,x,&
       energy,qsourcesplit,active)
  ! Add viscosity source in isotropic Newtonian fluids to w within ixO
  ! neglecting bulk viscosity
  ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw), wp(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in) :: energy,qsourcesplit
    logical, intent(inout) :: active

    double precision:: lambda(ixI^S,ndir,ndir),tmp(ixI^S),vlambda(ixI^S,ndir),qdtmu,divv23
    integer:: ix^L,ix^D

    if (viscInDiv) return

    if(qsourcesplit .eqv. vc_split) then
      active = .true.
      ! standard case, textbook viscosity
      ! Calculating viscosity sources
      if(.not.vc_4th_order) then
        ! involves second derivatives, two extra layers
        ix^L=ixO^L^LADD2;
        if({ ixImin^D>ixmin^D .or. ixImax^D<ixmax^D|.or.})&
          call mpistop("error for viscous source addition, 2 layers needed")
        ix^L=ixO^L^LADD1;
      else
        ! involves second derivatives, four extra layers
        ix^L=ixO^L^LADD4;
        if({ ixImin^D>ixmin^D .or. ixImax^D<ixmax^D|.or.})&
          call mpistop("error for viscous source addition"//&
          "requested fourth order gradients: 4 layers needed")
        ix^L=ixO^L^LADD2;
      end if
      qdtmu=qdt*vc_mu
      {^IFTHREED
      {do ix^DB=ixmin^DB,ixmax^DB\}
        ! idim=1, idir=1
        lambda(ix^D,1,1)=(wp(ix1+1,ix2,ix3,v1_)-wp(ix1-1,ix2,ix3,v1_))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))
        ! idim=1, idir=2
        lambda(ix^D,1,2)=(wp(ix1+1,ix2,ix3,v2_)-wp(ix1-1,ix2,ix3,v2_))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))
        ! idim=1, idir=3
        lambda(ix^D,1,3)=(wp(ix1+1,ix2,ix3,v3_)-wp(ix1-1,ix2,ix3,v3_))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))
        ! idim=2, idir=1
        lambda(ix^D,2,1)=(wp(ix1,ix2+1,ix3,v1_)-wp(ix1,ix2-1,ix3,v1_))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))
        ! idim=2, idir=2
        lambda(ix^D,2,2)=(wp(ix1,ix2+1,ix3,v2_)-wp(ix1,ix2-1,ix3,v2_))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))
        ! idim=2, idir=3
        lambda(ix^D,2,3)=(wp(ix1,ix2+1,ix3,v3_)-wp(ix1,ix2-1,ix3,v3_))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))
        ! idim=3, idir=1
        lambda(ix^D,3,1)=(wp(ix1,ix2,ix3+1,v1_)-wp(ix1,ix2,ix3-1,v1_))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))
        ! idim=3, idir=2
        lambda(ix^D,3,2)=(wp(ix1,ix2,ix3+1,v2_)-wp(ix1,ix2,ix3-1,v2_))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))
        ! idim=3, idir=3
        lambda(ix^D,3,3)=(wp(ix1,ix2,ix3+1,v3_)-wp(ix1,ix2,ix3-1,v3_))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))
        ! dv_i/d_j + dv_j/d_i
        lambda(ix^D,1,2)=lambda(ix^D,1,2)+lambda(ix^D,2,1)
        lambda(ix^D,2,1)=lambda(ix^D,1,2)
        lambda(ix^D,1,3)=lambda(ix^D,1,3)+lambda(ix^D,3,1)
        lambda(ix^D,3,1)=lambda(ix^D,1,3)
        lambda(ix^D,2,3)=lambda(ix^D,2,3)+lambda(ix^D,3,2)
        lambda(ix^D,3,2)=lambda(ix^D,2,3)
        divv23=two*third*(lambda(ix^D,1,1)+lambda(ix^D,2,2)+lambda(ix^D,3,3))
        lambda(ix^D,1,1)=two*lambda(ix^D,1,1)-divv23
        lambda(ix^D,2,2)=two*lambda(ix^D,2,2)-divv23
        lambda(ix^D,3,3)=two*lambda(ix^D,3,3)-divv23
        if(energy) then
          vlambda(ix^D,1)=wp(ix^D,v1_)*lambda(ix^D,1,1)+wp(ix^D,v2_)*lambda(ix^D,2,1)+wp(ix^D,v3_)*lambda(ix^D,3,1)
          vlambda(ix^D,2)=wp(ix^D,v1_)*lambda(ix^D,1,2)+wp(ix^D,v2_)*lambda(ix^D,2,2)+wp(ix^D,v3_)*lambda(ix^D,3,2)
          vlambda(ix^D,3)=wp(ix^D,v1_)*lambda(ix^D,1,3)+wp(ix^D,v2_)*lambda(ix^D,2,3)+wp(ix^D,v3_)*lambda(ix^D,3,3)
        end if
      {end do\}
      ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr)
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        w(ix^D,mom(1))=((lambda(ix1+1,ix2,ix3,1,1)-lambda(ix1-1,ix2,ix3,1,1))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))+&
                        (lambda(ix1,ix2+1,ix3,2,1)-lambda(ix1,ix2-1,ix3,2,1))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))+&
                        (lambda(ix1,ix2,ix3+1,3,1)-lambda(ix1,ix2,ix3-1,3,1))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3)))*qdtmu+w(ix^D,mom(1))
        w(ix^D,mom(2))=((lambda(ix1+1,ix2,ix3,1,2)-lambda(ix1-1,ix2,ix3,1,2))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))+&
                        (lambda(ix1,ix2+1,ix3,2,2)-lambda(ix1,ix2-1,ix3,2,2))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))+&
                        (lambda(ix1,ix2,ix3+1,3,2)-lambda(ix1,ix2,ix3-1,3,2))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3)))*qdtmu+w(ix^D,mom(2))
        w(ix^D,mom(3))=((lambda(ix1+1,ix2,ix3,1,3)-lambda(ix1-1,ix2,ix3,1,3))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))+&
                        (lambda(ix1,ix2+1,ix3,2,3)-lambda(ix1,ix2-1,ix3,2,3))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))+&
                        (lambda(ix1,ix2,ix3+1,3,3)-lambda(ix1,ix2,ix3-1,3,3))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3)))*qdtmu+w(ix^D,mom(3))
      {end do\}
      }
      {^IFTWOD
      {do ix^DB=ixmin^DB,ixmax^DB\}
        ! idim=1, idir=1
        lambda(ix^D,1,1)=(wp(ix1+1,ix2,v1_)-wp(ix1-1,ix2,v1_))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))
        ! idim=1, idir=2
        lambda(ix^D,1,2)=(wp(ix1+1,ix2,v2_)-wp(ix1-1,ix2,v2_))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))
        ! idim=2, idir=1
        lambda(ix^D,2,1)=(wp(ix1,ix2+1,v1_)-wp(ix1,ix2-1,v1_))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))
        ! idim=2, idir=2
        lambda(ix^D,2,2)=(wp(ix1,ix2+1,v2_)-wp(ix1,ix2-1,v2_))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))
        ! dv_i/d_j + dv_j/d_i
        lambda(ix^D,1,2)=lambda(ix^D,1,2)+lambda(ix^D,2,1)
        lambda(ix^D,2,1)=lambda(ix^D,1,2)
        divv23=two*third*(lambda(ix^D,1,1)+lambda(ix^D,2,2))
        lambda(ix^D,1,1)=two*lambda(ix^D,1,1)-divv23
        lambda(ix^D,2,2)=two*lambda(ix^D,2,2)-divv23
        if(ndir==3) then
          ! idim=1, idir=3
          lambda(ix1,ix2,1,3)=(wp(ix1+1,ix2,v3_)-wp(ix1-1,ix2,v3_))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))
          ! idim=2, idir=3
          lambda(ix1,ix2,2,3)=(wp(ix1,ix2+1,v3_)-wp(ix1,ix2-1,v3_))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))
          lambda(ix1,ix2,3,1)=lambda(ix1,ix2,1,3)
          lambda(ix1,ix2,3,2)=lambda(ix1,ix2,2,3)
          lambda(ix1,ix2,3,3)=-divv23
          if(energy) then
            vlambda(ix1,ix2,1)=wp(ix1,ix2,v1_)*lambda(ix1,ix2,1,1)+wp(ix1,ix2,v2_)*lambda(ix1,ix2,2,1)+wp(ix1,ix2,v3_)*lambda(ix1,ix2,3,1)
            vlambda(ix1,ix2,2)=wp(ix1,ix2,v1_)*lambda(ix1,ix2,1,2)+wp(ix1,ix2,v2_)*lambda(ix1,ix2,2,2)+wp(ix1,ix2,v3_)*lambda(ix1,ix2,3,2)
            vlambda(ix1,ix2,3)=wp(ix1,ix2,v1_)*lambda(ix1,ix2,1,3)+wp(ix1,ix2,v2_)*lambda(ix1,ix2,2,3)+wp(ix1,ix2,v3_)*lambda(ix1,ix2,3,3)
          end if
        else if(energy) then
          vlambda(ix1,ix2,1)=wp(ix1,ix2,v1_)*lambda(ix1,ix2,1,1)+wp(ix1,ix2,v2_)*lambda(ix1,ix2,2,1)
          vlambda(ix1,ix2,2)=wp(ix1,ix2,v1_)*lambda(ix1,ix2,1,2)+wp(ix1,ix2,v2_)*lambda(ix1,ix2,2,2)
        end if
      {end do\}
      ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr)
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        w(ix^D,mom(1))=((lambda(ix1+1,ix2,1,1)-lambda(ix1-1,ix2,1,1))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                        (lambda(ix1,ix2+1,2,1)-lambda(ix1,ix2-1,2,1))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2)))*qdtmu+w(ix^D,mom(1))
        w(ix^D,mom(2))=((lambda(ix1+1,ix2,1,2)-lambda(ix1-1,ix2,1,2))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                        (lambda(ix1,ix2+1,2,2)-lambda(ix1,ix2-1,2,2))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2)))*qdtmu+w(ix^D,mom(2))
        if(ndir==3) then
          w(ix1,ix2,mom(3))=((lambda(ix1+1,ix2,1,3)-lambda(ix1-1,ix2,1,3))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                             (lambda(ix1,ix2+1,2,3)-lambda(ix1,ix2-1,2,3))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2)))*qdtmu+w(ix1,ix2,mom(3))
        end if
      {end do\}
      }
      {^IFONED
      do ix1=ixmin1,ixmax1
        ! idim=1, idir=1
        lambda(ix1,1,1)=(wp(ix1+1,v1_)-wp(ix1-1,v1_))/(x(ix1+1,1)-x(ix1-1,1))
        ! dv_i/d_j + dv_j/d_i
        divv23=two*third*lambda(ix1,1,1)
        lambda(ix1,1,1)=two*lambda(ix1,1,1)-divv23
        if(ndir==2) then
          ! idim=1, idir=2
          lambda(ix1,1,2)=(wp(ix1+1,v2_)-wp(ix1-1,v2_))/(x(ix1+1,1)-x(ix1-1,1))
          ! dv_i/d_j + dv_j/d_i
          lambda(ix1,2,1)=lambda(ix1,1,2)
          lambda(ix1,2,2)=-divv23
          if(energy) then
            vlambda(ix1,1)=wp(ix1,v1_)*lambda(ix1,1,1)+wp(ix1,v2_)*lambda(ix1,2,1)
            vlambda(ix1,2)=wp(ix1,v1_)*lambda(ix1,1,2)+wp(ix1,v2_)*lambda(ix1,2,2)
          end if
        else if(ndir==3) then
          ! idim=1, idir=2
          lambda(ix1,1,2)=(wp(ix1+1,v2_)-wp(ix1-1,v2_))/(x(ix1+1,1)-x(ix1-1,1))
          ! dv_i/d_j + dv_j/d_i
          lambda(ix1,2,1)=lambda(ix1,1,2)
          lambda(ix1,2,2)=-divv23
          ! idim=1, idir=3
          lambda(ix1,1,3)=(wp(ix1+1,v3_)-wp(ix1-1,v3_))/(x(ix1+1,1)-x(ix1-1,1))
          lambda(ix1,3,1)=lambda(ix1,1,3)
          lambda(ix1,3,3)=-divv23
          if(energy) then
            vlambda(ix1,1)=wp(ix1,v1_)*lambda(ix1,1,1)+wp(ix1,v2_)*lambda(ix1,2,1)+wp(ix1,v3_)*lambda(ix1,3,1)
            vlambda(ix1,2)=wp(ix1,v1_)*lambda(ix1,1,2)+wp(ix1,v2_)*lambda(ix1,2,2)+wp(ix1,v3_)*lambda(ix1,3,2)
            vlambda(ix1,3)=wp(ix1,v1_)*lambda(ix1,1,3)+wp(ix1,v2_)*lambda(ix1,2,3)+wp(ix1,v3_)*lambda(ix1,3,3)
          end if
        else if(energy) then
          vlambda(ix1,1)=wp(ix1,v1_)*lambda(ix1,1,1)
        end if
      end do
      ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr)
      do ix1=ixOmin1,ixOmax1
        if(ndir==1) then
          w(ix1,mom(1))=(lambda(ix1+1,1,1)-lambda(ix1-1,1,1))/(x(ix1+1,1)-x(ix1-1,1))*qdtmu+w(ix1,mom(1))
        else if(ndir==2) then
          w(ix1,mom(1))=(lambda(ix1+1,1,1)-lambda(ix1-1,1,1))/(x(ix1+1,1)-x(ix1-1,1))*qdtmu+w(ix1,mom(1))
          w(ix1,mom(2))=(lambda(ix1+1,1,2)-lambda(ix1-1,1,2))/(x(ix1+1,1)-x(ix1-1,1))*qdtmu+w(ix1,mom(2))
        else
          w(ix1,mom(1))=(lambda(ix1+1,1,1)-lambda(ix1-1,1,1))/(x(ix1+1,1)-x(ix1-1,1))*qdtmu+w(ix1,mom(1))
          w(ix1,mom(2))=(lambda(ix1+1,1,2)-lambda(ix1-1,1,2))/(x(ix1+1,1)-x(ix1-1,1))*qdtmu+w(ix1,mom(2))
          w(ix1,mom(3))=(lambda(ix1+1,1,3)-lambda(ix1-1,1,3))/(x(ix1+1,1)-x(ix1-1,1))*qdtmu+w(ix1,mom(3))
        end if
      end do
      }
      if(energy) then
        ! de/dt= +div(v.dot.[mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v *kr])
        ! thus e=e+d_i v_j tensor_ji
        call divvector(vlambda,ixI^L,ixO^L,tmp)
        w(ixO^S,e_)=w(ixO^S,e_)+tmp(ixO^S)*qdtmu
      end if
    end if

  end subroutine viscosity_add_source_Cartesian

  subroutine viscosity_add_source_sphere(qdt,ixI^L,ixO^L,wCT,wp,w,x,&
       energy,qsourcesplit,active)
  ! Add viscosity source in isotropic Newtonian fluids to w within ixO
  ! neglecting bulk viscosity
  ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw), wp(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in) :: energy,qsourcesplit
    logical, intent(inout) :: active

    double precision :: lambda(ixI^S,ndir,ndir),vlambda(ixI^S,ndir),invr(ixI^S),tctan(ixI^S),invrsin(ixI^S)
    double precision :: qdtmu,tsin,tcos,divv23
    integer:: ix^L,ix^D

    if (viscInDiv) return

    if(qsourcesplit .eqv. vc_split) then
      active = .true.
      ! standard case, textbook viscosity
      ! Calculating viscosity sources
      if(.not.vc_4th_order) then
        ! involves second derivatives, two extra layers
        ix^L=ixO^L^LADD2;
        if({ ixImin^D>ixmin^D .or. ixImax^D<ixmax^D|.or.})&
          call mpistop("error for viscous source addition, 2 layers needed")
        ix^L=ixO^L^LADD1;
      else
        ! involves second derivatives, four extra layers
        ix^L=ixO^L^LADD4;
        if({ ixImin^D>ixmin^D .or. ixImax^D<ixmax^D|.or.})&
          call mpistop("error for viscous source addition"//&
          "requested fourth order gradients: 4 layers needed")
        ix^L=ixO^L^LADD2;
      end if

      ! construct lambda tensor: lambda_ij = gradv_ij + gradv_ji
      ! initialize
      qdtmu=qdt*vc_mu
      {^IFTHREED
      {do ix^DB=ixmin^DB,ixmax^DB\}
        invr(ix^D)=1.d0/x(ix^D,1)
        tcos=dcos(x(ix^D,2))
        ! idim=1, idir=1
        lambda(ix^D,1,1)=(wp(ix1+1,ix2,ix3,v1_)-wp(ix1-1,ix2,ix3,v1_))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))*two*qdtmu
        ! idim=1, idir=2
        lambda(ix^D,1,2)=(wp(ix1+1,ix2,ix3,v2_)-wp(ix1-1,ix2,ix3,v2_))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))
        ! idim=1, idir=3
        lambda(ix^D,1,3)=(wp(ix1+1,ix2,ix3,v3_)-wp(ix1-1,ix2,ix3,v3_))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))
        ! idim=2, idir=1
        lambda(ix^D,2,1)=((wp(ix1,ix2+1,ix3,v1_)-wp(ix1,ix2-1,ix3,v1_))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))-wp(ix^D,v2_))*invr(ix^D)
        ! idim=2, idir=2
        lambda(ix^D,2,2)=((wp(ix1,ix2+1,ix3,v2_)-wp(ix1,ix2-1,ix3,v2_))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))+wp(ix^D,v1_))*invr(ix^D)*&
                         two*qdtmu
        ! idim=2, idir=3
        lambda(ix^D,2,3)=((wp(ix1,ix2+1,ix3,v3_)-wp(ix1,ix2-1,ix3,v3_))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))+wp(ix^D,v2_)*tcos)*invr(ix^D)
        tsin=dsin(x(ix^D,2))
        invrsin(ix^D)=invr(ix^D)/tsin
        tctan(ix^D)=tcos/tsin
        ! idim=3, idir=1
        lambda(ix^D,3,1)=((wp(ix1,ix2,ix3+1,v1_)-wp(ix1,ix2,ix3-1,v1_))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))-wp(ix^D,v3_)*tsin)*invrsin(ix^D)
        ! idim=3, idir=2
        lambda(ix^D,3,2)=((wp(ix1,ix2,ix3+1,v2_)-wp(ix1,ix2,ix3-1,v2_))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))-wp(ix^D,v3_)*tcos)*invrsin(ix^D)
        ! idim=3, idir=3
        lambda(ix^D,3,3)=((wp(ix1,ix2,ix3+1,v3_)-wp(ix1,ix2,ix3-1,v3_))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))+wp(ix^D,v1_)*tsin+&
            wp(ix^D,v2_)*tcos)*invrsin(ix^D)*two*qdtmu
        ! dv_i/d_j + dv_j/d_i
        lambda(ix^D,1,2)=(lambda(ix^D,1,2)+lambda(ix^D,2,1))*qdtmu
        lambda(ix^D,2,1)=lambda(ix^D,1,2)
        lambda(ix^D,1,3)=(lambda(ix^D,1,3)+lambda(ix^D,3,1))*qdtmu
        lambda(ix^D,3,1)=lambda(ix^D,1,3)
        lambda(ix^D,2,3)=(lambda(ix^D,2,3)+lambda(ix^D,3,2))*qdtmu
        lambda(ix^D,3,2)=lambda(ix^D,2,3)
        divv23=third*(lambda(ix^D,1,1)+lambda(ix^D,2,2)+lambda(ix^D,3,3))
        lambda(ix^D,1,1)=lambda(ix^D,1,1)-divv23
        lambda(ix^D,2,2)=lambda(ix^D,2,2)-divv23
        lambda(ix^D,3,3)=lambda(ix^D,3,3)-divv23
        if(energy) then
          vlambda(ix^D,1)=wp(ix^D,v1_)*lambda(ix^D,1,1)+wp(ix^D,v2_)*lambda(ix^D,2,1)+wp(ix^D,v3_)*lambda(ix^D,3,1)
          vlambda(ix^D,2)=wp(ix^D,v1_)*lambda(ix^D,1,2)+wp(ix^D,v2_)*lambda(ix^D,2,2)+wp(ix^D,v3_)*lambda(ix^D,3,2)
          vlambda(ix^D,3)=wp(ix^D,v1_)*lambda(ix^D,1,3)+wp(ix^D,v2_)*lambda(ix^D,2,3)+wp(ix^D,v3_)*lambda(ix^D,3,3)
        end if
      {end do\}
      ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr)
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        w(ix^D,mom(1))=(lambda(ix1+1,ix2,ix3,1,1)-lambda(ix1-1,ix2,ix3,1,1))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))+&
                       (lambda(ix1,ix2+1,ix3,2,1)-lambda(ix1,ix2-1,ix3,2,1))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))*invr(ix^D)+&
                       (lambda(ix1,ix2,ix3+1,3,1)-lambda(ix1,ix2,ix3-1,3,1))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))*invrsin(ix^D)+&
                    two*lambda(ix^D,1,1)*invr(ix^D)+lambda(ix^D,2,1)*invr(ix^D)*tctan(ix^D)+w(ix^D,mom(1))
        w(ix^D,mom(2))=(lambda(ix1+1,ix2,ix3,1,2)-lambda(ix1-1,ix2,ix3,1,2))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))+&
                       (lambda(ix1,ix2+1,ix3,2,2)-lambda(ix1,ix2-1,ix3,2,2))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))*invr(ix^D)+&
                       (lambda(ix1,ix2,ix3+1,3,2)-lambda(ix1,ix2,ix3-1,3,2))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))*invrsin(ix^D)+&
                    two*lambda(ix^D,1,2)*invr(ix^D)+lambda(ix^D,2,2)*invr(ix^D)*tctan(ix^D)-lambda(ix^D,3,3)*tctan(ix^D)*invr(ix^D)**2+&
                    w(ix^D,mom(2))
        w(ix^D,mom(3))=(lambda(ix1+1,ix2,ix3,1,3)-lambda(ix1-1,ix2,ix3,1,3))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))+&
                       (lambda(ix1,ix2+1,ix3,2,3)-lambda(ix1,ix2-1,ix3,2,3))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))*invr(ix^D)+&
                       (lambda(ix1,ix2,ix3+1,3,3)-lambda(ix1,ix2,ix3-1,3,3))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))*invrsin(ix^D)+&
                    two*lambda(ix^D,1,3)*invr(ix^D)+lambda(ix^D,2,3)*(invr(ix^D)+invr(ix^D)**2)*tctan(ix^D)+w(ix^D,mom(3))
      {end do\}
      }
      {^IFTWOD
      {do ix^DB=ixmin^DB,ixmax^DB\}
        invr(ix1,ix2)=1.d0/x(ix1,ix2,1)
        ! idim=1, idir=1
        lambda(ix1,ix2,1,1)=(wp(ix1+1,ix2,v1_)-wp(ix1-1,ix2,v1_))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))*two*qdtmu
        ! idim=1, idir=2
        lambda(ix1,ix2,1,2)=(wp(ix1+1,ix2,v2_)-wp(ix1-1,ix2,v2_))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))
        ! idim=2, idir=1
        lambda(ix1,ix2,2,1)=((wp(ix1,ix2+1,v1_)-wp(ix1,ix2-1,v1_))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))-wp(ix1,ix2,v2_))*invr(ix1,ix2)
        ! idim=2, idir=2
        lambda(ix1,ix2,2,2)=((wp(ix1,ix2+1,v2_)-wp(ix1,ix2-1,v2_))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))+wp(ix1,ix2,v1_))*invr(ix1,ix2)*&
                         two*qdtmu
        ! dv_i/d_j + dv_j/d_i
        lambda(ix1,ix2,1,2)=(lambda(ix1,ix2,1,2)+lambda(ix1,ix2,2,1))*qdtmu
        lambda(ix1,ix2,2,1)=lambda(ix1,ix2,1,2)
        divv23=third*(lambda(ix1,ix2,1,1)+lambda(ix1,ix2,2,2))
        tcos=dcos(x(ix1,ix2,2))
        tsin=dsin(x(ix1,ix2,2))
        tctan(ix1,ix2)=tcos/tsin
        if(ndir==2) then
          lambda(ix1,ix2,1,1)=lambda(ix1,ix2,1,1)-divv23
          lambda(ix1,ix2,2,2)=lambda(ix1,ix2,2,2)-divv23
          if(energy) then
            vlambda(ix1,ix2,1)=wp(ix1,ix2,v1_)*lambda(ix1,ix2,1,1)+wp(ix1,ix2,v2_)*lambda(ix1,ix2,2,1)
            vlambda(ix1,ix2,2)=wp(ix1,ix2,v1_)*lambda(ix1,ix2,1,2)+wp(ix1,ix2,v2_)*lambda(ix1,ix2,2,2)
          end if
        else
          invrsin(ix1,ix2)=invr(ix1,ix2)/tsin
          ! idim=1, idir=3
          lambda(ix1,ix2,1,3)=(wp(ix1+1,ix2,v3_)-wp(ix1-1,ix2,v3_))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))
          ! idim=2, idir=3
          lambda(ix1,ix2,2,3)=((wp(ix1,ix2+1,v3_)-wp(ix1,ix2-1,v3_))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))+wp(ix1,ix2,v2_)*tcos)*invr(ix1,ix2)
          ! idim=3, idir=1
          lambda(ix1,ix2,3,1)=-wp(ix1,ix2,v3_)*tsin*invrsin(ix1,ix2)
          ! idim=3, idir=2
          lambda(ix1,ix2,3,2)=-wp(ix1,ix2,v3_)*tcos*invrsin(ix1,ix2)
          ! idim=3, idir=3
          lambda(ix1,ix2,3,3)=(wp(ix1,ix2,v1_)*tsin+wp(ix1,ix2,v2_)*tcos)*invrsin(ix1,ix2)*two*qdtmu
          lambda(ix1,ix2,1,3)=(lambda(ix1,ix2,1,3)+lambda(ix1,ix2,3,1))*qdtmu
          lambda(ix1,ix2,3,1)=lambda(ix1,ix2,1,3)
          lambda(ix1,ix2,2,3)=(lambda(ix1,ix2,2,3)+lambda(ix1,ix2,3,2))*qdtmu
          lambda(ix1,ix2,3,2)=lambda(ix1,ix2,2,3)
          lambda(ix1,ix2,1,1)=lambda(ix1,ix2,1,1)-divv23
          lambda(ix1,ix2,2,2)=lambda(ix1,ix2,2,2)-divv23
          lambda(ix1,ix2,3,3)=-divv23
          if(energy) then
            vlambda(ix1,ix2,1)=wp(ix1,ix2,v1_)*lambda(ix1,ix2,1,1)+wp(ix1,ix2,v2_)*lambda(ix1,ix2,2,1)+wp(ix1,ix2,v3_)*lambda(ix1,ix2,3,1)
            vlambda(ix1,ix2,2)=wp(ix1,ix2,v1_)*lambda(ix1,ix2,1,2)+wp(ix1,ix2,v2_)*lambda(ix1,ix2,2,2)+wp(ix1,ix2,v3_)*lambda(ix1,ix2,3,2)
            vlambda(ix1,ix2,3)=wp(ix1,ix2,v1_)*lambda(ix1,ix2,1,3)+wp(ix1,ix2,v2_)*lambda(ix1,ix2,2,3)+wp(ix1,ix2,v3_)*lambda(ix1,ix2,3,3)
          end if
        end if
      {end do\}
      ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr)
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        if(ndir==2) then
          w(ix1,ix2,mom(1))=(lambda(ix1+1,ix2,1,1)-lambda(ix1-1,ix2,1,1))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                            (lambda(ix1,ix2+1,2,1)-lambda(ix1,ix2-1,2,1))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))*invr(ix1,ix2)+&
                         two*lambda(ix1,ix2,1,1)*invr(ix1,ix2)+lambda(ix1,ix2,2,1)*invr(ix1,ix2)*tctan(ix1,ix2)+w(ix1,ix2,mom(1))
          w(ix1,ix2,mom(2))=(lambda(ix1+1,ix2,1,2)-lambda(ix1-1,ix2,1,2))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                         (lambda(ix1,ix2+1,2,2)-lambda(ix1,ix2-1,2,2))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))*invr(ix1,ix2)+&
                      two*lambda(ix1,ix2,1,2)*invr(ix1,ix2)+lambda(ix1,ix2,2,2)*invr(ix1,ix2)*tctan(ix1,ix2)+w(ix1,ix2,mom(2))
        else
          w(ix1,ix2,mom(1))=(lambda(ix1+1,ix2,1,1)-lambda(ix1-1,ix2,1,1))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                         (lambda(ix1,ix2+1,2,1)-lambda(ix1,ix2-1,2,1))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))*invr(ix1,ix2)+&
                      two*lambda(ix1,ix2,1,1)*invr(ix1,ix2)+lambda(ix1,ix2,2,1)*invr(ix1,ix2)*tctan(ix1,ix2)+w(ix1,ix2,mom(1))
          w(ix1,ix2,mom(2))=(lambda(ix1+1,ix2,1,2)-lambda(ix1-1,ix2,1,2))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                         (lambda(ix1,ix2+1,2,2)-lambda(ix1,ix2-1,2,2))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))*invr(ix1,ix2)+&
                      two*lambda(ix1,ix2,1,2)*invr(ix1,ix2)+lambda(ix1,ix2,2,2)*invr(ix1,ix2)*tctan(ix1,ix2)-&
                      lambda(ix1,ix2,3,3)*tctan(ix1,ix2)*invr(ix1,ix2)**2+w(ix1,ix2,mom(2))
          w(ix1,ix2,mom(3))=(lambda(ix1+1,ix2,1,3)-lambda(ix1-1,ix2,1,3))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                         (lambda(ix1,ix2+1,2,3)-lambda(ix1,ix2-1,2,3))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))*invr(ix1,ix2)+&
                      two*lambda(ix1,ix2,1,3)*invr(ix1,ix2)+lambda(ix1,ix2,2,3)*(invr(ix1,ix2)+invr(ix1,ix2)**2)*tctan(ix1,ix2)+w(ix1,ix2,mom(3))
        end if
      {end do\}
      }
      {^IFONED
      do ix1=ixmin1,ixmax1
        invr(ix1)=1.d0/x(ix1,1)
        ! idim=1, idir=1
        lambda(ix1,1,1)=(wp(ix1+1,v1_)-wp(ix1-1,v1_))/(x(ix1+1,1)-x(ix1-1,1))*two*qdtmu
        ! dv_i/d_j + dv_j/d_i
        divv23=third*lambda(ix1,1,1)
        lambda(ix1,1,1)=lambda(ix1,1,1)-divv23
        if(ndir==1) then
          if(energy) then
            vlambda(ix1,1)=wp(ix1,v1_)*lambda(ix1,1,1)
          end if
        else if(ndir==2) then
          ! idim=1, idir=2
          lambda(ix1,1,2)=(wp(ix1+1,v2_)-wp(ix1-1,v2_))/(x(ix1+1,1)-x(ix1-1,1))
          lambda(ix1,1,2)=lambda(ix1,1,2)*qdtmu
          lambda(ix1,2,1)=lambda(ix1,1,2)
          lambda(ix1,2,2)=-divv23
          if(energy) then
            vlambda(ix1,1)=wp(ix1,v1_)*lambda(ix1,1,1)+wp(ix1,v2_)*lambda(ix1,2,1)
            vlambda(ix1,2)=wp(ix1,v1_)*lambda(ix1,1,2)+wp(ix1,v2_)*lambda(ix1,2,2)
          end if
        else
          tcos=0.d0
          tsin=1.d0
          tctan(ix1)=tcos/tsin
          invrsin(ix1)=invr(ix1)/tsin
          ! idim=1, idir=3
          lambda(ix1,1,3)=(wp(ix1+1,v3_)-wp(ix1-1,v3_))/(x(ix1+1,1)-x(ix1-1,1))
          ! idim=3, idir=1
          lambda(ix1,3,1)=-wp(ix1,v3_)*tsin*invrsin(ix1)
          ! idim=3, idir=2
          lambda(ix1,3,2)=-wp(ix1,v3_)*tcos*invrsin(ix1)
          ! idim=3, idir=3
          lambda(ix1,3,3)=(wp(ix1,v1_)*tsin+wp(ix1,v2_)*tcos)*invrsin(ix1)*two*qdtmu
          lambda(ix1,1,3)=(lambda(ix1,1,3)+lambda(ix1,3,1))*qdtmu
          lambda(ix1,3,1)=lambda(ix1,1,3)
          lambda(ix1,3,2)=lambda(ix1,3,2)*qdtmu
          lambda(ix1,2,3)=lambda(ix1,3,2)
          lambda(ix1,3,3)=lambda(ix1,3,3)-divv23
          if(energy) then
            vlambda(ix1,1)=wp(ix1,v1_)*lambda(ix1,1,1)+wp(ix1,v2_)*lambda(ix1,2,1)+wp(ix1,v3_)*lambda(ix1,3,1)
            vlambda(ix1,2)=wp(ix1,v1_)*lambda(ix1,1,2)+wp(ix1,v2_)*lambda(ix1,2,2)+wp(ix1,v3_)*lambda(ix1,3,2)
            vlambda(ix1,3)=wp(ix1,v1_)*lambda(ix1,1,3)+wp(ix1,v2_)*lambda(ix1,2,3)+wp(ix1,v3_)*lambda(ix1,3,3)
          end if
        end if
      end do
      ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr)
      do ix1=ixOmin1,ixOmax1
        if(ndir==1) then
          w(ix1,mom(1))=(lambda(ix1+1,1,1)-lambda(ix1-1,1,1))/(x(ix1+1,1)-x(ix1-1,1))+&
                      two*lambda(ix1,1,1)*invr(ix1)+w(ix1,mom(1))
        else if(ndir==2) then
          w(ix1,mom(1))=(lambda(ix1+1,1,1)-lambda(ix1-1,1,1))/(x(ix1+1,1)-x(ix1-1,1))+&
                      two*lambda(ix1,1,1)*invr(ix1)+lambda(ix1,2,1)*invr(ix1)*tctan(ix1)+w(ix1,mom(1))
          w(ix1,mom(2))=(lambda(ix1+1,1,2)-lambda(ix1-1,1,2))/(x(ix1+1,1)-x(ix1-1,1))+&
                      two*lambda(ix1,1,2)*invr(ix1)+lambda(ix1,2,2)*invr(ix1)*tctan(ix1)+w(ix1,mom(2))
        else
          w(ix1,mom(1))=(lambda(ix1+1,1,1)-lambda(ix1-1,1,1))/(x(ix1+1,1)-x(ix1-1,1))+&
                      two*lambda(ix1,1,1)*invr(ix1)+lambda(ix1,2,1)*invr(ix1)*tctan(ix1)+w(ix1,mom(1))
          w(ix1,mom(2))=(lambda(ix1+1,1,2)-lambda(ix1-1,1,2))/(x(ix1+1,1)-x(ix1-1,1))+&
                      two*lambda(ix1,1,2)*invr(ix1)+lambda(ix1,2,2)*invr(ix1)*tctan(ix1)-lambda(ix1,3,3)*tctan(ix1)*invr(ix1)**2+&
                      w(ix1,mom(2))
          w(ix1,mom(3))=(lambda(ix1+1,1,3)-lambda(ix1-1,1,3))/(x(ix1+1,1)-x(ix1-1,1))+&
                      two*lambda(ix1,1,3)*invr(ix1)+lambda(ix1,2,3)*(invr(ix1)+invr(ix1)**2)*tctan(ix1)+w(ix1,mom(3))
        end if
      end do
      }
      if(energy) then
        ! de/dt= +div(v.dot.[mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v *kr])
        ! thus e=e+d_i v_j tensor_ji
        call divvector(vlambda,ixI^L,ixO^L,invr)
        w(ixO^S,e_)=w(ixO^S,e_)+invr(ixO^S)
      end if

    end if

  end subroutine viscosity_add_source_sphere

  subroutine viscosity_add_source_cylinder(qdt,ixI^L,ixO^L,wCT,wp,w,x,&
       energy,qsourcesplit,active)
  ! Add viscosity source in isotropic Newtonian fluids to w within ixO
  ! neglecting bulk viscosity
  ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw), wp(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in) :: energy,qsourcesplit
    logical, intent(inout) :: active

    double precision:: lambda(ixI^S,ndir,ndir),vlambda(ixI^S,ndir),invr(ixI^S)
    double precision :: qdtmu,divv23
    integer:: ix^L,ix^D

    if (viscInDiv) return

    if(qsourcesplit .eqv. vc_split) then
      active = .true.
      ! standard case, textbook viscosity
      ! Calculating viscosity sources
      if(.not.vc_4th_order) then
        ! involves second derivatives, two extra layers
        ix^L=ixO^L^LADD2;
        if({ ixImin^D>ixmin^D .or. ixImax^D<ixmax^D|.or.})&
          call mpistop("error for viscous source addition, 2 layers needed")
        ix^L=ixO^L^LADD1;
      else
        ! involves second derivatives, four extra layers
        ix^L=ixO^L^LADD4;
        if({ ixImin^D>ixmin^D .or. ixImax^D<ixmax^D|.or.})&
          call mpistop("error for viscous source addition"//&
          "requested fourth order gradients: 4 layers needed")
        ix^L=ixO^L^LADD2;
      end if

      ! construct lambda tensor: lambda_ij = gradv_ij + gradv_ji
      ! initialize
      qdtmu=qdt*vc_mu
      {^IFTHREED
      {do ix^DB=ixmin^DB,ixmax^DB\}
        invr(ix^D)=1.d0/x(ix^D,1)
        ! idim=1, idir=1
        lambda(ix^D,1,1)=(wp(ix1+1,ix2,ix3,v1_)-wp(ix1-1,ix2,ix3,v1_))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))*two*qdtmu
        ! idim=1, idir=2
        lambda(ix^D,1,2)=(wp(ix1+1,ix2,ix3,v2_)-wp(ix1-1,ix2,ix3,v2_))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))
        ! idim=1, idir=3
        lambda(ix^D,1,3)=(wp(ix1+1,ix2,ix3,v3_)-wp(ix1-1,ix2,ix3,v3_))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))
        ! idim=2, idir=1
        lambda(ix^D,2,1)=((wp(ix1,ix2+1,ix3,v1_)-wp(ix1,ix2-1,ix3,v1_))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))-wp(ix^D,v2_))*invr(ix^D)
        ! idim=2, idir=2
        lambda(ix^D,2,2)=((wp(ix1,ix2+1,ix3,v2_)-wp(ix1,ix2-1,ix3,v2_))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))+wp(ix^D,v1_))*invr(ix^D)*&
                         two*qdtmu
        ! idim=2, idir=3
        lambda(ix^D,2,3)=(wp(ix1,ix2+1,ix3,v3_)-wp(ix1,ix2-1,ix3,v3_))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))*invr(ix^D)
        ! idim=3, idir=1
        lambda(ix^D,3,1)=(wp(ix1,ix2,ix3+1,v1_)-wp(ix1,ix2,ix3-1,v1_))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))*invr(ix^D)
        ! idim=3, idir=2
        lambda(ix^D,3,2)=(wp(ix1,ix2,ix3+1,v2_)-wp(ix1,ix2,ix3-1,v2_))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))*invr(ix^D)
        ! idim=3, idir=3
        lambda(ix^D,3,3)=(wp(ix1,ix2,ix3+1,v3_)-wp(ix1,ix2,ix3-1,v3_))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))*invr(ix^D)*two*qdtmu
        ! dv_i/d_j + dv_j/d_i
        lambda(ix^D,1,2)=(lambda(ix^D,1,2)+lambda(ix^D,2,1))*qdtmu
        lambda(ix^D,2,1)=lambda(ix^D,1,2)
        lambda(ix^D,1,3)=(lambda(ix^D,1,3)+lambda(ix^D,3,1))*qdtmu
        lambda(ix^D,3,1)=lambda(ix^D,1,3)
        lambda(ix^D,2,3)=(lambda(ix^D,2,3)+lambda(ix^D,3,2))*qdtmu
        lambda(ix^D,3,2)=lambda(ix^D,2,3)
        divv23=third*(lambda(ix^D,1,1)+lambda(ix^D,2,2)+lambda(ix^D,3,3))
        lambda(ix^D,1,1)=lambda(ix^D,1,1)-divv23
        lambda(ix^D,2,2)=lambda(ix^D,2,2)-divv23
        lambda(ix^D,3,3)=lambda(ix^D,3,3)-divv23
        if(energy) then
          vlambda(ix^D,1)=wp(ix^D,v1_)*lambda(ix^D,1,1)+wp(ix^D,v2_)*lambda(ix^D,2,1)+wp(ix^D,v3_)*lambda(ix^D,3,1)
          vlambda(ix^D,2)=wp(ix^D,v1_)*lambda(ix^D,1,2)+wp(ix^D,v2_)*lambda(ix^D,2,2)+wp(ix^D,v3_)*lambda(ix^D,3,2)
          vlambda(ix^D,3)=wp(ix^D,v1_)*lambda(ix^D,1,3)+wp(ix^D,v2_)*lambda(ix^D,2,3)+wp(ix^D,v3_)*lambda(ix^D,3,3)
        end if
      {end do\}
      ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr)
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        w(ix^D,mom(1))=(lambda(ix1+1,ix2,ix3,1,1)-lambda(ix1-1,ix2,ix3,1,1))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))+&
                       (lambda(ix1,ix2+1,ix3,2,1)-lambda(ix1,ix2-1,ix3,2,1))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))*invr(ix^D)+&
                       (lambda(ix1,ix2,ix3+1,3,1)-lambda(ix1,ix2,ix3-1,3,1))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))+&
                       (lambda(ix^D,1,1)-lambda(ix^D,2,2))*invr(ix^D)+w(ix^D,mom(1))
        w(ix^D,mom(2))=(lambda(ix1+1,ix2,ix3,1,2)-lambda(ix1-1,ix2,ix3,1,2))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))+&
                       (lambda(ix1,ix2+1,ix3,2,2)-lambda(ix1,ix2-1,ix3,2,2))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))*invr(ix^D)+&
                       (lambda(ix1,ix2,ix3+1,3,2)-lambda(ix1,ix2,ix3-1,3,2))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))+&
                       two*lambda(ix^D,1,2)*invr(ix^D)+w(ix^D,mom(2))
        w(ix^D,mom(3))=(lambda(ix1+1,ix2,ix3,1,3)-lambda(ix1-1,ix2,ix3,1,3))/(x(ix1+1,ix2,ix3,1)-x(ix1-1,ix2,ix3,1))+&
                       (lambda(ix1,ix2+1,ix3,2,3)-lambda(ix1,ix2-1,ix3,2,3))/(x(ix1,ix2+1,ix3,2)-x(ix1,ix2-1,ix3,2))*invr(ix^D)+&
                       (lambda(ix1,ix2,ix3+1,3,3)-lambda(ix1,ix2,ix3-1,3,3))/(x(ix1,ix2,ix3+1,3)-x(ix1,ix2,ix3-1,3))+&
                        lambda(ix^D,1,3)*invr(ix^D)+w(ix^D,mom(3))
      {end do\}
      }
      {^IFTWOD
      {do ix^DB=ixmin^DB,ixmax^DB\}
        invr(ix^D)=1.d0/x(ix^D,1)
        ! idim=1, idir=1
        lambda(ix^D,1,1)=(wp(ix1+1,ix2,v1_)-wp(ix1-1,ix2,v1_))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))*two*qdtmu
        ! idim=1, idir=2
        lambda(ix^D,1,2)=(wp(ix1+1,ix2,v2_)-wp(ix1-1,ix2,v2_))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))
        if(phi_==2) then
          ! idim=2, idir=1
          lambda(ix^D,2,1)=((wp(ix1,ix2+1,v1_)-wp(ix1,ix2-1,v1_))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))-wp(ix^D,v2_))*invr(ix^D)
          ! idim=2, idir=2
          lambda(ix^D,2,2)=((wp(ix1,ix2+1,v2_)-wp(ix1,ix2-1,v2_))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))+wp(ix^D,v1_))*invr(ix^D)*&
                           two*qdtmu
          if(ndir==3) then
            ! idim=2, idir=3
            lambda(ix^D,2,3)=(wp(ix1,ix2+1,v3_)-wp(ix1,ix2-1,v3_))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))*invr(ix^D)
            ! idim=3, idir=1
            lambda(ix^D,3,1)=0.d0
            ! idim=3, idir=2
            lambda(ix^D,3,2)=0.d0
            ! idim=3, idir=3
            lambda(ix^D,3,3)=0.d0
          end if
        else
          ! idim=2, idir=1
          lambda(ix^D,2,1)=(wp(ix1,ix2+1,v1_)-wp(ix1,ix2-1,v1_))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))*invr(ix^D)
          ! idim=2, idir=2
          lambda(ix^D,2,2)=(wp(ix1,ix2+1,v2_)-wp(ix1,ix2-1,v2_))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))*invr(ix^D)*two*qdtmu
          if(ndir==3) then
            ! idim=2, idir=3
            lambda(ix^D,2,3)=(wp(ix1,ix2+1,v3_)-wp(ix1,ix2-1,v3_))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))*invr(ix^D)
            ! idim=3, idir=1
            lambda(ix^D,3,1)=-wp(ix^D,v3_)*invr(ix^D)
            ! idim=3, idir=2
            lambda(ix^D,3,2)=0.d0
            ! idim=3, idir=3
            lambda(ix^D,3,3)=wp(ix^D,v1_)*invr(ix^D)
          end if
        end if
        ! dv_i/d_j + dv_j/d_i
        lambda(ix^D,1,2)=(lambda(ix^D,1,2)+lambda(ix^D,2,1))*qdtmu
        lambda(ix^D,2,1)=lambda(ix^D,1,2)
        divv23=third*(lambda(ix^D,1,1)+lambda(ix^D,2,2)+lambda(ix^D,3,3))
        lambda(ix^D,1,1)=lambda(ix^D,1,1)-divv23
        lambda(ix^D,2,2)=lambda(ix^D,2,2)-divv23
        if(ndir==3) then
          ! idim=1, idir=3
          lambda(ix^D,1,3)=(wp(ix1+1,ix2,v3_)-wp(ix1-1,ix2,v3_))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))
          lambda(ix^D,1,3)=(lambda(ix^D,1,3)+lambda(ix^D,3,1))*qdtmu
          lambda(ix^D,3,1)=lambda(ix^D,1,3)
          lambda(ix^D,2,3)=(lambda(ix^D,2,3)+lambda(ix^D,3,2))*qdtmu
          lambda(ix^D,3,2)=lambda(ix^D,2,3)
          lambda(ix^D,3,3)=lambda(ix^D,3,3)-divv23
          if(energy) then
            vlambda(ix^D,1)=wp(ix^D,v1_)*lambda(ix^D,1,1)+wp(ix^D,v2_)*lambda(ix^D,2,1)+wp(ix^D,v3_)*lambda(ix^D,3,1)
            vlambda(ix^D,2)=wp(ix^D,v1_)*lambda(ix^D,1,2)+wp(ix^D,v2_)*lambda(ix^D,2,2)+wp(ix^D,v3_)*lambda(ix^D,3,2)
            vlambda(ix^D,3)=wp(ix^D,v1_)*lambda(ix^D,1,3)+wp(ix^D,v2_)*lambda(ix^D,2,3)+wp(ix^D,v3_)*lambda(ix^D,3,3)
          end if
        else if(energy) then
          vlambda(ix^D,1)=wp(ix^D,v1_)*lambda(ix^D,1,1)+wp(ix^D,v2_)*lambda(ix^D,2,1)
          vlambda(ix^D,2)=wp(ix^D,v1_)*lambda(ix^D,1,2)+wp(ix^D,v2_)*lambda(ix^D,2,2)
        end if
      {end do\}
      ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr)
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        if(ndir==2) then
          if(phi_==2) then
            w(ix^D,mom(1))=(lambda(ix1+1,ix2,1,1)-lambda(ix1-1,ix2,1,1))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                           (lambda(ix1,ix2+1,2,1)-lambda(ix1,ix2-1,2,1))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))*invr(ix^D)+&
                           (lambda(ix^D,1,1)-lambda(ix^D,2,2))*invr(ix^D)+w(ix^D,mom(1))
            w(ix^D,mom(2))=(lambda(ix1+1,ix2,1,2)-lambda(ix1-1,ix2,1,2))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                           (lambda(ix1,ix2+1,2,2)-lambda(ix1,ix2-1,2,2))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))*invr(ix^D)+&
                        two*lambda(ix^D,1,2)*invr(ix^D)+w(ix^D,mom(2))
          else
            w(ix^D,mom(1))=(lambda(ix1+1,ix2,1,1)-lambda(ix1-1,ix2,1,1))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                           (lambda(ix1,ix2+1,2,1)-lambda(ix1,ix2-1,2,1))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))+&
                            lambda(ix^D,1,1)*invr(ix^D)+w(ix^D,mom(1))
            w(ix^D,mom(2))=(lambda(ix1+1,ix2,1,2)-lambda(ix1-1,ix2,1,2))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                           (lambda(ix1,ix2+1,2,2)-lambda(ix1,ix2-1,2,2))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))+&
                            lambda(ix^D,1,2)*invr(ix^D)+w(ix^D,mom(2))
          end if
        else
          if(phi_==2) then
            w(ix^D,mom(1))=(lambda(ix1+1,ix2,1,1)-lambda(ix1-1,ix2,1,1))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                           (lambda(ix1,ix2+1,2,1)-lambda(ix1,ix2-1,2,1))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))*invr(ix^D)+&
                           (lambda(ix^D,1,1)-lambda(ix^D,2,2))*invr(ix^D)+w(ix^D,mom(1))
            w(ix^D,mom(2))=(lambda(ix1+1,ix2,1,2)-lambda(ix1-1,ix2,1,2))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                           (lambda(ix1,ix2+1,2,2)-lambda(ix1,ix2-1,2,2))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))*invr(ix^D)+&
                           two*lambda(ix^D,1,2)*invr(ix^D)+w(ix^D,mom(2))
            w(ix^D,mom(3))=(lambda(ix1+1,ix2,1,3)-lambda(ix1-1,ix2,1,3))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                           (lambda(ix1,ix2+1,2,3)-lambda(ix1,ix2-1,2,3))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))*invr(ix^D)+&
                            lambda(ix^D,1,3)*invr(ix^D)+w(ix^D,mom(3))
          else
            w(ix^D,mom(1))=(lambda(ix1+1,ix2,1,1)-lambda(ix1-1,ix2,1,1))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                           (lambda(ix1,ix2+1,2,1)-lambda(ix1,ix2-1,2,1))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))+&
                           (lambda(ix^D,1,1)-lambda(ix^D,2,2))*invr(ix^D)+w(ix^D,mom(1))
            w(ix^D,mom(2))=(lambda(ix1+1,ix2,1,3)-lambda(ix1-1,ix2,1,3))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                           (lambda(ix1,ix2+1,2,3)-lambda(ix1,ix2-1,2,3))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))+&
                           two*lambda(ix^D,1,3)*invr(ix^D)+w(ix^D,mom(2))
            w(ix^D,mom(3))=(lambda(ix1+1,ix2,1,2)-lambda(ix1-1,ix2,1,2))/(x(ix1+1,ix2,1)-x(ix1-1,ix2,1))+&
                           (lambda(ix1,ix2+1,2,2)-lambda(ix1,ix2-1,2,2))/(x(ix1,ix2+1,2)-x(ix1,ix2-1,2))+&
                            lambda(ix^D,1,2)*invr(ix^D)+w(ix^D,mom(3))
          end if
        end if
      {end do\}
      }
      {^IFONED
      {do ix^DB=ixmin^DB,ixmax^DB\}
        invr(ix^D)=1.d0/x(ix^D,1)
        if(ndir==1) then
          ! idim=1, idir=1
          lambda(ix^D,1,1)=(wp(ix1+1,v1_)-wp(ix1-1,v1_))/(x(ix1+1,1)-x(ix1-1,1))*two*qdtmu
          divv23=third*lambda(ix^D,1,1)
          lambda(ix^D,1,1)=lambda(ix^D,1,1)-divv23
          if(energy) then
            vlambda(ix^D,1)=wp(ix^D,v1_)*lambda(ix^D,1,1)
          end if
        else if(ndir==2) then
          ! idim=1, idir=1
          lambda(ix^D,1,1)=(wp(ix1+1,v1_)-wp(ix1-1,v1_))/(x(ix1+1,1)-x(ix1-1,1))*two*qdtmu
          ! idim=1, idir=2
          lambda(ix^D,1,2)=(wp(ix1+1,v2_)-wp(ix1-1,v2_))/(x(ix1+1,1)-x(ix1-1,1))
          if(phi_==2) then
            ! idim=2, idir=1
            lambda(ix^D,2,1)=-wp(ix^D,v2_)*invr(ix^D)
            ! idim=2, idir=2
            lambda(ix^D,2,2)=+wp(ix^D,v1_)*invr(ix^D)*two*qdtmu
          else
            ! idim=2, idir=1
            lambda(ix^D,2,1)=0.d0
            ! idim=2, idir=2
            lambda(ix^D,2,2)=0.d0
          end if
          ! dv_i/d_j + dv_j/d_i
          lambda(ix^D,1,2)=(lambda(ix^D,1,2)+lambda(ix^D,2,1))*qdtmu
          lambda(ix^D,2,1)=lambda(ix^D,1,2)
          divv23=third*(lambda(ix^D,1,1)+lambda(ix^D,2,2))
          lambda(ix^D,1,1)=lambda(ix^D,1,1)-divv23
          lambda(ix^D,2,2)=lambda(ix^D,2,2)-divv23
          if(energy) then
            vlambda(ix^D,1)=wp(ix^D,v1_)*lambda(ix^D,1,1)+wp(ix^D,v2_)*lambda(ix^D,2,1)
            vlambda(ix^D,2)=wp(ix^D,v1_)*lambda(ix^D,1,2)+wp(ix^D,v2_)*lambda(ix^D,2,2)
          end if
        else
          ! idim=1, idir=1
          lambda(ix^D,1,1)=(wp(ix1+1,v1_)-wp(ix1-1,v1_))/(x(ix1+1,1)-x(ix1-1,1))*two*qdtmu
          ! idim=1, idir=2
          lambda(ix^D,1,2)=(wp(ix1+1,v2_)-wp(ix1-1,v2_))/(x(ix1+1,1)-x(ix1-1,1))
          ! idim=1, idir=3
          lambda(ix^D,1,3)=(wp(ix1+1,v3_)-wp(ix1-1,v3_))/(x(ix1+1,1)-x(ix1-1,1))
          if(phi_==2) then
            ! idim=2, idir=1
            lambda(ix^D,2,1)=-wp(ix^D,v2_)*invr(ix^D)
            ! idim=2, idir=2
            lambda(ix^D,2,2)=+wp(ix^D,v1_)*invr(ix^D)*two*qdtmu
            ! idim=2, idir=3
            lambda(ix^D,2,3)=0.d0
            ! idim=3, idir=1
            lambda(ix^D,3,1)=0.d0
            ! idim=3, idir=2
            lambda(ix^D,3,2)=0.d0
            ! idim=3, idir=3
            lambda(ix^D,3,3)=0.d0
          else
            ! idim=2, idir=1
            lambda(ix^D,2,1)=0.d0
            ! idim=2, idir=2
            lambda(ix^D,2,2)=0.d0
            ! idim=2, idir=3
            lambda(ix^D,2,3)=0.d0
            ! idim=3, idir=1
            lambda(ix^D,3,1)=-wp(ix^D,v3_)*invr(ix^D)
            ! idim=3, idir=2
            lambda(ix^D,3,2)=0.d0
            ! idim=3, idir=3
            lambda(ix^D,3,3)=wp(ix^D,v1_)*invr(ix^D)*two*qdtmu
          end if
          ! dv_i/d_j + dv_j/d_i
          lambda(ix^D,1,2)=(lambda(ix^D,1,2)+lambda(ix^D,2,1))*qdtmu
          lambda(ix^D,2,1)=lambda(ix^D,1,2)
          lambda(ix^D,1,3)=(lambda(ix^D,1,3)+lambda(ix^D,3,1))*qdtmu
          lambda(ix^D,3,1)=lambda(ix^D,1,3)
          lambda(ix^D,2,3)=(lambda(ix^D,2,3)+lambda(ix^D,3,2))*qdtmu
          lambda(ix^D,3,2)=lambda(ix^D,2,3)
          divv23=third*(lambda(ix^D,1,1)+lambda(ix^D,2,2)+lambda(ix^D,3,3))
          lambda(ix^D,1,1)=lambda(ix^D,1,1)-divv23
          lambda(ix^D,2,2)=lambda(ix^D,2,2)-divv23
          lambda(ix^D,3,3)=lambda(ix^D,3,3)-divv23
          if(energy) then
            vlambda(ix^D,1)=wp(ix^D,v1_)*lambda(ix^D,1,1)+wp(ix^D,v2_)*lambda(ix^D,2,1)+wp(ix^D,v3_)*lambda(ix^D,3,1)
            vlambda(ix^D,2)=wp(ix^D,v1_)*lambda(ix^D,1,2)+wp(ix^D,v2_)*lambda(ix^D,2,2)+wp(ix^D,v3_)*lambda(ix^D,3,2)
            vlambda(ix^D,3)=wp(ix^D,v1_)*lambda(ix^D,1,3)+wp(ix^D,v2_)*lambda(ix^D,2,3)+wp(ix^D,v3_)*lambda(ix^D,3,3)
          end if
        end if
      {end do\}
      ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr)
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        if(ndir==1) then
          w(ix^D,mom(1))=(lambda(ix1+1,1,1)-lambda(ix1-1,1,1))/(x(ix1+1,1)-x(ix1-1,1))+w(ix^D,mom(1))
        else if(ndir==2) then
          if(phi_==2) then
            w(ix^D,mom(1))=(lambda(ix1+1,1,1)-lambda(ix1-1,1,1))/(x(ix1+1,1)-x(ix1-1,1))+&
                           (lambda(ix^D,1,1)-lambda(ix^D,2,2))*invr(ix^D)+w(ix^D,mom(1))
            w(ix^D,mom(2))=(lambda(ix1+1,1,2)-lambda(ix1-1,1,2))/(x(ix1+1,1)-x(ix1-1,1))+&
                           two*lambda(ix^D,1,2)*invr(ix^D)+w(ix^D,mom(2))
          else
            w(ix^D,mom(1))=(lambda(ix1+1,1,1)-lambda(ix1-1,1,1))/(x(ix1+1,1)-x(ix1-1,1))+&
                           (lambda(ix^D,1,1)-lambda(ix^D,2,2))*invr(ix^D)+w(ix^D,mom(1))
          end if
        else
          if(phi_==2) then
            w(ix^D,mom(1))=(lambda(ix1+1,1,1)-lambda(ix1-1,1,1))/(x(ix1+1,1)-x(ix1-1,1))+&
                           (lambda(ix^D,1,1)-lambda(ix^D,2,2))*invr(ix^D)+w(ix^D,mom(1))
            w(ix^D,mom(2))=(lambda(ix1+1,1,2)-lambda(ix1-1,1,2))/(x(ix1+1,1)-x(ix1-1,1))+&
                           two*lambda(ix^D,1,2)*invr(ix^D)+w(ix^D,mom(2))
            w(ix^D,mom(3))=(lambda(ix1+1,1,3)-lambda(ix1-1,1,3))/(x(ix1+1,1)-x(ix1-1,1))+&
                            lambda(ix^D,1,3)*invr(ix^D)+w(ix^D,mom(3))
          else
            w(ix^D,mom(1))=(lambda(ix1+1,1,1)-lambda(ix1-1,1,1))/(x(ix1+1,1)-x(ix1-1,1))+&
                           (lambda(ix^D,1,1)-lambda(ix^D,2,2))*invr(ix^D)+w(ix^D,mom(1))
            w(ix^D,mom(2))=(lambda(ix1+1,1,3)-lambda(ix1-1,1,3))/(x(ix1+1,1)-x(ix1-1,1))+&
                           two*lambda(ix^D,1,3)*invr(ix^D)+w(ix^D,mom(2))
            w(ix^D,mom(3))=(lambda(ix1+1,1,2)-lambda(ix1-1,1,2))/(x(ix1+1,1)-x(ix1-1,1))+&
                            lambda(ix^D,1,2)*invr(ix^D)+w(ix^D,mom(3))
          end if
        end if
      {end do\}
      }
      if(energy) then
        ! de/dt= +div(v.dot.[mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v *kr])
        ! thus e=e+d_i v_j tensor_ji
        call divvector(vlambda,ixI^L,ixO^L,invr)
        w(ixO^S,e_)=w(ixO^S,e_)+invr(ixO^S)
      end if

    end if

  end subroutine viscosity_add_source_cylinder

  subroutine viscosity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    ! Check diffusion time limit for dt < dtdiffpar * dx**2 / (mu/rho)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: dx^D, x(ixI^S,1:ndim)
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision, intent(inout) :: dtnew

    double precision :: tmp(ixI^S)
    double precision:: dtdiff_visc, dxinv2(1:ndim), max_mu
    integer:: idim

    ! Calculate the kinematic viscosity tmp=mu/rho
    ! here vc_mu must be non-zero!!!

    tmp(ixO^S)=vc_mu/w(ixO^S,rho_)

    if(slab_uniform)then
      ^D&dxinv2(^D)=one/dx^D**2;
      do idim=1,ndim
         dtdiff_visc=dtdiffpar/maxval(tmp(ixO^S)*dxinv2(idim))
         ! limit the time step
         dtnew=min(dtnew,dtdiff_visc)
      enddo
    else
      do idim=1,ndim
         max_mu=maxval(tmp(ixO^S)/block%ds(ixO^S,idim)**2)
         dtdiff_visc=dtdiffpar/max_mu
         ! limit the time step
         dtnew=min(dtnew,dtdiff_visc)
      enddo
    endif

  end subroutine viscosity_get_dt

  ! viscInDiv
  ! Get the viscous stress tensor terms in the idim direction
  ! Beware : a priori, won't work for ndir /= ndim
  ! Rq : we work with primitive w variables here
  ! Rq : ixO^L is already extended by 1 unit in the direction we work on

  subroutine visc_get_flux_prim(w, x, ixI^L, ixO^L, idim, f, energy)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: w(ixI^S, 1:nw), x(ixI^S, 1:^ND)
    double precision, intent(inout) :: f(ixI^S, nwflux)
    logical, intent(in) :: energy
    integer                         :: idir, i
    double precision :: v(ixI^S,1:ndir)

    double precision                :: divergence(ixI^S)

    double precision:: lambda(ixI^S,ndir) !, tmp(ixI^S) !gradV(ixI^S,ndir,ndir)

    if (.not. viscInDiv) return

    do i=1,ndir
     v(ixI^S,i)=w(ixI^S,i+1)
    enddo
    call divvector(v,ixI^L,ixO^L,divergence)

    call get_crossgrad(ixI^L,ixO^L,x,w,idim,lambda)
    lambda(ixO^S,idim) = lambda(ixO^S,idim) - (2.d0/3.d0) * divergence(ixO^S)

    ! Compute the idim-th row of the viscous stress tensor
    do idir = 1, ndir
      f(ixO^S, mom(idir)) = f(ixO^S, mom(idir)) - vc_mu * lambda(ixO^S,idir)
      if (energy) f(ixO^S, e_) = f(ixO^S, e_) - vc_mu * lambda(ixO^S,idir) * v(ixI^S,idir)
    enddo

  end subroutine visc_get_flux_prim

  ! Compute the cross term ( d_i v_j + d_j v_i in Cartesian BUT NOT IN
  ! CYLINDRICAL AND SPHERICAL )
  subroutine get_crossgrad(ixI^L,ixO^L,x,w,idim,cross)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: w(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    double precision, intent(out)   :: cross(ixI^S,ndir)

    double precision :: tmp(ixI^S), v(ixI^S)
    integer :: idir

    if (ndir/=ndim) call mpistop("This formula are probably wrong for ndim/=ndir")
    ! Beware also, we work w/ the angle as the 3rd component in cylindrical
    ! and the colatitude as the 2nd one in spherical
    cross(ixI^S,:)=zero
    tmp(ixI^S)=zero
    select case(coordinate)
    case (Cartesian,Cartesian_stretched)
      call cart_cross_grad(ixI^L,ixO^L,x,w,idim,cross)
    case (cylindrical)
      if (idim==1) then
        ! for rr and rz
        call cart_cross_grad(ixI^L,ixO^L,x,w,idim,cross)
        ! then we overwrite rth w/ the correct expression
        {^NOONED
        v(ixI^S)=w(ixI^S,mom(1)) ! v_r
        call gradient(v,ixI^L,ixO^L,2,tmp) ! d_th (rq : already contains 1/r)
        cross(ixI^S,2)=tmp(ixI^S)
        v(ixI^S)=w(ixI^S,mom(2))/x(ixI^S,1)  ! v_th / r
        call gradient(v,ixI^L,ixO^L,1,tmp) ! d_r
        cross(ixI^S,2)=cross(ixI^S,2)+tmp(ixI^S)*x(ixI^S,1)
        }
      elseif (idim==2) then
        ! thr (idem as above)
        v(ixI^S)=w(ixI^S,mom(1)) ! v_r
        {^NOONED
        call gradient(v,ixI^L,ixO^L,2,tmp) ! d_th
        cross(ixI^S,1)=tmp(ixI^S)
        v(ixI^S)=w(ixI^S,mom(2))/x(ixI^S,1)  ! v_th / r
        call gradient(v,ixI^L,ixO^L,1,tmp) ! d_r
        cross(ixI^S,1)=cross(ixI^S,1)+tmp(ixI^S)*x(ixI^S,1)
        ! thth
        v(ixI^S)=w(ixI^S,mom(2)) ! v_th
        call gradient(v,ixI^L,ixO^L,2,tmp) ! d_th
        cross(ixI^S,2)=two*tmp(ixI^S)
        v(ixI^S)=w(ixI^S,mom(1)) ! v_r
        cross(ixI^S,2)=cross(ixI^S,2)+two*v(ixI^S)/x(ixI^S,1) ! + 2 vr/r
        !thz
        v(ixI^S)=w(ixI^S,mom(3)) ! v_z
        call gradient(v,ixI^L,ixO^L,2,tmp) ! d_th
        }
        cross(ixI^S,3)=tmp(ixI^S)
        v(ixI^S)=w(ixI^S,mom(2))  ! v_th
        {^IFTHREED
        call gradient(v,ixI^L,ixO^L,3,tmp) ! d_z
        }
        cross(ixI^S,3)=cross(ixI^S,3)+tmp(ixI^S)
        {^IFTHREED
      elseif (idim==3) then
        ! for zz and rz
        call cart_cross_grad(ixI^L,ixO^L,x,w,idim,cross)
        ! then we overwrite zth w/ the correct expression
        !thz
        v(ixI^S)=w(ixI^S,mom(3)) ! v_z
        call gradient(v,ixI^L,ixO^L,2,tmp) ! d_th
        cross(ixI^S,2)=tmp(ixI^S)
        v(ixI^S)=w(ixI^S,mom(2))  ! v_th
        call gradient(v,ixI^L,ixO^L,3,tmp) ! d_z
        cross(ixI^S,2)=cross(ixI^S,2)+tmp(ixI^S)
        }
      endif
    case (spherical)
      if (idim==1) then
        ! rr (normal, simply 2 * dr vr)
        v(ixI^S)=w(ixI^S,mom(1)) ! v_r
        call gradient(v,ixI^L,ixO^L,1,tmp) ! d_r
        cross(ixI^S,1)=two*tmp(ixI^S)
        !rth
        v(ixI^S)=w(ixI^S,mom(1)) ! v_r
        {^NOONED
        call gradient(v,ixI^L,ixO^L,2,tmp) ! d_th (rq : already contains 1/r)
        cross(ixI^S,2)=tmp(ixI^S)
        v(ixI^S)=w(ixI^S,mom(2))/x(ixI^S,1)  ! v_th / r
        call gradient(v,ixI^L,ixO^L,1,tmp) ! d_r
        cross(ixI^S,2)=cross(ixI^S,2)+tmp(ixI^S)*x(ixI^S,1)
        }
        {^IFTHREED
        !rph
        v(ixI^S)=w(ixI^S,mom(1)) ! v_r
        call gradient(v,ixI^L,ixO^L,3,tmp) ! d_phi (rq : contains 1/rsin(th))
        cross(ixI^S,3)=tmp(ixI^S)
        v(ixI^S)=w(ixI^S,mom(3))/x(ixI^S,1) ! v_phi / r
        call gradient(v,ixI^L,ixO^L,1,tmp) ! d_r
        cross(ixI^S,3)=cross(ixI^S,3)+tmp(ixI^S)*x(ixI^S,1)
        }
      elseif (idim==2) then
        ! thr
        v(ixI^S)=w(ixI^S,mom(1)) ! v_r
        {^NOONED
        call gradient(v,ixI^L,ixO^L,2,tmp) ! d_th (rq : already contains 1/r)
        cross(ixI^S,1)=tmp(ixI^S)
        v(ixI^S)=w(ixI^S,mom(2))/x(ixI^S,1)  ! v_th / r
        call gradient(v,ixI^L,ixO^L,1,tmp) ! d_r
        cross(ixI^S,1)=cross(ixI^S,1)+tmp(ixI^S)*x(ixI^S,1)
        ! thth
        v(ixI^S)=w(ixI^S,mom(2)) ! v_th
        call gradient(v,ixI^L,ixO^L,2,tmp) ! d_th
        cross(ixI^S,2)=two*tmp(ixI^S)
        v(ixI^S)=w(ixI^S,mom(1)) ! v_r
        cross(ixI^S,2)=cross(ixI^S,2)+two*v(ixI^S)/x(ixI^S,1) ! + 2 vr/r
        }
        {^IFTHREED
        !thph
        v(ixI^S)=w(ixI^S,mom(2)) ! v_th
        call gradient(v,ixI^L,ixO^L,3,tmp) ! d_phi (rq : contains 1/rsin(th))
        cross(ixI^S,3)=tmp(ixI^S)
        v(ixI^S)=w(ixI^S,mom(3))/dsin(x(ixI^S,2)) ! v_ph / sin(th)
        call gradient(v,ixI^L,ixO^L,2,tmp) ! d_th
        cross(ixI^S,3)=cross(ixI^S,3)+tmp(ixI^S)*dsin(x(ixI^S,2))
        }
        {^IFTHREED
      elseif (idim==3) then
        !phr
        v(ixI^S)=w(ixI^S,mom(1)) ! v_r
        call gradient(v,ixI^L,ixO^L,3,tmp) ! d_phi
        cross(ixI^S,1)=tmp(ixI^S)
        v(ixI^S)=w(ixI^S,mom(3))/x(ixI^S,1) ! v_phi / r
        call gradient(v,ixI^L,ixO^L,1,tmp) ! d_r
        cross(ixI^S,1)=cross(ixI^S,1)+tmp(ixI^S)*x(ixI^S,1)
        !phth
        v(ixI^S)=w(ixI^S,mom(2)) ! v_th
        call gradient(v,ixI^L,ixO^L,3,tmp) ! d_phi
        cross(ixI^S,2)=tmp(ixI^S)
        v(ixI^S)=w(ixI^S,mom(3))/dsin(x(ixI^S,2)) ! v_ph / sin(th)
        call gradient(v,ixI^L,ixO^L,2,tmp) ! d_th
        cross(ixI^S,2)=cross(ixI^S,2)+tmp(ixI^S)*dsin(x(ixI^S,2))
        !phph
        v(ixI^S)=w(ixI^S,mom(3)) ! v_ph
        call gradient(v,ixI^L,ixO^L,3,tmp) ! d_phi
        cross(ixI^S,3)=two*tmp(ixI^S)
        v(ixI^S)=w(ixI^S,mom(1)) ! v_r
        cross(ixI^S,3)=cross(ixI^S,3)+two*v(ixI^S)/x(ixI^S,1) ! + 2 vr/r
        v(ixI^S)=w(ixI^S,mom(2)) ! v_th
        cross(ixI^S,3)=cross(ixI^S,3)+two*v(ixI^S)/(x(ixI^S,1)*dtan(x(ixI^S,2))) ! + 2 vth/(rtan(th))
        }
      endif
    case default
      call mpistop("Unknown geometry specified")
    end select

  end subroutine get_crossgrad

  !> yields d_i v_j + d_j v_i for a given i, OK in Cartesian and for some
  !> tensor terms in cylindrical (rr & rz) and in spherical (rr)
  subroutine cart_cross_grad(ixI^L,ixO^L,x,w,idim,cross)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: w(ixI^S, 1:nw), x(ixI^S, 1:^ND)
    double precision, intent(out)   :: cross(ixI^S,ndir)

    double precision :: tmp(ixI^S), v(ixI^S)
    integer :: idir

    v(ixI^S)=w(ixI^S,mom(idim))
    do idir=1,ndir
      call gradient(v,ixI^L,ixO^L,idir,tmp)
      cross(ixO^S,idir)=tmp(ixO^S)
    enddo
    do idir=1,ndir
      v(ixI^S)=w(ixI^S,mom(idir))
      call gradient(v,ixI^L,ixO^L,idim,tmp)
      cross(ixO^S,idir)=cross(ixO^S,idir)+tmp(ixO^S)
    enddo

  end subroutine cart_cross_grad

  subroutine visc_add_source_geom(qdt, ixI^L, ixO^L, wCT, w, x)
    use mod_global_parameters
    use mod_geometry
    ! w and wCT conservative variables here
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S, 1:ndim)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)
    ! to change and to set as a parameter in the parfile once the possibility to
    ! solve the equations in an angular momentum conserving form has been
    ! implemented (change tvdlf.t eg)
    double precision :: vv(ixI^S), divergence(ixI^S)
    double precision :: tmp(ixI^S),tmp1(ixI^S)
    integer          :: i

    if (.not. viscInDiv) return

    select case (coordinate)
    case (cylindrical)
      ! thth tensor term - - -
        ! 1st the cross grad term
{^NOONED
      call gradient(wCT(ixI^S,mom(2)),ixI^L,ixO^L,2,tmp1) ! d_th
      tmp(ixO^S)=two*(tmp1(ixO^S)+wCT(ixO^S,mom(1))/x(ixO^S,1)) ! 2 ( d_th v_th / r + vr/r )
        ! 2nd the divergence
      call divvector(wCT(ixI^S,mom(1:ndir)),ixI^L,ixO^L,divergence)
      tmp(ixO^S) = tmp(ixO^S) - (2.d0/3.d0) * divergence(ixO^S)
      ! s[mr]=-thth/radius
      w(ixO^S,mom(1))=w(ixO^S,mom(1))-qdt*vc_mu*tmp(ixO^S)/x(ixO^S,1)
      ! rth tensor term - - -
      call gradient(wCT(ixI^S,mom(1)),ixI^L,ixO^L,2,tmp) ! d_th
      vv(ixI^S)=wCT(ixI^S,mom(2))/x(ixI^S,1)  ! v_th / r
      call gradient(vv,ixI^L,ixO^L,1,tmp1) ! d_r
      tmp(ixO^S)=tmp(ixO^S)+tmp1(ixO^S)*x(ixO^S,1)
      ! s[mphi]=+rth/radius
      w(ixO^S,mom(2))=w(ixO^S,mom(2))+qdt*vc_mu*tmp(ixO^S)/x(ixO^S,1)
}
    case (spherical)
      ! thth tensor term - - -
      ! 1st the cross grad term
{^NOONED
      call gradient(wCT(ixI^S,mom(2)),ixI^L,ixO^L,2,tmp1) ! d_th
      tmp(ixO^S)=two*(tmp1(ixO^S)+wCT(ixO^S,mom(1))/x(ixO^S,1)) ! 2 ( 1/r * d_th v_th + vr/r )
      ! 2nd the divergence
      call divvector(wCT(ixI^S,mom(1:ndir)),ixI^L,ixO^L,divergence)
      tmp(ixO^S) = tmp(ixO^S) - (2.d0/3.d0) * divergence(ixO^S)
      ! s[mr]=-thth/radius
      w(ixO^S,mom(1))=w(ixO^S,mom(1))-qdt*vc_mu*tmp(ixO^S)/x(ixO^S,1)
}
      ! phiphi tensor term - - -
      ! 1st the cross grad term
{^IFTHREED
      call gradient(wCT(ixI^S,mom(3)),ixI^L,ixO^L,3,tmp1) ! d_phi
      tmp(ixO^S)=two*(tmp1(ixO^S)+wCT(ixO^S,mom(1))/x(ixO^S,1)+wCT(ixO^S,mom(2))/(x(ixO^S,1)*dtan(x(ixO^S,2)))) ! 2 ( 1/rsinth * d_ph v_ph + vr/r + vth/rtanth )
}
      ! 2nd the divergence
      tmp(ixO^S) = tmp(ixO^S) - (2.d0/3.d0) * divergence(ixO^S)
      ! s[mr]=-phiphi/radius
      w(ixO^S,mom(1))=w(ixO^S,mom(1))-qdt*vc_mu*tmp(ixO^S)/x(ixO^S,1)
      ! s[mth]=-cotanth*phiphi/radius
{^NOONED
      w(ixO^S,mom(2))=w(ixO^S,mom(2))-qdt*vc_mu*tmp(ixO^S)/(x(ixO^S,1)*dtan(x(ixO^S,2)))
}
      ! rth tensor term - - -
      call gradient(wCT(ixI^S,mom(1)),ixI^L,ixO^L,2,tmp) ! d_th (rq : already contains 1/r)
      vv(ixI^S)=wCT(ixI^S,mom(2))/x(ixI^S,1)  ! v_th / r
      call gradient(vv,ixI^L,ixO^L,1,tmp1) ! d_r
      tmp(ixO^S)=tmp(ixO^S)+tmp1(ixO^S)*x(ixO^S,1)
      ! s[mth]=+rth/radius
      w(ixO^S,mom(2))=w(ixO^S,mom(2))+qdt*vc_mu*tmp(ixO^S)/x(ixO^S,1)
      ! rphi tensor term - - -
{^IFTHREED
      call gradient(wCT(ixI^S,mom(1)),ixI^L,ixO^L,3,tmp) ! d_phi (rq : contains 1/rsin(th))
}
      vv(ixI^S)=wCT(ixI^S,mom(3))/x(ixI^S,1) ! v_phi / r
      call gradient(vv,ixI^L,ixO^L,1,tmp1) ! d_r
      tmp(ixO^S)=tmp(ixO^S)+tmp1(ixO^S)*x(ixO^S,1)
      ! s[mphi]=+rphi/radius
      w(ixO^S,mom(3))=w(ixO^S,mom(3))+qdt*vc_mu*tmp(ixO^S)/x(ixO^S,1)
      ! phith tensor term - - -
{^IFTHREED
      call gradient(wCT(ixI^S,mom(2)),ixI^L,ixO^L,3,tmp) ! d_phi
}
{^NOONED
      vv(ixI^S)=wCT(ixI^S,mom(3))/dsin(x(ixI^S,2)) ! v_ph / sin(th)
      call gradient(vv,ixI^L,ixO^L,2,tmp1) ! d_th
      tmp(ixO^S)=tmp(ixO^S)+tmp1(ixO^S)*dsin(x(ixO^S,2))
      ! s[mphi]=+cotanth*phith/radius
      w(ixO^S,mom(3))=w(ixO^S,mom(3))+qdt*vc_mu*tmp(ixO^S)/(x(ixO^S,1)*dtan(x(ixO^S,2)))
}
    end select

  end subroutine visc_add_source_geom

  subroutine sub_add_source(qdt,ixI^L,ixO^L,wCT,wp,w,x,&
       energy,qsourcesplit,active)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw), wp(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in) :: energy,qsourcesplit
    logical, intent(inout) :: active
  end subroutine sub_add_source

end module mod_viscosity
