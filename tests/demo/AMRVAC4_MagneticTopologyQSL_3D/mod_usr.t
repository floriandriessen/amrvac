module mod_usr
  use mod_mhd
  use mod_eos, only: eos
  use mod_global_parameters, only: par_files,mype,refine_max_level
  use mod_magnetic_topology, only: mt_params_read,mt_run_topology_task
  use mod_usr_methods, only: usr_special_convert,usr_refine_grid

  implicit none

  integer, parameter :: np0=200*2
  integer, parameter :: nb0=3
  double precision, parameter :: rho0=1.d0, p0=1.d0
  double precision, parameter :: a0=3.d0*0.75d0
  double precision, parameter :: F0=60.04d0*0.75d0
  double precision, parameter :: L_cha0=4.0d0*0.75d0
  double precision, parameter :: d_cha0=6.0d0*0.75d0
  double precision, parameter :: qc=800.d0*0.75d0
  double precision, parameter :: q_cha0=qc/sqrt(4.d0*dpi)/nb0
  double precision, parameter :: x_cha0=0.d0
  double precision :: xs0(np0,3)

contains

  subroutine usr_init()
    call set_coordinate_system("Cartesian_3D")

    usr_init_one_grid => initonegrid_usr
    usr_special_convert => topology_bench_special_convert
    usr_refine_grid => tdm_bench_refine_grid

    call calc_geom(xs0,np0)
    call mt_params_read(par_files)
    call mhd_activate()
  end subroutine usr_init

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer :: ixK^L
    double precision :: B1(ixI^S,1:ndir),B2(ixI^S,1:ndir)

    B1=zero
    B2=zero
    ixK^L=ixO^L^LADD1;

    call bipoB(ixI^L,ixO^L,x,L_cha0,d_cha0,q_cha0,x_cha0,nb0,B1)
    call rbslB(ixI^L,ixO^L,ixK^L,x,a0,F0,xs0,np0,B2)

    w(ixO^S,rho_)=rho0
    w(ixO^S,p_)=p0
    w(ixO^S,mom(:))=zero
    w(ixO^S,mag(1))=B1(ixO^S,1)+B2(ixO^S,1)
    w(ixO^S,mag(2))=B1(ixO^S,2)+B2(ixO^S,2)
    w(ixO^S,mag(3))=B1(ixO^S,3)+B2(ixO^S,3)

    call eos%to_conserved(ixI^L,ixO^L,w,x)
  end subroutine initonegrid_usr

  subroutine topology_bench_special_convert(qunitconvert)
    integer, intent(in) :: qunitconvert
    double precision :: t0,t1

    call topology_bench_print_amr_summary()
    call cpu_time(t0)
    call mt_run_topology_task()
    call cpu_time(t1)
    if (mype==0) write(*,'(a,es14.6)') 'TOPOLOGY_BENCH_TIME seconds=',t1-t0
  end subroutine topology_bench_special_convert

  subroutine topology_bench_print_amr_summary()
    use mod_forest, only: nleafs,nleafs_active,nleafs_level
    integer :: level

    if (mype==0) then
      write(*,'(a,i0,a,i0,a,i0)') 'TOPOLOGY_BENCH_AMR nleafs=',nleafs, &
        ' active=',nleafs_active,' refine_max_level=',refine_max_level
      do level=1,refine_max_level
        write(*,'(a,i0,a,i0)') 'TOPOLOGY_BENCH_AMR_LEVEL level=',level, &
          ' nleafs=',nleafs_level(level)
      end do
    end if
  end subroutine topology_bench_print_amr_summary

  subroutine tdm_bench_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    ! Deterministic benchmark-only AMR criterion. It refines a broad low/mid
    ! corona region on level 1 and a tighter rope/arcade region on level 2.
    ! This is for AMR-input tracing smoke tests, not a physics accuracy model.
    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout) :: refine, coarsen
    logical :: broad_region, rope_region

    broad_region = any(abs(x(ixO^S,1)) <= 8.5d0 .and. &
      abs(x(ixO^S,2)) <= 4.5d0 .and. &
      x(ixO^S,3) >= 0.d0 .and. x(ixO^S,3) <= 7.5d0)
    rope_region = any(abs(x(ixO^S,1)) <= 7.5d0 .and. &
      abs(x(ixO^S,2)) <= 2.5d0 .and. &
      x(ixO^S,3) >= 0.d0 .and. x(ixO^S,3) <= 5.5d0) .or. &
      any(abs(x(ixO^S,1)) <= 9.d0 .and. &
      abs(x(ixO^S,2)) <= 4.d0 .and. &
      x(ixO^S,3) >= 0.d0 .and. x(ixO^S,3) <= 2.d0)

    if (level==1 .and. broad_region) then
      refine = 1
      coarsen = -1
    else if (level==2 .and. rope_region) then
      refine = 1
      coarsen = -1
    end if
  end subroutine tdm_bench_refine_grid

  subroutine calc_geom(xs,npo)
    integer, intent(in) :: npo
    double precision, dimension(npo,3), intent(inout) :: xs

    integer :: i,np
    double precision, dimension(npo/2) :: xs2,zs2

    ! Copied from the p3 caseA_bench RBSL reference. The p3 setup mirrors
    ! these 200 points to make the full 400-point rope axis, then scales the
    ! full path by 0.75.
    np = npo / 2
    xs2 = (/ &
-8.903300e+00, -8.844232e+00, -8.784798e+00, -8.724057e+00, &
-8.661875e+00, -8.598199e+00, -8.533963e+00, -8.470245e+00, &
-8.405256e+00, -8.339208e+00, -8.271733e+00, -8.203644e+00, &
-8.135463e+00, -8.066539e+00, -7.996664e+00, -7.925153e+00, &
-7.853625e+00, -7.781541e+00, -7.708385e+00, -7.634195e+00, &
-7.559840e+00, -7.484918e+00, -7.408818e+00, -7.331528e+00, &
-7.254341e+00, -7.176209e+00, -7.097165e+00, -7.017736e+00, &
-6.937332e+00, -6.857022e+00, -6.774957e+00, -6.692981e+00, &
-6.610043e+00, -6.526851e+00, -6.442644e+00, -6.357574e+00, &
-6.272701e+00, -6.186711e+00, -6.099932e+00, -6.012841e+00, &
-5.925377e+00, -5.837106e+00, -5.748077e+00, -5.658759e+00, &
-5.568575e+00, -5.478455e+00, -5.387353e+00, -5.295916e+00, &
-5.204021e+00, -5.111341e+00, -5.018374e+00, -4.924821e+00, &
-4.830822e+00, -4.736469e+00, -4.641461e+00, -4.546297e+00, &
-4.450258e+00, -4.354206e+00, -4.257292e+00, -4.160448e+00, &
-4.062914e+00, -3.965270e+00, -3.866860e+00, -3.768346e+00, &
-3.669268e+00, -3.570153e+00, -3.470299e+00, -3.370301e+00, &
-3.270081e+00, -3.169385e+00, -3.068491e+00, -2.967168e+00, &
-2.865657e+00, -2.763959e+00, -2.661765e+00, -2.559429e+00, &
-2.456890e+00, -2.354073e+00, -2.251005e+00, -2.147808e+00, &
-2.044322e+00, -1.940575e+00, -1.836688e+00, -1.732687e+00, &
-1.628506e+00, -1.524181e+00, -1.419687e+00, -1.315040e+00, &
-1.210248e+00, -1.105341e+00, -1.000327e+00, -8.952107e-01, &
-7.900097e-01, -6.847337e-01, -5.794107e-01, -4.740308e-01, &
-3.686029e-01, -2.631342e-01, -1.576413e-01, -5.213195e-02, &
 5.338301e-02,  1.588924e-01,  2.643853e-01,  3.698535e-01, &
 4.752814e-01,  5.806603e-01,  6.859825e-01,  7.912589e-01, &
 8.964605e-01,  1.001573e+00,  1.106587e+00,  1.211493e+00, &
 1.316280e+00,  1.420932e+00,  1.525436e+00,  1.629779e+00, &
 1.733946e+00,  1.837921e+00,  1.941829e+00,  2.045562e+00, &
 2.149077e+00,  2.252359e+00,  2.355388e+00,  2.458116e+00, &
 2.560669e+00,  2.663100e+00,  2.765204e+00,  2.866989e+00, &
 2.968414e+00,  3.069671e+00,  3.170810e+00,  3.271419e+00, &
 3.371580e+00,  3.471536e+00,  3.571340e+00,  3.670754e+00, &
 3.769488e+00,  3.868233e+00,  3.966450e+00,  4.064170e+00, &
 4.161618e+00,  4.258974e+00,  4.355311e+00,  4.451638e+00, &
 4.547336e+00,  4.642486e+00,  4.737663e+00,  4.832211e+00, &
 4.926107e+00,  5.019607e+00,  5.112339e+00,  5.204946e+00, &
 5.296849e+00,  5.388728e+00,  5.479428e+00,  5.570431e+00, &
 5.659633e+00,  5.749492e+00,  5.837932e+00,  5.926145e+00, &
 6.014028e+00,  6.101186e+00,  6.187353e+00,  6.273376e+00, &
 6.359097e+00,  6.443304e+00,  6.527468e+00,  6.611354e+00, &
 6.693816e+00,  6.775979e+00,  6.857531e+00,  6.938433e+00, &
 7.018276e+00,  7.098536e+00,  7.177525e+00,  7.255111e+00, &
 7.332228e+00,  7.409378e+00,  7.485415e+00,  7.560391e+00, &
 7.634839e+00,  7.709347e+00,  7.782858e+00,  7.854275e+00, &
 7.925539e+00,  7.996866e+00,  8.067614e+00,  8.136335e+00, &
 8.203881e+00,  8.272254e+00,  8.339618e+00,  8.406000e+00, &
 8.470889e+00,  8.534105e+00,  8.598362e+00,  8.661923e+00, &
 8.724171e+00,  8.785269e+00,  8.844481e+00,  8.903300e+00 /)
    zs2 = (/ &
 0.000000e+00,  8.741635e-02,  1.745744e-01,  2.608161e-01, &
 3.460240e-01,  4.301367e-01,  5.138323e-01,  5.979164e-01, &
 6.810060e-01,  7.632576e-01,  8.443512e-01,  9.249464e-01, &
 1.005449e+00,  1.085312e+00,  1.164344e+00,  1.241918e+00, &
 1.319479e+00,  1.396511e+00,  1.472520e+00,  1.547530e+00, &
 1.622384e+00,  1.696670e+00,  1.769742e+00,  1.841566e+00, &
 1.913499e+00,  1.984399e+00,  2.054277e+00,  2.123722e+00, &
 2.192038e+00,  2.260473e+00,  2.326797e+00,  2.393221e+00, &
 2.458438e+00,  2.523336e+00,  2.586916e+00,  2.649330e+00, &
 2.712017e+00,  2.773166e+00,  2.833184e+00,  2.892752e+00, &
 2.951766e+00,  3.009575e+00,  3.066208e+00,  3.122378e+00, &
 3.177153e+00,  3.232033e+00,  3.285263e+00,  3.337920e+00, &
 3.389772e+00,  3.440210e+00,  3.490116e+00,  3.538907e+00, &
 3.586841e+00,  3.634064e+00,  3.679967e+00,  3.725543e+00, &
 3.769247e+00,  3.812923e+00,  3.854644e+00,  3.896533e+00, &
 3.936781e+00,  3.976772e+00,  4.014830e+00,  4.052624e+00, &
 4.088913e+00,  4.125095e+00,  4.159191e+00,  4.192854e+00, &
 4.225861e+00,  4.257364e+00,  4.288246e+00,  4.317690e+00, &
 4.346458e+00,  4.374581e+00,  4.400834e+00,  4.426521e+00, &
 4.451408e+00,  4.475109e+00,  4.497679e+00,  4.519657e+00, &
 4.540249e+00,  4.559479e+00,  4.577929e+00,  4.595729e+00, &
 4.612434e+00,  4.628223e+00,  4.642863e+00,  4.656365e+00, &
 4.668704e+00,  4.680014e+00,  4.690279e+00,  4.699447e+00, &
 4.707580e+00,  4.714676e+00,  4.721032e+00,  4.726362e+00, &
 4.730644e+00,  4.733764e+00,  4.735914e+00,  4.736984e+00, &
 4.736979e+00,  4.735902e+00,  4.733751e+00,  4.730609e+00, &
 4.726322e+00,  4.720969e+00,  4.714599e+00,  4.707509e+00, &
 4.699382e+00,  4.690167e+00,  4.679900e+00,  4.668583e+00, &
 4.656211e+00,  4.642745e+00,  4.628176e+00,  4.612496e+00, &
 4.595690e+00,  4.577737e+00,  4.559394e+00,  4.540091e+00, &
 4.519651e+00,  4.498064e+00,  4.475292e+00,  4.451200e+00, &
 4.426383e+00,  4.401061e+00,  4.374459e+00,  4.346665e+00, &
 4.317582e+00,  4.287922e+00,  4.257849e+00,  4.226049e+00, &
 4.192877e+00,  4.159085e+00,  4.124860e+00,  4.089513e+00, &
 4.052296e+00,  4.015117e+00,  3.976574e+00,  3.936795e+00, &
 3.896345e+00,  3.855667e+00,  3.812627e+00,  3.769579e+00, &
 3.725157e+00,  3.679579e+00,  3.634056e+00,  3.587227e+00, &
 3.539115e+00,  3.490227e+00,  3.439916e+00,  3.389378e+00, &
 3.337571e+00,  3.285719e+00,  3.231817e+00,  3.178426e+00, &
 3.122072e+00,  3.066774e+00,  3.009261e+00,  2.951402e+00, &
 2.893035e+00,  2.833605e+00,  2.772734e+00,  2.711670e+00, &
 2.650184e+00,  2.586629e+00,  2.523025e+00,  2.459061e+00, &
 2.393267e+00,  2.327087e+00,  2.260182e+00,  2.192494e+00, &
 2.123538e+00,  2.055066e+00,  1.985152e+00,  1.913676e+00, &
 1.841682e+00,  1.769743e+00,  1.696639e+00,  1.622440e+00, &
 1.547693e+00,  1.473019e+00,  1.397373e+00,  1.319731e+00, &
 1.241940e+00,  1.164229e+00,  1.085990e+00,  1.005970e+00, &
 9.249375e-01,  8.445923e-01,  7.634271e-01,  6.814618e-01, &
 5.983057e-01,  5.138446e-01,  4.301784e-01,  3.460002e-01, &
 2.608550e-01,  1.748806e-01,  8.757697e-02,  0.000000e+00 /)
    xs(1:np,1) = xs2(np:1:-1)
    xs(1:np,3) = zs2(np:1:-1)
    do i = 1, np
       xs(np+i,1) =  xs(np+1-i,1)
       xs(np+i,3) = -xs(np+1-i,3)
    end do
    xs(:,2) = zero
    xs=0.75d0*xs
  end subroutine calc_geom

  subroutine rbslB(ixI^L,ixO^L,ixK^L,x,a,F,xs,np,Bout)
    integer, intent(in) :: ixI^L,ixO^L,ixK^L,np
    double precision, intent(in) :: a,F
    double precision, dimension(np,3), intent(in) :: xs
    double precision, dimension(ixI^S,3), intent(in) :: x
    double precision, dimension(ixI^S,3), intent(inout) :: Bout
    integer :: ix^D,ixp
    double precision :: re_pi,I_cur,r_mag,KIr,KFr1,KFr2,Rdr
    double precision :: asr,or2,asrr
    double precision, dimension(3) :: Rpl,r_vec,Rcr
    double precision, dimension(ixI^S,3) :: BfrI,BfrF

    ! Direct RBSL magnetic-field formulation based on Titov et al. (2021)
    ! Appendix A/B. This demo computes the flux-rope magnetic field by direct
    ! B-field accumulation instead of vector-potential accumulation plus a
    ! finite-difference curl.
    re_pi=1.d0/dpi
    I_cur=-1.0d0*F*5.0d0*sqrt(2.0d0)/3.0d0/a
    BfrI=zero
    BfrF=zero

    {do ix^DB=ixOmin^DB,ixOmax^DB\}
    do ixp=1,np
      r_vec(1)=(x(ix^D,1)-xs(ixp,1))/a
      r_vec(2)=(x(ix^D,2)-xs(ixp,2))/a
      r_vec(3)=(x(ix^D,3)-xs(ixp,3))/a
      r_mag=sqrt(r_vec(1)**2+r_vec(2)**2+r_vec(3)**2)
      if(ixp .eq. 1) then
        Rpl(1)=0.5d0*(xs(ixp+1,1)-xs(np,1))
        Rpl(2)=0.5d0*(xs(ixp+1,2)-xs(np,2))
        Rpl(3)=0.5d0*(xs(ixp+1,3)-xs(np,3))
      else if(ixp .eq. np) then
        Rpl(1)=0.5d0*(xs(1,1)-xs(ixp-1,1))
        Rpl(2)=0.5d0*(xs(1,2)-xs(ixp-1,2))
        Rpl(3)=0.5d0*(xs(1,3)-xs(ixp-1,3))
      else
        Rpl(1)=0.5d0*(xs(ixp+1,1)-xs(ixp-1,1))
        Rpl(2)=0.5d0*(xs(ixp+1,2)-xs(ixp-1,2))
        Rpl(3)=0.5d0*(xs(ixp+1,3)-xs(ixp-1,3))
      end if
      Rcr(1)=Rpl(2)*r_vec(3)-Rpl(3)*r_vec(2)
      Rcr(2)=Rpl(3)*r_vec(1)-Rpl(1)*r_vec(3)
      Rcr(3)=Rpl(1)*r_vec(2)-Rpl(2)*r_vec(1)
      Rdr=Rpl(1)*r_vec(1)+Rpl(2)*r_vec(2)+Rpl(3)*r_vec(3)
      if(r_mag .le. 1.d-3) then
        KIr=1.698d0
        KFr1=3.902d0
        KFr2=1.73d0
      else if(r_mag .le. one) then
        asr=asin(r_mag)/r_mag
        or2=sqrt(one-r_mag**2)
        asrr=asin((one+two*r_mag*r_mag)/(5.d0-two*r_mag*r_mag))
        KIr=two*re_pi*((asr-or2)/r_mag**2+two*or2)
        KFr1=two*re_pi/r_mag**2*(or2-asr)+8.d0*re_pi*or2&
            +(5.d0-4.d0*r_mag**2)/sqrt(6.d0)*(one-two*re_pi*asrr)
        KFr2=two*re_pi/r_mag**4*(3.d0*asr-(3.d0+two*r_mag**2)*or2)&
            +two/sqrt(6.d0)*(one-two*re_pi*asrr)
      else
        KIr=one/r_mag**3
        KFr1=-one/r_mag**3
        KFr2=3.d0/r_mag**5
      end if
      BfrI(ix^D,:)=BfrI(ix^D,:)+I_cur*0.25d0*re_pi*KIr*Rcr(:)/a**2
      BfrF(ix^D,:)=BfrF(ix^D,:)+F*0.25d0*re_pi*(KFr1*Rpl(:)+KFr2*Rdr*r_vec(:))/a**3
    end do
    {end do\}
    Bout(ixO^S,1:3)=BfrI(ixO^S,1:3)+BfrF(ixO^S,1:3)
  end subroutine rbslB

  subroutine bipoB(ixI^L,ixO^L,x,L_cha,d_cha,q_cha,x_cha,nb,Bout)
    integer, intent(in) :: ixI^L,ixO^L,nb
    double precision, intent(in) :: L_cha,d_cha,q_cha,x_cha
    double precision, dimension(ixI^S,3), intent(in) :: x
    double precision, dimension(ixI^S,3), intent(inout) :: Bout
    integer :: i
    double precision, dimension(nb) :: xpos
    double precision, dimension(ixI^S) :: rpv,rmv
    double precision, dimension(ixI^S,1:3) :: rplus,rminu

    Bout=zero
    do i=1,nb
      xpos(i)=2.d0*(i-1)*x_cha/(nb-1)-x_cha
    end do
    do i=1,nb
      rplus(ixO^S,1)=x(ixO^S,1)-xpos(i)
      rminu(ixO^S,1)=x(ixO^S,1)-xpos(i)
      rplus(ixO^S,2)=x(ixO^S,2)-L_cha
      rminu(ixO^S,2)=x(ixO^S,2)+L_cha
      rplus(ixO^S,3)=x(ixO^S,3)+d_cha
      rminu(ixO^S,3)=x(ixO^S,3)+d_cha
      rpv(ixO^S)=dsqrt(^D&rplus(ixO^S,^D)**2+)
      rmv(ixO^S)=dsqrt(^D&rminu(ixO^S,^D)**2+)
      Bout(ixO^S,1)=Bout(ixO^S,1)+q_cha&
        *(rplus(ixO^S,1)/rpv(ixO^S)**3-rminu(ixO^S,1)/rmv(ixO^S)**3)
      Bout(ixO^S,2)=Bout(ixO^S,2)+q_cha&
        *(rplus(ixO^S,2)/rpv(ixO^S)**3-rminu(ixO^S,2)/rmv(ixO^S)**3)
      Bout(ixO^S,3)=Bout(ixO^S,3)+q_cha&
        *(rplus(ixO^S,3)/rpv(ixO^S)**3-rminu(ixO^S,3)/rmv(ixO^S)**3)
    end do
  endsubroutine bipoB

end module mod_usr
