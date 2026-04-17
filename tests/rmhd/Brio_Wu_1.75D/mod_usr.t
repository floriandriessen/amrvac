module mod_usr
  use mod_mhd
  implicit none

contains

  subroutine usr_init()
    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output

    call set_coordinate_system('Cartesian_1.75D')

    unit_length        = 1.001d10 !1.d8    ! cm
    unit_temperature   = 1.1d3  !1.d4   ! K
    unit_numberdensity = 1.001d10 !1.d9   ! cm^-3
    
    call mhd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
     mhd_gamma = 2.0d0

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    ! initialize one grid 
    integer, intent(in) :: ixG^L,ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision:: rholeft,rhoright,slocx^D
    double precision:: vleft(3),pleft,vright(3),pright
    double precision:: bx,byleft,bzleft,byright,bzright
    logical,save :: first=.true.


    bx=0.75d0 
    rholeft=1.0d0
    pleft=1.0d0 
    vleft=0.0d0
    byleft=1.0d0
    bzleft=zero
    
    rhoright=0.125d0
    pright=0.1d0
    vright=zero
    byright= -1.0d0
    bzright=0.0d0 

    if(first.and.mype==0) then
      print *,'BrioWu test'
      print *,'by=',byright,' bz=',bzright
      first=.false.
    endif

    slocx1=half*(xprobmax1+xprobmin1)
    where({^D&x(ixG^S,^D)<=slocx^D|.or.})
       w(ixG^S,rho_)     = rholeft
       w(ixG^S,mom(1))   = vleft(1)
       w(ixG^S,mom(2))   = vleft(2)
       w(ixG^S,mom(3))   = vleft(3)
       w(ixG^S,p_ )      = pleft
       w(ixG^S,mag(1) )  = bx
       w(ixG^S,mag(2) )  = byleft
       w(ixG^S,mag(3) )  = bzleft
       w(ixG^S,r_e) = const_rad_a*(unit_temperature)**4.d0/unit_pressure
    elsewhere
       w(ixG^S,rho_)     = rhoright
       w(ixG^S,mom(1))   = vright(1)
       w(ixG^S,mom(2))   = vright(2)
       w(ixG^S,mom(3))   = vright(3)
       w(ixG^S,p_  )     = pright
       w(ixG^S,mag(1) )  = bx
       w(ixG^S,mag(2) )  = byright
       w(ixG^S,mag(3) )  = bzright
       w(ixG^S,r_e) = const_rad_a*(0.8*unit_temperature)**4.d0/unit_pressure
    endwhere

    call mhd_to_conserved(ixG^L,ixG^L,w,x)
  
  end subroutine initonegrid_usr

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    use mod_fld
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision                   :: wlocal(ixI^S,nw)
    double precision :: lamb(ixI^S), R(ixI^S)

    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    call fld_get_fluxlimiter(wlocal,x,ixI^L,ixO^L,lamb,R,2)
    w(ixO^S,nw+1)=lamb(ixO^S)
    w(ixO^S,nw+2)=R(ixO^S)
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    character(len=*) :: varnames
    varnames='Lambda R' 

  end subroutine specialvarnames_output

end module mod_usr
