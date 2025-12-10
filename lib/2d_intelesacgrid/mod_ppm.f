module mod_ppm

  implicit none
  private

  public :: PPMlimiter
  public :: PPMlimitervar

contains

  subroutine PPMlimitervar(ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,&
     ixmax1,ixmax2,idims,q,qCT,qLC,qRC)

    ! references:
    ! Mignone et al 2005, ApJS 160, 199,
    ! Miller and Colella 2002, JCP 183, 26
    ! Fryxell et al. 2000 ApJ, 131, 273 (Flash)
    ! baciotti Phd (http://www.aei.mpg.de/~baiotti/Baiotti_PhD.pdf)
    ! old version : april 2009
    ! author: zakaria.meliani@wis.kuleuven.be
    ! current version : March 2023 by Chun Xia

    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, ixmin1,&
       ixmin2,ixmax1,ixmax2, idims
    double precision, intent(in)    :: q(ixImin1:ixImax1,ixImin2:ixImax2),&
       qCT(ixImin1:ixImax1,ixImin2:ixImax2)

    double precision, intent(inout) :: qRC(ixImin1:ixImax1,ixImin2:ixImax2),&
       qLC(ixImin1:ixImax1,ixImin2:ixImax2)

    double precision,dimension(ixImin1:ixImax1,ixImin2:ixImax2)  :: dqC,d2qC,&
       ldq
    double precision,dimension(ixImin1:ixImax1,ixImin2:ixImax2)  :: qMin,qMax,&
       tmp

    integer   :: lxCmin1,lxCmin2,lxCmax1,lxCmax2,lxRmin1,lxRmin2,lxRmax1,&
       lxRmax2,lxLmin1,lxLmin2,lxLmax1,lxLmax2
    integer   :: ixLmin1,ixLmin2,ixLmax1,ixLmax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,ixRmin1,ixRmin2,ixRmax1,ixRmax2,ixOLmin1,ixOLmin2,ixOLmax1,&
       ixOLmax2,ixORmin1,ixORmin2,ixORmax1,ixORmax2
    integer   :: hxLmin1,hxLmin2,hxLmax1,hxLmax2,hxCmin1,hxCmin2,hxCmax1,&
       hxCmax2,hxRmin1,hxRmin2,hxRmax1,hxRmax2
    integer   :: kxLmin1,kxLmin2,kxLmax1,kxLmax2,kxCmin1,kxCmin2,kxCmax1,&
       kxCmax2,kxRmin1,kxRmin2,kxRmax1,kxRmax2

    ixOmin1=ixmin1-kr(idims,1);ixOmin2=ixmin2-kr(idims,2)
    ixOmax1=ixmax1+kr(idims,1);ixOmax2=ixmax2+kr(idims,2); !ixO[ixMmin-1,ixMmax+1]
    ixOLmin1=ixOmin1-kr(idims,1);ixOLmin2=ixOmin2-kr(idims,2)
    ixOLmax1=ixOmax1;ixOLmax2=ixOmax2; !ixOL[ixMmin-2,ixMmax+1]
    ixORmin1=ixOLmin1+kr(idims,1);ixORmin2=ixOLmin2+kr(idims,2)
    ixORmax1=ixOLmax1+kr(idims,1);ixORmax2=ixOLmax2+kr(idims,2); !ixOR=[iMmin-1,ixMmax+2]
    ixLmin1=ixOmin1-kr(idims,1);ixLmin2=ixOmin2-kr(idims,2)
    ixLmax1=ixOmax1-kr(idims,1);ixLmax2=ixOmax2-kr(idims,2); !ixL[ixMmin-2,ixMmax]
    ixRmin1=ixOmin1+kr(idims,1);ixRmin2=ixOmin2+kr(idims,2)
    ixRmax1=ixOmax1+kr(idims,1);ixRmax2=ixOmax2+kr(idims,2); !ixR=[iMmin,ixMmax+2]

    hxCmin1=ixOmin1;hxCmin2=ixOmin2;hxCmax1=ixmax1;hxCmax2=ixmax2; !hxC = [ixMmin-1,ixMmax]
    hxLmin1=hxCmin1-kr(idims,1);hxLmin2=hxCmin2-kr(idims,2)
    hxLmax1=hxCmax1-kr(idims,1);hxLmax2=hxCmax2-kr(idims,2); !hxL = [ixMmin-2,ixMmax-1]
    hxRmin1=hxCmin1+kr(idims,1);hxRmin2=hxCmin2+kr(idims,2)
    hxRmax1=hxCmax1+kr(idims,1);hxRmax2=hxCmax2+kr(idims,2); !hxR = [ixMmin,ixMmax+1]

    kxCmin1=ixLmin1-1;kxCmin2=ixLmin2-1; kxCmax1=ixRmax1;kxCmax2=ixRmax2; !kxC=[iMmin-3,ixMmax+2]
    kxRmin1=kxCmin1+kr(idims,1);kxRmin2=kxCmin2+kr(idims,2)
    kxRmax1=kxCmax1+kr(idims,1);kxRmax2=kxCmax2+kr(idims,2); !kxR=[iMmin-2,ixMmax+3]

    lxCmin1=ixOmin1-kr(idims,1);lxCmin2=ixOmin2-kr(idims,2)
    lxCmax1=ixOmax1+kr(idims,1);lxCmax2=ixOmax2+kr(idims,2); !lxC=[iMmin-2,ixMmax+2]
    lxLmin1=lxCmin1-kr(idims,1);lxLmin2=lxCmin2-kr(idims,2)
    lxLmax1=lxCmax1-kr(idims,1);lxLmax2=lxCmax2-kr(idims,2); !lxL=[iMmin-3,ixMmax+1]
    lxRmin1=lxCmin1+kr(idims,1);lxRmin2=lxCmin2+kr(idims,2)
    lxRmax1=lxCmax1+kr(idims,1);lxRmax2=lxCmax2+kr(idims,2); !lxR=[iMmin-1,ixMmax+3]

    dqC(kxCmin1:kxCmax1,kxCmin2:kxCmax2)=q(kxRmin1:kxRmax1,&
       kxRmin2:kxRmax2)-q(kxCmin1:kxCmax1,kxCmin2:kxCmax2)
    ! Eq. 64,  Miller and Colella 2002, JCP 183, 26 
    d2qC(kxCmin1:kxCmax1,kxCmin2:kxCmax2)=half*(q(lxRmin1:lxRmax1,&
       lxRmin2:lxRmax2)-q(lxLmin1:lxLmax1,lxLmin2:lxLmax2))
    where(dqC(lxCmin1:lxCmax1,lxCmin2:lxCmax2)*dqC(lxLmin1:lxLmax1,&
       lxLmin2:lxLmax2)>zero)
       ! Store the sign of dwC in wMin
       qMin(kxCmin1:kxCmax1,kxCmin2:kxCmax2)= sign(one,d2qC(lxCmin1:lxCmax1,&
          lxCmin2:lxCmax2))
       ! Eq. 65,  Miller and Colella 2002, JCP 183, 26 
       ldq(lxCmin1:lxCmax1,lxCmin2:lxCmax2)= qMin(lxCmin1:lxCmax1,&
          lxCmin2:lxCmax2)*min(dabs(d2qC(lxCmin1:lxCmax1,lxCmin2:lxCmax2)),&
          2.0d0*dabs(dqC(lxLmin1:lxLmax1,lxLmin2:lxLmax2)),&
          2.0d0*dabs(dqC(lxCmin1:lxCmax1,lxCmin2:lxCmax2)))
    else where
       ldq(lxCmin1:lxCmax1,lxCmin2:lxCmax2)=zero
    end where

    ! Eq. 66, Miller and Colella 2002, JCP 183, 26
    qLC(ixOLmin1:ixOLmax1,ixOLmin2:ixOLmax2)=qLC(ixOLmin1:ixOLmax1,&
       ixOLmin2:ixOLmax2)+half*dqC(ixOLmin1:ixOLmax1,&
       ixOLmin2:ixOLmax2)+(ldq(ixOLmin1:ixOLmax1,&
       ixOLmin2:ixOLmax2)-ldq(ixORmin1:ixORmax1,ixORmin2:ixORmax2))/6.0d0
    qRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2)=qLC(ixLmin1:ixLmax1,ixLmin2:ixLmax2)

    ! make sure that min wCT(i)<wLC(i)<wCT(i+1) same for wRC(i)
    call extremaq(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,qCT,1,qMax,qMin)

    ! Eq. B8, page 217, Mignone et al 2005, ApJS
    qRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2)=max(qMin(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2),min(qMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
       qRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2)))
    qLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=max(qMin(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2),min(qMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
       qLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))

    ! Eq. B9, page 217, Mignone et al 2005, ApJS
    where((qRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2)-qCT(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))*(qCT(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)-qLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2))<=zero)
       qRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2)=qCT(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)
       qLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=qCT(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)
    end where

    qMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(qLC(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)-qRC(ixLmin1:ixLmax1,&
       ixLmin2:ixLmax2))*(qCT(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)-half*(qLC(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)+qRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2)))
    qMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(qLC(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)-qRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2))**2/6.0d0

    tmp(hxLmin1:hxLmax1,hxLmin2:hxLmax2)=qRC(hxLmin1:hxLmax1,hxLmin2:hxLmax2)
    ! Eq. B10, page 218, Mignone et al 2005, ApJS
    where(qMax(hxRmin1:hxRmax1,hxRmin2:hxRmax2)>qMin(hxRmin1:hxRmax1,&
       hxRmin2:hxRmax2))
       qRC(hxCmin1:hxCmax1,hxCmin2:hxCmax2)= 3.0d0*qCT(hxRmin1:hxRmax1,&
          hxRmin2:hxRmax2)-2.0d0*qLC(hxRmin1:hxRmax1,hxRmin2:hxRmax2)
    end where
    ! Eq. B11, page 218, Mignone et al 2005, ApJS
    where(qMax(hxCmin1:hxCmax1,hxCmin2:hxCmax2)<-qMin(hxCmin1:hxCmax1,&
       hxCmin2:hxCmax2))
       qLC(hxCmin1:hxCmax1,hxCmin2:hxCmax2)= 3.0d0*qCT(hxCmin1:hxCmax1,&
          hxCmin2:hxCmax2)-2.0d0*tmp(hxLmin1:hxLmax1,hxLmin2:hxLmax2)
    end where

  end subroutine PPMlimitervar

  subroutine PPMlimiter(ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
     ixmax2,idims,w,wCT,wLC,wRC)

    ! references:
    ! Mignone et al 2005, ApJS 160, 199, 
    ! Miller and Colella 2002, JCP 183, 26 
    ! Fryxell et al. 2000 ApJ, 131, 273 (Flash)
    ! baciotti Phd (http://www.aei.mpg.de/~baiotti/Baiotti_PhD.pdf)
    ! old version : april 2009
    ! author: zakaria.meliani@wis.kuleuven.be
    ! current version : March 2023 by Chun Xia
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, ixmin1,&
       ixmin2,ixmax1,ixmax2, idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
       wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision, intent(inout) :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw),wLC(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) 

    double precision,dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nwflux)  :: dwC,d2wC,ldw
    double precision,dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nwflux)  :: wMin,wMax,tmp
    double precision,dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: aa, ab, ac,&
        ki

    integer   :: lxCmin1,lxCmin2,lxCmax1,lxCmax2,lxLmin1,lxLmin2,lxLmax1,&
       lxLmax2,lxRmin1,lxRmin2,lxRmax1,lxRmax2
    integer   :: ixLLmin1,ixLLmin2,ixLLmax1,ixLLmax2,ixLmin1,ixLmin2,ixLmax1,&
       ixLmax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,ixRmin1,ixRmin2,ixRmax1,ixRmax2,&
       ixRRmin1,ixRRmin2,ixRRmax1,ixRRmax2,ixOLmin1,ixOLmin2,ixOLmax1,ixOLmax2,&
       ixORmin1,ixORmin2,ixORmax1,ixORmax2
    integer   :: hxLmin1,hxLmin2,hxLmax1,hxLmax2,hxCmin1,hxCmin2,hxCmax1,&
       hxCmax2,hxRmin1,hxRmin2,hxRmax1,hxRmax2
    integer   :: kxLLmin1,kxLLmin2,kxLLmax1,kxLLmax2,kxLmin1,kxLmin2,kxLmax1,&
       kxLmax2,kxCmin1,kxCmin2,kxCmax1,kxCmax2,kxRmin1,kxRmin2,kxRmax1,kxRmax2,&
       kxRRmin1,kxRRmin2,kxRRmax1,kxRRmax2
    integer   :: iw, idimss

    double precision, parameter :: betamin=0.75d0, betamax=0.85d0,Zmin=0.25d0,&
        Zmax=0.75d0,eta1=20.0d0,eta2=0.05d0,eps=0.01d0,kappa=0.1d0

    ixOmin1=ixmin1-kr(idims,1);ixOmin2=ixmin2-kr(idims,2)
    ixOmax1=ixmax1+kr(idims,1);ixOmax2=ixmax2+kr(idims,2); !ixO[ixMmin-1,ixMmax+1]
    ixOLmin1=ixOmin1-kr(idims,1);ixOLmin2=ixOmin2-kr(idims,2)
    ixOLmax1=ixOmax1;ixOLmax2=ixOmax2; !ixOL[ixMmin-2,ixMmax+1]
    ixORmin1=ixOLmin1+kr(idims,1);ixORmin2=ixOLmin2+kr(idims,2)
    ixORmax1=ixOLmax1+kr(idims,1);ixORmax2=ixOLmax2+kr(idims,2); !ixOR=[iMmin-1,ixMmax+2]
    ixLmin1=ixOmin1-kr(idims,1);ixLmin2=ixOmin2-kr(idims,2)
    ixLmax1=ixOmax1-kr(idims,1);ixLmax2=ixOmax2-kr(idims,2); !ixL[ixMmin-2,ixMmax]
    ixRmin1=ixOmin1+kr(idims,1);ixRmin2=ixOmin2+kr(idims,2)
    ixRmax1=ixOmax1+kr(idims,1);ixRmax2=ixOmax2+kr(idims,2); !ixR=[iMmin,ixMmax+2]

    hxCmin1=ixOmin1;hxCmin2=ixOmin2;hxCmax1=ixmax1;hxCmax2=ixmax2; !hxC = [ixMmin-1,ixMmax]
    hxLmin1=hxCmin1-kr(idims,1);hxLmin2=hxCmin2-kr(idims,2)
    hxLmax1=hxCmax1-kr(idims,1);hxLmax2=hxCmax2-kr(idims,2); !hxL = [ixMmin-2,ixMmax-1]
    hxRmin1=hxCmin1+kr(idims,1);hxRmin2=hxCmin2+kr(idims,2)
    hxRmax1=hxCmax1+kr(idims,1);hxRmax2=hxCmax2+kr(idims,2); !hxR = [ixMmin,ixMmax+1]

    kxCmin1=ixLmin1-1;kxCmin2=ixLmin2-1; kxCmax1=ixRmax1;kxCmax2=ixRmax2; !kxC=[iMmin-3,ixMmax+2]
    kxRmin1=kxCmin1+kr(idims,1);kxRmin2=kxCmin2+kr(idims,2)
    kxRmax1=kxCmax1+kr(idims,1);kxRmax2=kxCmax2+kr(idims,2); !kxR=[iMmin-2,ixMmax+3]

    lxCmin1=ixOmin1-kr(idims,1);lxCmin2=ixOmin2-kr(idims,2)
    lxCmax1=ixOmax1+kr(idims,1);lxCmax2=ixOmax2+kr(idims,2); !lxC=[iMmin-2,ixMmax+2]
    lxLmin1=lxCmin1-kr(idims,1);lxLmin2=lxCmin2-kr(idims,2)
    lxLmax1=lxCmax1-kr(idims,1);lxLmax2=lxCmax2-kr(idims,2); !lxL=[iMmin-3,ixMmax+1]
    lxRmin1=lxCmin1+kr(idims,1);lxRmin2=lxCmin2+kr(idims,2)
    lxRmax1=lxCmax1+kr(idims,1);lxRmax2=lxCmax2+kr(idims,2); !lxR=[iMmin-1,ixMmax+3]

    dwC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nwflux)=w(kxRmin1:kxRmax1,&
       kxRmin2:kxRmax2,1:nwflux)-w(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nwflux)
    ! Eq. 64,  Miller and Colella 2002, JCP 183, 26 
    d2wC(lxCmin1:lxCmax1,lxCmin2:lxCmax2,1:nwflux)=half*(w(lxRmin1:lxRmax1,&
       lxRmin2:lxRmax2,1:nwflux)-w(lxLmin1:lxLmax1,lxLmin2:lxLmax2,1:nwflux))
    where(dwC(lxCmin1:lxCmax1,lxCmin2:lxCmax2,1:nwflux)*dwC(lxLmin1:lxLmax1,&
       lxLmin2:lxLmax2,1:nwflux)>zero)
       ! Store the sign of dwC in wMin
       wMin(lxCmin1:lxCmax1,lxCmin2:lxCmax2,1:nwflux)= sign(one,&
          d2wC(lxCmin1:lxCmax1,lxCmin2:lxCmax2,1:nwflux))
       ! Eq. 65,  Miller and Colella 2002, JCP 183, 26 
       ldw(lxCmin1:lxCmax1,lxCmin2:lxCmax2,1:nwflux)= wMin(lxCmin1:lxCmax1,&
          lxCmin2:lxCmax2,1:nwflux)*min(dabs(d2wC(lxCmin1:lxCmax1,&
          lxCmin2:lxCmax2,1:nwflux)),2.0d0*dabs(dwC(lxLmin1:lxLmax1,&
          lxLmin2:lxLmax2,1:nwflux)),2.0d0*dabs(dwC(lxCmin1:lxCmax1,&
          lxCmin2:lxCmax2,1:nwflux)))
    else where
       ldw(lxCmin1:lxCmax1,lxCmin2:lxCmax2,1:nwflux)=zero
    end where

    ! Eq. 66,  Miller and Colella 2002, JCP 183, 26 
    wLC(ixOLmin1:ixOLmax1,ixOLmin2:ixOLmax2,1:nwflux)=wLC(ixOLmin1:ixOLmax1,&
       ixOLmin2:ixOLmax2,1:nwflux)+half*dwC(ixOLmin1:ixOLmax1,&
       ixOLmin2:ixOLmax2,1:nwflux)+(ldw(ixOLmin1:ixOLmax1,ixOLmin2:ixOLmax2,&
       1:nwflux)-ldw(ixORmin1:ixORmax1,ixORmin2:ixORmax2,1:nwflux))/6.0d0
    wRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,1:nwflux)=wLC(ixLmin1:ixLmax1,&
       ixLmin2:ixLmax2,1:nwflux)

    ! make sure that min w(i)<wLC(i)<w(i+1) same for wRC(i)
    call extremaw(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,w,1,wMax,wMin)

    ! Eq. B8, page 217, Mignone et al 2005, ApJS
    wRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,1:nwflux)=max(wMin(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1:nwflux),min(wMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1:nwflux),wRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,1:nwflux))) 
    wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)=max(wMin(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1:nwflux),min(wMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1:nwflux),wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)))

    ! Eq. B9, page 217, Mignone et al 2005, ApJS
    where((wRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,1:nwflux)-wCT(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1:nwflux))*(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1:nwflux)-wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux))<=zero)
       wRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,1:nwflux)=wCT(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,1:nwflux)
       wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)=wCT(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,1:nwflux)
    end where

    wMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)=(wLC(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1:nwflux)-wRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
       1:nwflux))*(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1:nwflux)-half*(wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1:nwflux)+wRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,1:nwflux)))
    wMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)=(wLC(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1:nwflux)-wRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
       1:nwflux))**2/6.0d0
    tmp(hxLmin1:hxLmax1,hxLmin2:hxLmax2,1:nwflux)=wRC(hxLmin1:hxLmax1,&
       hxLmin2:hxLmax2,1:nwflux)
    ! Eq. B10, page 218, Mignone et al 2005, ApJS
    where(wMax(hxRmin1:hxRmax1,hxRmin2:hxRmax2,1:nwflux)>wMin(hxRmin1:hxRmax1,&
       hxRmin2:hxRmax2,1:nwflux))
       wRC(hxCmin1:hxCmax1,hxCmin2:hxCmax2,&
          1:nwflux)= 3.0d0*wCT(hxRmin1:hxRmax1,hxRmin2:hxRmax2,&
          1:nwflux)-2.0d0*wLC(hxRmin1:hxRmax1,hxRmin2:hxRmax2,1:nwflux)
    end where
    ! Eq. B11, page 218, Mignone et al 2005, ApJS
    where(wMax(hxCmin1:hxCmax1,hxCmin2:hxCmax2,1:nwflux)<-wMin(hxCmin1:hxCmax1,&
       hxCmin2:hxCmax2,1:nwflux))
       wLC(hxCmin1:hxCmax1,hxCmin2:hxCmax2,&
          1:nwflux)= 3.0d0*wCT(hxCmin1:hxCmax1,hxCmin2:hxCmax2,&
          1:nwflux)-2.0d0*tmp(hxLmin1:hxLmax1,hxLmin2:hxLmax2,1:nwflux)
    end where

    ! flattening at the contact discontinuities
    if(flatcd)then
      ixRRmin1=ixRmin1+kr(idims,1);ixRRmin2=ixRmin2+kr(idims,2)
      ixRRmax1=ixRmax1+kr(idims,1);ixRRmax2=ixRmax2+kr(idims,2); !ixRR=[iMmin+1,ixMmax+3]
      kxLmin1=kxCmin1-kr(idims,1);kxLmin2=kxCmin2-kr(idims,2)
      kxLmax1=kxCmax1-kr(idims,1);kxLmax2=kxCmax2-kr(idims,2); !kxL=[iMmin-4,ixMmax+1]
      call ppm_flatcd(ixImin1,ixImin2,ixImax1,ixImax2,kxCmin1,kxCmin2,kxCmax1,&
         kxCmax2,kxLmin1,kxLmin2,kxLmax1,kxLmax2,kxRmin1,kxRmin2,kxRmax1,&
         kxRmax2,wCT,d2wC,aa,ab)
      if(any(kappa*aa(kxCmin1:kxCmax1,kxCmin2:kxCmax2)>=ab(kxCmin1:kxCmax1,&
         kxCmin2:kxCmax2)))then
        do iw=1,nwflux
          where(kappa*aa(kxCmin1:kxCmax1,kxCmin2:kxCmax2)>=ab(kxCmin1:kxCmax1,&
             kxCmin2:kxCmax2).and. dabs(dwC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,&
             iw))>smalldouble)
            wMax(kxCmin1:kxCmax1,kxCmin2:kxCmax2,iw) = wCT(kxRmin1:kxRmax1,&
               kxRmin2:kxRmax2,iw)-2.0d0*wCT(kxCmin1:kxCmax1,kxCmin2:kxCmax2,&
               iw)+wCT(kxLmin1:kxLmax1,kxLmin2:kxLmax2,iw)
          end where

          where(wMax(ixRmin1:ixRmax1,ixRmin2:ixRmax2,iw)*wMax(ixLmin1:ixLmax1,&
             ixLmin2:ixLmax2,iw)<zero .and.dabs(wCT(ixRmin1:ixRmax1,&
             ixRmin2:ixRmax2,iw)-wCT(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
             iw))-eps*min(dabs(wCT(ixRmin1:ixRmax1,ixRmin2:ixRmax2,iw)),&
             dabs(wCT(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
             iw)))>zero .and. kappa*aa(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)>=ab(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2).and. dabs(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             iw))>smalldouble)

            ac(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(wCT(ixLLmin1:ixLLmax1,&
               ixLLmin2:ixLLmax2,iw)-wCT(ixRRmin1:ixRRmax1,ixRRmin2:ixRRmax2,&
               iw)+4.0d0*dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               iw))/(12.0d0*dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw))
            wMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=max(zero,&
               min(eta1*(ac(ixOmin1:ixOmax1,ixOmin2:ixOmax2)-eta2),one))
          else where
            wMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=zero
          end where

          where(wMin(hxCmin1:hxCmax1,hxCmin2:hxCmax2,iw)>zero)
            wLC(hxCmin1:hxCmax1,hxCmin2:hxCmax2,iw) = wLC(hxCmin1:hxCmax1,&
               hxCmin2:hxCmax2,iw)*(one-wMin(hxCmin1:hxCmax1,hxCmin2:hxCmax2,&
               iw))+(wCT(hxCmin1:hxCmax1,hxCmin2:hxCmax2,&
               iw)+half*ldw(hxCmin1:hxCmax1,hxCmin2:hxCmax2,&
               iw))*wMin(hxCmin1:hxCmax1,hxCmin2:hxCmax2,iw)
          end where
          where(wMin(hxRmin1:hxRmax1,hxRmin2:hxRmax2,iw)>zero)
            wRC(hxCmin1:hxCmax1,hxCmin2:hxCmax2,iw) = wRC(hxCmin1:hxCmax1,&
               hxCmin2:hxCmax2,iw)*(one-wMin(hxRmin1:hxRmax1,hxRmin2:hxRmax2,&
               iw))+(wCT(hxRmin1:hxRmax1,hxRmin2:hxRmax2,&
               iw)-half*ldw(hxRmin1:hxRmax1,hxRmin2:hxRmax2,&
               iw))*wMin(hxRmin1:hxRmax1,hxRmin2:hxRmax2,iw)
          end where
        end do
      end if
    end if

    ! flattening at the shocks
    if(flatsh)then
       ! following MILLER and COLELLA 2002 JCP 183, 26
       kxCmin1=ixmin1-2;kxCmin2=ixmin2-2; kxCmax1=ixmax1+2;kxCmax2=ixmax2+2; !kxC=[ixMmin-2,ixMmax+2]
       ki=bigdouble
       do idimss=1,ndim
          kxLmin1=kxCmin1-kr(idimss,1);kxLmin2=kxCmin2-kr(idimss,2)
          kxLmax1=kxCmax1-kr(idimss,1);kxLmax2=kxCmax2-kr(idimss,2); !kxL=[ixMmin-3,ixMmax+1]
          kxRmin1=kxCmin1+kr(idimss,1);kxRmin2=kxCmin2+kr(idimss,2)
          kxRmax1=kxCmax1+kr(idimss,1);kxRmax2=kxCmax2+kr(idimss,2); !kxR=[ixMmin-1,ixMmax+3]
          kxLLmin1=kxLmin1-kr(idimss,1);kxLLmin2=kxLmin2-kr(idimss,2)
          kxLLmax1=kxLmax1-kr(idimss,1);kxLLmax2=kxLmax2-kr(idimss,2); !kxLL=[ixMmin-4,ixMmax]
          kxRRmin1=kxRmin1+kr(idimss,1);kxRRmin2=kxRmin2+kr(idimss,2)
          kxRRmax1=kxRmax1+kr(idimss,1);kxRRmax2=kxRmax2+kr(idimss,2); !kxRR=[ixMmin,ixMmax+4]

          call ppm_flatsh(ixImin1,ixImin2,ixImax1,ixImax2,kxCmin1,kxCmin2,&
             kxCmax1,kxCmax2,kxLLmin1,kxLLmin2,kxLLmax1,kxLLmax2,kxLmin1,&
             kxLmin2,kxLmax1,kxLmax2,kxRmin1,kxRmin2,kxRmax1,kxRmax2,kxRRmin1,&
             kxRRmin2,kxRRmax1,kxRRmax2,idimss,wCT,aa,ab)

          ! eq. B17, page 218, Mignone et al 2005, ApJS (Xi1 min)
          ac(kxCmin1:kxCmax1,kxCmin2:kxCmax2) = max(zero,min(one,&
             (betamax-aa(kxCmin1:kxCmax1,kxCmin2:kxCmax2))/(betamax-betamin)))
          ! eq. B18, page 218, Mignone et al 2005, ApJS (Xi1)
          ! recycling aa(ixL^S)
          where(wCT(kxRmin1:kxRmax1,kxRmin2:kxRmax2,&
             iw_mom(idimss))<wCT(kxLmin1:kxLmax1,kxLmin2:kxLmax2,&
             iw_mom(idimss)))
             aa(kxCmin1:kxCmax1,kxCmin2:kxCmax2) = max(ac(kxCmin1:kxCmax1,&
                kxCmin2:kxCmax2), min(one,(Zmax-ab(kxCmin1:kxCmax1,&
                kxCmin2:kxCmax2))/(Zmax-Zmin)))
          else where
             aa(kxCmin1:kxCmax1,kxCmin2:kxCmax2) = one
          end where
          ixLmin1=ixOmin1-kr(idimss,1);ixLmin2=ixOmin2-kr(idimss,2)
          ixLmax1=ixOmax1-kr(idimss,1);ixLmax2=ixOmax2-kr(idimss,2); !ixL=[ixMmin-2,ixMmax]
          ixRmin1=ixOmin1+kr(idimss,1);ixRmin2=ixOmin2+kr(idimss,2)
          ixRmax1=ixOmax1+kr(idimss,1);ixRmax2=ixOmax2+kr(idimss,2); !ixR=[ixMmin,ixMmax+2]
          ki(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=min(ki(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2),aa(ixLmin1:ixLmax1,ixLmin2:ixLmax2),&
             aa(ixOmin1:ixOmax1,ixOmin2:ixOmax2),aa(ixRmin1:ixRmax1,&
             ixRmin2:ixRmax2))
       end do
       ! recycling wMax
       do iw=1,nwflux
          where(dabs(ab(ixOmin1:ixOmax1,ixOmin2:ixOmax2)-one)>smalldouble)
             wMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                iw) = (one-ab(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)
          end where

          where(dabs(ab(hxCmin1:hxCmax1,hxCmin2:hxCmax2)-one)>smalldouble)
             wLC(hxCmin1:hxCmax1,hxCmin2:hxCmax2,iw) = ab(hxCmin1:hxCmax1,&
                hxCmin2:hxCmax2)*wLC(hxCmin1:hxCmax1,hxCmin2:hxCmax2,&
                iw)+wMax(hxCmin1:hxCmax1,hxCmin2:hxCmax2,iw)
          end where

          where(dabs(ab(hxRmin1:hxRmax1,hxRmin2:hxRmax2)-one)>smalldouble)
             wRC(hxCmin1:hxCmax1,hxCmin2:hxCmax2,iw) = ab(hxRmin1:hxRmax1,&
                hxRmin2:hxRmax2)*wRC(hxCmin1:hxCmax1,hxCmin2:hxCmax2,&
                iw)+wMax(hxRmin1:hxRmax1,hxRmin2:hxRmax2,iw)
          end where
       end do
    end if

  end subroutine PPMlimiter

  subroutine ppm_flatcd(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,ixLmin1,ixLmin2,ixLmax1,ixLmax2,ixRmin1,ixRmin2,ixRmax1,&
     ixRmax2,w,d2w,drho,dp)
    use mod_global_parameters
    use mod_physics, only: phys_gamma

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, ixLmin1,ixLmin2,ixLmax1,ixLmax2,&
        ixRmin1,ixRmin2,ixRmax1,ixRmax2
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw),&
        d2w(ixGlo1:ixGhi1,ixGlo2:ixGhi2, 1:nwflux)
    double precision, intent(inout) :: drho(ixGlo1:ixGhi1,ixGlo2:ixGhi2),&
        dp(ixGlo1:ixGhi1,ixGlo2:ixGhi2)

    drho(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = phys_gamma*abs(d2w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, iw_rho))/min(w(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
        iw_rho), w(ixRmin1:ixRmax1,ixRmin2:ixRmax2, iw_rho))
    dp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = abs(d2w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, iw_e))/min(w(ixLmin1:ixLmax1,ixLmin2:ixLmax2, iw_e),&
        w(ixRmin1:ixRmax1,ixRmin2:ixRmax2, iw_e))

  end subroutine ppm_flatcd

  subroutine ppm_flatsh(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,ixLLmin1,ixLLmin2,ixLLmax1,ixLLmax2,ixLmin1,ixLmin2,&
     ixLmax1,ixLmax2,ixRmin1,ixRmin2,ixRmax1,ixRmax2,ixRRmin1,ixRRmin2,&
     ixRRmax1,ixRRmax2,idims,w,drho,dp)
    use mod_global_parameters
    use mod_physics, only: phys_gamma

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, ixLLmin1,ixLLmin2,ixLLmax1,ixLLmax2,&
        ixLmin1,ixLmin2,ixLmax1,ixLmax2, ixRmin1,ixRmin2,ixRmax1,ixRmax2,&
        ixRRmin1,ixRRmin2,ixRRmax1,ixRRmax2
    integer, intent(in)             :: idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(inout) :: drho(ixImin1:ixImax1,ixImin2:ixImax2),&
        dp(ixImin1:ixImax1,ixImin2:ixImax2)

    ! eq. B15, page 218, Mignone and Bodo 2005, ApJS (beta1)
    where (abs(w(ixRRmin1:ixRRmax1,ixRRmin2:ixRRmax2,&
        iw_e)-w(ixLLmin1:ixLLmax1,ixLLmin2:ixLLmax2, iw_e))>smalldouble)
       drho(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = abs((w(ixRmin1:ixRmax1,&
          ixRmin2:ixRmax2, iw_e)-w(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
           iw_e))/(w(ixRRmin1:ixRRmax1,ixRRmin2:ixRRmax2,&
           iw_e)-w(ixLLmin1:ixLLmax1,ixLLmin2:ixLLmax2, iw_e)))
    else where
       drho(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = zero
    end where

    ! eq. 76, page 48, Miller and Colella 2002, JCoPh, adjusted
    dp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = abs(w(ixRmin1:ixRmax1,&
       ixRmin2:ixRmax2, iw_e)-w(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
        iw_e))/(phys_gamma*(w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
        iw_e)+w(ixLmin1:ixLmax1,ixLmin2:ixLmax2, iw_e)))

  end subroutine ppm_flatsh

  subroutine extremaq(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
     ixOmax2,q,nshift,qMax,qMin)

    use mod_global_parameters

    integer,intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: q(ixImin1:ixImax1,ixImin2:ixImax2)
    integer,intent(in)           :: nshift

    double precision, intent(out) :: qMax(ixImin1:ixImax1,ixImin2:ixImax2),&
       qMin(ixImin1:ixImax1,ixImin2:ixImax2)

    integer           :: ixsmin1,ixsmin2,ixsmax1,ixsmax2,ixsRmin1,ixsRmin2,&
       ixsRmax1,ixsRmax2,ixsLmin1,ixsLmin2,ixsLmax1,ixsLmax2,idims,jdims,kdims,&
       ishift,i,j

    do ishift=1,nshift
      idims=1
      ixsRmin1=ixOmin1+ishift*kr(idims,1);ixsRmin2=ixOmin2+ishift*kr(idims,2)
      ixsRmax1=ixOmax1+ishift*kr(idims,1);ixsRmax2=ixOmax2+ishift*kr(idims,2);
      ixsLmin1=ixOmin1-ishift*kr(idims,1);ixsLmin2=ixOmin2-ishift*kr(idims,2)
      ixsLmax1=ixOmax1-ishift*kr(idims,1);ixsLmax2=ixOmax2-ishift*kr(idims,2);
      if (ishift==1) then
        qMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=max(q(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2),q(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2),&
           q(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2))
        qMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=min(q(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2),q(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2),&
           q(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2))
      else
        qMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=max(qMax(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2),q(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2),&
           q(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2))
        qMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=min(qMin(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2),q(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2),&
           q(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2))
      end if
      
      idims=1
      jdims=idims+1
      do i=-1,1
        ixsmin1=ixOmin1+i*ishift*kr(idims,1)
        ixsmin2=ixOmin2+i*ishift*kr(idims,2)
        ixsmax1=ixOmax1+i*ishift*kr(idims,1)
        ixsmax2=ixOmax2+i*ishift*kr(idims,2);
        ixsRmin1=ixsmin1+ishift*kr(jdims,1)
        ixsRmin2=ixsmin2+ishift*kr(jdims,2)
        ixsRmax1=ixsmax1+ishift*kr(jdims,1)
        ixsRmax2=ixsmax2+ishift*kr(jdims,2);
        ixsLmin1=ixsmin1-ishift*kr(jdims,1)
        ixsLmin2=ixsmin2-ishift*kr(jdims,2)
        ixsLmax1=ixsmax1-ishift*kr(jdims,1)
        ixsLmax2=ixsmax2-ishift*kr(jdims,2);
        qMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=max(qMax(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2),q(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2),&
           q(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2))
        qMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=min(qMin(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2),q(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2),&
           q(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2))
      end do
     
      
    enddo

  end subroutine  extremaq

  subroutine extremaw(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
     ixOmax2,w,nshift,wMax,wMin)
    use mod_global_parameters

    integer,intent(in)            :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    integer,intent(in)            :: nshift

    double precision, intent(out) :: wMax(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nwflux),wMin(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux)

    integer          :: ixsmin1,ixsmin2,ixsmax1,ixsmax2,ixsRmin1,ixsRmin2,&
       ixsRmax1,ixsRmax2,ixsLmin1,ixsLmin2,ixsLmax1,ixsLmax2,idims,jdims,kdims,&
       ishift,i,j

    do ishift=1,nshift
      idims=1
      ixsRmin1=ixOmin1+ishift*kr(idims,1);ixsRmin2=ixOmin2+ishift*kr(idims,2)
      ixsRmax1=ixOmax1+ishift*kr(idims,1);ixsRmax2=ixOmax2+ishift*kr(idims,2);
      ixsLmin1=ixOmin1-ishift*kr(idims,1);ixsLmin2=ixOmin2-ishift*kr(idims,2)
      ixsLmax1=ixOmax1-ishift*kr(idims,1);ixsLmax2=ixOmax2-ishift*kr(idims,2);
      if (ishift==1) then
        wMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)= max(w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,1:nwflux),w(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2,&
           1:nwflux),w(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2,1:nwflux))
        wMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)= min(w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,1:nwflux),w(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2,&
           1:nwflux),w(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2,1:nwflux))
      else
        wMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1:nwflux)= max(wMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux),&
           w(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2,1:nwflux),w(ixsLmin1:ixsLmax1,&
           ixsLmin2:ixsLmax2,1:nwflux))
        wMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1:nwflux)= min(wMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux),&
           w(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2,1:nwflux),w(ixsLmin1:ixsLmax1,&
           ixsLmin2:ixsLmax2,1:nwflux))
      end if
      
      idims=1
      jdims=idims+1
      do i=-1,1
        ixsmin1=ixOmin1+i*ishift*kr(idims,1)
        ixsmin2=ixOmin2+i*ishift*kr(idims,2)
        ixsmax1=ixOmax1+i*ishift*kr(idims,1)
        ixsmax2=ixOmax2+i*ishift*kr(idims,2);
        ixsRmin1=ixsmin1+ishift*kr(jdims,1)
        ixsRmin2=ixsmin2+ishift*kr(jdims,2)
        ixsRmax1=ixsmax1+ishift*kr(jdims,1)
        ixsRmax2=ixsmax2+ishift*kr(jdims,2);
        ixsLmin1=ixsmin1-ishift*kr(jdims,1)
        ixsLmin2=ixsmin2-ishift*kr(jdims,2)
        ixsLmax1=ixsmax1-ishift*kr(jdims,1)
        ixsLmax2=ixsmax2-ishift*kr(jdims,2);
        wMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1:nwflux)= max(wMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux),&
           w(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2,1:nwflux),w(ixsLmin1:ixsLmax1,&
           ixsLmin2:ixsLmax2,1:nwflux))
        wMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1:nwflux)= min(wMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux),&
           w(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2,1:nwflux),w(ixsLmin1:ixsLmax1,&
           ixsLmin2:ixsLmax2,1:nwflux))
      end do
     
      
    enddo

  end subroutine  extremaw

end module mod_ppm
