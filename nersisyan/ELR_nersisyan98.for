      program ELR_magnetic
      implicit real*8 (a-h,o-z)
      real*8 n0_cmm3, n0_au, lam_D, kB
      real*8 KP_MIN, KP_MAX, KMAX, KMAX2
      complex*16 CFUN_KP,CRES
      character (24) filename
      character(len=:),allocatable:: str
      external CFUN_KP
      
      parameter (NVMX=29)

      common/ERROR/ERR_KP,ERR_KZ,ERR_EPS
      common/plasmaparams/wc_au,omegac_au,vt_au,lam_D
      common/velocity/u,xlam,theta
      common/besselbck/n

      dimension vp_au(NVMX)

c... Default velocity grid in a.u.
      data vp_au/   0.1D0,   0.2D0,    0.3D0,    0.4D0,    0.6D0,
     |              0.8D0,   1.0D0,    1.2D0,    1.4D0,    1.6D0,
     |              1.8D0,   2.0D0,    2.4D0,    2.8D0,    3.0D0,
     |              3.5d0,
     |              4.0D0,   4.5D0,    5.0D0,   6.32D0,   8.94D0,
     |           10.954D0, 12.65D0, 14.142D0, 16.733D0,  20.00D0,
     |            28.28D0, 47.00D0,  63.24D0/

      pi = 3.141592653589793
      ERR_KP = 1E-3
      ERR_KZ = 1E-3
      ERR_EPS = 1E-3

c... other parameters
      Zp = 1.d0                                                         ! charge
      Mp = 1.d0                                                         ! a.m.u.
      xme = 1.d0                                                        ! a.u.

c... Read input data 
      open(10,file='ELR_nersisyan98.inp',status='old')
      READ(10,*) vmin
      READ(10,*) vmax
      read(10,*) tev
      read(10,*) n0_cmm3 !cm^-3
      read(10,*) B0
      read(10,*) theta
      close(10)

c... Write output file
      write(*,*) 'Enter the filename ("_eq8_N98.dat" will be added)'
      read(*,*) filename

      str=trim(filename)//'_eq8_N98.dat'
      open(7,file=str,action='write',status='replace')
      close(7)

c... define constants and conversion factors
      kB = 1.380649E-23                                                 ! J/K
      a0 = 5.29177210903E-11                                            ! m
      a0_cm = 0.529177210903E-8                                         ! cm
      c_au = 137.0                                                      ! a.u.
      conv_amu_to_au = 1.82289E3                                        ! a.u. / amu
      conv_J_to_eV = 6.2414959614E18                                    ! eV / J
      conv_au_to_eV = 27.21132457                                       ! eV / a.u.
      conv_eV_to_au = 1.d0 / conv_au_to_eV                              ! a.u. / eV
      conv_J_to_au = conv_J_to_eV * conv_eV_to_au                       ! a.u. / J
      conv_cmm3_to_au = 1.0E6 * a0**3                                   ! a.u. * cm^3
      conv_Tesla_to_au = 1.d0 / 2.35E5                                  ! a.u. / T
      conv_Gauss_to_Tesla = 1.0E-4                                      ! T / G
      conv_Gauss_to_au = conv_Gauss_to_Tesla * conv_Tesla_to_au         ! a.u. / G
      conv_au_to_second = 2.418884326509E-17                            ! s / a.u.
      conv_factor = conv_au_to_eV / conv_au_to_second * 1E-6            ! for ELR in MeV / s

c... define constants and conversion factors
      T_K = tev * 11600                                                   ! K / eV 
      xMp_au = conv_amu_to_au * Mp                                      ! a.u.
      n0_au = n0_cmm3 * a0_cm ** 3                                      ! a.u.
      vt_au = dsqrt(kB * T_K * conv_J_to_au)                            ! a.u.
      wp_au = dsqrt(4.d0 * pi * n0_au)                                  ! a.u.
      lam_D = vt_au / wp_au                                             ! a.u.
      B0_au = B0 * conv_Gauss_to_au                                     ! a.u. / G
      omegac_au = B0_au / (xMp_au * c_au)
      wc_au = B0_au

      print*, "T_K=", T_K
      print*, "vt_au=", vt_au
      print*, "n0_au=", n0_au
      print*, "wp_au=", wp_au
      print*, "lam_D=", lam_D
      print*, "B0_au=", B0_au
      print*, "wc_au=", wc_au
      print*, "Omegac_au=", omegac_au
      print*,'B0=',B0
      print*,''

c... Loop over u (particle velocity)
      Nv = 20
      do 100 i = 1,Nv
         xlam = vp_au(i)
c.... add the minimum and maximum velocity given
         if ((xlam.lt.vmin).or.(xlam.gt.vmax)) go to 100
         if ((vp_au(iv+1).gt.vmin).and.(xlam.lt.vmin)) xlam = vmin
         if ((xlam.gt.vmax).and.(vp_au(iv-1).lt.vmax)) xlam = vmax

         u = xlam * vt_au
         B = wc_au / wp_au
         KMAX = B / lam_D
         KP_MIN = 0.d0
         KP_MAX = KMAX
         IMAX = 0
         IP = 0

c.... ELR - Eq. (8)
         cte = 2.d0 * Zp**2 / pi
         ELR_eq8 = 0
         do 110 n = 0,0
            call INTEG2(KP_MIN,KP_MAX,CFUN_KP,CRES,ERR_KP,IMAX,IP)
            ELR_eq8_part = cte * dreal(CRES)
            ELR_eq8 = ELR_eq8 + ELR_eq8_part
 110     continue
         open(unit=7,file=str,status='unknown',access='append')
         S_eq8 = ELR_eq8 / u
         write(7,*) xlam, ELR_eq8 * conv_factor, S_eq8
         close(7)
         print*,xlam, ELR_eq8 * conv_factor, S_eq8
 100  continue

      stop
      end
c
c-----------------------------------------------------------------------
c
c Argument of integral over k perpendicular, it goes from 0 to kmax
c
c-----------------------------------------------------------------------
      function CFUN_KP(KP)
      implicit real*8 (a-h,o-z)
      real*8 KP,KZ_MIN,KZ_MAX,KKP,KKZ,lam_D,x,KMAX
      complex*16 CFUN_KP,CFUN_KZ,CRES
      external CFUN_KZ

      common/ERROR/ERR_KP,ERR_KZ,ERR_EPS
      common/velocity/u,xlam,theta
      common/plasmaparams/wc_au,omegac_au,vt_au,lam_D
      common/besselbck/n
      common/KBCK/KKP,KKZ
      common/kmaxbck/KMAX

      KKP = KP
      v = u * dsin(theta)
      ac = v / omegac_au
      x = kp * ac
      KZ_MIN = 1e-6
      KZ_MAX = 20
      IMAX = 0
      IP = 0
      call INTEG(KZ_MIN,KZ_MAX,CFUN_KZ,CRES,ERR_KZ,IMAX,IP)
      CFUN_KP = KP * bessj(n,x) ** 2 * CRES

      return
      end
c
c-----------------------------------------------------------------------
c
c Argument of integral over k_z, it goes from 0 to infinity
c
c-----------------------------------------------------------------------
      function CFUN_KZ(kz)
      implicit real*8 (a-h,o-z)
      real*8 kz,kkz,kp,k2,wc_au
      complex*16 CFUN_KZ, cepsilon, ceps

      common/velocity/u,xlam,theta
      common/plasmaparams/wc_au,omegac_au,vt_au,lam_D
      common/besselbck/n
      common/KBCK/KP,KKZ

      KKZ = kz
      k2 = kz**2 + kp**2
      u0 = u * dcos(theta)
      omega = kz * u0 + n * omegac_au
      ceps = cepsilon(kz,kp,omega)
      elf = dimag(-1.d0 / ceps)
      CFUN_KZ = dcmplx(omega / k2 * elf, 0.d0)

      return
      end
c
c-----------------------------------------------------------------------
c
c Dielectric function given by Eq. (10) @ Nersisyan (1998)
c
c-----------------------------------------------------------------------
      function cepsilon(kz,kp,omega)

      real*8 kz,kp,KKZ,KKP,omg
      real*8 k,k2,p,omega,W,lam_D,vt_au,omegac_au,wc_au
      real*8 tmin, tmax, ERR_KP, ERR_KZ, ERR_EPS
      integer IMAX, IP
      complex*16 caux, cepsilon, carg_epsilon, CRESUL
      external INTEG, carg_epsilon

      common/plasmaparams/wc_au,omegac_au,vt_au,lam_D
      common/ERROR/ERR_KP,ERR_KZ,ERR_EPS
      common/KBCK/KKP,KKZ
      common/omegabck/omg

      KKP = kp
      KKZ = kz
      omg = omega
      k2 = kz**2 + kp**2
      p = omega / (dsqrt(2.d0*k2)*vt_au)
      tmin = 0
      tmax = 100
      IMAX = 0
      IP = 0
      call INTEG(tmin,tmax,carg_epsilon,CRESUL,ERR_EPS,IMAX,IP)
      caux = 1.d0 + dcmplx(0, 2.d0) * p * CRESUL
      cepsilon = 1.d0 + 1.d0 / (k2 * lam_D**2) * caux

      return
      end
c
c-----------------------------------------------------------------------
c
c
c
c-----------------------------------------------------------------------
c
      function carg_epsilon(t)

      real*8 t, u, v, xlam, theta, alpha, omega
      real*8 k2, KP, KZ, p, ac, omegac_au, vt_au, W, wc_au
      real*8 arg0, arg1
      complex*16 carg, carg_epsilon

      common/velocity/u,xlam,theta
      common/plasmaparams/wc_au,omegac_au,vt_au,lam_D
      common/KBCK/KP,KZ
      common/omegabck/omega

      v = u * dsin(theta)
      ac = v / omegac_au
      alpha = datan(kp / kz)
c      alpha = 0.d0
      k2 = kz**2 + kp**2
      p = omega / (dsqrt(2.d0 * k2) * vt_au)
      arg0 = dsqrt(k2) * vt_au
      arg1 = 1.d0 - dcos(dsqrt(2.d0) * wc_au * t / arg0)
c      W = (t * dcos(alpha))**2 + k2 * (ac * dsin(alpha * arg1))**2
      W = (t * dcos(alpha))**2 + k2 * (ac * dsin(alpha))**2 * arg1

      carg = dcmplx(0,2.d0 * p * t) - W
      carg_epsilon = zexp(carg)

      return
      end
c
c-----------------------------------------------------------------------
c
c Bessel function J of general integer order
c
c Returns the Bessel function Jn (x) for any real x and n ≥ 0
c USES bessj0,bessj1
c
c-----------------------------------------------------------------------
      FUNCTION bessj(n,x)
      INTEGER n,IACC
      REAL*8 bessj,x,BIGNO,BIGNI
      PARAMETER (IACC=40,BIGNO=1.e10,BIGNI=1.e-10)
      INTEGER j,jsum,m
      REAL*8 ax,bj,bjm,bjp,sum,tox,bessj0,bessj1

      ! if(n.lt.2) pause "bad argument n in bessj"
      if(n.eq.0) then 
         bessj=bessj0(x)
         return
      else if(n.eq.1) then
         bessj=bessj1(x)
         return
      endif

      ax=dabs(x)
      if(ax.eq.0.) then
         bessj=0.
      else if(ax.gt.float(n)) then            ! Upwards recurrence from J0 and J1 .
         tox=2./ax
         bjm=bessj0(ax)
         bj=bessj1(ax)
         do 11 j=1,n-1
            bjp=j*tox*bj-bjm
            bjm=bj
            bj=bjp
   11    continue
         bessj=bj
      else                                    ! Downwards recurrence from an even m here com
         tox=2./ax                            ! puted. Make IACC larger to increase accuracy.
         m=2*((n+int(dsqrt(float(IACC*n)*1.d0)))/2)
         bessj=0.
         jsum=0                               ! jsum will alternate between 0 and 1; when it is 1, we
         sum=0.                               ! accumulate in sum the even terms in (5.5.16).
         bjp=0.
         bj=1.
         do 12 j=m,1,-1                        ! The downward recurrence.
            bjm=j*tox*bj-bjp
            bjp=bj
            bj=bjm
            if(dabs(bj).gt.BIGNO) then               ! Renormalize to prevent overflows.
               bj=bj*BIGNI
               bjp=bjp*BIGNI
               bessj=bessj*BIGNI
               sum=sum*BIGNI
            endif
            if(jsum.ne.0) sum=sum+bj                 ! Accumulate the sum.
            jsum=1-jsum                             ! Change 0 to 1 or vice versa.
            if(j.eq.n) bessj=bjp                     ! Save the unnormalized answer.
   12    continue
         sum=2.*sum-bj                           ! Compute (5.5.16)
         bessj=bessj/sum                         ! and use it to normalize the answer.
      endif
      if(x.lt.0..and.mod(n,2).eq.1) bessj=-bessj

      return
      END
c
c-----------------------------------------------------------------------
c
c Bessel function J0
c
c Returns the Bessel function J0 (x) for any real x.
c
c-----------------------------------------------------------------------
      FUNCTION bessj0(x)
      REAL*8 bessj0,x
      REAL*8 ax,xx,z
      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,
     |  r5,r6,s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     |  s1,s2,s3,s4,s5,s6
      DATA p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,
     |   -.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/
     |   -.1562499995d-1,.1430488765d-3,-.6911147651d-5,.7621095161d-6,
     |   -.934945152d-7/
      DATA r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,
     |   651619640.7d0,-11214424.18d0,77392.33017d0,-184.9052456d0/,
     |   s1,s2,s3,s4,s5,s6/57568490411.d0,1029532985.d0,
     |   9494680.718d0,59272.64853d0,267.8532712d0,1.d0/

      if(dabs(x).lt.8.) then ! Direct rational function fit.
        y=x**2
        bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))
     |         /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
      else ! Fitting function (6.5.9).
        ax=dabs(x)
        z=8./ax
        y=z**2
        xx=ax-.785398164
        bessj0=dsqrt(.636619772/ax)*(dcos(xx)*(p1+y*(p2+y*(p3+y*(p4+y
     |         *p5))))-z*dsin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif

      return
      END
c
c-----------------------------------------------------------------------
c
c Bessel function J1
c
c Returns the Bessel function J1 (x) for any real x.
c
c-----------------------------------------------------------------------
      FUNCTION bessj1(x)
      REAL*8 bessj1,x
      REAL*8 ax,xx,z
      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,
     |   r5,r6,s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     |   s1,s2,s3,s4,s5,s6
      DATA r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0,
     |   242396853.1d0,-2972611.439d0,15704.48260d0,-30.16036606d0/,
     |   s1,s2,s3,s4,s5,s6/144725228442.d0,2300535178.d0,
     |   18583304.74d0,99447.43394d0,376.9991397d0,1.d0/
      DATA p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,
     |   .2457520174d-5,-.240337019d-6/, q1,q2,q3,q4,q5/.04687499995d0,
     |   -.2002690873d-3,.8449199096d-5,-.88228987d-6,.105787412d-6/

      if(dabs(x).lt.8.) then ! Direct rational approximation.
         y=x**2
         bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))
     |          /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
      else ! Fitting function (6.5.9).
         ax=dabs(x)
         z=8./ax
         y=z**2
         xx=ax-2.356194491
         bessj1=dsqrt(.636619772/ax)*(dcos(xx)*(p1+y*(p2+y*(p3+y*(p4+y
     |          *p5))))-z*dsin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
     |          *dsign(1.d0,x)
      endif

      return
      END
c
c-----------------------------------------------------------------------
c
C  Integration routine INTEG
c
C  CRESUL = INT. ON X [XL,+INFINITY] CFCT (X) dX
c
C  IT IS SUSPECTED THE MAIN CONTRIBUTION COMES FROM THE REGION
C  CFCT(XL) TO CFCT(XL+DX), AFTERWARDS MOVES FOR ADDING RANGES
C  OF DX
c
C  IT USES THE SUBROUTINE IDELZ
c
C      XL : LOWER LIMIT
C      ER : RELATIVE ERROR, SUGGESTED VALUES 1.E-3,1.E-2,
C    IMAX : NUMBER OF STEPS ('MOVES FOR '),DEFAULT:IF IMAX=0, IMAX=30
C    CFCT : COMPLEX FUNCTION TO INTEGRATE
C  CRESUL : RESULT
C      IP : NUMERO DE SUBDIVISIONES(DEBE SER<25),DEFAULT:IF IP=0, IP=24 
c
C NECESITA TENER ABIERTO UN ARCHIVO CON UNIT=3 PARA ESCRIBIR EN
C CASO DE PROBLEMAS
c
c-----------------------------------------------------------------------
      SUBROUTINE INTEG(XL,DX,CFCT,CRESUL,ER,IMAX,IP)
      
      IMPLICIT REAL*8 (A-B,D-H,O-Z)
      IMPLICIT DOUBLE COMPLEX (C)
      parameter(NIMX=60)
      
      EXTERNAL CFCT
      
      DATA C0/(0.D0,0.D0)/
      
      IF(IMAX.EQ.0)IIMAX=NIMX
      
    ! ------------------------------------------------------------------
    !   Primera integracion : [XL,XL+DX]
    ! ------------------------------------------------------------------
      XU=XL+DX
      CALL IDELZ(XL,XU,CFCT,CR0,ER,NDIM0,IP)
      CR=CR0
      R=ABS(CR)

    ! ------------------------------------------------------------------
    !   Inicializacion de las integraciones siguientes
    ! ------------------------------------------------------------------
      XUI=XU
      ERI=ER
      
    ! ------------------------------------------------------------------
    !  Sucesivos pasos de integracion
    ! ------------------------------------------------------------------
      DO 1 I=1,IIMAX
         CRI=C0
         XLI=XUI
         XUI=XLI+DX
         CALL IDELZ(XLI,XUI,CFCT,CRI,ERI,NDIMI,IP)
         CR=CR+CRI
         RI=ABS(CRI)
         R=ABS(CR)
         RATIO=RI/R
         IF(RATIO.LT.ER) GO TO 4
         IF(RATIO.LT.0.1) ERI=ER/RATIO
         IF(ERI.GT.0.1) ERI=0.1D0
    1 CONTINUE
      
      WRITE(3,200) RATIO
  200 FORMAT(1X,'IDELZ, I>IMAX, RATIO=',1PE10.3,/,80('-'))

    4 CONTINUE
      CRESUL=CR
      
      RETURN
      END
c
c-----------------------------------------------------------------------
c
c Adaptative integration routine IDELZ
c
c Integrates a function that is complex in a certain region by introducing 
c more points in a given interval.
c
c       XL,XU : lower and upper integration limits
c       CF    : function
c               (must be declared EXTERNAL in the calling routine)
c       CRES  : result
c       ERR   : integration error
c       N     : number of points 3+2*N
c       KMAX  : maximum number of partitions in a given interval
c               (must be less than 25). If KMAX=0, then KMAX=24
c
c Uses an 3-point Simpson rule integration routine.
c
c-----------------------------------------------------------------------
      SUBROUTINE IDELZ(XL,XU,CF,CRES,ERR,N,KMAX)	

      IMPLICIT REAL*8 (A-B,D-H,O-Z)
      IMPLICIT DOUBLE COMPLEX(C)
      parameter(NPMX=55)

      DIMENSION XXL(NPMX+5),XXU(NPMX+5),CINT(NPMX+5)
      DIMENSION CFXL(NPMX+5),CFXM(NPMX+5),CFXU(NPMX+5)

      DATA C0/(0.D0,0.D0)/

      CRES=C0

      IF(KMAX.EQ.0) KKMAX=NPMX

      IF(XL.GT.XU)GO TO 1000

    ! ------------------------------------------------------------------
    ! VALORES INICIALES
    ! ------------------------------------------------------------------

      N=0
      K=1

      XXL(1)=XL
      XXU(1)=XU
      XM=0.5D0*(XU+XL)

      CFXL(1)=CF(XL)
      CFXU(1)=CF(XU)
      CFXM(1)=CF(XM)

      CALL SIMP(XXL(1),XXU(1),XM,CFXL(1),CFXU(1),CFXM(1),CINT(1))

    ! ------------------------------------------------------------------
    ! COMIENZO DE LA ITERACION
    ! ------------------------------------------------------------------

 100  IF(K+1.GT.KKMAX) GO TO 900
 200  N=N+1

      CR=CINT(K)

      XXU(K+1)=XXU(K)
      CFXU(K+1)=CFXU(K)

      XXU(K)=0.5*(XXL(K)+XXU(K))
      CFXU(K)=CFXM(K)

      XXL(K+1)=XXU(K)
      CFXL(K+1)=CFXM(K)

      XM=0.5D0*( XXU(K)+XXL(K) )
      CFXM(K)=CF(XM)

      CALL SIMP(XXL(K),XXU(K),XM,CFXL(K),CFXU(K),CFXM(K),CINT(K))

      K=K+1

      XM=0.5D0*( XXU(K)+XXL(K) )
      CFXM(K)=CF(XM)

      CALL SIMP(XXL(K),XXU(K),XM,CFXL(K),CFXU(K),CFXM(K),CINT(K))

      CRR=CINT(K-1)+CINT(K)

      IF(ABS(CR-CRR).GT.(ERR*ABS(CRES+CRR)))GO TO 100        

    ! ------------------------------------------------------------------
    ! FIN DE LA ITERACION EN K
    ! ------------------------------------------------------------------

      CRES=CRES+CRR

      K=K-2
      IF(K.GT.0) GO TO 200

    ! ------------------------------------------------------------------
    ! FIN DE LA INTEGRACION
    ! ------------------------------------------------------------------

      RETURN

    ! ------------------------------------------------------------------
    ! MENSAJES DE ERROR
    ! ------------------------------------------------------------------

 900  ERROR=ABS(CR-CRR)/ABS(CRES+CRR)
      WRITE(3,10)ERROR,XXL(K),XXU(K)
      STOP 'ERROR IN INTEG1'

 1000 WRITE(3,11)

 10   FORMAT(1X,'INTEG1',/,1X,'ERROR=',1PE10.3,5X,       
     |       'XXL=',1PE10.3,5X,'XXU=',1PE10.3,/,80('-'))    
 11   FORMAT(1X,'INTEG1',/,1X,'ERROR : XL.GT.XU',/,80('-'))    
      STOP 'ERROR IN INTEG1'

      END
c
c-----------------------------------------------------------------------
c
c 3-point Simpson integration routine
c
c-----------------------------------------------------------------------
      SUBROUTINE SIMP(XL,XU,XM,CF0,CF2,CF1,CY)

      IMPLICIT REAL*8 (A-B,D-H,O-Z)
      IMPLICIT DOUBLE COMPLEX(C)

      W=0.5D0*(XU-XL)

      CY=(W/3.D0)*(CF0+4.D0*CF1+CF2)

      RETURN
      END

C **********************************************************************
C    SUBRUTINA DE INTEGRACION

C       Integra una funcion que es complicada solo en cierta region,
C       poniendo mas puntos en dicha zona.

C       XL,XU : limites inferior y superior de integracion.        
C       CF : funcion compleja a integrar
C            (Debe ser declarada EXTERNAL en el programa de llamamiento)
C       CRES : resultado
C       ERR : error de integracion pedido
C       N : el numero de puntos usados es 3+2*N
C       KMAX :numero maximo de particiones de un dado intervalo
C             Debe ser menor que 25. Si KMAX=0 toma KMAX=24

C       Usa una subrutina de integracion por Simpson con 3 puntos

C ----------------------------------------------------------------------
      SUBROUTINE INTEG2(XL,XU,CF,CRES,ERR,N,KMAX)	
      
      IMPLICIT REAL*8 (A-B,D-H,O-Z)
      IMPLICIT DOUBLE COMPLEX(C)
      parameter(NPMX=55)

      DIMENSION XXL(NPMX+5),XXU(NPMX+5),CINT(NPMX+5)
      DIMENSION CFXL(NPMX+5),CFXM(NPMX+5),CFXU(NPMX+5)

      DATA C0/(0.D0,0.D0)/

      CRES=C0

      IF(XL.GT.XU)GO TO 1000

      IF(KMAX.EQ.0)KKMAX=NPMX


C -----VALORES INICIALES

      N=0
      K=1

      XXL(1)=XL
      XXU(1)=XU
      XM=0.5D0*(XU+XL)

      CFXL(1)=CF(XL)
      CFXU(1)=CF(XU)
      CFXM(1)=CF(XM)

      CALL SIMP(XXL(1),XXU(1),XM,CFXL(1),CFXU(1),CFXM(1),CINT(1))

C -----COMIENZO DE LA ITERACION

 100  IF(K+1.GT.KKMAX) GO TO 900
 200  N=N+1

      CR=CINT(K)

      XXU(K+1)=XXU(K)
      CFXU(K+1)=CFXU(K)

      XXU(K)=0.5*(XXL(K)+XXU(K))
      CFXU(K)=CFXM(K)

      XXL(K+1)=XXU(K)
      CFXL(K+1)=CFXM(K)

      XM=0.5D0*( XXU(K)+XXL(K) )
      CFXM(K)=CF(XM)

      CALL SIMP(XXL(K),XXU(K),XM,CFXL(K),CFXU(K),CFXM(K),CINT(K))

      K=K+1

      XM=0.5D0*( XXU(K)+XXL(K) )
      CFXM(K)=CF(XM)

      CALL SIMP(XXL(K),XXU(K),XM,CFXL(K),CFXU(K),CFXM(K),CINT(K))

      CRR=CINT(K-1)+CINT(K)

      IF(ABS(CR-CRR).GT.(ERR*ABS(CRES+CRR)))GO TO 100        

C ------ FIN DE LA ITERACION EN K

      CRES=CRES+CRR

      K=K-2
      IF(K.GT.0) GO TO 200

C ------FIN DE LA INTEGRACION

      RETURN

C ----------------------------------------------------------------------
C ------MENSAJES DE ERROR
C ----------------------------------------------------------------------

 900  ERROR=ABS(CR-CRR)/ABS(CRES+CRR)
      WRITE(3,10)ERROR,XXL(K),XXU(K)
 10   FORMAT(1X,'INTEG3',/,1X,'ERROR=',1PE10.3,5X,       
     1       'XXL=',1PE10.3,5X,'XXU=',1PE10.3,/,80('-'))    
      STOP 'ERROR: INTEG3'

 1000 WRITE(3,11)
 11   FORMAT(1X,'INTEG3',/,1X,'ERROR : XL.GT.XU',/,80('-'))    
      STOP 'ERROR: INTEG3'

      END