      program ELR_strongB
      implicit real*8 (a-h,o-z)
      real*8 n0_cmm3, n0_au, lam_D, kB
      real*8 KP_MIN, KP_MAX, KMAX, KMAX2
      complex*16 CFUN_KP,CRES
      character (24) filename
      character(len=:),allocatable:: str
      external CFUN_KP

      parameter (NVMX=28)

      common/ERROR/ERR_KP,ERR_KZ
      common/plasmaparams/lam_D,vt_au
      common/velocity/u0,xlam

      dimension vp_au(NVMX)

c... Default velocity grid in a.u.
      data vp_au/   0.1D0,   0.2D0,    0.3D0,    0.4D0,    0.6D0,
     |              0.8D0,   1.0D0,    1.2D0,    1.4D0,    1.6D0,
     |              1.8D0,   2.0D0,    2.4D0,    2.8D0,    3.0D0,
     |              4.0D0,   4.5D0,    5.0D0,   6.32D0,   8.94D0,
     |           10.954D0, 12.65D0, 14.142D0, 16.733D0,  20.00D0,
     |            28.28D0, 47.00D0,  63.24D0/

      pi = 3.141592653589793
      ERR_KP = 1E-3
      ERR_KZ = 1E-3

c... other parameters
      Zp = 1.d0                                                         ! charge
      Mp = 1.d0                                                         ! a.m.u.
      xme = 1.d0                                                        ! a.u.

c... Read input data 
      open(10,file='numeric_strong_ELR.inp',status='old')
      READ(10,*) vmin
      READ(10,*) vmax
      read(10,*) tev
      read(10,*) n0_cmm3 !cm^-3
      read(10,*) B0
      close(10)

c... Write output file
      write(*,*) 'Enter the filename ("_eq9_N98.dat" will be added)'
      read(*,*) filename

      str=trim(filename)//'_eq9_N98.dat'
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
      print*,''

c... Loop over u0 (particle velocity || to strong magnetic field B0)
      Nv = 20
      do 100 i = 1, Nv
         xlam = vp_au(i)
c.... add the minimum and maximum velocity given
         if ((xlam.lt.vmin).or.(xlam.gt.vmax)) go to 100
         if ((vp_au(iv+1).gt.vmin).and.(xlam.lt.vmin)) xlam = vmin
         if ((xlam.gt.vmax).and.(vp_au(iv-1).lt.vmax)) xlam = vmax

         u0 = xlam * vt_au
         B = wc_au / wp_au
         KMAX = B / lam_D
         KMAX2 = u0 ** 2 + vt_au ** 2

c.... Numeric ELR for strong magnetic field - Eq. (9)
         cte = 2.d0 * Zp**2 * u0 / pi
         KP_MIN = 0.d0
         KP_MAX = KMAX
         IMAX = 0
         IP = 0
         call INTEG2(KP_MIN,KP_MAX,CFUN_KP,CRES,ERR_KP,IMAX,IP)
         ELR_eq9 = cte * dreal(CRES)
         S_eq9 = ELR_eq9 / u0
         L_eq9 = u0 * ELR_eq9 / (wp_au ** 2)

         open(unit=7,file=str,status='unknown',access='append')
         write(7,*) xlam, ELR_eq9 * conv_factor
         close(7)
         print*, xlam, ELR_eq9 * conv_factor

 100  continue
      close(7)

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
      real*8 KP,KZ_MIN,KZ_MAX,KKP,lam_D,kmax
      complex*16 CFUN_KP,CFUN_KZ,CRES
      external CFUN_KZ

      common/ERROR/ERR_KP,ERR_KZ
      common/KPBCK/KKP
      common/velocity/u0,xlam
      common/plasmaparams/lam_D,vt_au

      KKP = KP
      KZ_MIN = 1E-10
      KZ_MAX = 100
      IMAX = 0
      IP = 0
      call INTEG(KZ_MIN,KZ_MAX,CFUN_KZ,CRES,ERR_KZ,IMAX,IP)
      CFUN_KP = KP * CRES

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
      real*8 kz,kp,k2
      complex*16 CFUN_KZ, cepsilon

      common/KPBCK/kp
      common/velocity/u0,xlam

      k2 = kz**2 + kp**2
      elf = dimag(-1.d0 / cepsilon(kz,kp,xlam))
      CFUN_KZ = dcmplx(kz / k2 * elf, 0.d0)   

      return
      end
c
c-----------------------------------------------------------------------
c
c Dielectric function given by Eq. (10) @ Nersisyan (1998)
c
c-----------------------------------------------------------------------
      function cepsilon(kz,kp,xlam)

      real*8 kz,kp,k2,xlam,vt_au,X,Y,lam_D
      complex*16 caux, cepsilon

      common/plasmaparams/lam_D,vt_au

      k2 = kz**2 + kp**2
      caux = dcmplx(X(xlam), dabs(kz) / kz * Y(xlam))
      cepsilon = (1.d0,0d0) + (1.d0,0d0) / (k2 * lam_D**2) * caux

      return
      end
c
c-----------------------------------------------------------------------
c
c Compute Dawson integral F(x) = exp(−x^2) ∫_0^x exp(t^2)dt for any real x.
c This function was taken from Numerical Recipies in Fortran 77 by W. Press et al.
c
c-----------------------------------------------------------------------
      function dawson(x)
      integer NMAX
      real*8 dawson,x,H,A1,A2,A3,pi,const
      parameter (NMAX=6,H=0.4,A1=2./3.,A2=0.4,A3=2./7.)
      integer i,init,n0
      real*8 d1,d2,e1,e2,sum,x2,xp,xx,c(NMAX)
      save init,c
      data init/0/            ! Flag is 0 if we need to initialize, else 1.

      pi = 3.141592653589793
      const = 1.d0 / dsqrt(pi)
      if (init.eq.0) then
         init = 1
         do i = 1,NMAX
            c(i) = dexp(-((2.d0 * float(i) - 1.d0) * H) ** 2)
         enddo
      endif
      if (dabs(x).lt.0.2) then ! Use series expansion.
         x2 = x ** 2
         dawson = x * (1.d0 - A1 * x2 * (1.d0 - A2 * x2 * (1.d0-A3*x2)))
      else ! Use sampling theorem representation.
         xx = dabs(x)
         n0 = 2 * nint(0.5d0 * xx / H)
         xp = xx - float(n0) * H
         e1 = dexp(2.d0 * xp * H)
         e2 = e1 ** 2
         d1 = float(n0 + 1)
         d2 = d1 - 2.d0
         sum = 0.d0
         do i = 1,NMAX
            sum = sum + c(i) * (e1 / d1 + 1.d0 / (d2 * e1))
            d1 = d1 + 2.d0
            d2 = d2 - 2.d0
            e1 = e2 * e1
         enddo
         dawson = const * sign(dexp(-xp ** 2),x) * sum ! Constant is 1/√π.
      endif

      return
      end
c
c-----------------------------------------------------------------------
c
c Compute the real part of the plasma dispersion function W(z) 
c Ref. [23] @ Nersisyan (1998)
c
c-----------------------------------------------------------------------
      function X(z)
      real*8 X, z, dawson

      X = 1.d0 - dsqrt(2.d0) * z * dawson(z / dsqrt(2.d0))

      return
      end
c
c-----------------------------------------------------------------------
c
c Compute the imaginary part of the plasma dispersion function W(z) 
c Ref. [23] @ Nersisyan (1998)
c
c-----------------------------------------------------------------------
      function Y(z)
      real*8 Y, z

      pi = 3.141592653589793
      Y = dsqrt(pi / 2.d0) * z * dexp(-z*z / 2.d0)

      return
      end
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