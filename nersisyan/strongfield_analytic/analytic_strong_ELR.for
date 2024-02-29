      program ELR_strongB
      implicit real*8 (a-h,o-z)
      real*8 n0_cmm3, n0_au, lam_D, kB
      parameter (NVMX=28)

      common/plasmaparams/lam_D,vt_au
      common/velocity/u0

      dimension vp_au(NVMX)

c... Default velocity grid in a.u.
      data vp_au/   0.1D0,   0.2D0,    0.3D0,    0.4D0,    0.6D0,
     |              0.8D0,   1.0D0,    1.2D0,    1.4D0,    1.6D0,
     |              1.8D0,   2.0D0,    2.4D0,    2.8D0,    3.0D0,
     |              4.0D0,   4.5D0,    5.0D0,   6.32D0,   8.94D0,
     |           10.954D0, 12.65D0, 14.142D0, 16.733D0,  20.00D0,
     |            28.28D0, 47.00D0,  63.24D0/

      pi = 3.141592653589793

c... projectile
      Zp = 1.d0                                                         ! charge
      Mp = 1.d0                                                         ! a.m.u.

c... plasma electrons
      n0_cmm3 = 1.0E14                                                  ! cm^-3
      T = 10                                                            ! eV
      xme = 1.d0                                                        ! a.u.

c... magnetic field
c      B0 = 50E3                                                         ! G
      B0 = 80E3                                                         ! G

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
      T_K = T * 11600                                                   ! K / eV 
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
      print*,''


c... Loop over u0 (particle velocity || to strong magnetic field B0)
      Nv = 20
      do 100 i = 1, Nv
         xlam = vp_au(i)
         write(222, *) xlam, dawson(xlam)

         u0 = xlam * vt_au
         B = wc_au / wp_au

c.... ELR for large velocities - Eq. (19)
         S_eq19 = Zp**2 * wp_au ** 2 / (2.d0 * vt_au * xlam)
         S_eq19 = S_eq19 * conv_factor

c.... ELR for small velocities - Eq. (18)
         cte = Zp**2 * vt_au / (2.d0 * dsqrt(2*pi) * lam_D**2)
         fun1 = 2.d0 * xlam**2 * dlog(B)
         fun2 = xlam**4 * (dlog(B) - 1.d0 + pi/12.d0)
         G1 = xlam ** 2 - 3.d0/2.d0 * xlam ** 4
         G2 = xlam ** 4 * (pi / 12.d0 + 5.d0/4.d0) - 0.5d0 * xlam ** 2
         fun3 = G1 / B ** 2 + G2 / B ** 4
         S_eq18 = cte * (fun1 - fun2 + fun3) * conv_factor

c.... ELR for strong magnetic field - Eq. (14)
         cte2 = Zp**2 * vt_au / (2.d0 * pi * lam_D ** 2)
         YY = Y(xlam)
         XX = X(xlam)
         xn1 = YY**2 + (B**2 + XX)**2
         xd1 = YY**2 + XX**2
         F1 = YY / 2.d0 * dlog(xn1/xd1)
         atan_arg = (B**2 + XX)/YY
         F2 = B**2 * (pi/2.d0 - atan(atan_arg))
         F3 = atan(XX/YY) - atan(atan_arg)
         S_eq14 = cte2 * xlam * (F1 + F2 + XX * F3)
         S_eq14 = S_eq14 * conv_factor

         write(50,*) xlam, S_eq19, S_eq18, S_eq14
c         write(80,*) xlam, S_eq19, S_eq18, S_eq14

 100  continue

      stop
      end
c-----------------------------------------------------------------------
c
c Compute Dawson integral F(x) = exp(−x^2) ∫_0^x exp(t^2)dt for any real x.
c This function was taken from Numerical Recipies in Fortran 77 by W. Press et al.
c
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
c-----------------------------------------------------------------------
c
c Compute the real part of the plasma dispersion function W(z) Ref. [23] @ Nersisyan (1998)
c
      function X(z)
      real*8 X, z, dawson

      X = 1.d0 - dsqrt(2.d0) * z * dawson(z / dsqrt(2.d0))

      return
      end
c-----------------------------------------------------------------------
c
c Compute the imaginary part of the plasma dispersion function W(z) Ref. [23] @ Nersisyan (1998)
c
      function Y(z)
      real*8 Y, z

      pi = 3.141592653589793
      Y = dsqrt(pi / 2.d0) * z * dexp(-z*z / 2.d0)

      return
      end
