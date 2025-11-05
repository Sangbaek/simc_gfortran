      program call_maid
c
      real*8 siglab,sigcm
      real*8 EIMEV,EFMEV,THE,omega
      logical prod_in_cm
      common /cm_logical/  prod_in_cm
c
      prod_in_cm = .TRUE.
      EIMEV=4560.
      Q2 = 0.33e+6
      W = 1232.
      tarmass=938.27
      omega = (Q2 + W*W - tarmass*tarmass) / (2.*tarmass)
      EFMEV=EIMEV-omega
      pi = 3.14159
      THE = 180/pi*acos(1-Q2/2/EIMEV/EFMEV)
         call  get_xn_maid_07(Q2,W,EIMEV,EFMEV
     > ,0.436332,3.14159,siglab,sigcm)
         write(*,*) q2,w,eimev,efmev,the,siglab,sigcm
         
c
      end
