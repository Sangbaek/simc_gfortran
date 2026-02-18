	subroutine enerloss_new(len,dens,zeff,aeff,epart,mpart,typeflag,Eloss)

	implicit none

	real*8 thick,len,dens,zeff,aeff,epart,mpart,Eloss
	real*8 x,chsi,lambda,gauss1,Eloss_mp,gamma,Eloss_mp_pdg
	real*8 denscorr,CO,hnup,log10bg,I,beta,Eloss_mp_new
	real*8 ranlan1,GVAVIV
	real*8 FLAND
	real*8 grnd
	real*8 emax,kappa
	real*8 xmean,xlamx,xlamb,beta2
	real*8 IX(92)
	data IX /19.2, 41.8, 40.0, 63.7, 76.0, 78.0, 82.0 , 95.0,  115.0, 137.0,
     A 149.0, 156.0, 166.0, 173.0, 173.0, 180.0,  174.0, 188.0, 190.0, 191.0,
     B 216.0, 233.0, 245.0, 257.0,  272.0, 286.0, 297.0, 311.0, 322.0, 330.0,
     C 334.0, 350.0,  347.0, 348.0, 343.0, 352.0, 363.0, 366.0, 379.0, 393.0,
     D 417.0, 424.0, 428.0, 441.0, 449.0, 470.0, 470.0, 469.0,  488.0, 488.0, 
     E 487.0, 485.0, 491.0, 482.0, 488.0, 491.0,  501.0, 523.0, 535.0, 546.0, 
     F 560.0, 574.0, 580.0, 591.0,  614.0, 628.0, 650.0, 658.0, 674.0, 684.0, 
     G 694.0, 705.0,  718.0, 727.0, 736.0, 746.0, 757.0, 790.0, 790.0, 800.0,
     H 810.0, 823.0, 823.0, 830.0, 825.0, 794.0, 827.0, 826.0,  841.0, 847.0, 
     I 878.0, 890.0 /

	integer typeflag          !1=normal eloss (picked from distribution)
                                  !2=min eloss
	                          !3=max eloss
                                  !4=most probable eloss
        integer numerr
        data numerr /0/
	logical use_old /.FALSE./

	real*8 me
        parameter(me=0.51099906)
	
		integer ZI
		ZI = NINT(zeff)

	thick = len*dens
	gamma=epart/mpart
        beta = sqrt(1.-1./gamma**2)

C	if(zeff.eq.1) then	!Ionization potential in MeV
C	  I = 21.8e-06
C	else
C	  I = (16.*zeff**0.9)*1.0e-06
C	endif
	
	chsi = 0.307075/2.*zeff/aeff*thick/beta**2
	emax = 2*me*beta*beta*gamma*gamma/(1+2*gamma*me/mpart+(me/mpart)**2)
	kappa = chsi/emax

	hnup = 28.816e-06*sqrt(dens*zeff/aeff) !plasma frequency
	log10bg = log(beta*gamma)/log(10.)
	CO=log(hnup)-log(IX(ZI)*1e-6)+0.5

C DJG Get density effect correction (I got this from JV).

	if(log10bg.lt.0.) then
	  denscorr=0.
	elseif(log10bg.lt.3.) then
	  denscorr=CO+log(10.)*log10bg+abs(CO/27.)*(3.-log10bg)**3
	elseif(log10bg.lt.4.7) then
	  denscorr=CO+log(10.)*log10bg
	else
	  denscorr=CO+log(10.)*4.7
	endif
	   
	if (thick.le.0.) then
	  Eloss = 0.
	else
	  Eloss_mp = 0.1536e-03 * zeff/aeff * thick * ( 19.26 +
     &          log(thick/dens) )
	  Eloss_mp_new = 0.1536e-03 * zeff/aeff *thick/beta**2* (
     &          log(me/I**2) + 1.063 + 2.*log(gamma*beta) + 
     &		log(0.1536*zeff/aeff*thick/beta**2)-beta**2-denscorr)
	  Eloss_mp_pdg = chsi*(log(2.0*me*beta**2*gamma**2/(IX(ZI)*1e-6))
     & + log(chsi/(IX(ZI)*1e-6)) + 0.2 -beta**2 - denscorr)
c	  write(6,*) 'ELOSS',Eloss_mp,Eloss_mp_new 
! ........ convert to MeV, the unit of choice in THIS program
! ........ (cf. EVCOIN where GeV prevail)
C	  Eloss_mp = Eloss_mp_new*1000.
C	 Eloss = Eloss_mp_pdg

         if (use_old) then	  
	  x = grnd()
	  if(typeflag.eq.1)then
	    x=ranlan1(x)
	  elseif(typeflag.eq.2)then
	    x= -eloss_mp_pdg/chsi
	  elseif(typeflag.eq.3)then
	    x=25.
	  elseif(typeflag.eq.4)then
	    x= 0.
          endif
c
	  else  ! use new landau and vavilov
	  if (kappa < 0.01) then   
	     xmean = -beta**2-log(chsi)-0.422784
	     xlamx = FLAND(XMEAN)
 25          CALL GLANDG(XLAMB)
             IF(XLAMB.GT.XLAMX) GO TO 25
	  else
        	  x = grnd()
		  beta2= beta**2
		  xlamb = GVAVIV(kappa,beta2,x)
	  endif
	  if(typeflag.eq.1)then
	    x=xlamb
	  elseif(typeflag.eq.2)then
	    x= -eloss_mp_pdg/chsi
	  elseif(typeflag.eq.3)then
	    x=25.
	  elseif(typeflag.eq.4)then
	    x= 0.
          endif
c
	  endif
	  Eloss = x*chsi+eloss_mp_pdg
	endif

        if (eloss.gt.(epart-mpart)) then
	   numerr=numerr+1
c	   if (numerr.le.10) then
	      write(6,'(a)') 'Eloss>Total KE; forcing Eloss=KE'
c	      if (numerr.eq.10) write(6,*) '     FURTHER ELOSS ERRORS SUPPRESSED'
	      write(*,'(i8,a,i3,a,f10.5,a,f10.5)')  
     > numerr," typeflag = ",typeflag," part mass= "
     >,mpart," set eloss to Part kinematic energy = " ,epart-mpart
c	   endif
	   eloss=(epart-mpart)-0.0000001
        endif 

	return
	end
c
      FUNCTION FLAND(X)
	real*8 p1,p2,p3,p4,p5,p6
	real*8 X,fland
      PARAMETER (P1=.60715,P2=.67794,P3=.52382E-1,P4=.94753,
     +           P5=.74442,P6=1.1934)
C     FLAND(X) = P1+P6*X+(P2+P3*X)*EXP(P4*X+P5)
      FLAND = P1+P6*X+(P2+P3*X)*EXP(P4*X+P5)
C      RETURN
      END FUNCTION FLAND
c
      SUBROUTINE GLANDG(YRAN)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Copy of the CERN library routine GENLAN                       *
C.    *  Generation of LANDAU-distributed random numbers by 4-point    *
C.    *  interpolation in the previously-tabulated inverse cumulative  *
C.    *  distribution                                                  *
C.    *                                                                *
C.    *    ==>Called by : GLANDO                                       *
C.    *                                                                *
C.    ******************************************************************
C.
      implicit none
      real*8 RANGE1,RANGE2
      real*8 X1BOT,X2BOT,GP1INV,GP2INV,XCUM1,XCUM2
	real*8 X,YRAN,TABPO1,TABPO2,p,a,b
	real*8 grnd
	integer j
c
      PARAMETER (RANGE1=0.807069,RANGE2=0.994869)
      PARAMETER (X1BOT =0.0     ,X2BOT =0.791240)
      PARAMETER (GP1INV=118.836 ,GP2INV=486.178 )
C
      DIMENSION XCUM1(100) , XCUM2(100)
C
C====>  1ST TABLE OF INVERSE CUMULATIVE LANDAU POINTS. ( 0.0<P<0.833 )
      DATA XCUM1(  1) / -2.57000 /   ,   XCUM1(  2) / -2.15412 /
      DATA XCUM1(  3) / -1.94167 /   ,   XCUM1(  4) / -1.79583 /
      DATA XCUM1(  5) / -1.67975 /   ,   XCUM1(  6) / -1.58090 /
      DATA XCUM1(  7) / -1.49341 /   ,   XCUM1(  8) / -1.41397 /
      DATA XCUM1(  9) / -1.34057 /   ,   XCUM1( 10) / -1.27185 /
      DATA XCUM1( 11) / -1.20686 /   ,   XCUM1( 12) / -1.14490 /
      DATA XCUM1( 13) / -1.08545 /   ,   XCUM1( 14) / -1.02809 /
      DATA XCUM1( 15) /  -.97249 /   ,   XCUM1( 16) /  -.91839 /
      DATA XCUM1( 17) /  -.86557 /   ,   XCUM1( 18) /  -.81386 /
      DATA XCUM1( 19) /  -.76308 /   ,   XCUM1( 20) /  -.71311 /
      DATA XCUM1( 21) /  -.66384 /   ,   XCUM1( 22) /  -.61515 /
      DATA XCUM1( 23) /  -.56696 /   ,   XCUM1( 24) /  -.51919 /
      DATA XCUM1( 25) /  -.47176 /   ,   XCUM1( 26) /  -.42462 /
      DATA XCUM1( 27) /  -.37769 /   ,   XCUM1( 28) /  -.33092 /
      DATA XCUM1( 29) /  -.28426 /   ,   XCUM1( 30) /  -.23767 /
      DATA XCUM1( 31) /  -.19109 /   ,   XCUM1( 32) /  -.14447 /
      DATA XCUM1( 33) /  -.09779 /   ,   XCUM1( 34) /  -.05100 /
      DATA XCUM1( 35) /  -.00405 /   ,   XCUM1( 36) /   .04308 /
      DATA XCUM1( 37) /   .09044 /   ,   XCUM1( 38) /   .13806 /
      DATA XCUM1( 39) /   .18598 /   ,   XCUM1( 40) /   .23423 /
      DATA XCUM1( 41) /   .28285 /   ,   XCUM1( 42) /   .33187 /
      DATA XCUM1( 43) /   .38134 /   ,   XCUM1( 44) /   .43128 /
      DATA XCUM1( 45) /   .48174 /   ,   XCUM1( 46) /   .53275 /
      DATA XCUM1( 47) /   .58435 /   ,   XCUM1( 48) /   .63658 /
      DATA XCUM1( 49) /   .68948 /   ,   XCUM1( 50) /   .74309 /
      DATA XCUM1( 51) /   .79746 /   ,   XCUM1( 52) /   .85263 /
      DATA XCUM1( 53) /   .90865 /   ,   XCUM1( 54) /   .96556 /
      DATA XCUM1( 55) /  1.02342 /   ,   XCUM1( 56) /  1.08228 /
      DATA XCUM1( 57) /  1.14219 /   ,   XCUM1( 58) /  1.20322 /
      DATA XCUM1( 59) /  1.26542 /   ,   XCUM1( 60) /  1.32887 /
      DATA XCUM1( 61) /  1.39362 /   ,   XCUM1( 62) /  1.45976 /
      DATA XCUM1( 63) /  1.52736 /   ,   XCUM1( 64) /  1.59650 /
      DATA XCUM1( 65) /  1.66727 /   ,   XCUM1( 66) /  1.73976 /
      DATA XCUM1( 67) /  1.81407 /   ,   XCUM1( 68) /  1.89032 /
      DATA XCUM1( 69) /  1.96860 /   ,   XCUM1( 70) /  2.04905 /
      DATA XCUM1( 71) /  2.13180 /   ,   XCUM1( 72) /  2.21699 /
      DATA XCUM1( 73) /  2.30477 /   ,   XCUM1( 74) /  2.39531 /
      DATA XCUM1( 75) /  2.48878 /   ,   XCUM1( 76) /  2.58540 /
      DATA XCUM1( 77) /  2.68535 /   ,   XCUM1( 78) /  2.78889 /
      DATA XCUM1( 79) /  2.89626 /   ,   XCUM1( 80) /  3.00775 /
      DATA XCUM1( 81) /  3.12364 /   ,   XCUM1( 82) /  3.24429 /
      DATA XCUM1( 83) /  3.37005 /   ,   XCUM1( 84) /  3.50136 /
      DATA XCUM1( 85) /  3.63866 /   ,   XCUM1( 86) /  3.78246 /
      DATA XCUM1( 87) /  3.93334 /   ,   XCUM1( 88) /  4.09194 /
      DATA XCUM1( 89) /  4.25900 /   ,   XCUM1( 90) /  4.43533 /
      DATA XCUM1( 91) /  4.62186 /   ,   XCUM1( 92) /  4.81960 /
      DATA XCUM1( 93) /  5.02974 /   ,   XCUM1( 94) /  5.25368 /
      DATA XCUM1( 95) /  5.49312 /   ,   XCUM1( 96) /  5.74987 /
      DATA XCUM1( 97) /  6.02605 /   ,   XCUM1( 98) /  6.32428 /
      DATA XCUM1( 99) /  6.64773 /   ,   XCUM1(100) /  7.00000 /
C
C====>  2ND TABLE OF INVERSE CUMULATIVE LANDAU POINTS. ( 0.791<P<0.995 )
      DATA XCUM2(  1) /   5.50000 /   ,   XCUM2(  2) /   5.56120 /
      DATA XCUM2(  3) /   5.62347 /   ,   XCUM2(  4) /   5.68684 /
      DATA XCUM2(  5) /   5.75133 /   ,   XCUM2(  6) /   5.81699 /
      DATA XCUM2(  7) /   5.88383 /   ,   XCUM2(  8) /   5.95191 /
      DATA XCUM2(  9) /   6.02125 /   ,   XCUM2( 10) /   6.09190 /
      DATA XCUM2( 11) /   6.16391 /   ,   XCUM2( 12) /   6.23732 /
      DATA XCUM2( 13) /   6.31219 /   ,   XCUM2( 14) /   6.38855 /
      DATA XCUM2( 15) /   6.46646 /   ,   XCUM2( 16) /   6.54597 /
      DATA XCUM2( 17) /   6.62712 /   ,   XCUM2( 18) /   6.70998 /
      DATA XCUM2( 19) /   6.79460 /   ,   XCUM2( 20) /   6.88103 /
      DATA XCUM2( 21) /   6.96935 /   ,   XCUM2( 22) /   7.05962 /
      DATA XCUM2( 23) /   7.15192 /   ,   XCUM2( 24) /   7.24633 /
      DATA XCUM2( 25) /   7.34294 /   ,   XCUM2( 26) /   7.44182 /
      DATA XCUM2( 27) /   7.54306 /   ,   XCUM2( 28) /   7.64676 /
      DATA XCUM2( 29) /   7.75300 /   ,   XCUM2( 30) /   7.86188 /
      DATA XCUM2( 31) /   7.97351 /   ,   XCUM2( 32) /   8.08800 /
      DATA XCUM2( 33) /   8.20548 /   ,   XCUM2( 34) /   8.32610 /
      DATA XCUM2( 35) /   8.44997 /   ,   XCUM2( 36) /   8.57725 /
      DATA XCUM2( 37) /   8.70808 /   ,   XCUM2( 38) /   8.84262 /
      DATA XCUM2( 39) /   8.98103 /   ,   XCUM2( 40) /   9.12349 /
      DATA XCUM2( 41) /   9.27021 /   ,   XCUM2( 42) /   9.42141 /
      DATA XCUM2( 43) /   9.57730 /   ,   XCUM2( 44) /   9.73812 /
      DATA XCUM2( 45) /   9.90410 /   ,   XCUM2( 46) /  10.07552 /
      DATA XCUM2( 47) /  10.25265 /   ,   XCUM2( 48) /  10.43584 /
      DATA XCUM2( 49) /  10.62540 /   ,   XCUM2( 50) /  10.82169 /
      DATA XCUM2( 51) /  11.02508 /   ,   XCUM2( 52) /  11.23598 /
      DATA XCUM2( 53) /  11.45487 /   ,   XCUM2( 54) /  11.68222 /
      DATA XCUM2( 55) /  11.91855 /   ,   XCUM2( 56) /  12.16441 /
      DATA XCUM2( 57) /  12.42045 /   ,   XCUM2( 58) /  12.68734 /
      DATA XCUM2( 59) /  12.96580 /   ,   XCUM2( 60) /  13.25663 /
      DATA XCUM2( 61) /  13.56075 /   ,   XCUM2( 62) /  13.87912 /
      DATA XCUM2( 63) /  14.21280 /   ,   XCUM2( 64) /  14.56296 /
      DATA XCUM2( 65) /  14.93088 /   ,   XCUM2( 66) /  15.31799 /
      DATA XCUM2( 67) /  15.72593 /   ,   XCUM2( 68) /  16.15650 /
      DATA XCUM2( 69) /  16.61170 /   ,   XCUM2( 70) /  17.09378 /
      DATA XCUM2( 71) /  17.60523 /   ,   XCUM2( 72) /  18.14884 /
      DATA XCUM2( 73) /  18.72792 /   ,   XCUM2( 74) /  19.34622 /
      DATA XCUM2( 75) /  20.00796 /   ,   XCUM2( 76) /  20.71792 /
      DATA XCUM2( 77) /  21.48176 /   ,   XCUM2( 78) /  22.30618 /
      DATA XCUM2( 79) /  23.19862 /   ,   XCUM2( 80) /  24.16809 /
      DATA XCUM2( 81) /  25.22547 /   ,   XCUM2( 82) /  26.38320 /
      DATA XCUM2( 83) /  27.65658 /   ,   XCUM2( 84) /  29.06490 /
      DATA XCUM2( 85) /  30.63048 /   ,   XCUM2( 86) /  32.38211 /
      DATA XCUM2( 87) /  34.35555 /   ,   XCUM2( 88) /  36.59607 /
      DATA XCUM2( 89) /  39.16212 /   ,   XCUM2( 90) /  42.13195 /
      DATA XCUM2( 91) /  45.61204 /   ,   XCUM2( 92) /  49.74758 /
      DATA XCUM2( 93) /  54.74189 /   ,   XCUM2( 94) /  60.89935 /
      DATA XCUM2( 95) /  68.67370 /   ,   XCUM2( 96) /  78.81480 /
      DATA XCUM2( 97) /  92.61047 /   ,   XCUM2( 98) / 112.50807 /
      DATA XCUM2( 99) / 143.78539 /   ,   XCUM2(100) / 200.00000 /
C.
C.    ------------------------------------------------------------------
C.
	  x = grnd()
      IF(X .LT. 0.004) THEN
*
* *** Extreme left-hand tail
         YRAN = -SQRT(ABS(LOG(X)))
*
      ELSEIF(X.LE.RANGE1) THEN
*
* *** 4-point interpolation in the first cumulative table
         TABPO1 = (X-X1BOT)*GP1INV
         J = TABPO1 + 1
         J = MAX(J,2)
         J = MIN(J,98)
         P = TABPO1 - (J-1)
         A = (P+1.0) * XCUM1(J+2) - (P-2.0) * XCUM1(J-1)
         B = (P-1.0) * XCUM1(J)   - P       * XCUM1(J+1)
         YRAN = A*P*(P-1.0)*0.16666667+B*(P+1.0)*(P-2.0)*0.5
*
      ELSEIF(X.LE.RANGE2) THEN
*
* *** 4-point interpolation in the first cumulative table
         TABPO2 = (X-X2BOT)*GP2INV
         J = TABPO2 + 1
         J = MAX(J,2)
         J = MIN(J,98)
         P = TABPO2 - (J-1)
         A = (P+1.0) * XCUM2(J+2) - (P-2.0) * XCUM2(J-1)
         B = (P-1.0) * XCUM2(J)   - P       * XCUM2(J+1)
         YRAN = A*P*(P-1.0)*0.16666667+B*(P+1.0)*(P-2.0)*0.5
*
      ELSE
*
* *** 1/x**2 sampling for extreme Landau tail
         YRAN=200./grnd()
*
      ENDIF
      END
c
      function GVAVIV(RKAPPA,BETA2,RAN)
	implicit none
	real*8 GVAVIV,RKAPPA,BETA2,RAN
	real*8 AC,HC,H
	real*8 BKMNX1,BKMNY1,BKMNX2,BKMNY2,BKMNX3,BKMNY3,BKMXX1
	real*8 BKMXY1,BKMXX2,BKMXY2,BKMXX3,BKMXY3
	real*8 FBKX1,FBKX2,FBKX3,FBKY1,FBKY2,FBKY3
	real*8 EDGEC,FNINV,DRK,DSIGM,ALFA
	real*8 U1,U2,U3,U4,U5,U6,U7,U8
	real*8 V1,V2,V3,V4,V5,V6,V7,V8
	real*8 W1,W2,W3,W4,W5,W6,W8
	real*8 X,Y,XX,yy,x2,y2,x3,y3
	real*8 xy,p2,p3,q2,q3,pq
	real*8 t,rlam,fl,s,fu,s0
	integer itype,npt,j,n,fn,k
	real*8 wk
	real*8 glands,glande
	
C.
C.    ******************************************************************
C.    *                                                                *
C.    *                                                                *
C.    *  GVAVIV                                                        *
C.    *  ------                                                        *
C.    *                                                                *
C.    *  Initializes the parameters for a Vavilov distribution         *
C.    *                                                                *
C.    *  BETA2 is the particle velocity squared                        *
C.    *  RKAPPA the usual straggling parameter                         *
C.    *                                                                *
C.    *  This routine has been extracted from a proposed CERNLIB       *
C.    *  set of routines (VAVCOE,VAVFSM).                              *
C.    *  The authors of these routines have submitted an article       *
C.    *  to NIM:                                                       *
C.    *                                                                *
C.    *    Alberto Rotondi, Paolo Montagna                             *
C.    *    Fast calculation of Vavilov distribution NIMB 1990          *
C.    *                                                                *
C.    *       Author    K.S.Koelbig     ********                       *
C.    *                                                                *
C.    *   ==> Called by : GLANDO                                       *
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION AC(0:12),HC(0:8),H(9)
 
      PARAMETER
     1(BKMNX1 = 0.02, BKMNY1 = 0.05, BKMNX2 = 0.12, BKMNY2 = 0.05,
     2 BKMNX3 = 0.22, BKMNY3 = 0.05, BKMXX1 = 0.1 , BKMXY1 = 1   ,
     3 BKMXX2 = 0.2 , BKMXY2 = 1   , BKMXX3 = 0.3 , BKMXY3 = 1   )
      PARAMETER
     1(FBKX1 = 2/(BKMXX1-BKMNX1), FBKX2 = 2/(BKMXX2-BKMNX2),
     2 FBKX3 = 2/(BKMXX3-BKMNX3), FBKY1 = 2/(BKMXY1-BKMNY1),
     3 FBKY2 = 2/(BKMXY2-BKMNY2), FBKY3 = 2/(BKMXY3-BKMNY3))
 
      DIMENSION EDGEC(2:7),FNINV(5),DRK(5),DSIGM(5),ALFA(2:5)
      DIMENSION U1(13),U2(13),U3(13),U4(12),U5(13),U6(13),U7( 8),U8(13)
      DIMENSION V1(12),V2(12),V3(13),V4(12),V5(13),V6(13),V7(11),V8(11)
      DIMENSION W1(13),W2(11),W3(13),W4(13),W5(13),W6(13),       W8( 8)
 
      DATA FNINV /1, 0.5, 0.33333333, 0.25, 0.2/
 
      DATA (EDGEC(J),J=2,7)
     1/ 0.16666667E+0, 0.41666667E-1, 0.83333333E-2,
     2  0.13888889E-1, 0.69444444E-2, 0.77160493E-3/
 
      DATA (U1(K),K=1,13)
     1/ 0.25850868E+0,  0.32477982E-1, -0.59020496E-2,
     2  0.            , 0.24880692E-1,  0.47404356E-2,
     3 -0.74445130E-3,  0.73225731E-2,  0.           ,
     4  0.11668284E-2,  0.           , -0.15727318E-2,-0.11210142E-2/
 
      DATA (U2(K),K=1,13)
     1/ 0.43142611E+0,  0.40797543E-1, -0.91490215E-2,
     2  0.           ,  0.42127077E-1,  0.73167928E-2,
     3 -0.14026047E-2,  0.16195241E-1,  0.24714789E-2,
     4  0.20751278E-2,  0.           , -0.25141668E-2,-0.14064022E-2/
 
      DATA (U3(K),K=1,13)
     1/ 0.25225955E+0,  0.64820468E-1, -0.23615759E-1,
     2  0.           ,  0.23834176E-1,  0.21624675E-2,
     3 -0.26865597E-2, -0.54891384E-2,  0.39800522E-2,
     4  0.48447456E-2, -0.89439554E-2, -0.62756944E-2,-0.24655436E-2/
 
      DATA (U4(K),K=1,12)
     1/ 0.12593231E+1, -0.20374501E+0,  0.95055662E-1,
     2 -0.20771531E-1, -0.46865180E-1, -0.77222986E-2,
     3  0.32241039E-2,  0.89882920E-2, -0.67167236E-2,
     4 -0.13049241E-1,  0.18786468E-1,  0.14484097E-1/
 
      DATA (U5(K),K=1,13)
     1/-0.24864376E-1, -0.10368495E-2,  0.14330117E-2,
     2  0.20052730E-3,  0.18751903E-2,  0.12668869E-2,
     3  0.48736023E-3,  0.34850854E-2,  0.           ,
     4 -0.36597173E-3,  0.19372124E-2,  0.70761825E-3, 0.46898375E-3/
 
      DATA (U6(K),K=1,13)
     1/ 0.35855696E-1, -0.27542114E-1,  0.12631023E-1,
     2 -0.30188807E-2, -0.84479939E-3,  0.           ,
     3  0.45675843E-3, -0.69836141E-2,  0.39876546E-2,
     4 -0.36055679E-2,  0.           ,  0.15298434E-2, 0.19247256E-2/
 
      DATA (U7(K),K=1,8)
     1/ 0.10234691E+2, -0.35619655E+1,  0.69387764E+0,
     2 -0.14047599E+0, -0.19952390E+1, -0.45679694E+0,
     3  0.           ,  0.50505298E+0/
 
      DATA (U8(K),K=1,13)
     1/ 0.21487518E+2, -0.11825253E+2,  0.43133087E+1,
     2 -0.14500543E+1, -0.34343169E+1, -0.11063164E+1,
     3 -0.21000819E+0,  0.17891643E+1, -0.89601916E+0,
     4  0.39120793E+0,  0.73410606E+0,  0.           ,-0.32454506E+0/
 
      DATA (V1(K),K=1,12)
     1/ 0.27827257E+0, -0.14227603E-2,  0.24848327E-2,
     2  0.           ,  0.45091424E-1,  0.80559636E-2,
     3 -0.38974523E-2,  0.           , -0.30634124E-2,
     4  0.75633702E-3,  0.54730726E-2,  0.19792507E-2/
 
      DATA (V2(K),K=1,12)
     1/ 0.41421789E+0, -0.30061649E-1,  0.52249697E-2,
     2  0.           ,  0.12693873E+0,  0.22999801E-1,
     3 -0.86792801E-2,  0.31875584E-1, -0.61757928E-2,
     4  0.           ,  0.19716857E-1,  0.32596742E-2/
 
      DATA (V3(K),K=1,13)
     1/ 0.20191056E+0, -0.46831422E-1,  0.96777473E-2,
     2 -0.17995317E-2,  0.53921588E-1,  0.35068740E-2,
     3 -0.12621494E-1, -0.54996531E-2, -0.90029985E-2,
     4  0.34958743E-2,  0.18513506E-1,  0.68332334E-2,-0.12940502E-2/
 
      DATA (V4(K),K=1,12)
     1/ 0.13206081E+1,  0.10036618E+0, -0.22015201E-1,
     2  0.61667091E-2, -0.14986093E+0, -0.12720568E-1,
     3  0.24972042E-1, -0.97751962E-2,  0.26087455E-1,
     4 -0.11399062E-1, -0.48282515E-1, -0.98552378E-2/
 
      DATA (V5(K),K=1,13)
     1/ 0.16435243E-1,  0.36051400E-1,  0.23036520E-2,
     2 -0.61666343E-3, -0.10775802E-1,  0.51476061E-2,
     3  0.56856517E-2, -0.13438433E-1,  0.           ,
     4  0.           , -0.25421507E-2,  0.20169108E-2,-0.15144931E-2/
 
      DATA (V6(K),K=1,13)
     1/ 0.33432405E-1,  0.60583916E-2, -0.23381379E-2,
     2  0.83846081E-3, -0.13346861E-1, -0.17402116E-2,
     3  0.21052496E-2,  0.15528195E-2,  0.21900670E-2,
     4 -0.13202847E-2, -0.45124157E-2, -0.15629454E-2, 0.22499176E-3/
 
      DATA (V7(K),K=1,11)
     1/ 0.54529572E+1, -0.90906096E+0,  0.86122438E-1,
     2  0.           , -0.12218009E+1, -0.32324120E+0,
     3 -0.27373591E-1,  0.12173464E+0,  0.           ,
     4  0.           ,  0.40917471E-1/
 
      DATA (V8(K),K=1,11)
     1/ 0.93841352E+1, -0.16276904E+1,  0.16571423E+0,
     2  0.           , -0.18160479E+1, -0.50919193E+0,
     3 -0.51384654E-1,  0.21413992E+0,  0.           ,
     4  0.           ,  0.66596366E-1/
 
      DATA (W1(K),K=1,13)
     1/ 0.29712951E+0,  0.97572934E-2,  0.           ,
     2 -0.15291686E-2,  0.35707399E-1,  0.96221631E-2,
     3 -0.18402821E-2, -0.49821585E-2,  0.18831112E-2,
     4  0.43541673E-2,  0.20301312E-2, -0.18723311E-2,-0.73403108E-3/
 
      DATA (W2(K),K=1,11)
     1/ 0.40882635E+0,  0.14474912E-1,  0.25023704E-2,
     2 -0.37707379E-2,  0.18719727E+0,  0.56954987E-1,
     3  0.           ,  0.23020158E-1,  0.50574313E-2,
     4  0.94550140E-2,  0.19300232E-1/
 
      DATA (W3(K),K=1,13)
     1/ 0.16861629E+0,  0.           ,  0.36317285E-2,
     2 -0.43657818E-2,  0.30144338E-1,  0.13891826E-1,
     3 -0.58030495E-2, -0.38717547E-2,  0.85359607E-2,
     4  0.14507659E-1,  0.82387775E-2, -0.10116105E-1,-0.55135670E-2/
 
      DATA (W4(K),K=1,13)
     1/ 0.13493891E+1, -0.26863185E-2, -0.35216040E-2,
     2  0.24434909E-1, -0.83447911E-1, -0.48061360E-1,
     3  0.76473951E-2,  0.24494430E-1, -0.16209200E-1,
     4 -0.37768479E-1, -0.47890063E-1,  0.17778596E-1, 0.13179324E-1/
 
      DATA (W5(K),K=1,13)
     1/ 0.10264945E+0,  0.32738857E-1,  0.           ,
     2  0.43608779E-2, -0.43097757E-1, -0.22647176E-2,
     3  0.94531290E-2, -0.12442571E-1, -0.32283517E-2,
     4 -0.75640352E-2, -0.88293329E-2,  0.52537299E-2, 0.13340546E-2/
 
      DATA (W6(K),K=1,13)
     1/ 0.29568177E-1, -0.16300060E-2, -0.21119745E-3,
     2  0.23599053E-2, -0.48515387E-2, -0.40797531E-2,
     3  0.40403265E-3,  0.18200105E-2, -0.14346306E-2,
     4 -0.39165276E-2, -0.37432073E-2,  0.19950380E-2, 0.12222675E-2/
 
      DATA (W8(K),K=1,8)
     1/ 0.66184645E+1, -0.73866379E+0,  0.44693973E-1,
     2  0.           , -0.14540925E+1, -0.39529833E+0,
     3 -0.44293243E-1,  0.88741049E-1/
 
      GVAVIV=0
      IF(RKAPPA .LT. 0.01 .OR. RKAPPA .GT. 12) RETURN
      IF(RKAPPA .GE. 0.29) THEN
       ITYPE=1
       NPT=100
       WK=1/SQRT(RKAPPA)
       AC(0)=(-0.032227*BETA2-0.074275)*RKAPPA+
     1   (0.24533*BETA2+0.070152)*WK+(-0.55610*BETA2-3.1579)
       AC(8)=(-0.013483*BETA2-0.048801)*RKAPPA+
     1   (-1.6921*BETA2+8.3656)*WK+(-0.73275*BETA2-3.5226)
       DRK(1)=WK**2
       DSIGM(1)=SQRT(RKAPPA/(1-0.5*BETA2))
       DO 1 J = 1,4
       DRK(J+1)=DRK(1)*DRK(J)
       DSIGM(J+1)=DSIGM(1)*DSIGM(J)
    1  ALFA(J+1)=(FNINV(J)-BETA2*FNINV(J+1))*DRK(J)
 
       HC(0)=LOG(RKAPPA)+BETA2+0.42278434
       HC(1)=DSIGM(1)
       HC(2)=ALFA(3)*DSIGM(3)
       HC(3)=(3*ALFA(2)**2+ALFA(4))*DSIGM(4)-3
       HC(4)=(10*ALFA(2)*ALFA(3)+ALFA(5))*DSIGM(5)-10*HC(2)
       HC(5)=HC(2)**2
       HC(6)=HC(2)*HC(3)
       HC(7)=HC(2)*HC(5)
       DO 2 J = 2,7
    2  HC(J)=EDGEC(J)*HC(J)
       HC(8)=0.39894228*HC(1)
      ELSEIF(RKAPPA .GE. 0.22) THEN
       ITYPE=2
       NPT=150
       X=1+(RKAPPA-BKMXX3)*FBKX3
       Y=1+(SQRT(BETA2)-BKMXY3)*FBKY3
       XX=2*X
       YY=2*Y
       X2=XX*X-1
       X3=XX*X2-X
       Y2=YY*Y-1
       Y3=YY*Y2-Y
       XY=X*Y
       P2=X2*Y
       P3=X3*Y
       Q2=Y2*X
       Q3=Y3*X
       PQ=X2*Y2
       AC(1)=W1(1)+W1(2)*X+W1(4)*X3+W1(5)*Y+W1(6)*Y2+W1(7)*Y3+
     1  W1(8)*XY+W1(9)*P2+W1(10)*P3+W1(11)*Q2+W1(12)*Q3+W1(13)*PQ
       AC(2)=W2(1)+W2(2)*X+W2(3)*X2+W2(4)*X3+W2(5)*Y+W2(6)*Y2+
     1  W2(8)*XY+W2(9)*P2+W2(10)*P3+W2(11)*Q2
       AC(3)=W3(1)+W3(3)*X2+W3(4)*X3+W3(5)*Y+W3(6)*Y2+W3(7)*Y3+
     1  W3(8)*XY+W3(9)*P2+W3(10)*P3+W3(11)*Q2+W3(12)*Q3+W3(13)*PQ
       AC(4)=W4(1)+W4(2)*X+W4(3)*X2+W4(4)*X3+W4(5)*Y+W4(6)*Y2+W4(7)*Y3+
     1  W4(8)*XY+W4(9)*P2+W4(10)*P3+W4(11)*Q2+W4(12)*Q3+W4(13)*PQ
       AC(5)=W5(1)+W5(2)*X+W5(4)*X3+W5(5)*Y+W5(6)*Y2+W5(7)*Y3+
     1  W5(8)*XY+W5(9)*P2+W5(10)*P3+W5(11)*Q2+W5(12)*Q3+W5(13)*PQ
       AC(6)=W6(1)+W6(2)*X+W6(3)*X2+W6(4)*X3+W6(5)*Y+W6(6)*Y2+W6(7)*Y3+
     1  W6(8)*XY+W6(9)*P2+W6(10)*P3+W6(11)*Q2+W6(12)*Q3+W6(13)*PQ
       AC(8)=W8(1)+W8(2)*X+W8(3)*X2+W8(5)*Y+W8(6)*Y2+W8(7)*Y3+W8(8)*XY
       AC(0)=-3.05
      ELSEIF(RKAPPA .GE. 0.12) THEN
       ITYPE=3
       NPT=200
       X=1+(RKAPPA-BKMXX2)*FBKX2
       Y=1+(SQRT(BETA2)-BKMXY2)*FBKY2
       XX=2*X
       YY=2*Y
       X2=XX*X-1
       X3=XX*X2-X
       Y2=YY*Y-1
       Y3=YY*Y2-Y
       XY=X*Y
       P2=X2*Y
       P3=X3*Y
       Q2=Y2*X
       Q3=Y3*X
       PQ=X2*Y2
       AC(1)=V1(1)+V1(2)*X+V1(3)*X2+V1(5)*Y+V1(6)*Y2+V1(7)*Y3+
     1  V1(9)*P2+V1(10)*P3+V1(11)*Q2+V1(12)*Q3
       AC(2)=V2(1)+V2(2)*X+V2(3)*X2+V2(5)*Y+V2(6)*Y2+V2(7)*Y3+
     1  V2(8)*XY+V2(9)*P2+V2(11)*Q2+V2(12)*Q3
       AC(3)=V3(1)+V3(2)*X+V3(3)*X2+V3(4)*X3+V3(5)*Y+V3(6)*Y2+V3(7)*Y3+
     1  V3(8)*XY+V3(9)*P2+V3(10)*P3+V3(11)*Q2+V3(12)*Q3+V3(13)*PQ
       AC(4)=V4(1)+V4(2)*X+V4(3)*X2+V4(4)*X3+V4(5)*Y+V4(6)*Y2+V4(7)*Y3+
     1  V4(8)*XY+V4(9)*P2+V4(10)*P3+V4(11)*Q2+V4(12)*Q3
       AC(5)=V5(1)+V5(2)*X+V5(3)*X2+V5(4)*X3+V5(5)*Y+V5(6)*Y2+V5(7)*Y3+
     1  V5(8)*XY+V5(11)*Q2+V5(12)*Q3+V5(13)*PQ
       AC(6)=V6(1)+V6(2)*X+V6(3)*X2+V6(4)*X3+V6(5)*Y+V6(6)*Y2+V6(7)*Y3+
     1  V6(8)*XY+V6(9)*P2+V6(10)*P3+V6(11)*Q2+V6(12)*Q3+V6(13)*PQ
       AC(7)=V7(1)+V7(2)*X+V7(3)*X2+V7(5)*Y+V7(6)*Y2+V7(7)*Y3+
     1  V7(8)*XY+V7(11)*Q2
       AC(8)=V8(1)+V8(2)*X+V8(3)*X2+V8(5)*Y+V8(6)*Y2+V8(7)*Y3+
     1  V8(8)*XY+V8(11)*Q2
       AC(0)=-3.04
      ELSE
       ITYPE=4
       IF(RKAPPA .GE. 0.02) ITYPE=3
       NPT=200
       X=1+(RKAPPA-BKMXX1)*FBKX1
       Y=1+(SQRT(BETA2)-BKMXY1)*FBKY1
       XX=2*X
       YY=2*Y
       X2=XX*X-1
       X3=XX*X2-X
       Y2=YY*Y-1
       Y3=YY*Y2-Y
       XY=X*Y
       P2=X2*Y
       P3=X3*Y
       Q2=Y2*X
       Q3=Y3*X
       PQ=X2*Y2
       IF(ITYPE .EQ. 4) GO TO 4
       AC(1)=U1(1)+U1(2)*X+U1(3)*X2+U1(5)*Y+U1(6)*Y2+U1(7)*Y3+
     1  U1(8)*XY+U1(10)*P3+U1(12)*Q3+U1(13)*PQ
       AC(2)=U2(1)+U2(2)*X+U2(3)*X2+U2(5)*Y+U2(6)*Y2+U2(7)*Y3+
     1  U2(8)*XY+U2(9)*P2+U2(10)*P3+U2(12)*Q3+U2(13)*PQ
       AC(3)=U3(1)+U3(2)*X+U3(3)*X2+U3(5)*Y+U3(6)*Y2+U3(7)*Y3+
     1  U3(8)*XY+U3(9)*P2+U3(10)*P3+U3(11)*Q2+U3(12)*Q3+U3(13)*PQ
       AC(4)=U4(1)+U4(2)*X+U4(3)*X2+U4(4)*X3+U4(5)*Y+U4(6)*Y2+U4(7)*Y3+
     1  U4(8)*XY+U4(9)*P2+U4(10)*P3+U4(11)*Q2+U4(12)*Q3
       AC(5)=U5(1)+U5(2)*X+U5(3)*X2+U5(4)*X3+U5(5)*Y+U5(6)*Y2+U5(7)*Y3+
     1  U5(8)*XY+U5(10)*P3+U5(11)*Q2+U5(12)*Q3+U5(13)*PQ
       AC(6)=U6(1)+U6(2)*X+U6(3)*X2+U6(4)*X3+U6(5)*Y+U6(7)*Y3+
     1  U6(8)*XY+U6(9)*P2+U6(10)*P3+U6(12)*Q3+U6(13)*PQ
    4  AC(7)=U7(1)+U7(2)*X+U7(3)*X2+U7(4)*X3+U7(5)*Y+U7(6)*Y2+U7(8)*XY
       AC(8)=U8(1)+U8(2)*X+U8(3)*X2+U8(4)*X3+U8(5)*Y+U8(6)*Y2+U8(7)*Y3+
     1  U8(8)*XY+U8(9)*P2+U8(10)*P3+U8(11)*Q2+U8(13)*PQ
       AC(0)=-3.03
      ENDIF
      AC(9)=(AC(8)-AC(0))/NPT
      IF(ITYPE .EQ. 3) THEN
       X=(AC(7)-AC(8))/(AC(7)*AC(8))
       Y=1/LOG(AC(8)/AC(7))
       P2=AC(7)**2
       AC(11)=P2*(AC(1)*EXP(-AC(2)*(AC(7)+AC(5)*P2)-
     1            AC(3)*EXP(-AC(4)*(AC(7)+AC(6)*P2)))-0.045*Y/AC(7))/
     2            (1+X*Y*AC(7))
       AC(12)=(0.045+X*AC(11))*Y
      ENDIF
      IF(ITYPE .EQ. 4) AC(10)=0.995/GLANDS(AC(8))
 
      T=2*RAN/AC(9)
      RLAM=AC(0)
      FL=0
      S=0
      DO 21 N = 1,NPT
      RLAM=RLAM+AC(9)
      IF(ITYPE .EQ. 1) THEN
       FN=1
       X=(RLAM+HC(0))*HC(1)
       H(1)=X
       H(2)=X**2-1
       DO 31 K = 2,8
       FN=FN+1
   31  H(K+1)=X*H(K)-FN*H(K-1)
       Y=1+HC(7)*H(9)
       DO 32 K = 2,6
   32  Y=Y+HC(K)*H(K+1)
       FU=HC(8)*EXP(-0.5*X**2)*MAX(Y,0.)
      ELSEIF(ITYPE .EQ. 2) THEN
       X=RLAM**2
       FU=AC(1)*EXP(-AC(2)*(RLAM+AC(5)*X)-
     1    AC(3)*EXP(-AC(4)*(RLAM+AC(6)*X)))
      ELSEIF(ITYPE .EQ. 3) THEN
       IF(RLAM .LT. AC(7)) THEN
        X=RLAM**2
        FU=AC(1)*EXP(-AC(2)*(RLAM+AC(5)*X)-
     1     AC(3)*EXP(-AC(4)*(RLAM+AC(6)*X)))
       ELSE
        X=1/RLAM
        FU=(AC(11)*X+AC(12))*X
       ENDIF
      ELSE
       FU=AC(10)*GLANDE(RLAM)
      ENDIF
      S=S+FL+FU
      IF(S .GT. T) GO TO 22
   21 FL=FU
 
   22 S0=S-FL-FU
      GVAVIV=RLAM-AC(9)
      IF(S .GT. S0) GVAVIV=GVAVIV+AC(9)*(T-S0)/(S-S0)
      RETURN
      END
c
      FUNCTION GLANDE(X)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Copy of the CERN library routine DENLAN (G110)                *
C.    *                                                                *
C.    *    ==>Called by : GVACOE                                       *
C.    *                                                                *
C.    ******************************************************************
C.
	real*8 glande
	real*8 P1,P2,P3,P4,P5,P6
	real*8 Q1,Q2,Q3,Q4,Q5,Q6,a1,a2,v,u
      DIMENSION P1(0:4),P2(0:4),P3(0:4),P4(0:4),P5(0:4),P6(0:4)
      DIMENSION Q1(0:4),Q2(0:4),Q3(0:4),Q4(0:4),Q5(0:4),Q6(0:4)
      DIMENSION A1(1:3),A2(1:2)
C
      DATA (P1(I),I=0,4),(Q1(J),J=0,4)
     1/ 0.42598 94875E+0,-0.12497 62550E+0, 0.39842 43700E-1,
     2 -0.62982 87635E-2, 0.15111 62253E-2,
     3  1.0             ,-0.33882 60629E+0, 0.95943 93323E-1,
     4 -0.16080 42283E-1, 0.37789 42063E-2/
C
      DATA (P2(I),I=0,4),(Q2(J),J=0,4)
     1/ 0.17885 41609E+0, 0.11739 57403E+0, 0.14888 50518E-1,
     2 -0.13949 89411E-2, 0.12836 17211E-3,
     3  1.0             , 0.74287 95082E+0, 0.31539 32961E+0,
     4  0.66942 19548E-1, 0.87906 09714E-2/
C
      DATA (P3(I),I=0,4),(Q3(J),J=0,4)
     1/ 0.17885 44503E+0, 0.93591 61662E-1, 0.63253 87654E-2,
     2  0.66116 67319E-4,-0.20310 49101E-5,
     3  1.0             , 0.60978 09921E+0, 0.25606 16665E+0,
     4  0.47467 22384E-1, 0.69573 01675E-2/
C
      DATA (P4(I),I=0,4),(Q4(J),J=0,4)
     1/ 0.98740 54407E+0, 0.11867 23273E+3, 0.84927 94360E+3,
     2 -0.74377 92444E+3, 0.42702 62186E+3,
     3  1.0             , 0.10686 15961E+3, 0.33764 96214E+3,
     4  0.20167 12389E+4, 0.15970 63511E+4/
C
      DATA (P5(I),I=0,4),(Q5(J),J=0,4)
     1/ 0.10036 75074E+1, 0.16757 02434E+3, 0.47897 11289E+4,
     2  0.21217 86767E+5,-0.22324 94910E+5,
     3  1.0             , 0.15694 24537E+3, 0.37453 10488E+4,
     4  0.98346 98876E+4, 0.66924 28357E+5/
C
      DATA (P6(I),I=0,4),(Q6(J),J=0,4)
     1/ 0.10008 27619E+1, 0.66491 43136E+3, 0.62972 92665E+5,
     2  0.47555 46998E+6,-0.57436 09109E+7,
     3  1.0             , 0.65141 01098E+3, 0.56974 73333E+5,
     4  0.16591 74725E+6,-0.28157 59939E+7/
C
      DATA (A1(I),I=1,3)
     1/ 0.41666 66667E-1,-0.19965 27778E-1, 0.27095 38966E-1/
C
      DATA (A2(I),I=1,2)
     1/-0.18455 68670E+1,-0.42846 40743E+1/
C
      V=X
      IF(V .LT. -5.5) THEN
       U=EXP(V+1.0)
       GLANDE=0.3989422803*(EXP(-1.0/U)/SQRT(U))*
     1        (1.0+(A1(1)+(A1(2)+A1(3)*U)*U)*U)
      ELSE IF(V .LT. -1.0) THEN
       U=EXP(-V-1.0)
       GLANDE=EXP(-U)*SQRT(U)*
     1        (P1(0)+(P1(1)+(P1(2)+(P1(3)+P1(4)*V)*V)*V)*V)/
     2        (Q1(0)+(Q1(1)+(Q1(2)+(Q1(3)+Q1(4)*V)*V)*V)*V)
      ELSE IF(V .LT. 1.0) THEN
       GLANDE=(P2(0)+(P2(1)+(P2(2)+(P2(3)+P2(4)*V)*V)*V)*V)/
     1        (Q2(0)+(Q2(1)+(Q2(2)+(Q2(3)+Q2(4)*V)*V)*V)*V)
      ELSE IF(V .LT. 5.0) THEN
       GLANDE=(P3(0)+(P3(1)+(P3(2)+(P3(3)+P3(4)*V)*V)*V)*V)/
     1        (Q3(0)+(Q3(1)+(Q3(2)+(Q3(3)+Q3(4)*V)*V)*V)*V)
      ELSE IF(V .LT. 12.0) THEN
       U=1.0/V
       GLANDE=U**2*(P4(0)+(P4(1)+(P4(2)+(P4(3)+P4(4)*U)*U)*U)*U)/
     1             (Q4(0)+(Q4(1)+(Q4(2)+(Q4(3)+Q4(4)*U)*U)*U)*U)
      ELSE IF(V .LT. 50.0) THEN
       U=1.0/V
       GLANDE=U**2*(P5(0)+(P5(1)+(P5(2)+(P5(3)+P5(4)*U)*U)*U)*U)/
     1             (Q5(0)+(Q5(1)+(Q5(2)+(Q5(3)+Q5(4)*U)*U)*U)*U)
      ELSE IF(V .LT. 300.0) THEN
       U=1.0/V
       GLANDE=U**2*(P6(0)+(P6(1)+(P6(2)+(P6(3)+P6(4)*U)*U)*U)*U)/
     1             (Q6(0)+(Q6(1)+(Q6(2)+(Q6(3)+Q6(4)*U)*U)*U)*U)
      ELSE
       U=1.0/(V-V*LOG(V)/(V+1.0))
       GLANDE=U**2*(1.0+(A2(1)+A2(2)*U)*U)
      END IF
      RETURN
      END
c
c
      FUNCTION GLANDS(X)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Copy of the CERN library routine DSTLAN (G110)                *
C.    *                                                                *
C.    *    ==>Called by : GVACOE                                       *
C.    *                                                                *
C.    ******************************************************************
C.
	real*8 glands
	real*8 P1,P2,P3,P4,P5,P6
	real*8 Q1,Q2,Q3,Q4,Q5,Q6,a1,a2,v,u
      DIMENSION P1(0:4),P2(0:3),P3(0:3),P4(0:3),P5(0:3),P6(0:3)
      DIMENSION Q1(0:4),Q2(0:3),Q3(0:3),Q4(0:3),Q5(0:3),Q6(0:3)
      DIMENSION A1(1:3),A2(1:3)
C
      DATA (P1(I),I=0,4),(Q1(J),J=0,4)
     1/ 0.25140 91491E+0,-0.62505 80444E-1, 0.14583 81230E-1,
     2 -0.21088 17737E-2, 0.74112 47290E-3,
     3  1.0             ,-0.55711 75625E-2, 0.62253 10236E-1,
     4 -0.31373 78427E-2, 0.19314 96439E-2/
C
      DATA (P2(I),I=0,3),(Q2(J),J=0,3)
     1/ 0.28683 28584E+0, 0.35643 63231E+0, 0.15235 18695E+0,
     2  0.22513 04883E-1,
     3  1.0             , 0.61911 36137E+0, 0.17207 21448E+0,
     4  0.22785 94771E-1/
C
      DATA (P3(I),I=0,3),(Q3(J),J=0,3)
     1/ 0.28683 29066E+0, 0.30038 28436E+0, 0.99509 51941E-1,
     2  0.87338 27185E-2,
     3  1.0             , 0.42371 90502E+0, 0.10956 31512E+0,
     4  0.86938 51567E-2/
C
      DATA (P4(I),I=0,3),(Q4(J),J=0,3)
     1/ 0.10003 51630E+1, 0.45035 92498E+1, 0.10858 83880E+2,
     2  0.75360 52269E+1,
     3  1.0             , 0.55399 69678E+1, 0.19335 81111E+2,
     4  0.27213 21508E+2/
C
      DATA (P5(I),I=0,3),(Q5(J),J=0,3)
     1/ 0.10000 06517E+1, 0.49094 14111E+2, 0.85055 44753E+2,
     2  0.15321 53455E+3,
     3  1.0             , 0.50099 28881E+2, 0.13998 19104E+3,
     4  0.42000 02909E+3/
C
      DATA (P6(I),I=0,3),(Q6(J),J=0,3)
     1/ 0.10000 00983E+1, 0.13298 68456E+3, 0.91621 49244E+3,
     2 -0.96050 54274E+3,
     3  1.0             , 0.13398 87843E+3, 0.10559 90413E+4,
     4  0.55322 24619E+3/
C
      DATA (A1(I),I=1,3)
     1/-0.45833 33333E+0, 0.66753 47222E+0,-0.16417 41416E+1/
C
      DATA (A2(I),I=1,3)
     1/ 1.0             ,-0.42278 43351E+0,-0.20434 03138E+1/
C
      V=X
      IF(V .LT. -5.5) THEN
       U=EXP(V+1.0)
       GLANDS=0.3989422803*EXP(-1.0/U)*SQRT(U)*
     1        (1.0+(A1(1)+(A1(2)+A1(3)*U)*U)*U)
      ELSE IF(V .LT. -1.0) THEN
       U=EXP(-V-1.0)
       GLANDS=(EXP(-U)/SQRT(U))*
     1        (P1(0)+(P1(1)+(P1(2)+(P1(3)+P1(4)*V)*V)*V)*V)/
     2        (Q1(0)+(Q1(1)+(Q1(2)+(Q1(3)+Q1(4)*V)*V)*V)*V)
      ELSE IF(V .LT. 1.0) THEN
       GLANDS=(P2(0)+(P2(1)+(P2(2)+P2(3)*V)*V)*V)/
     1        (Q2(0)+(Q2(1)+(Q2(2)+Q2(3)*V)*V)*V)
      ELSE IF(V .LT. 4.0) THEN
       GLANDS=(P3(0)+(P3(1)+(P3(2)+P3(3)*V)*V)*V)/
     1        (Q3(0)+(Q3(1)+(Q3(2)+Q3(3)*V)*V)*V)
      ELSE IF(V .LT. 12.0) THEN
       U=1.0/V
       GLANDS=(P4(0)+(P4(1)+(P4(2)+P4(3)*U)*U)*U)/
     1        (Q4(0)+(Q4(1)+(Q4(2)+Q4(3)*U)*U)*U)
      ELSE IF(V .LT. 50.0) THEN
       U=1.0/V
       GLANDS=(P5(0)+(P5(1)+(P5(2)+P5(3)*U)*U)*U)/
     1        (Q5(0)+(Q5(1)+(Q5(2)+Q5(3)*U)*U)*U)
      ELSE IF(V .LT. 300.0) THEN
       U=1.0/V
       GLANDS=(P6(0)+(P6(1)+(P6(2)+P6(3)*U)*U)*U)/
     1        (Q6(0)+(Q6(1)+(Q6(2)+Q6(3)*U)*U)*U)
      ELSE
       U=1.0/(V-V*LOG(V)/(V+1.0))
       GLANDS=1.0-(A2(1)+(A2(2)+A2(3)*U)*U)*U
      END IF
      RETURN
      END
c	
c
C MEP:
C This is blatantly stealing the random landau function from cernlib 
C (I know cernlib is also loaded, but I am trying to slowly make simc 
C not cernlib dependant, so..)
	
      FUNCTION RANLAN1(X)
	  
	  implicit none
  	  
	  INTEGER I
	  REAL*8 RANLAN1
	  REAL*8 X,U,V
C      DIMENSION F(6:982)
	  REAL*8 F(6:982)
  
      DATA (F(I),I= 6,100)
     A/                                                  -2.244733,
     B -2.204365,-2.168163,-2.135219,-2.104898,-2.076740,-2.050397,
     C -2.025605,-2.002150,-1.979866,-1.958612,-1.938275,-1.918760,
     D -1.899984,-1.881879,-1.864385,-1.847451,-1.831030,-1.815083,
     E -1.799574,-1.784473,-1.769751,-1.755383,-1.741346,-1.727620,
     F -1.714187,-1.701029,-1.688130,-1.675477,-1.663057,-1.650858,
     G -1.638868,-1.627078,-1.615477,-1.604058,-1.592811,-1.581729,
     H -1.570806,-1.560034,-1.549407,-1.538919,-1.528565,-1.518339,
     I -1.508237,-1.498254,-1.488386,-1.478628,-1.468976,-1.459428,
     J -1.449979,-1.440626,-1.431365,-1.422195,-1.413111,-1.404112,
     K -1.395194,-1.386356,-1.377594,-1.368906,-1.360291,-1.351746,
     L -1.343269,-1.334859,-1.326512,-1.318229,-1.310006,-1.301843,
     M -1.293737,-1.285688,-1.277693,-1.269752,-1.261863,-1.254024,
     N -1.246235,-1.238494,-1.230800,-1.223153,-1.215550,-1.207990,
     O -1.200474,-1.192999,-1.185566,-1.178172,-1.170817,-1.163500,
     P -1.156220,-1.148977,-1.141770,-1.134598,-1.127459,-1.120354,
     Q -1.113282,-1.106242,-1.099233,-1.092255/
  
      DATA (F(I),I=101,200)
     A/-1.085306,-1.078388,-1.071498,-1.064636,-1.057802,-1.050996,
     B -1.044215,-1.037461,-1.030733,-1.024029,-1.017350,-1.010695,
     C -1.004064, -.997456, -.990871, -.984308, -.977767, -.971247,
     D  -.964749, -.958271, -.951813, -.945375, -.938957, -.932558,
     E  -.926178, -.919816, -.913472, -.907146, -.900838, -.894547,
     F  -.888272, -.882014, -.875773, -.869547, -.863337, -.857142,
     G  -.850963, -.844798, -.838648, -.832512, -.826390, -.820282,
     H  -.814187, -.808106, -.802038, -.795982, -.789940, -.783909,
     I  -.777891, -.771884, -.765889, -.759906, -.753934, -.747973,
     J  -.742023, -.736084, -.730155, -.724237, -.718328, -.712429,
     K  -.706541, -.700661, -.694791, -.688931, -.683079, -.677236,
     L  -.671402, -.665576, -.659759, -.653950, -.648149, -.642356,
     M  -.636570, -.630793, -.625022, -.619259, -.613503, -.607754,
     N  -.602012, -.596276, -.590548, -.584825, -.579109, -.573399,
     O  -.567695, -.561997, -.556305, -.550618, -.544937, -.539262,
     P  -.533592, -.527926, -.522266, -.516611, -.510961, -.505315,
     Q  -.499674, -.494037, -.488405, -.482777/
  
      DATA (F(I),I=201,300)
     A/ -.477153, -.471533, -.465917, -.460305, -.454697, -.449092,
     B  -.443491, -.437893, -.432299, -.426707, -.421119, -.415534,
     C  -.409951, -.404372, -.398795, -.393221, -.387649, -.382080,
     D  -.376513, -.370949, -.365387, -.359826, -.354268, -.348712,
     E  -.343157, -.337604, -.332053, -.326503, -.320955, -.315408,
     F  -.309863, -.304318, -.298775, -.293233, -.287692, -.282152,
     G  -.276613, -.271074, -.265536, -.259999, -.254462, -.248926,
     H  -.243389, -.237854, -.232318, -.226783, -.221247, -.215712,
     I  -.210176, -.204641, -.199105, -.193568, -.188032, -.182495,
     J  -.176957, -.171419, -.165880, -.160341, -.154800, -.149259,
     K  -.143717, -.138173, -.132629, -.127083, -.121537, -.115989,
     L  -.110439, -.104889, -.099336, -.093782, -.088227, -.082670,
     M  -.077111, -.071550, -.065987, -.060423, -.054856, -.049288,
     N  -.043717, -.038144, -.032569, -.026991, -.021411, -.015828,
     O  -.010243, -.004656,  .000934,  .006527,  .012123,  .017722,
     P   .023323,  .028928,  .034535,  .040146,  .045759,  .051376,
     Q   .056997,  .062620,  .068247,  .073877/
  
      DATA (F(I),I=301,400)
     A/  .079511,  .085149,  .090790,  .096435,  .102083,  .107736,
     B   .113392,  .119052,  .124716,  .130385,  .136057,  .141734,
     C   .147414,  .153100,  .158789,  .164483,  .170181,  .175884,
     D   .181592,  .187304,  .193021,  .198743,  .204469,  .210201,
     E   .215937,  .221678,  .227425,  .233177,  .238933,  .244696,
     F   .250463,  .256236,  .262014,  .267798,  .273587,  .279382,
     G   .285183,  .290989,  .296801,  .302619,  .308443,  .314273,
     H   .320109,  .325951,  .331799,  .337654,  .343515,  .349382,
     I   .355255,  .361135,  .367022,  .372915,  .378815,  .384721,
     J   .390634,  .396554,  .402481,  .408415,  .414356,  .420304,
     K   .426260,  .432222,  .438192,  .444169,  .450153,  .456145,
     L   .462144,  .468151,  .474166,  .480188,  .486218,  .492256,
     M   .498302,  .504356,  .510418,  .516488,  .522566,  .528653,
     N   .534747,  .540850,  .546962,  .553082,  .559210,  .565347,
     O   .571493,  .577648,  .583811,  .589983,  .596164,  .602355,
     P   .608554,  .614762,  .620980,  .627207,  .633444,  .639689,
     Q   .645945,  .652210,  .658484,  .664768/
  
      DATA (F(I),I=401,500)
     A/  .671062,  .677366,  .683680,  .690004,  .696338,  .702682,
     B   .709036,  .715400,  .721775,  .728160,  .734556,  .740963,
     C   .747379,  .753807,  .760246,  .766695,  .773155,  .779627,
     D   .786109,  .792603,  .799107,  .805624,  .812151,  .818690,
     E   .825241,  .831803,  .838377,  .844962,  .851560,  .858170,
     F   .864791,  .871425,  .878071,  .884729,  .891399,  .898082,
     G   .904778,  .911486,  .918206,  .924940,  .931686,  .938446,
     H   .945218,  .952003,  .958802,  .965614,  .972439,  .979278,
     I   .986130,  .992996,  .999875, 1.006769, 1.013676, 1.020597,
     J  1.027533, 1.034482, 1.041446, 1.048424, 1.055417, 1.062424,
     K  1.069446, 1.076482, 1.083534, 1.090600, 1.097681, 1.104778,
     L  1.111889, 1.119016, 1.126159, 1.133316, 1.140490, 1.147679,
     M  1.154884, 1.162105, 1.169342, 1.176595, 1.183864, 1.191149,
     N  1.198451, 1.205770, 1.213105, 1.220457, 1.227826, 1.235211,
     O  1.242614, 1.250034, 1.257471, 1.264926, 1.272398, 1.279888,
     P  1.287395, 1.294921, 1.302464, 1.310026, 1.317605, 1.325203,
     Q  1.332819, 1.340454, 1.348108, 1.355780/
  
      DATA (F(I),I=501,600)
     A/ 1.363472, 1.371182, 1.378912, 1.386660, 1.394429, 1.402216,
     B  1.410024, 1.417851, 1.425698, 1.433565, 1.441453, 1.449360,
     C  1.457288, 1.465237, 1.473206, 1.481196, 1.489208, 1.497240,
     D  1.505293, 1.513368, 1.521465, 1.529583, 1.537723, 1.545885,
     E  1.554068, 1.562275, 1.570503, 1.578754, 1.587028, 1.595325,
     F  1.603644, 1.611987, 1.620353, 1.628743, 1.637156, 1.645593,
     G  1.654053, 1.662538, 1.671047, 1.679581, 1.688139, 1.696721,
     H  1.705329, 1.713961, 1.722619, 1.731303, 1.740011, 1.748746,
     I  1.757506, 1.766293, 1.775106, 1.783945, 1.792810, 1.801703,
     J  1.810623, 1.819569, 1.828543, 1.837545, 1.846574, 1.855631,
     K  1.864717, 1.873830, 1.882972, 1.892143, 1.901343, 1.910572,
     L  1.919830, 1.929117, 1.938434, 1.947781, 1.957158, 1.966566,
     M  1.976004, 1.985473, 1.994972, 2.004503, 2.014065, 2.023659,
     N  2.033285, 2.042943, 2.052633, 2.062355, 2.072110, 2.081899,
     O  2.091720, 2.101575, 2.111464, 2.121386, 2.131343, 2.141334,
     P  2.151360, 2.161421, 2.171517, 2.181648, 2.191815, 2.202018,
     Q  2.212257, 2.222533, 2.232845, 2.243195/
  
      DATA (F(I),I=601,700)
     A/ 2.253582, 2.264006, 2.274468, 2.284968, 2.295507, 2.306084,
     B  2.316701, 2.327356, 2.338051, 2.348786, 2.359562, 2.370377,
     C  2.381234, 2.392131, 2.403070, 2.414051, 2.425073, 2.436138,
     D  2.447246, 2.458397, 2.469591, 2.480828, 2.492110, 2.503436,
     E  2.514807, 2.526222, 2.537684, 2.549190, 2.560743, 2.572343,
     F  2.583989, 2.595682, 2.607423, 2.619212, 2.631050, 2.642936,
     G  2.654871, 2.666855, 2.678890, 2.690975, 2.703110, 2.715297,
     H  2.727535, 2.739825, 2.752168, 2.764563, 2.777012, 2.789514,
     I  2.802070, 2.814681, 2.827347, 2.840069, 2.852846, 2.865680,
     J  2.878570, 2.891518, 2.904524, 2.917588, 2.930712, 2.943894,
     K  2.957136, 2.970439, 2.983802, 2.997227, 3.010714, 3.024263,
     L  3.037875, 3.051551, 3.065290, 3.079095, 3.092965, 3.106900,
     M  3.120902, 3.134971, 3.149107, 3.163312, 3.177585, 3.191928,
     N  3.206340, 3.220824, 3.235378, 3.250005, 3.264704, 3.279477,
     O  3.294323, 3.309244, 3.324240, 3.339312, 3.354461, 3.369687,
     P  3.384992, 3.400375, 3.415838, 3.431381, 3.447005, 3.462711,
     Q  3.478500, 3.494372, 3.510328, 3.526370/
  
      DATA (F(I),I=701,800)
     A/ 3.542497, 3.558711, 3.575012, 3.591402, 3.607881, 3.624450,
     B  3.641111, 3.657863, 3.674708, 3.691646, 3.708680, 3.725809,
     C  3.743034, 3.760357, 3.777779, 3.795300, 3.812921, 3.830645,
     D  3.848470, 3.866400, 3.884434, 3.902574, 3.920821, 3.939176,
     E  3.957640, 3.976215, 3.994901, 4.013699, 4.032612, 4.051639,
     F  4.070783, 4.090045, 4.109425, 4.128925, 4.148547, 4.168292,
     G  4.188160, 4.208154, 4.228275, 4.248524, 4.268903, 4.289413,
     H  4.310056, 4.330832, 4.351745, 4.372794, 4.393982, 4.415310,
     I  4.436781, 4.458395, 4.480154, 4.502060, 4.524114, 4.546319,
     J  4.568676, 4.591187, 4.613854, 4.636678, 4.659662, 4.682807,
     K  4.706116, 4.729590, 4.753231, 4.777041, 4.801024, 4.825179,
     L  4.849511, 4.874020, 4.898710, 4.923582, 4.948639, 4.973883,
     M  4.999316, 5.024942, 5.050761, 5.076778, 5.102993, 5.129411,
     N  5.156034, 5.182864, 5.209903, 5.237156, 5.264625, 5.292312,
     O  5.320220, 5.348354, 5.376714, 5.405306, 5.434131, 5.463193,
     P  5.492496, 5.522042, 5.551836, 5.581880, 5.612178, 5.642734,
     Q  5.673552, 5.704634, 5.735986, 5.767610/
  
      DATA (F(I),I=801,900)
     A/ 5.799512, 5.831694, 5.864161, 5.896918, 5.929968, 5.963316,
     B  5.996967, 6.030925, 6.065194, 6.099780, 6.134687, 6.169921,
     C  6.205486, 6.241387, 6.277630, 6.314220, 6.351163, 6.388465,
     D  6.426130, 6.464166, 6.502578, 6.541371, 6.580553, 6.620130,
     E  6.660109, 6.700495, 6.741297, 6.782520, 6.824173, 6.866262,
     F  6.908795, 6.951780, 6.995225, 7.039137, 7.083525, 7.128398,
     G  7.173764, 7.219632, 7.266011, 7.312910, 7.360339, 7.408308,
     H  7.456827, 7.505905, 7.555554, 7.605785, 7.656608, 7.708035,
     I  7.760077, 7.812747, 7.866057, 7.920019, 7.974647, 8.029953,
     J  8.085952, 8.142657, 8.200083, 8.258245, 8.317158, 8.376837,
     K  8.437300, 8.498562, 8.560641, 8.623554, 8.687319, 8.751955,
     L  8.817481, 8.883916, 8.951282, 9.019600, 9.088889, 9.159174,
     M  9.230477, 9.302822, 9.376233, 9.450735, 9.526355, 9.603118,
     N  9.681054, 9.760191, 9.840558, 9.922186,10.005107,10.089353,
     O 10.174959,10.261958,10.350389,10.440287,10.531693,10.624646,
     P 10.719188,10.815362,10.913214,11.012789,11.114137,11.217307,
     Q 11.322352,11.429325,11.538283,11.649285/
  
      DATA (F(I),I=901,982)
     A/11.762390,11.877664,11.995170,12.114979,12.237161,12.361791,
     B 12.488946,12.618708,12.751161,12.886394,13.024498,13.165570,
     C 13.309711,13.457026,13.607625,13.761625,13.919145,14.080314,
     D 14.245263,14.414134,14.587072,14.764233,14.945778,15.131877,
     E 15.322712,15.518470,15.719353,15.925570,16.137345,16.354912,
     F 16.578520,16.808433,17.044929,17.288305,17.538873,17.796967,
     G 18.062943,18.337176,18.620068,18.912049,19.213574,19.525133,
     H 19.847249,20.180480,20.525429,20.882738,21.253102,21.637266,
     I 22.036036,22.450278,22.880933,23.329017,23.795634,24.281981,
     J 24.789364,25.319207,25.873062,26.452634,27.059789,27.696581,
     K 28.365274,29.068370,29.808638,30.589157,31.413354,32.285060,
     L 33.208568,34.188705,35.230920,36.341388,37.527131,38.796172,
     M 40.157721,41.622399,43.202525,44.912465,46.769077,48.792279,
     N 51.005773,53.437996,56.123356,59.103894/
  
      U=1000*X
      I=U
      U=U-I
      IF(I .GE. 70 .AND. I .LE. 800) THEN
       RANLAN1=F(I)+U*(F(I+1)-F(I))
      ELSEIF(I .GE.  7 .AND. I .LE. 980) THEN
       RANLAN1=
     1  F(I)+U*(F(I+1)-F(I)-0.25*(1-U)*(F(I+2)-F(I+1)-F(I)+F(I-1)))
      ELSEIF(I. LT. 7) THEN
       V=LOG(X)
       U=1/V
       RANLAN1=((0.99858950+(3.45213058E1+1.70854528E1*U)*U)/
     1         (1         +(3.41760202E1+4.01244582  *U)*U))*
     2         (-LOG(-0.91893853-V)-1)
      ELSE
       U=1-X
       V=U**2
       IF(X .LE. 0.999) THEN
        RANLAN1=(1.00060006+2.63991156E2*U+4.37320068E3*V)/
     1        ((1         +2.57368075E2*U+3.41448018E3*V)*U)
       ELSE
        RANLAN1=(1.00001538+6.07514119E3*U+7.34266409E5*V)/
     1        ((1         +6.06511919E3*U+6.94021044E5*V)*U)
       ENDIF
      ENDIF
      RETURN
      END
