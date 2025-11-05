	subroutine musc(m2,p,rad_len,dth,dph)
C+_____________________________________________________________________
!
! MUSC - Simulate multiple scattering of any particle.
!
! ASSUMPTIONS: DTH and DPH given in milli-radians, RAD_LEN in radiation
!   lengths. The formula used is due to Rossi and Greisen (See the book
!   by Segre, NUCLEI AND PARTICLES, 1982, p. 48.) The formula assumes a
!   gaussian distribution for the scattering angle. It is further assumed
!   that the angles DTH and DPH are located at an angle with respect 
!   to the beam direction that is large enough so that DTH-DPH space can
!   be approximated as a cartesian space.
!
! D. Potterveld - Sept. 1985
!
! Add option for protons 
!   Precision not good enough for our thick targets. Latest values supplied:
!   Lynch and Dahl, NIM B58 (1991) 6.
!   Note, that there is a typo in Particle data booklet: eps = 0.088 instead
!   of 0.2! Precision: -5% deviation for H2, +8 for 238U at radiation 
!   lengths between ~ 0.01 and 1.
!
! H.J. Bulten - Aug. 1991
C-_____________________________________________________________________

	implicit none

	real*8 Es, epsilon
	parameter (Es = 13.6)		!MeV
	parameter (epsilon = 0.088)

	real*8 rad_len, dth, dph
	real*8 beta, theta_sigma
	real*8 m2, p
	real*8 B, sigr, signew, ampr, xtrans, powc, bp
	real*8 integral1, integral2, integral3, val, retval, sigv
	real*8 ymin, ymax, yval

	real*8 nsig_max
	parameter(nsig_max=99.0e0)      !max #/sigma for gaussian ran #s.

	real*8 gauss1, grnd

! Compute scattering angles, THETA_SCAT from a gaussian distribution,
! PHI_SCAT from uniform distribution.

	if (rad_len.eq.0) return
	if (p.lt.25.) write(6,*)
     >		'Momentum passed to musc.f should be in MeV, but p=',p

	beta = p / sqrt(m2+p*p)
c	theta_sigma = Es/p/beta * sqrt(rad_len) * (1+epsilon*log10(rad_len))
c Better form for beta .ne. 1
	theta_sigma = Es/p/beta * sqrt(rad_len) * (1+epsilon*log10(rad_len/beta**2))
	
c Now we are going to calculate the extended tails of the Moliere distribution to
c better represent the multiple scattering.  This is done by approximating the tail
c a 1/x^B distribution, and creating a continuous distribution from the standard 
c gaussian at < 1 sigma_original, a wider "patch" gaussian at sigma_orig < 2 and
c then the 1/x^B distribution.
c the B factor would typically be calculated using the A, Z, and thickness of
c material (not rad_len), but we will estimate it from the rad_len, A = Z = 1,
c and p = 1GeV.  Note that B is fairly insensitive to p, but at large Z (~50), the
c the provide estimation OVER estimates the B factor by about 50%.  This means
c materials with larger Z will have larger musc tails than in reality.

	B = -57.0121/log(rad_len) -3.13946 -328.837*rad_len + 8665.29*rad_len*rad_len
	if (B.GT.6) then
		B = 6
	end if
	
	B = B/sqrt(2.)

c with three "regions" each described by a different function, we must calculate the
c the integrals to determine probibilities.

	sigr =  (3.0435 -0.553706*B + 0.038984*B*B)
	signew = theta_sigma/sigr*sqrt(2.)
	ampr = 0.53768*sigr -0.463712

	xtrans =  (1.0/log(B)*4.9)*signew
	powc =   B*(2.11052 -0.313557*B + 0.0156546*B*B)


	bp = (exp(-0.5*xtrans**2/signew**2) + ampr*exp(-0.5*xtrans**2/(2.0*signew)**2))*xtrans**powc
c combined gaussian up to xtrans
	integral1 = sqrt(3.14159/2.0)*signew*(2.0*ampr*erf(xtrans/(2.0*sqrt(2.0)*signew)) + erf(xtrans/sqrt(2.0)/signew))

c power function integral from xtrans to infinity
	integral2 = bp/(powc -1.0)*xtrans**(1.0 - powc)

c wider gaussian integral, to determine when to throw according to this gaussian.
	integral3 = sqrt(3.14159/2.0)*signew*(2.0*ampr*erf(xtrans/(2.0*sqrt(2.0)*signew)))


c we choose which section to use based on which fraction we get:
	val = integral1/(integral1 + integral2);
	
	retval = 0
	sigv = signew
	
	if(grnd().LT.val) then
		if((grnd()*integral1).LT.integral3) then
			sigv = 2.0*signew
		end if
		retval = 2.0*xtrans
		do while (abs(retval).GT.xtrans)
			retval = sigv*gauss1(nsig_max)
		end do
	else
		ymin = 0.0;
		ymax = bp/(integral2*(powc-1.0))*xtrans**(1.0-powc);
		yval = ymin + grnd()*(ymax - ymin)
		retval = (yval*(powc-1.0)/bp*integral2)**(1.0/(1.0-powc));

		if(grnd() > 0.5) then
			retval = -retval;
		end if
	end if
	
	dth = dth + retval
	
c do again for phi

	retval = 0
	sigv = signew
	
	if(grnd().LT.val) then
		if((grnd()*integral1).LT.integral3) then
			sigv = 2.0*signew
		end if
		retval = 2.0*xtrans
		do while (abs(retval).GT.xtrans)
			retval = sigv*gauss1(nsig_max)
		end do
	else
		ymin = 0.0;
		ymax = bp/(integral2*(powc-1.0))*xtrans**(1.0-powc);
		yval = ymin + grnd()*(ymax - ymin)
		retval = (yval*(powc-1.0)/bp*integral2)**(1.0/(1.0-powc));

		if(grnd() > 0.5) then
			retval = -retval;
		end if
	end if
	
	dph = dph + retval

	return
	end
