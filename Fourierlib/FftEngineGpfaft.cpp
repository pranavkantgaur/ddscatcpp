#include "StdAfx.h"

#include "FftEngineGpfaft.h"

/* **
GPFAPACK - FORTRAN IMPLEMENTATION OF THE SELF-SORTING IN-PLACE GENERALIZED PRIME FACTOR (COMPLEX) FFT [GPFA]
WRITTEN BY CLIVE TEMPERTON RECHERCHE EN PREVISION NUMERIQUE / ECMWF
THE PACKAGE CONSISTS OF THE SETUP ROUTINE SETGPFA, TOGETHER WITH THE ROUTINES GPFA, GPFA2F, GPFA3F, GPFA5F
** */

void Debugga(Complex *c, int L)
{
	int i, j, n, ja, jb, jc;
    
    n = 96;
    FILE *file9 = fopen("aaa.txt", "w");
    for(i=0; i<L; ++i)
	{
	    ja = i % n;
	    jb = (i/n) % n;
		jc = i/(n*n);
		j = jc + n*(jb + ja*n);
		fprintf(file9, "%8d%8d%8d%8d%13.5f%13.5f\n", i, ja, jb, jc, c[i].re, c[i].im);
	}
	fclose(file9);
}

FftEngineGpfaft::FftEngineGpfaft(void)
{
	method = FftMethod_GPFAFT;
	trigX = trigY = trigZ = NULL;
	mxold = myold = mzold = 0;
}

FftEngineGpfaft::~FftEngineGpfaft(void)
{
	CleanDelete2(trigX);
	CleanDelete2(trigY);
	CleanDelete2(trigZ);
}

void FftEngineGpfaft::ElementaryInit(unsigned int &n, unsigned int N, Complex *&trig)
{
	if (n != N)
	{
		if (n)
		{
			CleanDelete2(trig);
		}
		trig = new Complex[N];
		Setgpfa(trig, N);
		n = N;
	}
}

void FftEngineGpfaft::Init(unsigned int nX, unsigned int nY, unsigned int nZ)
{
	ElementaryInit(nx, nX, trigX);
	ElementaryInit(ny, nY, trigY);
	ElementaryInit(nz, nZ, trigZ);
}

void FftEngineGpfaft::SayHello(char *Buffer)
{
	strcat(Buffer, " - using GPFA package from Clive Temperton");
}

//
// Setup routine for self-sorting in-place generalized prime factor (complex) FFT [gpfa]
bool FftEngineGpfaft::Setgpfa(Complex *trigs, unsigned int n)
{
/* **
!*********************************************************************
!                                                                    *
!     GPFAPACK - FORTRAN IMPLEMENTATION OF THE SELF-SORTING          *
!     IN-PLACE GENERALIZED PRIME FACTOR (COMPLEX) FFT [GPFA]         *
!                                                                    *
!     WRITTEN BY CLIVE TEMPERTON                                     *
!     RECHERCHE EN PREVISION NUMERIQUE / ECMWF                       *
!                                                                    *
!     THE PACKAGE CONSISTS OF THE SETUP ROUTINE SETGPFA, TOGETHER    *
!     WITH THE ROUTINES GPFA, GPFA2F, GPFA3F, GPFA5F                 *
!                                                                    *
!*********************************************************************

!        SUBROUTINE 'SETGPFA'
!        SETUP ROUTINE FOR SELF-SORTING IN-PLACE
!            GENERALIZED PRIME FACTOR (COMPLEX) FFT [GPFA]

!        CALL SETGPFA(TRIGS,N)

!        INPUT :
!        -----
!        N IS THE LENGTH OF THE TRANSFORMS. N MUST BE OF THE FORM:
!          -----------------------------------
!            N = (2**IP) * (3**IQ) * (5**IR)
!          -----------------------------------

!        OUTPUT:
!        ------
!        TRIGS IS A TABLE OF TWIDDLE FACTORS,
!          OF LENGTH 2*IPQR (REAL) WORDS, WHERE:
!          --------------------------------------
!            IPQR = (2**IP) + (3**IQ) + (5**IR)
!          --------------------------------------

!        WRITTEN BY CLIVE TEMPERTON 1990

!----------------------------------------------------------------------
** */

// Decompose n into factors 2,3,5
	unsigned int nj[3];
	bool bRes = PrimeDecompose(n, nj);

	if (!bRes)
	{
		printf(" *** WARNING!!!%10d IS NOT A LEGAL VALUE OF N ***\n", n);
		return false;
	}

// Compute list of rotated twiddle factors
	nj[0] = powint(2, nj[0]);
	nj[1] = powint(3, nj[1]);
	nj[2] = powint(5, nj[2]);

	real twopi = Pi + Pi;
	int i = 0;

	for(int ll=0; ll<3; ++ll)
	{
		unsigned int ni = nj[ll];
        if (ni == 1)
			continue;

		real del = twopi / ni;
		int irot = n / ni;
		int kink = irot % ni;
        unsigned int kk = 0;
        for(unsigned int k=0; k<ni; ++k)
		{
			real angle = kk * del;
			trigs[i].set(Cos(angle), Sin(angle));
			++i;
			kk += kink;
			if (kk >= ni) 
				kk -= ni;
		}
	}
	return true;
}

//
// Self-sorting in-place generalized prime factor (complex) FFT
bool FftEngineGpfaft::Gpfa(Complex *ab, Complex *trigs, int inc, int jump, unsigned int n, int lot, FftDirection isign)
{
/* **
       SUBROUTINE 'GPFA'
       SELF-SORTING IN-PLACE GENERALIZED PRIME FACTOR (COMPLEX) FFT

       *** THIS IS THE ALL-FORTRAN VERSION ***
           -------------------------------

       CALL GPFA(A,B,TRIGS,INC,JUMP,N,LOT,ISIGN)

       A IS FIRST REAL INPUT/OUTPUT VECTOR
       B IS FIRST IMAGINARY INPUT/OUTPUT VECTOR
       TRIGS IS A TABLE OF TWIDDLE FACTORS, PRECALCULATED
             BY CALLING SUBROUTINE 'SETGPFA'
       INC IS THE INCREMENT WITHIN EACH DATA VECTOR
       JUMP IS THE INCREMENT BETWEEN DATA VECTORS
       N IS THE LENGTH OF THE TRANSFORMS:
         -----------------------------------
           N = (2**IP) * (3**IQ) * (5**IR)
         -----------------------------------
       LOT IS THE NUMBER OF TRANSFORMS
       ISIGN = +1 FOR FORWARD TRANSFORM
             = -1 FOR INVERSE TRANSFORM

       WRITTEN BY CLIVE TEMPERTON
       RECHERCHE EN PREVISION NUMERIQUE
       ATMOSPHERIC ENVIRONMENT SERVICE, CANADA

       DEFINITION OF TRANSFORM
       -----------------------

       X(J) = SUM(K=0,...,N-1)(C(K)*EXP(ISIGN*2*I*J*K*PI/N))

       FOR A MATHEMATICAL DEVELOPMENT OF THE ALGORITHM USED,
       SEE:
        C TEMPERTON : "A GENERALIZED PRIME FACTOR FFT ALGORITHM
         FOR ANY N = (2**P)(3**Q)(5**R)",
         SIAM J. SCI. STAT. COMP., MAY 1992.
** */

//     DECOMPOSE N INTO FACTORS 2,3,5
//     ------------------------------

	unsigned int nj[3];
	bool bRes = PrimeDecompose(n, nj);

	if (!bRes)
	{
		printf(" *** WARNING!!!%10d IS NOT A LEGAL VALUE OF N ***\n", n);
		return false;
	}
//
// Compute the transform
	int i = 0;
	if (nj[0] > 0)
	{
		Gpfa2f(ab, trigs, inc, jump, n, nj[0], lot, isign);
		i += powint(2, nj[0]);
	}

	if (nj[1] > 0)
	{
		Gpfa3f(ab, trigs+i, inc, jump, n, nj[1], lot, isign);
		i += powint(3, nj[1]);
	}

	if (nj[2] > 0)
	{
		Gpfa5f(ab, trigs+i, inc, jump, n, nj[2], lot, isign);
	}
	return true;
}

void FftEngineGpfaft::Gpfa2f(Complex *ab, Complex *trigs, int inc, int jump, unsigned int n, unsigned int mm, int lot, FftDirection isign)
{
/* **
    fortran version of *gpfa2* -
    radix-2 section of self-sorting, in-place, generalized pfa
    central radix-2 and radix-8 passes included
    so that transform length can be any power of 2

    ***************************************************************
    *                                                             *
    *  N.B. LVR = LENGTH OF VECTOR REGISTERS, SET TO 128 FOR C90. *
    *  RESET TO 64 FOR OTHER CRAY MACHINES, OR TO ANY LARGE VALUE *
    *  (GREATER THAN OR EQUAL TO LOT) FOR A SCALAR COMPUTER.      *
    *                                                             *
    ***************************************************************
** */

	const real half_ = (real)0.5;
	const real onex_ = (real)1.;
	const int lvr = 64;

	unsigned int n2 = powint(2, mm);
	unsigned int inq = n / n2;
	int jstepx = (n2 - n) * inc;
	int ninc = n * inc;
	int ink = inc * inq;

	int m = 0;
	bool m2 = false;
	bool m8 = false;
	switch(mm % 4)
	{
	case 1:
		m = (mm - 1) / 2;
		m2 = true;
		break;

	case 3:
		m = (mm - 3) / 2;
		m8 = true;
		break;

	default:
		m = mm / 2;
		break;
	}
//
	int mh = (m + 1) / 2;
	int nblox = 1 + (lot - 1) / lvr;
	int left = lot;
	real s = (isign == FftForward) ? onex_ : -onex_;
	int istart = 0;

	Complex cs[7];
    int ipass, j, ja, jb, jc, jd, je, jf, jg, jh, ji, jj, jk, jl, jm, jn, jo, jp, jstep, jstepl;
    int k, l, la, mu, nvex, ll;
    unsigned int nu, jjj;
	real c1, c2, c3;
    unsigned int kk = 0;
//
// Loop on blocks of lvr transforms
	for(int nb=0; nb<nblox; ++nb)
	{
		if (left <= lvr)
			nvex = left;
		else
		{
			if (left < (2*lvr))
			{
				nvex = left / 2;
				nvex += nvex % 2;
			}
			else
				nvex = lvr;
		}
        left -= nvex;

        la = 1;
//
// Loop on type I radix-4 passes
		mu = inq % 4;
        if (isign == FftBackward) 
			mu = 4 - mu;
        real ss = onex_;
        if (mu == 3) 
			ss = -onex_;

        if (mh == 0) 
			goto l200;

		real t0, t1, t2, t3, u0, u1, u2, u3;
		for(ipass=0; ipass<mh; ++ipass)
		{
			jstep = (n * inc) / (4 * la);
			jstepl = jstep - ninc;
//
// k = 0 loop (no twiddle factors)
			for(jjj=0; jjj<=(n-1)*inc; jjj+=4*jstep)
			{
				ja = istart + jjj;
                for(nu=0; nu<inq; ++nu)					// "transverse" loop
				{
					II(jb, ja, jstepl)
					II(jc, jb, jstepl)
					II(jd, jc, jstepl)
					j = 0;
//
// Loop across transforms
// dir$ ivdep, shortloop
					for(l=0; l<nvex; ++l)
					{
						t0 =  ab[ja+j].re + ab[jc+j].re;
						t2 =  ab[ja+j].re - ab[jc+j].re;
						t1 =  ab[jb+j].re + ab[jd+j].re;
						t3 = (ab[jb+j].re - ab[jd+j].re) * ss;
						u0 =  ab[ja+j].im + ab[jc+j].im;
						u2 =  ab[ja+j].im - ab[jc+j].im;
						u1 =  ab[jb+j].im + ab[jd+j].im;
						u3 = (ab[jb+j].im - ab[jd+j].im) * ss;
						ab[ja+j].re = t0 + t1;
						ab[jc+j].re = t0 - t1;
						ab[ja+j].im = u0 + u1;
						ab[jc+j].im = u0 - u1;
						ab[jb+j].re = t2 - u3;
						ab[jd+j].re = t2 + u3;
						ab[jb+j].im = u2 + t3;
						ab[jd+j].im = u2 - t3;
						j += jump;
					}
					II(ja, ja, jstepx)
				}
			}
//
// Finished if n2 = 4
			if (n2 == 4)
				goto l490;
			kk = la;
//
// Loop on nonzero k
			for(k=ink; k<=jstep-ink; k+=ink)
			{
				LoadTrigs(cs, 3, trigs, kk, s);
//
// Loop along transform
				for(jjj=k; jjj<=(n-1)*inc; jjj+=4*jstep)
				{
					ja = istart + jjj;
					for(nu=0; nu<inq; ++nu)					// "transverse" loop
					{
						II(jb, ja, jstepl)
						II(jc, jb, jstepl)
						II(jd, jc, jstepl)
						j = 0;
//
// Loop across transforms
// dir$ ivdep,shortloop
						for(l=0; l<nvex; ++l)
						{
							t0 =  ab[ja+j].re + ab[jc+j].re;
							t2 =  ab[ja+j].re - ab[jc+j].re;
							t1 =  ab[jb+j].re + ab[jd+j].re;
							t3 = (ab[jb+j].re - ab[jd+j].re) * ss;
							u0 =  ab[ja+j].im + ab[jc+j].im;
							u2 =  ab[ja+j].im - ab[jc+j].im;
							u1 =  ab[jb+j].im + ab[jd+j].im;
							u3 = (ab[jb+j].im - ab[jd+j].im) * ss;
							ab[ja+j].re = t0 + t1;
							ab[ja+j].im = u0 + u1;
							ab[jb+j].re = cs[0].re * (t2-u3) - cs[0].im * (u2+t3);
							ab[jb+j].im = cs[0].im * (t2-u3) + cs[0].re * (u2+t3);
							ab[jc+j].re = cs[1].re * (t0-t1) - cs[1].im * (u0-u1);
							ab[jc+j].im = cs[1].im * (t0-t1) + cs[1].re * (u0-u1);
							ab[jd+j].re = cs[2].re * (t2+u3) - cs[2].im * (u2-t3);
							ab[jd+j].im = cs[2].im * (t2+u3) + cs[2].re * (u2-t3);
							j += jump;
						}								// -----( end of loop across transforms )
						II(ja, ja, jstepx)
					}
				}										// -----( end of loop along transforms )
				kk += la;
			}											// -----( end of loop on nonzero k )
			la *= 4;
		}												// -----( end of loop on type I radix-4 passes)
//
// Central radix-2 pass
l200:
		if (m2 == false) 
			goto l300;

		jstep = (n * inc) / (2 * la);
		jstepl = jstep - ninc;
//
// k=0 loop (no twiddle factors)
		for(jjj=0; jjj<=(n-1)*inc; jjj+=2*jstep)
		{
			ja = istart + jjj;
			for(nu=0; nu<inq; ++nu)						// "transverse" loop
			{
				II(jb, ja, jstepl)
				j = 0;
//
// Loop across transforms
// dir$ ivdep, shortloop
				for(l=0; l<nvex; ++l)
				{
					t0 = ab[ja+j].re - ab[jb+j].re;
					ab[ja+j].re += ab[jb+j].re;
					ab[jb+j].re = t0;
					u0 = ab[ja+j].im - ab[jb+j].im;
					ab[ja+j].im += ab[jb+j].im;
					ab[jb+j].im = u0;
					j += jump;
				}										// -----(end of loop across transforms)
				II(ja, ja, jstepx)
			}
		}
//
// Finished if n2=2
		if (n2 == 2) 
			goto l490;
		kk = la;
//
// Loop on nonzero k
		for(k=ink; k<=jstep-ink; k+=ink)
		{
			cs[0] = trigs[kk].MultIm(s);
//
// Loop along transforms
			for(jjj=k; jjj<=(n-1)*inc; jjj+=2*jstep)
			{
				ja = istart + jjj;
				for(nu=0; nu<inq; ++nu)					// "transverse" loop
				{
					II(jb, ja, jstepl)
					j = 0;
//
// Loop across transforms
					if (2*kk == n2/2)
					{
// dir$ ivdep, shortloop
						for(l=0; l<nvex; ++l)
						{
							t0          = (ab[ja+j].re - ab[jb+j].re) * ss;
							ab[ja+j].re += ab[jb+j].re;
							ab[jb+j].re = (ab[jb+j].im - ab[ja+j].im) * ss;
							ab[ja+j].im += ab[jb+j].im;
							ab[jb+j].im = t0;
							j += jump;
						}
					}
					else
					{
// dir$ ivdep, shortloop
						for(l=0; l<nvex; ++l)
						{
							t0 = ab[ja+j].re - ab[jb+j].re;
							ab[ja+j].re += ab[jb+j].re;
							u0 = ab[ja+j].im - ab[jb+j].im;
							ab[ja+j].im += ab[jb+j].im;
							ab[jb+j].re = cs[0].re * t0 - cs[0].im * u0;
							ab[jb+j].im = cs[0].im * t0 + cs[0].re * u0;
							j += jump;
						}
					}									// -----(end of loop across transforms)
					II(ja, ja, jstepx)
				}
			}											// -----(end of loop along transforms)
			kk += la;
		}												// -----(end of loop on nonzero k)
//
// -----(end of radix-2 pass)
		la = 2*la;
        goto l400;
//
// Central radix-8 pass
l300:
		if (m8 == false) 
			goto l400;
		jstep = (n * inc) / (8 * la);
		jstepl = jstep - ninc;
        mu = inq % 8;
        if (isign == FftBackward) 
			mu = 8 - mu;
        c1 = onex_;
        if ((mu == 3) || (mu == 7))
			c1 = -onex_;
        c2 = Sqrt(half_);
        if ((mu == 3) || (mu == 5))
			c2 = -c2;
        c3 = c1 * c2;
//
// Stage 1
		for(k=0; k<=jstep-ink; k+=ink)
		{
			for(jjj=k; jjj<=(n-1)*inc; jjj+=8*jstep)
			{
				ja = istart + jjj;
				for(nu=0; nu<inq; ++nu)					// "transverse" loop
				{
					II(jb, ja, jstepl)
					II(jc, jb, jstepl)
					II(jd, jc, jstepl)
					II(je, jd, jstepl)
					II(jf, je, jstepl)
					II(jg, jf, jstepl)
					II(jh, jg, jstepl)
					j = 0;
// dir$ ivdep, shortloop
					for(l=0; l<nvex; ++l)
					{
						t0 = ab[ja+j].re - ab[je+j].re;
						ab[ja+j].re += ab[je+j].re;
						t1 = c1*(ab[jc+j].re - ab[jg+j].re);
						ab[je+j].re = ab[jc+j].re + ab[jg+j].re;
						t2 = ab[jb+j].re - ab[jf+j].re;
						ab[jc+j].re = ab[jb+j].re + ab[jf+j].re;
						t3 = ab[jd+j].re - ab[jh+j].re;
						ab[jg+j].re = ab[jd+j].re + ab[jh+j].re;
						ab[jb+j].re = t0;
						ab[jf+j].re = t1;
						ab[jd+j].re = c2 * (t2-t3);
						ab[jh+j].re = c3 * (t2+t3);
						u0 = ab[ja+j].im - ab[je+j].im;
						ab[ja+j].im += ab[je+j].im;
						u1 = c1*(ab[jc+j].im - ab[jg+j].im);
						ab[je+j].im = ab[jc+j].im + ab[jg+j].im;
						u2 = ab[jb+j].im - ab[jf+j].im;
						ab[jc+j].im = ab[jb+j].im + ab[jf+j].im;
						u3 = ab[jd+j].im - ab[jh+j].im;
						ab[jg+j].im = ab[jd+j].im + ab[jh+j].im;
						ab[jb+j].im = u0;
						ab[jf+j].im = u1;
						ab[jd+j].im = c2*(u2 - u3);
						ab[jh+j].im = c3*(u2 + u3);
						j += jump;
					}
					II(ja, ja, jstepx)
				}
			}
		}
//
// Stage 2
//
// k=0 (no twiddle factors)
		for(jjj=0; jjj<=(n-1)*inc; jjj+=8*jstep)
		{
			ja = istart + jjj;
			for(nu=0; nu<inq; ++nu)					// "transverse" loop
			{
				II(jb, ja, jstepl)
				II(jc, jb, jstepl)
				II(jd, jc, jstepl)
				II(je, jd, jstepl)
				II(jf, je, jstepl)
				II(jg, jf, jstepl)
				II(jh, jg, jstepl)
				j = 0;
// dir$ ivdep, shortloop
				for(l=0; l<nvex; ++l)
				{
					t0 =  ab[ja+j].re + ab[je+j].re;
					t2 =  ab[ja+j].re - ab[je+j].re;
					t1 =  ab[jc+j].re + ab[jg+j].re;
					t3 = (ab[jc+j].re - ab[jg+j].re) * c1;
					u0 =  ab[ja+j].im + ab[je+j].im;
					u2 =  ab[ja+j].im - ab[je+j].im;
					u1 =  ab[jc+j].im + ab[jg+j].im;
					u3 = (ab[jc+j].im - ab[jg+j].im) * c1;
					ab[ja+j].re = t0 + t1;
					ab[je+j].re = t0 - t1;
					ab[ja+j].im = u0 + u1;
					ab[je+j].im = u0 - u1;
					ab[jc+j].re = t2 - u3;
					ab[jg+j].re = t2 + u3;
					ab[jc+j].im = u2 + t3;
					ab[jg+j].im = u2 - t3;
					t0 = ab[jb+j].re + ab[jd+j].re;
					t2 = ab[jb+j].re - ab[jd+j].re;
					t1 = ab[jf+j].re - ab[jh+j].re;
					t3 = ab[jf+j].re + ab[jh+j].re;
					u0 = ab[jb+j].im + ab[jd+j].im;
					u2 = ab[jb+j].im - ab[jd+j].im;
					u1 = ab[jf+j].im - ab[jh+j].im;
					u3 = ab[jf+j].im + ab[jh+j].im;
					ab[jb+j].re = t0 - u3;
					ab[jh+j].re = t0 + u3;
					ab[jb+j].im = u0 + t3;
					ab[jh+j].im = u0 - t3;
					ab[jd+j].re = t2 + u1;
					ab[jf+j].re = t2 - u1;
					ab[jd+j].im = u2 - t1;
					ab[jf+j].im = u2 + t1;
					j += jump;
				}
				II(ja, ja, jstepx)
			}
		}
		if (n2 == 8) 
			goto l490;
//
// Loop on nonzero k
		kk = la;
		for(k=ink; k<=jstep-ink; k+=ink)
		{
			LoadTrigs(cs, 7, trigs, kk, s);
//
			for(jjj=k; jjj<=(n-1)*inc; jjj+=8*jstep)
			{
				ja = istart + jjj;
				for(nu=0; nu<inq; ++nu)					// "transverse" loop
				{
					II(jb, ja, jstepl)
					II(jc, jb, jstepl)
					II(jd, jc, jstepl)
					II(je, jd, jstepl)
					II(jf, je, jstepl)
					II(jg, jf, jstepl)
					II(jh, jg, jstepl)
					j = 0;
// dir$ ivdep, shortloop
					for(l=0; l<nvex; ++l)
					{
						t0 =  ab[ja+j].re + ab[je+j].re;
						t2 =  ab[ja+j].re - ab[je+j].re;
						t1 =  ab[jc+j].re + ab[jg+j].re;
						t3 = (ab[jc+j].re - ab[jg+j].re) * c1;
						u0 =  ab[ja+j].im + ab[je+j].im;
						u2 =  ab[ja+j].im - ab[je+j].im;
						u1 =  ab[jc+j].im + ab[jg+j].im;
						u3 = (ab[jc+j].im - ab[jg+j].im) * c1;
						ab[ja+j].re = t0 + t1;
						ab[ja+j].im = u0 + u1;
						ab[je+j].re = cs[3].re * (t0-t1) - cs[3].im * (u0-u1);
						ab[je+j].im = cs[3].im * (t0-t1) + cs[3].re * (u0-u1);
						ab[jc+j].re = cs[1].re * (t2-u3) - cs[1].im * (u2+t3);
						ab[jc+j].im = cs[1].im * (t2-u3) + cs[1].re * (u2+t3);
						ab[jg+j].re = cs[5].re * (t2+u3) - cs[5].im * (u2-t3);
						ab[jg+j].im = cs[5].im * (t2+u3) + cs[5].re * (u2-t3);
						t0 = ab[jb+j].re + ab[jd+j].re;
						t2 = ab[jb+j].re - ab[jd+j].re;
						t1 = ab[jf+j].re - ab[jh+j].re;
						t3 = ab[jf+j].re + ab[jh+j].re;
						u0 = ab[jb+j].im + ab[jd+j].im;
						u2 = ab[jb+j].im - ab[jd+j].im;
						u1 = ab[jf+j].im - ab[jh+j].im;
						u3 = ab[jf+j].im + ab[jh+j].im;
						ab[jb+j].re = cs[0].re * (t0-u3) - cs[0].im * (u0+t3);
						ab[jb+j].im = cs[0].im * (t0-u3) + cs[0].re * (u0+t3);
						ab[jh+j].re = cs[6].re * (t0+u3) - cs[6].im * (u0-t3);
						ab[jh+j].im = cs[6].im * (t0+u3) + cs[6].re * (u0-t3);
						ab[jd+j].re = cs[2].re * (t2+u1) - cs[2].im * (u2-t1);
						ab[jd+j].im = cs[2].im * (t2+u1) + cs[2].re * (u2-t1);
						ab[jf+j].re = cs[4].re * (t2-u1) - cs[4].im * (u2+t1);
						ab[jf+j].im = cs[4].im * (t2-u1) + cs[4].re * (u2+t1);
						j += jump;
					}
					II(ja, ja, jstepx)
				}
			}
			kk += la;
		}

		la = 8*la;
//
// Loop on type II radix-4 passes
l400:
		int laincl;
		mu = inq % 4;
        if (isign == FftBackward) 
			mu = 4 - mu;
        ss = onex_;
        if (mu == 3) ss = -onex_;

		for(ipass=mh; ipass<m; ++ipass)
		{
			jstep = (n*inc)/(4*la);
			jstepl = jstep - ninc;
			laincl = la*ink - ninc;
//
// k=0 loop (no twiddle factors)
			for(ll=0; ll<=(la-1)*ink; ll+=4*jstep)
			{
				for(jjj=ll; jjj<=(n-1)*inc; jjj+=4*la*ink)
				{
					ja = istart + jjj;
					for(nu=0; nu<inq; ++nu)					// "transverse" loop
					{
						II(jb, ja, jstepl)
						II(jc, jb, jstepl)
						II(jd, jc, jstepl)
						II(je, ja, laincl)
						II(jf, je, jstepl)
						II(jg, jf, jstepl)
						II(jh, jg, jstepl)
						II(ji, je, laincl)
						II(jj, ji, jstepl)
						II(jk, jj, jstepl)
						II(jl, jk, jstepl)
						II(jm, ji, laincl)
						II(jn, jm, jstepl)
						II(jo, jn, jstepl)
						II(jp, jo, jstepl)
						j = 0;
//
// Loop across transforms
// dir$ ivdep, shortloop
						for(l=0; l<nvex; ++l)
						{
							t0 =  ab[ja+j].re + ab[jc+j].re;
							t2 =  ab[ja+j].re - ab[jc+j].re;
							t1 =  ab[jb+j].re + ab[jd+j].re;
							t3 = (ab[jb+j].re - ab[jd+j].re) * ss;
							ab[jc+j].re = ab[ji+j].re;
							u0 =  ab[ja+j].im + ab[jc+j].im;
							u2 =  ab[ja+j].im - ab[jc+j].im;
							u1 =  ab[jb+j].im + ab[jd+j].im;
							u3 = (ab[jb+j].im - ab[jd+j].im) * ss;
							ab[jb+j].re = ab[je+j].re;
							ab[ja+j].re = t0 + t1;
							ab[ji+j].re = t0 - t1;
							ab[ja+j].im = u0 + u1;
							ab[jc+j].im = u0 - u1;
							ab[jd+j].im = ab[jm+j].im;
							ab[je+j].re = t2 - u3;
							ab[jd+j].re = t2 + u3;
							ab[jb+j].im = u2 + t3;
							ab[jm+j].im = u2 - t3;
//
							t0 =  ab[jb+j].re + ab[jg+j].re;
							t2 =  ab[jb+j].re - ab[jg+j].re;
							t1 =  ab[jf+j].re + ab[jh+j].re;
							t3 = (ab[jf+j].re - ab[jh+j].re) * ss;
							ab[jg+j].re = ab[jj+j].re;
							u0 =  ab[je+j].im + ab[jg+j].im;
							u2 =  ab[je+j].im - ab[jg+j].im;
							u1 =  ab[jf+j].im + ab[jh+j].im;
							u3 = (ab[jf+j].im - ab[jh+j].im) * ss;
							ab[je+j].im = ab[jb+j].im;
							ab[jb+j].re = t0 + t1;
							ab[jj+j].re = t0 - t1;
							ab[jg+j].im = ab[jj+j].im;
							ab[jb+j].im = u0 + u1;
							ab[jj+j].im = u0 - u1;
							ab[jf+j].re = t2 - u3;
							ab[jh+j].re = t2 + u3;
							ab[jf+j].im = u2 + t3;
							ab[jh+j].im = u2 - t3;
//
							t0 =  ab[jc+j].re + ab[jk+j].re;
							t2 =  ab[jc+j].re - ab[jk+j].re;
							t1 =  ab[jg+j].re + ab[jl+j].re;
							t3 = (ab[jg+j].re - ab[jl+j].re) * ss;
							u0 =  ab[ji+j].im + ab[jk+j].im;
							u2 =  ab[ji+j].im - ab[jk+j].im;
							ab[jl+j].re = ab[jo+j].re;
							u1 =  ab[jg+j].im + ab[jl+j].im;
							u3 = (ab[jg+j].im - ab[jl+j].im) * ss;
							ab[ji+j].im = ab[jc+j].im;
							ab[jc+j].re = t0 + t1;
							ab[jk+j].re = t0 - t1;
							ab[jl+j].im = ab[jo+j].im;
							ab[jc+j].im = u0 + u1;
							ab[jk+j].im = u0 - u1;
							ab[jg+j].re = t2 - u3;
							ab[jo+j].re = t2 + u3;
							ab[jg+j].im = u2 + t3;
							ab[jo+j].im = u2 - t3;
//
							t0 =  ab[jm+j].re + ab[jl+j].re;
							t2 =  ab[jm+j].re - ab[jl+j].re;
							t1 =  ab[jn+j].re + ab[jp+j].re;
							t3 = (ab[jn+j].re - ab[jp+j].re) * ss;
							ab[jm+j].re = ab[jd+j].re;
							u0 =  ab[jd+j].im + ab[jl+j].im;
							u2 =  ab[jd+j].im - ab[jl+j].im;
							u1 =  ab[jn+j].im + ab[jp+j].im;
							u3 = (ab[jn+j].im - ab[jp+j].im) * ss;
							ab[jn+j].re = ab[jh+j].re;
							ab[jd+j].re = t0 + t1;
							ab[jl+j].re = t0 - t1;
							ab[jd+j].im = u0 + u1;
							ab[jl+j].im = u0 - u1;
							ab[jn+j].im = ab[jh+j].im;
							ab[jh+j].re = t2 - u3;
							ab[jp+j].re = t2 + u3;
							ab[jh+j].im = u2 + t3;
							ab[jp+j].im = u2 - t3;
							j += jump;
						}								// -----( end of loop across transforms )
						II(ja, ja, jstepx)
					}
				}
			}											// -----( end of double loop for k=0 )
//
// Finished if last pass
			if (ipass == m-1) 
				goto l490;

			kk = la;
//
// Loop on nonzero k
			for(k=ink; k<=jstep-ink; k+=ink)
			{
				LoadTrigs(cs, 3, trigs, kk, s);
//
// Double loop along first transform in block
				for(ll=k; ll<=(la-1)*ink; ll+=4*jstep)
				{
					for(jjj=ll; jjj<=(n-1)*inc; jjj+=4*la*ink)
					{
						ja = istart + jjj;
						for(nu=0; nu<inq; ++nu)					// "transverse" loop
						{
							II(jb, ja, jstepl)
							II(jc, jb, jstepl)
							II(jd, jc, jstepl)
							II(je, ja, laincl)
							II(jf, je, jstepl)
							II(jg, jf, jstepl)
							II(jh, jg, jstepl)
							II(ji, je, laincl)
							II(jj, ji, jstepl)
							II(jk, jj, jstepl)
							II(jl, jk, jstepl)
							II(jm, ji, laincl)
							II(jn, jm, jstepl)
							II(jo, jn, jstepl)
							II(jp, jo, jstepl)
							j = 0;
//
// Loop across transforms
// dir$ ivdep, shortloop
							for(l=0; l<nvex; ++l)
							{
								t0 =  ab[ja+j].re + ab[jc+j].re;
								t2 =  ab[ja+j].re - ab[jc+j].re;
								t1 =  ab[jb+j].re + ab[jd+j].re;
								t3 = (ab[jb+j].re - ab[jd+j].re) * ss;
								ab[jc+j].re = ab[ji+j].re;
								u0 =  ab[ja+j].im + ab[jc+j].im;
								u2 =  ab[ja+j].im - ab[jc+j].im;
								u1 =  ab[jb+j].im + ab[jd+j].im;
								u3 = (ab[jb+j].im - ab[jd+j].im) * ss;
								ab[jb+j].re = ab[je+j].re;
								ab[ja+j].re = t0 + t1;
								ab[ja+j].im = u0 + u1;
								ab[je+j].re = cs[0].re * (t2-u3) - cs[0].im * (u2+t3);
								ab[jb+j].im = cs[0].im * (t2-u3) + cs[0].re * (u2+t3);
								ab[jd+j].im = ab[jm+j].im;
								ab[ji+j].re = cs[1].re * (t0-t1) - cs[1].im * (u0-u1);
								ab[jc+j].im = cs[1].im * (t0-t1) + cs[1].re * (u0-u1);
								ab[jd+j].re = cs[2].re * (t2+u3) - cs[2].im * (u2-t3);
								ab[jm+j].im = cs[2].im * (t2+u3) + cs[2].re * (u2-t3);
//
								t0 =  ab[jb+j].re + ab[jg+j].re;
								t2 =  ab[jb+j].re - ab[jg+j].re;
								t1 =  ab[jf+j].re + ab[jh+j].re;
								t3 = (ab[jf+j].re - ab[jh+j].re) * ss;
								ab[jg+j].re = ab[jj+j].re;
								u0 =  ab[je+j].im + ab[jg+j].im;
								u2 =  ab[je+j].im - ab[jg+j].im;
								u1 =  ab[jf+j].im + ab[jh+j].im;
								u3 = (ab[jf+j].im - ab[jh+j].im) * ss;
								ab[je+j].im = ab[jb+j].im;
								ab[jb+j].re = t0 + t1;
								ab[jb+j].im = u0 + u1;
								ab[jg+j].im = ab[jj+j].im;
								ab[jf+j].re = cs[0].re * (t2-u3) - cs[0].im * (u2+t3);
								ab[jf+j].im = cs[0].im * (t2-u3) + cs[0].re * (u2+t3);
								ab[jj+j].re = cs[1].re * (t0-t1) - cs[1].im * (u0-u1);
								ab[jj+j].im = cs[1].im * (t0-t1) + cs[1].re * (u0-u1);
								ab[jh+j].re = cs[2].re * (t2+u3) - cs[2].im * (u2-t3);
								ab[jh+j].im = cs[2].im * (t2+u3) + cs[2].re * (u2-t3);
//
								t0 =  ab[jc+j].re + ab[jk+j].re;
								t2 =  ab[jc+j].re - ab[jk+j].re;
								t1 =  ab[jg+j].re + ab[jl+j].re;
								t3 = (ab[jg+j].re - ab[jl+j].re) * ss;
								u0 =  ab[ji+j].im + ab[jk+j].im;
								u2 =  ab[ji+j].im - ab[jk+j].im;
								ab[jl+j].re = ab[jo+j].re;
								u1 =  ab[jg+j].im + ab[jl+j].im;
								u3 = (ab[jg+j].im - ab[jl+j].im) * ss;
								ab[ji+j].im = ab[jc+j].im;
								ab[jc+j].re = t0 + t1;
								ab[jc+j].im = u0 + u1;
								ab[jl+j].im = ab[jo+j].im;
								ab[jg+j].re = cs[0].re * (t2-u3) - cs[0].im * (u2+t3);
								ab[jg+j].im = cs[0].im * (t2-u3) + cs[0].re * (u2+t3);
								ab[jk+j].re = cs[1].re * (t0-t1) - cs[1].im * (u0-u1);
								ab[jk+j].im = cs[1].im * (t0-t1) + cs[1].re * (u0-u1);
								ab[jo+j].re = cs[2].re * (t2+u3) - cs[2].im * (u2-t3);
								ab[jo+j].im = cs[2].im * (t2+u3) + cs[2].re * (u2-t3);
//
								t0 =  ab[jm+j].re + ab[jl+j].re;
								t2 =  ab[jm+j].re - ab[jl+j].re;
								t1 =  ab[jn+j].re + ab[jp+j].re;
								t3 = (ab[jn+j].re - ab[jp+j].re) * ss;
								ab[jm+j].re = ab[jd+j].re;
								u0 = ab[jd+j].im + ab[jl+j].im;
								u2 = ab[jd+j].im - ab[jl+j].im;
								ab[jn+j].re = ab[jh+j].re;
								u1 =  ab[jn+j].im + ab[jp+j].im;
								u3 = (ab[jn+j].im - ab[jp+j].im) * ss;
								ab[jn+j].im = ab[jh+j].im;
								ab[jd+j].re = t0 + t1;
								ab[jd+j].im = u0 + u1;
								ab[jh+j].re = cs[0].re * (t2-u3) - cs[0].im * (u2+t3);
								ab[jh+j].im = cs[0].im * (t2-u3) + cs[0].re * (u2+t3);
								ab[jl+j].re = cs[1].re * (t0-t1) - cs[1].im * (u0-u1);
								ab[jl+j].im = cs[1].im * (t0-t1) + cs[1].re * (u0-u1);
								ab[jp+j].re = cs[2].re * (t2+u3) - cs[2].im * (u2-t3);
								ab[jp+j].im = cs[2].im * (t2+u3) + cs[2].re * (u2-t3);
								j += jump;
							}							// -----(end of loop across transforms)
							II(ja, ja, jstepx)
						}
					}
				}										// -----( end of double loop for this k )
				kk += la;
			}											// -----( end of loop over values of k )
			la = 4*la;
		}												// -----( end of loop on type II radix-4 passes )
//
// -----( nvex transforms completed)
l490:
		istart += nvex*jump;
	}													// -----( end of loop on blocks of transforms )
}

void FftEngineGpfaft::Gpfa3f(Complex *ab, Complex *trigs, int inc, int jump, unsigned int n, unsigned int mm, int lot, FftDirection isign)
{
/* **
!     fortran version of *gpfa3* -
!     radix-3 section of self-sorting, in-place generalized PFA
!
!     ***************************************************************
!     *                                                             *
!     *  N.B. LVR = LENGTH OF VECTOR REGISTERS, SET TO 128 FOR C90. *
!     *  RESET TO 64 FOR OTHER CRAY MACHINES, OR TO ANY LARGE VALUE *
!     *  (GREATER THAN OR EQUAL TO LOT) FOR A SCALAR COMPUTER.      *
!     *                                                             *
!     ***************************************************************
** */

	const real half_ = (real)0.5;
	const real onex_ = (real)1.;
	const real sin60 = (real)0.866025403784437;
	const int lvr = 64;

//	FILE *file15 = fopen("ccc.txt", "w");

	unsigned int n3 = powint(3, mm);
	unsigned int inq = n / n3;
	int jstepx = (n3 - n) * inc;
	int ninc = n * inc;
	int ink = inc * inq;
	int mu = inq % 3;
	if (isign == FftBackward) 
		mu = 3 - mu;
	int m = mm;
	int mh = (m + 1) / 2;
	real s = (isign == FftForward) ? onex_ : -onex_;
	real c1 = sin60;
	if (mu == 2)
		c1 = -c1;

	int nblox = 1 + (lot - 1) / lvr;
	int left = lot;
	s = (isign == FftForward) ? onex_ : -onex_;
	int istart = 0;

	Complex cs[2];
	real t1, t2, t3, u1, u2, u3;
    int ipass, j, ja, jb, jc, jd, je, jf, jg, jh, ji, jstep, jstepl, laincl;
    int k, l, la, ll, nvex;
    unsigned int nu, jjj;
    unsigned int kk = 0;
//
// Loop on blocks of lvr transforms
	for(int nb=0; nb<nblox; ++nb)
	{
		if (left <= lvr)
			nvex = left;
        else
		{
			if (left < (2 * lvr))
			{
				nvex = left / 2;
				nvex += nvex % 2;
			}
			else
				nvex = lvr;
		}
        left -= nvex;
		la = 1;
//
// Loop on type I radix-3 passes
        for(ipass=0; ipass<mh; ++ipass)
		{
			jstep = (n*inc)/(3*la);
			jstepl = jstep - ninc;
//
// k = 0 loop (no twiddle factors)
			for(jjj=0; jjj<=(n-1)*inc; jjj+=3*jstep)
			{
				ja = istart + jjj;
				for(nu=0; nu<inq; ++nu)					// "transverse" loop
				{
					II(jb, ja, jstepl)
					II(jc, jb, jstepl)
					j = 0;
//
// Loop across transforms
// dir$ ivdep, shortloop
					for(l=0; l<nvex; ++l)
					{
						t1 =  ab[jb+j].re + ab[jc+j].re;
						t2 =  ab[ja+j].re - half_*t1;
						t3 = (ab[jb+j].re - ab[jc+j].re) * c1;
						u1 =  ab[jb+j].im + ab[jc+j].im;
						u2 =  ab[ja+j].im - half_*u1;
						u3 = (ab[jb+j].im - ab[jc+j].im) * c1;
						ab[ja+j].re = ab[ja+j].re + t1;
						ab[ja+j].im = ab[ja+j].im + u1;
						ab[jb+j].re = t2 - u3;
						ab[jb+j].im = u2 + t3;
						ab[jc+j].re = t2 + u3;
						ab[jc+j].im = u2 - t3;
						j += jump;
					}
					II(ja, ja, jstepx)
				}
			}
//
// Finished if n3 = 3
			if (n3 == 3) 
				goto l490;
			kk = la;
//
// Loop on nonzero k
			for(k=ink; k<=jstep-ink; k+=ink)
			{
				LoadTrigs(cs, 2, trigs, kk, s);
//
// Loop along transform
				for(jjj=k; jjj<=(n-1)*inc; jjj+=3*jstep)
				{
					ja = istart + jjj;
					for(nu=0; nu<inq; ++nu)					// "transverse" loop
					{
						II(jb, ja, jstepl)
						II(jc, jb, jstepl)
						j = 0;
// 
// Loop across transforms
// dir$ ivdep,shortloop
						for(l=0; l<nvex; ++l)
						{
							t1 =  ab[jb+j].re + ab[jc+j].re;
							t2 =  ab[ja+j].re - half_*t1;
							t3 = (ab[jb+j].re - ab[jc+j].re) * c1;
							u1 =  ab[jb+j].im + ab[jc+j].im;
							u2 =  ab[ja+j].im - half_*u1;
							u3 = (ab[jb+j].im - ab[jc+j].im) * c1;
							ab[ja+j].re = ab[ja+j].re + t1;
							ab[ja+j].im = ab[ja+j].im + u1;
							ab[jb+j].re = cs[0].re * (t2-u3) - cs[0].im * (u2+t3);
							ab[jb+j].im = cs[0].im * (t2-u3) + cs[0].re * (u2+t3);
							ab[jc+j].re = cs[1].re * (t2+u3) - cs[1].im * (u2-t3);
							ab[jc+j].im = cs[1].im * (t2+u3) + cs[1].re * (u2-t3);
							j += jump;
						}								// -----( end of loop across transforms )
						II(ja, ja, jstepx)
					}
				}										// -----( end of loop along transforms )
				kk += la;
			}											// -----( end of loop on nonzero k )
			la = 3*la;
		}												// -----( end of loop on type I radix-3 passes)
//
// Loop on type II radix-3 passes
		for(ipass=mh; ipass<m; ++ipass)
		{
			jstep = (n*inc)/(3*la);
			jstepl = jstep - ninc;
			laincl = la*ink - ninc;
//
// k=0 loop (no twiddle factors)
			for(ll=0; ll<=(la-1)*ink; ll+=3*jstep)
			{
				for(jjj=ll; jjj<=(n-1)*inc; jjj+=3*la*ink)
				{
					ja = istart + jjj;
					for(nu=0; nu<inq; ++nu)					// "transverse" loop
					{
						II(jb, ja, jstepl)
						II(jc, jb, jstepl)
						II(jd, ja, laincl)
						II(je, jd, jstepl)
						II(jf, je, jstepl)
						II(jg, jd, laincl)
						II(jh, jg, jstepl)
						II(ji, jh, jstepl)
						j = 0;
//
// Loop across transforms
// dir$ ivdep, shortloop
						for(l=0; l<nvex; ++l)
						{
							t1 =  ab[jb+j].re + ab[jc+j].re;
							t2 =  ab[ja+j].re - half_*t1;
							t3 = (ab[jb+j].re - ab[jc+j].re) * c1;
							ab[jb+j].re = ab[jd+j].re;
							u1 =  ab[jb+j].im + ab[jc+j].im;
							u2 =  ab[ja+j].im - half_*u1;
							u3 = (ab[jb+j].im - ab[jc+j].im) * c1;
							ab[jb+j].im = ab[jd+j].im;
							ab[ja+j].re = ab[ja+j].re + t1;
							ab[ja+j].im = ab[ja+j].im + u1;
							ab[jd+j].re = t2 - u3;
							ab[jd+j].im = u2 + t3;
							ab[jc+j].re = t2 + u3;
							ab[jc+j].im = u2 - t3;
//
							t1 =  ab[je+j].re + ab[jf+j].re;
							t2 =  ab[jb+j].re - half_*t1;
							t3 = (ab[je+j].re - ab[jf+j].re) * c1;
							ab[jf+j].re = ab[jh+j].re;
							u1 =  ab[je+j].im + ab[jf+j].im;
							u2 =  ab[jb+j].im - half_*u1;
							u3 = (ab[je+j].im - ab[jf+j].im) * c1;
							ab[jf+j].im = ab[jh+j].im;
							ab[jb+j].re = ab[jb+j].re + t1;
							ab[jb+j].im = ab[jb+j].im + u1;
							ab[je+j].re = t2 - u3;
							ab[je+j].im = u2 + t3;
							ab[jh+j].re = t2 + u3;
							ab[jh+j].im = u2 - t3;
//
							t1 =  ab[jf+j].re + ab[ji+j].re;
							t2 =  ab[jg+j].re - half_*t1;
							t3 = (ab[jf+j].re - ab[ji+j].re) * c1;
							t1 =  ab[jg+j].re + t1;
							ab[jg+j].re = ab[jc+j].re;
							u1 =  ab[jf+j].im + ab[ji+j].im;
							u2 =  ab[jg+j].im - half_*u1;
							u3 = (ab[jf+j].im - ab[ji+j].im) * c1;
							u1 =  ab[jg+j].im + u1;
							ab[jg+j].im = ab[jc+j].im;
							ab[jc+j].re = t1;
							ab[jc+j].im = u1;
							ab[jf+j].re = t2 - u3;
							ab[jf+j].im = u2 + t3;
							ab[ji+j].re = t2 + u3;
							ab[ji+j].im = u2 - t3;
							j += jump;
						}								// -----( end of loop across transforms )
						II(ja, ja, jstepx)
					}
				}
			}											// -----( end of double loop for k=0 )
//
// Finished if last pass
			if (ipass == m-1) 
				goto l490;
			kk = la;
//
// Loop on nonzero k
			for(k=ink; k<=jstep-ink; k+=ink)
			{
				LoadTrigs(cs, 2, trigs, kk, s);
//
// Double loop along first transform in block
				for(ll=k; ll<=(la-1)*ink; ll+=3*jstep)
				{
					for(jjj=ll; jjj<=(n-1)*inc; jjj+=3*la*ink)
					{
						ja = istart + jjj;
						for(nu=0; nu<inq; ++nu)					// "transverse" loop
						{
							II(jb, ja, jstepl)
							II(jc, jb, jstepl)
							II(jd, ja, laincl)
							II(je, jd, jstepl)
							II(jf, je, jstepl)
							II(jg, jd, laincl)
							II(jh, jg, jstepl)
							II(ji, jh, jstepl)
							j = 0;
//
// Loop across transforms
// dir$ ivdep, shortloop
							for(l=0; l<nvex; ++l)
							{
								t1 =  ab[jb+j].re + ab[jc+j].re;
								t2 =  ab[ja+j].re - half_*t1;
								t3 = (ab[jb+j].re - ab[jc+j].re) * c1;
								ab[jb+j].re = ab[jd+j].re;
								u1 =  ab[jb+j].im + ab[jc+j].im;
								u2 =  ab[ja+j].im - half_*u1;
								u3 = (ab[jb+j].im - ab[jc+j].im) * c1;
								ab[jb+j].im = ab[jd+j].im;
								ab[ja+j].re = ab[ja+j].re + t1;
								ab[ja+j].im = ab[ja+j].im + u1;
								ab[jd+j].re = cs[0].re * (t2-u3) - cs[0].im * (u2+t3);
								ab[jd+j].im = cs[0].im * (t2-u3) + cs[0].re * (u2+t3);
								ab[jc+j].re = cs[1].re * (t2+u3) - cs[1].im * (u2-t3);
								ab[jc+j].im = cs[1].im * (t2+u3) + cs[1].re * (u2-t3);
//
								t1 =  ab[je+j].re + ab[jf+j].re;
								t2 =  ab[jb+j].re - half_*t1;
								t3 = (ab[je+j].re - ab[jf+j].re) * c1;
								ab[jf+j].re = ab[jh+j].re;
								u1 =  ab[je+j].im + ab[jf+j].im;
								u2 =  ab[jb+j].im - half_*u1;
								u3 = (ab[je+j].im - ab[jf+j].im) * c1;
								ab[jf+j].im = ab[jh+j].im;
								ab[jb+j].re = ab[jb+j].re + t1;
								ab[jb+j].im = ab[jb+j].im + u1;
								ab[je+j].re = cs[0].re * (t2-u3) - cs[0].im * (u2+t3);
								ab[je+j].im = cs[0].im * (t2-u3) + cs[0].re * (u2+t3);
								ab[jh+j].re = cs[1].re * (t2+u3) - cs[1].im * (u2-t3);
								ab[jh+j].im = cs[1].im * (t2+u3) + cs[1].re * (u2-t3);
//
								t1 =  ab[jf+j].re + ab[ji+j].re;
								t2 =  ab[jg+j].re - half_*t1;
								t3 = (ab[jf+j].re - ab[ji+j].re) * c1;
								t1 += ab[jg+j].re;
								ab[jg+j].re = ab[jc+j].re;
								u1 =  ab[jf+j].im + ab[ji+j].im;
								u2 =  ab[jg+j].im - half_*u1;
								u3 = (ab[jf+j].im - ab[ji+j].im) * c1;
								u1 += ab[jg+j].im;
								ab[jg+j].im = ab[jc+j].im;
								ab[jc+j].re = t1;
								ab[jc+j].im = u1;
								ab[jf+j].re = cs[0].re * (t2-u3) - cs[0].im * (u2+t3);
								ab[jf+j].im = cs[0].im * (t2-u3) + cs[0].re * (u2+t3);
								ab[ji+j].re = cs[1].re * (t2+u3) - cs[1].im * (u2-t3);
								ab[ji+j].im = cs[1].im * (t2+u3) + cs[1].re * (u2-t3);
								j += jump;
							}							// -----(end of loop across transforms)
							II(ja, ja, jstepx)
						}
					}
				}										// -----( end of double loop for this k )
				kk += la;
			}											// -----( end of loop over values of k )
			la = 3*la;
		}												// -----( end of loop on type II radix-3 passes )
//
// -----( nvex transforms completed)
l490:
		istart += nvex*jump;
	}													// -----( end of loop on blocks of transforms )
//	fclose(file15);
}

void FftEngineGpfaft::Gpfa5f(Complex *ab, Complex *trigs, int inc, int jump, unsigned int n, unsigned int mm, int lot, FftDirection isign)
{
/* **
Fortran version of *gpfa5* - radix-5 section of self-sorting, in-place, generalized pfa

!     ***************************************************************
!     *                                                             *
!     *  N.B. LVR = LENGTH OF VECTOR REGISTERS, SET TO 128 FOR C90. *
!     *  RESET TO 64 FOR OTHER CRAY MACHINES, OR TO ANY LARGE VALUE *
!     *  (GREATER THAN OR EQUAL TO LOT) FOR A SCALAR COMPUTER.      *
!     *                                                             *
!     ***************************************************************
** */

	const int lvr = 64;
	const real sin36 = (real)0.587785252292473;
	const real sin72 = (real)0.951056516295154;
	const real qrt5  = (real)0.559016994374947;
	const real quat_ = (real)0.25;
	const real onex_ = (real)1.;

	unsigned int n5 = powint(5, mm);
	unsigned int inq = n / n5;
	int jstepx = (n5 - n) * inc;
	int ninc = n * inc;
	int ink = inc *  inq;
	int mu = inq % 5;
	if (isign == FftBackward) 
		mu = 5 - mu;

	int m = mm;
	int mh = (m + 1) / 2;
	real s = (isign == FftForward) ? onex_ : -onex_;
	real c1 = qrt5;
	real c2 = sin72;
	real c3 = sin36;
	if ((mu == 2) || (mu == 3))
	{
		c1 = -c1;
		c2 = sin36;
		c3 = sin72;
	}
	if ((mu == 3) || (mu == 4))
		c2 = -c2;
	if ((mu == 2) || (mu == 4))
		c3 = -c3;

	int nblox = 1 + (lot-1)/lvr;
	int left = lot;
	s = (isign == FftForward) ? onex_ : -onex_;
	int istart = 0;

	Complex cs[4];
	real t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11;
	real u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, ax, bx;
    int ipass, j, ja, jb, jc, jd, je, jf, jg, jh, ji, jj, jk, jl, jm, jn, jo, jp, jq, jr, js, jt, ju, jv, jw, jx, jy, jstep, jstepl, k, l, la, nvex, laincl, ll;
    unsigned int nu, jjj;
    unsigned int kk = 0;
//
// Loop on blocks of lvr transforms
	for(int nb=0; nb<nblox; ++nb)
	{
		if (left <= lvr)
			nvex = left;
        else
		{
			if (left < (2*lvr))
			{
				nvex = left/2;
				nvex += nvex%2;
			}
			else
				nvex = lvr;
		}
        left -= nvex;

		la = 1;
//
// Loop on type I radix-5 passes
		for(ipass=0; ipass<mh; ++ipass)
		{
			jstep = (n*inc)/(5*la);
			jstepl = jstep - ninc;
			kk = 0;
			for(k=0; k<=jstep-ink; k+=ink)
			{
				if (k > 0)
					LoadTrigs(cs, 4, trigs, kk, s);
//
// Loop along transform
				for(jjj=k; jjj<=(n-1)*inc; jjj+=5*jstep)
				{
					ja = istart + jjj;
					for(nu=0; nu<inq; ++nu)					// "transverse" loop
					{
                        II(jb, ja, jstepl)
						II(jc, jb, jstepl)
						II(jd, jc, jstepl)
						II(je, jd, jstepl)
						j = 0;
//
// Loop across transforms
						if (k == 0)
						{
// dir$ ivdep, shortloop
							for(l=0; l<nvex; ++l)
							{
								t1 = ab[jb+j].re + ab[je+j].re;
								t2 = ab[jc+j].re + ab[jd+j].re;
								t3 = ab[jb+j].re - ab[je+j].re;
								t4 = ab[jc+j].re - ab[jd+j].re;
								t5 = t1 + t2;
								t6 = c1 * (t1-t2);
								t7 = ab[ja+j].re - quat_*t5;
								ab[ja+j].re += t5;
								t8 = t7 + t6;
								t9 = t7 - t6;
								t10 = c3*t3 - c2*t4;
								t11 = c2*t3 + c3*t4;
								u1 = ab[jb+j].im + ab[je+j].im;
								u2 = ab[jc+j].im + ab[jd+j].im;
								u3 = ab[jb+j].im - ab[je+j].im;
								u4 = ab[jc+j].im - ab[jd+j].im;
								u5 = u1 + u2;
								u6 = c1 * (u1-u2);
								u7 = ab[ja+j].im - quat_*u5;
								ab[ja+j].im += u5;
								u8 = u7 + u6;
								u9 = u7 - u6;
								u10 = c3*u3 - c2*u4;
								u11 = c2*u3 + c3*u4;
								ab[jb+j].re = t8 - u11;
								ab[jb+j].im = u8 + t11;
								ab[je+j].re = t8 + u11;
								ab[je+j].im = u8 - t11;
								ab[jc+j].re = t9 - u10;
								ab[jc+j].im = u9 + t10;
								ab[jd+j].re = t9 + u10;
								ab[jd+j].im = u9 - t10;
								j += jump;
							}
						}
						else
						{
// dir$ ivdep,shortloop
							for(l=0; l<nvex; ++l)
							{
								t1 = ab[jb+j].re + ab[je+j].re;
								t2 = ab[jc+j].re + ab[jd+j].re;
								t3 = ab[jb+j].re - ab[je+j].re;
								t4 = ab[jc+j].re - ab[jd+j].re;
								t5 = t1 + t2;
								t6 = c1 * (t1-t2);
								t7 = ab[ja+j].re - quat_*t5;
								ab[ja+j].re += t5;
								t8 = t7 + t6;
								t9 = t7 - t6;
								t10 = c3*t3 - c2*t4;
								t11 = c2*t3 + c3*t4;
								u1 = ab[jb+j].im + ab[je+j].im;
								u2 = ab[jc+j].im + ab[jd+j].im;
								u3 = ab[jb+j].im - ab[je+j].im;
								u4 = ab[jc+j].im - ab[jd+j].im;
								u5 = u1 + u2;
								u6 = c1 * (u1-u2);
								u7 = ab[ja+j].im - quat_*u5;
								ab[ja+j].im += u5;
								u8 = u7 + u6;
								u9 = u7 - u6;
								u10 = c3*u3 - c2*u4;
								u11 = c2*u3 + c3*u4;
								ab[jb+j].re = cs[0].re * (t8-u11) - cs[0].im * (u8+t11);
								ab[jb+j].im = cs[0].im * (t8-u11) + cs[0].re * (u8+t11);
								ab[je+j].re = cs[3].re * (t8+u11) - cs[3].im * (u8-t11);
								ab[je+j].im = cs[3].im * (t8+u11) + cs[3].re * (u8-t11);
								ab[jc+j].re = cs[1].re * (t9-u10) - cs[1].im * (u9+t10);
								ab[jc+j].im = cs[1].im * (t9-u10) + cs[1].re * (u9+t10);
								ab[jd+j].re = cs[2].re * (t9+u10) - cs[2].im * (u9-t10);
								ab[jd+j].im = cs[2].im * (t9+u10) + cs[2].re * (u9-t10);
								j += jump;
							}
						}								// -----( end of loop across transforms )
						II(ja, ja, jstepx)
					}
				}										// -----( end of loop along transforms )
				kk += la;
			}											// -----( end of loop on nonzero k )
			la = 5 * la;
		}												// -----( end of loop on type I radix-5 passes)

		if (n != 5)
		{
//
// Loop on type II radix-5 passes
			for(ipass=mh; ipass<m; ++ipass)
			{
				jstep = (n*inc)/(5*la);
				jstepl = jstep - ninc;
				laincl = la*ink - ninc;
				kk = 0;
				for(k=0; k<=jstep-ink; k+=ink)
				{
					if (k > 0)
						LoadTrigs(cs, 4, trigs, kk, s);
//
// Double loop along first transform in block
					for(ll=k; ll<=(la-1)*ink; ll+=5*jstep)
					{
						for(jjj=ll; jjj<=(n-1)*inc; jjj+=5*la*ink)
						{
							ja = istart + jjj;
							for(nu=0; nu<inq; ++nu)					// "transverse" loop
							{
								II(jb, ja, jstepl)
								II(jc, jb, jstepl)
								II(jd, jc, jstepl)
								II(je, jd, jstepl)
								II(jf, ja, laincl)
								II(jg, jf, jstepl)
								II(jh, jg, jstepl)
								II(ji, jh, jstepl)
								II(jj, ji, jstepl)
								II(jk, jf, laincl)
								II(jl, jk, jstepl)
								II(jm, jl, jstepl)
								II(jn, jm, jstepl)
								II(jo, jn, jstepl)
								II(jp, jk, laincl)
								II(jq, jp, jstepl)
								II(jr, jq, jstepl)
								II(js, jr, jstepl)
								II(jt, js, jstepl)
								II(ju, jp, laincl)
								II(jv, ju, jstepl)
								II(jw, jv, jstepl)
								II(jx, jw, jstepl)
								II(jy, jx, jstepl)
								j = 0;
//
// Loop across transforms
								if (k == 0)
								{
// dir$ ivdep, shortloop
									for(l=0; l<nvex; ++l)
									{
										t1 = ab[jb+j].re + ab[je+j].re;
										t2 = ab[jc+j].re + ab[jd+j].re;
										t3 = ab[jb+j].re - ab[je+j].re;
										t4 = ab[jc+j].re - ab[jd+j].re;
										ab[jb+j].re = ab[jf+j].re;
										t5 = t1 + t2;
										t6 = c1 * (t1 - t2);
										t7 = ab[ja+j].re - quat_*t5;
										ab[ja+j].re += t5;
										t8 = t7 + t6;
										t9 = t7 - t6;
										ab[jc+j].re = ab[jk+j].re;
										t10 = c3*t3 - c2*t4;
										t11 = c2*t3 + c3*t4;
										u1 = ab[jb+j].im + ab[je+j].im;
										u2 = ab[jc+j].im + ab[jd+j].im;
										u3 = ab[jb+j].im - ab[je+j].im;
										u4 = ab[jc+j].im - ab[jd+j].im;
										ab[jb+j].im = ab[jf+j].im;
										u5 = u1 + u2;
										u6 = c1 * (u1 - u2);
										u7 = ab[ja+j].im - quat_*u5;
										ab[ja+j].im += u5;
										u8 = u7 + u6;
										u9 = u7 - u6;
										ab[jc+j].im = ab[jk+j].im;
										u10 = c3*u3 - c2*u4;
										u11 = c2*u3 + c3*u4;
										ab[jf+j].re = t8 - u11;
										ab[jf+j].im = u8 + t11;
										ab[je+j].re = t8 + u11;
										ab[je+j].im = u8 - t11;
										ab[jk+j].re = t9 - u10;
										ab[jk+j].im = u9 + t10;
										ab[jd+j].re = t9 + u10;
										ab[jd+j].im = u9 - t10;
//
										t1 = ab[jg+j].re + ab[jj+j].re;
										t2 = ab[jh+j].re + ab[ji+j].re;
										t3 = ab[jg+j].re - ab[jj+j].re;
										t4 = ab[jh+j].re - ab[ji+j].re;
										ab[jh+j].re = ab[jl+j].re;
										t5 = t1 + t2;
										t6 = c1 * (t1-t2);
										t7 = ab[jb+j].re - quat_*t5;
										ab[jb+j].re += t5;
										t8 = t7 + t6;
										t9 = t7 - t6;
										ab[ji+j].re = ab[jq+j].re;
										t10 = c3*t3 - c2*t4;
										t11 = c2*t3 + c3*t4;
										u1 = ab[jg+j].im + ab[jj+j].im;
										u2 = ab[jh+j].im + ab[ji+j].im;
										u3 = ab[jg+j].im - ab[jj+j].im;
										u4 = ab[jh+j].im - ab[ji+j].im;
										ab[jh+j].im = ab[jl+j].im;
										u5 = u1 + u2;
										u6 = c1 * (u1-u2);
										u7 = ab[jb+j].im - quat_*u5;
										ab[jb+j].im += u5;
										u8 = u7 + u6;
										u9 = u7 - u6;
										ab[ji+j].im = ab[jq+j].im;
										u10 = c3*u3 - c2*u4;
										u11 = c2*u3 + c3*u4;
										ab[jg+j].re = t8 - u11;
										ab[jg+j].im = u8 + t11;
										ab[jj+j].re = t8 + u11;
										ab[jj+j].im = u8 - t11;
										ab[jl+j].re = t9 - u10;
										ab[jl+j].im = u9 + t10;
										ab[jq+j].re = t9 + u10;
										ab[jq+j].im = u9 - t10;
//
										t1 = ab[jh+j].re + ab[jo+j].re;
										t2 = ab[jm+j].re + ab[jn+j].re;
										t3 = ab[jh+j].re - ab[jo+j].re;
										t4 = ab[jm+j].re - ab[jn+j].re;
										ab[jn+j].re = ab[jr+j].re;
										t5 = t1 + t2;
										t6 = c1 * (t1-t2);
										t7 = ab[jc+j].re - quat_*t5;
										ab[jc+j].re += t5;
										t8 = t7 + t6;
										t9 = t7 - t6;
										ab[jo+j].re = ab[jw+j].re;
										t10 = c3*t3 - c2*t4;
										t11 = c2*t3 + c3*t4;
										u1 = ab[jh+j].im + ab[jo+j].im;
										u2 = ab[jm+j].im + ab[jn+j].im;
										u3 = ab[jh+j].im - ab[jo+j].im;
										u4 = ab[jm+j].im - ab[jn+j].im;
										ab[jn+j].im = ab[jr+j].im;
										u5 = u1 + u2;
										u6 = c1 * (u1-u2);
										u7 = ab[jc+j].im - quat_*u5;
										ab[jc+j].im += u5;
										u8 = u7 + u6;
										u9 = u7 - u6;
										ab[jo+j].im = ab[jw+j].im;
										u10 = c3*u3 - c2*u4;
										u11 = c2*u3 + c3*u4;
										ab[jh+j].re = t8 - u11;
										ab[jh+j].im = u8 + t11;
										ab[jw+j].re = t8 + u11;
										ab[jw+j].im = u8 - t11;
										ab[jm+j].re = t9 - u10;
										ab[jm+j].im = u9 + t10;
										ab[jr+j].re = t9 + u10;
										ab[jr+j].im = u9 - t10;
//
										t1 = ab[ji+j].re + ab[jt+j].re;
										t2 = ab[jn+j].re + ab[js+j].re;
										t3 = ab[ji+j].re - ab[jt+j].re;
										t4 = ab[jn+j].re - ab[js+j].re;
										ab[jt+j].re = ab[jx+j].re;
										t5 = t1 + t2;
										t6 = c1 * (t1-t2);
										t7 = ab[jp+j].re - quat_*t5;
										ax = ab[jp+j].re + t5;
										t8 = t7 + t6;
										t9 = t7 - t6;
										ab[jp+j].re = ab[jd+j].re;
										t10 = c3*t3 - c2*t4;
										t11 = c2*t3 + c3*t4;
										ab[jd+j].re = ax;
										u1 = ab[ji+j].im + ab[jt+j].im;
										u2 = ab[jn+j].im + ab[js+j].im;
										u3 = ab[ji+j].im - ab[jt+j].im;
										u4 = ab[jn+j].im - ab[js+j].im;
										ab[jt+j].im = ab[jx+j].im;
										u5 = u1 + u2;
										u6 = c1 * (u1-u2);
										u7 = ab[jp+j].im - quat_*u5;
										bx = ab[jp+j].im + u5;
										u8 = u7 + u6;
										u9 = u7 - u6;
										ab[jp+j].im = ab[jd+j].im;
										u10 = c3*u3 - c2*u4;
										u11 = c2*u3 + c3*u4;
										ab[jd+j].im = bx;
										ab[ji+j].re = t8 - u11;
										ab[ji+j].im = u8 + t11;
										ab[jx+j].re = t8 + u11;
										ab[jx+j].im = u8 - t11;
										ab[jn+j].re = t9 - u10;
										ab[jn+j].im = u9 + t10;
										ab[js+j].re = t9 + u10;
										ab[js+j].im = u9 - t10;
//
										t1 = ab[jv+j].re + ab[jy+j].re;
										t2 = ab[jo+j].re + ab[jt+j].re;
										t3 = ab[jv+j].re - ab[jy+j].re;
										t4 = ab[jo+j].re - ab[jt+j].re;
										ab[jv+j].re = ab[jj+j].re;
										t5 = t1 + t2;
										t6 = c1 * (t1-t2);
										t7 = ab[ju+j].re - quat_*t5;
										ax = ab[ju+j].re + t5;
										t8 = t7 + t6;
										t9 = t7 - t6;
										ab[ju+j].re = ab[je+j].re;
										t10 = c3*t3 - c2*t4;
										t11 = c2*t3 + c3*t4;
										ab[je+j].re = ax;
										u1 = ab[jv+j].im + ab[jy+j].im;
										u2 = ab[jo+j].im + ab[jt+j].im;
										u3 = ab[jv+j].im - ab[jy+j].im;
										u4 = ab[jo+j].im - ab[jt+j].im;
										ab[jv+j].im = ab[jj+j].im;
										u5 = u1 + u2;
										u6 = c1 * (u1-u2);
										u7 = ab[ju+j].im - quat_*u5;
										bx = ab[ju+j].im + u5;
										u8 = u7 + u6;
										u9 = u7 - u6;
										ab[ju+j].im = ab[je+j].im;
										u10 = c3*u3 - c2*u4;
										u11 = c2*u3 + c3*u4;
										ab[je+j].im = bx;
										ab[jj+j].re = t8 - u11;
										ab[jj+j].im = u8 + t11;
										ab[jy+j].re = t8 + u11;
										ab[jy+j].im = u8 - t11;
										ab[jo+j].re = t9 - u10;
										ab[jo+j].im = u9 + t10;
										ab[jt+j].re = t9 + u10;
										ab[jt+j].im = u9 - t10;
										j += jump;
									}
								}
								else
								{
// dir$ ivdep, shortloop
									for(l=0; l<nvex; ++l)
									{
										t1 = ab[jb+j].re + ab[je+j].re;
										t2 = ab[jc+j].re + ab[jd+j].re;
										t3 = ab[jb+j].re - ab[je+j].re;
										t4 = ab[jc+j].re - ab[jd+j].re;
										ab[jb+j].re = ab[jf+j].re;
										t5 = t1 + t2;
										t6 = c1 * (t1-t2);
										t7 = ab[ja+j].re - quat_*t5;
										ab[ja+j].re += t5;
										t8 = t7 + t6;
										t9 = t7 - t6;
										ab[jc+j].re = ab[jk+j].re;
										t10 = c3*t3 - c2*t4;
										t11 = c2*t3 + c3*t4;
										u1 = ab[jb+j].im + ab[je+j].im;
										u2 = ab[jc+j].im + ab[jd+j].im;
										u3 = ab[jb+j].im - ab[je+j].im;
										u4 = ab[jc+j].im - ab[jd+j].im;
										ab[jb+j].im = ab[jf+j].im;
										u5 = u1 + u2;
										u6 = c1 * (u1-u2);
										u7 = ab[ja+j].im - quat_*u5;
										ab[ja+j].im += u5;
										u8 = u7 + u6;
										u9 = u7 - u6;
										ab[jc+j].im = ab[jk+j].im;
										u10 = c3*u3 - c2*u4;
										u11 = c2*u3 + c3*u4;
										ab[jf+j].re = cs[0].re * (t8-u11) - cs[0].im * (u8+t11);
										ab[jf+j].im = cs[0].im * (t8-u11) + cs[0].re * (u8+t11);
										ab[je+j].re = cs[3].re * (t8+u11) - cs[3].im * (u8-t11);
										ab[je+j].im = cs[3].im * (t8+u11) + cs[3].re * (u8-t11);
										ab[jk+j].re = cs[1].re * (t9-u10) - cs[1].im * (u9+t10);
										ab[jk+j].im = cs[1].im * (t9-u10) + cs[1].re * (u9+t10);
										ab[jd+j].re = cs[2].re * (t9+u10) - cs[2].im * (u9-t10);
										ab[jd+j].im = cs[2].im * (t9+u10) + cs[2].re * (u9-t10);
//
										t1 = ab[jg+j].re + ab[jj+j].re;
										t2 = ab[jh+j].re + ab[ji+j].re;
										t3 = ab[jg+j].re - ab[jj+j].re;
										t4 = ab[jh+j].re - ab[ji+j].re;
										ab[jh+j].re = ab[jl+j].re;
										t5 = t1 + t2;
										t6 = c1 * (t1-t2);
										t7 = ab[jb+j].re - quat_*t5;
										ab[jb+j].re += t5;
										t8 = t7 + t6;
										t9 = t7 - t6;
										ab[ji+j].re = ab[jq+j].re;
										t10 = c3*t3 - c2*t4;
										t11 = c2*t3 + c3*t4;
										u1 = ab[jg+j].im + ab[jj+j].im;
										u2 = ab[jh+j].im + ab[ji+j].im;
										u3 = ab[jg+j].im - ab[jj+j].im;
										u4 = ab[jh+j].im - ab[ji+j].im;
										ab[jh+j].im = ab[jl+j].im;
										u5 = u1 + u2;
										u6 = c1 * (u1-u2);
										u7 = ab[jb+j].im - quat_*u5;
										ab[jb+j].im += u5;
										u8 = u7 + u6;
										u9 = u7 - u6;
										ab[ji+j].im = ab[jq+j].im;
										u10 = c3*u3 - c2*u4;
										u11 = c2*u3 + c3*u4;
										ab[jg+j].re = cs[0].re * (t8-u11) - cs[0].im * (u8+t11);
										ab[jg+j].im = cs[0].im * (t8-u11) + cs[0].re * (u8+t11);
										ab[jj+j].re = cs[3].re * (t8+u11) - cs[3].im * (u8-t11);
										ab[jj+j].im = cs[3].im * (t8+u11) + cs[3].re * (u8-t11);
										ab[jl+j].re = cs[1].re * (t9-u10) - cs[1].im * (u9+t10);
										ab[jl+j].im = cs[1].im * (t9-u10) + cs[1].re * (u9+t10);
										ab[jq+j].re = cs[2].re * (t9+u10) - cs[2].im * (u9-t10);
										ab[jq+j].im = cs[2].im * (t9+u10) + cs[2].re * (u9-t10);
//
										t1 = ab[jh+j].re + ab[jo+j].re;
										t2 = ab[jm+j].re + ab[jn+j].re;
										t3 = ab[jh+j].re - ab[jo+j].re;
										t4 = ab[jm+j].re - ab[jn+j].re;
										ab[jn+j].re = ab[jr+j].re;
										t5 = t1 + t2;
										t6 = c1 * (t1-t2);
										t7 = ab[jc+j].re - quat_*t5;
										ab[jc+j].re += t5;
										t8 = t7 + t6;
										t9 = t7 - t6;
										ab[jo+j].re = ab[jw+j].re;
										t10 = c3*t3 - c2*t4;
										t11 = c2*t3 + c3*t4;
										u1 = ab[jh+j].im + ab[jo+j].im;
										u2 = ab[jm+j].im + ab[jn+j].im;
										u3 = ab[jh+j].im - ab[jo+j].im;
										u4 = ab[jm+j].im - ab[jn+j].im;
										ab[jn+j].im = ab[jr+j].im;
										u5 = u1 + u2;
										u6 = c1 * (u1-u2);
										u7 = ab[jc+j].im - quat_*u5;
										ab[jc+j].im += u5;
										u8 = u7 + u6;
										u9 = u7 - u6;
										ab[jo+j].im = ab[jw+j].im;
										u10 = c3*u3 - c2*u4;
										u11 = c2*u3 + c3*u4;
										ab[jh+j].re = cs[0].re * (t8-u11) - cs[0].im * (u8+t11);
										ab[jh+j].im = cs[0].im * (t8-u11) + cs[0].re * (u8+t11);
										ab[jw+j].re = cs[3].re * (t8+u11) - cs[3].im * (u8-t11);
										ab[jw+j].im = cs[3].im * (t8+u11) + cs[3].re * (u8-t11);
										ab[jm+j].re = cs[1].re * (t9-u10) - cs[1].im * (u9+t10);
										ab[jm+j].im = cs[1].im * (t9-u10) + cs[1].re * (u9+t10);
										ab[jr+j].re = cs[2].re * (t9+u10) - cs[2].im * (u9-t10);
										ab[jr+j].im = cs[2].im * (t9+u10) + cs[2].re * (u9-t10);
//
										t1 = ab[ji+j].re + ab[jt+j].re;
										t2 = ab[jn+j].re + ab[js+j].re;
										t3 = ab[ji+j].re - ab[jt+j].re;
										t4 = ab[jn+j].re - ab[js+j].re;
										ab[jt+j].re = ab[jx+j].re;
										t5 = t1 + t2;
										t6 = c1 * (t1-t2);
										t7 = ab[jp+j].re - quat_*t5;
										ax = ab[jp+j].re + t5;
										t8 = t7 + t6;
										t9 = t7 - t6;
										ab[jp+j].re = ab[jd+j].re;
										t10 = c3*t3 - c2*t4;
										t11 = c2*t3 + c3*t4;
										ab[jd+j].re = ax;
										u1 = ab[ji+j].im + ab[jt+j].im;
										u2 = ab[jn+j].im + ab[js+j].im;
										u3 = ab[ji+j].im - ab[jt+j].im;
										u4 = ab[jn+j].im - ab[js+j].im;
										ab[jt+j].im = ab[jx+j].im;
										u5 = u1 + u2;
										u6 = c1 * (u1-u2);
										u7 = ab[jp+j].im - quat_*u5;
										bx = ab[jp+j].im + u5;
										u8 = u7 + u6;
										u9 = u7 - u6;
										ab[jp+j].im = ab[jd+j].im;
										u10 = c3*u3 - c2*u4;
										u11 = c2*u3 + c3*u4;
										ab[jd+j].im = bx;
										ab[ji+j].re = cs[0].re * (t8-u11) - cs[0].im * (u8+t11);
										ab[ji+j].im = cs[0].im * (t8-u11) + cs[0].re * (u8+t11);
										ab[jx+j].re = cs[3].re * (t8+u11) - cs[3].im * (u8-t11);
										ab[jx+j].im = cs[3].im * (t8+u11) + cs[3].re * (u8-t11);
										ab[jn+j].re = cs[1].re * (t9-u10) - cs[1].im * (u9+t10);
										ab[jn+j].im = cs[1].im * (t9-u10) + cs[1].re * (u9+t10);
										ab[js+j].re = cs[2].re * (t9+u10) - cs[2].im * (u9-t10);
										ab[js+j].im = cs[2].im * (t9+u10) + cs[2].re * (u9-t10);
//
										t1 = ab[jv+j].re + ab[jy+j].re;
										t2 = ab[jo+j].re + ab[jt+j].re;
										t3 = ab[jv+j].re - ab[jy+j].re;
										t4 = ab[jo+j].re - ab[jt+j].re;
										ab[jv+j].re = ab[jj+j].re;
										t5 = t1 + t2;
										t6 = c1 * (t1-t2);
										t7 = ab[ju+j].re - quat_*t5;
										ax = ab[ju+j].re + t5;
										t8 = t7 + t6;
										t9 = t7 - t6;
										ab[ju+j].re = ab[je+j].re;
										t10 = c3*t3 - c2*t4;
										t11 = c2*t3 + c3*t4;
										ab[je+j].re = ax;
										u1 = ab[jv+j].im + ab[jy+j].im;
										u2 = ab[jo+j].im + ab[jt+j].im;
										u3 = ab[jv+j].im - ab[jy+j].im;
										u4 = ab[jo+j].im - ab[jt+j].im;
										ab[jv+j].im = ab[jj+j].im;
										u5 = u1 + u2;
										u6 = c1*(u1-u2);
										u7 = ab[ju+j].im - quat_*u5;
										bx = ab[ju+j].im + u5;
										u8 = u7 + u6;
										u9 = u7 - u6;
										ab[ju+j].im = ab[je+j].im;
										u10 = c3*u3 - c2*u4;
										u11 = c2*u3 + c3*u4;
										ab[je+j].im = bx;
										ab[jj+j].re = cs[0].re * (t8-u11) - cs[0].im * (u8+t11);
										ab[jj+j].im = cs[0].im * (t8-u11) + cs[0].re * (u8+t11);
										ab[jy+j].re = cs[3].re * (t8+u11) - cs[3].im * (u8-t11);
										ab[jy+j].im = cs[3].im * (t8+u11) + cs[3].re * (u8-t11);
										ab[jo+j].re = cs[1].re * (t9-u10) - cs[1].im * (u9+t10);
										ab[jo+j].im = cs[1].im * (t9-u10) + cs[1].re * (u9+t10);
										ab[jt+j].re = cs[2].re * (t9+u10) - cs[2].im * (u9-t10);
										ab[jt+j].im = cs[2].im * (t9+u10) + cs[2].re * (u9-t10);
										j += jump;
									}							// -----(end of loop across transforms)
								}
								II(ja, ja, jstepx)
							}
						}
					}										// -----( end of double loop for this k )
					kk += la;
				}											// -----( end of loop over values of k )
			la = 5 * la;
			}												// -----( end of loop on type II radix-5 passes )
//
// -----( nvex transforms completed)
		}
		istart += nvex*jump;
	}													// -----( end of loop on blocks of transforms )
}

void FftEngineGpfaft::LoadTrigs(Complex *cs, unsigned int num, Complex *trigs, unsigned int kk, real s)
{
	unsigned int za, zb;
	for(za=0, zb=kk; za<num; ++za, zb+=kk)
	{
		cs[za] = trigs[zb].MultIm(s);
	}
}

void FftEngineGpfaft::DoFFT(Complex *ccc, unsigned int mx, unsigned int my, unsigned int mz, FftDirection isign)
{
/* **
Subroutine CXFFT3N

Interface routine for "Generalized Prime Factor Algorithm" FFT code
of Temperton for computation of 3 dimensional FFTs

Calling program should have call of form

   CALL CXFFT3N(CX,MX,MY,MZ,ISIGN)

where
   CX(I,J,K)= Complex array of 3-vectors at locations (I,J,K), with
              I = 1,...,MX
              J = 1,...,MY
              K = 1,...,MZ

Upon return,

           Mx-1 My-1 Mz-1
CX(I,J,K)= sum  sum  sum CXin(u,v,w)*exp[ISIGN*2*pi*i*(u*I/Mx+v*J/My+w*K/Mz)
           u=0  v=0  w=0

Note that here C is defined as a real array C(*)
C(1) = Re[CX(1,1,1)]
C(2) = Im[CX(1,1,1)]
C(3) = Re[CX(2,1,1)]
C(4) = Im[CX(2,1,1)]
...

or, in general:

C(1+2*((I-1)+(J-1)*MY+(K-1)*MY*MZ) = Re[CX(I,J,K)]
C(2+2*((I-1)+(J-1)*MY+(K-1)*MY*MZ) = Im[CX(I,J,K)]

Interface written by P.J.Flatau

History:
07.06.30 (BTD) Increased MXTRIG from 1000 to 16384=4*4096
               (since current version of EXTEND has NF235 up to 4096).
08.06.05 (BTD) corrected comments
end history
// ChB: ccc is Complex vector, then skips are 1, not 2
** */

	Init(mx, my, mz);
//
// First dimension
	int inc = 1;
	int jump = mx;
	int lot = my * mz;
	Gpfa(ccc, trigX, inc, jump, mx, lot, isign);
//
// Second dimension
	inc = mx;
	jump = 1;
	lot = mx;
    for(unsigned int k=0; k<mz; ++k)						// One plane at a time
	{
        Gpfa(ccc + k * mx * my, trigY, inc, jump, my, lot, isign);
//		Debugga(ccc, nx*ny*nz);
	}
//
// Third dimension
	inc = mx * my;
	jump = 1;
	lot = mx * my;
	Gpfa(ccc, trigZ, inc, jump, mz, lot, isign);
}
