#include "StdAfx.h"

#include "FftEngineGpfaftNew.h"

FftEngineGpfaftNew::FftEngineGpfaftNew(void)
{
	method = FftMethod_GPFAFT;
	nx = ny = nz = 0;
	trigX = trigY = trigZ = NULL;
	mxold = myold = mzold = 0;
}

FftEngineGpfaftNew::~FftEngineGpfaftNew(void)
{
	CleanDelete2(trigX);
	CleanDelete2(trigY);
	CleanDelete2(trigZ);
}

void FftEngineGpfaftNew::LoadTrigs(Complex *cs, int num, Complex *trigs, int kk, real s)
{
	int za, zb;
	for(za=0, zb=kk; za<num; ++za, zb+=kk)
	{
		cs[za] = trigs[zb].MultIm(s);
	}
}

void FftEngineGpfaftNew::Gpfa5f(Complex *ab, Complex *trigs, int inc, int jump, int n, int mm, int lot, FftDirection isign)
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
	const real onex_ = (real)1.;

	int n5 = powint(5, mm);
	int inq = n / n5;
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

	Complex cs[4], tu[12];
	int ipass, j, ja, jb, jc, jd, je, jf, jg, jh, ji, jj, jk, jl, jm, jn, jo, jp, jq, jr, js, jt, ju, jv, jw, jx, jy, jjj, jstep, jstepl, k, l, la, nb, nvex, nu, laincl, ll;
	int kk = 0;
//
// Loop on blocks of lvr transforms
	for(nb=0; nb<nblox; ++nb)
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
								Gpfa5fStep1(tu, ab, jb+j, je+j, jc+j, jd+j, c1);
								Gpfa5fStep2Simple(tu, ab, ja+j, c2, c3);
								Gpfa5fStep3A(tu, ab, jb+j, je+j, jc+j, jd+j);
								j += jump;
							}
						}
						else
						{
// dir$ ivdep,shortloop
							for(l=0; l<nvex; ++l)
							{
								Gpfa5fStep1(tu, ab, jb+j, je+j, jc+j, jd+j, c1);
								Gpfa5fStep2Simple(tu, ab, ja+j, c2, c3);
								Gpfa5fStep3B(tu, ab, jb+j, je+j, jc+j, jd+j, cs);
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
										Gpfa5fStep1(tu, ab, jb+j, je+j, jc+j, jd+j, c1);
										Gpfa5fStep2(tu, ab, jb+j, jf+j, jc+j, jk+j, ja+j, ja+j, c2, c3);
										Gpfa5fStep3A(tu, ab, jf+j, je+j, jk+j, jd+j);
//
										Gpfa5fStep1(tu, ab, jg+j, jj+j, jh+j, ji+j, c1);
										Gpfa5fStep2(tu, ab, jh+j, jl+j, ji+j, jq+j, jb+j, jb+j, c2, c3);
										Gpfa5fStep3A(tu, ab, jg+j, jj+j, jl+j, jq+j);
//
										Gpfa5fStep1(tu, ab, jh+j, jo+j, jm+j, jn+j, c1);
										Gpfa5fStep2(tu, ab, jn+j, jr+j, jo+j, jw+j, jc+j, jc+j, c2, c3);
										Gpfa5fStep3A(tu, ab, jh+j, jw+j, jm+j, jr+j);
//
										Gpfa5fStep1(tu, ab, ji+j, jt+j, jn+j, js+j, c1);
										Gpfa5fStep2(tu, ab, jt+j, jx+j, jp+j, jd+j, jp+j, jd+j, c2, c3);
										Gpfa5fStep3A(tu, ab, ji+j, jx+j, jn+j, js+j);
//
										Gpfa5fStep1(tu, ab, jv+j, jy+j, jo+j, jt+j, c1);
										Gpfa5fStep2(tu, ab, jv+j, jj+j, ju+j, je+j, ju+j, je+j, c2, c3);
										Gpfa5fStep3A(tu, ab, jj+j, jy+j, jo+j, jt+j);
										j += jump;
									}
								}
								else
								{
// dir$ ivdep, shortloop
									for(l=0; l<nvex; ++l)
									{
										Gpfa5fStep1(tu, ab, jb+j, je+j, jc+j, jd+j, c1);
										Gpfa5fStep2(tu, ab, jb+j, jf+j, jc+j, jk+j, ja+j, ja+j, c2, c3);
										Gpfa5fStep3B(tu, ab, jf+j, je+j, jk+j, jd+j, cs);
//
										Gpfa5fStep1(tu, ab, jg+j, jj+j, jh+j, ji+j, c1);
										Gpfa5fStep2(tu, ab, jh+j, jl+j, ji+j, jq+j, jb+j, jb+j, c2, c3);
										Gpfa5fStep3B(tu, ab, jg+j, jj+j, jl+j, jq+j, cs);
//
										Gpfa5fStep1(tu, ab, jh+j, jo+j, jm+j, jn+j, c1);
										Gpfa5fStep2(tu, ab, jn+j, jr+j, jo+j, jw+j, jc+j, jc+j, c2, c3);
										Gpfa5fStep3B(tu, ab, jh+j, jw+j, jm+j, jr+j, cs);
//
										Gpfa5fStep1(tu, ab, ji+j, jt+j, jn+j, js+j, c1);
										Gpfa5fStep2(tu, ab, jt+j, jx+j, jp+j, jd+j, jp+j, jd+j, c2, c3);
										Gpfa5fStep3B(tu, ab, ji+j, jx+j, jn+j, js+j, cs);
//
										Gpfa5fStep1(tu, ab, jv+j, jy+j, jo+j, jt+j, c1);
										Gpfa5fStep2(tu, ab, jv+j, jj+j, ju+j, je+j, ju+j, je+j, c2, c3);
										Gpfa5fStep3B(tu, ab, jj+j, jy+j, jo+j, jt+j, cs);
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

void FftEngineGpfaftNew::Gpfa5fStep1(Complex *tu, Complex *ab, int j1, int j2, int j3, int j4, real c1)
{
	tu[0].Add(ab[j1], ab[j2]);
	tu[1].Add(ab[j3], ab[j4]);
	tu[2].Sub(ab[j1], ab[j2]);
	tu[3].Sub(ab[j3], ab[j4]);
	tu[4].Add(tu[0], tu[1]);
	tu[5].Sub(tu[0], tu[1], c1);
}

void FftEngineGpfaftNew::Gpfa5fStep2(Complex *tu, Complex *ab, int j1, int j5, int j3, int j7, int j6, int j0, real c2, real c3)
{
	ab[j1].Copy(ab[j5]);
	tu[6].Sub2(ab[j6], tu[4], (real)0.25);
	tu[11].Add(ab[j6], tu[4]);
	tu[7].Add(tu[6], tu[5]);
	tu[8].Sub(tu[6], tu[5]);
	ab[j3].Copy(ab[j7]);
	tu [9].Sub0(tu[2], c3, tu[3], c2);
	tu[10].Add0(tu[2], c2, tu[3], c3);
	ab[j0].Copy(tu[11]);
}

void FftEngineGpfaftNew::Gpfa5fStep3A(Complex *tu, Complex *ab, int j1, int j2, int j3, int j4)
{
	ab[j1].set(tu[7].re - tu[10].im, tu[7].im + tu[10].re); 
	ab[j2].set(tu[7].re + tu[10].im, tu[7].im - tu[10].re); 
	ab[j3].set(tu[8].re - tu[ 9].im, tu[8].im + tu[ 9].re); 
	ab[j4].set(tu[8].re + tu[ 9].im, tu[8].im - tu[ 9].re); 
}

void FftEngineGpfaftNew::Gpfa5fStep3B(Complex *tu, Complex *ab, int j1, int j2, int j3, int j4, Complex *cs)
{
	real a = tu[7].re - tu[10].im;
	real b = tu[7].im + tu[10].re;
	ab[j1].set(cs[0].re * a - cs[0].im * b, cs[0].im * a + cs[0].re * b);
	a = tu[7].re + tu[10].im;
	b = tu[7].im - tu[10].re;
	ab[j2].set(cs[3].re * a - cs[3].im * b, cs[3].im * a + cs[3].re * b);
	a = tu[8].re - tu[9].im;
	b = tu[8].im + tu[9].re;
	ab[j3].set(cs[1].re * a - cs[1].im * b, cs[1].im * a + cs[1].re * b);
	a = tu[8].re + tu[9].im;
	b = tu[8].im - tu[9].re;
	ab[j4].set(cs[2].re * a - cs[2].im * b, cs[2].im * a + cs[2].re * b);
}

void FftEngineGpfaftNew::Gpfa5fStep2Simple(Complex *tu, Complex *ab, int j1, real c2, real c3)
{
	tu[6].Sub(ab[j1], tu[4], (real)0.25);
	ab[j1].Add(ab[j1], tu[4]);
	tu[7].Add(tu[6], tu[5]);
	tu[8].Sub(tu[6], tu[5]);
	tu [9].Sub0(tu[2], c3, tu[3], c2);
	tu[10].Add0(tu[2], c2, tu[3], c3);
}

void FftEngineGpfaftNew::Gpfa3fStep1(Complex *tu, Complex *ab, int j1, int j2, int j3, real c1)
{
	tu[0].Add(ab[j1], ab[j2]);
	tu[1].Sub2(ab[j3], tu[0], (real)0.5);
	tu[2].Sub(ab[j1], ab[j2], c1);
}

void FftEngineGpfaftNew::Gpfa3fStep3A(Complex *tu, Complex *ab, int j1, int j2)
{
	ab[j1].set(tu[1].re - tu[2].im, tu[1].im + tu[2].re); 
	ab[j2].set(tu[1].re + tu[2].im, tu[1].im - tu[2].re); 
}

void FftEngineGpfaftNew::Gpfa3fStep3B(Complex *tu, Complex *ab, int j1, int j2, Complex *cs)
{
	real a = tu[1].re - tu[2].im;
	real b = tu[1].im + tu[2].re;
	ab[j1].set(cs[0].re * a - cs[0].im * b, cs[0].im * a + cs[0].re * b);
	a = tu[1].re + tu[2].im;
	b = tu[1].im - tu[2].re;
	ab[j2].set(cs[1].re * b - cs[1].im * a, cs[1].im * b + cs[1].re * a);
}

void FftEngineGpfaftNew::Gpfa3f(Complex *ab, Complex *trigs, int inc, int jump, int n, int mm, int lot, FftDirection isign)
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

	const real onex_ = (real)1.;
	const real sin60 = (real)0.866025403784437;
	const int lvr = 64;

//	FILE *file15 = fopen("ccc.txt", "w");

	int n3 = powint(3, mm);
	int inq = n / n3;
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

	Complex cs[2], tu[4];
    int ipass, j, ja, jb, jc, jd, je, jf, jg, jh, ji, jjj, jstep, jstepl, laincl;
	int k, l, la, ll, nb, nu, nvex;
	int kk = 0;
//
// Loop on blocks of lvr transforms
	for(nb=0; nb<nblox; ++nb)
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
						Gpfa3fStep1(tu, ab, jb+j, jc+j, ja+j, c1);
						ab[ja+j].Add(ab[ja+j], tu[0]);
						Gpfa3fStep3A(tu, ab, jb+j, jc+j);
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
							Gpfa3fStep1(tu, ab, jb+j, jc+j, ja+j, c1);
							ab[ja+j].Add(ab[ja+j], tu[0]);
							Gpfa3fStep3B(tu, ab, jb+j, jc+j, cs);
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
//			fprintf(file15, "Ipass1 = %d\n", ipass);
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
//
							Gpfa3fStep1(tu, ab, jb+j, jc+j, ja+j, c1);
                            ab[jb+j].Copy(ab[jd+j]);
							ab[ja+j].Add(ab[ja+j], tu[0]);
							Gpfa3fStep3A(tu, ab, jd+j, jc+j);
//
							Gpfa3fStep1(tu, ab, je+j, jf+j, jb+j, c1);
							ab[jf+j].Copy(ab[jh+j]);
							ab[jb+j].Add(ab[jb+j], tu[0]);
							Gpfa3fStep3A(tu, ab, je+j, jh+j);
//
							Gpfa3fStep1(tu, ab, jf+j, ji+j, jg+j, c1);
							tu[0].Add(tu[0], ab[jg+j]);
							ab[jg+j].Copy(ab[jc+j]);
							ab[jc+j].Copy(tu[0]);
							Gpfa3fStep3A(tu, ab, jf+j, ji+j);
//
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
//				fprintf(file15, "K %d %d %d\n", k, ink, jstep-ink);
				LoadTrigs(cs, 2, trigs, kk, s);
//
// Double loop along first transform in block
				for(ll=k; ll<=(la-1)*ink; ll+=3*jstep)
				{
					for(jjj=ll; jjj<=(n-1)*inc; jjj+=3*la*ink)
					{
//						fprintf(file15, "Jjj3 = %d %d %d %d\n", jjj, ll, (n-1)*inc, 3*la*ink);
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
//
								Gpfa3fStep1(tu, ab, jb+j, jc+j, ja+j, c1);
								ab[jb+j].Copy(ab[jd+j]);
								ab[ja+j].Add(ab[ja+j], tu[0]);
								Gpfa3fStep3B(tu, ab, jd+j, jc+j, cs);
//
								Gpfa3fStep1(tu, ab, je+j, jf+j, jb+j, c1);
								ab[jf+j].Copy(ab[jh+j]);
								ab[jb+j].Add(ab[jb+j], tu[0]);
								Gpfa3fStep3B(tu, ab, je+j, jh+j, cs);
//
								Gpfa3fStep1(tu, ab, jf+j, ji+j, jg+j, c1);
								tu[0].Add(tu[0], ab[jg+j]);
								ab[jg+j].Copy(ab[jc+j]);
								ab[jc+j].Copy(tu[0]);
								Gpfa3fStep3B(tu, ab, jf+j, ji+j, cs);
//
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

void FftEngineGpfaftNew::Gpfa2f(Complex *ab, Complex *trigs, int inc, int jump, int n, int mm, int lot, FftDirection isign)
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

	int n2 = powint(2, mm);
	int inq = n / n2;
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
	int ipass, j, ja, jb, jc, jd, je, jf, jg, jh, ji, jj, jk, jl, jm, jn, jo, jp, jjj, jstep, jstepl;
	int k, l, la, mu, nb, nu, nvex, ll;
	real c1, c2, c3;
	int kk = 0;
//
// Loop on blocks of lvr transforms
	for(nb=0; nb<nblox; ++nb)
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
        if (mu == 3 || mu == 7) 
			c1 = -onex_;
        c2 = Sqrt(half_);
        if (mu == 3 || mu == 5) 
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
