#include "StdAfx.h"

#include "GreenFunctionManager.h"
#include "AbstractFftEngine.h"
#include "DDscatCommons.h"
#include "Functions.h"
#include "Timeit.h"

void Ddebug(const char *Header, int nx, int ny, int nz, Array4Stacked<Complex> *C)
{
	fprintf(stderr, Header);
	fprintf(stderr, "\n");
	fprintf(stderr, "%d %d %d\n", nx, ny, nz);

	for(int jx=0; jx<nx; ++jx)
	{
		for(int jy=0; jy<ny; ++jy)
		{
			for(int jz=0; jz<nz; ++jz)
			{
				if ((C->Value(jx,jy,jz,0).re != (real)0.) || (C->Value(jx,jy,jz,0).im != (real)0.) || 
					(C->Value(jx,jy,jz,1).re != (real)0.) || (C->Value(jx,jy,jz,1).im != (real)0.) || 
					(C->Value(jx,jy,jz,2).re != (real)0.) || (C->Value(jx,jy,jz,2).im != (real)0.) ||
				    (C->Value(jx,jy,jz,3).re != (real)0.) || (C->Value(jx,jy,jz,3).im != (real)0.) || 
					(C->Value(jx,jy,jz,4).re != (real)0.) || (C->Value(jx,jy,jz,4).im != (real)0.) || 
					(C->Value(jx,jy,jz,5).re != (real)0.) || (C->Value(jx,jy,jz,5).im != (real)0.))
				fprintf(stderr, "%4d %4d %4d %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf\n", jx, jy, jz, 
					C->Value(jx,jy,jz,0).re, C->Value(jx,jy,jz,0).im, C->Value(jx,jy,jz,1).re, C->Value(jx,jy,jz,1).im,  
					C->Value(jx,jy,jz,2).re, C->Value(jx,jy,jz,2).im, C->Value(jx,jy,jz,3).re, C->Value(jx,jy,jz,3).im,  
					C->Value(jx,jy,jz,4).re, C->Value(jx,jy,jz,4).im, C->Value(jx,jy,jz,5).re, C->Value(jx,jy,jz,5).im);
			}
		}
	}
}

GreenFunctionManager *GreenFunctionManager::item = NULL;
GreenFunctionManager *GreenFunctionManager::GetInstance(void)
{
	if (!item)
	{
		item = new GreenFunctionManager;
		item->Init();
	}
	return item;
}

GreenFunctionManager::GreenFunctionManager(void) 
{ 
	subSys[SubsystemTypeElectric] = new SubsystemElectric;
	subSys[SubsystemTypeMagnetic] = new SubsystemMagnetic;
}

GreenFunctionManager::~GreenFunctionManager(void) 
{ 
	CleanDelete(subSys[SubsystemTypeElectric]);
	CleanDelete(subSys[SubsystemTypeMagnetic]);
}

void GreenFunctionManager::Kill(void)
{
	CleanDelete(item);
}

GreenFunctionManager::Subsystem::Subsystem(void) : mySize(0)
{
	myType = SubsystemTypeEnd;
	cxzcg = NULL;
	dcxsum = NULL;
	tmp6 = NULL;
	issym = NULL;
}

GreenFunctionManager::Subsystem::Subsystem(int mSize, SubsystemType mt) : mySize(mSize)
{
	myType = mt;
	cxzcg = NULL;
	dcxsum = new Complex[mySize];
	tmp6 = new Complex[mySize];
	issym = new int[3 * mySize];
}

GreenFunctionManager::Subsystem::~Subsystem(void)
{
	CleanDelete(cxzcg);
	CleanDelete2(dcxsum);
	CleanDelete2(tmp6);
	CleanDelete2(issym);
}

void GreenFunctionManager::Init(void)
{
	subSys[SubsystemTypeElectric]->SetDimension(0, 0, 0);
	subSys[SubsystemTypeMagnetic]->SetDimension(0, 0, 0);
}

void GreenFunctionManager::SetDimension(int nx, int ny, int nz)
{
	subSys[SubsystemTypeElectric]->SetDimension(nx, ny, nz);
	subSys[SubsystemTypeMagnetic]->SetDimension(nx, ny, nz);
}

void GreenFunctionManager::Subsystem::SetDimension(int nx, int ny, int nz)
{
	this->nx = nx;
	this->ny = ny;
	this->nz = nz;
}

void GreenFunctionManager::SetIpbc(bool value)
{
	subSys[SubsystemTypeElectric]->Ipbc() = value;
	subSys[SubsystemTypeMagnetic]->Ipbc() = value;
}

bool &GreenFunctionManager::Ipbc(void)
{
	return subSys[SubsystemTypeElectric]->Ipbc();			// TODO
}

GreenFunctionManager::Subsystem *GreenFunctionManager::GetSusystem(SubsystemType type)
{
	if (type != SubsystemTypeEnd)
		return subSys[type];
	else
		return NULL;
}

GreenFunctionManager::Subsystem *GreenFunctionManager::GetElectric(void)
{
	return subSys[SubsystemTypeElectric];
}

GreenFunctionManager::Subsystem *GreenFunctionManager::GetMagnetic(void)
{
	return subSys[SubsystemTypeMagnetic];
}

bool GreenFunctionManager::Subsystem::Allocate(void)
{
	if (cxzcg)
		return false;
	cxzcg = new Array4Stacked<Complex>;
	if (!cxzcg)
		return false;
	Complex *res = cxzcg->Dimension(2*nx, 2*ny, 2*nz, mySize);
	cxzcg->Clear();
	return (res != NULL);
}

void GreenFunctionManager::Subsystem::Deallocate(void)
{
	cxzcg->Deallocate();
}

void GreenFunctionManager::SubsystemElectric::Cisi(real x, real &ci, real &si)
{
// Given:
//     X = real argument
// Returns
//     CI = gamma + ln(x) + \int_0^x [(cos(t)-1)/t] dt = "cosine integral"
//     SI = \int_0^x [sin(t)/t] dt = "sine integral"
// 
// Code adapted from Numerical Recipes in FORTRAN (2e) 
// by Press, Teukolsky, Vetterling, and Flannery (1994)
//
	const real eps = (real)6.e-8;
	const real euler = (real)0.57721566;
	const real piby2 = (real)1.5707963;
	const real fpmin = (real)1.e-30;
	const real tmin = (real)2.;
	const real onex_ = (real)1.;
	const real zero_ = (real)0.;
	const int maxit = 100;
//
	real t = Fabs(x);
	if(t == zero_)
	{
		si = zero_;
		ci =-onex_ / fpmin;
		return;
	}
//
	bool bOk = true;
	Complex b, c, d, h, del;
	if (t > tmin) 
	{
		b.set(onex_, t);
		c.set(onex_ / fpmin, zero_);
		d = Complex(onex_, zero_) / b;
		h = d;
		bOk = true;
		for(int i=0; i<maxit; ++i)			// i=1...
		{
			real a = -(real)(i*i);
			b += 2.;
			d = Complex(onex_, zero_) / (d * a + b);
			c = b + Complex(a, zero_) / c;
			del = c * d;
			h *= del;
			if ((del - Complex(onex_, zero_)).Absc() < eps)
			{
				bOk = false;
				break;
			}
		}
		if (bOk)
			Errmsg("Fatal", "Cisi", "Fatal error: failed in cisi");
		h *= Complex(Cos(t), -Sin(t));
		ci = -h.re;
		si = piby2 + h.im;
	}
	else
	{
		real sumc, sums, su, sgn, fact, term, err;
		if(t < Sqrt(fpmin))
		{
			sumc = zero_;
			sums = t;
		}
		else
		{
			su = sums = sumc = zero_;
			sgn = fact = onex_;
			bool odd = true;
			bOk = true;
			for(int k=0; k<maxit; ++k)	// k=1...
			{
				fact *= (t / (k+1));
				term = fact / (k+1);
				su += sgn * term;
				err = term / Fabs(su);
				if(odd)
				{
					sgn = -sgn;
					sums = su;
					su = sumc;
				}
				else
				{
					sumc = su;
					su = sums;
				}
				if(err < eps)
				{
					bOk = false;
					break;
				}
				odd = !odd;
			}
			if (bOk)
				Errmsg("Fatal", "Cisi", "maxits exceeded in cisi");
		}
		si = sums;
		ci = sumc + log(t) + euler;
	}
	if (x < zero_) 
		si = -si;
}

void GreenFunctionManager::SubsystemElectric::DeepDirectCalc(real r2, real phasyz, real akd, real akd2, real gammakd4, real *x)
{
	real r = Sqrt(r2);
	real r3 = r * r2;
//
// PHASYZ = phase at (1,JPY*NPY+1,JPZ*NPZ+1) - phase at (1,1,1)
	real kr = akd * r;
	Complex cxikr((real)0., kr);
//
// Include artificial factor exp[-(gamma*kr)^4] to assist convergence
	if (Common0::GetInstance()->Get((unsigned int)myType)->Idipint() == 0)
	{
		Complex cxphas = (cxikr + Complex((real)0., phasyz) - gammakd4 * r2 * r2).exp() / r3;
		Complex cxfac = Complex((real)1. / r2, -kr / r2);
//
// II=1      -> M=1   a_1 (xx)
//      JJ=2      2   a_2 (xy)
//         3      3   a_3 (xz)
//    2           4   a_4 (yy)
//         3      5   a_5 (yz)
//    3           6   a_6 (zz)
		int m = -1;
		for(int ii=0; ii<3; ++ii)									// loop over ii
		{
			++m;
			real xx = x[ii] * x[ii];
			dcxsum[m] -= cxphas * (cxfac * (r2 - (real)3. * xx) + akd2 * (xx - r2));
			if (ii < 2)
			{
				for(int jj=ii+1; jj<3; ++jj)							// loop over jj
				{
					++m;
					xx = x[ii] * x[jj];
					dcxsum[m] += cxphas * (cxfac * (real)3. - akd2) * xx;
				}												// end loop over jj
			}
		}														// end loop over ii
	}
	else
	{
// Expressions for the Green coefficients correspond to those in 
// "A library for computing the filtered and non-filtered 3D Green's tensor 
//  associated with infinite homogeneous space and surfaces" 
// by P. Gay-Balmaz and O. Martin (2002; Computer Physics Communications,
// 144, 111-120) apart from a factor of 4*pi.

		real ciplus, ciminus, siplus, siminus;
		Complex cxexpikr = cxikr.exp();
		Complex cxphas = Complex(-gammakd4 * r2 * r2, phasyz).exp();
		real kfr = Pi * r;
		Cisi(kfr + kr, ciplus, siplus);
		Cisi(kfr - kr, ciminus, siminus);
		real coskr = Cos(kr);
		real sinkr = Sin(kr);
		real coskfr = Cos(kfr);
		real sinkfr = Sin(kfr);
// CXA0 = alpha(R) of Gay-Balmaz & Martin 2002
// CXA1 = (d/dR) alpha   = alpha'
// CXA2 = (d2/dR2) alpha = alpha"

		real cxa0 = sinkr * (ciplus - ciminus) + coskr * (Pi - siplus - siminus);
		real cxa1 = akd * sinkr * (-Pi + siplus + siminus) + akd * coskr * (ciplus - ciminus) - (real)2. * sinkfr / r;
		real cxa2 = akd2 * (sinkr * (ciminus - ciplus) + coskr * (siplus + siminus - Pi)) + (real)2. * (sinkfr - kfr * coskfr) / r2;

// Yurkin, Min & Hoekstra 2010 define G_ij such that G_ij*P_j = E field at i
// so that elements of G_ij are same as elements of our matrix CXZC
// notation:
// CXG0   = g_F(R)              of YMH2010
// CXG1   = (d/dR) g_F   = g_F' of YMH2010
// CXG2   = (d2/dR2) g_F = g_F" of YMH2010  
// CX2PIH = 2*pi*h(R)           of YMH2010
//
// CXG0 = exp(ikR)/R - alpha/(pi*R)
// CXG1 = exp(ikR)*(-1+ikR)/R^2 - alpha'/(pi*R) + alpha/(pi*R^2)
// CXG2 = exp(ikR)*(2-2ikR-(kR)^2)/R^3 - (2*alpha-2R*alpha'+R^2*alpha")/(pi*R^3)

		Complex cxg0 = (cxexpikr - cxa0 / Pi) / r;
		Complex cxg1 = (cxexpikr * (cxikr - (real)1.) + (cxa0 - cxa1 * r) / Pi) / r2;
		Complex cxg2 = (cxexpikr * (Complex((real)2., (real)0.) - cxikr * (real)2. - kr * kr) - ((real)2. * (cxa0 - r * cxa1) + r2 * cxa2) / Pi) / r3;
		real cx2pih = (sinkfr - kfr * coskfr) / (Pi * r3);
//
		int m = -1;
		for(int ii=0; ii<3; ++ii)									// loop over ii
		{
			++m;
//! diagonal elements 
//! II=1: M=1 xx
//! II=2: M=4 yy
//! II=3: M=6 zz
			real xx = x[ii] * x[ii];
			dcxsum[m] += cxphas * (cxg0 * akd2 + cxg1 / r + ((real)(2./3.)) * cx2pih + (cxg2 / r2 - cxg1 / r3) * xx);
			if (ii < 2)
			{
				for(int jj=ii+1; jj<3; ++jj)							// loop over jj
				{
					++m;
//! off-diagonal elements
//! II=1: M=2 xy
//! II=1: M=3 xz
//! II=2: M=5 yz
					xx = x[ii] * x[jj];
					dcxsum[m] += cxphas * (cxg2 / r2 - cxg1 / r3) * xx;
				}													// enddo JJ
			}
		}															// enddo II
	}
}

void GreenFunctionManager::SubsystemMagnetic::DeepDirectCalc(real r2, real phasyz, real akd, real akd2, real gammakd4, real *x)
{
//     k^2                    1
// B = --- * exp(ikr) * (1 - --- ) * r x p
//     r^2                   ikr 

// p = exp(i*phasyz) * p(1,1,1)

//     k^2                                    1
//   = --- * exp(ikr) * exp(i*phasyz) * (1 - --- ) * r x p(1,1,1)
//     r^2                                   ikr

// f(r) = (k/r)^2 * exp(ikr+i*phasyz) * [1 - 1/(ikr)]

// B = f(r) * r x p(1,1,1)     

// B_x =  c_1*p_y + c_2*p_z
// B_y = -c_1*p_x + c_3*p_z
// B_z = -c_2*p_x - c_3*p_y

// c_1 = -f * r_z
// c_2 =  f * r_y
// c_3 = -f * r_x
//
// Include artificial factor exp[-(gamma*kr)^4] to assist convergence
	real r = Sqrt(r2);
	Complex cxikr((real)0., akd * r);
	Complex cxphas = (cxikr + Complex((real)0., phasyz) - gammakd4 * r2 * r2).exp();
	Complex cxfac = Complex((real)1., (real)0.) - Complex((real)1., (real)0.) / cxikr;  
	Complex cxcoeff = cxfac * cxphas * akd2 / r2;
	dcxsum[0] -= cxcoeff * x[2];
	dcxsum[1] += cxcoeff * x[1];
	dcxsum[2] -= cxcoeff * x[0];
}

void GreenFunctionManager::Subsystem::DirectCalc(int ix, int iy, int iz, Vect3<real> &dx, Vect3<real> &ak, real akd, real akd2, real gamma, real pyddx, real pzddx)
{
/* **
calculates magnetic field at (I,J,K) due to dipole at (1,1,1) plus its replicas
given:
       IX,IY,IZ = 1, 1, 1 to do octant I>0, J>0, K>0 only
                  1, 1,-1              I>0, J>0, all K
                  1,-1, 1              I>0, K>0, all J
                  1,-1,-1              I>0, all J, all K
                 -1, 1, 1              all I, J>0, K>0
                 -1, 1,-1              all I, J>0, all K
                 -1,-1, 1              all I, all J, K>0
                 -1,-1,-1              all I, all J, all K
       IPBC     = 0 if only doing first octant (IX=IY=IZ=1)
                = 1 otherwise
                  N.B.: IPBC affects dimensioning of CXZG
       GAMMA    = factor used to assist convergence of sums over
                  replica dipoles when IPBC > 0
                  contribution from replica dipoles is suppressed
                  by factor exp(-gamma*(kr)^4)
                  typical value gamma = 0.005
                  The effective
                  range/d = 1/(gamma*k*d) = 400 if gamma=5e-3 and kd=0.5
                  range/lambda = 1/(2*pi*gamma) = 31.8
                  The sums are actually continued out to
                  r/d = 2* /(gamma*kd) = 800 if gamma=5e-3 , kd=0.5
                  [screening factor = exp(-16)=1.1e-7]
       NX,NY,NZ = size of first octant (I = 1 -> NX,
                                        J = 1 -> NY,
                                        K = 1 -> NZ)
       IPBC     = 0 for isolated target
                = 1 to use periodic boundary conditions
                    (y direction, z direction, or both y and z directions)
                  N.B.: IPBC affects dimensioning of CXZC
       GAMMA     = real coefficient used to assist convergence of sums
                   over replica dipoles by suppressing long-range
                   contributions with factor exp(-gamma*(kr)^2)
                   typical value gamma = 0.005
                   The effective range/d = 1/(gamma*k*d) = 400 if gamma=5e-3 and kd=0.5
                   range/lambda = 1/(2*pi*gamma) = 31.8
                   The sums are actually continued out to
                   r/d = 2* /(gamma*kd) = 800 if gamma=5e-3 , kd=0.5
                   [screening factor = exp(-16)=1.1e-7]
       PYDDX     = period of lattice in y direction/d
       PZDDX     = period of lattice in z direction/d
       DX(1-3)   = lattice spacing in x,y,z direction, in units of d = n**(-1/3).
       AK(1-3)   = k(1-3)*d, where k = k vector in vacuo
       AKD       = |k|d
       AKD2      = |kd|^2
       CXZC      = array with dimensions
                   (NX+1)*(NY+1)*(NZ+1)*3 if IPBC=0
                   (2*NX)*(2*NY)*(2*NZ)*3 if IPBC=1

returns:
       CXZG(I,J,K,M) = c_M to calculate magnetic field B at (I,J,K)
                       produced by dipole at (1,1,1) and replica dipoles if IPBC = 0
                       B_x =           c_1*P_y + c_2*P_z
                       B_y = -c_1*P_x          + c_3*P_z
                       B_z = -c_2*P_x -c_3*P_y

JPYM = maximum number of periods in Y direction rmax = JPYM*PYDDX
note that even for short-range interaction, need to extend sums
from JPY=-1 to JPY=+1 to include "edge" interactions with next TUC

PYD=0. for nonperiodic case
   >0. for periodic boundary conditions in y direction

Adapted from original subroutine ESELF written by Jeremy Goodman

History records of Fortran versions removed.
Copyright (C) 2013 B.T. Draine and P.J. Flatau
Copyright (C) 2013, C++ versions, Choliy V.

This code is covered by the GNU General Public License.
End history
** */
	const char *MyLabel = "DirectCalc";
	const real precision = (real)1.e-6;
#ifdef openmp
	int nthreads, tid;
#endif
	real gammakd4 = gamma * akd;
	real range = (real)2. / gammakd4;
	gammakd4 = gammakd4 * gammakd4;
	gammakd4 = gammakd4 * gammakd4;
	real range2 = range * range;
	int jpym = (pyddx <= (real)0.) ? 0 : 1 + nint_(range / pyddx);
//
// Compute 3 independent elements of 3x3 anti-symmetric matrix C_jk, where
// C_jk*P_k = magnetic field at location j due to dipole P_k at location
//
// C_jk = ( 0   c_1  c_2)
//        (-c_1  0   c_3)
//        (-c_2 -c_3  0 )_jk
//
// Initialize CXG(I,J,K,M) = c_M for magnetic field at (I,J,K)
//                          produced by a dipole at (1,1,1)
//                          and replica dipoles (if PYD or PYZ are nonzero).
//
// Need to calculate this for all (I,J,K) for one octant:
// I running from 1 to NX, J from 1 to NY, K from 1 to NZ
//
// X0,Y0,Z0 = X(I,J,K) - X(1,1,1) = vector from dipole location (1,1,1)
//                                  to point where B is to be calculated
// IX = +1 -> IMIN = 1
// IX = -1 -> IMIN = 2-NX   (1-IMIN = NX-1)

	int imin = (1 - nx) * (1 - ix) / 2;
	int jmin = (1 - ny) * (1 - iy) / 2;
	int kmin = (1 - nz) * (1 - iz) / 2;
//
// Determine elapsed cpu time for these sums
	static real t1 = (real)0., t2 = (real)0., t0 = Cpu_time();
//
#ifdef openmp
! fork a team of threads 
!$OMP PARALLEL               &
!$OMP&   PRIVATE(NTHREADS,TID)

      TID=OMP_GET_THREAD_NUM()
!      WRITE(0,*)'direct_calc ckpt 0.9: hello world from thread = ',TID

! only master thread does this:

      IF(TID.EQ.0)THEN
         NTHREADS=OMP_GET_NUM_THREADS()
         WRITE(CMSGNM,FMT='(A,I3)')'number of OpenMP threads =',NTHREADS
         CALL WRIMSG('DIRECT_CALCB',CMSGNM)
      ENDIF
!$OMP END PARALLEL
#endif

#ifdef openmp
! 2012.04.27 (BTD) added R to private variables
!                  removed FLUSH directive (should not be needed)
!$OMP PARALLEL                                   &
!$OMP&   PRIVATE(I,II,J,JPY,JPZ,K,M)             &
!$OMP&   PRIVATE(JCXZG,JPZM,KCXZG)               &
!$OMP&   PRIVATE(PHASY,PHASYZ,R,R2,RJPY,RJPZ) &
!$OMP&   PRIVATE(X,X0,X2,X2Y2,XX,Y0,Z0)       &
!$OMP&   PRIVATE(CXFAC,CXTERM,CXPHAS,DCXSUM)
!$OMP DO
#endif
//
    real x[3];
    int i, j, k;

	if (Common10::GetInstance()->Myid() == 0)		// For first dipole, determine time to sum over replicas
		t1 = Cpu_time();

	const real *dxData = dx.data;
	const real *akData = ak.data;

	for(k=kmin; k<nz; ++k)												// Loop over k
	{
		real z0 = dxData[2] * k;
		int kc = (k >= 0) ? k : 2*nz + k;
//
		for(j=jmin; j<ny; ++j)											// Loop over j
		{
			real y0 = dxData[1] * j;
			int jc = (j >= 0) ? j : 2*ny + j;
//
			for(i=imin; i<nx; ++i)										// Loop over i
			{
				real x0 = dxData[0] * i;
				int ic = (i >= 0) ? i : 2*nx + i;

				x[0] = x0;
				real x2 = x0 * x0;
				memset(dcxsum, 0, mySize * sizeof(Complex));
//
// JPY=0, JPZ=0 gives B field from dipole at (1,1,1)
// general JPY,JPZ gives B field from dipole at (1,JPY*NPY+1,JPZ*NPZ+1)
// replica dipoles have same magnitude of polarization as dipole (1,1,1), but different phase.
// PHASYZ = phase of replica dipole - phase of dipole (1,1,1)
				for(int jpy=-jpym; jpy<=jpym; ++jpy)										// loop over JPY
				{
					real rjpy = jpy * pyddx * dxData[1];
					x[1] = y0 - rjpy;
					real x2y2 = x2 + x[1] * x[1];
					real phasy = akData[1] * rjpy;
//
// PZDDX=0. for nonperiodic case
//      >0. for periodic boundary conditions in z direction
					int jpzm = (pzddx <= (real)0.) ? 0 : 1 + nint_(Sqrt(max(range2 - rjpy * rjpy, (real)0.)) / pzddx);
//
					for(int jpz=-jpzm; jpz<=jpzm; ++jpz)									// loop over JPZ
					{
						real rjpz = jpz * pzddx * dxData[2];
						x[2] = z0 - rjpz;
						real r2 = x2y2 + x[2] * x[2];
//
// Skip the self-interaction (R=0) case (I=J=K=1 and JPY=JPZ=0)
						if (r2 > precision)
						{
//
// PHASYZ = phase at (1,JPY*NPY+1,JPZ*NPZ+1) - phase at (1,1,1)
							real phasyz = phasy + akData[2] * rjpz;
							DeepDirectCalc(r2, phasyz, akd, akd2, gammakd4, x);
						}
					}																// end loop over jpz
				}																	// end loop over jpy
				cxzcg->InsertLastDim(ic, jc, kc, dcxsum);
//
				if((Common10::GetInstance()->Myid() == 0) && (i == imin) && (j == jmin) && (k == kmin))
				{

// NB: this call to CPU_TIME occurs within the first thread (I==IMIN, J==JMIN, K==KMIN)
					t2 = Cpu_time();
					t2 = (real)((nx-imin+1) * (ny-jmin+1) * (nz-kmin+1)) * (t2 - t1);
					PrintTimes(MyLabel, "Estimated", t2);
					PrintTimeForecasts(MyLabel, pyddx, pzddx);
				}
			}																	// End loop over i
		}																		// End loop over j
	}																			// End loop over k
#ifdef openmp
!$OMP END DO
!$OMP END PARALLEL
#endif

	t2 = Cpu_time();
	t2 = t2 - t0;
	if (Common10::GetInstance()->Myid() == 0)
		PrintTimes(MyLabel, "Actual", t2);
//
// If IPBC=1: set the elements with ICXZG=NX+1 or JCXZG=NY+1 or KCXZG=NZ+1 to zero
	if (ipbc == true)
		cxzcg->ClearIntoPlanes(nx, ny, nz);
}

void GreenFunctionManager::SubsystemElectric::LoadIssym(void)
{
	int iss[] = { 1, 1, 1,  -1, -1, 1,  -1, 1, -1,   1, 1, 1,   1, -1, -1,   1, 1, 1 };
	memcpy(issym, iss, 3 * mySize * sizeof(int));
}

void GreenFunctionManager::SubsystemMagnetic::LoadIssym(void)
{
	int iss[] = { 1, 1, -1,  1, -1, 1,  -1, 1, 1 };
	memcpy(issym, iss, 3 * mySize * sizeof(int));
}

void GreenFunctionManager::SubsystemElectric::PrepareTmp3(Complex *tmp3, int ign, int jgn, int kgn, unsigned int *pos2, Complex *cxzwData)
{
	if (ign * jgn < 0) tmp6[1] = -tmp6[1];
	if (ign * kgn < 0) tmp6[2] = -tmp6[2];
	if (jgn * kgn < 0) tmp6[4] = -tmp6[4];
	tmp3[0] = tmp6[0] * cxzwData[pos2[0]] + tmp6[1] * cxzwData[pos2[1]] + tmp6[2] * cxzwData[pos2[2]];
	tmp3[1] = tmp6[1] * cxzwData[pos2[0]] + tmp6[3] * cxzwData[pos2[1]] + tmp6[4] * cxzwData[pos2[2]];
	tmp3[2] = tmp6[2] * cxzwData[pos2[0]] + tmp6[4] * cxzwData[pos2[1]] + tmp6[5] * cxzwData[pos2[2]];
}

void GreenFunctionManager::SubsystemElectric::PrepareTmp3(Complex *tmp3, unsigned int *pos2, Complex *cxzwData)
{
	tmp3[0] = tmp6[0] * cxzwData[pos2[0]] + tmp6[1] * cxzwData[pos2[1]] + tmp6[2] * cxzwData[pos2[2]];
	tmp3[1] = tmp6[1] * cxzwData[pos2[0]] + tmp6[3] * cxzwData[pos2[1]] + tmp6[4] * cxzwData[pos2[2]];
	tmp3[2] = tmp6[2] * cxzwData[pos2[0]] + tmp6[4] * cxzwData[pos2[1]] + tmp6[5] * cxzwData[pos2[2]];
}

void GreenFunctionManager::SubsystemMagnetic::PrepareTmp3(Complex *tmp3, int ign, int jgn, int kgn, unsigned int *pos2, Complex *cxzwData)
{
	if (ign * jgn < 0) tmp6[0] = -tmp6[0];
	if (ign * kgn < 0) tmp6[1] = -tmp6[1];
	if (jgn * kgn < 0) tmp6[2] = -tmp6[2];
	tmp3[0] =  tmp6[0] * cxzwData[pos2[1]] + tmp6[1] * cxzwData[pos2[2]];
	tmp3[1] = -tmp6[0] * cxzwData[pos2[0]] + tmp6[2] * cxzwData[pos2[2]];
	tmp3[2] = -tmp6[1] * cxzwData[pos2[0]] - tmp6[2] * cxzwData[pos2[1]];
}

void GreenFunctionManager::SubsystemMagnetic::PrepareTmp3(Complex *tmp3, unsigned int *pos2, Complex *cxzwData)
{
	tmp3[0] = tmp6[0] * cxzwData[pos2[1]] + tmp6[1] * cxzwData[pos2[2]];
	tmp3[1] =-tmp6[1] * cxzwData[pos2[0]] + tmp6[2] * cxzwData[pos2[2]];
	tmp3[2] =-tmp6[2] * cxzwData[pos2[0]] - tmp6[2] * cxzwData[pos2[1]];
}

void GreenFunctionManager::Subsystem::Self(FftMethod cmethd, Complex *cxzp, real gamma, real pyd, real pzd, Vect3<real> &ak, real akd, Vect3<real> &dx, Array4Stacked<Complex> *cxzw, Complex *cxzeb)
{
/* **
Parameter GAMMA determines the range of the sums when periodic boundary conditions are employed.  
Dipole-dipole interactions are screened by a factor exp[-(gamma*kr)^4]
The effective range/d = 1/(gamma*k*d) = 2000 if gamma=1e-3 and kd=0.5

The sums are actually continued out to r/d = 2/(gamma*kd) = 4000 if gamma=1e-3 , kd=0.5 [screening factor = exp(-16)=1.1e-7]

NB: module DDCOMMON_0 must have previously set values of AK2OLD,AK3OLD,WOLD to be used by ESELF

Given the dipole moments, CXZP, at all points on a rectangular grid, oscillating at frequency AKD, compute the electric field 
amplitude, CXZE, at each point produced by all the other dipoles except the one at that point.
The relationship between the dipoles and the field values at the grid points can be expressed as a convolution,
and the convolution is efficiently evaluated using 3D FFTs.

options for computation of 3-dimensional FFT:

if CMETHD='GPFAFT':
          Use CXFFT3N interface to GPFA code of Temperton.
          Good points:
             -On CRAY the code is on average 30-40% faster in comparison to TMPRTN (and 10-15 faster in comparison to BRENNR)
             -On scalar machines the code is 2-8 faster in comparison to BRENNR or TMPRTN
             -Doesn't require additional storage
          Limitations:
              -Requires that NX,NY,NZ be of form (2**I)*(3**J)*(5**K); subroutine EXTEND takes care of choosing suitable NX,NY,
              -The choice of "lvr" variable in gpfa2f, gpfa3f, gpfa5f depends on machine. WARNING: on C90 use lvr=128, on
               all other CRAY's use lvr=64; wrong lvr will produce WRONG RESULTS. On scalar machines optimal lvr depends on cache
               length; sub-optimal choice degrades performance but still produces correct results.

   if CMETHD='FFTW21'
      Use CXFFTW interface to FFTW (Fastest Fourier Transform in the West) version 2.1.x from Frigo and Johnson.

   if CMETHD='FFTMKL':
      Use CXFFT3_MKL interface to Intel Math Kernel Library (MKL) FFT

INPUT:
       CXZP(I,J,K,L)   Lth cartesian component of the dipole moment at the grid point (I,J,K);
                       the DIMENSIONed length of CXZP in the calling routine is CXZP(NX,NY,NZ,3) [or CXZP(NX*NY*NZ,3) or CXZP(3*NX*NY*NZ)]
       NX,NY,NZ        Size of grid in x,y,z directions (INTEGER).
       PYD             (Period of lattice in y direction)/D(2)
       PZD             (Period of lattice in z direction)/D(3)
       DX(1-3)         Lattice spacing in x,y,z directions, in units of n**(-1./3.) .  Note that with this normalization we have DX(1)*DX(2)*DX(3)=1.
       AK(1-3)         k(1-3)*d, where k = k vector in vacuo, and d = effective lattice spacing = (dx*dy*dz)**(1/3
       AKD             Frequency of oscillation of dipoles and electric field; also absolute value of wave vector of incident wave [c=1]  (REAL*4).
       CXZC            (NX+1)*(NY+1)*(NZ+1)*6 array of Green function coefficients used internally by ESELF and
                       *************NOT TO BE OVERWRITTEN*********** between calls, because these coefficients are recomputed only if W has changed 
					   since the last call to ESELF.
       CXZW            Complex, scratch-space vector of length: 2*NX*2*NY*2*NY*3
                       See comment about FFT usage and CMETHD flag. Can be overwritten between calls to ESELF

OUTPUT:
       CXZE(I,J,K,L)   Lth component of dipole-generated electric field at grid point (I,J,K);
                       the DECLARED length of ZE in the calling program is CXZE(NX,NY,NZ,3)
                       [or CXZE(NX*NY*NZ,3) or CXZE(3*NX*NY*NZ)]

Originally written by Jeremy Goodman, Princeton Univ. Observatory, 90.
History:
Fortran versions history removed.

Copyright (C) 1993,1994,1996,1997,1999,2000,2003,2004,2005,2006,2007,
              2008,2011,2012 B.T. Draine and P.J. Flatau
This code is covered by the GNU General Public License.
** */
//
	const char *MyLabel = "Eself ";
	const real half_ = (real)0.5;
	const real precision = (real)1.e-6;
#ifdef openmp
	int nthreads, tid;
#endif
	AbstractFftEngine *fftEngine = AbstractFftEngine::GetEngine(cmethd);
	int nz2 = 2 * nz;
	int ny2 = 2 * ny;
	int nx2 = 2 * nx;
    int m;
//
// Check if we can skip recomputation of Green-function coefficients
	bool bSkip = Common0::GetInstance()->Get((unsigned int)myType)->Check(precision, akd, ak.data[1], ak.data[2], pyd, pzd);
//
// We have to recompute the Green-function coefficients giving components of field strength at a each grid point R
// produced by unit-valued component of dipole moment at point R', and then Fourier transform these components.
	if (!bSkip)
	{
		Common0::GetInstance()->Get((unsigned int)myType)->Set(ak.data[1], ak.data[2], akd);
		Common0::GetInstance()->Get((unsigned int)myType)->Ngrid() = nx2 * ny2 * nz2;
		real akd2 = akd * akd;
//
// We assume screening function exp(-(gamma*kr)^4) so range/d = 1/(gamma*kd) = 2000 if gamma=1e-3 and kd=0.5 although the sums are actually continued out to
// r/d = 2/(gamma*kd) = 4000 if gamma=1e-3, kd=0.5 [screening factor = exp(-16)=1.1e-7]
//
// PYD*d(2) = periodicity in Y direction
// PZD*d(3) = periodicity in Z direction
		if ((pyd > (real)0.) || (pzd > (real)0.))
		{
			char cmsgnm[80];
            sprintf(cmsgnm, "(PBC with PYD, PZD=%8.2lf%8.2lf, GAMMA=%9.2e", pyd, pzd, gamma);
		    Wrimsg(MyLabel, cmsgnm);
		}
//
// Compute the actual coefficients:
// Compute 6 independent elements of 3x3 symmetric matrix A_jk,where
// A_jk*P_k = -electric field at location j due to dipole P_k at location
// A_jk = (a_1  a_2  a_3)
//        (a_2  a_4  a_5)
//        (a_3  a_5  a_6)_jk
		if (ipbc == false)
		{
// initialize CXZC(I,J,K,M) = a_M for electric field at (I,J,K) produced by a dipole at (1,1,1) and replica dipoles (if PYD or PYZ are nonzero).
// need to calculate this for all (I,J,K) for one octant: I running from 1 to NX, J from 1 to NY, K from 1 to NZ
//
// Later obtain A_jk values for other octants by using symmetry
			DirectCalc(1, 1, 1, dx, ak, akd, akd2, gamma, pyd * dx.data[1], pzd * dx.data[2]);
//			Ddebug("Cxzc inside Self", nx, ny, nz, cxzcg);
//
// At this point, CXZC(I,J,K,1-6) contains the upper triangular part of t symmetric 3 x 3 matrix giving the electric field at grid point (i,j,k)
// produced by a dipole at (1,1,1)
// Fill out CXZC to twice the size in each grid dimension to accomodate negative lags [periodicity in each dimension is assumed, so (e.g.)
// nx < i <= 2*nx is equivalent to -nx < i <= 0], exploiting symmetries, and then Fourier transform.
// If PYDDX=0 and PZDDX=0 , need only do direct calculation of A matrix for first octant, since remaining octants can be obtained by symmetry.
// After calculating A matrix, store only the first octant of the transform, since the others can be obtained by symmetry.
//

//
// Extend a_1 = a_xx(x,y,z) : a -> +a for x -> -x   +a y -> -y   +a z -> -z
// Extend a_2 = a_xy(x,y,z) : a -> -a for x -> -x   -a y -> -y   +a z -> -z
// Extend a_3 = a_xz(x,y,z) : a -> -a for x -> -x   +a y -> -y   -a z -> -z
// Extend a_4 = a_yy(x,y,z) : a -> +a for x -> -x   +a y -> -y   +a z -> -z
// Extend a_5 = a_yz(x,y,z) : a -> +a for x -> -x   -a y -> -y   -a z -> -z
// Extend a_6 = a_zz(x,y,z) : a -> +a for x -> -x   +a y -> -y   +a z -> -z
			LoadIssym();
			for(int m=0; m<mySize; ++m)
			{
				Extnd(cxzcg->SubarrayData(m), issym + 3*m, cxzw->SubarrayData(0));
//				Ddebug("Cxzc inside Self after Extnd", nx, ny, nz, cxzcg);

				fftEngine->DoFFT(cxzw->SubarrayData(0), nx2, ny2, nz2, FftForward);
//				Ddebug("Cxzc inside Self after DoFFT", nx, ny, nz, cxzcg);

				Trim(cxzw->SubarrayData(0), cxzcg->SubarrayData(m));
//				Ddebug("Cxzc inside Self after Trim", nx, ny, nz, cxzcg);
			}
		}
		else
		{
// This point is reached when PYDDX or PZDDX are nonzero.
// When PBC are used for general direction of incident wave, all octants of A matrix require direct calculation: symmetries valid
// for single target no longer apply because of position-dependent phases of replica dipoles.
			DirectCalc(-1, -1, -1, dx, ak, akd, akd2, gamma, pyd * dx.data[1], pzd * dx.data[2]);

// The array CXZC(1-2*NX,1-2*NY,1-2*NZ,1-6) of A matrix coefficients now covers all octants.
// Fourier transform the A matrix CXZC:
			for(int m=0; m<mySize; ++m)
			{
				fftEngine->DoFFT(cxzcg->SubarrayData(m), nx2, ny2, nz2, FftForward);
			}
// CXZC now contains the full Fourier transform of the A convolution and should not be overwritten between calls to ESELF
		}
// End of recomputation of Green-function coefficients
	}
//
// Fourier transform the polarizations:
    for(m=0; m<3; ++m)
	{
		PadC(cxzp, nx, ny, nz, m, cxzw->SubarrayData(m));
		fftEngine->DoFFT(cxzw->SubarrayData(m), nx2, ny2, nz2, FftForward);
	}
//
// Multiply by F.t. of Green function.
	if (ipbc == false)
	{
//
// If IPBC=0, then only one octant of F.t. of Green function has been stored, but can recover others using symmetry.
#ifdef openmp
!$OMP PARALLEL DO                                            &
!$OMP&   PRIVATE(K,J,I,KSGN,KR,JSGN,JR,ISGN,IR)              &
!$OMP&   PRIVATE(CXXX,CXXY,CXXZ,CXYY,CXYZ,CXZZ,CXEX,CXEY,CXEZ)
#endif

		const int lastDimStep = nx2 * ny2 * nz2;
		Complex *cxzcData = cxzcg->GetData();				// ChB Those two arrays are stored in fortran-stacked format
		Complex *cxzwData = cxzw->GetData();
		Complex tmp3[3];
		unsigned int pos1, pos2[3];

		for(int k=0; k<nz2; ++k)
		{
			int ksgn = nint_(sign_(nz + half_ - k));
			int kr = min_(k, nz2 - k);
			for(int j=0; j<ny2; ++j)
			{
				int jsgn = nint_(sign_(ny + half_ - j));
				int jr = min_(j, ny2 - j);
				for(int i=0; i<nx2; ++i)
				{
					int isgn = nint_(sign_(nx + half_ - i));
					int ir = min_(i, nx2 - i);

					pos1 = ir + (jr + kr * ny2) * nx2;
					pos2[0] = i + (j + k * ny2) * nx2;
					pos2[1] = pos2[0] + lastDimStep;
					pos2[2] = pos2[1] + lastDimStep;
					for(int ii=0; ii<mySize; ++ii, pos1+=lastDimStep)
					{
						tmp6[ii] = cxzcData[pos1];
					}
					PrepareTmp3(tmp3, isgn, jsgn, ksgn, pos2, cxzwData);
					cxzwData[pos2[0]] = tmp3[0];
					cxzwData[pos2[1]] = tmp3[1];
					cxzwData[pos2[2]] = tmp3[2];
				}
			}
		}
#ifdef openmp
!$omp end parallel do
#endif
	}
	else
	{
//
// If IPBC=1, then the full F.t. of the Green function has been stored.
#ifdef openmp
!$OMP PARALLEL DO                                                  &
!$OMP&   PRIVATE(K,J,I,CXXX,CXXY,CXXZ,CXYY,CXYZ,CXZZ,CXEX,CXEY,CXEZ)
#endif

		const int lastDimStep = nx2 * ny2 * nz2;
		Complex tmp3[3];
		Complex *cxzcData = cxzcg->GetData();
		Complex *cxzwData = cxzw->GetData();
		unsigned int pos1, pos2[3];
		for(int k=0; k<nz2; ++k)
		{
			for(int j=0; j<ny2; ++j)
			{
				for(int i=0; i<nx2; ++i)
				{
					pos1 = pos2[0] = i + nx2 * (j + ny2 * k);
					pos2[1] = pos2[0] + lastDimStep;
					pos2[2] = pos2[1] + lastDimStep;
					for(int ii=0; ii<mySize; ++ii, pos1+=lastDimStep)
					{
						tmp6[ii] = cxzcData[pos1];
					}
					PrepareTmp3(tmp3, pos2, cxzwData);
					cxzwData[pos2[0]] = tmp3[0];
					cxzwData[pos2[1]] = tmp3[1];
					cxzwData[pos2[2]] = tmp3[2];
				}
			}
		}
#ifdef openmp
!$OMP END PARALLEL DO
#endif
	}
//      
// Inverse Fourier transform to obtain electric field:
	for(m=0; m<3; ++m)
	{
		fftEngine->DoFFT(cxzw->SubarrayData(m), nx2, ny2, nz2, FftBackward);
//
// Note: the Convex FFT routine already normalizes result.
//	For other FFT routines need to divide result by NGRID
		Complex fac((real)1. / (real)Common0::GetInstance()->Get((unsigned int)myType)->Ngrid(), (real)0.);

#ifdef openmp
!$OMP PARALLEL DO     &
!$OMP&   PRIVATE(I,J,K)
#endif
		ExportCorner(cxzw->SubarrayData(m), cxzeb, nx, ny, nz, 3, m, fac);
#ifdef openmp
!$omp end parallel do
#endif
	}
}

//
// Both cxa and cxb have identical dimensions (2*nx, 2*ny, 2*nz) and they are stored in Fortran manner
void GreenFunctionManager::Subsystem::Extnd(Complex *cxa, int *isym, Complex *cxb)
{
// Using symmetries, extend first octant of coefficient matrix to other octants.
//
// Originally written by Jeremy Goodman
// 080620 (ASL) modified to use OpenMP
// 080627 (BTD) eself_v4: further mods related to OpenMP
//              reordered several nested do loops
//              introduced variables J2 and K2 to speed computations
// end history
// Copyright (C) 1993, B.T. Draine and P.J. Flatau
// This code is covered by the GNU General Public License.

#ifdef openmp
	USE OMP_LIB   !Art for OpenMP function declarations
#endif

#ifdef openmp
!$OMP PARALLEL              & 
!$OMP&   PRIVATE(I,J,J2,K,K2)
!$OMP DO
#endif
	int sz[3] = { 2*nx, 2*ny, 2*nz };

	int i, j, k, index[4];
	for(k=0; k<nz; ++k)
	{
		index[2] = k * sz[1];
		for(j=0; j<ny; ++j)
		{
			index[1] = (index[2] + j) * sz[0];
			memcpy(cxb+index[1], cxa+index[1], nx*sizeof(cxa[0]));
		}
	}
#ifdef openmp
!$OMP ENDDO
#endif

// x -> -x

// !btd 080627 moved I=NX+1 out of loop
// ! the SINGLE directive specifies that enclosed code is to be executed
// ! by only one thread in the team
#ifdef openmp
!$OMP SINGLE
#endif

	for(k=0; k<nz; ++k)
	{
		index[2] = k * sz[1];
		for(j=0; j<ny; ++j)
		{
			index[1] = (index[2] + j) * sz[0];
			cxb[index[1] + nx].clear(); 
		}
	}

#ifdef openmp
!$OMP END SINGLE
!$OMP DO 
#endif

	for(k=0; k<nz; ++k)
	{
		index[0] = k * sz[1];
		for(j=0; j<ny; ++j)
		{
			index[1] = (index[0] + j) * sz[0];
			for(i=nx+1; i<sz[0]; ++i)
			{
				index[2] = index[1] + i;
				index[3] = index[1] + (sz[0] - i);
				cxb[index[2]] = cxa[index[3]] * isym[0];
			}
		}
	}

#ifdef openmp
!$OMP END DO
#endif

// y -> -y

// !btd 080627 moved J=NY+1 out of loop, switched order of loops I and J
// the SINGLE directive specifies that enclosed code is to be executed
// by only one thread in the team
#ifdef openmp
!$omp single
#endif

	for(k=0; k<nz; ++k)
	{
		index[0] = (k * sz[1] + ny) * sz[0];
		for(i=0; i<sz[0]; ++i)
		{
			cxb[index[0] + i].clear();
		}
	}

#ifdef openmp
!$OMP END SINGLE
!$OMP DO 
#endif

	for(k=0; k<nz; ++k)
	{
		index[0] = k * sz[1];
		for(j=ny+1; j<sz[1]; ++j)
		{
			int j2 = sz[1] - j;
			index[1] = (index[0] + j) * sz[0];
			index[2] = (index[0] + j2) * sz[0];
			for(i=0; i<sz[0]; ++i)
			{
				cxb[index[1] + i] = cxb[index[2] + i] * isym[1];
			}
		}
	}

#ifdef openmp
!$omp end do
#endif

// ! z -> -z

// !Art pulling the 3rd dimension to the outer most loop.
// !Art we'll do this expression in only the thread that has NZ+1
// !btd 080627 reorder loops: J,I,K -> K,J,I
// ! the SINGLE directive specifies that enclosed code is to be executed
// ! by only one thread in the team
#ifdef openmp
!$OMP SINGLE
#endif

	for(j=0; j<sz[1]; ++j)
	{
		index[0] = (nz * sz[1] + j) * sz[0];
		for(i=0; i<sz[0]; ++i)
		{
			cxb[index[0] + i].clear();
		}
	}

#ifdef openmp
!$OMP END SINGLE
!$OMP DO
#endif

	for(k=nz+1; k<sz[2]; ++k)
	{
		int k2 = sz[2] - k;
		index[0] = k * sz[1];
		index[1] = k2 * sz[1];
		for(j=0; j<sz[1]; ++j)
		{
			index[2] = (index[0] + j) * sz[0];
			index[3] = (index[1] + j) * sz[0];
			for(i=0; i<sz[0]; ++i)
			{
				cxb[index[2] + i] = cxb[index[3] + i] * isym[2];
			}
		}
	}

#ifdef openmp
!$OMP END DO
!$OMP END PARALLEL
#endif
}

void GreenFunctionManager::Subsystem::Trim(Complex *cxb, Complex *cxa)
{
// Copy the first octant of CXB into CXA
//
// Originally written by Jeremy Goodman
// History:
// 92.04.20 (BTD) removed ISYM from argument list of TRIM:
// 06.09.28 (BTD) eself v2.0
// 08.06.27 (BTD) eself_v4
//                
// Copyright (C) 1993,2006 B.T. Draine and P.J. Flatau
// This code is covered by the GNU General Public License.

#ifdef openmp
!$OMP PARALLEL             &
!$OMP&    PRIVATE(I,J,K)   &
!$OMP&    SHARED(NX,NY,NZ) &
!$OMP&    SHARED(CXA,CXB)
!$OMP DO
#endif

	int sz[3] = { 2*nx, 2*ny, 2*nz };
	int index[3];
	for(int k=0; k<=nz; ++k)
	{
		index[2] = k * sz[1];
		for(int j=0; j<=ny; ++j)
		{
			index[1] = (index[2] + j) * sz[0];
			memcpy(cxa+index[1], cxb+index[1], (nx+1) * sizeof(cxb[0]));
		}
	}

#ifdef openmp
!$OMP END DO
!$OMP END PARALLEL
#endif
}

void GreenFunctionManager::Subsystem::PadC(Complex *cxa, unsigned int nx, unsigned int ny, unsigned int nz, unsigned int m, Complex *cxb)
{
// Pad the array CXA with zeros out to twice its length in each dimension, and put the result in CXB.
//
// Originally written by Jeremy Goodman
//
// Copyright (C) 1993, B.T. Draine and P.J. Flatau
// This code is covered by the GNU General Public License.

#ifdef openmp
!$OMP PARALLEL        &
!$OMP&   PRIVATE(I,J,K)
!$OMP DO
#endif

	int size = 8 * nx * ny * nz;
	memset(cxb, 0, size * sizeof(Complex));
//	for(int i=0; i<size; ++i)
//	{
//		cxb[i] = czero;
//	}
	//bool bRes = cxb.Set(2*nx, 2*ny, 2*nz, czero);
	//if (!bRes)
	//{
	//	Errmsg("Fatal", "Pad", "Cannot zero cxb");
	//}

#ifdef openmp
!$OMP END DO
!$OMP DO
#endif

	ImportCorner(cxb, cxa, nx, ny, nz, 3, m);

#ifdef openmp
!$OMP END DO
!$OMP END PARALLEL
#endif
}

// Values from C-like array f are imported in this->data like Stacked, from 0 to (ax, ay, az) and only for this->lastDim == m
void GreenFunctionManager::Subsystem::ImportCorner(Complex *to, Complex *from, unsigned int ax, unsigned int ay, unsigned int az, unsigned int at, unsigned int m)
{
	unsigned int dim[3] = { 2*ax, 2*ay, 2*az };
	for(unsigned int k=0; k<az; ++k)
	{
		for(unsigned int j=0; j<ay; ++j)
		{
			for(unsigned int i=0; i<ax; ++i)
			{
				to[(k*dim[1] + j)*dim[0] + i] = from[((i*ay + j)*az + k)*at + m];
			}
		}
	}
}

// From is Fortran-like(2*ax,2*ay,2*az) array, to is C-like (ax,ay,az)
void GreenFunctionManager::Subsystem::ExportCorner(Complex *from, Complex *to, unsigned int ax, unsigned int ay, unsigned int az, unsigned int at, unsigned int m, const Complex &scale)
{
	unsigned int dim[3] = { 2*ax, 2*ay, 2*az };
	for(unsigned int k=0; k<az; ++k)
	{
		for(unsigned int j=0; j<ay; ++j)
		{
			for(unsigned int i=0; i<ax; ++i)
			{
				to[((i*ay + j)*az + k)*at + m] = from[(k*dim[1] + j)*dim[0] + i] * scale;
			}
		}
	}
}

void GreenFunctionManager::Subsystem::PrintTimes(const char *label, const char *prefix, real t2)
{
	char cmsgnm[256];
	if (t2 < (real)600.)
		sprintf(cmsgnm, "%s cputime to complete %s = %8.2f cpu-sec", prefix, label, t2);
	else
	{
		if ((t2 >= (real)600.) && (t2 < (real)3600.))
			sprintf(cmsgnm, "%s cputime to complete %s = %8.2f cpu-min", prefix, label, t2/60.);
		else
			sprintf(cmsgnm, "%s cputime to complete %s = %8.2f cpu-hr", prefix, label, t2/3600.);
	}
	Wrimsg(label, cmsgnm);
}

void GreenFunctionManager::Subsystem::PrintTimeForecasts(const char *label, real py, real pz)
{
	char cmsgnm[256];
	if (py * pz > (real)0.)
		sprintf(cmsgnm, "[cputime scales as 1/gamma^2]");
	else
	{
		if (py + pz == (real)0.)
			sprintf(cmsgnm, "[cputime is independent of gamma]");
		else
			sprintf(cmsgnm, "[cputime scales as 1/gamma]");
	}
	Wrimsg(label, cmsgnm);
}

bool GreenFunctionManager::AllocateMagneticBuffer(void)
{
	return subSys[SubsystemTypeMagnetic]->Allocate();
}

void GreenFunctionManager::DeallocateMagneticBuffer(void)
{
	subSys[SubsystemTypeMagnetic]->Deallocate();
}

bool GreenFunctionManager::AllocateElectricBuffer(void)
{
	return subSys[SubsystemTypeElectric]->Allocate();
}

void GreenFunctionManager::DeallocateElectricBuffer(void)
{
	subSys[SubsystemTypeElectric]->Deallocate();
}
