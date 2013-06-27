#include "StdAfx.h"

#include "Targspher.h"
#include "TargetManager.h"

/* **
Routine to construct irregular grain using spherical harmonics to define surface.

Input:
        AEFF  =effective radius/d (d=lattice spacing)
        PRINAX=0 to use (1,0,0) and (0,1,0) in TF as a1,a2
              =1 to use principal axes of moment of inertia tensor for a1 and a2
        BETA  = power law index
                if BETA = 0 : read ALM values from file CFLSHP
                if BETA > 1 : generate ALM values here using random number generator, with
                <|a_LM|^2> propto [1/(L+1)]*L**(-BETA) for L>1

        RLMAX = maximum L value to use in Y_LM expansion LMAX=NINT(RLMAX)
        RSEED = to choose seed for random number generator
              = 1, 2, 3, 4, ... ISEED=-NINT(RSEED) is passed to gaussian deviate program
        S2 = expectation value of <s^2>, where < > denotes average over grain surface

        CFLSHP=name of file containing spherical harmonic amplitudes (used only if BETA = 0)
        IOSHP=device number for "target.out" file =-1 to suppress printing of "target.out"
        MXNAT=dimensioning information (max number of atoms)

Output:
       A1(1-3)=unit vector (1,0,0) defining target axis 1
       A2(1-3)=unit vector (0,1,0) defining target axis 2
       NAT=number of atoms in target
       IXYZ(1-NAT,1-3)=(x-x0(1))/d,(y-x0(2))/d,(z-x0(3))/d for atoms of target
       X0(1-3)=(location/d) in Target Frame corresponding to dipole with
               IXYZ=(0,0,0).  This will be treated at the origin of physical coordinates in the TF.
               Here origin is set to be centroid of the target.
       CDESCR=string describing target (up to 67 chars.)

=======================================================================
                  Gaussian Sphere Geometry

The target is defined in polar coordinates by a radius function

                               Lmax  l
   R(theta,phi).= const * exp[ sum  sum  a_lm Y_lm(theta,phi) ]
                               l=1  m=-l

   We assume that the a_lm are gaussian random variables, with

                  4*pi   
   <|a_{1M}|^2> = ---- * f1* <s^2>
                   3     

                   4*pi     (1-f1)*<s^2>     1
   <|a_{LM}|^2> = ------- * ------------ * ------
                  (2*L+1)   zeta(beta)-1   L^beta

After generating the a_lm values according to above prescription for
L=1,...,NMAX, we calculate the actual <s^2> averaged over angles.
Because the a_lm have been chosen stochastically, this will not be
exactly equal to the input S2.
We then rescale all coefficients a_lm by a fixed factor so that
the resulting a_lm have angle-averaged <s^2> = S2 (summing over
L values up to NMAX).
For reference, we also calculate <s^2> with the sum restricted to
L values up to LMAX (as is used in the target generation).

For s to be real it is also required that

   a_{L,-M} = (-1)^m conjg(a_{L,M})    ---> Im[a_{L0}] = 0

Because of the above requirement, we do not calculate or store
the a_lm for m < 0, since all that we require to calculate the
target shape are the combination

   a_{L,m}Y_{L,m} + a_{L,-m}Y_{L,-m}

                 (2L+1)*(L-m)!
         = sqrt[ ------------- ] * P_Lm * 
                   4*pi*(L+m)!

             [ a_Lm e^{i*m*phi} + (-1)^m a_{L,-m}*e^{-i*m*phi) ]

where P_Lm(theta,phi) is the associated Legendre function.
With above requirement on a_{L,-M} it follows that

                 (2L+1)*(L-m)!
         = sqrt[ ------------- ] * P_Lm * 2 Re(a_Lm e^{i*m*phi})
                   4*pi*(L+m)!

History:
Fortran history records removed.
End history.

Copyright (C) 2002,2003,2004,2007,2008 B.T. Draine and P.J. Flatau
Copyright (c) 2012, C++ version, V.Choliy

This code is covered by the GNU General Public License.
** */

Targspher::Targspher(TargetManager *man) : LoadableTarget(man) 
{ 
	shortDescr = string("Targspher");
	longDescr = string("Gaussian sphere");
	nres = 20;
	lmax = nint_(shpar[3]);
}

void Targspher::Sizer(void)					// Sample surface to determine total extent of target in x,y,z directions
{
	Presizer();
//
	const real Twopi = twox_ * Pi;
	real xmin, xmax, ymin, ymax, zmin, zmax, rad, xx, yy, zz;
	xmin = xmax = ymin = ymax = zmin = zmax = zero_;
	for(int jth=0; jth<nres*lmax; ++jth)
	{
		real theta = Pi * ((real)jth - half_) / (real)(nres * lmax);
		real costh = Cos(theta);
		real sinth = Sin(theta);

		for(int jphi=0; jphi<nres*lmax; ++jphi)
		{
			real phi = Twopi * (real)jphi / (real)(nres*lmax);
			rad = shpar[0] * Rgspher(costh, phi, lmax, lmax, shpar[4], alm);
			xx = costh * rad;
			yy = sinth * Cos(phi) * rad;
			zz = sinth * Sin(phi) * rad;
			if (xx < xmin) xmin = xx;
			if (xx > xmax) xmax = xx;
			if (yy < ymin) ymin = yy;
			if (yy > ymax) ymax = yy;
			if (zz < zmin) zmin = zz;
			if (zz > zmax) zmax = zz;
		}
	}
//
	minJx = -nint_(xmin) - 1;
	minJy = -nint_(ymin) - 1;
	minJz = -nint_(zmin) - 1;
	maxJx =  nint_(xmax) + 1;
	maxJy =  nint_(ymax) + 1;
	maxJz =  nint_(zmax) + 1;
}

void Targspher::Descriptor(void)
{
	sprintf(freeDescr, " Gaussian sphere target of NAT=%7d dipoles", nat);
}

void Targspher::Vector(void)
{
	if (shpar[1] != 0)
		eigval = Prinaxis();				// Call PRINAXIS to compute principal axes (A1, A2 changed)
	else
		VectorA();							// Do not compute principal axes: simply set, a1=(1,0,0) and a2=(0,1,0) in target frame.
}

void Targspher::VectorX(void)					// Set X0 so that origin is at the centroid
{
	x0.Clear();
	for(int jx=0; jx<nat0; ++jx)
	{
		x0 += Vect3<real>((real)ixyz.Value(jx, 0), (real)ixyz.Value(jx, 1), (real)ixyz.Value(jx, 2));
	}
	x0 /= -nat;
}

void Targspher::Allocator(void)
{
	int curSize = nx * ny;
	ixyz.Dimension(curSize, 3);
	
	nat0 = 0;
	for(int jx=minJx; jx<=maxJx; ++jx)
	{
		real x = (real)jx;
		for(int jy=minJy; jy<=maxJy; ++jy)
		{
			real y = (real)jy;
			for(int jz=minJz; jz<=maxJz; ++jz)
			{
                real cosph, sinph, costh, phi;
				real z = (real)jz;
				real ryz2 = y*y + z*z;
				real ryz = Sqrt(ryz2);
				if (ryz2 > zero_)
				{
					cosph = y / ryz;
					sinph = z / ryz;
				}
				else
				{
					cosph = onex_;
					sinph = zero_;
				}
				real r = Sqrt(ryz2 + x*x);
				if (r > zero_)
					costh = x/r;
				else
					costh = onex_;
				if (cosph >= zero_)
					phi = Asin(sinph);
				else
					phi = Pi - Asin(sinph);
				real rad = shpar[0] * Rgspher(costh, phi, lmax, lmax, shpar[4], alm);
				if (rad >= r)
				{
					if (nat0 >= curSize)
					{
						ixyz.Extend(nz);
						curSize = ixyz.GetSize(0);
					}
                    ixyz.Fill3(nat0, jx, jy, jz);
					int index = GetLinearAddress(nat0);
					Composer(index, 0);									// Homogeneous target:
					iocc[index] = true;
					++nat0;
				}
			}
		}
	}
	ixyz.Close(nat0);
}

void Targspher::Presizer(void)
{
	const real Twopi = Pi + Pi;
	const real Fourpi = Twopi + Twopi;
	const real Roottwo = Sqrt(twox_);

	alm.Dimension(lmax, lmax+1);

	bool bOk = false;
	if (shpar[2] == zero_)
	{
		if (AlmReader() == false)
		{
			Errmsg("Fatal", "Targspher", " Almreader failed.");
		}
//
// Calculate S2
		real sum = zero_;
		for(int l=0; l<lmax; ++l)
		{
			sum += alm.Value(l, 0).modSquared();
			for(int m=0; m<l; ++m)
			{
				sum += twox_ * alm.Value(l, m).modSquared();
			}
		}
		s2lmax = sum / Fourpi;
		bOk = true;
	}
//
	if (shpar[2] > onex_)
	{
		int l, m;
//
// Set value of F1 = fraction of <s^2> contributed by L=1
// note that to leading order the L=1 terms produce merely a
// translation with no deformation, so we normally suppress the
// L=1 terms by setting F1=0.0
		f1 = zero_;
//
// Assign a_lm values using gaussian random deviates set seed
		int iseed = -nint_(shpar[5]);
//
// Calculate random numbers and store in appropriate elements of ALM
// Calculate them in appropriate order so that changes in NMAX will leave lower L values unaffected.  
// Specifying RSEED, F1, and BETA therefore serves to define a shape, with adjustment to LMAX allowing one to see important to the shape of high L components.
//
// Initialize GASDEV
		real ar = Gasdev(iseed);
		real ai = zero_;
		for(l=1; l<=lmax; ++l)
		{
			alm.Value(l, 0) = Complex(Gasdev(iseed), zero_);
			for(m=1; m<=l; ++m)
			{
				ar = Gasdev(iseed) / Roottwo;
				ai = Gasdev(iseed) / Roottwo;
				alm.Value(l, m) = Complex(ar, ai);
			}
		}
//
// ALM(L,M) now contains random Complex numbers a_lm , with
//    <a_lm> = 0 
//    <|a_lm|^2> = <[Re(a_lm)]^2> + <[Im(a_lm)]^2> = 1
//    Im[a_lm(L,0)] = 0
//    <[Re(a_lm)]^2> = <[Im(a_lm)]^2> = 1/2 for M > 0
//
// Now multiply the A(L,M) random deviates by appropriate factors to obtain intended <|a_lm|^2> values
//
// First treal L=1
// SIGLM = <|a_lm|^2>^{1/2} 
		real siglm = Sqrt(f1 * shpar[4] * Fourpi / (real)3.);
		alm.Value(1, 0) *= siglm;
		alm.Value(1, 1) *= siglm;
//
// set coefficients for L > 1:
// estimate SUM = [Riemann zeta function - 1]
		real sum = zero_;
		for(l=lmax; l>=2; --l)
		{
			sum += onex_ / Pow((real)l, shpar[2]);
		}
		real term = (onex_ - f1) * shpar[4] * Fourpi / sum;
//
// Note: we calculate coefficients up to NMAX
// We then use values only up to LMAX to generate shape, but we want to have correct amount of power in high order coefficients, even if they are later suppressed.
		for(l=2; l<=lmax; ++l)
		{
			siglm = Sqrt(term / ((real)(2*l+1) * Pow((real)l, shpar[2])));
			for(m=0; m<=l; ++m)
			{
				alm.Value(l, m) *= siglm;
			}
		}
//
// Compute actual <s^2> for this realization
		sum = zero_;
		for(l=0; l<=lmax; ++l)
		{
			sum += alm.Value(l, 0).modSquared();
			for(m=1; m<=l; ++m)
			{
				sum += alm.Value(l, m).modSquared() * twox_;
			}
		}
		sum /= Fourpi;
//
// At this point, we rescale the a_lm values so that we will have
// <s^2> = S2 exactly including all terms up to NMAX
		term = Sqrt(shpar[4] / sum);
		for(l=0; l<=lmax; ++l)
		{
			for(m=1; m<=l; ++m)
			{
				alm.Value(l, m) *= term;
			}
		}
//
// Now calculate S2LMAX = <s^2> contributed by terms up to LMAX
		sum = zero_;
		for(l=0; l<=lmax; ++l)
		{
			sum += alm.Value(l, 0).modSquared();
			for(m=1; m<=l; ++m)
			{
				sum += alm.Value(l, m).modSquared() * twox_;
			}
		}
		s2lmax = sum / Fourpi;
		bOk = true;
	}

	if (bOk == false)
		Errmsg("Fatal", "Targspher", " Invalid BETA");

//
// ALM values are now defined.  
// Proceed to generate target array, using a_lm values up to Lmax
	a1.Set(onex_, zero_, zero_);
	a2.Set(zero_, onex_, zero_);
//
// if IOSHP > 0 , calculate three cross-sectional slices
//
// Polar slice with phi=0
	if (manager->Ioshp())
	{
        int jth;
		Vect3<real> *phis  = new Vect3<real>[nres * lmax];
		Vect3<real> *rads  = new Vect3<real>[nres * lmax];
		Vect3<real> *thets = new Vect3<real>[nres * lmax];
        real theta, costh;
		real phi = zero_;
        for(jth=0; jth<=nres*lmax; ++jth)
		{
			theta = Twopi * (real)jth / ((real)(nres * lmax));
			costh = Cos(theta);
            thets[jth].data[0] = theta;
            phis[jth].data[0] = phi;
            real phiarg = phi;
            if (theta > Pi)
				phiarg = phi + Pi;
            rads[jth].data[0] = Rgspher(costh, phiarg, lmax, lmax, shpar[4], alm);
		}
//
// Polar slice with phi=pi/2
		phi = Pi / twox_;
		for(jth=0; jth<=nres*lmax; ++jth)
		{
			theta = Twopi * (real)jth / ((real)(nres * lmax));
			costh = Cos(theta);
            thets[jth].data[1] = theta;
            phis[jth].data[1] = phi;
            real phiarg = phi;
            if (theta > Pi)
				phiarg = phi + Pi;
            rads[jth].data[1] = Rgspher(costh, phiarg, lmax, lmax, shpar[4], alm);
		}
//
// Equatorial slice
		theta = Pi / twox_;
		costh = zero_;
		for(int jphi=0; jphi<nres*lmax; ++jphi)
		{
			phi = Twopi * (real)jphi / ((real)(nres * lmax));
			thets[jphi].data[2] = theta;
			costh = Cos(theta);
            phis[jphi].data[2] = phi;
            rads[jphi].data[2] = Rgspher(costh, phi, lmax, lmax, shpar[4], alm);
		}
//
		const char *Format7600 = "three cross sectional slices for spherical harmonic grain\n%3d = LMAX (highest order for Y_lm)\n------ slice 1 ------ ------ slice 2 ------ ------ slice 3 ------\n theta  phi   r/a_eff  theta  phi   r/a_eff  theta  phi   r/a_eff\n";
		const char *Format7700 = "%6.4lf%7.4lf%8.4lf%7.4lf%7.4lf%8.4lf%7.4lf%7.4lf%8.4lf";
		FILE *ioshpFile = fopen("spharm.out", "w");
		fprintf(ioshpFile, Format7600, lmax);
		for(jth=0; jth<nres*lmax; ++jth)
		{
			fprintf(ioshpFile, Format7700, thets[jth].data[0], phis[jth].data[0], rads[jth].data[0], 
				thets[jth].data[1], phis[jth].data[1], rads[jth].data[1], thets[jth].data[2], phis[jth].data[2], rads[jth].data[2]);
		}
		fclose(ioshpFile);

		delete [] phis;
		delete [] rads;
		delete [] thets;
	}
}

bool Targspher::AlmReader(void)
{
	char Buffer[256];

	FILE *ioshpFile = fopen(cflshp, "r");
	if (!ioshpFile)
		return false;
    char Format[256];
	fgets(Buffer, 255, ioshpFile);
    sprintf(Format, "%s%s", realFormat, realFormat);
    sscanf(Buffer, Format, &s2lmax, &shpar[4]);
	fgets(Buffer, 255, ioshpFile);
	fgets(Buffer, 255, ioshpFile);
	fgets(Buffer, 255, ioshpFile);
	fgets(Buffer, 255, ioshpFile);
	sscanf(Buffer, "%d", &lmax);
	fgets(Buffer, 255, ioshpFile);
//
// Read ALM values for nonnegative m
// LMAX = 1 : 2    (m=0 and 1)
// LMAX = 2 : 2+3 = 5
// LMAX = 3 : 2+3+4 = 9
// general: (LMAX^2+3*LMAX)/2
	real ar, ai;
	int l, m;
	int jth = (lmax + 3) * lmax / 2;
    sprintf(Format, "%%d%%d%s%s", realFormat, realFormat);
	for(int ja=0; ja<jth; ++ja)
	{
		fgets(Buffer, 255, ioshpFile);
        sscanf(Buffer, Format, &l, &m, &ar, &ai);
		alm.Value(l-1, m-1) = Complex(ar, ai);
	}
	fclose(ioshpFile);

	return true;
}

void Targspher::Printer(void)								// Technical output
{
	if (manager->Ioshp())
	{
		if (shpar[2] > 0)					// if this is a new shape, write the ALM value out to a file
		{
			FILE *ioshpFile2 = fopen("targspher.alm", "w");
			if (ioshpFile2)
			{
				const char *Format9100 = "%7.5lf%7.4lf%6.3lf%5.2lf%3d%4d = <s^2>_lmax, <s^2>, f_1, shpar[2], Lmax, ISEED for gaussian sphere\n%7d%8.4lf%8.4lf%8.4lf%8.4lf%8.4lf = NAT, a_eff/d, alpha_1-3\n";
                fprintf(ioshpFile2, Format9100, s2lmax, shpar[4], f1, shpar[2], lmax, nint_(shpar[5]), nat, shpar[0], Pow((real)(0.75 * nat / Pi), (real)(1./3.)), alpha.data[0], alpha.data[1], alpha.data[2]);
				a1.Fprintf(ioshpFile2, "%10.6lf", NULL, " = A_1 vector\n");
				a2.Fprintf(ioshpFile2, "%10.6lf", NULL, " = A_2 vector\n");
				fprintf(ioshpFile2, "%4d = Lmax\n L   M   Re(a_LM)  Im(a_LM)\n", lmax);
				for(int l=0; l<lmax; ++l)
				{
					for(int m=0; m<l; ++m)
					{
                        fprintf(ioshpFile2, "%3d%4d%10.6lf%10.6lf", l, m, alm.Value(l, m).re, alm.Value(l, m).im);
					}
				}
				fclose(ioshpFile2);
			}
			else
				Wrimsg("Targspher", "Cannot open beta - ioshp file.\n");
		}
//
		FILE *ioshpFile = fopen("target.out", "w");
		if (ioshpFile)
		{
			const char *Format9020 = "> TARGSPHER:%7.5lf%7.4lf%6.3lf%5.2lf%3d%4d = <s^2>_Lmax, <s^2>, f_1, shpar[2], Lmax, ISEED for gaussian sphere\n%10d%8.4lf%8.4lf%8.4lf%8.4lf%8.4lf = NAT, AEFF(in) a_eff/d, alpha_1-3\n";
			fprintf(ioshpFile, Format9020, s2lmax, shpar[4], f1, shpar[2], lmax, nint_(shpar[5]), nat, shpar[0], Pow((real)(0.75 * nat / Pi), (real)(1./3.)), alpha.data[0], alpha.data[1], alpha.data[2]); 
			OutShpar(ioshpFile);
			fprintf(ioshpFile, "%10d = NAT\n", nat0);
			OutVectorsAndArrays(ioshpFile);
			fclose(ioshpFile);
		}
		else
			fprintf(stderr, "Cannot open ioshpFile in Targspher::Printer for %s.", shortDescr.c_str());
	}
}

real Targspher::Rgspher(real costh, real phi, int lmax, int nmax, real s2, Array2F<Complex> &alm)
{
/* **
Given:
   COSTH = cos(theta)
   PHI   = phi (radians)
   LMAX = maximum order of spherical harmonics
   S2   = <s^2>
   ALM  = Complex array of spherical harmonic amplitudes a_lm
          L=1,...,LMAX , M=-0,..,L
          Note that we do not require the ALM for m<0, since these
          are obtained from the values for m>0 using
          a_{l,-m} = (-1)^m conjg(a_{l,m})

Returns
   RGSPHER = distance from "center" to surface for target with
             a_eff=1.
                   1             Lmax  l
           = ------------ * exp[ sum  sum a_lm Y_lm ]
             sqrt(1+<s^2>)       l=1  m=-l

Requires:
External function P_LM(L,M,X) to return associated Legendre polynomial at the moment this function resides in file tarspharm.f

Note: we assume that a_{l,-m} = (-1)^m conjg(a_{l,m}) so that sum over spherical harmonics gives real function

Fortran history records removed.

Copyright (C) 2002,2003 B.T. Draine and P.J. Flatau
Copyright (c) 2012, C++ version, V.Choliy

This code is covered by the GNU General Public License.
** */
	const real Fourpi = four_ * Pi;

	real sum = zero_;
	for(int l=0; l<lmax; ++l)
	{
		real fac = Sqrt((real)(2*l + 1) / Fourpi);
		real y_l0 = fac * P_lm(l+1, 0, costh);
		sum += alm.Value(1, 0).re * y_l0;

		real sign = onex_;
		for(int m=0; m<l; ++m)
		{
			fac /= Sqrt((real)((l+m) * (l+1-m)));
			sign = -sign;
			Complex cxterm = alm.Value(l, m) * Complex(Cos(m * phi), Sin(m * phi));
			sum += fac * P_lm(l+1, m+1, costh) * twox_ * cxterm.re;
		}
	}

	return (Exp(sum) / Sqrt(onex_ + s2));
}

real Targspher::P_lm(int l, int m, real x)
{
/* **
FUNCTION P_LM  = Associated Legendre polynomial P_{lm}(x=cos(theta))
given:
   L,M
   X
returns:
   P_LM = P_lm(x)
            associated Legendre polynomial
this routine follows closely the structure of function PLGNDR 
described in section 6.6 of "Numerical Recipes", by Press, W.H.,
Flannery, B.P., Teukolsky, S.A., and Vetterling, W.T. (1986, Cambridge
Univ. Press).

Modifications by B.T. Draine, Princeton Univ. Observatory
history
02.11.12 (BTD) first written
04.03.31 (BTD) changed WRITE(0 to CALL WRIMSG(
07.08.06 (BTD) converted to f90
end history
** */
	real pmm = onex_;
	if ((m > 0) && (m < l))
	{
		real sinth = Sqrt(onex_ - x*x);
		real fact = onex_;
		for(int i=0; i<m; ++i)
		{
			pmm *= (fact*sinth);
			fact += twox_;
		}
	}
	else
	{
		char Buffer[256];
		sprintf(Buffer, " fatal error in p_lm: m=%d < 0", m);
		Wrimsg(" p_lm ", Buffer);
		return 0; //exit(0);
	}
//
	if (l == m)
		return pmm;
	else
	{
		real pmmp1 = x * (real)(2*m + 1) * pmm;
		if (l == m+1)
			return pmmp1;
		else
		{
			real pll = pmmp1;
			for(int ll=m+2; ll<=l; ++ll)
			{
				pll = (x * (real)(2*ll - 1) * pmmp1 - (real)(ll + m - 1) * pmm) / (real)(ll - m);
				pmm = pmmp1;
				pmmp1 = pll;
			}
			return pll;
		}
	}
}

const char *TargetVerboseDescriptor_Gausssph(int num)
{

	return NULL;
}

REGISTER_TARGET(Gausssph,6,false,-1,0," -- gauss spheres -- ")
void Target_Gausssph::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "a_eff/d:\n");
	fprintf(stream, "%20.16lf\n", shpar[0]);
	fprintf(stream, "power law index beta > 1 (0 to input a_lm from file targspher.alm)\n");
	fprintf(stream, "%20.16lf\n", shpar[1]);
	if (shpar[2] == zero_)
	{
		fprintf(stream, "OK -- will read input params from targspher.alm\n");
		strcpy(cflshp, "targspher.alm");
	}
	else
	{
		if (shpar[2] > onex_)
		{
			fprintf(stream, "Lmax > 0 (e.g., 5)\n");
			fprintf(stream, "%20.16lf\n", shpar[3]);
			fprintf(stream, "<s^2> (e.g., 0.05) [this includes all a_lm up to NMAX (=100)\n");
			fprintf(stream, "%20.16lf\n", shpar[4]);
			fprintf(stream, "int seed > 0\n");
			fprintf(stream, "%20.16lf\n", shpar[5]);
		}
		else
		{
			fprintf(stream, "Invalid value for BETA\n");
		}
	}
	fprintf(stream, "PRINAX=0 for a_1=(1,0,0),a_2=(0,1,0) in TF, 1 for a_1,a_2 = principal axes\n");
	fprintf(stream, "%20.16lf\n", shpar[1]);
}

real Targspher::Gasdev(int idum)
{
/* **
Given:

   IDUM = seed for random number generator
          if IDUM < 0, reinitialize RAN1 with seed IDUM
          if IDUM .ge. 0, continue sequence from RAN1
Returns:

  GASDEV(IDUM) = gaussian random deviate with mean <gasdev>= 0 
                 and variance <gasdev**2>=1

This subroutine follows closely the structure of the routine gasdev
described in Numerical Recipes by Press, W.H., Flannery, B.P.,
Teukolsky, S.A., and Vetterling, W.T. (Cambridge Univ. Press).

NR version uses RAN1, but it was found that under RH7.1/pgf77 
gasdev does not have <gasdev> = 0 when we use RAN1.
Therefore have changed to use RAN3

Modifications by B.T. Draine, Princeton University Observatory,
2003.09.10
2007.08.06 (BTD) converted to f90
** */

	static bool iset = false;
	static real gset; 
//
	if (iset == false)
	{
		real v1, v2, rsq;
// when ISET = 0, we do not have an extra deviate on hand so generate two new deviates
		do						// generate random variates V1 and V2 on interval -1,1
		{
			v1 = twox_ * Ran3(idum) - onex_;
			v2 = twox_ * Ran3(idum) - onex_;
			rsq = v1 * v1 + v2 * v2;
		} while ((rsq >= onex_) || (rsq == zero_));
		real fac = Sqrt(-twox_ * Log(rsq) / rsq);
        gset = v1 * fac;
		iset = true;
        return v2 * fac;
	}
	else
	{
// when ISET = 1, we have an extra deviate on hand so use it, and change ISET back to 0
		iset = false;
		return gset;
	}
}
