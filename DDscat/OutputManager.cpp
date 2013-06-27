#include "StdAfx.h"

#include "OutputManager.h"

OutputManager *OutputManager::item = NULL;
OutputManager::OutputManager(void)
{
	currentTarget = NULL;
	param = NULL;
	jpbc = PeriodicEnd;
	betad = phid = thetad = (real)0.;
	file8 = file10 = file11 = file12 = NULL;
}

OutputManager::~OutputManager(void)
{
//	fclose(file8);
//	fclose(file10);
	fclose(file11);
//	fclose(file12);
}

OutputManager *OutputManager::GetInstance(void)
{
	if (!item)
		item = new OutputManager;

	return item;
}

void OutputManager::Kill(void)
{

}

void OutputManager::Init(AbstractTarget *currentTarget, DDscatParameters *param, DielectricManager *dielec, PeriodicBoundaryFlag jpbc)
{
	this->currentTarget = currentTarget;
	this->param = param;
	this->dielec = dielec;
	this->jpbc = jpbc;

	betmid = Degrad * param->GetOriData()->Betami();
	betmxd = Degrad * param->GetOriData()->Betamx();
	phimid = Degrad * param->GetOriData()->Phimin();
	phimxd = Degrad * param->GetOriData()->Phimax();
	thtmid = Degrad * param->GetOriData()->Thetmi();
	thtmxd = Degrad * param->GetOriData()->Thetmx();
}

bool OutputManager::Writefml(int iori, int irad, int iwav, int navg, int ncomp, char *cdescr, real aeff, real ak1, Vect3<real> &akr, 
	real wave, real xx, Vect3<Complex> &cxe01r, Vect3<Complex> &cxe02r, FourArray *cxfData)
{
/* **
Purpose of this module is to write out the Complex scattering amplitudes f_ml

Fortran versions history removed.

Copyright (C) 2007,2008,2011, B.T. Draine and P.J. Flatau
Copyright (C) 2013, C++ versions, Choliy V.
This code is covered by the GNU General Public License.
** */
//
// If IORI > 0, then write out scattering properties for this specific target orientation.
	std::string cframe = std::string(FrameEnumeratorVerbose(param->Cmdfrm()));
//
	int j;
	file8 = fopen(FileNamer::GetInstance()->GetCflfml(), "w");
	if (file8 == NULL)
	{
		Wrimsg("Writefml", "Cannot open file8, line 55");
		return false;
	}
	const char *cstamp = Version();
	const char *Format9030 = " DDSCAT --- %s\n TARGET ---%s\n %s --- method of solution \n %s --- prescription for polarizabilies\n%s --- shape\n%7d = NAT0 = number of dipoles\n";
	fprintf(file8, Format9030, cstamp, cdescr, SolEnumerator(param->Cmdsol()), AlphaEnumerator(param->Calpha()), param->TargetString(), currentTarget->Nat0());
	const char *Format9032 = "\n  AEFF=%11.5lf = effective radius (physical units)\n  WAVE=%11.5lf = vacuum wavelength (physical units)\nK*AEFF=%11.5lf = 2*pi*aeff/lambda\nNAMBIENT=%9.5lf = refractive index of ambient medium\n";
	fprintf(file8, Format9032, aeff, wave, xx, param->Nambient());
	for(j=0; j<ncomp; ++j)
	{
		const char *Format9031 = "n= (%7.4lf , %7.4lf),  eps.= (%8.4lf , %7.4lf)  |m|kd=%8.4lf for subs.%2d\n";
		real mkd = Pow((real)(4./3.) * Pi / currentTarget->Nat0(), (real)(1./3.)) * dielec->GetCxrfr(j).mod() * xx;
		fprintf(file8, Format9031, dielec->GetCxrfr(j).re, dielec->GetCxrfr(j).im, dielec->GetCxeps(j).re, dielec->GetCxeps(j).im, mkd, j);
	}
	const char *Format90331 = "   TOL=%10.3e = error tolerance for CCG method\n  NAVG=%6d = (theta,phi) values used in comp. of Qsca,g\n";
	fprintf(file8, Format90331, param->Tol(), navg);
	currentTarget->A1().Fprintf(file8, "%8.5lf", "(", ") = target axis A1 in Target Frame\n");
	currentTarget->A2().Fprintf(file8, "%8.5lf", "(", ") = target axis A2 in Target Frame\n");
	const char *Format9035 = " (%8.5lf%9.5lf%9.5lf) = k vector (latt. units) in TF\n (%8.5lf,%8.5lf)(%8.5lf,%8.5lf)(%8.5lf,%8.5lf)=inc.pol.vec. 1 in TF\n (%8.5lf,%8.5lf)(%8.5lf,%8.5lf)(%8.5lf,%8.5lf)=inc.pol.vec. 2 in TF\n";
	fprintf(file8, Format9035, akr.data[0], akr.data[1], akr.data[2], 
		cxe01r.data[0].re, cxe01r.data[0].im, cxe01r.data[1].re, cxe01r.data[1].im, cxe01r.data[2].re, cxe01r.data[2].im, 
		cxe02r.data[0].re, cxe02r.data[0].im, cxe02r.data[1].re, cxe02r.data[1].im, cxe02r.data[2].re, cxe02r.data[2].im);
	real rvar1 = ak1;
	real rvar2 = (real)0.;
	real rvar3 = (real)0.;
	const char *Format9037 = " (%8.5lf%9.5lf%9.5lf) = k vector (latt. units) in Lab Frame\n (%8.5lf,%8.5lf)(%8.5lf,%8.5lf)(%8.5lf,%8.5lf)=inc.pol.vec. 1 in LF\n (%8.5lf,%8.5lf)(%8.5lf,%8.5lf)(%8.5lf,%8.5lf)=inc.pol.vec. 2 in LF\n";
	fprintf(file8, Format9037, rvar1, rvar2, rvar3, 
		param->Cxe01_lf().data[0].re, param->Cxe01_lf().data[0].im, param->Cxe01_lf().data[1].re, param->Cxe01_lf().data[1].im, 
		param->Cxe01_lf().data[2].re, param->Cxe01_lf().data[2].im, param->Cxe02_lf().data[0].re, param->Cxe02_lf().data[0].im, 
		param->Cxe02_lf().data[1].re, param->Cxe02_lf().data[1].im, param->Cxe02_lf().data[2].re, param->Cxe02_lf().data[2].im);
   	const char *Format9040 = " BETA =%7.3lf = rotation of target around A1\n THETA=%7.3lf = angle between A1 and k\n  PHI =%7.3lf = rotation of A1 around k\n";
	fprintf(file8, Format9040, betad, thetad, phid);
//
	real pyd = currentTarget->Pyd();
	real pzd = currentTarget->Pzd();
	Complex cxfac;
	real sinalpha;
	switch(jpbc)
	{
	case PeriodicNo:
		{
			cxfac.unityRe();
			const char *Format9050 = "     Finite target:\n     e_m dot E(r) = i*exp(ikr)*f_ml*E_inc(0)/(kr)\n     m=1 in scatt. plane, m=2 perp to scatt. plane\n\n theta   phi  Re(f_11)   Im(f_11)   Re(f_21)   Im(f_21)   Re(f_12)   Im(f_12)   Re(f_22)   Im(f_22)\n";
			fprintf(file8, Format9050);
		}
		break;
//
// for periodic structures (JPC=1,2, or 3), it is assumed that 
// target axis a_1 = x_TF
//             a_2 = y_TF
//             a_3 = z_TF
	case PeriodicY:							// JPBC=1 : target is periodic in y direction
		{
			sinalpha = Sin(thetad / Degrad) * Cos(betad / Degrad);
			sinalpha = (real)1. - sinalpha * sinalpha;
			sinalpha = Sqrt(sinalpha);
			cxfac = (Complex((real)0., TwoPi / sinalpha)).sqrt() / (ak1 * pyd);
			rvar1 = ak1 * pyd;
			const char *Format9051 = "     Target periodic in y direction with L_y/d=%10.3e  kL_y = %10.3e\n     e_m dot E(r) = exp(ikr)*f_ml*[e_l dot E_inc(0)]/sqrt(kR)     [R = dist. from target]\n     m=1 in scatt. plane, m=2 perp to scatt. plane\n\n alpha  zeta  Re(f_11)   Im(f_11)   Re(f_21)   Im(f_21)   Re(f_12)   Im(f_12)   Re(f_22)   Im(f_22)\n";
			fprintf(file8, Format9051, pyd, rvar1);
		}
		break; 

	case PeriodicZ:							//      2                         z direction
		{
			sinalpha = Sin(thetad / Degrad) * Sin(betad / Degrad);
			sinalpha = (real)1. - sinalpha * sinalpha;
            sinalpha = Sqrt(sinalpha);
			cxfac = (Complex((real)0., TwoPi / sinalpha)).sqrt() /(ak1 * pzd);
			rvar1 = ak1 * pzd;
			const char *Format9052 = "     Target periodic in z direction with L_z/d=%10.3e  kL_z = %10.3e\n     e_m dot E(r) = exp(ikr)*f_ml*[e_l dot E_inc(0)]/sqrt(kR)     [R = dist. from target]\n     m=1 in scatt. plane, m=2 perp to scatt. plane\n\n alpha  zeta  Re(f_11)   Im(f_11)   Re(f_21)   Im(f_21)   Re(f_12)   Im(f_12)   Re(f_22)   Im(f_22)\n";
			fprintf(file8, Format9052, pzd, rvar1);
		}
		break; 

	case PeriodicBoth:							//      3                         y and z directions
		{
// 080630 BTD not certain how we want to set normalization factor CXFAC
//            THETAD = angle relative to normal
//            ALPHA  = 90-THETAD
//                   = angle between scattered vector and surface
			sinalpha = Cos(thetad / Degrad);
			cxfac.set((real)0., TwoPi / (ak1 * pyd * ak1 * pzd * sinalpha));
			rvar1 = ak1 * pyd;
			rvar2 = ak1 * pzd;
			const char *Format9053 = "     Target periodic in y,z directions with L_y/d,L_z/d=%10.3e%10.3e  kL_y=%10.3e kL_z=%10.3e\n      e_m dot E(r) = exp(ikr)*f_ml*[e_l dot E_inc(0)]\n     m=1 in scatt. plane, m=2 perp to scatt. plane\n\n alpha   phi  Re(f_11)   Im(f_11)   Re(f_21)   Im(f_21)   Re(f_12)   Im(f_12)   Re(f_22)   Im(f_22)\n";
			fprintf(file8, Format9053, pyd, pzd, rvar1, rvar2);
		}
		break;

	default:
		break;
	}

	for(int nd=0; nd<param->Nscat(); ++nd)
	{
		const char *Format9070 = "%6.1lf%6.1lf%11.3e%11.3e%11.3e%11.3e\n";
		const char *Format9071 = "%6.1lf%6.1lf%11.3e%11.3e%11.3e%11.3e%11.3e%11.3e%11.3e%11.3e\n";
		real phind = Degrad * param->Phin()[nd];
		real thetnd = Degrad * param->Thetan()[nd];
		if (param->Iorth() == 1)
		{
			fprintf(file8, Format9070, thetnd, phind, 
				(cxfac * cxfData->Cx11(nd)).re, (cxfac * cxfData->Cx11(nd)).im, (cxfac * cxfData->Cx21(nd)).re, (cxfac * cxfData->Cx21(nd)).im);
		}
		else
		{
			fprintf(file8, Format9071, thetnd, phind, 
				(cxfac * cxfData->Cx11(nd)).re, (cxfac * cxfData->Cx11(nd)).im, (cxfac * cxfData->Cx21(nd)).re, (cxfac * cxfData->Cx21(nd)).im, 
				(cxfac * cxfData->Cx12(nd)).re, (cxfac * cxfData->Cx12(nd)).im, (cxfac * cxfData->Cx22(nd)).re, (cxfac * cxfData->Cx22(nd)).im);
		}
	}
	fclose(file8);
	return true;
}

bool OutputManager::Writesca(int itheta, int ibeta, int iphi, bool wantIobin, int iori, int irad, int iwav, int navg, int *itnum, int ncomp, int nori, 
	char *cdescr, real aeff, real ak1, Vect3<real> &ak_tf, real wave, real xx, SumPackage &sumPackage, real *s1111, real *s2121, real *sm, real *smori, 
	Complex *cx1121, Vect3<Complex> &cxe01_tf, Vect3<Complex> &cxe02_tf, FourArray *cxfData, real pyddx, real pzddx, Vect6<real> &xMinmax)
{
/* **
Purpose of WRITESCA is to collect a lot of the output code into a single module.

Given:
    JPBC   = 0 for isolated finite target
           = 1 for target periodic in y direction in Target Frame
           = 2 for target periodic in z direction in Target Frame
           = 3 for target periodic in y and z directions in TF

If IORI = 0, then print orientational average
If IORI > 0, print output for specific orientation

History:
Fortran versions history removed.

Copyright (C) 1996,1998,2003,2004,2005,2006,2007,2008
              B.T. Draine and P.J. Flatau
This code is covered by the GNU General Public License.
** */

	const real half_ = (real)0.5;
	const real zero_ = (real)0.;
	const real threeFourth = (real)0.75;

	const char *Format9020 = " DDSCAT --- %s\n TARGET ---%s\n %s --- DDA method\n %s --- CCG method\n %s --- shape \n%8d     = NAT0 = number of dipoles\n%12.8lf = d/aeff for this target [d=dipole spacing]\n%12.6lf = d (physical units)";
	const char *Format9031 = "n= (%7.4lf , %7.4lf),  eps.= (%8.4lf , %7.4lf)  |m|kd=%8.4lf for subs.%2d\n";
	const char *Format9032 = "\n  AEFF=%14.6lf = effective radius (physical units)\n  WAVE=%14.6lf = wavelength (in vacuo, physical units)\nK*AEFF=%14.6lf = 2*pi*aeff/lambda\nNABBIENT=%12.6lf = refractive index of ambient medium\n";
	const char *Format9034 = "  NAVG=%6d = (theta,phi) values used in comp. of Qsca,g\n";
	const char *Format9037 = "(%8.5lf%9.5lf%9.5lf ) = k vector (latt. units) in Lab Frame\n(%8.5lf,%8.5lf)(%8.5lf,%8.5lf)(%8.5lf,%8.5lf )=inc.pol.vec. 1 in LF\n(%8.5lf,%8.5lf)(%8.5lf,%8.5lf)(%8.5lf,%8.5lf )=inc.pol.vec. 2 in LF\n";
	const char *Format9041 = "\n%7.4lf = ETASCA = param. controlling # of scatt. dirs used to calculate <cos> etc.\n";
	const char *Format9050 = "          Qext         Qabs         Qsca       g(1)=<cos>  <cos^2>     Qbk       Qpha\n JO=1: %12.4e%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e\n";
	const char *Format9150 = "         Qsca*g(1)   Qsca*g(2)   Qsca*g(3)   iter  mxiter  Nsca\n JO=1: %12.4e%12.4e%12.4e%7d%7d%7d\n";
	const char *Format9250 = "        Qtrqab[0]  Qtrqab[1]  Qtrqab[2]  Qtrqsc[0]  Qtrqsc[1]  Qtrqsc[2]\n JO=1: %12.4e%12.4e%12.4e%12.4e%12.4e%12.4e\n";
	const char *Format9060 = "**** Selected scattering directions  [note: incident pol state 1 is FIXED!]\n ND THETA   PHI <|f11|^2> <|f21|^2> Re<f11*f21> Im<f11*f21>\n";
	const char *Format9070 = "%3d%6.1lf%6.1lf%10.3lf%10.3lf%11.3lf%11.3lf";
	const char *Format9055 = " JO=2: %12.4e%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e\n mean: %12.4e%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e\n Qpol= %12.4e                                                  dQpha=%12.4e\n";
	const char *Format9155 = " JO=2: %12.4e%12.4e%12.4e%7d%7d%7d\n mean: %12.4e%12.4e%12.4e\n";
	const char *Format9056 = "%10.4e%11.4e%11.4e%11.4e%11.4e%12.4e%11.4e%11.4e%6d\n";
	const char *Format9057 = "%10.4e%11.4e%12.4e%12.4e%12.4e\n";
	const char *Format9080 = "            Mueller matrix elements for selected scattering directions in %s\n";
	const char *Format9090_ = "%6.2lf%7.2lf%9.5lf";
	const char *Format9090  = "%12.4e";

	const char *Format90xx = " theta    phi    Pol.    ";
	const char *Format91xx = " alpha   zeta    Pol.    ";
	const char *Format91yy = "S_%1d%1d        ";
	const char *Format9255 = " Jd=2:%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e\n mean:%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e\n";

	const char *cstamp = Version();
//
// Compute DAEFF=d/aeff and DPHYS=d (physical units)
	real daeff = Pow(FourPi / ((real)3. * currentTarget->Nat0()), (real)(1./3.));
	real dphys = wave * ak1 / TwoPi;
	real qav[18];
//
// If IORI > 0, then write out scattering properties for this specific target orientation.
	std::string cframe = FrameEnumeratorVerbose(param->Cmdfrm());
//
	Vect3<real> a1_lf, a2_lf;
	real g[2], g2[2];
	int i, j;
	if (iori >= 0)								// Compute G=<cos(theta)> and G2=<cos^2(theta)>
	{
		for(j=0; j<param->Iorth(); ++j)
		{
			if (jpbc == PeriodicNo)
			{
				g[j]  = sumPackage.Qscag().Q()[j].data[0] / sumPackage.Qsca().Q()[j];
				g2[j] = sumPackage.Qscag2().Q()[j] / sumPackage.Qsca().Q()[j];
			}
			else								// QSCAT and QSCAG are not available when JPBC > 0
			{
				g[j] = g2[j] = zero_;
			}
		}
//
		real absco[2], rvar1, rvar2, rvar3;
		if (param->Iwrksc() == 1)
		{
			const char *Format9030 = "\n----- physical extent of target volume in Target Frame ------\n%14.6lf%14.6lf = xmin,xmax (physical units)\n%14.6lf%14.6lf = ymin,ymax (physical units)\n%14.6lf%14.6lf = zmin,zmax (physical units)";
			const char *Format9035 = "(%8.5lf%9.5lf%9.5lf ) = k vector (latt. units) in TF\n(%8.5lf,%8.5lf)(%8.5lf,%8.5lf)(%8.5lf,%8.5lf )=inc.pol.vec. 1 in TF\n(%8.5lf,%8.5lf)(%8.5lf,%8.5lf)(%8.5lf,%8.5lf )=inc.pol.vec. 2 in TF";

			FileNamer::GetInstance()->Namer(iwav, irad, iori);
			file8 = fopen(FileNamer::GetInstance()->GetCflsca(), "w");
			if (file8 == NULL)
			{
				Wrimsg("Writesca", "Cannot open file8, line 116");
				return false;
			}
			fprintf(file8, Format9020, cstamp, cdescr, AlphaEnumerator(param->Calpha()), SolEnumerator(param->Cmdsol()), param->TargetString(), currentTarget->Nat0(), daeff, dphys);
			fprintf(file8, Format9030, xMinmax.data[0]*dphys, xMinmax.data[1]*dphys, xMinmax.data[2]*dphys, xMinmax.data[3]*dphys, xMinmax.data[4]*dphys, xMinmax.data[5]*dphys);
//
// Following code to be enabled for noncubic treatment:
			fprintf(file8, Format9032, aeff, wave, xx, param->Nambient());
			for(j=0; j<ncomp; ++j)
			{
				real mkd = daeff * dielec->GetCxrfr(j).mod() * xx;
				fprintf(file8, Format9031, dielec->GetCxrfr(j).re, dielec->GetCxrfr(j).im, dielec->GetCxeps(j).re, dielec->GetCxeps(j).im, mkd, j);
			}
			fprintf(file8, "   TOL=%10.3e = error tolerance for CCG method\n", param->Tol());
			currentTarget->A1().Fprintf(file8, "%8.5lf", "(", ") = target axis A1 in Target Frame\n");
			currentTarget->A2().Fprintf(file8, "%8.5lf", "(", ") = target axis A2 in Target Frame\n");
			if (jpbc == PeriodicNo)
				fprintf(file8, Format9034, navg);
			fprintf(file8, Format9035, ak_tf.data[0], ak_tf.data[1], ak_tf.data[2], 
				cxe01_tf.data[0].re, cxe01_tf.data[0].im, cxe01_tf.data[1].re, cxe01_tf.data[1].im, cxe01_tf.data[2].re, cxe01_tf.data[2].im, 
				cxe02_tf.data[0].re, cxe02_tf.data[0].im, cxe02_tf.data[1].re, cxe02_tf.data[1].im, cxe02_tf.data[2].re, cxe02_tf.data[2].im);
// 
// Calculate A1_LF = target axis A1 in Lab Frame A2_LF = target axis A2 in Lab Frame
			const char *Format9036 = "\n(%8.5lf%9.5lf%9.5lf ) = target axis A1 in Lab Frame\n(%8.5lf%9.5lf%9.5lf ) = target axis A2 in Lab Frame\n";
			const char *Format9040 = " BETA =%7.3lf = rotation of target around A1\n THETA=%7.3lf = angle between A1 and k\n  PHI =%7.3lf = rotation of A1 around k\n";
			real bb = param->GetOriData()->Beta(ibeta);
			real tt = param->GetOriData()->Theta(itheta);
			real pp = param->GetOriData()->Phi(iphi);
			a1_lf.Set(	Cos(tt), Sin(tt) * Cos(pp), Sin(tt) * Sin(pp));
			a2_lf.Set( -Sin(tt) * Cos(bb), Cos(tt) * Cos(bb) * Cos(pp) - Sin(bb) * Sin(pp), Cos(tt) * Cos(bb) * Sin(pp) + Sin(bb) * Cos(pp));
			fprintf(file8, Format9036, a1_lf.data[0], a1_lf.data[1], a1_lf.data[2], a2_lf.data[0], a2_lf.data[1], a2_lf.data[2]);
			rvar1 = ak1;
			rvar2 = zero_;
			rvar3 = zero_;
			fprintf(file8, Format9037, rvar1, rvar2, rvar3, 
				param->Cxe01_lf().data[0].re, param->Cxe01_lf().data[0].im, param->Cxe01_lf().data[1].re, param->Cxe01_lf().data[1].im, 
				param->Cxe01_lf().data[2].re, param->Cxe01_lf().data[2].im, param->Cxe02_lf().data[0].re, param->Cxe02_lf().data[0].im, 
				param->Cxe02_lf().data[1].re, param->Cxe02_lf().data[1].im, param->Cxe02_lf().data[2].re, param->Cxe02_lf().data[2].im);
			fprintf(file8, Format9040, betad, thetad, phid);

			switch(jpbc)
			{
			case PeriodicNo:
				fprintf(file8, Format9041, param->Etasca());
				fprintf(file8, Format9050, sumPackage.Qext().Q()[0], sumPackage.Qabs().Q()[0], sumPackage.Qsca().Q()[0], g[0], g2[0], sumPackage.Qbksca().Q()[0], sumPackage.Qpha().Q()[0]);
				break; 
			
			case PeriodicY:
			case PeriodicZ:
				fprintf(file8, "*** need to figure out what we want to print for total cross section/unit length for JPBC=1,2 ***\n");
				break;

			case PeriodicBoth:					// Abs.coeff = 1-R-T = [Qabs*pi/cos(theta_i)]*(d/L_y)*(d/L_z)*(3*N/4*pi)^{2/3}
				absco[0] = Pi * (ak1 / Fabs(ak_tf.data[0])) * sumPackage.Qabs().Q()[0] * Pow((threeFourth * (real)currentTarget->Nat0() / Pi), (real)(2./3.)) / (pyddx * pzddx);
				fprintf(file8, " absorption coeff.  iter  mxiter\n JO=1: %12.4e %6d %6d\n", absco[0], itnum[0], param->Mxiter());
				break;

			default:
				break;
			}
		}

		const real nnn = (real)999.;
		if (param->Iorth() == 1)
		{
//
// Assign "missing value" value if there is only one polarization
			for(i=0; i<18; ++i)
				qav[i] = -nnn;
			sumPackage.Qext().Q()[1] = -nnn;
			sumPackage.Qabs().Q()[1] = -nnn;
			sumPackage.Qsca().Q()[1] = -nnn;
			sumPackage.Qbksca().Q()[1] = -nnn;
			sumPackage.Qpha().Q()[1] = -nnn;
			g[1] = -nnn;
			g2[1] = -nnn;
			sumPackage.Qscag().Q()[1].Set(-nnn, -nnn, -nnn);
			sumPackage.Qscag2().Q()[1] = -nnn;
		}

		if (param->Iorth() == 2)
		{
			if (jpbc == PeriodicNo)
			{
				qav[0] = half_ * sumPackage.Qext().GetSumQ();
				qav[1] = half_ * sumPackage.Qabs().GetSumQ();
				qav[2] = half_ * sumPackage.Qsca().GetSumQ();
				qav[3] = half_ * (sumPackage.Qscag().Q()[0].data[0] + sumPackage.Qscag().Q()[1].data[0]) / qav[2];
				qav[4] = half_ * sumPackage.Qscag2().GetSumQ() / qav[2];
				qav[5] = half_ * sumPackage.Qbksca().GetSumQ();
				qav[6] = half_ * sumPackage.Qpha().GetSumQ();
				qav[7] = sumPackage.Qext().Q()[0] - sumPackage.Qext().Q()[1];
				qav[8] = sumPackage.Qpha().GetDifQ();
				for(j=0; j<3; ++j)
				{
					qav[9+j] = half_ * (sumPackage.Qscag().Q()[0].data[j] + sumPackage.Qscag().Q()[1].data[j]);
				}
				if (param->Cmdtrq() == TorqMethod_DOTORQ)
				{
					for(j=0; j<3; ++j)
					{
						qav[12+j] = half_ * (sumPackage.Qtrqab().Q()[0].data[j] + sumPackage.Qtrqab().Q()[1].data[j]);
						qav[15+j] = half_ * (sumPackage.Qtrqsc().Q()[0].data[j] + sumPackage.Qtrqsc().Q()[1].data[j]);
					}
				}
			}
			else
			{
				qav[0] = half_ * sumPackage.Qext().GetSumQ();
				qav[1] = half_ * sumPackage.Qabs().GetSumQ();
				for(j=2; j<6; ++j)
				{
					qav[j] = zero_;
				}
				qav[6] = half_ * sumPackage.Qpha().GetSumQ();
				qav[7] = sumPackage.Qext().GetDifQ();
				qav[8] = sumPackage.Qpha().GetDifQ();
				for(j=9; j<18; ++j)
				{
					qav[j] = zero_;
				}
			}
			if (param->Iwrksc() == 1)
			{
				switch(jpbc)
				{
				case PeriodicNo:
					fprintf(file8, Format9055, sumPackage.Qext().Q()[1], sumPackage.Qabs().Q()[1], sumPackage.Qsca().Q()[1], g[1], g2[1], sumPackage.Qbksca().Q()[1], sumPackage.Qpha().Q()[1], 
						qav[0], qav[1], qav[2], qav[3], qav[4], qav[5], qav[6], qav[7], qav[8]);
					break;

				case PeriodicY:
				case PeriodicZ:
					fprintf(file8, "*** need to figure out what we want to print for total cross section/unit length for JPBC=1,2 ***\n");
					break;

				case PeriodicBoth:
					absco[1] = Pi * (ak1 / Fabs(ak_tf.data[0])) * sumPackage.Qabs().Q()[1] * Pow((threeFourth * (real)currentTarget->Nat0() / Pi), (real)(2./3.)) / (pyddx * pzddx);
					qav[1] = half_*(absco[0] + absco[1]);
					fprintf(file8, " JO=2: %12.4e%7d%7d\n mean: %12.4e\n", absco[1], itnum[1], param->Mxiter(), qav[1]);
					break;

				default:
					break;
				}
			}
		}

		if (param->Iwrksc() == 1)
		{
			if (jpbc == PeriodicNo)
			{
				fprintf(file8, Format9150, sumPackage.Qscag().Q()[0].data[0], sumPackage.Qscag().Q()[0].data[1], sumPackage.Qscag().Q()[0].data[2], itnum[0], param->Mxiter(), navg);
				if (param->Iorth() == 2)
					fprintf(file8, Format9155, sumPackage.Qscag().Q()[1].data[0], sumPackage.Qscag().Q()[1].data[1], sumPackage.Qscag().Q()[1].data[2], itnum[1], param->Mxiter(), navg, qav[9], qav[10], qav[11]);
			}
		}

		if (param->Cmdtrq() == TorqMethod_DOTORQ)
		{
			if (param->Iwrksc() == 1)
			{
				fprintf(file8, Format9250, sumPackage.Qtrqab().Q()[0].data[0], sumPackage.Qtrqab().Q()[0].data[1], sumPackage.Qtrqab().Q()[0].data[2], 
					sumPackage.Qtrqsc().Q()[0].data[0], sumPackage.Qtrqsc().Q()[0].data[1], sumPackage.Qtrqsc().Q()[0].data[2]);
				if (param->Iorth() == 2)
					fprintf(file8, Format9255, sumPackage.Qtrqab().Q()[1].data[0], sumPackage.Qtrqab().Q()[1].data[1], sumPackage.Qtrqab().Q()[1].data[2], 
						sumPackage.Qtrqsc().Q()[1].data[0], sumPackage.Qtrqsc().Q()[1].data[1], sumPackage.Qtrqsc().Q()[1].data[2], qav[12], qav[13], qav[14], qav[15], qav[16], qav[17]); 
			}
		}
		else									// Add missing value code (set to -999.)
		{
			for(j=0; j<2; ++j)
			{
				sumPackage.Qtrqab().Q()[j].Set(-nnn, -nnn, -nnn);
				sumPackage.Qtrqsc().Q()[j].Set(-nnn, -nnn, -nnn);
			}
		}

		int nd;
		if (param->Iorth() == 1)
		{
			if (param->Iwrksc() == 1)
				fprintf(file8, Format9060);
			for(nd=0; nd<param->Nscat(); ++nd)					// Convert scattering angles to degrees
			{
				real phind = Degrad * param->Phin()[nd];
				real thetnd = Degrad * param->Thetan()[nd];
				rvar1 = (cxfData->Cx11(nd).conjg() * cxfData->Cx11(nd)).re;
				rvar2 = (cxfData->Cx21(nd).conjg() * cxfData->Cx21(nd)).re;
				Complex cxvar1 = cxfData->Cx11(nd).conjg() * cxfData->Cx21(nd);
				if (param->Iwrksc() == 1) 
					fprintf(file8, Format9070, nd, thetnd, phind, rvar1, rvar2, cxvar1.re, cxvar1.im); 
			}
		}
		else
		{
			if (param->Iwrksc() == 1)
			{
				int nsmelts = param->Nsmelts();
				fprintf(file8, Format9080, cframe.c_str());												// Write file header information
				if (nsmelts > 0)
				{
					if ((jpbc == PeriodicNo) || (jpbc == PeriodicBoth))
						fprintf(file8, Format90xx);
					else
						fprintf(file8, Format91xx);
					for(j=0; j<nsmelts; ++j)
						fprintf(file8, Format91yy, param->Smind1()[j], param->Smind2()[j]);
					fprintf(file8, "\n");
				}
				for(nd=0; nd<param->Nscat(); ++nd)
				{
					int index = 16 * nd;
// Convert scattering angles to degrees
					real phind = Degrad * param->Phin()[nd];
					real thetnd = Degrad * param->Thetan()[nd];
// Compute PPOL = degree of linear polarization of scattered light for incident unpolarized light
					real ppol = Sqrt(sm[index + 4] * sm[index + 4] + sm[index + 8] * sm[index + 8]) / sm[index];
					fprintf(file8, Format9090_, thetnd, phind, ppol);
					for(j=0; j<nsmelts; ++j)
						fprintf(file8, Format9090, sm[index + 4*(param->Smind1()[j]-1) + (param->Smind2()[j]-1)]);
					fprintf(file8, "\n");
				}
			}
		}
		if (param->Iwrksc() == 1)
			fclose(file8);
	}
	else
	{
		bool bWhatSum = sumPackage.UseWhatSum();
		const char *Format9042 = "%8.3lf%8.3lf = beta_min, beta_max ;  NBETA =%2d\n%8.3lf%8.3lf = theta_min, theta_max; NTHETA=%2d\n%8.3lf%8.3lf = phi_min, phi_max   ;   NPHI =%2d\n";
		const char *Format9043 = " Results averaged over %4d target orientations\n                   and %4d incident polarizations\n";
//
// This module writes out results from orientational averaging. (arrive here when called with IORI=0)
		file8 = fopen(FileNamer::GetInstance()->GetCflavg(), "w");
		if (file8 == NULL)
		{
			Wrimsg("Writesca", "Cannot open file8, line 372");
		}
		fprintf(file8, Format9020, cstamp, cdescr, AlphaEnumerator(param->Calpha()), SolEnumerator(param->Cmdsol()), param->TargetString(), currentTarget->Nat0(), daeff, dphys);
//
// Following code to be enabled for noncubic treatment: 
		//WRITE(8,FMT=9020)CSTAMP,CDESCR,CMDSOL,CALPHA,CSHAPE,NAT0,DX
		fprintf(file8, Format9032, aeff, wave, xx, param->Nambient());
		for(j=0; j<ncomp; ++j)
		{
			real mkd = Pow((real)(Pi / (threeFourth * currentTarget->Nat0())), (real)(1./3.)) * dielec->GetCxrfr(j).mod() * xx;	
			fprintf(file8, Format9031, dielec->GetCxrfr(j).re, dielec->GetCxrfr(j).im, dielec->GetCxeps(j).re, dielec->GetCxrfr(j).im, mkd, j);
		}
		fprintf(file8, "   TOL=%10.3e = error tolerance for CCG method\n", param->Tol());
		currentTarget->A1().Fprintf(file8, "%8.5lf", "(", ") = target axis A1 in Target Frame\n");
		currentTarget->A2().Fprintf(file8, "%8.5lf", "(", ") = target axis A2 in Target Frame\n");
		if (jpbc == PeriodicNo)
			fprintf(file8, Format9034, navg);
		real rvar1 = ak1;
		real rvar2 = zero_;
		real rvar3 = zero_;
		fprintf(file8, Format9037, rvar1, rvar2, rvar3, 
			param->Cxe01_lf().data[0].re, param->Cxe01_lf().data[0].im, param->Cxe01_lf().data[1].re, param->Cxe01_lf().data[1].im, 
			param->Cxe01_lf().data[2].re, param->Cxe01_lf().data[2].im, param->Cxe02_lf().data[0].re, param->Cxe02_lf().data[0].im, 
			param->Cxe02_lf().data[1].re, param->Cxe02_lf().data[1].im, param->Cxe02_lf().data[2].re, param->Cxe02_lf().data[2].im);
		fprintf(file8, Format9042, betmid, betmxd, param->GetOriData()->GetNbeta(), thtmid, thtmxd, param->GetOriData()->GetNtheta(), phimid, phimxd, param->GetOriData()->GetNphi());
		if (jpbc == PeriodicNo)
			fprintf(file8, Format9041, param->Etasca());
		fprintf(file8, Format9043, nori, param->Iorth());
		if (jpbc == PeriodicNo)
		{
			rvar1 = sumPackage.Qscag().Qsum(bWhatSum)[0].data[0] / sumPackage.Qsca().Qsum(bWhatSum)[0];
			rvar2 = sumPackage.Qscag2().Qsum(bWhatSum)[0] / sumPackage.Qsca().Qsum(bWhatSum)[0];
			fprintf(file8, Format9050, sumPackage.Qext().Qsum(bWhatSum)[0], sumPackage.Qabs().Qsum(bWhatSum)[0], sumPackage.Qsca().Qsum(bWhatSum)[0], rvar1, rvar2, 
				sumPackage.Qbksca().Qsum(bWhatSum)[0], sumPackage.Qpha().Qsum(bWhatSum)[0]);
		}
		else
		{
			const char *Format9051 = "          Qext       Qabs       Qsca\n JO=1: %12.4e%12.4e%12.4e\n";
			sumPackage.Qsca().Qsum(bWhatSum)[0] = sumPackage.Qext().Qsum(bWhatSum)[0] - sumPackage.Qabs().Qsum(bWhatSum)[0];
			fprintf(file8, Format9051, sumPackage.Qext().Qsum(bWhatSum)[0], sumPackage.Qabs().Qsum(bWhatSum)[0], sumPackage.Qsca().Qsum(bWhatSum)[0]);
		}

		int nd;
		if (param->Iorth() == 1)
		{
			if (jpbc == PeriodicNo)
			{
				fprintf(file8, Format9150, sumPackage.Qscag().Qsum(bWhatSum)[0].data[0], sumPackage.Qscag().Qsum(bWhatSum)[0].data[1], sumPackage.Qscag().Qsum(bWhatSum)[0].data[2], 
					itnum[0], param->Mxiter(), navg);
				if (param->Cmdtrq() == TorqMethod_DOTORQ)
					fprintf(file8, Format9250, sumPackage.Qtrqab().Qsum(bWhatSum)[0].data[0], sumPackage.Qtrqab().Qsum(bWhatSum)[0].data[1], sumPackage.Qtrqab().Qsum(bWhatSum)[0].data[2], 
						sumPackage.Qtrqsc().Qsum(bWhatSum)[0].data[0], sumPackage.Qtrqsc().Qsum(bWhatSum)[0].data[1], sumPackage.Qtrqsc().Qsum(bWhatSum)[0].data[2]);
			}
			else
			{
				const char *Format9151 = "                                            iter  mxiter\n";
				fprintf(file8, Format9151);
			}
			fprintf(file8, Format9060);
			for(nd=0; nd<param->Nscat(); ++nd)						// ! Convert scattering angles to degrees
			{
				real phind = Degrad * param->Phin()[nd];
				real thetnd = Degrad * param->Thetan()[nd];
                fprintf(file8, Format9070, nd, thetnd, phind, s1111[nd], s2121[nd], cx1121[nd].re, cx1121[nd].im);
			}
// Write orientationally-averaged Q values (except Qpha) to 'qtable':
			file10 = fopen("qtable", "a+");
			if (!file10)
			{
				Wrimsg("Writesca", "Cannot open qtable file 10\n");
				return false;
			}
			fseek(file10, 0L, SEEK_END);
            if (jpbc == PeriodicNo)
			{
				rvar1 = sumPackage.Qscag().Qsum(bWhatSum)[0].data[0] / sumPackage.Qsca().Qsum(bWhatSum)[0];
				rvar2 = sumPackage.Qscag2().Qsum(bWhatSum)[0] / sumPackage.Qsca().Qsum(bWhatSum)[0];
				fprintf(file10, Format9056, aeff, wave, sumPackage.Qext().Qsum(bWhatSum)[0], sumPackage.Qabs().Qsum(bWhatSum)[0], sumPackage.Qsca().Qsum(bWhatSum)[0], rvar1, rvar2, 
					sumPackage.Qscag2().Qsum(bWhatSum)[0], navg);
			}
			else
			{
				fprintf(file10, Format9057, aeff, wave, sumPackage.Qext().Qsum(bWhatSum)[0], sumPackage.Qabs().Qsum(bWhatSum)[0], sumPackage.Qsca().Qsum(bWhatSum)[0]);
			}
            fclose(file10);

			file12 = fopen("qtable2", "a+");
			if (!file12)
			{
				Wrimsg("Writesca", "Cannot open qtable2 file 12\n");
				return false;
			}
			fseek(file12, 0L, SEEK_END);
			fprintf(file12, Format9057, aeff, wave, sumPackage.Qpha().Qsum(bWhatSum)[0]);
			fclose(file12);
		}
		else
		{
			if (jpbc == PeriodicNo)
			{
				rvar1 = sumPackage.Qscag().Qsum(bWhatSum)[1].data[0] / sumPackage.Qsca().Qsum(bWhatSum)[1];
				qav[0] = half_ * (sumPackage.Qext().Qsum(bWhatSum)[0] + sumPackage.Qext().Qsum(bWhatSum)[1]);
				qav[1] = half_ * (sumPackage.Qabs().Qsum(bWhatSum)[0] + sumPackage.Qabs().Qsum(bWhatSum)[1]);
				qav[2] = half_ * (sumPackage.Qsca().Qsum(bWhatSum)[0] + sumPackage.Qsca().Qsum(bWhatSum)[1]);
				qav[3] = half_ * (sumPackage.Qscag().Qsum(bWhatSum)[0].data[0] + sumPackage.Qscag().Qsum(bWhatSum)[1].data[0]) / qav[2];
				qav[4] = half_ * (sumPackage.Qscag2().Qsum(bWhatSum)[0] + sumPackage.Qscag2().Qsum(bWhatSum)[1]) / qav[2];
				qav[5] = half_ * (sumPackage.Qbksca().Qsum(bWhatSum)[0] + sumPackage.Qbksca().Qsum(bWhatSum)[1]);
				qav[6] = half_ * (sumPackage.Qpha().Qsum(bWhatSum)[0] + sumPackage.Qpha().Qsum(bWhatSum)[1]);
				qav[7] = sumPackage.Qext().Qsum(bWhatSum)[0] - sumPackage.Qext().Qsum(bWhatSum)[1];
				qav[8] = sumPackage.Qpha().Qsum(bWhatSum)[0] - sumPackage.Qpha().Qsum(bWhatSum)[1];
				for(j=0; j<3; ++j)
				{
					qav[9+j] = half_ * (sumPackage.Qscag().Qsum(bWhatSum)[0].data[j] + sumPackage.Qscag().Qsum(bWhatSum)[1].data[j]);
				}
				if (param->Cmdtrq() == TorqMethod_DOTORQ)
				{
					for(j=0; j<3; ++j)
					{
						qav[12+j] = half_ * (sumPackage.Qtrqab().Qsum(bWhatSum)[0].data[j] + sumPackage.Qtrqab().Qsum(bWhatSum)[1].data[j]);
						qav[15+j] = half_ * (sumPackage.Qtrqsc().Qsum(bWhatSum)[0].data[j] + sumPackage.Qtrqsc().Qsum(bWhatSum)[1].data[j]);
					}
				}
				rvar2 = half_ * sumPackage.Qscag2().GetSumQ() / qav[2];
			}
			else
			{
				sumPackage.Qsca().Qsum(bWhatSum)[1] = sumPackage.Qext().Qsum(bWhatSum)[1] - sumPackage.Qabs().Qsum(bWhatSum)[1];
				rvar1 = rvar2 = zero_;
				qav[0] = half_ * (sumPackage.Qext().Qsum(bWhatSum)[0] + sumPackage.Qext().Qsum(bWhatSum)[1]);
				qav[1] = half_ * (sumPackage.Qabs().Qsum(bWhatSum)[0] + sumPackage.Qabs().Qsum(bWhatSum)[1]);
				qav[2] = half_ * (sumPackage.Qsca().Qsum(bWhatSum)[0] + sumPackage.Qsca().Qsum(bWhatSum)[1]);
				for(j=3; j<6; ++j)
					qav[j] = zero_;
				qav[6] = half_ * (sumPackage.Qpha().Qsum(bWhatSum)[0] + sumPackage.Qpha().Qsum(bWhatSum)[1]);
				qav[7] = sumPackage.Qext().Qsum(bWhatSum)[0] - sumPackage.Qext().Qsum(bWhatSum)[1];
				qav[8] = sumPackage.Qpha().Qsum(bWhatSum)[0] - sumPackage.Qpha().Qsum(bWhatSum)[1];
				for(j=9; j<12; ++j)
					qav[j] = zero_;
			}
//
			if (jpbc == PeriodicNo)
			{
				fprintf(file8, Format9055, sumPackage.Qext().Qsum(bWhatSum)[1], sumPackage.Qabs().Qsum(bWhatSum)[1], sumPackage.Qsca().Qsum(bWhatSum)[1], rvar1, rvar2, 
					sumPackage.Qbksca().Qsum(bWhatSum)[1], sumPackage.Qpha().Qsum(bWhatSum)[1], qav[0], qav[1], qav[2], qav[3], qav[4], qav[5], qav[6], qav[7], qav[8]);
				fprintf(file8, Format9150, sumPackage.Qscag().Qsum(bWhatSum)[0].data[0], sumPackage.Qscag().Qsum(bWhatSum)[0].data[1], sumPackage.Qscag().Qsum(bWhatSum)[0].data[2], 
					itnum[0], param->Mxiter(), navg);
				fprintf(file8, Format9155, sumPackage.Qscag().Qsum(bWhatSum)[1].data[0], sumPackage.Qscag().Qsum(bWhatSum)[1].data[1], sumPackage.Qscag().Qsum(bWhatSum)[1].data[2], 
					itnum[1], param->Mxiter(), navg, qav[9], qav[10], qav[11]);
				if (param->Cmdtrq() == TorqMethod_DOTORQ)
				{
					fprintf(file8, Format9250, sumPackage.Qtrqab().Qsum(bWhatSum)[0].data[0], sumPackage.Qtrqab().Qsum(bWhatSum)[0].data[1], sumPackage.Qtrqab().Qsum(bWhatSum)[0].data[2], 
						sumPackage.Qtrqsc().Qsum(bWhatSum)[0].data[0], sumPackage.Qtrqsc().Qsum(bWhatSum)[0].data[1], sumPackage.Qtrqsc().Qsum(bWhatSum)[0].data[2]);
					fprintf(file8, Format9255, sumPackage.Qtrqab().Qsum(bWhatSum)[1].data[0], sumPackage.Qtrqab().Qsum(bWhatSum)[1].data[1], sumPackage.Qtrqab().Qsum(bWhatSum)[1].data[2], 
						sumPackage.Qtrqsc().Qsum(bWhatSum)[1].data[0], sumPackage.Qtrqsc().Qsum(bWhatSum)[1].data[1], sumPackage.Qtrqsc().Qsum(bWhatSum)[1].data[2], qav[12], qav[13], qav[14], qav[15], qav[16], qav[17]);
				}
			}
            else
			{
				const char *Format9054 = " JO=2: %12.4e%12.4e%12.4e\n mean: %12.4e%12.4e%12.4e\n Qpol= %12.4e\n";
				fprintf(file8, Format9054, sumPackage.Qext().Qsum(bWhatSum)[1], sumPackage.Qabs().Qsum(bWhatSum)[1], sumPackage.Qsca().Qsum(bWhatSum)[1], qav[0], qav[1], qav[2], qav[7]);
			}
//
// Write orientationally-averaged Q values (except Qpha) to 'qtable':
			file10 = fopen("qtable", "a+");
			if (!file10)
			{
				Wrimsg("Writesca", "Cannot open qtable (a+) file 10\n");
				return false;
			}
			fseek(file10, 0L, SEEK_END);
			fprintf(file10, Format9056, aeff, wave, qav[0], qav[1], qav[2], qav[3], qav[4], qav[5], navg);
			fclose(file10);

			file12 = fopen("qtable2", "a+");
			if (!file12)
			{
				Wrimsg("Writesca", "Cannot open qtable2 (a+) file 12\n");
				return false;
			}
			fseek(file12, 0L, SEEK_END);
			fprintf(file12, Format9057, aeff, wave, qav[6], qav[7], qav[8]);
			fclose(file12);
//	
			fprintf(file8, Format9080, cframe.c_str());
			fprintf(file8, Format90xx);
			int nsmelts = param->Nsmelts();
			for(j=0; j<nsmelts; ++j)
				fprintf(file8, Format91yy, param->Smind1()[j], param->Smind2()[j]);
			fprintf(file8, "\n");
			for(nd=0; nd<param->Nscat(); ++nd)						// Convert scattering angles to degrees
			{
				int index = 16 * nd;
				real phind = Degrad * param->Phin()[nd];
				real thetnd = Degrad * param->Thetan()[nd];
//
// Compute PPOL = degree of linear polarization of scattered light for incident unpolarized light (Bohren & Huffman p. 382)
				real ppol = Sqrt(smori[index + 4] * smori[index + 4] + smori[index + 8] * smori[index + 8]) / smori[index];
				fprintf(file8, Format9090_, thetnd, phind, ppol);
				for(j=0; j<nsmelts; ++j)
					fprintf(file8, Format9090, smori[index + 4*(param->Smind1()[j]-1) + (param->Smind2()[j]-1)]);
				fprintf(file8, "\n");
			}
		}
		fclose(file8);
	}
//
// Write FORTRAN unformatted file
	if (param->Cbinflag() != BinflagMethod_NOTBIN)
	{
		Writebin(wantIobin, cdescr, navg, ncomp, nori, iwav, irad, iori, aeff, wave, xx, ak1, ak_tf, cxe01_tf, cxe02_tf, sumPackage, g, qav, sm);
	}
	return true;
}

void OutputManager::Writepol(real aeff, Vect3<real> &akr, real wave, Vect3<Complex> &cxe0r, Matrix *theMatrix, Complex *cxpol, const char *cflpol)
{
/* **
Given:
   NRWORD   = length (bytes) of real word [int]
   MXNAT    = dimensioning parameter
   NX       = (1/d)*(maximum extent of target in x-direction) (TF)
   NY       = (1/d)*(maximum extent of target in y-direction) (TF)
   NZ       = (1/d)*(maximum extent of target in z-direction) (TF)
   NAT0     = number of occupied dipole sites
   IANISO   = 0 for isotropic materials
              1 for anisotropic materials with optical axes aligned
                with Target Frame (i.e., DF=TF)
              2 for general anisotropic material, with arbitrary
                orientation of DF at each site relative to TF
   ICOMP(IA,JD)=composition of site # IA for principal axis JD=1,2,3 in "dielectric frame" where dielectric tensor is diagonalized
   IXYZ0(IA,JD)=x/d,y/d,z/d for physical dipole # IA, JD=1-3 where (x,y,z) are in Target Frame
   JPBC
   PYD      = periodicity/d of target replication in y_TF direction
            = 0 if no replication (single target)
   PZD      = periodicity/d of target replication in z_TF direction
            = 0 if no replication (single target)
   GAMMA
   NAMBIENT = (real) refractive index of ambient medium
   A1
   A2
   AKR(1-3) = k_x*d,k_y*d,k_z*d where (k_x,k_y,k_z)=k in ambient medium in Target Frame
   DX(1-3)  = lattice spacing/d in x,y,z directions, with DX(1)*DX(2)*DX(3)=1 normally have DX(1)=DX(2)=DX(3)=1
   X0(1-3)  = location/(DX*d) in TF of lattice site IX=0,IY=0,IZ=0
   WAVE     = wavelength in ambient medium (physical units)
   CXE0R(1-3)=E_x,E_y,E_z for incident wave in ambient medium in Target Frame

!??? BTD 110819 I don't think this is correct
!    CXADIA(3i-2) =alphainv_xx at location i in TF
!          (3i-1) =alphainv_yy at location i in TF
!          (3i)   =alphainv_zz at location i in TF
!          where alphainv = inverse of (Complex polarizability tensor)/d
!    CXAOFF(3i-2) =alpha_yz=alphainv_zy at location i in TF
!          (3i-1) =alpha_zx=alphainv_xz at location i in TF
!          (3i)   =alpha_xy=alphainv_yx at location i in TF
!          where alphainv = inverse of (Complex polarizability matrix/d^
!    CXPOL(i)        =P_x at location i in TF
!         (NAT0+i)   =P_y at location i
!         (2*NAT0+i) =P_z at location i
!    CFLPOL = filename

writes NAT0,IXYZ,CXPOL to unformatted file CFLPOL

Note: IXYZ0, CXPOL, CXADIA, CXAOFF are "reduced" arrays
      IXYZ0 is produced by subroutine EXTEND
      CXPOL is reduced by subroutine REDUCE called from GETFML

      CXADIA is reduced by REDUCE called from GETFLM
      CXAOFF is reduced by REDUCE called from GETFLM

      where "reduced" arrays are limited to values for occupied sites

Note: if IANISO=0 or 1: do not store BETADF,THETADF,PHIDF
      if IANISO=2       do store BETADF,THETADF,PHIDF

B.T. Draine, Princeton Univ. Observatory, 2006.04.08

History
Fortran versions history removed.

Copyright (C) 2006,2007,2008,2011, B.T. Draine and P.J. Flatau
This code is covered by the GNU General Public License.
** */
	int j;

	int nrword = sizeof(real);
	int na = currentTarget->Nat0();
	Vect6<int> ixyzMinmax = currentTarget->IxyzMinmax();
#ifdef _WIN32
	int modex = O_WRONLY|O_BINARY|O_CREAT;
#else
	int modex = O_WRONLY|O_CREAT;
#endif	
	int file32 = open(cflpol, modex, S_IWRITE);
	if (file32 == -1)
	{
		fprintf(stderr, "Error opening pol file in Writepol.\n");
		return;
	}
	write(file32, &nrword, sizeof(int));
	write(file32, &na, sizeof(int));
	write(file32, &currentTarget->Ianiso(), sizeof(int));
	ixyzMinmax.Write(file32);
	write(file32, &currentTarget->Nx(), sizeof(int));
	write(file32, &currentTarget->Ny(), sizeof(int));
	write(file32, &currentTarget->Nz(), sizeof(int));
//
	for(j=0; j<na; ++j)
	{
		currentTarget->Ixyz().WriteLine(file32, j);
	}
	j = (int)jpbc;
	write(file32, &j, sizeof(int));
	write(file32, &currentTarget->Pyd(), sizeof(real));
	write(file32, &currentTarget->Pzd(), sizeof(real));
	write(file32, &param->Gamma(), sizeof(real));
	currentTarget->A1().Write(file32);
	currentTarget->A2().Write(file32);
	currentTarget->Dx().Write(file32);
	currentTarget->X0().Write(file32);
//
// Reduce composition vector to only 3*NAT0 elements
	int nxyz = currentTarget->Nx() * currentTarget->Ny() * currentTarget->Nz();
	for(j=0; j<nxyz; ++j)
	{
		if(currentTarget->Icomp().Value(j, 0) > 0)
		{
			currentTarget->Icomp().WriteLine(file32, j);
		}
	}
//
// Write refractive index of ambient medium
	write(file32, &param->Nambient(), sizeof(real));
//
// Polarization vector CXPOL has already been reduced
	write(file32, &wave, sizeof(real));
	write(file32, &aeff, sizeof(real));
	akr.Write(file32);
	cxe0r.Write(file32);
	for(j=0; j<3*na; ++j)
	{
		cxpol[j].Write(file32);
	}
//
// Store elements of A matrix for future use
	theMatrix->WriteDiagonal(file32);
	theMatrix->WriteOffDiagonal(file32);
	if (currentTarget->Ianiso() == TargetIsAnisotropicDisoriented)
	{
		AbstractDFData *data = ((LoadableTarget *)currentTarget)->GetDFdata();
		real t[3];
		for(j=0; j<na; ++j)
		{
			t[0] = data->GetBetadf(j);
			t[1] = data->GetThetadf(j);
			t[2] = data->GetPhidf(j);
			write(file32, t, 3 * sizeof(real));
		}
	}
	close(file32);
}

bool OutputManager::Writebin(bool wantIobin, char *cdescr, int navg, int ncomp, int nori, int iwav, int irad, int iori, real aeff, real wave, real xx, real ak1,  
	Vect3<real>&akr, Vect3<Complex>&cxe01r, Vect3<Complex>&cxe02r, SumPackage &sumPackage, real *g, real *qav, real *sm)
{
/* **
Writes FORTRAN unformatted binary file to file cbinfile ('dd.bin')

Original version created by P.J.Flatau, University of California, SD
History:
Fortran versions history removed.

Copyright (C) 1996,1998,1999,2003,2007 B.T. Draine and P.J. Flatau
This code is covered by the GNU General Public License.
** */
	const char *cbinfile = "dd.bin";
	static int ientry = 0;

	int i;
	if (ientry == 0)
	{
#ifdef _WIN32
		int modex = O_WRONLY|O_BINARY|O_CREAT;
#else
		int modex = O_WRONLY|O_CREAT;
#endif		
		int iobinFile = open(cbinfile, modex, S_IWRITE);
		if (iobinFile == -1)
		{
			Wrimsg("Writebin", "Cannot open binary file during first call to binary write");
			return false;
		}
		else
		{
			Wrimsg("Writebin", "First call to binary write Ok");
		}
		ientry = 1;
// ! General run information:
		const char *cstamp = Version();
		write(iobinFile, cstamp, (int)(strlen(cstamp)+1));
		write(iobinFile, cdescr, 67);
		write(iobinFile, FftEnumerator(param->Cmdfft()), 6);
		write(iobinFile, AlphaEnumerator(param->Calpha()), 6);
		write(iobinFile, param->TargetString(), 9);
		write(iobinFile, TorqEnumerator(param->Cmdtrq()), 6);
// ! mostly dimensions:
		write(iobinFile, &param->Iorth(), sizeof(int));
		write(iobinFile, &currentTarget->Nx(), sizeof(int));
		write(iobinFile, &currentTarget->Ny(), sizeof(int));
		write(iobinFile, &currentTarget->Nz(), sizeof(int));
		write(iobinFile, &currentTarget->Nat0(), sizeof(int));
		write(iobinFile, &currentTarget->Nat(), sizeof(int));		// ChB: binary file format changed!!!
		write(iobinFile, &param->Nscat(), sizeof(int));
		write(iobinFile, &ncomp, sizeof(int));
		write(iobinFile, &nori, sizeof(int));
		write(iobinFile, &param->Nwav(), sizeof(int));
		write(iobinFile, &param->Nrad(), sizeof(int));
		int ntt = param->GetOriData()->GetNtheta();
		int npp = param->GetOriData()->GetNphi();
		int nbb = param->GetOriData()->GetNbeta();
		write(iobinFile, &ntt, sizeof(int));
		write(iobinFile, &npp, sizeof(int));
		write(iobinFile, &nbb, sizeof(int));
		i = TimerManager::GetInstance()->GetTimerNumber();
		write(iobinFile, &i, sizeof(int));
// ! maximum dimension (not really needed)
//		write(iobin, &mxnx, sizeof(int));					// ChB: binary file format changed!!!
//		write(iobin, &mxny, sizeof(int));					// ChB: binary file format changed!!!
//		write(iobin, &mxnz, sizeof(int));					// ChB: binary file format changed!!!
//		write(iobin, &mxcomp, sizeof(int));					// ChB: binary file format changed!!!
//		write(iobin, &mxsca, sizeof(int));					// ChB: binary file format changed!!!
		write(iobinFile, &betmid, sizeof(real));
		write(iobinFile, &betmxd, sizeof(real));
		write(iobinFile, &thtmid, sizeof(real));
		write(iobinFile, &thtmxd, sizeof(real));
		write(iobinFile, &phimid, sizeof(real));
		write(iobinFile, &phimxd, sizeof(real));
// ! dump scattering angles
		write(iobinFile, param->Phin(), param->Nscat() * sizeof(real));
		write(iobinFile, param->Thetan(), param->Nscat() * sizeof(real));
// ! misc run information
		write(iobinFile, &param->Tol(), sizeof(real));
		write(iobinFile, &navg, sizeof(int));
		write(iobinFile, param->Wavea(), param->Nwav() * sizeof(real));
		write(iobinFile, param->Aeffa(), param->Nrad() * sizeof(real));
		write(iobinFile, param->GetOriData()->Theta(), param->GetOriData()->GetNtheta() * sizeof(real));
		write(iobinFile, param->GetOriData()->Beta(), param->GetOriData()->GetNbeta() * sizeof(real));
		write(iobinFile, param->GetOriData()->Phi(), param->GetOriData()->GetNphi() * sizeof(real));
		currentTarget->A1().Write(iobinFile);
		currentTarget->A2().Write(iobinFile);
		for(i=0; i<currentTarget->Nat(); ++i)
		{
			currentTarget->Icomp().WriteLine(iobinFile, i);
		}
		for(i=0; i<currentTarget->Nat(); ++i)
		{
			currentTarget->Ixyz().WriteLine(iobinFile, i);
		}
		currentTarget->Dx().Write(iobinFile);
		close(iobinFile);
	}
//
// Binary write of  for particular "iwav, irad, itheta, ibeta iphi"
// Let us try to keep only  one version of the binary file. Thus, in
// case of IORTH=1 or 'DOTORQ' some values are "missing"
// flatau -->draine  do we need smori for "ppol" do we need cxf_ij ?

	if ((iori != -1) && (DDscatParameters::GetInstance()->Cbinflag() == BinflagMethod_ALLBIN))
	{
		int iu;
#ifdef _WIN32
		int modex = O_WRONLY|O_BINARY|O_CREAT;
#else
		int modex = O_WRONLY|O_CREAT;
#endif				
		int iobinFile = open(cbinfile, modex);
		if (iobinFile == -1)
		{
			Wrimsg("Writebin", "Cannot open binary file during second or next call to binary write");
			return false;
		}
		write(iobinFile, &iwav, sizeof(int));
		write(iobinFile, &irad, sizeof(int));
		write(iobinFile, &iori, sizeof(int));
		for(iu=0; iu<TimerManager::GetInstance()->GetTimerNumber(); ++iu)
		{
			real tt = TimerManager::GetInstance()->GetValue(iu);
			write(iobinFile, &tt, sizeof(real));
		}
        write(iobinFile, &aeff, sizeof(real));
        write(iobinFile, &wave, sizeof(real));
        write(iobinFile, &xx, sizeof(real));
        write(iobinFile, &ak1, sizeof(real));
        write(iobinFile, &betad, sizeof(real));
        write(iobinFile, &thetad, sizeof(real));
        write(iobinFile, &phid, sizeof(real));
        akr.Write(iobinFile);
        cxe01r.Write(iobinFile);
        cxe02r.Write(iobinFile);
		write(iobinFile, dielec->GetCxrfr(), ncomp * sizeof(Complex));
		write(iobinFile, dielec->GetCxeps(), ncomp * sizeof(Complex));
        write(iobinFile, sumPackage.Qext().Q(), 2 * sizeof(real));
        write(iobinFile, sumPackage.Qabs().Q(), 2 * sizeof(real));
        write(iobinFile, sumPackage.Qsca().Q(), 2 * sizeof(real));
        write(iobinFile, g, 2*sizeof(real));
        write(iobinFile, sumPackage.Qbksca().Q(), 2 * sizeof(real));
        write(iobinFile, sumPackage.Qpha().Q(), 2 * sizeof(real));
        write(iobinFile, qav, 17*sizeof(real));
		sumPackage.Qtrqab().Q()[0].Write(iobinFile);
		sumPackage.Qtrqab().Q()[1].Write(iobinFile);
		sumPackage.Qtrqsc().Q()[0].Write(iobinFile);
		sumPackage.Qtrqsc().Q()[1].Write(iobinFile);
		write(iobinFile, sm, 16 * param->Nscat() * sizeof(real));				// ! write sm as 16 independent matrices (IDL handles 7 subscripts only)
		close(iobinFile);
// ! endif "if(iori. ne. 0 .and. cbinflag.eq.'ALLCDF') then"
	}
	return true;
}

bool OutputManager::WriteHeadings(int nori)
{
	const char *cstamp = Version();
	file10 = fopen("qtable", "w+");
	if (file10 == NULL)
	{
		Wrimsg("DDscat", "Cannot open file10, line 1034");
		return false;
	}
	file11 = fopen("mtable", "w+");
	if (file11 == NULL)
	{
		Wrimsg("DDscat", "Cannot open file11, line 1039");
		return false;
	}
	file12 = fopen("qtable2", "w+");
	if (file12 == NULL)
	{
		Wrimsg("DDscat", "Cannot open file12, line 1044");
		return false;
	}
//
	int ll;
	const char *Format9xxx = "%8.3lf%8.3lf = beta_min, beta_max ;  NBETA =%2d\n%8.3lf%8.3lf = theta_min, theta_max; NTHETA=%2d\n%8.3lf%8.3lf = phi_min, phi_max   ;   NPHI =%2d\n";
	const char *Format9043 = "%7.4lf = ETASCA (param. controlling # of scatt. dirs used to calculate <cos> etc.\n Results averaged over %4d target orientations\n                   and %2d incident polarizations\n   aeff       wave       Q_ext     Q_abs      Q_sca     g(1)=<cos>  <cos^2>    Q_bk       Nsca\n";
	const char *Format9045 = " Results averaged over %4d target orientations\n                   and %2d incident polarizations\n   aeff       wave       Q_pha\n";
	const char *Format9046 = " Results averaged over %4d target orientations\n                   and %2d incident polarizations\n   aeff       wave       Q_pha       Q_pol       Q_cpol\n";
	const char *Format9021 = " %s\n";
	const char *Format9020 = " DDSCAT --- %s\n TARGET ---%s\n %s --- method of solution \n %s --- prescription for polarizabilities\n %s --- shape \n%7d = NAT0 = number of dipoles\n";
	fprintf(file10, Format9020, cstamp, currentTarget->GetFreeDescription(), FftEnumerator(param->Cmdfft()), AlphaEnumerator(param->Calpha()), param->TargetString(), currentTarget->Nat0());
	fprintf(file11, Format9020, cstamp, currentTarget->GetFreeDescription(), SolEnumerator(param->Cmdsol()), AlphaEnumerator(param->Calpha()), param->TargetString(), currentTarget->Nat0());
	fprintf(file12, Format9020, cstamp, currentTarget->GetFreeDescription(), FftEnumerator(param->Cmdfft()), AlphaEnumerator(param->Calpha()), param->TargetString(), currentTarget->Nat0());
	for(ll=0; ll<currentTarget->Ncomp(); ++ll)
	{
		fprintf(file10, Format9021, dielec->GetFileName(ll));
		fprintf(file11, Format9021, dielec->GetFileName(ll));
		fprintf(file12, Format9021, dielec->GetFileName(ll));
	}
	OriData *oriData = param->GetOriData();
	fprintf(file11, " wave(um)   f(cm-1)     Re(m)     Im(m)    Re(eps)   Im(eps)\n");
	fprintf(file10, Format9xxx, betmid, betmxd, oriData->GetNbeta(), thtmid, thtmxd, oriData->GetNtheta(), phimid, phimxd, oriData->GetNphi());
	fprintf(file12, Format9xxx, betmid, betmxd, oriData->GetNbeta(), thtmid, thtmxd, oriData->GetNtheta(), phimid, phimxd, oriData->GetNphi());
	if(param->Iorth() == 1)
	{
		fprintf(file10, Format9043, param->Etasca(), nori, param->Iorth());
		fprintf(file12, Format9045, nori, param->Iorth());
	}
	else
	{
		fprintf(file10, Format9043, param->Etasca(), nori, param->Iorth());
		fprintf(file12, Format9046, nori, param->Iorth());
	}
	fclose(file10);
	fclose(file12);

	return true;
}

void OutputManager::WriteDielec(real wave)
{
	const char *Format9058 = "%10.4e%11.4e%10.3e%10.3e%10.3e%10.3e\n";
	real freq = (real)1.e4 / wave;
	Complex cxen = (dielec->GetCxeps(0)).sqrt();
	fprintf(file11, Format9058, wave, freq, cxen.re, cxen.im, dielec->GetCxeps(0).re, dielec->GetCxeps(0).im);
}