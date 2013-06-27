#include "StdAfx.h"

#include "DDscatParameters.h"
#include "DDscatCommons.h"
#include "TargetManager.h"
#include "DDscatMain.h"
#include "Enumerator.h"
#include "Functions.h"
#include "DielectricManager.h"
#include "AbstractFftEngine.h"
#include "FileNamer.h"

DDscatParameters *DDscatParameters::item = NULL;

DDscatParameters::DDscatParameters(void)
{
	Init();
}

DDscatParameters::~DDscatParameters(void)
{
	Destroy();
}

DDscatParameters *DDscatParameters::GetInstance()
{
	if (!item)
		item = new DDscatParameters;

	return item;
}

void DDscatParameters::Kill()
{
	CleanDelete(item);
}

// 
// Handles the reading of input parameters from the text file as well 
// as elementart processing with those input parameters to generate arrays
// 
// Given:
//	CFLPAR	= name of file (normally 'ddscat.par')
//	NRFLD	= 0 on first call
//			= 1 if previous call set NRFLD=1 to prepare to do nearfield calculation
// 
// Returns:
//			 0 is all file read and processed Ok
//			-1 if file open error occured snumber of lines read otherwise
//
int DDscatParameters::Load(const char *cflpar, NearfieldMethod &nrfld)
{
	const char *ReaparLabel = "Reapar";
	char bbb[256];
//
// The first line
	bool bTemp = reader.Open(cflshp);
	if (!bTemp)
	{
		Wrimsg(ReaparLabel, cflpar);
		Errmsg("DDscatParameters::Load", "Cannot open parameters file.", "");
		return -1;
	}
//
// Specify whether torques are to be calculated:
//	CMDTRQ*6 = 'DOTORQ' -  calculate torque on grain
//			 = 'NOTORQ' -  do not calculate torque on grain
	reader.Read(ReaparLabel);
	int test = reader.ExtractFirstWord('\'');
	if (test == -1)
	{
		Wrimsg(ReaparLabel, " wrong definition of cmdtrq");
		return reader.GetLineCounter();
	}
	bTemp = ProcessCmtorq(reader.GetBuffer(), ReaparLabel);
	if (!bTemp)
		return reader.GetLineCounter();
//
// Define method used for iterative solution of Complex linear equations
//	CMDSOL*6	= 'PETRKP' -  Petravic & Kuo-Petravic method
//				= 'PBCGST' -  PIM BiConjugate Gradient with Stabilization
//				= 'PBCGS2' -  M.A Botcheve implementation of BiCGstab enhanced to improve convergence properties with finite precision arithmetic
//				= 'GPBICG' -  Tang et al conjugate gradient solver implemented by P.C. Chaumet & A. Rahmani
//				= 'QMRCCG' -  PIM interface to QMR solver implemented by P.C. Chaumet and A. Rahmani
//				= 'SBIGCG' -  Simplified Bi Conjugate Gradient Method, implemented by Sarkar, Yang, & Arvas 1988
	reader.Read();
	test = reader.ExtractFirstWord('\'');
	if (test == -1)
	{
		Wrimsg(ReaparLabel, " wrong definition of CMDSOL");
		return reader.GetLineCounter();
	}
	bTemp = ProcessCmdsol(reader.GetBuffer(), ReaparLabel);
	if (!bTemp)
		return reader.GetLineCounter();
//
// Define FFT method:
//	CMDFFT*6	= 'GPFAFT' -  GPFA code of Temperton
//				= 'FFTW21' -  FFTW 2.1.x code of Frigo & Johnson
//				= 'FFTMKL' -  Use DFTI from Intel MKL
//      IDIPINT	= 0 for point dipole interaction method
//				= 1 for filtered coupled dipole (FCD) interaction method
//       NRFLDB	= 0 : to skip near-field calculation of B
//				= 1 : to perform near-field calculation of B
	reader.Read();
	test = reader.ExtractFirstWord('\'');
	if (test == -1)
	{
		Wrimsg(ReaparLabel, " wrong definition of CMDFFT");
		return reader.GetLineCounter();
	}
	bTemp = ProcessCmdFFT(reader.GetBuffer(), ReaparLabel);
	if (!bTemp)
		return reader.GetLineCounter();
//
// Define prescription for computing polarizabilities:
//	CALPHA*6	= 'LATTDR' - Lattice Dispersion Relation of Draine & Goodman (1993)
//				= 'GKDLDR' - Lattice Dispersion Relation of Gutkowicz-Krusin & Draine (2004)
//              = 'FLTRCD' - Filtered Discrete Dipole treatment of Piller & Martin (1998) as corrected by
//							 Gay-Balmaz & Martin (2002) and discussed by Yurkin, Min, & Hoekstra (2010)
	reader.Read();
	test = reader.ExtractFirstWord('\'');
	if (test == -1)
	{
		Wrimsg(ReaparLabel, " wrong definition of CALPHA");
		return reader.GetLineCounter();
	}
	bTemp = ProcessCalpha(reader.GetBuffer(), ReaparLabel);
	if (!bTemp)
		return reader.GetLineCounter();
//
// Binary file flag
	reader.Read();
	test = reader.ExtractFirstWord('\'');
	if (test == -1)
	{
		Wrimsg(ReaparLabel, " wrong definition of CBINFLAG");
		return reader.GetLineCounter();
	}
	bTemp = ProcessCbinflag(reader.GetBuffer(), ReaparLabel);
	if (!bTemp)
		return reader.GetLineCounter();
//
// Read upper bound on target extent
	Wrimsg(ReaparLabel, "Use Cpp version of par file, max memory dimensions read but ignored.");
//
// Define shape:
	reader.Read();
// 
//	CSHAPE*9	= 'FROM_FILE' shape and composition data will later be read from file CFLSHP for dielectric tensors that are diagonal in Target Frame
//				= 'ANIFRMFIL' read shape and composition data from file, for general anisotropic dielectric tensors with orientation angles 
//							  THETADF,PHIDF,BETAD relative to Target Frame.
//				= 'FRMFILPBC' shape and composition data for TUC will later be read from file CFLSHP for dielectric tensors that are diagonal in Target Frame 
//							  PYD and PZD are input via ddscat.par
//				= 'ANIFILPBC' shape and composition data for TUC will later be read from file CFLSHP for general anisotropic dielectricc tensors with orientation 
//							  ngles THETADF,PHIDF,BETAD relative to Target Frame. PYD and PZD are input via ddscat.par
//				= 'ANIELLIPS' ellipsoid of anisotropic material
//				= 'ANIRCTNGL' homogeneous anisotropic rectangular target
//				= 'ANI_ELL_2' two touching anisotropic ellipsoids of materials 1-6
//				= 'ANI_ELL_3' three touching anisotropic ellipsoids of materials 1-9
//				= 'BISLINPBC' bilayer slab with periodic grid of lines parallel to z on top, with y-period/d=PYD [or, if PYD=0, a single line]
//				= 'CONELLIPS' two concentric ellipsoids of materials 1,2
//				= 'CYLINDER1' homogeneous finite cylinder
//				= 'CYLNDRCAP' homogeneous cylinder with hemispherical endcaps
//				= 'CYLNDRPBC' 1-d or 2-d array of finite cylinders
//				= 'DSKBLYPBC' 1-d or 2-d array of disk on bilayer rect. slab
//				= 'DSKRCTNGL' single disk on rectangular slab
//				= 'DSKRCTPBC' 1-d or 2-d array of disk on rectangular slab
//				= 'DW1996TAR' 13-cube target used by Draine & Weingartner 1996
//				= 'EL_IN_RCT' ellipsoid embedded in rectangular block
//				= 'ELLIPSOID' ellipsoid (homogeneous and isotropic)
//				= 'ELLIPSPBC' 1-d or 2-d array of ellipsoids
//				= 'ELLIPSO_2' two touching isotropic ellipsoids of materials 1 and 2
//				= 'ELLIPSO_3'  three touching isotropic ellipsoids of materials 1,2,3
//				= 'GAUSS_SPH'  gaussian sphere target
//				= 'HEXGONPBC'  1-d or 2-d array of hexagonal prisms
//				= 'HEX_PRISM'  homogeneous hexagonal prism
//				= 'LYRD_SLAB'  layered slab target, with up to 4 separate material layers
//				= 'LYRSLBPBC'  1-d or 2-d array of layered rect. slab targets, with up to 4 material layers
//				= 'MLTBLOCKS'  collection of cubic blocks defined by data in file 'blocks.par'
//				= 'RCTGLBLK3'  isolated target: 3 rectangular blocks with centers on x-axis
//				= 'RCTGLPRSM'  homogeneous rectangular prism
//				= 'RCTGL_PBC'  1-d or 2-d array of rectangular targets
//				= 'SLAB_HOLE'  rectangular block with cylindrical hole
//				= 'SLBHOLPBC'  1-d or 2-d array of rectangular blocks with cylindrical hole
//				= 'SPHERES_N'  multisphere target = union of N spheres
//				= 'SPHRN_PBC'  1-d or 2-d array of multisphere target
//				= 'SPHROID_2'  two touching spheroids with symmetry axes at specified angle!
//				= 'SPH_ANI_N'  multisphere target, with spheres that can have different, anisotropic, composition
//				= 'TETRAHDRN'  regular tetrahedron
//				= 'TRILYRPBC'  periodic target: 3 layer rectangular structure
//				= 'TRNGLPRSM'  triangular prism (homogeneous and isotropic)
//				= 'UNIAXICYL'  cylinder of unixaxial material

	reader.Read();
	test = reader.ExtractFirstWord('\'');
	if (test == -1)
	{
		strcat(reader.GetBuffer(), " Unrecognized shape directive");
		Wrimsg(ReaparLabel, reader.GetBuffer());
		return reader.GetLineCounter();
	}
	reader.RemoveSymbols('_', ' ');
	strcpy(targetString, reader.GetBuffer());
	bTemp = ProcessTarget(targetString, ReaparLabel);
	if (!bTemp)
		return reader.GetLineCounter();
//
// following code disabled 03.01.29
// (retain for use in future noncubic version)
// 
//   Obtain lattice anisotropy parameters DX(1-3)
//   For cubic lattice, DX(1)=DX(2)=DX(3)
//   Note that we do not require here that DX(1)*DX(2)*DX(3)=1 upon
//   input; DX is renormalized here before being returned to DDSCAT
// 
//      READ(IOPAR,FMT=*,ERR=99)DX(1),DX(2),DX(3)
//      WRITE(CMSGNM,FMT='(F8.3,F8.3,F8.3,A)')DX,' = relative lattice spacings dx,dy,dz'
//      CALL WRIMSG('REAPAR',CMSGNM)
// 
//      DELTA=(DX(1)*DX(2)*DX(3))**(1./3.)
//      DX(1)=DX(1)/DELTA
//      DX(2)=DX(2)/DELTA
//      DX(3)=DX(3)/DELTA
//      WRITE(CMSGNM,FMT='(F8.3,F8.3,F8.3,A)')DX,' = normalized lattice spacings dx,dy,dz'
//      CALL WRIMSG('REAPAR',CMSGNM)
// and replaced by following:
	dx.Set((real)1., (real)1., (real)1.);
//
// Obtain names of file(s) containing dielectric function(s)
//	NCOMP = number of different dielectric functions
//	CFLEPS(1-NCOMP) = names of files containing dielectric functions
	reader.Read();
	int j = reader.ScanInt();
	sprintf(bbb, "NCOMP=%d", j);
	Wrimsg(ReaparLabel, bbb);
//
// *** Check that NCOMP=2 if CSHAPE=UNIAXICYL *******************************
//                      3           ANIELLIPS
//                      2           ELLIPSO_2
//                      3           ELLIPSO_3
//                      6           ANI_ELL_2
//                      9           ANI_ELL_3
//                      2           CONELLIPS
//                      2           SPHROID_2
//                      2-4         LYRD_SLAB
//                      1-4         LYRSLBPBC
	if (!TargetManager::GetInstance()->AnalyseNcomp(targetString, j))
		return reader.GetLineCounter();
//
	char curFormat[64];
	DielectricManager::GetInstance()->InitFilePool(TargetManager::GetInstance()->Ncomp());
	for(j=0; j<TargetManager::GetInstance()->Ncomp(); ++j)
	{
		reader.Read();
		if ((reader.GetBuffer()[0] == '=') || (reader.GetBuffer()[1] == '='))				// just copy all loaded file names cyclically
		{
			int res = DielectricManager::GetInstance()->CycleFileNames();
			if (res)
			{
				sprintf(bbb, "Wrong equality sign command found in par file");
				Wrimsg(ReaparLabel, bbb);
				return reader.GetLineCounter();
			}
			else
				break;
		}
		else
		{
			reader.ExtractFirstWord('\'');
			DielectricManager::GetInstance()->AddDataFromFile(j, reader.GetBuffer());
		}
		sprintf(bbb, "%d %s", j, DielectricManager::GetInstance()->GetFileName(j));
		Wrimsg(ReaparLabel, bbb);
	}
//
// Debug dielectric table
//	DielectricManager::GetInstance()->Debug();
//
// Specify whether NEARFIELD calculation is to be done
// and specify fractional expansion of computational volume
//
	reader.Read(ReaparLabel);
	char aa = reader.GetFirstNonemptyCharacter();
	if ((aa == '\'') || isalpha(aa))
	{
		Wrimsg(ReaparLabel, "Possibly wrong Icomp or too uch composition files");
		return reader.GetLineCounter();
	}
	j = reader.ScanInt();
	nrfldb = NearfieldBMethodElectric;
	switch(j)
	{
	case 0:
		sprintf(bbb, "%2d = nrfld : nearfield calculation not desired", j);
		Wrimsg(ReaparLabel, bbb);
		break;

	case 1:
		sprintf(bbb, "%2d = nrfld : calculate nearfield E in specified volume", j);
		Wrimsg(ReaparLabel, bbb);
		break;

	case 2:
		sprintf(bbb, "%2d = nrfld : calculate nearfield E and B in specified volume", j);
		Wrimsg(ReaparLabel, bbb);
		j = 1;
		nrfldb = NearfieldBMethodBoth;
		break;

	default:
		sprintf(bbb, "%6d = nrfld but nrfld should be 0,1, or 2: fatal error in ddscat.par", j);
		Wrimsg(ReaparLabel, bbb);
		return reader.GetLineCounter();
		break;
	}
	nrfld = (NearfieldMethod)((int)nrfld + j);
//
// Upon return to DDSCAT:
// NRFLD	= 0: no interest in nearfield calculation
//			= 1: prepare to do nearfield calculation on next pass
//          = 2: and NRFLDB=0: perform nearfield calculation of E only 
//          = 2: and NRFLDB=1: perform nearfield calculation of both E and B
	reader.Read();
	if(nrfld != NearfieldMethodNo)
		extendxyz.Load(reader.GetBuffer(), "%lf"); 
	else
		extendxyz.Clear();

	Wrimsg(ReaparLabel, " fractional extension in -x,+x,-y,+y,-z,+z directions");
	extendxyz.Sprintf(bbb, "%6.3lf");
	Wrimsg(ReaparLabel, bbb);
//
//     Define INITIALIZATION:
//     INIT = 0 to start with |X0>=0
//            1 to obtain |X0> from 4th-order expansion in polarizability
//            2 to read |X0> from file solvp.in
//  disabled 08.03.12 v7.0.5 since we always use INIT=0
//  skip line:
//       READ(IOPAR,FMT=*,ERR=99)CLINE
//       CALL WRIMSG(' ',CLINE)
//       CWHERE='error reading INIT in ddscat.par'
//       READ(IOPAR,FMT=*,ERR=99)INIT
//       IF(INIT==0)THEN
//          CALL WRIMSG('REAPAR','INIT=0 to start with |X0>=0 (CCG method)')
//       ELSEIF(INIT==1)THEN
//          CALL WRIMSG('REAPAR','INIT=1 to obtain |X0> from 4th-order expansion (CCG)')
//       ELSEIF(INIT==2)THEN
//          CALL WRIMSG('REAPAR','INIT=2 to read |X0> from file solvp.in')
//       ELSE
//          CALL ERRMSG('FATAL','REAPAR',' Wrong value of INIT')
//       ENDIF
// !***********************************************************************

//     Define error tolerance:
//     TOL= maximum acceptable value of |Ax-E|/|E|
	const real tolLowerLimit = (real)1.e-5;
	reader.Read();
	*strchr(reader.GetBuffer(), ' ') = '\0';
	tol = reader.ScanReal(realFormat);
	if (tol < tolLowerLimit)
	{
		tol = tolLowerLimit;
		Wrimsg(ReaparLabel, "TOL corrected to be on lower limit.\n");
	}
	sprintf(bbb, "%10.3e = TOL = max. acceptable normalized residual |Ax-E|/|E|", tol);
	Wrimsg(ReaparLabel, bbb);
//
	reader.Read();
	mxiter = reader.ScanInt();
	sprintf(bbb, "%10d = MXITER", mxiter);
	Wrimsg(ReaparLabel, bbb);
//
// 2009.08.27 (BTD) add code to catch problem with
// we expect TOL to be between 1e-10 and <1.
// if outside this range, there has probably been an error reading TOL
	if ((tol < 1.e-10) || (tol > (real)1.))
	{
// Note: if the number of diel.fn. files is less than NCOMP or greater than NCOMP, this will lead to subsequent errors
//       in reading ddscat.par, which will probably affect reading of TOL
		sprintf(bbb, "%10.3le = TOL", tol);
		Wrimsg(ReaparLabel, bbb);
		Wrimsg(ReaparLabel, "Appears that there has been an error reading TOL");
		Wrimsg(ReaparLabel, "Check whether there are NCOMP diel.fn. files");
		Wrimsg(ReaparLabel," error reading ddscat.par file");
        return reader.GetLineCounter();
	}
//	
	if ((tol < tolLowerLimit) || (tol > 1.e-1))
	{
		Wrimsg(ReaparLabel, " strange value of tol ");
		return reader.GetLineCounter();
	}
//
// Define summation limit parameter GAMMA for PBC calculations summations will be carried out to distance r=2/(k*alpha)
// with suppression factor exp[-(alpha*kr)^4]=exp[-16] GAMMA is used only when JPBC=1,2, or 3
	reader.Read();
	gamma = reader.ScanReal(realFormat);
//
	if (TargetManager::GetInstance()->Jpbc() != PeriodicNo)
	{
		sprintf(bbb, "%10.3lf = GAMMA = replica dipole summation limiter for PBC", gamma);
		Wrimsg(ReaparLabel, bbb);
		if ((gamma > 1.e-1) || (gamma < 1.e-4))
		{
			Wrimsg(ReaparLabel, " strange value of gamma ");
			return reader.GetLineCounter();
		}
	}
	else
		Wrimsg(ReaparLabel, "[GAMMA is not used in present  (non-PBC) calculation]");
//
// Define angular resolution used for calculating radiation force,
//     <cos>, <cos^2>, and radiation torque (if CMTORQ='DOTORQ')
//      ETASCA = parameter controlling number of scattering angles
//             = 1 is a good choice, but can reduce this if higher accuracy is required for angular averages
//      number of scattering angles will be approximately 6*[(3+x)/ETASCA]^2
	reader.Read();
	etasca = reader.ScanReal(realFormat);
	if (etasca > 1.e-3)
		sprintf(bbb, "%10.3lf = ETASCA (parameter controlling number of scattering angles", etasca);
	else
		sprintf(bbb, "%10.3lf  is not an appropriate value for ETASCA: check ddscat.par", etasca);
	Wrimsg(ReaparLabel, bbb);
//
// Define WAVELENGTH :
//	WAVINI	= first wavelength (physical units)
//	WAVEND	= last wavelength (physical units)
//	NWAV	= number of wavelengths
//	CDIVID	= 'LIN', 'LOG', or 'INV' for equal increments in wavelength, log(wavelength), or frequency
//			= 'TAB' to read wavelengths from file 'wave.tab'
	real wavini, wavend;
	reader.Read();
	sprintf(curFormat, "%s%s%%d", realFormat, realFormat);
	sscanf(reader.GetBuffer(), curFormat, &wavini, &wavend, &nwav);
	char *ia = strchr(reader.GetBuffer(), '\'');
	char cdivid[4];
	strncpy(cdivid, ia + 1, 3);
	cdivid[3] = '\0';
	wavDivid = CdividEnumerator(cdivid);

	switch(wavDivid)
	{
	case CdividLin:
	case CdividInv:
	case CdividLog:
		wavea = (real *)malloc(nwav * sizeof(real));
		sprintf(bbb, "%3d wavelengths from  %7.4lf to %7.4lf", nwav, wavini, wavend);
		Wrimsg(ReaparLabel, bbb);
		Divide(wavDivid, wavini, wavend, nwav, wavea);
		break;
		
	case CdividTab:
		{
			FILE *file28 = fopen("wave.tab", "r");
			if (!file28)
			{
				Wrimsg("DDscatParameters", "Cannot open wave.tab file 28\n");
				return -1;
			}
			dgets(bbb, 255, file28);

			int curSize = 16;
			int deltaSize = 16;
			wavea = (real *)malloc(curSize * sizeof(real));

			nwav = 0;
			while(1)
			{
				char *ia = dgets(bbb, 255, file28);
				if (!ia)
					break;
			    wavini = reader.ScanReal(realFormat);
				wavea[nwav++] = wavini;
				if (nwav >= curSize)
				{
					wavea = (real *)realloc(wavea, (curSize + deltaSize) * sizeof(real));
					curSize += deltaSize;
				}
			}
			fclose(file28);
			if (curSize != nwav)
				wavea = (real *)realloc(wavea, nwav * sizeof(real));

			sprintf(bbb, "%d  wavelengths from %lf  to %lf", nwav, wavea[0], wavea[nwav]);
			Wrimsg(ReaparLabel, bbb);
		}
		break;

	default:
		Wrimsg(ReaparLabel, " CDIVID for wavelengths must be LIN,LOG,INV, or TAB");
		return reader.GetLineCounter();
		break;
	}
//
// Define NAMBIENT = refractive index of ambient medium
	real aefini, aefend; 
	reader.Read();
	nambient = reader.ScanInt();
	sprintf(bbb, "%7.4f = NAMBIENT = refractive index of ambient medium", nambient);
	Wrimsg(ReaparLabel, bbb);
//
// Define	AEFF	= effective radius a_eff 
//					= radius of equal-solid-volume sphere (physical units)
//			AEFINI	= first a_eff (physical units)
//			AEFEND	= last a_eff (physical units)
//			NRAD	= number of a_eff
//			CDIVID	= 'LIN', 'LOG', or 'INV' for equal increments in a_eff, log(a_eff), or 1/a_eff
//					= 'TAB' to read aeff from file 'aeff.tab'

	reader.Read();
	sscanf(reader.GetBuffer(), curFormat, &aefini, &aefend, &nrad);
	aeffa = (real *)malloc(nrad * sizeof(real));
	ia = strchr(reader.GetBuffer(), '\'');
	strncpy(cdivid, ia + 1, 3);
	cdivid[3] = '\0';
	aefDivid = CdividEnumerator(cdivid);

	switch(aefDivid)
	{
	case CdividLin:
	case CdividInv:
	case CdividLog:
		sprintf(bbb, "%d eff. radii from %lf to %lf", nrad, aefini, aefend);
		Wrimsg(ReaparLabel, bbb);
		Divide(wavDivid, aefini, aefend, nrad, aeffa);
		break;

	case CdividTab:
		{
			FILE *file28 = fopen("aeff.tab", "r");
			if (!file28)
			{
				Wrimsg("DDscatParameters", "Cannot open aeff.tab file 28\n");
				return -1;
			}
			dgets(bbb, 255, file28);					// skip one header line

			nrad = 0;
			while(1)
			{
				char *ia = dgets(bbb, 255, file28);
				if (!ia)
					break;
				nrad++;
				aefini = reader.ScanReal(realFormat);
				aeffa[nrad] = aefini;
			}
			fclose(file28);
			sprintf(bbb, "%3d' eff. radii from %7.4lf to %7.4lf", nrad, aeffa[0], aeffa[nrad-1]);
			Wrimsg(ReaparLabel, bbb);
		}
		break;

	default:
		Wrimsg(ReaparLabel, " CDIVID for radii must be LIN,LOG,INV or TAB");
		return reader.GetLineCounter();
		break;
	}
//
// Define incident polarizations (in Lab Frame)
// It is assumed that incident radiation is along x axis in Lab Frame
//
// Read Complex polarization vector CXE01_LF=e01 (normalize if necessary).
	reader.Read();
	sprintf(curFormat, "(%s,%s) (%s,%s) (%s,%s)", realFormat, realFormat, realFormat, realFormat, realFormat, realFormat);
	sscanf(reader.GetBuffer(), curFormat, &cxe01_lf.data[0].re, &cxe01_lf.data[0].im, &cxe01_lf.data[1].re, &cxe01_lf.data[1].im, &cxe01_lf.data[2].re, &cxe01_lf.data[2].im);
	real e1 = cxe01_lf.data[0].mod();
	if (e1 != (real)0.)
	{
		Wrimsg(ReaparLabel, " cxe01_lf(1) must be zero!");
		return reader.GetLineCounter();
	}
// Normalize:
	e1 = Sqrt(cxe01_lf.data[1].modSquared() + cxe01_lf.data[2].modSquared());
	cxe01_lf.data[1] /= e1;
	cxe01_lf.data[2] /= e1;
// Construct orthogonal normalized polarization vector CXE02=e02 using xhat cross e01* = e02
	cxe02_lf.data[0].clear();
	cxe02_lf.data[1] = -(cxe01_lf.data[2].conjg());
	cxe02_lf.data[2] =   cxe01_lf.data[1].conjg();
// IORTH	= 1 to calculate only for single polarization
//			= 2 to also calculate for orthogonal polarization
	reader.Read();
	iorth = reader.ScanInt();
	sprintf(bbb, "IORTH=%2d", iorth);
	Wrimsg(ReaparLabel, bbb);
//
// Specify whether or not to write ".sca" files  skip line:
// IWRKSC	= 0 to NOT write ".sca" file for each target orientation
//			= 1 to write ".sca" file for each target orientation
	reader.Read();
	iwrksc = reader.ScanInt(); 
	sprintf(bbb, "IWRKSC=%2d", iwrksc);
	Wrimsg(ReaparLabel, bbb);
//
// Specify whether or not to write ".pol" files
// 12.01.23 (BTD) * eliminate this -- no longer needed since we now include support for fast nearfield calculations, and program call readE to read the nearfield output files.
//                * now initialize IWRPOL=0, although this can then be over-ridden to set IWRPOL=1 when NRFLD=1
//
// IWRPOL	= 0 to NOT write ".pol1" and ".pol2" file for each target orientation (beta,theta)
//			= 1 to write ".pol1" and ".pol2" file for each target orientation (beta,theta)
//
//	strcpy(cwhere, "error reading IWRPOL in ddscat.par");
//	fgets(Buffer, 255, );
//	sscanf(Buffer, "%d", &iwrpol);
//	sprintf(cmsgnm, "IWRPOL=%2d", iwrpol);
//	Wrimsg(ReaparLabel, cmsgnm);
	iwrpol = 0;
//
// In the event that user set NRFLD=1 but IWRPOL=0, set IWRPOL to 1
	if((nrfld == NearfieldMethodPrepare) && (iwrpol <= 0))
	{
		iwrpol = 1;
		Wrimsg(ReaparLabel, "set iwrpol=1 because nrfld=1");
	}
//
// Read information determining target rotations
	{
		real a, b;
		int nb, nt, np;
		sprintf(curFormat, "%s%s%%d", realFormat, realFormat);
		reader.Read(ReaparLabel);
		sscanf(reader.GetBuffer(), curFormat, &a, &b, &nb);
		oridata->SetBetas(a, b, nb);
		reader.Read(ReaparLabel);
		sscanf(reader.GetBuffer(), curFormat, &a, &b, &nt);
		oridata->SetThetas(a, b, nt);
		reader.Read(ReaparLabel);
   		sscanf(reader.GetBuffer(), curFormat, &a, &b, &np);
		oridata->SetPhis(a, b, np);
		oridata->SetDimensions(nb, nt, np);
		FileNamer::GetInstance()->Init(nb * nt * np);
	}
//
// check that user has not requested more than 1000 orientations and IWRKSC=1
	if (iwrksc > 0 && (oridata->GetNbeta() * oridata->GetNtheta() * oridata->GetNphi()) > 1000)
	{
		Wrimsg(ReaparLabel, "error: if iwrksc=1, nbeta*ntheta*nphi must be .le. 1000");
		return reader.GetLineCounter();
	}
	reader.Read(ReaparLabel);
	sscanf(reader.GetBuffer(), "%d%d%d", &iwav0, &irad0, &iori0);
//
// Check that IWAV0,IRAD0,IORI0 are OK:
	if (iwav0+1 > nwav)
	{
		Wrimsg(ReaparLabel, "iwav0+1 > nwav");
		return reader.GetLineCounter();
	}
	if (irad0+1 > nrad) 
	{
		Wrimsg(ReaparLabel, "irad0+1 > nrad");
		return reader.GetLineCounter();
	}
	if (iori0+1 > oridata->GetNbeta() * oridata->GetNtheta() * oridata->GetNphi()) 
	{
		Wrimsg(ReaparLabel, "iori0+1 > nbeta * ntheta * nphi");
		return reader.GetLineCounter();
	}
//
// If NPHI>1, then set IORTH=2 regardless of value input.
	if ((iorth == 1) && (oridata->GetNphi() > 1))
	{
		iorth = 2;
		Wrimsg(ReaparLabel, "set iorth=2 since nphi>1 ");
	}
//
	switch(iorth)
	{
	case 1:
		Wrimsg(ReaparLabel, "Calculate only for single polarization ");
		break;

	case 2:
		Wrimsg(ReaparLabel, "Do orthogonal polarization for each case ");
		break;

	default:
		Wrimsg(ReaparLabel, " WRONG VALUE OF IORTH ");
		return reader.GetLineCounter();
		break;
	}
	sprintf(bbb, "%7.2lf%7.2lf  Range of BETA values ; NBETA =%d", oridata->Betami(), oridata->Betamx(), oridata->GetNbeta());
	Wrimsg(ReaparLabel, bbb);
	sprintf(bbb, "%7.2lf%7.2lf  Range of THETA values; NTHETA=%d", oridata->Thetmi(), oridata->Thetmx(), oridata->GetNtheta());
    Wrimsg(ReaparLabel, bbb);
	sprintf(bbb, "%7.2lf%7.2lf  Range of PHI values ;   NPHI =%d", oridata->Phimin(), oridata->Phimax(), oridata->GetNphi());
	Wrimsg(ReaparLabel, bbb);
//
// Convert from degrees to radians
	oridata->ToRadians();
//
// Specify elements of scattering matrix to be printed out
	reader.Read(ReaparLabel);
	nsmelts = reader.ScanInt();
//
// NSMELTS = number of elements of scattering matrix element to be printed out. 
// Must be no greater than 9 if NSMELTS is zero, then 6 "default" elements are output (ChB, no any more)
	if (nsmelts <= 0)
	{
		nsmelts = 6;
		smind1 = new int[nsmelts];
		smind2 = new int[nsmelts];
		smind1[0] = 11;
		smind1[1] = 21;
		smind1[2] = 31;
		smind1[3] = 41;
		smind1[4] = 12;
		smind1[5] = 13;
		reader.Read();
	}
	else
	{
		smind1 = new int[nsmelts];
		smind2 = new int[nsmelts];
		reader.Read();
		ia = strtok(reader.GetBuffer(), " \t");
		for(j=0; j<nsmelts; ++j)
		{
			sscanf(ia, "%d", &smind1[j]);
			ia = strtok(NULL, " \t");
		}
	}
	for(j=0; j<nsmelts; ++j)
	{
		smind2[j] = smind1[j] % 10;
		smind1[j] /= 10;
	}
//
// Specify scattering directions to be calculated
// Two options:
// CMDFRM	= 'LFRAME': specify scattering directions n in Lab Frame (where incident beam is in x-direction)
//			  THETAN,PHIN = Direction of n from n0 in calculation of the scattering matrix; cos(THETAN) is n0 \dot n
//			  and PHIN is azimuthal angle of n from Lab xy pla 
// Note: THETA1, THETA2, and PHI for each scatterin plane are entered in degrees, and immediately converted to radians.
// Arbitrary number of scattering planes may be considered.
//
// CMDFRM	= 'TFRAME': specify scattering directions in Target Frame, defined by target axes a1,a2,a3
// If JPBC = 0:
//			THETAN,PHIN = Direction of n relative to a1,a2,a3: 
//			THETAN = angle between n and a1
//			PHIN   = angle between a1-n plane and a1-a2 plane
//    JPBC = 1:
//			THETAN = Diffraction order along y_TF axis
//			PHIN   = azimuthal angle around y_TF
//	  JPBC = 2:
//			THETAN = Diffraction order along z_TF axis
//			PHIN   = azimuthal angle around z_TF
//	  JPBC = 3:
//			THETAN = Diffraction order along y_TF
//			PHIN   = Diffraction order along z_TF
//where we first run through transmitted waves 1 -> NSCAT/2 
//      and then run through reflected waves NSCAT/2+1 -> NSCAT
// 
// Three cases:
//	JPBC = 0: single isolated target.
//				specify
//				phi for scattering plane
//				thetamin, thetamax, dtheta for scattering plane
//	JPBC = 1,2: periodic in one dimension.
//				specify
//				diffraction order in direction of target periodicity
//				phimin, phimax, dphi for scattering cone
//	JPBC=3:		periodic in two dimensions
//				specify
//				difraction order in y direction and order in z direction
	real phi1, theta1, theta2, dtheta, delta;
	int nplanes, jplane, nsca0;
// 
	reader.Read();
	test = reader.ExtractFirstWord('\'');
	if (test == -1)
	{
		Wrimsg(ReaparLabel, " Error reading FRAME code");
		return reader.GetLineCounter();
	}
	bTemp = ProcessCmdfrm(reader.GetBuffer(), ReaparLabel);
	if (!bTemp)
		return reader.GetLineCounter();
//	
	reader.Read();
	nplanes = reader.ScanInt();
	switch(TargetManager::GetInstance()->Jpbc())
	{
	case PeriodicNo:
		sprintf(bbb, "%4d = number of scattering planes", nplanes);
		break;

	case PeriodicY:
	case PeriodicZ:
		sprintf(bbb, "%4d = number of scattering cones", nplanes);
		break;

	case PeriodicBoth:
		sprintf(bbb, "%4d = number of diffraction orders for transmission", nplanes);
		break;

	default:
		break;
	}
	Wrimsg(ReaparLabel, bbb);
//
	nscat = 0;
	if (nplanes > 0)
	{
		if (TargetManager::GetInstance()->Jpbc() != PeriodicBoth)
		{
			sprintf(curFormat, "%s%s%s%s", realFormat, realFormat, realFormat, realFormat);
			for(jplane=0; jplane<nplanes; ++jplane)
			{
				reader.Read();
				sscanf(reader.GetBuffer(), curFormat, &phi1, &theta1, &theta2, &dtheta);
				if (TargetManager::GetInstance()->Jpbc() == PeriodicNo)
					sprintf(bbb, "%7.1lf%7.1lf%7.1f = phi, theta_min, theta_max for scatt. plane%4d", phi1, theta1, theta2, jplane);
				else
					sprintf(bbb, "%7.1lf%7.1lf%7.1f = order, zeta_min, zeta_max for scattering cone", phi1, theta1, theta2);
				Wrimsg(ReaparLabel, bbb);
				if (TargetManager::GetInstance()->Jpbc() == PeriodicNo)
				{
					if((Fabs(theta1-theta2) > (real)0.) && (dtheta == (real)0.))
					{
						Wrimsg(ReaparLabel, "DTHETA=0 in ddscat.par!");
						return reader.GetLineCounter();
					}
				}
// ! Convert to radians
				phi1 /= Degrad;
				theta1 /= Degrad;
				theta2 /= Degrad;
				dtheta /= Degrad;
// ! Allow for possibility that user did not enter correct sign for
// ! DTHETA (if I did it, others will too...)
				if (theta2 < theta1)
					dtheta = -Fabs(dtheta);
// ! compute theta values for this scattering plane/cone
				nsca0 = 1;
				delta = theta2 - theta1;
				if (delta)
					nsca0 = 1 + nint_(delta/dtheta);
// ChB: ! realloc to have sufficient space
				thetan = (real *)realloc(thetan, (nscat + nsca0)*sizeof(real));
				phin = (real *)realloc(phin, (nscat + nsca0)*sizeof(real));
				thetan[nscat] = theta1;
				if(nsca0 > 2)
				{
					for(j=1; j<nsca0-1; ++j)
					{
						thetan[nscat + j] = theta1 + j*dtheta;
					}
				}
				thetan[nscat + nsca0-1] = theta2;
				for(j=0; j<nsca0; ++j)
				{
					phin[nscat + j] = phi1;
				}
				nscat += nsca0;
				sprintf(bbb, "%4d = number of scattering angles in this scattering %s", nsca0, ((TargetManager::GetInstance()->Jpbc() == PeriodicNo) ? "plane" : "cone"));
				Wrimsg(ReaparLabel, bbb);
			}
		}
		else
		{
			real aa, bb;
			sprintf(curFormat, "%s%s", realFormat, realFormat);
			phin = (real *)malloc(2 * nplanes * sizeof(real));
			thetan = (real *)malloc(2 * nplanes * sizeof(real));
			for(jplane=0; jplane<nplanes; ++jplane)
			{
				reader.Read();
				sscanf(reader.GetBuffer(), curFormat, &aa, &bb);
				phin[nscat] = aa;
				thetan[nscat] = bb;
				++nscat;
			}
			for(j=0; j<nscat; ++j)
			{
				phin[j + nscat] = phin[j];
				thetan[j + nscat] = thetan[j];
			}
			nscat = nscat + nscat;
		}
	}
	reader.Close();

	return 0;
}

void DDscatParameters::Init(void)
{
	cmdtrq = TorqMethod_End;
	cmdsol = SolMethod_End;
	cmdfft = FftMethod_End;
	fftEngine = NULL;
	calpha = AlphaMethod_End;
	cbinflag = BinflagMethod_End;
	cshape = TargetType_End;
	cmdfrm = FrameCode_End;
	wavDivid = CdividEnd;
	aefDivid = CdividEnd;
	cflshp = (char *)malloc(10);
	strcpy(cflshp, "shape.dat");
	numShpar = 0;
	wavea = NULL;
	thetan = NULL;
	phin = NULL;
	aeffa = NULL;
	smind1 = smind2 = NULL;
	oridata = new OriData;
}

void DDscatParameters::Clear2(void)
{
	wavea = NULL;
	thetan = NULL;
	phin = NULL;
	aeffa = NULL;
}

void DDscatParameters::Destroy(void)
{
	if (cflshp)
		free(cflshp);
	if (wavea)
		free(wavea);
	if (thetan)
		free(thetan);
	if (phin)
		free(phin);
	if (aeffa)
		free(aeffa);
	Clear2();
	CleanDelete2(smind1);
	CleanDelete2(smind2);
	delete oridata;
}

void DDscatParameters::OutParam(const char *Label)
{
	char Buffer[256];
	for(int j=0; j<numShpar; ++j)
	{
		sprintf(Buffer, "%lf = sphar[%d]", shpar[j], j);
		Wrimsg(Label, Buffer);
	}
}

void DDscatParameters::Orient()
{
	oridata->Orient();
}

void DDscatParameters::Scavec(Vect3<real> *ensc, Vect3<real> *em1, Vect3<real> *em2)
{
/* **
Given:
      MXSCA=dimension information for ENSC,PHIN,THETAN
      NSCAT=number of scattering directions
      THETAN(1-NSCAT)=scattering angles theta
      PHIN(1-NSCAT)=scattering angles phi
      CXE01(1-3)=Complex polarization vector 1 (phi=0 direction is defined by x,y plane (Lab Frame), where incident radiation propagates along the x axis.

Returns:
      ENSC(1-3,1-NSCAT)=scattering vectors in Lab Frame
      EM1(1-3,1-NSCAT)=scattered pol vectors parallel to scat. plane in Lab Frame
      EM2(1-3,1-NSCAT)=scattered pol vectors perp. to scat. plane in Lab Frame
 It is assumed that incident propagation vector is in x-direction in Lab Frame

 History:
 96.11.06 (BTD): Changed definition of scattering angle phi
                 Previously, phi was measured from plane containing
                 incident k vector (i.e., Lab x-axis) and Re(CXE01)
                 Henceforth, phi is measured from Lab x,y plane.
 10.01.30 (BTD): cosmetic changes
 end history
Copyright (C) 1993,1996,2010 B.T. Draine and P.J. Flatau
This code is covered by the GNU General Public License.
** */

	for(int i=0; i<nscat; ++i)
	{
		real cosphi = Cos(phin[i]);
		real sinphi = Sin(phin[i]);
		real costhe = Cos(thetan[i]);
		real sinthe = Sin(thetan[i]);
		ensc[i].Set(costhe, sinthe*cosphi, sinthe*sinphi);
		em1[i].Set(-sinthe, costhe*cosphi, costhe*sinphi);
		em2[i].Set((real)0., -sinphi, cosphi);
	}
}

int DDscatParameters::LoadXml(const char *fileName, NearfieldMethod &nrfld)
{
	const real onex_ = (real)1.;
	const int maxDeepLevel = 4;
	const int rootNodeLevel = 0;
	const char *ReaparLabel = "XmlReapar";
	char bbb[256];
//
	xmlDocPtr doc = xmlParseFile(fileName);
	if(doc == NULL)
	{
		fprintf(stderr, "Cannot open xml parameters file.\n");
		return -1;
	}
	else
		fprintf(stdout, "Using Xml parameter file.\n");
//
	xmlNodePtr theNodes[maxDeepLevel];
	theNodes[rootNodeLevel] = xmlDocGetRootElement(doc);
	if (!theNodes[rootNodeLevel])
	{
		Wrimsg(ReaparLabel, "Cannot find root node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	sprintf(bbb, " ========= Parameter file for v%s", xmlGetProp(theNodes[rootNodeLevel], (const xmlChar *)"ver"));
	Wrimsg(ReaparLabel, bbb);
	
	theNodes[rootNodeLevel + 1] = XmlFindChild(theNodes[rootNodeLevel], (const xmlChar *)"Preliminaries");
	if (!theNodes[rootNodeLevel + 1])
	{
		Wrimsg(ReaparLabel, "Cannot find root->Preliminaries node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	Wrimsg(ReaparLabel, "**** Preliminaries ****\n\n");
//	
	theNodes[rootNodeLevel + 2] = XmlFindChild(theNodes[rootNodeLevel + 1], (const xmlChar *)"Cmtorq");
	if (!theNodes[rootNodeLevel + 2])
	{
		Wrimsg(ReaparLabel, "Cannot find root->Preliminaries->Cmtorq node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	strcpy(bbb, (char *)xmlGetProp(theNodes[rootNodeLevel + 2], (const xmlChar *)"Value"));
	bool bTemp = ProcessCmtorq(bbb, ReaparLabel);
	if (!bTemp)
	{
		xmlFreeDoc(doc);
		return -1;
	}
//
// Define method used for iterative solution of Complex linear equations
//	CMDSOL*6	= 'PETRKP' -  Petravic & Kuo-Petravic method
//				= 'PBCGST' -  PIM BiConjugate Gradient with Stabilization
//				= 'PBCGS2' -  M.A Botcheve implementation of BiCGstab enhanced to improve convergence properties with finite precision arithmetic
//				= 'GPBICG' -  Tang et al conjugate gradient solver implemented by P.C. Chaumet & A. Rahmani
//				= 'QMRCCG' -  PIM interface to QMR solver implemented by P.C. Chaumet and A. Rahmani
	theNodes[rootNodeLevel + 2] = XmlFindChild(theNodes[rootNodeLevel + 1], (const xmlChar *)"Cmdsol");
	if (!theNodes[rootNodeLevel + 2])
	{
		Wrimsg(ReaparLabel, "Cannot find root->Preliminaries->Cmdsol node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	strcpy(bbb, (char *)xmlGetProp(theNodes[rootNodeLevel + 2], (const xmlChar *)"Value"));
	bTemp = ProcessCmdsol(bbb, ReaparLabel);
	if (!bTemp)
	{
		xmlFreeDoc(doc);
		return -1;
	}		
//
// Define FFT method:
//	CMDFFT*6	= 'GPFAFT' -  GPFA code of Temperton
//				= 'FFTW21' -  FFTW 2.1.x code of Frigo & Johnson
//				= 'FFTMKL' -  Use DFTI from Intel MKL
	theNodes[rootNodeLevel + 2] = XmlFindChild(theNodes[rootNodeLevel + 1], (const xmlChar *)"CmdFFT");
	if (!theNodes[rootNodeLevel + 2])
	{
		Wrimsg(ReaparLabel, "Cannot find root->Preliminaries->CmdFFT node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	strcpy(bbb, (char *)xmlGetProp(theNodes[rootNodeLevel + 2], (const xmlChar *)"Value"));
	bTemp = ProcessCmdFFT(bbb, ReaparLabel);
	if (!bTemp)
	{
		xmlFreeDoc(doc);
		return -1;
	}
//
// Define prescription for computing polarizabilities:
//	CALPHA*6	= 'LATTDR' - Lattice Dispersion Relation of Draine & Goodman (1993)
//				= 'GKDLDR' - Lattice Dispersion Relation of Gutkowicz-Krusin & Draine (2004)
	theNodes[rootNodeLevel + 2] = XmlFindChild(theNodes[rootNodeLevel + 1], (const xmlChar *)"Calpha");
	if (!theNodes[rootNodeLevel + 2])
	{
		Wrimsg(ReaparLabel, "Cannot find root->Preliminaries->Calpha node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	strcpy(bbb, (char *)xmlGetProp(theNodes[rootNodeLevel + 2], (const xmlChar *)"Value"));
	bTemp = ProcessCalpha(bbb, ReaparLabel);
	if (!bTemp)
	{
		xmlFreeDoc(doc);
		return -1;
	}
//
// Binary file flag
	theNodes[rootNodeLevel + 2] = XmlFindChild(theNodes[rootNodeLevel + 1], (const xmlChar *)"Cbinflag");
	if (!theNodes[rootNodeLevel + 2])
	{
		Wrimsg(ReaparLabel, "Cannot find root->Preliminaries->Cbinflag node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	strcpy(bbb, (char *)xmlGetProp(theNodes[rootNodeLevel + 2], (const xmlChar *)"Value"));
	bTemp = ProcessCbinflag(bbb, ReaparLabel);
	if (!bTemp)
	{
		xmlFreeDoc(doc);
		return -1;
	}
//
// Read upper bound on target extent
	Wrimsg(ReaparLabel, "Use Cpp xml version of par file, max memory dimensions are ignored.");
// 
//	CSHAPE*9	= 'FROM_FILE' shape and composition data will later be read from file CFLSHP for dielectric tensors that are diagonal in Target Frame
//				= 'ANIFRMFIL' read shape and composition data from file, for general anisotropic dielectric tensors with orientation angles 
//							  THETADF,PHIDF,BETAD relative to Target Frame.
//				= 'FRMFILPBC' shape and composition data for TUC will later be read from file CFLSHP for dielectric tensors that are diagonal in Target Frame 
//							  PYD and PZD are input via ddscat.par
//				= 'ANIFILPBC' shape and composition data for TUC will later be read from file CFLSHP for general anisotropic dielectricc tensors with orientation 
//							  ngles THETADF,PHIDF,BETAD relative to Target Frame. PYD and PZD are input via ddscat.par
//				= 'ANIELLIPS' ellipsoid of anisotropic material
//				= 'ANIRCTNGL' homogeneous anisotropic rectangular target
//				= 'ANI_ELL_2' two touching anisotropic ellipsoids of materials 1-6
//				= 'ANI_ELL_3' three touching anisotropic ellipsoids of materials 1-9
//				= 'BISLINPBC' bilayer slab with periodic grid of lines parallel to z on top, with y-period/d=PYD [or, if PYD=0, a single line]
//				= 'CONELLIPS' two concentric ellipsoids of materials 1,2
//				= 'CYLINDER1' homogeneous finite cylinder
//				= 'CYLNDRCAP' homogeneous cylinder with hemispherical endcaps
//				= 'CYLNDRPBC' 1-d or 2-d array of finite cylinders
//				= 'DSKBLYPBC' 1-d or 2-d array of disk on bilayer rect. slab
//				= 'DSKRCTNGL' single disk on rectangular slab
//				= 'DSKRCTPBC' 1-d or 2-d array of disk on rectangular slab
//				= 'DW1996TAR' 13-cube target used by Draine & Weingartner 1996
//				= 'EL_IN_RCT' ellipsoid embedded in rectangular block
//				= 'ELLIPSOID' ellipsoid (homogeneous and isotropic)
//				= 'ELLIPSPBC' 1-d or 2-d array of ellipsoids
//				= 'ELLIPSO_2' two touching isotropic ellipsoids of materials 1 and 2
//				= 'ELLIPSO_3'  three touching isotropic ellipsoids of materials 1,2,3
//				= 'GAUSS_SPH'  gaussian sphere target
//				= 'HEXGONPBC'  1-d or 2-d array of hexagonal prisms
//				= 'HEX_PRISM'  homogeneous hexagonal prism
//				= 'LYRD_SLAB'  layered slab target, with up to 4 separate material layers
//				= 'LYRSLBPBC'  1-d or 2-d array of layered rect. slab targets, with up to 4 material layers
//				= 'MLTBLOCKS'  collection of cubic blocks defined by data in file 'blocks.par'
//				= 'RCTGLBLK3'  isolated target: 3 rectangular blocks with centers on x-axis
//				= 'RCTGLPRSM'  homogeneous rectangular prism
//				= 'RCTGL_PBC'  1-d or 2-d array of rectangular targets
//				= 'SLAB_HOLE'  rectangular block with cylindrical hole
//				= 'SLBHOLPBC'  1-d or 2-d array of rectangular blocks with cylindrical hole
//				= 'SPHERES_N'  multisphere target = union of N spheres
//				= 'SPHRN_PBC'  1-d or 2-d array of multisphere target
//				= 'SPHROID_2'  two touching spheroids with symmetry axes at specified angle!
//				= 'SPH_ANI_N'  multisphere target, with spheres that can have different, anisotropic, composition
//				= 'TETRAHDRN'  regular tetrahedron
//				= 'TRILYRPBC'  periodic target: 3 layer rectangular structure
//				= 'TRNGLPRSM'  triangular prism (homogeneous and isotropic)
//				= 'UNIAXICYL'  cylinder of unixaxial material
	theNodes[rootNodeLevel + 1] = XmlFindChild(theNodes[rootNodeLevel], (const xmlChar *)"TargetGeometryAndComposition");
	if (!theNodes[rootNodeLevel + 1])
	{
		Wrimsg(ReaparLabel, "Cannot find root->TargetGeometryAndComposition node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	Wrimsg(ReaparLabel, "**** Target Geometry and Composition ****\n\n");
//	
	theNodes[rootNodeLevel + 2] = XmlFindChild(theNodes[rootNodeLevel + 1], (const xmlChar *)"Cshape");	
	if (!theNodes[rootNodeLevel + 2])
	{
		Wrimsg(ReaparLabel, "Cannot find root->TargetGeometryAndComposition->Cshape node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	strcpy(bbb, (char *)xmlGetProp(theNodes[rootNodeLevel + 2], (const xmlChar *)"Name"));	
	bTemp = ProcessTarget(bbb, ReaparLabel);
	if (!bTemp)
	{
		xmlFreeDoc(doc);
		return -1;
	}
//
// Read shape parameters
	theNodes[rootNodeLevel + 3] = theNodes[rootNodeLevel + 2]->children;
	numShpar = LoadShparFromXml(theNodes[rootNodeLevel + 3]);
//
// following code disabled 03.01.29
// (retain for use in future noncubic version)
// 
//   Obtain lattice anisotropy parameters DX(1-3)
//   For cubic lattice, DX(1)=DX(2)=DX(3)
//   Note that we do not require here that DX(1)*DX(2)*DX(3)=1 upon
//   input; DX is renormalized here before being returned to DDSCAT
// 
//      READ(IOPAR,FMT=*,ERR=99)DX(1),DX(2),DX(3)
//      WRITE(CMSGNM,FMT='(F8.3,F8.3,F8.3,A)')DX,' = relative lattice spacings dx,dy,dz'
//      CALL WRIMSG('REAPAR',CMSGNM)
// 
//      DELTA=(DX(1)*DX(2)*DX(3))**(1./3.)
//      DX(1)=DX(1)/DELTA
//      DX(2)=DX(2)/DELTA
//      DX(3)=DX(3)/DELTA
//      WRITE(CMSGNM,FMT='(F8.3,F8.3,F8.3,A)')DX,' = normalized lattice spacings dx,dy,dz'
//      CALL WRIMSG('REAPAR',CMSGNM)
// and replaced by following:
	dx.Set(onex_, onex_, onex_);
//
// Obtain names of file(s) containing dielectric function(s)
//	NCOMP = number of different dielectric functions
//	CFLEPS(1-NCOMP) = names of files containing dielectric functions
	theNodes[rootNodeLevel + 2] = XmlFindChild(theNodes[rootNodeLevel + 1], (const xmlChar *)"Ncomp");	
	if (!theNodes[rootNodeLevel + 2])
	{
		Wrimsg(ReaparLabel, "Cannot find root->TargetGeometryAndComposition->Ncomp node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	int j = XmlGetIntProperty(theNodes[rootNodeLevel + 2], (const xmlChar *)"Amount");
	sprintf(bbb, "NCOMP=%d", j);
	Wrimsg(ReaparLabel, bbb);
//
// *** Check that NCOMP=2 if CSHAPE=UNIAXICYL *******************************
//                      3           ANIELLIPS
//                      2           ELLIPSO_2
//                      3           ELLIPSO_3
//                      6           ANI_ELL_2
//                      9           ANI_ELL_3
//                      2           CONELLIPS
//                      2           SPHROID_2
//                      2-4         LYRD_SLAB
//                      1-4         LYRSLBPBC
	if (!TargetManager::GetInstance()->AnalyseNcomp(targetString, j))
		return reader.GetLineCounter();	
//
	char curFormat[64];
	theNodes[rootNodeLevel + 3] = theNodes[rootNodeLevel + 2]->children;
	DielectricManager::GetInstance()->InitFilePool(TargetManager::GetInstance()->Ncomp());
	for(j=0; j<TargetManager::GetInstance()->Ncomp(); ++j)
	{
		if (!xmlStrcmp(theNodes[rootNodeLevel + 3]->name, (const xmlChar *)"Dielec"))
		{
			xmlChar *uri = xmlGetProp(theNodes[rootNodeLevel + 3], (const xmlChar *)"File");
			if ((uri[0] == '=') || (uri[1] == '='))						// just copy all loaded file names cyclically
			{
				int res = DielectricManager::GetInstance()->CycleFileNames();
				if (res)
				{
					sprintf(bbb, "Wrong equality sign command found in par file");
					Wrimsg(ReaparLabel, bbb);
					return reader.GetLineCounter();
				}
				else
					break;
			}
			else
			{
				DielectricManager::GetInstance()->AddDataFromFile(j, (char *)uri);
			}
		}
		else
			theNodes[rootNodeLevel + 3] = theNodes[rootNodeLevel + 3]->next;
		sprintf(bbb, "%d %s", j, DielectricManager::GetInstance()->GetFileName(j));
		Wrimsg(ReaparLabel, bbb);
	}
//
// Specify whether NEARFIELD calculation is to be done
// and specify fractional expansion of computational volume
//
	theNodes[rootNodeLevel + 1] = XmlFindChild(theNodes[rootNodeLevel], (const xmlChar *)"NearfieldCalculation");	
	if (!theNodes[rootNodeLevel + 1])
	{
		Wrimsg(ReaparLabel, "Cannot find root->NearfieldCalculation node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	j = XmlGetIntProperty(theNodes[rootNodeLevel + 1], (const xmlChar *)"Nrfld");
	nrfldb = NearfieldBMethodElectric;
	switch(j)
	{
	case 0:
		sprintf(bbb, "%2d = nrfld : nearfield calculation not desired", j);
		Wrimsg(ReaparLabel, bbb);
		break;

	case 1:
		sprintf(bbb, "%2d = nrfld : calculate nearfield E in specified volume", j);
		Wrimsg(ReaparLabel, bbb);
		break;

	case 2:
		sprintf(bbb, "%2d = nrfld : calculate nearfield E and B in specified volume", j);
		Wrimsg(ReaparLabel, bbb);
		j = 1;
		nrfldb = NearfieldBMethodBoth;
		break;

	default:
		sprintf(bbb, "%6d = nrfld but nrfld should be 0,1, or 2: fatal error in ddscat.par", j);
		Wrimsg(ReaparLabel, bbb);
		return reader.GetLineCounter();
		break;
	}
	nrfld = (NearfieldMethod)((int)nrfld + j);
//
// Upon return to DDSCAT:
// NRFLD	= 0: no interest in nearfield calculation
//			= 1: prepare to do nearfield calculation on next pass
//			= 2: perform nearfield calculation
	if(nrfld != NearfieldMethodNo)
	{
		const char *aa[] = { "Xm", "Xp", "Ym", "Yp", "Zm", "Zp" };
		theNodes[rootNodeLevel + 2] = XmlFindChild(theNodes[rootNodeLevel + 1], (const xmlChar *)"Extendxyz");
		for(int k=0; k<6; ++k)
		{
			extendxyz.data[k] = XmlGetDoubleProperty(theNodes[rootNodeLevel + 2], (const xmlChar *)aa[k]);
		}
	}
	else
		extendxyz.Clear();
	Wrimsg(ReaparLabel, " fractional extension in -x,+x,-y,+y,-z,+z directions");
	extendxyz.Sprintf(bbb, "%6.3lf");
	Wrimsg(ReaparLabel, bbb);
//
//     Define INITIALIZATION:
//     INIT = 0 to start with |X0>=0
//            1 to obtain |X0> from 4th-order expansion in polarizability
//            2 to read |X0> from file solvp.in
//  disabled 08.03.12 v7.0.5 since we always use INIT=0
//  skip line:
//       READ(IOPAR,FMT=*,ERR=99)CLINE
//       CALL WRIMSG(' ',CLINE)
//       CWHERE='error reading INIT in ddscat.par'
//       READ(IOPAR,FMT=*,ERR=99)INIT
//       IF(INIT==0)THEN
//          CALL WRIMSG('REAPAR','INIT=0 to start with |X0>=0 (CCG method)')
//       ELSEIF(INIT==1)THEN
//          CALL WRIMSG('REAPAR','INIT=1 to obtain |X0> from 4th-order expansion (CCG)')
//       ELSEIF(INIT==2)THEN
//          CALL WRIMSG('REAPAR','INIT=2 to read |X0> from file solvp.in')
//       ELSE
//          CALL ERRMSG('FATAL','REAPAR',' Wrong value of INIT')
//       ENDIF
// !***********************************************************************

//     Define error tolerance:
//     TOL= maximum acceptable value of |Ax-E|/|E|
	const real tolLowerLimit = (real)1.e-5;
	Wrimsg(ReaparLabel, "**** Error Tolerance ****\n");
	theNodes[rootNodeLevel + 1] = XmlFindChild(theNodes[rootNodeLevel], (const xmlChar *)"Tol");
	if (!theNodes[rootNodeLevel + 1])
	{
		Wrimsg(ReaparLabel, "Cannot find root->Tol node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	tol = XmlGetDoubleProperty(theNodes[rootNodeLevel + 1], (const xmlChar *)"Value");
	if (tol < tolLowerLimit)
	{
		tol = tolLowerLimit;
		Wrimsg(ReaparLabel, "TOL corrected to be on lower limit.\n");
	}
	sprintf(bbb, "%10.3e = TOL = max. acceptable normalized residual |Ax-E|/|E|", tol);
	Wrimsg(ReaparLabel, bbb);
//
	Wrimsg(ReaparLabel, "**** maximum number of iterations allowed ****\n");
	theNodes[rootNodeLevel + 1] = XmlFindChild(theNodes[rootNodeLevel], (const xmlChar *)"Mxiter");
	if (!theNodes[rootNodeLevel + 1])
	{
		Wrimsg(ReaparLabel, "Cannot find root->Mxiter node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	mxiter = XmlGetIntProperty(theNodes[rootNodeLevel + 1], (const xmlChar *)"Value");
	sprintf(bbb, "%10d = MXITER\n", mxiter);
	Wrimsg(ReaparLabel, bbb);
//
// 2009.08.27 (BTD) add code to catch problem with
// we expect TOL to be between 1e-10 and <1.
// if outside this range, there has probably been an error reading TOL
	if(tol < 1.e-10 || tol > onex_)
	{
// Note: if the number of diel.fn. files is less than NCOMP or greater than NCOMP, this will lead to subsequent errors
//       in reading ddscat.par, which will probably affect reading of TOL
		sprintf(bbb, "%10.3le = TOL", tol);
		Wrimsg(ReaparLabel, bbb);
		Wrimsg(ReaparLabel, "Appears that there has been an error reading TOL");
		Wrimsg(ReaparLabel, "Check whether there are NCOMP diel.fn. files");
		Wrimsg(ReaparLabel," error reading ddscat.par file");
		xmlFreeDoc(doc);
        return -1;
	}
//	
	sprintf(bbb, "%10.3le  = TOL = max. acceptable normalized residue |Ax-E|/|E|", tol);
	Wrimsg(ReaparLabel, bbb);
	if (tol < tolLowerLimit || tol > 1.e-1)
	{
		Wrimsg(ReaparLabel, "strange value of tol ");
		return -1;
	}
//
// Define summation limit parameter GAMMA for PBC calculations summations will be carried out to distance r=2/(k*alpha)
// with suppression factor exp[-(alpha*kr)^4]=exp[-16] GAMMA is used only when JPBC=1,2, or 3 skip line
	theNodes[rootNodeLevel + 1] = XmlFindChild(theNodes[rootNodeLevel], (const xmlChar *)"Gamma");
	if (!theNodes[rootNodeLevel + 1])
	{
		Wrimsg(ReaparLabel, "Cannot find root->Gamma node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	gamma = XmlGetDoubleProperty(theNodes[rootNodeLevel + 1], (const xmlChar *)"Value");
	Wrimsg(ReaparLabel, "**** Interaction cutoff parameter for PBC calculations ****\n");
	if (TargetManager::GetInstance()->Jpbc() != PeriodicNo)
	{
		sprintf(bbb, "%10.3lf = GAMMA = replica dipole summation limiter for PBC", gamma);
		Wrimsg(ReaparLabel, bbb);
		if ((gamma > 1.e-1) || (gamma < 1.e-4))
		{
			Wrimsg(ReaparLabel, " strange value of gamma ");
			return reader.GetLineCounter();
		}
	}
	else
		Wrimsg(ReaparLabel, "[GAMMA is not used in present  (non-PBC) calculation]");
//
// Define angular resolution used for calculating radiation force,
//     <cos>, <cos^2>, and radiation torque (if CMTORQ='DOTORQ')
//      ETASCA = parameter controlling number of scattering angles
//             = 1 is a good choice, but can reduce this if higher accuracy is required for angular averages
//      number of scattering angles will be approximately 6*[(3+x)/ETASCA]^2
	theNodes[rootNodeLevel + 1] = XmlFindChild(theNodes[rootNodeLevel], (const xmlChar *)"Etasca");
	if (!theNodes[rootNodeLevel + 1])
	{
		Wrimsg(ReaparLabel, "Cannot find root->Etasca node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	etasca = XmlGetDoubleProperty(theNodes[rootNodeLevel + 1], (const xmlChar *)"Value");
	Wrimsg(ReaparLabel, "**** Angular resolution for calculation of <cos>, etc. ****\n");	
	if (etasca > 1.e-3)
		sprintf(bbb, "%10.3lf = ETASCA (parameter controlling number of scattering angles", etasca);
	else
		sprintf(bbb, "%10.3lf  is not an appropriate value for ETASCA: check ddscat.par", etasca);
	Wrimsg(ReaparLabel, bbb);
//
// Define WAVELENGTH :
//	WAVINI	= first wavelength (physical units)
//	WAVEND	= last wavelength (physical units)
//	NWAV	= number of wavelengths
//	CDIVID	= 'LIN', 'LOG', or 'INV' for equal increments in wavelength, log(wavelength), or frequency
//			= 'TAB' to read wavelengths from file 'wave.tab'
	char cdivid[4];
	theNodes[rootNodeLevel + 1] = XmlFindChild(theNodes[rootNodeLevel], (const xmlChar *)"VacuumWavelengths");
	if (!theNodes[rootNodeLevel + 1])
	{
		Wrimsg(ReaparLabel, "Cannot find root->VacuumWavelengths node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	Wrimsg(ReaparLabel, "**** Vacuum wavelengths (micron) ****\n");
	real wavini = XmlGetDoubleProperty(theNodes[rootNodeLevel + 1], (const xmlChar *)"First");
	real wavend = XmlGetDoubleProperty(theNodes[rootNodeLevel + 1], (const xmlChar *)"Last");
	nwav = XmlGetIntProperty(theNodes[rootNodeLevel + 1], (const xmlChar *)"HowMany");
	strcpy(cdivid, (char *)xmlGetProp(theNodes[rootNodeLevel + 1], (const xmlChar *)"How"));
	wavDivid = CdividEnumerator(cdivid);
	switch(wavDivid)
	{
	case CdividLin:
	case CdividInv:
	case CdividLog:
		wavea = (real *)malloc(nwav * sizeof(real));
		sprintf(bbb, "%3d wavelengths from  %7.4lf to %7.4lf", nwav, wavini, wavend);
		Wrimsg(ReaparLabel, bbb);
		Divide(wavDivid, wavini, wavend, nwav, wavea);
		break;
		
	case CdividTab:
		{
			FILE *file28 = fopen("wave.tab", "r");
			if (!file28)
			{
				Wrimsg("DDscatParameters", "Cannot open wave.tab file 28\n");
				return -1;
			}
			dgets(bbb, 255, file28);

			int curSize = 16;
			int deltaSize = 16;
			wavea = (real *)malloc(curSize * sizeof(real));

			nwav = 0;
			while(1)
			{
				char *ia = dgets(bbb, 255, file28);
				if (!ia)
					break;
			    wavini = reader.ScanReal(realFormat);
				wavea[nwav++] = wavini;
				if (nwav >= curSize)
				{
					wavea = (real *)realloc(wavea, (curSize + deltaSize) * sizeof(real));
					curSize += deltaSize;
				}
			}
			fclose(file28);
			if (curSize != nwav)
				wavea = (real *)realloc(wavea, nwav * sizeof(real));

			sprintf(bbb, "%d  wavelengths from %lf  to %lf", nwav, wavea[0], wavea[nwav]);
			Wrimsg(ReaparLabel, bbb);
		}
		break;

	default:
		Wrimsg(ReaparLabel, " CDIVID for wavelengths must be LIN,LOG,INV, or TAB");
		return reader.GetLineCounter();
		break;
	}
//
// Define NAMBIENT = refractive index of ambient medium
	theNodes[rootNodeLevel + 1] = XmlFindChild(theNodes[rootNodeLevel], (const xmlChar *)"Nambient");
	if (!theNodes[rootNodeLevel + 1])
	{
		Wrimsg(ReaparLabel, "Cannot find root->Nambient node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	nambient = XmlGetDoubleProperty(theNodes[rootNodeLevel + 1], (const xmlChar *)"Value");
	Wrimsg(ReaparLabel, "**** Refractive index of ambient medium\n");
	sprintf(bbb, "%7.4f = NAMBIENT = refractive index of ambient medium\n", nambient);
	Wrimsg(ReaparLabel, bbb);
//
// Define	AEFF	= effective radius a_eff 
//					= radius of equal-solid-volume sphere (physical units)
//			AEFINI	= first a_eff (physical units)
//			AEFEND	= last a_eff (physical units)
//			NRAD	= number of a_eff
//			CDIVID	= 'LIN', 'LOG', or 'INV' for equal increments in a_eff, log(a_eff), or 1/a_eff
//					= 'TAB' to read aeff from file 'aeff.tab'
	theNodes[rootNodeLevel + 1] = XmlFindChild(theNodes[rootNodeLevel], (const xmlChar *)"Aeff");
	if (!theNodes[rootNodeLevel + 1])
	{
		Wrimsg(ReaparLabel, "Cannot find root->Aeff node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	Wrimsg(ReaparLabel, "**** Effective Radii (micron) ****\n");
	real aefini = XmlGetDoubleProperty(theNodes[rootNodeLevel + 1], (const xmlChar *)"First");
	real aefend = XmlGetDoubleProperty(theNodes[rootNodeLevel + 1], (const xmlChar *)"Last");
	nrad = XmlGetIntProperty(theNodes[rootNodeLevel + 1], (const xmlChar *)"HowMany");
	strcpy(cdivid, (char *)xmlGetProp(theNodes[rootNodeLevel + 1], (const xmlChar *)"How"));
	aefDivid = CdividEnumerator(cdivid);
//
	switch(aefDivid)
	{
	case CdividLin:
	case CdividInv:
	case CdividLog:
		sprintf(bbb, "%d eff. radii from %lf to %lf", nrad, aefini, aefend);
		Wrimsg(ReaparLabel, bbb);
		Divide(wavDivid, aefini, aefend, nrad, aeffa);
		break;

	case CdividTab:
		{
			FILE *file28 = fopen("aeff.tab", "r");
			if (!file28)
			{
				Wrimsg("DDscatParameters", "Cannot open aeff.tab file 28\n");
				return -1;
			}
			dgets(bbb, 255, file28);					// skip one header line

			nrad = 0;
			while(1)
			{
				char *ia = dgets(bbb, 255, file28);
				if (!ia)
					break;
				nrad++;
				aefini = reader.ScanReal(realFormat);
				aeffa[nrad] = aefini;
			}
			fclose(file28);
			sprintf(bbb, "%3d' eff. radii from %7.4lf to %7.4lf", nrad, aeffa[0], aeffa[nrad-1]);
			Wrimsg(ReaparLabel, bbb);
		}
		break;

	default:
		Wrimsg(ReaparLabel, " CDIVID for radii must be LIN,LOG,INV or TAB");
		return reader.GetLineCounter();
		break;
	}	
//
// Define incident polarizations (in Lab Frame)
// It is assumed that incident radiation is along x axis in Lab Frame
	theNodes[rootNodeLevel + 1] = XmlFindChild(theNodes[rootNodeLevel], (const xmlChar *)"IncidentPolarization");
	if (!theNodes[rootNodeLevel + 1])
	{
		Wrimsg(ReaparLabel, "Cannot find root->IncidentPolarization node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	Wrimsg(ReaparLabel, "**** Define Incident Polarizations ****\n");
//
// Read Complex polarization vector CXE01_LF=e01 (normalize if necessary).
	theNodes[rootNodeLevel + 2] = XmlFindChild(theNodes[rootNodeLevel + 1], (const xmlChar *)"PolarizationState");
	if (!theNodes[rootNodeLevel + 2])
	{
		Wrimsg(ReaparLabel, "Cannot find root->IncidentPolarization->PolarizationState node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	theNodes[rootNodeLevel + 3] = XmlFindChild(theNodes[rootNodeLevel + 2], (const xmlChar *)"X");	
	if (!theNodes[rootNodeLevel + 3])
	{
		Wrimsg(ReaparLabel, "Cannot find root->IncidentPolarization->PolarizationState->X node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	cxe01_lf.data[0].re = XmlGetDoubleProperty(theNodes[rootNodeLevel + 3], (const xmlChar *)"Re");
	cxe01_lf.data[0].im = XmlGetDoubleProperty(theNodes[rootNodeLevel + 3], (const xmlChar *)"Im");
	theNodes[rootNodeLevel + 3] = XmlFindChild(theNodes[rootNodeLevel + 2], (const xmlChar *)"Y");	
	if (!theNodes[rootNodeLevel + 3])
	{
		Wrimsg(ReaparLabel, "Cannot find root->IncidentPolarization->PolarizationState->Y node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	cxe01_lf.data[1].re = XmlGetDoubleProperty(theNodes[rootNodeLevel + 3], (const xmlChar *)"Re");
	cxe01_lf.data[1].im = XmlGetDoubleProperty(theNodes[rootNodeLevel + 3], (const xmlChar *)"Im");
	theNodes[rootNodeLevel + 3] = XmlFindChild(theNodes[rootNodeLevel + 2], (const xmlChar *)"Z");	
	if (!theNodes[rootNodeLevel + 3])
	{
		Wrimsg(ReaparLabel, "Cannot find root->IncidentPolarization->PolarizationState->Z node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	cxe01_lf.data[2].re = XmlGetDoubleProperty(theNodes[rootNodeLevel + 3], (const xmlChar *)"Re");
	cxe01_lf.data[2].im = XmlGetDoubleProperty(theNodes[rootNodeLevel + 3], (const xmlChar *)"Im");
//
	real e1 = cxe01_lf.data[0].mod();
	if (e1 != (real)0.)
	{
		Wrimsg(ReaparLabel, " cxe01_lf(1) must be zero!");
		xmlFreeDoc(doc);
		return -1;
	}
// Normalize:
	e1 = Sqrt(cxe01_lf.data[1].modSquared() + cxe01_lf.data[2].modSquared());
	cxe01_lf.data[1] /= e1;
	cxe01_lf.data[2] /= e1;
// Construct orthogonal normalized polarization vector CXE02=e02 using xhat cross e01* = e02
	cxe02_lf.data[0].clear();
	cxe02_lf.data[1] = -(cxe01_lf.data[2].conjg());
	cxe02_lf.data[2] =   cxe01_lf.data[1].conjg();
//	
// IORTH	= 1 to calculate only for single polarization
//			= 2 to also calculate for orthogonal polarization
	iorth = XmlGetIntProperty(theNodes[rootNodeLevel + 1], (const xmlChar *)"Iorth");
	sprintf(bbb, "IORTH=%2d\n", iorth);
	Wrimsg(ReaparLabel, bbb);
//
// Specify whether or not to write ".sca" files  skip line:
	Wrimsg(ReaparLabel, "**** Specify which output files to write ****\n");
//
// IWRKSC	= 0 to NOT write ".sca" file for each target orientation
//			= 1 to write ".sca" file for each target orientation
	theNodes[rootNodeLevel + 1] = XmlFindChild(theNodes[rootNodeLevel], (const xmlChar *)"Iwrksc");
	if (!theNodes[rootNodeLevel + 1])
	{
		Wrimsg(ReaparLabel, "Cannot find root->Iwrksc node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	iwrksc = XmlGetIntProperty(	theNodes[rootNodeLevel + 1], (const xmlChar *)"Value");
	sprintf(bbb, "IWRKSC=%2d\n", iwrksc);
	Wrimsg(ReaparLabel, bbb);
//
// Specify whether or not to write ".pol" files
// 12.01.23 (BTD) * eliminate this -- no longer needed since we now include support for fast nearfield calculations, and program call readE to read the nearfield output files.
//                * now initialize IWRPOL=0, although this can then be over-ridden to set IWRPOL=1 when NRFLD=1
//
// IWRPOL	= 0 to NOT write ".pol1" and ".pol2" file for each target orientation (beta,theta)
//			= 1 to write ".pol1" and ".pol2" file for each target orientation (beta,theta)
//
//	strcpy(cwhere, "error reading IWRPOL in ddscat.par");
//	fgets(Buffer, 255, iopar);
//	sscanf(Buffer, "%d", &iwrpol);
//	sprintf(cmsgnm, "IWRPOL=%2d", iwrpol);
//	Wrimsg(ReaparLabel, cmsgnm);
	iwrpol = 0;
//
// In the event that user set NRFLD=1 but IWRPOL=0, set IWRPOL to 1
	if((nrfld == NearfieldMethodPrepare) && (iwrpol <= 0))
	{
		iwrpol = 1;
		Wrimsg(ReaparLabel, "set iwrpol=1 because nrfld=1");
	}
//
// Read information determining target rotations
	real a, b;
	int nb, nt, np;	
	Wrimsg(ReaparLabel, "**** Prescribe Target Rotations ****\n");
	theNodes[rootNodeLevel + 1] = XmlFindChild(theNodes[rootNodeLevel], (const xmlChar *)"PrescribeTargetRotations");
	if (!theNodes[rootNodeLevel + 1])
	{
		Wrimsg(ReaparLabel, "Cannot find root->PrescribeTargetRotations node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	theNodes[rootNodeLevel + 2] = XmlFindChild(theNodes[rootNodeLevel + 1], (const xmlChar *)"Beta");
	if (!theNodes[rootNodeLevel + 2])
	{
		Wrimsg(ReaparLabel, "Cannot find root->PrescribeTargetRotations->Beta node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	a = XmlGetDoubleProperty(theNodes[rootNodeLevel + 2], (const xmlChar *)"Min");
	b = XmlGetDoubleProperty(theNodes[rootNodeLevel + 2], (const xmlChar *)"Max");
	nb = XmlGetIntProperty(theNodes[rootNodeLevel + 2], (const xmlChar *)"Number");
	oridata->SetBetas(a, b, nb);
	theNodes[rootNodeLevel + 2] = XmlFindChild(theNodes[rootNodeLevel + 1], (const xmlChar *)"Theta");
	if (!theNodes[rootNodeLevel + 2])
	{
		Wrimsg(ReaparLabel, "Cannot find root->PrescribeTargetRotations->Theta node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	a = XmlGetDoubleProperty(theNodes[rootNodeLevel + 2], (const xmlChar *)"Min");
	b = XmlGetDoubleProperty(theNodes[rootNodeLevel + 2], (const xmlChar *)"Max");
	nt = XmlGetIntProperty(theNodes[rootNodeLevel + 2], (const xmlChar *)"Number");
	oridata->SetThetas(a, b, nt);
	theNodes[rootNodeLevel + 2] = XmlFindChild(theNodes[rootNodeLevel + 1], (const xmlChar *)"Phi");
	if (!theNodes[rootNodeLevel + 2])
	{
		Wrimsg(ReaparLabel, "Cannot find root->PrescribeTargetRotations->Phi node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	a = XmlGetDoubleProperty(theNodes[rootNodeLevel + 2], (const xmlChar *)"Min");
	b = XmlGetDoubleProperty(theNodes[rootNodeLevel + 2], (const xmlChar *)"Max");
	np = XmlGetIntProperty(theNodes[rootNodeLevel + 2], (const xmlChar *)"Number");
	oridata->SetPhis(a, b, np);
	oridata->SetDimensions(nb, nt, np);
	FileNamer::GetInstance()->Init(nb * nt * np);
//
// check that user has not requested more than 1000 orientations and IWRKSC=1
	if (iwrksc > 0 && (oridata->GetNbeta() * oridata->GetNtheta() * oridata->GetNphi()) > 1000)
	{
		Wrimsg(ReaparLabel, "error: if iwrksc=1, nbeta*ntheta*nphi must be .le. 1000");
		xmlFreeDoc(doc);
		return -1;
	}
//	
	Wrimsg(ReaparLabel, "**** Specify first IWAV, IRAD, IORI (normally 0 0 0) ****\n");	
	theNodes[rootNodeLevel + 1] = XmlFindChild(theNodes[rootNodeLevel], (const xmlChar *)"SpecifyFirst");	
	if (!theNodes[rootNodeLevel + 1])
	{
		Wrimsg(ReaparLabel, "Cannot find root->SpecifyFirst node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	iwav0 = XmlGetIntProperty(theNodes[rootNodeLevel + 1], (const xmlChar *)"Iwav");
	irad0 = XmlGetIntProperty(theNodes[rootNodeLevel + 1], (const xmlChar *)"Irad");
	iori0 = XmlGetIntProperty(theNodes[rootNodeLevel + 1], (const xmlChar *)"Iori");
//
// Check that IWAV0,IRAD0,IORI0 are OK:
	if (iwav0+1 > nwav)
	{
		Wrimsg(ReaparLabel, "iwav0+1 > nwav");
		xmlFreeDoc(doc);
		return -1;
	}
	if (irad0+1 > nrad) 
	{
		Wrimsg(ReaparLabel, "irad0+1 > nrad");
		xmlFreeDoc(doc);
		return -1;
	}
	if (iori0+1 > oridata->GetNbeta() * oridata->GetNtheta() * oridata->GetNphi()) 
	{
		Wrimsg(ReaparLabel, "iori0+1 > nbeta * ntheta * nphi");
		xmlFreeDoc(doc);
		return -1;
	}
//
// If NPHI>1, then set IORTH=2 regardless of value input.
	if ((iorth == 1) && (oridata->GetNphi() > 1))
	{
		iorth = 2;
		Wrimsg(ReaparLabel, "set iorth=2 since nphi>1 ");
	}
	switch(iorth)
	{
	case 1:
		Wrimsg(ReaparLabel, "Calculate only for single polarization ");
		break;

	case 2:
		Wrimsg(ReaparLabel, "Do orthogonal polarization for each case ");
		break;

	default:
		Wrimsg(ReaparLabel, " WRONG VALUE OF IORTH ");
		xmlFreeDoc(doc);
		return -1;
		break;
	}
	sprintf(bbb, "%7.2lf%7.2lf  Range of BETA values ; NBETA =%d", oridata->Betami(), oridata->Betamx(), oridata->GetNbeta());
	Wrimsg(ReaparLabel, bbb);
	sprintf(bbb, "%7.2lf%7.2lf  Range of THETA values; NTHETA=%d", oridata->Thetmi(), oridata->Thetmx(), oridata->GetNtheta());
    Wrimsg(ReaparLabel, bbb);
	sprintf(bbb, "%7.2lf%7.2lf  Range of PHI values ;   NPHI =%d", oridata->Phimin(), oridata->Phimax(), oridata->GetNphi());
	Wrimsg(ReaparLabel, bbb);
//
// Convert from degrees to radians
	oridata->ToRadians();
//
// Specify elements of scattering matrix to be printed out
	Wrimsg(ReaparLabel, "**** Select Elements of S_ij Matrix to Print ****\n");
	theNodes[rootNodeLevel + 1] = XmlFindChild(theNodes[rootNodeLevel], (const xmlChar *)"S_ijMatrix");	
	if (!theNodes[rootNodeLevel + 1])
	{
		Wrimsg(ReaparLabel, "Cannot find root->S_ijMatrix node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	nsmelts = XmlGetIntProperty(theNodes[rootNodeLevel + 1], (const xmlChar *)"Nsmelts");
//
// NSMELTS = number of elements of scattering matrix element to be printed out. 
// Must be no greater than 9 if NSMELTS is zero, then 6 "default" elements are output (ChB, no any more)
	if (nsmelts <= 0)
	{
		nsmelts = 6;
		smind1 = new int[nsmelts];
		smind2 = new int[nsmelts];
		smind1[0] = 11;
		smind1[1] = 21;
		smind1[2] = 31;
		smind1[3] = 41;
		smind1[4] = 12;
		smind1[5] = 13;
		reader.Read();
	}
	else
	{
		smind1 = new int[nsmelts];
		smind2 = new int[nsmelts];
		theNodes[rootNodeLevel + 2] = XmlFindChild(theNodes[rootNodeLevel + 1], (const xmlChar *)"ij");	
		if (!theNodes[rootNodeLevel + 2])
		{
			Wrimsg(ReaparLabel, "Cannot find root->S_ijMatrix->ij node.\n");
			xmlFreeDoc(doc);
			return -1;
		}
		strcpy(bbb, (char *)xmlGetProp(theNodes[rootNodeLevel + 2], (const xmlChar *)"Value"));
		char *ia = strtok(bbb, " \t");
		for(j=0; j<nsmelts; ++j)
		{
			sscanf(ia, "%d", &smind1[j]);
			ia = strtok(NULL, " \t");
		}
	}
	for(j=0; j<nsmelts; ++j)
	{
		smind2[j] = smind1[j] % 10;
		smind1[j] /= 10;
	}	
//
// Specify scattering directions to be calculated
// Two options:
// CMDFRM	= 'LFRAME': specify scattering directions n in Lab Frame (where incident beam is in x-direction)
//			  THETAN,PHIN = Direction of n from n0 in calculation of the scattering matrix; cos(THETAN) is n0 \dot n
//			  and PHIN is azimuthal angle of n from Lab xy pla 
// Note: THETA1, THETA2, and PHI for each scatterin plane are entered in degrees, and immediately converted to radians.
// Arbitrary number of scattering planes may be considered.
//
// CMDFRM	= 'TFRAME': specify scattering directions in Target Frame, defined by target axes a1,a2,a3
// If JPBC = 0:
//			THETAN,PHIN = Direction of n relative to a1,a2,a3: 
//			THETAN = angle between n and a1
//			PHIN   = angle between a1-n plane and a1-a2 plane
//    JPBC = 1:
//			THETAN = Diffraction order along y_TF axis
//			PHIN   = azimuthal angle around y_TF
//	  JPBC = 2:
//			THETAN = Diffraction order along z_TF axis
//			PHIN   = azimuthal angle around z_TF
//	  JPBC = 3:
//			THETAN = Diffraction order along y_TF
//			PHIN   = Diffraction order along z_TF
//where we first run through transmitted waves 1 -> NSCAT/2 
//      and then run through reflected waves NSCAT/2+1 -> NSCAT
// 
// Three cases:
//	JPBC = 0: single isolated target.
//				specify
//				phi for scattering plane
//				thetamin, thetamax, dtheta for scattering plane
//	JPBC = 1,2: periodic in one dimension.
//				specify
//				diffraction order in direction of target periodicity
//				phimin, phimax, dphi for scattering cone
//	JPBC=3:		periodic in two dimensions
//				specify
//				difraction order in y direction and order in z direction
	theNodes[rootNodeLevel + 1] = XmlFindChild(theNodes[rootNodeLevel], (const xmlChar *)"ScatteredDirections");	
	if (!theNodes[rootNodeLevel + 2])
	{
		Wrimsg(ReaparLabel, "Cannot find root->ScatteredDirections node.\n");
		xmlFreeDoc(doc);
		return -1;
	}
	strcpy(bbb, (char *)xmlGetProp(theNodes[rootNodeLevel + 1], (const xmlChar *)"Cmdfrm"));
	bTemp = ProcessCmdfrm(bbb, ReaparLabel);
	if (!bTemp)
	{
		xmlFreeDoc(doc);
		return -1;	
	}
//	
	int nplanes = XmlGetIntProperty(theNodes[rootNodeLevel + 1], (const xmlChar *)"Nplanes");
	switch(TargetManager::GetInstance()->Jpbc())
	{
	case PeriodicNo:
		sprintf(bbb, "%4d = number of scattering planes", nplanes);
		break;

	case PeriodicY:
	case PeriodicZ:
		sprintf(bbb, "%4d = number of scattering cones", nplanes);
		break;

	case PeriodicBoth:
		sprintf(bbb, "%4d = number of diffraction orders for transmission", nplanes);
		break;

	default:
		break;
	}
	Wrimsg(ReaparLabel, bbb);
//	
	nscat = 0;
	if (nplanes > 0)
	{
		if (TargetManager::GetInstance()->Jpbc() != PeriodicBoth)
		{
			int nsca0;
			theNodes[rootNodeLevel + 2] = theNodes[rootNodeLevel + 1]->children;
			for(int jplane=0; jplane<nplanes; ++jplane)
			{
				while(xmlStrcmp(theNodes[rootNodeLevel + 2]->name, (const xmlChar *)"Plane"))
					theNodes[rootNodeLevel + 2] = theNodes[rootNodeLevel + 2]->next;
				real phi1 = XmlGetDoubleProperty(theNodes[rootNodeLevel + 2], (const xmlChar *)"Phi");
				real theta1 = XmlGetDoubleProperty(theNodes[rootNodeLevel + 2], (const xmlChar *)"MinTheta");
				real theta2 = XmlGetDoubleProperty(theNodes[rootNodeLevel + 2], (const xmlChar *)"MaxTheta");
				real dtheta = XmlGetDoubleProperty(theNodes[rootNodeLevel + 2], (const xmlChar *)"Dtheta");
#pragma message ("next variable pp needs TODO")
// TODO				int pp = XmlGetIntProperty(theNodes[rootNodeLevel + 2], (const xmlChar *)"N");
				if (TargetManager::GetInstance()->Jpbc() == PeriodicNo)
					sprintf(bbb, "%7.1lf%7.1lf%7.1f = phi, theta_min, theta_max for scatt. plane%4d", phi1, theta1, theta2, jplane);
				else
					sprintf(bbb, "%7.1lf%7.1lf%7.1f = order, zeta_min, zeta_max for scattering cone", phi1, theta1, theta2);
				Wrimsg(ReaparLabel, bbb);
				if (TargetManager::GetInstance()->Jpbc() == PeriodicNo)
				{
					if((Fabs(theta1-theta2) > (real)0.) && (dtheta == (real)0.))
					{
						Wrimsg(ReaparLabel, "DTHETA=0 in ddscat.par!");
						xmlFreeDoc(doc);
						return -1;
					}
				}
// ! Convert to radians
				phi1 /= Degrad;
				theta1 /= Degrad;
				theta2 /= Degrad;
				dtheta /= Degrad;
// ! Allow for possibility that user did not enter correct sign for
// ! DTHETA (if I did it, others will too...)
				if (theta2 < theta1)
					dtheta = -Fabs(dtheta);
// ! compute theta values for this scattering plane/cone
				nsca0 = 1;
				real delta = theta2 - theta1;
				if (delta)
					nsca0 = 1 + nint_(delta/dtheta);
// ChB: ! realloc to have sufficient space
				thetan = (real *)realloc(thetan, (nscat + nsca0)*sizeof(real));
				phin = (real *)realloc(phin, (nscat + nsca0)*sizeof(real));
				thetan[nscat] = theta1;
				if(nsca0 > 2)
				{
					for(j=1; j<nsca0-1; ++j)
					{
						thetan[nscat + j] = theta1 + j*dtheta;
					}
				}
				thetan[nscat + nsca0-1] = theta2;
				for(j=0; j<nsca0; ++j)
				{
					phin[nscat + j] = phi1;
				}
				nscat += nsca0;
				sprintf(bbb, "%4d = number of scattering angles in this scattering %s", nsca0, ((TargetManager::GetInstance()->Jpbc() == PeriodicNo) ? "plane" : "cone"));
				Wrimsg(ReaparLabel, bbb);
			}
		}
		else
		{
			sprintf(curFormat, "%s%s", realFormat, realFormat);
			phin = (real *)malloc(2 * nplanes * sizeof(real));
			thetan = (real *)malloc(2 * nplanes * sizeof(real));
			for(int jplane=0; jplane<nplanes; ++jplane)
			{
				while(xmlStrcmp(theNodes[rootNodeLevel + 2]->name, (const xmlChar *)"Order"))
					theNodes[rootNodeLevel + 2] = theNodes[rootNodeLevel + 2]->next;
				real aa = XmlGetDoubleProperty(theNodes[rootNodeLevel + 2], (const xmlChar *)"Phi");
				real bb = XmlGetDoubleProperty(theNodes[rootNodeLevel + 2], (const xmlChar *)"Theta");
				phin[nscat] = aa;
				thetan[nscat] = bb;
				++nscat;
			}
			for(j=0; j<nscat; ++j)
			{
				phin[j + nscat] = phin[j];
				thetan[j + nscat] = thetan[j];
			}
			nscat = nscat + nscat;
		}
	}
	xmlFreeDoc(doc);

	return 0;
}

void OriData::SetDimensions(int nb, int nt, int np)
{
	nbeta = nb;
	ntheta = nt;
	nphi = np;
//
	if (nbeta > 0)
	{
		beta = (real *)malloc(nbeta * sizeof(real));
		wgtb = (real *)malloc(nbeta * sizeof(real));
	}
	if (ntheta > 0)
	{
		theta = (real *)malloc(ntheta * sizeof(real));
		if (nphi > 0)
		{
			wgta.Dimension(ntheta, nphi);
			phi = (real *)malloc(nphi * sizeof(real));
		}
	}
}

void OriData::Orient()
// Orient(real *beta, real *theta, real *phi, Array2F<real> &wgta, real *wgtb)
{
/* **
 Given:
        BETAMI=minimum value of beta (radians)
        BETAMX=maximum value of beta (radians)
        THETMI=minimum value of theta (radians)
        THETMX=maximum value of theta (radians)
        PHIMIN=minimum value of phi (radians)
        PHIMAX=maximum value of phi (radians)
        MXBETA,MXTHET,MXPHI=parameters for dimensioning of arrays BETA,THETA,PHI
        NBETA=desired number of values of beta
        NTHETA=desired number of values of theta
        NPHI=desired number of values of PHI
 Returns:
        BETA(1-NBETA)=beta values (radians)
        THETA(1-NTHETA)=theta values (radians)
        PHI(1-NPHI)=phi values (radians)
        WGTA(1-NTHETA,1-NPHI)=weighting of each orientation of target axis a1 in Lab Frame (sum of weights = 1)
        WGTB(1-NBETA)=weighting of each rotation of target around a1 (sum of weights = 1)

 Note: it is assumed that target orientation weight function can be factored into WGTA*WGTB -- i.e., that rotation around a1 are decoupled from orientation of a1.

 Purpose: to generate list of desired target orientations

 Definitions: beta=angle of rotation of target around target axis a1 (beta=0 is defined to be such that a2 lies in a1-x plane, with a2 in direction of increasing theta)
              theta=angle between target axis a1 and x axis (0.le.theta.le.pi)
              phi=angle of rotation of axis a1 around axis (phi=0 is defined to be such that a1 lies in xy plane)

 Present version assumes:
        beta to be uniformly distributed between BETAMI and BETAMX
        cos(theta) to be uniformly distributed between cos(THETMI) and cos(THETMX)
        phi to be uniformly distributed between PHIMIN and PHIMAX

 Values assigned to BETA and PHI are midpoints of uniform intervals in beta and phi
 Note that if NBETA=1, BETA is set to midpoint of range in phi
 Likewise, if NPHI=1, PHI is set to midpoint of range in phi

 If NTHETA=1: specify a single value of THETA cos(theta)=0.5*(cos(thetami)+cos(thetamx))

 If NTHETA>1:
    If NTHETA is even:
       divide range of cos(theta) into NTHETA equal intervals
       set THETA to midpoints of these intervals
       give intervals equal weights
    If NTHETA is odd:
       divide range of cos(theta) into NTHETA-1 equal intervals
       set THETA to endpoints of these intervals; first value of
       THETA is THETAMI and last value of THETA is THETAMX

 Note: values chosen for beta, phi are at midpoints of uniform intervals

 B.T.Draine, Princeton Univ. Obs., 89.11.20
 History:
 92.04.01 (BTD) Modify to allow Simpson's rule or trapezoidal rule
                for integration over cos(theta).
 99.03.01 (BTD) Modified to properly assign THETA values when NPHI>1
                Code up to this time was not consistent with
                description in UserGuide.  This problem was brought
                to my attention by Miroslav Kocifaj
                (Miroslav.Kocifaj@swh.sk).
 99.03.05 (BTD) Corrected errors in assignment of THETA values.
                Thanks to Miroslav Kocifaj for detecting these.
 99.04.30 (BTD) Corrected errors in assignment of BETA and
                PHI values when BETAMI.ne.0 or PHIMIN.ne.0
                These errors (inadvertent omission of BETAMI and
                PHIMIN from expressions assigning BETA(J) and PHI(J))
                were introduced in 99.03.05 revision.
                Thanks to Henriette Lemke for calling attention to these errors.
 Copyright (C) 1993,1999 B.T. Draine and P.J. Flatau
 This code is covered by the GNU General Public License.
** */

	const real half_ = (real)0.5;
	const real onex_ = (real)1.;
//
// Assign values to BETA:
	int i, j;
//
	real delta = (betamx - betami) / (real)nbeta;
	for(j=0; j<nbeta; ++j)
	{
		beta[j] = betami + delta * (j + half_);
	}
//
// Assign values to THETA:
	if (!(ntheta % 2))
	{
		delta = (Cos(thetmx) - Cos(thetmi)) / (real)ntheta;										// THETA is even:
		for(j=0; j<ntheta; ++j)
		{
			theta[j] = Acos(Cos(thetmi) + delta * (j + half_));
		}
	}
	else
	{
		if (ntheta == 1)																		// THETA is odd:
		{
			theta[0] = Acos(half_ * (Cos(thetmi) + Cos(thetmx)));
		}
		else
		{
			delta = (Cos(thetmx) - Cos(thetmi)) / (real)(ntheta - 1);
			theta[0] = thetmi;
			theta[ntheta - 1] = thetmx;
			for(j=1; j<ntheta-1; ++j)
			{
				theta[j] = Acos(Cos(thetmi) + delta * j);
			}
		}
	}
//
// Assign values to PHI:
	real wg;
	delta = (phimax - phimin) / (real)nphi;
	for(j=0; j<nphi; ++j)
	{
		phi[j] = phimin + delta * (j + half_);
	}
//
// Specify weight factors WGTA, WGTB
//    (Note: weight function WGTA = 4*pi*P/(NTHETA*NPHI), where
//     P=(probability/solid angle) of orientation in direction THETA,PHI .
//     The orientational averaging program automatically samples uniformly
//     in cos(theta) to allow for d(solid angle)=sin(theta)d(theta)d(phi).
//
//    Present version assumes random orientations.
//
//    When NTHETA is >1 and even, we use trapezoidal integration
//    When NTHETA is >1 and odd, we use Simpson's rule integration.
//
	if ((ntheta == 1) || !(ntheta % 2))												// either 1 or even number of theta values: midpoints of intervals
	{
		wg = onex_ / (real)(ntheta * nphi);
		for(i=0; i<ntheta; ++i)
		{
			for(j=0; j<nphi; ++j)
			{
				wgta.Value(i, j) = wg;
			}
		}
	}
	else																			// odd number >1 of theta values: use Simpson's rule weighting
	{
		wg = onex_ / (real)(3 * (ntheta - 1) * nphi);
		for(j=0; j<nphi; ++j)
		{
			wgta.Value(0, j) = wg;
			wgta.Value(ntheta-1, j) = wg;
		}
        wg = (real)4. / (real)(3 * (ntheta - 1) * nphi);
        for(i=1; i<ntheta; i+=2)
		{
			for(j=0; j<nphi; ++j)
			{
				wgta.Value(i, j) = wg;
			}
		}
        if (ntheta >= 5)
		{
			wg = (real)2. / (real)(3 * (ntheta - 1) * nphi);
			for(i=2; i<ntheta-1; i+=2)
			{
				for(j=0; j<nphi; ++j)
				{
					wgta.Value(i, j) = wg;
				}
			}
		}
	}
//
// Assign weights for rotation angle beta:
	wg = onex_ / (real)nbeta;
	for(j=0; j<nbeta; ++j)
	{
		wgtb[j] = wg;
	}
}

bool DDscatParameters::Divide(CdividMethod wDivid, real x1, real x2, int nx, real *ary)
{
/* **
!***********************************************************************
! Given:
!       CDIVID='LIN','INV', or 'LOG'
!       X1 = lower limit to interval
!       X2 = upper limit to interval
!       NX = number of elements
!       MXARY=dimensioning information for array ARY

! Returns:
!       ARY(1-NX)=vector of points with
!                 ARY(1)=X1
!                 ARY(NX)=X2
!                 ARY(J) values spaced either
!                        linearly (if CDIVID.EQ.'LIN')
!                        linearly in 1/X (if CDIVID.EQ.'INV')
!                        logarithmically (if CDIVID.EQ.'LOG')

! Copyright (C) 1993, B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
** */
	const char *DivideLabel = "Divide";
	const real onex_ = (real)1.;

//	if (nx > mxary)
//		Errmsg("Fatal", DivideLabel, "  NX .GT. MXARY ");
//!      IF(X1.GT.X2)THEN
//!         CALL ERRMSG('FATAL','DIVIDE','  X1 .GT. X2 ')
//!      ENDIF
	if (nx < 1)
	{
		Wrimsg(DivideLabel, "  NX .LT. 1 ");
		return false;
	}
	if (nx == 1)
	{
		Wrimsg(DivideLabel, " Only one element initialized ");
		ary[0] = x1;
	}
	else
	{
		int i;
		real delta = (real)0.;
		switch(wDivid)
		{
		case CdividLin:
			delta = (x2 - x1) / (real)(nx-1);
			for(i=0; i<nx; ++i)
				ary[i] = x1 + delta * (real)(i-1);
			break;

		case CdividInv:
			delta = (onex_ / x2 - onex_ / x1) / (real)(nx-1);
			for(i=0; i<nx; ++i)
				ary[i] = onex_ / (onex_ / x1 + delta * (real)(i-1));
			break;

		case CdividLog:
			delta = Log(x2 / x1) / (real)(nx-1);
			for(i=0; i<nx; ++i)
				ary[i] = Exp(Log(x1) + delta * (real)(i-1));
			break;

		default:
			break;
		}
	}
	return true;
}

char *DDscatParameters::dgets(char *Buffer, int num, FILE *file)
{
	char *ia = fgets(Buffer, num, file);
	while((Buffer[strlen(Buffer) - 1] == 0x0d) || (Buffer[strlen(Buffer) - 1] == 0x0a))
	{
		Buffer[strlen(Buffer) - 1] = '\0';
	}
	return ia;
}

bool DDscatParameters::ProcessCmtorq(char *Buffer, const char *Label)
{
	cmdtrq = TorqEnumerator(Buffer);
	switch(cmdtrq)
	{
	case TorqMethod_DOTORQ:
		strcat(Buffer, " - compute torques ");
		break;

	case TorqMethod_NOTORQ:
		strcat(Buffer, " - do not compute torques ");
		break;

	default:
		break;
	}
	if (cmdtrq == TorqMethod_End)
	{
		Wrimsg(Label, " wrong definition of cmdtrq");
		return false;
	}
	else
	{
		Wrimsg(Label, reader.GetBuffer());
		return true;
	}
}

bool DDscatParameters::ProcessCmdsol(char *Buffer, const char *Label)
{
	cmdsol = SolEnumerator(Buffer);
	if (cmdsol == SolMethod_End)
	{
		Wrimsg(Label, " wrong definition of CMDSOL");
		return false;
	}
	else
	{
		strcat(reader.GetBuffer(), " - CCG Method  ");
		Wrimsg(Label, reader.GetBuffer());
		return true;
	}
}

bool DDscatParameters::ProcessCmdFFT(char *Buffer, const char *Label)
{
	cmdfft = FftEnumerator(Buffer);
	if (cmdfft == FftMethod_End)
	{
		Wrimsg(Label, "DDSCAT 7.2 only supports FFT options FFTW21, GPFAFT, and FFTMKL");
		Wrimsg(Label, " wrong definition of cmdfft");
		return false;
	}
	else
	{
		fftEngine = AbstractFftEngine::GetEngine(cmdfft);
		fftEngine->SayHello(reader.GetBuffer());
		Wrimsg(Label, reader.GetBuffer());
		return true;
	}
}

bool DDscatParameters::ProcessCalpha(char *Buffer, const char *Label)
{
	calpha = AlphaEnumerator(Buffer);
	switch(calpha)
	{
	case AlphaMethod_LATTDR:
		Wrimsg(Label, "LATTDR - Draine & Goodman (1993) LDR for alpha");
		Common0::GetInstance()->Get((unsigned int)0)->Idipint() = 0;
		Common0::GetInstance()->Get((unsigned int)1)->Idipint() = 0;
		break;

	case AlphaMethod_GKDLDR:
		Wrimsg(Label, "GKDLDR - Gutkowicz-Krusin & Draine (2004) LDR for alpha");
		Common0::GetInstance()->Get((unsigned int)0)->Idipint() = 0;
		Common0::GetInstance()->Get((unsigned int)1)->Idipint() = 0;
		break;

	case AlphaMethod_FLTRCD:
		Wrimsg(Label, "FLTRCD - Filtered Coupled Dipole (Piller & Martin 1998)");
		Common0::GetInstance()->Get((unsigned int)0)->Idipint() = 1;
		Common0::GetInstance()->Get((unsigned int)1)->Idipint() = 1;
		break;

	default:
		Wrimsg(Label, " wrong definition of CALPHA");
		return false;
	}
	return true;
}

bool DDscatParameters::ProcessCbinflag(char *Buffer, const char *Label)
{
	cbinflag = BinflagEnumerator(Buffer);
	if (cbinflag == BinflagMethod_End)
	{
		Wrimsg(Label, " Wrong definition of CBINFLAG");
		return false;
	}
	else
	{
		strcat(reader.GetBuffer(), " - Unformatted binary dump option");
		Wrimsg(Label, reader.GetBuffer());
		return true;
	}
}

bool DDscatParameters::ProcessTarget(char *Buffer, const char *Label)
{
	cshape = TargetManager::GetInstance()->TargetTypeFromString(Buffer);
	strcpy(reader.GetBuffer(), Buffer);
	if (cshape != TargetType_End)
	{
		strcat(reader.GetBuffer(), " - Shape definition ");
		Wrimsg(Label, reader.GetBuffer());
	}
	else
	{
		strcat(reader.GetBuffer(), " - Unrecognized shape directive");
		Wrimsg(Label, reader.GetBuffer());
		return false;
	}
//
// Read shape parameters
	if (cshape != TargetType_End)
	{
		reader.Read();
		TargetManager::GetInstance()->PreloadTarget(targetString, reader.GetBuffer());		// jpbc is calculated here too
	}
	return true;
}

bool DDscatParameters::ProcessCmdfrm(char *Buffer, const char *Label)
{
	cmdfrm = FrameEnumerator(Buffer);
	if (cmdfrm == FrameCode_End)
	{
		Wrimsg(Label, " Error reading ddscat.par file");
		return false;
	}
	if (TargetManager::GetInstance()->Jpbc() != PeriodicNo)
	{
		if (cmdfrm == FrameCode_LFRAME)
		{
			Wrimsg(Label, " cannot use LFRAME when JPBC != PeriodicNo");
			return false;
		}
	}
	if (TargetManager::GetInstance()->Jpbc() == PeriodicNo)
	{
		char bbb[256];
		sprintf(bbb, "cmdfrm=%s : scattering directions given in %s frame", FrameEnumerator(cmdfrm), (cmdfrm == FrameCode_LFRAME) ? "lab" : "target");
		Wrimsg(Label, bbb);
	}
	return true;
}

int DDscatParameters::LoadShparFromXml(xmlNodePtr nod1)
{
	int num = 0;
	while(1)
	{
		nod1 = nod1->next;
		if (!xmlStrcmp(nod1->name, (const xmlChar *) "Shpar")) 
		{
			int pos = XmlGetIntProperty(nod1, (const xmlChar *)"Pos");
			real value = XmlGetDoubleProperty(nod1, (const xmlChar *)"Value");
			shpar[pos-1] = value;
			++num;
		}
	}
	return num;
}
