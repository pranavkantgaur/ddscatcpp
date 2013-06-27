#include "StdAfx.h"

#define OLDVERSION

#include "Definitions.h"
#include "Complex.h"
#include "Vect3.h"
#include "DDscatCommons.h"
#include "DDscatMain.h"
#include "DDscatParameters.h"
#include "TargetManager.h"
#include "DielectricManager.h"
#include "GreenFunctionManager.h"
#include "ArrayF.h"
#include "DipoleData.h"
#include "SumHolder.h"
#include "FileNamer.h"
#include "OutputManager.h"
#ifndef OLDVERSION
	#include "DDscatEngine.h"
#endif

void DDscat(NearfieldMethod &nfrld);

/* **
Each call to subroutine DDSCAT must
1. Supply CFLPAR to subroutine DDSCAT
2. Reset values of AK2OLD,AK3OLD,WOLD between calls to DDSCAT
   to force recalculation of A_ij by subroutine ESELF
   (If not done, earlier A_ij values might be inadvertently reused).
This is accomplished by setting values here, which are communicated
through MODULE DDCOMMON_0
3. Communicate NGRID calculated in ESELF to other routines
** */
#ifndef OLDVERSION
int main(int argc, const char *argv[])
{
	DDscatEngine *theEngine = new DDscatEngine;
	if (argc == 2)
		theEngine->SetParameterFile(argv[1]);
	else
		theEngine->SetParameterFile();

	Common0::GetInstance()->InitReals();
	theEngine->SayHello();
	theEngine->Run();
	theEngine->SayBye();

	if (theEngine->GetNrfld() != NearfieldMethodNo)
	{
		Common0::GetInstance()->InitReals();
		theEngine->SayHello();
		theEngine->Run();
		theEngine->SayBye();
	}

	delete theEngine;
	TargetManager::Kill();

	return 0;
}

#else
int main(int argc, const char *argv[])
{
	const char *Cflpar_default = "ddscat.par";
	char cmsgnm[256];
	NearfieldMethod nrfld;
//
// initialize NEARFIELD to 0 before first call to DDSCAT
	strcpy(Common0::GetInstance()->Cflpar(), Cflpar_default);
	if (argc == 2)
	{
		strcpy(Common0::GetInstance()->Cflpar(), argv[1]);
	}
	sprintf(cmsgnm, " using parameter file=%s", Common0::GetInstance()->Cflpar());
	Wrimsg("DDscat", cmsgnm);

	nrfld = NearfieldMethodNo;
	Common0::GetInstance()->InitReals();

	DDscat(nrfld);

	if(nrfld == NearfieldMethodNo)
	{
		sprintf(cmsgnm, "normal termination, with no nearfield calculation");
		Wrimsg("DDscat", cmsgnm);
	}
	else
	{
		Common0::GetInstance()->InitReals();

		sprintf(cmsgnm, "use calculated polarization to do nearfield calculations");
		Wrimsg("DDscat", cmsgnm);
     
		DDscat(nrfld);
		sprintf(cmsgnm, "normal termination after nearfield calculation");
		Wrimsg("DDscat", cmsgnm);
	}

	Common0::Kill();
	Common1::Kill();
	Common4::Kill();
	Common10::Kill();

    TargetManager::Kill();
	return 0;
}
#endif

#ifdef OLDVERSION
//
// Main function of DDscat application
//
// Original version created by Choliy V., Kyiv Shevchenko University 
// Conforms to fotran version of DDscat 7.2.0 by 
// P.J.Flatau, Colorado State Univ. and B.T.Draine, Princeton Univ. Obs.
//
// History:
// 12.02.27 (ChB): Initial version.
// end history
//
// Copyright (C) 2012 V.Choliy
// Copyright (C) All fortran versions of DDscat, 1993-2012, B.T. Draine and P.J. Flatau
// This code is covered by the GNU General Public License.
//
void DDscat(NearfieldMethod &nrfld)
{
	const real zero_ = (real)0.;
	const real onex_ = (real)1.;
	const real half_ = (real)0.5;
	const char *DDscatLabel = "DDscat";
/* **
           Instructions for Enabling or Disabling MPI
           ==========================================

 To enable MPI use (this requires that you have installed the "full" version of DDSCAT):
                      uncomment INCLUDE 'mpif.h' statement
                      comment out INTEGER MPI_COMM_WORLD statement
                      compile mpi_subs.f

 For non-MPI version (either the "plain" version of DDSCAT, or the "full" version with MPI support disabled:
                      comment out INCLUDE 'mpif.h' statement
                      uncomment INTEGER MPI_COMM_WORLD statement
                      compile mpi_fake.f
** */
#ifdef mpi
	#include "mpif.h"
#else
 	int MPI_Comm_world = 0;
#endif

/* **
    Adjustable Parameters:
    MXNX   = max. extent of target in x direction
    MXNY   =                          y
    MXNZ   =                          z
    MXPBC  = 0 if not planning to use PBC option
           = 1 for best memory use with PBC option
    MXCOMP = max. number of different dielectric functions
    MXTHET = max. number of target rotation angles THETA
    MXBETA = max. number of target rotation angles BETA
    MXPHI  = max. number of target rotation angles PHI
    MXSCA  = max. number of scattered directions
    MXRAD  = max. number of radii
    MXWAV  = max. number of wavelengths
    MXWAVT = max. number of wavelengths in dielectric tables
** */
// !12.04.16 to support OpenMP call
#ifdef openmp
	int nthreads, tid, omp_get_num_threads, omp_get_thread_num;
#endif
//
// *** Set parameter MXPBC
//	const int mxpbc = 1;
//
// Set parameters MXCOMP,MXTHET,MXBETA,MXPHI,MXRAD,MXWAV,MXWAVT,MXSCA			// ChB: now all arrays are allocated dynamically
//
// Derived Parameters:
	int mxnat, mxn3, mxcxsc;
//
// Local variables:
	char cmsgnm[1024];								// ChB: write buffer
	int navg, nrword, posaeff, posak_tf, poscxe0;
	real aeff, ak1, ak3, cosbet, cosphi, costhe, cword, daeff, pyddx, pzddx, rword, sinbet, sinphi, sinthe, wave, xx;

/* **
Target Properties:
AEFF=aeff=effective radius of target (physical units)
PIA2=pi*(aeff/d)**2=pi*(3*NAT/4*pi)**(2/3)

Incident Wave Properties:
AK1=(k*d), where k=2*pi/wave=propagation vector and d=dipole spacing
AK3=(k*d)**3
AKE2=(k*d \dot E0), where E0 is incident polarization vector
WAVE=wavelength (in vacuo) in physical units
XX=k*aeff=2*pi*aeff/wave

Scattering Properties:
G=<cos(theta)> for scattered radiation
QSCGSUM(1-3,1-2)=running (weighted) sum of g(1-3)*Q_sca over target orientations, for incident polarizations 1-2, where
           g(1)=<cos(theta)> for scattered radiation
           g(2)=<sin(theta)cos(phi)> for scattered radiation
           g(3)=<sin(theta)sin(phi)> for scattered radiation
           theta is measured relative to incident direction
           phi is measured relative to x,y plane (Lab Frame)
QPHA(1-2)=(phaseshift cross section)/(pi*aeff**2)
QEXSUM(1-2)=running sum of QEXT over target orientations
QABSUM(1-2)=running sum of QABS over target orientations
QBKSUM(1-2)=running sum of QBKSCA over target orientations
QSCSUM(1-2)=running sum of QSCA over target orientations
QSCG2SUM(1-2)=running sum of cos^2*Q_sca
QTRQSCSUM(1-3,1-2)=running sum of QTRQSC(1-3,1-2) over target orientat
THETND=scattering angle theta in degrees
PHIND=scattering angle phi in degrees
SINPHI=sin(phi)
DEPOLR=depolarization ratio for scattered of initially polarized light into direction theta,phi
PPOL=percentage polarization for scattering of initially unpolarized light into direction theta,phi
** */

	int ibeta, ibeta1, ibeth, ibeth1, idir, iori, iori1, iphi;
    int iphi1, irad, irad1, itask, itheta, itheta1, iwav, iwav1, j, lace, laxi, lclm, lgi, lpi, lqi;
	int lsc0, nbeth, nori, numprocs, itnum[2], itnummx[2];

/* **
MYID     = MPI process identifier that runs from 0-(NUMPROCS-1). 0 is master process.
NUMPROCS = number of MPI parallel processes (including the master)
IERR     = MPI error code, which we ignore
NAT=number of dipoles in target
NAT3=3*NAT
NBETA=number of different beta values for target orientation
NBETH=NBETA*NTHETA=number of target orientations for outer orientation
NTHETA=number of different theta values for target orientation
NPHI=number of different phi values for target orientation
NORI=NBETA*NTHETA*NPHI=total number of target orientations
NSCAT=number of different scattering directions for each orientation
NWAV=number of different wavelengths
NRAD=number of different target radii
NSMELTS=number of scattering matrix elements to print out
INIT=0,1,2 for choice of initial vector |x0> in Complex conjugate grad method
IORTH=1 or 2 to do 1 or 2 incident polarizations
IWRKSC=0 or 1 to suppress or generate "wAArBBkCCC.sca" files for each target orientation
NX=x-length of extended target
NY=y-length of extended target
NZ=z-length of extended target

SMIND1(1-NSMELTS) = index 1 of scattering matrix elements to be printed out (e.g., 1 for element S_{13})
SMIND2(1-NSMELTS) = index 2 of scattering matrix elements to be printed out (e.g., 3 for element S_{13})
** */
//	char cflavg[16], cflpol2[18], cflsca[18], cfle1[16], cfle2[16], cfleb1[17], cfleb2[17];

	Vect3<Complex> cxe01_tf, cxe02_tf;
	FourArray *cxfData = NULL;

/* **
Arrays describing target properties:
CXADIA(1-3*NAT)=diagonal elements of A matrix
CXEPS(1-MXCOMP)=dielectric const. for each of MXCOMP materials
CXREFR(1-3)=Complex refractive index at current wavelength
CXALPH(1-3*NAT)=diagonal elements of dipole polarizability tensor alpha
CXALOF(1-3*NAT)=off-diagonal elements of dipole polarizability tensor

Arrays describing incident radiation:
CXE01_LF(1-3)=incident polarization vector e01 prior to target rotation
CXE02_LF(1-3)=incident polarization vector e02 prior to target rotation
CXE01R_TF(1-3)=rotated incident polarization vector e01 in Target Frame
CXE02R_TF(1-3)=rotated incident polarization vector e02 in Target Frame
CXE0R(1-3)=incident polarization vector (=CXE01R or CXE02R) being used
CXE_TF(1-3*NAT)=incident electric field vector in Target Frame

Array describing polarization of target for one incident wave:
CXXI(1-3*NAT)=polarization vector at each dipole

Arrays describing scattering properties of target:
CXFF11(1-NSCAT)=scattering matrix element f11 for each of NSCAT directions
CXFF12(1-NSCAT)=                          f12
CXFF21(1-NSCAT)=                          f21
CXFF22(1-NSCAT)=                          f22
CX1121(1-NSCAT)=sum of f11* \times f21 over orientations

Scratch arrays:
CXSC(1-MXCXSC)=scratch array used by DDACCG
** */

	TargetManager *manager = TargetManager::GetInstance();
	manager->Ioshp() = true;
	AbstractTarget *currentTarget = manager->GetCurrentTarget();
	DDscatParameters *paramDDscat = DDscatParameters::GetInstance();
	OriData *oridata = paramDDscat->GetOriData();

/* **
IXYZ0(1-NAT0,1-3)=[x-X0(1)]/d,
                  [y-X0(2)]/d,
                  [z-X0(3)]/d for dipoles 1-NAT0 in real target
                  where offset vector X0(1-3) is set by routine TARGET
                  X0(1-3)*DX(1-3) = location of lattice site (0,0,0)
Scratch arrays:
ISCR1,ISCR2,ISCR3,ISCR4 (used by EXTEND)
** */
	
	SumPackage sumPackage;
	Vect3<real> en0_tf;

	real *orderm = NULL;
	real *ordern = NULL;
	Complex *cx1121 = NULL;
	Complex *cx1121_1 = NULL;
	Vect3<real> *aks_tf = NULL;
	Vect3<real> *em1_lf = NULL;
	Vect3<real> *em1_tf = NULL;
	Vect3<real> *em2_tf = NULL;
	Vect3<real> *em2_lf = NULL;
	Vect3<real> *ensc_lf = NULL;
	Vect3<real> *ensc_tf = NULL;

	real *s1111 = NULL;
	real *s1111_1 = NULL;
	real *s2121 = NULL;
	real *s2121_1 = NULL;

	real *sm = NULL;
	real *smori = NULL;
	real *smori_1 = NULL;

/* **
IANISO=0 for isotropic target
       1 for anisotropic target, but with principal axes of
             dielectric tensor || x_TF,y_TF,z_TF
             (i.e., for every dipole, Dielectric Frame = Target Frame)
       2 for general anisotropic target, where orientation of
             Dielectric Frame is specified for every dipole
                                 
Arrays describing target properties:
A1(1-3)=direction of target axis a1 in Target Frame
       =initial direction of target axis a1 in Lab Frame
A2(1-3)=direction of target axis a2 in Target Frame
       =initial direction of target axis a2 in Lab Frame
A3(1-3)=direction of target axis a3=a1xa2 in Target Frame
       =initial direction of target axis a3 in Lab Frame
AEFFA(1-NRAD)=values of effective radius of target (physical units)
SHPAR(1-10)=up to 10 parameters for defining target geometry
X0(1-3)=offset for IXYZ array, such that
        x/d = IXYZ(J,1) - X0(1)
        y/d = IXYZ(J,2) - X0(2)
        z/d = IXYZ(J,3) - X0(3)
Arrays describing orientation of "Dielectric Frame" for anisotropic
material relative to Target Frame, for each dipole location:
BETADF(1-NAT) = angle BETADF (radians)
PHIDF(1-NAT) = angle PHIDF (radians)
THETADF(1-NAT) = angle THETADF (radians)

Arrays describing target rotation:
THETA(1-NTHETA)=desired values of target rotation angle theta
BETA(1-NBETA)=desired values of target rotation angle beta
PHI(1-NPHI)=desired values of target rotation angle phi
WGTA(1-NTHETA,1-NPHI)=weight factor for each orientation of target axis a1 in Lab Frame
WGTB(1-NBETA)=weight factor for rotation of target around a1

Arrays describing incident wave:
WAVEA(1-NWAV)=values of wavelength (physical units)
EN0_TF(3)=unit vector in direction of incidence in (rotated) Target Frame
AK_TF(3)=k*d vector in (rotated) Target Frame
XLR(3)=unit vector xlab in target frame
YLR(3)=unit vector ylab in target frame
ZLR(3)=unit vector zlab in target frame
Arrays describing scattering properties:
ENSC_LF(1-3,1-NSCAT)=unit scattering vectors in Lab Frame
ENSC_TF(1-3,1-NSCAT)=unit scattering vectors in Target Frame
AKS_TF(1-3,1-NSCAT)=scattering k-vectors in Target Frame
EM1_LF(1-3,1-NSCAT)=scattering polarization state 1 in Lab Frame
EM2_LF(1-3,1-NSCAT)=scattering polarization state 2 in Lab Frame
EM1_TF(1-3,1-NSCAT)=scattering polarization state 1 in Target Frame
EM2_TF(1-3,1-NSCAT)=scattering polarization state 2 in Target Frame
THETAN(1-NSCAT)=values of theta for NSCAT scattering directions nhat at which scattering matrix fml is to be calculated
For finite target (JPBC=0), THETAN is specified in the Lab Frame
PHIN(1-NSCAT)=values of phi for scattering NSCAT directions nhat at which scattering matrix fml is to be calculated
For finite target (JPBC=0), PHIN is specified in the Lab Frame.
S1111(1-NSCAT)=sum of |f11|**2 over target orientations
S2121(1-NSCAT)=       |f21|**2
QABS(1-2)=(absorption cross section)/(pi*aeff**2) for 2 polarizations
QBKSCA(1-2)=4*pi*(diff.scatt.cross section for theta=pi)/(pi*aeff**2) for 2 polarizations
QEXT(1-2)=(extinction cross section)/(pi*aeff**2)
QSCA(1-2)=(scattering cross section)/(pi*aeff**2)=QEXT-QABS
QSCAT(1-2)=(scattering cross section)/(pi*aeff**2) estimated by SCAT

Scratch array:
SCRRS1(1-NAT,3)=scratch array used by DDACCG and SCAT
SCRRS2(1-NAT)=scratch array used by SCAT
** */

// Set device numbers:
// IDVERR = device number for error messages
// IOPAR = device number for reading ddscat.par file
// IOSHP = device number for writing target.out (IOSHP<0 to suppress)
// IDVSHP = device number for reading shape file

// IDVOUT = device number for running I/O
//          on Solaris IDVOUT=0 produces unbuffered output to standard output
//                     IDVOUT=6 produces buffered output to standard output
// Initialize IDVOUT (cannot be done in DATA statement because IDVOUT
// is in COMMON/M6/

// 
// MPI environment setup:
	int ierr;
	MPI_Init(ierr);
	ierr = MPI_Comm_rank(MPI_Comm_world, &(Common10::GetInstance()->Myid()));
	ierr = MPI_Comm_size(MPI_Comm_world, &numprocs);

	if (!Common10::GetInstance()->Myid())
	{
		sprintf(cmsgnm, "numproc=%4d", numprocs);
		Wrimsg(DDscatLabel, cmsgnm);
	}

/* **
MPI_INIT is called by each process to initialize the MPI environment
MPI_COMM_RANK is called by each process to determine its identifier which is stored in MYID
MPI_COMM_SIZE is called by each process to determine the number of parallel processes, stored in NUMPROCS
If not using MPI, these calls can be turned into dummies by including
subroutines in mpi_fake.f (also need to comment out INCLUDE 'mpif.h'
and uncomment declaration of int variable MPI_COMM_WORLD as noted above)
** */
//
// Initialize variables to pass to SHARE1:
//	int mxpbc_sh = mxpbc;
//
// Call subroutine NAMID to obtain unique filename for output file written by this process.
	char cfllog[32];
	Namid(Common10::GetInstance()->Myid(), cfllog);
//
/* **
!=======================================================================
When following statement is enabled, output is buffered.
If using MPI, each MPI process will write to a unique file named
ddscat.log_xxx , where xxx=000,001,002, etc.
If it is necessary to have unbuffered output for debugging purposes,
comment out the following OPEN statement, in which case running
diagnostic information will be written to standard output.
Note that in this case, if using MPI, all MPI processes will write
to a single output stream.

//       OPEN(UNIT=IDVOUT,FILE=CFLLOG)
!=======================================================================
** */
	const char *cstamp = Version();
	const int versnum = VersionNum();

	if(!Common10::GetInstance()->Myid())								// ! begin if(myid==0)... #1
	{
		const char *Format9000 = " --- %s";
		sprintf(cmsgnm, Format9000, cstamp);
		Wrimsg(DDscatLabel, cmsgnm);
		sprintf(cmsgnm, "Versnum=%d", versnum);
		Wrimsg(DDscatLabel, cmsgnm);
#ifdef openmp
		sprintf(cmsgnmt, Format9000, "compiled with OpenMP enabled");
		Wrimsg(DDscatLabel, cmsgnm);
// ! fork a team of threads to ascertain number of OpenMP threads
!$OMP PARALLEL PRIVATE(NTHREADS,TID)
		tid = omp_get_thread_num();
//! only master thread does this:
		if(tid == 0)
		{
			nthreads = omp_get_num_threads();
			sprintf(cmsgnm, "Number of openmp threads = %3d\n", nthreads);
			Wrimsg(DDscatLabel, cmsgnm);
		}
// ! all threads join master thread and disband
!$OMP END PARALLEL
#endif

#ifdef DOUBLEPRECISION
		strcpy(cmsgnm, "Double-precision version");
#else
		strcpy(cmsgnm, "Single-precision version");
#endif
		nrword = sizeof(real);
		rword = (real)sizeof(real);					// NRWORD,RWORD = length (bytes) of REAL word
		cword = (real)2. * rword;					// cword = length (bytes) of Complex word
		Wrimsg(DDscatLabel, cmsgnm);
	}												// --- endif(myid==0)... #1

/* **
!**** Read program control parameters
CMDFFT=choice of solution methods
CMDFRM=specify scattering directions in either Lab Frame ('LFRAME') or target frame ('TFRAME')
CSHAPE=description of shape
SHPAR(1-10)= up to 10 parameters for defining target geometry
INIT=selection of |x0> in CCG method
TOL=error tolerance in CCG method
MXWAV=dimensioning information for array WAVEA
NWAV=number of wavelengths to be used
WAVEA(1-NWAV)=wavelengths (physical units)
MXRAD=dimensioning information for array AEFFA
NRAD=number of radii to be used
AEFFA(1-NRAD)=effective radii (physical units)
IORTH=1 or 2 for 1 or 2 incident polarization states to be used
NBETA=number of beta values to be used for target orientation
BETAMI=minimum value of beta (radians)
BETAMX=maximum value of beta (radians)
NTHETA=number of theta values to be used for target orientation
THETMI=minimum value of theta (radians)
THETMX=maximum value of theta (radians)
NPHI=number of phi values to be used for target orientation
PHIMIN=minimum value of phi (radians)
PHIMAX=maximum value of phi (radians)
DEGRAD=degrees/radian (conversion factor)
MXSCA=dimensioning information for arrays PHIN, THETAN
NSCAT=number of scattered directions for computation of fml
THETAN(1-NSCAT)=theta values (radians) for scattered directions for finite targets (JPBC=0), THETAN is in Lab Frame
PHIN(1-NSCAT)=phi values (radians) for scattered directions for finite targets (JPBC=0), PHIN is in Lab Frame
for JPBC=1 or 2 (infinite target with periodicity in 1-d): 
	PHIN(1-NSCAT) specifies scattering order 
	THETAN(1-NSCAT) specifies azimuthal angle of scattering direction around periodicity axis
for JPBC=3 (infinite target with periodicity in y_TF and z_TF)
	PHIN(1-NSCAT) specifies scattering order M
	THETAN(1-NSCAT) specifies scattering order N
** */

	Vect6<int> iMinmax, jMinmax, ipnf;
	Vect6<real> xMinmax;
	Vect3<real> a3;
	long file22IxyzPosition, poswave;
	int file22, file23, paramMxnx = 0, paramMxny = 0, paramMxnz = 0;
#ifdef _WIN32
	int modex = O_RDONLY|O_BINARY;
#else
	int modex = O_RDONLY;
#endif	
	if (!Common10::GetInstance()->Myid())	// begin if(myid==0)... #2
	{
 		int res = paramDDscat->Load(Common0::GetInstance()->Cflpar(), nrfld);
		if (res > 0)
		{
			fprintf(stderr, "Something went wrong in Load parameters, line = %d.\n", res);
			return;
		}
//
// Prepare to do nearfield calculation: geometric preliminaries
		if (nrfld == NearfieldMethodDo)
		{
			int i, delq[6];
			FileNamer::GetInstance()->Namer(paramDDscat->Iwav0(), paramDDscat->Irad0(), paramDDscat->Iori0());
			file22 = open(FileNamer::GetInstance()->GetCflpol1(), modex);
			if (file22 == -1)
			{
				Errmsg(DDscatLabel, "Error opeinf file22 in DDscat.", "");
			}
			read(file22, &nrword, sizeof(int));
			read(file22, &currentTarget->Nat0(), sizeof(int));
			read(file22, &currentTarget->Ianiso(), sizeof(int));
			iMinmax.Read(file22);
			read(file22, &currentTarget->Nx(), sizeof(int));			// NX,NY,NZ = size of original computational volume
			read(file22, &currentTarget->Ny(), sizeof(int));
			read(file22, &currentTarget->Nz(), sizeof(int));
///fprintf(stderr, "i: %5d %5d %5d %5d %5d %5d\n", iMinmax.data[0], iMinmax.data[1], iMinmax.data[2], iMinmax.data[3], iMinmax.data[4], iMinmax.data[5]);
//
// Here determine the actual size of the original target
// we wish to extend volume by 
// fraction SHPAR(1) in -xlab direction, fraction SHPAR(2) in +xlab direction
// fraction SHPAR(3) in -ylab direction, fraction SHPAR(4) in +ylab direction
// fraction SHPAR(5) in -zlab direction, fraction SHPAR(6) in +zlab direction
			delq[0] = nint_((iMinmax.data[1] - iMinmax.data[0] + 1) * paramDDscat->Extendxyz().data[0]);
			delq[1] = nint_((iMinmax.data[1] - iMinmax.data[0] + 1) * paramDDscat->Extendxyz().data[1]);
			delq[2] = nint_((iMinmax.data[3] - iMinmax.data[2] + 1) * paramDDscat->Extendxyz().data[2]);
			delq[3] = nint_((iMinmax.data[3] - iMinmax.data[2] + 1) * paramDDscat->Extendxyz().data[3]);
			delq[4] = nint_((iMinmax.data[5] - iMinmax.data[4] + 1) * paramDDscat->Extendxyz().data[4]);
			delq[5] = nint_((iMinmax.data[5] - iMinmax.data[4] + 1) * paramDDscat->Extendxyz().data[5]);
//fprintf(stderr, "dq: %5d %5d %5d %5d %5d %5d\n", delq[0], delq[1], delq[2], delq[3], delq[4], delq[5]);
//
			jMinmax.data[0] = iMinmax.data[0] - delq[0];
			jMinmax.data[1] = iMinmax.data[1] + delq[1];
			jMinmax.data[2] = iMinmax.data[2] - delq[2];
			jMinmax.data[3] = iMinmax.data[3] + delq[3];
			jMinmax.data[4] = iMinmax.data[4] - delq[4];
			jMinmax.data[5] = iMinmax.data[5] + delq[5];
//
// Define new computational volume to allow near-field calculation
			paramMxnx = jMinmax.data[1] - jMinmax.data[0] + 1;
			paramMxny = jMinmax.data[3] - jMinmax.data[2] + 1;
			paramMxnz = jMinmax.data[5] - jMinmax.data[4] + 1;
//
// If CMETHD = GPFAFT, need to allow for requirement that extended NX,NY,NZ be factorizable by 2,3,5
			if (paramDDscat->Cmdfft() == FftMethod_GPFAFT)
			{
				paramMxnx = GetNearestFactor(paramMxnx);
				paramMxny = GetNearestFactor(paramMxny);
				paramMxnz = GetNearestFactor(paramMxnz);
//
// In case MXNX,MXNY,MXNZ have changed, update JXMAX,JYMAX,JZMAX
				jMinmax.data[1] = paramMxnx + jMinmax.data[0] - 1;
				jMinmax.data[3] = paramMxny + jMinmax.data[2] - 1;
				jMinmax.data[5] = paramMxnz + jMinmax.data[4] - 1;
			}
// 
// Do not close file, because we will do more reading below
//
// Continue reading file 22, We will not allocate memory until call to Target
//
// Prepare to arrays allocation, set anchol in file 22
#ifdef _WIN32
			file22IxyzPosition = tell(file22);
#else
			file22IxyzPosition = lseek(file22, (long)0, SEEK_CUR);
#endif			
//
// Read ixyz in dummy variables to determine ipnf only
			Vect3<int> temp;
			iMinmax.Set(ihuge_, -ihuge_, ihuge_, -ihuge_, ihuge_, -ihuge_);
			for(i=0; i<currentTarget->Nat0(); ++i)
			{
				temp.Read(file22);
				if (temp.data[0] < iMinmax.data[0]) iMinmax.data[0] = temp.data[0];
				if (temp.data[0] > iMinmax.data[1]) iMinmax.data[1] = temp.data[0];
				if (temp.data[1] < iMinmax.data[2]) iMinmax.data[2] = temp.data[1];
				if (temp.data[1] > iMinmax.data[3]) iMinmax.data[3] = temp.data[1];
				if (temp.data[2] < iMinmax.data[4]) iMinmax.data[4] = temp.data[2];
				if (temp.data[2] > iMinmax.data[5]) iMinmax.data[5] = temp.data[2];
			}
// 
// IPNF(1) = padding in -x direction to create near-field zone, IPNF(2) = padding in +x direction to create near-field zone
// IPNF(3) = padding in -y direction to create near-field zone, IPNF(4) = padding in +y direction to create near-field zone
// IPNF(5) = padding in -z direction to create near-field zone, IPNF(6) = padding in +z direction to create near-field zone
			ipnf.data[0] = iMinmax.data[0] - jMinmax.data[0];
			ipnf.data[1] = jMinmax.data[1] - iMinmax.data[1];
			ipnf.data[2] = iMinmax.data[2] - jMinmax.data[2];
			ipnf.data[3] = jMinmax.data[3] - iMinmax.data[3];
			ipnf.data[4] = iMinmax.data[4] - jMinmax.data[4];
			ipnf.data[5] = jMinmax.data[5] - iMinmax.data[5];
///fprintf(stderr, "p: %5d %5d %5d %5d %5d %5d\n", ipnf.data[0], ipnf.data[1], ipnf.data[2], ipnf.data[3], ipnf.data[4], ipnf.data[5]);
// 
// Return to the last fixed position in file 22
			lseek(file22, file22IxyzPosition, SEEK_SET);
//
// Now we can allocate all arrays in target parameters
			int sizeX = iMinmax.data[1] - iMinmax.data[0] + 1 + ipnf.data[0] + ipnf.data[1]; 
			int sizeY = iMinmax.data[3] - iMinmax.data[2] + 1 + ipnf.data[2] + ipnf.data[3]; 
			int sizeZ = iMinmax.data[5] - iMinmax.data[4] + 1 + ipnf.data[4] + ipnf.data[5]; 
			if (paramDDscat->Cmdfft() == FftMethod_GPFAFT)
			{
				sizeX = GetNearestFactor(sizeX); 
				sizeY = GetNearestFactor(sizeY); 
				sizeZ = GetNearestFactor(sizeZ); 
			}
			currentTarget->AllocateArrays(sizeX, sizeY, sizeZ, true);			// Dimensions are quite precise
			currentTarget->SetMin(iMinmax);
//
// Continue loading target
//
///fprintf(stderr, "j: %5d %5d %5d %5d %5d %5d\n", jMinmax.data[0], jMinmax.data[1], jMinmax.data[2], jMinmax.data[3], jMinmax.data[4], jMinmax.data[5]);
///fprintf(stderr, "i: %5d %5d %5d %5d %5d %5d\n", iMinmax.data[0], iMinmax.data[1], iMinmax.data[2], iMinmax.data[3], iMinmax.data[4], iMinmax.data[5]);
//fprintf(stderr, "s: %5d %5d %5d\n", shiftt.data[0], shiftt.data[1], shiftt.data[2]);
//
// Read more data stored in file 22 and obtain X0
			int *absoluteAddresses = new int[currentTarget->Nat0()];
			real delz[3];
			delz[0] = (1 - iMinmax.data[0] + ipnf.data[0]);
			delz[1] = (1 - iMinmax.data[2] + ipnf.data[2]);
			delz[2] = (1 - iMinmax.data[4] + ipnf.data[4]);
///fprintf(stderr, "i: %5d %5d %5d\n", (1 - iMinmax.data[0] + ipnf.data[0]), (1 - iMinmax.data[2] + ipnf.data[2]), (1 - iMinmax.data[4] + ipnf.data[4]));
///fprintf(stderr, "i: %5d %5d %5d %5d %5d %5d\n", ipnf.data[0], ipnf.data[1], ipnf.data[2], ipnf.data[3], ipnf.data[4], ipnf.data[5]);
///fprintf(stderr, "d: %lf %lf %lf\n", delz[0], delz[1], delz[2]);
///fprintf(stderr, "u: %lf %lf %lf\n", 1 + currentTarget->GetMinJx() - jMinmax.data[0], 1 + currentTarget->GetMinJy() - jMinmax.data[2], 1 + currentTarget->GetMinJz() - jMinmax.data[4]);
			for(i=0; i<currentTarget->Nat0(); ++i)
			{
				temp.Read(file22);
				temp.data[0] += delz[0];
				temp.data[1] += delz[1];
				temp.data[2] += delz[2];
				currentTarget->Ixyz().Fill3(i, temp.data[0], temp.data[1], temp.data[2]);
				absoluteAddresses[i] = currentTarget->GetLinearAddress(i);
			}
///fprintf(stderr, "\n\n");
///fprintf(stderr, "%d %d %d %d %d\n", currentTarget->Nat0(), currentTarget->Nat(), currentTarget->Icomp().GetTotalSize(), currentTarget->Icomp().GetSize(0), currentTarget->Icomp().GetSize(1));
///fprintf(stderr, "%d %d %d\n", currentTarget->Nx(), currentTarget->Ny(), currentTarget->Nz());
///fprintf(stderr, "\n\n");
///fflush(stderr);
			read(file22, &i, sizeof(int));
			manager->Jpbc() = (PeriodicBoundaryFlag)i;
			read(file22, &currentTarget->Pyd(), sizeof(real));
			read(file22, &currentTarget->Pzd(), sizeof(real));
			read(file22, &paramDDscat->Gamma(), sizeof(real));
			currentTarget->A1().Read(file22);
			currentTarget->A2().Read(file22);
			currentTarget->Dx().Read(file22);
			currentTarget->X0().Read(file22);
//currentTarget->X0().Fprintf(stderr, "%lf");
			for(i=0; i<3; ++i)
			{
				currentTarget->X0().data[i] -= delz[i];
			}
//currentTarget->X0().Fprintf(stderr, "%lf");
			Vect3<short> shortTemp;
			for(i=0; i<currentTarget->Nat0(); ++i)
			{
				shortTemp.Read(file22);
///if (i < 32)
///fprintf(stderr, "%6d %4d %4d %4d %6d\n", i, shortTemp.data[0], shortTemp.data[1], shortTemp.data[2], absoluteAddresses[i]);
				currentTarget->Icomp().Fill3(absoluteAddresses[i], shortTemp.data[0], shortTemp.data[1], shortTemp.data[2]);
				currentTarget->Iocc(absoluteAddresses[i]) = true;
///if (i == 32) fflush(stderr);
			}
			read(file22, &paramDDscat->Nambient(), sizeof(real));
#ifdef _WIN32			
			poswave = tell(file22);
#else
			poswave = lseek(file22, (long)0, SEEK_CUR);
#endif			
			close(file22);
			delete [] absoluteAddresses;
		}
		else
		{
//
// Load target data from Target manager, allocate arrays
			manager->CashedDx() = paramDDscat->Dx();
			currentTarget = manager->LoadTarget(paramDDscat->TargetString(), paramDDscat->Cflshp());
// 
// Return from TARGET with following arrays defined:
//      NAT0   = number of occupied lattice sites
// IXYZ0(J,1-3)= (IX,IY,IZ) for occupied sites J=1-NAT0
//     X0(1-3) = specifies location in TF of lattice site (IX,IY,IZ)=(0,0,0) :
//       x_TF  = (IX+X0(1))*d*DX(1)
//       y_TF  = (IY+X0(2))*d*DX(2)
//       z_TF  = (IZ+X0(3))*d*DX(3)
//     A1(1-3) = target axis A1 in TF
//     A2(1-3) = target axis A2 in TF
//       PYD   = periodicity/[d*DX(2)] in y_TF direction (for JPBC=1 or 3)
//       PZD   = periodicity/[d*DX(3)] in z_TF direction (for JPBC=2 or 3)
// ICOMP(J,1-3)= composition at lattice site J, for directions 1-3 in local "Dielectric Frame" where dielectric tensor is diagonalized
// 	             dimension ICOMP(3*MXNAT0)
//
//               J = 1 - NAT0                     =comp for 1x,2x,3x,..,NAT0x
//               J = NAT0+1 - MXNAT0              =0
//               J = MXNAT0+1 - MXNAT0+NAT0       =comp for 1y,2y,3y,..,NAT0y
//               J = MXNAT0+NAT0+1 - 2*MXNAT0     =0
//               J = 2*MXNAT0+1 - 2*MXNAT0+NAT0   =comp for 1z,2z,3z,...,NAT0z
//               J = 2*MXNAT0+NAT0+1 - 3*MXNAT0   =0
// 
// for targets composed of general anisotropic materials (IANISO=2), also define:
//    BETADF(J) = rotation angle beta defining orientation in TF of "dielectric frame" (DF) for site J
//     PHIDF(J) = rot. angle phi defining orientation in TF of DF for site J
//   THETADF(J) = rot. angle theta defining orientation in TF of DF for site J
			sprintf(cmsgnm, "     pyd=%lf     pzd=%lf", currentTarget->Pyd(), currentTarget->Pzd());
			Wrimsg("DDscat", cmsgnm);
//
// Check to see if NCOMP_NEED is compatible with NCOMP
			int ncomp_need = manager->IcompNeed();
			if (ncomp_need > manager->Ncomp())
			{
				sprintf(cmsgnm, "error:ncomp_need=%4d >%4d = ncomp\n", ncomp_need, manager->Ncomp());
				Wrimsg("DDdscat", cmsgnm);
				sprintf(cmsgnm, "ddscat.par does not list enough different refractive index filenames\n");
				Errmsg("Fatal", DDscatLabel, "ncomp is too small\n");
			}
		}
		iMinmax = currentTarget->IxyzMinmax();							// Get min and max values from param->ixyz
		int iu;
		for(iu=0; iu<3; ++iu)
		{
			xMinmax.data[2*iu]   = (real)iMinmax.data[2*iu]   - half_ + currentTarget->X0().data[iu];
			xMinmax.data[2*iu+1] = (real)iMinmax.data[2*iu+1] + half_ + currentTarget->X0().data[iu];
		}
//
// Print minimax information
		Wrimsg(DDscatLabel, "dimensions of physical target:");
		for(iu=0; iu<3; ++iu)
		{
			sprintf(cmsgnm, "%6d%6d = min, max values of j%c", iMinmax.data[2*iu], iMinmax.data[2*iu+1], 'x'+iu);
			Wrimsg(DDscatLabel, cmsgnm);
		}
		currentTarget->X0().Sprintf(cmsgnm, "%10.5lf", NULL, " = x0(1-3)");
		Wrimsg(DDscatLabel, cmsgnm);
		Wrimsg(DDscatLabel, " --- physical extent of target volume (occupied sites) ---");
		for(iu=0; iu<3; ++iu)
		{
			sprintf(cmsgnm, "%12.5lf%12.5lf = min, max values of %c/d", xMinmax.data[2*iu], xMinmax.data[2*iu+1], 'x'+iu);
			Wrimsg(DDscatLabel, cmsgnm);
		}
//
//  Calculate target axis A3 in Target Frame
		a3 = currentTarget->A1() * currentTarget->A2();
//
		pyddx = currentTarget->Pyd() * currentTarget->Dx().data[1];
		pzddx = currentTarget->Pzd() * currentTarget->Dx().data[2];
//
// Calculate DAEFF = d/aeff for this target (d=lattice spacing)
		daeff = Pow((FourPi / (3 * currentTarget->Nat0())), (real)(1./3.));
//
// We don't need Extend any more, anyway the following information is valuable, as now all memory is allocated by Target or during read from file 22
//
// EXTEND reorders the occupied lattice sites, returns the extents NX,NY,NZ for the extended target.
// 
// if necessary, shifts JX,JY,JZ for each site to that JX,JY,JZ run from 1-NX,1-NY,1-NZ
// IXYZ0(JA,1-3) = coordinates JX,JY,JZ of occupied lattice sites JA=1-NAT0
// location JX,JY,JZ in extended target corresponds to 
// JE=NX*NY*(JZ-1)+NX*(JY-1)+JX, where JE runs from 1 to NX*NY*NZ
//    ICOMP(JE,1-3)   = composition for each site in extended target
//                    = 0 if site is unoccupied
//       IOCC(JE)     = 1 if occupied, 0 if not
// BETADF(JA=1-NAT0)  = rotation angle beta for DF in reordered physical target for JA=1-NAT0
//  PHIDF(JA=1-NAT0)  = rotation angle phi for DF in reordered physical target
// THETADF(JA=1-NAT0) = rotation angle theta for DF in reordered physical target
//    X0(1-3)         = location in TF of lattice site JX,JY,JZ=0,0,0
// 
// Memory reallocation inside TargetParameters is not necessary
//		paramDDscat->Mxnx() = currentTarget->Nx();
//		paramDDscat->Mxny() = currentTarget->Ny();
//		paramDDscat->Mxnz() = currentTarget->Nz();
		mxnat = currentTarget->Nx() * currentTarget->Ny() * currentTarget->Nz();
		mxn3 = 3 * mxnat;
//
// MXCXSC specifies required Complex scratch array CXSC
// increase from 7*MXN3 required by PETR
//           to 10*MXN3 required by PIMCBICGSTAB method
//                 (PIM Complex BiConjugate Gradient with Stabilization)
// ver7.2:
// increase to  12*MXN3 required by GPBICG
//                 (Tang et al CCG method implemented by Chaumet & Rahmani)

		if (paramDDscat->Cmdsol() == SolMethod_End)
		{
            sprintf(cmsgnm, "unrecognized CMDSOL=%s", SolEnumerator(paramDDscat->Cmdsol()));
			Errmsg("Fatal", DDscatLabel, cmsgnm);
		}
		mxcxsc = 2 * mxn3;
// 
// Pseudo names of scratch array  positions:
		lace = 0;
		laxi = 1 * mxn3;
		lclm = 2 * mxn3;
		lgi  = 3 * mxn3;
		lpi  = 4 * mxn3;
		lqi  = 5 * mxn3;
		lsc0 = 6 * mxn3;
	}								// !--- endif(myid==)... #2
//
// Allocate arrays according to their dimensions given by Target
	if (paramDDscat->Nscat() > 0)
	{
		int size = paramDDscat->Nscat();
		cxfData = new FourArray;
		cxfData->Alloc(size);
		orderm = new real[size];
		ordern = new real[size];
		s1111 = new real[size];
		s1111_1 = new real[size];
		s2121 = new real[size];
		s2121_1 = new real[size];
		cx1121 = new Complex[size];
		cx1121_1 = new Complex[size];
		aks_tf = new Vect3<real>[size];
		em1_lf = new Vect3<real>[size];
		em1_tf = new Vect3<real>[size];
		em2_lf = new Vect3<real>[size];
		em2_tf = new Vect3<real>[size];
		ensc_lf = new Vect3<real>[size];
		ensc_tf = new Vect3<real>[size];
		smori_1 = new real[4 * 4 * size];
		sm = new real[4 * 4 * size];
		smori = new real[4 * 4 * size];
	}
	if (manager->Ncomp() > 0)
	{
		DielectricManager::GetInstance()->Allocate(manager->Ncomp());
	}
	DielectricManager *dielec = DielectricManager::GetInstance();
// 
// This MPI_Barrier is required to ensure that LACE, etc have been evaluated by MYID=0 before being SHAREd by other threads
	ierr = MPI_Barrier(MPI_Comm_world);
//
// Share info between computers      
	Share0(lace, laxi, lclm, lgi, lpi, lqi, mxn3, mxnat, paramMxnx, paramMxny, paramMxnz, (manager->Jpbc() != PeriodicNo) ? 1 : 0, mxcxsc, Common10::GetInstance()->Myid(), currentTarget->Nat0());
//
// Begin allocate memory for scratch etc. arrays
// Cxzw is allocated and zerofied during first call to Getfml, accessible via Common4
	real *scrrs2 = NULL;
	DipoleData *theDipoleData = NULL;
	Complex *cxscr1 = NULL;
	Complex *cxsc = NULL;
	Matrix *theMatrix = NULL;
	Matrix *theTensor = NULL;
	if (!Common10::GetInstance()->Myid())
	{
		scrrs2 = (real *)malloc(currentTarget->Nat() * sizeof(real));
		cxsc = new Complex[mxcxsc];
		Common4::GetInstance()->AllocateCxzw(2 * currentTarget->Nx(), 2 * currentTarget->Ny(), 2 * currentTarget->Nz(), 3);
		theDipoleData = new DipoleData;
		theDipoleData->Allocate(3 * currentTarget->Nat());
//fprintf(stderr, " --- CCC\n");
//theDipoleData->Debug(stderr);
		cxscr1 = new Complex[3 * currentTarget->Nat()];
		theMatrix = new Matrix;
		theMatrix->Allocate(3 * currentTarget->Nat());
		paramDDscat->SetMatrix(theMatrix);
		theTensor = new Matrix;
		theTensor->Allocate(3 * currentTarget->Nat());
	}
// End allocate memory

//
// Instantiate Output manager
	OutputManager *outputer = OutputManager::GetInstance();
	outputer->Init(currentTarget, paramDDscat, dielec, manager->Jpbc());
//
	if (!Common10::GetInstance()->Myid())															// !--- begin if(myid==0)... #23
	{
		if(manager->Jpbc() != PeriodicNo)
		{
			memcpy(orderm, paramDDscat->Phin(), paramDDscat->Nscat() * sizeof(real));
			memcpy(ordern, paramDDscat->Thetan(), paramDDscat->Nscat() * sizeof(real));
		}
// 
// JPBC = 0 if PBC not used
//      = 1 for PBC in y direction only
//      = 2 for PBC in z direction only
//      = 3 for PBC in y and z directions
//
// IPBC controls memory use.
// At this time use same approach for JPBC=1 and 2 as for JPBC=3
// May change this in future.
		switch(manager->Jpbc())
		{
		case PeriodicNo:
			GreenFunctionManager::GetInstance()->SetIpbc(false);
			break;

		case PeriodicY:
		case PeriodicZ:
			if (paramDDscat->Cmdfrm() == FrameCode_LFRAME)
				Errmsg("Fatal", DDscatLabel, " lframe = invalid option for periodic target (jpbc != PeriodicNo)\n");
            GreenFunctionManager::GetInstance()->SetIpbc(true);
			break;

		case PeriodicBoth:
            if (paramDDscat->Cmdfrm() == FrameCode_LFRAME)
				Errmsg("Fatal", DDscatLabel, " lframe = invalid option for periodic target (jpbc != PeriodicNo)\n");
            GreenFunctionManager::GetInstance()->SetIpbc(true);
			break;

		default:
			break;
		}
		GreenFunctionManager::GetInstance()->SetDimension(currentTarget->Nx(), currentTarget->Ny(), currentTarget->Nz());
		GreenFunctionManager::GetInstance()->AllocateElectricBuffer();
//
		nbeth = oridata->GetNtheta() * oridata->GetNbeta();
		nori  = nbeth * oridata->GetNphi();
//
// If JPBC=0, obtain scattering vectors and scattering pol. vectors in Lab Frame
		if (manager->Jpbc() == PeriodicNo)
		{
			paramDDscat->Scavec(ensc_lf, em1_lf, em2_lf);
// if JPBC = 0:
// ENSC_LF(1-3,1-NSCAT) = scattering vectors in Lab Frame
// EM1_LF(1-3,1-NSCAT)  = scattering polarization state 1 in Lab Frame
// EM2_LF(1-3,1-NSCAT)  = scattering polarization state 2 in Lab Frame
		}
		else
		{
// if JPBC > 0:
// initialize ENSC_LF,EM1_LF,EM2_LF to ensure that there are no divisions by zero
// in ROTATE.  Note that for JPBC > 0 these vectors are not actually used
// AKS_TF,EM1_TF,EM2_TF will be calculated later by subroutine PBCSCAVEC
			for(j=0; j<paramDDscat->Nscat(); ++j)
			{
				ensc_lf[j].Set(onex_, zero_, zero_);
				em1_lf[j].Set(zero_, onex_, zero_);
				em2_lf[j].Set(zero_, zero_, onex_);
			}
		}
//
// Note: IXYZ0(1-NAT0,1-3) contain coordinates of occupied lattice sites
//       [reordered to standard order used for arrays produced by REDUCE]
//       BETADF(1-NAT0),PHIDF(1-NAT0),THETADF(1-NAT0) contain orientation
//                         angles for anisotropic material at occupied
//                         lattice sites (not used for isotropic material
		const char *Format9011 = "       %7d = NAT  = number of dipoles in extended target\n       %4d%4d%4d = x,y,z length of extended target (Targ. Frame)";
		sprintf(cmsgnm, Format9011, currentTarget->Nat(), currentTarget->Nx(), currentTarget->Ny(), currentTarget->Nz());
		Wrimsg(DDscatLabel, cmsgnm);
//
// Specify target orientations
// MXBETA=dimensioning information for array BETA of beta values
// MXTHET=dimensioning information for array THETA of theta values
// MXPHI=dimensioning information for array PHI of phi values
// BETAMI,BETAMX=minimum,maximum beta values (radians)
// THETMI,THETMX=minimum,maximum theta values (radians)
// PHIMIN,PHIMAX=minimum,maximum phi values (radians)
// BETA(1-NBETA)=beta values for target orientations (radians)
// THETA(1-NTHETA)=theta values for target orientations (radians)
// PHI(1-NPHI)=phi values for target orientations (radians)
// WGTA(1-NTHETA,1-NPHI)=weight for each orientation of axis a1 in LF
// WGTB(1-NBETA)=weight for each rotation of target around a1
//
		paramDDscat->Orient();
//
// Open file 'qtable' for summary of Q values and
//           'mtable' for summary of refractive indices
// except if NRFLD=2, in which case these files have already been written during previous pass with NRFLD=1
//
// Write out header for 'qtable', 'qtable2', and 'mtable' files
		if (nrfld != NearfieldMethodDo)
		{
			outputer->WriteHeadings(nori);
		}
	}									// !--- end if(myid==0)... #23
// The master process needs to share information from the above section with the slave processes. All processes call this subroutine, which does nothing if MPI is not being used.
	ierr = MPI_Barrier(MPI_Comm_world);
//
// Why not use cbinfile ? and open it in writesca ?
//	int iobin = -1;
//	if (paramDDscat->Cbinflag() == BinflagMethod_ALLBIN || paramDDscat->Cbinflag() == BinflagMethod_ORIBIN)
//	{
//		iobin = open("dd.bin", O_CREAT|O_WRONLY|O_BINARY, S_IWRITE);
//		if (iobin == -1)
//		{
//			Errmsg(DDscatLabel, "Error opening iobin file in DDscat.", "");
//		}
//	}
//
// ChB 7.2.1. nambient is a parameter
	Share1(currentTarget, paramDDscat, dielec->GetFileNames(), 
		daeff, ensc_lf, em1_lf, em2_lf, manager->Ioshp(), GreenFunctionManager::GetInstance()->Ipbc(), 
		manager->Jpbc(), mxn3, mxnat, Common10::GetInstance()->Myid(), nbeth, manager->Ncomp(), nori, orderm, ordern, 
		currentTarget->Pyd(), pyddx, currentTarget->Pzd(), pzddx, xMinmax);
//
// IWAV0,IRAD0,IORI0 were obtained from REAPAR
	iwav1 = paramDDscat->Iwav0();
	irad1 = paramDDscat->Irad0();
	iori1 = paramDDscat->Iori0();
	itheta1 = paramDDscat->Iori0() / (oridata->GetNbeta() * oridata->GetNphi());
	ibeta1  = paramDDscat->Iori0() /  oridata->GetNphi() - itheta1 * oridata->GetNbeta();
	ibeth1  = paramDDscat->Iori0() /  oridata->GetNphi();
	iphi1 = iori1 - (itheta1 * oridata->GetNbeta() + ibeta1) * oridata->GetNphi();
//
// Loop over wavelengths:
	for(iwav=iwav1; iwav<paramDDscat->Nwav(); ++iwav)
	{
		if(nrfld == NearfieldMethodDo)
		{
			FileNamer::GetInstance()->Namer(iwav, irad1, iori1);
			int file22 = open(FileNamer::GetInstance()->GetCflpol1(), modex);
			if (file22 == -1)
				Errmsg(DDscatLabel, "Error opening file22 in DDscat.", "");
			lseek(file22, (long)poswave, SEEK_SET);
			read(file22, &wave, sizeof(real));
#ifdef _WIN32
			posaeff = tell(file22);
#else
			posaeff = lseek(file22, (long)0, SEEK_CUR);
#endif					
            close(file22);
		}
		else
			wave = paramDDscat->Wavea(iwav);
//
// Determine dielectric properties of target
		if (!Common10::GetInstance()->Myid())
			dielec->PrepareForWave(wave);
// 
// Again the master process needs to share information with the slaves
// 080723 BTD add mpi_barrier to ensure that E1A,E2A have been read by MYID=0 before being SHAREd by other threads.
		ierr = MPI_Barrier(MPI_Comm_world);
		Share2(dielec);
//
// Obtain Complex refractive index for the constituent materials
		const char *Format9031 = " n= (%7.4lf ,%7.4lf),  epsilon= (%8.4lf , %7.4lf)  for material %2d";
 		for(j=0; j<manager->Ncomp(); ++j)
		{
			Complex tr = dielec->GetCxrfr(j);
			Complex te = dielec->GetCxeps(j);
            if (!Common10::GetInstance()->Myid())
			{
				sprintf(cmsgnm, Format9031, tr.re, tr.im, te.re, te.im, j);
				Wrimsg(DDscatLabel, cmsgnm);
			}
		}
//
// Now convert to relative refractive indices
		if(Fabs(paramDDscat->Nambient() - onex_) > (real)0.000001)
			dielec->UseNambient(paramDDscat->Nambient());
//
// Loop over irad
		for(irad=irad1; irad<paramDDscat->Nrad(); ++irad)
		{
			if(nrfld == NearfieldMethodDo)
			{
				FileNamer::GetInstance()->Namer(iwav, irad, iori1);
				int file22 = open(FileNamer::GetInstance()->GetCflpol1(), modex);
				if (file22 == -1)
					Errmsg(DDscatLabel, "Error opening file22 in DDscat.", "");
				lseek(file22, (long)posaeff, SEEK_SET);
				read(file22, &aeff, sizeof(aeff));
#ifdef _WIN32				
	            posak_tf = tell(file22);
#else
	            posak_tf = lseek(file22, (long)0, SEEK_CUR);
#endif				
				close(file22);
			}
            else
				aeff = paramDDscat->Aeffa()[irad];
            xx = TwoPi * aeff / (wave / paramDDscat->Nambient());
//
// Compute AK1=length of k vector in "natural" units=k*d (natural unit of length = d = lattice spacing)
// Remember that NAT = number of dipoles including "vacuum" sites in extended target NAT0= number of dipoles in "real" target
			ak1 = xx * daeff;
            if(!Common10::GetInstance()->Myid())
			{
				const char *Format9032 = "%12.6lf = AEFF = effective radius (physical units)\n >DDSCAT %12.6lf = d/aeff for this target\n >DDSCAT %12.6lf = WAVE = wavelength in vacuo (physical units)\n >DDSCAT %12.6lf = wavelength in ambient medium (phys. units)\n >DDSCAT %12.6lf = k_amb*aeff = 2*pi*aeff/lambda_amb\n >DDSCAT %12.6lf = k_amb*d";
				sprintf(cmsgnm, Format9032, aeff, daeff, wave, wave/paramDDscat->Nambient(), xx, ak1);
				Wrimsg(DDscatLabel, cmsgnm);
			}
			ak3 = ak1 * ak1 * ak1;
//
// PIA2=pi*(aeff/d)**2
//			pia2 = Pi * Pow((real)0.75 * currentTarget->Nat0() / Pi, (real)(2./3.));
//
// Initialize various sums over target orientation.
			memset(s1111, 0, paramDDscat->Nscat() * sizeof(real));
			memset(s2121, 0, paramDDscat->Nscat() * sizeof(real));
			memset(s1111_1, 0, paramDDscat->Nscat() * sizeof(real));
			memset(s2121_1, 0, paramDDscat->Nscat() * sizeof(real));
			memset(cx1121, 0, paramDDscat->Nscat() * sizeof(Complex));
			memset(cx1121_1, 0, paramDDscat->Nscat() * sizeof(Complex));

			memset(smori, 0, 16 * paramDDscat->Nscat() * sizeof(real));
			memset(smori_1, 0, 16 * paramDDscat->Nscat() * sizeof(real));
//
// Loop over target orientations (target rotations in Lab Frame) (actually, loop over lab rotations in Target Frame)
// The outer orientation loop over IBETH=1,NBETH is divided up among the parallel processes.
			iori = paramDDscat->Iori0();
//
// ITNUMMX will store the maximum number of iterations used for any of the orientations for current size and wavelength, so that this can be written to the waarbbkccc.sca output file
			itnummx[0] = itnummx[1] = 0;
//
// Loop over ibeth
			for(ibeth=Common10::GetInstance()->Myid()+ibeth1; ibeth<nbeth; ibeth+=numprocs)
			{
				itheta = (int)(ibeth / oridata->GetNbeta());
				outputer->SetThetad(oridata->Theta(itheta) * Degrad);
				ibeta = ibeth % oridata->GetNbeta();
				outputer->SetBetad(oridata->Beta(ibeta) * Degrad);
				sinbet = Sin(oridata->Beta(ibeta));
				cosbet = Cos(oridata->Beta(ibeta));
				sinthe = Sin(oridata->Theta(itheta));
				costhe = Cos(oridata->Theta(itheta));
//
// Loop over phi
				for(iphi=iphi1; iphi<oridata->GetNphi(); ++iphi)
				{
					outputer->SetPhid(oridata->Phi(iphi) * Degrad);
					iori = (itheta * oridata->GetNbeta() + ibeta) * oridata->GetNphi() + iphi;
//
					if(nrfld == NearfieldMethodDo)
					{
						FileNamer::GetInstance()->Namer(iwav, irad, iori);
						int file22 = open(FileNamer::GetInstance()->GetCflpol1(), modex);
						if (file22 == -1)
						{
							Errmsg(DDscatLabel, "Error opening file22 in DDscat.", "");
						}
						lseek(file22, (long)posak_tf, SEEK_SET);
						Common1::GetInstance()->Ak_tf().Read(file22);
#ifdef _WIN32				
	            poscxe0 = tell(file22);
#else
	            poscxe0 = lseek(file22, (long)0, SEEK_CUR);
#endif							
						cxe01_tf.Read(file22);
						if(paramDDscat->Iorth() == 2)
						{
							file23 = open(FileNamer::GetInstance()->GetCflpol2(), modex);
							if (file23 == -1)
								Errmsg(DDscatLabel, "Error opening file23 in DDscat.", "");
							lseek(file23, (long)poscxe0, SEEK_SET);
                            cxe02_tf.Read(file23);
						}
					}										
//        
// For specified target rotations BETA,THETA,PHI in the Lab Frame, and predefined 
// 
//     ENSC_LF(1-3,1-NSCAT) scattering directions in the Lab Frame 
//     EM1_LF(1-3,1-NSCAT) scattered pol. vectors 1 in the Lab Frame
//     EM2_LF(1-3,1-NSCAT) scattered pol. vectors 2 in the Lab Frame
// 
//     subroutine ROTATE computes
// 
//     AKS_TF(1-3,1-NSCAT) propagation vectors in the Target Frame
//     EM1_TF(1-3,1-NSCAT) scattered pol vectors 1 in Target Frame
//     EM2_TF(1-3,1-NSCAT) scattered pol vectors 2 in Target Frame
// 
//     Note: subroutine ROTATE assumes that the angles THETAN,PHIN define
//     scattering directions in the Lab Frame -- if these values are
//     actually for the Target Frame, values of AKS_TF,EM1_TF,EM2_TF returned
//     by ROTATE will later be replaced by correct values.
					Rotate(currentTarget->A1(), currentTarget->A2(), ak1, oridata->Beta(ibeta), oridata->Theta(itheta), oridata->Phi(iphi), 
						en0_tf, cxe01_tf, cxe02_tf, paramDDscat->Nscat(), ensc_lf, em1_lf, em2_lf, aks_tf, em1_tf, em2_tf);
// fprintf(stderr, "Rotate: %lf %lf %lf %lf %lf %lf\n", cxe01_tf.data[0].re, cxe01_tf.data[0].im, cxe01_tf.data[1].re, cxe01_tf.data[1].im, cxe01_tf.data[2].re, cxe01_tf.data[2].im);
// fprintf(stderr, "Rotate: %lf %lf %lf %lf %lf %lf\n", cxe02_tf.data[0].re, cxe02_tf.data[0].im, cxe02_tf.data[1].re, cxe02_tf.data[1].im, cxe02_tf.data[2].re, cxe02_tf.data[2].im);
//
// Calculate ENSC_TF = unit scattering vector in Target Frame
					for(idir=0; idir<paramDDscat->Nscat(); ++idir)
					{
						ensc_tf[idir].Copy(aks_tf[idir]);
						ensc_tf[idir] /= ak1;
					}
//
// Calculate XLR,YLR,ZLR = lab unit vectors XL,YL,ZL in Target Frame
					sinphi = Sin(oridata->Phi(iphi));
					cosphi = Cos(oridata->Phi(iphi));
//					
// If CMDFRM = 'TFRAME' then angles THETAN,PHIN define scattering directions in the Target Frame, so ignore vectors AKS_TF,EM1_TF,EM2_TF calculated by ROTATE, and instead simply use ENSC_LF,EM1_LF,EM2_LF
					if (paramDDscat->Cmdfrm() == FrameCode_TFRAME)
					{
						if(manager->Jpbc() == PeriodicNo)
						{
							for(idir=0; idir<paramDDscat->Nscat(); ++idir)
							{
								aks_tf[idir].Copy(ensc_lf[idir]);
								aks_tf[idir] *= ak1;
								em1_tf[idir].Copy(em1_lf[idir]);
								em2_tf[idir].Copy(em2_lf[idir]);
							}
						}
						else
						{
							Vect3<real>xlr = currentTarget->A1() * costhe          - currentTarget->A2() * sinthe * cosbet + a3 * sinthe * sinbet;
							Vect3<real>ylr = currentTarget->A1() * sinthe * cosphi + currentTarget->A2() * (costhe * cosbet * cosphi - sinbet * sinphi) - a3 * (costhe * sinbet * cosphi + cosbet * sinphi);
							Vect3<real>zlr = currentTarget->A1() * sinthe * sinphi + currentTarget->A2() * (costhe * cosbet * sinphi + sinbet * cosphi) - a3 * (costhe * sinbet * sinphi - cosbet * cosphi);
//
// JPBC= 1, 2, or 3: call PBCSCAVEC to calculate scattering vectors AKS_TF,EM1_TF,EM2_TF in TF, ENSC_LF,EM1_LF,EM2_LF in LF
							Pbcscavec(manager->Jpbc(), paramDDscat->Nscat(), pyddx, pzddx, currentTarget->A1(), currentTarget->A2(), 
								oridata->Theta(itheta), oridata->Beta(ibeta), 
								xlr, ylr, zlr, ak1, en0_tf, cxe01_tf, orderm, ordern, paramDDscat->Thetan(), paramDDscat->Phin(), aks_tf, em1_tf, em2_tf, ensc_lf, em1_lf, em2_lf);
						}
					}
//
					Common1::GetInstance()->Ak_tf().Copy(en0_tf);
					Common1::GetInstance()->Ak_tf() *= ak1;
					FileNamer::GetInstance()->Namer(iwav, irad, (paramDDscat->Iwrksc() == 1) ? iori : 1);
// initialize CXXI:
					theDipoleData->ClearCxxi();
//fprintf(stderr, " --- BBB ---\n");
//theDipoleData->Debug(stderr);
					if(nrfld == NearfieldMethodDo)
					{
// read CXPOL1 and CXADIA (in reduced form...)
						for(int iu=0; iu < 3*currentTarget->Nat0(); ++iu)
						{
							cxsc[iu].Read(file22);
						}
						theMatrix->ReadDiagonal(file22);
						if(paramDDscat->Iorth() == 2)
						{
// read CXPOL2
							for(int iu=0; iu < 3*currentTarget->Nat0(); ++iu)
							{
								read(file23, cxscr1+lace+iu, sizeof(Complex));
							}
						}
						Nearfield(currentTarget, paramDDscat, FileNamer::GetInstance()->GetCfle1(), FileNamer::GetInstance()->GetCfle2(), 
							FileNamer::GetInstance()->GetCfleb1(), FileNamer::GetInstance()->GetCfleb2(), Common10::GetInstance()->Myid(), ak1, Common1::GetInstance()->Ak_tf(), 
							manager->Ncomp(), aeff, wave, theMatrix->Diagonal(), dielec, cxsc, cxscr1+lace, cxe01_tf, cxe02_tf, Common4::GetInstance()->Cxzw());
					}
					else
					{
						Getfml(currentTarget, paramDDscat, Common1::GetInstance()->Ak_tf(), ak3, aks_tf, theMatrix, theTensor, theDipoleData, 
							cxe01_tf, cxe02_tf, dielec, cxfData, cxsc, cxscr1, em1_tf, em2_tf, ibeth, ibeth1, iphi, iphi1, itask, itnum, manager->Jpbc(), 
							lace, Common10::GetInstance()->Myid(), navg, sumPackage);
//
// if IORTH=1, GETFML returns reduced polarization array for orientation IPHI=IPHI1 and JO=1 (incident wave CXE01_TF) in CXSCR1, and scattering amplitude factors CXF11,CXF21 for NSCAT directions
// if IORTH=2, GETFML also returns reduced polarization array for JO=2 (incident wave CXE02_TF) in CXSC(LACE) and scattering amplitude factors CXF12,CXF22 
						if (paramDDscat->Iwrksc() == 1)
						{
							outputer->Writefml(iori, irad, iwav, navg, manager->Ncomp(), currentTarget->GetFreeDescription(), aeff, ak1, 
								Common1::GetInstance()->Ak_tf(), wave, xx, cxe01_tf, cxe02_tf, cxfData);
						}
	
						for(j=0; j<2; ++j)
						{
							if(itnum[j] > itnummx[j])
								itnummx[j] = itnum[j];
						}
//
// Compute Mueller scattering matrix
						Getmueller(paramDDscat, ibeta, iphi, itheta, manager->Jpbc(), orderm, ordern, ak1, aks_tf, ensc_lf, ensc_tf, pyddx, pzddx, 
							sm, smori_1, s1111_1, s2121_1, cx1121_1, cxfData, sumPackage, em1_lf, em2_lf, em1_tf, em2_tf);
//
// Choose whether or not to write out scattering properties for current orientation; conditional may be replaced if desired.
// If IWRKSC=1, write scattering properties of current orientation:
						if (paramDDscat->Iwrksc() == 1)
						{
							sumPackage.UseWhatSum() = true;
							outputer->Writesca(itheta, ibeta, iphi, manager->Ioshp(), iori, irad, iwav, navg, itnum, manager->Ncomp(), nori, currentTarget->GetFreeDescription(), aeff, 
								ak1, Common1::GetInstance()->Ak_tf(), wave, xx, sumPackage, s1111_1, s2121_1, sm, smori_1, cx1121, cxe01_tf, cxe02_tf, cxfData, pyddx, pzddx, xMinmax);
						}
// Only write out polarization array if IPHI=1 (for IPHI>1, can reconstruct polarization array by linear combination of CFPOL1 and CFPOL2)
						if ((paramDDscat->Iwrpol() == 1) && (iphi == 0))
						{
							outputer->Writepol(aeff, Common1::GetInstance()->Ak_tf(), wave, cxe01_tf, theMatrix, cxscr1, FileNamer::GetInstance()->GetCflpol1());
							if (paramDDscat->Iorth() == 2)
							{
								outputer->Writepol(aeff, Common1::GetInstance()->Ak_tf(), wave, cxe02_tf, theMatrix, cxsc + lace, FileNamer::GetInstance()->GetCflpol2());
							}
						}
					}							// end loop over IPHI
				}
				iphi1 = 0;
			}									// end loop over IBETH

			if (nrfld != NearfieldMethodDo)
			{
				ibeth1 = 0;
// Collect the partial sums over target orientation from the parallel processes. This routine is a dummy if MPI is not being used.
				Colsum(paramDDscat, Common10::GetInstance()->Myid(), sumPackage, s1111, s1111_1, s2121, s2121_1, cx1121, cx1121_1, smori, smori_1);
//
// Write dielectric information to 'mtable'
				if (!Common10::GetInstance()->Myid())
				{
					outputer->WriteDielec(wave);
				}
//
// Set IORI=0 prior to calling WRITESCA in order to print orientational averages
				iori = -1;
//
// Now print out orientational averages
				if (!Common10::GetInstance()->Myid())
				{
					sumPackage.UseWhatSum() = false;
					outputer->Writesca(itheta, ibeta, iphi, manager->Ioshp(), iori, irad, iwav, navg, itnummx, manager->Ncomp(), nori, currentTarget->GetFreeDescription(),
						aeff, ak1, Common1::GetInstance()->Ak_tf(), wave, xx, sumPackage, s1111, s2121, sm, smori, cx1121, cxe01_tf, cxe02_tf, cxfData, pyddx, pzddx, xMinmax);
				}
			}
		}											// ! end loop over IRAD
		irad1 = 0;
	}												// ! end loop over IWAV
//
// Close binary file (this perhaps should be moved to writesca)
	if (!Common10::GetInstance()->Myid())
	{
//		if (paramDDscat->Cbinflag() == BinflagMethod_NOTBIN) 
//		{
//			Wrimsg("ddscat", " close dd.bin");
//		}
	}
	Wrimsg(DDscatLabel, "Return from DDSCAT to calling program");
//
// Deallocate all
	delete [] em1_lf;
	delete [] em1_tf;
	delete [] em2_lf;
	delete [] em2_tf;
	delete [] ensc_lf;
	delete [] ensc_tf;
	delete [] aks_tf;
	delete [] orderm;
	delete [] ordern;
	delete [] s1111;
	delete [] s1111_1;
	delete [] s2121;
	delete [] s2121_1;
	delete [] cxsc;
	delete cxfData;
	delete [] cx1121;
	delete [] cx1121_1;
	delete [] sm;
	delete [] smori;
	delete [] smori_1;
	free(scrrs2);

	delete theMatrix;
	delete theTensor;
	delete theDipoleData;
	delete [] cxscr1;
//
// MPI environment shutdown:
	ierr = MPI_Finalize();

	DDscatParameters::Kill();
	OutputManager::Kill();
	DielectricManager::Kill();
	GreenFunctionManager::Kill();
	FileNamer::Kill();
}

#endif
