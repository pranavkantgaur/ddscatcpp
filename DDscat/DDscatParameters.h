#ifndef __DDSCATPARAMETERS_H__
#define __DDSCATPARAMETERS_H__

#include "Vect3.h"
#include "Vect6.h"
#include "Complex.h"
#include "TargetDefinitions.h"
#include "Definitions.h"
#include "AbstractFftEngine.h"
#include "Matrix.h"
#include "ArrayF.h"
#include "XmlHelper.h"
#include "Functions.h"

//
// Singleton class to keep program parameters read from initialization file
//
// Original version created by Choliy V., Kyiv Shevchenko University 
// Conforms to fotran version of DDscat 7.3.0 by 
// P.J.Flatau, Colorado State Univ. and B.T.Draine, Princeton Univ. Obs.
//
// History:
// 12.02.26 (ChB): Initial version.
// end history
//
// Copyright (C) 2012,2013 C++ versions, Choliy V.
// Copyright (C) All fortran versions of DDscat, 1993-2013, B.T. Draine and P.J. Flatau
// This code is covered by the GNU General Public License.
//

class OriData
{
protected:
	Array2F<real> wgta;
	real betami, betamx, thetmi, thetmx, phimin, phimax;
	int nbeta, ntheta, nphi;
	real *beta, *theta, *phi, *wgtb;

public:
	OriData() { beta = theta = phi = wgtb = NULL; }
	~OriData() 
	{
		if (beta) free(beta);
		if (theta) free(theta);
		if (phi) free(phi);
		if (wgtb) free(wgtb);
	}

public:
	inline real &Betami() { return betami; }
	inline real &Betamx() { return betamx; }
	inline real &Thetmi() { return thetmi; }
	inline real &Thetmx() { return thetmx; }
	inline real &Phimin() { return phimin; }
	inline real &Phimax() { return phimax; }
	inline int GetNbeta() { return nbeta; }
	inline int GetNtheta() { return ntheta; }
	inline int GetNphi() { return nphi; }
	inline real *Beta() { return beta; }
	inline real &Beta(int index) { return beta[index]; }
	inline real *Theta() { return theta; }
	inline real &Theta(int index) { return theta[index]; }
	inline real *Phi() { return phi; }
	inline real &Phi(int index) { return phi[index]; }
	inline real GetWgtb(int index) { return wgtb[index]; }
	inline real GetWgtaValue(int r, int c) { return wgta.Value(r, c); }

public:
	void Orient();
	void SetDimensions(int nb, int nt, int np);
	void SetBetas(real mi, real mx, int n)
	{
		betami = mi; betamx = mx; nbeta = n;
	}
	void SetThetas(real mi, real mx, int n)
	{
		thetmi = mi; thetmx = mx; ntheta = n;
	}
	void SetPhis(real mi, real mx, int n)
	{
		phimin = mi; phimax = mx; nphi = n;
	}
	void ToRadians()
	{
		betami /= Degrad;
		betamx /= Degrad;
		thetmi /= Degrad;
		thetmx /= Degrad;
		phimin /= Degrad;
		phimax /= Degrad;
	}
	void ToDegrees()
	{
		betami *= Degrad;
		betamx *= Degrad;
		thetmi *= Degrad;
		thetmx *= Degrad;
		phimin *= Degrad;
		phimax *= Degrad;
	}
};

class Reader
{
protected:
	FILE *file;
	int lineCounter, lineSize;
	char buffer[256];

public:
	Reader() { lineCounter = 0; file = NULL; }
	~Reader() { }
	bool Open(const char *fileName)
	{
		lineCounter = 0;
		file = fopen(fileName, "r");
		return (file == NULL) ? false : true;
	}
	void Close(void)
	{
		fclose(file);
		file = NULL;
	}
	char *Read(const char *Label = NULL)					// Read the next line skipping comment lines, "\' ", "\'***", "!"
	{
		while(1)
		{
			char *ia = ReadNextLine(Label);
			if (!ia)
				return NULL;
			if (buffer[0] == '!')
				continue;
			if (!memcmp(buffer, "\' ", 2))
				continue;
			if (!memcmp(buffer, "\'***", 4))
				continue;
			return ia;
		}
	}
	char *ReadNextLine(const char *Label = NULL)			// Read the next line
	{
		++lineCounter;
		char *ia = fgets(buffer, 255, file);
		if (feof(file))
		{
			buffer[0] = '\0';
			return NULL;
		}
		lineSize = (int)(strlen(buffer) - 1);
		while(lineSize)
		{
			if ((buffer[lineSize] == 0x0d) || (buffer[lineSize] == 0x0a) || (buffer[lineSize] == ' '))
			{
				buffer[lineSize] = '\0';
				lineSize--;
			}
			else
				break;
		}
		if (Label)
			fprintf(stdout, " >%s %s\n", Label, buffer); 
		return ia;
	}
	int ExtractFirstWord(char fc)
	{
		char *ia = NULL;
		while(buffer[0] == fc)
			ShiftLeft();
		ia = strchr(buffer, fc);
		if (ia)
			*ia = '\0';
		strupr(buffer);
		return (ia - buffer);
	}
	void RemoveSymbols(char sym, char toSymbol)
	{
		int i = 0;
		while((buffer[i] != toSymbol) && (buffer[i] != '\0'))
		{
			if (buffer[i] == sym)
			{
				ShiftLeft(i);
				continue;
			}
			else
				++i;
		}
	}
	int ScanInt(void)
	{
		int res;
		sscanf(buffer, "%d", &res);
		return res;
	}
	real ScanReal(const char *format)
	{
		real res;
		sscanf(buffer, format, &res);
		return res;
	}
	char GetFirstNonemptyCharacter(void)
	{
		for(int i=0; i<lineSize; ++i)
		{
			if (buffer[i] == ' ')
				continue;
			return buffer[i];
		}
		return 0x0a;
	}

protected:
	void ShiftLeft(int shift = 0)
	{
		char *ia = buffer + shift;
		do
		{
			*ia = *(ia + 1);
			++ia;
		} while(*ia);
	}

public:
	inline char *GetBuffer(void) { return buffer; }
	inline int GetLineCounter(void) { return lineCounter; }
};

class DDscatParameters
{
protected:
	static DDscatParameters *item;
	Reader reader;

public:
	static DDscatParameters *GetInstance();
	static void Kill();

protected:
	DDscatParameters(void);
	virtual ~DDscatParameters(void);
	void Init(void);
	void Destroy(void);
	void OutParam(const char *Label);
	bool Divide(CdividMethod wDivid, real x1, real x2, int nx, real *ary);
	char *dgets(char *Buffer, int num, FILE *file);

protected:
	OriData *oridata;
    real shpar[12], tol, gamma, etasca, nambient;
	real *thetan, *phin;
	real *aeffa, *wavea;
	int *smind1, *smind2;
	int numShpar, mxiter, nwav, iorth, iwrksc, iwrpol, nscat;
	char *cflshp;
	int nrad, iwav0, irad0, iori0, nsmelts;
	NearfieldBMethod nrfldb;
	Vect3<real> dx;
	Vect6<real> extendxyz;
	Vect3<Complex> cxe01_lf, cxe02_lf;
	FftMethod cmdfft;
	SolMethod cmdsol;
	TorqMethod cmdtrq;
	AlphaMethod calpha; 
	BinflagMethod cbinflag;
	char targetString[16];
	TargetType cshape;
	FrameCode cmdfrm;
	CdividMethod wavDivid, aefDivid;
	AbstractFftEngine *fftEngine;
	Matrix *theMatrix;

public:
	int Load(const char *cflpar, NearfieldMethod &nrfld);
	int LoadXml(const char *fileName, NearfieldMethod &nrfld);
	void Orient(void);
	void Scavec(Vect3<real> *ensc, Vect3<real> *em1, Vect3<real> *em2);

public:
	inline Matrix *GetMatrix() { return theMatrix; }
	inline void SetMatrix(Matrix *op) { theMatrix = op; }
	inline OriData *GetOriData() { return oridata; }
	inline TorqMethod Cmdtrq() { return cmdtrq; }
	inline SolMethod Cmdsol() { return cmdsol; }
	inline FftMethod Cmdfft() { return cmdfft; }
	inline AbstractFftEngine *GetFftEngine() { return fftEngine; }
	inline AlphaMethod Calpha() { return calpha; }
	inline BinflagMethod Cbinflag() { return cbinflag; }
	inline TargetType Cshape() { return cshape; }
	inline char *TargetString() { return targetString; }
	inline char *Cflshp() { return cflshp; }
	inline FrameCode Cmdfrm() { return cmdfrm; }
	inline int &Mxiter() { return mxiter; }
	inline int &Nwav() { return nwav; }
	inline int &Nrad() { return nrad; }
	inline int &Nscat() { return nscat; }
	inline int &Iorth() { return iorth; }
	inline int &Iwrksc() { return iwrksc; }
	inline int &Iwrpol() { return iwrpol; }
	inline int &Iwav0() { return iwav0; }
	inline int &Irad0() { return irad0; }
	inline int &Iori0() { return iori0; }
	inline int &Nsmelts() { return nsmelts; }
	inline NearfieldBMethod &Nrfldb() { return nrfldb; }
	inline int NumShpar() { return numShpar; }
	inline int *Smind1() { return smind1; }
	inline int *Smind2() { return smind2; }
	inline real *Shpar() { return shpar; }
	inline real *Aeffa() { return aeffa; }
	inline real &Tol() { return tol; }
	inline real *Thetan() { return thetan; }
	inline real &Gamma() { return gamma; }
	inline real &Etasca() { return etasca; }
	inline real &Nambient() { return nambient; }
	inline real *Phin() { return phin; }
	inline real *Wavea() { return wavea; }
	inline real &Wavea(int index) { return wavea[index]; }
	inline Vect3<Complex> &Cxe01_lf() { return cxe01_lf; }
	inline Vect3<Complex> &Cxe02_lf() { return cxe02_lf; }
	inline Vect6<real> &Extendxyz() { return extendxyz; }
	inline Vect3<real> &Dx() { return dx; }
	
protected:
	void Clear2(void);
	bool ProcessCmtorq(char *Buffer, const char *Label);
	bool ProcessCmdsol(char *Buffer, const char *Label);
	bool ProcessCmdFFT(char *Buffer, const char *Label);
	bool ProcessCalpha(char *Buffer, const char *Label);
	bool ProcessCbinflag(char *Buffer, const char *Label);
	bool ProcessTarget(char *Buffer, const char *Label);
	bool ProcessCmdfrm(char *Buffer, const char *Label);
	int LoadShparFromXml(xmlNodePtr nod1);
};

/* **
Subroutine REAPAR handles the reading of input parameters from the
"ddscat.par" file, as well as elementary processing with those input
parameters to generate arrays.

Given:
      CFLPAR = name of file (normally 'ddscat.par')
      NRFLD    = 0 on first call
               = 1 if previous call set NRFLD=1 to prepare to do nearfield calculation
Returns:
      CMDTRQ   = 'NOTORQ' or 'DOTORQ'
      CMDSOL   = 'PBCGST' or 'PBCGS2' or 'PETRKP' or 'GPBICG' or 'QMRCCG'
      CMDFFT   = 'GPFAFT' or 'FFTW21' or 'FFTMKL'
      CALPHA   = 'LATTDR' or 'GKDLDR' or 'FLTRCD'
      CBINFLAG = 'ALLBIN' or 'ORIBIN' or 'NOTBIN'
      CSHAPE   = 'RCTGLPRSM' or 'ELLIPSOID' or ... other valid shape
      JPBC     = 0 :isolated finite targets
               = 1 :target periodic in y direction, finite in z
               = 2 :target finite in y direction, periodic in z
               = 3 :target periodic in both y and z directions
      NCOMP    = number of compositions
      CFLEPS(1-NCOMP) = filenames for dielectric functions
      NRFLD   = 0 : skip near-field calculation
              = 1 : to prepare to do additional near-field calculation
              = 2 : to perform near-field calculation itself
      IDIPINT = 0 for point dipole interaction method
			  = 1 for filtered coupled dipole (FCD) interaction method
      NRFLDB   = 0 : to skip near-field calculation of B
               = 1 : to perform near-field calculation of B
      EXTNDXYZ(1-6) = fractional expansion of calculational volume
                 in -x,+x,-y,+y,-z,+z directions (used only if NRFLD=1) 
      INIT     = 0 (no longer used)
      TOL      = error tolerance
      GAMMA    = parameter controlling integration limit for PBC
      ETASCA   = parameter controlling accuracy of calculation of <cos(theta)>
                 1 is OK, 0.5 gives high accuracy
      NWAV     = number of wavelengths
      WAVEA(1-NWAV) = wavelengths (in vacuo, physical units)
      NAMBIENT = (real) refractive index of ambient medium
      NRAD     = number of radii
      AEFFA(1-NRAD) = target effective radii (physical units)
      CXE01_LF(1-3)= (Complex) polarization state 1
      CXE02_LF(1-3)= Complex polarization state 2 (orthogonal to 1)
      IORTH    = 1 or 2 (number of incident polarization states to use
      IWRKSC   = 0 or 1 (not write/write ".sca" files)
      IWRPOL   = 0 or 1 (not write/write ".pol" files)
      NBETA    = number of beta values for target rotation
      BETAMI   = minimum beta value (rad) [input from file in deg]
      BETAMX   = maximum beta value (rad) [input from file in deg]
      NTHETA   = number of theta values for target rotation
      THETAMI  = minimum theta value (rad) [input from file in deg]
      THETAMX  = maximum theta value (rad) [input from file in deg]
      NPHI     = number of PHI values for target rotation
      PHIMIN   = minimum PHI value (rad) [input from file in deg]
      PHIMAX   = maximum PHI value (rad) [input from file in deg]
      NSMELTS  = number of elements of S matrix to calculate (less than or equal to 9)
      SMINDI(J) = indices of S matrix elements to be calculated (e.g., 11 , 21 , 31 , 41 , 12 , 13)
      CMDFRM   = 'LFRAME' or 'TFRAME' (Lab frame or Target frame)
      NSCA     = number of scattering directions
If JPBC=0:
      THETAN(1-NSCA) = theta for scattering directions (rad)
      PHIN(1-NSCA) = phi for scattering directions (rad)
If JPBC=1 or 2:
      THETAN(1-NSCA) = order_y or order_z for scattering directions
      PHIN(1-NSCA)   = psi (radians) for scattering directions, where psi = rotation around axis of
                       target periodicity (y_TF if JPBC=1, z_TF if JPB
If JPBC=3: periodic in both y and z directions
      THETAN(1-NSCA) = order_y for diffraction
      PHIN(1-NSCA)   = order_z for diffraction
** */

#endif // __DDSCATPARAMETERS_H__