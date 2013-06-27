#ifndef __PROCESSHELPER_H__
#define __PROCESSHELPER_H__

#include "Processlib.h"
#include "ArrayF.h"
#include "Vect3.h"
#include "Complex.h"
#include "Vtr.h"
#include "Linia.h"
#include "Enumerator.h"

#include <string>
using namespace std;

class PROCESSLIB_API ProcessHelper
{
protected:
	Vect3<Complex> cxe0_tf, cxb0_tf;
	Vect3<real> x0, ak_tf, s_inc;
	Array2F<short> icomp;
	string stringVersion;
	Complex *cxbsca, *cxeps, *cxpol, *cxesca, *cxadia, *cxeinc, *cxbinc;
	real aeff, wave, nambient, snorm, *s;
	real dphys, xmin, xmax, ymin, ymax, zmin, zmax, sumerr2;
	int nrword_nf, nxyz, nat0, nat3, ncomp, nx, ny, nz, intVersion;
	NearfieldBMethod nrfldb;

public:
	ProcessHelper(void);
	virtual ~ProcessHelper(void);

public:
	inline int GetNrwordNf(void) { return nrword_nf; }
	inline const char *GetStringVersion(void) { return stringVersion.c_str(); }
	inline int GetIntVersion(void) { return intVersion; }
	inline real GetNambient(void) { return nambient; }
	inline real GetSnorm(void) { return snorm; }
	inline real GetSumerr2(void) { return sumerr2; }
	inline real GetDphys(void) { return dphys; }
	inline real GetAeff(void) { return aeff; }
	inline real GetWave(void) { return wave; }
	inline int GetNat0(void) { return nat0; }
	inline int GetNx(void) { return nx; }
	inline int GetNy(void) { return ny; }
	inline int GetNz(void) { return nz; }
	inline int GetNcomp(void) { return ncomp; }
	inline real GetMinX(void) { return xmin; }
	inline real GetMaxX(void) { return xmax; }
	inline real GetMinY(void) { return ymin; }
	inline real GetMaxY(void) { return ymax; }
	inline real GetMinZ(void) { return zmin; }
	inline real GetMaxZ(void) { return zmax; }

public:
	bool Load(const char *cflename, int nrword);
	void PrepareNormsAndPointing();
	int SanityCheck(void);
	bool PrepareVTKFile(const char *filename, int ivtr);
	bool IsLiniaOk(const Linia &op);
	bool WriteLiniaData(const char *filename, const Linia &op);
	void OutMinMax(void);
};

#endif