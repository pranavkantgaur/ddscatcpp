#pragma once

#include "ArrayF.h"
#include "Vect3.h"
#include "Complex.h"

class ReadnfHolder
{
protected:
	Complex *cxbsca, *cxeps, *cxpol, *cxesca, *cxadia, *cxeinc, *cxbinc;
	Array2F<short> *icomp;
	int nrfldb, nx, ny, nz, nxyz, nat0, nat3, ncomp, nambient;
	real xmin, xmax, ymin, ymax, zmin, zmax;
	real aeff, wave, dphys;
	Vect3<real> x0, ak_tf;
	Vect3<Complex> cxe0_tf, cxb0_tf;

public:
	ReadnfHolder(void);
	virtual ~ReadnfHolder(void);

public:
	bool Readnf(const char *cflename);
	void AllocateEArrays();
	void AllocateBArrays();
	bool IsPointInside(real xa, real xb, real ya, real yb, real za, real zb);
	void Output(FILE *outputFile, const Vect3<real> &s_inc, real snorm);
	void Evaluate(FILE *outputFile, real xa, real xb, real ya, real yb, real za, real zb, int nab, Array3F <Vect3<real> > &s);

public:
	inline int GetNrfldb() { return nrfldb; }
	inline int GetNx() { return nx; }
	inline int GetNy() { return ny; }
	inline int GetNz() { return nz; }
	inline int GetNxyz() { return nxyz; }
	inline int GetNcomp() { return ncomp; }
	inline Complex GetCxeinc(int index) { return cxeinc[index]; }
	inline Complex GetCxbinc(int index) { return cxbinc[index]; }
	inline Complex GetCxesca(int index) { return cxesca[index]; }
	inline Complex GetCxbsca(int index) { return cxbsca[index]; }
	inline Vect3<Complex> &Cxe0_tf() { return cxe0_tf; }
	inline Vect3<Complex> &Cxb0_tf() { return cxb0_tf; }
	inline Vect3<real> &GetX0() { return x0; }
	inline Vect3<real> &GetAk_tf() { return ak_tf; }
	inline real GetXmin() { return xmin; }
	inline real GetXmax() { return xmax; }
	inline real GetYmin() { return ymin; }
	inline real GetYmax() { return ymax; }
	inline real GetZmin() { return zmin; }
	inline real GetZmax() { return zmax; }
	inline real GetDphys() { return dphys; }
	inline short GetIcomp(int jx, int jy) { return icomp->Value(jx, jy); }
};
