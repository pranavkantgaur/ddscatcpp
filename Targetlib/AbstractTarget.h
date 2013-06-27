#ifndef __ABSTRACTTARGET_H__
#define __ABSTRACTTARGET_H__

#include <vector>
#include "Targetlib.h"
#include "Definitions.h"
#include "Enumerator.h"
#include "TargetDefinitions.h"
#include "Vect3.h"
#include "Vect6.h"
#include "Matc3.h"
#include "ArrayF.h"
#include "Complex.h"

extern const real quat_;
extern const real half_;
extern const real zero_;
extern const real onex_;
extern const real twox_;
extern const real four_;

#define REGISTER_TARGET(x,y,z,aa,bb,cc) \
AbstractTarget *Create##x() \
{ \
	AbstractTarget *result = new Target_##x(TargetManager::GetInstance()); \
	result->Build(); \
	return result; \
} \
const bool reg##x = TargetManager::GetInstance()->RegisterTarget(TargetType_##x, #x, Create##x, y, z, aa, bb, cc, TargetVerboseDescriptor_##x);

class TargetManager;

class TARGETLIB_API AbstractTarget
{
protected:
	TargetManager *manager;
	AbstractTarget(void);
	string shortDescr, longDescr;
	char freeDescr[256];
	Vect3<real> a1, a2, dx, x0;
	Array2F<int> ixyz;
	Array2F<short> icomp;
	int ncomp;
	real *shpar;
	real pyd, pzd;
	bool *iocc;
	IsotropicFlag ianiso; 
	int nat, nat0, nx, ny, nz;
	int minJx, maxJx, minJy, maxJy, minJz, maxJz;

public:
	AbstractTarget(TargetManager *man);
	virtual ~AbstractTarget(void);
	void Init(void);
	void OutVectorsAndArrays(FILE *file);
	void Debug(const char *fileName);
	Vect3<real> Prinaxis(void);
	void Dsyevj3(Matc3<real> &a, Matc3<real> &q, Vect3<real> &w);
	void AllocateArrays(int sizeX, int sizeY, int sizeZ, bool bWithIxyz = false);
	Vect6<int> IxyzMinmax(void);
	int GetLinearAddress(int index);
	int GetRelativeLinearAddress(int ax, int ay, int az, int mX = 0, int mY = 0, int mZ = 0);
	real DipoleScalar(const Vect3<real> &op, int absIndex);
	void DipoleVector(int absIndex, real *res);
	void SetMin(const Vect6<int> &op);
	bool IsDipolePhysicallyAnisotropic(unsigned int index);

public:
	inline Vect3<real> &A1(void) { return a1; }
	inline Vect3<real> &A2(void) { return a2; }
	inline Vect3<real> &Dx(void) { return dx; }
	inline Vect3<real> &X0(void) { return x0; }
	inline Array2F<int> &Ixyz(void) { return ixyz; }
	inline Array2F<short> &Icomp(void) { return icomp; }
	inline int &Ncomp(void) { return ncomp; }
	inline bool *Iocc(void) { return iocc; }
	inline bool &Iocc(int index) { return iocc[index]; }
	inline IsotropicFlag &Ianiso(void) { return ianiso; }
	inline int &Nat(void) { return nat; }
	inline int &Nat0(void) { return nat0; }
	inline int &Nx(void) { return nx; }
	inline int &Ny(void) { return ny; }
	inline int &Nz(void) { return nz; }
	inline char *GetFreeDescription() { return freeDescr; }
	inline real &Pyd(void) { return pyd; }
	inline real &Pzd(void) { return pzd; }
	inline int GetMinJx(void) { return minJx; }
	inline int GetMinJy(void) { return minJy; }
	inline int GetMinJz(void) { return minJz; }

public:
	virtual bool IsLoadableTarget(void) { return false; }
	virtual bool IsMultiBlockTarget(void) { return false; }
	virtual void Build(void);
	virtual void Sizer(void) { }
	virtual void Vector(void);
	virtual void Descriptor(void) { }
	virtual void Allocator(void) { }
	virtual void Printer(void);
	virtual void Composer(int index, int item);
	virtual void SayHello(void);
	virtual void OutShpar(FILE *file);
	virtual void PreparePyzd(void);
	virtual void PrepareIaniso(void);

public:
	void Evale(Vect3<Complex> &cxe00, Vect3<real> &akd, Complex *cxe);
	void Reduce(Complex *cxv);
	void Unreduce(Complex *cxv);

protected:
	virtual void VectorA(void);
	virtual void VectorX(void) { }
	virtual void ShiftDipolesAndX(void);
	int Analyze(int &number);
	void InternalMinMax(int x, int y, int z);
	virtual void PrepareIanisoSpecial(void);
	void IsotropicComposer(int index, int item);
	void AnisotropicComposer(int index, int item);
};

#endif
