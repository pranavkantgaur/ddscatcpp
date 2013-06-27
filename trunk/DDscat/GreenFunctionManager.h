#ifndef __GREENFUNCTIONMANAGER_H__
#define __GREENFUNCTIONMANAGER_H__

#include "Enumerator.h"
#include "Complex.h"
#include "ArrayF.h"
#include "Vect3.h"

typedef enum _SubsystemType
{
	SubsystemTypeElectric, SubsystemTypeMagnetic, SubsystemTypeEnd
} SubsystemType;

class GreenFunctionManager
{
	class Subsystem
	{
	protected:
		const int mySize;
		int nx, ny, nz, *issym;
		Array4Stacked<Complex> *cxzcg;
		Complex *dcxsum, *tmp6;
		bool ipbc;
		SubsystemType myType;
		Subsystem(void);

	public:
		Subsystem(int mSize, SubsystemType mt);
		~Subsystem(void);

	public:
		inline bool &Ipbc(void) { return ipbc; }
		inline Array4Stacked<Complex> *Cxzcg(void) { return cxzcg; }

	public:
		void SetDimension(int nx, int ny, int nz);
		bool Allocate(void);
		void Deallocate(void);
		void PrintTimes(const char *label, const char *prefix, real t2);
		void PrintTimeForecasts(const char *label, real py, real pz);
		void Extnd(Complex *cxa, int *isym, Complex *cxb);
		void Trim(Complex *cxb, Complex *cxa);
		void PadC(Complex *cxa, unsigned int nx, unsigned int ny, unsigned int nz, unsigned int m, Complex *cxb);
		void ImportCorner(Complex *to, Complex *from, unsigned int ax, unsigned int ay, unsigned int az, unsigned int at, unsigned int m);
		void ExportCorner(Complex *from, Complex *to, unsigned int ax, unsigned int ay, unsigned int az, unsigned int at, unsigned int m, const Complex &scale);
		void DirectCalc(int ix, int iy, int iz, Vect3<real> &dx, Vect3<real> &ak, real akd, real akd2, real gamma, real pyddx, real pzddx);
		void Self(FftMethod cmethd, Complex *cxzp, real gamma, real pyd, real pzd, Vect3<real> &ak, real akd, Vect3<real> &dx, Array4Stacked<Complex> *cxzw, Complex *cxze);

	protected:
		virtual void LoadIssym(void) = 0;
		virtual void PrepareTmp3(Complex *tmp3, int ign, int jgn, int kgn, unsigned int *pos, Complex *cxzwData) = 0;
		virtual void PrepareTmp3(Complex *tmp3, unsigned int *pos, Complex *cxzwData) = 0;
		virtual void DeepDirectCalc(real r2, real phasyz, real akd, real akd2, real gammakd4, real *x) = 0;
	};

	class SubsystemElectric : public Subsystem
	{
	public:
		SubsystemElectric(void) : Subsystem(6, SubsystemTypeElectric) { }
		~SubsystemElectric(void) { }

	protected:
		void DeepDirectCalc(real r2, real phasyz, real akd, real akd2, real gammakd4, real *x);
		void Cisi(real x, real &ci, real &si);
		void LoadIssym(void);
		void PrepareTmp3(Complex *tmp3, int ign, int jgn, int kgn, unsigned int *pos, Complex *cxzwData);
		void PrepareTmp3(Complex *tmp3, unsigned int *pos, Complex *cxzwData);
	};

	class SubsystemMagnetic : public Subsystem
	{
	public:
		SubsystemMagnetic(void) : Subsystem(3, SubsystemTypeMagnetic) { }
		~SubsystemMagnetic(void) { }

	protected:
		void DeepDirectCalc(real r2, real phasyz, real akd, real akd2, real gammakd4, real *x);
		void LoadIssym(void);
		void PrepareTmp3(Complex *tmp3, int ign, int jgn, int kgn, unsigned int *pos, Complex *cxzwData);
		void PrepareTmp3(Complex *tmp3, unsigned int *pos, Complex *cxzwData);
	};

protected:
	static GreenFunctionManager *item;
	GreenFunctionManager(void);
	virtual ~GreenFunctionManager(void);
	Subsystem *subSys[SubsystemTypeEnd];

public:
	static GreenFunctionManager *GetInstance();
	static void Kill(void);
	void Init(void);

public:
	Subsystem *GetSusystem(SubsystemType type);
	Subsystem *GetElectric(void);
	Subsystem *GetMagnetic(void);
	void SetDimension(int nx, int ny, int nz);
	bool AllocateElectricBuffer(void);
	void DeallocateElectricBuffer(void);
	bool AllocateMagneticBuffer(void);
	void DeallocateMagneticBuffer(void);
	void SetIpbc(bool value);
	bool &Ipbc(void);
};

#endif // __GREENFUNCTIONMANAGER_H__