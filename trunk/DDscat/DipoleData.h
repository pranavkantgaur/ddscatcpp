#ifndef __DIPOLEDATA_H__
#define __DIPOLEDATA_H__

#include "Complex.h"

class DipoleData
{
protected:
	Complex *cxe_tf, *cxxi;
	int curSize;

public:
	DipoleData();
	~DipoleData();

public:
	inline Complex &Cxe_tf(int index) { return cxe_tf[index]; }
	inline Complex &Cxxi  (int index) { return cxxi[index]; }

	inline Complex *Cxe_tf() { return cxe_tf; }
	inline Complex *Cxxi  () { return cxxi; }

public:
	void Allocate(int size);
	void ClearCxxi();
	void ClearCxeTf();
	void Debug(FILE *file, int from = 0, int toto = 0);

protected:
	void Delete();
};

#endif // __SCRATCHARRAYS_H__