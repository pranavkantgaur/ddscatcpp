#ifndef __FOURARRAY_H__
#define __FOURARRAY_H__

#include "General.h"
#include "Complex.h"

class GENERALLIB_API FourArray
{
protected:
	Complex *cxs1, *cxs2, *cxs3, *cxs4;
	int curSize;

public:
	FourArray();
	~FourArray();
	void Alloc(int newSize);

public:
	inline Complex &Cxs1(int index) { return cxs1[index]; }
	inline Complex &Cxs2(int index) { return cxs2[index]; }
	inline Complex &Cxs3(int index) { return cxs3[index]; }
	inline Complex &Cxs4(int index) { return cxs4[index]; }
	inline Complex &Cx11(int index) { return cxs1[index]; }
	inline Complex &Cx12(int index) { return cxs2[index]; }
	inline Complex &Cx21(int index) { return cxs3[index]; }
	inline Complex &Cx22(int index) { return cxs4[index]; }

	inline Complex *Cx11() { return cxs1; }
	inline Complex *Cx12() { return cxs2; }
	inline Complex *Cx21() { return cxs3; }
	inline Complex *Cx22() { return cxs4; }
	inline Complex *Cxs1() { return cxs1; }
	inline Complex *Cxs2() { return cxs2; }
	inline Complex *Cxs3() { return cxs3; }
	inline Complex *Cxs4() { return cxs4; }

protected:
	void Delete();
};

#endif