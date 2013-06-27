#ifndef __MATRIX_H__
#define __MATRIX_H__

#include "Vect3.h"
#include "Complex.h"
#include "DipoleData.h"
#include "LoadableTarget.h"

//
// The data in cxadia, cxaoff is placed according to C rules!
class Matrix
{
protected:
	Complex *cxadia, *cxaoff;
	unsigned int curSize;

public:
	Matrix(void);
	virtual ~Matrix(void);
	void Allocate(unsigned int newSize);

public:
	inline Complex *Diagonal() { return cxadia; }
	inline Complex *OffDiagonal() { return cxaoff; }
	inline Complex &Diagonal(int index) { return cxadia[index]; }
	inline Complex &OffDiagonal(int index) { return cxaoff[index]; }
	inline unsigned int GetSize() { return curSize; }

public:
	int WriteDiagonal(int file);
	int WriteOffDiagonal(int file);
	int ReadDiagonal(int file);
	int ReadOffDiagonal(int file);
	void Debug(FILE *file, int from = 0, int toto = 0);
	void Evala(Matrix *theTensor, int nat);
	void Evalq(DipoleData *theDipoleData, Vect3<real> &ak, int nat3, real e02, real &cabs, real &cext, real &cpha, int imethd);
	void AbstractAlphadiag(AbstractTarget *currentTarget, Complex *dielec, real b1, real b2, real *b3, const Complex &cxrr);
	void LoadableAlphadiag(LoadableTarget *loadableTarget);

protected:
	void Delete();
};

#endif // __MATRIX_H__