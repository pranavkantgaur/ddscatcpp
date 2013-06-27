#ifndef __TARGSPHER_H__
#define __TARGSPHER_H__

#include "LoadableTarget.h"
#include "Targetlib.h"
#include "ArrayF.h"

class Targspher : public LoadableTarget
{
protected:
	Array2F<Complex> alm;
	Vect3<real> eigval;
	real s2lmax, f1;
	int lmax, nres;
	Targspher(void) {  nres = 20;   lmax = nint_(shpar[3]);  }

public:
	Targspher(TargetManager *man);
	virtual ~Targspher(void) 
	{
		alm.Deallocate();
	}

public:
	virtual void Sizer(void);
	virtual void Descriptor(void);
	virtual void Vector(void);
	virtual void VectorX(void);
	virtual void Allocator(void);
	virtual void Printer(void);

protected:
	real Gasdev(int idum);
	void Presizer(void);
	bool AlmReader(void);
	real Rgspher(real costh, real phi, int lmax, int nmax, real s2, Array2F<Complex> &alm);
	real P_lm(int l, int m, real x);
};

class Target_Gausssph : public Targspher 
{ 
protected: 
	Target_Gausssph(void) {}

public: 
	Target_Gausssph(TargetManager *man) : Targspher(man) {} 
	virtual ~Target_Gausssph(void) {} 
	virtual void SayHello(FILE *stream);
}; 

#endif // __TARGSPHER_H__