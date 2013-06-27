#ifndef __TARPRSM_H__
#define __TARPRSM_H__

#include "AbstractTarget.h"
#include "Targetlib.h"

class Tarprsm : public AbstractTarget
{
protected:
	real cotbeta, cotgamma;
	Tarprsm(void) { }

public:
	Tarprsm(TargetManager *man);
	virtual ~Tarprsm(void) { }

public:
	virtual void Sizer(void);
	virtual void Descriptor(void);
	virtual void Vector(void);
	virtual void VectorX(void);
	virtual void Allocator(void);
	virtual void OutShpar(FILE *file);
};

class Target_Trnglprsm : public Tarprsm 
{ 
protected: 
	Target_Trnglprsm(void) {}

public: 
	Target_Trnglprsm(TargetManager *man) : Tarprsm(man) {} 
	virtual ~Target_Trnglprsm(void) {}
	virtual void SayHello(FILE *stream);
};

#endif // __TARPRSM_H__