#ifndef __TARCYLCAP_H__
#define __TARCYLCAP_H__

#include "AbstractTarget.h"
#include "Targetlib.h"

class Tarcylcap : public AbstractTarget
{
protected:
	real xcm, ycm, zcm;
	int jxu, jxl;
	Tarcylcap(void) { }

public:
	Tarcylcap(TargetManager *man);
	virtual ~Tarcylcap(void) { }

public:
	virtual void Sizer(void);
	virtual void Descriptor(void);
	virtual void Vector(void);
	virtual void VectorX(void);
	virtual void Allocator(void);
	virtual void OutShpar(FILE *file);
};

class Target_Cylndrcap : public Tarcylcap 
{ 
protected: 
	Target_Cylndrcap(void) {}

public:
	Target_Cylndrcap(TargetManager *man) : Tarcylcap(man) {}
	virtual ~Target_Cylndrcap(void) {}
	virtual void SayHello(FILE *stream);
};

#endif // __TARCYLCAP_H__