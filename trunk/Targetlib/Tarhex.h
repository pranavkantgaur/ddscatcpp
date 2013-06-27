#ifndef __TARHEX_H__
#define __TARHEX_H__

#include "AbstractTarget.h"
#include "Targetlib.h"

class Tarhex : public AbstractTarget
{
protected:
	int iori;
	real acm, bcm, ccm;

	Tarhex(void) { }

public:
	Tarhex(TargetManager *man);
	virtual ~Tarhex(void) { }

public:
	virtual void Sizer(void);
	virtual void Descriptor(void);
	virtual void Vector(void);
	virtual void VectorX(void);
	virtual void Allocator(void);
	virtual void OutShpar(FILE *file);

protected:
	bool Testhex(real x, real y, real z);
};

class Target_Hexgonpbc : public Tarhex 
{ 
protected: 
	Target_Hexgonpbc(void) {}

public: 
	Target_Hexgonpbc(TargetManager *man) : Tarhex(man) {} 
	virtual ~Target_Hexgonpbc(void) {}
	virtual void SayHello(FILE *stream);
}; 

class Target_HexPrism : public Tarhex 
{ 
protected: 
	Target_HexPrism(void) {}

public: 
	Target_HexPrism(TargetManager *man) : Tarhex(man) {} 
	virtual ~Target_HexPrism(void) {} 
	virtual void SayHello(FILE *stream);
}; 

#endif // __TARHEX_H__