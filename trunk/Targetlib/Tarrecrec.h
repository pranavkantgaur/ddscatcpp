#ifndef __TARRECREC_H__
#define __TARRECREC_H__

#include "AbstractTarget.h"
#include "Targetlib.h"

class Tarrecrec : public AbstractTarget
{
protected:
	int jx1min, jx1max, jx2min, jx2max;
	int jy1min, jy1max, jy2min, jy2max;
	int jz1min, jz1max, jz2min, jz2max;
	Tarrecrec(void) { }

public:
	Tarrecrec(TargetManager *man);
	virtual ~Tarrecrec(void) { }

public:
	virtual void Sizer(void);
	virtual void Descriptor(void);
	virtual void Vector(void);
	virtual void VectorX(void);
	virtual void Allocator(void);
	virtual void OutShpar(FILE *file);
};

class Target_Rctgrctg : public Tarrecrec 
{ 
protected: 
	Target_Rctgrctg(void) {}

public: 
	Target_Rctgrctg(TargetManager *man) : Tarrecrec(man) {} 
	virtual ~Target_Rctgrctg(void) {}
	virtual void SayHello(FILE *stream);
}; 

class Target_Recrecpbc : public Tarrecrec 
{ 
protected: 
	Target_Recrecpbc(void) {}

public: 
	Target_Recrecpbc(TargetManager *man) : Tarrecrec(man) {} 
	virtual ~Target_Recrecpbc(void) {} 
	virtual void SayHello(FILE *stream);
};

#endif // __TARRECREC_H__
