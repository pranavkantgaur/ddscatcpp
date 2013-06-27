#ifndef __TARRCTELL_H__
#define __TARRCTELL_H__

#include "LoadableTarget.h"
#include "Targetlib.h"

class Tarrctell : public AbstractTarget
{
protected:
	int nin;
	real xoff, yoff, zoff;
	Tarrctell(void) { }

public:
	Tarrctell(TargetManager *man);
	virtual ~Tarrctell(void){ }

public:
	virtual void Sizer(void);
	virtual void Descriptor(void);
	virtual void Vector(void);
	virtual void VectorX(void);
	virtual void Allocator(void);
};

class Target_Elinrct : public Tarrctell 
{ 
protected: 

public: 
	Target_Elinrct(void) {}
	Target_Elinrct(TargetManager *man) : Tarrctell(man) {} 
	virtual ~Target_Elinrct(void) {} 
	virtual void SayHello(FILE *stream);
};

#endif // __TARRCTELL_H__