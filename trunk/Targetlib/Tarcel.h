#ifndef __TARCELL_H__
#define __TARCELL_H__

#include "AbstractTarget.h"
#include "Targetlib.h"

class Tarcel : public AbstractTarget
{
protected:
	int nin;
	Tarcel(void) { }

public:
	Tarcel(TargetManager *man);
	virtual ~Tarcel(void) { }

public:
	virtual void Sizer(void);
	virtual void Descriptor(void);
	virtual void Vector(void);
	virtual void VectorX(void);
	virtual void Allocator(void);
	virtual void OutShpar(FILE *file);
};

class Target_Conellips : public Tarcel 
{ 
protected: 
	Target_Conellips(void) {}

public: 
	Target_Conellips(TargetManager *man) : Tarcel(man) {} 
	virtual ~Target_Conellips(void) {} 
	void SayHello(FILE *stream);
};

#endif // __TARCELL_H__