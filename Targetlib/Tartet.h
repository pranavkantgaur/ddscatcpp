#ifndef __TARTET_H__
#define __TARTET_H__

#include "AbstractTarget.h"
#include "Targetlib.h"

class Tartet : public AbstractTarget
{
protected:
	real xoff, yoff, zoff;
	Tartet(void) { }

public:
	Tartet(TargetManager *man);
	virtual ~Tartet(void) { }

public:
	virtual void Sizer(void);
	virtual void Descriptor(void);
	virtual void Vector(void);
	virtual void VectorX(void);
	virtual void Allocator(void);
	virtual void OutShpar(FILE *file);
};

class Target_Tetrahdrn : public Tartet 
{ 
protected: 
	Target_Tetrahdrn(void) {}

public: 
	Target_Tetrahdrn(TargetManager *man) : Tartet(man) {} 
	virtual ~Target_Tetrahdrn(void) {} 
	virtual void SayHello(FILE *stream);
};

#endif // __TARTET_H__