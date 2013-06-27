#ifndef __TARSLABHOLE_H__
#define __TARSLABHOLE_H__

#include "AbstractTarget.h"
#include "Targetlib.h"

class Tarslabhole : public AbstractTarget
{
protected:
	Tarslabhole(void) { }

public:
	Tarslabhole(TargetManager *man);
	virtual ~Tarslabhole(void) { }

public:
	virtual void Sizer(void);
	virtual void Descriptor(void);
	virtual void Vector(void);
	virtual void VectorX(void);
	virtual void Allocator(void);
	virtual void OutShpar(FILE *file);
};

class Target_Slabhole : public Tarslabhole 
{ 
protected: 
	Target_Slabhole(void) {}

public: 
	Target_Slabhole(TargetManager *man) : Tarslabhole(man) {}
	virtual ~Target_Slabhole(void) {}
	virtual void SayHello(FILE *stream);
};

class Target_Slbholpbc : public Tarslabhole 
{ 
protected: 
	Target_Slbholpbc(void) {}

public: 
	Target_Slbholpbc(TargetManager *man) : Tarslabhole(man) {} 
	virtual ~Target_Slbholpbc(void) {}
	virtual void SayHello(FILE *stream);
};

#endif // __TARSLABHOLE_H__