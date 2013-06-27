#ifndef __TARLYRSLAB_H__
#define __TARLYRSLAB_H__

#include "LoadableTarget.h"
#include "Targetlib.h"

class Tarlyrslab : public AbstractTarget
{
protected:
	Tarlyrslab(void) { }

public:
	Tarlyrslab(TargetManager *man);
	virtual ~Tarlyrslab(void) { }

public:
	virtual void Sizer(void);
	virtual void Descriptor(void);
	virtual void Vector(void);
	virtual void VectorX(void);
	virtual void Allocator(void);
};

class Target_Layrdslab : public Tarlyrslab 
{ 
protected: 
	Target_Layrdslab(void) {}

public: 
	Target_Layrdslab(TargetManager *man) : Tarlyrslab(man) {} 
	virtual ~Target_Layrdslab(void) {} 
	virtual void SayHello(FILE *stream);
}; 

class Target_Lyrslbpbc : public Tarlyrslab 
{ 
protected: 
	Target_Lyrslbpbc(void) {}

public: 
	Target_Lyrslbpbc(TargetManager *man) : Tarlyrslab(man) {} 
	virtual ~Target_Lyrslbpbc(void) {} 
	virtual void SayHello(FILE *stream);
};

#endif // __TARLYRSLAB_H__