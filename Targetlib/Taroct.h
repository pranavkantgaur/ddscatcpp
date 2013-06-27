#ifndef __TAROCT_H__
#define __TAROCT_H__

#include "AbstractTarget.h"
#include "Targetlib.h"

class TARGETLIB_API Taroct : public AbstractTarget
{
protected:
	Taroct(void) { }
	real xoff, yoff, zoff;

public:
	Taroct(TargetManager *man);
	virtual ~Taroct(void) { }

public:
	virtual void Sizer(void);
	virtual void Descriptor(void);
	virtual void Vector(void);
	virtual void VectorX(void);
	virtual void Allocator(void);
	virtual void OutShpar(FILE *file);
};

class TARGETLIB_API Target_OctPrism : public Taroct
{ 
protected: 
	Target_OctPrism(void) {}

public: 
	Target_OctPrism(TargetManager *man) : Taroct(man) {} 
	virtual ~Target_OctPrism(void) {} 
	virtual void SayHello(FILE *stream);
}; 

#endif // __TAROCT_H__