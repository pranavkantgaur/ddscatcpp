#ifndef __TARNSP_H__
#define __TARNSP_H__

#include "LoadableTarget.h"
#include "Targetlib.h"

class TARGETLIB_API Tarnsp : public LoadableTarget
{
protected:
	Tarnsp(void) { }

public:
	Tarnsp(TargetManager *man);
	virtual ~Tarnsp(void) { }

public:
	virtual void Sizer(void);
	virtual void Descriptor(void);
	virtual void Vector(void);
	virtual void VectorX(void);
	virtual void Allocator(void);
	virtual void Reader(void);
};

class TARGETLIB_API Target_SpheresN : public Tarnsp 
{ 
protected: 
	Target_SpheresN(void) {}

public: 
	Target_SpheresN(TargetManager *man) : Tarnsp(man) {} 
	virtual ~Target_SpheresN(void) {}
	void SayHello(FILE *stream);
}; 

#endif // __TARNSP_H__