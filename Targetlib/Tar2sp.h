#ifndef __TAR2SP_H__
#define __TAR2SP_H__

#include "AbstractTarget.h"
#include "Targetlib.h"

class TARGETLIB_API Tar2sp : public AbstractTarget
{
protected:
	Tar2sp(void) { }
	real xc1, yc1, zc1;
	int nat1, nat2;

public:
	Tar2sp(TargetManager *man);
	virtual ~Tar2sp(void) { }

public:
	virtual void Sizer(void);
	virtual void Descriptor(void);
	virtual void Vector(void);
	virtual void VectorX(void);
	virtual void Allocator(void);
	virtual void OutShpar(FILE *file);
};

class Target_Sphroid2 : public Tar2sp 
{ 
protected: 
	Target_Sphroid2(void) {}

public: 
	Target_Sphroid2(TargetManager *man) : Tar2sp(man) {} 
	virtual ~Target_Sphroid2(void) {} 
	virtual void SayHello(FILE *stream);
};

#endif // __TAR2SP_H__