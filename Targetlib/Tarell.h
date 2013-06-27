#ifndef __TARELL_H__
#define __TARELL_H__

#include "AbstractTarget.h"
#include "Targetlib.h"

class TARGETLIB_API Tarell : public AbstractTarget
{
protected:
	Tarell(void) { }

public:
	Tarell(TargetManager *man);
	virtual ~Tarell(void) { }

public:
	virtual void Sizer(void);
	virtual void Descriptor(void);
	virtual void Vector(void);
	virtual void VectorX(void);
	virtual void Allocator(void);
	virtual void OutShpar(FILE *file);
};

class TARGETLIB_API Target_Aniellips : public Tarell
{
protected:
	Target_Aniellips(void) {}

public:
	Target_Aniellips(TargetManager *man) : Tarell(man) {}
	virtual ~Target_Aniellips(void) { }
	virtual void SayHello(FILE *stream);
	virtual void Composer(int index, int item);
	void PrepareIaniso();
};

class TARGETLIB_API Target_Ellipsoid : public Tarell 
{ 
protected: 
	Target_Ellipsoid(void) {}

public:
	Target_Ellipsoid(TargetManager *man) : Tarell(man) {}
	virtual ~Target_Ellipsoid(void) {} 
	virtual void SayHello();
};

#endif // __TARELL_H__