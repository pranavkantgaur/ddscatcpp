#ifndef __TARNEL_H__
#define __TARNEL_H__

#include "AbstractTarget.h"
#include "Targetlib.h"

class TARGETLIB_API TarNel : public AbstractTarget
{
protected:
	int mySize;						// amount of identical ellipsoids in the chain
	int delta;						// distance between their surfaces along x
	TarNel(void);

public:
	TarNel(TargetManager *man);
	virtual ~TarNel(void);

public:
	virtual void Sizer(void);
	virtual void Allocator(void);
	virtual void Vector(void);
	virtual void VectorX(void);
	virtual void OutShpar(FILE *file);
	virtual void Descriptor(void);
};

class TARGETLIB_API Target_AniEllN : public TarNel 
{ 
protected: 
	Target_AniEllN(void) {}

public:
	Target_AniEllN(TargetManager *man) : TarNel(man) {}
	virtual ~Target_AniEllN(void) {}
	virtual void SayHello(FILE *stream);
	virtual void Composer(int index, int item);
	void PrepareIaniso();
};

class TARGETLIB_API Target_EllipsoN : public TarNel 
{ 
protected: 
	Target_EllipsoN(void) {}

public:
	Target_EllipsoN(TargetManager *man) : TarNel(man) {}
	virtual ~Target_EllipsoN(void) {} 
	virtual void SayHello(FILE *stream);
	virtual void Composer(int index, int item);
}; 

#endif //__TARNEL_H__