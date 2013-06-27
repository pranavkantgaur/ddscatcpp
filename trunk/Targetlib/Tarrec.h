#ifndef __TARREC_H__
#define __TARREC_H__

#include "AbstractTarget.h"
#include "Targetlib.h"

class TARGETLIB_API Tarrec : public AbstractTarget
{
protected:
	Tarrec(void) { }

public:
	Tarrec(TargetManager *man);
	virtual ~Tarrec(void) { }

public:
	virtual void Sizer(void);
	virtual void Descriptor(void);
	virtual void Vector(void);
	virtual void VectorX(void);
	virtual void Allocator(void);
	virtual void OutShpar(FILE *file);
};

class TARGETLIB_API Target_Rctglpbc : public Tarrec 
{ 
protected: 
	Target_Rctglpbc(void) {}

public: 
	Target_Rctglpbc(TargetManager *man) : Tarrec(man) {} 
	virtual ~Target_Rctglpbc(void) {} 
	virtual void SayHello(FILE *stream);
	virtual void PreparePyzd();
};

class TARGETLIB_API Target_Rctglprsm : public Tarrec 
{ 
protected: 
	Target_Rctglprsm(void) { }

public: 
	Target_Rctglprsm(TargetManager *man) : Tarrec(man) { }
	virtual ~Target_Rctglprsm(void) { }
	virtual void SayHello(FILE *stream);
};

#endif // __TARREC_H__