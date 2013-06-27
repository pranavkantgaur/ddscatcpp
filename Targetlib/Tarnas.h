#ifndef __TARNAS_H__
#define __TARNAS_H__

#include "TargetManager.h"
#include "LoadableTarget.h"
#include "Targetlib.h"

class TARGETLIB_API Tarnas : public LoadableTarget
{
protected:
	Tarnas(void) { }

public:
	Tarnas(TargetManager *man);
	virtual ~Tarnas(void) { }

public:
	virtual void Sizer(void);
	virtual void Descriptor(void);
	virtual void Vector(void);
	virtual void VectorX(void);
	virtual void Allocator(void);
	virtual void Reader(void);
};

class TARGETLIB_API Target_Sphrn_pbc : public Tarnas 
{ 
protected: 
	Target_Sphrn_pbc(void) {}

public: 
	Target_Sphrn_pbc(TargetManager *man) : Tarnas(man) 
	{ 
		man->ExtendShpar(1);
		shpar = man->GetShpar();
		shpar[3] = shpar[2]; shpar[2] = shpar[1]; shpar[1] = (real)0.; 
	} 
	virtual ~Target_Sphrn_pbc(void) {} 
	void SayHello(FILE *stream);
	virtual void PreparePyzd();
}; 

class TARGETLIB_API Target_SphAniN : public Tarnas 
{ 
protected: 
	Target_SphAniN(void) {}

public: 
	Target_SphAniN(TargetManager *man) : Tarnas(man) 
	{ 
		man->ExtendShpar(2);
		shpar = man->GetShpar();
		shpar[2] = shpar[3] = (real)0.; 
	} 
	virtual ~Target_SphAniN(void) {}
	void SayHello(FILE *stream);
	void PrepareIaniso();
}; 

#endif // __TARNAS_H__