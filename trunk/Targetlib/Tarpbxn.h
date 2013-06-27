#ifndef __TARPBXN_H__
#define __TARPBXN_H__

#include "AbstractTarget.h"
#include "Targetlib.h"
#include "TargetManager.h"

class Tarpbxn : public AbstractTarget
{
protected:
	Tarpbxn(void) { }

public:
	Tarpbxn(TargetManager *man);
	virtual ~Tarpbxn(void) { }

public:
	virtual void Sizer(void);
	virtual void Descriptor(void);
	virtual void Vector(void);
	virtual void VectorX(void);
	virtual void Allocator(void);
	virtual void OutShpar(FILE *file);
	virtual void ShiftDipolesAndX(void);

protected:
	void InitSlab(int jxmin, int jxmax, int element, int &curSize);
	void InitDisk(int jxmin, int jxmax, int element, int &curSize);
};

class Target_Dskblypbc : public Tarpbxn 
{ 
protected: 
	Target_Dskblypbc(void) {}

public: 
	Target_Dskblypbc(TargetManager *man) : Tarpbxn(man) {} 
	virtual ~Target_Dskblypbc(void) {} 
	virtual void SayHello(FILE *stream);
}; 

class Target_Dskrctngl : public Tarpbxn 
{ 
protected: 
	Target_Dskrctngl(void) {}

public:
	Target_Dskrctngl(TargetManager *man) : Tarpbxn(man)
	{
		man->ExtendShpar(1);
		shpar = man->GetShpar();
		for(int i=5; i>3; --i)
			shpar[i] = shpar[i-1];
		shpar[3] = (real)0.;
	}
	virtual ~Target_Dskrctngl(void) {}
	virtual void SayHello(FILE *stream);
};

class Target_Dskrctpbc : public Tarpbxn 
{ 
protected: 
	Target_Dskrctpbc(void) {}

public: 
	Target_Dskrctpbc(TargetManager *man) : Tarpbxn(man) 
	{
		man->ExtendShpar(1);
		shpar = man->GetShpar();
		for(int i=7; i>3; --i)
			shpar[i] = shpar[i-1];
		shpar[3] = (real)0.;
	} 
	virtual ~Target_Dskrctpbc(void) {} 
	virtual void SayHello(FILE *stream);
	virtual void PreparePyzd();
};

#endif // __TARPBXN_H__