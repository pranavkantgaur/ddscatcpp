#ifndef __TARCYL_H__
#define __TARCYL_H__

#include "AbstractTarget.h"
#include "Targetlib.h"

class Tarcyl : public AbstractTarget
{
protected:
	real xcm, ycm, zcm;
	Tarcyl(void) { }

public:
	Tarcyl(TargetManager *man);
	virtual ~Tarcyl(void) { }

public:
	virtual void Sizer(void);
	virtual void Descriptor(void);
	virtual void Vector(void);
	virtual void VectorX(void);
	virtual void Allocator(void);
	virtual void OutShpar(FILE *file);

protected:
	void PreSizer(int jx, int jy, int jz);
};

class Target_Cylinder1 : public Tarcyl 
{ 
protected: 
	Target_Cylinder1(void) {}

public: 
	Target_Cylinder1(TargetManager *man) : Tarcyl(man) {} 
	virtual ~Target_Cylinder1(void) {} 
	virtual void SayHello(FILE *stream);
};

class Target_Cylndrpbc : public Tarcyl 
{ 
protected: 
	Target_Cylndrpbc(void) {}

public:
	Target_Cylndrpbc(TargetManager *man) : Tarcyl(man) {}
	virtual ~Target_Cylndrpbc(void) {}
	virtual void SayHello(FILE *stream);
	virtual void PreparePyzd();
};

class Target_Uniaxicyl : public Tarcyl 
{ 
protected: 
	Target_Uniaxicyl(void) {}

public: 
	Target_Uniaxicyl(TargetManager *man) : Tarcyl(man) {} 
	virtual ~Target_Uniaxicyl(void) {} 
	virtual void SayHello(FILE *stream);
	void PrepareIaniso();
	void Composer(int index, int item);
};

#endif // __TARCYL_H__