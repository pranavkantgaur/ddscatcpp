#ifndef __TARBLOCKS_H__
#define __TARBLOCKS_H__

#include "AbstractTarget.h"
#include "Targetlib.h"

class TARGETLIB_API Tarblocks : public AbstractTarget
{
protected:
	vector<Vect3<real> *> xyzb;
	Vect3<real> eigval;
	int blksiz, nblocks, iprinax;
	Tarblocks(void) { }

public:
	Tarblocks(TargetManager *man);									// shpar == (int nblocks, int blksiz, int iprinax)
	virtual ~Tarblocks(void) { }

public:
	void Loader(real &miX, real &maX, real &miY, real &maY, real &miZ, real &maZ);
	virtual void Sizer(void);
	virtual void Descriptor(void);
	virtual void Vector(void);
	virtual void VectorX(void);
	virtual void Allocator(void);

public:
	virtual bool IsMultiBlockTarget(void) { return true; }
//	virtual void SetVector(const Vect3<real> &op) { xyzb.push_back(op); }
};

class Target_Dw1996tar : public Tarblocks 
{ 
protected: 
	Target_Dw1996tar(void) {}

public: 
	Target_Dw1996tar(TargetManager *man) : Tarblocks(man) {} 
	virtual ~Target_Dw1996tar(void) {} 
	virtual void SayHello(FILE *stream);
};

class Target_Mltblocks : public Tarblocks 
{ 
protected: 
	Target_Mltblocks(void) {}

public:
	Target_Mltblocks(TargetManager *man) : Tarblocks(man) {}
	virtual ~Target_Mltblocks(void) {}
	virtual void SayHello(FILE *stream);
}; 

#endif // __TARBLOCKS_H__