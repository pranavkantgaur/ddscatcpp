#ifndef __TARRCTBLK3_H__
#define __TARRCTBLK3_H__

#include "LoadableTarget.h"
#include "Targetlib.h"

class Tarrctblk3 : public AbstractTarget
{
protected:
	Tarrctblk3(void) { }

public:
	Tarrctblk3(TargetManager *man);
	virtual ~Tarrctblk3(void) { }

public:
	virtual void Sizer(void);
	virtual void Descriptor(void);
	virtual void Vector(void);
	virtual void Allocator(void);
	virtual void OutShpar(FILE *file);

protected:
	int ElementaryBlock(int jxmin, int jxmax, int jymin, int jymax, int jzmin, int jzmax, int comp, int &curSize);
};

class Target_Bislinpbc : public Tarrctblk3 
{ 
protected:
	Target_Bislinpbc(void) {}

public: 
	Target_Bislinpbc(TargetManager *man) : Tarrctblk3(man) {}
	virtual ~Target_Bislinpbc(void) {} 
	virtual void SayHello(FILE *stream);
};

class Target_Rctglblk3 : public Tarrctblk3 
{ 
protected: 
	Target_Rctglblk3(void) {}

public:
	Target_Rctglblk3(TargetManager *man) : Tarrctblk3(man) {}
	virtual ~Target_Rctglblk3(void) {}
	virtual void SayHello(FILE *stream);
};

class Target_Trilyrpbc : public Tarrctblk3 
{ 
protected: 
	Target_Trilyrpbc(void) {}

public: 
	Target_Trilyrpbc(TargetManager *man) : Tarrctblk3(man) {}  
	virtual ~Target_Trilyrpbc(void) {}
	virtual void SayHello(FILE *stream);
};

#endif // __TARRCTBLK3_H__