#ifndef __TARANIREC_H__
#define __TARANIREC_H__

#include "AbstractTarget.h"
#include "Targetlib.h"

class TARGETLIB_API Taranirec : public AbstractTarget
{
protected:
	Taranirec(void) { }

public:
	Taranirec(TargetManager *man);
	virtual ~Taranirec(void) { }

public:
	virtual void Sizer(void);
	virtual void Descriptor(void);
	virtual void Vector(void);
	virtual void VectorX(void);
	virtual void Allocator(void);
	virtual void OutShpar(FILE *file);
};

class TARGETLIB_API Target_Anirctngl : public Taranirec 
{ 
protected: 
	Target_Anirctngl(void) {}

public: 
	Target_Anirctngl(TargetManager *man) : Taranirec(man) {}
	virtual ~Target_Anirctngl(void) {}
	void SayHello(FILE *stream);
	void PrepareIaniso();
}; 

#endif // __TARANIREC_H__