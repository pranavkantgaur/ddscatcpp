#ifndef __TAR2EL_H__
#define __TAR2EL_H__

#include "TarNel.h"
#include "Targetlib.h"

class TARGETLIB_API Tar2el : public TarNel
{
protected:
	Tar2el(void) { mySize = 2; }

public:
	Tar2el(TargetManager *man);
	virtual ~Tar2el(void) { }
};

class TARGETLIB_API Target_AniEll2 : public Tar2el 
{ 
protected: 
	Target_AniEll2(void) {}

public:
	Target_AniEll2(TargetManager *man) : Tar2el(man) {}
	virtual ~Target_AniEll2(void) {}
	virtual void SayHello(FILE *stream);
	virtual void Composer(int index, int item);
	void PrepareIaniso();
};

class TARGETLIB_API Target_Ellipso2 : public Tar2el 
{ 
protected: 
	Target_Ellipso2(void) {}

public:
	Target_Ellipso2(TargetManager *man) : Tar2el(man) {}
	virtual ~Target_Ellipso2(void) {} 
	virtual void SayHello(FILE *stream);
	virtual void Composer(int index, int item);
}; 

#endif //__TAR2EL_H__