#ifndef __TAR3EL_H__
#define __TAR3EL_H__

#include "TarNel.h"
#include "Targetlib.h"

class TARGETLIB_API Tar3el : public TarNel
{
protected:
	Tar3el(void) { mySize = 3; }

public:
	Tar3el(TargetManager *man);
	virtual ~Tar3el(void) { }
};

class TARGETLIB_API Target_AniEll3 : public Tar3el 
{ 
protected:
	Target_AniEll3(void) {}

public:
	Target_AniEll3(TargetManager *man) : Tar3el(man) {}
	virtual ~Target_AniEll3(void) {}
	virtual void SayHello(FILE *stream);
	virtual void Composer(int index, int item);
	void PrepareIaniso();
}; 

class TARGETLIB_API Target_Ellipso3 : public Tar3el 
{ 
protected: 
	Target_Ellipso3(void) {}

public:
	Target_Ellipso3(TargetManager *man) : Tar3el(man) {}
	virtual ~Target_Ellipso3(void) {}
	virtual void SayHello(FILE *stream);
	virtual void Composer(int index, int item);
}; 

#endif // __TAR3EL_H__