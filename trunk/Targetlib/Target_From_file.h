#ifndef __TARGETFROMFILE_H__
#define __TARGETFROMFILE_H__

#include "LoadableTarget.h"
#include "Targetlib.h"

class TARGETLIB_API Target_From_file : public LoadableTarget
{ 
protected: 

public: 
	Target_From_file(void) {}
	Target_From_file(TargetManager *man) : LoadableTarget(man) {} 
	virtual ~Target_From_file(void) {}

public:
	virtual void Sizer(void);
	virtual void Descriptor(void);
	virtual void Vector(void);
	virtual void VectorX(void);
	void PrepareIaniso();
}; 

class TARGETLIB_API Target_Frmfilpbc : public LoadableTarget
{ 
protected: 
	Target_Frmfilpbc(void) {}

public: 
	Target_Frmfilpbc(TargetManager *man) : LoadableTarget(man) {} 
	virtual ~Target_Frmfilpbc(void) {} 
	void SayHello(FILE *stream);
	virtual void PreparePyzd(void);
	virtual void Sizer(void);
	void PrepareIaniso();
}; 

#endif // __TARGETFROMFILE_H__