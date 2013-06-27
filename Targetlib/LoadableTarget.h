#ifndef __LOADABLETARGET_H__
#define __LOADABLETARGET_H__

#include "AbstractTarget.h"
#include "Targetlib.h"
#include "Definitions.h"
#include "TargetDefinitions.h"
#include "Vect3.h"
#include "Vect6.h"
#include "AbstractDFData.h"

/* **
The class to manage the targets, which parameters are read from file.

Copyright (C) 2012,2013, C++ versions, Choliy V.
This code is covered by the GNU General Public License.
** */
class TARGETLIB_API LoadableTarget : public AbstractTarget
{
protected:
	AbstractDFData *dfdata;
	char *cflshp;
	int icomp_need;
	Vect3<real> alpha;
	vector<Item0 *> *elements;
	int nsph;
	LoadableTarget(void) { Init(NULL); 	dfdata = NULL; elements = NULL; }

public:
	LoadableTarget(TargetManager *man) : AbstractTarget(man) { Init("shape.dat"); dfdata = NULL; elements = NULL; }
	virtual ~LoadableTarget(void);

public:
	inline AbstractDFData *GetDFdata() { return dfdata; }
	inline int IcompNeed(void) { return icomp_need; }

public:
	virtual void Allocator(void);
	virtual bool IsLoadableTarget(void) { return true; }
	void SetFileName(const char *fileName) { cflshp = new char[strlen(fileName)+1]; strcpy(cflshp, fileName); }
	void Printer(void);
	virtual void Reader(void);

protected:
	void Init(const char *c)
	{
		nsph = 0;
		cflshp = NULL;
		if (c)
		{
			cflshp = new char [strlen(c) + 1];
			strcpy(cflshp, c);
		}
	}
    void OutVectors(FILE *file);
	virtual void PrepareIanisoSpecial(void);
};

class Target_Anifilpbc : public LoadableTarget
{ 
protected: 

public: 
	Target_Anifilpbc(void) {}
	Target_Anifilpbc(TargetManager *man) : LoadableTarget(man) {} 
	virtual ~Target_Anifilpbc(void) {}
	void PrepareIaniso();
}; 

class Target_Anifrmfil : public LoadableTarget 
{ 
protected: 

public: 
	Target_Anifrmfil(void) {}
	Target_Anifrmfil(TargetManager *man) : LoadableTarget(man) {} 
	virtual ~Target_Anifrmfil(void) {} 
	void PrepareIaniso();
}; 

#endif // __LOADABLETARGET_H__