#ifndef __TARGETMANAGER_H__
#define __TARGETMANAGER_H__

#ifdef MSVC
    #pragma warning(disable:4251)
#endif

#include <string>
#include <map>
#include <vector>

#include "Targetlib.h"
#include "TargetDefinitions.h"
#include "AbstractTarget.h"
#include "Functions.h"
#include "Enumerator.h"
#include "Vect6.h"

typedef AbstractTarget *(*TargetCreator)();
typedef const char *(*TargetDescriptor)(int);

/* **
The class to manipulate the targets with.

Copyright (c) 2012, C++ version, V.Choliy
This code is covered by the GNU General Public License.
** */
class TARGETLIB_API TargetManager
{
protected:
	class Item
	{
	public:
		int numShpar, numJpbc, numNcomp;
		bool useFile;
	public:
		Item(void) { numShpar = 0; useFile = false; numJpbc = -1; numNcomp = 0;}
		Item(int ns, bool bf, int nj, int nc) { numShpar = ns; useFile = bf; numJpbc = nj; numNcomp = nc; }
		~Item(void) { }
	};
protected:
	AbstractTarget *currentTarget;
	TargetType currentType;
	Vect3<real> cashedDx;
	int icomp_need;
	bool ioshp;
	real *shpar;
	int shparNum, ncomp;
	char *cashedFilename;
	PeriodicBoundaryFlag jpbc;

protected:
	TargetManager(void);
	virtual ~TargetManager(void);

protected:
	static TargetManager *manager;
	map<TargetType, pair<string, TargetCreator> > *factory;
	map<string, TargetType> *helper;
	map<TargetType, Item *> *shparer;
	map<TargetType, string> *hello;
	map<TargetType, TargetDescriptor> *descriptor;

public:
	static TargetManager *GetInstance(void);
	static void Kill(void);
	TargetType TargetTypeFromString(const char *key);
	const char *StringFromTargetType(TargetType type);

public:
	bool RegisterTarget(TargetType Type, const char *Name, TargetCreator Creator, int numShpar, bool bFile, int numJpbc, int numNcomp, const char *longText, TargetDescriptor Descriptor);
	AbstractTarget *LoadTarget(const char *key, const char *fileName);
	void PreloadTarget(const char *key, const char *Buffer = NULL);
	bool AnalyseNcomp(const char *key, int ncomp);
	void SayHello(void);
	void Renew(void);
	void ExtendShpar(int cnt);

public:
	inline AbstractTarget *GetCurrentTarget(void) { return currentTarget; }
	inline char *GetCashedFileName() { return cashedFilename; }
	inline Vect3<real> &CashedDx(void) { return cashedDx; }
	inline int IcompNeed(void) { return icomp_need; }
	inline bool &Ioshp(void) { return ioshp; }
	inline real GetShpar(int index) { return shpar[index]; }
	inline real *GetShpar(void) { return shpar; }
	inline PeriodicBoundaryFlag &Jpbc(void) { return jpbc; }
	inline int &Ncomp(void) { return ncomp; }

protected:
	const char *Extract(const char *Buffer, int num, real *shpar);
	void OutParam(const char *Label);
	void Clear(void);
	void PrepareKey(const char *key, char *out);
};

#endif
