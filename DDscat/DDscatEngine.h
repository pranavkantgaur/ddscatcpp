#ifndef __DDSCAT_ENGINE_H__
#define __DDSCAT_ENGINE_H__

#include "DDscatCommons.h"
#include "TargetManager.h"
#include "DDscatMain.h"
#include "FileNamer.h"
#include "OutputManager.h"
#include "GreenFunctionManager.h"

class DDscatEngine
{
protected:
	char Message[256], DDscatLabel[7];
	string cflpar_default, cflpar;
	NearfieldMethod nrfld;
	TargetManager *manager;
	DDscatParameters *paramDDscat;
	OutputManager *outputer;
	GreenFunctionManager *greener;

public:
	DDscatEngine(void);
	virtual ~DDscatEngine(void);

public:
	void Run(void);
	void SetParameterFile(const char *arg = NULL);
	void SayHello(void);
	void SayBye(void);

public:
	inline NearfieldMethod GetNrfld(void) { return nrfld; }
};

#endif