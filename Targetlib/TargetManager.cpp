#include "StdAfx.h"
#include "TargetManager.h"
#include "LoadableTarget.h"

#ifndef _WIN32
char *strlwr(char *Buffer)
{
	const int diff = ('A' - 'a');
	size_t len = strlen(Buffer);
	for(size_t i=0; i<len; ++i)
	{
		if ((Buffer[i] >= 'A') && (Buffer[i] <= 'Z'))
			Buffer[i] -= diff;
	}
	return Buffer;
}
#endif

TargetManager *TargetManager::manager = NULL;

TargetManager::TargetManager(void) : currentTarget(NULL)
{
	factory = new map<TargetType, pair<string, TargetCreator> >;
	helper = new map<string, TargetType>;
	shparer = new map<TargetType, Item *>;
	hello = new map<TargetType, string>;
	descriptor = new map<TargetType, TargetDescriptor>;
	Clear();
}

TargetManager::~TargetManager(void)
{
	CleanDelete(factory);
	CleanDelete(helper);
	CleanDelete(shparer);
	CleanDelete(hello);
	CleanDelete(descriptor);
	if (shpar)
		free(shpar);
	shpar = NULL;
	if (cashedFilename)
		free(cashedFilename);
	cashedFilename = NULL;
}

void TargetManager::Clear(void)
{
	currentType = TargetType_End;
	icomp_need = 1;
	ioshp = false;
	cashedDx.Clear();
	jpbc = PeriodicNo;
	ncomp = 0;
	shpar = NULL;
	shparNum = 0;
	cashedFilename = NULL;
	currentTarget = NULL;
}

void TargetManager::Renew(void)
{
	if (shpar)
		free(shpar);
	if (cashedFilename)
		delete [] cashedFilename;
	if (currentTarget)
		delete currentTarget;
	Clear();
}

bool TargetManager::RegisterTarget(TargetType Type, const char *Name, TargetCreator Creator, int numShpar, bool bFile, int numJpbc, int numNcomp, const char *longText, TargetDescriptor Descriptor)
{
	char nn[32];
	PrepareKey(Name, nn);
	pair<string, TargetCreator> pax = pair<string, TargetCreator>(string(nn), Creator);
	factory->insert(pair<TargetType, pair<string, TargetCreator> >(Type, pax));
	helper->insert(pair<string, TargetType>(string(nn), Type));
	Item *item = new Item(numShpar, bFile, numJpbc, numNcomp);
	shparer->insert(pair<TargetType, Item *>(Type, item));
	hello->insert(pair<TargetType, string>(Type, string(longText)));
	descriptor->insert(pair<TargetType, TargetDescriptor>(Type, Descriptor));

	return true;
}

TargetManager *TargetManager::GetInstance(void)
{
	if (!manager)
		manager = new TargetManager;

	return manager;
}

void TargetManager::Kill(void)
{
	if (manager)
		delete manager;
}

void TargetManager::PrepareKey(const char *key, char *out)
{
	strcpy(out, key);
	strlwr(out);
	out[0] += ('A' - 'a');
}

TargetType TargetManager::TargetTypeFromString(const char *key)
{
	char temp[32];
	PrepareKey(key, temp);
	map<string, TargetType>::iterator ia = helper->find(string(temp));
	if (ia == helper->end())
		return TargetType_End;
	else
		return ia->second;
}

const char *TargetManager::StringFromTargetType(TargetType type)
{
	map<TargetType, pair<string, TargetCreator> >::iterator ia = factory->find(type);
	if (ia == factory->end())
		return NULL;
	else
		return ia->second.first.c_str();
}

void TargetManager::PreloadTarget(const char *key, const char *Buffer)
{
	TargetType type = TargetTypeFromString(key);
	map<TargetType, Item *>::iterator ia = shparer->find(type);
	if (ia != shparer->end())
	{
		shparNum = ia->second->numShpar;
		if (shparNum)
			shpar = (real *)malloc(sizeof(real) * shparNum);
		if (Buffer)
		{
			const char *finis = Buffer;
			if (shparNum)
				finis = Extract(Buffer, shparNum, shpar);
			if (ia->second->useFile)
			{
				if (cashedFilename)
					free(cashedFilename);
				cashedFilename = (char *)malloc(strlen(finis) + 1);
				strcpy(cashedFilename, finis + 1);
				cashedFilename[strlen(cashedFilename) - 1] = '\0';
			}
		}
		else
		{
			TargetDescriptor desc = descriptor->find(type)->second;
			fprintf(stdout, "The target %s uses %d parameters. Please, provide them one per line.\n", key, shparNum);
			for(int i=0; i<shparNum; ++i)
			{
				char Buffer[32];
				real a;
				fprintf(stdout, "%2d.%s\n", i+1, desc(i));
				fflush(stdin);
				fgets(Buffer, 32, stdin);
				if (strchr(Buffer, '.'))
				{
                    sscanf(Buffer, realFormat, &a);
				}
				else
				{
					int ja;
					sscanf(Buffer, "%d", &ja);
					a = (real)ja;
				}
				shpar[i] = a;
			}
			if (ia->second->useFile)
			{
				char Bf[256];
				fprintf(stdout, "The target needs file name. Please, enter.\n");
				gets(Bf);
				if (cashedFilename)
					free(cashedFilename);
				cashedFilename = (char *)malloc(strlen(Bf) + 1);
				strcpy(cashedFilename, Bf);
			}
		}
//
		if (ia->second->useFile)
		{
			Wrimsg("Reapar", "shape file=");
			Wrimsg("Reapar", cashedFilename);
			OutParam("Reapar");
		}
//	
// Need to set JPBC for PBC targets
// JPBC	= 0 if PBC are not used
//		= 1 if PBC in y direction only
//		= 2 if PBC in z direction only
//		= 3 if PBC in y and z directions
		jpbc = PeriodicNo;
		int jpbctmp = ia->second->numJpbc;
		if (jpbctmp > -1)
		{
			if (shpar[jpbctmp] > (real)0.)
				jpbc = PeriodicY;
			if (jpbctmp + 1 < shparNum)
			{
				if (shpar[jpbctmp + 1] > (real)0.)
				{
					if (jpbc == PeriodicNo)
						jpbc = PeriodicZ;
					else
						jpbc = PeriodicBoth;
				}
			}
			else
			{
				if (jpbc == PeriodicNo)
					jpbc = PeriodicZ;
				else
					jpbc = PeriodicBoth;
			}
		}
	}
}

bool TargetManager::AnalyseNcomp(const char *key, int ncomp)
{
	this->ncomp = ncomp;
	TargetType type = TargetTypeFromString(key);
	map<TargetType, Item *>::iterator ia = shparer->find(type);
	if (ia != shparer->end())
	{
		int ncomptmp = ia->second->numNcomp;
		if (ncomptmp != 0)
		{
			char out[256];
			if (ncomptmp > 0)
			{
				if (ncomp != ncomptmp)
				{
					sprintf(out, "Ncomp must be %d for target %s\n", ncomptmp, key);
					Wrimsg("Reapar", out);
					return false;
				}
			}
			else
			{
				int strt = abs(ncomptmp);
				int j = 4;
				while(j > 1)
				{
					if (shpar[strt] <= (real)0.) 
					{
						--j;
						--strt;
					}
					else
						break;
				}
				if (ncomp < j) 
				{
					sprintf(out, "%s: need ncomp = %d", key, j);
					Wrimsg("Reapar", out);
					return false;
				}
			}
		}
		return true;
	}
	else
		return false;
}

AbstractTarget *TargetManager::LoadTarget(const char *key, const char *fileName)
{
	TargetType type = TargetTypeFromString(key);
	if (type != TargetType_End)
	{
		if (currentType != type)
		{
			delete currentTarget;
			map<TargetType, pair<string, TargetCreator> >::iterator ia = factory->find(type);
			if (ia == factory->end())
				Clear();
			else
			{
				currentTarget = (ia->second).second();
				currentType = type;
				if (currentTarget->IsLoadableTarget())
				{
					if (fileName)
						((LoadableTarget *)currentTarget)->SetFileName(fileName);
					icomp_need = ((LoadableTarget *)currentTarget)->IcompNeed();
				}
				currentTarget->Ncomp() = ncomp;
			}
		}
	}
	else
		Clear();

	return currentTarget;
}

const char *TargetManager::Extract(const char *Input, int num, real *shpar)
{
	char Buffer[256];
	strcpy(Buffer, Input);

	char *ia = Buffer;
	char *ib = Buffer;
	for(int i=0; i<num; ++i)
	{
		ia = strtok(ib, " \t");
		if (strchr(ia, '.'))
			shpar[i] = (real)atof(ia);
		else
			shpar[i] = (real)atoi(ia);
		ib = NULL;
	}
	return strtok(NULL, " \t");
}

void TargetManager::OutParam(const char *Label)
{
	char Buffer[256];
	for(int j=0; j<shparNum; ++j)
	{
		sprintf(Buffer, "sphar[%2d] = %lf", j, shpar[j]);
		Wrimsg(Label, Buffer);
	}
}

void TargetManager::SayHello(void)
{
	map<TargetType, string>::iterator ia, ib;
	for(ia=hello->begin(), ib=hello->end(); ia!=ib; ++ia)
	{
		map<TargetType, pair<string, TargetCreator> >:: iterator ic = factory->find(ia->first);
		fprintf(stdout, "%s - %s\n", ic->second.first.c_str(), ia->second.c_str());
	}
}

void TargetManager::ExtendShpar(int cnt)
{
	if (cnt > 0)
	{
		shpar = (real *)realloc(shpar, (cnt + shparNum) * sizeof(real));
		shparNum += cnt;
	}
}