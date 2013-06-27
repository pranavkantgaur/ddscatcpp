#ifndef __LINIA_H__
#define __LINIA_H__

#include "Processlib.h"
#include "Definitions.h"

class PROCESSLIB_API Linia
{
protected:
	real data[6];
	int nab;

public:
	Linia(void);
	Linia(real xxa, real xxb, real yya, real yyb, real zza, real zzb, int nnab);
	virtual ~Linia(void);

public:
	void Scanf(char *Buffer);
	void Debug(void);
	void ParametricPoint(real zeta, real *xtf) const;

public:
	inline real GetData(int index) const { return data[index]; }
	inline int GetNab(void) const { return nab; }
};

#endif