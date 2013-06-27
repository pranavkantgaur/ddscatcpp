#include "StdAfx.h"

#include "Linia.h"

Linia::Linia(void)
{
	memset(data, 0, 6*sizeof(real));
	nab = 0;
}

Linia::~Linia(void)
{

}

Linia::	Linia(real xxa, real yya, real zza, real xxb, real yyb, real zzb, int nnab)
{
	data[0] = xxa;
	data[1] = yya;
	data[2] = zza;
	data[3] = xxb;
	data[4] = yyb;
	data[5] = zzb;
	nab = nnab;
}

void Linia::Scanf(char *Buffer)
{
	char *ia = strtok(Buffer, " \t");
	for(int i=0; i<6; ++i)
	{
		data[i] = (real)atof(ia);
		ia = strtok(NULL, " \t");
	}
	nab = atoi(ia);
}

void Linia::Debug()
{
	fprintf(stderr, "xa,ya,za = %13.2e%13.2e%13.2e\n", data[0], data[1], data[2]);
	fprintf(stderr, "xb,yb,zb = %13.2e%13.2e%13.2e\n", data[3], data[4], data[5]);
}

void Linia::ParametricPoint(real zeta, real *xtf) const
{
	xtf[0] = data[0] + (data[3] - data[0]) * zeta;
	xtf[1] = data[1] + (data[4] - data[1]) * zeta;
	xtf[2] = data[2] + (data[5] - data[2]) * zeta;
}