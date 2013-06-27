#include "StdAfx.h"
#include "DipoleData.h"

DipoleData::DipoleData() : cxe_tf(NULL), cxxi(NULL), curSize(0)
{ 

}

DipoleData::~DipoleData() 
{ 
	Delete(); 
}

void DipoleData::Allocate(int size)
{
	if (curSize != size)
	{
		Delete();
	}
	cxe_tf = new Complex[size];
	cxxi = new Complex[size];
	curSize = size;
}

void DipoleData::Delete()
{
	CleanDelete2(cxe_tf);
	CleanDelete2(cxxi);
	curSize = 0;
}

void DipoleData::ClearCxxi()
{
	memset(cxxi, 0, curSize * sizeof(Complex));
}

void DipoleData::ClearCxeTf()
{
	memset(cxe_tf, 0, curSize * sizeof(Complex));
}

void DipoleData::Debug(FILE *file, int from, int toto)
{
	fprintf(file, "Size = %d\n", curSize);
	if ((from == 0) && (toto == 0))
		toto = curSize;
	for(int i=from; i<toto; ++i)
	{
		if (cxxi[i].re != (real)0. || cxxi[i].im != (real)0. || cxe_tf[i].re != (real)0. || cxe_tf[i].im != (real)0.)
			fprintf(file, "%3d %20.8e %20.8e %20.8e %20.8e\n", i, cxxi[i].re, cxxi[i].im, cxe_tf[i].re, cxe_tf[i].im);
	}
}