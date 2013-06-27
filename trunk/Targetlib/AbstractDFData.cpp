#include "StdAfx.h"
#include "AbstractDFData.h"

const real zero_ = (real)0.;

AbstractDFData::AbstractDFData(void)
{
	elems = NULL;
	mySize = 0;
	js = NULL;
}

AbstractDFData::AbstractDFData(vector<Item0 *> *op)
{
	elems = op;
	mySize = 0;
	js = NULL;
}

AbstractDFData::~AbstractDFData(void)
{
	CleanDelete2(js);
}

void AbstractDFData::Allocate(int size)
{
	js = (int *)realloc(js, size * sizeof(int));
	for(int i=0; i<size; ++i)
		js[i] = -1;
	mySize = size;
}

IsotropicFlag AbstractDFData::GetIaniso(void)
{
	for(int j=0; j<mySize; ++j)
	{
		int jjs = js[j];
		if (jjs == -1)
			continue;
		Item *item = (Item *)elems->at(jjs);
		if(item->Be() != zero_) return TargetIsAnisotropicDisoriented;
		if(item->Ph() != zero_) return TargetIsAnisotropicDisoriented;
		if(item->Th() != zero_) return TargetIsAnisotropicDisoriented;
	}
	return TargetIsAnisotropic;
}

bool AbstractDFData::AreAnglesZero(int index)
{
	if (!IsAllocated())
		return true;
	int jjs = js[index];
	if (jjs == -1)
		return true;
	Item *item = (Item *)elems->at(jjs);
	if(item->Be() != zero_) return false;
	if(item->Ph() != zero_) return false;
	if(item->Th() != zero_) return false;
	return true;
}

void AbstractDFData::SetAngles(int index, int jjs)
{
	js[index] = jjs;
}

real AbstractDFData::GetBetadf(int index)
{
	int jjs = js[index];
	return (jjs == -1) ? zero_ : ((Item *)elems->at(jjs))->Be();
}

real AbstractDFData::GetPhidf(int index)
{
	int jjs = js[index];
	return (jjs == -1) ? zero_ : ((Item *)elems->at(jjs))->Ph();
}

real AbstractDFData::GetThetadf(int index)
{
	int jjs = js[index];
	return (jjs == -1) ? zero_ : ((Item *)elems->at(jjs))->Th();
}
