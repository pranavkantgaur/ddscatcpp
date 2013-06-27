#include "StdAfx.h"

#include "FourArray.h"

FourArray::FourArray() 
{ 
	cxs1 = cxs2 = cxs3 = cxs4 = NULL; 
	curSize = 0; 
}

FourArray::~FourArray() 
{
	Delete();
}

void FourArray::Alloc(int newSize)
{
	if (curSize != newSize)
	{
		if (curSize)
			Delete();
		cxs1 = new Complex[newSize];
		cxs2 = new Complex[newSize];
		cxs3 = new Complex[newSize];
		cxs4 = new Complex[newSize];
		curSize = newSize;
	}
}

void FourArray::Delete()
{
	CleanDelete2(cxs1);
	CleanDelete2(cxs2);
	CleanDelete2(cxs3);
	CleanDelete2(cxs4);
	curSize = 0;
}
