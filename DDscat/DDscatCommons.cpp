#include "StdAfx.h"

#include "DDscatCommons.h"

Common0 *Common0::item = NULL;
Common0 *Common0::GetInstance()
{
	if (!item)
	{
		item = new Common0;
		item->Init();
	}
	return item;
}

void Common0::Kill()
{
	CleanDelete(item);
}

void Common0::InitReals(void)
{
	data[0].Set(-(real)999.);
	data[1].Set(-(real)999.);
}

void Common0::Init(void)
{
	InitReals();
}

Common1 *Common1::item = NULL;
Common1 *Common1::GetInstance()
{
	if (!item)
	{
		item = new Common1;
		item->Init();
	}
	return item;
}

void Common1::Kill()
{
	CleanDelete(item);
}

void Common1::Init()
{
	ak_tf.Clear();
}

Common4 *Common4::item = NULL;
Common4 *Common4::GetInstance()
{
	if (!item)
	{
		item = new Common4;
		item->Init();
	}
	return item;
}

void Common4::Kill()
{
	CleanDelete(item);
}

void Common4::Init()
{
	cxzw = new Array4Stacked<Complex>;
}

Common4::~Common4()
{ 
	CleanDelete(cxzw); 
}

void Common4::AllocateCxzw(int sX, int sY, int sZ, int sA)
{
	cxzw->Dimension(sX, sY, sZ, sA);
}

Common10 *Common10::item = NULL;
Common10 *Common10::GetInstance()
{
	if (!item)
	{
		item = new Common10;
		item->Init();
	}
	return item;
}

void Common10::Kill()
{
	CleanDelete(item);
}

void Common10::Init()
{
	myid = 0;
}
