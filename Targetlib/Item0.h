#pragma once

#include "Definitions.h"

class Item0
{
protected:
	real as, as2, xs, ys, zs;
	
public:
	Item0(void);
	~Item0(void);

public:
	inline real &As() { return as; }
	inline real &As2() { return as2; }
	inline real &Xs() { return xs; }
	inline real &Ys() { return ys; }
	inline real &Zs() { return zs; }
	inline real GetAsCube() { return as * as * as;} 

public:
	void Sscanf(const char *Buffer, const char *Format);
	real DistSquared(real x, real y, real z);
};

class Item : public Item0
{
protected:
	real be, ph, th;
	int ic1, ic2, ic3;

public:
	Item(void);
	~Item(void);

public:
	inline real &Be() { return be; }
	inline real &Ph() { return ph; }
	inline real &Th() { return th; }
	inline int &Ic1() { return ic1; }
	inline int &Ic2() { return ic2; }
	inline int &Ic3() { return ic3; }

public:
	void Sscanf(const char *Buffer, const char *Format);
	real Dist2(Item *op);
	bool Overlap(Item *op);
	bool IsDisoriented(Item *op);
	bool IsAniso();
};