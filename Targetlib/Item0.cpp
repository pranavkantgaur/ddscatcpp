#include "StdAfx.h"
#include "Item0.h"

Item0::Item0(void)
{
	as = as2 = xs = ys = zs = (real)0.;
}

Item0::~Item0(void) 
{ 

}

void Item0::Sscanf(const char *Buffer, const char *Format)
{
	sscanf(Buffer, Format, &xs, &ys, &zs, &as);
	as2 = (real)0.;
}

real Item0::DistSquared(real x, real y, real z)
{
	return (x - xs) * (x - xs) + (y - ys) * (y - ys) + (z - zs) * (z - zs);
}

Item::Item(void) 
{ 
	be = ph = th = (real)0.; 
	ic1 = ic2 = ic3 = 0;
}

Item::~Item(void) 
{ 

}

void Item::Sscanf(const char *Buffer, const char *Format)
{
	sscanf(Buffer, Format, &xs, &ys, &zs, &as, &ic1, &ic2, &ic3, &th, &ph, &be);
	as2 = (real)0.;
}

real Item::Dist2(Item *op)
{
	return (xs - op->xs) * (xs - op->xs) + (ys - op->ys) * (ys - op->ys) + (zs - op->zs) * (zs - op->zs) - (as + op->as) * (as + op->as);
}

bool Item::Overlap(Item *op)
{
	return ((ic1 != op->ic1) || (ic2 != op->ic2) || (ic3 != op->ic3));
}

bool Item::IsDisoriented(Item *op)
{
	return ((be != op->be) || (ph != op->ph) || (th != op->th));
}

bool Item::IsAniso()
{
	return ((ic1 != ic2) || (ic1 != ic3));
}
