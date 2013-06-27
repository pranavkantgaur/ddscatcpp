#pragma once

#include "Targetlib.h"
#include "Definitions.h"
#include "Enumerator.h"
#include "Item0.h"

class TARGETLIB_API AbstractDFData
{
protected:
	vector<Item0 *> *elems;
	int mySize, *js;
	AbstractDFData(void);

public:
	AbstractDFData(vector<Item0 *> *op);
	virtual ~AbstractDFData(void);

public:
	void Allocate(int size);
	IsotropicFlag GetIaniso(void);
	bool AreAnglesZero(int index);
	void SetAngles(int index, int jjs);
	real GetBetadf(int index);
	real GetPhidf(int index);
	real GetThetadf(int index);

public:
	inline bool IsAllocated() { return (mySize != 0); }
};