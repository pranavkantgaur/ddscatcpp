#ifndef __SUMHOLDER_H__
#define __SUMHOLDER_H__

#include "Definitions.h"
#include "Vect3.h"

template <typename T>
class SumHolder
{
protected:
	T q[2], qsum[2], qsum1[2];

public:
	SumHolder(void)
	{
		memset(q, 0, 2*sizeof(T));
		memset(qsum, 0, 2*sizeof(T));
		memset(qsum1, 0, 2*sizeof(T));
	}
	virtual ~SumHolder(void) { }

public:
	inline T *Q() { return q; }
	inline T *Qsum() { return qsum; }
	inline T *Qsum1() { return qsum1; }
	inline T *Qsum(bool bWhat) { return bWhat ? qsum1 : qsum; }
	inline T GetSumQ() { return q[0] + q[1]; }
	inline T GetDifQ() { return q[0] - q[1]; }
};

class SumPackage
{
protected:
	bool bWhatSum;
	SumHolder<real> qabs, qext, qsca, qbksca, qpha, qscag2;
	SumHolder<Vect3<real> > qscag, qtrqsc, qtrqab;

public:
	SumPackage() { bWhatSum = false; }
	virtual ~SumPackage() { }

public:
	inline bool &UseWhatSum() { return bWhatSum; }
	inline SumHolder<real> &Qabs() { return qabs; }
	inline SumHolder<real> &Qext() { return qext; }
	inline SumHolder<real> &Qsca() { return qsca; }
	inline SumHolder<real> &Qbksca() { return qbksca; }
	inline SumHolder<real> &Qpha() { return qpha; }
	inline SumHolder<real> &Qscag2() { return qscag2; }
	inline SumHolder<Vect3<real> > &Qscag() { return qscag; }
	inline SumHolder<Vect3<real> > &Qtrqsc() { return qtrqsc; }
	inline SumHolder<Vect3<real> > &Qtrqab() { return qtrqab; }

public:
	void GeneralSummate(int iorth)
	{
		for(int nd=0; nd<iorth; ++nd)
		{
			qsca.Qsum()[nd]   += qsca.Qsum1()[nd];
			qabs.Qsum()[nd]   += qabs.Qsum1()[nd];
			qext.Qsum()[nd]   += qext.Qsum1()[nd];
			qbksca.Qsum()[nd] += qbksca.Qsum1()[nd];
			qpha.Qsum()[nd]   += qpha.Qsum1()[nd];
			qscag2.Qsum()[nd] += qscag2.Qsum1()[nd];
			qscag.Qsum()[nd]  += qscag.Qsum1()[nd];
			qtrqab.Qsum()[nd] += qtrqab.Qsum1()[nd];
			qtrqsc.Qsum()[nd] += qtrqsc.Qsum1()[nd];
		}
	}
};

#endif // __SUMHOLDER_H__