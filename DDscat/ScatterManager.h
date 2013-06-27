#ifndef __SCATTERMANAGER_H__
#define __SCATTERMANAGER_H__

#include "AbstractTarget.h"
#include "Complex.h"
#include "Enumerator.h"
#include "DipoleData.h"

class ScatterManager
{
protected:
	AbstractTarget *currentTarget;
	real *scrrs1, *scrrs2;
	TorqMethod torqMethod;
	Complex *cscr1, *cscr2, *cscr3, *cscr4;

public:
	ScatterManager(void);
	ScatterManager(AbstractTarget *target);
	virtual ~ScatterManager(void);

public:
	inline void SetTorqMethod(TorqMethod method) { torqMethod = method; }

public:
    void Scat(Vect3<real> &ak_tf, Vect3<real> *aks_tf, Vect3<real> *em1_tf, Vect3<real> *em2_tf, real e02, real etasca,
		real &cbksca, real &csca, Vect3<real> &cscag, real &cscag2, Vect3<real> &ctrqin, Vect3<real> &ctrqsc, DipoleData *theDipoleData,
		Vect3<Complex> &cxe01_tf, Complex *cxf1l, Complex *cxf2l, int myid, PeriodicBoundaryFlag jpbc, int &navg, int ndir);

protected:
	void Init();
};

#endif // __SCATTERMANAGER_H__