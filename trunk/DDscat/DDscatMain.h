#ifndef __DDSCAT_MAIN_H__
#define __DDSCAT_MAIN_H__

#include "DDscatMainNew.h"

void Namid(int myid, char *cfllog);
const char *Version();
int VersionNum();
void Getfml(AbstractTarget *currentTarget, DDscatParameters *param, Vect3<real> &ak_tf, real ak3, Vect3<real> *aks_tf, Matrix *theMatrix, Matrix *theTensor, DipoleData *theDipoleData, 
	Vect3<Complex> &cxe01_tf, Vect3<Complex> &cxe02_tf, DielectricManager *dielec, FourArray *cxfData, Complex *cxsc, Complex *cxscr1, 
	Vect3<real> *em1_tf, Vect3<real> *em2_tf, int ibeth, int ibeth1, int iphi, int iphi1, int &itask, 
	int *itnum, PeriodicBoundaryFlag jpbc, int lace, int myid, int &navg, SumPackage &sumPackage);
void Getmueller(DDscatParameters *param, int ibeta, int iphi, int itheta, PeriodicBoundaryFlag jpbc, real *orderm, real *ordern, real ak1, Vect3<real> *aks_tf, 
	Vect3<real> *ensc_lf, Vect3<real> *ensc_tf, real pyddx, real pzddx, real *sm, real *smori, real *s1111, real *s2121, Complex *cx1121, 
	FourArray *cxfData, SumPackage &sumPackage, Vect3<real> *em1_lf,Vect3<real> *em2_lf, Vect3<real> *em1_tf, Vect3<real> *em2_tf);
void Pbcscavec(PeriodicBoundaryFlag jpbc, int nscat, real pyddx, real pzddx, Vect3<real> &a1, Vect3<real> &a2, real theta, real beta, Vect3<real> &xlr, Vect3<real> &ylr, Vect3<real> &zlr, 
	real ak1, Vect3<real> &en0r, Vect3<Complex> &cxe01r, real *orderm, real *ordern, real *thetan, real *phin, Vect3<real> *aksr, Vect3<real> *em1r, Vect3<real> *em2r, 
	Vect3<real> *ensc, Vect3<real> *em1, Vect3<real> *em2);
void Rotate(Vect3<real> &a1, Vect3<real> &a2, real ak1, real beta, real theta, real phi, 
	Vect3<real> &en0r, Vect3<Complex> &cxe01r, Vect3<Complex> &cxe02r, int nscat, Vect3<real> *ensc, 
	Vect3<real> *em1, Vect3<real> *em2, Vect3<real> *aksr, Vect3<real> *em1r, Vect3<real> *em2r);
void Alphadiag(AbstractTarget *currentTarget, Vect3<real> &akr, Matrix *theTensor, Vect3<Complex> &cxe0r, DielectricManager *dielec, int myid);
void Nuller(Complex *cxvec, bool *iocc, int nat);
void Prod3c(Matc3<real> &rm, Vect3<Complex> &cvin, Vect3<Complex> &cvout);
void Prod3v(const Matc3<real> &rm, Vect3<real> *vin, Vect3<real> *vout, int nv);
void Rot2(const Vect3<real> &a, real theta, Matc3<real> &rm);
void Nearfield(AbstractTarget *currentTarget, DDscatParameters *param, const char *cfle1, const char *cfle2, const char *cfleb1, const char *cfleb2, int myid, real akd, Vect3<real> &akr,  
	int ncomp, real aeff, real wave, Complex *cxadia, DielectricManager *dielec, Complex *cxpol1, Complex *cxpol2, Vect3<Complex> &cxe01r, Vect3<Complex> &cxe02r, Array4Stacked<Complex> *cxzw);

#include "Timeit.h"

#endif // __DDSCAT_MAIN_H__
