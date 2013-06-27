#ifndef __CPROD_H__
#define __CPROD_H__

void Diagl(Complex *u, Complex *v, int *);
void Cdiv(Complex *u, Complex *v, Complex *cxa, int n);
void Matvec(Complex *cxx, Complex *cxy, int *);
void Cmatvec(Complex *cxx, Complex *cxy, int *);
void Cprod(Vect3<real> &akr, real gamma, Vect3<real> &dx, FftMethod cmethd, char cwhat, Complex *cxadia, Complex *cxaoff, 
	Complex *cxx, Complex *cxy, Array4Stacked<Complex> *cxzw, bool *iocc, int nat, int nat0, int nat3, real pyd, real pzd);

#endif