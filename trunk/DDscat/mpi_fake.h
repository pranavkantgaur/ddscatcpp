#ifndef __MPI_FAKE_H__
#define __MPI_FAKE_H__

#include "TargetDefinitions.h"
#include "Complex.h"
#include "Vect3.h"
#include "Vect6.h"
#include "Enumerator.h"
#include "DielectricManager.h"
#include "ArrayF.h"
#include "SumHolder.h"
#include "AbstractTarget.h"
#include "DDscatParameters.h"

typedef int MPI_Comm;

// mpi_fake.f contains dummy subroutines that take the place
// ! of
// !    SUBROUTINE COLSUM
// !    SUBROUTINE SHARE0
// !    SUBROUTINE SHARE1
// !    SUBROUTINE SHARE2
// !    SUBROUTINE MPI_INIT
// !    SUBROUTINE MPI_COMM_RANK
// !    SUBROUTINE MPI_COMM_SIZE
// !    SUBROUTINE MPI_FINALIZE
// !    SUBROUTINE MPI_BARRIER
// ! which are called when MPI is being used.
// ! These dummy routines do nothing when MPI is not being used.

void Colsum(DDscatParameters *param, int myid, SumPackage &sumPackage, real *s1111, real *s1111_1, real *s2121, real *s2121_1, 
	Complex *cx1121, Complex *cx1121_1, real *smori, real *smori_1);

void Share0(int lace, int laxi, int lclm, int lgi, int lpi, int lqi,
	int mxn3, int mxnat, int mxnx, int mxny, int mxnz, int mxpbc, int mxcxsc, int myid, int nat0);

void Share1(AbstractTarget *currentTarget, DDscatParameters *param, 
    const vector<string> &cfleps, real &daeff, Vect3<real> *ensc, Vect3<real> *em1, Vect3<real> *em2,
	bool &WantIobin, bool &ipbc, PeriodicBoundaryFlag &jpbc, int &mxn3, int &mxnat, 
	int &myid, int &nbeth, int &ncomp, int &nori, real *orderm, real *ordern,  
	real &pyd, real &pyddx, real &pzd, real &pzddx, Vect6<real> &iMinmax);

void Share2(DielectricManager *dielec);

// =======================================================================
void MPI_Init(int &nn);
int MPI_Comm_rank(MPI_Comm MPI_Comm_world, int *myid);
int MPI_Comm_size(MPI_Comm MPI_Comm_world, int *numprocs);
int MPI_Finalize();
int MPI_Barrier(MPI_Comm MPI_Comm_world);

#endif // __MPI_FAKE_H__
