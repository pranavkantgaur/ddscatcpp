#ifndef __DDSCAT_MAINNEW_H__
#define __DDSCAT_MAINNEW_H__

#include "ArrayF.h"
#include "Complex.h"
#include "Vect3.h"
#include "Vect6.h"
#include "Matc3.h"
#include "Matrix.h"
#include "DipoleData.h"
#include "Functions.h"
#include "Enumerator.h"
#include "AbstractTarget.h"
#include "SumHolder.h"
#include "DDscatParameters.h"
#include "FourArray.h"
#include "CleanDelete.h"

#ifdef mpi
	#include "mpif.h"
#else
	#include "mpi_fake.h"
	extern int MPI_Comm_world;
#endif

#ifdef MKL

#else
	#include "Mkl_fake.h"
#endif

#include "Timeit.h"

#endif // __DDSCAT_MAIN_H__
