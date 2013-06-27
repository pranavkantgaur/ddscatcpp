#include "StdAfx.h"

#include "Functions.h"

#ifndef _WIN32
char *strupr(char *Buffer)
{
	const int diff = ('A' - 'a');
	size_t len = strlen(Buffer);
	for(size_t i=0; i<len; ++i)
	{
		if ((Buffer[i] >= 'a') && (Buffer[i] <= 'z'))
			Buffer[i] += diff;
	}
	return Buffer;
}
#endif

// Extracts num real numbers from Buffer according to Format
void ExtractFromInt(const char *Buffer, const char *Format, int num, real *shpar)
{
	int a;
	const char *ia = Buffer;
	for(int i=0; i<num; ++i)
	{
		sscanf(ia, Format, &a);
		shpar[i] = (real)a;
		ia = strchr(ia+1, ' ');
	}
}

// Extracts num real numbers from Buffer according to Format
void ExtractFromReal(const char *Buffer, const char *Format, int num, real *shpar)
{
	real a;
	const char *ia = Buffer;
	for(int i=0; i<num; ++i)
	{
		sscanf(ia, Format, &a);
		shpar[i] = a;
		ia = strchr(ia+1, ' ');
	}
}

// Get nearest 2*3*5 factorizable number
int GetNearestFactor(int num)
{
	const int nf235[] = 
	{
		   1,    2,    3,    4,    5,    6,    8,    9,   10,   12,   15,   16,   18,   20,
		  24,   25,   27,   30,   32,   36,   40,   45,   48,   50,   54,   60,   64,   72,
		  75,   80,   81,   90,   96,  100,  108,  120,  125,  128,  135,  144,  150,  160,
		 162,  180,  192,  200,  216,  225,  240,  243,  250,  256,  270,  288,  300,  320,
		 324,  360,  375,  384,  400,  405,  432,  450,  480,  486,  500,  512,  540,  576,  
		 600,  625,  640,  648,  675,  720,  729,  750,  768,  800,  810,  864,  900,  960,  
		 972, 1000, 1024, 1080, 1125, 1152, 1200, 1215, 1250, 1280, 1296, 1350, 1440, 1458, 
		1500, 1536, 1600, 1620, 1728, 1800, 1875, 1920, 1944, 2000, 2025, 2048, 2160, 2187, 
		2250, 2304, 2400, 2430, 2500, 2560, 2592, 2700, 2880, 2916, 3000, 3072, 3125, 3200, 
		3240, 3375, 3456, 3600, 3645, 3750, 3840, 3888, 4000, 4050, 4096
	};
	const int size = sizeof(nf235) / sizeof(nf235[0]);

	for(int j=0; j<size; ++j)
	{
		if (num <= nf235[j]) 
			return nf235[j];
	}
	return -1;
}

Vect3<int> GetNearestFactoredVector(const Vect3<int> &op)
{
	return Vect3<int>(GetNearestFactor(op.data[0]), GetNearestFactor(op.data[1]), GetNearestFactor(op.data[2]));
}

void Wrimsg(const char *csubrt, const char *cmsgnm)
{
/* **
! Standard procedure for writing messages
! History:
! 96.11.14 (PJF) Remove "getset" and hardwire" ioout"
! 96.11.20 (BTD) change IOOUT to IDVOUT

! Copyright (C) 1993,1996 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
** */

// ! Note: on some systems (e.g., Solaris) IDVOUT=0 generates unbuffered
// ! output to "standard output".  On other systems IDVOUT=0 may not be
// ! valid; then set IDVOUT=6 to get "standard output", probably buffered.

	const char *Format9000 = " >%s %s\n";
	fprintf(stdout, Format9000, csubrt, cmsgnm);
}

void Errmsg(const char *cstats, const char *csubrt, const char *cmsgnm)
{
/* **
Given:
      CSTATS = 'WARNING' or 'FATAL'
      CSUBRT = name of subroutine
      CMESGN = message
Prints a warning message in a standardized way, and STOPs if CSTATS='FATAL'

Copyright (C) 1993,1996,2004,2010 B.T. Draine and P.J. Flatau
This code is covered by the GNU General Public License.

History:
96.11.14 (PJF) Remove "getset" and hardwire "ioerr"
04.05.23 (BTD) cleanup
10.05.08 (BTD) modify output
end history
** */
	const char *Format9000 = "\n >>>>> FATAL ERROR IN PROCEDURE: %s\n";
	const char *Format9010 = " >>>>> %s\n";
	const char *Format9020 = " >>>>> EXECUTION ABORTED \n";
	const char *Format9030 = "\n >>>>> WARNING IN PROCEDURE: %s\n";

	if (!strcmp(cstats, "FATAL") || !strcmp(cstats, "fatal"))
	{
		fprintf(stderr, Format9000, csubrt);
        fprintf(stderr, Format9010, cmsgnm);
        fprintf(stderr, Format9020);
        return;					// exit(0);
	}
	else
	{
		fprintf(stderr, Format9030, csubrt);
		fprintf(stderr, Format9010, cmsgnm);
	}
}

const char *NumberInWords(int number)
{
	const char *text[] = { "One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine", "Ten", "A lot of" };
	if ((number < 1) || (number > 11))
		number = 10;
	return text[number];
}

int ExtractFirstWord(char *Buffer, char fc)
{
	char *ia = NULL;
	while(Buffer[0] == fc)
		ShiftLeft(Buffer);
	ia = strchr(Buffer, fc);
	if (ia)
		*ia = '\0';
	strupr(Buffer);
	return (ia - Buffer);
}

void ShiftLeft(char *Buffer)
{
	char *ia = Buffer;
	do
	{
		*ia = *(ia + 1);
		++ia;
	} while(*ia);
}

void RemoveSymbols(char *Buffer, char sym, char toSymbol)
{
	int i = 0;
	while((Buffer[i] != toSymbol) && (Buffer[i] != '\0'))
	{
		if (Buffer[i] == sym)
		{
			ShiftLeft(Buffer+i);
			continue;
		}
		else
			++i;
	}
}

Complex Scalar(Complex *c, real *d)
{
	return c[0] * d[0] + c[1] * d[1] + c[2] * d[2];
}

Complex Scalar(Complex *c, Complex *d)
{
	return c[0] * d[0] + c[1] * d[1] + c[2] * d[2];
}

void Vector(Complex *c, real *d, Complex *res)
{
	res[0] = c[1] * d[2] - c[2] * d[1]; 
	res[1] = c[2] * d[0] - c[0] * d[2];
	res[2] = c[0] * d[1] - c[1] * d[0];
}

void Vector(Complex *c, Complex *d, Complex *res)
{
	res[0] = c[1] * d[2] - c[2] * d[1]; 
	res[1] = c[2] * d[0] - c[0] * d[2];
	res[2] = c[0] * d[1] - c[1] * d[0];
}

real Ran3(int idum)
{
/* **
Given:
    IDUM = integer seed (always used on first call, but
           disregarded on subsequent calls unless IDUM < 0)
Returns:
    RAN3 = uniform random deviate between 0.0 and 1.0.

Usage:
    On first call, will be initialized with whatever seed
    IDUM is provided.
    On subsequent calls, result does not depend on IDUM 
    unless IDUM<0, in which case IDUM is used to reinitialize
    the sequence.

This routine is based on Knuth's recommendation for a portable
random number generator.  Present Fortran 90 implementation follows
closely the f77 FUNCTION RAN3 described in section 7.1 of the book 
"Numerical Recipes" by W.H. Press, B.P. Flannery, S.A. Teukolsky, and 
W.T. Vetterling (1986; Cambridge University Press).

B.T. Draine, Princeton University Observatory
** */
	const int mbig = 1000000000;
	const int mseed = 161803398;
	const int mz = 0;
	const real fac = 1. / mbig;

	static bool iff = false; 
	static int inext, inextp;
	static int ma[55];

// Initialization:
	if ((idum < 0) || (iff == false))
	{
		iff = true;
		int mj = mseed - abs(idum);
		mj = mj % mbig;
		ma[54] = mj;
		int mk = 1;
		for(int i=0; i<54; ++i)
		{
			int ii = (21 * (i+1)) % 55;
			ma[ii] = mk;
			mk = mj - mk;
			if (mk < mz) mk += mbig;
			mj = ma[ii];
		}
		for(int k=0; k<4; ++k)
		{
			for(int i=0; i<55; ++i)
			{
				ma[i] -= ma[(i+31) % 55];
				if (ma[i] < mz)
					ma[i] += mbig;
			}
		}
		inext = 0;
		inextp = 31;
		idum = 1;
	}
//
// Under normal use (second and subsequent calls), this is the first executed statement:
	++inext;
	if (inext == 55)
		inext = 0;
	++inextp;
	if (inextp == 55)
		inextp = 0;
	int mj = ma[inext] - ma[inextp];
	if (mj < mz)
		mj += mbig;
	ma[inext] = mj;

	return mj * fac;
}
