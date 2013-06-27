#include "StdAfx.h"

#include "Definitions.h"
#include "Functions.h"
#include "Timeit.h"

real Cpu_time()
{
	real result = (real)0.;
#ifdef _WIN32
	result = (real)clock() / CLOCKS_PER_SEC;
#else
	result = (real)0.;
#endif
	return result;
}

real Timeit(const char *cmsgtm)
{
/* **
Subroutine TIMEIT

Given:
      CMSGTM = string

Returns:
      If odd-numbered call (first, third, etc.): no input/output:
       
      If even-numbered call:
        * Call WRIMSG to print out CMSGTM and elapsed cpu time since previous call
        * Return
          DTIME = elapsed cputime (sec) on
  
This version of timeit uses the system call "etime" available under
Linux, Solaris, and some other bsd-like Unix systems (e.g., Convex unix).

Fortran version history removed.

Copyright (C) 1993,1994,1995,2004,2007 B.T. Draine and P.J. Flatau
Copyright (c) C++ version, 2012, V.Choliy
This code is covered by the GNU General Public License.
** */

	static bool isOdd = true;
	static real t1 = (real)0.;

	real dtime = (real)0.;
	if (isOdd)
	{
		t1 = Cpu_time();
		dtime = t1;
	}
	else
	{
		real t2 = Cpu_time();
		dtime = t2 - t1;	
		if (strcmp(cmsgtm, "NOPRINT") && strcmp(cmsgtm, "noprint"))
		{
			char cmsgnm[256];
			sprintf(cmsgnm, " Timing results for: %s", cmsgtm);
			Wrimsg("Timeit", cmsgnm);
			if (dtime > 1.e4)
				sprintf(cmsgnm, "%9.0f = CPU time (sec)", dtime);
			else
				sprintf(cmsgnm, "%9.3f = CPU time (sec)", dtime);
            Wrimsg("Timeit", cmsgnm);
		}
	}
	isOdd = !isOdd;

	return dtime;
}