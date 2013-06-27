#include "StdAfx.h"

#include "Complex.h"

void Unreduce(Complex *cxv, bool *iocc, int nat, int nat0)
{
/* **
  Given:
   CXV(1:MXN3)   = complex vector in "reduced" form
   IOCC(1:MXNAT) = 0/1 if site is vacant/occupied
   MXN3          = dimensioning information
   MXNAT         = dimensioning information
   NAT           = number of lattice sites
   NAT0          = number of occupied lattice sites

 Returns:
   CXV(1:MXN3)   = complex vector in "natural" form

 history:
 11.08.16 (BTD) created for use in version 7.2.1
 end history
 Copyright (C) 2011
               B.T. Draine and P.J. Flatau
 This code is covered by the GNU General Public License.
** */

	int jocc = 3*(nat0 - 1);
	for(int j=nat-1; j>=0; j--)
	{
		if(iocc[j] == true)
		{
			cxv[3*j  ] = cxv[jocc  ];		
			cxv[jocc  ].clear();
			cxv[3*j+1] = cxv[jocc+1];		
			cxv[jocc+1].clear();
			cxv[3*j+2] = cxv[jocc+2];		
			cxv[jocc+2].clear();
			jocc -= 3;
		}
	}
}