#include "StdAfx.h"

//
// History of major changes to DDSCAT package.
// Fortran version history removed.
// Copyright (C) 1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,
//               2005,2006,2007,2008,2009,2010,2011,2012,2013 
//               B.T. Draine and P.J. Flatau
// This code is covered by the GNU General Public License.
// Copyright (C) 2012,2013 C++ version, Choliy V.

const char *Version()
{
	static const char *version = "DDscat.C++ 7.3.0++ [01 May 2013].";
	return version;
}

int VersionNum()
{
	static const int versionnum = 730;
	return versionnum;
}