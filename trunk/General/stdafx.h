// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#ifndef __STDAFXGENERAL_H__
#define __STDAFXGENERAL_H__

#include <cmath>
using namespace std;

#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers

// Windows Header Files:
#ifdef _WIN32
	#include <windows.h>
	#include <io.h>
	#include <fcntl.h>
	#include <sys/stat.h>
#else
	#include <unistd.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "CleanDelete.h"

// TODO: reference additional headers your program requires here

#endif // __STDAFXGENERAL_H__