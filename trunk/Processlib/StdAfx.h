// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once


#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers
// Windows Header Files:
#ifdef _WIN32
	#include <windows.h>
	#include <io.h>
#else
	#include <unistd.h>
#endif	
#include <fcntl.h>
#include <sys/stat.h>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include <string>

using namespace std;

#include "CleanDelete.h"

// TODO: reference additional headers your program requires here
