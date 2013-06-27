#ifndef __SCANER_H__
#define __SCANER_H__

#include "General.h"

class GENERALLIB_API Scaner
{
public:
	Scaner() { }
	virtual ~Scaner() { }

public:
	static int ScanInt(const char *Buf);
	static int ScanInt(const char *Buf, int len);
	static int ScanIntHex(const char *Buf);
	static int ScanIntHex(const char *buf, int len);
	static real ScanReal(const char *buf);
	static real ScanReal(const char *buf, int len);
	static real ScanExponent(const char *buf);
	static real ScanExponent(const char *buf, int len);
	static real ScanDouble(const char *buf);
	static real ScanDouble(const char *buf, int len);
	static int Count(const char *buf, const char key);
};

#endif // __SCANER_H__