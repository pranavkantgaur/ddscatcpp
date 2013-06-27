#include "StdAfx.h"

#include "Definitions.h"
#include "Scaner.h"

int Scaner::ScanInt(const char *Buf)
{
	return ScanInt(Buf, strlen(Buf));
}

int Scaner::ScanInt(const char *Buf, int len)
{
	if (len < 1) 
		return 0;
	bool bNeg = false;
	int res = 0;
	for(int i=0; i<len; ++i)
	{
		if (Buf[i] == ' ') continue;
		if (Buf[i] == '-') { bNeg = true; continue; }
		if (Buf[i] == '+') continue;
		if (Buf[i] == '\t') throw ("Tab symbol found in Scaner::ScanInt");
		if ((Buf[i] < '0') || (Buf[i] > '9')) break;
		res = 10*res + (Buf[i] - '0');
	}
	return bNeg ? -res : res;
}

int Scaner::ScanIntHex(const char *Buf)
{
	return ScanIntHex(Buf, strlen(Buf));
}

int Scaner::ScanIntHex(const char *buf, int len)
{
	int i, res;
	int bS = 1;
	for(res = i = 0; i<len; i++)
	{
		if (buf[i] == '+') continue;
		if (buf[i] == '.') continue;
		if (buf[i] == '-') { bS = 0; continue; }
		char ch = toupper(buf[i]);
		if (ch == ' ') ch = '0';
		if (!isdigit(ch) && ((ch < 'A') || (ch > 'F'))) break;
		res = 16 * res + ((ch >= 'A') && (ch <= 'F') ? (ch - 'A' + 10) : (ch - '0'));
	}
	return bS ? res : -res;
}

real Scaner::ScanReal(const char *Buf)		// Automatic decision if there is exponent or now
{
	return ScanReal(Buf, strlen(Buf));
}

real Scaner::ScanReal(const char *Buf, int len)
{
	int ia = strcspn(Buf, "EeDd");
	if (ia)
		return ScanExponent(Buf, len);
	else
		return ScanDouble(Buf, len);
}

real Scaner::ScanExponent(const char *Buf)
{
	return ScanExponent(Buf, strlen(Buf));
}

real Scaner::ScanExponent(const char *Buf, int len)
{
	int i; 
	real p = (real)1.;
	real a = (real)0.;
	bool bX = true;
	bool bS = true;
	for(i=0; i<len; i++)
	{
		if((toupper(Buf[i]) == 'E') || (toupper(Buf[i]) == 'D')) break; 
		if(Buf[i] == '+') continue;
		if(Buf[i] == '-') {  bS = false;  continue;  }
		if(Buf[i] == '.') {  bX = false;  continue;  }
		if(Buf[i] == '\t') throw ("Tab symbol found during ScanExponent");
		a = (real)10. * a + ((Buf[i] == ' ') ? 0 : (Buf[i] - '0'));
		if (!bX) p *= 10;
	}
	a /= (real)p;

	bX = true;
	p = 0;
	for(i++; i<len; i++)
	{
		if(Buf[i] == '+') continue;
		if(Buf[i] == '-') { bX = false; continue; }
		if(Buf[i] == '\t') throw ("Tab symbol found during ScanExponent");
		p = 10 * p + ((Buf[i] == ' ') ? 0 : (Buf[i] - '0'));
	}
	p = Pow((real)10., p);
	a = bX ? a*p : a/p;

	return bS ? a : -a;
}

real Scaner::ScanDouble(const char *Buf)
{
	return ScanDouble(Buf, strlen(Buf));
}

real Scaner::ScanDouble(const char *Buf, int len)
{
	bool bX = true;
	bool bS = true;
	real p = (real)1.;
	real a = (real)0.;
	for(int i=0; i<len; i++)
	{
		if(Buf[i] == '+') continue;
		if(Buf[i] == '-') {  bS = false;  continue;  }
		if(Buf[i] == '.') {  bX = false;  continue;  }
		if(Buf[i] == '\t') throw ("Tab symbol found during ScanDouble");
		a = (real)10. * a + ((Buf[i] == ' ') ? 0 : (Buf[i] - '0'));
		if (!bX) p *= 10;
	}
	a /= (real)p;
	return bS ? a : -a;
}

int Scaner::Count(const char *buf, const char key)
{
	int res = 0;
	
	const char *ia = buf;
	while(*ia)
	{
		if (*ia == key)
			++res;
		++ia;
	}
	return res;
}
