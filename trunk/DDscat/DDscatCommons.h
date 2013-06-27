#ifndef __DDSCAT_COMMONS_H__
#define __DDSCAT_COMMONS_H__

#include "Definitions.h"
#include "Complex.h"
#include "Vect3.h"
#include "ArrayF.h"

class Common0
{
	class Datax
	{
	protected:
		real ak2, ak3, w;
		int ngrid, idipint;

	public:
		Datax() { ak2 = ak3 = w = (real)0.; ngrid = idipint = 0; }
		~Datax() { }

	public:
		inline real &Ak2old() { return ak2; }
		inline real &Ak3old() { return ak3; }
		inline real &Wold() { return w; }
		inline int &Ngrid() { return ngrid; }
		inline int &Idipint() { return idipint; }
		inline void Set(real x) { ak2 = ak3 = w = x; }
		inline void Set(real a2, real a3, real ww) { ak2 = a2; ak3 = a3; w = ww; }
//
		bool Check(real prec, real akd, real ak2, real ak3, real pyd, real pzd)
		{
			real precakd = prec * akd;
			bool bSkip = false;
			bool bCheck = (Fabs(w - akd) < precakd);
			if (bCheck)
			{
				if (pyd == (real)0.)
				{
					if (pzd == (real)0.)
						bSkip = true;
					else
						if (Fabs(ak3 - this->ak3) < precakd)
							bSkip = true;
				}
				else
				{
					if (pzd == (real)0.)
					{
						if (Fabs(ak2 - this->ak2) < precakd)
							bSkip = true;
					}
					else
					{
						if ((Fabs(ak2 - this->ak2) < precakd) && (Fabs(ak3 - this->ak3) < precakd))
							bSkip = true;
					}
				}
			}
			return bSkip;
		}
	};

protected:
	static Common0 *item;
	Common0() {}
	virtual ~Common0() { }
	Datax data[2];
	char cflpar[80];

public:
	static Common0 *GetInstance();
	static void Kill();
	void Init(void);
	void InitReals(void);

public:
	inline Datax *Get(unsigned int i) { return ((i < 2) ? data+i : NULL); }
	inline char *Cflpar(void) { return cflpar; }
};

class Common1
{
protected:
	static Common1 *item;
	Common1() { }
	virtual ~Common1() { }

	Vect3<real> ak_tf;

public:
	static Common1 *GetInstance();
	static void Kill();
	void Init();

public:
	inline Vect3<real> &Ak_tf() { return ak_tf; }
};

class Common4
{
protected:
	static Common4 *item;
	Common4() { }
	virtual ~Common4();

	Array4Stacked<Complex> *cxzw;

public:
	static Common4 *GetInstance();
	static void Kill();
	void Init();

public:
	inline Array4Stacked<Complex> *Cxzw() { return cxzw; }
	void AllocateCxzw(int sX, int sY, int sZ, int sA);
};

class Common10
{
protected:
	static Common10 *item;
	Common10() { }
	virtual ~Common10() { }

	int myid;

public:
	static Common10 *GetInstance();
	static void Kill();
	void Init();

public:
	inline int &Myid() { return myid; }
};

#endif // __DDSCAT_COMMONS_H__