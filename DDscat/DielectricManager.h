#ifndef __DIELECTRICMANAGER_H__
#define __DIELECTRICMANAGER_H__

#include <string>
#include <vector>
#include <map>

#include "Definitions.h"
#include "Complex.h"

class DielectricManager
{
	class Table
	{
	protected:
        int ieps1;
        unsigned int i2;
		real *wva, *e1a, *e2a;
        unsigned int curSize;

	public:
		Table() { wva = e1a = e2a = NULL; curSize = 0; ieps1 = 0; i2 = 0; } 
		virtual ~Table()
		{
			CleanDelete2(wva);
			CleanDelete2(e1a);
			CleanDelete2(e2a);
		}
		bool Load(const char *fileName);
		bool Interpolate(real wave, real &e1, real &e2);
		inline int GetIeps1() { return ieps1; }
		inline real *Wva() { return wva; }
		inline real *E1a() { return e1a; }
		inline real *E2a() { return e2a; }
		inline real &Wva(int index) { return wva[index]; }
		inline real &E1a(int index) { return e1a[index]; }
		inline real &E2a(int index) { return e2a[index]; }

	protected:
		real Parab3(real x, real x1, real x2, real x3, real y1, real y2, real y3);
		real Parab4(real x, real x1, real x2, real x3, real x4, real y1, real y2, real y3, real y4);
	};

	class TableRef
	{
	protected:
		real data[3]; //wl, rn, cn;

	public:
		TableRef(void) { memset(data, 0, 3*sizeof(real)); }
		~TableRef(void) { }

	public:
		inline real GetWl(void) { return data[0]; }
		inline real GetRn(void) { return data[1]; }
		inline real GetCn(void) { return data[2]; }
		inline real GetWlt(void) { return data[0]; }
		inline real GetTabre(void) { return data[1]; }
		inline real GetTabim(void) { return data[2]; }

	public:
		char *FromString(char *Buffer)
		{
			char *ia = strtok(Buffer, "\t ");
			for(int index=0; index<2; ++index)
			{
				data[index] = atof(ia);
				ia = strtok(NULL, "\t ");
			}
			return ia;
		}
	};

	class TableRefExtra : public TableRef
	{
	protected:
		real extraData[6];		// data[0], data[1], extraData[0-2], data[2], extraData[3-5]

	public:
		TableRefExtra(void) { memset(extraData, 0, 6*sizeof(real)); }
		~TableRefExtra(void) { }

	public:
		inline real GetTabret(int index) { return (!index) ? data[1] : extraData[index-1]; }
		inline real GetTabimt(int index) { return (!index) ? data[2] : extraData[index+2]; }

	public:
		char *FromString(char *Buffer)
		{
			int id = 0;
			int ied = 0;
			char *ia = strtok(Buffer, "\t ");
			for(int index=0; index<9; ++index)
			{
				real a = atof(ia);
				if ((index == 0) || (index == 1) || (index == 5))
					data[id++] = a;
				else
					extraData[ied++] = a;	
				ia = strtok(NULL, "\t ");
			}
			return ia;
		}
	};

protected:
	vector<pair<string, int> > *fileNames;
	vector<Table *> *data;
	vector<TableRef *> *waterData, *iceData;
	vector<TableRefExtra *> *iceExtraData;
	Complex *cxeps, *cxrfr;
	real curWave;
    unsigned int curSize;

protected:
	DielectricManager(void);
	virtual ~DielectricManager(void);

public:
	static DielectricManager *item;

public:
	static DielectricManager *GetInstance();
	static void Kill();

public:
    void Allocate(int size);
	void PrepareForWave(real wave);
    int GetFileNameLength(unsigned int index);
	const char *GetFileName(int index);
	vector<string> GetFileNames();
	void AddDataFromFile(int index, const char *filename);
	int CycleFileNames(void);
	void InitFilePool(int size);
    Complex GetCxrfr(unsigned int index);
    Complex GetCxeps(unsigned int index);
	void UseNambient(real namb);
	void Debug(void);
	bool Refwat(int iunit, real xlam, real t, real &rn, real &cn, real &absind, real &abscof);
	bool Refice(int iunit, real xlam, real t, real &rn, real &cn, real &absind, real &abscof);
	bool LoadWaterData(const char *fileName);
	bool LoadIceData(const char *fileName);

public:
	inline int GetNcomp() const { return data->size(); }
	inline real GetCurWave() const { return curWave; }
	inline Complex *GetCxrfr() const { return cxrfr; }
	inline Complex *GetCxeps() const { return cxeps; }

protected:
	real Sumq(real wl, real wlcen, real bet, real del, real gam);
	void DeleteDataTable(vector<TableRef *> *op);
	void DeleteDataTable(vector<TableRefExtra *> *op);
};

#endif // __DIELECTRICMANAGER_H__
