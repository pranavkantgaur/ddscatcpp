#include "StdAfx.h"

#include "DielectricManager.h"
#include "Functions.h"

const char *RefIce = "Refice";
const char *RefWat = "Refwat";

class DeleteTableElement
{
public:
	template <typename T>
	void operator() (const T *ptr) const
	{
		delete ptr;
	}
};

DielectricManager *DielectricManager::item = NULL;
DielectricManager::DielectricManager(void)
{
	fileNames = NULL;
	data = NULL;
	cxeps = cxrfr = NULL;
	curSize = 0;
	waterData = NULL;
	iceData = NULL;
	iceExtraData = NULL;
}

DielectricManager::~DielectricManager(void)
{
	if (fileNames)
		delete fileNames;
	if (data)
		delete data;
	fileNames = NULL;
	data = NULL;
	CleanDelete2(cxeps);
	CleanDelete2(cxrfr);
	DeleteDataTable(iceData);
	DeleteDataTable(iceExtraData);
	DeleteDataTable(waterData);
}

DielectricManager *DielectricManager::GetInstance()
{
	if (!item)
		item = new DielectricManager;
	return item;
}

void DielectricManager::Kill()
{
	CleanDelete(item);
}

void DielectricManager::Allocate(int size)
{
	curSize = size;
	cxeps = new Complex[curSize];
	cxrfr = new Complex[curSize];
}

void DielectricManager::InitFilePool(int size)
{
	fileNames = new vector<pair<string, int> >;
	fileNames->resize(size);
	data = new vector<Table *>;
}

void DielectricManager::Debug(void)
{
    for(unsigned int i=0; i<fileNames->size(); ++i)
	{
		pair<string, int> pair = fileNames->at(i);
		fprintf(stdout, "%3d %3d %s", i, pair.second, pair.first.c_str());
		if (pair.second < 0)
		{
			pair = fileNames->at(-pair.second - 1);
			fprintf(stdout, "\t%3d %s", pair.second, pair.first.c_str());
		}
		fprintf(stdout, "\n");
	}
}

void DielectricManager::AddDataFromFile(int index, const char *filename)
{
	const char *MyLabel = "Dielec";
	char Buffer[256];

	int pos = -1;
    for(unsigned int i=0; i<data->size(); ++i)
	{
		if(!strcmp(fileNames->at(i).first.c_str(), filename))
		{
			pos = i;
			break;
		}
	}
	if (pos != -1)
	{
		data->push_back(NULL);
		fileNames->at(index).first = string(" ");
		fileNames->at(index).second = -pos - 1;
		sprintf(Buffer, "File <%s> reassigned into position %d.", filename, pos);
	}
    else
	{
		Table *table = new Table;
		bool bOk = table->Load(filename);
		if (!bOk)
		{
			fprintf(stderr, "Cannot open dielec file %s.\n", filename);
			return;	
		}
		data->push_back(table);
		fileNames->at(index).first = string(filename);
		fileNames->at(index).second = (int)data->size();
		sprintf(Buffer, "File <%s> succesfully loaded.", filename);
	}
	Wrimsg(MyLabel, Buffer);
}

int DielectricManager::CycleFileNames(void)
{
	int indexFrom = 0;
	int indexToto = data->size();
	while(data->size() < curSize)
	{
		data->push_back(NULL);
		fileNames->at(indexToto).first = string(" ");
		if (fileNames->at(indexFrom).second > 0)
			fileNames->at(indexToto).second = -indexFrom - 1;
		else
			fileNames->at(indexToto).second = fileNames->at(indexFrom).second;
		++indexFrom;
		++indexToto;
	}
	return 0;
}

bool DielectricManager::Table::Load(const char *fileName)
{
	FILE *file = fopen(fileName, "r");
	if (!file)
		return false;
//
// Read header line:
	char Buffer[256];
	fgets(Buffer, 255, file);
	while((Buffer[strlen(Buffer)-1] == 0x0d) || (Buffer[strlen(Buffer)-1] == 0x0a))
		Buffer[strlen(Buffer)-1] = '\0';
	Wrimsg("Header:", Buffer);
//
//  Read line specifying columns for wavelength,Re(n),Im(n),Re(eps),Im(eps
	int iwv, iren, iimn, ieps2;
	fgets(Buffer, 255, file);
	sscanf(Buffer, "%d%d%d%d%d", &iwv, &iren, &iimn, &ieps1, &ieps2);
//
	int icol = iwv;
	if (iren > icol) icol = iren;
	if (iimn > icol) icol = iimn;
	if (ieps1 > icol) icol = ieps1;
	if (ieps2 > icol) icol = ieps2;
//
// Skip header line:
	char *ia = NULL;
	real xin[10];
	fgets(Buffer, 255, file);

	curSize = 16;
	int stepSize = 16;
	wva = (real *)malloc(curSize * sizeof(real));
	e1a = (real *)malloc(curSize * sizeof(real));
	e2a = (real *)malloc(curSize * sizeof(real));

    unsigned int nwavt = 0;
	while(!feof(file))
	{
		ia = fgets(Buffer, 255, file);
		if (!ia)
			break;
		ia = Buffer;
		for(int jj=0; jj<icol; ++jj)
		{
			ia = strtok(ia, " \t");
			xin[jj] = (real)atof(ia);
			ia = NULL;
		}
		wva[nwavt] = xin[iwv-1];
		if (ieps1 > 0)
		{
			e1a[nwavt] = xin[ieps1-1];
			e2a[nwavt] = xin[ieps2-1];
		}
		else
		{
			e1a[nwavt] = xin[iren-1];
			e2a[nwavt] = xin[iimn-1];
		}
		nwavt++;
		if (nwavt >= curSize)
		{
			curSize += stepSize;
			wva = (real *)realloc(wva, curSize * sizeof(real));
			e1a = (real *)realloc(e1a, curSize * sizeof(real));
			e2a = (real *)realloc(e2a, curSize * sizeof(real));
		}
	}
//
// Everything is read, fix size to only necessary number of lines
	wva = (real *)realloc(wva, nwavt * sizeof(real));
	e1a = (real *)realloc(e1a, nwavt * sizeof(real));
	e2a = (real *)realloc(e2a, nwavt * sizeof(real));
	curSize = nwavt;
	fclose(file);
//
	return true;
}

void DielectricManager::PrepareForWave(real wave)
{
	real e1, e2;
	vector<Table *>::iterator ia, ib;
	int i = 0;
	for(ia=data->begin(), ib=data->end(); ia!=ib; ++ia, ++i)
	{
		if (*ia)
		{
			bool bRes = (*ia)->Interpolate(wave, e1, e2);
			if (!bRes)
				Wrimsg("Dielec interpolator", "Bad result");
			if ((*ia)->GetIeps1() > 0)
				cxeps[i] = Complex(e1, e2);
			else
				cxeps[i] = Complex(e1 * e1 - e2 * e2, (real)2. * e1 * e2);
		}
		else
		{
			cxeps[i] = cxeps[fileNames->at(-fileNames->at(i).second - 1).second - 1];
		}
		cxrfr[i] = cxeps[i].sqrt();
	}
	curWave = wave;
}

Complex DielectricManager::GetCxrfr(unsigned int index)
{
	if (index < curSize)
		return cxrfr[index];
	else
		return Complex();
}

Complex DielectricManager::GetCxeps(unsigned int index)
{
	if (index < curSize)
		return cxeps[index];
	else
		return Complex();
}

const char *DielectricManager::GetFileName(int index) 
{ 
    if (index < (int)fileNames->size())
	{
		if (fileNames->at(index).second > 0)
			return fileNames->at(index).first.c_str();
		else
			return fileNames->at(-fileNames->at(index).second - 1).first.c_str();
	}
	else
		return NULL;
}

int DielectricManager::GetFileNameLength(unsigned int index)
{
    if (index < fileNames->size())
	{
		if (fileNames->at(index).second > 0)
			return fileNames->at(index).second;
		else
			return fileNames->at(-fileNames->at(index).second - 1).second;
	}
	else
		return ihuge_;
}

void DielectricManager::UseNambient(real namb)
{
    for(unsigned int i=0; i<curSize; ++i)
	{
		cxrfr[i] /= namb;
		cxeps[i] /= (namb * namb);
	}
}

vector<string> DielectricManager::GetFileNames(void)
{
	vector<string> res;
	vector<pair<string, int> >::iterator ia, ib;
	for(ia=fileNames->begin(), ib=fileNames->end(); ia!=ib; ++ia)
	{
		if (ia->second > 0)
			res.push_back(ia->first);
		else
			res.push_back(fileNames->at(-ia->second - 1).first);
	}
	return res;
}

bool DielectricManager::Table::Interpolate(real x, real &y, real &z)
{
// void Interp(real x, real &y, real &z, real *xa -> wva, real *ya -> e1a, real *za -> e2a, int &ntab)
/* **
Given array of tabulated values XA(1-NTAB), YA(1-NTAB), ZA(1-NTAB), and given independent variable X, this routine interpolates for two dependent variables: Y and Z.
Uses only parabolic interpolation for safety. 
When called with X outside range of tabulated values XA, returns extrapolated values Y and Z but issues warning statement.
Note: XA values must be either monotonically increasing OR monotonically decreasing.
B.T.Draine, Princeton Univ. Observatory

History records of Fortran versions removed.

Copyright (C) 1993, B.T. Draine and P.J. Flatau
Copyright (C) 2013, C++ versions, V.Choliy
This code is covered by the GNU General Public License.
** */
	real sgn, x1, x2, x3, x4, y1, y2, y3, y4;
	int i1, i3, i4;

	if  (i2 >= curSize - 1) 
		i2 = curSize - 2;
	int inc = 0;
	sgn = (real)1.;
//
// Check whether X is increasing or decreasing:
	if (wva[0] >= wva[curSize - 1]) 
		sgn = -(real)1.;
//
// Check whether outside table limits
	const char *Format6990 = "Warning from INTERP: outside table limits for X=%12.5lf\n";
	if ((sgn*(x - wva[0])) < (real)0.)
	{
		i2 = 1;
		fprintf(stderr, Format6990, x);
        goto l4600;
	}
	if (sgn*(x - wva[curSize - 1]) > (real)0.)
	{
		i2 = curSize - 2;
		fprintf(stderr, Format6990, x);
        goto l4600;
	}
//
// X is within table limits.  Find I2
l1000:
	if (sgn*(wva[i2] - x) > (real)0.) goto l4000;
	if (sgn*(wva[i2] - x) < (real)0.) goto l2000;
	y = e1a[i2];
	z = e2a[i2];
	return true;

l2000:
	if (inc >= (real)0.) goto l2200;
	if (i2 + 3 - curSize <= (real)0.) goto l4700;
	if (i2 + 3 - curSize  > (real)0.) goto l2500;

l2200:
	inc = 1;
	++i2;
	if (i2 <= curSize - 2) goto l1000;

l2500:
	i2 = curSize - 2;
	goto l4600;

l4000:  
	if (inc <= (real)0.) goto l4200;
	--i2;
	if (i2 - 1 < (real)0.) 
		goto l4500;
	else 
		goto l4700;

l4200:
	inc = -1;
	--i2;
	if (i2 >= 1) 
		goto l1000;

l4500:
	i2 = 1;

l4600:
	i1 = i2 - 1;
	i3 = i2 + 1;
	x1 = wva[i1];
	x2 = wva[i2];
	x3 = wva[i3];
	y1 = e1a[i1];
	y2 = e1a[i2];
	y3 = e1a[i3];
	y = Parab3(x, x1, x2, x3, y1, y2, y3);
	y1 = e2a[i1];
	y2 = e2a[i2];
	y3 = e2a[i3];
	z = Parab3(x, x1, x2, x3, y1, y2, y3);
	return true;

l4700:
	i1 = i2 - 1;
	i3 = i2 + 1;
	i4 = i2 + 2;
	x1 = wva[i1];
	x2 = wva[i2];
	x3 = wva[i3];
	x4 = wva[i4];
	y1 = e1a[i1];
	y2 = e1a[i2];
	y3 = e1a[i3];
	y4 = e1a[i4];
	y = Parab4(x, x1, x2, x3, x4, y1, y2, y3, y4);
	y1 = e2a[i1];
	y2 = e2a[i2];
	y3 = e2a[i3];
	y4 = e2a[i4];
	z = Parab4(x, x1, x2, x3, x4, y1, y2, y3, y4);
	return true;
}

real DielectricManager::Table::Parab3(real x, real x1, real x2, real x3, real y1, real y2, real y3)
{
/* **
Subroutine PARAB3 does parabolic interpolation, with parabola constrained to fit (x1,y1),(x2,y2),(x3,y3) exactly.
B.T.Draine, Institute for Advanced Study, March 1980.

Copyright (C) 1993, B.T. Draine and P.J. Flatau
This code is covered by the GNU General Public License.
** */
	real a = (y3 - y2 - (x3 - x2) * (y1 - y2) / (x1 - x2)) / ((x3 - x2) * (x3 - x1));
	real b = (y1 - y2) / (x1 - x2) - a * (x1 - x2);

	return (a * (x - x2) + b) * (x - x2) + y2;
}

real DielectricManager::Table::Parab4(real x, real x1, real x2, real x3, real x4, real y1, real y2, real y3, real y4)
{
/* **
Subroutine PARAB4 is designed to do parabolic interpolation, with parabola constrained to match (x2,y2) and (x3,y3) exactly, and to minimize sum of squared deviations from (x1,y1) and (x4,y4).
It is assumed that x1.lt.x2.le.x.lt.x3.lt.x4 
Fit parabola Y=A*(X-X2)**2+B*(X-X2)+Y2
B.T.Draine, Institute for Advanced Study, March 1980.

Copyright (C) 1993, B.T. Draine and P.J. Flatau
This code is covered by the GNU General Public License.
** */
	real a = ((x1-x2)*(y3-y2)/(x3-x2)+y2-y1)*(x1-x2)*(x1-x3) + ((x4-x2)*(y3-y2)/(x3-x2)+y2-y4)*(x4-x2)*(x4-x3);
	a = -a / (((x1-x2)*(x1-x3))*((x1-x2)*(x1-x3)) + ((x4-x2)*(x4-x3))*((x4-x2)*(x4-x3)));
	real b = (y3-y2)/(x3-x2) - a*(x3-x2);

	return (a*(x-x2)+b) * (x-x2) + y2;
}

void DielectricManager::DeleteDataTable(vector<TableRef *> *op)
{
	if (op)
		for_each(op->begin(), op->end(), DeleteTableElement());
	CleanDelete(op);
}

void DielectricManager::DeleteDataTable(vector<TableRefExtra *> *op)
{
	if (op)
		for_each(op->begin(), op->end(), DeleteTableElement());
	CleanDelete(op);
}

bool DielectricManager::LoadWaterData(const char *fileName)
{
	if (!waterData)
	{
		FILE *file = fopen(fileName, "r");
		if (!file)
			return false;

		int num;
		char Buffer[256];
		fgets(Buffer, 255, file);
		sscanf(Buffer, "%d", &num);
		waterData = new vector<TableRef *>;
		waterData->reserve(num);
		for(int i=0; i<num; ++i)
		{
			fgets(Buffer, 255, file);
			TableRef *tableref = new TableRef;
			tableref->FromString(Buffer);
			waterData->push_back(tableref);
		}
		fclose(file);
	}
	return true;
}

bool DielectricManager::LoadIceData(const char *fileName)
{
	if (!iceData)
	{
		FILE *file = fopen(fileName, "r");
		if (!file)
			return false;

		int i, num;
		char Buffer[256];
		fgets(Buffer, 255, file);
		sscanf(Buffer, "%d", &num);
		iceData = new vector<TableRef *>;
		iceData->reserve(num);
		for(i=0; i<num; ++i)
		{
			fgets(Buffer, 255, file);
			TableRef *tableref = new TableRef;
			tableref->FromString(Buffer);
			iceData->push_back(tableref);
		}
		fgets(Buffer, 255, file);
		sscanf(Buffer, "%d", &num);
		iceExtraData = new vector<TableRefExtra *>;
		iceExtraData->reserve(num);
		for(i=0; i<num; ++i)
		{
			fgets(Buffer, 255, file);
			TableRefExtra *tableref = new TableRefExtra;
			tableref->FromString(Buffer);
			iceExtraData->push_back(tableref);
		}
		fclose(file);
	}
	return true;
}

// FUNCTION FOR TREATING ABSORPTION BANDS NOT CONSIDERED IN THE DEBYE THEORY
real DielectricManager::Sumq(real wl, real wlcen, real bet, real del, real gam)
{
	return bet * Exp(Pow(-Fabs(Log10(wl / wlcen) / del), gam));
}

bool DielectricManager::Refwat(int iunit, real xlam, real t, real &rn, real &cn, real &absind, real &abscof)
{
/* **
     DEFINES WAVELENGTH DEPENDENT COMPLEX INDEX OF REFRACTION FOR WATER
     ALLOWABLE WAVELENGTH RANGE EXTENDS FROM .2 MICRONS TO 10 CM
     TEMPERATURE DEPENDENCE ONLY CONSIDERED BEYOND 0.1 CM

     ERIC A. SMITH
     DEPT OF ATMOSPHERIC SCIENCE
     COLORADO STATE UNIVERSITY
     FORT COLLINS,CO  80523
     TEL   303-491-8533

     REFERENCES
     0.2 UM - 0.69 UM
     HALE,G., AND M. QUERRY,1972.
     OPTICAL CONSTANTS OF WATER IN THE 200 NM TO 200 UM WAVELENGTH REGI
     APPLIED OPTICS,12,3,555-563.

     0.69 UM - 2.0 UM
     PALMER,K.F., AND D. WILLIAMS,1974.
     OPTICAL PROPERTIES OF WATER IN THE NEAR INFRARED.
     JOURNAL OF THE OPTICAL SOCIETY OF AMERICA,64,8,1107-1110.

     2.0 UM - 1000.0 UM
     DOWNING,H.D., AND D. WILLIAMS,1975.
     OPTICAL CONSTANTS OF WATER IN THE INFRARED.
     JOURNAL OF GEOPHYSICAL REVIEW,80,12,1656-1661.

     1.0 MM - 10.0 CM
     RAY,P.S.,1972.
     BROADBAND COMPLEX REFRACTIVE INDICES OF ICE AND WATER.
     APPLIED OPTICS,11,8,1836-1844.

     INPUT PARAMETERS
     IUNIT = 0 FOR WAVELENGTH SPECIFIED IN MICRONS
           = 1 FOR WAVELENGTH SPECIFIED IN MILLIMETERS
           = 2 FOR WAVELENGTH SPECIFIED IN CENTIMETERS
           = 3 FOR WAVELENGTH SPECIFIED IN INVERSE CENTIMETERS ( WAVE N
     XLAM = WAVELENGTH ( MICRONS OR MM OR CM OR CM**-1 )
     T = TEMPERATURE ( DEGREES KELVIN )

     OUTPUT PARAMETERS

     RN = REAL PORTION ( SCATTERING )
     CN = COMPLEX PORTION ( ABSORPTION )
     ABSIND = ABSORPTIVE INDEX ( CN/RN )
     ABSCOF = ABORPTION COEFFICIENT ( 4*PI*CN/XLAM )
** */
	const real cutwat = (real)1000.;
	const int numwat = waterData->size();

// Zero parameters
	rn = cn = absind = abscof = (real)0.;

// Convert wavelength to microns
	real wl = xlam;
	switch(iunit)
	{
	case 1:
		wl *= (real)1000.;
		break;

	case 2:
		wl *= (real)10000.;
		break;

	case 3:
		wl = (real)10000. / wl;
		break;

	default:
		return false;
	}

// Region from 0.2 micron to 1000.0 micron  -  table lookup
	if (wl <= cutwat)
	{
		int i1 = numwat - 1;
		for(int i=1; i<numwat; ++i)
		{
			if (wl > waterData->at(i)->GetWl())
				continue;
			i1 = i - 1;
			break;
		}
		TableRef *t1 = waterData->at(i1);
		TableRef *t2 = waterData->at(i1 + 1);
        real fac = (wl - t1->GetWl()) / (t2->GetWl() - t1->GetWl());
		rn = t1->GetRn() + fac * (t2->GetRn() - t1->GetRn());
		cn = t1->GetCn() + fac * (t2->GetCn() - t1->GetCn());
		absind = cn / rn;
		abscof = (real)4. * Pi * cn / wl;
		return true;
	}
//
//  Region from 0.1 cm to 10 cm
//  Extension of debye theorey based on the work of
//  Cole,K.S.,and R.H.Cole,1941.Jour.Chem.Phys.,9,p 341.
//  Define temperature terms and wavelength in cm
	real tc =  t - (real)273.15;
	real t1 = tc + (real)273.0;
	real t2 = tc - (real)25.0;
	real xl = wl / (real)10000.0;
//
//  Define frequency independent conductivity(sigma) and  spread parameter(alpha)
//  in classical debye theory these terms are zero
//  sigma given by Saxton,J.A.,1949.Wireless Engineer,26,P 288.
//  alpha given by Ray ( equation 7b )
	real sigma = (real)12.5664e8;
	real alpha = (real)(-16.8129 / t1 + 0.0609265);
//
//  Define static dielectric constant(es) - ray eqn 4
//  High frequency dielectric constant(e00) - ray eqn 7a
//  relaxtion wavelength in cm(xlams) - ray eqn 7c
//  temperature dependence of es given by
//  Wyman,J.,and E.N.Ingalls,1938.Jour.Am.Chem.Soc.,60,P 1182.
	real es = 78.54 * (1.0 + (-4.579e-3 + (1.19e-5 - 2.8e-8 * t2) * t2) * t2);
	real e00 = 5.27137 + (0.0216474 - 0.00131198 * tc) * tc;
	real xlams = 0.00033836 * Exp(2513.98 / t1);
//
//  Calculate expressions used for dielectric constant
	real term = Pi * alpha / 2;
	real sint = Sin(term);
	real cost = Cos(term);
	real xlrat = xlams / xl;
	real powtrm = Pow(xlrat, 1.-alpha);
	real denom = 1.0 + 2 * powtrm * sint + Pow(xlrat, 2.0*(1-alpha));
//
//  Calculation of dielectric constant
//
	real er = e00 + (es-e00) * (1.0 + powtrm * sint) / denom;
//
//  Imaginary part or loss term - ray eqn 6
	real ei = (sigma * xl / 18.8496e10) + (es-e00) * powtrm * cost / denom;
//
//  Complex permittivity
	Complex e(er, -ei);
// 
//  Complex index of refraction - ray eqn 1
	Complex m = e.sqrt();
	rn =  m.re;
	cn = -m.im;
//
//  Correction to imaginary index to account for the
//  remaining absorption bands - ray eqn 8(table 2)
	if (wl <= (real)3000.)
		cn = cn + Sumq(wl, 17.0,0.39,0.45,1.3) + Sumq(wl, 62.0,0.41,0.35,1.7) + Sumq(wl,300.0,0.25,0.47,3.0);
//
// Absorptive quanities
	absind = cn / rn;
	abscof = (real)4. * Pi * cn / wl;
	return true;
}

bool DielectricManager::Refice(int iunit, real xlam, real t, real &rn, real &cn, real &absind, real &abscof)
{
/* **
     DEFINES WAVELENGTH DEPENDENT COMPLEX INDEX OF REFRACTION FOR ICE.
     ALLOWABLE WAVELENGTH RANGE EXTENDS FROM 0.045 MICRONS TO 8.6 METER
     TEMPERATURE DEPENDENCE ONLY CONSIDERED BEYOND 167 MICRONS.

     INTERPOLATION IS DONE     RN  VS. LOG(XLAM)
                               RN  VS.        T
                           LOG(CN) VS. LOG(XLAM)
                           LOG(CN) VS.        T
                           
     STEPHEN G. WARREN - 1983
     DEPT. OF ATMOSPHERIC SCIENCES
     UNIVERSITY OF WASHINGTON
     SEATTLE, WA  98195

     BASED ON
        WARREN,S.G.,1984.
        OPTICAL CONSTANTS OF ICE FROM THE ULTRAVIOLET TO THE MICROWAVE.
        APPLIED OPTICS,23,1206-1225

     INPUT PARAMETERS
     IUNIT = 0 FOR WAVELENGTH SPECIFIED IN MICRONS
           = 1 FOR WAVELENGTH SPECIFIED IN MILLIMETERS
           = 2 FOR WAVELENGTH SPECIFIED IN CENTIMETERS
           = 3 FOR WAVELENGTH SPECIFIED IN INVERSE CENTIMETERS ( WAVE N
     XLAM = WAVELENGTH ( MICRONS OR MM OR CM OR CM**-1 )
     T = TEMPERATURE ( DEGREES KELVIN )

     OUTPUT PARAMETERS
     RN = REAL PORTION ( SCATTERING )
     CN = COMPLEX PORTION ( ABSORPTION )
     ABSIND = ABSORPTIVE INDEX ( CN/RN )
     ABSCOF = ABORPTION COEFFICIENT ( 4*PI*CN/XLAM )

     REFERENCE TEMPERATURES ARE -1.0,-5.0,-20.0, AND -60.0 DEG CENTIGRA
** */
	const real temref[4] = { (real)272.16, (real)268.16, (real)253.16, (real)213.16 };
	const real wlmin = (real)0.045;
	const real wlmax = (real)8.6e6;
	const real cutice = (real)167.;
	const int nwl = iceData->size();
	const int nwlt = iceExtraData->size();
//
// Zero parameters
	rn = cn = absind = abscof = (real)0.;
//
// Convert wavelength to microns
	real alam = xlam;
	switch(iunit)
	{
	case 1:
		alam *= (real)1000.;
		break;

	case 2:
		alam *= (real)10000.;
		break;

	case 3:
		alam = (real)10000. / alam;
		break;
	}
	if ((alam < wlmin) || (alam > wlmax))
		return false;
//
	int j;
	int i = -1;
	if (alam < cutice)
	{
// Region from 0.045 microns to 167.0 microns - no temperature depend
		i = -1;
		for(j=1; j<nwl; ++j)
		{
			i = j;
			if (alam < iceData->at(j)->GetWl())
				break;
		}
		real x1 = Log(iceData->at(i-1)->GetWl());
		real x2 = Log(iceData->at(i)->GetWl());
		real y1 = iceData->at(i-1)->GetTabre();
		real y2 = iceData->at(i)->GetTabre();
		real x = Log(alam);
		real y = ((x-x1)*(y2-y1) / (x2-x1)) + y1;
		rn = y;
		y1 = Log(Fabs(iceData->at(i-1)->GetTabim()));
		y2 = Log(Fabs(iceData->at(i)->GetTabim()));
		y = ((x-x1)*(y2-y1) / (x2-x1)) + y1;
		cn = Exp(y);
		absind = cn / rn;
		abscof = (real)4. * Pi * cn / alam;
		return true;
	}
//
// Region from 167.0 microns to 8.6 meters - temperature dependence
	real tk = t;
	if (tk > temref[0]) tk = temref[0];
	if (tk < temref[3]) tk = temref[3];
	i = 0;
	for(j=1; j<4; ++j)
	{
		i = j;
		if (tk >= temref[j])
			break;
	}
	int lt1 = i;
	int lt2 = i - 1;
	i = 0;
	for(j=1; j<nwlt; ++j)
	{
		i = j;
		if (alam <= iceExtraData->at(j)->GetWlt())
			break;
	}
	real x1 = Log(iceExtraData->at(i-1)->GetWlt());
	real x2 = Log(iceExtraData->at(i)->GetWlt());
	real y1 = iceExtraData->at(i-1)->GetTabret(lt1);
	real y2 = iceExtraData->at(i)->GetTabret(lt1);
	real x = Log(alam);
	real ylo = ((x-x1) * (y2-y1) / (x2-x1)) + y1;
	y1 = iceExtraData->at(i-1)->GetTabret(lt2);
	y2 = iceExtraData->at(i)->GetTabret(lt2);
	real yhi = ((x-x1)*(y2-y1) / (x2-x1)) + y1;
	real t1 = temref[lt1];
	real t2 = temref[lt2];
	real y = ((tk-t1)*(yhi-ylo) / (t2-t1)) + ylo;
	rn = y;
	y1 = Log(Fabs(iceExtraData->at(i-1)->GetTabimt(lt1)));
	y2 = Log(Fabs(iceExtraData->at(i)->GetTabimt(lt1)));
	ylo = ((x-x1)*(y2-y1) / (x2-x1)) + y1;
	y1 = Log(Fabs(iceExtraData->at(i-1)->GetTabimt(lt2)));
	y2 = Log(Fabs(iceExtraData->at(i)->GetTabimt(lt2)));
	yhi = ((x-x1)*(y2-y1) / (x2-x1)) + y1;
	y = ((tk-t1)*(yhi-ylo) / (t2-t1)) + ylo;
	cn = Exp(y);
//
// Absorptive quanities
	absind = cn / rn;
	abscof = (real)4. * Pi * cn / alam;
	return true;
}
