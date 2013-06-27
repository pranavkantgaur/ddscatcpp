// TestDDscat.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "ArrayF.h"
#include "Scaner.h"
#include "AbstractFftEngine.h"

typedef bool (*TestFunction) ();
// bool TestGpfaftFft();
bool TestGpfaftFftEngine();

void PrintTest(FILE *file990, const char *Label, Array3F<Complex> &ar, int nx, int ny, int nz);

int main(int argc, const char *argv[])
{
	TestFunction tests[] = { TestGpfaftFftEngine };
	int testNumber = sizeof(tests) / sizeof(tests[0]);

	for(int i=0; i<testNumber; ++i)
	{
		bool bRes = tests[i]();
		printf("Test #%2d finished %s\n", i, bRes ? "Ok" : "Bad");
	}

	return 0;
}

bool TestGpfaftFftEngine()
{
	int nx = 16;
	int ny = 32;
	int nz = 32;
	int num = 8*nx*ny*nz;
	Array3F<Complex> data;
	data.Dimension(2*nx, 2*ny, 2*nz);

	char Buffer[256];
	FILE *file9 = fopen("FftTestdata.txt", "r");
	if (!file9)
		return false;
	for(int a=0; a<num; ++a)
	{
		int xx, yy, zz;
		fgets(Buffer, 255, file9);
		sscanf(Buffer, "%4d%4d%4d", &xx, &yy, &zz);
		real re = Scaner::ScanExponent(Buffer+12, 25);
		real im = Scaner::ScanExponent(Buffer+37, 25);
		data.Value(xx, yy, zz) = Complex(re, im);
	}
	fclose(file9);

	file9 = fopen("in.txt", "w");
	if (!file9)
	{
		fprintf(stderr, "TestGpfaftFftEngine - Cannot open in.txt file.\n");
		return false;
	}
	PrintTest(file9, "Input", data, 2*nx, 2*ny, 2*nz);
	fclose(file9);

	AbstractFftEngine *engine = AbstractFftEngine::GetEngine(FftMethod_GPFAFT);
	engine->DoFFT(data.GetData(), 2*nx, 2*ny, 2*nz, FftForward);

	file9 = fopen("ou.txt", "w");
	if (!file9)
	{
		fprintf(stderr, "TestGpfaftFftEngine - Cannot open on.txt file.\n");
		return false;
	}
	PrintTest(file9, "Output", data, 2*nx, 2*ny, 2*nz);
	fclose(file9);

	return true;
}

void PrintTest(FILE *file990, const char *Label, Array3F<Complex> &ar, int nx, int ny, int nz)
{
	fprintf(file990, "%s\n", Label);
	for(int x=0; x<nx; ++x)
	{
		for(int y=0; y<ny; ++y)
		{
			for(int z=0; z<nz; ++z)
			{
				fprintf(file990, "%4d%4d%4d%25.16e%25.16e\n", x, y, z, ar.Value(x, y, z).re, ar.Value(x, y, z).im); 
			}
		}
	}
	fprintf(file990, "\n");
}

/* **
bool TestGpfaftFft()
{
	int nx = 16;
	int ny = 32;
	int nz = 32;
	int num = 8*nx*ny*nz;
	Array3F<Complex> data;
	data.Allocate(2*nx, 2*ny, 2*nz);

	char Buffer[256];
	FILE *file9 = fopen("FftTestdata.txt", "r");
	if (!file9)
		return false;
	for(int a=0; a<num; ++a)
	{
		int xx, yy, zz;
		fgets(Buffer, 255, file9);
		sscanf(Buffer, "%4d%4d%4d", &xx, &yy, &zz);
		real re = ScanExponent(Buffer+12, 25);
		real im = ScanExponent(Buffer+37, 25);
		data.Value(xx, yy, zz) = Complex(re, im);
	}
	fclose(file9);

	file9 = fopen("in.txt", "w");
	PrintTest(file9, "Input", data, 2*nx, 2*ny, 2*nz);
	fclose(file9);

	Cxfft3n(data, 2*nx, 2*ny, 2*nz, 1);

	file9 = fopen("ou.txt", "w");
	PrintTest(file9, "Output", data, 2*nx, 2*ny, 2*nz);
	fclose(file9);

	return true;
}
** */