#ifndef __VTR_H__
#define __VTR_H__

#include <cstdio>
#include <cstring>
#include "Vtrlib.h"
#include "Definitions.h"
#include "ArrayF.h"
#include "Complex.h"

class VTRLIB_API Vtr
{
	static const char *outFormat;
	static const char *endDataArray;
	static const char *dataArrayFormat;
	static const char *xmlHeader;

	class VTRLIB_API VtrFileHandle
	{
	protected:
		FILE *handle;
		char *prefix;
        unsigned int currentTab;
		char *tabs;

	public:
        int unit, counter, restart;
        unsigned int nxx, nyy, nzz;
		bool first, isOpened;

	public:
		VtrFileHandle();
		virtual ~VtrFileHandle();

	public:
		inline char *GetPrefix(void) const { return prefix; }
		inline FILE *GetHandle(void) const { return handle; }
		inline void SetHandle(FILE *hd) { handle = hd; }
		
	public:
		void CloseHandle();
		bool OpenHandle(const char *f, const char *mode);
		void SetPrefix(const char *data);
		void Clear(void);
		void Fprintf(const char *Format, ...);

	protected:
		void Tabify(void);
		void Tabify(bool dir);
	};

protected:
	VtrFileHandle fd;
	int iproc, nb_procs;

public:
	Vtr(void);
	virtual ~Vtr(void);

protected:
	void HandleError(const char *name, const char *message)
	{
		fprintf(stderr, "\n    *** Error *** %s: %s\n", name, message);
	}
	void HandleWarning(const char *name, const char *message)
	{
		fprintf(stderr, "\n   *** Warning *** %s: %s\n", name, message);
	}
	void HandleInfo(const char *name, const char *message)
	{
		fprintf(stderr, "\n   *** Info *** %s: %s\n", name, message);
	}
	void WriteDataArray(real *x, int ni, const char *Name, int size);
	void WriteDataVector(real **x, const char *Prefix, const char *Name);
	void WriteDataVector(real ***x, const char *Prefix, const char *Name);
	void WriteDataArray(int *x, int ni, const char *Name, int size);
	void WriteDataVector(int **x, const char *Prefix, const char *Name);
	void WriteDataVector(int ***x, const char *Prefix, const char *Name);
	void PrepareFormat(Array3F<real> &a, char *Format, int num);
	void PrepareFormat(Array3F<int> &a, char *Format, int num);
	void PrepareFormat(real a, char *Format, int num);
	void PrepareFormat(int a, char *Format, int num);

public:
	void CloseFile();
	void CollectFile();
	void OpenFile(const char *prefix, int proc_rank = -1, int num_procs = -1, int restart = -1);

public:
template <typename T>
    void WriteMesh2d(T *x, T *y, unsigned int xSize, unsigned int ySize);
template <typename T>
    void WriteMesh3d(T *x, T *y, T *z, unsigned int xSize, unsigned int ySize, unsigned int zSize);
template <typename T>
    void WriteVector2d(const char *name, T **vx, T **vy, unsigned int xSize, unsigned int ySize);
template <typename T>
    void WriteVector3d(const char *name, T ***vx, T ***vy, T ***vz, unsigned int xSize, unsigned int ySize, unsigned int zSize);
template <typename T>
    void WriteVector3d(const char *name, Array3F<T> &vx, Array3F<T> &vy, Array3F<T> &vz, unsigned int xSize, unsigned int ySize, unsigned int zSize);
template <typename T>
    void WriteScalar2d(const char *name, T **field, unsigned int xSize, unsigned int ySize);
template <typename T>
    void WriteScalar3d(const char *name, T ***field, unsigned int xSize, unsigned int ySize, unsigned int zSize);
template <typename T>
	void WriteScalar3d(const char *name, Array3F<T> &array);

public:
	inline VtrFileHandle &FileHandle() { return fd; }

// TODO
//	void WriteMesh2d(Complex *x, Complex *y, int xSize, int ySize);
//	void WriteMesh3d(Complex *x, Complex *y, Complex *z, int xSize, int ySize, int zSize);
//	void WriteVector2d(char *name, Complex **vx, int xSize1, int ySize1, Complex **vy, int xSize2, int ySize2);
//	void WriteVector3d(char *name, Complex ***vx, int xSize1, int ySize1, int zSize1, Complex ***vy, int xSize2, int ySize2, int zSize2, Complex ***vz, int xSize3, int ySize3, int zSize3);
//	void WriteScalar2d(char *name, Complex **field, int xSize, int ySize);
//	void WriteScalar3d(char *name, Complex ***field, int xSize, int ySize, int zSize);
//	void WriteVector3d(char *name, Array3F<Complex> &array);
};

/* **
VTR_mod.f90 --- XML VTR ASCII data file

Auteur          : Jalel Chergui (LIMSI-CNRS) <Jalel.Chergui@limsi.fr>
Date            : Wed Jul 26 14:36:52 2006
Dern. mod. par  : Jalel Chergui (LIMSI-CNRS) <Jalel.Chergui@limsi.fr>
Dern. mod. le   : Wed Sep 16 14:36:29 2009
** */

template <typename T>
void Vtr::WriteMesh2d(T *x, T *y, unsigned int xSize, unsigned int ySize)
{
	real z[] = { (real)0. };

	fd.nxx = xSize;
	fd.nyy = ySize;
    fd.Fprintf("<RectilinearGrid WholeExtent=\"1 %d 1 %d 1 1 \">\n", fd.nxx, fd.nyy);
    fd.Fprintf("<Piece Extent=\"1 %d 1 %d 1 1 \">\n", fd.nxx, fd.nyy);
    fd.Fprintf("<Coordinates>\n");
	WriteDataArray(x, fd.nxx, "X_COORDINATES", 1);
	WriteDataArray(y, fd.nyy, "Y_COORDINATES", 1);
	WriteDataArray(z,      1, "Z_COORDINATES", 1);
    fd.Fprintf("</Coordinates>\n");
    fd.Fprintf("<PointData>\n");
}

template <typename T>
void Vtr::WriteMesh3d(T *x, T *y, T *z, unsigned int xSize, unsigned int ySize, unsigned int zSize)
{
    fd.nxx = xSize;
	fd.nyy = ySize;
	fd.nzz = zSize;
    fd.Fprintf("<RectilinearGrid WholeExtent=\"1 %d 1 %d 1 %d\">\n", fd.nxx, fd.nyy, fd.nzz);
    fd.Fprintf("<Piece Extent=\"1 %d 1 %d 1 %d\">\n", fd.nxx, fd.nyy, fd.nzz);
    fd.Fprintf("<Coordinates>\n");
	WriteDataArray(x, fd.nxx, "X_COORDINATES", 1);
	WriteDataArray(y, fd.nyy, "Y_COORDINATES", 1);
	WriteDataArray(z, fd.nzz, "Z_COORDINATES", 1);
    fd.Fprintf("</Coordinates>\n");
	fd.Fprintf("<PointData>\n");
}

//
// It is implied that both arrays have identical dimensions
template <typename T>
void Vtr::WriteVector2d(const char *name, T **vx, T **vy, unsigned int xSize, unsigned int ySize)
{
	if ((xSize != fd.nxx) || (ySize != fd.nyy)) 
		HandleWarning("Vtr::WriteVector2d", "Incompatible component and mesh sizes.");
//
	char Format[256];
    PrepareFormat(*vx, Format, 3);
    fd.Fprintf(dataArrayFormat, name, 3);
    for(unsigned int j=0; j<fd.nyy; ++j)
        for(unsigned int i=0; i<fd.nxx; ++i)
			fd.Fprintf(Format, vx[i][j], vy[i][j], (real)0.);
	fd.Fprintf(endDataArray);
//
	WriteDataVector(vx, "X", name);
	WriteDataVector(vy, "Y", name);
}

//
// It is implied that three arrays have identical dimensions
template <typename T>
void Vtr::WriteVector3d(const char *name, T ***vx, T ***vy, T ***vz, unsigned int xSize, unsigned int ySize, unsigned int zSize)
{
    if ((xSize != fd.nxx) || (ySize != fd.nyy) || (zSize != fd.nzz)) 
		HandleWarning("Vtr::WriteVector3d", "Incompatible component and mesh sizes.");
//    
	char Format[256];
	PrepareFormat(**vx, Format, 3);
    fd.Fprintf(dataArrayFormat, name, 3);
    for(unsigned int k=0; k<fd.nzz; ++k)
        for(unsigned int j=0; j<fd.nyy; ++j)
            for(unsigned int i=0; i<fd.nxx; ++i)
				fd.Fprintf(Format, vx[i][j][k], vy[i][j][k], vz[i][j][k]);
	fd.Fprintf(endDataArray);
//
	WriteDataVector(vx, "X", name);
	WriteDataVector(vy, "Y", name);
	WriteDataVector(vz, "Z", name);
}

template <typename T>
void Vtr::WriteVector3d(const char *name, Array3F<T> &vx, Array3F<T> &vy, Array3F<T> &vz, unsigned int xSize, unsigned int ySize, unsigned int zSize)
{
    if ((xSize != fd.nxx) || (ySize != fd.nyy) || (zSize != fd.nzz)) 
		HandleWarning("Vtr::WriteVector3d", "Incompatible component and mesh sizes.");
//    
	char Format[256];
    unsigned int i, j, k;
	PrepareFormat(vx, Format, 3);
    fd.Fprintf(dataArrayFormat, name, 3);
    for(k=0; k<fd.nzz; ++k)
        for(j=0; j<fd.nyy; ++j)
            for(i=0; i<fd.nxx; ++i)
				fd.Fprintf(Format, vx.Value(i, j, k), vy.Value(i, j, k), vz.Value(i, j, k));
	fd.Fprintf(endDataArray);
//
	sprintf(Format, "X_%s", name);
    fd.Fprintf(dataArrayFormat, Format, 1);
    for(k=0; k<fd.nzz; ++k)
        for(j=0; j<fd.nyy; ++j)
            for(i=0; i<fd.nxx; ++i)
				fd.Fprintf(outFormat, vx.Value(i, j, k));
	fd.Fprintf(endDataArray);
//
	sprintf(Format, "Y_%s", name);
    fd.Fprintf(dataArrayFormat, Format, 1);
    for(k=0; k<fd.nzz; ++k)
        for(j=0; j<fd.nyy; ++j)
            for(i=0; i<fd.nxx; ++i)
				fd.Fprintf(outFormat, vy.Value(i, j, k));
	fd.Fprintf(endDataArray);
//
	sprintf(Format, "Z_%s", name);
    fd.Fprintf(dataArrayFormat, Format, 1);
    for(k=0; k<fd.nzz; ++k)
        for(j=0; j<fd.nyy; ++j)
            for(i=0; i<fd.nxx; ++i)
				fd.Fprintf(outFormat, vz.Value(i, j, k));
	fd.Fprintf(endDataArray);
}

template <typename T>
void Vtr::WriteScalar2d(const char *name, T **field, unsigned int xSize, unsigned int ySize)
{
	if ((xSize != fd.nxx) || (ySize != fd.nyy))
		HandleWarning("Vtr::WriteScalar2d", "Incompatible FIELD and MESH sizes.");

    fd.Fprintf(dataArrayFormat, name, 1);
    for(unsigned int j=0; j<fd.nyy; ++j)
        for(unsigned int i=0; i<fd.nxx; ++i)
			fd.Fprintf(outFormat, field[i][j]);
	fd.Fprintf(endDataArray);
}

template <typename T>
void Vtr::WriteScalar3d(const char *name, T ***field, unsigned int xSize, unsigned int ySize, unsigned int zSize)
{
	if ((xSize != fd.nxx) || (ySize != fd.nyy) || (zSize != fd.nzz))
		HandleWarning("Vtr::WriteScalar3d", "Incompatible FIELD and MESH sizes.");

    fd.Fprintf(dataArrayFormat, name, 1);
    for(unsigned int k=0; k<fd.nzz; ++k)
        for(unsigned int j=0; j<fd.nyy; ++j)
            for(unsigned int i=0; i<fd.nxx; ++i)
				fd.Fprintf(outFormat, field[i][j][k]);
	fd.Fprintf(endDataArray);
}

template <typename T>
void Vtr::WriteScalar3d(const char *name, Array3F<T> &array)
{
	if ((array.GetSize(0) != fd.nxx) || (array.GetSize(1) != fd.nyy) || (array.GetSize(2) != fd.nzz)) 
		HandleWarning("Vtr::WriteVector3d", "Incompatible component and mesh sizes.");

	char Format[256];
	PrepareFormat(array, Format, 1);
    fd.Fprintf(dataArrayFormat, name, 1);
    for(unsigned int k=0; k<fd.nzz; ++k)
        for(unsigned int j=0; j<fd.nyy; ++j)
            for(unsigned int i=0; i<fd.nxx; ++i)
				fd.Fprintf(Format, array.Value(i, j, k));
	fd.Fprintf(endDataArray);
}

#endif // __VTR_H__
