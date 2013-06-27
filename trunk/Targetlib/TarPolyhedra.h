#ifndef __TARPOLYHEDRA_H__
#define __TARPOLYHEDRA_H__

#include "AbstractTarget.h"
#include "Targetlib.h"

class TARGETLIB_API TarPolyhedra : public AbstractTarget
{
protected:
	int numNormals, numVertices, nlong;
	Vect3<real> *normals, *vertices;
	TarPolyhedra(void) { numNormals = 0; numVertices = 0; normals = NULL; vertices = NULL; }

public:
	TarPolyhedra(TargetManager *man) : AbstractTarget(man) { }
	virtual ~TarPolyhedra(void) 
	{
		if (numNormals)
			CleanDelete2(normals);
		if (numVertices)
			CleanDelete2(vertices);
	}

public:
	virtual void Sizer(void) { }
	virtual void Descriptor(void);
	virtual void Vector(void);
	virtual void VectorX(void);
	virtual void Allocator(void);
	virtual void OutShpar(FILE *file);

protected:
	virtual bool Check(real x, real y, real z) = 0;
	virtual void GenerateVertices(void) = 0;
	virtual void GenerateNormals(void) = 0;
};

class TARGETLIB_API Target_Octahedron : public TarPolyhedra
{
public:
	Target_Octahedron(void) { numNormals = 8; numVertices = 6; }

public:
	Target_Octahedron(TargetManager *man);
	virtual ~Target_Octahedron(void) { }
	virtual void SayHello(FILE *stream);
	virtual void Sizer(void);
	virtual void Descriptor(void);
	virtual void Allocator(void);

protected:
	virtual bool Check(real x, real y, real z);
	virtual void GenerateVertices(void);
	virtual void GenerateNormals(void);
};

class TARGETLIB_API Target_Icosahedron : public TarPolyhedra
{
public:
	Target_Icosahedron(void) { numNormals = 20; numVertices = 12; }

public:
	Target_Icosahedron(TargetManager *man);
	virtual ~Target_Icosahedron(void) { }
	virtual void SayHello(FILE *stream);
	virtual void Sizer(void);
	virtual void Descriptor(void);

protected:
	virtual bool Check(real x, real y, real z);
	virtual void GenerateVertices(void);
	virtual void GenerateNormals(void);
};

class TARGETLIB_API Target_Dodecahedron : public TarPolyhedra
{
public:
	Target_Dodecahedron(void) { numNormals = 12; numVertices = 20; }

public:
	Target_Dodecahedron(TargetManager *man);
	virtual ~Target_Dodecahedron(void) { }
	virtual void SayHello(FILE *stream);
	virtual void Sizer(void);
	virtual void Descriptor(void);

protected:
	virtual bool Check(real x, real y, real z);
	virtual void GenerateVertices(void);
	virtual void GenerateNormals(void);
};

#endif // __TARPOLYHEDRA_H__