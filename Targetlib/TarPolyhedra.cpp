#include "StdAfx.h"

#include "TarPolyhedra.h"
#include "TargetManager.h"

void TarPolyhedra::Descriptor(void)
{
	sprintf(freeDescr, " Unspecified platonic shape, definitevly the error");
}

void TarPolyhedra::Vector(void)
{
	VectorA();
}

void TarPolyhedra::VectorX(void)
{
//
// Set X0 = location in TF of IXYZ=0 0 0
// We assume that origin of TF is located at centroid of hex prism
// Set composition
	x0.Clear();
	for(int i=0; i<nat0; ++i)
	{
		x0 += Vect3<real>((real)ixyz.Value(i, 0), (real)ixyz.Value(i, 1), (real)ixyz.Value(i, 2));
	}
	x0 /= -(real)nat0;
}

void TarPolyhedra::Allocator(void)
{
//
// Current version of TARPLATONIC is restricted to cubic lattices
	if ((dx.data[0] != onex_) || (dx.data[1] != onex_))
		Errmsg("Fatal", shortDescr.c_str(), " TargetIcosahedron does not support noncubic lattice");
//
	int curSize = nx * ny;
	ixyz.Dimension(curSize, 3);
//
	nat0 = 0;
	for(int jx=minJx; jx<=maxJx; ++jx)
	{
		real x = (real)jx;
		for(int jy=minJy; jy<=maxJy; ++jy)
		{
			real y = (real)jy;
			for(int jz=minJz; jz<=maxJz; ++jz)
			{
				real z = (real)jz;
				if (Check(x, y, z))
				{
					if (nat0 >= curSize)
					{
						ixyz.Extend(nz);
						curSize = ixyz.GetSize(0);
					}
					ixyz.Fill3(nat0, jx, jy, jz);
					int index = GetLinearAddress(nat0);
					Composer(index, 0);
					iocc[index] = true;
					++nat0;
				}
			}
		}
	}
	ixyz.Close(nat0);
}

void TarPolyhedra::OutShpar(FILE *file)
{
	fprintf(file, " Nlong = %8.4lf\n", shpar[0]); 
}

// => Target_Octahedron

const char *TargetVerboseDescriptor_Octahedron(int num)
{
	static const char *descr[1] = 
	{
		"distance between opposite vertices"
	};
	if (num < 1)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(Octahedron,1,false,-1,0,"Normal octahedron")
void Target_Octahedron::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "distance between opposite vertices:\n");
	fprintf(stream, "%20.16lf\n", manager->GetShpar(0));
}

Target_Octahedron::Target_Octahedron(TargetManager *man) : TarPolyhedra(man) 
{ 
	numNormals = 8;
	numVertices = 6;
	shortDescr = string("Octahedron");
	longDescr = string("Norman octahedron");
	normals = new Vect3<real>[numNormals];
	vertices = new Vect3<real>[numVertices];
}

void Target_Octahedron::Descriptor(void)
{
	sprintf(freeDescr, " Single octahedron particle NAT=%7d dipoles", nat0);
}

void Target_Octahedron::Sizer(void)
{
	dx = manager->CashedDx();

	minJx = minJy = minJz =  ihuge_;
	maxJx = maxJy = maxJz = -ihuge_;

	GenerateVertices();
	GenerateNormals();

	int nlong = (int)shpar[0];
	real nlongHalf = nlong / (real)2.;
//
	for(int jx=0; jx<=nlong; ++jx)
	{
		real x = (real)jx - nlongHalf;
		for(int jy=0; jy<=nlong; ++jy)
		{
			real y = (real)jy - nlongHalf;
			for(int jz=0; jz<=nlong; ++jz)
			{
				real z = (real)jz - nlongHalf;
				if (Check(x, y, z))
				{
					InternalMinMax(jx, jy, jz);
				}
			}
		}
	}

	AllocateArrays(maxJx - minJx + 1, maxJy - minJy + 1, maxJz - minJz + 1);
}

void Target_Octahedron::Allocator(void)
{
//
// Current version of TARPLATONIC is restricted to cubic lattices
	if ((dx.data[0] != onex_) || (dx.data[1] != onex_))
		Errmsg("Fatal", shortDescr.c_str(), " TargetOctahedron does not support noncubic lattice");
//
	int nlong = (int)shpar[0];
	real nlongHalf = nlong / (real)2.;
//
	int curSize = nx * ny;
	ixyz.Dimension(curSize, 3);
//
	nat0 = 0;
	for(int jx=minJx; jx<=maxJx; ++jx)
	{
		real x = (real)jx - nlongHalf;
		for(int jy=minJy; jy<=maxJy; ++jy)
		{
			real y = (real)jy - nlongHalf;
			for(int jz=minJz; jz<=maxJz; ++jz)
			{
				real z = (real)jz - nlongHalf;
				if (Check(x, y, z))
				{
					if (nat0 >= curSize)
					{
						ixyz.Extend(nz);
						curSize = ixyz.GetSize(0);
					}
					ixyz.Fill3(nat0, jx, jy, jz);
					int index = GetLinearAddress(nat0);
					Composer(index, 0);
					iocc[index] = true;
					++nat0;
				}
			}
		}
	}
	ixyz.Close(nat0);
}

bool Target_Octahedron::Check(real x, real y, real z)
{
	const real Sn = ((int)shpar[0] + 1) / (real)6.;		// the point at the center of face
	real tmp[3];
	bool res = true;
	for(int i=0; i<numNormals; ++i)
	{
		tmp[0] = (i & 1) ? (x + Sn) : (x - Sn);
		tmp[1] = (i & 2) ? (y + Sn) : (y - Sn);
		tmp[2] = (i & 4) ? (z + Sn) : (z - Sn);
		res = res && (normals[i].Scalar(tmp) <= zero_);
	}
	return res;
}

void Target_Octahedron::GenerateVertices(void)
{
	const real sqrt2 = (real)1. / Sqrt((real)2.);
	const real al = sqrt2 * shpar[0];
 	vertices[0].Set(zero_, zero_,    al);
	vertices[1].Set(   al, zero_, zero_);
	vertices[2].Set(zero_,    al, zero_);
	vertices[3].Set(  -al, zero_, zero_);
	vertices[4].Set(zero_,   -al, zero_);
	vertices[5].Set(zero_, zero_,   -al);
}

void Target_Octahedron::GenerateNormals(void)
{
	const int ind[8][4] =			// [numNormals][4]
	{
		{ 1, 0, 2, 0 }, { 2, 0, 3, 1 }, { 3, 0, 4, 2 }, { 4, 0, 1, 3 }, 
		{ 2, 5, 1, 4 }, { 3, 5, 2, 5 }, { 4, 5, 3, 6 }, { 1, 5, 4, 7 }
	};
	Vect3<real> vdirxyz1, vdirxyz2;
	for(int i=0; i<numNormals; ++i)
	{
		vdirxyz1 = vertices[ind[i][0]] - vertices[ind[i][1]];
		vdirxyz1.Ortify();
		vdirxyz2 = vertices[ind[i][2]] - vertices[ind[i][1]];
		vdirxyz2.Ortify();
		normals[ind[i][3]] = vdirxyz1 * vdirxyz2;
		normals[ind[i][3]].Ortify();
	}
}

// => Target_Icosahedron

const char *TargetVerboseDescriptor_Icosahedron(int num)
{
	static const char *descr[1] = 
	{
		"length of the edge"
	};
	if (num < 1)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(Icosahedron,1,false,-1,0,"Normal icosahedron")
void Target_Icosahedron::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "length of the edge:\n");
	fprintf(stream, "%20.16lf\n", manager->GetShpar(0));
}

Target_Icosahedron::Target_Icosahedron(TargetManager *man) : TarPolyhedra(man) 
{
	numNormals = 20;
	numVertices = 12;
	shortDescr = string("Icosahedron");
	longDescr = string("Norman icosahedron");
	normals = new Vect3<real>[numNormals];
	vertices = new Vect3<real>[numVertices];
}

void Target_Icosahedron::Descriptor(void)
{
	sprintf(freeDescr, " Single icosahedron particle NAT=%7d dipoles", nat0);
}

void Target_Icosahedron::Sizer(void)
{
	dx = manager->CashedDx();

	minJx = minJy = minJz =  ihuge_;
	maxJx = maxJy = maxJz = -ihuge_;

	GenerateVertices();
	GenerateNormals();

	int na  = (int)vertices[ 1].data[0];
	int nb  = (int)vertices[10].data[1];
	int nzz = (int)vertices[ 0].data[2];

	for(int ix=-na; ix<=na; ++ix)
	{
		real x = (real)ix;
		for(int iy=-nb; iy<=nb; ++iy)
		{
			real y = (real)iy;
			for(int iz=-nzz; iz<=nzz; ++iz)
			{
				real z = (real)iz;
				bool bOk = Check(x, y, z);
				if (bOk)
				{
					InternalMinMax(ix, iy, iz);
				}
			}
		}
	}
//
	AllocateArrays(maxJx - minJx + 1, maxJy - minJy + 1, maxJz - minJz + 1);
}

bool Target_Icosahedron::Check(real x, real y, real z)
{
	const int indd[20] = { 0, 0, 0, 0, 0, 5, 6, 7, 8, 8, 9, 9, 9, 10, 10, 11, 11, 11, 11, 11 };

	Vect3<real> tmp;
	bool res = true;
	for(int i=0; i<numNormals; ++i)
	{
		tmp.Set(x, y, z); 
		res = res && ((normals[i].Scalar(tmp - vertices[indd[i]])) <= (real)0.);
	}
	return res;
}

void Target_Icosahedron::GenerateVertices(void)
{
	const real phia = (real)(Pi * 26.56505 / 180.); 
	const real ang72 = (real)(Pi * 72. / 180.);
	const real ang36 = (real)(Pi * 36. / 180.);

	real factor = shpar[0] / (real)1.0514622;			// radius of enscribing sphere
	vertices[0].Set(zero_, zero_, factor);
	for(int i=0; i<5; ++i)
	{
		real ang1 = ang72 * i;
		real ang2 = ang1 - ang36;
		vertices[i+1].Set(factor * Cos(-ang1) * Cos( phia), factor * Sin(-ang1) * Cos( phia), factor * Sin( phia));
		vertices[i+6].Set(factor * Cos(-ang2) * Cos(-phia), factor * Sin(-ang2) * Cos(-phia), factor * Sin(-phia));
	}
	vertices[11].Set(zero_, zero_, -factor);
}

void Target_Icosahedron::GenerateNormals(void)
{
	const int ind[20][4] =			// [numNormals][4]
	{
			{ 2,  0, 1,  0 }, { 3,  0, 2,  1 }, { 4,  0, 3,  2 }, { 5,  0,  4,  3 }, {  1,  0,  5,  4 }, 
			{ 1,  5, 6,  5 }, { 1,  7, 2,  7 }, { 2,  8, 3,  9 }, { 3,  9,  4, 11 }, {  4, 10,  5, 13 }, 
			{ 7,  1, 6,  6 }, { 8,  2, 7,  8 }, { 9,  3, 8, 10 }, {10,  4,  9, 12 }, {  6,  5, 10, 14 }, 
			{ 6, 11, 7, 15 }, { 7, 11, 8, 16 }, { 8, 11, 9, 17 }, { 9, 11, 10, 18 }, { 10, 11,  6, 19 }
	};
	Vect3<real> vdirxyz1, vdirxyz2;
	for(int i=0; i<numNormals; ++i)
	{
		vdirxyz1 = vertices[ind[i][0]] - vertices[ind[i][1]];
		vdirxyz1.Ortify();
		vdirxyz2 = vertices[ind[i][2]] - vertices[ind[i][1]];
		vdirxyz2.Ortify();
		normals[ind[i][3]] = vdirxyz1 * vdirxyz2;
		normals[ind[i][3]].Ortify();
	}
}

// => Target_Dodecahedron

const char *TargetVerboseDescriptor_Dodecahedron(int num)
{
	static const char *descr[1] = 
	{
		"length of the edge"
	};
	if (num < 1)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(Dodecahedron,1,false,-1,0,"Normal dodecahedron")
void Target_Dodecahedron::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "length of the edge:\n");
	fprintf(stream, "%20.16lf\n", manager->GetShpar(0));
}

Target_Dodecahedron::Target_Dodecahedron(TargetManager *man) : TarPolyhedra(man) 
{
	numNormals = 12;
	numVertices = 20;
	shortDescr = string("Icosahedron");
	longDescr = string("Norman icosahedron");
	normals = new Vect3<real>[numNormals];
	vertices = new Vect3<real>[numVertices];
}

void Target_Dodecahedron::Descriptor(void)
{
	sprintf(freeDescr, " Single dodecahedron particle NAT=%7d dipoles", nat0);
}

void Target_Dodecahedron::GenerateVertices(void)
{
	const real phia = (real)(Pi * 52.62263590 / 180.); 
	const real phib = (real)(Pi * 10.81231754 / 180.);
	const real ang72 = (real)(Pi * 72. / 180.);
	const real ang36 = (real)(Pi * 36. / 180.);

	real factor = shpar[0] / (real)0.713644;			// radius of enscribing sphere

	for(int i=0; i<5; ++i)
	{
		real ang1 = ang72 * i;
		real ang2 = ang1 + ang36;
		vertices[i   ].Set(factor * Cos(ang1) * Cos( phia), factor * Sin(ang1) * Cos( phia), factor * Sin( phia));
		vertices[i +5].Set(factor * Cos(ang1) * Cos( phib), factor * Sin(ang1) * Cos( phib), factor * Sin( phib));
		vertices[i+10].Set(factor * Cos(ang2) * Cos(-phib), factor * Sin(ang2) * Cos(-phib), factor * Sin(-phib));
		vertices[i+15].Set(factor * Cos(ang2) * Cos(-phia), factor * Sin(ang2) * Cos(-phia), factor * Sin(-phia));
	}
}

void Target_Dodecahedron::GenerateNormals(void)
{
	const int ind[12][4] =			// [numNormals][4]
	{
			{  1,  0,  2, 0 }, {  6,  0, 1,  1 }, {  7,  1,  2,  2 }, {  8,  2,  3,  3 }, 
			{  9,  3,  4, 4 }, {  4,  0, 5,  5 }, { 16, 15, 11,  6 }, { 17, 16, 12,  7 }, 
			{ 18, 17, 13, 8 }, { 19, 18, 14, 9 }, { 15, 19, 10, 10 }, { 17, 15, 16, 11 }
	};
	Vect3<real> vdirxyz1, vdirxyz2;
	for(int i=0; i<numNormals; ++i)
	{
		vdirxyz1 = vertices[ind[i][0]] - vertices[ind[i][1]];
		vdirxyz1.Ortify();
		vdirxyz2 = vertices[ind[i][2]] - vertices[ind[i][1]];
		vdirxyz2.Ortify();
		normals[ind[i][3]] = vdirxyz1 * vdirxyz2;
		normals[ind[i][3]].Ortify();
	}
}

void Target_Dodecahedron::Sizer(void)
{
	dx = manager->CashedDx();

	minJx = minJy = minJz =  ihuge_;
	maxJx = maxJy = maxJz = -ihuge_;

	GenerateVertices();
	GenerateNormals();

	int na  = (int)vertices[5].data[0];
	int nb  = (int)vertices[6].data[1];
	int nzz = (int)vertices[0].data[2];

	for(int ix=-na; ix<=na; ++ix)
	{
		real x = (real)ix;
		for(int iy=-nb; iy<=nb; ++iy)
		{
			real y = (real)iy;
			for(int iz=-nzz; iz<=nzz; ++iz)
			{
				real z = (real)iz;
				bool bOk = Check(x, y, z);
				if (bOk)
				{
					InternalMinMax(ix, iy, iz);
				}
			}
		}
	}
//
	AllocateArrays(maxJx - minJx + 1, maxJy - minJy + 1, maxJz - minJz + 1);
}

bool Target_Dodecahedron::Check(real x, real y, real z)
{
	const int indd[12] = { 0, 0, 1, 2, 3, 4, 15, 16, 17, 18, 19, 15 };

	Vect3<real> tmp;
	bool res = true;
	for(int i=0; i<numNormals; ++i)
	{
		tmp.Set(x, y, z);
		res = res && ((normals[i].Scalar(tmp - vertices[indd[i]])) <= (real)0.);
	}
	return res;
}