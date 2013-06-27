#include "StdAfx.h"

#include <algorithm>

#include "LoadableTarget.h"
#include "TargetManager.h"
#include "Line.h"

class DeleteTableElement
{
public:
	template <typename T>
	void operator() (const T *ptr) const
	{
		delete ptr;
	}
};

LoadableTarget::~LoadableTarget(void) 
{ 
    CleanDelete2(cflshp);
	CleanDelete(dfdata);
	if (elements)
		for_each(elements->begin(), elements->end(), DeleteTableElement());
	CleanDelete(elements);
}

void LoadableTarget::OutVectors(FILE *file)
{
	a1.Fprintf(file, "%10.6lf", NULL, " = A_1 vector\n");
	a2.Fprintf(file, "%10.6lf", NULL, " = A_2 vector\n");
	dx.Fprintf(file, "%10.6lf", NULL, " = lattice spacings (d_x,d_y,d_z)/d\n");		
	x0.Fprintf(file, "%10.6lf", NULL, " = lattice offset x0(1-3) = (x_TF,y_TF,z_TF)/d ");		 
}

void LoadableTarget::Printer(void)
{
	if (manager->Ioshp())
	{
		FILE *ioshpFile = fopen("target.out", "w");
		if (ioshpFile)
		{
			fprintf(ioshpFile, " >%s %s;", shortDescr.c_str(), longDescr.c_str()); 
			fprintf(ioshpFile, "%10d = NAT", nat0);
			if (shpar)
				fprintf(ioshpFile, ", Diamx = %8.4lf\n", shpar[0]);
			else
				fprintf(ioshpFile, "\n");
			fprintf(ioshpFile, "%8.4lf%8.4lf%8.4lf = alpha_1-3\n", alpha.data[0], alpha.data[1], alpha.data[2]);
			OutVectors(ioshpFile);
			fprintf(ioshpFile, "for dipole 0 0 0\n     JA  IX  IY  IZ ICOMP(x,y,z)\n");
			for(int ja=0; ja<nat0; ++ja)
			{
				fprintf(ioshpFile, "%7d", ja);
				ixyz.PrintLine(ioshpFile, "%4d", ja);
				int index = GetLinearAddress(ja);
				icomp.PrintLine(ioshpFile, "%2d", index);
				if (dfdata && dfdata->IsAllocated())
					fprintf(ioshpFile, "%10.6lf%10.6lf%10.6lf", dfdata->GetBetadf(index), dfdata->GetPhidf(index), dfdata->GetThetadf(index));
				fprintf(ioshpFile, "\n");
			}
			fclose(ioshpFile);
		}
		else
			fprintf(stderr, "Cannot open ioshpFile in AbstractFile::Printer for %s.", shortDescr.c_str());
	}
}

void LoadableTarget::Reader(void)
{
	fprintf(stderr, ">FROM_FILE open file=%s\n", cflshp);

	FILE *idvshpFile = fopen(cflshp, "r");
	if (idvshpFile)
	{
		char Buffer[256], wrk[256];
		int i, ja, ix, iy, iz, icc1, icc2, icc3;
		fgets(Buffer, 255, idvshpFile);
//
		char *ia = strchr(Buffer, '>');
		if (!ia) 
			ia = Buffer;
		ia++;
		char *ib = strchr(ia, ' ');
		if (!ib)
			ib = Buffer + strlen(Buffer) - 1;
		memcpy(wrk, ia, ib-ia);
		wrk[ib-ia] = '\0';
		shortDescr = string(wrk);
//
		ia = ib+1;
		while(*ia == ' ') ++ia;
		ib = strchr(Buffer, ';');
		if (!ib)
			ib = Buffer + strlen(Buffer) - 1;
		memcpy(wrk, ia, ib-ia);
		wrk[ib-ia] = '\0';
        longDescr = string(wrk);
//
		int nnat;
		fgets(Buffer, 255, idvshpFile);
		sscanf(Buffer, "%d", &nnat);
		fgets(Buffer, 255, idvshpFile);
		a1.Load(Buffer, realFormat);
		fgets(Buffer, 255, idvshpFile);
		a2.Load(Buffer, realFormat);
		fgets(Buffer, 255, idvshpFile);
		dx.Load(Buffer, realFormat);
		fgets(Buffer, 255, idvshpFile);
		x0.Load(Buffer, realFormat);
//
		fgets(Buffer, 255, idvshpFile);		// skip
		minJx = ihuge_;		maxJx = -ihuge_;
		minJy = ihuge_;		maxJy = -ihuge_;
		minJz = ihuge_;		maxJz = -ihuge_;
		for(i=0; i<nnat; ++i)
		{
			fgets(Buffer, 255, idvshpFile);
			sscanf(Buffer, "%d%d%d%d%d%d%d", &ja, &ix, &iy, &iz, &icc1, &icc2, &icc3);
			InternalMinMax(ix, iy, iz);
		}
		fclose(idvshpFile);
		fprintf(stderr, ">FROM_FILE close file=%s\n", cflshp);
	}
	else
	{
		Wrimsg("From File", "Error: Cannot open shape file");
	}
}

bool sortPredicate(const Line *op1, const Line *op2)
{
	return *op1 < *op2;
}

void LoadableTarget::Allocator(void)
{
	FILE *idvshpFile = fopen(cflshp, "r");
	if (idvshpFile)
	{
		int i, ja;
		char Buffer[256];
		for(i=0; i<7; ++i)
		{
			fgets(Buffer, 255, idvshpFile);
			if (i == 1)
			{
				sscanf(Buffer, "%d", &nat0);
				ixyz.Dimension(nat0, 3);
			}
		}
//
		Line::nx = nx;
		Line::ny = ny;
		Line::nz = nz;
		vector<Line *> *temp = new vector<Line *>;
		for(i=0; i<nat0; ++i)
		{
			fgets(Buffer, 255, idvshpFile);
			Line *line = new Line;
			sscanf(Buffer, "%d%d%d%d%d%d%d", &ja, &line->ix, &line->iy, &line->iz, &line->icc1, &line->icc2, &line->icc3);
			temp->push_back(line);
		}
		fclose(idvshpFile);
//
		std::sort(temp->begin(), temp->end(), sortPredicate);
		for(i=0; i<nat0; ++i)
		{
			Line *line = temp->at(i);
			ixyz.Fill3(i, line->ix, line->iy, line->iz);
			int index = GetLinearAddress(i);
			icomp.Fill3(index, line->icc1, line->icc2, line->icc3);
			iocc[index] = true;
		}
		if (temp)
			for_each(temp->begin(), temp->end(), DeleteTableElement());
		CleanDelete(temp);
	}
	else
		Wrimsg("From File::Allocator", "Error: Cannot open shape file");
}

void LoadableTarget::PrepareIanisoSpecial(void)
{
	ianiso = TargetIsAnisotropic;
	if (dfdata && dfdata->IsAllocated())
	{
		ianiso = dfdata->GetIaniso();
	}
}

// Loadable targets below
const char *TargetVerboseDescriptor_Anifilpbc(int num)
{

	return NULL;
}

REGISTER_TARGET(Anifilpbc,2,true,0,0," -- anifilpbc -- ")

void Target_Anifilpbc::PrepareIaniso(void)
{
	PrepareIanisoSpecial();
}

const char *TargetVerboseDescriptor_Anifrmfil(int num)
{

	return NULL;
}

REGISTER_TARGET(Anifrmfil,0,false,-1,0," -- anifrmfil -- ")

void Target_Anifrmfil::PrepareIaniso(void)
{
	PrepareIanisoSpecial();
}
