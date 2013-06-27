#include "StdAfx.h"

#include "Target_From_file.h"
#include "TargetManager.h"

const char *TargetVerboseDescriptor_From_file(int num)
{

	return NULL;
}

REGISTER_TARGET(From_file,0,false,-1,0," -- from file -- ")
void Target_From_file::Sizer(void)
{
	Reader();
	AllocateArrays(maxJx - minJx + 1, maxJy - minJy + 1, maxJz - minJz + 1);
}

void Target_From_file::Descriptor(void)
{
	sprintf(freeDescr, " From file target, %s: %s", shortDescr.c_str(), longDescr.c_str());
}

void Target_From_file::Vector(void)
{
	VectorA();
}

void Target_From_file::VectorX(void)
{

}

void Target_From_file::PrepareIaniso()
{
	ianiso = TargetIsAnisotropic;
}

const char *TargetVerboseDescriptor_Frmfilpbc(int num)
{
// TODO
	return NULL;
}

REGISTER_TARGET(Frmfilpbc,2,true,0,0," -- frmfilpbc -- ")
void Target_Frmfilpbc::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "periodicity in y and z directions\n");
	fprintf(stream, "%20.16lf %20.16lf\n", shpar[0], shpar[1]);
	fprintf(stream, "Other parameters and shape itself are read from shape.dat file\n");
}

void Target_Frmfilpbc::PreparePyzd(void)
{
	pyd = shpar[0];
	pzd = shpar[1];
}

void Target_Frmfilpbc::Sizer(void)
{
	Reader();
	AllocateArrays(maxJx - minJx + 1, maxJy - minJy + 1, maxJz - minJz + 1);
}

void Target_Frmfilpbc::PrepareIaniso()
{
	ianiso = TargetIsAnisotropic;
}
