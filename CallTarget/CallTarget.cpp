#include "StdAfx.h"

#include "Definitions.h"
#include "Enumerator.h"
#include "Functions.h"
#include "Vect3.h"
#include "TargetManager.h"

char *mgets(char *Buffer, int num, FILE *file);

//
// Program CALLTARGET is used to generate target arrays using target generation routines employed by DDSCAT.
//
// Original version created by Choliy V., Kyiv Shevchenko University 
// based upon fotran version of DDscat 7.2.0 by
// P.J.Flatau, Colorado State Univ. and B.T.Draine, Princeton Univ. Obs.
//
// History:
// All history records of fortran versions removed. 
// 10.04.26 (ChB): Initial version.
// end history
//
// Copyright (C) 2012, V.Choliy
// Copyright (C) All fortran versions of DDscat, 1993-2012, B.T. Draine and P.J. Flatau
// This code is covered by the GNU General Public License.
bool CallTarget(int argc, const char *argv[])
{
	TargetManager *manager = TargetManager::GetInstance();
	manager->CashedDx() = Vect3<real>((real)1., (real)1., (real)1.);					// restrict to cubic lattice:
	manager->Ioshp() = true;

	char cshape[32], *shapeParam = NULL;
	if (argc >= 2)											// file name in command line -> single line with all file in it except target name
	{
		FILE *file = fopen(argv[1], "r");
		if (!file)
		{
			fprintf(stderr, "Cannot open parameter file. Stop.\n");
			TargetManager::Kill();
			return -1;
		}
		char Buffer[80];
		mgets(Buffer, 80, file);
		strcpy(cshape, Buffer);

		long pos1 = ftell(file);
		fseek(file, 0L, SEEK_END);
		long pos2 = ftell(file);
		fseek(file, pos1, SEEK_SET);
		shapeParam = new char[pos2-pos1+16];
		shapeParam[0] = '\0';
		while(!feof(file))
		{
			mgets(Buffer, 80, file);
			strcat(shapeParam, Buffer);
			strcat(shapeParam, " ");
		}
		fclose(file);
	}
	else
	{
		manager->SayHello();
		fprintf(stdout, "What shape? (Enter choice) Choices:\n");
		fflush(stdin);
		mgets(cshape, 32, stdin);
		strupr(cshape);
	}

//	TargetType shapeNum = manager->TargetTypeFromString(cshape);
    manager->PreloadTarget(cshape, shapeParam);
	AbstractTarget *target = manager->LoadTarget(cshape, manager->GetCashedFileName());
	bool bOk = (target != NULL);
	if (bOk)
		manager->GetCurrentTarget()->SayHello();
	else
		printf("Cannot process input/load target. Stop.\n");

	delete [] shapeParam;
	TargetManager::Kill();

	return true;
}

char *mgets(char *Buffer, int num, FILE *file)
{
	char *res = fgets(Buffer, num, file);
	while ((Buffer[strlen(Buffer)-1] == 0x0a) || (Buffer[strlen(Buffer)-1] == 0x0d)) 
		Buffer[strlen(Buffer)-1] = '\0';

	return res;
}
