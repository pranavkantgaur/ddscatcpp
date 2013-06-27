#include "StdAfx.h"

#include "Vtr.h"

#ifdef _WIN32
	const char slash = '\\';
#else
	const char slash = '/';
#endif

Vtr::VtrFileHandle::VtrFileHandle(void)
{
	prefix = NULL;
	counter = restart = 0;
	Clear();
}

Vtr::VtrFileHandle::~VtrFileHandle(void)
{
	CleanDelete2(prefix);
	if (tabs)
		free(tabs);
	tabs = NULL;
}

void Vtr::VtrFileHandle::Clear(void)
{
	handle = NULL;
	unit = nxx = nyy = nzz = 0;
	first = true;
	isOpened = false;
	currentTab = 0;
	tabs = (char *)malloc(32 * sizeof(char));
	memset(tabs, ' ', 32 * sizeof(char));
}

void Vtr::VtrFileHandle::SetPrefix(const char *data)
{
	size_t len = prefix ? strlen(prefix) : 0;
	size_t lenData = strlen(data);
	if ((len != lenData) && prefix)
		delete [] prefix;
	prefix = new char [lenData + 1];
	strcpy(prefix, data);
}

bool Vtr::VtrFileHandle::OpenHandle(const char *f, const char *mode)
{
	handle = fopen(f, mode);
	return (handle != NULL) ? true : false;
}

void Vtr::VtrFileHandle::CloseHandle()
{
	fclose(handle);
	handle = NULL;
}

void Vtr::VtrFileHandle::Tabify(void)
{
	fwrite(tabs, sizeof(char), currentTab, handle); 
}

void Vtr::VtrFileHandle::Tabify(bool dir)
{
	if (dir)
	{
		Tabify();
		currentTab += 3;
		if (currentTab > sizeof(tabs))
		{
			tabs = (char *)realloc(tabs, currentTab * sizeof(char));
			memset(tabs, ' ', currentTab * sizeof(char));
		}
	}
	else
	{
		currentTab -= 3;
		Tabify();
	}
}

void Vtr::VtrFileHandle::Fprintf(const char *Format, ...)
{
	bool bDelayed = false;
	if (Format[0] == '<')
	{
		switch(Format[1])
		{
		case '/':
			Tabify(false);
			break;

		case '?':
			Tabify();
			break;

		default:
			Tabify(true);
			size_t len = strlen(Format) - 1;
			while((Format[len] == 0x0a) || (Format[len] == 0x0d))
				--len;
			if ((Format[len] == '>') && (Format[len-1] == '/'))
				bDelayed = true;
			break;
		}
	}
	else
		Tabify();

	va_list list;
	va_start(list, Format);
	vfprintf(handle, Format, list);
	va_end(list);

	if (Format[0] != '<')
		fprintf(handle, "\n");
	else
	{
		if (bDelayed)
			currentTab -= 3;
	}
}

const char *Vtr::outFormat = realFormat;
const char *Vtr::xmlHeader = "<?xml version=\"1.0\"?>\n";
const char *Vtr::endDataArray = "</DataArray>\n";
const char *Vtr::dataArrayFormat = "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"ascii\">\n";

Vtr::Vtr(void)
{
	iproc = 0;
	nb_procs = 1;
	outFormat = realFormat;
}

Vtr::~Vtr(void)
{

}

void Vtr::OpenFile(const char *prefix, int proc_rank, int num_procs, int restart)
{
	fd.SetPrefix(prefix);
	if ((proc_rank != -1) && (num_procs != -1))
	{
		iproc = proc_rank;
		nb_procs = num_procs;
	}
//
	if (fd.first && (restart != -1))
	{
		fd.restart = restart;
		fd.counter = restart;
		fd.first = false;
	}
	++fd.counter;
//
	char f[256];
	if (proc_rank != -1)
		sprintf(f, "%s_%d_%d.vtr", fd.GetPrefix(), iproc, fd.counter);
	else
		sprintf(f, "%s_%d.vtr", fd.GetPrefix(), fd.counter);
	bool bOp = fd.OpenHandle(f, "w+");
	if (!bOp) 
		fprintf(stderr, "Problem creating file %s.\n", f);

	fd.Fprintf(xmlHeader);
	fd.Fprintf("<VTKFile type=\"RectilinearGrid\" version=\"0.1\" format=\"ascii\">\n");
	fd.isOpened = true;
}

void Vtr::CloseFile()
{
	if (fd.isOpened)
	{
		fd.Fprintf("</PointData>\n");
		fd.Fprintf("</Piece>\n");
		fd.Fprintf("</RectilinearGrid>\n");
		fd.Fprintf("</VTKFile>\n");
		fd.CloseHandle();
		fd.Clear();
	}
	else
		HandleWarning("VTR_close_file", "No such file to close. Please, check file descriptor.");
}

void Vtr::WriteDataArray(real *x, int nx, const char *Name, int size)
{
	fd.Fprintf(dataArrayFormat, Name, size);
	for(int i=0; i<nx; ++i)
		fd.Fprintf(outFormat, x[i]);
	fd.Fprintf(endDataArray);
}

void Vtr::WriteDataVector(real **x, const char *Prefix, const char *Name)
{
	char Buffer[256];
	sprintf(Buffer, "%s_%s", Prefix, Name);
	fd.Fprintf(dataArrayFormat, Buffer, 1);
    for(unsigned int j=0; j<fd.nyy; ++j)
        for(unsigned int i=0; i<fd.nxx; ++i)
			fd.Fprintf(outFormat, x[i][j]);
	fd.Fprintf(endDataArray);
}

void Vtr::WriteDataVector(real ***x, const char *Prefix, const char *Name)
{
	char Buffer[256];
	sprintf(Buffer, "%s_%s", Prefix, Name);
    fd.Fprintf(dataArrayFormat, Buffer, 1);
    for(unsigned int k=0; k<fd.nzz; ++k)
        for(unsigned int j=0; j<fd.nyy; ++j)
            for(unsigned int i=0; i<fd.nxx; ++i)
				fd.Fprintf(outFormat, x[i][j][k]);
	fd.Fprintf(endDataArray);
}

void Vtr::WriteDataArray(int *x, int nx, const char *Name, int size)
{
	fd.Fprintf(dataArrayFormat, Name, size);
	for(int i=0; i<nx; ++i)
		fd.Fprintf("%d", x[i]);
	fd.Fprintf(endDataArray);
}

void Vtr::WriteDataVector(int **x, const char *Prefix, const char *Name)
{
	char Buffer[256];
	sprintf(Buffer, "%s_%s", Prefix, Name);
    fd.Fprintf(dataArrayFormat, Buffer, 1);
    for(unsigned int j=0; j<fd.nyy; ++j)
        for(unsigned int i=0; i<fd.nxx; ++i)
			    fd.Fprintf("%d", x[i][j]);
	fd.Fprintf(endDataArray);
}

void Vtr::WriteDataVector(int ***x, const char *Prefix, const char *Name)
{
	char Buffer[256];
	sprintf(Buffer, "%s_%s", Prefix, Name);
    fd.Fprintf(dataArrayFormat, Buffer, 1);
    for(unsigned int k=0; k<fd.nzz; ++k)
        for(unsigned int j=0; j<fd.nyy; ++j)
            for(unsigned int i=0; i<fd.nxx; ++i)
				fd.Fprintf("%d", x[i][j][k]);
	fd.Fprintf(endDataArray);
}

void Vtr::PrepareFormat(Array3F<real> &a, char *Format, int num)
{
	strcpy(Format, outFormat);
	for(int i=1; i<num; ++i)
	{
		strcat(Format, " ");
		strcat(Format, outFormat);
	}
}

void Vtr::PrepareFormat(Array3F<int> &a, char *Format, int num)
{
	strcpy(Format, "%d");
	for(int i=1; i<num; ++i)
	{
		strcat(Format, " %d");
	}
}

void Vtr::PrepareFormat(real a, char *Format, int num)
{
	strcpy(Format, outFormat);
	for(int i=1; i<num; ++i)
	{
		strcat(Format, " ");
		strcat(Format, outFormat);
	}
}

void Vtr::PrepareFormat(int a, char *Format, int num)
{
	strcpy(Format, "%d");
	for(int i=1; i<num; ++i)
	{
		strcat(Format, " %d");
	}
}

void Vtr::CollectFile()
{
	if (iproc == 0)
	{
		char filename[256];
		sprintf(filename, "%s.pvd", fd.GetPrefix());
		fd.SetHandle(fopen(filename, "w"));
		if (fd.GetHandle() == NULL)
			HandleError("VTR_collect_file: Error, problem creating file\n", " ");
		fd.Fprintf(xmlHeader);
		fd.Fprintf("<VTKFile type=\"Collection\" version=\"0.1\" format=\"ascii\">\n");
	    fd.Fprintf("<Collection>\n");

		char *vtrfile = strrchr(fd.GetPrefix(), slash);
		if (!vtrfile && strlen(fd.GetPrefix()))
			vtrfile = fd.GetPrefix();
		if (nb_procs == 1)
		{
			for(int shot=0; shot<fd.counter; ++shot)
				fd.Fprintf("<DataSet timestep=\"%d\" part=\"0\" file=\"%s_%d.vtr\"/>\n", shot, vtrfile, shot+1);
		}
		else
		{
			for(int k=0; k<nb_procs; ++k)
				for(int shot=0; shot<fd.counter; ++shot)
					fd.Fprintf("<DataSet timestep=\"%d\" part=\"%d\" file=\"%s_%d_%d.vtr\"/>\n", shot, k, vtrfile, k, shot+1);
		}
		fd.Fprintf("</Collection>\n");
		fd.Fprintf("</VTKFile>\n");
		fd.CloseHandle();
	}
	fd.counter = 0;
	fd.restart = 0;
	fd.first = true;
	iproc = 0;
	nb_procs = 1;
}

//void Vtr::WriteVector3d(char *name, Array3F<Complex> &array)
//{
// TODO:
//}
