#ifndef __VECT6_H__
#define __VECT6_H__

#include <stdlib.h>

template <typename T>
class Vect6
{
public:
	T data[6];

public:
	Vect6() { Clear(); }
	~Vect6() { }
	inline void Clear() { memset(data, 0, 6*sizeof(T)); }
	inline void Set(T a, T b, T c, T d, T e, T f) { data[0] = a; data[1] = b; data[2] = c; data[3] = d; data[4] = e; data[5] = f;}
	void Load(char *Buffer, const char *Format)
	{
        char *ia = Buffer;
		T x;
		for(int i=0; i<6; ++i)
		{
			char *ib = strtok(ia, " ");
            x = (real)atof(ib);
			data[i] = x;
			if (!i)
				ia = NULL;
		}
	}
	void Fprintf(FILE *file, const char *Format, const char *Prefix = NULL, const char *Suffix = NULL)
	{
		if (Prefix)
			fprintf(file, Prefix);
		for(int i=0; i<6; ++i)
			fprintf(file, Format, data[i]);
		if (Suffix)
			fprintf(file, Suffix);
	}
	void Sprintf(char *Buffer, const char *Format, const char *Prefix = NULL, const char *Suffix = NULL)
	{
		char bx[256];
		if (Prefix)
			strcpy(Buffer, Prefix);
		else
			Buffer[0] = '\0';
		for(int i=0; i<6; ++i)
		{
			sprintf(bx, Format, data[i]);
			strcat(Buffer, bx);
		}
		if (Suffix)
			strcat(Buffer, Suffix);
	}
	void Write(int file)
	{
		write(file, data, 6*sizeof(T));
	}
	void Read(int file)
	{
		read(file, data, 6*sizeof(T));
	}
};

#endif
