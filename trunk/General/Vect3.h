#ifndef __VECT3_H__
#define __VECT3_H__

#include <stdio.h>
#include <string.h>
#include "General.h"
#include "Definitions.h"

template <typename T>
class Vect3
{
public:
	T data[3];

public:
	Vect3(void) { Clear(); }
	Vect3(T a, T b, T c) { Set(a, b, c); }
	Vect3(T *a) { Set(a[0], a[1], a[2]); }
	~Vect3(void) { }
	inline void Clear(void) { memset(data, 0, 3*sizeof(T)); }
	inline void Set(T a, T b, T c) { data[0] = a; data[1] = b; data[2] = c; }
	void Load(char *Buffer, const char *Format)
	{
		char *ia = Buffer;
		T x;
		for(int i=0; i<3; ++i)
		{
			char *ib = strtok(ia, " ");
			sscanf(ib, Format, &x);
			data[i] = x;
			if (!i)
				ia = NULL;
		}
	}

public:
	inline Vect3 operator+=(const Vect3 &op)
	{
		data[0] += op.data[0];
		data[1] += op.data[1];
		data[2] += op.data[2];
		return *this;
	}
	inline void Ortify(void)
	{
		real m = Mod();
		data[0] /= m;
		data[1] /= m;
		data[2] /= m;
	}
	Vect3 Magnify(const Vect3 &factor)
	{
		return Vect3(data[0] * factor.data[0], data[1] * factor.data[1], data[2] * factor.data[2]);
	}
	inline Vect3 operator*(const Vect3 &op) const
	{
		return Vect3(data[1]*op.data[2] - data[2]*op.data[1], data[2]*op.data[0] - data[0]*op.data[2], data[0]*op.data[1] - data[1]*op.data[0]);
	}
	template <typename U>
	inline void Vector(U *op, U *res) const
	{
		res[0] = op[2] * data[1] - op[1] * data[2];
		res[1] = op[0] * data[2] - op[2] * data[0];
		res[2] = op[1] * data[0] - op[0] * data[1];
	}
	inline Vect3 operator+(const Vect3 &op) const
	{
		return Vect3(data[0] + op.data[0], data[1] + op.data[1], data[2] + op.data[2]);
	}
	inline Vect3 Add(T *op)
	{
		return Vect3(data[0] + op[0], data[1] + op[1], data[2] + op[2]);
	}
	inline Vect3 operator-(const Vect3 &op) const
	{
		return Vect3(data[0] - op.data[0], data[1] - op.data[1], data[2] - op.data[2]);
	}
	inline Vect3 Sub(T *op)
	{
		return Vect3(data[0] - op[0], data[1] - op[1], data[2] - op[2]);
	}
	inline Vect3 operator/(T op)
	{
		return Vect3(data[0]/op, data[1]/op, data[2]/op);
	}
	inline Vect3 operator*(T op)
	{
		return Vect3(data[0]*op, data[1]*op, data[2]*op);
	}
	inline Vect3 operator-=(const Vect3 &op)
	{
		data[0] -= op.data[0];
		data[1] -= op.data[1];
		data[2] -= op.data[2];
		return *this;
	}
	inline Vect3 operator*=(T op)
	{
		data[0] *= op;
		data[1] *= op;
		data[2] *= op;
		return *this;
	}
	inline Vect3 operator/=(T op)
	{
		data[0] /= op;
		data[1] /= op;
		data[2] /= op;
		return *this;
	}
	inline Vect3 operator-(void) const
	{
		return Vect3(-data[0], -data[1], -data[2]);
	}
	inline T Scalar(const Vect3<T> &op) const 
	{
		return data[0]*op.data[0] + data[1]*op.data[1] + data[2]*op.data[2];
	}
	template <typename U>
	inline U Scalar(U *op)
	{
		return op[0]*data[0] + op[1]*data[1] + op[2]*data[2];
	}
	inline real Mod(void) const 
	{ 
		return sqrt(ModSquared());
	}
	inline real ModSquared(void) const
	{
		return data[0]*data[0] + data[1]*data[1] + data[2]*data[2];
	}
	inline void Copy(const Vect3<T> &src)
	{
		memcpy(data, src.data, 3*sizeof(T));
	}
	void Fprintf(FILE *file, const char *Format, const char *Prefix = NULL, const char *Suffix = NULL)
	{
		if (Prefix)
			fprintf(file, Prefix);
		for(int i=0; i<3; ++i)
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
		for(int i=0; i<3; ++i)
		{
			sprintf(bx, Format, data[i]);
			strcat(Buffer, bx);
		}
		if (Suffix)
			strcat(Buffer, Suffix);
	}
	void Write(int file)
	{
		write(file, data, 3*sizeof(T));
	}
	void Read(int file)
	{
		read(file, data, 3*sizeof(T));
	}
};

#endif
