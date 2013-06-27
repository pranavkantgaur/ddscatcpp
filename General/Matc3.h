#ifndef __MATC3_H__
#define __MATC3_H__

#include "Vect3.h"

template <typename T> 
class Matc3
{
protected:
	T data[3][3];

public:
	Matc3() 
	{
		Clear(); 
	}
	virtual ~Matc3() 
	{ 
	}
	void Clear()
	{
		for(int i=0; i<3; ++i)
		{
			for(int j=0; j<3; ++j)
			{
				data[i][j] = T();
			}
		}
	}
	void Unity()
	{
		for(int i=0; i<3; ++i)
		{
			for(int j=0; j<3; ++j)
			{
				data[i][j] = T();
			}
			data[i][i] = T(1.);
		}
	}
	T &Data(int i, int j) 
	{ 
		return data[i][j]; 
	}
	void Set(real a11, real a12, real a13, real a21, real a22, real a23, real a31, real a32, real a33)
	{
		data[0][0] = a11;  data[0][1] = a12;  data[0][2] = a13; 
		data[1][0] = a21;  data[1][1] = a22;  data[1][2] = a23; 
		data[2][0] = a31;  data[2][1] = a32;  data[2][2] = a33; 
	}

public:
	Matc3<T> operator *(const Matc3<T> &op) const
	{
		Matc3<T> res;

		for(int i=0; i<3; ++i)
		{
			for(int j=0; j<3; ++j)
			{
				T a = T();
				for(int k=0; k<3; ++k)
				{
					a += data[i][k] * op.data[k][j]; 
				}
				res.data[i][j] = a;
			}
		}
		return res;
	}

	Vect3<T> operator *(Vect3<T> &op) const
	{
		Vect3<T> res;

		for(int i=0; i<3; ++i)
		{
			T a = T();
			for(int j=0; j<3; ++j)
			{
				a += data[i][j] * op.data[j];
			}
			res.data[i] = a;
		}
		return res;
	}

	Matc3 operator/=(T op)
	{
		for(int i=0; i<3; ++i)
		{
			for(int j=0; j<3; ++j)
			{
				data[i][j] /= op;
			}
		}
		return *this;
	}

//
// The function builds the rotation matrix (3 x 3) for the given axis (any length) and the angle
// theResultVector = theMatric * the Vector
// Reference: Leubner, C., 1977,  Coordinate-free  rotation operator,Am. J. Phys. 47, 727---729.
// Copyright (C) 1993, B.T. Draine and P.J. Flatau
// Copyright (C) 2013, C++ version, Choliy V.
// This code is covered by the GNU General Public License.
	Matc3 Rot2(const Vect3<real> &a, real theta)
	{
		real ct = Cos(theta);	
		real st = Sin(theta);	
		real anorm = a.Mod();	
		Vect3<real> ahat(a);	
		ahat /= anorm;
//
// cos(\theta) {\bf 1}  term:
		Set(ct, -st*ahat.data[2], st*ahat.data[1], st*ahat.data[2], ct, -st*ahat.data[0], -st*ahat.data[1], st*ahat.data[0], ct);
//
// aa-dyadic (outer-product):	
		for(int i=0; i<3; ++i)	
		{
			for(int j=0; j<3; ++j)
			{
				Data(i,j) += ((real)1. - ct) * ahat.data[j] * ahat.data[i];
			}
		}
	}
};

#endif