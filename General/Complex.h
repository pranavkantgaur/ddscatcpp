#ifndef __COMPLEX_H__
#define __COMPLEX_H__

#ifdef _WIN32
	#include <io.h>
#else
	#include <unistd.h>
#endif	
#include "General.h"
#include "Definitions.h"

class Complex
{
public:
	real re, im;

public:
	Complex() : re((real)0.), im((real)0.) { }
	Complex(real r, real i)  : re(r), im(i) { }
	~Complex() { }

public:
	inline void clear() { re = im = (real)0.; }
	inline void unityIm() { re = (real)0.; im = (real)1.; }
	inline void unityRe() { re = (real)1.; im = (real)0.; }
	inline void set(real r, real i) { re = r; im = i; }

public:
	inline Complex conjg()
	{
		return Complex(re, -im);
	}

	inline Complex operator*(const Complex &op) const
	{
		return Complex(re*op.re - im*op.im, re*op.im + im*op.re); 
	}
	inline Complex operator*=(const Complex &op)
	{
		real x = re*op.re - im*op.im;
		real y = re*op.im + im*op.re;

		re = x;
		im = y;

 		return *this; 
	}
	inline Complex operator*(real op) const
	{ 
		return Complex(re*op, im*op); 
	}
	inline Complex operator*=(real op)
	{
		re *= op; 
		im *= op; 
		return *this; 
	}
	inline Complex operator/(const Complex &op) const
	{
		return Complex(re*op.re + im*op.im, op.re*im - re*op.im) / op.modSquared();
	}
	inline Complex operator/=(const Complex &op)
	{
		real x = re*op.re + im*op.im;
		real y = op.re*im - re*op.im;
		real m = op.modSquared();

		re = x / m;
		im = y / m;

		return *this;
	}
	inline Complex operator/(real op) const
	{
		return Complex(re/op, im/op);
	}
	inline Complex operator/=(real op)
	{ 
		re /= op; 
		im /= op; 
		return *this; 
	}
	inline Complex operator+(const Complex &op) const
	{ 
		return Complex(re + op.re, im + op.im); 
	}
	inline Complex operator+=(const Complex &op)
	{ 
		re += op.re; im += op.im; return *this;
	}
	inline Complex operator+(real op) const
	{ 
		return Complex(re + op, im); 
	}
	inline Complex operator+=(real op)
	{ 
		re += op; return *this;
	}
	inline Complex operator+() const
	{
		return *this;
	}
	inline Complex operator-(const Complex &op) const
	{ 
		return Complex(re - op.re, im - op.im); 
	}
	inline Complex operator-=(const Complex &op)
	{ 
		re -= op.re; im -= op.im; return *this; 
	}
	inline Complex operator-(real op) const
	{ 
		return Complex(re - op, im); 
	}
	inline Complex operator-=(real op)
	{
		re -= op;
		return * this;
	}
	inline Complex operator-() const
	{ 
		return Complex(-re, -im); 
	}
	inline bool operator==(const Complex &op) const
	{
		if (re != op.re)
			return false;
		if (im != op.im)
			return false;
		return true;
	}

public:
	inline Complex Add(const Complex &opa, const Complex &opb)
	{
		re = opa.re + opb.re;
		im = opa.im + opb.im;
		return *this;
	}
	inline Complex Add(const Complex &opa, const Complex &opb, real op)
	{
		re = (opa.re + opb.re) * op;
		im = (opa.im + opb.im) * op;
		return *this;
	}
	inline Complex Add2(const Complex &opa, const Complex &opb, real op)
	{
		re = opa.re + opb.re * op;
		im = opa.im + opb.im * op;
		return *this;
	}
	inline Complex Add0(const Complex &opa, real ra, const Complex &opb, real rb)
	{
		re = opa.re * ra + opb.re * rb;
		im = opa.im * ra + opb.im * rb;
		return *this;
	}
	inline Complex Sub(const Complex &opa, const Complex &opb)
	{
		re = opa.re - opb.re;
		im = opa.im - opb.im;
		return *this;
	}
	inline Complex Sub(const Complex &opa, const Complex &opb, real op)
	{
		re = (opa.re - opb.re) * op;
		im = (opa.im - opb.im) * op;
		return *this;
	}
	inline Complex Sub0(const Complex &opa, real ra, const Complex &opb, real rb)
	{
		re = opa.re * ra - opb.re * rb;
		im = opa.im * ra - opb.im * rb;
		return *this;
	}
	inline Complex Sub2(const Complex &opa, const Complex &opb, real op)
	{
		re = opa.re - opb.re * op;
		im = opa.im - opb.im * op;
		return *this;
	}
	inline Complex Copy(const Complex &op)
	{
		re = op.re;
		im = op.im;
		return *this;
	}

public:
	inline real arg() const
	{
		return ::Atan2(im, re);
	}
	inline real mod() const
	{
		return ::Sqrt(modSquared());
	}
	inline real modSquared() const
	{
		return re*re + im*im; 
	}
	Complex sqrt() const
	{ 
		real m = ::Sqrt(mod());
		real a = arg() / (real)2.;
		return Complex(m * Cos(a), m * Sin(a));
	}
	inline real abs() const
	{
		return mod();
	}
	Complex exp() const
	{
		return Complex(::Cos(im), ::Sin(im)) * ::Exp(re);
	}
	void Write(int file)
	{
		write(file, &re, sizeof(real));
		write(file, &im, sizeof(real));
	}
	void Read(int file)
	{
		read(file, &re, sizeof(real));
		read(file, &im, sizeof(real));
	}

//	Complex AddWithIm(const Complex &op) const		// (a + bi) + i*(c + di)
//	{
//		return Complex(re - op.im, im + op.re);
//	}
//      Complex SubWithIm(const Complex &op) const		// (a + bi) - i*(c + di)
//	{
//		return Complex(re + op.im, im - op.re);
//	}
        Complex MultRe(real op)
	{
		return Complex(re*op, im);
	}
        Complex MultIm(real op)
	{
		return Complex(re, im*op);
	}
	real Absc()
	{
		return Fabs(re) + Fabs(im); 
	}
};

#endif
