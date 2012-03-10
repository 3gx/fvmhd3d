#ifndef __PFLOAT_H__
#define __PFLOAT_H__

#include <cassert>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include "vector3.h"

#define LARGEU 0xffffffffu

template <int ch> 
struct pfloat {
  unsigned int uval;
  static double xscale_i2f,  xscale_f2i, scale;
  
  static void set_scale(const double _scale)
	{
		scale = _scale;
		assert(scale > 0.0);
    xscale_i2f = scale / double(LARGEU);
    xscale_f2i = 1.0   / xscale_i2f;
  }
	pfloat() {}
  pfloat(const double val) 
	{
#if 0
    const double tmp = std::min(std::max(0.0, val), scale) * xscale_f2i;
#else
    const double tmp = val * xscale_f2i;
    if      (val <  0.0  ) uval =                                  0;
		else if (val >= scale) uval =                                  0xffffffffu;
		else                   uval = static_cast<unsigned int>(tmp) & 0xffffffffu;
#endif
  }
	pfloat(const unsigned int _uval) : uval(_uval) {}

	const double to_double() const {return xscale_i2f * static_cast<double>(uval);}
  double operator-(const pfloat rhs) const 
	{
    const unsigned int s = (uval - rhs.uval);
    return xscale_i2f * static_cast<double>((int)s);
  }

	friend const pfloat add_half(const pfloat x, const pfloat y)
	{
		unsigned long long uval = (const unsigned long long)x.uval + (const unsigned long long)y.uval;
		uval = uval >> 1;
		return pfloat(static_cast<unsigned int>(uval));	
	}	
	pfloat& dec() 
	{
#if 1
		const unsigned int lhs = uval - 0x00000002u;
    uval = (lhs <= uval) ? lhs : 0;
#endif
    return *this;
	}
 	pfloat& half() 
	{
    const unsigned int lhs = uval >> 1; 
    uval = (lhs <= uval) ? lhs : 0;
    return *this;
  };
  pfloat& twice() 
	{
    const unsigned int lhs = uval << 1;
    uval = (lhs >= uval) ? lhs : 0xffffffffu;
    return *this;
  };
  pfloat& add(const pfloat rhs) 
	{
    const unsigned int lhs = uval + rhs.uval;
    uval = (lhs >= uval) ? lhs : 0xffffffffu;
    return *this;
  }
  pfloat& sub(const pfloat rhs) 
	{
    const unsigned int lhs = uval - rhs.uval;
    uval = (lhs <= uval) ? lhs : 0;
    return *this;
  }
  pfloat& padd(const pfloat rhs) 
	{
    uval += rhs.uval;
    return *this;
  }
  pfloat& psub(const pfloat rhs) 
	{
    uval -= rhs.uval;
    return *this;
  }
  pfloat& div(const unsigned int rhs) 
	{
    const unsigned int lhs = uval/rhs;
    uval = (lhs <= uval) ? lhs : 0;
    return *this;
  }


  const bool aleq(const pfloat rhs) const 
	{
    const int           s = uval;
    const unsigned int  p = abs(s); // & -2; ///? may not work here
    return p <= rhs.uval;
  }
  bool ale(const pfloat rhs) const 
	{
    const int           s = uval;
    const unsigned int  p = abs(s); // & -2; ///? may not work here
    return p < rhs.uval;
  }
  bool leq(const pfloat rhs) const 
	{
    return uval <= rhs.uval;
  }

#if 0	
  pfloat& min(const pfloat rhs) 
	{
    uval = std::min(uval, rhs.uval);
    return *this;
  }
  pfloat& max(const pfloat rhs) {
    uval = std::max(uval, rhs.uval);
    return *this;
  }
#else
  friend const pfloat min(const pfloat lhs, const pfloat rhs) 
	{
    return pfloat(std::min(lhs.uval, rhs.uval));
	}
  friend const pfloat max(const pfloat lhs, const pfloat rhs) 
	{
    return pfloat(std::max(lhs.uval, rhs.uval));
	}
#endif


#if 0
  pfloat& add(const double x) {
    double xh = x / 2.0;
    while (xh < xmin) xh += (xmax - xmin);
    while (xh > xmax) xh -= (xmax - xmin);
    const double tmp = (xh - xmin) * xscale_f2i;
    const unsigned int itmp = static_cast<unsigned int>(floor(tmp));
    uval += (itmp << 1);
    return *this;
  }
#endif

  friend void pfloat_merge(
			const pfloat x1, const pfloat h1, 
			const pfloat x2, const pfloat h2,
			pfloat &xc, pfloat &hs) 
	{
		unsigned long long_shift = 0x100000000lu;
		//     long_shift += 0xu;
		//     long_shift = 0;

		const unsigned long ux1 = x1.uval + long_shift;
		const unsigned long uh1 = h1.uval;
		const unsigned long ux2 = x2.uval + long_shift;
		const unsigned long uh2 = h2.uval;

		//     fprintf(stderr, "long_shift= 0x%.16lx\n", long_shift);

		const unsigned long min1 = ux1 - uh1;
		const unsigned long max1 = ux1 + uh1;
		const unsigned long min2 = ux2 - uh2;
		const unsigned long max2 = ux2 + uh2;

		//     fprintf(stderr, "min1= 0x%.16lx\n", min1);
		//     fprintf(stderr, "max1= 0x%.16lx\n", max1);
		//     fprintf(stderr, "min2= 0x%.16lx\n", min2);
		//     fprintf(stderr, "max2= 0x%.16lx\n", max2);

		unsigned long minc = std::min(min1, min2);
		unsigned long maxc = std::max(max1, max2);

		//     fprintf(stderr, "minc= 0x%.16lx\n", minc);
		//     fprintf(stderr, "maxc= 0x%.16lx\n", maxc);


		minc = minc >> 1;
		maxc = maxc >> 1;

		const unsigned long uxc = maxc + minc - long_shift;
		const unsigned long uhs = maxc - minc;


		xc.uval = uxc;
		hs.uval = uhs;

		//     xc.uval = x2.uval;
		//     hs.uval = h2.uval;


	}

};

template <int ch> double pfloat<ch>::scale      = 1.0;
template <int ch> double pfloat<ch>::xscale_i2f = 1.0 / double(LARGEU);
template <int ch> double pfloat<ch>::xscale_f2i = 1.0 * double(LARGEU);

struct pfloat3 {
	pfloat<0> x;
	pfloat<1> y;
	pfloat<2> z;

	pfloat3(const dvec3 &p) : 
		x(pfloat<0>(p.x)), y(pfloat<1>(p.y)), z(pfloat<2>(p.z)) {}

	const dvec3 to_dvec3() const 
	{
		return dvec3(x.to_double(), y.to_double(), z.to_double());
	}


	pfloat3() {}
	~pfloat3() {}
};

#endif // _PERIODIC_H_

