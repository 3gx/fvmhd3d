#ifndef __SIMD_H__
#define __SIMD_H__

namespace simd {
	enum {vlen = 4};

	typedef float     _v4sf __attribute__ ((vector_size(16)));
	typedef int       _v4si __attribute__ ((vector_size(16)));
	typedef long long _v2di __attribute__ ((vector_size(16)));
	typedef double    _v2df __attribute__ ((vector_size(16)));

	struct vfloat {
		_v4sf val;
		vfloat() {};
		vfloat(_v4sf _val) : val(_val) {}
		vfloat(const float f) {
			val = (_v4sf){f, f, f, f};
		}
		vfloat(const float a, const float b, const float c, const float d) {
			val = (_v4sf){a, b, c, d};
		}
		vfloat(const float *p, const _v4si idx) {
			val = vfloat(
					p[__builtin_ia32_vec_ext_v4si(idx, 0)],
					p[__builtin_ia32_vec_ext_v4si(idx, 1)],
					p[__builtin_ia32_vec_ext_v4si(idx, 2)],
					p[__builtin_ia32_vec_ext_v4si(idx, 3)]);
		}
		vfloat(const float *p, const bool aligned = true) {
			if (aligned) {
				val = *(_v4sf*)p;
			} else {
				val = __builtin_ia32_loadups(p);
			}
		}
		vfloat(const double *p, const bool aligned = true) {
			if (aligned) {
				val = __builtin_ia32_movlhps(
						__builtin_ia32_cvtpd2ps(*(_v2df*)(p + 0)),
						__builtin_ia32_cvtpd2ps(*(_v2df*)(p + 2)));
			} else {
				val = (_v4sf){p[0], p[1], p[2], p[3]};
			}
		}

		const vfloat& to_vfloat() const {return *this;} 

		const vfloat& operator =(const vfloat a) {val  = a.val; return *this;}
		const vfloat& operator+=(const vfloat a) {val += a.val; return *this;}
		const vfloat& operator-=(const vfloat a) {val -= a.val; return *this;}
		const vfloat& operator*=(const vfloat a) {val *= a.val; return *this;}
		const vfloat& operator/=(const vfloat a) {val /= a.val; return *this;}
		const vfloat operator -()  const {return vfloat(-val);}
		const vfloat operator +(const vfloat a) const {return vfloat(val + a.val);}
		const vfloat operator -(const vfloat a) const {return vfloat(val - a.val);}
		const vfloat operator *(const vfloat a) const {return vfloat(val * a.val);}
		const vfloat operator /(const vfloat a) const {return vfloat(val / a.val);}

		const vfloat operator!=(const vfloat a) const {
			return vfloat((_v4sf)__builtin_ia32_cmpneqps(val, a.val));
		}
		const vfloat operator==(const vfloat a) const {
			return vfloat((_v4sf)__builtin_ia32_cmpeqps(val, a.val));
		}

		const vfloat operator<(const vfloat &rhs) const{
			return (_v4sf)__builtin_ia32_cmpltps(val, rhs.val);
		}
		const vfloat operator<=(const vfloat &rhs) const{
			return (_v4sf)__builtin_ia32_cmpleps(val, rhs.val);
		}
		const vfloat operator>(const vfloat &rhs) const{
			return (_v4sf)__builtin_ia32_cmpgtps(val, rhs.val);
		}
		const vfloat operator>=(const vfloat &rhs) const{
			return (_v4sf)__builtin_ia32_cmpgeps(val, rhs.val);
		}
		const vfloat operator|(const vfloat &rhs) const{
			return __builtin_ia32_orps(val, rhs.val);
		}
		const vfloat operator&(const vfloat &rhs) const{
			return __builtin_ia32_andps(val, rhs.val);
		}

//		const vfloat operator &(const vfloat a) const {
//			return vfloat(__builtin_ia32_andps(val, a.val));
//		}

		const float operator[] (const int i) const {
			switch(i) {
				case 0:	return __builtin_ia32_vec_ext_v4sf(val, 0);
				case 1:	return __builtin_ia32_vec_ext_v4sf(val, 1);
				case 2:	return __builtin_ia32_vec_ext_v4sf(val, 2);
				case 3:	return __builtin_ia32_vec_ext_v4sf(val, 3);
				default: return 0.0f;	
			}
		}

		const float sum() const {
			const vfloat tmp  = __builtin_ia32_haddps(val, val);
			const vfloat tmp1 = __builtin_ia32_haddps(tmp, tmp);
			return __builtin_ia32_vec_ext_v4sf(tmp1, 0);
		}

		const float min() const {
			const float v0 = (*this)[0];
			const float v1 = (*this)[1];
			const float v2 = (*this)[2];
			const float v3 = (*this)[3];
			const float c01 = (v0 < v1) ? v0 : v1;
			const float c23 = (v2 < v3) ? v2 : v3;
			return (c01 < c23) ? c01 : c23;
		}

		const float max() const {
			const float v0 = (*this)[0];
			const float v1 = (*this)[1];
			const float v2 = (*this)[2];
			const float v3 = (*this)[3];
			const float c01 = (v0 > v1) ? v0 : v1;
			const float c23 = (v2 > v3) ? v2 : v3;
			return (c01 > c23) ? c01 : c23;
		}


		const vfloat rcp() const {
			return vfloat(1.0f)/(*this);
		}

#if 0
		const vfloat rcp_NR() const {
			_v4sf xm0 = val;
			_v4sf xm1 = __builtin_ia32_rcpps(xm0);
			xm0 = __builtin_ia32_mulps(xm0, xm1);
			xm0 = __builtin_ia32_mulps(xm0, xm1);
			xm1 = __builtin_ia32_addps(xm1, xm1);
			xm1 = __builtin_ia32_subps(xm1, xm0);
			return vfloat(xm1);
		}
#endif
#if 0
		const vfloat rcp_safe() const {
			_v4sf xm0 = val;
			_v4sf xm1 = __builtin_ia32_rcpps(xm0);
			xm0 = __builtin_ia32_mulps(xm0, xm1);
			xm0 = __builtin_ia32_mulps(xm0, xm1);
			xm1 = __builtin_ia32_addps(xm1, xm1);
			xm1 = __builtin_ia32_subps(xm1, xm0);
			return vfloat(xm1);
		}
#endif
		const vfloat sqrt() const {
			return vfloat(__builtin_ia32_sqrtps(val));
		}

#if 0
		const vfloat rsqrt() const {
			const vfloat x = val;
			const vfloat y = __builtin_ia32_rsqrtps(x);
			return (vfloat(-0.5f) * y) * (x*y*y + vfloat(-3.0f));
		}
#endif

		operator _v4sf() const{
			return val;
		}
		const _v2df as_v2df() const { // reinterpret cast
			return (_v2df)val;
		}
		const _v4si as_v4si() const{
			return (_v4si)val;
		}

		const _v2df dsum() const{
			_v2df vl = __builtin_ia32_cvtps2pd(val);
			_v2df vh = __builtin_ia32_cvtps2pd(__builtin_ia32_movhlps(val, val));
			return vl+vh;
		}

		static vfloat min(const vfloat x, const vfloat y) {
			return vfloat(__builtin_ia32_minps(x, y));
		}
		static vfloat max(const vfloat x, const vfloat y) {
			return vfloat(__builtin_ia32_maxps(x, y));
		}
		static vfloat abs(const vfloat x) {
			const _v4si mask = {0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff};
			return vfloat(__builtin_ia32_andps(x, (_v4sf)mask));
		}

#if 0
		const vfloat operator<(const vfloat rhs) const {
			return (_v4sf)__builtin_ia32_cmpltps(val, rhs.val);
		}
		const vfloat operator>(const vfloat rhs) const {
			return (_v4sf)__builtin_ia32_cmpgtps(val, rhs.val);
		}
#endif
		const vfloat select(const vfloat a, const vfloat b) const {
			const _v4sf tmp1 = __builtin_ia32_andps (val, a);
			const _v4sf tmp2 = __builtin_ia32_andnps(val, b);
			return vfloat((_v4sf)__builtin_ia32_orps(tmp1, tmp2));
		}
	};

	struct vint;

	struct vdouble {
		_v2df lo, hi;
		vdouble() {}
		//			vdouble(const vdouble val) : lo(val.lo), hi(val.hi)  {}
		vdouble(const _v2df _lo, const _v2df _hi) : lo(_lo), hi(_hi) {}
		vdouble(const double f)  {
			lo = hi = (_v2df){f, f};
		}
		vdouble(const double a, const double b, const double c, const double d) {
			lo = (_v2df){a,b};
			hi = (_v2df){c,d};
		}
		vdouble(const double *p, const bool aligned = true) {
			if (aligned) {
				lo = *(_v2df*)(p+0);
				hi = *(_v2df*)(p+2);
			} else {
				lo = __builtin_ia32_loadupd(p+0);
				hi = __builtin_ia32_loadupd(p+2);
			}
		}

		vdouble(const double *p, const vint);

		vdouble(const vfloat a) {
			lo = __builtin_ia32_cvtps2pd(_v4sf(a));
			hi = __builtin_ia32_cvtps2pd(__builtin_ia32_movhlps(_v4sf(a), _v4sf(a)));
		}

		const vdouble& operator =(const vdouble a) {lo  = a.lo; hi  = a.hi; return *this;}
		const vdouble& operator+=(const vdouble a) {lo += a.lo; hi += a.hi; return *this;}
		const vdouble& operator-=(const vdouble a) {lo -= a.lo; hi -= a.hi; return *this;}
		const vdouble& operator*=(const vdouble a) {lo *= a.lo; hi *= a.hi; return *this;}
		const vdouble& operator/=(const vdouble a) {lo /= a.lo; hi /= a.hi; return *this;}
		const vdouble operator -()  const {return vdouble(-lo, -hi);}
		const vdouble operator +(const vdouble a) const {return vdouble(lo + a.lo, hi + a.hi);}
		const vdouble operator -(const vdouble a) const {return vdouble(lo - a.lo, hi - a.hi);}
		const vdouble operator *(const vdouble a) const {return vdouble(lo * a.lo, hi * a.hi);}
		const vdouble operator /(const vdouble a) const {return vdouble(lo / a.lo, hi /a.hi);}

		const vdouble operator!=(const vdouble a) const {
			return vdouble(
					(_v2df)__builtin_ia32_cmpneqpd(lo, a.lo),
					(_v2df)__builtin_ia32_cmpneqpd(hi, a.hi));
		}
		const vdouble operator==(const vdouble a) const {
			return vdouble(
					(_v2df)__builtin_ia32_cmpeqpd(lo, a.lo),
					(_v2df)__builtin_ia32_cmpeqpd(hi, a.hi));
		}
		const vdouble operator &(const vdouble a) const {
			return vdouble(
					(_v2df)__builtin_ia32_andpd(lo, a.lo),
					(_v2df)__builtin_ia32_andpd(hi, a.hi));
		}

		const vdouble operator<(const vdouble &rhs) const{
			return vdouble(
					(_v2df)__builtin_ia32_cmpltpd(lo, rhs.lo),
					(_v2df)__builtin_ia32_cmpltpd(hi, rhs.hi));
		}
		const vdouble operator<=(const vdouble &rhs) const{
			return vdouble(
					(_v2df)__builtin_ia32_cmplepd(lo, rhs.lo),
					(_v2df)__builtin_ia32_cmplepd(hi, rhs.hi));
		}
		const vdouble operator>(const vdouble &rhs) const{
			return vdouble(
					(_v2df)__builtin_ia32_cmpgtpd(lo, rhs.lo),
					(_v2df)__builtin_ia32_cmpgtpd(hi, rhs.hi));
		}
		const vdouble operator>=(const vdouble &rhs) const{
			return vdouble(
					(_v2df)__builtin_ia32_cmpgepd(lo, rhs.lo),
					(_v2df)__builtin_ia32_cmpgepd(hi, rhs.hi));
		}
		const vdouble operator|(const vdouble &rhs) const{
			return vdouble(
					(_v2df)__builtin_ia32_orpd(lo, rhs.lo),
					(_v2df)__builtin_ia32_orpd(hi, rhs.hi));
		}
		const vdouble operator&(const vdouble &rhs) const{
			return vdouble(
					(_v2df)__builtin_ia32_andpd(lo, rhs.lo),
					(_v2df)__builtin_ia32_andpd(hi, rhs.hi));
		}
		const double operator[] (const int i) const {
			switch(i) {
				case 0:	return __builtin_ia32_vec_ext_v2df(lo, 0);
				case 1:	return __builtin_ia32_vec_ext_v2df(lo, 1);
				case 2:	return __builtin_ia32_vec_ext_v2df(hi, 0);
				case 3:	return __builtin_ia32_vec_ext_v2df(hi, 1);
				default: return 0.0;	
			}
		}

		static const _v2df select(const _v2df val, const _v2df a, const _v2df b) {
			const _v2df tmp1 = __builtin_ia32_andpd (val, a);
			const _v2df tmp2 = __builtin_ia32_andnpd(val, b);
			return (_v2df)__builtin_ia32_orpd(tmp1, tmp2);
		}
		const vdouble select(const vdouble a, const vdouble b) const {
			return vdouble(select(lo, a.lo, b.lo), select(hi, a.hi, b.hi));
		}

		const double sum() const {
			const _v2df tmp = hi + lo;
			const _v2df tmp1 = __builtin_ia32_haddpd(tmp, tmp);
			return __builtin_ia32_vec_ext_v2df(tmp1,0);
		}
		const double min() const {
			const _v2df tmp = __builtin_ia32_minpd(hi, lo);
			const double lhs = __builtin_ia32_vec_ext_v2df(tmp, 0);
			const double rhs = __builtin_ia32_vec_ext_v2df(tmp, 1);
		  return (lhs < rhs) ? lhs : rhs;	
		}

		const double max() const {
			const _v2df tmp = __builtin_ia32_maxpd(hi, lo);
			const double lhs = __builtin_ia32_vec_ext_v2df(tmp, 0);
			const double rhs = __builtin_ia32_vec_ext_v2df(tmp, 1);
		  return (lhs > rhs) ? lhs : rhs;	
		}
		const vdouble rcp() const {
			return vdouble(1.0)/(*this);
		}
#if 0
		static const _v2df rcp_NR(const _v2df val) {
			_v2df xm0 = val;
			_v2df xm1 = __builtin_ia32_rcppd(xm0);
			xm0 = __builtin_ia32_mulpd(xm0, xm1);
			xm0 = __builtin_ia32_mulpd(xm0, xm1);
			xm1 = __builtin_ia32_addpd(xm1, xm1);
			xm1 = __builtin_ia32_subpd(xm1, xm0);i
				return xm1;
		}
#else
		static const _v2df rcp_NR(const _v2df val) {
			return (_v2df){1.0, 1.0}/val;
		}
#endif
		const vdouble rcp_NR() const {
			return vdouble(rcp_NR(lo), rcp_NR(hi));
		}

		const vdouble sqrt() const {
			return vdouble(
					__builtin_ia32_sqrtpd(lo),
					__builtin_ia32_sqrtpd(hi));
		}

#if 0
		const vdouble rsqrt() const {
			const vdouble x = *this;
			const vdouble y = vdouble(__buildin_ia32_rsqrtpd(x.lo), __buildin_ia32_rsqrtpd(x.hi));
			return (vdouble(-0.5) * y) * (x*y*y + vdouble(-3.0));
		}
#else
		const vdouble rsqrt() const {
			return vdouble(1.0)/sqrt();
		}
#endif
		const _v4sf lo_as_v2df() const { // reinterpret cast
			return (_v4sf)lo;
		}
		const _v4sf hi_as_v2df() const { // reinterpret cast
			return (_v4sf)hi;
		}
		const _v4si lo_as_v4si() const{
			return (_v4si)lo;
		}
		const _v4si hi_as_v4si() const{
			return (_v4si)hi;
		}
		const vfloat to_vfloat() const {
			return	
				vfloat(
						__builtin_ia32_movlhps(
							__builtin_ia32_cvtpd2ps(lo),
							__builtin_ia32_cvtpd2ps(hi)));
		}
#if 0
		const _v2df dsum() const{
			_v2df vl = __builtin_ia32_cvtps2pd(val);
			_v2df vh = __builtin_ia32_cvtps2pd(__builtin_ia32_movhlps(val, val));
			return vl+vh;
		}
#endif
		static const vdouble min(const vdouble x, const vdouble y) {
			return vdouble(__builtin_ia32_minpd(x.lo, y.lo), __builtin_ia32_minpd(x.hi, y.hi));
		}
		static const vdouble max(const vdouble x, const vdouble y) {
			return vdouble(__builtin_ia32_maxpd(x.lo, y.lo), __builtin_ia32_maxpd(x.hi, y.hi));
		}
		static const vdouble abs(const vdouble x) {
			const _v2di mask = {0x7fffffffffffffff, 0x7fffffffffffffff};
			return vdouble(__builtin_ia32_andpd(x.lo, (_v2df)mask),	__builtin_ia32_andpd(x.hi, (_v2df)mask));
		}
	};

	struct vint {
		_v4si val;
		vint() {};
		vint(_v4si _val) : val(_val) {}
		vint(const int f) {
			val = (_v4si){f, f, f, f};
		}
		vint(const int a, const int b, const int c, const int d) {
			val = (_v4si){a, b, c, d};
		}
		vint(const int *p, const _v4si idx) {
			val = vint(
					p[__builtin_ia32_vec_ext_v4si(idx, 0)],
					p[__builtin_ia32_vec_ext_v4si(idx, 1)],
					p[__builtin_ia32_vec_ext_v4si(idx, 2)],
					p[__builtin_ia32_vec_ext_v4si(idx, 3)]);
		}
		vint(const int *p, const bool aligned = true) {
			if (aligned) {
				val = *(_v4si*)p;
			} else {
				val = (_v4si)__builtin_ia32_loadups((float*)p);
			}
		}
		vint(const long long p[]) {
			val = (_v4si){p[0], p[1], p[2], p[3]};
		}

		const vint& operator =(const vint a) {val  = a.val; return *this;}
		const vint& operator+=(const vint a) {val += a.val; return *this;}
		const vint& operator-=(const vint a) {val -= a.val; return *this;}
		const vint& operator*=(const vint a) {val *= a.val; return *this;}
		const vint& operator/=(const vint a) {val /= a.val; return *this;}
		const vint operator -()  const {return vint(-val);}
		const vint operator +(const vint a) const {return vint(val + a.val);}
		const vint operator -(const vint a) const {return vint(val - a.val);}
		const vint operator *(const vint a) const {return vint(val * a.val);}
		const vint operator /(const vint a) const {return vint(val / a.val);}

#if 0
		const vint operator!=(const vint a) const {
			return vint((_v4si)__builtin_ia32_pcmpeqd128(val, a.val));
		}
		const vint operator =(const vint a) const {
			return vint((_v4si)__builtin_ia32_pcmpeqd128(val, a.val));
		}
		const vint operator &(const vint a) const {
			return vint(__builtin_ia32_andps(val, a.val));
		}
#endif

		const int operator[] (const int i) const {
			switch(i) {
				case 0:	return __builtin_ia32_vec_ext_v4si(val, 0);
				case 1:	return __builtin_ia32_vec_ext_v4si(val, 1);
				case 2:	return __builtin_ia32_vec_ext_v4si(val, 2);
				case 3:	return __builtin_ia32_vec_ext_v4si(val, 3);
				default: return 0.0f;	
			}
		}

		const int sum() const {
			const int v0 = (*this)[0];
			const int v1 = (*this)[1];
			const int v2 = (*this)[2];
			const int v3 = (*this)[3];
			const int s01 = v0 + v1;
			const int s23 = v2 + v3;
			return s01 + s23;
		}

		const int min() const {
			const int v0 = (*this)[0];
			const int v1 = (*this)[1];
			const int v2 = (*this)[2];
			const int v3 = (*this)[3];
			const int c01 = (v0 < v1) ? v0 : v1;
			const int c23 = (v2 < v3) ? v2 : v3;
			return (c01 < c23) ? c01 : c23;
		}

		const int max() const {
			const int v0 = (*this)[0];
			const int v1 = (*this)[1];
			const int v2 = (*this)[2];
			const int v3 = (*this)[3];
			const int c01 = (v0 > v1) ? v0 : v1;
			const int c23 = (v2 > v3) ? v2 : v3;
			return (c01 > c23) ? c01 : c23;
		}

		operator _v4si() const{
			return val;
		}

#if 0
		static vint min(const vint x, const vint y) {
			return vint(__builtin_ia32_minps(x, y));
		}
		static vint max(const vint x, const vint y) {
			return vint(__builtin_ia32_maxps(x, y));
		}
		static vint abs(const vint x) {
			const _v4si mask = {0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff};
			return vint(__builtin_ia32_andps(x, (_v4si)mask));
		}
#endif
	};

	inline void transpose(
			const vfloat v0, const vfloat v1, const vfloat v2, const vfloat v3,
			vfloat &x, vfloat &y, vfloat &z, vfloat &w) 
	{
		_v4sf t0 = __builtin_ia32_unpcklps(v0, v2);
		_v4sf t1 = __builtin_ia32_unpckhps(v0, v2);
		_v4sf t2 = __builtin_ia32_unpcklps(v1, v3);
		_v4sf t3 = __builtin_ia32_unpckhps(v1, v3);

		x = __builtin_ia32_unpcklps(t0, t2);
		y = __builtin_ia32_unpckhps(t0, t2);
		z = __builtin_ia32_unpcklps(t1, t3);
		w = __builtin_ia32_unpckhps(t1, t3);
	}

	inline void transpose(const vfloat *vin, vfloat *vout) {
		_v4sf v0 = vin[0].val;
		_v4sf v1 = vin[1].val;
		_v4sf v2 = vin[2].val;
		_v4sf v3 = vin[3].val;
		_v4sf t0 = __builtin_ia32_unpcklps(v0, v2);
		_v4sf t1 = __builtin_ia32_unpckhps(v0, v2);
		_v4sf t2 = __builtin_ia32_unpcklps(v1, v3);
		_v4sf t3 = __builtin_ia32_unpckhps(v1, v3);

		_v4sf &x = vout[0].val;
		_v4sf &y = vout[1].val;
		_v4sf &z = vout[2].val;
		_v4sf &w = vout[3].val;

		x = __builtin_ia32_unpcklps(t0, t2);
		y = __builtin_ia32_unpckhps(t0, t2);
		z = __builtin_ia32_unpcklps(t1, t3);
		w = __builtin_ia32_unpckhps(t1, t3);
	}

	inline void transpose(vfloat *data) {
		const _v4sf v0 = data[0].val;
		const _v4sf v1 = data[1].val;
		const _v4sf v2 = data[2].val;
		const _v4sf v3 = data[3].val;
		_v4sf t0 = __builtin_ia32_unpcklps(v0, v2);
		_v4sf t1 = __builtin_ia32_unpckhps(v0, v2);
		_v4sf t2 = __builtin_ia32_unpcklps(v1, v3);
		_v4sf t3 = __builtin_ia32_unpckhps(v1, v3);

		_v4sf &x = data[0].val;
		_v4sf &y = data[1].val;
		_v4sf &z = data[2].val;
		_v4sf &w = data[3].val;

		x = __builtin_ia32_unpcklps(t0, t2);
		y = __builtin_ia32_unpckhps(t0, t2);
		z = __builtin_ia32_unpcklps(t1, t3);
		w = __builtin_ia32_unpckhps(t1, t3);
	}

	inline void transpose(const _v2df v0, const _v2df v1, _v2df &t0, _v2df &t1) {
		t0 = __builtin_ia32_unpcklpd(v0, v1);
		t1 = __builtin_ia32_unpckhpd(v0, v1);
	}

	inline void transpose(const vdouble *vin, vdouble *vout) {
		transpose(vin[0].lo, vin[1].lo, vout[0].lo, vout[1].lo);
		transpose(vin[0].hi, vin[1].hi, vout[2].lo, vout[3].lo);
		transpose(vin[2].lo, vin[3].lo, vout[0].hi, vout[1].hi);
		transpose(vin[2].hi, vin[3].hi, vout[2].hi, vout[3].hi);
#if 0	
		fprintf(stderr, "  %10.0f   %10.0f   %10.0f   %10.0f   \n", vin[0][0], vin[0][1], vin[0][2], vin[0][3]);
		fprintf(stderr, "  %10.0f   %10.0f   %10.0f   %10.0f   \n", vin[1][0], vin[1][1], vin[1][2], vin[1][3]);
		fprintf(stderr, "  %10.0f   %10.0f   %10.0f   %10.0f   \n", vin[2][0], vin[2][1], vin[2][2], vin[2][3]);
		fprintf(stderr, "  %10.0f   %10.0f   %10.0f   %10.0f   \n", vin[3][0], vin[3][1], vin[3][2], vin[3][3]);
		fprintf(stderr, " ---------------------\n");
		fprintf(stderr, "  %10.0f   %10.0f   %10.0f   %10.0f   \n", vout[0][0], vout[0][1], vout[0][2], vout[0][3]);
		fprintf(stderr, "  %10.0f   %10.0f   %10.0f   %10.0f   \n", vout[1][0], vout[1][1], vout[1][2], vout[1][3]);
		fprintf(stderr, "  %10.0f   %10.0f   %10.0f   %10.0f   \n", vout[2][0], vout[2][1], vout[2][2], vout[2][3]);
		fprintf(stderr, "  %10.0f   %10.0f   %10.0f   %10.0f   \n", vout[3][0], vout[3][1], vout[3][2], vout[3][3]);
		//		assert(false);
#endif
	}

	inline void transpose(vdouble *data) {
#if 0
		fprintf(stderr, "  %10.0f   %10.0f   %10.0f   %10.0f   \n", data[0][0], data[0][1], data[0][2], data[0][3]);
		fprintf(stderr, "  %10.0f   %10.0f   %10.0f   %10.0f   \n", data[1][0], data[1][1], data[1][2], data[1][3]);
		fprintf(stderr, "  %10.0f   %10.0f   %10.0f   %10.0f   \n", data[2][0], data[2][1], data[2][2], data[2][3]);
		fprintf(stderr, "  %10.0f   %10.0f   %10.0f   %10.0f   \n", data[3][0], data[3][1], data[3][2], data[3][3]);
		fprintf(stderr, " -----+++++++++++-----\n");
#endif
		_v2df vout2lo, vout3lo;
		transpose(data[0].lo, data[1].lo, data[0].lo, data[1].lo);
		transpose(data[0].hi, data[1].hi, vout2lo, vout3lo);
		transpose(data[2].lo, data[3].lo, data[0].hi, data[1].hi);
		transpose(data[2].hi, data[3].hi, data[2].hi, data[3].hi);
		data[2].lo = vout2lo;
		data[3].lo = vout3lo;
#if 0	
		fprintf(stderr, "  %10.0f   %10.0f   %10.0f   %10.0f   \n", data[0][0], data[0][1], data[0][2], data[0][3]);
		fprintf(stderr, "  %10.0f   %10.0f   %10.0f   %10.0f   \n", data[1][0], data[1][1], data[1][2], data[1][3]);
		fprintf(stderr, "  %10.0f   %10.0f   %10.0f   %10.0f   \n", data[2][0], data[2][1], data[2][2], data[2][3]);
		fprintf(stderr, "  %10.0f   %10.0f   %10.0f   %10.0f   \n", data[3][0], data[3][1], data[3][2], data[3][3]);
		fprintf(stderr, " ---------------------\n");
		//		assert(false);
#endif
	}

	inline void broadcast(const vfloat vin, vfloat *vout) {
		_v4sf &x = vout[0].val;
		_v4sf &y = vout[1].val;
		_v4sf &z = vout[2].val;
		_v4sf &w = vout[3].val;

		x = __builtin_ia32_shufps(vin, vin, 0x00);
		y = __builtin_ia32_shufps(vin, vin, 0x55);
		z = __builtin_ia32_shufps(vin, vin, 0xaa);
		w = __builtin_ia32_shufps(vin, vin, 0xff);

	}

	inline const simd::vfloat rcp_safe(const simd::vfloat x) {
#ifndef _FPESIG_ENABLE_
		return (x != simd::vfloat(0.0f)).select(x.rcp(), simd::vfloat(0.0f));
#else
		const vfloat tmp = (x != vfloat(0.0f)).select(x, vfloat(HUGE));
		return tmp.rcp();
#endif
	}


	inline const simd::vdouble rcp_safe(const simd::vdouble x) {
#ifndef _FPESIG_ENABLE_
		return (x != simd::vdouble(0.0)).select(x.rcp(), simd::vdouble(0.0));
#else
		const vdouble tmp = (x != vdouble(0.0)).select(x, vfloat(HUGE));
		return tmp.rcp();
#endif
	}
	inline vdouble::vdouble(const double *p, const vint i) {
		lo = (_v2df){p[i[0]], p[i[1]]};
		hi = (_v2df){p[i[2]], p[i[3]]};
	}

	inline const simd::vfloat min(const simd::vfloat x, const simd::vfloat y) {
		return vfloat::min(x,y);
	}

	inline const simd::vfloat max(const simd::vfloat x, const simd::vfloat y) {
		return vfloat::max(x,y);
	}

	inline const simd::vdouble min(const simd::vdouble x, const simd::vdouble y) {
		return vdouble::min(x,y);
	}

	inline const simd::vdouble max(const simd::vdouble x, const simd::vdouble y) {
		return vdouble::max(x,y);
	}

	inline const vdouble abs(const vdouble x) {
		const _v2di mask = {0x7fffffffffffffff, 0x7fffffffffffffff};
		return vdouble(__builtin_ia32_andpd(x.lo, (_v2df)mask),	__builtin_ia32_andpd(x.hi, (_v2df)mask));
	}
	inline const vfloat abs(const vfloat x) {
		const _v4si mask = {0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff};
		return vfloat(__builtin_ia32_andps(x, (_v4sf)mask));
	}

	inline const vdouble sqrt(const vdouble x) {
		return vdouble(
				__builtin_ia32_sqrtpd(x.lo),
				__builtin_ia32_sqrtpd(x.hi));
	}

	inline const vfloat sqrt(const vfloat x) {
		return vfloat(__builtin_ia32_sqrtps(x.val));
	}

	inline const vdouble rcp(const vdouble x) {
		return vdouble(1.0)/x;
	}

	inline const vfloat rcp(const vfloat x) {
		return vfloat(1.0)/x;
	}

	inline const vdouble rsqrt(const vdouble x) {
		return rcp(sqrt(x));
	}

	inline void swap (vdouble &v0, vdouble &v1) {
		const vdouble tmp = v0;
		v0 = v1;
		v1 = tmp;
	}

};

#endif // __SIMD_H__

