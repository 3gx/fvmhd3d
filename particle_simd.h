#ifndef __PARTICLE_H__
#define __PARTICLE_H__

template<class DOUBLE, class INT, int VLEN>
struct ParticleBase {
	enum {ACTIVE_BIT = (1LLU << 31)};
	enum {REFINE_BIT = (1LLU << 30), DEREFINE_BIT = (1LLU << 29)};

	DOUBLE  x,  y,  z;           // 3
	DOUBLE vx, vy, vz;           // 6
	DOUBLE ax, ay, az, tmp;    
	DOUBLE rmax;                 // 7
	INT    idx[VLEN];            // 8
	INT    rung[VLEN];           // 9
	DOUBLE tlast;                // 10
	DOUBLE vrel2max;             // 11
	DOUBLE Volume;               // 12
	// Total: 512 = 128x4 bits
	ParticleBase() {}; 
	ParticleBase(const DOUBLE &_x, const DOUBLE &_y, const DOUBLE &_z) : x(_x), y(_y), z(_z) {}
	ParticleBase(
			const DOUBLE & _x, const DOUBLE & _y, const DOUBLE & _z,
			const DOUBLE &_vx, const DOUBLE &_vy, const DOUBLE &_vz) :
		x(_x), y(_y), z(_z), vx(_vx), vy(_vy), vz(_vz)  {}
	ParticleBase(
			const DOUBLE & _x, const DOUBLE & _y, const DOUBLE & _z,
			const DOUBLE &_vx, const DOUBLE &_vy, const DOUBLE &_vz,
			const INT &_idx) :
		x(_x), y(_y), z(_z), vx(_vx), vy(_vy), vz(_vz)  {
			for (int ch = 0; ch < VLEN; ch++)
				idx[ch] = _idx;
		}
	inline const bool is_active() const;
	inline void set_inactive();
	inline void set_active();
	
	inline const bool is_derefine() const;
	inline void unset_derefine();
	inline void set_derefine();
	~ParticleBase() {};
};


template<>
inline const bool ParticleBase<double, unsigned long long, 1>::is_active() const 
{
	return ((idx[0] & ACTIVE_BIT) == ACTIVE_BIT);
}
	template<>
inline void ParticleBase<double, unsigned long long, 1>::set_inactive() 
{
	idx[0] &= ~ACTIVE_BIT;
}
	template<>
inline void ParticleBase<double, unsigned long long, 1>::set_active() 
{
	idx[0] |= ACTIVE_BIT;
}

template<>
inline const bool ParticleBase<double, unsigned long long, 1>::is_derefine() const 
{
	return ((idx[0] & DEREFINE_BIT) == DEREFINE_BIT);
}
	template<>
inline void ParticleBase<double, unsigned long long, 1>::unset_derefine() 
{
	idx[0] &= ~DEREFINE_BIT;
}
	template<>
inline void ParticleBase<double, unsigned long long, 1>::set_derefine() 
{
	idx[0] |= DEREFINE_BIT;
}


struct ParticleSimd : ParticleBase<simd::vdouble, unsigned long long, simd::vlen> {
	ParticleSimd() {};
	~ParticleSimd() {};

	ParticleSimd(const ParticleBase<double, unsigned long long, 1> &p) {
		x    = simd::vdouble(p.x );
		y    = simd::vdouble(p.y );
		z    = simd::vdouble(p.z );
		vx   = simd::vdouble(p.vx);
		vy   = simd::vdouble(p.vy);
		vz   = simd::vdouble(p.vz);
		ax   = simd::vdouble(p.ax);
		ay   = simd::vdouble(p.ay);
		az   = simd::vdouble(p.az);
		rmax = simd::vdouble(p.rmax);
		for (int ch = 0; ch < simd::vlen; ch++)
		{
			idx [ch] = p.idx [0];
			rung[ch] = p.rung[0];
		}
		tlast    = simd::vdouble(p.tlast   );
		vrel2max = simd::vdouble(p.vrel2max);
		Volume     = simd::vdouble(p.Volume);
	};

	ParticleSimd(const ParticleBase<double, unsigned long long, 1> *p) {
		asm("#vload_beg");
		simd::vdouble *vin  = (simd::vdouble *)p;
		simd::vdouble *vout = (simd::vdouble *)this;

		const int stride = sizeof(ParticleBase<double, unsigned long long, 1>);
		const int nparts = stride/sizeof(simd::vdouble);

		assert(nparts * sizeof(simd::vdouble) == stride);

		static simd::vdouble el[simd::vlen];

		asm("#test1");
		for (int i = 0, j = 0; i < nparts; i++, j += simd::vlen) {
			asm("#innloop-beg");
			for (int ch = 0, offset = 0; ch < simd::vlen; ch++, offset += nparts) {
				el[ch] = vin[i + offset];
			}
			asm("#innloop-end");
			asm("#transpose1");
			simd::transpose(el, &vout[j]);
			asm("#transpose2");
		}
		asm("#test2");

	};

#if 1
	ParticleSimd(const ParticleBase<double, unsigned long long, 1> *p, const simd::vint addr) {
		asm("#ParticleSimd_vgather_beg");
		simd::vdouble *in  = (simd::vdouble *)p;
		simd::vdouble *out = (simd::vdouble *)this;

		const int stride = sizeof(ParticleBase<double, unsigned long long, 1>);
		const int nparts = stride/sizeof(simd::vdouble);
		assert(nparts * sizeof(simd::vdouble) == stride);

		simd::vint offset = addr * simd::vint(nparts);
		static simd::vdouble tmp[simd::vlen];

		for (int i = 0, j = 0; i < nparts; i++, j += simd::vlen, offset += simd::vint(1)) {
#if 1
			for (int ch = 0; ch < simd::vlen; ch++) {
#if 0
				if (offset[ch] < 0) {
					fprintf(stderr, " ch= %d  offset[ch]= %d [ %g ]\n",
							ch, offset[ch], offset[ch] * (1.0/nparts));
					assert(offset[ch] < 0);
				}
#endif
				tmp[ch] = in[offset[ch]];
			}
#else
			tmp[0] = in[offset[0]];
			tmp[1] = in[offset[1]];
			tmp[2] = in[offset[2]];
			tmp[3] = in[offset[3]];
#endif
			simd::transpose(tmp, &out[j]);
		}
		asm("#ParticleSimd_vgather_end");

	};
#else
	ParticleSimd(const ParticleBase<double, unsigned long long, 1> *p, const simd::vint addr) {
		asm("#ParticleSimd_vgather_beg");

		const int stride = sizeof(ParticleBase<double, unsigned long long, 1>)/sizeof(double);
		const int nparts = stride/simd::vlen;
		assert(nparts * simd::vlen == stride);

		simd::vdouble  *in = (simd::vdouble*)p;
		simd::vdouble *out = (simd::vdouble*)this;

		const simd::vint voffset = addr * simd::vint(nparts);
		const int *offset = (const int*)&voffset;

		static simd::vdouble  tmp[stride];
		switch(simd::vlen) {
			case 4:
				asm("#read_tmp_data");
				for (int ch = 0, i = 0; ch < simd::vlen; ch++, i += 2) {
					tmp[i + 0] = in[offset[ch] + 0];
					tmp[i + 1] = in[offset[ch] + 1];
				}

				asm("#shuffle_tmp_data");
				// shuffle data in memory
				for (int i = 0, j = 0; i < stride; i += simd::vlen, j++) {
					out[i + 0] = tmp[j + 0*nparts];
					out[i + 1] = tmp[j + 1*nparts];
					out[i + 2] = tmp[j + 2*nparts];
					out[i + 3] = tmp[j + 3*nparts];
				}
				break;

			default:
				asm("#read_tmp_data");
				for (int j = 0; j < nparts; j++) {
					for (int ch = 0, i = 0; ch < simd::vlen; ch++, i += nparts) {
						tmp[i + j] = in[offset[ch] + j];
					}
				}

				asm("#shuffle_tmp_data");
				// shuffle data in memory
				for (int i = 0, j = 0; i < stride; i += simd::vlen, j++) {
					for (int ch = 0, offset = 0; ch < simd::vlen; ch++, offset += nparts)
						out[i + ch] = tmp[j + offset];
				}
		}


		// transpose temp array

		for (int i = 0; i < stride; i += simd::vlen) {
			simd::transpose(&out[i]);
		}

		asm("#ParticleSimd_vgather_end");
	}

#endif

};

typedef ParticleBase<double, unsigned long long, 1> Particle;
typedef ParticleSimd         Particle_simd;
#endif // __PARTICLE_H__

