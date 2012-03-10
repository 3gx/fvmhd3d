#ifndef __FLUID__H__
#define __FLUID__H__

template<int FLUIDNSCALARS, class REAL>
struct FluidBase {
		enum {NMHD   = 10};
		enum {NSCALARS = FLUIDNSCALARS};
		enum {NFLUID = NMHD + NSCALARS};
		enum {
			MASS  = 0,
			MOMX  = 1,
			MOMY  = 2,
			MOMZ  = 3,
			ENER  = 4,
			WBX   = 5,
			WBY   = 6,
			WBZ   = 7,
			MPSI  = 8,
			MENTR = 9
		};
		enum {
			DENS = 0,
			VELX = 1,
			VELY = 2,
			VELZ = 3,
			ETHM = 4,
			BX   = 5,
			BY   = 6,
			BZ   = 7,
			PSI  = 8,
			ENTR = 9
		};
		REAL  data[NFLUID];

		FluidBase(const REAL v) {
			for (int k = 0; k < NFLUID; k++)
				data[k] = v; 
		}
		FluidBase(const FluidBase &v) {
			for (int k = 0; k < NFLUID; k++)
				data[k] = v.data[k];
		}
		FluidBase() {}
		~FluidBase() {}

    const vec3 get_vel() const {return vec3(data[VELX], data[VELY], data[VELZ]);}
    const vec3 get_B()   const {return vec3(data[  BX], data[  BY], data[  BZ]);}

		const REAL& operator[](const int i) const {return data[i];}
		      REAL& operator[](const int i)       {return data[i];}

		const REAL& scal(const int i) const {return data[NMHD + i];}
		      REAL& scal(const int i)       {return data[NMHD + i];}

	  const FluidBase to_primitive(const REAL volume) const {
			FluidBase prim;

			const REAL iM = REAL(1.0)/data[MASS];
			const REAL iV = REAL(1.0)/volume;

			prim[DENS] = data[MASS] * iV;
			prim[VELX] = data[MOMX] * iM;
			prim[VELY] = data[MOMY] * iM;
			prim[VELZ] = data[MOMZ] * iM;
			prim[BX  ] = data[WBX ] * iV;
			prim[BY  ] = data[WBY ] * iV;
			prim[BZ  ] = data[WBZ ] * iV;
#if 1
			prim[PSI ] = data[MPSI] * iM;
#else
			prim[PSI ] = data[MPSI];
#endif
			prim[ENTR] = data[MENTR] * iM;

#ifdef __ENER_UB__
			prim[ETHM] = data[ENER] * iV;
#ifndef __ENER_U__
			prim[ETHM] -= REAL(0.5) * (sqr(prim[BX]) + sqr(prim[BY]) + sqr(prim[BZ]));
#endif
#else
			prim[ETHM] = 
				data[ENER] * iV - 
				REAL(0.5) * (sqr(prim[  BX]) + sqr(prim[  BY]) + sqr(prim[  BZ])) -
				REAL(0.5) * (sqr(prim[VELX]) + sqr(prim[VELY]) + sqr(prim[VELZ])) * prim[DENS];

#endif

			for (int i = 0; i < NSCALARS; i++) 
				prim[NMHD + i] = data[NMHD + i] * iM;

			return prim;
		}

		const FluidBase to_conservative(const REAL volume) const {
			FluidBase cons;

			const REAL V = volume;
			const REAL M = data[DENS] * volume;

			cons[MASS] = data[DENS] * V;
			cons[MOMX] = data[VELX] * M;
			cons[MOMY] = data[VELY] * M;
			cons[MOMZ] = data[VELZ] * M;
			cons[ WBX] = data[  BX] * V;
			cons[ WBY] = data[  BY] * V;
			cons[ WBZ] = data[  BZ] * V;
#if 1
			cons[MPSI] = data[ PSI] * M;
#else
			cons[MPSI] = data[ PSI] ;
#endif
			cons[MENTR] = data[ENTR] * M;

#ifdef __ENER_UB__
			cons[ENER] = data[ETHM] * V;
#ifndef __ENER_U__
			cons[ENER] += REAL(0.5) * (sqr(data[BX]) + sqr(data[BY]) + sqr(data[BZ])) * V;
#endif
#else
			cons[ENER] = data[ETHM] * V +
				REAL(0.5) * (sqr(data[  BX]) + sqr(data[  BY]) + sqr(data[  BZ])) * V +
				REAL(0.5) * (sqr(data[VELX]) + sqr(data[VELY]) + sqr(data[VELZ])) * M;
#endif

			for (int i = 0; i < NSCALARS; i++)
				cons[NMHD + i] = data[NMHD + i] * M;

			return cons;
		}
};


typedef FluidBase<NEXTRASCALARS, real > Fluid;

struct Fluid_st
{
	int   bnd;
	vec3  pos, vel;
  real  dt;
	Fluid w;    
#if 0
  vec3 J;
#endif
	Fluid_st() {}
};

struct Fluid_rec
{ 
	int bnd;
	real tlast;      // time of the last update, negative if inactive
	vec3 pos, vel;
#if 0
	vec3 J, Jdot; // current
	vec3 c0;      // reconstruction centre
#endif

	Fluid w, x, y, z, t;
	Fluid_rec() {}
  Fluid_rec(const Fluid &_w) : w(_w), x(0.0), y(0.0), z(0.0), t(0.0) {}
	const bool is_active() const {return tlast >= 0;}
};

struct FluidD
{
	Fluid U;
	vec3  gradPsi;
	real  divB;
	FluidD() {};
	FluidD(const real v) : U(v), gradPsi(v), divB(v) {}
	FluidD(const Fluid &_U, const vec3 &_gradPsi, const real _divB) :
		U(_U), gradPsi(_gradPsi), divB(_divB) {}
};

struct FluidExtra
{
	vec3 J;
};


#endif // __FLUID__H__

