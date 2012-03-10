#ifndef __FVMHD3D3D_H__
#define __FVMHD3D3D_H__

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <map>
#include <cmath>
#include <cassert>

#define NMAXPROC 256

#include "timerlib.h"
#include "mytimer.h"
#include "Distribute.h"
#include "myMPI.h"
// #include "simd.h"
#include "vector3.h"
#include "Boundary.h"
#include "Octree.h"
#include "bOctree.h"
#include "Scheduler.h"

#if 1
template<class T>
inline void clear_vec(std::vector<T> &vec) {std::vector<T>().swap(vec);}
#else
template<class T>
inline void clear_vec(std::vector<T> &vec) {vec.clear();}
#endif

template<class T>
inline void fit_vec(std::vector<T> &vec) {std::vector<T>(vec).swap(vec);}

#if 0
	template<class T>
inline void fit_reserve_vec(std::vector<T> &vec, const int nres)
{
	const int size = vec.size();
	vec.resize(std::max(size, nres));
	fit_vec(vec);
	vec.resize(size);
}
#endif

typedef unsigned int uint128_t __attribute__((mode(TI)));

template<int NBIT>
struct BitOps
{
	unsigned long long key[NBIT >> 6];     // NBIT/128
	BitOps() {
		assert(NBIT % sizeof(unsigned long long) == 0);
		clear();
	}
	void clear()
	{
		const unsigned int Nel = NBIT >> 6;
		for (unsigned int i = 0; i < Nel; i++)
			key[i] = 0;
	}
	void operator|=(const int bit)
	{
		assert(bit < NBIT);
		key[bit >> 6] |= (1ULL << (bit & 0x3F));
	}
	void operator&=(const int bit)
	{
		assert(bit < NBIT);
		key[bit >> 6] &= (1ULL << (bit & 0x3F));
	}
	void operator!=(const int bit)
	{
		assert(bit < NBIT);
		key[bit >> 6] != (1ULL << (bit & 0x3F));
	}
	void operator^=(const int bit)
	{
		assert(bit < NBIT);
		key[bit >> 6] ^= (1ULL << (bit & 0x3F));
	}
	void clear(const int bit)
	{
		assert(bit < NBIT);
		key[bit >> 6] &= ~(1ULL << (bit & 0x3F));
	}
	const bool is_set(const int bit) const
	{
		assert(bit < NBIT);
		return ((key[bit >> 6] & (1ULL << (bit & 0x3F))) == (1ULL << (bit & 0x3F)));
	}
	void set(const int bit)
	{
		assert(bit < NBIT);
		key[bit >> 6] |= (1ULL << (bit & 0x3F));
	}
};

namespace fvmhd3d 
{
  enum {NEXTRASCALARS = 0};

	typedef double real;
	typedef vector3<real> vec3;
	typedef Boundary<vec3> boundary;

	inline real sqr(const real x) {return x*x;}

#include "Particle.h"
#include "Fluid.h"

	struct ParticleFluidStruct
	{
		Particle p;
		Fluid  U;
		FluidD dU;
    Fluid_rec Wrec;
		ParticleFluidStruct() {}
		ParticleFluidStruct(const Particle &_p, const Fluid &_U, const FluidD &_dU, const Fluid_rec &_Wrec) :
			p(_p), U(_U), dU(_dU), Wrec(_Wrec) {}
	};
	
  struct ParticleFluidStructLite
	{
		Particle p;
		Fluid  U;
		FluidD dU;
		ParticleFluidStructLite() {}
		ParticleFluidStructLite(const Particle &_p, const Fluid &_U, const FluidD &_dU) :
			p(_p), U(_U), dU(_dU) {}
	};


#if 0
#include "Cell.h"
#include "Face.h"
#else
	struct SiteImport
	{
		protected:
			vec3  _pos;
			int   _proc;
			int   _id;
			real  _rmax;
		public:
			SiteImport(const vec3 &__pos, const int &__proc, const int &__id, const real &__rmax) :
				_pos(__pos), _proc(__proc), _id(__id), _rmax(__rmax) {
					assert((_id & 0x30000000) == 0);
				}
			const int id() const {return _id & 0x0FFFFFFF;}
			const int proc() const {return _proc;}
			const bool is_active() const {return (_id & 0x80000000) == 0x80000000;}
      void   set_active() {_id |=  (0x80000000);}
      void unset_active() {_id &= ~(0x80000000);}
			const bool is_ngb   () const {return (_id & 0x40000000) == 0x40000000;}
			void   set_local_ngb() {_id |=  (0x20000000);}
			void unset_local_ngb() {_id &= ~(0x20000000);}
			const bool is_local_ngb() const {return (_id & 0x20000000) == 0x20000000;}
			void set_hasngb() {_id |= 0x10000000;}
			void unset_hasngb() {_id &= ~(0x10000000);}
			const bool is_hasngb() const {return (_id & 0x10000000) == 0x10000000;}
			const vec3& pos() const {return _pos;}
			vec3& pos() {return _pos;}
			const real& rmax() const {return _rmax;}
			real& rmax() {return _rmax;}
	};

	struct Face
	{
		protected:
		public:
			int s1, s2;
			vec3 centroid;
			vec3 n;
			Face() {}
			template<bool DEBUG>
				int ngb(const int i) const {
					if (DEBUG)
					{
						if      (i == s1) return s2;
						else if (i == s2) return s1;

						assert(false);
						return -1;	
					}
					else
						return i == s1 ? s2 : s1;
				};
			const real area()   const {return n.abs();}
	};

	struct Ngb
	{
		int proc, remote_id, id;
		Ngb() {}
		Ngb(const int _proc, const int _id) :
			proc(_proc), remote_id(_id), id(-1) {};
	};

	template<typename T>
	struct Neighbours
	{
		protected:
			std::vector<T> list;
		public:
			Neighbours()                     {list.reserve(32);}
			void push_back(const T v)        {list.push_back(v);}
			size_t size()  const             {return list.size();}
			T  operator[](const int i) const {return list[i];}
			T& operator[](const int i)       {return list[i];}
#if 0
			void clear()                     {std::vector<T>().swap(list); list.reserve(64);}
			~Neighbours()                    {std::vector<T>().swap(list);}
#else
			void clear()                     {list.clear();}
			~Neighbours()                    {list.clear();}
#endif
	};

	struct Cell
	{
		real Volume;
		vec3 centroid;
		Neighbours<int> ngb;
		Cell() : Volume(-1.0) {}
		Cell(const real v, const vec3 &c) : Volume(v), centroid(c) {}
		const Neighbours<int>& faces() const 
    {
#if 0
      assert(ngb.size() > 0);
      assert(ngb[0] >= 0);
#endif
      return ngb;
    }
		Neighbours<int>& faces() 
    {
#if 0
      assert(ngb.size() > 0);
      assert(ngb[0] >= 0);
#endif
      return ngb;
    }
    const bool is_ngb() const 
    {
      assert(ngb.size() > 0);
      return (ngb[0] < 0);
    }
	};

#endif

#include "Energy.h"

	typedef double TREAL;
	typedef vector3<TREAL> TVEC3;
	typedef Boundary<TVEC3> TBOUNDARY;


	struct LightSite;
	struct Site;
	struct ImportSite;

	struct LightSite
	{
		TVEC3 pos;
		int  id;
		LightSite(const TVEC3 &_pos, const int _id) : pos(_pos), id(_id) {}
	};

	struct Site
	{
		TVEC3 pos;
		int id;
		TREAL rmax;
		bool failed;
		Site() {}
		Site(const TVEC3 &_pos, const int _id, const TREAL _rmax, const bool _failed = true) :
			pos(_pos), id(_id), rmax(_rmax), failed(_failed) {}
	};

	struct ImportSite
	{
		int proc;
		int id;
		int oct;
		TVEC3 pos;
		ImportSite(const int _proc, const int _id, const int _oct, const TVEC3 &_pos) :
			proc(_proc), id(_id), oct(_oct), pos(_pos) {};
		ImportSite() : proc(-1) {}
	};


	inline int int_abs(const int i) {return i < 0 ? -1-i : i;}
	inline int int_map(const int i) {return -1-i;}


	struct system
	{

		enum {NLEAF_LOC = 128, NLEAF_GLB = 8};

	
		///////////////////////
		
		std::vector<Particle> ptcl_local, ptcl_import;
		std::vector<FluidExtra> Wextra_local, Wextra_import;
		std::vector<float> divBi;

		std::vector<Fluid    >    U_local,    U_import;
		std::vector<FluidD   >   dU_local,   dU_import;
		std::vector<Fluid_rec> Wrec_local, Wrec_import;
    std::vector<Fluid_st > Wst_import;
		
		///////////////////////

		std::vector<SiteImport> site_import;
			
		std::vector<int> site_active_list;
		std::vector<int> site_with_ngb_active_list;

		std::vector< Neighbours<Ngb> > ngb_list;
		std::vector<Cell>  cell_local, cell_list;
		std::vector<Face>  face_list;
		std::vector<Face*> face_active_list;
		int                nface_active;


		//////////////////

		int sample_freq;
		int myproc, nproc;
		unsigned long long local_n, global_n, virtual_n, boundary_n, import_n, active_n;
		unsigned long long nactive_glb;
		bool all_active;

		std::vector<boundary> proc_tiles;
    std::vector<int     > proc_procs;
		bOctree::Tree<NLEAF_GLB, 4> proc_tree;

		boundary global_domain;
		vec3     global_domain_size;
		Distribute distribute_glb;
		Octree::Tree<NLEAF_LOC, 128> local_tree;   // permits up to 128 inserts
		Octree::Tree<NLEAF_LOC, 128> import_tree;  // permits up to 128 inserts
		bool distribute_data_flag;

		std::vector<int       >  sites_added_local;   // local sites that are in DT
		std::vector<ImportSite>  import_site_list;    // list of remote sites
		int nimport_site_list;
		std::vector<bool      > site_ngb_used;
#if 0
		std::vector<int       > active_sites;
		std::vector<Site      > active_site_list;
#endif
		int  nactive_box;


		bool mpi_debug_flag;


		int iteration;
//		bool isoeos_flag;
		real t_global, dt_global, t_prev;
		real gamma_gas, courant_no;

#if 0
		std::vector<Fluid    > U, dU;
		std::vector<Fluid_st > Wst;
		std::vector<Fluid_rec> Wrec;
		std::vector<real>      gradPsi;
		std::vector<float>     divBi;
#endif

		int  n_restart;
		int  di_log;
		real dt_dump, dt_restart;
		real t_end;

//		std::vector<int> active_faces, active_faces_with_ngb;
		std::vector<int> active_ptcl;
		Scheduler scheduler;

		enum PROFILING_TYPE {
			ITERATE = 0,
			ACTIVEP,
			ACTIVEF,
			PREDICT,
			DISTRIBUTE,
			COLLECT,
			DISTW,
			BATCH,
			MESH,
			FLUID,
			UPDATE,
			GRADIENTS,
			RECONSTRCT,
			PVEL,
			TIMESTEP,
			MISC,
			MPICOM,
			NPROFILING_TYPE};
		double tprofile[NPROFILING_TYPE + 5];

		void determine_sampling_freq();
		void collect_sample_data(std::vector<vec3>&, const std::vector<vec3>&);
		void distribute_data(const bool, const bool, const bool);
		void distribute_work(const int nbatch = 1);
    void distribute_particle_fluid_data();
		void collect_results();
    void distribute_fluid_update_data();
		void collect_reconstruction();
		void get_site_list_batch(
				std::vector<             int  > &active_sites,
				std::vector< std::vector<int> > &active_list,
				const int batch = 0);
		void compute_proc_domain(const std::vector<boundary>&, const std::vector<int>&);

		void sort_local_data();		
		void exchange_ptcl();
		void exchange_primitive();
		void exchange_Wstate();
		void exchange_reconstruction();
		void compute_fluid_grad(const bool);
		void exchange_pvel();

		//
		// >>>>  iterate.cpp
		//
		void dump_binary(const char*);
		void read_binary(const char*, const int);

		//
		// >>>>  iterate.cpp
		//
		void get_active_ptcl(const bool);
		void get_active_faces();
		void fluid_update(const bool);
		void iterate_step();
		void iterate();
		void dump_profile_info();
		const Energy get_energy() const;


		//////////////////    GEOMETRY CALLS
		//

		void clear_mesh(const bool);
		void build_mesh_active(
				std::vector<             int  > &active_sites,
				std::vector< std::vector<int> > &active_list,
				std::vector<             int  > &active_sites_ngb);

		void extract_ngb_from_mesh();
		int  build_mesh_global();
		int  build_mesh(std::vector<Site>&, int nbox);
		int  build_voronoi_cells(
				std::vector<Site> &site_list, 
				int nbox);
		int build_mesh(
				std::vector< std::vector<int> > site_list,
				const Octree::Tree<NLEAF_LOC, 128> &tree,
				const bool store_ngb);

		void get_sites2add(
				const std::vector<LightSite>&,
				const TBOUNDARY&,
				const TVEC3&,
				const TVEC3&,
				const int,
				std::vector<ImportSite>&);
		void relax_mesh(const int);
		bool get_cell_vtx_list(const Site &si, std::vector<TVEC3> &vtx_list);
		bool remove_site(const int, std::vector< std::pair<int, std::pair<real, real> > >&);
		bool insert_site(const vec3&, const int, std::vector< std::pair<int, std::pair<real, real> > >&);


		/////////////////// PROBLEM CALLS
		//
		void set_geometry(const bool);
		void set_problem (const bool);
    bool problem_force_distribute();
		real compute_pressure(const Fluid&) const;
		real compute_entropy_from_ethm(const Fluid&) const;
		real compute_ethm_from_entropy(const Fluid&) const;
		real compute_ethm_update(const Fluid&, const int i) const;
		bool compute_update_prob(Fluid &W, int i);
    void problem_set_boundary(const int, const real);
//		real compute_pressure_gradient(const Fluid&, const Fluid&) const;
		void predict(const int i);
    const real extra_timestep(const int i);
    void correct(const int i);
		bool compute_pvel1();

		bool refine(const int);
		bool derefine(const int);

		void do_refinement();
		void do_derefinement();



		///////////////////   FLUID CALLS
		//	
		void compute_update();
		void compute_update0();
    void compute_update_first_order(const std::vector<int>&);
		inline void compute_flux(
				const vec3&,
				const Fluid&, 
				const Fluid&, 
				const int, 
				const vec3&, 
				real&, 
				real&, 
				Fluid&);

		inline void riemann_solver(
				Fluid &flux,
				const real Bn_ij,
				const real wn_ij,
				const real dens_L, const real pres_L, const real ethm_L,
				const real velx_L, const real vely_L, const real velz_L, const real By_L, const real Bz_L,
				const real dens_R, const real pres_R, const real ethm_R,
				const real velx_R, const real vely_R, const real velz_R, const real By_R, const real Bz_R
				) const;

		//////////////////
		//
		void compute_pvel();
		void compute_reconstruction();
    void slope_limiter();
		void compute_timesteps(const bool shared = false);
		real compute_timestep (const int);
		void adjust_timesteps();

		system(const int _myproc, const int _nproc) : myproc(_myproc), nproc(_nproc)
		{
			iteration = 0;
			courant_no = 0.4;
			t_global = dt_global = t_prev = 0.0;
			mpi_debug_flag = false;
			local_n = global_n =  virtual_n = 0;
			boundary_n = 0;
			gamma_gas = 5.0/3.0;

			n_restart  = 2;
			di_log     = 1;
			dt_dump    = HUGE;
			dt_restart = HUGE;
			t_end      = 1.0;
		}
		~system() {myMPI::free_type();}

		vec3 periodic(const vec3 &pos)
		{
			real x = pos.x;
			real y = pos.y;
			real z = pos.z;
			while (x <  global_domain.get_rmin().x) {x += global_domain_size.x;}
			while (x >= global_domain.get_rmax().x) {x -= global_domain_size.x;}
			while (y <  global_domain.get_rmin().y) {y += global_domain_size.y;}
			while (y >= global_domain.get_rmax().y) {y -= global_domain_size.y;}
			while (z <  global_domain.get_rmin().z) {z += global_domain_size.z;}
			while (z >= global_domain.get_rmax().z) {z -= global_domain_size.z;}
      assert(x >= global_domain.get_rmin().x);
      assert(y >= global_domain.get_rmin().y);
      assert(z >= global_domain.get_rmin().z);
      assert(x <  global_domain.get_rmax().x);
      assert(y <  global_domain.get_rmax().y);
      assert(z <  global_domain.get_rmax().z);
			return vec3(x,y,z);
		}

	};	


}
#endif // __FVMHD3D3D_H__
