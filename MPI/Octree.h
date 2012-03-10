#ifndef __OCTREE_H__
#define __OCTREE_H__

#include <vector>
#include <cassert>
#include <memory.h>
#include "vector3.h"
#include "Boundary.h"
#include "memory_pool.h"
#include "mytimer.h"

#if 0
#define PREFETCH(x) __builtin_prefetch(x);
#else
#define PREFETCH(x)
#endif

// #define __FP32__

namespace Octree {

#define myassert(x) assert(x);
	enum {VOROSKIPFLAG = (1 << 31)};

	typedef vector3<float>  fvec3;
	typedef vector3<double> dvec3;
	struct float4{
		typedef float  v4sf __attribute__ ((vector_size(16)));
		typedef double v2df __attribute__ ((vector_size(16)));
		static v4sf v4sf_abs(v4sf x){
			typedef int v4si __attribute__ ((vector_size(16)));
			v4si mask = {0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff};
			return __builtin_ia32_andps(x, (v4sf)mask);
		}
		union{
			v4sf v;
			struct{
				float x, y, z, w;
			};
		};
		float4() : v((v4sf){0.f, 0.f, 0.f, 0.f}) {}
		float4(float x, float y, float z, float w) : v((v4sf){x, y, z, w}) {}
		float4(const fvec3 &f) : v((v4sf){f.x, f.y, f.z, 0}) {}
		float4(const dvec3 &f) : v((v4sf){f.x, f.y, f.z, 0}) {}
		float4(v4sf _v) : v(_v) {}
		float4 abs(){
			typedef int v4si __attribute__ ((vector_size(16)));
			v4si mask = {0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff};
			return float4(__builtin_ia32_andps(v, (v4sf)mask));
		}
		void dump(){
			std::cerr << x << " "
				<< y << " "
				<< z << " "
				<< w << std::endl;
		}
		const fvec3 xyz() {return fvec3(x,y,z);}
		v4sf operator=(const float4 &rhs){
			v = rhs.v;
			return v;
		}
		float4(const float4 &rhs){
			v = rhs.v;
		}
	};
#ifndef __FP32__
	typedef double          real;
	typedef vector3<double> vec3;
	typedef long long       rint;
#else
	typedef float           real;
	typedef vector3<float>  vec3;
	typedef int             rint;
#endif

	typedef Boundary<vec3> rBoundary;


	struct Particle {
		vec3 pos;
		int  idx;
		Particle() {}
		Particle(const vec3 &_pos, const int _idx) : pos(_pos), idx(_idx) {}
	};

	struct Body {
		Body     *next;
		vec3     pos;
		int      id;

		Body() {};
		Body(Particle &p, const int i) {
			next = NULL;
			pos  = vec3(p.pos.x, p.pos.y, p.pos.z);
			id   = i;
		}
		~Body() {};

		void update(const vec3 _pos, const int _id)
		{
			pos = _pos;
			id  = _id;
		}

	};

	template<typename T>
		inline const T sign(const T x) {return (x < 0) ? -1 : +1;}

	template<bool PERIODIC>	
		inline const vec3 periodic(const vec3 &p, const vec3 &Lbox)
		{
			if (PERIODIC)
				return vec3(
						p.x - ((std::abs(p.x) < 0.5*Lbox.x) ? 0.0 : sign(p.x) * Lbox.x),
						p.y - ((std::abs(p.y) < 0.5*Lbox.y) ? 0.0 : sign(p.y) * Lbox.y),
						p.z - ((std::abs(p.z) < 0.5*Lbox.z) ? 0.0 : sign(p.z) * Lbox.z)
						);
			else
				return p;
		}


	template<int NLEAF>	struct Node {
		fvec3 centre;
		float hsize;
		rBoundary inner; //, outer;
		int nparticle;
		bool touch;

		union {
			Node *child;
			Body *pfirst;
		};
		Node *parent;

		const bool isleaf()    const {return nparticle <= NLEAF; }
		const bool isempty()   const {return nparticle == 0; }
		const bool istouched() const {return touch;}

#if 1
		Node() : nparticle(0), touch(false), child(NULL), parent(NULL) {}
#else
		Node() {}
#endif

#if 1
		void clear() {
			nparticle = 0;
			touch = false;
			child = NULL;
			parent = NULL;
		}
#endif

		Node(Node &parent, const int ic) {
			hsize = parent.hsize * 0.5f;
			centre.x = parent.centre.x + hsize * ((ic & 1) ? 1.0f : -1.0f);
			centre.y = parent.centre.y + hsize * ((ic & 2) ? 1.0f : -1.0f);
			centre.z = parent.centre.z + hsize * ((ic & 4) ? 1.0f : -1.0f);

			nparticle    = 0;
			touch        = false;
			child        = NULL;
			this->parent = &parent;
		}
		~Node() {
			//			clear();
		};

		void assign_root(const rBoundary &bnd) {
			centre = bnd.centre();
			hsize  = 1.0f;

			const fvec3 hsize3 = bnd.hsize();
			const float hsize0 = std::max(hsize3.x, std::max(hsize3.y, hsize3.z));
			while (hsize > hsize0) hsize *= 0.5f;
			while (hsize < hsize0) hsize *= 2.0f;

			child  = NULL;
			parent = NULL;

		}

		void push_body(Body &p) {
			p.next = pfirst;
			pfirst = &p;
			nparticle++;
			touch = true;

#if 0
      const rBoundary bnd = rBoundary(centre-1.0001*hsize, centre+1.0001*hsize);
      if (!(bnd.isinbox(p.pos)))
      {
        fprintf(stderr, " pos= %g %g %g \n",
            p.pos.x,
            p.pos.y,
            p.pos.z);
        bnd.dump(stderr, true);
      }
      assert(bnd.isinbox(p.pos));
#endif
    }

#if 1
    static inline int octant_index(const fvec3 &lhs, const fvec3 &rhs) {
      return
        (((lhs.x <= rhs.x) ? 1 : 0) +  
         ((lhs.y <= rhs.y) ? 2 : 0) + 
         ((lhs.z <= rhs.z) ? 4 : 0));
    }
#else
    static inline int octant_index(const float4 &center, const float4 &pos){
      int i = __builtin_ia32_movmskps(
          (float4::v4sf)__builtin_ia32_cmpltps(center.v, pos.v));
      return 7 & i;
    }
#endif


    void divide_node(memory_pool<Node> &pool) {
      static Body *bodyarray[NLEAF + 1];
      assert(nparticle <= NLEAF + 1);
      Body *p = pfirst;

      for (int i = 0; i < nparticle; i++) {
        bodyarray[i] = p;
        p = p->next;
      }
      pfirst = NULL;

      assert(child == NULL);
      child = pool.get(8);

      for (int ic = 0; ic < 8; ic++) {
        child[ic] = Node(*this, ic);
      }

      for (int i = 0; i < nparticle; i++) {
        Body &p = *bodyarray[i];
        const int ic = octant_index(centre, p.pos);
        child[ic].push_body(p);
      }

      for (int ic = 0; ic < 8; ic++) {
        if (!child[ic].isleaf()) {
          child[ic].divide_node(pool);
        }
      }

    }

    static void insert_body(Node &node, Body &p, memory_pool<Node> &pool) {
      Node *current = &node;
      while (true) {
        if (current->isleaf()) {
          current->push_body(p);
          if (!current->isleaf())
            current->divide_node(pool);
          return;
        } else {
          current->nparticle++;
          current->touch = true;
          const int ic = octant_index(current->centre, p.pos);
          current = &current->child[ic];
          continue;
        }
      }
    }

    ///////
    ///////
    ///////

    void sanity_check() const 
    {
//      const rBoundary bnd = rBoundary(centre, hsize);
      const rBoundary bnd = rBoundary(centre-1.0001*hsize, centre+1.0001*hsize);
      //      assert(bnd.isinbox(inner));
      if (isleaf()) 
      {
        for (Body *bp = pfirst; bp != NULL; bp = bp->next) 
          if (bp->id >= 0)
          {
            if (!(bnd.isinbox(bp->pos)))
            {
              fprintf(stderr, " pos= %g %g %g \n",
                  bp->pos.x,
                  bp->pos.y,
                  bp->pos.z);
              bnd.dump(stderr, true);
            }
            assert(bnd.isinbox(bp->pos));
          }
      } 
      else 
      {
        for (int ic = 0; ic < 8; ic++)
          child[ic].sanity_check();
      }
    }

    rBoundary& inner_boundary() 
    {
      inner = rBoundary();

      //			if (!touch   ) return inner;
      if (isempty()) return inner;

      if (isleaf()) 
      {
        for (Body *bp = pfirst; bp != NULL; bp = bp->next)
          if (bp->id >= 0)
            inner.merge(rBoundary(bp->pos));
      }
      else 
      {
        for (int ic = 0; ic < 8; ic++) 
          inner.merge(child[ic].inner_boundary());
      }
      return inner;
    }

    template<class Site, typename T>	
      void walk_boundary(const rBoundary &bnd, std::vector<Site> &site_list, 
          const vector3<T> &global_domain_size) const 
      {
        if (isempty()) return;
        if (!bnd.overlap(inner)) 
          if (!bnd.overlap(inner, global_domain_size)) return;

        if (isleaf()) 
        {
          for (Body *bp = pfirst; bp != NULL; bp = bp->next) 
            if (bp->id >= 0)
              if (bnd.overlap(bp->pos, global_domain_size))
                site_list.push_back(Site(bp->pos, bp->id));
        } 
        else 
        {
          for (int ic = 0; ic < 8; ic++)
            child[ic].walk_boundary(bnd, site_list, global_domain_size);
        }
      }
    
    template<class Site, typename T>	
      bool walk_boundary_flag(const rBoundary &bnd, std::vector<Site> &site_list, 
          const vector3<T> &global_domain_size, bool flag = false) const 
      {
        if (isempty()) return flag;
        if (!bnd.overlap(inner)) 
          if (!bnd.overlap(inner, global_domain_size))
            return flag;

        if (isleaf()) 
        {
          for (Body *bp = pfirst; bp != NULL; bp = bp->next) 
            if (bnd.overlap(bp->pos, global_domain_size, flag))
              site_list.push_back(Site(bp->pos, bp->id));
        } 
        else 
          for (int ic = 0; ic < 8; ic++)
            flag |= child[ic].walk_boundary_flag(bnd, site_list, global_domain_size, flag);
        return flag;
      }


    void extract_leaves(std::vector<Node*> &leaf_list) 
    {
      if (isempty()) return;

      if (isleaf()) 
        leaf_list.push_back(this);
      else
        for (int ic = 0; ic < 8; ic++)
          child[ic].extract_leaves(leaf_list);
    }

    void get_rmax(std::vector< std::pair<int, real> > &rmax_list) 
    {
      if (isempty()) return;

      if (isleaf()) 
      {
        for (Body *bp = pfirst; bp != NULL; bp = bp->next) 
          if (bp->id >= 0)
            rmax_list.push_back(std::make_pair(bp->id, hsize));
      }
      else 
        for (int ic = 0; ic < 8; ic++)
          child[ic].get_rmax(rmax_list);
    }

    void dump_particles_in_Zorder(std::vector<int> &particle_list) 
    {
      if (isempty()) return;

      if (isleaf()) 
        for (Body *bp = pfirst; bp != NULL; bp = bp->next) 
        {
          assert(bp->id >= 0);
          particle_list.push_back(bp->id);
        }
      else
        for (int ic = 0; ic < 8; ic++)
          child[ic].dump_particles_in_Zorder(particle_list);
    }


  };

  template<int NLEAF, int NIMPMAX> struct Tree {
    int n_bodies, n_nodes, n_leaves;

    memory_pool< Node<NLEAF> > node_lists[NIMPMAX];
    std::vector< Body >        body_lists[NIMPMAX];
    int nimp;

    Node<NLEAF> root;
    std::vector< Node<NLEAF>* > leaf_list;

    Tree() : n_bodies(0), n_nodes(0), nimp(0) {}
    ~Tree() {}

    void clear() {
      n_bodies = n_nodes = 0;
      for (int i = 0; i < nimp; i++) 
      {
        body_lists[i].clear();
        node_lists[i].reset();
      }
      root.clear();
      nimp = 0;
    }

    void update_bodies()
    {
      for (int i = 0; i < nimp; i++)
        for (int j = 0; j < (const int)body_lists[i].size(); j++)
          body_lists[i][j].update();
    }

    const size_t get_leaves() {
      leaf_list.clear();
      leaf_list.reserve(128);
      root.extract_leaves(leaf_list);
      n_leaves = leaf_list.size();
      return leaf_list.size();
    }


    const Node<NLEAF>* leaf(const int i) const {return leaf_list[i];}

    void set_domain(const rBoundary &domain)  {
      root.assign_root(domain);
    }

    void insert(
        Particle *p, 
        const int np, 
        const int offs,
        const int nmax_nodes) 
    {
      assert(nimp < NIMPMAX);
      memory_pool< Node<NLEAF> > &node_list = node_lists[nimp];
      std::vector< Body        > &body_list = body_lists[nimp];

      body_list.resize(np);
      node_list.resize(nmax_nodes);

      for (int i = 0; i < np; i++) {
        body_list[i] = Body(p[i], offs + i);
      }

      n_bodies += np;

      for (int i = 0; i < np; i++) {
        Node<NLEAF>::insert_body(root, body_list[i], node_list);
      }
      n_nodes += node_list.used;

      assert(node_list.free_unused());

      nimp++;
    }


  };

}
#endif // __OCTREE_H__
