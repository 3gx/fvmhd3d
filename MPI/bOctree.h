#ifndef __BOCTREE_H__
#define __BOCTREE_H__

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

namespace bOctree {

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


	struct Body {
		rBoundary pp;
		Body     *next;
		vec3     xcache;
		int      id;

		Body() {};
		Body(rBoundary &p, const int i) : pp(p)
		{
			next = NULL;
			xcache = p.centre();
			id   = i;
		}
		~Body() {};

		void update(const rBoundary &p)
		{
			pp = p;
			xcache = p.centre();
		}

	};

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
		}

#if 0
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
				const int ic = octant_index(centre, p.xcache);
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
					const int ic = octant_index(current->centre, p.xcache);
					current = &current->child[ic];
					continue;
				}
			}
		}

		///////
		///////
		///////
		
		rBoundary& inner_boundary() 
		{
			inner = rBoundary();

			//			if (!touch   ) return inner;
			if (isempty()) return inner;

			if (isleaf()) 
			{
				for (Body *bp = pfirst; bp != NULL; bp = bp->next)
						inner.merge(bp->pp);
			}
			else 
			{
				for (int ic = 0; ic < 8; ic++) 
					inner.merge(child[ic].inner_boundary());
			}
			return inner;
		}

		template<typename T>	
			void walk_boundary(const rBoundary &bnd, std::vector<int> &site_list, 
					const vector3<T> &global_domain_size) const 
			{
				if (isempty()) return;
				if (!bnd.overlap(inner)) 
					if (!bnd.overlap(inner, global_domain_size)) return;

				if (isleaf()) 
				{
					for (Body *bp = pfirst; bp != NULL; bp = bp->next) 
							if (bnd.overlap(bp->pp, global_domain_size))
								site_list.push_back(bp->id);
				} 
				else 
				{
					for (int ic = 0; ic < 8; ic++)
						child[ic].walk_boundary(bnd, site_list, global_domain_size);
				}
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

	};

	template<int NLEAF, int NIMPMAX> struct Tree {
		int n_bodies, n_nodes;

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
			return leaf_list.size();
		}


		const Node<NLEAF>* leaf(const int i) const {return leaf_list[i];}

		void set_domain(const rBoundary &domain)  {
			root.assign_root(domain);
		}

		void insert(
				rBoundary *p, 
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
