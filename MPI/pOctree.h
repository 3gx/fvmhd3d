#ifndef _POCTREEBND_H_
#define _POCTREEBND_H_

#include "pboundary.h"
#include "memory_pool.h"

namespace pOctree {

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
	struct Body {
		pBoundary *pp;
		Body   *next;
		pfloat3 xcache;
		int     id;

		Body() {};
		~Body() {};
		Body(pBoundary &p, const int i) 
		{
			pp      = &p;
			next    = NULL;
			xcache  = p.pcentre;
			id      = i;
		};


	};


	////////

	template<int NLEAF>
		struct Node {
			float4 centre;     // x,y,z, half-length      // 4
			pBoundary inner;                              // 6
			int nparticle;                                // 7
			int touch;                                    // 8

			union {
				Node *child;                                
				Body *pfirst;                                    
			};                                            // 10

			bool isleaf() const {return nparticle <= NLEAF;}
			bool isempty() const {return nparticle == 0;}
			bool istouched() const {return touch;}

			Node() {clear();}
			void clear()   {nparticle = 0; touch = 0; child = NULL;}
			Node(const Node &parent, const int ic) 
			{
				centre.w = parent.centre.w*0.5f;
				centre.x = parent.centre.x + centre.w * ((ic & 1) ? 1.0f : -1.0f);
				centre.y = parent.centre.y + centre.w * ((ic & 2) ? 1.0f : -1.0f);
				centre.z = parent.centre.z + centre.w * ((ic & 4) ? 1.0f : -1.0f);

				child     = NULL;
				nparticle = 0;
				touch     = 0;
			};
			~Node() {};

			////

			void assign_root(const pBoundary &bnd) 
			{
				const dvec3 c(bnd.pcentre.to_dvec3());
				const dvec3 s(bnd.phsize .to_dvec3());

				centre.x = c.x;
				centre.y = c.y;
				centre.z = c.z;
				centre.w = fmax(s.x, fmax(s.y,s.z));

				float hsize = 1.0f;
				while(hsize > centre.w) hsize *= 0.5f;
				while(hsize < centre.w) hsize *= 2.0f;

				centre.w = hsize;
				child = NULL;

			};

			void push_body(Body &p) 
			{
				p.next = pfirst;
				pfirst = &p;
			}

			///////

#if 0
			static inline int octant_index(const pfloat3 &lhs, const pfloat3 &rhs) {
				return (
						((lhs.x.leq(rhs.x)) ? 1 : 0) +  
						((lhs.y.leq(rhs.y)) ? 2 : 0) + 
						((lhs.z.leq(rhs.z)) ? 4 : 0));
			}
#else
			static inline int octant_index(float4 &center, float4 &pos){
				int i = __builtin_ia32_movmskps(
						(float4::v4sf)__builtin_ia32_cmpltps(center.v, pos.v));
				return 7 & i;
			}
#endif

			//////

			void divide_treenode(memory_pool<Node> &pool) 
			{
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

#if 0      
					pfloat3 c;
					c.x.set(centre.x);
					c.y.set(centre.y);
					c.z.set(centre.z);
					int ic = octant_index(c, p.xcache);
#else
					float4 x4;
					x4.x = p.xcache.x.to_double();
					x4.y = p.xcache.y.to_double();
					x4.z = p.xcache.z.to_double();
					x4.w = 0.0f;
					int ic = octant_index(centre, x4);
#endif

					child[ic].push_body(p);
					child[ic].nparticle++;
					child[ic].touch = 1;
				}
				for (int ic = 0; ic < 8; ic++) {
					if (!child[ic].isleaf()) {
						child[ic].divide_treenode(pool);
					}
				}
			}

			//////


			static void insert_body(Node &node, Body &p, memory_pool<Node> &pool) 
			{
				Node *current = &node;
				while(true) {

					if (current->isleaf()) {
						for (Body *bp = current->pfirst; bp != NULL; bp = bp->next) {
							assert( !(
										bp->xcache.x.uval == p.xcache.x.uval &&
										bp->xcache.y.uval == p.xcache.y.uval &&
										bp->xcache.z.uval == p.xcache.z.uval) );
						}
						current->push_body(p);
						current->nparticle++;
						current->touch = 1;
						if (!current->isleaf()) 
							current->divide_treenode(pool);
						return;
					} else {
						current->nparticle++;
						current->touch = 1;

#if 0	
						pfloat3 c;
						c.x.set(current->centre.x);
						c.y.set(current->centre.y);
						c.z.set(current->centre.z);
						int ic = octant_index(c, p.xcache);
#else
						float4 x4;
						x4.x = p.xcache.x.to_double();
						x4.y = p.xcache.y.to_double();
						x4.z = p.xcache.z.to_double();
						x4.w = 0.0f;
						int ic = octant_index(current->centre, x4);
#endif

						current = &current->child[ic];
						continue;
					}

				}

			}

			////

			void sanity_check() {
				if (isleaf()) {

					pBoundary bnd = pBoundary(centre);
					for (Body *bp = pfirst; bp != NULL; bp = bp->next)  {
						if (!bnd.isinbox(bp->pp->pcentre)) {
							pfloat3 pos = bp->pp->pcentre;
							fprintf(stderr, "pos= %20.16lg %20.16lg %20.16lg\n",
									pos.x.to_double(),
									pos.y.to_double(),
									pos.z.to_double());
							fprintf(stderr, "bnd= "); bnd.dump(stderr, true);
							assert(bnd.isinbox(bp->pp->pcentre));
						}
					}

				} else {

					for (int ic = 0; ic < 8; ic++)
						child[ic].sanity_check();

				}

			}


			////////

			pBoundary& inner_boundary() 
			{
				inner = pBoundary();

//				if (!touch   ) return inner;
				if (isempty()) return inner;

				if (isleaf()) {

					for (Body *bp = pfirst; bp != NULL; bp = bp->next) 
						inner.merge(*bp->pp);

				} else {

					for (int ic = 0; ic < 8; ic++) 
						inner.merge(child[ic].inner_boundary());

				}
				return inner;
			}

			void walk_pBoundary(
					const pBoundary  &bnd, 
					std::vector<int> &export_list) const 
			{
				if (isempty()) return;
				if (!bnd.overlap(inner)) return;

				if (isleaf()) {

					for (Body *bp = pfirst; bp != NULL; bp = bp->next) 
						if (bnd.overlap(*bp->pp))
							export_list.push_back(bp->id);

				} else {

					for (int ic = 0; ic < 8; ic++) 
						child[ic].walk_pBoundary(bnd, export_list);

				}
			}

			void walk_pfloat3(
					const pfloat3           &ppos, 
					std::vector<int>        &export_list) const 
			{
				if (isempty()) return;
				if (!inner.isinbox(ppos)) return;

				if (isleaf()) {

					for (Body *bp = pfirst; bp != NULL; bp = bp->next) 
						if (bp->pp->isinbox(ppos))
							export_list.push_back(bp->id);

				} else {

					for (int ic = 0; ic < 8; ic++) 
						child[ic].walk_pfloat3(ppos, export_list);

				}
			}

			///////////

			void extract_leaves(std::vector<Node*> &leaf_list) {
				if (isempty()) return;

				if (isleaf()) 
					leaf_list.push_back(this);
				else
					for (int ic = 0; ic < 8; ic++)
						child[ic].extract_leaves(leaf_list);
			}

#if 0
			void dump_particles_in_Zorder(std::vector<particle*> &particle_list) {
				if (isempty()) return;

				if (isleaf()) {
					for (poctbody *bp = pfirst; bp != NULL; bp = bp->next) 
						if (bp->islocal())
							particle_list.push_back(bp->pp);
				} else {
					for (int ic = 0; ic < 8; ic++)
						child[ic].dump_particles_in_Zorder(particle_list);

				}
			}
#endif

		};

	template<int NLEAF, int NIMP>
		struct Tree {
			int n_bodies, n_nodes;
			memory_pool< Node<NLEAF> > node_lists[NIMP];
			std::vector< Body        > body_lists[NIMP];
			int nimp;
			Node<NLEAF> root;
			std::vector<Node<NLEAF>*> leaf_list;

			Tree() : nimp(0) {clear();};
			~Tree() {};

			void clear() {
				n_bodies = n_nodes = 0;
				for (int i = 0; i < nimp; i++) 
				{
//					erase_vec(body_lists[i]);
					body_lists[i].clear();
					node_lists[i].reset();
				}
				root.clear();
				nimp = 0;
			}

			void get_leaves() {
				leaf_list.clear();
				leaf_list.reserve(128);
				root.extract_leaves(leaf_list);
			}

			void set_domain(const pBoundary &domain)  {
				root.assign_root(domain);
			}

			void insert(
					pBoundary *p, 
					const int np, 
					const int offs,
					const int nmax_nodes) 
			{
				assert(nimp < NIMP);
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

#if 0
			void insert(particle *p, const int np) {
				const int n_bodies_old = n_bodies;
				for (int i = 0; i < np; i++) {
					assert(n_bodies < n_bodies_max);
					body_list[n_bodies++] = poctbody(p[i]);
				}

				n_bodies = n_bodies_old;
				for (int i = 0; i < np; i++) {
					poctnode<NLEAF>::insert_poctbody(root, body_list[n_bodies++], node_list);
				}
				n_nodes = node_list.used;
				assert(n_nodes < n_nodes_max);

			}
#endif

		};

}

#endif // _OCTREE_H_
