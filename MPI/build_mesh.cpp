#include "fvmhd3d.h"
#define SMALLDIFF1 1.0e-10

#if 0
#define _DEBUG_PRINT_
#endif

#ifdef _DEBUG_PRINT_
#define MY_MPIBARRIER() MPI_Barrier(MPI_COMM_WORLD)
#define MY_TIMER() mytimer::get_wtime()
#else
#define MY_MPIBARRIER()
#define MY_TIMER() 0.0
#endif

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h> 

#if 0
typedef CGAL::Exact_predicates_exact_constructions_kernel K; 
#define _EXACT_CONSTRUCTIONS_
#else
typedef CGAL::Exact_predicates_inexact_constructions_kernel K; 
#endif

typedef CGAL::Triangulation_vertex_base_with_info_3<int, K> Vb; 
typedef CGAL::Triangulation_data_structure_3<Vb>            Tds; 
typedef CGAL::Delaunay_triangulation_3<K, Tds> DT;

typedef DT::Vertex_iterator TVertex_iterator;
typedef DT::Vertex_handle   TVertex_handle;
typedef DT::Point           TPoint;
typedef DT::Edge            TEdge;
typedef DT::Cell_handle     TCell_handle;
typedef DT::Cell            TCell;
typedef DT::Cell_circulator TCell_circulator;
typedef DT::Facet_circulator TFacet_circulator;
			
static int nsend0[NMAXPROC], nrecv0[NMAXPROC];

template<class Triangulation> 
struct Traits_for_spatial_sort:public Triangulation::Geom_traits
{ 
	typedef typename Triangulation::Geom_traits Gt; 
	typedef std::pair<const typename Triangulation::Point*,int> Point_3; 

	struct Less_x_3
	{ 
		bool operator()(const Point_3& p,const Point_3& q) const 
		{ 
			return typename Gt::Less_x_3()(*(p.first),*(q.first)); 
		} 
	}; 

	struct Less_y_3
	{ 
		bool operator()(const Point_3& p,const Point_3& q) const 
		{ 
			return typename Gt::Less_y_3()(*(p.first),*(q.first)); 
		} 
	}; 

	struct Less_z_3
	{ 
		bool operator()(const Point_3& p,const Point_3& q) const 
		{ 
			return typename Gt::Less_z_3()(*(p.first),*(q.first)); 
		} 
	};   

	Less_x_3 less_x_3_object () const {return Less_x_3();} 
	Less_y_3 less_y_3_object () const {return Less_y_3();} 
	Less_z_3 less_z_3_object () const {return Less_z_3();} 
}; 

namespace fvmhd3d
{
	struct Edge
	{
		int s1, s2;
		Edge(const int _s1, const int _s2) : s1(_s1 < _s2 ? _s1 : _s2), s2(_s1 < _s2 ? _s2 : _s1) {}
		Edge() {};
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
	};

	struct int0
	{
		int oct;
		int0() : oct(0) {}
	};


	DT                                  *T_ptr;                   // Delaunay triagnulation
	TVertex_handle                       hint;                // Vertex_handle
  std::vector<TVertex_handle>          sites_in_DT;         // vertex_handles of local sites
	std::map< std::pair<int, int>, int0> sites_added_remote;  // local sites that are exported
//	static std::vector<int>                     sites_added_local;   // local sites that are in DT
//	static std::vector<ImportSite>              import_site_list;    // list of remote sites
	int nsites_in_DT;                                         // number of sites in DT
	std::map< std::pair<int, int>, int0> return_site_list_map;
//	static std::vector<int>                     sites_added_local_list;
	std::vector<int       > active_sites;

//	static std::vector<bool> site_ngb_used;
//	static std::vector<bool> active_sites;


	// Use Lloyd relaxation method
	void system::relax_mesh(const int niter)
	{
		if (myproc == 0)
			fprintf(stderr, " --- relaxing mesh --- \n");
		
		clear_mesh(false);

		for (int iter = 0; iter < niter; iter++)
		{
			const double t00 = mytimer::get_wtime();
			distribute_data(false, false, false);

			const double t10 = mytimer::get_wtime();
			int nattempt = build_mesh_global();
			clear_mesh(false);

			const double t20 = mytimer::get_wtime();
			for (int i = 0; i < (int)local_n; i++)
				ptcl_local[i].set_pos(periodic(cell_local[i].centroid));

			const double t30 = mytimer::get_wtime();

			double dt10 = t10 - t00;
			double dt20 = t20 - t10;
			double dt30 = t30 - t00;

			double dt10max, dt20max, dt30max;
			MPI_Allreduce(&dt10, &dt10max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
			MPI_Allreduce(&dt20, &dt20max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
			MPI_Allreduce(&dt30, &dt30max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

			if (myproc == 0)
				fprintf(stderr, "relaxing mesh: iteration %d out of %d  [ %g / %g sec :: %g / %g cells/s/proc/thread ]\n",
						iter, niter, 
						dt20max, dt30max,
						global_n/nproc/dt20max,
						global_n/nproc/dt30max);

			double volume_loc = 0.0;
			{
				std::vector<TREAL> v(local_n);
				for (int i = 0; i < (int)local_n; i++)
					v[i] = cell_local[i].Volume;
				std::sort(v.begin(), v.end());  // sort volumes from low to high, to avoid roundoff errors
				for (int i = 0; i < (int)local_n; i++)
					volume_loc += v[i];
			}


			double volume_glob = 0.0;	
			int    nattempt_max, nattempt_min;
			MPI_Allreduce(&volume_loc, &volume_glob,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&nattempt,   &nattempt_max, 1, MPI_INT,    MPI_MAX, MPI_COMM_WORLD);
			MPI_Allreduce(&nattempt,   &nattempt_min, 1, MPI_INT,    MPI_MIN, MPI_COMM_WORLD);

			const double volume_exact = global_domain_size.x*global_domain_size.y*global_domain_size.z;
			if (myproc == 0)
			{
				fprintf(stderr, "   computed_volume= %g  exact_volume= %g diff= %g [ %g ]  nattempt= %d %d \n",
						volume_glob, volume_exact, 
						volume_glob - volume_exact,	(volume_glob - volume_exact)/volume_exact,
						nattempt_min, nattempt_max);
			}
		}
	}

	void system::clear_mesh(const bool flush_mem)
	{
		hint = DT::Vertex_handle();

		sites_added_local.resize(local_n);
		sites_added_remote.clear();

		nsites_in_DT = 0;		

//		if (flush_mem)
		{
			clear_vec(sites_added_local);
			clear_vec(active_sites);
			clear_vec(sites_in_DT);
//			clear_vec(cell_local);
			clear_vec(import_site_list);
			clear_vec(site_ngb_used);
			cell_local.clear();
		}
	}

	void system::extract_ngb_from_mesh()
  {
    std::vector< Neighbours<Ngb> > ngb_list;

		ngb_list.resize(local_n);
		for (int i = 0; i < (int)local_n; i++)
		{
			const Cell &ci = cell_local[i];
			const int nj = ci.ngb.size(); 
			ngb_list[i].clear();
			for (int jp = 0; jp < nj; jp++)   // loop over neighbours
			{
				const int jDT = ci.ngb[jp];
				const ImportSite &s = import_site_list[jDT];
				ngb_list[i].push_back(Ngb(s.proc, s.id));
			}
		}
    this->ngb_list.swap(ngb_list);
	}

	int system::build_mesh_global()
	{
    DT T;
    T_ptr = &T;
		sites_added_local.resize(local_n);
		sites_in_DT.resize(local_n);

    clear_vec(cell_local);
		cell_local.resize(local_n);

		active_sites.resize(local_n);
		for (int i = 0; i < (int)local_n; i++)
			active_sites[i] = 0;

		import_site_list.resize(2*local_n);
		site_ngb_used.resize(2*local_n);
		for (int i = 0; i < 2*(int)local_n; i++)
			site_ngb_used[i] = false;
		nimport_site_list = 0;

		for (int i = 0; i < (int)local_n; i++)
			sites_added_local[i] = 0;

		nsites_in_DT = 0;
		hint = DT::Vertex_handle();

		const int nleaves = local_tree.n_leaves;
		std::vector<Site> active_site_list(nleaves * NLEAF_LOC);

		nactive_box = 0;
		for (int leaf = 0; leaf < nleaves; leaf++)
		{
			int el = nactive_box * NLEAF_LOC;
			for (Octree::Body *bp = local_tree.leaf(leaf)->pfirst; bp != NULL; bp = bp->next)
				if (bp->id >= 0)
				{
					assert(!ptcl_local[bp->id].is_virtual());
					active_site_list[el++] = Site(bp->pos, bp->id, ptcl_local[bp->id].rmax);
				}

			if (el > nactive_box * NLEAF_LOC) 
			{
				active_site_list[el - 1].id = int_map(active_site_list[el-1].id); // -1-active_site_list[el-1].id;
				nactive_box++;
			}
		}

		const int nattempt = build_mesh(active_site_list, nactive_box);

		return nattempt;
	}

	void system::get_sites2add(
			const std::vector<LightSite> &site_list,
			const TBOUNDARY &outer_box,
			const TVEC3     &outer_box_centre,
			const TVEC3     &outer_box_hsize,
			const int       proc,
			std::vector<ImportSite> &sites2add)
	{
		const int nsite = site_list.size();
		for (int site = 0; site < nsite; site++)
		{
			asm("#get_sites2add-inner_loop-beg");
			const LightSite &s = site_list[site];
			assert(s.id >= 0);

			TVEC3 dr = s.pos - outer_box_centre;

			int oct = 0;
			if      (dr.x >  outer_box_hsize.x) {dr.x -= global_domain_size.x; oct += 1;}
			else if (dr.x < -outer_box_hsize.x) {dr.x += global_domain_size.x; oct += 1;}
			if      (dr.y >  outer_box_hsize.y) {dr.y -= global_domain_size.y; oct += 2;}
			else if (dr.y < -outer_box_hsize.y) {dr.y += global_domain_size.y; oct += 2;}
			if      (dr.z >  outer_box_hsize.z) {dr.z -= global_domain_size.z; oct += 4;}
			else if (dr.z < -outer_box_hsize.z) {dr.z += global_domain_size.z; oct += 4;}

#if 0      // sanity check
			{
				const TVEC3 drt = dr.abseach() - outer_box_hsize;
				if (drt.x > 0.0) fprintf(stderr, " site= %d  drt.x= %g \n", site, drt.x);
				if (drt.y > 0.0) fprintf(stderr, " site= %d  drt.y= %g \n", site, drt.y);
				if (drt.z > 0.0) fprintf(stderr, " site= %d  drt.z= %g \n", site, drt.z);
				assert(drt.x <= 0.0);
				assert(drt.y <= 0.0);
				assert(drt.z <= 0.0);
			}
#endif

			oct = 1 << oct;

			int *value;
			if (proc == myproc) 
			{
				value = &sites_added_local [s.id];
#if 0
				if (*value == -1)
				{
					*value = 0;
					sites_added_local_list.push_back(s.id);
				}
#endif
			}
			else	value = &sites_added_remote[std::make_pair(proc, s.id)].oct;


			if ((*value & oct) == oct) continue;
			*value |= oct;

			sites2add.push_back(ImportSite(proc, s.id, oct, dr + outer_box_centre));
			asm("#get_sites2add-inner_loop-end");
		}
	}

	int system::build_mesh(
			std::vector<Site> &site_list, 
			int nbox)

	{
    DT &T = *T_ptr;
		// extract number of sites in each of the box
		//
		std::vector<TBOUNDARY> local_outer_box(nbox);
		std::vector<int>   nbox_size(nbox);
		for (int i = 0; i < nbox; i++)
		{
			int el = i * NLEAF_LOC;
			int j;
			for (j = 0; j < NLEAF_LOC; j++)
				if (site_list[el + j].id < 0)
				{
					site_list[el+j].id = int_map(site_list[el+j].id); 
					break;
				}
			assert(j < NLEAF_LOC);
			nbox_size[i] = j + 1;
		}

		for (int ibox = 0; ibox < nbox; ibox++)
			for (int isite = ibox*NLEAF_LOC; isite < ibox*NLEAF_LOC + nbox_size[ibox]; isite++)
			{
				const int id = site_list[isite].id;
				assert(id >= 0);
				active_sites[id] |= 1;
				cell_local       [id].ngb.clear();
			}

		int  nfailed = 1;     // number of sites that failed to be built in a pass
		int nattempt = 0;     // number of passes required to triangulate the system
		int nngb_tot = 0;     // total number of neighbours of all sites
    int nngb_min = local_n;
    int nngb_max = 0;

		// main loop
		//
		while (nfailed > 0)
		{
			// compute outer_box 
			//
			for (int ibox = 0; ibox < nbox; ibox++)
			{
				for (int isite = ibox*NLEAF_LOC; isite < ibox*NLEAF_LOC + nbox_size[ibox]; isite++)
					if (site_list[isite].failed)
						local_outer_box[ibox].merge(TBOUNDARY(
									site_list[isite].pos - site_list[isite].rmax,
									site_list[isite].pos + site_list[isite].rmax));
#if 1       // box size must be smaller than the domain size
				const vec3 sz = local_outer_box[ibox].hsize() * 2.0;
        if (sz.x >= global_domain_size.x ||
            sz.y >= global_domain_size.y ||
            sz.z >= global_domain_size.z)
        {
          fprintf(stderr, "myproc= %d :: %g %g %g [ %g %g %g ] nattempt= %d\n",
              myproc, sz.x, sz.y, sz.x,
              global_domain_size.x,
              global_domain_size.y,
              global_domain_size.z,
              nattempt);
        }
				assert(sz.x < global_domain_size.x);
				assert(sz.y < global_domain_size.y);
				assert(sz.z < global_domain_size.z);
#endif
			}

			// add sites that overlap with the outer_box to DT
			// the domain is assumed to be periodic ...
			//
			{
				std::vector<ImportSite> sites2add;
				sites2add.reserve(16);
				{
					std::vector< LightSite> sites_lst;
					sites_lst.reserve(16);

					// first add local sites
					//
					for (int ibox = 0; ibox < nbox; ibox++)
					{
						sites_lst.clear();
						local_tree.root.walk_boundary(local_outer_box[ibox], sites_lst, global_domain_size);
						get_sites2add(
								sites_lst,
								local_outer_box[ibox], 
								local_outer_box[ibox].centre(),
								local_outer_box[ibox].hsize(),
								myproc,
								sites2add);
					}
				}

				// now add remote sites
				//
				{
					std::vector<TBOUNDARY> bnd_send[NMAXPROC];
					std::vector<TBOUNDARY> bnd_recv[NMAXPROC];
					for (int p = 0; p < nproc; p++)
						bnd_send[p].reserve(16);

					std::vector<int > remote_tiles;
					std::vector<bool> exported(nproc);
					remote_tiles.reserve(16);

					for (int ibox = 0; ibox < nbox; ibox++)
					{
						remote_tiles.clear();	
						proc_tree.root.walk_boundary(local_outer_box[ibox], remote_tiles, global_domain_size);
						const int ntile = remote_tiles.size();

						for (int p = 0; p < nproc; p++) exported[p] = false;

						for (int itile = 0; itile < ntile; itile++)
						{
							const int tile = remote_tiles[itile];
							const int proc = proc_procs[tile]; // distribute_glb.tiles2proc[tile];
              assert(proc  >= 0);
              assert(proc < nproc);

							if (proc == myproc || exported[proc]) continue;

							bnd_send[proc].push_back(local_outer_box[ibox]);
							exported[proc] = true;
						}
					}

#if 0
					myMPI::all2all(bnd_send, bnd_recv, myproc, nproc, mpi_debug_flag);
#else
					myMPI::all2all<true>(bnd_send, bnd_recv, myproc, nproc, 1, nsend0, nrecv0);
#endif

					std::vector<double> site_send[NMAXPROC];
					std::vector<double> site_recv[NMAXPROC];
					for (int p = 0; p < nproc; p++)
						site_send[p].reserve(16);

					std::vector<ImportSite> sites2add_tmp;
					std::vector< LightSite> sites_lst;
					sites2add_tmp.reserve(16);
					sites_lst.reserve(16);
					for (int p = 0; p < nproc; p++)
						for (int q = 0; q < (const int)bnd_recv[p].size(); q++)
						{
							const int proc = p;
							assert(proc != myproc);
							const TBOUNDARY &outer_box = bnd_recv[p][q];

							sites_lst.clear();
							sites2add_tmp.clear();
							local_tree.root.walk_boundary(outer_box, sites_lst, global_domain_size);
							get_sites2add(
									sites_lst,
									outer_box, 
									outer_box.centre(),
									outer_box.hsize(),
									proc,
									sites2add_tmp);
							for (int j = 0; j < (const int)sites2add_tmp.size(); j++)
							{
								const TVEC3              &dr = sites2add_tmp[j].pos;
								const int jid = sites2add_tmp[j].id;
								const int id[] = 
								{
									((active_sites[jid] & 1) == 1) ? jid : int_map(jid),
									sites2add_tmp[j].oct
								};
								const double *dbl_id = (double*)id;

								site_send[proc].push_back((double)dr.x);
								site_send[proc].push_back((double)dr.y);
								site_send[proc].push_back((double)dr.z);
								site_send[proc].push_back(dbl_id[0]);
							}
						}

#if 0
					myMPI::all2all(site_send, site_recv, myproc, nproc, mpi_debug_flag);
#else
					myMPI::all2all<true>(site_send, site_recv, myproc, nproc, 1, nsend0, nrecv0);
#endif

					for (int p = 0; p < nproc; p++)
						for (int q = 0; q < (const int)site_recv[p].size(); q += 4)
						{
							const TREAL xj = (TREAL)site_recv[p][q + 0];
							const TREAL yj = (TREAL)site_recv[p][q + 1];
							const TREAL zj = (TREAL)site_recv[p][q + 2];
							const int *lid = (int*)&site_recv[p][q + 3];
							const int  id  = lid[0];    // remote id of the particles on proc 'p'
							const int  oct = lid[1];

							sites2add.push_back(ImportSite(p, id, oct, TVEC3(xj, yj, zj)));
						}
				}


				{
					// CGAL optimized insertion into DT
					const int np = sites2add.size();
					std::vector< TPoint >                        pnts(np);
					std::vector< std::pair<const TPoint*, int> > points(np);
					for (int ip = 0; ip < np; ip++)
					{
						pnts  [ip]        = TPoint(sites2add[ip].pos.x, sites2add[ip].pos.y, sites2add[ip].pos.z);
						points[ip].first  = &pnts[ip];
						points[ip].second = ip;
					}
					spatial_sort(points.begin(), points.end(), Traits_for_spatial_sort<DT>());

					for (int ip = 0; ip < np; ip++)
					{
						const int  i = points[ip].second;

						int id = sites2add[i].id;
						if (sites2add[i].oct != 1 || sites2add[i].proc != myproc)
						{
							id = local_n + nimport_site_list;
							import_site_list[id] = sites2add[i];
							nimport_site_list++;
							if (id + 1 >= (int)import_site_list.size())
							{
								import_site_list.resize(id + 1 + local_n);
								site_ngb_used.resize(id + 1 + local_n);
								for (int i = 0; i < (int)local_n; i++)
									site_ngb_used[id + 1 + i] = false;
							}
						}
						else
							import_site_list[id] = sites2add[i];

						hint = T.insert(*points[ip].first, hint);
						hint->info() = id;

						if (id < (int)local_n) sites_in_DT[id] = hint;

						nsites_in_DT++;
					}
				}
			}

			// now attempt to construct Voronoi cells of input sites
			// 
			{
				nfailed = 0;
				std::vector<TEdge> edges;
				std::vector<int  > site_ngb_list;
				std::vector<TCell_handle> Tcells;
				const int NMAXEDGE = 1024;
				std::vector<TVEC3> vertex_list[NMAXEDGE];
				Edge               edge_list  [NMAXEDGE];

				for (int ibox = 0; ibox < nbox; ibox++)
				{
					Site *sites2process = &site_list[ibox * NLEAF_LOC];
					const int nsite = nbox_size[ibox]; //sites2process.size();
					assert(nsite > 0 && nsite <= NLEAF_LOC);

					for (int isite = 0; isite < nsite; isite++)
					{
						if (!sites2process[isite].failed) continue;

						const int id = sites2process[isite].id;

						const TVertex_handle &vi = sites_in_DT[id];
						assert(vi != DT::Vertex_handle());
						assert(vi->info() == id);

						edges.clear();
#if 0       // Use native CGAL edge extractor
						T.finite_incident_edges(vi, std::back_inserter(edges));
#else
						{  // Optimized edge extractor

							Tcells.clear();
							T.incident_cells(vi, std::back_inserter(Tcells));

							site_ngb_list.clear();

							const int ncells = Tcells.size();
							for (int icell = 0; icell < ncells; icell++)
							{
								const TCell_handle &ci = Tcells[icell];

								const bool is_infinite = T.is_infinite(ci);

								int idx = -1;
								for (int iv = 0; iv < 4; iv++)
									if (ci->vertex(iv) == vi)
										idx = iv;

								int iadd = 0;
								for (int iv = 0; iv < 4; iv++)
								{
									if (iv == idx) continue;

									const TVertex_handle &v = ci->vertex(iv);
									if (is_infinite)
										if (T.is_infinite(v)) continue;

									const int id_orig = v->info();
									const int id = int_abs(id_orig);
									if (site_ngb_used[id]) continue;

									iadd++;
									site_ngb_used[id] = true;
									site_ngb_list.push_back(id_orig);
									edges.push_back(TEdge(ci, idx, iv));
								}
								assert(iadd < 4);
							}

							const int nngb = site_ngb_list.size();

							for (int j = 0; j < nngb; j++)
								site_ngb_used[int_abs(site_ngb_list[j])] = false;
						}
#endif

						bool failed_flag = false;
						TREAL r2max = 0.0;
						int  ninf  = 0;

						assert((int)edges.size() <= NMAXEDGE);
						int nedge = 0;

						for (std::vector<TEdge>::iterator edge_it = edges.begin(); edge_it != edges.end(); edge_it++)
						{
							const TCell_circulator cc_end = T.incident_cells(*edge_it);
							TCell_circulator cc(cc_end);

							vertex_list[nedge].clear();

							edge_list[nedge] = Edge(
									int_abs(edge_it->get<0>()->vertex(edge_it->get<1>())->info()),
									int_abs(edge_it->get<0>()->vertex(edge_it->get<2>())->info()) );
							assert(edge_list[nedge].s1 != edge_list[nedge].s2);
							assert(edge_list[nedge].s1 <  edge_list[nedge].s2);

							do
							{
								if (T.is_infinite(cc))
								{
									ninf++; 
									failed_flag = true;
								}
								else
								{
									const TPoint c = T.dual(cc);
									const TVEC3 centre = TVEC3(
#ifndef _EXACT_CONSTRUCTIONS_
											c.x(), c.y(), c.z()
#else
											to_double(c.x()), to_double(c.y()), to_double(c.z())
#endif
											);
									vertex_list[nedge].push_back(centre);
									r2max = std::max(r2max, (centre - sites2process[isite].pos).norm2());
								}

								cc++;
							} while (cc != cc_end);

							nedge++;
						}

						TREAL rmax = 2.000001 * std::sqrt(r2max);

						// do not let search radius grow faster than 25% of the previous size
						rmax = std::min(rmax, 1.25 * sites2process[isite].rmax);

						// no mirror particles at this stage..., rmax < half-box-size
						assert(rmax < 0.5*global_domain_size.x);
						assert(rmax < 0.5*global_domain_size.y);
						assert(rmax < 0.5*global_domain_size.z);

						// store slightly higher rmax for caching
						sites2process[isite].rmax = 1.1*rmax;

						const TBOUNDARY bi(
								sites2process[isite].pos - rmax, 
								sites2process[isite].pos + rmax);

						if (!local_outer_box[ibox].isinbox(bi)) failed_flag = true;

						sites2process[isite].failed = failed_flag;

						if (failed_flag) 
						{
							nfailed++;
							continue;
						}

						// if cell is complete, compute its volume & faces area
						//
						TREAL cell_volume   = 0.0;
						TVEC3 cell_centroid = 0.0;
            int nngb_i = 0;
						for (int edge = 0; edge < nedge; edge++)
						{
							const int nvtx = vertex_list[edge].size();

							TVEC3 c = 0.0;
							for (int j = 0; j < nvtx; j++)
								c += vertex_list[edge][j];
							c *= 1.0/(TREAL)nvtx;

							vec3 norm(0.0);
							const vec3 &centroid = c;
							TVEC3 v1 = vertex_list[edge].back() - c;
							real area1 = 0.0;
							for (int j = 0; j < nvtx; j++)
							{
								const TVEC3 v2 = vertex_list[edge][j] - c;
								const TVEC3 norm3 = v1.cross(v2);
								const TREAL area3 = norm3.abs();
								norm  +=      norm3;
								area1 +=      area3;
								v1 = v2;
							}

							const real area0 = area1; // face.n.abs();
							const real L1 = std::sqrt(area0);
							const real L2 = (centroid - sites2process[isite].pos).abs();
							const real area = (L1 < SMALLDIFF1*L2) ? 0.0 : area0;

							const real nabs = norm.abs();
							if (area > 0.0 && nabs > 0.0)
							{
								const int jid_orig = edge_list[edge].ngb<false>(id);
								const int jid = int_abs(jid_orig);

								cell_local[id].ngb.push_back(jid);
								assert(nabs > 0.0);
							}
							else
								norm = 0.0;

							if (norm.abs() == 0.0)  continue;

							nngb_tot++;
              nngb_i++;

							const TVEC3 cv  = sites2process[isite].pos - centroid;
							v1    = vertex_list[edge].back() - centroid;
							const TREAL fourth = 1.0/4.0;
							for (int j = 0; j < nvtx; j++)
							{
								const TVEC3 v2 = vertex_list[edge][j] - centroid;
								const TVEC3 c4   = centroid + (v1 + v2 + cv) * fourth;
								const TREAL vol4 = std::abs(v1.cross(v2) * cv);
								cell_volume   +=      vol4;
								cell_centroid += c4 * vol4;
								v1 = v2;
							}

						}  // for edge < nedge
            nngb_min = std::min(nngb_min, nngb_i);
            nngb_max = std::max(nngb_max, nngb_i);

						cell_centroid *= 1.0/cell_volume;
						cell_volume   *= 1.0/6.0;

						assert(id < (int)local_n);

						cell_local[id].Volume   = cell_volume;
						cell_local[id].centroid = cell_centroid;
						ptcl_local[id].rmax     = sites2process[isite].rmax;
					}  // for site < nsite
				}
			}  // construction of Voronoi...
			int nfailed_loc = nfailed;
			MPI_Allreduce(&nfailed_loc, &nfailed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

#ifdef _DEBUG_PRINT_
			fprintf(stderr, "      proc= %d  nfailed= %d nsites_in_DT= %d <ngb>= %g (min= %d  max= %d) \n",
					myproc, nfailed_loc, nsites_in_DT, (double)nngb_tot/(double)local_n, nngb_min, nngb_max);
#endif
			nattempt++;
		} // main loop while (nfailed > 0)

		{
			unsigned long long nngb_loc = nngb_tot;
			unsigned long long nngb_glb;
      int nngb_min_glb;
      int nngb_max_glb;
			int nattempt_min, nattempt_max;
			MPI_Allreduce(&nngb_loc, &nngb_glb, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&nattempt, &nattempt_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
			MPI_Allreduce(&nattempt, &nattempt_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
			MPI_Allreduce(&nngb_min, &nngb_min_glb, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
			MPI_Allreduce(&nngb_max, &nngb_max_glb, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
			if (myproc == 0)
				fprintf(stderr, " <ngb>= %g  (min= %d max= %d) attempt= %d %d \n",
						(double)nngb_glb/(double)global_n, nngb_min_glb, nngb_max_glb,  nattempt_min, nattempt_max);
		}

		for (int ibox = 0; ibox < nbox; ibox++)
			for (int isite = ibox*NLEAF_LOC; isite < ibox*NLEAF_LOC + nbox_size[ibox]; isite++)
				active_sites[site_list[isite].id] = 0;

		return nattempt;
	} // build_mesh

#if 0
	bool system::get_cell_vtx_list(
			const Site &si,
			std::vector<TVEC3> &vtx_list)
	{
		assert(si.id < local_n);

		const TVertex_handle &vi = sites_in_DT[si.id];
		assert(vi != DT::Vertex_handle());
		assert(vi->info() == si.id);

		static std::vector<TCell_handle> Tcells;
		Tcells.clear();
		T.incident_cells(vi, std::back_inserter(Tcells));

		const int ncells = Tcells.size();
		bool success_flag = true;
		for (int icell = 0; icell < ncells; icell++)
		{
			const TCell_handle &ci = Tcells[icell];
			const bool is_infinite = T.is_infinite(ci);

			assert(!is_infinite);
			if (is_infinite) {success_flag = false; continue;}

			const real volume =
				std::abs(CGAL::Tetrahedron_3<K>( 
							ci->vertex(0)->point(), 
							ci->vertex(1)->point(), 
							ci->vertex(2)->point(), 
							ci->vertex(3)->point()).volume());
			if (volume == 0.0) {success_flag = false; continue;}

			bool local_success_flag = true;
			for (int iv = 0; iv < 4; iv++)
			{
				const TVertex_handle &v = ci->vertex(iv);
				if (is_infinite)
					if (T.is_infinite(v)) continue;
				if (v->info() >= local_n) local_success_flag = false;
				if (!local_success_flag) break;
			}
			success_flag &= local_success_flag;

#if 1
			if (!local_success_flag) 
				continue;
#endif

			const TPoint tc = T.dual(ci);
#ifndef _EXACT_CONSTRUCTIONS_
			const TVEC3 centre = TVEC3(tc.x(), tc.y(), tc.z());
#else
			const TVEC3 centre = TVEC3(to_double(tc.x()), to_double(tc.y()), to_double(tc.z()));
#endif
			vtx_list.push_back(centre);
		}

		return success_flag;
	}
#endif

#if 0
	bool system::remove_site(
			const Site &si,
			std::vector< std::pair<int, std::pair<real, real> > > &ngb_volume_list)   // oldVol, newVol
	{
		assert(si.id < local_n);

		const TVertex_handle &vi = sites_in_DT[si.id];
		assert(vi != DT::Vertex_handle());
		assert(vi->info() == si.id);

#if 0
		if (site_ngb_used.size() < import_site_list.size())
		{
			const int n1 = site_ngb_used.size();
			site_ngb_used.resize(import_site_list.size());
			const int n2 = site_ngb_used.size();
			for (int i = n1; i < n2; i++)
				site_ngb_used[i] = false;
		}
#endif

		static std::vector<TCell_handle> Tcells;
		Tcells.clear();
		T.incident_cells(vi, std::back_inserter(Tcells));

		static std::vector<TVertex_handle> site_ngb_list;
		site_ngb_list.clear();

		const int ncells = Tcells.size();
		bool success_flag = true;
		for (int icell = 0; icell < ncells; icell++)
		{
			const TCell_handle &ci = Tcells[icell];
			const bool is_infinite = T.is_infinite(ci);

			int idx = -1;
			for (int iv = 0; iv < 4; iv++)
				if (ci->vertex(iv) == vi)
					idx = iv;

			int iadd = 0;
			for (int iv = 0; iv < 4; iv++)
			{
				if (iv == idx) continue;

				const TVertex_handle &v = ci->vertex(iv);
				if (is_infinite)
					if (T.is_infinite(v)) continue;

				const int id = v->info();
				if (site_ngb_used[id]) continue;

				iadd++;
				site_ngb_used[id] = true;
				site_ngb_list.push_back(v);
				if (id >= local_n) success_flag = false;
				if (!success_flag) break;
			}
			assert(iadd < 4);
		}

		const int nngb = site_ngb_list.size();
		for (int j = 0; j < nngb; j++)
			site_ngb_used[site_ngb_list[j]->info()] = false;

		if (!success_flag) return false;

		////////////// COMPUTE OLD VOLUMES
		//

		static std::vector<TREAL> volume_old;
		volume_old.resize(nngb);
		for (int j = 0; j < nngb; j++)
		{
			const TVertex_handle &vi = site_ngb_list[j];

			// now compute new volumes of each of the neighbours

			static std::vector<TCell_handle> Tcells;
			Tcells.clear();
			T.incident_cells(vi, std::back_inserter(Tcells));

			static std::vector<int> site_ngb_list1;
			static std::vector<TEdge> edges;

			site_ngb_list1.clear();
			edges.clear();

			const int ncells = Tcells.size();
			for (int icell = 0; icell < ncells; icell++)
			{
				const TCell_handle &ci = Tcells[icell];
				const bool is_infinite = T.is_infinite(ci);

				int idx = -1;
				for (int iv = 0; iv < 4; iv++)
					if (ci->vertex(iv) == vi)
						idx = iv;

				int iadd = 0;
				for (int iv = 0; iv < 4; iv++)
				{
					if (iv == idx) continue;

					const TVertex_handle &v = ci->vertex(iv);
					if (is_infinite)
						if (T.is_infinite(v)) continue;

					const int id = v->info();
					if (site_ngb_used[id]) continue;

					iadd++;
					site_ngb_used[id] = true;
					site_ngb_list1.push_back(id);
					edges.push_back(TEdge(ci, idx, iv));
				}
				assert(iadd < 4);
			}

			const int nngb1	= site_ngb_list1.size();	
			for (int j1 = 0; j1 < nngb1; j1++)
				site_ngb_used[site_ngb_list1[j1]] = false;

			int nj = 0;
			TREAL volume_sj = 0.0;
			for (std::vector<TEdge>::iterator edge_it = edges.begin(); edge_it != edges.end(); edge_it++)
			{
				const TCell_circulator cc_end = T.incident_cells(*edge_it);
				TCell_circulator cc(cc_end);

				const int iv1 = edge_it->get<0>()->vertex(edge_it->get<1>())->info();
				const int iv2 = edge_it->get<0>()->vertex(edge_it->get<2>())->info();
				const int nj_id = (iv1 == vi->info()) ? iv2 : iv1;
				assert(nj_id == site_ngb_list1[nj++]);

				static std::vector<TVEC3> vertex_list;
				vertex_list.clear();
				TVEC3 c(0.0);
				do
				{
					assert(!T.is_infinite(cc));
					const TPoint tc = T.dual(cc);
#ifndef _EXACT_CONSTRUCTIONS_
					const TVEC3 centre = TVEC3(tc.x(), tc.y(), tc.z());
#else
					const TVEC3 centre = TVEC3(to_double(tc.x()), to_double(tc.y()), to_double(tc.z()));
#endif
					vertex_list.push_back(centre);
					c += centre;

					cc++;
				} while (cc != cc_end);

				const int nvtx = vertex_list.size();
				c *= 1.0/(TREAL)nvtx;

				TVEC3 normal(0.0);
				TVEC3 v1 = vertex_list.back() - c;
				for (int j = 0; j < nvtx; j++)
				{
					const TVEC3 v2 = vertex_list[j] - c;
					const TVEC3 norm3 = v1.cross(v2);
					normal += norm3;
					v1 = v2;
				}
				TREAL area = normal.abs();
				if (area == 0.0) continue;

				const TPoint pi = vi->point();
				const TVEC3 posj(pi.x(), pi.y(), pi.z());
				const TREAL volume = std::abs(normal * (posj - c));
				volume_sj += volume;
			}	
			volume_sj *= 1.0/6.0;
			volume_old[j] = volume_sj;
		}

		////////////// SANITY CHECK
		//

#if 1     // SANITY CHECK ...  \sum_j (Vnew_j - Vold_j) = V
#define _SANITY_CHECK_ENABLED_
		TREAL si_vol = 0.0;
		{
			static std::vector<int  > site_ngb_list;
			static std::vector<TEdge> edges;
			{
				static std::vector<TCell_handle> cells;
				cells.clear();
				T.incident_cells(vi, std::back_inserter(cells));

				site_ngb_list.clear();
				edges.clear();

				const int ncells = cells.size();
				for (int icell = 0; icell < ncells; icell++)
				{
					const TCell_handle &ci = cells[icell];

					const bool is_infinite = T.is_infinite(ci);

					int idx = -1;
					for (int iv = 0; iv < 4; iv++)
						if (ci->vertex(iv) == vi)
							idx = iv;

					int iadd = 0;
					for (int iv = 0; iv < 4; iv++)
					{
						if (iv == idx) continue;

						const TVertex_handle &v = ci->vertex(iv);
						const int id = v->info();
						if (site_ngb_used[id]) continue;
						if (is_infinite)
							if (T.is_infinite(v)) continue;

						iadd++;
						site_ngb_used[id] = true;
						site_ngb_list.push_back(id);
						edges.push_back(TEdge(ci, idx, iv));
					}
					assert(iadd < 4);
				}

				const int nngb = site_ngb_list.size();

				for (int i = 0; i < nngb; i++)
					site_ngb_used[site_ngb_list[i]] = false;
			}

			static std::vector<TVEC3> vertex_list;
			vertex_list.reserve(32);


			for (std::vector<TEdge>::iterator edge_it = edges.begin(); edge_it != edges.end(); edge_it++)
			{
				const TCell_circulator cc_end = T.incident_cells(*edge_it);
				TCell_circulator cc(cc_end);

				vertex_list.clear();

				TVEC3 c(0.0);
				do
				{
					assert(!T.is_infinite(cc));

					const TPoint tc = T.dual(cc);
					const TVEC3 centre = TVEC3(
#ifndef _EXACT_CONSTRUCTIONS_
							tc.x(), tc.y(), tc.z()
#else
							to_double(tc.x()), to_double(tc.y()), to_double(tc.z())
#endif
							);
					vertex_list.push_back(centre);
					c += centre;

					cc++;
				} while (cc != cc_end);

				const int nvtx = vertex_list.size();
				c *= 1.0/(TREAL)nvtx;

				TREAL area = 0.0;
				TVEC3 normal(0.0);
				TVEC3 v1 = vertex_list.back() - c;
				for (int j = 0; j < nvtx; j++)
				{
					const TVEC3 v2 = vertex_list[j] - c;
					const TVEC3 norm3 = v1.cross(v2);
					const TREAL area3 = norm3.abs();
					area   += area3;
					normal += norm3;
					v1 = v2;
				}
				if (area == 0.0) continue;
				assert(std::abs(normal.abs() -  area) <= SMALLDIFF1*area);

				si_vol += std::abs(normal * (si.pos - c));

			}

			si_vol *= 1.0/6.0;

		}
#endif


		////////////// REMOVE SITE
		//

		T.remove(vi);      // Remove site from DT

		////////////// COMPUTE NEW VOLUMES
		//

		TREAL si_vol_check = 0.0;
		for (int j = 0; j < nngb; j++)
		{
			const TVertex_handle &vi = site_ngb_list[j];

			// now compute new volumes of each of the neighbours

			static std::vector<TCell_handle> Tcells;
			Tcells.clear();
			T.incident_cells(vi, std::back_inserter(Tcells));

			static std::vector<int> site_ngb_list1;
			static std::vector<TEdge> edges;

			site_ngb_list1.clear();
			edges.clear();

			const int ncells = Tcells.size();
			for (int icell = 0; icell < ncells; icell++)
			{
				const TCell_handle &ci = Tcells[icell];
				const bool is_infinite = T.is_infinite(ci);

				int idx = -1;
				for (int iv = 0; iv < 4; iv++)
					if (ci->vertex(iv) == vi)
						idx = iv;

				int iadd = 0;
				for (int iv = 0; iv < 4; iv++)
				{
					if (iv == idx) continue;

					const TVertex_handle &v = ci->vertex(iv);
					if (is_infinite)
						if (T.is_infinite(v)) continue;

					const int id = v->info();
					if (site_ngb_used[id]) continue;
					iadd++;
					site_ngb_used[id] = true;
					site_ngb_list1.push_back(id);
					edges.push_back(TEdge(ci, idx, iv));
				}
				assert(iadd < 4);
			}

			const int nngb1	= site_ngb_list1.size();	
			for (int j1 = 0; j1 < nngb1; j1++)
				site_ngb_used[site_ngb_list1[j1]] = false;

			int nj = 0;
			TREAL volume_sj = 0.0;
			for (std::vector<TEdge>::iterator edge_it = edges.begin(); edge_it != edges.end(); edge_it++)
			{
				const TCell_circulator cc_end = T.incident_cells(*edge_it);
				TCell_circulator cc(cc_end);

				const int iv1 = edge_it->get<0>()->vertex(edge_it->get<1>())->info();
				const int iv2 = edge_it->get<0>()->vertex(edge_it->get<2>())->info();
				const int nj_id = (iv1 == vi->info()) ? iv2 : iv1;
				assert(nj_id == site_ngb_list1[nj++]);

				static std::vector<TVEC3> vertex_list;
				vertex_list.clear();
				TVEC3 c(0.0);
				do
				{
					assert(!T.is_infinite(cc));
					const TPoint tc = T.dual(cc);
#ifndef _EXACT_CONSTRUCTIONS_
					const TVEC3 centre = TVEC3(tc.x(), tc.y(), tc.z());
#else
					const TVEC3 centre = TVEC3(to_double(tc.x()), to_double(tc.y()), to_double(tc.z()));
#endif
					vertex_list.push_back(centre);
					c += centre;

					cc++;
				} while (cc != cc_end);

				const int nvtx = vertex_list.size();
				c *= 1.0/(TREAL)nvtx;

				TVEC3 normal(0.0);
				TVEC3 v1 = vertex_list.back() - c;
				for (int j = 0; j < nvtx; j++)
				{
					const TVEC3 v2 = vertex_list[j] - c;
					const TVEC3 norm3 = v1.cross(v2);
					normal += norm3;
					v1 = v2;
				}
				TREAL area = normal.abs();
				if (area == 0.0) continue;

				const TPoint pi = vi->point();
#ifndef _EXACT_CONSTRUCTIONS_
				const TVEC3 posj(pi.x(), pi.y(), pi.z());
#else
				const TVEC3 posj(to_double(pi.x()), to_double(pi.y()), to_double(pi.z()));
#endif
				const TREAL volume = std::abs(normal * (posj - c));
				volume_sj += volume;
			}	
			volume_sj *= 1.0/6.0;
			ngb_volume_list.push_back(std::make_pair(
						vi->info(), 
						std::make_pair((real)volume_old[j], (real)volume_sj)
						));
			const TREAL dv = volume_sj - volume_old[j];
			assert(dv >= -SMALLDIFF1*volume_old[j]);
			si_vol_check += dv;
		}

#ifdef _SANITY_CHECK_ENABLED_
		if (!(std::abs(si_vol - si_vol_check) <= SMALLDIFF1*si_vol))
		{
			fprintf(stderr, " si_vol= %g   si_vol_check= %g diff= %g [ %g ] \n",
					si_vol, si_vol_check, si_vol - si_vol_check, 
					(si_vol - si_vol_check)/si_vol);
		}
		assert(std::abs(si_vol - si_vol_check) <= SMALLDIFF1*si_vol);
#endif

		return true;
	}
#endif

#if 0
	bool system::insert_site(
			const Site &si,
			std::vector< std::pair<int, std::pair<real, real> > > &ngb_volume_list)   // oldVol, newVol
	{
		assert(si.id < local_n);
		if (si.id >= (int)sites_in_DT.size())
			sites_in_DT.resize(si.id+1);

		////////////// INSERT SITE
		//

		hint = T.insert(TPoint(si.pos.x, si.pos.y, si.pos.z), hint);
		hint->info() = si.id;
		sites_in_DT[si.id] = hint;


		////////////// ATTEMPT TO REMOVE IT, COMPUTES VOLUEMS & CHECKS FOR REMOTE NGB
		//

		static std::vector< std::pair<int, std::pair<real, real> > > volume;
		volume.clear();

		const bool success_flag = remove_site(si, volume);

		if (!success_flag)  // if failed to remove, means one of the ngb are from remote proc
		{
			T.remove(hint);    // cannot insert site, forcefully remove it
			return false;
		}

		////////////// INSERT SITE AGAIN
		//

		hint = T.insert(TPoint(si.pos.x, si.pos.y, si.pos.z), hint);
		hint->info() = si.id;
		sites_in_DT[si.id] = hint;

		////////////// STORE VOLUMES
		//

		const int nngb = volume.size();
		for (int j = 0; j < nngb; j++)
			ngb_volume_list.push_back(
					std::make_pair(volume[j].first,
						std::make_pair(volume[j].second.second, volume[j].second.first)
						)
					);

		return true;

	}
#endif

}

