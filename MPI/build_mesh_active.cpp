#include "fvmhd3d.h"
#define SMALLDIFF 1.0e-16

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h> 

#include <CGAL/Spatial_sort_traits_adapter_3.h>
#include <boost/iterator/counting_iterator.hpp>

// #include <CGAL/hilbert_sort.h>
#include <CGAL/spatial_sort.h>

#if 0
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
#endif

namespace fvmhd3d
{
	typedef CGAL::Exact_predicates_inexact_constructions_kernel K; 

	typedef CGAL::Triangulation_vertex_base_with_info_3<std::pair<int, int>, K> Vb; 
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
	
	typedef CGAL::Spatial_sort_traits_adapter_3<K, TPoint*> Search_traits_3;

#if 0
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
#endif

#if 1
	DT *T1_ptr;
#else
	DT T1;
#endif
  unsigned long long n_in_DT1;
	std::vector<TVertex_handle> sites_in_DT1;
//	TVertex_handle hint1;
	std::vector<bool> sites_ngb_used;
	std::vector<unsigned long long> sites_added_list;

	double dt_tree, dt_ngb, dt_prep, dt_sort, dt_tria, dt_voro;

  int system::build_mesh(
      std::vector< std::vector<int> > active_list,
      const Octree::Tree<NLEAF_LOC, 128> &tree,
      const bool store_ngb)
  {
#if 1
		DT &T1 = *T1_ptr;
#endif
    const int nimport = site_import.size();

    int nattempt = 0;
    while (active_list.size() > 0)
    {
      std::vector< std::pair<TPoint, long long> >points2insert;
      std::vector<boundary> outer_box_vec(active_list.size());

			double t0 = mytimer::get_wtime();
      for (int box = 0; box < (const int)active_list.size(); box++)
      {
        const std::vector<int> &list = active_list[box];
        const int nlist = list.size();
        assert(nlist <= NLEAF_LOC); 
        boundary &outer_box = outer_box_vec[box];
        for (int ip = 0; ip < nlist; ip++)
        {
          const int i = list[ip];
          const vec3 ipos = site_import[i].pos();
          const real rmax = site_import[i].rmax()*1.01;
          assert(rmax > 0.0);
          outer_box.merge(boundary(ipos - rmax, ipos + rmax));
        }

        const vec3 outer_box_centre = outer_box.centre();
        const vec3 outer_box_hsize  = outer_box.hsize();
#if 0
        if (
            !(outer_box_hsize.x < 0.5*global_domain_size.x) ||
            !(outer_box_hsize.y < 0.5*global_domain_size.y) ||
            !(outer_box_hsize.z < 0.5*global_domain_size.z) )
          fprintf(stderr, " myproc= %d: attempt= %d\n", myproc, nattempt);

        assert(outer_box_hsize.x < 0.5*global_domain_size.x);
        assert(outer_box_hsize.y < 0.5*global_domain_size.y);
        assert(outer_box_hsize.z < 0.5*global_domain_size.z);
#endif

        static std::vector<LightSite> sites_lst;
        sites_lst.clear();
#if 0
        const bool pflag = tree.root.walk_boundary_flag(outer_box, sites_lst, global_domain_size);
#else
        const bool pflag = true;
        tree.root.walk_boundary(outer_box, sites_lst, global_domain_size);
#endif

        const int nsites_lst = sites_lst.size();
        for (int jsite = 0; jsite < nsites_lst; jsite++)
        {
          const LightSite &s = sites_lst[jsite];
          assert(s.id >= 0);
          assert(s.id < nimport);

          TVEC3 dr = s.pos - outer_box_centre;
          unsigned long long oct = 0;
          
          if (pflag)
          {
            if      (dr.x >  outer_box_hsize.x) {dr.x -= global_domain_size.x; oct +=  1;}
            else if (dr.x < -outer_box_hsize.x) {dr.x += global_domain_size.x; oct +=  2;}
            if      (dr.y >  outer_box_hsize.y) {dr.y -= global_domain_size.y; oct +=  4;}
            else if (dr.y < -outer_box_hsize.y) {dr.y += global_domain_size.y; oct +=  8;}
            if      (dr.z >  outer_box_hsize.z) {dr.z -= global_domain_size.z; oct += 16;}
            else if (dr.z < -outer_box_hsize.z) {dr.z += global_domain_size.z; oct += 32;}
          }
          const int oct0 = oct;
          oct = 1LLU << oct;

          assert(s.id >= 0);
          assert(s.id < (int)sites_added_list.size());
          unsigned long long &value = sites_added_list[s.id];
          if ((value & oct) == oct) continue;
          value |= oct;

          const vec3 jpos = dr + outer_box_centre;
          points2insert.push_back(std::make_pair(TPoint(jpos.x, jpos.y, jpos.z), oct0 == 0 ? s.id : -1-s.id));
        }
      }
			dt_tree += mytimer::get_wtime() - t0;

			t0 = mytimer::get_wtime();
      {
        const size_t n0 = sites_ngb_used.size();
        sites_ngb_used.resize(n_in_DT1 + points2insert.size());
        const size_t n1 = sites_ngb_used.size();
        for (size_t i = n0; i < n1; i++)
          sites_ngb_used[i] = false;
      }
			dt_ngb += mytimer::get_wtime() - t0;

			t0 = mytimer::get_wtime();
      const int np = points2insert.size();
      std::vector< TPoint >                        pnts  (np);
			std::vector<  int   >                        pnts_idx(np);
      std::vector< std::pair<const TPoint*, int> > points(np);
      for (int ip = 0; ip < np; ip++)
      {
        pnts  [ip]        = points2insert[ip].first;
        points[ip].first  = &pnts[ip];
        points[ip].second = points2insert[ip].second;
				pnts_idx[ip] = points2insert[ip].second;
      }
			dt_prep += mytimer::get_wtime() - t0;

#if 0
			t0 = mytimer::get_wtime();
      spatial_sort(points.begin(), points.end(), Traits_for_spatial_sort<DT>());
			dt_sort += mytimer::get_wtime() - t0;
#endif
#if 0
			t0 = mytimer::get_wtime();
//			CGAL::hilbert_sort(points.begin(), points.end(), CGAL::Hilbert_sort_median_policy(), Traits_for_spatial_sort<DT>());
			CGAL::spatial_sort(points.begin(), points.end(), Traits_for_spatial_sort<DT>(), CGAL::Hilbert_sort_median_policy());
			dt_sort += mytimer::get_wtime() - t0;
#endif

#if 1
			t0 = mytimer::get_wtime();
			std::vector<std::ptrdiff_t> indices;
			indices.reserve(pnts.size());
			std::copy(
					boost::counting_iterator<std::ptrdiff_t>(0),
					boost::counting_iterator<std::ptrdiff_t>(pnts.size()),
					std::back_inserter(indices));
			CGAL::spatial_sort(indices.begin(), indices.end(), Search_traits_3(&(pnts[0])), CGAL::Hilbert_sort_median_policy());
			dt_sort += mytimer::get_wtime() - t0;


			DT::Vertex_handle hint1;
			t0 = mytimer::get_wtime();
			for (std::vector<std::ptrdiff_t>::iterator it = indices.begin(); it != indices.end(); it++)
      {
        hint1        = T1.insert(pnts[*it], hint1);
        assert(hint1 != TVertex_handle());
        const int jid = pnts_idx[*it] < 0 ? -1-pnts_idx[*it] : pnts_idx[*it];
        hint1->info() = std::make_pair(n_in_DT1++, jid);

        if (pnts_idx[*it] >= 0)
        { 
          assert(jid < nimport);
          sites_in_DT1[jid] = hint1;
        }
      }
			dt_tria += mytimer::get_wtime() - t0;
#endif

#if 0
			DT::Vertex_handle hint1;
			t0 = mytimer::get_wtime();
      for (int ip = 0; ip < np; ip++)
      {
        hint1        = T1.insert(*points[ip].first, hint1);
        assert(hint1 != TVertex_handle());
        const int jid = points[ip].second < 0 ? -1-points[ip].second : points[ip].second;
        hint1->info() = std::make_pair(n_in_DT1++, jid);

        if (points[ip].second >= 0)
        { 
          assert(jid < nimport);
          sites_in_DT1[jid] = hint1;
        }
      }
			dt_tria += mytimer::get_wtime() - t0;
#endif

      std::vector< std::vector<int> > active_list_new;

      std::vector<TEdge> Tedges;
      std::vector<TCell_handle> Tcells;
      std::vector<int> site_ngb_list;

			t0 = mytimer::get_wtime();
      for (int box = 0; box < (const int)active_list.size(); box++)
      {
        const std::vector<int> &list = active_list[box];
        const int nlist = list.size();
        std::vector<int> list_new;

        const boundary &outer_box = outer_box_vec[box];

        for (int ip = 0; ip < nlist; ip++)
        {
          const int isite = list[ip];
          const TVertex_handle vi = sites_in_DT1[isite];
          assert(vi != TVertex_handle());
          assert(vi->info().second == isite);
          const vec3 ipos(vi->point().x(), vi->point().y(), vi->point().z());
          Tedges.clear();
          Tcells.clear();
          site_ngb_list.clear();
          T1.incident_cells(vi, std::back_inserter(Tcells));

          const int ncells = Tcells.size();
          for (int icell = 0; icell < ncells; icell++)
          {
            const TCell_handle &ci = Tcells[icell];
            const bool is_infinite = T1.is_infinite(ci);

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
                if (T1.is_infinite(v)) continue;
              assert (v != TVertex_handle());

              const int id = v->info().first;
              assert(id >= 0);
              assert(id < (int)n_in_DT1);
              if (sites_ngb_used[id]) continue;

              iadd++;
              sites_ngb_used[id] = true;
              site_ngb_list.push_back(id);
              Tedges.push_back(TEdge(ci, idx, iv));
            }
            assert(iadd < 4);
          }
          const int nngb = site_ngb_list.size();
          for (int j = 0; j < nngb; j++)
            sites_ngb_used[site_ngb_list[j]] = false;

          TREAL r2max = 0.0;

          const int NMAXEDGE = 1024;
          static std::vector<TVEC3> vertex_list[NMAXEDGE];

          bool fail = false;
          int nedge = 0;
          // for each Delaunay edge trace it dual Voronoi face
          int cnt0 = 0;
          for (std::vector<TEdge>::iterator edge_it = Tedges.begin(); edge_it != Tedges.end(); edge_it++)
          {
            const TCell_circulator cc_end = T1.incident_cells(*edge_it);
            TCell_circulator cc(cc_end);
            vertex_list[nedge].clear();

            // compute face vertices
            do
            {
              cnt0++;
              if (T1.is_infinite(cc))
              {
                fail = true;
                r2max = HUGE;
              }
              else
              {
                const TPoint c = T1.dual(cc);
                const TVEC3 centre = TVEC3(c.x(), c.y(), c.z());
                r2max = std::max(r2max, (centre - ipos).norm2());

                vertex_list[nedge].push_back(centre);
              }

              cc++;
            } while (cc != cc_end);
            nedge++;
            assert(nedge < NMAXEDGE);
          }  // for edge < nedge

          if (cnt0 == 0) r2max = HUGE;
          assert(r2max > 0.0);
          real rmax = 2.01*std::sqrt(r2max);
#if 1
          rmax = std::min(rmax, 
              std::min(std::min(global_domain_size.x, global_domain_size.y),
                global_domain_size.z)/4.0);
#else
          assert(site_import[isite].rmax() > 0.0);
          rmax = std::min(rmax, 1.25 * site_import[isite].rmax());
#endif

          // no mirror particles at this stage..., rmax < half-box-size
          assert(rmax < 0.5*global_domain_size.x);
          assert(rmax < 0.5*global_domain_size.y);
          assert(rmax < 0.5*global_domain_size.z);

          assert(isite < nimport);
          site_import[isite].rmax() = rmax*1.03;

          const boundary bi(ipos - rmax, ipos + rmax);
          if (!outer_box.isinbox(bi)) fail = true;

          if (fail)
          {
            if (nattempt > 20)
              site_import[isite].rmax() = 
                std::min(std::min(global_domain_size.x, global_domain_size.y), global_domain_size.z)/4.0;
            list_new.push_back(isite);
            continue;
          }

          TREAL volume = 0.0;
          TVEC3 centroid = 0.0;
          vec3 sum_area  = 0.0;
          int edge = 0;
          int nface_built = 0.0;
          for (std::vector<TEdge>::iterator edge_it = Tedges.begin(); edge_it != Tedges.end(); edge_it++, edge++)
          {
            assert(edge < nedge);
            const int nvtx = vertex_list[edge].size();
            const int i1 = edge_it->get<0>()->vertex(edge_it->get<1>())->info().second;
            const int i2 = edge_it->get<0>()->vertex(edge_it->get<2>())->info().second;
            assert(i1 >= 0); 
            assert(i1 < nimport);
            assert(i2 >= 0); 
            assert(i2 < nimport);
            Face face;
            if (i1 == isite)
            {
              face.s1 = isite;
              face.s2 = i2;
            }
            else
            {
              face.s1 = i2;
              face.s2 = i1;
            }
            assert(face.s1 == isite);
            assert(face.s2 != isite);



            int face_id = -1;
            const int nj = cell_list[isite].ngb.size();

            for (int j = 0; j < nj; j++)
              if (face_list[cell_list[isite].faces()[j]].ngb<false>(isite) == face.s2)
              {
                face_id = cell_list[isite].faces()[j];
                break;
              }

            if (face_id != -1)        // if face found, get it from the list
            {
              assert(face_id >= 0 && face_id < (int)face_list.size());
              face = face_list[face_id];
            } 
            else                      // else compute new face
            {
              TVEC3 c = 0.0;
              for (int j = 0; j < nvtx; j++)
                c += vertex_list[edge][j];
              c *= 1.0/(TREAL)nvtx;

              face.n        = 0.0;
              face.centroid = 0.0;
              TVEC3 v1 = vertex_list[edge].back() - c;
              const TREAL third = 1.0/3.0;
              real area1 = 0.0;
              for (int j = 0; j < nvtx; j++)
              {
                const TVEC3 v2 = vertex_list[edge][j] - c;

                const TVEC3 norm3 = v1.cross(v2);
                const TREAL area3 = norm3.abs();
                const TVEC3 c3    = c + (v1 + v2) * third;

                face.n        +=      norm3;
                face.centroid += c3 * area3;
                area1         +=      area3;

                v1 = v2;
              }

              const TREAL area0 = area1;
              const TREAL L1 = std::sqrt(area0);
              const TREAL L2 = (face.centroid - ipos).abs();
              const TREAL area = (L1 < SMALLDIFF*L2) ? 0.0 : area0;

              const real nabs = face.n.abs();
              if (area > 0.0 && nabs > 0.0)
              {
                assert(nabs > 0.0);
                const TREAL iarea = 1.0/area;
                face.centroid *= iarea;
                face.n        *= 0.5;

                if ((face.centroid - ipos)*face.n < 0.0)
                  face.n *= -1.0;

                assert(face.s2 >= 0);
                const int jid = face.s2;
                assert(jid >= 0);
                assert(jid < nimport);

                if (store_ngb && !site_import[jid].is_active())
                {
                  site_import[jid].set_local_ngb();
                  site_import[isite].set_hasngb();
                }
              }
              else
              {
                face.n = 0.0;
              }
            }

            if (face.area() == 0.0)  continue;

            nface_built++;

            if (face_id == -1)
            {
              assert(face.s1 == isite);
              assert(face.s2 < nimport);
              face_list.push_back(face);

              const vec3 dri = face.centroid - site_import[face.s1].pos(); //Wrec_import[face.s1].pos;
              assert(std::abs(dri.x) < 0.5*global_domain_size.x);
              assert(std::abs(dri.y) < 0.5*global_domain_size.y);
              assert(std::abs(dri.z) < 0.5*global_domain_size.z);
              assert(dri*face.n > 0.0);


#if 0
              if (store_ngb)
                assert(ptcl_import[face.s1].boundary == 0);
#endif

              cell_list[face.s1].ngb.push_back(face_list.size() - 1);
              cell_list[face.s2].ngb.push_back(face_list.size() - 1);
            }

            TVEC3 cv  = ipos - face.centroid;

#if 1	
            if      (cv.x >  0.5*global_domain_size.x) face.centroid.x += global_domain_size.x;
            else if (cv.x < -0.5*global_domain_size.x) face.centroid.x -= global_domain_size.x;
            if      (cv.y >  0.5*global_domain_size.y) face.centroid.y += global_domain_size.y;
            else if (cv.y < -0.5*global_domain_size.y) face.centroid.y -= global_domain_size.y;
            if      (cv.z >  0.5*global_domain_size.z) face.centroid.z += global_domain_size.z;
            else if (cv.z < -0.5*global_domain_size.z) face.centroid.z -= global_domain_size.z;

            cv  = ipos - face.centroid;
            assert(std::abs(cv.x) <= global_domain_size.x);
            assert(std::abs(cv.y) <= global_domain_size.y);
            assert(std::abs(cv.z) <= global_domain_size.z);
#endif

            TVEC3 v1  = vertex_list[edge].back() - face.centroid;
            const TREAL fourth = 1.0/4.0;
            for (int j = 0; j < nvtx; j++)
            {
              const TVEC3 v2 = vertex_list[edge][j] - face.centroid;

              const TVEC3 c4   = face.centroid + (v1 + v2 + cv) * fourth;
              const TREAL vol4 = std::abs(v1.cross(v2) * cv);

              volume   +=      vol4;
              centroid += c4 * vol4;

              v1 = v2;
            }

          }  // for edge < nedge
          assert(edge == nedge);

          if (volume == 0.0)
            assert(nface_built == 0.0);
          assert(volume != 0.0);
          assert(volume > 0.0);
          centroid *= 1.0/volume;
          volume   *= 1.0/6.0;

          cell_list[isite].Volume = volume;
          cell_list[isite].centroid = centroid;
        }

        if (!list_new.empty())
          active_list_new.push_back(list_new);
      }
			dt_voro += mytimer::get_wtime() - t0;

      active_list.swap(active_list_new);

      nattempt++;
      if (!(nattempt < 36))
      {
        fprintf(stderr, " n_in_DT= %lld \n", n_in_DT1);
        fprintf(stderr, " nactive= %d nimport= %d\n", (int)active_ptcl.size(), (int)site_import.size());
        fprintf(stderr, "active_list.size= %d \n", (int)active_list.size());
        for (int box = 0; box < (const int)active_list.size(); box++)
        {
          const std::vector<int> &list = active_list[box];
          const int nlist = list.size();
          assert(nlist <= NLEAF_LOC); 
          boundary outer_box;
          for (int ip = 0; ip < nlist; ip++)
          {
            const int i = list[ip];
            const vec3 ipos = site_import[i].pos();
            const real rmax = site_import[i].rmax()*1.03; //*1.5; //0323;
            assert(rmax > 0.0);
            outer_box.merge(boundary(ipos - rmax, ipos + rmax));
            fprintf(stderr, " i= %d [ %d ] : pos= %g %g %g [ %g ]  rmax= %g \n",
                i, nlist, ipos.x, ipos.y, ipos.z, ipos.abs(), rmax);
          }
          std::vector<LightSite> sites_lst, sites_lst1;
          tree.root.walk_boundary(outer_box, sites_lst, global_domain_size);
          tree.root.walk_boundary(global_domain, sites_lst1, global_domain_size);
          fprintf(stderr, "n1= %d  n2= %d \n",
              (int)sites_lst.size(), (int)sites_lst1.size());
          fprintf(stderr, " points2insert= %d \n", (int)points2insert.size());

          const vec3 outer_box_centre = outer_box.centre();
          const vec3 outer_box_hsize  = outer_box.hsize();
          fprintf(stderr, " box= %d : c= %g %g %g [ %g ]  hsize= %g %g %g [ %g ]\n",
              box, 
              outer_box_centre.x,
              outer_box_centre.y,
              outer_box_centre.z,
              outer_box_centre.abs(),
              outer_box_hsize.x,
              outer_box_hsize.y,
              outer_box_hsize.z,
              outer_box_hsize.abs());
        }
      }
      assert(nattempt < 36);
    }

    return nattempt;
  }

  void system::build_mesh_active(
      std::vector<             int  > &active_sites,
      std::vector< std::vector<int> > &active_list,
      std::vector<             int  > &active_sites_ngb)
  {
    const int nimport = site_import.size();

#if 1
		DT T1;
		T1_ptr = &T1;
#endif


    const double t00 = mytimer::get_wtime();

    //    cell_list.clear();
    assert(cell_list.empty());
    assert(face_list.empty());
    cell_list.resize(nimport);
    //    face_list.clear();

    const double t05 = mytimer::get_wtime();

    //    T1.clear();
    n_in_DT1 = 0;
    sites_in_DT1.resize(nimport);
//    hint1 = TVertex_handle();
    sites_added_list.resize(nimport);
    for (int i = 0; i < nimport; i++)
      sites_added_list[i] = 0;

		dt_tree = dt_ngb = dt_prep = dt_sort = dt_tria = dt_voro = 0.0;

    const double t10 = mytimer::get_wtime();
    const int nattempt_act = build_mesh(active_list, import_tree, true);


    const double t20 = mytimer::get_wtime();
#if 0
    MPI_Barrier(MPI_COMM_WORLD);
    if (myproc == 0)
      fprintf(stderr, " ngb -- \n");
    MPI_Barrier(MPI_COMM_WORLD);
#endif

#if 0 
    T1.clear();
    n_in_DT1 = 0;
    site_mirror_map.resize(nimport);
    for (int i = 0; i < nimport; i++)
      site_mirror_map[i] = i;

    sites_in_DT1.resize(nimport);
    hint1 = TVertex_handle();
    sites_added_list.resize(nimport);
    for (int i = 0; i < nimport; i++)
      sites_added_list[i] = 0;
#endif


#if 0
    std::vector< int > local_ngb_list;
    std::vector< std::vector<int> > active_list_ngb;
    const int nleaves = import_tree.n_leaves;
    for (int leaf = 0; leaf < nleaves; leaf++)
    {
      boundary outer_box;
      std::vector<int> list;
      for (Octree::Body *bp = import_tree.leaf(leaf)->pfirst; bp != NULL; bp = bp->next)
        if (site_import[bp->id].is_local_ngb())
        {
          list.push_back(bp->id);
          active_sites_ngb.push_back(bp->id);
          local_ngb_list.push_back(bp->id);
          const int i = bp->id;
          const vec3 ipos = site_import[i].pos();
          const real rmax = site_import[i].rmax()*1.0323;
          assert(rmax > 0.0);
          outer_box.merge(boundary(ipos - rmax, ipos + rmax));
        }
      const vec3 outer_box_centre = outer_box.centre();
      const vec3 outer_box_hsize  = outer_box.hsize();

      if (
          !(outer_box_hsize.x < 0.5*global_domain_size.x) ||
          !(outer_box_hsize.y < 0.5*global_domain_size.y) ||
          !(outer_box_hsize.z < 0.5*global_domain_size.z) )
      {
        fprintf(stderr, " myproc= %d: attempt= %d\n", myproc, -1);
        int pc = 0;
        int pc1 = 0;
        for (Octree::Body *bp = import_tree.leaf(leaf)->pfirst; bp != NULL; bp = bp->next)
        {
          if (site_import[bp->id].is_local_ngb())
          {
            const int i = bp->id;
            const vec3 ipos = site_import[i].pos();
            const real rmax = site_import[i].rmax()*1.0323;
            assert(rmax > 0.0);
            fprintf(stderr, "j= %d %d : ipos= %g %g %g [ %g ] rmax = %g \n",
                pc++, pc1, ipos.x, ipos.y, ipos.z, ipos.abs(), rmax);
          }
          pc1++;
        }
        fprintf(stderr, "pc= %d  ppc1= %d \n", pc, pc1);
        outer_box.dump(stderr, true);
        global_domain.dump(stderr, true);
        Boundary<fvec3>(
            import_tree.leaf(leaf)->centre-1.0001*import_tree.leaf(leaf)->hsize,
            import_tree.leaf(leaf)->centre+1.0001*import_tree.leaf(leaf)->hsize).dump(stderr,true);
        fprintf(stderr,"centre= %g  %g %g  hsize= %g\n",
            import_tree.leaf(leaf)->centre.x,
            import_tree.leaf(leaf)->centre.y,
            import_tree.leaf(leaf)->centre.z,
            import_tree.leaf(leaf)->hsize);
        import_tree.leaf(leaf)->inner.dump(stderr,true);

      }

#if 1
      assert(outer_box_hsize.x < 0.5*global_domain_size.x);
      assert(outer_box_hsize.y < 0.5*global_domain_size.y);
      assert(outer_box_hsize.z < 0.5*global_domain_size.z);
#endif

      assert(list.size() <= NLEAF_LOC);
      if (!list.empty())
        active_list_ngb.push_back(list);
    }


    const int nattempt_ngb = build_mesh(active_list_ngb, import_tree, false);

    for (size_t i = 0; i < local_ngb_list.size(); i++)
      site_import[local_ngb_list[i]].unset_local_ngb();

#if 1
    fprintf(stderr, " proc= %d : nactive= %d  ngb_list= %d [ %d ] n_in_DT= %lld   nimport= %d  nattempt= %d %d\n", 
        myproc, (int)active_sites.size(), (int)local_ngb_list.size(), (int)active_sites_ngb.size(), n_in_DT1, nimport,
        nattempt_act, nattempt_ngb);
#endif
#endif
#if 0
    fprintf(stderr, " proc= %d : nactive= %d  ngb_list= %d [ %d ] n_in_DT= %lld   nimport= %d  nattempt= %d %d\n", 
        myproc, (int)active_sites.size(), 0, 0, n_in_DT1, nimport,
        nattempt_act, 0);
#endif

    const double t30 = mytimer::get_wtime();
#if 0
		T1.clear();
#endif
    clear_vec(sites_in_DT1);
    clear_vec(sites_added_list);
    clear_vec(sites_ngb_used);
    const double t35 = mytimer::get_wtime();

#if 0
		if (nactive_glb > 0.2*global_n)
		{
			fprintf(stderr, "myproc= %d:  nattempt= %d n= %d  dn= %d n_in_DT= %d:: tall= %g [ res= %g fill= %g mesh=%g  mesh_ngb= %g  clear= %g ] \n",
					myproc, nattempt_act, (int)active_ptcl.size(), (int)(site_import.size() - active_ptcl.size()), (int)n_in_DT1, t35- t00, t05-t00, t10-t05, t20-t10, t30-t20, t35-t30);
		}
#endif

#if 0
		if (nactive_glb > 0.2*global_n)
		{
			fprintf(stderr, "myproc= %d:  nattempt= %d n= %d  dn= %d n_in_DT= %d:: tall= %g [ mesh= %g :: tree= %g  ngb= %g prep= %g sort= %g tria= %g voro= %g] \n",
					myproc, nattempt_act, (int)active_ptcl.size(), (int)(site_import.size() - active_ptcl.size()), (int)n_in_DT1,  t35-t00, t20-t10, dt_tree, dt_ngb, dt_prep, dt_sort, dt_tria, dt_voro);
		}
#endif

#if 0
		if (active_sites.size() == 1)
		{
			fprintf(stderr, " myproc= %d: ngb_list_size= %d \n", myproc, (int)cell_list[active_sites[0]].ngb.size());
		}
#endif

		return;
	}


	////////////////////////////////////////
	////////////////////////////////////////
	////////////////////////////////////////
	////////////////////////////////////////

#if 0
	bool system::remove_site(
			const int id,
			std::vector< std::pair<int, std::pair<real, real> > > &ngb_volume_list) // j, (oldVol, newVol)
			{
				const int nimport = site_import.size();
				assert(id >= 0);
				assert(id < nimport);

				assert(site_import[id].is_active());

				const TVertex_handle &vi = sites_in_DT1[id];
				assert(vi != DT::Vertex_handle());
				assert(vi->info() == id);

				static std::vector<TEdge> edges_i;
				edges_i.clear();
				T1.incident_edges(vi, std::back_inserter(edges_i));

				T1.remove(vi);      // Remove site from DT

				////////////// COMPUTE NEW VOLUMES

				const TVertex_handle &vi0 = vi;
				real tot_vol_check = 0.0;
				for (std::vector<TEdge>::iterator eij = edges_i.begin(); eij != edges_i.end(); eij++)
				{
					const TVertex_handle &vi1 = eij->get<0>()->vertex(eij->get<1>());
					const TVertex_handle &vi2 = eij->get<0>()->vertex(eij->get<2>());
					assert(vi1 == vi0 || vi2 == vi0);

					const TVertex_handle &vi = (vi0 == vi2) ? vi1 : vi2;
					assert(vi != DT::Vertex_handle());
					assert(vi->info() >= 0);
					assert(vi->info() < (int)site_mirror_map.size()); 

					const int j = site_map(vi->info());
					assert(j >= 0);
					assert(j <  nimport);

					cell_list[j] = Cell();

					// now compute new volumes of each of the neighbours

					static std::vector<TEdge> edges;
					edges.clear();
					T1.incident_edges(vi, std::back_inserter(edges));

					TREAL volume_sj = 0.0;
					for (std::vector<TEdge>::iterator edge_it = edges.begin(); edge_it != edges.end(); edge_it++)
					{
						const TCell_circulator cc_end = T1.incident_cells(*edge_it);
						TCell_circulator cc(cc_end);

						const TVertex_handle &vi1 = edge_it->get<0>()->vertex(edge_it->get<1>());
						const TVertex_handle &vi2 = edge_it->get<0>()->vertex(edge_it->get<2>());
						assert(vi1 == vi || vi2 == vi);

						const TVertex_handle &vj = (vi == vi2) ? vi1 : vi2;
						assert(vj != DT::Vertex_handle());
						assert(vj->info() >= 0);
						assert(vj->info() < (int)site_mirror_map.size()); 

						const int j1 = site_map(vj->info());
						assert(j1 >= 0);
						assert(j1 <  nimport);

						cell_list[j].ngb.push_back(-1-j1);     // negative, means it is ngb, else face index

						static std::vector<TVEC3> vertex_list;
						vertex_list.clear();
						TVEC3 c(0.0);
						do
						{
							assert(!T1.is_infinite(cc));
							const TPoint tc = T1.dual(cc);
							const TVEC3 centre = TVEC3(tc.x(), tc.y(), tc.z());
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
					ngb_volume_list.push_back(std::make_pair(
								vi->info(), 
								std::make_pair(cell_list[j].Volume, volume_sj)
								));
					const TREAL dv = volume_sj - cell_list[j].Volume;
					assert(dv >= -SMALLDIFF*cell_list[j].Volume);
					tot_vol_check += dv;
					cell_list[j].Volume = volume_sj;
				}

				if (!(std::abs(cell_list[id].Volume - tot_vol_check) <= SMALLDIFF*cell_list[id].Volume))
				{
					const real si_vol = cell_list[id].Volume;
					const real si_vol_check = tot_vol_check;
					fprintf(stderr, " si_vol= %g   si_vol_check= %g diff= %g [ %g ] \n",
							si_vol, si_vol_check, si_vol - si_vol_check, 
							(si_vol - si_vol_check)/si_vol);
					assert(std::abs(si_vol - si_vol_check) <= SMALLDIFF*si_vol);
				}

				return true;
			}
#endif

#if 0

	bool system::insert_site(
			const vec3 &pos, const int id,
			std::vector< std::pair<int, std::pair<real, real> > > &ngb_volume_list)   // oldVol, newVol
	{
		const int nimport = site_import.size();

		hint1 = T1.insert(TPoint(pos.x, pos.y, pos.z), hint1);
		hint1->info() = id;

		// compute new volumes

		const TVertex_handle &vi = hint1;

		static std::vector<TEdge> edges_i;
		edges_i.clear();
		T1.incident_edges(vi, std::back_inserter(edges_i));

		cell_list[id] = Cell();

		////////////// COMPUTE NEW i-VOLUMES

		TREAL volume_si = 0.0;
		for (std::vector<TEdge>::iterator eij = edges_i.begin(); eij != edges_i.end(); eij++)
		{
			const TCell_circulator cc_end = T1.incident_cells(*eij);
			TCell_circulator cc(cc_end);

			const TVertex_handle &vi1 = eij->get<0>()->vertex(eij->get<1>());
			const TVertex_handle &vi2 = eij->get<0>()->vertex(eij->get<2>());
			assert(vi1 == vi || vi2 == vi);

			const TVertex_handle &vj = (vi == vi2) ? vi1 : vi2;
			assert(vj != DT::Vertex_handle());
			assert(vj->info() >= 0);
			assert(vj->info() < (int)site_mirror_map.size()); 

			const int j = site_map(vj->info());
			assert(j >= 0);
			assert(j <  nimport);

			cell_list[id].ngb.push_back(-1-j);     // negative, means it is ngb, else face index

			static std::vector<TVEC3> vertex_list;
			vertex_list.clear();
			TVEC3 c(0.0);
			do
			{
				assert(!T1.is_infinite(cc));
				const TPoint tc = T1.dual(cc);
				const TVEC3 centre = TVEC3(tc.x(), tc.y(), tc.z());
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
			volume_si += volume;
		}	
		volume_si *= 1.0/6.0;
		cell_list[id].Volume = volume_si;

		////////////// COMPUTE NEW j-VOLUMES


		const TVertex_handle &vi0 = vi;
		real tot_vol_check = 0.0;
		for (std::vector<TEdge>::iterator eij = edges_i.begin(); eij != edges_i.end(); eij++)
		{
			const TVertex_handle &vi1 = eij->get<0>()->vertex(eij->get<1>());
			const TVertex_handle &vi2 = eij->get<0>()->vertex(eij->get<2>());
			assert(vi1 == vi0 || vi2 == vi0);

			const TVertex_handle &vi = (vi0 == vi2) ? vi1 : vi2;
			assert(vi != DT::Vertex_handle());
			assert(vi->info() >= 0);
			assert(vi->info() < (int)site_mirror_map.size()); 

			const int j = site_map(vi->info());
			assert(j >= 0);
			assert(j <  nimport);

			cell_list[j] = Cell();

			// now compute new volumes of each of the neighbours

			static std::vector<TEdge> edges;
			edges.clear();
			T1.incident_edges(vi, std::back_inserter(edges));

			TREAL volume_sj = 0.0;
			for (std::vector<TEdge>::iterator edge_it = edges.begin(); edge_it != edges.end(); edge_it++)
			{
				const TCell_circulator cc_end = T1.incident_cells(*edge_it);
				TCell_circulator cc(cc_end);

				const TVertex_handle &vi1 = edge_it->get<0>()->vertex(edge_it->get<1>());
				const TVertex_handle &vi2 = edge_it->get<0>()->vertex(edge_it->get<2>());
				assert(vi1 == vi || vi2 == vi);

				const TVertex_handle &vj = (vi == vi2) ? vi1 : vi2;
				assert(vj != DT::Vertex_handle());
				assert(vj->info() >= 0);
				assert(vj->info() < (int)site_mirror_map.size()); 

				const int j1 = site_map(vj->info());
				assert(j1 >= 0);
				assert(j1 <  nimport);

				cell_list[j].ngb.push_back(-1-j1);     // negative, means it is ngb, else face index

				static std::vector<TVEC3> vertex_list;
				vertex_list.clear();
				TVEC3 c(0.0);
				do
				{
					assert(!T1.is_infinite(cc));
					const TPoint tc = T1.dual(cc);
					const TVEC3 centre = TVEC3(tc.x(), tc.y(), tc.z());
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
			ngb_volume_list.push_back(std::make_pair(
						vi->info(), 
						std::make_pair(cell_list[j].Volume, volume_sj)
						));
			const TREAL dv = volume_sj - cell_list[j].Volume;
			assert(dv >= -SMALLDIFF*cell_list[j].Volume);
			tot_vol_check += dv;
			cell_list[j].Volume = volume_sj;
		}

		if (!(std::abs(cell_list[id].Volume - tot_vol_check) <= SMALLDIFF*cell_list[id].Volume))
		{
			const real si_vol = cell_list[id].Volume;
			const real si_vol_check = tot_vol_check;
			fprintf(stderr, " si_vol= %g   si_vol_check= %g diff= %g [ %g ] \n",
					si_vol, si_vol_check, si_vol - si_vol_check, 
					(si_vol - si_vol_check)/si_vol);
			assert(std::abs(si_vol - si_vol_check) <= SMALLDIFF*si_vol);
		}

		return true;
	}
#endif



}

// #include "build_mesh.cpp"

