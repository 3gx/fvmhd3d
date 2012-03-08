#include "fvmhd3d.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h> 
#include <CGAL/Spatial_sort_traits_adapter_3.h>
#include <boost/iterator/counting_iterator.hpp>
#include <CGAL/spatial_sort.h>


namespace fvmhd3d
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K; 
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
  typedef CGAL::Spatial_sort_traits_adapter_3<K, TPoint*> Search_traits_3;
  
  void System::globalMesh_build(const bool relax, CkCallback &cb)
  {
    globalMesh_completeCb = cb;
    globalMesh_relax_flag = relax;

#if 0
    assert(globalMesh_T_ptr == NULL);
    assert(globalMesh_Tvtx_list_ptr = NULL);
#endif

    globalMesh_T_ptr         = (void*)new DT;
    globalMesh_Tvtx_list_ptr = (void*)(new std::vector<TVertex_handle>(local_n));
   
    cell_list.resize(local_n); 
		assert((int)ngb_list.size() == local_n);

    assert(globalMesh_failed_ptcl.empty());
  
    assert(nSend_cntr == 0);

    DT &T = *(DT*)globalMesh_T_ptr;

    std::vector<TVertex_handle> &Tvtx_list = *(std::vector<TVertex_handle>*)globalMesh_Tvtx_list_ptr;

//		assert((int)ptcl_list.size() == local_n);

    {
      const int np = local_n;
      std::vector<TPoint> Tpoints(np);
      std::vector< int  > Tpoints_idx(np);

      for (int i = 0; i < np; i++)
      {
        const vec3 &pos = ptcl_list[i].get_pos();
        Tpoints    [i] = TPoint(pos.x, pos.y, pos.z);
        Tpoints_idx[i] = i;
      }

      std::vector<std::ptrdiff_t> indices;
      indices.reserve(np);
      std::copy(
          boost::counting_iterator<std::ptrdiff_t>(0),
          boost::counting_iterator<std::ptrdiff_t>(np),
          std::back_inserter(indices));
      CGAL::spatial_sort(indices.begin(), indices.end(), Search_traits_3(&(Tpoints[0])), CGAL::Hilbert_sort_median_policy());

      DT::Vertex_handle hint1;
      for (std::vector<std::ptrdiff_t>::iterator it = indices.begin(); it != indices.end(); it++)
      {
        hint1         = T.insert(Tpoints[*it], hint1);
        assert(hint1 != TVertex_handle());
        hint1->info() = Tpoints_idx[*it];
        assert(hint1->info() >= 0);
        assert(hint1->info() <  np);
        Tvtx_list[hint1->info()] = hint1;
      }


      std::vector<TEdge> Tedges;
      std::vector<TCell_handle> Tcells;
      std::vector<bool> sites_ngb_used(local_n, false);
      std::vector<int> site_ngb_list;

      for (int i = 0; i < np; i++)
      {
        const TVertex_handle &vi = Tvtx_list[i];
        assert(vi != TVertex_handle());
        assert(vi->info() == i);

        const vec3 &ipos = ptcl_list[i].get_pos();

        Tedges.clear();
        Tcells.clear();
        site_ngb_list.clear();
        T.incident_cells(vi, std::back_inserter(Tcells));

        const int ncells = Tcells.size();
        for (int icell = 0; icell < ncells; icell++)
        {
          const TCell_handle &ci = Tcells[icell];

          int idx = -1;
          for (int iv = 0; iv < 4; iv++)
            if (ci->vertex(iv) == vi)
              idx = iv;

          int iadd = 0;
          for (int iv = 0; iv < 4; iv++)
          {
            if (iv == idx) continue;

            const TVertex_handle &v = ci->vertex(iv);
            if (T.is_infinite(v)) 
              continue;
            assert (v != TVertex_handle());

            const int id = v->info();
            assert(id >= 0);
            assert(id < (int)local_n);
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

        real r2max = (Tedges.empty()) ? HUGE : 0.0;

        const int NMAXEDGE = 1024;
        static std::vector<vec3> vertex_list[NMAXEDGE];

        int nedge = 0;
        // for each Delaunay edge trace it dual Voronoi face
        for (std::vector<TEdge>::iterator edge_it = Tedges.begin(); edge_it != Tedges.end(); edge_it++)
        {
          const TCell_circulator cc_end = T.incident_cells(*edge_it);
          TCell_circulator cc(cc_end);
          vertex_list[nedge].clear();

          // compute face vertices
          do
          {
            if (T.is_infinite(cc))
            {
              r2max = HUGE;
            }
            else
            {
              const TPoint c = T.dual(cc);
              const vec3 centre(c.x(), c.y(), c.z());
              r2max = std::max(r2max, (centre - ipos).norm2());

              vertex_list[nedge].push_back(centre);
            }

            cc++;
          } while (cc != cc_end);
          nedge++;
          assert(nedge < NMAXEDGE);
        }  // for edge < nedge

#if 0
        const vec3 domain_hsize = global_domain_size/2;
#else
        const vec3 domain_hsize = domains->proc_domains[thisIndex].hsize() * 2.0;
#endif
        real lmax = std::min(std::min(domain_hsize.x, domain_hsize.y), domain_hsize.z);
        if (ptcl_list[i].get_rmax() < 1e9)
          lmax = 2.0*ptcl_list[i].get_rmax();

        real rmax = lmax;
        if (r2max < HUGE)
          rmax = 2.01*std::sqrt(r2max);

        rmax = std::min(rmax, lmax);

        ptcl_list[i].set_rmax(rmax);

        const boundary bi(ipos - rmax, ipos + rmax);
        if (!domains->proc_domains[thisIndex].isinbox(bi)) 
        {
          globalMesh_failed_ptcl.push_back(std::make_pair(i, rmax));
          continue;
        }

        ngb_list[i].clear();

        real cell_volume   = 0.0;
        vec3 cell_centroid = 0.0;
        for (int edge = 0; edge < nedge; edge++)
        {
          const int nvtx = vertex_list[edge].size();

          vec3 c = 0.0;
          for (int j = 0; j < nvtx; j++)
            c += vertex_list[edge][j];
          c *= 1.0/(real)nvtx;

          vec3 norm(0.0);
          const vec3 &centroid = c;
          vec3 v1 = vertex_list[edge].back() - c;
          real area1 = 0.0;
          for (int j = 0; j < nvtx; j++)
          {
            const vec3 v2 = vertex_list[edge][j] - c;
            const vec3 norm3 = v1.cross(v2);
            const real area3 = norm3.abs();
            norm  +=      norm3;
            area1 +=      area3;
            v1 = v2;
          }

          const real SMALLDIFF1 = 1.0e-10;
          const real area0 = area1; // face.n.abs();
          const real L1 = std::sqrt(area0);
          const real L2 = (centroid - ipos).abs();
          const real area = (L1 < SMALLDIFF1*L2) ? 0.0 : area0;

          if (area > 0.0 && norm.abs() > 0.0)
            ngb_list[i].push_back(std::make_pair(thisIndex, site_ngb_list[edge]));
          else
            continue;

          const vec3 cv  = ipos - centroid;
          v1    = vertex_list[edge].back() - centroid;
          const real fourth = 1.0/4.0;
          for (int j = 0; j < nvtx; j++)
          {
            const vec3 v2   = vertex_list[edge][j] - centroid;
            const vec3 c4   = centroid + (v1 + v2 + cv) * fourth;
            const real vol4 = std::abs(v1.cross(v2) * cv);
            cell_volume    +=      vol4;
            cell_centroid  += c4 * vol4;
            v1 = v2;
          }
        }  // for edge < nedge
        cell_centroid *= 1.0/cell_volume;
        cell_volume   *= 1.0/6.0;

        cell_list[i].centroid = cell_centroid;
        cell_list[i].Volume   = cell_volume;

      }
    }

    {
      std::vector< std::pair<int, int> > export_list;
      std::vector<int> remote_tiles;
      const int nfailed = globalMesh_failed_ptcl.size();
      for (int i = 0; i < nfailed; i++)
      {
        const Particle &p = ptcl_list[globalMesh_failed_ptcl[i].first];
        const real   rmax =           globalMesh_failed_ptcl[i].second;
        remote_tiles.clear();	
        domains->proc_tree.root.walk_boundary(boundary(p.get_pos()-rmax, p.get_pos()+rmax), remote_tiles, global_domain_size);
        assert(remote_tiles.size() > 0);
        const int ntiles = remote_tiles.size();
        for (int j = 0; j < ntiles; j++)
          if (remote_tiles[j] != thisIndex)
            export_list.push_back(std::make_pair(remote_tiles[j], p.id()));
      }

      std::sort(export_list.begin(), export_list.end(), std_pair_first_sort());

      const int nexport = export_list.size();
      CkVec<Particle> export_ptcl;
      export_list.push_back(std::make_pair(-1,-1));

      for (int i = 0; i < nexport; i++)
      {
        const int iElement = export_list[i].first;
        const int iId      = export_list[i].second;
        export_ptcl.push_back(ptcl_list[iId]);
        assert(iElement >= 0);
        assert(iElement < numElements);
        assert(iElement != thisIndex);
        if (iElement != export_list[i+1].first && export_ptcl.size() > 0)
        {
          nSend_cntr++;
          systemProxy[iElement].globalMesh_recvPtcl(export_ptcl, thisIndex);
          export_ptcl.clear();
        }
      }
    }

    nSend_cntr++;
    globalMesh_recvPtcl_ticket();
  }
  
  void System::globalMesh_recvPtcl_ticket()
  {
    nSend_cntr--;
    assert(nSend_cntr >= 0);
    if (nSend_cntr == 0)
      contribute(CkCallback(CkIndex_System::globalMesh_complete(), thisProxy));
  }
  
  void System::globalMesh_recvPtcl(const CkVec<Particle> &ptcl_in, const int recvIndex)
  {
    assert(recvIndex != thisIndex);
    const int n     = ptcl_in  .size();

    const vec3 box_centre = domains->proc_domains[thisIndex].centre();
    const vec3 box_hsize  = domains->proc_domains[thisIndex].hsize();
    const vec3 c = global_domain.centre();
    const vec3 s = global_domain_size;

    vec3 ppos[8]; 
    for (int i = 0; i < n; i++)
    {
      int jmin = 0;
      ppos[0] = ptcl_in[i].get_pos();
      vec3 dr = (ppos[jmin] - box_centre).abseach() - box_hsize;
      dr = (dr + dr.abseach())*0.5;
      vec3 dr_min = dr;

      const bool flagx = ppos[0].x < c.x;
      const bool flagy = ppos[0].y < c.y;
      const bool flagz = ppos[0].z < c.z;
      for (int j = 1; j < 8; j++)
      {
        ppos[j] = ppos[0];
        vec3 &p = ppos[j];
        p.x = ((j&1) == 1) ? (flagx ? p.x + s.x : p.x - s.x) : p.x;
        p.y = ((j&2) == 2) ? (flagy ? p.y + s.y : p.y - s.y) : p.y;
        p.z = ((j&4) == 4) ? (flagz ? p.z + s.z : p.z - s.z) : p.z;
        dr = (p - box_centre).abseach() - box_hsize;
        dr = (dr + dr.abseach())*0.5;
        if (dr.norm2() < dr_min.norm2())
        {
          jmin = j;
          dr_min = dr;
        }
      }

      ptcl_list.push_back(ptcl_in[i]);
      ptcl_list.back() = ppos[jmin];
    }

    systemProxy[recvIndex].globalMesh_recvPtcl_ticket();
  }
  
  void System::globalMesh_complete()
  {
    DT &T = *(DT*)globalMesh_T_ptr;

    std::vector<TVertex_handle>   &Tvtx_list = *(std::vector<TVertex_handle>*)globalMesh_Tvtx_list_ptr;

		assert((int)cell_list.size() == local_n);

    { 
      const int np = ptcl_list.size() - local_n;
      std::vector<TPoint> Tpoints    (np);
      std::vector< int  > Tpoints_idx(np);

      for (int i = 0; i < np; i++)
      {
        const vec3 &pos = ptcl_list[local_n + i].get_pos();
        Tpoints    [i] = TPoint(pos.x, pos.y, pos.z);
        Tpoints_idx[i] = local_n + i;
      }

      std::vector<std::ptrdiff_t> indices;
      indices.reserve(np);
      std::copy(
          boost::counting_iterator<std::ptrdiff_t>(0),
          boost::counting_iterator<std::ptrdiff_t>(np),
          std::back_inserter(indices));
      CGAL::spatial_sort(indices.begin(), indices.end(), Search_traits_3(&(Tpoints[0])), CGAL::Hilbert_sort_median_policy());

      DT::Vertex_handle hint1;
      for (std::vector<std::ptrdiff_t>::iterator it = indices.begin(); it != indices.end(); it++)
      {
        hint1         = T.insert(Tpoints[*it], hint1);
        assert(hint1 != TVertex_handle());
        hint1->info() = Tpoints_idx[*it];
        assert(hint1->info() >= local_n);
        assert(hint1->info() <  local_n+np);
      }
    }

    std::vector<TEdge> Tedges;
    std::vector<TCell_handle> Tcells;
    std::vector<bool> sites_ngb_used(ptcl_list.size(), false);
    std::vector<int> site_ngb_list;

    const int np = globalMesh_failed_ptcl.size();
    for (int idx = 0; idx < np; idx++)
    {
      const int i = globalMesh_failed_ptcl[idx].first;

      const TVertex_handle &vi = Tvtx_list[i];
      assert(vi != TVertex_handle());
      assert(vi->info() == i);

      const vec3 &ipos = ptcl_list[i].get_pos();

      Tedges.clear();
      Tcells.clear();
      site_ngb_list.clear();
      T.incident_cells(vi, std::back_inserter(Tcells));

      const int ncells = Tcells.size();
      for (int icell = 0; icell < ncells; icell++)
      {
        const TCell_handle &ci = Tcells[icell];

        int idx = -1;
        for (int iv = 0; iv < 4; iv++)
          if (ci->vertex(iv) == vi)
            idx = iv;

        int iadd = 0;
        for (int iv = 0; iv < 4; iv++)
        {
          if (iv == idx) continue;

          const TVertex_handle &v = ci->vertex(iv);
          assert(!T.is_infinite(v));
          assert (v != TVertex_handle());

          const int id = v->info();
          assert(id >= 0);
          assert(id < (int)ptcl_list.size());
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

      real r2max = (Tedges.empty()) ? HUGE : 0.0;

      const int NMAXEDGE = 1024;
      static std::vector<vec3> vertex_list[NMAXEDGE];

      int nedge = 0;
      // for each Delaunay edge trace it dual Voronoi face
      for (std::vector<TEdge>::iterator edge_it = Tedges.begin(); edge_it != Tedges.end(); edge_it++)
      {
        const TCell_circulator cc_end = T.incident_cells(*edge_it);
        TCell_circulator cc(cc_end);
        vertex_list[nedge].clear();

        // compute face vertices
        do
        {
          if (T.is_infinite(cc))
          {
            assert(false);
          }
          else
          {
            const TPoint c = T.dual(cc);
            const vec3 centre(c.x(), c.y(), c.z());
            r2max = std::max(r2max, (centre - ipos).norm2());

            vertex_list[nedge].push_back(centre);
          }

          cc++;
        } while (cc != cc_end);
        nedge++;
        assert(nedge < NMAXEDGE);
      }  // for edge < nedge

      const real rmax = 2.01*std::sqrt(r2max);
      ptcl_list[i].set_rmax(rmax);
      assert(ptcl_list[i].get_rmax() < 1e9);

      ngb_list[i].clear();

      real cell_volume   = 0.0;
      vec3 cell_centroid = 0.0;
      for (int edge = 0; edge < nedge; edge++)
      {
        const int nvtx = vertex_list[edge].size();

        vec3 c = 0.0;
        for (int j = 0; j < nvtx; j++)
          c += vertex_list[edge][j];
        c *= 1.0/(real)nvtx;

        vec3 norm(0.0);
        const vec3 &centroid = c;
        vec3 v1 = vertex_list[edge].back() - c;
        real area1 = 0.0;
        for (int j = 0; j < nvtx; j++)
        {
          const vec3 v2 = vertex_list[edge][j] - c;
          const vec3 norm3 = v1.cross(v2);
          const real area3 = norm3.abs();
          norm  +=      norm3;
          area1 +=      area3;
          v1 = v2;
        }

        const real SMALLDIFF1 = 1.0e-10;
        const real area0 = area1; // face.n.abs();
        const real L1 = std::sqrt(area0);
        const real L2 = (centroid - ipos).abs();
        const real area = (L1 < SMALLDIFF1*L2) ? 0.0 : area0;

        if (area > 0.0 && norm.abs() > 0.0)
          ngb_list[i].push_back(ptcl_list[site_ngb_list[edge]].chare_id());
        else
          continue;

        const vec3 cv  = ipos - centroid;
        v1    = vertex_list[edge].back() - centroid;
        const real fourth = 1.0/4.0;
        for (int j = 0; j < nvtx; j++)
        {
          const vec3 v2   = vertex_list[edge][j] - centroid;
          const vec3 c4   = centroid + (v1 + v2 + cv) * fourth;
          const real vol4 = std::abs(v1.cross(v2) * cv);
          cell_volume    +=      vol4;
          cell_centroid  += c4 * vol4;
          v1 = v2;
        }
      }  // for edge < nedge
      cell_centroid *= 1.0/cell_volume;
      cell_volume   *= 1.0/6.0;

      cell_list[i].centroid = cell_centroid;
      cell_list[i].Volume   = cell_volume;

    }

    ptcl_list.resize(local_n);

    double volume_loc = 0.0;
    for (int i = 0; i < local_n; i++)
      volume_loc  += cell_list[i].Volume;

    if (globalMesh_relax_flag)
      for (int i = 0; i < local_n; i++)
        ptcl_list[i] = cell_list[i].centroid;


    delete (DT*)globalMesh_T_ptr;
    delete (std::vector<TVertex_handle>*)globalMesh_Tvtx_list_ptr;

    globalMesh_T_ptr = NULL;
    globalMesh_Tvtx_list_ptr = NULL;
    clear_vec(cell_list);
    clear_vec(globalMesh_failed_ptcl);

    contribute(sizeof(double), &volume_loc, CkReduction::sum_double, globalMesh_completeCb);
  }

}
