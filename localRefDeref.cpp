#include "fvmhd3d.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h> 
#include <CGAL/Spatial_sort_traits_adapter_3.h>
#include <boost/iterator/counting_iterator.hpp>
#include <CGAL/spatial_sort.h>


namespace fvmhd3d
{

#if 0
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

  void System::localRefDeref(CkCallbkac &cb)
  {
    MainCB = cb;
    
    std::vector<int> active = active_list;
    active_list.clear();
    const int nactive = active.size();     // this time active_list means ref/deref
    active_list.clear();
    for (int i = 0; i < nactive; i++)
    {
      const int id = 
      ptcl_
    }
    
    ptcl_act.clear();
    ptcl_import_ptr = new std::vector<Particle>();

    const int ngbPass = 2;
    systemProxy[thisIndex].localMesh_import(ngbPass, CkCallback(CkIndex_System::IterateIII(), systemProxy[thisIndex]));
  }
  
  void System::localRefDeref_build()
  {
    assert(T_ptr == NULL);
    DT  *To = new T();
    T_ptr = (void*)To;
    DT &T = *(DT*)T_ptr;

    assert(Tvxt_ptr == 0.0);
    std::vector<TVertex_handle> *Tvtx_o = new std:vector<TVertex_handler>(ptclact.size());
    Tvtx_ptr = (void*)Tvtx_o;
    std::vector<TVertex_handle> &Tvtx_list = *(std::vector<TVertex_handle>*)Tvtx_ptr;

    assert(n_in_DTloc == 0);

    std::vector<TPoint> Tpoints;
    std::vector< int  > Tpoints_idx;
    Tpoints    .resize(ptcl_act.size());
    Tpoints_idx.resize(ptcl_act.size());

    cell_list.clear();
    cell_list.resize(ptcl_act.size());
    for (int i = 0; i < (const int)ptcl_act.size(); i++)
    {
      const vec3 &pos = ptcl_act[i]->get_pos();
      Tpoints    [i] = TPoint(pos.x, pos.y, pos.z);
      Tpoints_idx[i] = i;
      assert(cell_list[i].ngb.empty());
    }

    std::vector<std::ptrdiff_t> indices;
    indices.reserve(Tpoints.size());
    std::copy(
        boost::counting_iterator<std::ptrdiff_t>(0),
        boost::counting_iterator<std::ptrdiff_t>(Tpoints.size()),
        std::back_inserter(indices));
    CGAL::spatial_sort(indices.begin(), indices.end(), Search_traits_3(&(Tpoints[0])), CGAL::Hilbert_sort_median_policy());

    Tvtx_list.resize(ptcl_act.size());

    DT::Vertex_handle hint1;
    for (std::vector<std::ptrdiff_t>::iterator it = indices.begin(); it != indices.end(); it++)
    {
      hint1         = T.insert(Tpoints[*it], hint1);
      assert(hint1 != TVertex_handle());
      hint1->info() = Tpoints_idx[*it];
      assert(hint1->info() >= 0);
      assert(hint1->info() <  (int)ptcl_act.size());
      Tvtx_list[hint1->info()] = hint1;
      n_in_DT++;
    }
    assert(n_in_DT == (int)ptcl_act.size());
  }

  struct dVolumeS
  {
    int  Id;
    real oldVol;
    real newVol;
  };

  int System::localRefDeref_deref()
  {
  }

  bool System::localRefDeref_remove(
      const int vtx_id,
      std::vector< dVolumeS > &ngb_volumes)
  {
    assert(T_ptr == NULL);
    DT  *To = new T();
    T_ptr = (void*)To;
    DT &T = *(DT*)T_ptr;

    assert(Tvxt_ptr == 0.0);
    std::vector<TVertex_handle> *Tvtx_o = new std:vector<TVertex_handler>(ptclact.size());
    Tvtx_ptr = (void*)Tvtx_o;
    std::vector<TVertex_handle> &Tvtx_list = *(std::vector<TVertex_handle>*)Tvtx_ptr;

    const TVertex_handle &vi = Tvtx_list[vtx_id];
    assert(vi != DT::Vertex_handle());
    assert(vi->info() > 0.0);
    assert(vi->info() < (int)active_list.size());
    assert(vi->info() == vtx_id);

    std::vector<TEdge>        Tedges;
    std::vector<TCell_handle> Tcells;
    std::vector<bool> sites_ngb_used(n_in_DTloc, false);
    std::vector<int>  site_ngb_list;

    const vec3 &ipos = ptcl_act[i]->get_pos();

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
        assert(id < n_in_DTloc);
        if (sites_ngb_used[id]) continue;

        iadd++;
        sites_ngb_used[id] = true;
        site_ngb_list.push_back(v->info());
        if (id >= active_list.size())
        {
          success_flag = false;
        }
        Tedges.push_back(TEdge(ci, idx, iv));
      }
      assert(iadd < 4);
    }

    const int nngb = site_ngb_list.size();
    for (int j = 0; j < nngb; j++)
      sites_ngb_used[site_ngb_list[j]] = false;

    if (!success_flag) return false;

    // this cell is suitable for derefinement, so proceed as planned

    ////////////// COMPUTE OLD VOLUMES of ngb cells
    //

    static std::vector<real> oldVol
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
          if (id >= active_list.size()) 
          {
            success_flag = false;
          }
          site_ngb_list1.push_back(id);
          edges.push_back(TEdge(ci, idx, iv));
        }
        assert(iadd < 4);
      }

      const int nngb1	= site_ngb_list1.size();	
      for (int j1 = 0; j1 < nngb1; j1++)
        site_ngb_used[site_ngb_list1[j1]] = false;

      if (!success_flag) return false;

      int nj = 0;
      real volume_sj = 0.0;
      for (std::vector<TEdge>::iterator edge_it = edges.begin(); edge_it != edges.end(); edge_it++)
      {
        const TCell_circulator cc_end = T.incident_cells(*edge_it);
        TCell_circulator cc(cc_end);

        const int iv1 = edge_it->get<0>()->vertex(edge_it->get<1>())->info();
        const int iv2 = edge_it->get<0>()->vertex(edge_it->get<2>())->info();
        const int nj_id = (iv1 == vi->info()) ? iv2 : iv1;
        assert(nj_id == site_ngb_list1[nj++]);

        static std::vector<vec3> vertex_list;
        vertex_list.clear();
        vec3 c(0.0);
        do
        {
          assert(!T.is_infinite(cc));
          const TPoint tc = T.dual(cc);
#ifndef _EXACT_CONSTRUCTIONS_
          const vec3 centre = vec3(tc.x(), tc.y(), tc.z());
#else
          const vec3 centre = vec3(to_double(tc.x()), to_double(tc.y()), to_double(tc.z()));
#endif
          vertex_list.push_back(centre);
          c += centre;

          cc++;
        } while (cc != cc_end);

        const int nvtx = vertex_list.size();
        c *= 1.0/(real)nvtx;

        vec3 normal(0.0);
        vec3 v1 = vertex_list.back() - c;
        for (int j = 0; j < nvtx; j++)
        {
          const vec3 v2 = vertex_list[j] - c;
          const vec3 norm3 = v1.cross(v2);
          normal += norm3;
          v1 = v2;
        }
        real area = normal.abs();
        if (area == 0.0) continue;

        const TPoint pi = vi->point();
        const vec3 posj(pi.x(), pi.y(), pi.z());
        const real volume = std::abs(normal * (posj - c));
        volume_sj += volume;
      }	
      volume_sj *= 1.0/6.0;
      volume_old[j] = volume_sj;
    }

    ////////////// SANITY CHECK
    //

#if 1     // SANITY CHECK ...  \sum_j (Vnew_j - Vold_j) = V
#define _SANITY_CHECK_ENABLED_
    real si_vol = 0.0;
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

      static std::vector<vec3> vertex_list;
      vertex_list.reserve(32);


      for (std::vector<TEdge>::iterator edge_it = edges.begin(); edge_it != edges.end(); edge_it++)
      {
        const TCell_circulator cc_end = T.incident_cells(*edge_it);
        TCell_circulator cc(cc_end);

        vertex_list.clear();

        vec3 c(0.0);
        do
        {
          assert(!T.is_infinite(cc));

          const TPoint tc = T.dual(cc);
          const vec3 centre = vec3(tc.x(), tc.y(), tc.z());
          vertex_list.push_back(centre);
          c += centre;

          cc++;
        } while (cc != cc_end);

        const int nvtx = vertex_list.size();
        c *= 1.0/(real)nvtx;

        real area = 0.0;
        vec3 normal(0.0);
        vec3 v1 = vertex_list.back() - c;
        for (int j = 0; j < nvtx; j++)
        {
          const vec3 v2 = vertex_list[j] - c;
          const vec3 norm3 = v1.cross(v2);
          const real area3 = norm3.abs();
          area   += area3;
          normal += norm3;
          v1 = v2;
        }
        if (area == 0.0) continue;
        assert(std::abs(normal.abs() -  area) <= SMALLDIFF*area);

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

    real si_vol_check = 0.0;
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
      real volume_sj = 0.0;
      for (std::vector<TEdge>::iterator edge_it = edges.begin(); edge_it != edges.end(); edge_it++)
      {
        const TCell_circulator cc_end = T.incident_cells(*edge_it);
        TCell_circulator cc(cc_end);

        const int iv1 = edge_it->get<0>()->vertex(edge_it->get<1>())->info();
        const int iv2 = edge_it->get<0>()->vertex(edge_it->get<2>())->info();
        const int nj_id = (iv1 == vi->info()) ? iv2 : iv1;
        assert(nj_id == site_ngb_list1[nj++]);

        static std::vector<vec3> vertex_list;
        vertex_list.clear();
        vec3 c(0.0);
        do
        {
          assert(!T.is_infinite(cc));
          const TPoint tc = T.dual(cc);
          const vec3 centre = vec3(tc.x(), tc.y(), tc.z());
          vertex_list.push_back(centre);
          c += centre;

          cc++;
        } while (cc != cc_end);

        const int nvtx = vertex_list.size();
        c *= 1.0/(real)nvtx;

        vec3 normal(0.0);
        vec3 v1 = vertex_list.back() - c;
        for (int j = 0; j < nvtx; j++)
        {
          const vec3 v2 = vertex_list[j] - c;
          const vec3 norm3 = v1.cross(v2);
          normal += norm3;
          v1 = v2;
        }
        real area = normal.abs();
        if (area == 0.0) continue;

        const TPoint pi = vi->point();
        const vec3 posj(pi.x(), pi.y(), pi.z());
        const real volume = std::abs(normal * (posj - c));
        volume_sj += volume;
      }	
      volume_sj *= 1.0/6.0;
      ngb_volume_list.push_back(std::make_pair(
            vi->info(), 
            std::make_pair((real)volume_old[j], (real)volume_sj)
            ));
      const real dv = volume_sj - volume_old[j];
      assert(dv >= -SMALLDIFF*volume_old[j]);
      si_vol_check += dv;
    }

#ifdef _SANITY_CHECK_ENABLED_
    if (!(std::abs(si_vol - si_vol_check) <= SMALLDIFF*si_vol))
    {
      fprintf(stderr, " si_vol= %g   si_vol_check= %g diff= %g [ %g ] \n",
          si_vol, si_vol_check, si_vol - si_vol_check, 
          (si_vol - si_vol_check)/si_vol);
    }
    assert(std::abs(si_vol - si_vol_check) <= SMALLDIFF*si_vol);
#endif

    return true;

  } /* end System::localRefDeref_remove(..) */

  void System::localRefDeref_compute_total_volume(CkCallback &cb)
  {
    double volume_loc = 0;
    for (int i = 0; i < (const int)active_list.size(); i++)
      volume_loc += cell_list[i].Volume;

    contribute(sizeof(double), &volume_loc, CkReduction::sum_double, cb);
  }
#endif

  void System::ImportPtclAct()
  {
    const int nremote     = nimport_glb;
    const int iremote_end = ptcl_act.size();
    const int iremote_beg = iremote_end - nremote;

    std::vector< std::pair<int, int> > request_list;
    request_list.reserve(nremote);
    for (int i = iremote_beg; i < iremote_end; i++)
      if (ptcl_act[i]->is_active())
      {
        assert(ptcl_act[i]->chare() != thisIndex);
        request_list.push_back(std::make_pair(ptcl_act[i]->chare(), i));
      }

    std::sort(request_list.begin(), request_list.end(), std_pair_first_sort());

    const int nrequest = request_list.size();
    CkVec< pair<int, int> > sites2request;
    sites2request.reserve(nrequest);
    request_list.push_back(std::make_pair(-1,-1));

    for (int i = 0; i < nrequest; i++)
    {
      const int iElement = request_list[i].first;
      const int iId      = request_list[i].second;
      sites2request.push_back(std::make_pair(ptcl_act[iId]->id(), iId));
      assert(iElement >= 0);
      assert(iElement < numElements);
      assert(iElement != thisIndex);
      if (iElement != request_list[i+1].first && sites2request.size() > 0)
      {
        ImportPtclAct_nRequested++;
        systemProxy[iElement].ImportPtclAct_request(sites2request, thisIndex);
        sites2request.clear();
      }
    }

    if (ImportPtclAct_nRequested == 0)
      ImportPtclActCb.send();
  }

  void System::ImportPtclAct_request(const CkVec< pair<int,int> > &reqData, const int recvIndex)
  {
    assert(thisIndex != recvIndex);
    const int nrecv = reqData.size();
    assert(nrecv > 0);
    CkVec< pair<int, Particle> > data2send;
    data2send.reserve(nrecv);
    for (int i = 0; i < nrecv; i++)
    {
      const int  local_id = reqData[i].first;
      const int remote_id = reqData[i].second;
      assert(local_id >= 0);
      assert(local_id < local_n);
      assert(ptcl_list[local_id].is_active());
      if (ptcl_list[i].is_refine() || ptcl_list[i].is_derefine())
        data2send.push_back(std::make_pair(remote_id, ptcl_list[local_id]));
    }

    systemProxy[recvIndex].ImportPtclAct_recv(data2send);
  }

  void System::ImportPtclAct_recv(const CkVec< pair<int, Particle> > &recvUpdates)
  {
    const int nrecv = recvUpdates.size();
    for (int i = 0; i < nrecv; i++)
    {
      const int   iId      = recvUpdates[i].first;
      const Particle &pi   = recvUpdates[i].second;
      assert(iId >= (int)(nactive_loc + nimport_loc));
      assert(iId <  (int) ptcl_act.size());
      assert(ptcl_act[iId]->is_active());
      assert(ptcl_act[iId]->id() == pi.id());
      assert(ptcl_act[iId]->chare() == pi.chare());
      const vec3 pos = ptcl_act[iId]->get_pos();
      ptcl_act[iId]->set_volume(pi.get_volume());
      ptcl_act[iId]->set_status(pi.get_status());
      ptcl_act[iId]->set_rmax(pi.get_rmax());
    }
    ImportPtclAct_nRequested--;
    assert(ImportPtclAct_nRequested >= 0);

    if (ImportPtclAct_nRequested == 0)
      ImportPtclActCb.send();
  }


}

