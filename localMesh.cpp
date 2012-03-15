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

  void System::localMesh_import_pass(std::vector<int> sites2process, const int ngbPass)
  {
    const int n2process = sites2process.size();
    assert(nSend_cntr == 0);

    int isite = 0;
    while (isite < n2process)
    {
      const int nsend     = sites2process[isite++];
      const int sendIndex = sites2process[isite++];
      assert(sendIndex >= 0);
      assert(sendIndex <  numElements);

      CkVec<Particle> return_list;
      CkVec<  int   > site_in;
      return_list.reserve(nsend);
      site_in    .reserve(nsend);
      assert(isite + nsend <= n2process);
      for (int i = 0; i < nsend; i++)
      {
        const int iId = sites2process[isite++];
        assert(iId >= 0);
        assert(iId <  local_n);
        if (ngb_list_hash[iId].use(sendIndex))
        {
          ngb_list_hash_used.push_back(iId);
          site_in.push_back(iId);

          Problem_predict_meshpoint_position(iId);
          ptcl_list[iId] = mesh_pnts[iId].pos;
          slopeLimiter_all[iId] = 1.0;

          if (!(sendIndex == thisIndex && ptcl_list[iId].is_active()))
            return_list.push_back(ptcl_list[iId]);
        }
      }

      const int nrecv = site_in.size();

      if (nrecv > 0)
      {
        if (sendIndex != thisIndex)
        {
          nRecv_cntr++;
          systemProxy[sendIndex].localMesh_insertPtcl(return_list, thisIndex);
        }
        else
          localMesh_insertPtcl(return_list, thisIndex);
      }

      if (nrecv > 0 && localMesh_ngbPass + 1 < ngbPass)
      {
        std::vector< pair<int, int> > export_list;
        export_list.reserve(nrecv);

        std::map< std::pair<int, int>, bool> remote_map;

        for (int irecv = 0; irecv < nrecv; irecv++)
        {
          const int i = site_in[irecv];
          assert(i >= 0);
          assert(i < local_n);
          const Neighbours< pair<int, int> > &ngb = ngb_list[i];
          const int nj = ngb.size();
          for (int j = 0; j < nj; j++)
          {
            const int jIndex = ngb[j].first;
            const int jId    = ngb[j].second;
            if (jIndex == thisIndex)
            {
              assert(jId >= 0);
              assert(jId <  local_n);
              if (!ngb_list_hash[jId].is_used(sendIndex))
                export_list.push_back(ngb[j]);
            }
            else
            {
              const int size0 = remote_map.size();
              remote_map[ngb[j].make_pair()] = true;
              if ((int)remote_map.size() == size0+1)
                export_list.push_back(ngb[j]);
              else
                assert((int)remote_map.size() == size0);
            }
          }
        }

        std::sort(export_list.begin(), export_list.end(), std_pair_first_sort());

        const int nexport = export_list.size();
        CkVec<int> sites2export;
        sites2export.reserve(nexport);
        export_list.push_back(std::make_pair(-1, -1));
        for (int i = 0; i < nexport; i++)
        {
          const int iElement =   export_list[i].first;
          assert(export_list[i].second >= 0);
          sites2export.push_back(export_list[i].second);
          assert(iElement >= 0);
          assert(iElement < numElements);
          if (iElement != export_list[i+1].first && sites2export.size() > 0)
          {
            if (thisIndex != iElement)
            {
              nSend_cntr++;
              thisProxy[iElement].localMesh_import_new(sites2export, std::make_pair(sendIndex, thisIndex));
            }
            else
              localMesh_import_new(sites2export, std::make_pair(sendIndex, thisIndex));
            sites2export.clear();
          }
        }
      }

    }  /* while isite < n2process */
  }

  void System::localMesh_import_new(const CkVec<int> &recvData, const pair<int, int> recvPair)
  {
    const int sendIndex = recvPair.first;
    const int recvIndex = recvPair.second;


    const int nsend = recvData.size();
    localMesh_sites2process.push_back(nsend);
    localMesh_sites2process.push_back(sendIndex);
    for (int i = 0; i < nsend; i++)
      localMesh_sites2process.push_back(recvData[i]);

    if (recvIndex != thisIndex)
      systemProxy[recvIndex].localMesh_import_new_ticket();
  };

  void System::localMesh_insertPtcl(const CkVec<Particle> &ptcl_in, const int recvIndex)
  {
    if (recvIndex != thisIndex)
    {
      const vec3 box_centre = domains->proc_domains[thisIndex].centre();
      const vec3 box_hsize  = domains->proc_domains[thisIndex].hsize();
      const vec3 c = global_domain.centre();
      const vec3 s = global_domain_size;

      vec3 ppos[8]; 
      const int nrecv    = ptcl_in.size();
      for (int i = 0; i < nrecv; i++)
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

        ptcl_import_ptr->push_back(ptcl_in[i]);
        ptcl_import_ptr->back() = ppos[jmin];
      }
      
      systemProxy[recvIndex].localMesh_insertPtcl_ticket();
    }
    else
    {
      const int nrecv = ptcl_in.size();
      for (int i = 0; i < nrecv; i++)
      {
        assert(ptcl_in[i].chare() == thisIndex);
        ptcl_act.push_back(&ptcl_list[ptcl_in[i].id()]);
      }
    }

  }; 

  void System::localMesh_import_complete()
  {
    nactive_loc = active_list.size();
    nimport_loc = ptcl_act.size();
    nimport_glb = ptcl_import_ptr->size();

    std::vector<Particle*> ptcl1 = ptcl_act;
    ptcl_act.clear();


    for (std::vector<int>::iterator it = active_list.begin(); it != active_list.end(); it++)
      ptcl_act.push_back(&ptcl_list[*it]);
    for (std::vector<Particle*>::iterator it = ptcl1.begin(); it != ptcl1.end(); it++)
      ptcl_act.push_back(*it);
    for (int i = 0; i < nimport_glb; i++)
      ptcl_act.push_back(&(*ptcl_import_ptr)[i]);

    /* set locally active particles as active */
    for (int i = 0; i < nactive_loc; i++)
      assert(ptcl_act[i]->is_active());

    for (int i = nactive_loc; i < nactive_loc + nimport_loc; i++)
      assert(!ptcl_act[i]->is_active());

    assert((int)ptcl_list.size() == local_n);

    for (int i = 0; i < (const int)ngb_list_hash_used.size(); i++)
      ngb_list_hash[ngb_list_hash_used[i]] = Hash<int>();
    ngb_list_hash_used.clear();

    localMesh_sites2process.clear();

#if 1
    {
      assert(mesh_import_ptr == NULL);
      assert(Wrec_import_ptr == NULL);
      assert(Wrec_minmax_import_ptr == NULL);
      assert(  dU_import_ptr == NULL);
      assert(slopeLimiter_import_ptr == NULL);

      mesh_act.resize(ptcl_act.size());
      Wrec_act.resize(ptcl_act.size());
      Wrec_minmax_act.resize(ptcl_act.size());
      dU_act  .resize(ptcl_act.size());
      U_act   .resize(ptcl_act.size());
      slopeLimiter_act.resize(ptcl_act.size());

      mesh_import_ptr         = new std::vector<MeshPoint>(nimport_glb);
      Wrec_minmax_import_ptr  = new std::vector<pair<Fluid_flt,Fluid_flt> >(nimport_glb);
      Wrec_import_ptr         = new std::vector<Fluid_rec>(nimport_glb + nimport_loc);
      dU_import_ptr           = new std::vector<FluidD   >(nimport_glb + nimport_loc, 0.0);
      slopeLimiter_import_ptr = new std::vector<Fluid_flt>(nimport_glb, 1.0);

      for (int i = 0; i < nactive_loc; i++)
      {
        const int Id = active_list[i];

        mesh_act        [i] = &mesh_pnts       [Id];
        Wrec_act        [i] = &Wrec_list       [Id];
        Wrec_minmax_act [i] = &Wrec_minmax_list[Id];
        dU_act          [i] = &  dU_list       [Id];
        U_act           [i] = &   U_list       [Id];
        slopeLimiter_act[i] = &slopeLimiter_all[Id];
      }

      for (int i = nactive_loc; i < nactive_loc + nimport_loc; i++)
      {
        const int Id = ptcl_act[i]->id();
        assert(Id >= 0);
        assert(Id <  local_n);
        assert(ptcl_act[i]->chare() == thisIndex);

        mesh_act        [i] = &mesh_pnts       [Id];
        slopeLimiter_act[i] = &slopeLimiter_all[Id];
        Wrec_minmax_act [i] = &Wrec_minmax_list[Id];
        U_act           [i] = &   U_list       [Id];

        (*Wrec_import_ptr)[i-nactive_loc] = Wrec_list[Id];
        Wrec_act[i] = &(*Wrec_import_ptr)[i-nactive_loc];
        dU_act  [i] = &(*  dU_import_ptr)[i-nactive_loc];
      }

      for (int i = 0; i < nimport_glb; i++)
      {
        const int j = i + ptcl_act.size() - nimport_glb;

        mesh_act        [j] = &(*        mesh_import_ptr)[i];
        Wrec_minmax_act [j] = &(* Wrec_minmax_import_ptr)[i];
        Wrec_act        [j] = &(*        Wrec_import_ptr)[i + nimport_loc];
        dU_act          [j] = &(*          dU_import_ptr)[i + nimport_loc];
        slopeLimiter_act[j] = &(*slopeLimiter_import_ptr)[i];
      }

      for (int i = 0; i < nactive_loc + nimport_loc; i++)
      {
        Wrec_act[i]->pos  = mesh_act[i]->pos;
        Wrec_act[i]->bnd  = mesh_act[i]->boundary;
      }
    }
#endif

#if 1
    localMesh_completeCb.send();
#else
    contribute(localMesh_completeCb);
#endif
  }

  void System::localMesh_destroy()
  {
    delete (DT*)T_ptr;
    delete (std::vector<TVertex_handle>*)Tvtx_ptr;
    T_ptr      = NULL;
    Tvtx_ptr   = NULL;
    n_in_DTloc = 0;
  }

  void System::localMesh_build()
  {
    cell_list.resize(active_list.size());

    assert(T_ptr == NULL);
    DT  *To = new DT();
    T_ptr = (void*)To;
    DT &T = *(DT*)T_ptr;

    assert(Tvtx_ptr == NULL);
    std::vector<TVertex_handle> *Tvtx_o = new std::vector<TVertex_handle>(ptcl_act.size());
    Tvtx_ptr = (void*)Tvtx_o;
    std::vector<TVertex_handle> &Tvtx_list = *(std::vector<TVertex_handle>*)Tvtx_ptr;

    if (active_list.size() == 0)
      assert((int)ptcl_list.size() == local_n);

    std::vector<TPoint> Tpoints;
    std::vector< int  > Tpoints_idx;
    Tpoints    .resize(ptcl_act.size());
    Tpoints_idx.resize(ptcl_act.size());

    assert(n_in_DTloc == 0);
    T.clear();
    face_list.clear();

    for (int i = 0; i < (const int)ptcl_act.size(); i++)
    {
      const vec3 &pos = ptcl_act[i]->get_pos();
      Tpoints    [i] = TPoint(pos.x, pos.y, pos.z);
      Tpoints_idx[i] = i;
      if (i < (int)active_list.size())
        cell_list[i].ngb.clear();
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
      n_in_DTloc++;
    }
    assert(n_in_DTloc == (int)ptcl_act.size());

    std::vector<TEdge> Tedges;
    std::vector<TCell_handle> Tcells;
    std::vector<bool> sites_ngb_used(n_in_DTloc, false);
    std::vector<int>  site_ngb_list;

    const int nactive = active_list.size();
    for (int i = 0; i < nactive; i++)
    {
      const TVertex_handle &vi = Tvtx_list[i];
      assert(vi != TVertex_handle());
      assert(vi->info() == i);

      const vec3 &ipos = ptcl_act[i]->get_pos();

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
          assert(id < n_in_DTloc);
          if (sites_ngb_used[id]) continue;

          iadd++;
          sites_ngb_used[id] = true;
          site_ngb_list.push_back(v->info());
          Tedges.push_back(TEdge(ci, idx, iv));
        }
        assert(iadd < 4);
      }

      const int nngb = site_ngb_list.size();
      for (int j = 0; j < nngb; j++)
        sites_ngb_used[site_ngb_list[j]] = false;

      assert(!Tedges.empty());
      real r2max = 0.0;

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
      ptcl_act[i]->set_rmax(rmax);

      real cell_volume   = 0.0;
      vec3 cell_centroid = 0.0;
      int edge = 0;
      for (std::vector<TEdge>::iterator edge_it = Tedges.begin(); edge_it != Tedges.end(); edge_it++, edge++)
      {
        assert(edge < nedge);
        const int nvtx = vertex_list[edge].size();

        const int i1 = edge_it->get<0>()->vertex(edge_it->get<1>())->info();
        const int i2 = edge_it->get<0>()->vertex(edge_it->get<2>())->info();

        Face face;
        face.s1 = (i1 == i) ? i1 : i2;
        face.s2 = (i1 == i) ? i2 : i1;
        assert(face.s1 == i);
        assert(face.s2 != i);

        int face_id = -1;
        const int nj = cell_list[i].ngb.size();
        for (int j = 0; j < nj; j++)
          if (face_list[cell_list[i].ngb[j]].ngb<false>(i) == face.s2)
          {
            face_id = cell_list[i].ngb[j];
            break;
          }

        if (face_id != -1)
        {
          assert(face_id >= 0);
          assert(face_id < (const int)face_list.size());
          face = face_list[face_id];
        }
        else
        {
          vec3 c = 0.0;
          for (int j = 0; j < nvtx; j++)
            c += vertex_list[edge][j];
          c *= 1.0/(real)nvtx;

          face.n        = 0.0;
          face.centroid = 0.0;
          real area1 = 0.0;
          const real third = 1.0/3.0;
          vec3 v1 = vertex_list[edge].back() - c;
          for (int j = 0; j < nvtx; j++)
          {
            const vec3 v2 = vertex_list[edge][j] - c;

            const vec3 norm3 = v1.cross(v2);
            const real area3 = norm3.abs();
            const vec3 c3    = c + (v1 + v2) * third;

            face.n        +=      norm3;
            face.centroid += c3 * area3;
            area1         +=      area3;

            v1 = v2;
          }

          const real SMALLDIFF1 = 1.0e-10;
          const real area0 = area1;
          const real L1 = std::sqrt(area0);
          const real L2 = (face.centroid - ipos).abs();
          const real area = (L1 < SMALLDIFF1*L2) ? 0.0 : area0;

          const real nabs = face.n.abs();
          if (area > 0.0 && nabs > 0.0)
          {
            const real iarea = 1.0/area;
            face.centroid *= iarea;
            face.n        *= 0.5;

            if ((face.centroid - ipos)*face.n < 0.0)
              face.n *= -1.0;

            const int jid = face.s2;
            assert(jid >= 0);
            assert(jid < (const int)ptcl_act.size());


            if (jid >= (const int)active_list.size())
              ptcl_act[jid]->set_ngb();

#if 0
            if (!site_import[jid].is_active())
            {
              site_import[jid].set_local_ngb();
              site_import[isite].set_hasngb();
            }
#endif
          }
          else
          {
            face.n = 0.0;
          }
        }

        if (face.area() == 0.0)  continue;

        if (face_id == -1)
        {
          face_list.push_back(face);
          cell_list[face.s1].ngb.push_back(face_list.size() - 1);
          if (face.s2 < (const int)active_list.size())
            cell_list[face.s2].ngb.push_back(face_list.size() - 1);
        }

        const vec3 cv  = ipos - face.centroid;
        vec3 v1 = vertex_list[edge].back() - face.centroid;
        const real fourth = 1.0/4.0;
        for (int j = 0; j < nvtx; j++)
        {
          const vec3 v2   = vertex_list[edge][j] - face.centroid;
          const vec3 c4   = face.centroid + (v1 + v2 + cv) * fourth;
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

  } /* end System::localMeshbuild(..) */

  void System::localMesh_compute_total_volume(CkCallback &cb)
  {
    double volume_loc = 0;
    for (int i = 0; i < (const int)active_list.size(); i++)
      volume_loc += cell_list[i].Volume;

    contribute(sizeof(double), &volume_loc, CkReduction::sum_double, cb);
  }
}

