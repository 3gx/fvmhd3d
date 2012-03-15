#include "fvmhd3d.h"

namespace fvmhd3d
{
#if 1
  void System::GetActiveList(CkCallback &cb)
  {
    MainCB = cb;
    int max_rung_loc = scheduler.get_max_rung();
    contribute(sizeof(int), &max_rung_loc, CkReduction::max_int, CkCallback(CkIndex_System::GetActiveListII(0), thisProxy));
  }
  void System::GetActiveListII(CkReductionMsg *msg)
  {
    scheduler.set_max_rung(*(int*)msg->getData());
    delete msg;

    active_list.clear();
    t_previous = scheduler.get_tsys();
    assert(t_previous >= 0.0);
    assert(t_previous == t_global);
    dt_global = scheduler.pull_active_list(active_list);
    t_global  = scheduler.get_tsys();

    for (std::vector<int>::iterator it = active_list.begin(); it != active_list.end(); it++)
    {
      ptcl_list[*it].set_active();
      if (!(mesh_pnts[*it].tend == t_global))
      {
        fprintf(stderr ," thisIndex= %d  tbeg= %g tend= %g  t_global= %g \n",
            thisIndex, mesh_pnts[*it].tbeg, mesh_pnts[*it].tend, t_global);
      }
      assert(mesh_pnts[*it].tend == t_global);
    }

    int nactive = active_list.size();
    assert(nactive <= local_n);

    contribute(MainCB);
  }

  /*******************************/

  void System::BuildLocalMesh(CkCallback &cb)
  {
    MainCB = cb;
    
    assert(ptcl_import_ptr == NULL);
    assert(ptcl_act.empty());

    ptcl_import_ptr = new std::vector<Particle>();

    const int ngbPass = 2;
    systemProxy[thisIndex].localMesh_import(ngbPass, CkCallback(CkIndex_System::BuildLocalMeshII(), systemProxy[thisIndex]));
  }
  void System::BuildLocalMeshII()
  {
    localMesh_build();
#if 1
    for (int i = 0; i < nactive_loc; i++)
    {
      const int id = ptcl_act[i]->id();
      ptcl_act[i]->set_volume(cell_list[i].Volume);

      if (Problem_meshpoint_refine(id)) 
        ptcl_act[i]->set_refine();
      if (Problem_meshpoint_derefine(id))
        ptcl_act[i]->set_derefine();
    }
    ImportPtclActCb = CkCallback(CkIndex_System::BuildLocalMeshIII(), systemProxy[thisIndex]);
    contribute(CkCallback(CkIndex_System::ImportPtclAct(), thisProxy));
#else
    contribute(MainCB);
#endif
  } 
  void System::BuildLocalMeshIII()
  {
    contribute(MainCB);
  }
  
  /*******************************/

  void System::RefineDerefine(CkCallback &cb)
  {
    MainCB = cb;
    contribute(cb);

  }
 
 
  /*******************************/

  void System::ComputeFluidUpdate(CkCallback &cb)
  {
    MainCB = cb;
    get_active_faces();
   
    /***** predict non-active local particles to current time *****/ 

    for (int i = nactive_loc; i < nactive_loc+nimport_loc; i++)
    {
      assert(!ptcl_act[i]->is_active());
      const real  dt = t_global - mesh_act[i]->tbeg;
      assert(dt > 0.0);
      for (int k = 0; k < Fluid::NFLUID; k++)
        Wrec_act[i]->w[k] += Wrec_act[i]->t[k]*dt;
    }

    computePrimitives(false);
    computePvel();
    
    for (int i = 0; i < nactive_loc+nimport_loc; i++)
    {
      assert(ptcl_act[i]->chare() == thisIndex);
      Wrec_act[i]->vel = mesh_act[i]->vel;
    }
    
    ImportFluidPrimitivesCb = CkCallback(CkIndex_System::ComputeFluidUpdateII(), systemProxy[thisIndex]);
    contribute(CkCallback(CkIndex_System::CkIndex_System::ImportFluidPrimitives(), thisProxy));
  }
  void System::ComputeFluidUpdateII()
  {
    computeTimestep();
    ImportNewDtCb = CkCallback(CkIndex_System::ComputeFluidUpdateIII(), systemProxy[thisIndex]);
    contribute(CkCallback(CkIndex_System::CkIndex_System::ImportNewDt(), thisProxy));
  }
  void System::ComputeFluidUpdateIII() 
  {
    computeTimestepLimiter();
    exchangeNewDtAct_returnCb = CkCallback(CkIndex_System::ComputeFluidUpdateIV(), thisProxy);
    contribute(CkCallback(CkIndex_System::exchangeNewDtAct(), thisProxy));
  }
  void System::ComputeFluidUpdateIV()
  {
    for (int i = 0; i < (int)ptcl_act.size(); i++)
    {
      Wrec_act[i]->tend = mesh_act[i]->tend;
      Wrec_act[i]->vel  = mesh_act[i]->vel;
      Wrec_act[i]->pos  = mesh_act[i]->pos;
      Wrec_act[i]->bnd  = mesh_act[i]->boundary;
      Wrec_act[i]->etaJ = mesh_act[i]->etaJ;
      Wrec_act[i]->acc  = mesh_act[i]->acc1;
    }

    computePredictor();

    for (int i = 0; i < nactive_loc; i++)
    {
      const int Id = ptcl_act[i]->id();

      Problem_correct_meshpoint_position(Id);
      ptcl_list[Id] = mesh_pnts[Id].pos_orig;
      mesh_pnts[Id].rung = scheduler.get_rung(mesh_pnts[Id].dt_new);
      mesh_pnts[Id].tbeg = t_global;
      mesh_pnts[Id].tend = t_global + scheduler.get_dt(mesh_pnts[Id].rung);
      assert(mesh_pnts[Id].tend > mesh_pnts[Id].tbeg);

      const int nface = cell_list[i].ngb.size();
      ngb_list[Id] = Neighbours< pair<int, int> >();

      for (int iface = 0; iface < nface; iface++)
      {
        const Face &face = face_list[cell_list[i].ngb[iface]];
        const int j      = face.ngb<true>(i);
        ngb_list[Id].push_back(ptcl_act[j]->chare_id());
      }

      scheduler.push_particle(Id, mesh_pnts[Id].rung);
    }

    for (int i = 0; i < (int)ptcl_act.size(); i++)
    {
      Wrec_act[i]->tend = mesh_act[i]->tend;
      Wrec_act[i]->vel  = mesh_act[i]->vel;
      Wrec_act[i]->pos  = mesh_act[i]->pos;
      Wrec_act[i]->bnd  = mesh_act[i]->boundary;
      Wrec_act[i]->etaJ = mesh_act[i]->etaJ;
      Wrec_act[i]->acc  = mesh_act[i]->acc1;
    }

    computeReconstruction();

    ImportFluidData_completeCb = CkCallback(CkIndex_System::ComputeFluidUpdateV(), systemProxy[thisIndex]);
    contribute(CkCallback(CkIndex_System::CkIndex_System::ImportFluidData(), thisProxy));
  }

  void System::ComputeFluidUpdateV()  
  {
    for (int i = 0; i < (int)ptcl_act.size(); i++)
    {
      Wrec_act[i]->tend = mesh_act[i]->tend;
      Wrec_act[i]->vel  = mesh_act[i]->vel;
      Wrec_act[i]->pos  = mesh_act[i]->pos;
      Wrec_act[i]->bnd  = mesh_act[i]->boundary;
      Wrec_act[i]->etaJ = mesh_act[i]->etaJ;
      Wrec_act[i]->acc  = mesh_act[i]->acc1;
    }
#if 0
    systemProxy[thisIndex].computeFluidUpdate_main(CkCallback(CkIndex_System::Iterate_complete(), thisProxy));
#else
    systemProxy[thisIndex].computeFluidUpdate_main(MainCB);
#endif
  }

  void System::CompleteIteration(CkCallback &cb)
  {
    ptcl_list.resize(local_n);
    assert(mesh_pnts.size() == ptcl_list.size());
    assert(   U_list.size() == ptcl_list.size());
    assert(  dU_list.size() == ptcl_list.size());
    assert(Wrec_list.size() == ptcl_list.size());
    assert(Wrec_minmax_list.size() == ptcl_list.size());

    clear_vec(cell_list);
    clear_vec(face_list);
    clear_vec(face_active_list);

    clear_vec(ptcl_act);
    clear_vec(mesh_act);
    clear_vec(dU_act);
    clear_vec(U_act);
    clear_vec(Wrec_act);
    clear_vec(Wrec_minmax_act);
    clear_vec(slopeLimiter_act);

    delete ptcl_import_ptr;
    delete mesh_import_ptr;
    delete Wrec_import_ptr;
    delete Wrec_minmax_import_ptr;
    delete   dU_import_ptr;
    delete slopeLimiter_import_ptr;

    ptcl_import_ptr = NULL;
    mesh_import_ptr = NULL;
    Wrec_import_ptr = NULL;
    Wrec_minmax_import_ptr = NULL;
    dU_import_ptr   = NULL;
    slopeLimiter_import_ptr = NULL;

    localMesh_destroy();

    for (std::vector<int>::iterator it = active_list.begin(); it != active_list.end(); it++)
      ptcl_list[*it].unset_active();

    iteration++;
    {
      const int Nel = 5;
      std::vector<double> data(Nel, 0.0);
      if (thisIndex == 0)
      {
        data[0] = t_global;
        data[1] = dt_global;
        data[2] = iteration;
      }
      data[3] = active_list.size(); 
      data[4] = local_n;
      contribute(Nel*sizeof(double), &data[0], CkReduction::sum_double, cb);
    }
  }
#endif

#if 1
  void System::get_active_list(CkCallback &cb)
  {
    get_active_list_cb = cb;
    int max_rung_loc = scheduler.get_max_rung();
    contribute(sizeof(int), &max_rung_loc, CkReduction::max_int, CkCallback(CkIndex_System::get_active_listII(0), thisProxy));
  }

  void System::get_active_listII(CkReductionMsg *msg)
  {
    scheduler.set_max_rung(*(int*)msg->getData());
    delete msg;

    active_list.clear();
    t_previous = scheduler.get_tsys();
    assert(t_previous >= 0.0);
    assert(t_previous == t_global);
    dt_global = scheduler.pull_active_list(active_list);
    t_global  = scheduler.get_tsys();

    for (std::vector<int>::iterator it = active_list.begin(); it != active_list.end(); it++)
    {
      ptcl_list[*it].set_active();
      if (!(mesh_pnts[*it].tend == t_global))
      {
        fprintf(stderr ," thisIndex= %d  tbeg= %g tend= %g  t_global= %g \n",
            thisIndex, mesh_pnts[*it].tbeg, mesh_pnts[*it].tend, t_global);
      }
      assert(mesh_pnts[*it].tend == t_global);
    }

    int nactive = active_list.size();
    assert(nactive <= local_n);

    get_active_list_cb.send();
  }

  void System::get_active_faces()
  {
    face_active_list.clear();
    const int nactive = active_list.size();
    for (int i = 0; i < nactive; i++)
    {
      if (mesh_pnts[active_list[i]].is_boundary()) continue;
      assert(ptcl_list[active_list[i]].is_active());
      const Cell &ci = cell_list[i];

      const int nface = ci.ngb.size();
      for (int iface = 0; iface < nface; iface++)
      {
        Face &face = face_list[ci.ngb[iface]];
        if (face.s1 < 0) continue;
        face_active_list.push_back(&face);
        face.s1 = -1-face.s1;
      }
    }

    for (std::vector<Face*>::iterator it = face_active_list.begin(); it != face_active_list.end(); it++)
    {
      Face &face = **it;
      assert(face.s1 < 0);
      face.s1 = -1-face.s1;
    }
  }

  void System::Iterate(CkCallback &cb)
  {
    Iterate_completeCb = cb;
    assert(ptcl_import_ptr == NULL);
    assert(ptcl_act.empty());
    ptcl_import_ptr = NULL;
    systemProxy[thisIndex].get_active_list(CkCallback(CkIndex_System::IterateII(), systemProxy[thisIndex]));
  }

  void System::IterateII()
  {
    assert(thisIndex >= 0);
    assert(thisIndex <  numElements);
    assert(ptcl_import_ptr == NULL);
    if (!ptcl_act.empty())
    {
      fprintf(stderr, "thisIndex= %d  ptcl_act.size= %d \n",
          thisIndex, (int)ptcl_act.size());
    }
    assert(ptcl_act.empty());

    ptcl_act.clear();
    ptcl_import_ptr = new std::vector<Particle>();

    const int ngbPass = 2;
#if 0
    systemProxy[thisIndex].localMesh_import(ngbPass, CkCallback(CkIndex_System::IterateIII(), systemProxy[thisIndex]));
#else
    localMesh_import(ngbPass, CkCallback(CkIndex_System::IterateIII(), systemProxy[thisIndex]));
#endif
  }

  void System::IterateIII()
  {
    localMesh_build ();
    get_active_faces();

    /////////////////

#if 1
    for (int i = nactive_loc; i < nactive_loc+nimport_loc; i++)
    {
      assert(!ptcl_act[i]->is_active());
      const real  dt = t_global - mesh_act[i]->tbeg;
      assert(dt > 0.0);
      for (int k = 0; k < Fluid::NFLUID; k++)
        Wrec_act[i]->w[k] += Wrec_act[i]->t[k]*dt;
    }
#endif

    computePrimitives(true);
    computePrimitives(false);
    computePvel();

    for (int i = 0; i < nactive_loc+nimport_loc; i++)
    {
      assert(ptcl_act[i]->chare() == thisIndex);
      Wrec_act[i]->vel = mesh_act[i]->vel;
    }

    /////////////////

    ImportFluidPrimitivesCb = CkCallback(CkIndex_System::IterateIV(), systemProxy[thisIndex]);
    contribute(CkCallback(CkIndex_System::CkIndex_System::ImportFluidPrimitives(), thisProxy));
  }

  void System::IterateIV()
  {
    computeTimestep();

#if 0
    exchangeNewDt_returnCb = CkCallback(CkIndex_System::IterateVII(), thisProxy);
    contribute(CkCallback(CkIndex_System::exchangeNewDt(), thisProxy));
#else
    ImportNewDtCb = CkCallback(CkIndex_System::IterateVII(), systemProxy[thisIndex]);
    contribute(CkCallback(CkIndex_System::CkIndex_System::ImportNewDt(), thisProxy));
#endif
  }
  void System::IterateVII() 
  {
    computeTimestepLimiter();
#if 1
    exchangeNewDtAct_returnCb = CkCallback(CkIndex_System::IterateVI(), thisProxy);
    contribute(CkCallback(CkIndex_System::exchangeNewDtAct(), thisProxy));
#else
    exchangeNewDt_returnCb = CkCallback(CkIndex_System::IterateVI(), thisProxy);
    contribute(CkCallback(CkIndex_System::exchangeNewDt(), thisProxy));
#endif
  }

  void System::IterateVI()  
  {
    for (int i = 0; i < (int)ptcl_act.size(); i++)
    {
      Wrec_act[i]->tend = mesh_act[i]->tend;
      Wrec_act[i]->vel  = mesh_act[i]->vel;
      Wrec_act[i]->pos  = mesh_act[i]->pos;
      Wrec_act[i]->bnd  = mesh_act[i]->boundary;
      Wrec_act[i]->etaJ = mesh_act[i]->etaJ;
      Wrec_act[i]->acc  = mesh_act[i]->acc1;
    }

    computePredictor();

    for (int i = 0; i < nactive_loc; i++)
    {
      const int Id = ptcl_act[i]->id();

      Problem_correct_meshpoint_position(Id);
      ptcl_list[Id] = mesh_pnts[Id].pos_orig;
      mesh_pnts[Id].rung = scheduler.get_rung(mesh_pnts[Id].dt_new);
      mesh_pnts[Id].tbeg = t_global;
      mesh_pnts[Id].tend = t_global + scheduler.get_dt(mesh_pnts[Id].rung);
      assert(mesh_pnts[Id].tend > mesh_pnts[Id].tbeg);

      const int nface = cell_list[i].ngb.size();
      ngb_list[Id] = Neighbours< pair<int, int> >();

      for (int iface = 0; iface < nface; iface++)
      {
        const Face &face = face_list[cell_list[i].ngb[iface]];
        const int j      = face.ngb<true>(i);
        ngb_list[Id].push_back(ptcl_act[j]->chare_id());
      }

      scheduler.push_particle(Id, mesh_pnts[Id].rung);
    }

    for (int i = 0; i < (int)ptcl_act.size(); i++)
    {
      Wrec_act[i]->tend = mesh_act[i]->tend;
      Wrec_act[i]->vel  = mesh_act[i]->vel;
      Wrec_act[i]->pos  = mesh_act[i]->pos;
      Wrec_act[i]->bnd  = mesh_act[i]->boundary;
      Wrec_act[i]->etaJ = mesh_act[i]->etaJ;
      Wrec_act[i]->acc  = mesh_act[i]->acc1;
    }

    computeReconstruction();

    ImportFluidData_completeCb = CkCallback(CkIndex_System::IterateV(), systemProxy[thisIndex]);
    contribute(CkCallback(CkIndex_System::CkIndex_System::ImportFluidData(), thisProxy));
  }

  void System::IterateV()  
  {
    for (int i = 0; i < (int)ptcl_act.size(); i++)
    {
      Wrec_act[i]->tend = mesh_act[i]->tend;
      Wrec_act[i]->vel  = mesh_act[i]->vel;
      Wrec_act[i]->pos  = mesh_act[i]->pos;
      Wrec_act[i]->bnd  = mesh_act[i]->boundary;
      Wrec_act[i]->etaJ = mesh_act[i]->etaJ;
      Wrec_act[i]->acc  = mesh_act[i]->acc1;
    }
    systemProxy[thisIndex].computeFluidUpdate_main(CkCallback(CkIndex_System::Iterate_complete(), thisProxy));
  }

  void System::IterateVIII()  {}
  void System::IterateIX()  {}
  void System::IterateX()  {}

  void System::Iterate_complete()
  {
    ptcl_list.resize(local_n);
    assert(mesh_pnts.size() == ptcl_list.size());
    assert(   U_list.size() == ptcl_list.size());
    assert(  dU_list.size() == ptcl_list.size());
    assert(Wrec_list.size() == ptcl_list.size());
    assert(Wrec_minmax_list.size() == ptcl_list.size());

    clear_vec(cell_list);
    clear_vec(face_list);
    clear_vec(face_active_list);

    clear_vec(ptcl_act);
    clear_vec(mesh_act);
    clear_vec(dU_act);
    clear_vec(U_act);
    clear_vec(Wrec_act);
    clear_vec(Wrec_minmax_act);
    clear_vec(slopeLimiter_act);

    delete ptcl_import_ptr;
    delete mesh_import_ptr;
    delete Wrec_import_ptr;
    delete Wrec_minmax_import_ptr;
    delete   dU_import_ptr;
    delete slopeLimiter_import_ptr;

    ptcl_import_ptr = NULL;
    mesh_import_ptr = NULL;
    Wrec_import_ptr = NULL;
    Wrec_minmax_import_ptr = NULL;
    dU_import_ptr   = NULL;
    slopeLimiter_import_ptr = NULL;

    localMesh_destroy();

    for (std::vector<int>::iterator it = active_list.begin(); it != active_list.end(); it++)
      ptcl_list[*it].unset_active();

    iteration++;
    {
      const int Nel = 5;
      std::vector<double> data(Nel, 0.0);
      if (thisIndex == 0)
      {
        data[0] = t_global;
        data[1] = dt_global;
        data[2] = iteration;
      }
      data[3] = active_list.size(); 
      data[4] = local_n;
      contribute(Nel*sizeof(double), &data[0], CkReduction::sum_double, Iterate_completeCb);
    }
  }
#endif

}
