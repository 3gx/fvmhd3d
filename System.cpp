#include "fvmhd3d.h"
#include "distributenew.h"

namespace fvmhd3d
{
  System::System() 
  {
    usesAtSync = CmiTrue;
    init();
  }
  System::System(CkMigrateMessage *msg) 
  {
    init();
    delete msg;
  }
  
#if 0
  void System::ckAboutToMigrate()
  {
  }


  void System::ckJustMigrated()
  {
  }
#endif

  void System::init()
  {
    nSend_list = std::vector<int>(numElements, 0);
    nRecv_list = std::vector<int>(numElements, 0);
		nSend_cntr = 0;
		nRecv_cntr = 0;
    slopeLimiter_nRequestedUpdates = 0;
    ImportFluidData_nRequestedUpdates = 0;
    ImportFluidPrimitives_nRequested = 0;
    ImportNewDt_nRequested = 0;
    ImportPtclAct_nRequested = 0;

    ptcl_import_ptr = NULL;
    mesh_import_ptr = NULL;
    Wrec_import_ptr = NULL;
    U_import_ptr = NULL;
    dU_import_ptr = NULL;
    slopeLimiter_import_ptr = NULL;
    Wrec_minmax_import_ptr = NULL;
    
    globalMesh_T_ptr = NULL;
    globalMesh_Tvtx_list_ptr = NULL;

    T_ptr      = NULL;
    Tvtx_ptr   = NULL;
    n_in_DTloc = 0;

    ptcl_act.clear();
    mesh_act.clear();
    Wrec_act.clear();
    U_act.clear();
    dU_act.clear();
    slopeLimiter_act.clear();
    Wrec_minmax_act.clear();
  }

  void Main::domainDecomposition(CkCallback &cb)
  {
    domainDecomposition_completeCb = cb;
#if 0
#if 0
    const int NMAXSAMPLE = 1000000;
    const int nsample = std::min(NMAXSAMPLE, int(numPtcl/numElements));
#else
    const int nsample = numPtcl/numElements;
#endif
    const int  sample_freq = (numPtcl + nsample - 1)/nsample;
#else
    const int sample_freq = 10;
#endif

    sample_pos_nRecv = 0;
    sample_pos.clear();
    systemProxy.sample_particles(sample_freq);
  }

  void System::sample_particles(const int sample_freq)
  {
    moveParticles_flag = false;
   
    for (int i = 0; i < local_n; i++)
    {
      const vec3 pos = periodic(ptcl_list[i].get_pos());
      ptcl_list[i] = pos;
      mesh_pnts[i].set_pos(pos);
    }

    CkVec<vec3> ptcl_sample;
    const int n = ptcl_list.size();
    ptcl_sample.reserve(n);
    for (int i = 0; i < n; i += sample_freq)
      ptcl_sample.push_back(ptcl_list[i].get_pos());
    mainProxy.domainDecompositionII(ptcl_sample);
  }


  void Main::domainDecompositionII(const CkVec<vec3> &list)
  {
    sample_pos_nRecv++;
    const int nrecv = list.size();
    for (int i = 0; i < nrecv; i++)
      sample_pos.push_back(list[i]);

    if (sample_pos_nRecv == numElements)
    {
#if 0
      CkPrintf(" [main::domainDecompositionII]  recvd all \n");
#endif
      DistributeNew<real, vec3, boundary> distribute(global_domain);
      distribute.determine_division(sample_pos, numElements);
      proc_domains.resize(numElements);
      double v1 = 0;
      for (int i = 0; i < numElements; i++)
      {
        proc_domains[i] = distribute.tiles[i];
        v1 += 
          proc_domains[i].hsize().x*
          proc_domains[i].hsize().y*
          proc_domains[i].hsize().z;
      }
      v1 *= 8.0;
#if 0
      const real v0 = global_domain_size.x*global_domain_size.y*global_domain_size.z;
      CkPrintf(" [main::domainDecompositionII]  volume: exact= %g  compute= %g  diff= %g [ %g ] \n",
          v0, v1, v1 - v0, (v1-v0)/v0);
#endif
      domainDecomposition_completeCb.send();
    }
  }

  const vec3 System::periodic(const vec3 &pos) const
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

  void System::generateGeometry(const int param, CkCallback &cb)
  {
    Problem_generate_geometry(param);

    mesh_pnts.resize(local_n);
    U_list   .resize(local_n);
    dU_list  .resize(local_n);
    Wrec_list.resize(local_n);
    Wrec_minmax_list.resize(local_n);
    divBi.resize(local_n);
    for (int i = 0; i < local_n; i++)
      mesh_pnts[i].rung = 0;

    contribute(sizeof(int), &generateGeometry_nRelax, CkReduction::max_int, cb);
  }
  void System::generateGeometryII() {}
  void System::generateGeometryIII() {}
  void System::generateGeometryIV() {}

  void System::generateIC(const int param, CkCallback &cb)
  {
    generateIC_completeCb = cb;
    generateIC_param      = param;
    active_list.resize(local_n);
    for (int i = 0; i < local_n; i++)
    {
      active_list[i] = i;
      ptcl_list  [i].set_active();
    }

    Problem_generate_IC(generateIC_param);
    for (int i = 0; i < local_n; i++)
    {
      Wrec_minmax_list[i] = std::make_pair(0.0, 0.0);
      mesh_pnts[i].set_pos(ptcl_list[i].get_pos());
      mesh_pnts[i].set_vel(0.0);
      mesh_pnts[i].acc0  = 0.0;
      mesh_pnts[i].acc1  = 0.0;
      mesh_pnts[i].tbeg   = t_global;
      mesh_pnts[i].tend   = t_global;
      mesh_pnts[i].rung     = 0;
      mesh_pnts[i].status   = 0;
    }
    
    assert(ptcl_import_ptr == NULL);
    assert(ptcl_act.empty());
    
    ptcl_act.clear();
    ptcl_import_ptr = new std::vector<Particle>();

#if 1
    systemProxy[thisIndex].localMesh_import(2, CkCallback(CkIndex_System::generateIC_II(), systemProxy[thisIndex]));
#else
    systemProxy[thisIndex].localMesh_import(2, CkCallback(CkIndex_System::generateIC_II(), thisProxy));
#endif
  }

  void System::generateIC_II()
  {

    localMesh_build();
    get_active_faces();

    for (int i = 0; i < nactive_loc; i++) 
    {
      assert(thisIndex == ptcl_act[i]->chare());
      assert(i == ptcl_act[i]->id());
      mesh_act[i]->Volume = cell_list[i].Volume;
      *U_act[i] = Wrec_act[i]->w.to_conservative(mesh_act[i]->Volume);
      *dU_act[i] = 0.0;
    }

    computePrimitives(false);
    computePvel();
    for (int i = 0; i < nactive_loc; i++)
    {
      Wrec_act[i]->vel = mesh_act[i]->vel;
    }

    assert(nimport_loc == 0);

    for (int i = 0; i < nactive_loc + nimport_loc; i++)
    {
      Wrec_act[i]->tend = mesh_act[i]->tend;
      Wrec_act[i]->vel  = mesh_act[i]->vel;
      Wrec_act[i]->pos  = mesh_act[i]->pos;
      Wrec_act[i]->bnd  = mesh_act[i]->boundary;
      Wrec_act[i]->etaJ = mesh_act[i]->etaJ;

      assert(Wrec_act[i]->w[Fluid::DENS] > 0.0);
      assert(Wrec_act[i]->w[Fluid::ETHM] > 0.0);
    }

    ImportFluidPrimitivesCb = CkCallback(CkIndex_System::generateIC_III(), systemProxy[thisIndex]);
    contribute(CkCallback(CkIndex_System::CkIndex_System::ImportFluidPrimitives(), thisProxy));

  }

  void System::generateIC_III()
  {
    computeTimestep();

    double dt_min = HUGE;
    for (int i = 0 ; i < nactive_loc; i++)
      dt_min = std::min(dt_min, mesh_act[i]->dt_new);
    contribute(sizeof(double), &dt_min, CkReduction::min_double, CkCallback(CkIndex_System::generateIC_complete(0), thisProxy));
  }
  
  void System::generateIC_IV()
  {
  }


  void System::generateIC_complete(CkReductionMsg *msg)
  {
    double dt_min = *(double*)msg->getData();
    delete msg;

    computeTimestep();

    scheduler.flush_list();
    for (int i = 0; i < local_n; i++)
    {
      Wrec_list[i] = Fluid_rec(Wrec_list[i].w);
      dU_list  [i] = 0.0;
      mesh_pnts[i].dt_new = dt_min;
      const real dt_new = mesh_pnts[i].dt_new;
      assert(dt_new > 0.0);
      const int  rung   = scheduler.get_rung(dt_new);
      mesh_pnts[i].rung = rung;
      mesh_pnts[i].tend = mesh_pnts[i].tbeg + scheduler.get_dt(rung);
//      mesh_pnts[i].tlast = mesh_pnts[i].tend;
      assert(mesh_pnts[i].tend > mesh_pnts[i].tbeg);
      scheduler.push_particle(i, rung);
    }

    for (int i = 0; i < local_n; i++)
    {
      mesh_pnts[i].set_pos(mesh_pnts[i].pos);
      mesh_pnts[i].set_vel(mesh_pnts[i].vel);
      ptcl_list[i].unset_active();
      assert(!ptcl_list[i].is_active());
    }

    cell_list.clear();

    ptcl_act.clear();
    mesh_act.clear();
    dU_act.clear();
    Wrec_act.clear();
    Wrec_minmax_act.clear();
    slopeLimiter_act.size();

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

    localMesh_compute_total_volume(generateIC_completeCb);

    localMesh_destroy();
  }

  //////////////

  void System::computeEnergy(CkCallback &cb)
  {
    Energy E(mesh_pnts, U_list);
    contribute(Energy::NVAR*sizeof(double), E.data, CkReduction::sum_double, cb);
  } 

  Energy Main::computeEnergy()
  {
    CkReductionMsg *msg;
    systemProxy.computeEnergy(CkCallbackResumeThread((void*&)msg));
    const Energy E((double*)msg->getData());
    delete msg;
    return E;
  }

  void System::get_info(CkCallback &cb)
  {
    const int Nel = 6;
    std::vector<double> data(Nel, 0.0);
    if (thisIndex == 0)
    {
      data[0] = t_global;
      data[1] = t_end;
      data[2] = dt_snap;
      data[3] = dt_restart;
      data[4] = iteration;
    }
    data[5] = local_n;
    contribute(Nel*sizeof(double), &data[0], CkReduction::sum_double, cb);
  }

}
