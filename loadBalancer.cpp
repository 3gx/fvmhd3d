#include "fvmhd3d.h"

namespace fvmhd3d
{

  void System::loadBalancer(CkCallback &cb)
  {
    loadBalancer_completeCb = cb;
    AtSync();
  }

  void System::ResumeFromSync()
  {
#if 1
    contribute(loadBalancer_completeCb);
#endif
  }


  void System::pup(PUP::er &p)
  {
    CBase_System::pup(p);
    __sdag_pup(p);
  

		p|loadBalancer_completeCb;
#if 1	
		p|scheduler;
    p|local_n;

    if (p.isUnpacking())
    {
      ptcl_list.resize(local_n);
      mesh_pnts.resize(local_n);
      U_list   .resize(local_n);
      dU_list  .resize(local_n);
      Wrec_list.resize(local_n);
      Wrec_minmax_list.resize(local_n);
      ngb_list.resize(local_n);
    }
    PUParray(p, &ptcl_list[0], local_n);
    PUParray(p, &mesh_pnts[0], local_n);
    PUParray(p, &   U_list[0], local_n);
    PUParray(p, &  dU_list[0], local_n);
    PUParray(p, &Wrec_list[0], local_n);
    PUParray(p, &Wrec_minmax_list[0], local_n);
    PUParray(p, &ngb_list[0], local_n);

    p|iteration;
    p|t_end;
    p|dt_snap;
    p|dt_restart;
    p|t_global;
    p|dt_global;
    p|gamma_gas;
    p|courant_no;

    if (p.isUnpacking())
    {
      domains = (globalDomains*)CkLocalNodeBranch(globalDomainsProxy);
      assert(domains != NULL);

      slopeLimiter_all.resize(local_n);
      divBi.resize(local_n);

      scheduler.flush_list();
      for (int i = 0; i < local_n; i++)
      {
        assert(ptcl_list[i].   id() == i);
        assert(ptcl_list[i].chare() == thisIndex);
        scheduler.push_particle(i, mesh_pnts[i].rung);
      }

      clear_vec(cell_list);
      clear_vec(active_list);

      ngb_list_hash          = std::vector< Hash<int> >(local_n);
      ngb_list_hash_used     = std::vector<int>();
    }
#endif
  
  }

};
