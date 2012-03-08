#include <sys/stat.h>
#include "fvmhd3d.h"

#ifndef __MACOSX_
#define __LINUX__
#endif

#ifdef __MACOSX__
#include <Accelerate/Accelerate.h>
#include <xmmintrin.h>
inline void fpe_catch() {
  _mm_setcsr( _MM_MASK_MASK &~
      (_MM_MASK_OVERFLOW|_MM_MASK_INVALID|_MM_MASK_DIV_ZERO) );
}
#elif defined __LINUX__
#include <fenv.h>
void fpe_catch(void) {
  /* Enable some exceptions. At startup all exceptions are masked. */
  feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}
#else
crap
void fpe_catch(void) {}
#endif

namespace fvmhd3d
{

  CProxy_Main                    mainProxy;
  CProxy_System                systemProxy;
  CProxy_globalDomains  globalDomainsProxy;
  int numElements;
  char     problem_string[256];
  char     problem_path  [256];
  boundary global_domain;
	double   dtWall_doLB;
  double   dtWall_doCheckPoint;
  vec3     global_domain_size;
 
#if 0 
  CProxy_CkCacheManager cachePtcl;
  CProxy_CkCacheManager cacheFluid;
#endif

  Main::Main(CkMigrateMessage* msg)
  {
    IsRestarting = true;
  }

  Main::Main(CkArgMsg *m)
  {
    IsRestarting = false;

    fpe_catch();

    const int numElements_min = 8;

    if (m->argc <= 2)
    {
      CkPrintf("Usage: %s numElements output_path \n", m->argv[0]);
      CkExit();
    }
    assert(m->argc > 2);
    numElements = atoi(m->argv[1]);
    sprintf(problem_path, "%s", m->argv[2]);
    CkPrintf(" [main] numElements= %d \n", numElements);
    CkPrintf(" [main] problem_path= %s \n", problem_path);
    if (numElements < numElements_min)
    {
      CkPrintf("FATAL: numElements must be >= %d \n", numElements_min);
      CkExit();
    }
    assert(numElements >= numElements_min);


		dtWall_doLB         = m->argc > 3 ? atof(m->argv[3]) : 120.0;
		dtWall_doCheckPoint = m->argc > 4 ? atof(m->argv[4]) : 420.0;
		CkPrintf(" [main] loadBalancer will be called every %g sec \n", dtWall_doLB);
		CkPrintf(" [main] checkPoint   will be called every %g sec \n", dtWall_doCheckPoint);
    CkPrintf(" [main] execline: ");
    for (int i = 0; i < m->argc; i++)
      CkPrintf("%s ", m->argv[i]);
    CkPrintf("\n");

    Problem_set_global_domain();


    CProxy_globalDomains _globalDomainsProxy = CProxy_globalDomains::ckNew();
    globalDomainsProxy = _globalDomainsProxy;

    systemProxy = CProxy_System::ckNew(numElements);

#if 0
    const int cacheSize = 128;
    cachePtcl  = CProxy_CkCacheManager::ckNew(cacheSize, systemProxy.ckLocMgr()->getGroupID());
    cacheFluid = CProxy_CkCacheManager::ckNew(cacheSize, systemProxy.ckLocMgr()->getGroupID());
#endif


    thisProxy.startSimulation();
  }

  void Main::startSimulation()
  {

    /**************** generate Geometry ***************/

    CkPrintf(" [main::startSimulation] generate_geometry\n");
    int nRelax = 0;
    {
      CkReductionMsg *msg;
      systemProxy.generateGeometry(1, CkCallbackResumeThread((void*&)msg));
      nRelax = std::max(nRelax, *(int*)msg->getData());
      delete msg;
    }
    CkPrintf(" [main::startSimulation] relaxing geometry in %d Lloyd's iterations \n", nRelax);

    /**************** relax Geometry ***************/

    for (int iter = 0; iter < nRelax + 1; iter++)
    {
      CkPrintf(" [main::startSimulation] domainDecomposition\n");
      mainProxy.domainDecomposition(CkCallbackResumeThread());
      CkPrintf(" [main::startSimulation] setDomains\n");
      globalDomainsProxy.setDomains(proc_domains, CkCallbackResumeThread());

      CkPrintf(" [main::startSimulation] move particles around \n");
      systemProxy.moveParticles(CkCallbackResumeThread());

      const double t0 = CkWallTimer();
      CkPrintf(" [main::startSimulation] globalMesh_build \n");
      CkReductionMsg *msg;
      systemProxy.globalMesh_build(iter < nRelax, CkCallbackResumeThread((void*&)msg));
      const real v0 = global_domain_size.x*global_domain_size.y*global_domain_size.z;
      const real v1 = *(double*)msg->getData();
      delete msg;
			CkPrintf(" [main::startSimulation] volume: exact= %g  compute= %g  diff= %g [ %g ] \n",
					v0, v1, v1 - v0, (v1-v0)/v0);
			CkPrintf(" [main::startSimulation] done in %g sec \n", CkWallTimer() - t0);
			{
				const double t0 = CkWallTimer();
				CkPrintf(" *** System::loadBalancer() call \n");
				systemProxy.loadBalancer(CkCallbackResumeThread());
				CkPrintf(" *** System::loadBalancer() done in %g sec \n", CkWallTimer() - t0);
			}

		}

		/**************** generate IC ***************/

		{ 
			CkPrintf(" [main::startSimulation] generate_IC\n");
			const double t0 = CkWallTimer();
			CkReductionMsg *msg;
			systemProxy.generateIC(1, CkCallbackResumeThread((void*&)msg));
			CkPrintf(" [main::startSimulation] 1 done in %g sec :: \n", CkWallTimer() - t0);
			const real v0 = global_domain_size.x*global_domain_size.y*global_domain_size.z;
			const real v1 = (*(double*)msg->getData());
			delete msg;
			CkPrintf(" [main::startSimulation]1  volume: exact= %g  compute= %g  diff= %g [ %g ] \n \n",
					v0, v1, v1 - v0, (v1-v0)/v0);
		}

		/************** dump IC into a file *************/

		E0 = computeEnergy();
		E0.print_mass    (" M1: ");
		E0.print_energy  (" E0: ");
		E0.print_momentum(" M0: ");

    {
      double dt_restart;
      {
        CkReductionMsg *msg;
        systemProxy.get_info(CkCallbackResumeThread((void*&)msg));
        t_global   = ((double*)msg->getData())[0];
        t_end      = ((double*)msg->getData())[1];
        dt_snap    = ((double*)msg->getData())[2];
        dt_restart = ((double*)msg->getData())[3];
        iteration  = (int)((double*)msg->getData())[4];
        global_n   = (int)((double*)msg->getData())[5];
        delete msg;
      }

      {
        CkVec<char> filename(256);
        char filepath[256];
        sprintf(filepath, "%s/init", problem_path);
        mkdir(filepath, S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IWOTH);
        sprintf(&filename[0], "%s/", filepath);
        systemProxy.dump_binary(filename, global_n, CkCallbackResumeThread());
      }

      CkPrintf(" [main::startSimulation] :: start_iteration= %d  global_n= %d  t_global= %g  t_end= %g  dt_snap= %g  dt_restart = %g \n",
          iteration, global_n, t_global, t_end, dt_snap, dt_restart);
      t_snap = t_global + dt_snap;
    }


    nmoved = 0;
    nfac = 5.2;


    dtWall_lastLB = -1.0;
    checkpoint_flag = true;
    tWall0 = CkWallTimer();
    thisProxy.doIterations();
  }

  void Main::doIterations()
  {
    if (dtWall_lastLB == -1.0)
    {
      const double t0 = CkWallTimer();
      CkPrintf(" *** System::loadBalancer() call \n");
      systemProxy.loadBalancer(CkCallbackResumeThread());
      CkPrintf(" *** System::loadBalancer() done in %g sec \n", CkWallTimer() - t0);
      dtWall_lastLB = 0.0;
    }

    /**************** start Simulation **********/

    double dtWall_lastCheckPoint  = 0.0;
    while (t_global < t_end)
    {
      const double t0 = CkWallTimer();
      CkReductionMsg *msg;


#if 0
      systemProxy.Iterate(CkCallbackResumeThread((void*&)msg));
#else
      const double t00 = CkWallTimer();
      systemProxy.GetActiveList (CkCallbackResumeThread());
      const double t10 = CkWallTimer();
      systemProxy.BuildLocalMesh(CkCallbackResumeThread());
      const double t20 = CkWallTimer();
      systemProxy.RefineDerefine(CkCallbackResumeThread());
      const double t30 = CkWallTimer();
      systemProxy.ComputeFluidUpdate(CkCallbackResumeThread());
      const double t40 = CkWallTimer();
      systemProxy.CompleteIteration(CkCallbackResumeThread((void*&)msg));
      const double t50 = CkWallTimer();
#if 0
      CkPrintf(" dt_all= %g :: Active= %g  Mesh= %g  RefDeref= %g  Fluid= %g  Complete= %g \n",
          t50 - t00, t10-t00, t20-t10, t30-t20, t40-t30, t50-t40);
#endif
#endif

      t_global  = ((double*)msg->getData())[0];
      iteration = ((double*)msg->getData())[2];
      const double dt_global = ((double*)msg->getData())[1];
      const int     n_active = ((double*)msg->getData())[3];
      const int     n_global = ((double*)msg->getData())[4];
      global_n = n_global;
      delete msg;
      const double dtWall = CkWallTimer() - t0;
      dtWall_lastLB += dtWall;
      dtWall_lastCheckPoint += dtWall;
      CkPrintf("iter= %d :: t= %g dt= %g Nact= %d  Nglb= %d [ %g ]  done in %g sec :: runtime= %g h [ %g  %g  %g ]   \n", 
          iteration,  t_global, dt_global, n_active, n_global, 1.0*n_active/n_global,
          dtWall, (CkWallTimer() - tWall0)/3600.0,
          n_global/dtWall, n_active/dtWall, n_active/dtWall/CkNumPes());
      assert(n_active > 0);
      assert(n_global > 0);

      nmoved += n_active;

      if (n_active == n_global)
      {
        const Energy E1 = computeEnergy();
        const Energy dE = (E1 - E0)/E0.abs();
        E1.print_energy  (" E1: ");
        dE.print_energy  (" dE: ");
        E1.print_mass    (" M1: ");
        dE.print_mass    (" dM: "); 
        E1.print_momentum(" Mom1: ");
        dE.print_momentum(" dMom: ");

      }

      //      if (n_active == n_global)
      if (t_global >= t_snap)
      {
        assert(n_active == n_global);
        CkVec<char> filename(256);
        char filepath[256];
        sprintf(filepath, "%s/iter%.6d", problem_path, (int)((t_snap+0.001*dt_snap)/dt_snap));
        mkdir(filepath, S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IWOTH);
        sprintf(&filename[0], "%s/", filepath);
        CkPrintf(" >>> Writing snapshot to %s \n", filepath);
        systemProxy.dump_binary(filename, n_global, CkCallbackResumeThread());
        t_snap += dt_snap;
      }


      if ((nmoved > (unsigned int)(nfac*n_global)) && (n_active == n_global))
        //      if(false)
        //      if (t_global >= t_snap)
      {  
        nmoved = 0;
        CkPrintf(" [main::startSimulation] domainDecomposition\n");
        mainProxy.domainDecomposition(CkCallbackResumeThread());
        CkPrintf(" [main::startSimulation] setDomains\n");
        globalDomainsProxy.setDomains(proc_domains, CkCallbackResumeThread());

        CkPrintf(" [main::startSimulation] move particles around \n");
        systemProxy.moveParticles(CkCallbackResumeThread());

        const double t0 = CkWallTimer();
        CkPrintf(" [main::startSimulation] globalMesh_build \n");
        CkReductionMsg *msg;
        systemProxy.globalMesh_build(false, CkCallbackResumeThread((void*&)msg));
        const real v0 = global_domain_size.x*global_domain_size.y*global_domain_size.z;
        const real v1 = *(double*)msg->getData();
        delete msg;
        CkPrintf(" [main::startSimulation] volume: exact= %g  compute= %g  diff= %g [ %g ] \n",
            v0, v1, v1 - v0, (v1-v0)/v0);
        CkPrintf(" [main::startSimulation] done in %g sec \n", CkWallTimer() - t0);

      }

      if (dtWall_lastLB > dtWall_doLB)
      {
        dtWall_lastLB = 0.0;
        const double t0 = CkWallTimer();
        CkPrintf(" *** System::loadBalancer() call \n");
        systemProxy.loadBalancer(CkCallbackResumeThread());
        CkPrintf(" *** System::loadBalancer() done in %g sec \n", CkWallTimer() - t0);
      }

      if (dtWall_lastCheckPoint > dtWall_doCheckPoint)
      {
        char filepath[256];
        if (checkpoint_flag)
        {
          sprintf(filepath, "%s/checkpoint1", problem_path);
          checkpoint_flag = false;
        }
        else
        {
          sprintf(filepath, "%s/checkpoint2", problem_path);
          checkpoint_flag = true;
        }
        CkStartCheckpoint(filepath, CkCallback(CkIndex_Main::restart(), thisProxy));
        return;
      }

    }

    CkExit();
  }

  void Main::restart()
  {
    if (IsRestarting)
    {
      CkPrintf(" --- restarting --- \n");
      dtWall_lastLB = -1.0;
      IsRestarting = false;

      CkVec<char> filename(256);
      char filepath[256];
      sprintf(filepath, "%s/restart", problem_path);
      mkdir(filepath, S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IWOTH);
      sprintf(&filename[0], "%s/", filepath);
      CkPrintf(" >>> Writing snapshot to %s \n", filepath);
      systemProxy.dump_binary(filename, global_n, CkCallbackResumeThread());
      tWall0 = CkWallTimer();
    }
    thisProxy.doIterations();
  }

  void Main::pup(PUP::er &p)
  {
    CBase_Main::pup(p);

    p|proc_domains;
    p|dtWall_lastLB;
    p|t_global;
    p|t_end;
    p|t_snap;
    p|dt_snap;
    p|iteration;
    p|nfac;
    p|nmoved;
    p|E0;
    p|checkpoint_flag;
    p|global_n;
  }


}

#include "fvmhd3d.def.h"

