#ifndef __FVMHD3D_H__
#define __FVMHD3D_H__

#include <pup.h>
#include <vector>
#include <algorithm>
#include <cassert>
#include <charm++.h>
#if 0
#include "CkCache.h"
#endif

#include "vector3.h"
#include "Boundary.h"
#include "bOctree.h"
#include "Scheduler.h"

#define NMAXELEMENTS 65536

template<class T>
inline void clear_vec(std::vector<T> &vec) {std::vector<T>().swap(vec);}
template<class T>
inline void fit_vec(std::vector<T> &vec) {std::vector<T>(vec).swap(vec);}


struct Rand48
{
  double drand(){
    update();
    return (stat&0xFFFFFFFFFFFF)*(1.0/281474976710656.0);
  }
  long lrand(){
    update();
    return (long)(stat>>17)&0x7FFFFFFF;
  }
  long mrand(){
    update();
    return(long)(stat>>16)&0xFFFFFFFF;
  }
  void srand(const long seed){
    stat = (seed<<16)+0x330E;
  }
  private:
  long long stat;
  void update(){
    stat = stat*0x5DEECE66D + 0xB;
  }
};

namespace fvmhd3d
{
  enum {NEXTRASCALARS = 0};

  typedef double real;
  typedef vector3<real> vec3;
  typedef Boundary<vec3> boundary;

  inline real sqr(const real x) {return x*x;}
  struct Particle
  {
    private:
      int    _id, _chare;
      vec3   _pos;
      float  _rmax;
      float  _volume;
      int    _status;

    public:

      void pup(PUP::er &p)
      {
        p|_id;
        p|_chare;
        p|_pos;
        p|_rmax;
        p|_volume;
        p|_status;
      }

      Particle() {}
      Particle(const int id, const int chare, const vec3 &pos, const float rmax = 1e10)   :
        _id(id), _chare(chare), _pos(pos), _rmax(rmax), _volume(-1.0), _status(0)   {}


      const float& get_volume() const { return _volume;}
      void set_volume(const float vol) { _volume = vol;}
      void operator=(const std::pair<int, int> in) {_chare = in.first; _id = in.second;};
      void operator=(const vec3 &pos)              {_pos = pos;}

      const int    id() const {return _id & 0x0FFFFFFF;}
      const int chare() const {return _chare;}
      const std::pair<int, int> chare_id() const {return std::make_pair(_chare, id());}
      const vec3& get_pos() const {return _pos;}

      void      set_active()       {        _id |=  (0x80000000);}
      void    unset_active()       {        _id &= ~(0x80000000);}
      const bool is_active() const {return (_id &    0x80000000) == 0x80000000;}

      void      set_ngb()       {        _id |=  (0x40000000);}
      void    unset_ngb()       {        _id &= ~(0x40000000);}
      const bool is_ngb() const {return (_id &    0x40000000) == 0x40000000;}
      
      void      set_refine()       {        _status |=  (0x40000000);}
      void    unset_refine()       {        _status &= ~(0x40000000);}
      const bool is_refine() const {return (_status &    0x40000000) == 0x40000000;}
      
      void      set_derefine()       {        _status |=  (0x20000000);}
      void    unset_derefine()       {        _status &= ~(0x20000000);}
      const bool is_derefine() const {return (_status &    0x20000000) == 0x20000000;}

      int get_status() const {return _status;}
      void set_status(const int st) {_status = st;}

      void set_rmax(const float rmax) {_rmax = rmax;}
      const float get_rmax() const  {return _rmax;}
  };

  template <class T1, class T2> 
    struct pair
    {
      typedef T1 first_type;
      typedef T2 second_type;

      T1 first;
      T2 second;
      pair() : first(T1()), second(T2()) {}
      pair(const T1& x, const T2& y) : first(x), second(y) {}
      template <class U, class V>
        pair (const pair<U,V> &p) : first(p.first), second(p.second) { }
      template <class U, class V>
        pair (const std::pair<U,V> &p) : first(p.first), second(p.second) { }

      const std::pair<T1, T2> make_pair() const 
      {
        return std::make_pair(first, second);
      }

      void pup(PUP::er &p)
      {
        p|first;
        p|second;
      }
    };

  template<typename T>
    struct Neighbours
    {
      protected:
        CkVec<T> list;
      public:
        Neighbours()                     {list.reserve(32);}
        void push_back(const T v)        {list.push_back(v);}
        void resize (const int  n)       {list.resize(n);}
        size_t size()  const             {return list.size();}
        T  operator[](const int i) const {return list[i];}
        T& operator[](const int i)       {return list[i];}
        void clear()                     {list.clear();}
        ~Neighbours()                    {list.clear();}

#if 1
        void pup(PUP::er &p)
        {
          p|list;
        }
#endif
    };

  template<class T>
    struct CellT
    {
      real Volume;
      vec3 centroid;
      Neighbours<T> ngb;
      CellT() {}
      CellT(const real v, const vec3 &c) : Volume(v), centroid(c) {}

      void pup(PUP::er &p)
      {
        p|Volume;
        p|centroid;
#if 0
        int n = ngb.size();
        p|n;
        if (p.isUnpacking)
          ngb.resize(n);
        PUParray(p, &ngb[0], n);
#else
        p|ngb;
#endif
      }
    };

#include "MeshPoint.h"
#include "Fluid.h"
}

#include "fvmhd3d.decl.h"

namespace fvmhd3d
{
#include "Energy.h"

  extern /* readonly */ CProxy_Main                    mainProxy;
  extern /* readonly */ CProxy_System                systemProxy;
  extern /* readonly */ CProxy_globalDomains  globalDomainsProxy;

  extern /* readonly */ int numElements;

  extern /* readonly */ char     problem_string[256];
  extern /* readonly */ char     problem_path  [256];
  extern /* readonly */ boundary global_domain;
  extern /* readonly */  vec3    global_domain_size;

	extern /* readonly */ double   dtWall_doLB;

#if 0
  extern /* readonly */ CProxy_CkCacheManager cachePtcl;
  extern /* readonly */ CProxy_CkCacheManager cacheFluid;
#endif


  struct std_pair_first_sort
  {
    bool operator() (const pair<int, int> &lhs, const pair<int, int> &rhs)
    {
      return lhs.first < rhs.first;
    }
  };

  struct std_pair_second_sort
  {
    bool operator() (const pair<int, int> &lhs, const pair<int, int> &rhs)
    {
      return lhs.second < rhs.second;
    }
  };

  class MeshFluidGrad_msg : public CMessage_MeshFluidGrad_msg
  {
    public:
      int n;
      int       *  id_list;
      MeshPoint *mesh_pnts;
      Fluid_rec *Wrec_list;
      pair<Fluid_flt, Fluid_flt> *Wrec_minmax;
  };

  class MeshFluidGradComplete_msg : public CMessage_MeshFluidGradComplete_msg
  {
    public:
      int n;
      Particle  *ptcl_list;
      MeshPoint *mesh_pnts;
      Fluid     *   U_list;
      FluidD    *  dU_list;
      Fluid_rec *Wrec_list;
      pair<Fluid_flt, Fluid_flt> *Wrec_minmax;
  };


  struct Face
  {
    protected:
    public:
      int s1, s2;
      vec3 centroid;
      vec3 n;
      Face() {}
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
      const real area()   const {return n.abs();}
  };


  template<class T>
    struct Hash
    {
      CkVec<T> list;
      Hash()  {list.reserve(16);}
      ~Hash() {}

      void pup(PUP::er &p)
      {
        p|list;
      }

      const bool use(const T p)
      {
        const int n = list.size();
        for (int i = 0; i < n; i++)
          if (list[i] == p)
            return false;
        list.push_back(p);
        return true;
      }

      const bool is_used(const T p) const
      {
        const int n = list.size();
        for (int i = 0; i < n; i++)
          if (list[i] == p)
            return true;
        return false;
      }
    };

  class globalDomains : public CBase_globalDomains
  {
    enum {NLEAF = 8, NINS = 4};
    public:
 
    std::vector<boundary> proc_domains;
    bOctree::Tree<NLEAF, NINS> proc_tree;
  

    globalDomains();
    globalDomains(CkMigrateMessage*);

    void setDomains(const CkVec<boundary>&, CkCallback &cb);
    void pup(PUP::er &p);
  };

  class System : public CBase_System
  {
    private:
      
      typedef CellT<int> Cell;

      /****** migratable objects ******/

      Scheduler scheduler;

      int local_n;
      std::vector<Particle > ptcl_list;
      std::vector<MeshPoint> mesh_pnts;
      std::vector<Fluid    >    U_list;
      std::vector<FluidD   >   dU_list;
      std::vector<Fluid_rec> Wrec_list;
      std::vector<pair<Fluid_flt, Fluid_flt> > Wrec_minmax_list;

      std::vector< Neighbours< pair<int, int> >  > ngb_list;

      int  iteration;
      real t_end, dt_snap, dt_restart;
      real t_previous, t_global, dt_global;
      real gamma_gas;
      real courant_no;

      CkCallback loadBalancer_completeCb; 

      /****** objects that are reconstructed every timestep ******/

      globalDomains *domains;      /* pointer to NodeGroup that store global domain data */

      std::vector<float     > divBi;
      std::vector<FluidExtra> Wextra_act;
      std::vector<int> active_list;

      std::vector<  Cell  > cell_list;
      std::vector<  Face  > face_list;
      std::vector<  Face* > face_active_list;

      int nactive_loc;
      int nimport_loc, nimport_glb;
      std::vector<Particle >           *ptcl_import_ptr;
      std::vector<MeshPoint>           *mesh_import_ptr;
      std::vector<Fluid_rec>           *Wrec_import_ptr; 
      std::vector<pair<Fluid_flt, Fluid_flt> > *Wrec_minmax_import_ptr;

      std::vector<Fluid    > *   U_import_ptr; 
      std::vector<FluidD   > *  dU_import_ptr; 

      std::vector<Particle* > ptcl_act;
      std::vector<MeshPoint*> mesh_act;
      std::vector<Fluid_rec*> Wrec_act;
      std::vector<pair<Fluid_flt, Fluid_flt>* > Wrec_minmax_act;
      std::vector<FluidD*   > dU_act;
      std::vector<Fluid*    >  U_act;

      std::vector<Fluid_flt >  slopeLimiter_all;
      std::vector<Fluid_flt > *slopeLimiter_import_ptr;
      std::vector<Fluid_flt*>  slopeLimiter_act;


      std::vector< Hash<int> > ngb_list_hash;
      std::vector<      int  > ngb_list_hash_used;


      std::vector<int> nSend_list;
      std::vector<int> nRecv_list;
      int              nSend_cntr;
      int              nRecv_cntr;

      CkCallback MainCB;

      /****************** globalMesh ****************/

    public:

      void globalMesh_build(const bool, CkCallback&);
      void globalMesh_recvPtcl_ticket();
      void globalMesh_recvPtcl(const CkVec<Particle>&, const int);
      void globalMesh_complete();

    private:

      CkCallback globalMesh_completeCb;
      bool  globalMesh_relax_flag;
      void *globalMesh_T_ptr;
      void *globalMesh_Tvtx_list_ptr;
      std::vector< std::pair<int, real> >  globalMesh_failed_ptcl;


      /**********************************************/

      void init();

    public:
      System();
      System(CkMigrateMessage*);

      /**** Load balancer ****/

      void pup(PUP::er&);
      void ResumeFromSync();
      void loadBalancer(CkCallback&);

      /**** Problem (user-defined) dependent methods ****/

      int generateGeometry_nRelax;
      int generateGeometry_param;
      void generateGeometry(const int param, CkCallback&);
      void generateGeometryII();
      void generateGeometryIII();
      void generateGeometryIV();
      CkCallback generateIC_completeCb;
      int        generateIC_param;
      void generateIC         (const int param, CkCallback&);
      void generateIC_II();
      void generateIC_III();
      void generateIC_IV();
      void generateIC_complete(CkReductionMsg*);

      void Problem_generate_geometry(const int);
      void Problem_generate_IC      (const int);
      bool Problem_compute_update(Fluid&, const int);
      real Problem_compute_pressure(const Fluid&);
      real Problem_compute_ethm_update(const Fluid&, const int);
      real Problem_compute_ethm_from_entropy(const Fluid&);
      real Problem_compute_entropy_from_ethm(const Fluid&);
      real Problem_extra_timestep_criterion(const int);
      bool Problem_computePvel();
			real Problem_enforce_limiter(const int);
      void Problem_set_boundary(const int);

      void Problem_predict_meshpoint_position(const int);
      void Problem_correct_meshpoint_position(const int);
  
      bool Problem_meshpoint_refine(const int);
      bool Problem_meshpoint_derefine(const int);

      void computeEnergy(CkCallback&);

      void get_info(CkCallback&);


      /**** Domain decomposition methods ****/

      void sample_particles(const int);

      CkCallback moveParticles_completeCb;
      bool moveParticles_flag;
      void moveParticles     (CkCallback&);
      void moveParticles_recv(MeshFluidGradComplete_msg*);
      void moveParticles_reduction(CkReductionMsg*);
      void moveParticles_complete();
      void sort_local_data();

      /**** Local mesh construction mehods *****/

      int localMesh_ngbPass;
      CkCallback localMesh_completeCb;
      std::vector<int> localMesh_sites2process;
      void *T_ptr;
      void *Tvtx_ptr;
      int   n_in_DTloc;
      void localMesh_import    (const int, CkCallback&);
      void localMesh_importII  (const std::vector<int>&);
      void localMesh_importIII ();
      void localMesh_import_reduction(CkReductionMsg*);
      void localMesh_import_add2process(const CkVec<int>&, const pair<int, int>);
      void localMesh_import_recvTicket();
      void localMesh_insertPtcl(const CkVec<Particle>&, const int);
      void localMesh_import_complete  ();
      void localMesh_build();
      void localMesh_compute_total_volume(CkCallback&);
      void localMesh_destroy();

      /**** Fluid dynamics methods ****/

      CkCallback computeFluidUpdate_completeCb;
      void computeFluidUpdate_main(CkCallback &cb);

      int slopeLimiter_nRequestedUpdates;
      CkCallback slopeLimiter_completeCb;
      void slopeLimiter(CkCallback&);
      void slopeLimiter_recvLimiter(const CkVec< pair<int, Fluid_flt> >&, const int);
      void slopeLimiter_recvTicket();
      void slopeLimiter_exchange();
      void slopeLimiter_requestLimiter(const CkVec< pair<int, int> >&, const int);
      void slopeLimiter_recvNewLimiter(const CkVec< pair<int, Fluid_flt> >&);

      void computeFluidUpdate         ();
      void computeFluidUpdateI        ();
      void computeFluidUpdateII       ();
      void computeFluidUpdate_complete();

      void computeFluidUpdate_sendIncUpdate();
      void computeFluidUpdate_recvIncUpdate(const CkVec< pair<int, FluidD> >&, const int);
      void computeFluidUpdate_recvIncUpdateTicket();

      CkCallback exchangeNewDt_returnCb;
      void exchangeNewDt();
      void exchangeNewDt_recv(const CkVec< pair<int, real> >&, const int);
      void exchangeNewDt_recvTicket();

      CkCallback exchangeNewDtAct_returnCb;
      void exchangeNewDtAct();
      void exchangeNewDtAct_recv(const CkVec< pair<int, real> >&, const int);
      void exchangeNewDtAct_recvTicket();

      void computePrimitives     (const bool);
      void computePvel           ();
      void computeTimestep       ();
      void computeTimestepLimiter();
      void computeReconstruction ();
      void computePredictor      ();

      void computeFlux(
          const vec3  &wij,
          const Fluid &Wi,
          const Fluid &Wj,
          const int    jbnd,
          const vec3  &normal,
          real &psi_ij,
          real &Bn_ij,
          Fluid &flux);

      void riemannSolver(
          Fluid &flux,
          const real Bx,
          const real w,
          const real dens_L, const real pres_L, const real ethm_L,
          const real velx_L, const real vely_L, const real velz_L, const real By_L, const real Bz_L,
          const real dens_R, const real pres_R, const real ethm_R,
          const real velx_R, const real vely_R, const real velz_R, const real By_R, const real Bz_R
          ) const;
      void riemannSolver_HLLD(
          Fluid &flux,
          const real Bx,
          const real w,
          const real dens_L, const real pres_L, const real ethm_L,
          const real velx_L, const real vely_L, const real velz_L, const real By_L, const real Bz_L,
          const real dens_R, const real pres_R, const real ethm_R,
          const real velx_R, const real vely_R, const real velz_R, const real By_R, const real Bz_R
          ) const;
      void riemannSolver_HLLE(
          Fluid &flux,
          const real Bx,
          const real w,
          const real dens_L, const real pres_L, const real ethm_L,
          const real velx_L, const real vely_L, const real velz_L, const real By_L, const real Bz_L,
          const real dens_R, const real pres_R, const real ethm_R,
          const real velx_R, const real vely_R, const real velz_R, const real By_R, const real Bz_R
          ) const;
      
      /***************************/
      /***************************/
      /***************************/

      void GetActiveList(CkCallback&);
      void GetActiveListII(CkReductionMsg*);
  
      void BuildLocalMesh(CkCallback&);
      void BuildLocalMeshII();
      void BuildLocalMeshIII();
  
      void ComputeFluidUpdate(CkCallback&);
      void ComputeFluidUpdateII();
      void ComputeFluidUpdateIII();
      void ComputeFluidUpdateIV();
      void ComputeFluidUpdateV();
  
      void CompleteIteration(CkCallback&);

      void RefineDerefine(CkCallback&);

      CkCallback ImportPtclActCb;
      int ImportPtclAct_nRequested;
      void ImportPtclAct();
      void ImportPtclAct_request(const CkVec< pair<int,int> >&, const int);
      void ImportPtclAct_recv(const CkVec< pair<int, Particle> >&);

      
      /***************************/
      /***************************/
      /***************************/


      CkCallback get_active_list_cb;
      void get_active_list  (CkCallback&);
      void get_active_listII(CkReductionMsg*);

      void get_active_faces();

      CkCallback Iterate_completeCb;
      void Iterate(CkCallback&);
      void IterateII();
      void IterateIII();
      void IterateIV();
      void IterateV();
      void IterateVI();
      void IterateVII();
      void IterateVIII();
      void IterateIX();
      void IterateX();
      void Iterate_complete();

      /************ ImportFluidData.cpp **********/

      int ImportFluidData_nRequestedUpdates;
      CkCallback ImportFluidData_completeCb;
      void ImportFluidData();
      void ImportFluidData_request(const CkVec< pair<int,int> >&, const int);
      void ImportFluidData_recv(MeshFluidGrad_msg*);

      int ImportFluidPrimitives_nRequested;
      CkCallback ImportFluidPrimitivesCb;
      void ImportFluidPrimitives();
      void ImportFluidPrimitives_request(const CkVec< pair<int,int> >&, const int);
      void ImportFluidPrimitives_recv   (const CkVec< pair<int, pair<Fluid, MeshPoint> > >&);
      
      int ImportNewDt_nRequested;
      CkCallback ImportNewDtCb;
      void ImportNewDt();
      void ImportNewDt_request(const CkVec< pair<int,int > >&, const int);
      void ImportNewDt_recv   (const CkVec< pair<int,real> >&);
      /************ IpmortFluidData.cpp **********/

      const vec3 periodic(const vec3&) const;

      /**************** I/O **************/

      void dump_binary(const CkVec<char>&, const int, CkCallback&);
      void read_binary(const CkVec<char>&);

  };

  class Main : public CBase_Main
  {
    private:

      int               sample_pos_nRecv;
      std::vector<vec3> sample_pos;
      CkVec<boundary>   proc_domains;

      Energy E0;

      bool IsRestarting;
      unsigned long long global_n;
      double t_global;
      double t_end;
      double t_snap;
      double dt_snap;
      int iteration;
      int nfac;
      unsigned long long nmoved;
      double dtWall_lastLB;
      double tWall0;
      bool checkpoint_flag;

    public:

      Main(CkArgMsg*);
      Main(CkMigrateMessage*);

      void Problem_set_global_domain();

      void startSimulation();
      void doIterations();
      CkCallback domainDecomposition_completeCb;
      void domainDecomposition(CkCallback&cb);
      void domainDecompositionII(const CkVec<vec3>&);
      void done();
      void restart();

      Energy computeEnergy();
      void pup(PUP::er &p);
  };
}

#endif /* __FVMHD3D_H__ */
