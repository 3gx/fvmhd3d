#include "fvmhd3d.h"

/****** Problem specific methods ******/

namespace fvmhd3d
{

#define _CYLINDER_

#define GM               35.0
#define GM_EPS           0.005

#define RoutBND 14.0

#if 1
#define BND_RADIUS       0.01
#define BND_RADIUS1      0.02
#endif

#if 0
#define BND_RADIUS       0.003
#define BND_RADIUS0      0.0005
#define BND_RADIUS1      0.01
#endif

#define DENSvac          0.01
// #define TEMPvac          100.0

#define RCLOUD  3.5
#define DCLOUD  1.0e4
// #define TCLOUD  100.0
#define MUCLOUD 2.0

#define _V01_

#define HoR        0.1
#define MAGNETISATION 0.1

#ifdef _C01_
#define prob_string   "C01"
#define XCL						(5.0)
#define YCL           (2.0)
#define VX0           (-120.0)
#define VY0           0.0
#endif

#ifdef _I01_
#define prob_string   "I01"
#define XCL						(5.0)
#define YCL           (1.0)
#define VX0           (-50.0)
#define VY0           0.0
#endif

#ifdef _I02_
#define prob_string   "I02"
#define XCL						(5.0)
#define YCL           (2.0)
#define VX0           (-50.0)
#define VY0           0.0
#endif

#ifdef _V01_
#define prob_string   "V01"
#define XCL						(5.0)
#define YCL           (3.0)
#define VX0           (-30.0)
#define VY0           0.0
#endif

#ifdef _V02_
#define prob_string   "V02"
#define XCL						(5.0)
#define YCL           (3.0)
#define VX0           (-50.0)
#define VY0           0.0
#endif

#ifdef _V03_
#define prob_string   "V03"
#define XCL						(5.0)
#define YCL           (3.0)
#define VX0           (-80.0)
#define VY0           0.0
#endif

#ifdef _V04_
#define prob_string   "V04"
#define XCL						(5.0)
#define YCL           (3.0)
#define VX0           (-100.0)
#define VY0           0.0
#endif

#define MSUN     (1.9891e33)
#define PC       (3.08568024849531e18)
#define kB       (1.380658e-16)
#define M_U      (1.66054e-24)
#define GCONST   (6.67259e-8)
#define SIGMASB  (6.67051e-5)

#define MSCALE   (1.0e5*MSUN)
#define RSCALE   (1.0 * PC)

#define DUNIT    (MSCALE/(RSCALE*RSCALE*RSCALE*M_U*MUCLOUD))
#define VUNIT    (std::sqrt(GCONST*MSCALE/RSCALE)/1e5)
#define TUNIT    (3.0*kB/M_U/1e10)
#define EUNIT    (GCONST*MSCALE/RSCALE)
#define TIMEUNIT (RSCALE/1e5/VUNIT/3.15569088e7/1e6)

  
  inline real get_r2(const vec3 &pos) 
  {
#ifndef _CYLINDER_
    return pos.norm2();
#else
    return sqr(pos.x) + sqr(pos.y);
#endif
  }
  
  inline real get_cs2(const vec3 r)
  {
    return (GM/std::sqrt(get_r2(r)))*sqr(HoR);
  }


  void Main::Problem_set_global_domain()
  {
    const char string[256] = "Orszag-Tang vortex\n";
    sprintf(problem_string, "%s", string);
    
    const double hLbox = 16;
    const vec3 rmin(-hLbox);
    const vec3 rmax(+hLbox);

    global_domain = boundary(rmin, rmax);
    global_domain_size = global_domain.hsize() * 2.0;

  };

  bool System::Problem_computePvel()
  {
    for (int i = 0; i < nactive_loc; i++)
    {
      MeshPoint &p = *mesh_act[i];
      p.vel = 0.0;

      const Fluid W = U_act[i]->to_primitive(p.Volume);
      p.vel = W.get_vel();

#if 0
      const real D  = std::log(W[Fluid::DENS]);
      const real DV = std::log(DVAC);
      const real dD = DV * 0.01/2;
      const real r  = (D-DV)/std::abs(dD);
      const real vfac = (-r > 20) ? 0.0 : 1.0/(1.0 + std::exp(-r));
      assert(vfac >= 0.0);
      assert(vfac <= 1.0);
      p.vel *= vfac;
#endif


#if 1
      const real B2   = W.get_B().norm2();
      const real pres = Problem_compute_pressure(W);
      const real cs2  = (gamma_gas * pres + B2)/W[Fluid::DENS];
      const real vel2 = W.get_vel().norm2();

      const real vabs = std::sqrt(cs2 + vel2);

      const vec3 centroid = cell_list[i].centroid - p.pos;
      const real d = centroid.abs();
      if (d == 0.0) continue;

      const real eta = 0.25f;
      const real ki  = 1.0f;

      const real f1  = 0.9f;
      const real f2  = 1.1f;

      const real R   = std::pow(cell_list[i].Volume * (3.0/(4.0*M_PI)), 1.0/3.0);
      const real fac = d/(eta*R);

      real f;
      if      (fac < f1) f = 0.0;
      else if (fac < f2) f = (d - f1*eta*R)/((f2 - f1)*eta*R);
      else               f = 1.0;

      real tau = d / vabs;

      f *= ki/tau;
      p.vel += centroid*f;

      if (p.is_boundary()) p.vel = 0.0;
      if (p.pos.abs() > RoutBND) p.vel = 0.0;
#endif
      
      p.vel = p.vel - (p.vel*p.pos)*p.pos/p.pos.norm2();  // subtract radial component

#if 0
      p.vel = 0.0;
#endif
    }

    return true;
  }

  void System::Problem_generate_geometry(const int param)
  {

    const double dt_max = 1.0/128.0;
    scheduler = Scheduler(dt_max);

    t_end      = 5.0 + 1.0/65536;
    dt_restart = 1.0/64;
    dt_snap    = 1.0/16;

    dt_restart = std::max(dt_restart, dt_max);
    dt_snap    = std::max(dt_snap,    dt_max);

    int Npnts_glb = 2e5;

#if 1
    Npnts_glb  = 2e6;
#endif


#if 0
    Npnts_glb  = 2e7;
#endif


    const int Npnts =  Npnts_glb/numElements;

    ptcl_list.clear();
    ptcl_list.reserve(Npnts);

    Rand48 rnd;
    rnd.srand(123 + 123*thisIndex);

    {
      const real pow  = 3.0/2.0;

      int pc = 0;

      const real Rmax = global_domain_size.abs();
      const real eps  = BND_RADIUS1;
      assert(pow > 1.0);
      const real r0   = eps/std::sqrt(pow - 1.0);
      const real n3rmax = r0*r0/std::pow(r0*r0 + eps*eps, pow);

      while (pc < Npnts)
      {
        bool pick = false;
        real r = eps;
        while (!pick)
        {
          r   = Rmax*drand48();
          const real n3  = 1.0/std::pow(r*r + eps*eps, pow);
          const real n3r = r*r*n3;
          pick = (n3rmax*drand48() < n3r);
        }
#if 0
        real theta = 0;
        pick = false;
        while (!pick)
        {
          theta = drand48()*M_PI;
          const real sigma = (M_PI/2.0) * 0.075;
          const real f= 0.05;
          const real p = (f + (1.0 - f)*exp(-sqr(theta - M_PI/2.0)/2/sqr(sigma)))*sin(theta);
          pick = (drand48() < p);
        }
#else
        const real theta = acos(1.0 - 2.0*drand48());
#endif
        const real phi   = 2*M_PI*drand48();
        const vec3 pos(r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), r*cos(theta));
        if (pos.x < global_domain.get_rmin().x) continue;
        if (pos.x > global_domain.get_rmax().x) continue;
        if (pos.y < global_domain.get_rmin().y) continue;
        if (pos.y > global_domain.get_rmax().y) continue;
        if (pos.z < global_domain.get_rmin().z) continue;
        if (pos.z > global_domain.get_rmax().z) continue;
        ptcl_list.push_back(Particle(ptcl_list.size(), thisIndex, pos));
        pc++;
      }
    }

    local_n = ptcl_list.size();

    generateGeometry_nRelax = 5;
  }

  void System::Problem_generate_IC(const int param)
  {
    if (thisIndex == 0)
    {
      CkPrintf(" ********* Cloud capture ************* \n");
      CkPrintf(" **** MAGNETISATION=%g     \n", MAGNETISATION);
      CkPrintf(" **** H/R          =%g     \n", HoR  );
      CkPrintf(" **** GRAVITY_MASS= %g     \n", GM);
      CkPrintf(" **** GRAVITY_EPS=  %g     \n", GM_EPS);
      CkPrintf(" ---  \n");	
    }

    gamma_gas  = 1.0;
    courant_no = 0.8;

    t_global  = 0;
    iteration = 0;

    const real xcl = XCL;
    const real ycl = YCL;

    const real vx_cl = VX0/VUNIT;
    const real vy_cl = VY0/VUNIT;

    const real vorb  = std::sqrt(sqr(vx_cl) + sqr(vy_cl));
    const real Rinit = std::sqrt(sqr(xcl)   + sqr(ycl)  );

    const real tinfall = Rinit/vorb;

    const real dcl = (DCLOUD/DUNIT);
    const real cs2 = get_cs2(vec3(xcl, ycl, 0.0)); //TCLOUD*Tunit/sqr(vunit);

    if (thisIndex == 0)
    {
      CkPrintf("Ro= %g pc, x= %g  y= %g; vx= %g vy= %g vt= %g  tinfall= %g Myr [%g]\n",
          Rinit,
          xcl, ycl,
          vx_cl, vy_cl,
          vorb,
          tinfall * TIMEUNIT, tinfall);
    }


    for (int i = 0; i < local_n; i++) 
    {
      const Particle &pi = ptcl_list[i];

      const vec3 &pos = pi.get_pos();

      if (pos.abs() < BND_RADIUS || pos.abs() > RoutBND)
        mesh_pnts[i].boundary = MeshPoint::DIOD;
      else  
        mesh_pnts[i].boundary = MeshPoint::NO_BOUNDARY;

      const real Rdist = (pos - vec3(xcl, ycl, 0.0)).abs();
      const real inv_beta = MAGNETISATION;

      real dens = dcl;
      real pres = dcl * cs2;
      real b0   = std::sqrt(2.0*pres * inv_beta);

      real bx(0), by(0), bz(0);

#if 0
      bx = by = b0/sqrt(2.0);
#else
      bx = by = bz = b0/sqrt(3.0);
#endif

      real vx(vx_cl), vy(vy_cl), vz(0.0);

      real scalar = 1.0;
      if (Rdist > RCLOUD)
      {
        dens = DENSvac/DUNIT;
        const real csig = std::sqrt(get_cs2(pos));
        vx = (1 - 2.0*drand48()) * csig;
        vy = (1 - 2.0*drand48()) * csig;
        vz = (1 - 2.0*drand48()) * csig;
        scalar = -1.0;
      }
      Fluid m;

      m[Fluid::DENS] = dens;
      m[Fluid::ETHM] = get_cs2(pos)*dens;
      m[Fluid::VELX] = vx;
      m[Fluid::VELY] = vy;
      m[Fluid::VELZ] = vz;
      m[Fluid::BX  ] = bx;
      m[Fluid::BY  ] = by;
      m[Fluid::BZ  ] = bz;
      m[Fluid::PSI ] = 0.0;
      m[Fluid::ENTR] = 1.0;

      Wrec_list[i] = Fluid_rec(m);

      mesh_pnts[i].idx  = thisIndex*1000000 + i+1;

    }
  }

  const std::pair<vec3, real> gacc(const vec3 &pos)
  {
    const real ds2  = pos.norm2();
    const real ds   = std::sqrt(ds2);
    const real ids  = (ds > 0.0) ? 1.0/ds : 0.0;
    const real ids2 = ids*ids;

    const real gravity_mass = GM;
    const real gravity_eps  = GM_EPS;
    assert(gravity_eps > 0.0);
    const real h  = gravity_eps;
    const real q  = ds/h;
    const real q2 = q*q;
    const real q3 = q2*q;
    const real q4 = q2*q2;
    const real q5 = q3*q2;

    real pot;
    real acc;
    if (q < 1) 
    {
      acc = -1.0/sqr(h) * (4.0/3*q  - 6.0/5*q3 + 0.5f*q4);
      pot = +1.0/    h  * (2.0/3*q2 - 0.3f *q4 + 0.1f*q5 - 7.0/5);
    } 
    else if (q < 2) 
    {
      acc = -1.0/sqr(h) * (8.0/3*q - 3.0*q2 + 6.0/5*q3 - 1.0/6*q4 - 1.0/15/q2);
      pot = +1.0/    h  * (4.0/3*q2 - q3 + 0.3*q4 - 1.0/30*q5 - 8.0/5 + 1.0/15/q);
    } 
    else 
    {
      acc = -ids2;
      pot = -ids;
    }

    acc *= gravity_mass*ids;
    pot *= gravity_mass;

    return std::make_pair(acc * pos, pot);
  }
  const vec3 project(const vec3 &r, const vec3 &v)
  {
    return (v - r*((r*v)/r.norm2()));
  }
  const std::pair<vec3, vec3> rotate(const vec3 &r, const vec3 &vt, const real dt)
  {
    const real th = vt.abs()*dt/r.abs();

    const real sinth = sin(th);
    const real costh = cos(th);

    const vec3 rv = r.cross(vt);
    const vec3 n = rv/rv.abs();
    const vec3 n2(sqr(n.x), sqr(n.y), sqr(n.z));

    const real Axx = n2.x + (n2.y + n2.z)*costh;
    const real Ayy = n2.y + (n2.x + n2.z)*costh;
    const real Azz = n2.z + (n2.x + n2.y)*costh;
    const real Axy = n.x*n.y*(1-costh) - n.z*sinth;
    const real Axz = n.x*n.z*(1-costh) + n.y*sinth; 
    const real Ayx = n.x*n.y*(1-costh) + n.z*sinth;
    const real Ayz = n.y*n.z*(1-costh) - n.x*sinth;
    const real Azx = n.x*n.z*(1-costh) - n.y*sinth;
    const real Azy = n.y*n.z*(1-costh) + n.x*sinth;

    return std::make_pair(
        vec3(
          vec3(Axx, Axy, Axz)*r,
          vec3(Ayx, Ayy, Ayz)*r,
          vec3(Azx, Azy, Azz)*r),
        vec3(
          vec3(Axx, Axy, Axz)*vt,
          vec3(Ayx, Ayy, Ayz)*vt,
          vec3(Azx, Azy, Azz)*vt)
        );
  }

  void System::Problem_predict_meshpoint_position(const int Id)
  {
    MeshPoint &p = mesh_pnts[Id];

    const std::pair<vec3, real> f0 = gacc(p.pos_orig);
    p.acc0 = f0.first;

    const real dt     = t_global - p.tbeg;
    p.pos = p.pos_orig;
    p.vel = p.vel_orig; 

    if (p.pos.cross(p.vel).norm2() > 0.0)
    {
      const std::pair<vec3,vec3> rv1 = rotate(p.pos, p.vel, dt);
      p.pos = rv1.first;
    }

    const std::pair<vec3, real> f1 = gacc(p.pos);
    p.acc1 = f1.first;
  }

  void System::Problem_correct_meshpoint_position(const int Id)
  {
    MeshPoint &p = mesh_pnts[Id];
    p.pos_orig = p.pos;
    p.vel_orig = p.vel;
  }

  bool System::Problem_compute_update(Fluid &Uc, const int Id)
  {
    return false;
  }

  real System::Problem_extra_timestep_criterion(const int Id)
  {
    return HUGE;
  }

  real System::Problem_compute_ethm_update(const Fluid &W, const int i)
  {
    //		return Problem_compute_ethm_from_entropy(W);
    return get_cs2(mesh_act[i]->pos)*W[Fluid::DENS];
  }

  real System::Problem_compute_pressure(const Fluid &W)
  {
    //    return gamma_gas > 1.0 ? (gamma_gas - 1.0) * W[Fluid::ETHM] : W[Fluid::ETHM];
    return W[Fluid::ETHM];
  }

  real System::Problem_compute_entropy_from_ethm(const Fluid &W)
  {
    return 1.0;
  }

  real System::Problem_compute_ethm_from_entropy(const Fluid &W)
  {
    assert(false);
    return -1.0;
  }






}
