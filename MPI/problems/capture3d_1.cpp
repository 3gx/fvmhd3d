#include "fvmhd3d.h"

namespace fvmhd3d 
{

#define GM               35.0
#define GM_EPS           0.003

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

#define _C01_

#define HoverR        0.03
#define MAGNETISATION 0.0

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


  bool system::compute_pvel1()
  {
    const int nactive_site = site_active_list.size();
    for (int isite = 0; isite < nactive_site; isite++)
    {
      const int i = site_active_list[isite];
      assert(ptcl_import[i].boundary == 0);

      Particle &p = ptcl_import[i];

      const Fluid &W = W_st[i].w;
      
      p.vel = vec3(W[Fluid::VELX], W[Fluid::VELY], W[Fluid::VELZ]);

      const real B2   = sqr(W[Fluid::BX]) + sqr(W[Fluid::BY]) + sqr(W[Fluid::BZ]);
      const real pres = compute_pressure(W);

      const real cs2  = (gamma_gas * pres + B2)/W[Fluid::DENS];
      const real vel2 = vec3(W[Fluid::VELX], W[Fluid::VELY], W[Fluid::VELZ]).norm2();
      const real vabs = std::sqrt(cs2) + std::sqrt(vel2);

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

      assert(p.idx < 0);

      {
        const vec3 R(p.pos.x, p.pos.y, 0.0);
        const vec3 V(p.vel.x, p.vel.y, 0.0);
        p.vel = V - (V*R)*R/R.norm2();  // subtract radial component
      }

      //      p.vel = 0.0;

      W_st[i].vel = p.vel;
    }
    return true;
  }

  const vec3 project(const vec3 &r, const vec3 &v)
  {
    const vec3 V(v.x, v.y, 0.0);
    const vec3 R(r.x, r.y, 0.0);
    return (V - R*((R*V)/R.norm2()));
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

  void system::predict(const int i)
  {
    const real dth = 0.5*(ptcl_local[i].tend - ptcl_local[i].tlast);
    if (ptcl_local[i].idx < 0) // particles move on a sphere of const R
    {
      const vec3 r0(ptcl_local[i].orig_pos);
      const vec3 v0(ptcl_local[i].orig_vel);

      if (r0.cross(v0).norm2() == 0.0) 
      {
        ptcl_local[i].pos = ptcl_local[i].orig_pos;
        ptcl_local[i].pos = periodic(ptcl_local[i].pos);
        ptcl_local[i].vel = project(ptcl_local[i].pos, ptcl_local[i].orig_vel);
      }
      else
      {
        ptcl_local[i].vel = project(r0, v0);
        const std::pair<vec3,vec3> rv1 = rotate(r0, ptcl_local[i].vel, dth);
        ptcl_local[i].pos = periodic(rv1.first);
      }
    }
    else
    {
      fprintf(stderr, "i= %d : ptcl_local[i].idx= %d \n",
          i, ptcl_local[i].idx);
      assert(false);
      ptcl_local[i].pos = periodic(ptcl_local[i].orig_pos + ptcl_local[i].orig_vel * dth);
      ptcl_local[i].vel = ptcl_local[i].orig_vel;
    }
  }

  void system::correct(const int i)
  {
    const real dth_pre = 0.5*(ptcl_local[i].tend - ptcl_local[i].tlast);
    assert(ptcl_local[i].tend > ptcl_local[i].tlast);
    assert(scheduler.get_tsys() > ptcl_local[i].tlast);
    assert(scheduler.get_tsys() == ptcl_local[i].tend);
    ptcl_local[i].tlast = scheduler.get_tsys();
    ptcl_local[i].tend  = scheduler.get_tsys() + scheduler.get_dt(ptcl_local[i].rung);

    if (ptcl_local[i].idx < 0) // particles move on a sphere of const R
    {
      const vec3 r0(ptcl_local[i].pos);
      const vec3 v0(ptcl_local[i].vel);

      if (r0.cross(v0).norm2() == 0.0) 
      {
        ptcl_local[i].orig_pos = ptcl_local[i].pos;
        ptcl_local[i].orig_pos = periodic(ptcl_local[i].orig_pos);
        ptcl_local[i].orig_vel = project(ptcl_local[i].orig_pos, ptcl_local[i].vel);
      }
      else
      {
        const vec3 va = project(r0, v0);
        const std::pair<vec3,vec3> rv1 = rotate(r0, va, dth_pre);
        ptcl_local[i].orig_pos = periodic(rv1.first);
        ptcl_local[i].orig_vel = rv1.second;

        const vec3 r0 = ptcl_local[i].orig_pos;
        const vec3 v0 = ptcl_local[i].orig_vel;
      }
    }
    else
    {
      assert(false);
      if (ptcl_local[i].pos.abs() < 1.0*BND_RADIUS)
      {
        ptcl_local[i].orig_pos = ptcl_local[i].pos;
        ptcl_local[i].orig_vel = ptcl_local[i].vel = 0.0;
        if (ptcl_local[i].boundary == 0) boundary_n++;
        ptcl_local[i].boundary = Particle::DIOD;
      }
      else
      {
        ptcl_local[i].orig_vel =          ptcl_local[i].vel;
        ptcl_local[i].orig_pos = periodic(ptcl_local[i].pos + ptcl_local[i].vel * dth_pre);
      }
    }
  }


  bool system::compute_update_prob(Fluid &U, const int i0)
  {
    const int i = (i0 < 0) ? -1-i0 : i0;

#if 0
    Fluid W = (i0 >= 0) ? U.to_primitive(cell_list[i].Volume) : U;
#endif

    if (i0 < 0)
    {
      const Fluid &W = U;

      const real DFLOOR = 3.33*DENSvac/DUNIT;
      if (
          W[Fluid::DENS] <= 1.001*DFLOOR || 
          ptcl_import[i].pos.abs() > 10.0 ||
          ptcl_import[i].pos.abs() < 1.5*BND_RADIUS)
      {
        W_rec[i].t = 0.0;
        W_rec[i].x = 0.0;
        W_rec[i].y = 0.0;
        W_rec[i].z = 0.0;
      }
    }

#if 0
    U = (i0 >= 0) ? W.to_conservative(cell_list[i].Volume) : W;
#endif

    return false;
  }

  const vec3 system::gpot_acc(const int i, real &gpot) const
  {
    const vec3 pos  = ptcl_import[i].pos;
    const real ds2  = pos.abs();
    const real ds   = std::sqrt(ds2);
    const real ids  = (ds > 0.0) ? 1.0/ds : 0.0;
    const real ids2 = ids*ids;

    const real gravity_eps  = GM_EPS;
    const real gravity_mass = GM;
    assert(gravity_eps > 0.0f);
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

    gpot = pot;
    return acc * pos;
  }
  inline real get_r2(const vec3 &r)
  {
#if 1
    return r.norm2();
#else
    return sqr(r.x) + sqr(r.y);
#endif
  }
  inline real get_cs2(const vec3 r)
  {
    return (GM/std::sqrt(get_r2(r)))*sqr(HoverR);
  }
  real system::compute_ethm_update(const Fluid &W, const int i) const
  {
    return get_cs2(ptcl_import[i].pos)*W[Fluid::DENS];
  }
  real system::compute_pressure(const Fluid& W) const
  {
    return W[Fluid::ETHM];
  }
  real system::compute_entropy_from_ethm(const Fluid& W) const
  {
    return 1.0;
  }
  real system::compute_ethm_from_entropy(const Fluid& W) const
  {
    assert(false);
    return -1.0;
  }

  bool system::refine(const int i) 
  {
    return false;
  }
  bool system::derefine(const int i)
  {
    return false;
  }

  void system::set_geometry(const bool init) 
  {
    const double dt_max = 1.0/64;
    scheduler = Scheduler(dt_max);

    t_end   = 2.5 + 1.0/128;

    n_restart  = 2;
    dt_restart = dt_max;

    dt_dump = 1.0/16;

    di_log = 100;

    global_n = local_n = 0;


    const double Lbox = 32;
    const vec3 rmin(-Lbox/2.0, -Lbox/2.0, -Lbox/2.0);
    const vec3 rmax(+Lbox/2.0, +Lbox/2.0, +Lbox/2.0);
    global_domain = boundary(rmin, rmax);
    global_domain_size = global_domain.hsize() * 2.0;

    const vec3 Len3 = global_domain.hsize() * 2.0;
    pfloat<0>::set_scale(Len3.x);
    pfloat<1>::set_scale(Len3.y);
    pfloat<2>::set_scale(Len3.z);

    Distribute::int3 nt(1, 1, 1);
    switch(nproc) 
    {
      case 1: break;
      case 2: nt.x = 2; nt.y = 1; nt.z = 1; break;
      case 4: nt.x = 2; nt.y = 2; nt.z = 1; break;
      case 6: nt.x = 3; nt.y = 2; nt.z = 1; break;
      case 8: nt.x = 4; nt.y = 2; nt.z = 1; break;
      case 10: nt.x = 5; nt.y = 2; nt.z = 1; break;
      case 12: nt.x = 4; nt.y = 3; nt.z = 1; break;
      case 14: nt.x = 7; nt.y = 2; nt.z = 1; break;
      case 16: nt.x = 4; nt.y = 4; nt.z = 1; break;
      case 32: nt.x = 8; nt.y = 4; nt.z = 1; break;
      case 64: nt.x = 8; nt.y = 8; nt.z = 1; break;
      case 128: nt.x = 8; nt.y = 8; nt.z = 2; break;
      default: assert(false);
    }

    const Distribute::int3 nt_glb(nt);
    const pBoundary pglobal_domain(pfloat3(0.0), pfloat3(Len3));
    distribute_glb.set(nproc, nt, pglobal_domain);

    if (!init) return;

    if (myproc == 0) 
    {
      int Namb = 1e5;
      int Ncld = 1e5;

#if 1
      Namb = 0.5e6;
      Ncld = 1.5e6;
#endif

#if 0
      Namb = 2.0e5;
      Ncld = 1;
#endif

#if 0
      Namb = 2.0e6;
      Ncld = 1;
#endif

      ptcl_local.clear();
      ptcl_local.reserve(128);

      fprintf(stderr, "Ncloud= %d \n", Ncld);
      const real dr_amb = std::pow(Lbox*Lbox*Lbox/Namb,      1.0/3.0);
      const real dr_cld = 2.0*RCLOUD*std::pow(1.0/Ncld, 1.0/3.0);
      fprintf(stderr, "dr_amb=  %g   dr_cloud= %g\n", dr_amb, dr_cld);
      fprintf(stderr, "rmin= %g %g %g \n", 
          global_domain.rmin.x,
          global_domain.rmin.y,
          global_domain.rmin.z);
      fprintf(stderr, "rmax= %g %g %g \n", 
          global_domain.rmax.x,
          global_domain.rmax.y,
          global_domain.rmax.z);
      fprintf(stderr, "\n -------------------------- \n\n");

#if 0
      int Nbnd = 32;
      {
        int pc = 0;
        vec3 pos(HUGE/100);
        while (pc < Nbnd)
        {
          while(pos.abs() > BND_RADIUS)
          {
            pos.x = (1.0 - 2.0*drand48()) * BND_RADIUS;	
            pos.y = (1.0 - 2.0*drand48()) * BND_RADIUS;	
            pos.z = (1.0 - 2.0*drand48()) * BND_RADIUS;	
          }
          const vec3 vel(0.0, 0.0, 0.0);
          Particle p;
          p.set_pos(pos);
          p.vel = vel;
          p.orig_vel = p.vel;
          p.boundary = 0;
          p.idx = -(pc + 1) ;
          const real rmax = dr_cld;
          p.rmax = rmax;

          ptcl_local.push_back(p);
          pc++;
          if (pc%100000 == 0)
            fprintf(stderr, "bnd: i= %d \n", pc);
        }
      }
#endif

#if 0
      {
        int pc = 0;
        vec3 pos(HUGE/100);
        while (pc < Ncld)
        {
          while(pos.abs() > RCLOUD)
          {
            pos.x = (1.0 - 2.0*drand48()) * RCLOUD;	
            pos.y = (1.0 - 2.0*drand48()) * RCLOUD;	
            pos.z = (1.0 - 2.0*drand48()) * RCLOUD;	
          }
          pos.x += XCL;
          pos.y += YCL;
          const vec3 vel(0.0, 0.0, 0.0);
          Particle p;
          p.set_pos(pos);
          p.vel = vel;
          p.orig_vel = p.vel;
          p.boundary = 0;
          p.idx = pc + 1;
          const real rmax = dr_cld;
          p.rmax = rmax;

          ptcl_local.push_back(p);
          pc++;
          if (pc%100000 == 0)
            fprintf(stderr, "cloud: i= %d \n", pc);
        }
      }
#endif

#if 1
      {
        int pc = 0;

        const real Rmax = Len3.abs();
        const real eps  = BND_RADIUS1;
        const real pow  = 3.0/2.0;
        assert(pow > 1.0);
        const real r0   = eps/std::sqrt(pow - 1.0);
        const real n3rmax = r0*r0/std::pow(r0*r0 + eps*eps, pow);

        while (pc < Namb)
        {
          bool pick = false;
          real r = eps;
          while (!pick)
          {
            r   = Rmax*drand48();
            //            if (r < BND_RADIUS0) continue;
            const real n3  = 1.0/std::pow(r*r + eps*eps, pow);
            const real n3r = r*r*n3;
            pick = (n3rmax*drand48() < n3r);
          }
#if 0
          const real theta = acos(1.0 - 2.0*drand48());
#else
          real theta = 0;
          pick = false;
          while (!pick)
          {
            theta = drand48()*M_PI;
            const real sigma = (M_PI/2.0) * HoverR*2.0;
            const real f= 0.05;
            const real p = (f + (1.0 - f)*exp(-sqr(theta - M_PI/2.0)/2/sqr(sigma)))*sin(theta);
            pick = (drand48() < p);
          }
#endif
          const real phi   = 2*M_PI*drand48();
          const vec3 pos(r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), r*cos(theta));
          if (pos.x < global_domain.rmin.x) continue;
          if (pos.x > global_domain.rmax.x) continue;
          if (pos.y < global_domain.rmin.y) continue;
          if (pos.y > global_domain.rmax.y) continue;
          if (pos.z < global_domain.rmin.z) continue;
          if (pos.z > global_domain.rmax.z) continue;
          const vec3 vel(0.0, 0.0, 0.0);
          Particle p;
          p.set_pos(pos);
          p.vel = vel;
          p.orig_vel = p.vel;
          p.boundary = 0;
          p.idx = -(Ncld*10 + ptcl_local.size() + pc + 1);
          const real rmax = 0.1*r;
          p.rmax = rmax;

          ptcl_local.push_back(p);
          pc++;
          if (pc%100000 == 0)
            fprintf(stderr, "ambient: i= %d \n", pc);
        }
      }
#endif



      local_n  = ptcl_local.size();
      global_n = local_n;

      fprintf(stderr, "  *** proc= %d : local_n= %llu  global_n= %llu \n", myproc, local_n, global_n);
    } // myproc == 0

    MPI_Bcast(&global_n,  1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    if (myproc == 0)
      fprintf(stderr, " ***  Distrubiting data \n");

    all_active = true;

    for (int k = 0; k < 5; k++)
      distribute_data(false,false);

    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(stderr, " *** proc= %d : local_n= %llu  global_n= %llu \n", myproc, local_n, global_n);
    MPI_Barrier(MPI_COMM_WORLD);

    if (myproc == 0)
      fprintf(stderr, " --- relax -- \n");
    relax_mesh(3);
    if (myproc == 0)
      fprintf(stderr, " ---- done --- \n");
    {
      distribute_data(false, false);
      const double t10 = mytimer::get_wtime();
      clear_mesh(false);
      int nattempt = build_mesh_global();
      double dt10 = mytimer::get_wtime() - t10;

      double volume_loc = 0.0;
      {
        std::vector<TREAL> v(local_n);
        for (int i = 0; i < (int)local_n; i++)
          v[i] = cell_local[i].Volume;
        std::sort(v.begin(), v.end());  // sort volumes from low to high, to avoid roundoff errors
        for (int i = 0; i < (int)local_n; i++)
          volume_loc += v[i];
      }

      double dt10max;
      MPI_Allreduce(&dt10, &dt10max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      double volume_glob = 0.0;	
      int    nattempt_max, nattempt_min;
      MPI_Allreduce(&volume_loc, &volume_glob,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&nattempt,   &nattempt_max, 1, MPI_INT,    MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&nattempt,   &nattempt_min, 1, MPI_INT,    MPI_MIN, MPI_COMM_WORLD);

      const double volume_exact = global_domain_size.x*global_domain_size.y*global_domain_size.z;
      if (myproc == 0)
      {
        fprintf(stderr, "first call build_mesh:[ %g  sec ::  %g cells/s/proc/thread ]\n",
            dt10max,
            global_n/nproc/dt10max);
        fprintf(stderr, "   computed_volume= %g  exact_volume= %g diff= %g [ %g ]  nattempt= %d %d \n",
            volume_glob, volume_exact, 
            volume_glob - volume_exact,	(volume_glob - volume_exact)/volume_exact,
            nattempt_min, nattempt_max);
      }
    }


    extract_ngb_from_mesh();
  }

  void system::set_problem(const bool init) 
  {

    if (myproc == 0) 
    {
      fprintf(stderr, " ********* Cloud capture ************* \n");
      fprintf(stderr, " **** MAGNETISATION=%g     \n", MAGNETISATION);
      //			fprintf(stderr, " **** COOLINGRATE  =%g     \n", COOLINGRATE  );
      fprintf(stderr, " **** HoverR       =%g     \n", HoverR  );
      fprintf(stderr, " **** GRAVITY_MASS= %g     \n", GM);
      //			fprintf(stderr, " **** GRAVITY_RACC= %g     \n", GRAVITY_RACC);
      fprintf(stderr, " **** GRAVITY_EPS=  %g     \n", GM_EPS);
      //			fprintf(stderr, " **** TCLOUD=       %g     \n", TCLOUD);
      //			fprintf(stderr, " **** DEREF_SIZE=   %g     \n", DEREF_SIZE);
      fprintf(stderr, " ---  \n");	
    }



    gamma_gas  = 1.0;
    courant_no = 0.8;

    const real xcl = XCL;
    const real ycl = YCL;

    const real vx_cl = VX0/VUNIT;
    const real vy_cl = VY0/VUNIT;

    const real vorb  = std::sqrt(sqr(vx_cl) + sqr(vy_cl));
    const real Rinit = std::sqrt(sqr(xcl)   + sqr(ycl)  );

    const real tinfall = Rinit/vorb;

    const real dcl = (DCLOUD/DUNIT);
    const real cs2 = get_cs2(vec3(xcl, ycl, 0.0)); //TCLOUD*Tunit/sqr(vunit);

    if (myproc == 0)
      fprintf(stderr, "Ro= %g pc, x= %g  y= %g; vx= %g vy= %g vt= %g  tinfall= %g Myr [%g]\n",
          Rinit,
          xcl, ycl,
          vx_cl, vy_cl,
          vorb,
          tinfall * TIMEUNIT, tinfall);

    if (!init) return;

    W_st.resize(local_n);
    W_rec.resize(local_n);

    for (int i = 0; i < (int)local_n; i++) 
    {
      Particle &pi = ptcl_local[i];

      const real R  = pi.pos.abs();
      assert(R > 0.0);

      if (pi.pos.abs() < BND_RADIUS)
        pi.boundary = Particle::DIOD;
      else
        pi.boundary = 0;

      const real Rdist = (pi.pos - vec3(xcl, ycl, 0.0)).abs();
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
        vx = vy = vz = 1.0e-4*std::sqrt(get_cs2(pi.pos));
        scalar = -1.0;
      }


      Fluid m;

      m[Fluid::DENS] = dens;
      m[Fluid::ETHM] = get_cs2(pi.pos)*dens;
      m[Fluid::VELX] = vx;
      m[Fluid::VELY] = vy;
      m[Fluid::VELZ] = vz;
      m[Fluid::BX  ] = bx;
      m[Fluid::BY  ] = by;
      m[Fluid::BZ  ] = bz;
      m[Fluid::PSI ] = 0.0;
      m[Fluid::ENTR] = 1.0;
      for (int k = 0 ; k < Fluid::NSCALARS; k++)
        m.scal(k) = scalar;

      W_st[i].w  = m;
      W_rec[i].w = m;
      Ulocal [i] = m.to_conservative(cell_local[i].Volume);
      dUlocal[i] = 0.0;
      ptcl_local[i].volume = cell_local[i].Volume;
      W_st [i].bnd = pi.boundary;
      W_rec[i].bnd = pi.boundary;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (myproc == 0)
      fprintf(stderr , " pvel ... \n");

    get_active_ptcl(true);

    MPI_Barrier(MPI_COMM_WORLD);
    if (myproc == 0)
      fprintf(stderr , " pvel ... \n");

    cell_list.swap(cell_local);
    ptcl_import.swap(ptcl_local);
    site_active_list.swap(active_ptcl);

    for (int i = 0; i < (int)local_n; i++)
    {
      real gpot;
      ptcl_import[i].acc = gpot_acc(i, gpot);
    }


    compute_pvel();
    compute_timesteps(true);

    cell_list.swap(cell_local);
    ptcl_import.swap(ptcl_local);
    site_active_list.swap(active_ptcl);

    for (int i = 0; i < (int)local_n; i++)
    {
      ptcl_local[i].rung += 3;
      ptcl_local[i].tend  = 0.0 + scheduler.get_dt(ptcl_local[i].rung);
      ptcl_local[i].vel = project(ptcl_local[i].orig_pos, ptcl_local[i].vel);
      ptcl_local[i].orig_vel = ptcl_local[i].vel;
    }
    all_active = true;
    scheduler.flush_list();
    boundary_n = 0;
    for (int i = 0; i < (int)local_n; i++)
    {
      if (ptcl_local[i].boundary == 0)
        scheduler.push_particle(i, (int)ptcl_local[i].rung);
      else
        boundary_n++;
      ptcl_local[i].unset_active();
    }
    

    unsigned long long boundary_glb;
    MPI_Allreduce(&boundary_n, &boundary_glb, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

    if (myproc == 0)
      fprintf(stderr, "boundary_glb= %lld\n", boundary_glb);

    clear_mesh(true);

    MPI_Barrier(MPI_COMM_WORLD);
    if (myproc == 0) fprintf(stderr, " proc= %d: complete problem setup \n", myproc);
  }
}

