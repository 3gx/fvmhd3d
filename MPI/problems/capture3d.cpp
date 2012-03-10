#include "fvmhd3d.h"

namespace fvmhd3d 
{

#define GM               35.0
#define GM_EPS           0.005

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

#define HoverR        0.1
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



  bool system::compute_pvel1()
  {
    const int nactive_site = site_active_list.size();
    for (int isite = 0; isite < nactive_site; isite++)
    {
      const int i = site_active_list[isite];
      
      Particle &p = ptcl_import[i];

			const Fluid W = U_import[i].to_primitive(ptcl_import[i].volume);
      
      p.vel = vec3(W[Fluid::VELX], W[Fluid::VELY], W[Fluid::VELZ]);

      const real B2   = sqr(W[Fluid::BX]) + sqr(W[Fluid::BY]) + sqr(W[Fluid::BZ]);
      const real pres = compute_pressure(W);

      const real cs2  = (gamma_gas * pres + B2)/W[Fluid::DENS];
      const real vel2 = vec3(W[Fluid::VELX], W[Fluid::VELY], W[Fluid::VELZ]).norm2();
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
      
      if (p.pos.abs() < BND_RADIUS) p.vel = 0.0;
      if (p.is_boundary())          p.vel = 0.0;
		
      if (p.idx > 0 && p.is_boundary() && p.pos.abs() > 0.25*BND_RADIUS && dt_global > 0.0)
      {
        const real dt = 32.0*dt_global;
        const real h = std::pow(p.volume/(4.0*M_PI/3.0), 1.0/3.0);
        p.vel = -p.pos * (h/dt * 1.0/p.pos.abs())*drand48();
      }

    }
    return true;
  }

  std::pair<vec3, real> gacc(const vec3 &pos)
  {
    const real ds2  = pos.norm2();
    const real ds   = std::sqrt(ds2);
    const real ids  = (ds > 0.0) ? 1.0/ds : 0.0;
    const real ids2 = ids*ids;

    const real gravity_eps  = GM_EPS;
    const real gravity_mass = GM;
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

  void system::predict(const int i)
  {
    assert(i >= 0);
    Particle &p = ptcl_local[i];
#if 1
    const real dth = 0.5*(p.tend - p.tlast);
    const std::pair<vec3, real> fh = gacc(p.orig_pos + p.orig_vel*dth);
    p.acc0 = p.acc1 = fh.first;

    const real dt     = t_global - p.tlast;
    p.vel = p.orig_vel; 
    p.pos = p.orig_pos + p.orig_vel*dt;
#else
    const std::pair<vec3, real> f0 = gacc(p.orig_pos);
    p.acc0 = f0.first;

    const real dt  = t_global - p.tlast;
    const real fac = p.is_boundary() ? 0.0 : 0.5;
    p.pos = p.orig_pos + p.orig_vel*dt + p.acc0*dt*dt*fac;
    p.vel = p.orig_vel; 

    const std::pair<vec3, real> f1 = gacc(p.pos);
    p.acc1 = f1.first;
#endif

    if (p.is_active() && (p.pos.abs() < BND_RADIUS))
    {
      if (!p.is_boundary()) boundary_n++;
      p.boundary = Particle::DIOD;
    }
  }

  void system::correct(const int i0)
  {
    const int i = i0 < 0 ? -1-i0 : i0;
    Particle &p = ptcl_local[i];
    {
      if (p.is_boundary())
      {
        p.tlast = scheduler.get_tsys();
        p.tend  = scheduler.get_tsys() + scheduler.get_dt(p.rung);

        p.unset_remove();
        p.orig_pos = p.pos;
        p.orig_vel = p.vel;
        if (p.idx > 0 && p.pos.abs() < 0.975*BND_RADIUS) p.set_remove();
        problem_set_boundary(i, p.volume);
        return;
      }
      else
        assert(p.pos.abs() >= BND_RADIUS);

#if 1
      const real dth = 0.5*(t_global - p.tlast);
      const vec3 ipos = p.orig_pos + (p.vel + p.orig_vel)*dth;
      const real r0 = (ipos - p.pos).abs();
      const real h = std::pow(p.volume/(4.0*M_PI/3.0), 1.0/3.0);
      if (r0 < 0.25*h) p.orig_pos = ipos;
      else             p.orig_pos = p.pos;
#else
      p.orig_pos = p.pos;
#endif
      p.orig_vel = p.vel;
      const std::pair<vec3, real> f = gacc(p.orig_pos);
      p.acc0 = p.acc1 = f.first;
      p.pot = f.second;

      p.tlast = scheduler.get_tsys();
      p.tend  = scheduler.get_tsys() + scheduler.get_dt(p.rung);
    }
  }

  void system::problem_set_boundary(const int i, const real v)
  {
    Particle &p = ptcl_local[i];
    if (!p.is_boundary()) return;
    assert(i == p.local_id);
    p.volume = v;
    U_local[p.local_id][Fluid::MOMX] = 0.0;
    U_local[p.local_id][Fluid::MOMY] = 0.0;
    U_local[p.local_id][Fluid::MOMZ] = 0.0;
    U_local[p.local_id][Fluid::WBX] = 0.0;
    U_local[p.local_id][Fluid::WBY] = 0.0;
    U_local[p.local_id][Fluid::WBZ] = 0.0;
    U_local[p.local_id][Fluid::DENS] = DENSvac/DUNIT/100.0*p.volume;
    U_local[p.local_id][Fluid::ETHM] = DENSvac/DUNIT/100.0*p.volume;
    return;
  }


  bool system::compute_update_prob(Fluid &U, const int i0)
  {
#if 1
    const int i = (i0 < 0) ? -1-i0 : i0;

#if 1
    Fluid W = (i0 >= 0) ? U.to_primitive(cell_list[i].Volume) : U;
#endif

    if (i0 < 0)
    {
      const Fluid &W = U;

      const real DFLOOR = 3.33*DENSvac/DUNIT;
      if (
          W[Fluid::DENS] <= 1.001*DFLOOR
          || ptcl_import[i].pos.abs() > 10.0 
          || ptcl_import[i].pos.abs() < 1.3*BND_RADIUS
         )
      {
        Wrec_import[i].t = 0.0;
        Wrec_import[i].x = 0.0;
        Wrec_import[i].y = 0.0;
        Wrec_import[i].z = 0.0;
      }
    }

#if 0
    U = (i0 >= 0) ? W.to_conservative(cell_list[i].Volume) : W;
#endif
#endif

    return false;
  }

  const real system::extra_timestep(const int i)
  {
    const Particle &p = ptcl_import[i];
    if (p.is_boundary()) return HUGE;
#if 0
    const std::pair<vec3, real> f0 = gacc(ptcl_import[i].pos);
    return 
      std::min(
          ptcl_import[i].pos.abs()/(1.0/HUGE + ptcl_import[i].vel.abs()),
          std::sqrt(2.0*ptcl_import[i].pos.abs()/(1.0/HUGE + f0.first.abs())));
#else
    const real r  = p.pos.abs();
    //    const real vr = std::abs(p.pos*p.vel)/(1.0/HUGE + r);
    return courant_no*r/(1.0/HUGE + p.vel.abs());
#endif
  }

  inline real get_r2(const vec3 &r)
  {
#if 1
    return r.norm2();
#else
    return sqr(r.x) + sqr(r.y) + sqr(r.y);
#endif
  }
  inline real get_cs2(const vec3 &r)
  {
#if 0
    return (GM/std::sqrt(get_r2(r)))*sqr(HoverR);
#else
    return std::abs(gacc(r).second)*sqr(HoverR);
#endif
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

  bool system::problem_force_distribute()
  {
#if 1
    unsigned long long boundary_glb;
    MPI_Allreduce(&boundary_n, &boundary_glb, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    return (boundary_glb > 20000);
#else
    return false;
#endif
  }

  void system::set_geometry(const bool init) 
  {
    const double dt_max = 1.0/128;
    scheduler = Scheduler(dt_max);

    t_end   = 5.0 + 1.0/65536.0;

    n_restart  = 1;
    dt_restart = 1.0/32;

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
//      case 6: nt.x = 3; nt.y = 2; nt.z = 1; break;
      case 8: nt.x = 4; nt.y = 2; nt.z = 1; break;
//      case 10: nt.x = 5; nt.y = 2; nt.z = 1; break;
//      case 12: nt.x = 4; nt.y = 3; nt.z = 1; break;
//      case 14: nt.x = 7; nt.y = 2; nt.z = 1; break;
      case 16: nt.x = 4; nt.y = 4; nt.z = 1; break;
      case 32: nt.x = 8; nt.y = 4; nt.z = 1; break;
      case 64: nt.x = 8; nt.y = 8; nt.z = 1; break;
      case 128: nt.x = 8; nt.y = 8; nt.z = 2; break;
      case 256: nt.x = 8; nt.y = 8; nt.z = 4; break;
      case 512: nt.x = 16; nt.y = 8; nt.z = 4; break;
      default: assert(false);
    }

    const Distribute::int3 nt_glb(nt);
    const pBoundary pglobal_domain(pfloat3(0.0), pfloat3(Len3));
    distribute_glb.set(nproc, nt, pglobal_domain);

    if (!init) return;

    if (myproc == 0) 
    {

      int Ncld = 1.5e5;
      int Namb = 0.5e5;

#if 1
      Ncld = 1.5e6;
      Namb = 0.5e6;
#endif

#if 0
      Ncld = 1.5e7;
      Namb = 0.5e7;
#endif

      int Nbnd = 512;


      ptcl_local.clear();
      ptcl_local.reserve(128);

      fprintf(stderr, "Ncloud= %d \n", Ncld);
      const real dr_amb = std::pow(Lbox*Lbox*Lbox/Namb,      1.0/3.0);
      const real dr_cld = 2.0*RCLOUD*std::pow(1.0/Ncld, 1.0/3.0);
      fprintf(stderr, "dr_amb=  %g   dr_cloud= %g\n", dr_amb, dr_cld);
      fprintf(stderr, "rmin= %g %g %g \n", 
          global_domain.get_rmin().x,
          global_domain.get_rmin().y,
          global_domain.get_rmin().z);
      fprintf(stderr, "rmax= %g %g %g \n", 
          global_domain.get_rmax().x,
          global_domain.get_rmax().y,
          global_domain.get_rmax().z);
      fprintf(stderr, "\n -------------------------- \n\n");

#if 1
      {
        int pc = 0;
        vec3 pos(HUGE/100);
        while (pc < Nbnd)
        {
          pos.x = (1.0 - 2.0*drand48()) * BND_RADIUS;
          pos.y = (1.0 - 2.0*drand48()) * BND_RADIUS;
          pos.z = (1.0 - 2.0*drand48()) * BND_RADIUS;
          const vec3 vel(0.0, 0.0, 0.0);
          Particle p;
          p.set_pos(pos);
          p.vel = vel;
          p.orig_vel = p.vel;
          p.boundary = 0;
          p.idx = pc+1;
          const real rmax = 1.5*BND_RADIUS/std::pow(Nbnd, 1.0/3.0);
          p.rmax = rmax;

          ptcl_local.push_back(p);
          pc++;
          if (pc%100000 == 0)
            fprintf(stderr, "bnd: i= %d \n", pc);
        }
      }
#endif

#if 1
      {
        int pc = 0;
        vec3 pos(HUGE/100);
        while (pc < Namb)
        {
          pos.x = (1.0 - 2.0*drand48()) * Lbox/2.0;
          pos.y = (1.0 - 2.0*drand48()) * Lbox/2.0;
          pos.z = (1.0 - 2.0*drand48()) * Lbox/2.0;
          const vec3 vel(0.0, 0.0, 0.0);
          Particle p;
          p.set_pos(pos);
          p.vel = vel;
          p.orig_vel = p.vel;
          p.boundary = 0;
          p.idx = 1000000000 + pc+1;
          const real rmax = 1.5*Lbox/std::pow(Namb, 1.0/3.0);
          p.rmax = rmax;

          ptcl_local.push_back(p);
          pc++;
          if (pc%100000 == 0)
            fprintf(stderr, "amb: i= %d \n", pc);
        }
      }
#endif

#if 1
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
          p.idx = 2000000000 + pc+1;
          const real rmax = dr_cld;
          p.rmax = rmax;

          ptcl_local.push_back(p);
          pc++;
          if (pc%100000 == 0)
            fprintf(stderr, "cld: i= %d \n", pc);
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
      distribute_data(false,false,false);
		fit_vec(ptcl_local);

    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(stderr, " *** proc= %d : local_n= %llu  global_n= %llu \n", myproc, local_n, global_n);
    MPI_Barrier(MPI_COMM_WORLD);

    if (myproc == 0)
      fprintf(stderr, " --- relax -- \n");
    relax_mesh(3);
    if (myproc == 0)
      fprintf(stderr, " ---- done --- \n");
    {
      distribute_data(false, false, false);
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

    U_local.resize(local_n);
    dU_local.resize(local_n);
    Wrec_local.resize(local_n);

    for (int i = 0; i < (int)local_n; i++) 
    {
      Particle &pi = ptcl_local[i];

      const real R  = pi.pos.abs();
      assert(R > 0.0);

      if (pi.pos.abs() < BND_RADIUS)
      {
        pi.idx = 0;
        pi.boundary = Particle::DIOD;
      }
      else
      {
        pi.boundary = 0;
      }

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
        const real csig = std::sqrt(get_cs2(pi.pos));
        vx = (1 - 2.0*drand48()) * csig;
        vy = (1 - 2.0*drand48()) * csig;
        vz = (1 - 2.0*drand48()) * csig;
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

      Wrec_local[i]        = Fluid_rec(m);
      U_local  [i]       = m.to_conservative(cell_local[i].Volume);
      dU_local[i] = 0.0;
      ptcl_local[i].volume = cell_local[i].Volume;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (myproc == 0)
      fprintf(stderr , " pvel ... \n");

    get_active_ptcl(true);

    MPI_Barrier(MPI_COMM_WORLD);
    if (myproc == 0)
      fprintf(stderr , " pvel ... \n");


    cell_list.swap(cell_local);
    U_import.swap(U_local);
    ptcl_import.swap(ptcl_local);
    site_active_list.swap(active_ptcl);

    compute_pvel();
    compute_timesteps(true);

    cell_list.swap(cell_local);
    U_import.swap(U_local);
    ptcl_import.swap(ptcl_local);
    site_active_list.swap(active_ptcl);

    for (int i = 0; i < (int)local_n; i++)
    {
      ptcl_local[i].rung += 3;
      ptcl_local[i].tend  = 0.0 + scheduler.get_dt(ptcl_local[i].rung);
      ptcl_local[i].orig_vel = ptcl_local[i].vel;
    }
    all_active = true;
    scheduler.flush_list();
    boundary_n = 0;
    for (int i = 0; i < (int)local_n; i++)
    {
      scheduler.push_particle(i, (int)ptcl_local[i].rung);
      if (ptcl_local[i].is_boundary())
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

