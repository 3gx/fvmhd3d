#include "fvmhd3d.h"

namespace fvmhd3d 
{

#define GM 1.0
#define GM_EPS 0.01
#define BND_RADIUS 0.01

/******* disk mass & geometry ******/

#define Mf 0.1
#define R0 1.0
#define H0 0.05

/******* disk cs, beta & vphi ********/

#define beta 100
#define Vf 0.95
#define HoR 0.1

#define VOLDISK (4.0*sqr(M_PI)*R0*sqr(H0))

/******** the rest ********/

#define D0     (Mf*GM/VOLDISK)
#define DFLOOR (D0*1.0e-7)
#define DVAC   (D0*1.0e-4)

#if 1
#define iso_cs (HoR*std::sqrt(GM/R0))
#define B0     (std::sqrt(2.0*D0*sqr(iso_cs)/beta))
// #define B0     0.0
#endif

#ifndef iso_cs
  real iso_cs;
#endif


  bool system::compute_pvel1()
  {
    const int nactive_site = site_active_list.size();
    for (int isite = 0; isite < nactive_site; isite++)
    {
      const int i = site_active_list[isite];

      Particle &p = ptcl_import[i];

      const Fluid &W = W_st[i].w;
      
      p.vel = vec3(W[Fluid::VELX], W[Fluid::VELY], W[Fluid::VELZ]);

      // if density below DVAC, particle are fixed in space (modulo Lloyd regularization)
      const real D  = std::log(W[Fluid::DENS]);
      const real DV = std::log(DVAC);
      const real dD = DV * 0.01/2;
      const real r  = (D-DV)/std::abs(dD);
      const real vfac = (-r > 20) ? 0.0 : 1.0/(1.0 + std::exp(-r));
      assert(vfac >= 0.0);
      assert(vfac <= 1.0);
      p.vel *= vfac;



      const real B2   = sqr(W[Fluid::BX]) + sqr(W[Fluid::BY]) + sqr(W[Fluid::BZ]);
      const real pres = compute_pressure(W);

      const real cs2  = (gamma_gas * pres + B2)/W[Fluid::DENS];
      const real vel2 = vec3(W[Fluid::VELX], W[Fluid::VELY], W[Fluid::VELZ]).norm2();
      const real vabs = std::sqrt(cs2 + 1.0*vel2);

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

//      p.vel *= vfac;

      W_st[i].vel = p.vel;
    }
    return true;
  }

  void system::predict(const int i)
  {
    const real dt = t_global - ptcl_local[i].tlast;
    ptcl_local[i].vel = ptcl_local[i].orig_vel + 1.0*ptcl_local[i].acc*dt;
    ptcl_local[i].pos = 
      ptcl_local[i].orig_pos  +
      (ptcl_local[i].orig_vel +
       ptcl_local[i].acc*(dt*0.5))*dt;
    ptcl_local[i].pos = periodic(ptcl_local[i].pos);
  }

  bool system::compute_update_prob(Fluid &U, const int i)
  {
    if (ptcl_import[i].pos.abs() < BND_RADIUS)
    {
      const real V = cell_list[i].Volume;
      Fluid W = U.to_primitive(V);

      const real d      = DFLOOR*10;
      const real e   = sqr(iso_cs)*DFLOOR;

      real vx,vy,vz;
      real bx,by,bz;

      vx = vy = vz = 0;
      bx = by = bz = 0;

      W[Fluid::DENS] = d;
      W[Fluid::ETHM] = e;
      W[Fluid::VELX] = vx;
      W[Fluid::VELY] = vy;
      W[Fluid::VELZ] = vz;
      W[Fluid::BX  ] = bx;
      W[Fluid::BY  ] = by;
      W[Fluid::BZ  ] = bz;
      W[Fluid::PSI ] = 0;
      W[Fluid::ENTR] = 1.0;

      ptcl_import[i].vel = 0.0;
      ptcl_import[i].acc = 0.0;

      U = W.to_conservative(V);

      return true;
    }

    const real V = cell_list[i].Volume;
    Fluid W = U.to_primitive(V);
    W[Fluid::DENS] = std::max(W[Fluid::DENS], DFLOOR);
    W[Fluid::ETHM] = sqr(iso_cs) * W[Fluid::DENS];
    U = W.to_conservative(V);

    return false;
  }

  const vec3 system::gpot_acc(const int i, real &gpot) const
  {
    const vec3 pos  = ptcl_import[i].pos;
    const real ds2  = pos.norm2();
    const real ds   = std::sqrt(ds2);
    if (ds < BND_RADIUS)
    {
      gpot = 0.0;
      return vec3(0.0);
    }
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

  real system::compute_ethm_update(const Fluid &W, const int i) const
  {
    return sqr(iso_cs)*W[Fluid::DENS];
  }
  real system::compute_pressure(const Fluid& W) const
  {
    return sqr(iso_cs) * W[Fluid::DENS];
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
  real system::compute_pressure_gradient(const Fluid &W, const Fluid &dW) const
  {
    const real dpdethm = 0.0;
    const real dpddens = sqr(iso_cs);
    return dpddens * dW[Fluid::DENS] + dpdethm * dW[Fluid::ETHM];
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

    t_end   = 50.0 + dt_max;

    n_restart = 2;
    dt_restart = dt_max;

    dt_dump = dt_max / 64;
    dt_dump = 1.0/16;

    di_log = 100;

    global_n = local_n = 0;


    const double Lbox = 8;
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

      int Namb  = 1e5;
      int Ndisk = 1e5;

      ptcl_local.clear();
      ptcl_local.reserve(128);

      const real H1 = 2.0*H0;

      fprintf(stderr, "Namb=  %d \n", Namb);
      fprintf(stderr, "Ndisk= %d \n", Ndisk);
      const real dr_amb  = std::pow(Lbox*Lbox*Lbox/Namb,      1.0/3.0);
      const real dr_disk = std::pow(0.5*VOLDISK*sqr(H1/H0)/Ndisk, 1.0/3.0);
      fprintf(stderr, "dr_amb=  %g \n", dr_amb );
      fprintf(stderr, "dr_disk= %g \n", dr_disk);
      fprintf(stderr, "rmin= %g %g %g \n", 
          global_domain.rmin.x,
          global_domain.rmin.y,
          global_domain.rmin.z);
      fprintf(stderr, "rmax= %g %g %g \n", 
          global_domain.rmax.x,
          global_domain.rmax.y,
          global_domain.rmax.z);
      fprintf(stderr, "\n -------------------------- \n\n");

      for (int i = 0; i < Namb; i++)
      {
        const vec3 pos = global_domain.rmin + dvec3(
            drand48()*Len3.x,
            drand48()*Len3.y,
            drand48()*Len3.z);
        
        dvec3 vel(0.0, 0.0, 0.0);
        Particle p;
        p.set_pos(pos);
        p.vel = vel;
        p.orig_vel = p.vel;
        p.boundary = 0;
        p.idx = i + 1;
        p.rmax = 2.0*dr_amb;
        ptcl_local.push_back(p);
      }

      for (int i = 0; i < Ndisk; i++)
      {
        const real phi   = 2*M_PI*drand48();
        const real theta = 2*M_PI*drand48();
        const real s     = H1*std::sqrt(drand48());

        const real R = R0 + s*sin(theta);
        const real x = R * cos(phi);
        const real y = R * sin(phi);
        const real z = s * cos(theta);


        const vec3 pos(x,y,z);
        
        dvec3 vel(0.0, 0.0, 0.0);
        Particle p;
        p.set_pos(pos);
        p.vel = vel;
        p.boundary = 0;
        p.idx = 10*Namb + i + 1;
        p.rmax = 2.0*dr_disk;
        ptcl_local.push_back(p);
      }
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
    relax_mesh(5);
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

#ifndef iso_cs
    const real B0 = 0.2*std::sqrt(2.0*D0)*std::sqrt(GM/R0)*HoR*R0;
    iso_cs = std::sqrt(0.5*beta * sqr(B0)/D0);
#endif
    if (myproc == 0)
    {
      fprintf(stderr, " ********* Setting up tori_collapse problem ************* \n");
      fprintf(stderr, " GM_EPS= %g \n", GM_EPS);
      fprintf(stderr, " beta= %g \n", beta);
      fprintf(stderr, " Mf=  %g \n",Mf);
      fprintf(stderr, " Vf=  %g \n", Vf);
      fprintf(stderr, " HoR= %g \n", HoR);
      fprintf(stderr, " R0= %g \n", R0);
      fprintf(stderr, " H0= %g \n", H0);
      fprintf(stderr, " CS= %g \n", iso_cs);
      fprintf(stderr, " D0= %g \n", D0);
      fprintf(stderr, " B0= %g \n", B0);
    }


    gamma_gas  = 1.0;
    courant_no = 0.4;

    if (!init) return;

    W_st.resize(local_n);

    int n1 = 0;
    int n2 = 0;

    for (int i = 0; i < (int)local_n; i++) 
    {
      const Particle &pi = ptcl_local[i];

      real x = pi.pos.x;
      real y = pi.pos.y;
      real z = pi.pos.z;

      const real r = std::sqrt(x*x + y*y);
      assert(r > 0.0);

      const real s2 = sqr(r - R0) + z*z;
      const real q  = s2/(2.0*sqr(H0));
      
      const real dens = std::max(D0*std::exp(-q), 10*DFLOOR);
      const real bmag = (q > 10.0) ? 0.0 : B0*std::exp(-0.5*q);
      const real vphi = Vf * std::sqrt(GM/r);

      const real bx = -bmag*y/r;
      const real by = +bmag*x/r;
      const real vx = -vphi*y/r;
      const real vy = +vphi*x/r;

      const real bz = 0.0;
      const real vz = 0.0;

      real scalar = -1;
      n2++;
      if (s2 < sqr(2.0*H0))
      {
        n2--;
        scalar = 1;
        n1++;
      }


      Fluid m;

      m[Fluid::DENS] = dens;
      m[Fluid::ETHM] = sqr(iso_cs)*dens;
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
      Ulocal [i] = m.to_conservative(cell_local[i].Volume);
      dUlocal[i] = 0.0;
      ptcl_local[i].volume = cell_local[i].Volume;
#if 0
      real gpot;
      ptcl_local[i].acc = gpot_acc(i, gpot);
#endif

#if 0
      const double L      = std::pow(cell_local[i].Volume, 1.0/3);
      const double cs_est = std::sqrt((p*gamma_gas + (bx*bx+by*by+bz*bz))/d);
      const double v_est  = std::sqrt(vx*vx + vy*vy + vz*vz);
      const double dt_est = 0.1 * courant_no * L/(cs_est + v_est);

      ptcl_local[i].tlast = 0.0;
      ptcl_local[i].rung = scheduler.get_rung(dt_est);

      dt_min = std::min(dt_min, dt_est);
#endif
    }
    fprintf(stderr, "myproc= %d: n1= %d  n2= %d  local_n= %llu\n", myproc,
        n1, n2, local_n);
    assert(n1+n2 == (int)local_n);

#if 0
    double dt_min_glob;
    MPI_Allreduce(&dt_min, &dt_min_glob, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif

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
      ptcl_local[i].rung += 4;
      ptcl_local[i].tend  = 0.0 + scheduler.get_dt(ptcl_local[i].rung);
      ptcl_local[i].orig_vel = ptcl_local[i].vel;
    }
    all_active = true;
    scheduler.flush_list();
    for (int i = 0; i < (int)local_n; i++)
      scheduler.push_particle(i, (int)ptcl_local[i].rung);

    clear_mesh(true);

    MPI_Barrier(MPI_COMM_WORLD);
    if (myproc == 0) fprintf(stderr, " proc= %d: complete problem setup \n", myproc);
  }
}

