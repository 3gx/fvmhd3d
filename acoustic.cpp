#include "fvmhd3d.h"

/****** Problem specific methods ******/

namespace fvmhd3d
{

#if 1
#define NX 16
#define NY 32
#define NZ 16
#endif


  void Main::Problem_set_global_domain()
  {
    const char string[256] = "Acoustic wave";
    sprintf(problem_string, "%s", string);
    
		const double Lx = 1.0;
		const vec3 rmin(0.0);
		const vec3 rmax((Lx/NY)*NX, Lx, (Lx/NY)*NZ);

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
#endif

#if 0
      p.vel = 0.0;
#endif
    }

    return true;
  }

  void System::Problem_generate_geometry(const int param)
  {

#if 0
    const double dt_max = 1.0/32;
#else
    const double dt_max = 1.0/64;
#endif
    scheduler = Scheduler(dt_max);

    t_end      = 5.0;
    dt_restart = 1.0/64;
    dt_snap    = dt_max;

    dt_restart = std::max(dt_restart, dt_max);
    dt_snap    = std::max(dt_snap,    dt_max);
      

    ptcl_list.clear();

    if (thisIndex == 0)
    {

      const int nglob = NX*NY*NZ;
      ptcl_list.reserve(nglob);

      const vec3 Len3 = global_domain.hsize() * 2.0;
      const dvec3 dr = dvec3(Len3.x/NX, Len3.y/NY, Len3.z/NZ);
      for (int k = 0; k < NZ; k++) 
        for (int j = 0; j < NY; j++) 
          for (int i = 0; i < NX; i++) 
          {
            dvec3 pos = global_domain.get_rmin() + dvec3(i*dr.x, j*dr.y, k*dr.z) + 0.5*dr;
#if 1
            {
              const real f = 1.0e-6;
              pos += vec3(drand48()*dr.x*f, drand48()*dr.y*f, drand48()*dr.z*f);
            }
#endif


#if 1
            pos = global_domain.get_rmin() + dvec3(
                drand48()*Len3.x,
                drand48()*Len3.y,
                drand48()*Len3.z);
#else
#define _UNIFORM_MESH_
#endif
            ptcl_list.push_back(Particle(ptcl_list.size(), thisIndex, pos));

#if 0
#ifndef _UNIFORM_MESH_
            pos = global_domain.get_rmin() + vec3(0.0, Len3.y/8, 0.0) + dvec3(
                drand48()*Len3.x,
                drand48()*Len3.y/4.0,
                drand48()*Len3.z);
            ptcl_list.push_back(Particle(ptcl_list.size(), thisIndex, pos));
#endif
#endif
          }
    }

    local_n = ptcl_list.size();

#ifndef _UNIFORM_MESH_
    generateGeometry_nRelax = 7;
#else
    generateGeometry_nRelax = 1;
#endif
  }

  const real cs2 = 1.0;

  void System::Problem_generate_IC(const int param)
  {
    if (thisIndex == 0)
    {
      CkPrintf(" ********* Setting up %s  ************* \n", problem_string);
      CkPrintf(" NX= %d  NY= %d  NZ= %d \n", NX, NY, NZ);
    }

    gamma_gas  = 1.0;
    courant_no = 0.4;

    t_global  = 0;
    iteration = 0;

    for (int i = 0; i < local_n; i++) 
    {
      const Particle &pi = ptcl_list[i];

      const vec3 &pos = pi.get_pos();

      const real y = pos.y;
			
      real d, p, vx, vy, vz, bx, by, bz;
      real scal = 1.0;

			vx = 0;
			vy = 0;
			vz = 0;
			bx = by = bz = 0;
			p = 1;
			d = 1;
			
      const real ampl = 0.0001;
			const real kwave = 2.0 * 2*M_PI;
			p  = 1.0; 
			d  = 1 + ampl * sin(kwave*y);
			vy = -gamma_gas * ampl * sin(kwave*y);

			scal = d;


      Fluid m;

      m[Fluid::DENS] = d;
      m[Fluid::ETHM] = cs2*d;
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

  void System::Problem_predict_meshpoint_position(const int Id)
  {
    MeshPoint &p = mesh_pnts[Id];

    const real dt     = t_global - p.tbeg;
    p.vel = p.vel_orig; 
    p.pos = p.pos_orig + p.vel_orig*dt;
    
    p.acc0 = p.acc1 = 0.0;
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
    return cs2*W[Fluid::DENS];
  }

  real System::Problem_compute_pressure(const Fluid &W)
  {
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

  real System::Problem_enforce_limiter(const int i)
  {
    return 1.0;
  }

  void System::Problem_set_boundary(const int i)
  {
  }






}
