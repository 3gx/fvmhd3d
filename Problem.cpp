#include "fvmhd3d.h"

/****** Problem specific methods ******/
#if 0
#include "mri_disk.cpp"
#endif

#if 0
#include "mri_cyl.cpp"
#endif

#if 1
#include "alfven.cpp"
#endif

#if 0
#include "acoustic.cpp"
#endif

#if 0
#include "capture3d.cpp"
#endif

#if 0
#include "blast.cpp"
#endif


#if 0
namespace fvmhd3d
{


#if 0
#define __ADVECT_PULSE_TEST__
#endif

#if 1
#define NX 64
#define NY 64
#define NZ 16
#endif

#if 0
#define NX 128
#define NY 128
#define NZ 32
#endif

  void Main::Problem_set_global_domain()
  {
    const char string[256] = "Orszag-Tang vortex\n";
    sprintf(problem_string, "%s", string);

    const double Lx = 1.0;
    const vec3 rmin(0.0);
    const vec3 rmax(Lx, (Lx/NX)*NY, (Lx/NX)*NZ);

    global_domain = boundary(rmin, rmax);
    global_domain_size = global_domain.hsize() * 2.0;

  };

  bool System::Problem_computePvel()
  {
#if 1
    return false;
#else
    for (int i = 0; i < nactive_loc; i++)
      mesh_act[i]->vel = vec3(0.0);
    return true;
#endif
  }

  void System::Problem_generate_geometry(const int param)
  {

    const double dt_max = 1.0/64;
    scheduler = Scheduler(dt_max);

    t_end      = 2.5;
    dt_restart = 1.0/64;
    dt_snap    = 1.0/16;

    dt_restart = std::max(dt_restart, dt_max);
    dt_snap    = std::max(dt_snap,    dt_max);

#if 1
    dt_snap = dt_max;
#endif


    const int global_n = NX*NY*NZ;
    std::vector<int> n_per_element(numElements, 0);

    for (int i = 0; i < global_n; i++)
      n_per_element[i % numElements]++;

    local_n = n_per_element[thisIndex];

    ptcl_list.clear();
    ptcl_list.reserve(local_n);

    Rand48 rnd;
    rnd.srand(123 + 123*thisIndex);
    for (int i = 0; i < local_n; i++)
    {
      vec3 pos;
      pos.x = global_domain.get_rmin().x + rnd.drand() * global_domain_size.x;
      pos.y = global_domain.get_rmin().y + rnd.drand() * global_domain_size.y;
      pos.z = global_domain.get_rmin().z + rnd.drand() * global_domain_size.z;
      ptcl_list.push_back(Particle(i, thisIndex, pos));
    }

    generateGeometry_nRelax = 2;
  }

  void System::Problem_generate_IC(const int param)
  {
    const real b0 = 1.0/sqrt(4.0*M_PI);
    const real d0 = 25.0/(36.0*M_PI);
    const real v0 = 1.0;
    const real p0 = 5.0/(12*M_PI);

    gamma_gas     = 5.0/3;
    courant_no    = 0.4;

    t_global  = 0;
    iteration = 0;

    real adv = 0.0;

    for (int i = 0; i < local_n; i++) 
    {
      const Particle &pi = ptcl_list[i];

      const real x = pi.get_pos().x;
      const real y = pi.get_pos().y;

      real d, p, vx, vy, vz, bx, by, bz;

      vx = -v0 * sin(2.0*M_PI*y) + adv;
      vy = +v0 * sin(2.0*M_PI*x) + adv;
      vz =  0.0;

      bx = -b0*sin(2*M_PI*y);
      by = +b0*sin(4*M_PI*x);
      bz = 0.0;

      d = d0;
      p = p0;
      real scal = 1;

#ifdef __ADVECT_PULSE_TEST__
      vx = 0;
      vy = 0;
      vz = 0;
      bx = by = bz = 0;
      p = 1.0;
      d = 1.0;

      vx = 0; vy = 1;
      if (y > 0.25 && y < 0.75) 
      {
        d = 10;
      }
#endif


      Fluid m;

      m[Fluid::DENS] = d ;
      m[Fluid::ETHM] = p/(gamma_gas - 1.0);
      m[Fluid::VELX] = vx;
      m[Fluid::VELY] = vy;
      m[Fluid::VELZ] = vz;
      m[Fluid::BX  ] = bx;
      m[Fluid::BY  ] = by;
      m[Fluid::BZ  ] = bz;
      m[Fluid::PSI ] = 0.0;
      m[Fluid::ENTR] = Problem_compute_entropy_from_ethm(m);
      for (int k = 0 ; k < Fluid::NSCALARS; k++)
        m.scal(k) = scal;

      Wrec_list[i] = Fluid_rec(m);
      
      mesh_pnts[i].idx      = thisIndex*(local_n << 1)+ i + 1;
      mesh_pnts[i].boundary = MeshPoint::NO_BOUNDARY;
      
    }
  }

  void System::Problem_predict_meshpoint_position(const int Id)
  {
    MeshPoint &p = mesh_pnts[Id];
    assert(p.idx >= 0);
    p.acc0 = p.acc1 = 0.0;

    const real dt = t_global - p.tbeg;
    p.vel = p.vel_orig; 
    p.pos = p.pos_orig + p.vel_orig*dt;
  }

  void System::Problem_correct_meshpoint_position(const int Id)
  {
    MeshPoint &p = mesh_pnts[Id];
    const real dth = 0.5*(t_global - p.tbeg);

    const vec3 ipos = p.pos_orig + (p.vel + p.vel_orig)*dth;
    const real r0 = (ipos - p.pos).abs();
    const real h = std::pow(p.Volume/(4.0*M_PI/3.0), 1.0/3.0);
    if (r0 < 0.25*h) p.pos_orig = ipos;
    else             p.pos_orig = p.pos;

    p.set_vel(p.vel);
  }

  bool System::Problem_compute_update(Fluid &Uc, const int Id)
  {
    return false;
  }

  real System::Problem_extra_timestep_criterion(const int Id)
  {
    return HUGE;
  }

  real System::Problem_compute_ethm_update(const Fluid &W, const int Id)
  {
    //		return Problem_compute_ethm_from_entropy(W);
    return W[Fluid::ETHM];
  }

  real System::Problem_compute_pressure(const Fluid &W)
  {
    return gamma_gas > 1.0 ? (gamma_gas - 1.0) * W[Fluid::ETHM] : W[Fluid::ETHM];
  }

  real System::Problem_compute_entropy_from_ethm(const Fluid &W)
  {
    assert(gamma_gas > 1.0);
    return (gamma_gas - 1.0) * W[Fluid::ETHM]/std::pow(W[Fluid::DENS], gamma_gas);
  }

  real System::Problem_compute_ethm_from_entropy(const Fluid &W)
  {
    assert(gamma_gas > 1.0);
    return W[Fluid::ENTR] * std::pow(W[Fluid::DENS], gamma_gas)/(gamma_gas - 1.0);
  }




}
#endif
