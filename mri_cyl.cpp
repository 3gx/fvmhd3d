#include "fvmhd3d.h"

/****** Problem specific methods ******/

namespace fvmhd3d
{

#define GM     1.0
#define GM_EPS 0.1

#define BND_RADIUS0 1.0
#define BND_RADIUS1 8.0

#define RinBND  1.0
#define RoutBND 8.0
  
  inline real get_r2(const vec3 &pos) 
  {
    return sqr(pos.x) + sqr(pos.y);
  }
  inline real get_cs2(const vec3 r)
  {
#if 0
    return (GM/std::sqrt(get_r2(r)))*sqr(0.1);
#else
    return 0.01; //(GM/2.0)*sqr(0.1);
#endif
  }

  void Main::Problem_set_global_domain()
  {
    const char string[256] = "MRI in Cylinder";
    sprintf(problem_string, "%s", string);

    vec3 hLbox(8.3, 8.3, 0.5);
    const vec3 rmin(-hLbox);
    const vec3 rmax( hLbox);

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
      const real vel2 = (1.0*W.get_vel()).norm2();

      const real vabs = std::sqrt(cs2 + vel2);

      const vec3 centroid = cell_list[i].centroid - p.pos;
      const real d = centroid.abs();
      if (d == 0.0) continue;

      const real eta = 0.25f;
      const real ki  = 0.5f;

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
      const vec3 dv = centroid*f;
      p.vel += dv;

      if (p.is_boundary()) p.vel = 0.0;
#endif

      const vec3 Rc(p.pos.x, p.pos.y, 0.0);

      p.vel.z = 0.0;
      p.vel = p.vel - (p.vel*Rc)*Rc/Rc.norm2();  // subtract radial component

#if 0
      p.vel = 0.0;
#endif
    }

    return true;
  }

  void System::Problem_generate_geometry(const int param)
  {

    const double dt_max = 1.0;
    scheduler = Scheduler(dt_max);

    t_end      = 200;
    dt_restart = 1.0/64;
    dt_snap    = 1.0;

    dt_restart = std::max(dt_restart, dt_max);
    dt_snap    = std::max(dt_snap,    dt_max);


    int Nz = 24;

#if 0
    Nz = 64;
#endif

    const real HoR = 1.0;
    const int Ndisk_glb = sqr(4.0/HoR)*Nz*Nz*Nz;


    const int Ndisk = Ndisk_glb/numElements;
    const int Namb  = Ndisk/4;

    ptcl_list.clear();
    ptcl_list.reserve(Ndisk);

    Rand48 rnd;
    rnd.srand(123 + 123*thisIndex);

    {
      int pc = 0;
      vec3 pos(HUGE/100);
      const vec3 hLbox = global_domain_size/2.0;
      while (pc < Namb)
      {
        pos.x = (1.0 - 2.0*rnd.drand()) * hLbox.x;
        pos.y = (1.0 - 2.0*rnd.drand()) * hLbox.y;
        pos.z = (1.0 - 2.0*rnd.drand()) * hLbox.z;

        ptcl_list.push_back(Particle(ptcl_list.size(), thisIndex, pos));
        pc++;
      }
    }

    {
      int pc = 0;
      vec3 pos(HUGE/100);
      const vec3 hLbox = global_domain_size/2.0;
      while (pc < Ndisk)
      {
        pos.x = (1.0 - 2.0*rnd.drand()) * 4.0;
        pos.y = (1.0 - 2.0*rnd.drand()) * 4.0;
        const real R = std::sqrt(sqr(pos.x) + sqr(pos.y));
        if (R < 1.0 || R > 4.0) continue;

        pos.z = (1.0 - 2.0*rnd.drand()) * hLbox.z;

        ptcl_list.push_back(Particle(ptcl_list.size(), thisIndex, pos));
        pc++;
      }
    }

    local_n = ptcl_list.size();
    if (thisIndex == 0)
      CkPrintf("local_n= %d \n" ,local_n);

    generateGeometry_nRelax = 3;
  }

  void System::Problem_generate_IC(const int param)
  {
    if (thisIndex == 0)
    {
      fprintf(stderr, " ********* Setting up %s problem ************* \n", problem_string);
    }

    gamma_gas  = 5.0/3.0;
    courant_no = 0.8;

    t_global  = 0;
    iteration = 0;

    Rand48 rnd;
    rnd.srand(1023 + 123*thisIndex);

    for (int i = 0; i < local_n; i++) 
    {
      const Particle &pi = ptcl_list[i];

      const vec3 &pos = pi.get_pos();

      const real x = pos.x;
      const real y = pos.y;

      const real R = std::sqrt(get_r2(pos)); 
      assert(R > 0.0);

      if (R < BND_RADIUS0)
      {
        mesh_pnts[i].boundary = MeshPoint::OUTFLOW;
      }
      else if (R > RoutBND)
      { 
        mesh_pnts[i].boundary = MeshPoint::OUTFLOW;
      }
      else  
        mesh_pnts[i].boundary = MeshPoint::NO_BOUNDARY;

      const real vphi = std::sqrt(GM/R);
      real vx = -vphi * y/R;
      real vy = +vphi * x/R;
      real vz = 0.0;
#if 0
      if (R > RoutBND)
        vx = vy = vz = 0;
#endif

#if 1
      if (!mesh_pnts[i].is_boundary())
      {
        const real dv = vphi*1.0e-4;
        const real dvz = (1.0 - 2*rnd.drand()) * dv;
        const real dvr = (1.0 - 2*rnd.drand()) * dv;
        vx += dvr*x/R;
        vy += dvr*y/R;
        vz += dvz;
      }
#endif

      real bx = 0.0;
      real by = 0.0;
      real bz = 0.0;

      real B0 = 0.05513;
      //      B0 = 0.052848;
      real n = 2;
#if 1
      if (R > 2.0 && R < 4)
        bz = B0/n * sin(2*M_PI*(R-2.0));
#else
      // bz  = (tanh((R-1.5)/0.01) - tanh((R-3.5)/0.01))*0.5;
      bz  = (tanh((R-1.2)/0.1) - tanh((R-3.8)/0.1))*0.5;
      bz *= B0/n;
#endif
      if (mesh_pnts[i].is_boundary())
        bx = by = bz = 0.0;

      real dens = 1;
      real pres = get_cs2(pos)*dens/gamma_gas;


      Fluid m;

      m[Fluid::DENS] = dens;
      m[Fluid::ETHM] = pres;
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
    const real ds2  = get_r2(pos);
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

    return std::make_pair(acc * vec3(pos.x, pos.y, 0.0), pot);
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

#if 1
    const std::pair<vec3, real> f0 = gacc(p.pos_orig);
    p.acc0 = f0.first;

    const real dt = t_global - p.tbeg;
    p.pos = p.pos_orig;
    p.vel = p.vel_orig; 

    if (p.pos.cross(p.vel).norm2() > 0.0)
    {
      assert(p.vel.z == 0.0);
      const std::pair<vec3,vec3> rv1 = rotate(vec3(p.pos.x, p.pos.y, 0.0), p.vel, dt);
      assert(rv1.first.z == 0.0);
      p.pos = vec3(rv1.first.x, rv1.first.y, p.pos.z);
      //      p.vel += p.acc0*dt;
    }

    const std::pair<vec3, real> f1 = gacc(p.pos);
    p.acc1 = f1.first;
    //    p.acc0 = p.acc1 = (p.acc0 + p.acc1)*0.5;
#endif

#if 0
    const std::pair<vec3, real> f0 = gacc(p.pos_orig);
    p.acc0 = f0.first;

    const real dt     = t_global - p.tbeg;
    p.pos = p.pos_orig;
    p.vel = p.vel_orig; 
    if (!p.is_boundary())
    {
#if 0
      p.pos = p.pos_orig + p.vel_orig * dt;
#else
      p.pos = p.pos_orig + p.vel_orig * dt + p.acc0 * dt*dt*0.5;
      p.vel = p.vel_orig + p.acc0 * dt;
#endif
    }

    const std::pair<vec3, real> f1 = gacc(p.pos);
    p.acc1 = f1.first;
    p.acc0 = p.acc1 = (p.acc0 + p.acc1)*0.5;
#endif
  }

  void System::Problem_correct_meshpoint_position(const int Id)
  {
    MeshPoint &p = mesh_pnts[Id];
    p.pos_orig = p.pos;
    p.vel_orig = p.vel;

#if 0
    if (!p.is_boundary())
    {
      const real dt     = t_global - p.tbeg;
      p.pos_orig = p.pos + p.acc0*dt*dt*0.5;
    }
#endif
  }

  bool System::Problem_compute_update(Fluid &Uc, const int Id)
  {
    return false;
  }

  real System::Problem_extra_timestep_criterion(const int Id)
  {
    return HUGE;
  }

#if 1
#endif
  real System::Problem_compute_ethm_update(const Fluid &W, const int i)
  {
#if 0
    return W[Fluid::ETHM];
    //    return Problem_compute_ethm_from_entropy(W);
#else
    return get_cs2(mesh_act[i]->pos)*W[Fluid::DENS]/gamma_gas;
#endif
  }

  real System::Problem_compute_pressure(const Fluid &W)
  {
#if 0
    return gamma_gas > 1.0 ? (gamma_gas - 1.0) * W[Fluid::ETHM] : W[Fluid::ETHM];
#else
    return W[Fluid::ETHM];
#endif
  }

  real System::Problem_compute_entropy_from_ethm(const Fluid &W)
  {
#if 0
    assert(gamma_gas > 1.0);
    return (gamma_gas - 1.0) * W[Fluid::ETHM]/std::pow(W[Fluid::DENS], gamma_gas);
#else
    return 1.0;
#endif
  }

  real System::Problem_compute_ethm_from_entropy(const Fluid &W)
  {
#if 0
    assert(gamma_gas > 1.0);
    return W[Fluid::ENTR] * std::pow(W[Fluid::DENS], gamma_gas)/(gamma_gas - 1.0);
#else
    assert(false);
    return -1.0;
#endif
  }

  real System::Problem_enforce_limiter(const int i)
  {
    const vec3 &pos = mesh_act[i]->pos;

    const real R = std::sqrt(get_r2(pos)); 
    if (R < 1.1) return 0.0;
    else         return 1.0;
  }

  void System::Problem_set_boundary(const int i)
  {
    const vec3 &pos = mesh_act[i]->pos;

    const real x = pos.x;
    const real y = pos.y;

    const real R = std::sqrt(get_r2(pos)); 
    const real vphi = std::sqrt(GM/R);
    real vx = -vphi * y/R;
    real vy = +vphi * x/R;
#if 0
    const real vr = R <= 2 ? -std::sqrt(get_cs2(pos)) : 0.0;
    vx += vr * x/R;
    vy += vr * y/R;
#endif
    Wrec_act[i]->w[Fluid::VELX] = vx;
    Wrec_act[i]->w[Fluid::VELY] = vy;
    Wrec_act[i]->w[Fluid::VELZ] = 0.0;
  }

  bool System::Problem_meshpoint_refine(const int i)
  {
    return false;
  }

  bool System::Problem_meshpoint_derefine(const int i)
  {
    return false;
  }




}
