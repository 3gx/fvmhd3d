#include "fvmhd3d.h"

/****** Problem specific methods ******/

namespace fvmhd3d
{

#define GM     1.0
#define GM_EPS 0.1
#define HoR    0.1

#if 1
#define BND_RADIUS0 0.5
#define BND_RADIUS1 1.0
#define Rin  1.0
#define B_Rin  2.0
#define B_Rout 3.0
#define BND_RADIUS 1.0
#else
#define BND_RADIUS0 0.125
#define BND_RADIUS1 0.75
#define Rin  0.25
#define B_Rin  2.0
#define B_Rout 3.0
#define BND_RADIUS 0.25
#endif

#define Rout 12.0
#define RoutBND 14.0

#define DVAC 1.0e-8
#define _CYLINDER_
  
  inline real get_r2(const vec3 &pos) 
  {
#ifndef _CYLINDER_
    return pos.norm2();
#else
    return sqr(pos.x) + sqr(pos.y);
#endif
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

#if 1
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

#if 1      
      p.vel = p.vel - (p.vel*p.pos)*p.pos/p.pos.norm2();  // subtract radial component
#else
      const vec3 Rc(p.pos.x, p.pos.y, 0.0);
      p.vel = p.vel - (p.vel*Rc)*Rc/Rc.norm2();  // subtract radial component
			if (Rc.abs() < 1.1) p.vel.z = 0.0;
#endif

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

    t_end      = 2500;
    dt_restart = 1.0/64;
    dt_snap    = 1.0;

    dt_restart = std::max(dt_restart, dt_max);
    dt_snap    = std::max(dt_snap,    dt_max);
      
    int Namb_glb  = 0.5e5;
    int Ndisk_glb = 1.5e5;

#if 0
    Namb_glb  = 0.5e6;
    Ndisk_glb = 1.5e6;
#endif

#if 0
    Namb_glb  = 0.5e7;
    Ndisk_glb = 1.5e7;
#endif

//    Namb_glb = 2e5;


    const int  Namb =  Namb_glb/numElements;
    const int Ndisk = Ndisk_glb/numElements;

    ptcl_list.clear();
    ptcl_list.reserve(Namb + Ndisk);

    Rand48 rnd;
    rnd.srand(123 + 123*thisIndex);

    // sample disk
#if 1
    {
      const int N_per_ann = 100;
      const int Nann = Ndisk/N_per_ann;
      const real rout = 0.75*Rout;
      const real eps = 
        std::pow(rout/Rin, 
            1.0/(real)Nann);

      for (int j = 0; j < Nann; j++)
      {
        const real R0 = Rin * std::pow(eps, (real)(j+0));
        const real R1 = Rin * std::pow(eps, (real)(j+1));
        const real dR = R1 - R0;
        const real sig_slope = +1.99f;
        const real mass_ratio = 
          (std::pow(R1,   2.0 - sig_slope) - std::pow(R0,  2.0 - sig_slope))/
          (std::pow(rout, 2.0 - sig_slope) - std::pow(Rin, 2.0 - sig_slope));
        const int Nj = (int)(mass_ratio * Ndisk); 

        for (int i = 0; i < Nj; i++)
        {
          const real ecc = 0.0;
          real sma, phi, R;
          bool flag = true;

          while (flag)
          {

            // SEMI-MAJOR AXIS
            //
            sma  = R0 + dR * rnd.drand();
            phi = 2*M_PI * rnd.drand();

            // CONVERT SEMI-MAJOR AXIS TO RADIUS
            R = sma * ((1.0-sqr(ecc))/(1.0 + (ecc*cos(phi))));

            const real area_scale = R / (sma*(1.0+ecc));
            const real scale    = area_scale;
            const real fv       = rnd.drand();

            if (fv <= scale) flag = false;
          }

          const real H = 4.0 *  HoR * R;
          real z;
          flag = true;
          while (flag)
          {
            z = (2.0*rnd.drand() - 1.0)*R;
            const real fz = rnd.drand();
            const real scale = exp(-sqr(z)/2.0/sqr(H));
            if (fz <= scale) flag = false;
          }

          const vec3 pos(R*cos(phi), R*sin(phi), z);
          ptcl_list.push_back(Particle(ptcl_list.size(), thisIndex, pos));


        }
      }
    }
#endif

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

    local_n = ptcl_list.size();

    generateGeometry_nRelax = 5;
  }

  void System::Problem_generate_IC(const int param)
  {
    if (thisIndex == 0)
    {
      fprintf(stderr, " ********* Setting up mri_disk3d problem ************* \n");
      fprintf(stderr, " GM= %g \n", GM);
      fprintf(stderr, " GM_EPS= %g \n", GM_EPS);
      fprintf(stderr, " HoR= %g \n", HoR);
      fprintf(stderr, " Rin= %g \n", Rin);
      fprintf(stderr, " Rout= %g \n", Rout);
    }
    
    gamma_gas  = 1.0;
    courant_no = 0.8;

    t_global  = 0;
    iteration = 0;

    const real D0       = 1.0;
    const real DENS_MIN = DVAC;


    for (int i = 0; i < local_n; i++) 
    {
      const Particle &pi = ptcl_list[i];

      const vec3 &pos = pi.get_pos();

      const real x = pos.x;
      const real y = pos.y;
      const real z = pos.z;
      
      const real R  = std::sqrt(get_r2(pos)); 
      assert(R > 0.0);
#if 1
      if (pos.abs() < BND_RADIUS)
			{
				mesh_pnts[i].boundary = MeshPoint::OUTFLOW;
			}
			else if (pos.abs() > RoutBND)
			{
        mesh_pnts[i].boundary = MeshPoint::DIOD;
			}
      else  
#endif
			{
        mesh_pnts[i].boundary = MeshPoint::NO_BOUNDARY;
			}

      const real H = HoR * R;

      const real vphi = std::sqrt(GM/R * (1 - sqr(H/R)));
      real vx = -vphi * y/R;
      real vy = +vphi * x/R;
      real vz = 0.0;
      if (R > Rout)
        vx = vy = vz = 0;



      const real dv = vphi * 1.0e-4;
      vx += drand48() * dv;
      vy += drand48() * dv;

      real bx = 0.0;
      real by = 0.0;
      real bz = 0.0;

      real dmid = D0*std::pow(x*x+y*y, -3.0/4.0);
      if (R > Rout)
        dmid = DENS_MIN;

      const real  beta = 25.0;
      const real Hbeta = H*sqrt(1.0 + 1.0/beta);

      real dens = exp(-sqr(z)/2/sqr(Hbeta)) * dmid;
      if (dens < DENS_MIN)
      {
        vx = drand48() * dv;
        vy = drand48() * dv;
        dens =  DENS_MIN;
      }
      
      real pres = sqr(vphi*HoR)*dens;
     
      const real bf = 1.0*std::sqrt(2.0*pres/beta);

      bx = -bf * y/R;
      by = +bf * x/R;
      bz = 0.0;
      if (R < Rin*1.1 || R > 0.9*Rout || std::abs(z) > 16.0*HoR*R*sqrt(1.0+1.0/beta))
        bx = by = bz = 0;

			bx = by = bz = 0.0;
//      if (mesh_pnts[i].boundary != MeshPoint::NO_BOUNDARY)
      if (pos.abs() > RoutBND)
      {
        vx = vy = vz = 0.0;
        bx = by = bz = 0.0;
        dens = DVAC;
        pres = DVAC;
      }

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
#if 1
      const std::pair<vec3,vec3> rv1 = rotate(p.pos, p.vel, dt);
      p.pos = rv1.first;
#else
      const std::pair<vec3,vec3> rv1 = rotate(vec3(p.pos.x, p.pos.y, 0.0), vec3(p.vel.x, p.vel.y,0.0), dt);
      p.pos = vec3(rv1.first.x, rv1.first.y, p.pos.z + p.vel.z * dt);
#endif
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

	real System::Problem_enforce_limiter(const int i)
	{
#if 0
		if (mesh_act[i]->pos.abs() < 1.1*BND_RADIUS) return 0.0;
		else                                         return 1.0;
#else
    if (mesh_act[i]->is_boundary()) return 0.0;
    else                         		return 1.0;
#endif
	}

  real System::Problem_extra_timestep_criterion(const int Id)
  {
    return HUGE;
  }

  inline real get_cs2(const vec3 r)
  {
    return (GM/std::sqrt(get_r2(r)))*sqr(HoR);
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
  
  void System::Problem_set_boundary(const int i)
  {
    assert(mesh_act[i]->is_boundary());
    assert(ptcl_act[i]->is_active());
    const vec3 &pos = mesh_act[i]->pos;

    const real x = pos.x;
    const real y = pos.y;
    const real z = pos.z;

    const real R = std::sqrt(get_r2(pos)); 

    if (R < 2.0*BND_RADIUS)
    {
      const real vphi = std::sqrt(GM/R);
      real vx = -vphi * y/R;
      real vy = +vphi * x/R;
      real vz = 0.0;
      vx = vy = vz = 0;
      const real r = pos.abs();
      const real vr = -1*std::sqrt(get_cs2(pos));
      vx += vr * x/r;
      vy += vr * y/r;
      vz += vr * z/r;
      Wrec_act[i]->w[Fluid::VELX] = vx;
      Wrec_act[i]->w[Fluid::VELY] = vy;
      Wrec_act[i]->w[Fluid::VELZ] = vz;
      Wrec_act[i]->w[Fluid::BX  ] = 0.0;
      Wrec_act[i]->w[Fluid::BY  ] = 0.0;
      Wrec_act[i]->w[Fluid::BZ  ] = 0.0;
      Wrec_act[i]->w[Fluid::DENS] = DVAC;
      Wrec_act[i]->w[Fluid::ETHM] = get_cs2(pos)*DVAC; //get_cs2(pos)(Wrec_act[i]
    }
  }
}
