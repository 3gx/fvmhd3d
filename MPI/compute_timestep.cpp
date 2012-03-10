#include "fvmhd3d.h"


namespace fvmhd3d {

	real system::compute_timestep(const int i0)
	{
    const int i = i0 < 0 ? -1-i0 : i0;
		const Particle &pi = ptcl_import[i];
		const Fluid     Wi = U_import   [i].to_primitive(pi.volume);

    real csi = 0.0, vrel = 0.0;
    if (!pi.is_boundary())
    {
      const real Pi  = compute_pressure(Wi) * gamma_gas;
      const real B2i = Wi.get_B().norm2();
      csi  = std::sqrt( (Pi + B2i) / Wi[Fluid::DENS] );
      vrel = (Wi.get_vel() - (i0 < 0 ? 0.0 : 1.0) * pi.vel).abs();
#if 0
      vrel = Wi.get_vel().abs();
#endif
    }

    const real Volume = pi.volume;
    assert(Volume > 0.0);
    const real h = std::pow(Volume/(4.0*M_PI/3.0), 1.0/3.0);
    real dti = h/(csi + vrel + 1.0/HUGE);
    if (i0 >= 0)
    {
      assert(site_import[i].is_active());
      const real fac = pi.is_boundary() ? 0.0 : 1.0;
      const vec3 vi  = Wi.get_vel();
      const vec3 vip = pi.vel;
      const Cell     &ci = cell_list  [i];
      const int nface = ci.faces().size();
      for (int iface = 0; iface < nface; iface++)
      {
        const int j = face_list[ci.faces()[iface]].ngb<false>(i);

        const vec3 vjp = ptcl_import[j].vel;
        const vec3 wij = (vjp + vip)*0.5;
        const vec3 dwi = (vjp - vip);

        dti = std::min(dti, h/std::max(fac*(csi + (vi - wij).abs()), 1.0/HUGE + dwi.abs()));
      }
    }

    assert(dti > 0.0);
    dti *= courant_no;
    return dti;
  }

  std::vector<real> *dti_list_ptr;
  void system::compute_timesteps(bool shared) 
  {

    if (shared)
    {
      real dt_loc(HUGE);

      for (int i = 0; i < (int)local_n; i++)
      {
        real dti = std::min(extra_timestep(i), compute_timestep(-1-i));
        ptcl_import[i].rung  = scheduler.get_rung(dti);
        dt_loc = std::min(dt_loc, dti);
      }

      double dt_glob;
      MPI_Allreduce(&dt_loc, &dt_glob, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

      const int rung = scheduler.get_rung(dt_glob);

      for (int i = 0; i < (int)local_n; i++)
        ptcl_import[i].rung = rung;
    } 
    else 
    {
      const int nimport = site_import.size();

      const int nactive_site = site_active_list.size();
      std::vector<real> dti_list(nimport, HUGE);

      const real dtfac = -1.0;

      for (int i = 0; i < nimport; i++)
        if (site_import[i].is_active() || site_import[i].is_local_ngb())
          dti_list[i] = (ptcl_import[i].tend - ptcl_import[i].tlast) * (ptcl_import[i].is_boundary() ? dtfac :  +1.0);

      double dt_min_loc = HUGE;;
      for (int isite = 0; isite < nactive_site; isite++)
      {
        const int i = site_active_list[isite];
        dti_list[i] = std::min(extra_timestep(i), compute_timestep(i));
        dt_min_loc = std::min(dt_min_loc, (double)dti_list[i]);
        if (ptcl_import[i].is_boundary()) dti_list[i] = dtfac*dti_list[i];
      }
      dti_list_ptr = &dti_list;

#if 1
      double dt_min_glb;
      MPI_Allreduce(&dt_min_loc, &dt_min_glb, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

      for (std::vector<int>::iterator it = site_active_list.begin(); it < site_active_list.end(); it++)
        if (dti_list[*it] >= 8.0*dt_min_glb) dti_list[*it] *= 0.5;
#endif



      for (int k = 0; k < 1; k++)	
        adjust_timesteps();
    }
  }

  // make sure that ngb steps are at most factor of two different
  void system::adjust_timesteps()
  {
    std::vector<real> &dt_list = *dti_list_ptr;
    const int nimport = site_import.size();
    std::vector<bool> ngb_hash(nimport, false);

    int nj = 0, jfailed1 = 0, jfailed2 = 0, jfailed3 = 0;
    const int nactive_site = site_active_list.size();
    for (int isite = 0; isite < nactive_site; isite++)
    {
      const int i = site_active_list[isite];
      if (ptcl_import[i].is_boundary()) continue;

      const Cell &ci = cell_list[i];

      assert(dt_list[i] > 0.0);

      const real scale_fac = 2.0;
      const int nface = ci.faces().size();
      for (int iface = 0; iface < nface; iface++)
      {
        const int j = face_list[ci.faces()[iface]].ngb<false>(i);
        if (dt_list[j] < 0.0) continue;
        ngb_hash[j] = true;
        dt_list[i] = std::min(dt_list[i], scale_fac*dt_list[j]);
        dt_list[j] = std::min(dt_list[j], scale_fac*dt_list[i]);
      }
    }

    for (int i = 0; i < nimport; i++)
    {
      if (!(site_import[i].is_active() || site_import[i].is_local_ngb())) continue;

      if (ptcl_import[i].is_active())
      {
        ptcl_import[i].rung  = scheduler.get_rung(std::abs(dt_list[i]));
        assert(ptcl_import[i].tend == t_global);
      }
      else if (ngb_hash[i])
      {
        assert(dt_list[i] > 0.0);
#if 1
        assert(!ptcl_import[i].is_boundary());
#endif
        assert(ptcl_import[i].tend > t_global);
        nj++;

        if (Wrec_import[i].bnd <= -128 && !ptcl_import[i].is_boundary())
        {
#if 1
          assert(false);
#endif
          jfailed1++;
          const int irung = scheduler.get_rung(dt_list[-128-Wrec_import[i].bnd]);
          const int jrung = scheduler.get_rung(compute_timestep(-1-i));
          ptcl_import[i].rung = std::max(ptcl_import[i].rung, std::max(jrung, irung));
          ptcl_import[i].tend = t_global;
        }
        else
        {
          const real tlast = ptcl_import[i].tlast;
          const real tend  = ptcl_import[i].tend;
          if ((tend - tlast) > dt_list[i])
          {
            assert(ptcl_import[i].tlast < t_global);
            assert(ptcl_import[i].tend  > t_global);
            const int rung = std::max(Scheduler::exp_of(scheduler.dtmax) - Scheduler::exp_of(std::abs(0.5*dt_list[i])), 0);
            const real dt  = scheduler.get_dt(rung);
            if (ptcl_import[i].tlast + dt <= t_global)
            {
              jfailed2++;
              ptcl_import[i].rung = scheduler.get_rung(dt);
              ptcl_import[i].tend = t_global + scheduler.get_dt(ptcl_import[i].rung);
              ptcl_import[i].rung = -1-ptcl_import[i].rung;
            }
            else
            {
              jfailed3++;
              assert(rung > ptcl_import[i].rung);
              ptcl_import[i].rung = rung;
              ptcl_import[i].tend = ptcl_import[i].tlast + dt;
            }
          }
        }
      }
    }
#if 0
    const int jfailed = jfailed1+jfailed2+jfailed3;
    if (jfailed != 0)
      fprintf(stderr, " myproc= %d :: nj= %d  jfailed= %d [ %g ] ( %d %d %d ) \n", myproc,
          nj, jfailed,  nj> 0 ? 1.0*jfailed/nj : 0.0, jfailed1, jfailed2, jfailed3);
#endif
  }

};
