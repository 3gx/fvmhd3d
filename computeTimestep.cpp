#include "fvmhd3d.h"

namespace fvmhd3d
{

  void System::computeTimestep()
  {
    const int nactive = active_list.size();
    for (int i = 0; i < nactive; i++)
    {
      const MeshPoint &pi = *mesh_act[i];

      const real Volume = pi.Volume;
      assert(Volume > 0.0);

      const Fluid  Wi = U_act[i]->to_primitive(Volume);
      const real    h = std::pow(Volume/(4.0*M_PI/3.0), 1.0/3.0);

      real csi = 0.0, vrel = 0.0;
      if (!pi.is_boundary())
      {
        const real Pi  = Problem_compute_pressure(Wi) * gamma_gas;
        const real B2i = Wi.get_B().norm2();
        csi  = std::sqrt( (Pi + B2i) / Wi[Fluid::DENS] );
        vrel = (Wi.get_vel() - pi.vel).abs();
#if 0
        vrel = std::max(vrel, pi.vel.abs());
#endif
      }

      real dti = h/(csi + vrel + 1.0/HUGE);
      const real fac = pi.is_boundary() ? 0.0 : 1.0;
      const vec3 vi  = Wi.get_vel();
      const vec3 wi  = pi.vel;
      const Cell &ci = cell_list[i];
      const int nface = ci.ngb.size();
      for (int iface = 0; iface < nface; iface++)
      {
        const int j = face_list[ci.ngb[iface]].ngb<true>(i);
        const Fluid_rec &Wrec_j = *Wrec_act[j];

        const vec3 vj  = Wrec_j.w.get_vel();
        const vec3 vij = (vj + vi)*0.5;
        const vec3 dvi = (vj - vi);

        const vec3 wj  = Wrec_j.vel;
        const vec3 wij = (wj + wi)*0.5;
        const vec3 dwi = (wj - wi);

        const real dvmax = std::sqrt(dwi.norm2() + fac*dvi.norm2() + 1.0/HUGE);


        dti = std::min(dti, h/(csi + dvmax));
      }

      assert(dti > 0.0);

      dti = std::min(dti, Problem_extra_timestep_criterion(i));

      mesh_act[i]->dt_new = dti * courant_no;
    }
  }
  
  void System::computeTimestepLimiter()
  {
    const int nactive = active_list.size();
    const int      np = ptcl_act.size();

    std::vector<real> dti_list(np, -HUGE);

    for (int i = 0; i < nactive; i++)
    {
      assert(ptcl_act[i]->is_active());
      dti_list[i] = mesh_act[i]->dt_new;
    }

    for (int i = nactive; i < np; i++)
      if (ptcl_act[i]->is_ngb())
      {
        if (ptcl_act[i]->is_active())
          dti_list[i] = mesh_act[i]->dt_new;
        else
          dti_list[i] = mesh_act[i]->tend - mesh_act[i]->tbeg;
      }

    const real scale_fac = 2.0;

    for (int i = 0; i < nactive; i++)
    {
      const Cell &ci = cell_list[i];
      const int nj = ci.ngb.size();
      for (int k = 0; k < nj; k++)
      {
        const int j = face_list[ci.ngb[k]].ngb<false>(i);
        assert(dti_list[i] > 0.0);
        assert(dti_list[j] > 0.0);
        dti_list[j] = std::min(dti_list[j], scale_fac*dti_list[i]);
        dti_list[i] = std::min(dti_list[i], scale_fac*dti_list[j]);
      }
    }

    for (int i = 0; i < np; i++)
      if (i < nactive_loc || ptcl_act[i]->is_ngb())
      {
        assert(dti_list[i] > 0);
        if (ptcl_act[i]->is_active())
          assert(dti_list[i] <= mesh_act[i]->dt_new);
        mesh_act[i]->dt_new = dti_list[i];
      }
  }

}
