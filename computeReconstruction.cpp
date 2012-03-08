#include "fvmhd3d.h"

namespace fvmhd3d
{
  inline real Phi(const real x, const real eps1)
  {
#if 1
    return 1.0;
#endif
    assert(x >= 0.0);
#if 1
    return std::min(1.0, x);
#else
		const real eps = 1.0e-16;
    return (x*x + 2.0*x + eps)/(x*x + x + 2.0 + eps);
#endif
  }
  
  inline real limiter(const real dWi, const real dWj_min, const real dWj_max, const real h3, const real Wi)
  {
    if (std::abs(dWi) <= 1.0e-16*std::max(std::abs(dWj_max), std::abs(dWj_min)))  return 1.0;

    const real invWi = 1.0/dWi;
#if 0
    const real K   = 0.1;
    const real eps = K*K*K*h3* sqr(invWi);
    if (dWi < 0.0) return Phi(dWj_min * invWi, eps);
    else           return Phi(dWj_max * invWi, eps);
#else
    const real K   = 0.01;
    const real dWj = std::min(-dWj_min, dWj_max);
    assert(dWj >= 0.0);
    const real dWi_max = K*std::max(std::abs(Wi), dWj);
    const real eps = sqr(dWi_max*std::abs(invWi));
    
    if (dWi < 0.0) return Phi(dWj_min * invWi, eps);
    else           return Phi(dWj_max * invWi, eps);
#endif

  }

  void System::computeReconstruction()
  {
    const int nactive = active_list.size();

    std::vector<Fluid_st> Wst_active(ptcl_act.size());
    const int np     = ptcl_act.size();
    for (int i = 0; i < np; i++)
    {
      Wst_active[i].w   = Wrec_act[i]->w; 
      Wst_active[i].bnd = Wrec_act[i]->bnd;
      Wst_active[i].pos = Wrec_act[i]->pos;
      Wst_active[i].vel = Wrec_act[i]->vel;
      Wst_active[i].tend = Wrec_act[i]->tend;
    }

    std::vector<vec3> centroid_list;
    std::vector<real> dtij_list;
    for (int i = 0; i < nactive; i++)
    {
      const Cell     &ci   =  cell_list[i];
      const Fluid_st &Wst  = Wst_active[i];
      const Fluid    &Wi   = Wst.w;
      const vec3     &ipos = Wst.pos;

      if (mesh_act[i]->is_boundary())
      {
        Wrec_act[i]->w = Wi;
        Wrec_act[i]->x = 0.0;
        Wrec_act[i]->y = 0.0;
        Wrec_act[i]->z = 0.0;
        Wrec_act[i]->t = 0.0;
        continue;
      }


      Fluid Wx(0.0);
      Fluid Wy(0.0);
      Fluid Wz(0.0);
      Fluid_flt dWj_min(0.0), dWj_max(0.0);

      centroid_list.clear();
      dtij_list.clear();
      const int nface = ci.ngb.size();
      bool zero_flag = false;

      vec3 sum_centroid(0.0);
      const vec3 ivel = Wi.get_vel();

      const real dtI = Wst.tend - t_global;
      assert(dtI > 0.0);

      for (int iface = 0; iface < nface; iface++)
      {
        const Face &face = face_list[ci.ngb[iface]];

        vec3 dri = face.centroid - ipos;
        if      (dri.x >  0.5*global_domain_size.x) dri.x -= global_domain_size.x;
        else if (dri.x < -0.5*global_domain_size.x) dri.x += global_domain_size.x;
        if      (dri.y >  0.5*global_domain_size.y) dri.y -= global_domain_size.y;
        else if (dri.y < -0.5*global_domain_size.y) dri.y += global_domain_size.y;
        if      (dri.z >  0.5*global_domain_size.z) dri.z -= global_domain_size.z;
        else if (dri.z < -0.5*global_domain_size.z) dri.z += global_domain_size.z;
        assert(std::abs(dri.x) < 0.5*global_domain_size.x);
        assert(std::abs(dri.y) < 0.5*global_domain_size.y);
        assert(std::abs(dri.z) < 0.5*global_domain_size.z);

        const vec3 centroid = dri;
        centroid_list.push_back(centroid);
        sum_centroid += centroid;


        const real area     = face.area();
        assert(area > 0.0);
        const vec3 normal   = face.n * ((centroid * face.n < 0.0) ? (-1.0/area) : (1.0/area));
        const real dsh  = centroid * normal;
        assert(dsh > 0.0);

        const real ids  = (dsh != 0.0) ? 0.5/dsh : 0.0;
        const real idsA = area * ids;

        const vec3 drh = normal   * dsh;
        const vec3 fij = centroid - drh;

        const int j = face.ngb<false>(i);
        assert(j >= 0);
        assert(j < (int)ptcl_act.size());

        if (Wst_active[j].bnd != MeshPoint::NO_BOUNDARY && Wst_active[j].bnd != MeshPoint::INFLOW)
        {
          zero_flag = true;
          break;
        }

        const real dtJ = Wst_active[j].tend - t_global;
//        assert(dtJ > 0.0);
        const real dth = 0.5*std::min(dtI, dtJ);
        dtij_list.push_back(dtI - dth);
        const Fluid &Wj = Wst_active[j].w;
        for (int k = 0; k < Fluid::NFLUID; k++) 
        {
          const real  sum = Wj[k] + Wi[k];
          const real diff = Wj[k] - Wi[k];
          dWj_min[k] = std::min((real)dWj_min[k], diff);
          dWj_max[k] = std::max((real)dWj_max[k], diff);
          Wx[k] += (sum*drh.x + diff*fij.x)*idsA;
          Wy[k] += (sum*drh.y + diff*fij.y)*idsA;
          Wz[k] += (sum*drh.z + diff*fij.z)*idsA;
        }
      }

      if (zero_flag)
      {
        Wrec_act[i]->w = Wi;
        Wrec_act[i]->x = 0.0;
        Wrec_act[i]->y = 0.0;
        Wrec_act[i]->z = 0.0;
        Wrec_act[i]->t = 0.0;
        continue;
      }

      assert(cell_list[i].Volume > 0.0);
      const real invV = 1.0/cell_list[i].Volume;
      for (int k = 0; k < Fluid::NFLUID; k++) 
      {
        Wx[k] *= invV;
        Wy[k] *= invV;
        Wz[k] *= invV;
      }

      const vec3 dxB = Wx.get_B();
      const vec3 dyB = Wy.get_B();
      const vec3 dzB = Wz.get_B();

      Wrec_act[i]->J.x = dyB.z - dzB.y;
      Wrec_act[i]->J.y = dzB.x - dxB.z;
      Wrec_act[i]->J.z = dxB.y - dyB.x;


      //
      // Limit gradients in primitive variables
      //

      Fluid limiter_min(1.0);
      real sum_dt = 0.0;
      for (int iface = 0; iface < nface; iface++)
        for (int k = 0; k < Fluid::NFLUID; k++) 
        {
          const real dWx = vec3(Wx[k], Wy[k], Wz[k])*centroid_list[iface];
          const real dWt = 0.0*Wrec_act[i]->t[k]*dtij_list[iface];
          assert(dtij_list[iface] > 0.0);
          const real dWi = dWx + dWt;
          sum_dt += dtij_list[iface];
          limiter_min[k] = std::min(limiter_min[k], limiter(dWi, dWj_min[k], dWj_max[k], mesh_act[i]->Volume, Wi[k]));
#if 0
          if (k == Fluid::ETHM)
          {
            assert(Wi[k] + limiter_min[k]*dWi >= 0.0);
//            assert(Wi[k] + dWi >= 0.0);
          }
          if (k == Fluid::DENS)
          {
            assert(Wi[k] + limiter_min[k]*dWi >= 0.0);
//            assert(Wi[k] + dWi >= 0.0);
          }
#endif
        }

      *Wrec_minmax_act[i] = std::make_pair(dWj_min, dWj_max);



			const real f0 = Problem_enforce_limiter(i);
      for (int k = 0; k < Fluid::NFLUID; k++) 
      {
        Wrec_minmax_act[i]->first [k] += Wi[k];
        Wrec_minmax_act[i]->second[k] += Wi[k];
        real tau = f0*limiter_min[k];

#if 0
        if (k == Fluid::BX) tau = 1;
        if (k == Fluid::BY) tau = 1;
        if (k == Fluid::BZ) tau = 1;
#endif

        assert(tau == 1.0);
        Wx[k] *= tau;
        Wy[k] *= tau;
        Wz[k] *= tau;

#if 1
        Wrec_act[i]->t[k] *= tau;
#endif

        Wrec_act[i]->w[k] = Wi[k];
        Wrec_act[i]->x[k] = Wx[k];
        Wrec_act[i]->y[k] = Wy[k];
        Wrec_act[i]->z[k] = Wz[k];
      }

#if 0
      // Positivity for the *-state, Berthon, JCoP 2006
      {
        const int k = Fluid::DENS;
        const real dWx = vec3(Wx[k], Wy[k], Wz[k])*sum_centroid;
        const real dWp = dWx + Wrec_act[i]->t[k]*sum_dt;
        const real dW  = std::max(dWx, dWp);

        if (dW > 0.0)
        {
          const real alpha = std::min(1.0, 0.9*Wi[k]/dW);
          assert(alpha >= 0.0);
          assert(alpha <= 1.0);
          Wrec_act[i]->x[k] *= alpha;
          Wrec_act[i]->y[k] *= alpha;
          Wrec_act[i]->z[k] *= alpha;
          Wrec_act[i]->t[k] *= alpha;
        }
      }
#endif
    }


  }  /** end System::computeReconstruction() **/


}
