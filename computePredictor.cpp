#include "fvmhd3d.h"


#define DEDNER
//#define DEDNER_ENERGY
//#define REMOVE_MONPOLE


namespace fvmhd3d
{
  void System::computePredictor()
  {
    const int nactive = active_list.size();

#if 0
    /***** first order update *****/
    for (int i = 0; i < nactive; i++)
      Wrec_active[i].t = 0.0;
    return;
#endif
    
    std::vector<Fluid_st> Wst(ptcl_act.size());
    const int np     = ptcl_act.size();
    for (int i = 0; i < np; i++)
    {
      Wst[i].w   = Wrec_act[i]->w; 
      if (i < nactive_loc || ptcl_act[i]->is_ngb())
      {
        assert(Wst[i].w[Fluid::DENS] > 0.0);
        assert(Wst[i].w[Fluid::ETHM] > 0.0);
      }
      Wst[i].bnd = Wrec_act[i]->bnd;
      Wst[i].pos = Wrec_act[i]->pos;
      Wst[i].vel = Wrec_act[i]->vel;
      Wst[i].J   = Wrec_act[i]->J;
      Wst[i].etaJ = Wrec_act[i]->etaJ;
    }


    std::vector<FluidD> dUdt(nactive, 0.0);
    std::vector<real  > dVdt(nactive, 0.0);

    const int nactive_face = face_active_list.size();
    for (int iface = 0; iface < nactive_face; iface++)
    {
      const Face &face = *face_active_list[iface];

      int i  = face.s1;
      int j  = face.s2; 

      //// j <-| i |-> j  ....

      const real area(face.area());
      assert(area > 0.0);
      vec3 normal(face.n * (1.0/area));

      Fluid_st Wst_i = Wst[i];
      Fluid_st Wst_j = Wst[j];

      const vec3 &ipos1 = Wst_i.pos;
      vec3 dri = face.centroid - ipos1;

      assert(std::abs(dri.x) < 0.5*global_domain_size.x);
      assert(std::abs(dri.y) < 0.5*global_domain_size.y);
      assert(std::abs(dri.z) < 0.5*global_domain_size.z);


      if (Wst_i.bnd != MeshPoint::NO_BOUNDARY)
      {
        std::swap(Wst_i, Wst_j);
        std::swap(i, j);
        normal *= -1.0;
        dri    *= -1.0;
      }

      assert (Wst_i.bnd == MeshPoint::NO_BOUNDARY);

      const real dsh = dri * normal;
      const real dsl = 2.0 * dsh;
      assert(dsl != 0.0);
      assert(dsl  > 0.0);

      const vec3 drj = dri - normal * dsl;

      const Fluid &Wi =  Wst_i.w;
      const Fluid &Wj = (Wst_j.bnd == MeshPoint::NO_BOUNDARY || Wst_j.bnd == MeshPoint::INFLOW) ? Wst_j.w : Wst_i.w;

      assert(Wi[Fluid::DENS] > 0.0);
      assert(Wi[Fluid::ETHM] > 0.0);
      assert(Wj[Fluid::DENS] > 0.0);
      assert(Wj[Fluid::ETHM] > 0.0);

      const vec3 &vi = Wst_i.vel;
      const vec3 &vj = Wst_j.vel;

      const vec3  vij = vj + vi;
      const vec3 dvij = vj - vi;

      const vec3  fij = dri - normal * dsh;
#if 0
      const real  fdv = dvij * fij;
      const real ids2 = -0.5*fdv * ((dsh > 0.0) ? 1.0/sqr(dsh) : 0.0);
      const vec3  wij = vij*0.5 +  ids2 * dri;
#else
      const real fdv  = dvij*fij;
      const vec3 rij  = dsl*normal;
      const vec3 dwij = fdv * rij/sqr(dsl);
      const vec3 wij  = (vij*0.5 - dwij);
#endif

      __builtin_prefetch(&dUdt[i]);
      __builtin_prefetch(&dUdt[j]);


      real psi, Bn;
      Fluid flux;
      computeFlux(wij, Wi, Wj, Wst_j.bnd, normal, psi, Bn, flux);

      for (int k = 0; k < Fluid::NFLUID; k++)
        dUdt[i].U[k] -= flux[k] * area;
      dUdt[i].gradPsi += normal * (psi * area);
      dUdt[i].divB    += Bn * area;

      const real dV = area * (rij*wij)/dsl;
      dVdt[i] += dV;

      if (j >= nactive || Wst_j.bnd != MeshPoint::NO_BOUNDARY)
        continue;
      
      dVdt[j] -= dV;

      for (int k = 0; k < Fluid::NFLUID; k++)
        dUdt[j].U[k] += flux[k] * area;
      dUdt[j].gradPsi -= normal * (psi * area);
      dUdt[j].divB    -= Bn * area;
    }

    int ethm_failed = 0;
    int adjust_dt_cnt_dens = 0;
    int adjust_dt_cnt_ethm = 0;
    for (int i = 0; i < nactive; i++)
    {
      if (mesh_act[i]->is_boundary()) 
        continue;

      const Fluid_flt &dotUi = dUdt[i].U;
      const float    divB = dUdt[i].divB;
      const vec3   &dPsi = dUdt[i].gradPsi;

      const Fluid &U0 = *U_act [i];
      const Fluid &W0 =  U0.to_primitive(cell_list[i].Volume);

#if 1
      assert(U0[Fluid::MASS] > 0.0);
      bool flag_cnt = false;
      while(U0[Fluid::MASS] + dotUi[Fluid::MASS]*mesh_act[i]->dt_new <= 0.0)
      {
        if (!flag_cnt)
        { 
          adjust_dt_cnt_dens++;
          flag_cnt = true;
        }
        mesh_act[i]->dt_new *= 0.5;
      }

      flag_cnt = false;
      bool flag_ener = true;
      Fluid W1;
      flag_ener = false;
      while(flag_ener)
      {
        const real dt = mesh_act[i]->dt_new;
        Fluid Uc;
        for (int k = 0; k < Fluid::NFLUID; k++)
          Uc[k] = U0[k] + dotUi[k]*dt;

        assert(Uc[Fluid::MASS] > 0.0);

        if (!Problem_compute_update(Uc, i))
        {
          assert(Uc[Fluid::MASS] > 1.0e-32);

          // subtract divB terms from the induction equation, to make it Gallilean invariant
          Uc[Fluid::WBX] -= W0[Fluid::VELX] * divB*dt;
          Uc[Fluid::WBY] -= W0[Fluid::VELY] * divB*dt;
          Uc[Fluid::WBZ] -= W0[Fluid::VELZ] * divB*dt;

#ifdef REMOVE_MONPOLE
          // subtract force due to monopoles (divB), this breaks conservation!
          const real vB = 
            W0[Fluid::VELX]*W0[Fluid::BX] + 
            W0[Fluid::VELY]*W0[Fluid::BY] + 
            W0[Fluid::VELZ]*W0[Fluid::BZ];

          Uc[Fluid::MOMX] -= W0[Fluid::BX] * divB*dt;
          Uc[Fluid::MOMY] -= W0[Fluid::BY] * divB*dt;
          Uc[Fluid::MOMZ] -= W0[Fluid::BZ] * divB*dt;
          Uc[Fluid::ENER] -= vB     * divB*dt;
#endif

#ifdef DEDNER
          // Hyperbolic divB cleaning (Dedner et al 2002)
          assert(W0[Fluid::ETHM] > 0.0);
          const real pres = Problem_compute_pressure(W0);
          const real B2   = vec3(W0[Fluid::BX], W0[Fluid::BY], W0[Fluid::BZ]).norm2();
          const real dcs2 = gamma_gas * pres + B2;
          const real ch_sig = std::sqrt(dcs2/W0[Fluid::DENS]);

          Uc[Fluid::WBX ] -= dPsi.x*dt;
          Uc[Fluid::WBY ] -= dPsi.y*dt;
          Uc[Fluid::WBZ ] -= dPsi.z*dt;
          Uc[Fluid::MPSI] -= dcs2 * divB*dt;
#ifdef DEDNER_ENERGY
          Uc[Fluid::ENER] -= vec3(W0[Fluid::BX], W0[Fluid::BY], W0[Fluid::BZ]) * dPsi * dt;
#endif
#endif

#ifdef DEDNER
          const double cr = 0.5;
          const real l_min = std::pow(cell_list[i].Volume/(4.0*M_PI/3.0), 1.0/3.0);
          Uc[Fluid::MPSI] *= std::exp(-cr*ch_sig/l_min * dt);
#endif


          assert(Uc[Fluid::MASS] > 0.0);

          // GRAVITY
#if 1
          const Fluid &Ui = U0;
          const vec3 dmom = 
            (Ui[Fluid::MASS]*mesh_act[i]->acc1 + Uc[Fluid::MASS]*mesh_act[i]->acc1)*0.5*dt;
          Uc[Fluid::MOMX] += dmom.x;
          Uc[Fluid::MOMY] += dmom.y;
          Uc[Fluid::MOMZ] += dmom.z;

          const vec3 mom0(Ui[Fluid::MOMX], Ui[Fluid::MOMY], Ui[Fluid::MOMZ]);
          const vec3 mom1(Uc[Fluid::MOMX], Uc[Fluid::MOMY], Uc[Fluid::MOMZ]);
          Uc[Fluid::ENER] += (mom0*mesh_act[i]->acc1 + mom1*mesh_act[i]->acc1)*dt*0.5;
#endif
        }

        W1 = Uc.to_primitive(cell_list[i].Volume*std::exp(dVdt[i]*dt/cell_list[i].Volume));

        assert(W1[Fluid::DENS] > 0.0);

        W1[Fluid::ETHM] = Problem_compute_ethm_update(W1, i);

        if (W1[Fluid::ETHM] <= 0.0)
        {
#if 0
          W1[Fluid::ETHM] = Problem_compute_ethm_from_entropy(W1);
          ethm_failed++;
          W1[Fluid::ETHM] = W0[Fluid::ETHM];
#else
          mesh_act[i]->dt_new *= 0.5;
          if (!flag_cnt)
          {
            adjust_dt_cnt_ethm++;
            flag_cnt = true;
          }
#endif
        }
        else
        {
          W1[Fluid::ENTR] = Problem_compute_entropy_from_ethm(W1);
          flag_ener = false;
        }
      }
#endif
      
      const real dt = scheduler.get_dt(scheduler.get_rung(mesh_act[i]->dt_new))/128;
      Fluid Uc;
      for (int k = 0; k < Fluid::NFLUID; k++)
        Uc[k] = U0[k] + dotUi[k]*dt;

      assert(Uc[Fluid::MASS] > 0.0);

      if (!Problem_compute_update(Uc, i))
      {
        assert(Uc[Fluid::MASS] > 1.0e-32);

        // subtract divB terms from the induction equation, to make it Gallilean invariant
        Uc[Fluid::WBX] -= W0[Fluid::VELX] * divB*dt;
        Uc[Fluid::WBY] -= W0[Fluid::VELY] * divB*dt;
        Uc[Fluid::WBZ] -= W0[Fluid::VELZ] * divB*dt;

#ifdef REMOVE_MONPOLE
        // subtract force due to monopoles (divB), this breaks conservation!
        const real vB = 
          W0[Fluid::VELX]*W0[Fluid::BX] + 
          W0[Fluid::VELY]*W0[Fluid::BY] + 
          W0[Fluid::VELZ]*W0[Fluid::BZ];

        Uc[Fluid::MOMX] -= W0[Fluid::BX] * divB*dt;
        Uc[Fluid::MOMY] -= W0[Fluid::BY] * divB*dt;
        Uc[Fluid::MOMZ] -= W0[Fluid::BZ] * divB*dt;
        Uc[Fluid::ENER] -= vB     * divB*dt;
#endif

#ifdef DEDNER
        // Hyperbolic divB cleaning (Dedner et al 2002)
        assert(W0[Fluid::ETHM] > 0.0);
        const real pres = Problem_compute_pressure(W0);
        const real B2   = vec3(W0[Fluid::BX], W0[Fluid::BY], W0[Fluid::BZ]).norm2();
        const real dcs2 = gamma_gas * pres + B2;
        const real ch_sig = std::sqrt(dcs2/W0[Fluid::DENS]);

        Uc[Fluid::WBX ] -= dPsi.x*dt;
        Uc[Fluid::WBY ] -= dPsi.y*dt;
        Uc[Fluid::WBZ ] -= dPsi.z*dt;
        Uc[Fluid::MPSI] -= dcs2 * divB*dt;
#ifdef DEDNER_ENERGY
        Uc[Fluid::ENER] -= vec3(W0[Fluid::BX], W0[Fluid::BY], W0[Fluid::BZ]) * dPsi * dt;
#endif
#endif

#ifdef DEDNER
        const double cr = 0.5;
        const real l_min = std::pow(cell_list[i].Volume/(4.0*M_PI/3.0), 1.0/3.0);
        Uc[Fluid::MPSI] *= std::exp(-cr*ch_sig/l_min * dt);
#endif


        assert(Uc[Fluid::MASS] > 0.0);

        // GRAVITY
#if 0
        const Fluid &Ui = U0;
        const vec3 dmom = 
          (Ui[Fluid::MASS]*mesh_act[i]->acc1 + Uc[Fluid::MASS]*mesh_act[i]->acc1)*0.5*dt;
        Uc[Fluid::MOMX] += dmom.x;
        Uc[Fluid::MOMY] += dmom.y;
        Uc[Fluid::MOMZ] += dmom.z;

        const vec3 mom0(Ui[Fluid::MOMX], Ui[Fluid::MOMY], Ui[Fluid::MOMZ]);
        const vec3 mom1(Uc[Fluid::MOMX], Uc[Fluid::MOMY], Uc[Fluid::MOMZ]);
        Uc[Fluid::ENER] += (mom0*mesh_act[i]->acc1 + mom1*mesh_act[i]->acc1)*dt*0.5;
#endif
      }

      W1 = Uc.to_primitive(cell_list[i].Volume*std::exp(dVdt[i]*dt/cell_list[i].Volume));

      assert(W1[Fluid::DENS] > 0.0);

      W1[Fluid::ETHM] = Problem_compute_ethm_update(W1, i);

      if (W1[Fluid::ETHM] <= 0.0)
      {
#if 1
        W1[Fluid::ETHM] = Problem_compute_ethm_from_entropy(W1);
        ethm_failed++;
//        W1[Fluid::ETHM] = W0[Fluid::ETHM];
#else
        mesh_act[i]->dt_new *= 0.5;
        if (!flag_cnt)
        {
          adjust_dt_cnt_ethm++;
          flag_cnt = true;
        }
#endif
      }
      else
      {
        W1[Fluid::ENTR] = Problem_compute_entropy_from_ethm(W1);
        flag_ener = false;
      }

      assert(W1[Fluid::ETHM] > 0.0);

      assert(dt > 0.0);
      const real idt = 1.0/dt;
      for (int k = 0; k < Fluid::NFLUID; k++)
        Wrec_act[i]->t[k] = (W1[k] - W0[k]) * idt;
#if 0
      Wrec_act[i]->t[Fluid::VELX] += mesh_act[i]->acc1.x;
      Wrec_act[i]->t[Fluid::VELY] += mesh_act[i]->acc1.y;
      Wrec_act[i]->t[Fluid::VELZ] += mesh_act[i]->acc1.z;
#endif
    }

#if 0
    if (adjust_dt_cnt_dens > 0)
      CkPrintf(" System::computePredictor(): thisIndex= %d, adjust_dt_cnt_dens= %d [ %g ] \n",
          thisIndex, adjust_dt_cnt_dens, nactive > 0 ? 1.0*adjust_dt_cnt_dens/nactive: 0.0);
    if (adjust_dt_cnt_ethm > 0)
      CkPrintf(" System::computePredictor(): thisIndex= %d, adjust_dt_cnt_ethm= %d [ %g ] \n",
          thisIndex, adjust_dt_cnt_ethm, nactive > 0 ? 1.0*adjust_dt_cnt_ethm/nactive: 0.0);
#endif
//    assert(ethm_failed == 0);

  }  /** end System::computePredictor() **/

}
