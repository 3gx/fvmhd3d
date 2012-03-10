#include "fvmhd3d.h"

/** ENABLE FOR CLOCKCYCLE TIMINGS ***/
// #define TICKPROFILE

#ifdef TICKPROFILE
#define RDTSCL0(x) rdtscl(x)
#else
#define RDTSCL0(x)
#endif

#define HLLD
#define DEDNER
//#define DEDNER_ENERGY
//#define REMOVE_MONPOLE

namespace fvmhd3d 
{

  enum failed_type {FAILED_MASS = 1, FAILED_ETHM = 2};

	void system::compute_update()
	{
#if 0
    {
      const int nactive_site = site_active_list.size();
      for (int isite = 0; isite < nactive_site; isite++)
      {
        const int i = site_active_list[isite];


        const Cell     &ci     = cell_list[i];
        const vec3     &ipos   = ptcl_import[i].pos;
        vec3 sum_area = 0.0;
        const int nface = ci.faces().size();	
        for (int iface = 0; iface < nface; iface++)
        {
          const Face &face = face_list[ci.faces()[iface]];
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
          const real area     = face.area();
          assert(area > 0.0);
          const vec3 normal   = face.n * ((centroid * face.n < 0.0) ? (-1.0/area) : (1.0/area));

          const real dsh  = centroid * normal;
          assert(dsh > 0.0);

          sum_area += normal * area;
        }
        if (!(sum_area.abs() < 1.0e-10*std::pow(cell_list[i].Volume, 2.0/3.0)))
        {
          fprintf(stderr, " nactive_ste= %d   with_ngb= %d\n",
              nactive_site, nactive_site);
          fprintf(stderr, " isite= %d  active= %d, area= %g %g %g   %g   vol= %g pos= %g   a= %g \n",
              isite,ptcl_import[i].is_active(),
              sum_area.x, sum_area.y, sum_area.z,
              sum_area.abs(),
              cell_list[i].Volume, ptcl_import[i].pos.abs(),
              std::pow(cell_list[i].Volume, 2.0/3.0));

          for (int iface = 0; iface < nface; iface++)
          {
            Face &face = face_list[ci.faces()[iface]];
            int s1 = face.s1 < 0 ? -1-face.s1 : face.s1;
            fprintf(stderr, "iface= %d i= %d nimport= %d :  %d:%d  %d:%d bnd= %d:%d [%d;%d] area= %g \n",
                iface, i, (int)site_import.size(),
                s1, face.s2,
                ptcl_import[s1].is_active(),
                ptcl_import[face.s2].is_active(),
                ptcl_import[s1].boundary,
                ptcl_import[face.s2].boundary,
                Wrec_import[s1].bnd,
                Wrec_import[face.s2].bnd,
                face.area());
          }
          fprintf(stderr, "volume= %g \n", cell_list[i].Volume);
        }
        assert(sum_area.abs() < 1.0e-10*std::pow(cell_list[i].Volume, 2.0/3.0));
      }

    }
#endif

    std::vector< std::pair<int, failed_type> > failed_ptcl;
    std::vector<vec3> area_test(site_import.size(), 0.0);

		const int nimport = site_import.size(); 
		const std::vector<FluidD> dU0(dU_import);
    for (int i = 0; i < nimport; i++)
      if (site_import[i].is_active() || site_import[i].is_local_ngb())
      {
        dU_import[i] = 0.0;
//        compute_update_prob(Wrec_import[i].w, -1-i);
      }

    const int nactive_face = nface_active;
    for (int iface = 0; iface < nactive_face; iface++)
    {
      const Face &face = *face_active_list[iface];

      int i  = face.s1;
      int j  = face.s2;

#if 1
      assert(i != j);
      assert(i >= 0);
      assert(i < nimport);
      assert(j >= 0);
      assert(j < nimport);
#endif

#if 1
      assert(ptcl_import[i].boundary != -1289);
      assert(ptcl_import[j].boundary != -1289);
#endif

      real area  (face.area());
      vec3 normal(face.n * (1.0/area));

      //// j <-| i |-> j  ....

      Fluid_rec Wrec_i = Wrec_import[i];
      Fluid_rec Wrec_j = Wrec_import[j];

      const vec3 ipos = Wrec_i.pos;
      vec3 dri = face.centroid - ipos;

      assert(std::abs(dri.x) < 0.5*global_domain_size.x);
      assert(std::abs(dri.y) < 0.5*global_domain_size.y);
      assert(std::abs(dri.z) < 0.5*global_domain_size.z);

      if (Wrec_i.bnd != Particle::NO_BOUNDARY)
      {
        std::swap(Wrec_i, Wrec_j);
        std::swap(i, j);
        normal *= -1.0;
        dri    *= -1.0;
      }


#if 1
      assert(Wrec_i.bnd == Particle::NO_BOUNDARY);
#else
      if (Wrec_i.bnd != Particle::NO_BOUNDARY) continue;
#endif

      const real dtI = t_global - Wrec_i.tlast;
      const real dtJ = t_global - Wrec_j.tlast;
      const real dt  = std::min(dtI, dtJ);
      assert(dt > 0.0);
      const real dti = dtI - 0.5*dt;
      const real dtj = dtJ - 0.5*dt;

      const real dsh = dri * normal;
      const real dsl = 2.0 * dsh;
#if 1
      assert(dsl != 0.0);
      assert(dsl  > 0.0);
#else
      assert(dsl > 0.0);
#endif

      const vec3 drj = dri - normal * dsl;

      Fluid Wi(Wrec_i.w), Wj(Wrec_j.w);


      if (Wrec_j.bnd != Particle::NO_BOUNDARY)
      {
        for (int k = 0; k < Fluid::NFLUID; k++) 
          Wi[k] = Wj[k] = Wrec_i.w[k];
      }
      else
      {
        assert(Wi[Fluid::DENS] > 0.0);
        assert(Wi[Fluid::ETHM] > 0.0);
        assert(Wj[Fluid::DENS] > 0.0);
        assert(Wj[Fluid::ETHM] > 0.0);

#if 0
        for (int k = 0; k < Fluid::NFLUID; k++) 
        {
          Wi[k] += vec3(Wrec_i.x[k], Wrec_i.y[k], Wrec_i.z[k]) * dri;
          Wj[k] += vec3(Wrec_j.x[k], Wrec_j.y[k], Wrec_j.z[k]) * drj;
        }

        assert(Wi[Fluid::DENS] > 0.0);
        assert(Wi[Fluid::ETHM] > 0.0);
        assert(Wj[Fluid::DENS] > 0.0);
        assert(Wj[Fluid::ETHM] > 0.0);

        for (int k = 0; k < Fluid::NFLUID; k++) 
        {
          Wi[k] += Wrec_i.t[k]*dti;
          Wj[k] += Wrec_j.t[k]*dtj;
        }
#else
        for (int k = 0; k < Fluid::NFLUID; k++) 
        {
          Wi[k] += vec3(Wrec_i.x[k], Wrec_i.y[k], Wrec_i.z[k]) * dri + Wrec_i.t[k]*dti;
          Wj[k] += vec3(Wrec_j.x[k], Wrec_j.y[k], Wrec_j.z[k]) * drj + Wrec_j.t[k]*dtj;
        }
#endif
      }

      assert(Wi[Fluid::DENS] > 0.0);
      assert(Wi[Fluid::ETHM] > 0.0);
      assert(Wj[Fluid::DENS] > 0.0);
      assert(Wj[Fluid::ETHM] > 0.0);

      const vec3 &vi = Wrec_i.vel;
      const vec3 &vj = Wrec_j.vel;

      const vec3  vij = vj + vi;
      const vec3 dvij = vj - vi;

      const vec3  fij = dri - normal * dsh;
      const real  fdv = dvij * fij;
      const real ids2 = -0.5*fdv * ((dsh > 0.0) ? 1.0/sqr(dsh) : 0.0);
      const vec3  wij = vij*0.5 +  ids2 * dri;

      real psi, Bn;
      Fluid flux;
      compute_flux(wij, Wi, Wj, Wrec_j.bnd, normal, psi, Bn, flux);

#if 0
      // resistivity
      const vec3 ResFlux = normal * ((Wrec_i.J + Wrec_i.Jdot*dti) + (Wrec_j.J + Wrec_j.Jdot*dtj))*0.5;
      flux[Fluid::BX] += ResFlux.x; 
      flux[Fluid::BY] += ResFlux.y; 
      flux[Fluid::BZ] += ResFlux.z; 
#endif

      // GRAVITY flux
#if 0
      const vec3 aij     = (Wrec_i.acc + Wrec_j.acc)*0.5;
      flux[Fluid::ENER] += flux[Fluid::MASS]*(aij*dri);
#endif

      const real tarea = dt * area;
      const vec3 dPsi  = normal * (psi * tarea);


      for (int k = 0; k < Fluid::NFLUID; k++)
        dU_import[i].U[k]  -= flux[k] * tarea;
      dU_import[i].gradPsi += dPsi;
      dU_import[i].divB    += Bn * tarea;   

      area_test[i] += normal*area;

      if (Wrec_j.bnd != Particle::NO_BOUNDARY) continue;     // if jsite is a boundary, skip its update

      for (int k = 0; k < Fluid::NFLUID; k++)
        dU_import[j].U[k]  += flux[k] * tarea;
      dU_import[j].gradPsi -= dPsi;
      dU_import[j].divB    -= Bn * tarea;   

      area_test[j] -= normal*area;
    }


    int ethm_failed = 0;
    //    int failed = 0;
    const int nactive_site = site_active_list.size();
    for (int isite = 0; isite < nactive_site; isite++)
    {
      const int i = site_active_list[isite];
#if 0
      const real v0 = ptcl_import[i].volume;
      const real v1 = cell_list[i].Volume;
      if (std::abs(v0 - v1) > 1.0e-10*std::abs(v0))
      {
        fprintf(stderr, "isite= %d  i= %d id= %d  v0= %g  v1= %g  df= %g [ %g ] :: pos= %g %g %g   vel= %g %g %g \n",
            isite, i, ptcl_import[i].local_id,  v0, v1, v1 - v0, (v1 - v0)/v0, 
            ptcl_import[i].pos.x,
            ptcl_import[i].pos.y,
            ptcl_import[i].pos.z,
            ptcl_import[i].vel.x,
            ptcl_import[i].vel.y,
            ptcl_import[i].vel.z);
        failed++;
      }
      if (failed > 10) assert(false);
#endif


#if 0 
      ptcl_import[i].volume = cell_list[i].Volume;
      continue;
#endif
      if (ptcl_import[i].is_boundary())
      { 
        ptcl_import[i].volume = cell_list[i].Volume;
        continue;
      }

      {
        const vec3 areai = area_test[i];
        if (!(areai.abs() < 1.0e-10*std::pow(cell_list[i].Volume, 2.0/3.0)))
        {
          fprintf(stderr, "myproc= %d::  [%lld; %d] isite= %d  active= %d, area= %g %g %g   %g   vol= %g pos= %g   a= %g \n",
              myproc, ptcl_import[i].idx, ptcl_import[i].boundary,
              isite,ptcl_import[i].is_active(),
              areai.x, areai.y, areai.z,
              areai.abs(),
              cell_list[i].Volume, ptcl_import[i].pos.abs(),
              std::pow(cell_list[i].Volume, 2.0/3.0));
        }
        assert(areai.abs() < 1.0e-10*std::pow(cell_list[i].Volume, 2.0/3.0));
      }

      Fluid Uc(U_import[i]);
      for (int k = 0; k < Fluid::NFLUID; k++)
        Uc[k] += dU0[i].U[k] + dU_import[i].U[k];

      const real divB = dU0[i].divB    + dU_import[i].divB;
      const vec3 dPsi = dU0[i].gradPsi + dU_import[i].gradPsi;

      const real dt = t_global - Wrec_import[i].tlast;
      dU_import[i].U    = 0.0;
      dU_import[i].divB = dt > 0.0 ? divB/dt/cell_list[i].Volume : divB;
      dU_import[i].gradPsi = 0.0;

      if (Uc[Fluid::MASS] <= 0.0)
      {
        fprintf(stderr, "update: proc= %d idx= %lld i= %d : dt= %g [ %g ] pos= %g  vel= %g \n",
            myproc, ptcl_import[i].idx, i, 
            ptcl_import[i].tend - ptcl_import[i].tlast,
            scheduler.get_dt(ptcl_import[i].rung),
            ptcl_import[i].pos.abs(), ptcl_import[i].vel.abs());
        failed_ptcl.push_back(std::make_pair(i, FAILED_MASS));
        continue;
      }
      assert(Uc[Fluid::MASS] > 0.0);    // check for NaN

      if (!compute_update_prob(Uc, i))
      {
        // --- interpolate in primitives ---
        //
        const real dth = 0.5 * dt;
        Fluid Wh(Wrec_import[i].w);
        for (int k = 0; k < Fluid::NFLUID; k++)
          Wh[k] += Wrec_import[i].t[k] * dth;

        // subtract divB terms from the induction equation, to make it Gallilean invariant
        Uc[Fluid::WBX] -= Wh[Fluid::VELX] * divB;
        Uc[Fluid::WBY] -= Wh[Fluid::VELY] * divB;
        Uc[Fluid::WBZ] -= Wh[Fluid::VELZ] * divB;

#ifdef REMOVE_MONPOLE
        // subtract force due to monopoles (divB), this breaks conservation!
        const real vB = 
          Wh[Fluid::VELX]*Wh[Fluid::BX] + 
          Wh[Fluid::VELY]*Wh[Fluid::BY] + 
          Wh[Fluid::VELZ]*Wh[Fluid::BZ];

        Uc[Fluid::MOMX] -= Wh[Fluid::BX] * divB;
        Uc[Fluid::MOMY] -= Wh[Fluid::BY] * divB;
        Uc[Fluid::MOMZ] -= Wh[Fluid::BZ] * divB;
        Uc[Fluid::ENER] -= vB     * divB;
#endif

#ifdef DEDNER
        /* Hyperbolic divB cleaning (Dedner et al 2002) */
        assert(Wh[Fluid::ETHM] > 0.0);
        const real pres = compute_pressure(Wh);
        const real B2   = vec3(Wh[Fluid::BX], Wh[Fluid::BY], Wh[Fluid::BZ]).norm2();
        const real dcs2 = gamma_gas * pres + B2;
        const real ch_sig = std::sqrt(dcs2/Wh[Fluid::DENS]);
        Uc[Fluid::WBX ] -= dPsi.x;
        Uc[Fluid::WBY ] -= dPsi.y;
        Uc[Fluid::WBZ ] -= dPsi.z;
        Uc[Fluid::MPSI] -= dcs2 * divB;
#ifdef DEDNER_ENERGY
        Uc[Fluid::ENER] -= vec3(Wh[Fluid::BX], Wh[Fluid::BY], Wh[Fluid::BZ]) * dPsi;
#endif
#endif

        // GRAVITY 2nd
#if 1
        const Fluid Ui = U_import[i];
        const vec3 dmom = 
          (Ui[Fluid::MASS]*ptcl_import[i].acc0 + Uc[Fluid::MASS]*ptcl_import[i].acc1)*0.5*dt;
        Uc[Fluid::MOMX] += dmom.x;
        Uc[Fluid::MOMY] += dmom.y;
        Uc[Fluid::MOMZ] += dmom.z;

        const vec3 mom0(Ui[Fluid::MOMX], Ui[Fluid::MOMY], Ui[Fluid::MOMZ]);
        const vec3 mom1(Uc[Fluid::MOMX], Uc[Fluid::MOMY], Uc[Fluid::MOMZ]);
        Uc[Fluid::ENER] += (mom0*ptcl_import[i].acc0 + mom1*ptcl_import[i].acc1)*dt*0.5;
#endif

        const double cr = 0.5;
        const real l_min = std::pow(cell_list[i].Volume/(4.0*M_PI/3.0), 1.0/3.0);
        Uc[Fluid::MPSI] *= std::exp(-cr * ch_sig/l_min * dt);
      }
      ptcl_import[i].volume = cell_list[i].Volume;

      assert(Uc[Fluid::MASS] > 0.0);
      Fluid W1 = Uc.to_primitive(cell_list[i].Volume);

      assert(W1[Fluid::DENS] > 0.0);
      W1[Fluid::ETHM] = compute_ethm_update(W1, i);

      if (W1[Fluid::ETHM] <= 0.0)
      {
        W1[Fluid::ETHM] = compute_ethm_from_entropy(W1);
        ethm_failed++;
        assert(false);
#if 0
        failed_ptcl.push_back(std::make_pair(i, FAILED_ETHM));
        continue;
#endif
      }
      else
        W1[Fluid::ENTR] = compute_entropy_from_ethm(W1);

      assert(W1[Fluid::ETHM] > 0.0);

      U_import [i]      = W1.to_conservative(cell_list[i].Volume);

      const Fluid Wi = U_import[i].to_primitive(ptcl_import[i].volume);
      assert(Wi[Fluid::ETHM] > 0);
    }


    assert(ethm_failed == 0);
    for (int ifailed = 0; ifailed < (const int)failed_ptcl.size(); ifailed++)
    {
      const int i = failed_ptcl[ifailed].first;
      const real dt = t_global - Wrec_import[i].tlast;
      Fluid W1;
      for (int k = 0; k < Fluid::NFLUID; k++)
        W1[k] = Wrec_import[i].w[k] + Wrec_import[i].t[k]*dt;
      W1 = (W1.to_conservative(ptcl_import[i].volume)).to_primitive(cell_list[i].Volume);

      assert(W1[Fluid::MASS] > 0.0);

      ptcl_import[i].volume = cell_list[i].Volume;

      assert(W1[Fluid::DENS] > 0.0);
      W1[Fluid::ETHM] = compute_ethm_update(W1, i);

      if (W1[Fluid::ETHM] <= 0.0)
      {
        W1[Fluid::ETHM] = compute_ethm_from_entropy(W1);
        ethm_failed++;
        assert(false);
      }
      else
        W1[Fluid::ENTR] = compute_entropy_from_ethm(W1);

      assert(W1[Fluid::ETHM] > 0.0);

      U_import [i]   = W1.to_conservative(cell_list[i].Volume);
      const Fluid Wi = U_import[i].to_primitive(ptcl_import[i].volume);
      assert(Wi[Fluid::ETHM] > 0);

#if 0
      ////////////// FAILED_MASS
      const Cell     &ci     = cell_list[i];
      const int nface = ci.faces().size();	
      for (int iface = 0; iface < nface; iface++)
      {
        const Face &face = face_list[ci.faces()[iface]];
        const int j = face.ngb<false>(i);
        Wrec_import[j].bnd = -128 - i;
        dU_import[j].U = 0.0;
        dU_import[j].divB = 0.0;
        dU_import[j].gradPsi = 0.0;

        const real dt = t_global - Wrec_import[j].tlast;
        Fluid W1;
        for (int k = 0; k < Fluid::NFLUID; k++)
          W1[k] = Wrec_import[j].w[k] + Wrec_import[j].t[k]*dt;
        assert(W1[Fluid::MASS] > 0.0);
        assert(W1[Fluid::ETHM] > 0.0);

        U_import [j]   = W1.to_conservative(ptcl_import[j].volume);

        const Fluid Wi = U_import[j].to_primitive(ptcl_import[j].volume);
        assert(Wi[Fluid::ETHM] > 0);
      }
#endif



    }
  }
  /////////////////////////////////
  /////////////////////////////////
  /////////////////////////////////

  void system::compute_update_first_order(const std::vector<int> &ilist)
  {
    for (int isite = 0; isite < (const int)ilist.size(); isite++)
    {
      const int i = ilist[isite];

      FluidD dUdt(0.0); 

      const Cell     &ci   = cell_list[i];
      const int nface = ci.faces().size();	
      for (int iface = 0; iface < nface; iface++)
      {
        const Face &face = face_list[ci.faces()[iface]];

        int i  = face.s1;
        int j  = face.s2; 

        real area  (face.area());
        vec3 normal(face.n * (1.0/area));

        Fluid_rec Wrec_i = Wrec_import[i];
        Fluid_rec Wrec_j = Wrec_import[j];

        const vec3 ipos = Wrec_i.pos;
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


        if (Wrec_i.bnd != Particle::NO_BOUNDARY)
        {
          std::swap(Wrec_i, Wrec_j);
          std::swap(i, j);
          normal *= -1.0;
          dri    *= -1.0;
        }
        assert(Wrec_i.bnd == Particle::NO_BOUNDARY);

        const real dsh = dri * normal;
        const real dsl = 2.0 * dsh;
#if 1
				assert(dsl != 0.0);
				assert(dsl  > 0.0);
#else
				assert(dsl > 0.0);
#endif

				const Fluid &Wi =  Wrec_i.w;
				const Fluid &Wj = (Wrec_j.bnd == Particle::NO_BOUNDARY) ? Wrec_j.w : Wrec_i.w;

				assert(Wi[Fluid::DENS] > 0.0);
				assert(Wi[Fluid::ETHM] > 0.0);
				assert(Wj[Fluid::DENS] > 0.0);
				assert(Wj[Fluid::ETHM] > 0.0);

				const vec3 &vi = Wrec_i.vel;
				const vec3 &vj = Wrec_j.vel;

				const vec3  vij = vj + vi;
				const vec3 dvij = vj - vi;

				const vec3  fij = dri - normal * dsh;
				const real  fdv = dvij * fij;
				const real ids2 = -0.5*fdv * ((dsh > 0.0) ? 1.0/sqr(dsh) : 0.0);
				const vec3  wij = vij*0.5 +  ids2 * dri;

				real psi, Bn;
				Fluid flux;
				compute_flux(wij, Wi, Wj, Wrec_j.bnd, normal, psi, Bn, flux);

				for (int k = 0; k < Fluid::NFLUID; k++)
					dUdt.U[k] -= flux[k] * area;
				dUdt.gradPsi += normal * (psi * area);
				dUdt.divB    += Bn * area;
			}

			const Fluid dotUi = dUdt.U;
			const real divB   = dUdt.divB;
			const vec3 dPsi   = dUdt.gradPsi;

			const Fluid U0(U_import[i]);
			assert(U0[Fluid::MASS] > 0.0);

			const Fluid W0 = U0.to_primitive(cell_list[i].Volume);

			const real dt = scheduler.get_dt(ptcl_import[i].rung);

			Fluid Uc;
			for (int k = 0; k < Fluid::NFLUID; k++)
				Uc[k] = U0[k] + dotUi[k]*dt;

			assert(Uc[Fluid::MASS] > 0.0);

			if (!compute_update_prob(Uc, i))
			{
				assert(Uc[Fluid::MASS] > 0.0);

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
				const real pres = compute_pressure(W0);
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

				const double cr = 0.5;
				const real l_min = std::pow(cell_list[i].Volume/(4.0*M_PI/3.0), 1.0/3.0);
				Uc[Fluid::MPSI] *= std::exp(-cr*ch_sig/l_min * dt);


				assert(Uc[Fluid::MASS] > 0.0);

				// GRAVITY
#if 1
				const Fluid &Ui = U0;
				const vec3 dmom = (Ui[Fluid::MASS]*ptcl_import[i].acc0 + Uc[Fluid::MASS]*ptcl_import[i].acc1)*0.5*dt;
				Uc[Fluid::MOMX] += dmom.x;
				Uc[Fluid::MOMY] += dmom.y;
				Uc[Fluid::MOMZ] += dmom.z;

				const vec3 mom0(Ui[Fluid::MOMX], Ui[Fluid::MOMY], Ui[Fluid::MOMZ]);
				const vec3 mom1(Uc[Fluid::MOMX], Uc[Fluid::MOMY], Uc[Fluid::MOMZ]);
				Uc[Fluid::ENER] += (mom0*ptcl_import[i].acc0 + mom1*ptcl_import[i].acc1)*dt*0.5;
#endif
			}

			Fluid W1 = Uc.to_primitive(cell_list[i].Volume);

			assert(W1[Fluid::DENS] > 0.0);

			W1[Fluid::ETHM] = compute_ethm_update(W1, i);

			if (W1[Fluid::ETHM] <= 0.0)
			{
				W1[Fluid::ETHM] = compute_ethm_from_entropy(W1);
				assert(false);
				W1[Fluid::ETHM] = W0[Fluid::ETHM];
			}
			else
				W1[Fluid::ENTR] = compute_entropy_from_ethm(W1);

			assert(W1[Fluid::ETHM] > 0.0);

			U_import [i]      = W1.to_conservative(cell_list[i].Volume);
			dU_import[i].U    = 0.0;
			dU_import[i].divB = dt > 0.0 ? divB/dt/cell_list[i].Volume : divB;
			dU_import[i].gradPsi = 0.0;

			const Fluid Wi = U_import[i].to_primitive(ptcl_import[i].volume);
			assert(Wi[Fluid::ETHM] > 0);
		}
	}

	/////////////////////////////////
	/////////////////////////////////
	/////////////////////////////////

	void system::compute_update0()
	{
#if 0
		{
			int nactive_site = site_active_list.size();
			for (int isite = 0; isite < nactive_site; isite++)
			{
				const int i = site_active_list[isite];
				Wrec_import[i].t = 0.0;
			}
			return;
		}
#endif

		const int nimport = site_import.size();
		std::vector<FluidD> dUdt(nimport, 0.0);

		const int nactive_face = nface_active;
		for (int iface = 0; iface < nactive_face; iface++)
		{
			const Face &face = *face_active_list[iface];

			int i  = face.s1;
			int j  = face.s2; 

#if 1
			assert(i != j);
			assert(i >= 0);
			assert(i < nimport);
			assert(j >= 0);
			assert(j < nimport);
#endif

#if 1
			assert(ptcl_import[i].boundary != -1289);
			assert(ptcl_import[j].boundary != -1289);
#endif
			//// j <-| i |-> j  ....

			const real area(face.area());
			assert(area > 0.0);
			vec3 normal(face.n * (1.0/area));

			Fluid_st Wst_i = Wst_import[i];
			Fluid_st Wst_j = Wst_import[j];

			const vec3 &ipos1 = Wst_i.pos;
			vec3 dri = face.centroid - ipos1;

			assert(std::abs(dri.x) < 0.5*global_domain_size.x);
			assert(std::abs(dri.y) < 0.5*global_domain_size.y);
			assert(std::abs(dri.z) < 0.5*global_domain_size.z);


			if (Wst_i.bnd != Particle::NO_BOUNDARY)
			{
				std::swap(Wst_i, Wst_j);
				std::swap(i, j);
				normal *= -1.0;
				dri    *= -1.0;
			}

			if (Wst_i.bnd != Particle::NO_BOUNDARY)
				continue;

			const real dsh = dri * normal;
			const real dsl = 2.0 * dsh;
#if 1
			assert(dsl != 0.0);
			assert(dsl  > 0.0);
#else
			assert(dsl > 0.0);
#endif

			const vec3 drj = dri - normal * dsl;

			const Fluid &Wi =  Wst_i.w;
			const Fluid &Wj = (Wst_j.bnd == Particle::NO_BOUNDARY) ? Wst_j.w : Wst_i.w;

			assert(Wi[Fluid::DENS] > 0.0);
			assert(Wi[Fluid::ETHM] > 0.0);
			assert(Wj[Fluid::DENS] > 0.0);
			assert(Wj[Fluid::ETHM] > 0.0);

			const vec3 &vi = Wst_i.vel;
			const vec3 &vj = Wst_j.vel;

			const vec3  vij = vj + vi;
			const vec3 dvij = vj - vi;

			const vec3  fij = dri - normal * dsh;
			const real  fdv = dvij * fij;
			const real ids2 = -0.5*fdv * ((dsh > 0.0) ? 1.0/sqr(dsh) : 0.0);
			const vec3  wij = vij*0.5 +  ids2 * dri;

			real psi, Bn;
			Fluid flux;
			compute_flux(wij, Wi, Wj, Wst_j.bnd, normal, psi, Bn, flux);

#if 0
			// resistivity
			const vec3 ResFlux = normal * (Wst_i.J + Wst_j.J)*0.5;
			flux[Fluid::BY] += ResFlux.y; 
			flux[Fluid::BZ] += ResFlux.z; 
#endif

			for (int k = 0; k < Fluid::NFLUID; k++)
			{
				dUdt[i].U[k] -= flux[k] * area;
				dUdt[j].U[k] += flux[k] * area;
			}
			dUdt[i].gradPsi += normal * (psi * area);
			dUdt[j].gradPsi -= normal * (psi * area);
			dUdt[i].divB    += Bn * area;
			dUdt[j].divB    -= Bn * area;
		}

		int nactive_site = site_active_list.size();
		int ethm_failed = 0;
		int adjust_dt_cnt = 0;
		for (int isite = 0; isite < nactive_site; isite++)
		{
			const int i = site_active_list[isite];
			if (Wst_import[i].bnd != Particle::NO_BOUNDARY) continue;

			const Fluid dotUi = dUdt[i].U;
			const real divB   = dUdt[i].divB;
			const vec3 dPsi   = dUdt[i].gradPsi;

			const Fluid U0(U_import[i]);
			assert(cell_list[i].Volume > 0.0);
			const Fluid W0 = U0.to_primitive(cell_list[i].Volume);

			assert(U0[Fluid::MASS] > 0.0);
			bool flag_cnt = false;
			//      while(std::abs(dotUi[Fluid::MASS]*scheduler.get_dt(ptcl_import[i].rung)) >= 0.5*U0[Fluid::MASS])
			while(U0[Fluid::MASS] + dotUi[Fluid::MASS]*scheduler.get_dt(ptcl_import[i].rung) <= 0.0*U0[Fluid::MASS])
			{
				if (!flag_cnt)
				{ 
					adjust_dt_cnt++;
					flag_cnt = true;
				}
				ptcl_import[i].rung++;
			}
			assert(ptcl_import[i].rung < Scheduler::RUNGMAX);

			const real dt = scheduler.get_dt(ptcl_import[i].rung);

			Wrec_import[i].tlast = -dt;

			Fluid Uc;
			for (int k = 0; k < Fluid::NFLUID; k++)
				Uc[k] = U0[k] + dotUi[k]*dt;

			assert(Uc[Fluid::MASS] > 1.0e-32);


			if (!compute_update_prob(Uc, i))
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
				const real pres = compute_pressure(W0);
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

				const double cr = 0.5;
				const real l_min = std::pow(cell_list[i].Volume/(4.0*M_PI/3.0), 1.0/3.0);
				Uc[Fluid::MPSI] *= std::exp(-cr*ch_sig/l_min * dt);


				assert(Uc[Fluid::MASS] > 0.0);

				// GRAVITY
#if 1
				const Fluid &Ui = U0;
				const vec3 dmom = (Ui[Fluid::MASS]*ptcl_import[i].acc0 + Uc[Fluid::MASS]*ptcl_import[i].acc1)*0.5*dt;
				Uc[Fluid::MOMX] += dmom.x;
				Uc[Fluid::MOMY] += dmom.y;
				Uc[Fluid::MOMZ] += dmom.z;

				const vec3 mom0(Ui[Fluid::MOMX], Ui[Fluid::MOMY], Ui[Fluid::MOMZ]);
				const vec3 mom1(Uc[Fluid::MOMX], Uc[Fluid::MOMY], Uc[Fluid::MOMZ]);
				Uc[Fluid::ENER] += (mom0*ptcl_import[i].acc0 + mom1*ptcl_import[i].acc1)*dt*0.5;
#endif
			}

			Fluid W1 = Uc.to_primitive(cell_list[i].Volume);

			assert(W1[Fluid::DENS] > 0.0);

			W1[Fluid::ETHM] = compute_ethm_update(W1, i);

			if (W1[Fluid::ETHM] <= 0.0)
			{
				W1[Fluid::ETHM] = compute_ethm_from_entropy(W1);
				ethm_failed++;
				//assert(false);
				W1[Fluid::ETHM] = W0[Fluid::ETHM];
			}
			else
				W1[Fluid::ENTR] = compute_entropy_from_ethm(W1);

			assert(W1[Fluid::ETHM] > 0.0);


#if 0
			{
				assert(dt > 0.0);
				Uc = W1.to_conservative(cell_list[i].Volume);
				const real idt = 1.0/(dt * cell_list[i].Volume);
				Fluid dUdt;
				for (int k = 0; k < Fluid::NFLUID; k++)
					dUdt[k] = (Uc[k] - U0[k]) * idt;

				Wrec_import[i].t[Fluid::DENS] = dUdt[Fluid::DENS];
				const vec3 vdot = dUdt.get_vel()*(1.0/U0[Fluid::DENS]) - U0.get_vel()*(dUdt[Fluid::DENS]/sqr(U0[Fluid::DENS]));
				Wrec_import[i].t[Fluid::VELX] = vdot.x;
				Wrec_import[i].t[Fluid::VELY] = vdot.y;
				Wrec_import[i].t[Fluid::VELZ] = vdot.z;
				Wrec_import[i].t[Fluid::ETHM] =  dUdt[Fluid::ENER]
					- (U0.get_vel()*dUdt.get_vel())/U0[Fluid::DENS]
					+  U0.get_vel().norm2()/(2.0*sqr(U0[Fluid::DENS]))*dUdt[Fluid::DENS]
					-  U0.get_B() * dUdt.get_B();
				Wrec_import[i].t[Fluid::BX]   = dUdt[Fluid::BX];
				Wrec_import[i].t[Fluid::BY]   = dUdt[Fluid::BY];
				Wrec_import[i].t[Fluid::BZ]   = dUdt[Fluid::BZ];
				Wrec_import[i].t[Fluid::PSI]  = dUdt[Fluid::PSI ]/U0[Fluid::DENS] - U0[Fluid::PSI ]/sqr(U0[Fluid::DENS])*dUdt[Fluid::DENS];
				Wrec_import[i].t[Fluid::ENTR] = dUdt[Fluid::ENTR]/U0[Fluid::DENS] - U0[Fluid::ENTR]/sqr(U0[Fluid::DENS])*dUdt[Fluid::DENS];
				for (int k = 0; k < Fluid::NSCALARS; k++) 
					Wrec_import[i].t.scal(k) = dUdt.scal(k)/U0[Fluid::DENS] - U0.scal(k)/sqr(U0[Fluid::DENS])*dUdt[Fluid::DENS];
			}
#else
			assert(dt > 0.0);
			const real idt = 1.0/dt;
			for (int k = 0; k < Fluid::NFLUID; k++)
				Wrec_import[i].t[k] = (W1[k] - W0[k]) * idt;
#endif
		}
#if 1
		if (adjust_dt_cnt > 0)
			fprintf(stderr, " update0: myproc= %d, adjust_dt_cnt= %d [ %g ] \n",
					myproc, adjust_dt_cnt, nactive_site > 0 ? 1.0*adjust_dt_cnt/nactive_site : 0.0);
#endif
	}


	inline void system::compute_flux(
			const vec3  &wij,
			const Fluid &Wi,
			const Fluid &Wj,
			const int    jbnd,
			const vec3  &normal,
			real &psi_ij,
			real &Bn_ij,
			Fluid &flux)
	{
		asm("#COMPUTEFLUX-CALL");
		// dpflop 24
		// spflop 120
		const real  ex = normal.x;
		const real  ey = normal.y;
		const real  ez = normal.z;
		const real  ds = std::sqrt(ex*ex + ey*ey);
		const real ids = (ds != 0.0) ? 1.0/ds : 0.0;
		//spflop 3 + 1*20 + 1*10
		const real cosph = (ids == 0.0) ? 1.0 : ex * ids;
		const real sinph = (ids == 0.0) ? 0.0 : ey * ids;
		const real costh = ez;
		const real sinth = ds;
		//spflop 4
		const real Axx =  cosph*sinth;
		const real Axy =  sinth*sinph;
		const real Axz =  costh;
		const real Ayx = -sinph;
		const real Ayy =  cosph;
		const real Ayz =  0.0;
		const real Azx = -costh*cosph;
		const real Azy = -costh*sinph;
		const real Azz =  sinth;
		//spflop 7
		real dens_L =     Wi[Fluid::DENS];
		real ethm_L =     Wi[Fluid::ETHM];
		real velx_L = Axx*Wi[Fluid::VELX] + Axy*Wi[Fluid::VELY] + Axz*Wi[Fluid::VELZ];
		real vely_L = Ayx*Wi[Fluid::VELX] + Ayy*Wi[Fluid::VELY] + Ayz*Wi[Fluid::VELZ];
		real velz_L = Azx*Wi[Fluid::VELX] + Azy*Wi[Fluid::VELY] + Azz*Wi[Fluid::VELZ];
		real Bx_L   = Axx*Wi[Fluid::BX]   + Axy*Wi[Fluid::BY]   + Axz*Wi[Fluid::BZ];
		real By_L   = Ayx*Wi[Fluid::BX]   + Ayy*Wi[Fluid::BY]   + Ayz*Wi[Fluid::BZ];
		real Bz_L   = Azx*Wi[Fluid::BX]   + Azy*Wi[Fluid::BY]   + Azz*Wi[Fluid::BZ];
		real psiL   =     Wi[Fluid::PSI];
		real entrL  =     Wi[Fluid::ENTR];
		//dpflop 30
		real dens_R =     Wj[Fluid::DENS];
		real ethm_R =     Wj[Fluid::ETHM];
		real velx_R = Axx*Wj[Fluid::VELX] + Axy*Wj[Fluid::VELY] + Axz*Wj[Fluid::VELZ];
		real vely_R = Ayx*Wj[Fluid::VELX] + Ayy*Wj[Fluid::VELY] + Ayz*Wj[Fluid::VELZ];
		real velz_R = Azx*Wj[Fluid::VELX] + Azy*Wj[Fluid::VELY] + Azz*Wj[Fluid::VELZ];
		real Bx_R   = Axx*Wj[Fluid::BX]   + Axy*Wj[Fluid::BY]   + Axz*Wj[Fluid::BZ];
		real By_R   = Ayx*Wj[Fluid::BX]   + Ayy*Wj[Fluid::BY]   + Ayz*Wj[Fluid::BZ];
		real Bz_R   = Azx*Wj[Fluid::BX]   + Azy*Wj[Fluid::BY]   + Azz*Wj[Fluid::BZ];
		real psiR   =     Wj[Fluid::PSI];
		real entrR  =     Wj[Fluid::ENTR];
		//dpflop 30

		real pres_L = compute_pressure(Wi);
		real pres_R = compute_pressure(Wj);

		const real wn_ij = vec3(Axx, Axy, Axz)*wij;
		switch(jbnd)
		{
			case Particle::REFLECTING:
				velx_R = -velx_R;
				break;
#if 1
			case Particle::DIOD:
				velx_R = std::max(velx_R - 1.2*wn_ij, 0.0) + 1.2*wn_ij;
				break;
#endif
		};


		//dpflop  4
		const real cfl2 = (gamma_gas*pres_L + sqr(Bx_L) + sqr(By_L) + sqr(Bz_L))/dens_L;
		const real cfr2 = (gamma_gas*pres_R + sqr(Bx_R) + sqr(By_R) + sqr(Bz_R))/dens_R;
		const real chL = std::sqrt(cfl2 + 0.0*sqr(velx_L - wn_ij));
		const real chR = std::sqrt(cfr2 + 0.0*sqr(velx_R - wn_ij));

		Bn_ij  = (Bx_L*chL + Bx_R*chR +         (psiL - psiR))/(chL + chR);
		psi_ij = (psiL*chR + psiR*chL + chR*chL*(Bx_L - Bx_R))/(chL + chR);
		//spflop 15 + 1*10


		//dpflop  5
		Fluid F;
		asm("#CALL_RP-BEG");
		riemann_solver(
				F, 
				Bn_ij, wn_ij,
				dens_L, pres_L, ethm_L, velx_L, vely_L, velz_L, By_L, Bz_L,
				dens_R, pres_R, ethm_R, velx_R, vely_R, velz_R, By_R, Bz_R);
		asm("#CALL_RP-END");

		const real iAxx =  cosph*sinth;
		const real iAxy = -sinph;
		const real iAxz = -costh*cosph;
		const real iAyx =  sinth*sinph;
		const real iAyy =  cosph;
		const real iAyz = -costh*sinph;
		const real iAzx =  costh;
		const real iAzy =  real(0.0);
		const real iAzz =  sinth;
		//sflop 7
		flux[Fluid::MASS] =      F[Fluid::MASS];
		flux[Fluid::ENER] =      F[Fluid::ENER];
		flux[Fluid::MOMX] = iAxx*F[Fluid::MOMX] + iAxy*F[Fluid::MOMY] + iAxz*F[Fluid::MOMZ];
		flux[Fluid::MOMY] = iAyx*F[Fluid::MOMX] + iAyy*F[Fluid::MOMY] + iAyz*F[Fluid::MOMZ];
		flux[Fluid::MOMZ] = iAzx*F[Fluid::MOMX] + iAzy*F[Fluid::MOMY] + iAzz*F[Fluid::MOMZ];
		flux[Fluid::WBX ] = iAxx*F[Fluid::WBX ] + iAxy*F[Fluid::WBY ] + iAxz*F[Fluid::WBZ ];
		flux[Fluid::WBY ] = iAyx*F[Fluid::WBX ] + iAyy*F[Fluid::WBY ] + iAyz*F[Fluid::WBZ ];
		flux[Fluid::WBZ ] = iAzx*F[Fluid::WBX ] + iAzy*F[Fluid::WBY ] + iAzz*F[Fluid::WBZ ];
		flux[Fluid::MPSI] =      F[Fluid::MASS] * (F[Fluid::MASS] > 0.0 ? psiL  : psiR);
		//    flux[Fluid::MPSI] = F[Fluid::MASS] * psi_ij;
		flux[Fluid::MENTR] =     F[Fluid::MASS] * (F[Fluid::MASS] > 0.0 ? entrL : entrR);
		//dpflop 32

		for (int k = 0; k < Fluid::NSCALARS; k++) 
			flux.scal(k) = F[Fluid::MASS] * (F[Fluid::MASS] > 0.0 ? Wi.scal(k) : Wj.scal(k));
		//dpflop 6
		asm("#COMPUTEFLUX-EXIT");
	}

	inline void system::riemann_solver(
			Fluid &flux,
			const real Bx,
			const real w,
			const real dens_L, const real pres_L, const real ethm_L,
			const real velx_L, const real vely_L, const real velz_L, const real By_L, const real Bz_L,
			const real dens_R, const real pres_R, const real ethm_R,
			const real velx_R, const real vely_R, const real velz_R, const real By_R, const real Bz_R
			) const
	{
		asm("#RP_BEG");

#ifdef HLLD
		{
			const real signBx = Bx == 0.0 ? 0.0 : (Bx > 0.0 ? +1.0 : -1.0);

			const real momx_L = dens_L*velx_L;
			const real momy_L = dens_L*vely_L;
			const real momz_L = dens_L*velz_L;

			const real momx_R = dens_R*velx_R;
			const real momy_R = dens_R*vely_R;
			const real momz_R = dens_R*velz_R;
			//flop 6
			const real chalf = 0.5;

			const real B2_L   = sqr(Bx)     + sqr(By_L)   + sqr(Bz_L);
			const real v2_L   = sqr(velx_L) + sqr(vely_L) + sqr(velz_L);
			const real etot_L = ethm_L      + chalf*(dens_L*v2_L + B2_L);
			//flop 14
			const real B2_R   = sqr(Bx)     + sqr(By_R)   + sqr(Bz_R);
			const real v2_R   = sqr(velx_R) + sqr(vely_R) + sqr(velz_R);
			const real etot_R = ethm_R      + chalf*(dens_R*v2_R + B2_R);
			//flop 14
			const real gpl  = gamma_gas * pres_L;
			const real gpr  = gamma_gas * pres_R;
			const real gpbl = gpl + B2_L;
			const real gpbr = gpr + B2_R;
			//flop 4
			//flop//flop//flop//flop//flop

			const real cfl2  = gpbl / dens_L;
			const real cfr2  = gpbr / dens_R;
#if 1
			const real cfmax = std::sqrt(std::max(cfl2, cfr2));
#else
			const real cfmax = std::max(std::sqrt(std::max(cfl2, cfr2)),
					std::max(std::abs(velx_L), std::abs(velx_R)));
#endif

			//flop 6 + 1*10 + 2*10
			const real S_L = std::min(velx_L, velx_R) - cfmax;
			const real S_R = std::max(velx_L, velx_R) + cfmax;

			//flop 4
			//flop//flop//flop//flop//flop

			const real pT_L = pres_L + chalf * B2_L;
			const real pT_R = pres_R + chalf * B2_R;
			//flop 6
			const real iSM = 1.0/((S_R - velx_R)*dens_R - (S_L - velx_L)*dens_L);
			const real S_M  =   iSM * ((S_R - velx_R)*momx_R - (S_L - velx_L)*momx_L - pT_R + pT_L);

			//flop 13 + 1*10
			const real ipTs = 1.0/((S_R - velx_R)*dens_R - (S_L - velx_L)*dens_L);
			const real pT_s = ipTs * ((S_R - velx_R)*dens_R*pT_L - (S_L - velx_L)*dens_L*pT_R +
					dens_L*dens_R*(S_R - velx_R)*(S_L - velx_L)*(velx_R - velx_L));
			//flop 22 + 1*10
			const real velx_L_s  = S_M;
			const real velx_L_ss = S_M;
			const real velx_R_s  = S_M;
			const real velx_R_ss = S_M;
			const real B2x       = Bx*Bx;
			//flop 24
			const real iSLmSM   = 1.0/(S_L - S_M);
			const real iSRmSM   = 1.0/(S_R - S_M);
			const real dens_L_s = dens_L * (S_L - velx_L) * iSLmSM;
			const real dens_R_s = dens_R * (S_R - velx_R) * iSRmSM;
			const real divL     = dens_L * (S_L - velx_L)*(S_L - S_M) - B2x;
			const real divR     = dens_R * (S_R - velx_R)*(S_R - S_M) - B2x;
			const real idivL    = (divL != 0.0) ? 1.0/divL : 0.0;
			const real idivR    = (divR != 0.0) ? 1.0/divR : 0.0;
			const real vely_L_s = vely_L - Bx*By_L*(S_M - velx_L) * idivL;
			const real velz_L_s = velz_L - Bx*Bz_L*(S_M - velx_L) * idivL;
			const real   By_L_s = By_L * (dens_L*sqr(S_L - velx_L) - B2x) * idivL;
			const real   Bz_L_s = Bz_L * (dens_L*sqr(S_L - velx_L) - B2x) * idivL;
			//flop 42 + 2*10 + 2*20
			const real vely_R_s = vely_R - Bx*By_R*(S_M - velx_R) * idivR;
			const real velz_R_s = velz_R - Bx*Bz_R*(S_M - velx_R) * idivR;
			const real   By_R_s = By_R * (dens_R*sqr(S_R - velx_R) - B2x) * idivR;
			const real   Bz_R_s = Bz_R * (dens_R*sqr(S_R - velx_R) - B2x) * idivR;
			//flop 22
			const real   vB_L   = velx_L  *Bx + vely_L  *By_L   + velz_L  *Bz_L;
			const real   vB_L_s = velx_L_s*Bx + vely_L_s*By_L_s + velz_L_s*Bz_L_s;
			const real etot_L_s = ((S_L - velx_L)*etot_L - pT_L*velx_L + pT_s*S_M + Bx*(vB_L - vB_L_s)) * iSLmSM;
			//flop 20
			const real   vB_R   = velx_R  *Bx + vely_R  *By_R   + velz_R  *Bz_R;
			const real   vB_R_s = velx_R_s*Bx + vely_R_s*By_R_s + velz_R_s*Bz_R_s;
			const real etot_R_s = ((S_R - velx_R)*etot_R - pT_R*velx_R + pT_s*S_M + Bx*(vB_R - vB_R_s)) * iSRmSM;
			//flop 20
			const real dens_L_ss = dens_L_s;
			const real dens_R_ss = dens_R_s;
			const real sDens_L_s = std::sqrt(dens_L_s);
			const real sDens_R_s = std::sqrt(dens_R_s);
			//flop 2*10
			const real    S_L_s  = S_M - std::abs(Bx/sDens_L_s);
			const real    S_R_s  = S_M + std::abs(Bx/sDens_R_s);
			//flop 6 + 2*10
			const real idsqroot  = 1.0/(sDens_L_s + sDens_R_s);
			const real  vely_ss = idsqroot*(sDens_L_s*vely_L_s + sDens_R_s*vely_R_s + (By_R_s - By_L_s)*signBx);
			const real  velz_ss = idsqroot*(sDens_L_s*velz_L_s + sDens_R_s*velz_R_s + (Bz_R_s - Bz_L_s)*signBx);
			//flop 15 + 1*10
			const real By_ss = idsqroot*(sDens_L_s*By_R_s + sDens_R_s*By_L_s + sDens_L_s*sDens_R_s*(vely_R_s - vely_L_s)*signBx);
			const real Bz_ss = idsqroot*(sDens_L_s*Bz_R_s + sDens_R_s*Bz_L_s + sDens_L_s*sDens_R_s*(velz_R_s - velz_L_s)*signBx);
			//flop 18
			const real vely_L_ss = vely_ss;
			const real velz_L_ss = velz_ss;
			const real   By_L_ss = By_ss;
			const real   Bz_L_ss = Bz_ss;

			const real vely_R_ss = vely_ss;
			const real velz_R_ss = velz_ss;
			const real   By_R_ss = By_ss;
			const real   Bz_R_ss = Bz_ss;

			const real vB_L_ss   = velx_L_ss*Bx + vely_L_ss*By_L_ss + velz_L_ss*Bz_L_ss;
			const real etot_L_ss = etot_L_s - sDens_L_s*(vB_L_s - vB_L_ss)*signBx;
			//flop 8
			const real vB_R_ss   = velx_R_ss*Bx + vely_R_ss*By_R_ss + velz_R_ss*Bz_R_ss;
			const real etot_R_ss = etot_R_s + sDens_R_s*(vB_R_s - vB_R_ss)*signBx;
			//flop 8
			const real Fdens_L = dens_L*velx_L;
			const real Fmomx_L = momx_L*velx_L + pT_L - B2x;
			const real Fmomy_L = momy_L*velx_L        - Bx*By_L;
			const real Fmomz_L = momz_L*velx_L        - Bx*Bz_L;
			const real Fetot_L = etot_L*velx_L + pT_L*velx_L - Bx*vB_L; 
			//flop 15
			const real Fdens_R = dens_R*velx_R;
			const real Fmomx_R = momx_R*velx_R + pT_R - B2x;
			const real Fmomy_R = momy_R*velx_R        - Bx*By_R;
			const real Fmomz_R = momz_R*velx_R        - Bx*Bz_R;
			const real Fetot_R = etot_R*velx_R + pT_R*velx_R - Bx*vB_R;
			//flop 15
			const real momx_L_s  = dens_L_s *velx_L_s;
			const real momy_L_s  = dens_L_s *vely_L_s;
			const real momz_L_s  = dens_L_s *velz_L_s;

			const real momx_L_ss = dens_L_ss*velx_L_ss;
			const real momy_L_ss = dens_L_ss*vely_L_ss;
			const real momz_L_ss = dens_L_ss*velz_L_ss;

			const real momx_R_s  = dens_R_s *velx_R_s;
			const real momy_R_s  = dens_R_s *vely_R_s;
			const real momz_R_s  = dens_R_s *velz_R_s;

			const real momx_R_ss = dens_R_ss*velx_R_ss;
			const real momy_R_ss = dens_R_ss*vely_R_ss;
			const real momz_R_ss = dens_R_ss*velz_R_ss;

			const real Fby_L  = By_L*velx_L - Bx * vely_L;
			const real Fbz_L  = Bz_L*velx_L - Bx * velz_L;

			const real Fby_R  = By_R*velx_R - Bx * vely_R;
			const real Fbz_R  = Bz_R*velx_R - Bx * velz_R;
			//flop 24
			const bool w_lt_SL  = w <= S_L;
			const bool w_le_SLs = w <= S_L_s;
			const bool w_le_SM  = w <= S_M;
			const bool w_le_SRs = w <= S_R_s;
			const bool w_le_SR  = w <= S_R;
			//flop 5
			const real fdens = (w_le_SM) ? Fdens_L : Fdens_R;
			const real fetot = (w_le_SM) ? Fetot_L : Fetot_R;
			const real fmomx = (w_le_SM) ? Fmomx_L : Fmomx_R;
			const real fmomy = (w_le_SM) ? Fmomy_L : Fmomy_R;
			const real fmomz = (w_le_SM) ? Fmomz_L : Fmomz_R;
			const real fby   = (w_le_SM) ? Fby_L   :   Fby_R;
			const real fbz   = (w_le_SM) ? Fbz_L   :   Fbz_R;
			//flop 7
			const real cnull(0.0);

			const real a = 
				(w_le_SLs ? 	cnull : 
				 (w_le_SM ? S_L_s - w : 
					(w_le_SRs ? S_R_s - w :	cnull)
				 ));

			//flop 5
			const real b  = 
				(w_lt_SL ? cnull :
				 (w_le_SLs ? +S_L - w :
					(w_le_SM  ? -S_L_s + S_L :
					 (w_le_SRs ? -S_R_s + S_R :
						(w_le_SR  ? +S_R - w : cnull)
					 ))));

			//flop 11
			const real c = 
				(w_lt_SL ? -w : 
				 (w_le_SLs ? -S_L : 
					(w_le_SM  ? -S_L : 
					 (w_le_SRs ? -S_R : 
						(w_le_SR  ?	-S_R : -w)
					 ))));

			//flop 11
			const real dens    = (w_le_SM ? dens_L :    dens_R   );
			const real dens_s  = (w_le_SM ? dens_L_s :  dens_R_s );
			const real dens_ss = (w_le_SM ? dens_L_ss : dens_R_ss);

			const real etot    = (w_le_SM ? etot_L :    etot_R   );
			const real etot_s  = (w_le_SM ? etot_L_s :  etot_R_s );
			const real etot_ss = (w_le_SM ? etot_L_ss : etot_R_ss);

			const real momx    = (w_le_SM ? momx_L :    momx_R   );
			const real momx_s  = (w_le_SM ? momx_L_s :  momx_R_s );
			const real momx_ss = (w_le_SM ? momx_L_ss : momx_R_ss);

			const real momy    = (w_le_SM ? momy_L :    momy_R   );
			const real momy_s  = (w_le_SM ? momy_L_s :  momy_R_s );
			const real momy_ss = (w_le_SM ? momy_L_ss : momy_R_ss);

			const real momz    = (w_le_SM ? momz_L :    momz_R   );
			const real momz_s  = (w_le_SM ? momz_L_s :  momz_R_s );
			const real momz_ss = (w_le_SM ? momz_L_ss : momz_R_ss);

			const real by    = (w_le_SM ? By_L :    By_R   );
			const real by_s  = (w_le_SM ? By_L_s :  By_R_s );
			const real by_ss = (w_le_SM ? By_L_ss : By_R_ss);

			const real bz    = (w_le_SM ? Bz_L :    Bz_R   );
			const real bz_s  = (w_le_SM ? Bz_L_s :  Bz_R_s );
			const real bz_ss = (w_le_SM ? Bz_L_ss : Bz_R_ss);

			//flop 21
			flux[Fluid::MASS] = fdens + a * dens_ss + b * dens_s + c * dens;
			flux[Fluid::ENER] = fetot + a * etot_ss + b * etot_s + c * etot;
			flux[Fluid::MOMX] = fmomx + a * momx_ss + b * momx_s + c * momx;
			flux[Fluid::MOMY] = fmomy + a * momy_ss + b * momy_s + c * momy;
			flux[Fluid::MOMZ] = fmomz + a * momz_ss + b * momz_s + c * momz;
			flux[Fluid::WBY ] = fby   + a * by_ss   + b * by_s   + c * by;
			flux[Fluid::WBZ ] = fbz   + a * bz_ss   + b * bz_s   + c * bz;
			flux[Fluid::WBX ] = -w * Bx;

		}
#else // HLLE
		{

			const real momx_L = dens_L*velx_L;
			const real momy_L = dens_L*vely_L;
			const real momz_L = dens_L*velz_L;

			const real momx_R = dens_R*velx_R;
			const real momy_R = dens_R*vely_R;
			const real momz_R = dens_R*velz_R;

			const real chalf = real(0.5);

			const real B2_L   = sqr(Bx)     + sqr(By_L)   + sqr(Bz_L);
			const real v2_L   = sqr(velx_L) + sqr(vely_L) + sqr(velz_L);
			const real etot_L = ethm_L      + chalf*(dens_L*v2_L + B2_L);
			//flop 14
			const real B2_R   = sqr(Bx)     + sqr(By_R)   + sqr(Bz_R);
			const real v2_R   = sqr(velx_R) + sqr(vely_R) + sqr(velz_R);
			const real etot_R = ethm_R      + chalf*(dens_R*v2_R + B2_R);

			const real pT_L = pres_L + chalf * B2_L;
			const real pT_R = pres_R + chalf * B2_R;

#if 1
			const real gpl  = real(gamma_gas) * pres_L;
			const real gpr  = real(gamma_gas) * pres_R;
			const real gpbl = gpl + B2_L;
			const real gpbr = gpr + B2_R;

			const real cfl2  = gpbl/ dens_L;
			const real cfr2  = gpbr/ dens_R;
			const real cfmax = std::max(std::sqrt(std::max(cfl2, cfr2)),
					std::max(std::abs(velx_L), std::abs(velx_R)));
#else
			const real cf2   = (std::max(pres_L, pres_R) + std::max(B2_L, B2_R)*0.5)/std::min(dens_L, dens_R);
			const real cfmax = std::sqrt(cf2);
#endif


			//flop 6 + 1*10 + 2*10
#if 1
			const real S_L = std::min(velx_L, velx_R) - cfmax;
			const real S_R = std::max(velx_L, velx_R) + cfmax;
#else
			const real S_R = std::max(std::abs(velx_L), std::abs(velx_R)) + cfmax;
			const real S_L = -S_R;
#endif

			real vB_L = velx_L*Bx + vely_L*By_L   + velz_L*Bz_L;
			real vB_R = velx_R*Bx + vely_R*By_R   + velz_R*Bz_R;

			real Fdens_L = dens_L*velx_L;
			real Fmomx_L = momx_L*velx_L + pT_L - Bx*Bx;
			real Fmomy_L = momy_L*velx_L        - Bx*By_L; 
			real Fmomz_L = momz_L*velx_L        - Bx*Bz_L;
			real Fetot_L = etot_L*velx_L + pT_L*velx_L - Bx*vB_L;
			real Fby_L   = By_L  *velx_L - Bx*vely_L;
			real Fbz_L   = Bz_L  *velx_L - Bx*velz_L;

			real Fdens_R = dens_R*velx_R;
			real Fmomx_R = momx_R*velx_R + pT_R - Bx*Bx;
			real Fmomy_R = momy_R*velx_R        - Bx*By_R; 
			real Fmomz_R = momz_R*velx_R        - Bx*Bz_R; 
			real Fetot_R = etot_R*velx_R + pT_R*velx_R - Bx*vB_R;
			real Fby_R   = By_R  *velx_R        - Bx*vely_R;
			real Fbz_R   = Bz_R  *velx_R        - Bx*velz_R;


			real U_dens = (S_R*dens_R - S_L*dens_L + Fdens_L - Fdens_R)/(S_R - S_L);
			real U_momx = (S_R*momx_R - S_L*momx_L + Fmomx_L - Fmomx_R)/(S_R - S_L);
			real U_momy = (S_R*momy_R - S_L*momy_L + Fmomy_L - Fmomy_R)/(S_R - S_L);
			real U_momz = (S_R*momz_R - S_L*momz_L + Fmomz_L - Fmomz_R)/(S_R - S_L);
			real U_etot = (S_R*etot_R - S_L*etot_L + Fetot_L - Fetot_R)/(S_R - S_L);

			flux[Fluid::MASS] = (S_R*Fdens_L - S_L*Fdens_R + S_L*S_R*(dens_R - dens_L))/(S_R - S_L);
			flux[Fluid::MOMX] = (S_R*Fmomx_L - S_L*Fmomx_R + S_L*S_R*(momx_R - momx_L))/(S_R - S_L);
			flux[Fluid::MOMY] = (S_R*Fmomy_L - S_L*Fmomy_R + S_L*S_R*(momy_R - momy_L))/(S_R - S_L);
			flux[Fluid::MOMZ] = (S_R*Fmomz_L - S_L*Fmomz_R + S_L*S_R*(momz_R - momz_L))/(S_R - S_L);
			flux[Fluid::ENER] = (S_R*Fetot_L - S_L*Fetot_R + S_L*S_R*(etot_R - etot_L))/(S_R - S_L);

			real U_by = (S_R*By_R - S_L*By_L + Fby_L - Fby_R)/(S_R - S_L);
			real U_bz = (S_R*Bz_R - S_L*Bz_L + Fbz_L - Fbz_R)/(S_R - S_L);

			flux[Fluid::WBY] = (S_R*Fby_L - S_L*Fby_R + S_L*S_R*(By_R - By_L))/(S_R - S_L);
			flux[Fluid::WBZ] = (S_R*Fbz_L - S_L*Fbz_R + S_L*S_R*(Bz_R - Bz_L))/(S_R - S_L);

			flux[Fluid::MASS] -= w*U_dens;
			flux[Fluid::MOMX] -= w*U_momx;
			flux[Fluid::MOMY] -= w*U_momy;
			flux[Fluid::MOMZ] -= w*U_momz;
			flux[Fluid::ENER] -= w*U_etot;
			flux[Fluid::WBY]  -= w*U_by;
			flux[Fluid::WBZ]  -= w*U_bz;

			flux[Fluid::WBX]   = -w*Bx;


		}
#endif

		//flop 45
		asm("#RP_END");
	}

};
