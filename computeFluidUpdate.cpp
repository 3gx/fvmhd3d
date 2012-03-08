#include "fvmhd3d.h"

#define DEDNER
//#define DEDNER_ENERGY
//#define REMOVE_MONPOLE

namespace fvmhd3d
{
  struct FluidUpdateStruct
  {
    Fluid_rec Wrec;
    FluidD    dU;
  };

  void System::computeFluidUpdate_main(CkCallback &cb)
  {
    computeFluidUpdate_completeCb = cb;
#if 1
    systemProxy[thisIndex].slopeLimiter(CkCallback(CkIndex_System::computeFluidUpdateI(), systemProxy[thisIndex]));
#else
    computeFluidUpdateII();
#endif
  }

  void System::computeFluidUpdateI()
  {
#if 1
    for (int i = 0; i < nactive_loc; i++)
    {
			const real f0 = Problem_enforce_limiter(i);
      for (int k = 0; k < Fluid::NFLUID; k++)
      {
        real f = f0*(*slopeLimiter_act[i])[k];
        assert(f >= 0.0);
        assert(f <= 1.0);
        //        assert(f > 0.9);

        assert(f == 1.0);
        Wrec_act[i]->x[k] *= f;
        Wrec_act[i]->y[k] *= f;
        Wrec_act[i]->z[k] *= f;
        Wrec_act[i]->t[k] *= f;
      }
#if 0
      const vec3 acc = mesh_act[i]->acc1;
      Wrec_act[i]->t[Fluid::VELX] += acc.x;
      Wrec_act[i]->t[Fluid::VELY] += acc.y;
      Wrec_act[i]->t[Fluid::VELZ] += acc.z;
#endif
    }
#endif
    for (int i = nactive_loc; i < (const int)ptcl_act.size(); i++)
      if (ptcl_act[i]->is_ngb())
      {
				const real f0 = Problem_enforce_limiter(i);
				for (int k = 0; k < Fluid::NFLUID; k++)
				{
					real f = f0*(*slopeLimiter_act[i])[k];
					assert(f >= 0.0);
					assert(f <= 1.0);
					Wrec_act[i]->x[k] *= f;
					Wrec_act[i]->y[k] *= f;
					Wrec_act[i]->z[k] *= f;
					Wrec_act[i]->t[k] *= f;
				}
#if 0
				const vec3 acc = mesh_act[i]->acc1;
				Wrec_act[i]->t[Fluid::VELX] += acc.x;
				Wrec_act[i]->t[Fluid::VELY] += acc.y;
				Wrec_act[i]->t[Fluid::VELZ] += acc.z;
#endif
			}
		computeFluidUpdateII();

	}

	void System::computeFluidUpdateII()
	{
		computeFluidUpdate();
		computeFluidUpdate_sendIncUpdate();
	}

	void System::computeFluidUpdate_complete()
	{
		assert(false);
	}

	void System::computeFluidUpdate()
	{
		std::vector<FluidUpdateStruct> Wact(ptcl_act.size());
		for (int i = 0; i < (const int)ptcl_act.size(); i++)
		{
			Wact[i].Wrec = *Wrec_act[i];
			Wact[i].dU   = 0.0;
//      assert(mesh_act[i]->tlast == Wrec_act[i]->tend);
		}

		const int nactive_face = face_active_list.size();
		for (int iface = 0; iface < nactive_face; iface++)
		{
			const Face &face = *face_active_list[iface];

			int i  = face.s1;
			int j  = face.s2;

			real area  (face.area());
			vec3 normal(face.n * (1.0/area));

			//// j <-| i |-> j  ....

			Fluid_rec Wrec_i = Wact[i].Wrec;
			Fluid_rec Wrec_j = Wact[j].Wrec;

			__builtin_prefetch(face_active_list[iface+1]);

			const vec3 ipos = Wrec_i.pos;
			vec3 dri = face.centroid - ipos;

			if (Wrec_i.bnd != MeshPoint::NO_BOUNDARY)
			{
				std::swap(Wrec_i, Wrec_j);
				std::swap(i, j);
				normal *= -1.0;
				dri    *= -1.0;
			}

			assert(Wrec_i.bnd == MeshPoint::NO_BOUNDARY);

			const real dtI = Wrec_i.tend - t_global;
			const real dtJ = Wrec_j.tend - t_global;
			const real dt  = std::min(dtI, dtJ);
			assert(dt > 0.0);
#if 1
			const real dti = dtI - 0.5*dt;
			const real dtj = dtJ - 0.5*dt;
#else
			const real dti = Wrec_i.tend - mesh_act[i]->tbeg - 0.5*dt;
			const real dtj = Wrec_j.tend - mesh_act[j]->tbeg - 0.5*dt;
#endif

			const real dsh = dri * normal;
			const real dsl = 2.0 * dsh;
			assert(dsl != 0.0);
			assert(dsl  > 0.0);

			const vec3 drj = dri - normal * dsl;

			Fluid &Wi = Wrec_i.w;
			Fluid &Wj = Wrec_j.w;



#if 0
			if ((Wrec_j.bnd != MeshPoint::NO_BOUNDARY) && (Wrec_j.bnd != MeshPoint::INFLOW))
			{
				for (int k = 0; k < Fluid::NFLUID; k++) 
					Wi[k] = Wj[k] = Wrec_i.w[k];
			}
			else
#endif
        if (Wrec_j.bnd == MeshPoint::NO_BOUNDARY)
			{
				for (int k = 0; k < Fluid::NFLUID; k++) 
				{
#if 1
					Wi[k] += vec3(Wrec_i.x[k], Wrec_i.y[k], Wrec_i.z[k]) * dri + Wrec_i.t[k]*dti;
					Wj[k] += vec3(Wrec_j.x[k], Wrec_j.y[k], Wrec_j.z[k]) * drj + Wrec_j.t[k]*dtj;
#endif
				}
#if 0
				const vec3 acc = (Wrec_i.acc + Wrec_j.acc)*0.5;
				Wi[Fluid::VELX] += acc.x*dti;
				Wi[Fluid::VELY] += acc.y*dti;
				Wi[Fluid::VELZ] += acc.z*dti;
				Wj[Fluid::VELX] += acc.x*dtj;
				Wj[Fluid::VELY] += acc.y*dtj;
				Wj[Fluid::VELZ] += acc.z*dtj;
#endif
        assert(Wi[Fluid::DENS] > 0.0);
        assert(Wi[Fluid::ETHM] > 0.0);
        assert(Wj[Fluid::DENS] > 0.0);
        assert(Wj[Fluid::ETHM] > 0.0);
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
#if 0
      const real  fdv = dvij * fij;
      const real ids2 = -0.5*fdv * ((dsh > 0.0) ? 1.0/sqr(dsh) : 0.0);
      const vec3  wij = vij*0.5 +  ids2 * dri;
#else
      const real fdv  = dvij*fij;
      const vec3 rji  = dsl*normal;
      const vec3 dwij = fdv * rji/sqr(dsl);
      const vec3 wij  = vij*0.5 - dwij;
#endif

      real psi, Bn;
      Fluid flux;
      computeFlux(wij, Wi, Wj, Wrec_j.bnd, normal, psi, Bn, flux);

#if 0
      // resistivity
      const vec3 ResFlux = normal.cross((Wrec_i.etaJ*Wrec_i.J + Wrec_j.etaJ*Wrec_j.J)*0.5);
      flux[Fluid::BX] += ResFlux.x; 
      flux[Fluid::BY] += ResFlux.y; 
      flux[Fluid::BZ] += ResFlux.z; 
#endif

      const real tarea = dt * area;
      const vec3 dPsi  = normal * (psi * tarea);


      for (int k = 0; k < Fluid::NFLUID; k++)
        Wact[i].dU.U[k]  -= flux[k] * tarea;
      Wact[i].dU.gradPsi += dPsi;
      Wact[i].dU.divB    += Bn * tarea;   


      for (int k = 0; k < Fluid::NFLUID; k++)
        Wact[j].dU.U[k]  += flux[k] * tarea;
      Wact[j].dU.gradPsi -= dPsi;
      Wact[j].dU.divB    -= Bn * tarea;   
    }

    for (int i = 0; i < (const int)ptcl_act.size(); i++)
    {
      //      if (i >= nactive_loc && i < nactive_loc + nimport_loc) continue;
      //      if (i >= nactive_loc + nimport_loc) continue;
      dU_act[i]->divB    += Wact[i].dU.divB;
      dU_act[i]->gradPsi += Wact[i].dU.gradPsi;
      for (int k = 0; k < Fluid::NFLUID; k++)
        dU_act[i]->U[k] += Wact[i].dU.U[k];
    }

  }  /* end System::computeFluidUpdate */

  ///////////////////////////
  ///////////////////////////
  ///////////////////////////

  void System::computeFluidUpdate_sendIncUpdate()
  {
#if 0
    contribute(computeFluidUpdate_completeCb);
    return;
#endif
    assert(nSend_cntr == 0);

    const int ibeg = active_list.size();
    const int iend = ptcl_act.size();

    std::vector< std::pair<int, int> > request_list;
    request_list.reserve(iend - ibeg);
    for (int i = ibeg; i < iend; i++)
      if (ptcl_act[i]->is_ngb())
        if (!ptcl_act[i]->is_active())
          request_list.push_back(std::make_pair(ptcl_act[i]->chare(), i));

    std::sort(request_list.begin(), request_list.end(), std_pair_first_sort());

    const int nrequest = request_list.size();
    CkVec< pair<int, FluidD> > data2send;
    data2send.reserve(nrequest);
    request_list.push_back(std::make_pair(-1,-1));

    for (int i = 0; i < nrequest; i++)
    {
      const int iElement = request_list[i].first;
      const int iId      = request_list[i].second;
      //			dU_act[iId]->data[0] = mesh_act[iId]->tlast;
      assert(!ptcl_act[iId]->is_active());
      data2send.push_back(std::make_pair(ptcl_act[iId]->id(), *dU_act[iId]));
      assert(iElement >= 0);
      assert(iElement < numElements);
      if (iElement != request_list[i+1].first && data2send.size() > 0)
      {
        if (iElement != thisIndex)
        {
          nSend_cntr++;
          systemProxy[iElement].computeFluidUpdate_recvIncUpdate(data2send, thisIndex);
        }
        else
        {
          computeFluidUpdate_recvIncUpdate(data2send, thisIndex);
        }
        data2send.clear();
      }
    }

    nSend_cntr++;
    computeFluidUpdate_recvIncUpdateTicket();
  }

  void System::computeFluidUpdate_recvIncUpdateTicket()
  {
    nSend_cntr--;
    assert(nSend_cntr >= 0);
    if (nSend_cntr == 0)
    {
#if 0
      contribute(CkCallback(CkIndex_System::computeFluidUpdate_complete(), thisProxy));
#else
      contribute(computeFluidUpdate_completeCb);
#endif
    }
  }

  void System::computeFluidUpdate_recvIncUpdate(const CkVec< pair<int, FluidD> > &dataRecv, const int recvIndex)
  {
    const int nrecv = dataRecv.size();
    for (int i = 0; i < nrecv; i++)
    {
      const int local_id = dataRecv[i].first;
      const FluidD   &dU = dataRecv[i].second;
      assert(local_id >= 0);
      assert(local_id <  local_n);
      assert(!ptcl_list[local_id].is_active());

      //      mesh_pnts[local_id].tlast = dU.data[0];


      if (!mesh_pnts[local_id].is_boundary())
      {
        for (int k = 0; k < Fluid::NFLUID; k++)
          dU_list[local_id].U[k]  += dU.U[k];
        dU_list[local_id].divB    += dU.divB;
        dU_list[local_id].gradPsi += dU.gradPsi;
      }

    }
    if (thisIndex != recvIndex)
      systemProxy[recvIndex].computeFluidUpdate_recvIncUpdateTicket();
  }

  void System::computePrimitives(const bool predictor_flag)
  {
    int ethm_failed = 0;
    std::vector<std::pair<int, int> > failed_fluid;

    for (int i = 0; i < nactive_loc; i++)
    {
      if (mesh_act[i]->is_boundary())
      { 
        Wrec_act[i]->w = U_act[i]->to_primitive(mesh_act[i]->Volume);
        Problem_set_boundary(i);
        Wrec_act[i]->x = 0.0;
        Wrec_act[i]->y = 0.0;
        Wrec_act[i]->z = 0.0;
        Wrec_act[i]->t = 0.0;
        mesh_act[i]->Volume = cell_list[i].Volume;
        *U_act[i] = Wrec_act[i]->w.to_conservative(mesh_act[i]->Volume);
        continue;
      }

      Fluid Uc;
      for (int k = 0; k < Fluid::NFLUID; k++)
        Uc[k] = (*U_act[i])[k] + dU_act[i]->U[k];

      const real divB = dU_act[i]->divB;
      const vec3 dPsi = dU_act[i]->gradPsi;

      const real dt = t_global - mesh_act[i]->tbeg;
      divBi[ptcl_act[i]->id()] = dt > 0.0 ? divB/dt/cell_list[i].Volume : divB;

      if (Uc[Fluid::MASS] <= 0.0)
      {
        fprintf(stderr, "update: dens: thisIndex= %d idx= %d i= %d : dt= %g [ %g ] pos= %g  vel= %g \n",
            thisIndex, mesh_act[i]->idx, i, 
            mesh_act[i]->tend - mesh_act[i]->tbeg,
            scheduler.get_dt(mesh_act[i]->rung),
            mesh_act[i]->pos.abs(), mesh_act[i]->vel.abs());
        failed_fluid.push_back(std::make_pair(i, 0)); // FAILED_MASS));
        continue;
      }
      assert(Uc[Fluid::MASS] > 0.0);    // check for NaN

      if (!Problem_compute_update(Uc, i))
      {
        // --- interpolate in primitives ---
        //
        const real dth = 0.5 * dt;
        Fluid Wh; 
        for (int k = 0; k < Fluid::NFLUID; k++)
          Wh[k] = Wrec_act[i]->w[k] + Wrec_act[i]->t[k] * dth;

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
        const real pres = Problem_compute_pressure(Wh);
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
        if (!predictor_flag)
        {
          const Fluid &Ui = *U_act[i];
          const vec3 dmom = 
            (Ui[Fluid::MASS]*mesh_act[i]->acc0 + Uc[Fluid::MASS]*mesh_act[i]->acc1)*0.5*dt;
          Uc[Fluid::MOMX] += dmom.x;
          Uc[Fluid::MOMY] += dmom.y;
          Uc[Fluid::MOMZ] += dmom.z;

          const vec3 mom0(Ui[Fluid::MOMX], Ui[Fluid::MOMY], Ui[Fluid::MOMZ]);
          const vec3 mom1(Uc[Fluid::MOMX], Uc[Fluid::MOMY], Uc[Fluid::MOMZ]);
          Uc[Fluid::ENER] += (mom0*mesh_act[i]->acc0 + mom1*mesh_act[i]->acc1)*dt*0.5;
        }
#endif

#ifdef DEDNER
        const double cr = 0.5;
        const real l_min = std::pow(cell_list[i].Volume/(4.0*M_PI/3.0), 1.0/3.0);
        Uc[Fluid::MPSI] *= std::exp(-cr * ch_sig/l_min * dt);
#endif
      }

      assert(Uc[Fluid::MASS] > 0.0);
      Fluid W1 = Uc.to_primitive(cell_list[i].Volume);

      assert(W1[Fluid::DENS] > 0.0);
      W1[Fluid::ETHM] = Problem_compute_ethm_update(W1, i);

      if (W1[Fluid::ETHM] <= 0.0)
      {
        assert(W1[Fluid::ENTR] > 0.0);
        W1[Fluid::ETHM] = Problem_compute_ethm_from_entropy(W1);
        ethm_failed++;
#if 0
        failed_fluid.push_back(std::make_pair(i, 0)); // FAILED_MASS));
        fprintf(stderr, "update: ethm: thisIndex= %d idx= %d i= %d : dt= %g [ %g ] pos= %g  vel= %g \n",
            thisIndex, mesh_act[i]->idx, i, 
            mesh_act[i]->tend - mesh_act[i]->tbeg,
            scheduler.get_dt(mesh_act[i]->rung),
            mesh_act[i]->pos.abs(), mesh_act[i]->vel.abs());
#endif
#if 0
        failed_fluid.push_back(std::make_pair(i, 0)); // FAILED_MASS));
        continue;
#endif
      }
      else
        W1[Fluid::ENTR] = Problem_compute_entropy_from_ethm(W1);


      assert(W1[Fluid::ETHM] > 0.0);
#if 1
      if (predictor_flag)
      {
        if (dt > 0.0)
          for (int k = 0; k < Fluid::NFLUID; k++)
            Wrec_act[i]->t[k] = (W1[k] - Wrec_act[i]->w[k])/dt;
        else
          Wrec_act[i]->t = 0.0;
      }
#endif

      if (!predictor_flag)
      {
        mesh_act[i]->Volume = cell_list[i].Volume;
        Wrec_act[i]->w = W1;
        *U_act  [i]    = W1.to_conservative(cell_list[i].Volume);
        *dU_act [i]    = 0.0;

        const Fluid Wi = U_act[i]->to_primitive(cell_list[i].Volume);
        assert(Wi[Fluid::ETHM] > 0);
      }
    }

    for (int ifailed = 0; ifailed < (const int)failed_fluid.size(); ifailed++)
    {
      const int i = failed_fluid[ifailed].first;
      const real dt = t_global - mesh_act[i]->tbeg;
      Fluid W1;
      for (int k = 0; k < Fluid::NFLUID; k++)
        W1[k] = Wrec_act[i]->w[k] + Wrec_act[i]->t[k]*dt;
      W1 = (W1.to_conservative(mesh_act[i]->Volume)).to_primitive(cell_list[i].Volume);

      assert(W1[Fluid::MASS] > 0.0);

      mesh_act[i]->Volume = cell_list[i].Volume;

      assert(W1[Fluid::DENS] > 0.0);
      W1[Fluid::ETHM] = Problem_compute_ethm_update(W1, i);

      if (W1[Fluid::ETHM] <= 0.0)
      {
        W1[Fluid::ETHM] = Problem_compute_ethm_from_entropy(W1);
        //        ethm_failed++;
        assert(false);
      }
      else
        W1[Fluid::ENTR] = Problem_compute_entropy_from_ethm(W1);

      assert(W1[Fluid::ETHM] > 0.0);

      Wrec_act[i]->w = W1;
      *U_act[i]      = W1.to_conservative(cell_list[i].Volume);
      *dU_act[i]     = 0.0;

      const Fluid Wi = U_act[i]->to_primitive(mesh_act[i]->Volume);
      assert(Wi[Fluid::ETHM] > 0);
    }
  }

  ///////////////////
  /////////////////////
  ///////////////////

  void System::exchangeNewDt()
  {
    assert(nSend_cntr == 0);

    const int ibeg = active_list.size();
    const int iend = ptcl_act.size();

    std::vector< std::pair<int, int> > request_list;
    request_list.reserve(iend - ibeg);
    for (int i = ibeg; i < iend; i++)
      if (ptcl_act[i]->is_ngb())
        if (!ptcl_act[i]->is_active())
          request_list.push_back(std::make_pair(ptcl_act[i]->chare(), i));

    std::sort(request_list.begin(), request_list.end(), std_pair_first_sort());

    const int nrequest = request_list.size();
    CkVec< pair<int, real> > data2send;
    data2send.reserve(nrequest);
    request_list.push_back(std::make_pair(-1,-1));

    for (int i = 0; i < nrequest; i++)
    {
      const int iElement = request_list[i].first;
      const int iId      = request_list[i].second;
      data2send.push_back(std::make_pair(ptcl_act[iId]->id(), mesh_act[iId]->dt_new));
      assert(iElement >= 0);
      assert(iElement < numElements);
      if (iElement != request_list[i+1].first && data2send.size() > 0)
      {
        if (thisIndex != iElement)
        {
          nSend_cntr++;
          systemProxy[iElement].exchangeNewDt_recv(data2send, thisIndex);
        }
        else
        {
          exchangeNewDt_recv(data2send, thisIndex);
        }
        data2send.clear();
      }
    }

    nSend_cntr++;
    exchangeNewDt_recvTicket();
  }

  void System::exchangeNewDt_recvTicket()
  {
    nSend_cntr--;
    assert(nSend_cntr >= 0);
    if (nSend_cntr == 0)
      contribute(exchangeNewDt_returnCb);
  }

  void System::exchangeNewDt_recv(const CkVec< pair<int, real> > &dataRecv, const int recvIndex)
  {
    const int nrecv = dataRecv.size();
    for (int i = 0; i < nrecv; i++)
    {
      const int local_id = dataRecv[i].first;
      assert(local_id >= 0);
      assert(local_id <  local_n);
      assert(!ptcl_list[local_id].is_active());

      const real dt_new = dataRecv[i].second;
      const real tbeg   = mesh_pnts[local_id].tbeg;
      const real tend   = mesh_pnts[local_id].tend;
      if (tend - tbeg > dt_new)
      {
        const int old_rung = mesh_pnts[local_id].rung;
        mesh_pnts[local_id].rung = scheduler.get_rung(0.5*dt_new);
        mesh_pnts[local_id].tend = t_global + scheduler.get_dt(mesh_pnts[local_id].rung);

        const int new_rung = mesh_pnts[local_id].rung;
        scheduler.move<true>(old_rung, new_rung, local_id);
      }
    }
    if (thisIndex != recvIndex)
      systemProxy[recvIndex].exchangeNewDt_recvTicket();
  }

  ///////////////////////

  void System::exchangeNewDtAct()
  {
    assert(nSend_cntr == 0);

    const int ibeg = active_list.size();
    const int iend = ptcl_act.size();

    std::vector< std::pair<int, int> > request_list;
    request_list.reserve(iend - ibeg);
    for (int i = ibeg; i < iend; i++)
      if (ptcl_act[i]->is_ngb())
        request_list.push_back(std::make_pair(ptcl_act[i]->chare(), i));

    std::sort(request_list.begin(), request_list.end(), std_pair_first_sort());

    const int nrequest = request_list.size();
    CkVec< pair<int, real> > data2send;
    data2send.reserve(nrequest);
    request_list.push_back(std::make_pair(-1,-1));

    for (int i = 0; i < nrequest; i++)
    {
      const int iElement = request_list[i].first;
      const int iId      = request_list[i].second;
      data2send.push_back(std::make_pair(ptcl_act[iId]->id(), mesh_act[iId]->dt_new));
      assert(iElement >= 0);
      assert(iElement < numElements);
      if (iElement != request_list[i+1].first && data2send.size() > 0)
      {
        if (thisIndex != iElement)
        {
          nSend_cntr++;
          systemProxy[iElement].exchangeNewDtAct_recv(data2send, thisIndex);
        }
        else
        {
          exchangeNewDtAct_recv(data2send, thisIndex);
        }
        data2send.clear();
      }
    }

    nSend_cntr++;
    exchangeNewDtAct_recvTicket();
  }

  void System::exchangeNewDtAct_recvTicket()
  {
    nSend_cntr--;
    assert(nSend_cntr >= 0);
    if (nSend_cntr == 0)
      contribute(exchangeNewDtAct_returnCb);
  }

  void System::exchangeNewDtAct_recv(const CkVec< pair<int, real> > &dataRecv, const int recvIndex)
  {
    const int nrecv = dataRecv.size();
    for (int i = 0; i < nrecv; i++)
    {
      const int local_id = dataRecv[i].first;
      assert(local_id >= 0);
      assert(local_id <  local_n);

      const real dt_new = dataRecv[i].second;
      if (ptcl_list[local_id].is_active())
      {
        if (dt_new < mesh_pnts[local_id].dt_new)
          mesh_pnts[local_id].dt_new = dt_new;
      }
      else
      {
#if 0
        const real tbeg   = mesh_pnts[local_id].tbeg;
        const real tend   = mesh_pnts[local_id].tend;
        if (tend - tbeg > dt_new)
        {
          const int old_rung = mesh_pnts[local_id].rung;
          mesh_pnts[local_id].rung = scheduler.get_rung(0.5*dt_new);
          mesh_pnts[local_id].tend = t_global + scheduler.get_dt(mesh_pnts[local_id].rung);

          const int new_rung = mesh_pnts[local_id].rung;
          scheduler.move<true>(old_rung, new_rung, local_id);
        }
#endif
      }
    }
    if (thisIndex != recvIndex)
      systemProxy[recvIndex].exchangeNewDtAct_recvTicket();
  }

}
