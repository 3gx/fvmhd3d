#include "fvmhd3d.h"

namespace fvmhd3d {

  inline real Phi(const real x, const real eps)
  {
    assert(x >= 0.0);
#if 0
    return std::min(1.0, x);
#else
    const real fac  = 0.2;
    const real eps1 = 1.0e-16;
    return (x*x + 2.0*x + fac*eps1)/(x*x + x + 2.0 + fac*eps1);
#endif
  }
  
  inline real limiter(const real dWi, const real dWj_min, const real dWj_max, const real eps)
  {
    if (std::abs(dWi) <= 1.0e-16*std::max(std::abs(dWj_max), std::abs(dWj_min)))  return 1.0;
    else if (dWi < 0.0) return Phi(dWj_min/dWi, eps);
    else                return Phi(dWj_max/dWi, eps);
  }

  void system::compute_reconstruction()
  {
#if 0
    {
      const int nactive_site = site_active_list.size();
      for (int isite = 0; isite < nactive_site; isite++)
      {
        const int i = site_active_list[isite];
        assert(!ptcl_import[i].is_boundary());
        Wrec_import[i] = Fluid_rec(Wst_import[i].w);
      }
      return;
    }
#endif
    const int nactive_site = site_active_list.size();
      
    std::vector<vec3 > centroid_list;
    for (int isite = 0; isite < nactive_site; isite++)
    {
      const int i = site_active_list[isite];
      if (ptcl_import[i].is_boundary())
      {
        Wrec_import[i] = Fluid_rec(Wst_import[i].w);
        continue;
      }

      const Cell     &ci   = cell_list[i];
      const Fluid_st &Wst  = Wst_import[i];
      const Fluid    &Wi   = Wst.w;
      const vec3     &ipos = Wst.pos;

      Fluid Wx(0.0);
      Fluid Wy(0.0);
      Fluid Wz(0.0);
      Fluid dWj_min(0.0), dWj_max(0.0);

      centroid_list.clear();

      //
      // compute_gradient ...
      //
      const int nface = ci.faces().size();	
      bool zero_flag = false;
      assert(nface > 0);
      vec3 sum_centroid(0.0);
      const vec3 ivel(Wi.get_vel());
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
        centroid_list.push_back(centroid);
        sum_centroid += centroid;
        const real area     = face.area();
        assert(area > 0.0);
        const vec3 normal   = face.n * ((centroid * face.n < 0.0) ? (-1.0/area) : (1.0/area));
        const real dsh  = centroid * normal;
        assert(dsh > 0.0);

        const real ids  = (dsh != 0.0) ? 0.5/dsh : 0.0;
        const real idsA = area * ids;

        const vec3 drh = normal * dsh;
        const vec3 fij = centroid - drh;

        const int j = face.ngb<false>(i);
        assert(j >= 0);

        if (Wst_import[j].bnd != Particle::NO_BOUNDARY)
        {
          zero_flag = true;
          break;
        }
        const Fluid &Wj = Wst_import[j].w;

        for (int k = 0; k < Fluid::NFLUID; k++) 
        {
          const real  sum = Wj[k] + Wi[k];
          const real diff = Wj[k] - Wi[k];
          dWj_min[k] = std::min(dWj_min[k], diff);
          dWj_max[k] = std::max(dWj_max[k], diff);
          Wx[k] += (sum*drh.x + diff*fij.x)*idsA;
          Wy[k] += (sum*drh.y + diff*fij.y)*idsA;
          Wz[k] += (sum*drh.z + diff*fij.z)*idsA;
        }
      }
      if (zero_flag)
      {
        Wrec_import[i].w = Wi;
        Wrec_import[i].x = 0.0;
        Wrec_import[i].y = 0.0;
        Wrec_import[i].z = 0.0;
        Wrec_import[i].t = 0.0;
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
      
      //
      // Limit gradients in primitive variables
      //

      Fluid limiter_min(1.0);
      const real sum_dt = scheduler.get_dt(ptcl_import[i].rung);
      for (int iface = 0; iface < nface; iface++)
        for (int k = 0; k < Fluid::NFLUID; k++) 
        {
          const real dWx = vec3(Wx[k], Wy[k], Wz[k])*centroid_list[iface];
          limiter_min[k] = std::min(limiter_min[k], limiter(dWx, dWj_min[k], dWj_max[k], ptcl_import[i].volume));
#if 0
          const real dWt = Wrec_import[i].t[k]*sum_dt;
          const real dW  = dWx + dWt;
          limiter_min[k] = std::min(limiter_min[k], limiter(dW , dWj_min[k], dWj_max[k], ptcl_import[i].volume));
#endif
#if 0
          limiter_min[k] = std::min(limiter_min[k], limiter(dWt, dWj_min[k], dWj_max[k], ptcl_import[i].volume));
#endif
        }
      

      for (int k = 0; k < Fluid::NFLUID; k++) 
      {
        const real tau = limiter_min[k];
  
        Wx[k] *= tau;
        Wy[k] *= tau;
        Wz[k] *= tau;

        
        Wrec_import[i].t[k] *= tau;

        Wrec_import[i].w[k] = Wi[k];
        Wrec_import[i].x[k] = Wx[k];
        Wrec_import[i].y[k] = Wy[k];
        Wrec_import[i].z[k] = Wz[k];
      }

#if 1
      // Positivity for the *-state, Berthon, JCoP 2006 & Waagan, JCoP 2009
      {
        int k = Fluid::DENS;
        real dWx = vec3(Wx[k], Wy[k], Wz[k])*sum_centroid;
        real dWp = dWx + Wrec_import[i].t[k]*sum_dt;
        real dW  = std::max(dWx, dWp);
        
        if (dW > 0.0)
        {
          real alpha = std::min(1.0, 0.9*Wi[k]/dW);
          assert(alpha >= 0.0);
          Wrec_import[i].x[k] *= alpha;
          Wrec_import[i].y[k] *= alpha;
          Wrec_import[i].z[k] *= alpha;
          Wrec_import[i].t[k] *= alpha;
        }
      }
#endif
    }
  }

  void system::slope_limiter()
  {
#if 1
    const int nimport = site_import.size();
    std::vector<std::pair<Fluid, Fluid> > dWj_minmax(nimport);
    std::vector<Fluid> min_limiter(nimport);

    std::vector<Fluid_st> Wst(nimport);

    for (int i = 0; i < nimport; i++)
      if (site_import[i].is_active() || site_import[i].is_local_ngb())
      {
        Wst[i].dt = t_global - Wrec_import[i].tlast;
        Wst[i].bnd = Wrec_import[i].bnd;      
        for (int k = 0; k < Fluid::NFLUID; k++)
        {
          min_limiter[i][k] = 1.0;
          dWj_minmax [i].first[k] = dWj_minmax[i].second[k] = 0.0;
          Wst        [i].w[k]     = Wrec_import[i].w[k];
        }
      }
      else
        Wst[i].dt = -1.0;

    const int nactive_face = nface_active;
    std::vector< std::pair<real, real> > dtij_list(nactive_face);
    for (int iface = 0; iface < nactive_face; iface++)
    {
      const Face &face = *face_active_list[iface];

      int i  = face.s1;
      int j  = face.s2; 

#if 1
      assert(ptcl_import[i].boundary != -1289);
      assert(ptcl_import[j].boundary != -1289);
#endif

      const Fluid_st &Wi = Wst[i];
      const Fluid_st &Wj = Wst[j];

      if (Wi.bnd != Particle::NO_BOUNDARY || Wj.bnd != Particle::NO_BOUNDARY)
        continue;
      
      __builtin_prefetch(face_active_list[iface+1]);
      __builtin_prefetch(&dWj_minmax[i]);
      __builtin_prefetch(&dWj_minmax[j]);

      const real dtI = Wi.dt;
      const real dtJ = Wj.dt;
      const real dth = 0.5*std::min(dtI, dtJ);
      assert(dth > 0.0);
      dtij_list[iface] = std::make_pair(dtI - dth, dtJ - dth);


      for (int k = 0; k < Fluid::NFLUID; k++)
      {
        const real diff = Wj.w[k] - Wi.w[k];
        dWj_minmax[i].first [k] = std::min(dWj_minmax[i].first [k],  diff);
        dWj_minmax[i].second[k] = std::max(dWj_minmax[i].second[k],  diff);
        dWj_minmax[j].first [k] = std::min(dWj_minmax[j].first [k], -diff);
        dWj_minmax[j].second[k] = std::max(dWj_minmax[j].second[k], -diff);
      }
    }

    static std::vector<Fluid> dWj_min_local, dWj_max_local;
    dWj_min_local.resize(local_n);
    dWj_max_local.resize(local_n);
    std::vector< int > idx_send[NMAXPROC];
    std::vector< int > idx_recv[NMAXPROC];
    std::vector<Fluid> fluid_send[NMAXPROC];
    std::vector<Fluid> fluid_recv[NMAXPROC];
    static int nsend_out[NMAXPROC], nrecv_out[NMAXPROC];
    static int nsend_in [NMAXPROC], nrecv_in [NMAXPROC];
#if 1
    {
      const int nimport = site_import.size();

      for (int i = 0; i < nimport; i++)
      {
        const SiteImport &s = site_import[i];
        if (!(s.is_local_ngb() || s.is_hasngb())) continue;
        assert(site_import[i].is_active() || site_import[i].is_local_ngb());

        const int proc = s.proc();
        const int id   = s.id();
        ptcl_import[i].local_id = id;
        assert(proc < nproc);

        const int  local_id = id;
        const int remote_id = i;

				if (proc != myproc)
				{
					idx_send  [proc].push_back( local_id);
					idx_send  [proc].push_back(remote_id);
					fluid_send[proc].push_back(dWj_minmax[i].first);
					fluid_send[proc].push_back(dWj_minmax[i].second);
				}
				else
					for (int k = 0; k < Fluid::NFLUID; k++)
						dWj_min_local[local_id][k] = dWj_max_local[local_id][k] = 0.0;
			}

			myMPI::all2all<true >(  idx_send,   idx_recv, myproc, nproc, 2, nsend_out, nrecv_out);
			myMPI::all2all<false>(fluid_send, fluid_recv, myproc, nproc, 2, nsend_out, nrecv_out);

			for (int p = 0; p < nproc; p++)
				for (int q = 0; q < (const int)idx_recv[p].size(); q += 2)
				{
					const int local_id = idx_recv[p][q + 0];
					assert(local_id >= 0.0);
					assert(local_id < (int)local_n);

					for (int k = 0; k < Fluid::NFLUID; k++)
						dWj_min_local[local_id][k] = dWj_max_local[local_id][k] = 0.0;
				}

			for (int p = 0; p < nproc; p++)
			{
				for (int q = 0; q < (const int)idx_recv[p].size(); q += 2)
				{
					assert (p != myproc);
					const int local_id =   idx_recv[p][q + 0];
					const Fluid &minWj = fluid_recv[p][q + 0];
					const Fluid &maxWj = fluid_recv[p][q + 1];
					for (int k = 0; k < Fluid::NFLUID; k++)
					{
						dWj_min_local[local_id][k] = std::min(dWj_min_local[local_id][k], minWj[k]);
						dWj_max_local[local_id][k] = std::max(dWj_max_local[local_id][k], maxWj[k]);
					}
				}
				idx_send  [p].clear();
				fluid_send[p].clear();
			}
     
      for (int i = 0; i < nimport; i++)
      {
        const SiteImport &s = site_import[i];
        if (!(s.is_local_ngb() || s.is_hasngb())) continue;

        const int proc = s.proc();
        const int id   = s.id();

        const int  local_id = id;
        const int remote_id = i;

				if (proc != myproc) continue;
				const Fluid &minWj = dWj_minmax[remote_id].first;
				const Fluid &maxWj = dWj_minmax[remote_id].second;
				for (int k = 0; k < Fluid::NFLUID; k++)
				{
					dWj_min_local[local_id][k] = std::min(dWj_min_local[local_id ][k], minWj[k]);
					dWj_max_local[local_id][k] = std::max(dWj_max_local[local_id ][k], maxWj[k]);
					dWj_minmax[remote_id].first [k] = std::min(dWj_minmax[remote_id].first [k], dWj_min_local[local_id][k]);
					dWj_minmax[remote_id].second[k] = std::max(dWj_minmax[remote_id].second[k], dWj_max_local[local_id][k]);
				}
			}

			for (int p = 0; p < nproc; p++)
				for (int q = 0; q < (const int)idx_recv[p].size(); q += 2)
				{
					const int local_id  = idx_recv[p][q + 0];
					const int remote_id = idx_recv[p][q + 1];

					idx_send  [p].push_back( local_id);
					idx_send  [p].push_back(remote_id);
					fluid_send[p].push_back(dWj_min_local[local_id]);
					fluid_send[p].push_back(dWj_max_local[local_id]);
				}

			myMPI::all2all<true >(  idx_send,   idx_recv, myproc, nproc, 2,  nsend_in, nrecv_in);
			myMPI::all2all<false>(fluid_send, fluid_recv, myproc, nproc, 2, nsend_in, nrecv_in);

			for (int p = 0; p < nproc; p++)
			{
				for (int q = 0; q < (const int)idx_recv[p].size(); q += 2)
				{
					assert (p != myproc);
					const int local_id  = idx_recv[p][q + 0];
					const int remote_id = idx_recv[p][q + 1];
					assert(remote_id >= 0);
					assert(remote_id < nimport);
					assert(local_id == ptcl_import[remote_id].local_id);

					const Fluid &minWj = fluid_recv[p][q + 0];
					const Fluid &maxWj = fluid_recv[p][q + 1];

					for (int k = 0; k < Fluid::NFLUID; k++)
					{
						dWj_minmax[remote_id].first [k] = std::min(dWj_minmax[remote_id].first [k], minWj[k]);
						dWj_minmax[remote_id].second[k] = std::max(dWj_minmax[remote_id].second[k], maxWj[k]);
					}
				}
				idx_send  [p].clear();
				idx_recv  [p].clear();
				fluid_send[p].clear();
				fluid_recv[p].clear();
			}
		}
#endif

		for (int iface = 0; iface < nactive_face; iface++)
		{
			const Face &face = *face_active_list[iface];

			int i  = face.s1;
			int j  = face.s2;

			const Fluid_rec &Wrec_i = Wrec_import[i];
			const Fluid_rec &Wrec_j = Wrec_import[j];

#if 1
			assert(ptcl_import[i].boundary != -1289);
			assert(ptcl_import[j].boundary != -1289);
#endif

			if (Wrec_i.bnd != Particle::NO_BOUNDARY || Wrec_j.bnd != Particle::NO_BOUNDARY)
			{
				for (int k = 0; k < Fluid::NFLUID; k++)
					min_limiter[i][k] = min_limiter[j][k] = 0;
				continue;
			}

			__builtin_prefetch(face_active_list[iface+1]);
			__builtin_prefetch(&min_limiter[i]);
			__builtin_prefetch(&min_limiter[j]);


			const vec3 dri = face.centroid - Wrec_import[i].pos;
#if 1
			assert(std::abs(dri.x) < 0.5*global_domain_size.x);
			assert(std::abs(dri.y) < 0.5*global_domain_size.y);
			assert(std::abs(dri.z) < 0.5*global_domain_size.z);
#endif


			const real area  (face.area());
			const vec3 normal(face.n * (1.0/area));
			const real dsh = dri * normal;
			const real dsl = 2.0 * dsh;
			assert(dsl > 0.0);
			const vec3 drj = dri - normal * dsl;

			for (int k = 0; k < Fluid::NFLUID; k++)
			{
				const real dWit = Wrec_i.t[k]*dtij_list[iface].first;
				const real dWix = vec3(Wrec_i.x[k], Wrec_i.y[k], Wrec_i.z[k]) * dri;
				const real dWi  = dWix + dWit;
				min_limiter[i][k] = std::min(min_limiter[i][k], limiter(dWi,  dWj_minmax[i].first[k], dWj_minmax[i].second[k], 0.0));
				
        const real dWjt = Wrec_j.t[k]*dtij_list[iface].second;
				const real dWjx = vec3(Wrec_j.x[k], Wrec_j.y[k], Wrec_j.z[k]) * drj;
				const real dWj  = dWjx + dWjt;
				min_limiter[j][k] = std::min(min_limiter[j][k], limiter(dWj,  dWj_minmax[j].first[k], dWj_minmax[j].second[k], 0.0));
			}
#if 0
			min_limiter[i][k] = std::min(min_limiter[i][k], limiter(dWix, dWj_min[i][k], dWj_max[i][k], ptcl_import[i].volume));
			min_limiter[j][k] = std::min(min_limiter[j][k], limiter(dWjx, dWj_min[j][k], dWj_max[j][k], ptcl_import[j].volume));
			min_limiter[i][k] = std::min(min_limiter[i][k], limiter(dWit, dWj_min[i][k], dWj_max[i][k], ptcl_import[i].volume));
			min_limiter[j][k] = std::min(min_limiter[j][k], limiter(dWjt, dWj_min[j][k], dWj_max[j][k], ptcl_import[j].volume));
#endif
		}


#if 1
		{
			std::vector<Fluid> &min_limiter_local = dWj_min_local;

			const int nimport = site_import.size();

			for (int i = 0; i < nimport; i++)
			{
				const SiteImport &s = site_import[i];
				const int proc = s.proc();
				const int id   = s.id();
				ptcl_import[i].local_id = id;
				assert(proc < nproc);

				//      if (!(site_import[i].is_active() || site_import[i].is_local_ngb()))  continue;
				//        if (!(s.is_active() || s.is_local_ngb()))  continue;
				if (!(s.is_local_ngb() || s.is_hasngb())) continue;
				assert(site_import[i].is_active() || site_import[i].is_local_ngb());

				const int  local_id = id;
				const int remote_id = i;

				if (proc != myproc)
				{
					idx_send  [proc].push_back(local_id);
					idx_send  [proc].push_back(remote_id);
					fluid_send[proc].push_back(min_limiter[i]);
				}
				else
					for (int k = 0; k < Fluid::NFLUID; k++)
						min_limiter_local[local_id][k] = 1.0;
			}

			myMPI::all2all<false>(  idx_send,   idx_recv, myproc, nproc, 2, nsend_out, nrecv_out);
			myMPI::all2all<false>(fluid_send, fluid_recv, myproc, nproc, 1, nsend_out, nrecv_out);

			for (int p = 0; p < nproc; p++)
				for (int q = 0; q < (const int)idx_recv[p].size(); q += 2)
				{
					assert(p!=myproc);
					const int local_id = idx_recv[p][q + 0];
					assert(local_id >= 0.0);
					assert(local_id < (int)local_n);
					for (int k = 0; k < Fluid::NFLUID; k++)
						min_limiter_local[local_id][k] = 1.0;
				}


			for (int p = 0; p < nproc; p++)
			{
				for (int q = 0; q < (const int)idx_recv[p].size(); q += 2)
				{
					const int local_id = idx_recv[p][q + 0];

					const Fluid &min_limit = fluid_recv[p][q >> 1];

					for (int k = 0; k < Fluid::NFLUID; k++)
						min_limiter_local[local_id][k] = std::min(min_limiter_local[local_id][k], min_limit[k]);
				}
				idx_send   [p].clear();
				fluid_send [p].clear();
			}

			for (int i = 0; i < nimport; i++)
			{
				const SiteImport &s = site_import[i];
				if (!(s.is_local_ngb() || s.is_hasngb())) continue;

				const int proc = s.proc();
				const int id   = s.id();

				const int  local_id = id;
				const int remote_id = i;

				if (proc != myproc) continue;
				const Fluid &min_limit = min_limiter[remote_id];
				for (int k = 0; k < Fluid::NFLUID; k++)
				{
					min_limiter_local[local_id][k] = std::min(min_limiter_local[ local_id][k], min_limit[k]);
					min_limiter[remote_id][k] = std::min(min_limiter[remote_id][k], min_limiter_local[local_id][k]);
				}
			}

			for (int p = 0; p < nproc; p++)
				for (int q = 0; q < (const int)idx_recv[p].size(); q += 2)
				{
					assert (p != myproc);
					const int local_id  = idx_recv[p][q + 0];
					const int remote_id = idx_recv[p][q + 1];
					assert(local_id >= 0);
					assert(local_id < (int)local_n);

					idx_send  [p].push_back(local_id);
					idx_send  [p].push_back(remote_id);
					fluid_send[p].push_back(min_limiter_local[local_id]);
				}

			myMPI::all2all<false>(  idx_send,   idx_recv, myproc, nproc, 2, nsend_in, nrecv_in);
			myMPI::all2all<false>(fluid_send, fluid_recv, myproc, nproc, 1, nsend_in, nrecv_in);

			for (int p = 0; p < nproc; p++)
				for (int q = 0; q < (const int)idx_recv[p].size(); q += 2)
				{
					assert(p!=myproc);
					const int local_id  = idx_recv[p][q + 0];
					const int remote_id = idx_recv[p][q + 1];
					assert(remote_id >= 0);
					assert(remote_id < nimport);
					assert(local_id == ptcl_import[remote_id].local_id);

					const Fluid &min_limit = fluid_recv[p][q >> 1];

					for (int k = 0; k < Fluid::NFLUID; k++)
						min_limiter[remote_id][k] = std::min(min_limiter[remote_id][k], min_limit[k]);
				}
		}
#endif 

		for (int i = 0; i < nimport; i++)
		{
			if (!(site_import[i].is_active() || site_import[i].is_local_ngb()))  continue;
			for (int k = 0; k < Fluid::NFLUID; k++) 
			{
				const real tau = min_limiter[i][k];

				assert(tau >= 0.0 && tau <= 1.0);

				Wrec_import[i].x[k] *= tau;
				Wrec_import[i].y[k] *= tau;
				Wrec_import[i].z[k] *= tau;
				Wrec_import[i].t[k] *= tau;
			}
		}


#endif // slope_limiter
	}

}
