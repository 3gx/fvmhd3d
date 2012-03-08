#include "fvmhd3d.h"

namespace fvmhd3d
{
  inline real Phi(const real x, const real eps)
  {
#if 1
    return 1.0;
#endif
    assert(x >= 0.0);
#if 1
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

  void System::slopeLimiter(CkCallback &cb)
  {
    assert(nSend_cntr == 0);

    slopeLimiter_completeCb = cb;

    std::vector<vec3> area_test(ptcl_act.size(), 0.0);
#if 0
    Wextra_act = std::vector<FluidExtra>(ptcl_act.size(), 0.0);
#endif
    
    const int np = ptcl_act.size();
    for (std::vector<Face*>::iterator it = face_active_list.begin(); it != face_active_list.end(); it++)
    {
      const Face &face = **it;

      int i  = face.s1;
      int j  = face.s2;

      assert(i >= 0);
      assert(i < np);
      assert(j >= 0);
      assert(j < np);

      const Fluid_rec &Wrec_i = *Wrec_act[i];
      const Fluid_rec &Wrec_j = *Wrec_act[j];


      const real dtI = Wrec_i.tend - t_global;
      const real dtJ = Wrec_j.tend - t_global;
      const real dth = 0.5*std::min(dtI, dtJ);
      assert(dtI > 0.0);
      assert(dtJ > 0.0);
      assert(dth > 0.0);

      __builtin_prefetch(*(it+1));
      __builtin_prefetch(slopeLimiter_act[i]);
      __builtin_prefetch(slopeLimiter_act[j]);
      __builtin_prefetch(&Wextra_act[i]);
      __builtin_prefetch(&Wextra_act[j]);


      const vec3 dri = face.centroid - Wrec_act[i]->pos;
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

      area_test[i] += area*normal;
      area_test[j] -= area*normal;
      
      if (Wrec_i.bnd != MeshPoint::NO_BOUNDARY || Wrec_j.bnd != MeshPoint::NO_BOUNDARY)
      {
        for (int k = 0; k < Fluid::NFLUID; k++)
          (*slopeLimiter_act[i])[k] = (*slopeLimiter_act[j])[k] = 0;
        continue;
      }
    
#if 0 
      { 
        const real idsA = area * 0.5/dsh;
        vec3 drh = normal * dsh;
        vec3 fij = dri    - drh;
        drh *= idsA;
        fij *= idsA;

        vec3 Bi   = Wrec_i.w.get_B() + Wrec_i.t.get_B()*dtI;
        vec3 Bj   = Wrec_j.w.get_B() + Wrec_j.t.get_B()*dtJ;
        if       (Wrec_i.bnd != MeshPoint::NO_BOUNDARY) Bj = Bi;
        else if  (Wrec_j.bnd != MeshPoint::NO_BOUNDARY) Bi = Bj;

        const vec3 sum  = Bj + Bi;
        const vec3 diff = Bj - Bi;
        const real dxBy = (drh.x*sum.y + fij.x*diff.y);
        const real dxBz = (drh.x*sum.z + fij.x*diff.z);
        const real dyBx = (drh.y*sum.x + fij.y*diff.x);
        const real dyBz = (drh.y*sum.z + fij.y*diff.z);
        const real dzBy = (drh.z*sum.y + fij.z*diff.y);
        const real dzBx = (drh.z*sum.x + fij.z*diff.x);

        const fvec3 dJ(dyBz - dzBy, dzBx - dxBz, dxBy - dyBx);

        Wextra_act[i].J += dJ;
        Wextra_act[j].J -= dJ;
      }
#endif






      for (int k = 0; k < Fluid::NFLUID; k++)
      {
#if 1
        const real dWit = Wrec_i.t[k]*(dtI - dth);
        const real dWjt = Wrec_j.t[k]*(dtJ - dth);
#else
        const real dWit = Wrec_i.t[k]*(Wrec_i.tend - mesh_act[i]->tbeg - dth);
        const real dWjt = Wrec_j.t[k]*(Wrec_j.tend - mesh_act[j]->tbeg - dth);
#endif
        const real dWix = vec3(Wrec_i.x[k], Wrec_i.y[k], Wrec_i.z[k]) * dri;
        const real dWi  = dWix + dWit;
        (*slopeLimiter_act[i])[k] = std::min((real)(*slopeLimiter_act[i])[k], 
            limiter(dWi, 
              std::min(Wrec_minmax_act[i]->first [k] - Wrec_i.w[k], 0.0),
              std::max(Wrec_minmax_act[i]->second[k] - Wrec_i.w[k], 0.0),
              mesh_act[i]->Volume));

        const real dWjx = vec3(Wrec_j.x[k], Wrec_j.y[k], Wrec_j.z[k]) * drj;
        const real dWj  = dWjx + dWjt;
        (*slopeLimiter_act[j])[k] = std::min((real)(*slopeLimiter_act[j])[k], 
            limiter(dWj,
              std::min(Wrec_minmax_act[j]->first [k] - Wrec_j.w[k], 0.0),
              std::max(Wrec_minmax_act[j]->second[k] - Wrec_j.w[k], 0.0),
              mesh_act[j]->Volume));
#if 0
        if (k == Fluid::ETHM)
        {
          if (i < nactive_loc)
          {
            assert(Wrec_i.w[k] + (*slopeLimiter_act[i])[k]*dWi > 0.0);
            assert(Wrec_i.w[k] + dWi > 0.0);
          }
          if (j < nactive_loc)
          {
            assert(Wrec_j.w[k] + (*slopeLimiter_act[j])[k]*dWj > 0.0);
            assert(Wrec_j.w[k] + dWj > 0.0);
          }
        }
        if (k == Fluid::DENS)
        {
          if (i < nactive_loc)
          {
            assert(Wrec_i.w[k] + (*slopeLimiter_act[i])[k]*dWi > 0.0);
            assert(Wrec_i.w[k] + dWi > 0.0);
          }
          if (j < nactive_loc)
          {
            assert(Wrec_j.w[k] + (*slopeLimiter_act[j])[k]*dWj > 0.0);
            assert(Wrec_j.w[k] + dWj > 0.0);
          }
        }
        if (k == Fluid::ETHM)
        {
          if (i >= nactive_loc)
          {
            assert(Wrec_i.w[k] > 0.0);
            assert(Wrec_i.w[k] + (*slopeLimiter_act[i])[k]*dWi > 0.0);
          }
          if (j > nactive_loc)
          {
            assert(Wrec_j.w[k] > 0.0);
            assert(Wrec_j.w[k] + (*slopeLimiter_act[j])[k]*dWj > 0.0);
          }
        }
        if (k == Fluid::DENS)
        {
          if (i > nactive_loc)
          {
            assert(Wrec_i.w[k] > 0.0);
            assert(Wrec_i.w[k] + (*slopeLimiter_act[i])[k]*dWi > 0.0);
          }
          if (j > nactive_loc)
          {
            assert(Wrec_j.w[k] > 0.0);
            if (!(Wrec_j.w[k] + (*slopeLimiter_act[j])[k]*dWj > 0.0))
            {
              fprintf(stderr, "w= %g  dWj= %g  slope= %g %g  minmax= %g %g \n",
                  Wrec_j.w[k], dWj,
                  (*slopeLimiter_act[j])[k],
                  (*slopeLimiter_act[j])[k]*dWj,
                  Wrec_minmax_act[j]->first[k], 
                  Wrec_minmax_act[j]->second[k]);
            }
            assert(Wrec_j.w[k] + (*slopeLimiter_act[j])[k]*dWj > 0.0);
          }
        }
#endif
      }
    }

    for (int i = 0; i < nactive_loc; i++)
    {
      if (mesh_act[i]->is_boundary())
        continue;

      {
        const vec3 areai = area_test[i];
        if (!(areai.abs() < 1.0e-10*std::pow(cell_list[i].Volume, 2.0/3.0)))
        {
          fprintf(stderr, "thisIndex= %d::  [%d; %d] isite= %d  active= %d, area= %g %g %g   %g   vol= %g pos= %g   a= %g \n",
              thisIndex, mesh_act[i]->idx, mesh_act[i]->boundary,
              i, ptcl_act[i]->is_active(),
              areai.x, areai.y, areai.z,
              areai.abs(),
              cell_list[i].Volume, mesh_act[i]->pos.abs(),
              std::pow(cell_list[i].Volume, 2.0/3.0));
        }
        assert(areai.abs() < 1.0e-10*std::pow(cell_list[i].Volume, 2.0/3.0));
      }
    }

    /* send limiting slopes of remote ngb cells */
    {
      const int nremote     = nimport_glb;
      const int iremote_end = ptcl_act.size();
      const int iremote_beg = iremote_end - nremote;

      std::vector< std::pair<int, int> > request_list;
      request_list.reserve(nremote);
      for (int i = iremote_beg; i < iremote_end; i++)
        if (ptcl_act[i]->is_ngb())
          request_list.push_back(std::make_pair(ptcl_act[i]->chare(), i));

      std::sort(request_list.begin(), request_list.end(), std_pair_first_sort());

      const int nrequest = request_list.size();
      CkVec< pair<int, Fluid_flt> > data2send;
      data2send.reserve(nrequest);
      request_list.push_back(std::make_pair(-1,-1));

      for (int i = 0; i < nrequest; i++)
      {
        const int iElement = request_list[i].first;
        const int iId      = request_list[i].second;
        data2send.push_back(std::make_pair(ptcl_act[iId]->id(), *slopeLimiter_act[iId]));
        assert(iElement >= 0);
        assert(iElement < numElements);
        assert(iElement != thisIndex);
        if (iElement != request_list[i+1].first && data2send.size() > 0)
        {
          nSend_cntr++;
          systemProxy[iElement].slopeLimiter_recvLimiter(data2send, thisIndex);
          data2send.clear();
        }
      }
    }

    nSend_cntr++;
    slopeLimiter_recvTicket();
  }

  void System::slopeLimiter_recvLimiter(const CkVec< pair<int, Fluid_flt> > &recvData, const int recvIndex)
  {
    assert(thisIndex != recvIndex);

    const int nrecv = recvData.size();
    for (int i = 0; i < nrecv; i++)
    {
      const int  local_id = recvData[i].first;
      const Fluid_flt  &limit = recvData[i].second;
      assert(local_id >= 0);
      assert(local_id <  local_n);
      for (int k = 0; k < Fluid::NFLUID; k++)
        slopeLimiter_all[local_id][k] = std::min(slopeLimiter_all[local_id][k], limit[k]);
    }

    systemProxy[recvIndex].slopeLimiter_recvTicket();
  }

  void System::slopeLimiter_recvTicket()
  {
    nSend_cntr--;
    assert(nSend_cntr >= 0);
    if (nSend_cntr == 0)
      contribute(CkCallback(CkIndex_System::slopeLimiter_exchange(), thisProxy));
  }

  void System::slopeLimiter_exchange()
  {
    assert(slopeLimiter_nRequestedUpdates == 0);
    /* send limiting slopes of remote ngb cells */
    const int nremote     = nimport_glb;
    const int iremote_end = ptcl_act.size();
    const int iremote_beg = iremote_end - nremote;


    std::vector< std::pair<int, int> > request_list;
    request_list.reserve(nremote);
    for (int i = iremote_beg; i < iremote_end; i++)
      if (ptcl_act[i]->is_ngb())
        request_list.push_back(std::make_pair(ptcl_act[i]->chare(), i));

    std::sort(request_list.begin(), request_list.end(), std_pair_first_sort());

    const int nrequest = request_list.size();
    CkVec< pair<int, int> > data2send;
    data2send.reserve(nrequest);
    request_list.push_back(std::make_pair(-1,-1));

    for (int i = 0; i < nrequest; i++)
    {
      const int iElement = request_list[i].first;
      const int iId      = request_list[i].second;
      data2send.push_back(std::make_pair(ptcl_act[iId]->id(), iId));
      assert(iElement >= 0);
      assert(iElement < numElements);
      assert(iElement != thisIndex);
      if (iElement != request_list[i+1].first && data2send.size() > 0)
      {
        slopeLimiter_nRequestedUpdates++;
        systemProxy[iElement].slopeLimiter_requestLimiter(data2send, thisIndex);
        data2send.clear();
      }
    }

    if (slopeLimiter_nRequestedUpdates == 0)
    {
#if 1
      slopeLimiter_completeCb.send();
#else
      contribute(slopeLimiter_completeCb);
#endif
    }
  }

  void System::slopeLimiter_requestLimiter(const CkVec< pair<int, int> > &dataRequest, const int recvIndex)
  {
    const int nrecv = dataRequest.size();
    CkVec< pair<int, Fluid_flt> > sendLimiter(nrecv);
    for (int i = 0; i < nrecv; i++)
    {
      const int  local_id = dataRequest[i].first;
      const int remote_id = dataRequest[i].second;
      assert(local_id >= 0);
      assert(local_id <  local_n);
      sendLimiter[i] = std::make_pair(remote_id, slopeLimiter_all[local_id]);
    }

    systemProxy[recvIndex].slopeLimiter_recvNewLimiter(sendLimiter);
  }

  void System::slopeLimiter_recvNewLimiter(const CkVec< pair<int, Fluid_flt> > &recvLimiter)
  {
    const int nrecv = recvLimiter.size();
    for (int i = 0; i < nrecv; i++)
    {
      const int       Id = recvLimiter[i].first;
      assert(Id >= (int)(nactive_loc + nimport_loc));
      assert(Id <  (int) ptcl_act.size());
      const Fluid_flt &limit = recvLimiter[i].second;
      for (int k = 0; k < Fluid::NFLUID; k++)
        (*slopeLimiter_act[Id])[k] = std::min((*slopeLimiter_act[Id])[k], limit[k]);
    }
    slopeLimiter_nRequestedUpdates--;
    assert(slopeLimiter_nRequestedUpdates >= 0);

    if (slopeLimiter_nRequestedUpdates == 0)
    {
#if 1
      slopeLimiter_completeCb.send();
#else
      contribute(slopeLimiter_completeCb);
#endif
    }
  }

}
