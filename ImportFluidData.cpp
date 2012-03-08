#include "fvmhd3d.h"

namespace fvmhd3d
{

#if 0
  void System::ImportFluidPrimitives(CkCallback &cb)
//  {
    ImportFluidPrimitivesCb = cb;
#else
  void System::ImportFluidPrimitives()
  {
#endif

    assert(ImportFluidPrimitives_nRequested == 0);

    /********* Request results for active mesh-cells from remote chares *********/

    const int nremote     = nimport_glb;
    const int iremote_end = ptcl_act.size();
    const int iremote_beg = iremote_end - nremote;

    std::vector< std::pair<int, int> > request_list;
    request_list.reserve(nremote);
    for (int i = iremote_beg; i < iremote_end; i++)
      if (ptcl_act[i]->is_ngb())
      {
        assert(ptcl_act[i]->chare() != thisIndex);
        request_list.push_back(std::make_pair(ptcl_act[i]->chare(), i));
      }

    std::sort(request_list.begin(), request_list.end(), std_pair_first_sort());

    const int nrequest = request_list.size();
    CkVec< pair<int, int> > sites2request;
    sites2request.reserve(nrequest);
    request_list.push_back(std::make_pair(-1,-1));

    for (int i = 0; i < nrequest; i++)
    {
      const int iElement = request_list[i].first;
      const int iId      = request_list[i].second;
      sites2request.push_back(std::make_pair(ptcl_act[iId]->id(), iId));
      assert(iElement >= 0);
      assert(iElement < numElements);
      assert(iElement != thisIndex);
      if (iElement != request_list[i+1].first && sites2request.size() > 0)
      {
        ImportFluidPrimitives_nRequested++;
        systemProxy[iElement].ImportFluidPrimitives_request(sites2request, thisIndex);
        sites2request.clear();
      }
    }

    if (ImportFluidPrimitives_nRequested == 0)
      ImportFluidPrimitivesCb.send();
  }

  void System::ImportFluidPrimitives_request(const CkVec< pair<int,int> > &reqData, const int recvIndex)
  {
    assert(thisIndex != recvIndex);
    const int nrecv = reqData.size();
    assert(nrecv > 0);
    CkVec< pair<int, pair<Fluid, MeshPoint> > > data2send(nrecv);
    for (int i = 0; i < nrecv; i++)
    {
      const int  local_id = reqData[i].first;
      const int remote_id = reqData[i].second;
      assert(local_id >= 0);
      assert(local_id < local_n);
      data2send[i].first  = remote_id;
      Fluid Wi = Wrec_list[local_id].w;
#if 1
      if (!ptcl_list[local_id].is_active())
      {
        Fluid Wi = Wrec_list[local_id].w;
        const real dt = t_global - mesh_pnts[local_id].tbeg;
        for (int k = 0; k < Fluid::NFLUID; k++)
          Wi[k] += Wrec_list[local_id].t[k]*dt;
      }
#endif
      data2send[i].second = std::make_pair(Wi, mesh_pnts[local_id]);
    }

    systemProxy[recvIndex].ImportFluidPrimitives_recv(data2send);
  }

  void System::ImportFluidPrimitives_recv(const CkVec< pair<int, pair<Fluid, MeshPoint> > > &recvUpdates)
  {
    const int nrecv = recvUpdates.size();
    for (int i = 0; i < nrecv; i++)
    {
      const int   iId      = recvUpdates[i].first;
      assert(iId >= (int)(nactive_loc + nimport_loc));
      assert(iId <  (int) ptcl_act.size());
      Wrec_act [iId]->w   = recvUpdates[i].second.first;
      *mesh_act[iId]      = recvUpdates[i].second.second;

      Wrec_act[iId]->vel = mesh_act[iId]->vel;
      assert(Wrec_act[iId]->w[Fluid::DENS] > 0.0);
      assert(Wrec_act[iId]->w[Fluid::ETHM] > 0.0);
    }
    ImportFluidPrimitives_nRequested--;
    assert(ImportFluidPrimitives_nRequested >= 0);

    if (ImportFluidPrimitives_nRequested == 0)
      ImportFluidPrimitivesCb.send();
  }

  //////////////////
  //////////////////
  //////////////////

#if 0
  void System::ImportFluidData(CkCallback &cb)
//  {
    ImportFluidData_completeCb = cb;
#else
  void System::ImportFluidData()
  {
#endif
    assert(ImportFluidData_nRequestedUpdates == 0);

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
    CkVec< pair<int, int> > sites2request;
    sites2request.reserve(nrequest);
    request_list.push_back(std::make_pair(-1,-1));

    for (int i = 0; i < nrequest; i++)
    {
      const int iElement = request_list[i].first;
      const int iId      = request_list[i].second;
      sites2request.push_back(std::make_pair(ptcl_act[iId]->id(), iId));
      assert(iElement >= 0);
      assert(iElement < numElements);
      assert(iElement != thisIndex);
      if (iElement != request_list[i+1].first && sites2request.size() > 0)
      {
        ImportFluidData_nRequestedUpdates++;
        systemProxy[iElement].ImportFluidData_request(sites2request, thisIndex);
        sites2request.clear();
      }
    }


    if (ImportFluidData_nRequestedUpdates == 0)
    {
#if 1
      ImportFluidData_completeCb.send();
#else
      contribute(ImportFluidData_completeCb);
#endif
    }
  }

  void System::ImportFluidData_request(const CkVec< pair<int,int> > &reqData, const int recvIndex)
  {
    assert(thisIndex != recvIndex);
    const int nrecv = reqData.size();
    assert(nrecv > 0);
    MeshFluidGrad_msg *msg = new (nrecv, nrecv, nrecv, nrecv) MeshFluidGrad_msg;
    msg->n = nrecv;
    for (int i = 0; i < nrecv; i++)
    {
      const int  local_id = reqData[i].first;
      const int remote_id = reqData[i].second;
      assert(local_id >= 0);
      assert(local_id < local_n);
      msg->  id_list[i]   = remote_id;
      msg->mesh_pnts[i]   = mesh_pnts[local_id];
      msg->Wrec_list[i]   = Wrec_list[local_id];
    
#if 1
      if (!ptcl_list[local_id].is_active())
      {
        Fluid Wi = Wrec_list[local_id].w;
        const real dt = t_global - mesh_pnts[local_id].tbeg;
        for (int k = 0; k < Fluid::NFLUID; k++)
          Wi[k] += Wrec_list[local_id].t[k]*dt;
        msg->Wrec_list[i].w = Wi;
      }
#endif

      msg->Wrec_minmax[i] = Wrec_minmax_list[local_id];
    }

    systemProxy[recvIndex].ImportFluidData_recv(msg);
  }

  void System::ImportFluidData_recv(MeshFluidGrad_msg *msg)
  {
    const int nrecv = msg->n;
    assert(nrecv > 0);
    for (int i = 0; i < nrecv; i++)
    {
      const int   id  = msg->  id_list[i];
      assert(id >= (int)(ptcl_act.size() - nimport_glb) );
      assert(id <  (int) ptcl_act.size() );

      *mesh_act[id] = msg->mesh_pnts[i];
      *Wrec_act[id] = msg->Wrec_list[i];
      *Wrec_minmax_act[id] = msg->Wrec_minmax[i];

      Wrec_act[id]->tend = mesh_act[id]->tend;
      Wrec_act[id]->vel  = mesh_act[id]->vel;
      Wrec_act[id]->pos  = mesh_act[id]->pos;
      Wrec_act[id]->bnd  = mesh_act[id]->boundary;
      Wrec_act[id]->etaJ = mesh_act[id]->etaJ;
      Wrec_act[id]->acc  = mesh_act[id]->acc1;

      assert(Wrec_act[id]->w[Fluid::DENS] > 0.0);
      assert(Wrec_act[id]->w[Fluid::ETHM] > 0.0);
    }

    ImportFluidData_nRequestedUpdates--;
    assert(ImportFluidData_nRequestedUpdates >= 0);

    if (ImportFluidData_nRequestedUpdates == 0)
    {
#if 1
      ImportFluidData_completeCb.send();
#else
      contribute(ImportFluidData_completeCb);
#endif
    }
    delete msg;
  }

  ///////////////////
  ///////////////////
  ///////////////////
  
  void System::ImportNewDt()
  {

    assert(ImportNewDt_nRequested == 0);

    /********* Request results for active mesh-cells from remote chares *********/

    const int nremote     = nimport_glb;
    const int iremote_end = ptcl_act.size();
    const int iremote_beg = iremote_end - nremote;

    std::vector< std::pair<int, int> > request_list;
    request_list.reserve(nremote);
    for (int i = iremote_beg; i < iremote_end; i++)
      if (ptcl_act[i]->is_ngb())
        if (ptcl_act[i]->is_active())
        {
          assert(ptcl_act[i]->chare() != thisIndex);
          request_list.push_back(std::make_pair(ptcl_act[i]->chare(), i));
        }

    std::sort(request_list.begin(), request_list.end(), std_pair_first_sort());

    const int nrequest = request_list.size();
    CkVec< pair<int, int> > sites2request;
    sites2request.reserve(nrequest);
    request_list.push_back(std::make_pair(-1,-1));

    for (int i = 0; i < nrequest; i++)
    {
      const int iElement = request_list[i].first;
      const int iId      = request_list[i].second;
      sites2request.push_back(std::make_pair(ptcl_act[iId]->id(), iId));
      assert(iElement >= 0);
      assert(iElement < numElements);
      assert(iElement != thisIndex);
      if (iElement != request_list[i+1].first && sites2request.size() > 0)
      {
        ImportNewDt_nRequested++;
        systemProxy[iElement].ImportNewDt_request(sites2request, thisIndex);
        sites2request.clear();
      }
    }

    if (ImportNewDt_nRequested == 0)
      ImportNewDtCb.send();
  }

  void System::ImportNewDt_request(const CkVec< pair<int,int> > &reqData, const int recvIndex)
  {
    assert(thisIndex != recvIndex);
    const int nrecv = reqData.size();
    assert(nrecv > 0);
    CkVec< pair<int, real> > data2send;
    data2send.reserve(nrecv);
    for (int i = 0; i < nrecv; i++)
    {
      const int  local_id = reqData[i].first;
      const int remote_id = reqData[i].second;
      assert(local_id >= 0);
      assert(local_id < local_n);
      assert(ptcl_list[local_id].is_active());
      data2send.push_back(std::make_pair(remote_id, mesh_pnts[local_id].dt_new));
    }

    systemProxy[recvIndex].ImportNewDt_recv(data2send);
  }

  void System::ImportNewDt_recv(const CkVec< pair<int, real> > &recvUpdates)
  {
    const int nrecv = recvUpdates.size();
    for (int i = 0; i < nrecv; i++)
    {
      const int   iId      = recvUpdates[i].first;
      assert(iId >= (int)(nactive_loc + nimport_loc));
      assert(iId <  (int) ptcl_act.size());
      assert(ptcl_act[iId]->is_active());
      mesh_act[iId]->dt_new = recvUpdates[i].second;
    }
    ImportNewDt_nRequested--;
    assert(ImportNewDt_nRequested >= 0);

    if (ImportNewDt_nRequested == 0)
      ImportNewDtCb.send();
  }


}
