#include "fvmhd3d.h"
#include <algorithm>

namespace fvmhd3d
{
  
  void System::moveParticles(CkCallback &cb)
  {
    moveParticles_completeCb = cb;

    domains = (globalDomains*)CkLocalNodeBranch(globalDomainsProxy);
    assert(domains != NULL);

    const int n = ptcl_list.size();
    assert(n == local_n);

#if 1
    for (int i = 0; i < local_n; i++)
    {
      const vec3 pos = periodic(mesh_pnts[i].pos);
      mesh_pnts[i].set_pos(pos);
      ptcl_list[i] = pos;
      assert(global_domain.isinbox(ptcl_list[i].get_pos()));
      assert(global_domain.isinbox(mesh_pnts[i].pos));
      assert(ptcl_list[i].get_pos().norm2() == mesh_pnts[i].pos.norm2());
    }
#endif

    int iloc = 0;
    std::vector< std::pair<int, int> > to_which_elements(n);
    std::vector<int> remote_tiles;
    for (int i = 0; i < n; i++)
    {
      remote_tiles.clear();	
      domains->proc_tree.root.walk_boundary(boundary(ptcl_list[i].get_pos()), remote_tiles, global_domain_size);
      assert(remote_tiles.size() > 0);
      const int send2Element = remote_tiles[0];

      to_which_elements[i] = std::make_pair(send2Element, -1);
      if (send2Element == thisIndex)
      {
        std::swap(to_which_elements[i], to_which_elements[iloc]);
        std::swap(        ptcl_list[i],         ptcl_list[iloc]);
        std::swap(        mesh_pnts[i],         mesh_pnts[iloc]);
        std::swap(           U_list[i],            U_list[iloc]);
        std::swap(          dU_list[i],           dU_list[iloc]);
        std::swap(        Wrec_list[i],         Wrec_list[iloc]);
        std::swap( Wrec_minmax_list[i],  Wrec_minmax_list[iloc]);
        iloc++;
      }
    }
    for (int i = 0; i < iloc; i++)
      ptcl_list[i] = std::make_pair(thisIndex, i);

    for (int i = iloc; i < n; i++)
    {
      assert(to_which_elements[i].first != thisIndex);
      to_which_elements[i].second = i;
    }

    std::sort(to_which_elements.begin()+iloc, to_which_elements.end(), std_pair_first_sort());

    to_which_elements.push_back(std::make_pair(-1, -1));
    std::vector<int> how_many_to_send;
    int n2send = 0;
    for (int i = iloc; i < n; i++)
    {
      n2send++;
      if (to_which_elements[i].first != to_which_elements[i+1].first)
      {
        how_many_to_send.push_back(n2send);
        n2send = 0;
      }
      assert(to_which_elements[i].first != thisIndex);
    }


    std::vector<int> elementsList;

    assert(iloc <= local_n);
    if (iloc < local_n)
    {
      int cnt = 0; 
      int nsend = how_many_to_send[cnt++];
      assert(nsend > 0);
      MeshFluidGradComplete_msg *msg = new (nsend, nsend, nsend, nsend, nsend, nsend) MeshFluidGradComplete_msg;
      msg->n = nsend;
      for (int i = iloc; i < n; i++)
      {
        const int iElement = to_which_elements[i].first;
        const int iId      = to_which_elements[i].second;
        assert(iElement != thisIndex);
        nsend--;
        assert(nsend >= 0);
        msg->ptcl_list[nsend] = ptcl_list[iId];
        msg->mesh_pnts[nsend] = mesh_pnts[iId];
        msg->   U_list[nsend] =    U_list[iId];
        msg->  dU_list[nsend] =   dU_list[iId];
        msg->Wrec_list[nsend] = Wrec_list[iId];
        msg->Wrec_minmax[nsend] = Wrec_minmax_list[iId];
        if (to_which_elements[i+1].first != iElement)
        {
          nSend_list[iElement] += msg->n;
          elementsList.push_back(iElement);
          systemProxy[iElement].moveParticles_recv(msg);
          assert(cnt <= (int)how_many_to_send.size());
          if (cnt < (int)how_many_to_send.size())
          {
            nsend = how_many_to_send[cnt++];
            msg = new (nsend, nsend, nsend, nsend, nsend, nsend) MeshFluidGradComplete_msg;
            msg->n = nsend;
          }
          else
            assert(i + 1 == n);
        }
      }
    }

    ptcl_list.resize(iloc);
    mesh_pnts.resize(iloc);
    U_list   .resize(iloc);
    dU_list  .resize(iloc);
    Wrec_list.resize(iloc);
    Wrec_minmax_list.resize(iloc);

    //		moveParticles_flag = true;

    contribute(numElements*sizeof(int), &nSend_list[0], CkReduction::sum_int, 
        CkCallback(CkIndex_System::moveParticles_reduction(NULL), thisProxy));

    for (std::vector<int>::iterator it = elementsList.begin(); it != elementsList.end(); it++)
      nSend_list[*it] = 0;
  }

  void System::moveParticles_reduction(CkReductionMsg *msg)
  {
    nRecv_list[thisIndex] -= ((int*)msg->getData())[thisIndex];
    if (nRecv_list[thisIndex] == 0)
      contribute(CkCallback(CkIndex_System::moveParticles_complete(), thisProxy));
    delete msg;
  }

  void System::moveParticles_recv(MeshFluidGradComplete_msg *msg)
  {
    moveParticles_flag = true;
    if (ptcl_import_ptr == NULL)
    {
      assert(mesh_import_ptr == NULL);
      assert(Wrec_import_ptr == NULL);
      assert(   U_import_ptr == NULL);
      assert(  dU_import_ptr == NULL);
      assert(Wrec_minmax_import_ptr == NULL);
      ptcl_import_ptr = new std::vector<Particle >();
      mesh_import_ptr = new std::vector<MeshPoint>();
      Wrec_import_ptr = new std::vector<Fluid_rec>();
      U_import_ptr    = new std::vector<Fluid    >();
      dU_import_ptr   = new std::vector<FluidD   >();
      Wrec_minmax_import_ptr = new std::vector< pair<Fluid_flt, Fluid_flt>   >();
    }

    const int nrecv = msg->n;
    for (int i = 0; i < nrecv; i++)
    {
      ptcl_import_ptr->push_back(msg->ptcl_list[i]);
      mesh_import_ptr->push_back(msg->mesh_pnts[i]);
      U_import_ptr   ->push_back(msg->   U_list[i]);
      dU_import_ptr  ->push_back(msg->  dU_list[i]);
      Wrec_import_ptr->push_back(msg->Wrec_list[i]);
      Wrec_minmax_import_ptr->push_back(msg->Wrec_minmax[i]);
    }
    nRecv_list[thisIndex] += msg->n;
    if (nRecv_list[thisIndex] == 0)
      contribute(CkCallback(CkIndex_System::moveParticles_complete(), thisProxy));
    delete msg;
  }

  void System::moveParticles_complete()
  {
    assert(nRecv_list[thisIndex] == 0);

    const int nrecv = ptcl_import_ptr == NULL ? 0 : ptcl_import_ptr->size();
    for (int i = 0; i < nrecv; i++)
    {
      ptcl_list.push_back((*ptcl_import_ptr)[i]);
      mesh_pnts.push_back((*mesh_import_ptr)[i]);
      U_list   .push_back((*   U_import_ptr)[i]);
      dU_list  .push_back((*  dU_import_ptr)[i]);
      Wrec_list.push_back((*Wrec_import_ptr)[i]);
      Wrec_minmax_list.push_back((*Wrec_minmax_import_ptr)[i]);
    }

    local_n = ptcl_list.size();
    slopeLimiter_all.resize(local_n);
    divBi.resize(local_n);


    fit_vec(ptcl_list);
    fit_vec(mesh_pnts);
    fit_vec(U_list);
    fit_vec(dU_list);
    fit_vec(Wrec_list);
    fit_vec(Wrec_minmax_list);
    fit_vec(slopeLimiter_all);
    fit_vec(divBi);


    moveParticles_flag = false;

    if (nrecv > 0)
    {
      assert(mesh_import_ptr != NULL);
      assert(Wrec_import_ptr != NULL);
      assert(   U_import_ptr != NULL);
      assert(  dU_import_ptr != NULL);
      assert(Wrec_minmax_import_ptr != NULL);
      delete ptcl_import_ptr;
      delete mesh_import_ptr;
      delete    U_import_ptr;
      delete   dU_import_ptr;
      delete Wrec_import_ptr;
      delete Wrec_minmax_import_ptr;
    }
    else
    {
      assert(mesh_import_ptr == NULL);
      assert(Wrec_import_ptr == NULL);
      assert(   U_import_ptr == NULL);
      assert(  dU_import_ptr == NULL);
      assert(Wrec_minmax_import_ptr == NULL);
    }

    ptcl_import_ptr = NULL;
    mesh_import_ptr = NULL;
    Wrec_import_ptr = NULL;
    U_import_ptr    = NULL;
    dU_import_ptr   = NULL;
    Wrec_minmax_import_ptr = NULL;
    
    clear_vec(ptcl_act);
    clear_vec(mesh_act);
    clear_vec(dU_act);
    clear_vec(U_act);
    clear_vec(Wrec_act);
    clear_vec(Wrec_minmax_act);
    clear_vec(slopeLimiter_act);

    sort_local_data();
    
    scheduler.flush_list();
    for (int i = 0; i < local_n; i++)
    {
      ptcl_list[i] = std::make_pair(thisIndex, i);
      scheduler.push_particle(i, mesh_pnts[i].rung);
    }


    clear_vec(cell_list);
    clear_vec(active_list);

#if 0
    T.clear();
    n_in_DT = 0;
    Tvtx_list.clear();
#endif


		clear_vec(ngb_list);
    ngb_list               = std::vector< Neighbours< pair<int, int> > >(local_n);
    ngb_list_hash          = std::vector< Hash<int> >(local_n);
    ngb_list_hash_used     = std::vector<int>();

    contribute(moveParticles_completeCb);
  }

}
