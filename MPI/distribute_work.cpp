#include "fvmhd3d.h"

struct cmp_proc_key
{
	bool operator() (const std::pair<int, int> &lhs, const std::pair<int, int> &rhs){
		return lhs.first < rhs.first;
	}
};

namespace fvmhd3d 
{
#if 0
  std::vector< std::pair<int, int> > sites2send[NMAXPROC];
#endif

  std::vector<int> send2remote_fluid_data_idx;
  std::vector<int> recv2remote_fluid_data_idx;

  static int nsend                  [NMAXPROC], nrecv                 [NMAXPROC];
  static int nsend2local_active_ngb [NMAXPROC], nrecv2local_active_ngb[NMAXPROC];
  static int nsend2remote_active_ngb[NMAXPROC], nrecv2remote_active_ngb[NMAXPROC];

	void system::distribute_work(const int nbatch)
	{
		{

			/** let other proc know how many particles this proc does **/
			int nactive_loc = active_ptcl.size();
      std::vector< std::pair<int, int> > active_list;
  
      static int nactive_proc [NMAXPROC];
      MPI_Allgather(&nactive_loc, 1, MPI_INT, nactive_proc, 1, MPI_INT, MPI_COMM_WORLD);

      int nactive_glb = 0;
			int np_min_loc = global_n;
			int np_max_loc = 0;
      for (int p = 0; p < nproc; p++)
      {
        nactive_glb += nactive_proc[p];
				np_min_loc = std::min(np_min_loc, nactive_proc[p]);
				np_max_loc = std::max(np_max_loc, nactive_proc[p]);
      }

			const int np_mean     = nactive_glb/nproc;
			const int dn_max_glb  = np_max_loc - np_mean;
      const int dn_min_crit = 3072; //(int)((global_n/nproc) * 0.01);

			const bool no_work_distribute = true|| all_active || dn_max_glb < dn_min_crit;
#if 0
			MPI_Barrier(MPI_COMM_WORLD);
      for (int p = 0; p < nproc; p++)
        if (myproc == 0)
          fprintf(stderr ," p= %d : nact= %d  dn_max_glb= %d %d  work_dist= %s \n",
              p, nactive_proc[p], dn_max_glb, dn_min_crit, no_work_distribute ? "false" : "true");
			MPI_Barrier(MPI_COMM_WORLD);
#endif

      if (no_work_distribute)
        for (int ip = 0; ip < nactive_loc; ip++)
          active_list.push_back(std::make_pair(myproc, active_ptcl[ip]));
      else
			{
				std::vector< std::pair<int, int> > ptcl_list;
				ptcl_list.reserve(nproc);
				for (int p = 0; p < nproc; p++)
					ptcl_list.push_back(std::make_pair(nactive_proc[p], p));
				
				std::vector<int> active_p = active_ptcl;
				std::sort(active_p.begin(), active_p.end());

				int ni = 0;
				const int q = 3;
				const int nsend_min = 1024;
				
				std::sort(ptcl_list.begin(), ptcl_list.end(), cmp_proc_key());

				for (int i = 0; i < nproc; i++)
				{
					const int idonor = nproc - 1 - i;
					const int irecvr_beg = q*i;
					const int irecvr_end = q*(i+1);
					if (irecvr_end >= idonor)
						break;

					int nmean = ptcl_list[idonor].first;
					for (int irecvr = irecvr_beg; irecvr < irecvr_end; irecvr++)
						nmean += ptcl_list[irecvr].first;
					nmean = nmean / (q + 1);

					int nsend_tot = 0;
					for (int irecvr = irecvr_beg, j=0; irecvr < irecvr_end; irecvr++,j++)
					{
						int nrecv = nmean - ptcl_list[irecvr].first;
						if (nrecv < nsend_min) continue;
						const int nsend = nrecv;
						nsend_tot += nsend;

						if (ptcl_list[idonor].second == myproc)
						{
							for (int ip = ni; ip < ni + nsend; ip++)
								active_list.push_back(std::make_pair(ptcl_list[irecvr].second, active_p[ip]));
							ni += nsend;
						}
					}
					nactive_proc[ptcl_list[idonor].second] -= nsend_tot;
					assert(nactive_proc[ptcl_list[idonor].second] > 0);
				}

				for (int ip = ni; ip < ni + nactive_proc[myproc]; ip++)
					active_list.push_back(std::make_pair(myproc, active_p[ip]));
				ni += nactive_proc[myproc];
				assert(ni == nactive_loc);
			}


#if 1
			std::vector< std::pair<int, int> > sites2send[NMAXPROC];
#else
			for (int p = 0; p < nproc; p++)
				sites2send[p].clear();
#endif

			/** now copy surface sites to this proc **/
			static std::vector< BitOps<NMAXPROC> > used_sites;
			used_sites.resize(local_n);

			std::vector<int > idx_used;
			std::vector<Ngb*> ngb_stored;


			const int ni = active_list.size();
			for (int ip = 0; ip < ni; ip++)
			{
				const int proc = active_list[ip].first;
				const int i    = active_list[ip].second;
				assert(proc < nproc);
				assert(i >= 0);
				assert(i  < (int)local_n);
				assert(!used_sites[i].is_set(proc));
				used_sites[i].set(proc);
				sites2send[proc].push_back(std::make_pair(myproc, i | (1 << 31)));
				idx_used.push_back(i);
			}


			const int ingb_max = 2;
			for (int ingb = 0; ingb < ingb_max; ingb++)
			{
				std::vector<int> site_send[NMAXPROC];
				std::vector<int> site_recv[NMAXPROC];

				std::vector< std::pair<int, int> > active_list_new;

				const int ni = active_list.size();
				for (int ip = 0; ip < ni; ip++)
				{
					const int proc = active_list[ip].first;
					const int i    = active_list[ip].second;
					Neighbours<Ngb> &ngb = ngb_list[i];
					const int nj = ngb.size(); 
					for (int jp = 0; jp < nj; jp++)   // loop over neighbours
					{
						Ngb &ngbi = ngb[jp];

						if (ngbi.id == -1)
						{
							if (ngbi.proc == myproc) 
							{
								ngbi.id = ngbi.remote_id;
								assert(ngbi.id >= 0);
								assert(ngbi.id  < (int)local_n);
							}
							else                    
							{
								used_sites.push_back(BitOps<NMAXPROC>());
								ngbi.id = used_sites.size() - 1;
							}
						}

						const int &j = ngbi.id;

						if (used_sites[j].is_set(proc)) continue;

						ngb_stored.push_back(&ngbi);
						idx_used.push_back(ngbi.id);

						used_sites[j].set(proc);

						if (j < (int)local_n)
						{
							assert(ngb[jp].proc == myproc);
							active_list_new.push_back(std::make_pair(proc, j));
						}
						else
						{
							assert(ngb[jp].proc != myproc);
							site_send[ngb[jp].proc].push_back(proc);
							site_send[ngb[jp].proc].push_back(ngb[jp].remote_id);
						}
					}
				}

				myMPI::all2all<true>(site_send, site_recv, myproc, nproc, 1, nsend, nrecv);

				int nrecv = 0;
				for (int p = 0; p < nproc; p++)
					for (int q = 0; q < (const int)site_recv[p].size(); q++)
						nrecv++;
				assert(nrecv % 2 == 0);

				for (int p = 0; p < nproc; p++)
					for (int q = 0; q < (const int)site_recv[p].size(); q += 2)
					{
						assert(p != myproc);
						const int proc = site_recv[p][q + 0];
						const int j    = site_recv[p][q + 1];
						assert(j < (int)local_n);
						if (used_sites[j].is_set(proc)) continue;
						idx_used.push_back(j);
						used_sites[j].set(proc);
						active_list_new.push_back(std::make_pair(proc, j));
					}

				for (size_t i = 0; i < active_list_new.size(); i++)
					sites2send[active_list_new[i].first].push_back(
							std::make_pair(myproc, active_list_new[i].second | (ingb == 0) << 30));

				active_list.swap(active_list_new);
			}

			// clear used_sites  & ngb;

			for (size_t i = 0; i < idx_used.size(); i++)
				used_sites[idx_used[i]].clear();

			for (size_t i = 0; i < ngb_stored.size(); i++)
				ngb_stored[i]->id = -1;



			/*** now send particles to other proc ***/
			{
				site_import.clear();

				std::vector<double> site_send[NMAXPROC];
				std::vector<double> site_recv[NMAXPROC];

				int nsend_loc = 0;
				int nsend_rem = 0;
				for (int p = 0; p < nproc; p++)
				{
					site_send[p].reserve(sites2send[p].size());
					for (int q = 0; q < (const int)sites2send[p].size(); q++)
					{
						const int proc = sites2send[p][q].first;
						assert(proc == myproc);
						const int i0   = sites2send[p][q].second;
						const int i    = i0 & 0x3FFFFFFF;
						assert(i >= 0);
						assert(i < (int)local_n);

						const int   id[]  = {proc, i0};

						predict(i);

						const vec3   pos  = ptcl_local[i].pos;
						const double rmax = ptcl_local[i].rmax;
						if (p == myproc)
						{
							nsend_loc++;
							assert(proc == p);
							site_import.push_back(SiteImport(pos, id[0], id[1], rmax));
							assert(site_import.back().proc() == myproc);
							assert(site_import.back().id() < (int)local_n);
						}
						else
						{
							nsend_rem++;
							const double *dbl_id = (double*)id;
							site_send[p].push_back(pos.x);
							site_send[p].push_back(pos.y);
							site_send[p].push_back(pos.z);
							site_send[p].push_back(dbl_id[0]);
							site_send[p].push_back(rmax);

#if 1
							{	
								//sanity_check
								const int *lid = (int*)&dbl_id[0];
								const int proc  = lid[0];
								const int id    = lid[1];
								const SiteImport s(pos, proc, id, rmax);
								assert(s.proc() == myproc);
								assert(s.id() < (int)local_n);
							}
#endif
						}
					}
				}

				myMPI::all2all<true>(site_send,  site_recv,  myproc, nproc, 5, nsend, nrecv);

				int nrecv = 0;	
				for (int p = 0; p < nproc; p++)
					for (int q = 0; q < (const int)site_recv[p].size(); q++)
						nrecv++;
				assert(nrecv % 5 == 0);

				int nr = 0;
				for (int p = 0; p < nproc; p++)
					for (int q = 0; q < (const int)site_recv[p].size(); q += 5)
					{
						nr++;
						assert(p != myproc);
						const real x = site_recv[p][q + 0];
						const real y = site_recv[p][q + 1];
						const real z = site_recv[p][q + 2];
						const int *lid = (int*)&site_recv[p][q + 3];
						const real rmax = site_recv[p][q + 4];
						const int proc  = lid[0];
						const int id    = lid[1];
						site_import.push_back(SiteImport(vec3(x,y,z), proc, id, rmax));
						assert(proc == p);
					}
			}
		}

#if 0
		{
			MPI_Barrier(MPI_COMM_WORLD);
			int nact = 0;
			for (int i = 0; i < (const int)site_import.size(); i++)
				if (site_import[i].is_active()) 
					nact++;
			fprintf(stderr, " proc= %d:  act= %d\n", myproc, nact);
			MPI_Barrier(MPI_COMM_WORLD);
		}
#endif

		{
			/*** Now build a tree of imported particles ***/
			const int nimport = site_import.size();
			import_n = nimport;


			import_tree.clear();
			import_tree.set_domain(
					boundary(
						global_domain.centre() - global_domain.hsize()*1.5,
						global_domain.centre() + global_domain.hsize()*1.5));
			std::vector<Octree::Particle> tree_ptcl(nimport);
			for (int i = 0; i < nimport; i++)
				tree_ptcl[i] = Octree::Particle(site_import[i].pos(), i);
			import_tree.insert(&tree_ptcl[0], nimport, 0, std::max(2*nimport, (int)local_n));
			import_tree.root.inner_boundary();
			import_tree.root.sanity_check();
			import_tree.get_leaves();


#if 0
			fprintf(stderr, " myproc= %d : import_tree_build= %g sec \n", myproc, t20 - t10);
#endif
		}

	}

	void system::distribute_particle_fluid_data()
	{
		const int nimport = site_import.size();

		assert(ptcl_import.empty());
		assert(U_import.empty());
		assert(dU_import.empty());
		assert(Wrec_import.empty());

		ptcl_import.resize(nimport);
		U_import.resize(nimport);
		dU_import.resize(nimport);
		Wrec_import.resize(nimport);

		for (int i = 0; i < nimport; i++)
			ptcl_import[i].boundary = -1289;


		std::vector< ParticleFluidStruct > fluid_send[NMAXPROC];
		std::vector< ParticleFluidStruct > fluid_recv[NMAXPROC];

		std::vector<int> idx_send[NMAXPROC];
		std::vector<int> idx_recv[NMAXPROC];

		send2remote_fluid_data_idx.clear();
		recv2remote_fluid_data_idx.clear();
		
		int nsend_cnt = 0;
		int nmyid = 0;
		send2remote_fluid_data_idx.push_back(0);
		recv2remote_fluid_data_idx.push_back(0);
		for (int i = 0; i < nimport; i++)
			if (site_import[i].is_active() || site_import[i].is_local_ngb())
			{
				const int      proc = site_import[i].proc();
				const int  local_id = site_import[i].id();
				const int remote_id = i;
				if (proc != myproc)
				{
					idx_send[proc].push_back( local_id);
					idx_send[proc].push_back(remote_id);
				}
				else
				{
					send2remote_fluid_data_idx.push_back(local_id);
					recv2remote_fluid_data_idx.push_back(remote_id);
					ptcl_import[remote_id] = ptcl_local[local_id];
					U_import   [remote_id] =    U_local[local_id];
					dU_import  [remote_id] =   dU_local[local_id];
					Wrec_import[remote_id] = Wrec_local[local_id];
					Wrec_import[remote_id].bnd      = remote_id;
					ptcl_import[remote_id].local_id = local_id;
					nmyid++;
				}
				nsend_cnt++;
			}
		send2remote_fluid_data_idx[0] = nmyid;
		recv2remote_fluid_data_idx[0] = nmyid;

		myMPI::all2all<true>(idx_send, idx_recv, myproc, nproc, 2, nsend2local_active_ngb, nrecv2local_active_ngb);

		for (int p = 0; p < nproc; p++)
		{
			if (p == myproc) continue; 
			fluid_send[p].reserve(idx_recv[p].size());
			send2remote_fluid_data_idx.push_back(idx_recv[p].size() >> 1);
			for (int q = 0; q < (const int)(idx_recv[p].size()); q += 2)
			{
				const int  local_id = idx_recv[p][q + 0];
				const int remote_id = idx_recv[p][q + 1];

				send2remote_fluid_data_idx.push_back(local_id);

				fluid_send[p].push_back(
						ParticleFluidStruct(
							ptcl_local[local_id], U_local[local_id], dU_local[local_id], Wrec_local[local_id]));

				fluid_send[p].back().p.local_id = local_id;
				fluid_send[p].back().Wrec.bnd   = remote_id;
			}
		}

		myMPI::all2all<true>(fluid_send, fluid_recv, myproc, nproc, 1, nsend2remote_active_ngb, nrecv2remote_active_ngb);

		int nrecv_cnt = recv2remote_fluid_data_idx[0];
		for (int p = 0; p < nproc; p++)
		{
			if (p == myproc) continue;
			const int nq = fluid_recv[p].size();
			recv2remote_fluid_data_idx.push_back(nq);
			for (int q = 0; q < nq; q++)
			{
				nrecv_cnt++;
				const int  local_id = fluid_recv[p][q].p.local_id;
				const int remote_id = fluid_recv[p][q].Wrec.bnd;
				assert(remote_id >= 0);
				assert(remote_id <  nimport);
				assert(site_import[remote_id].id()   == local_id);
				assert(site_import[remote_id].proc() == p);
				recv2remote_fluid_data_idx.push_back(remote_id);

				assert(ptcl_import[remote_id].boundary == -1289);

				ptcl_import[remote_id] = fluid_recv[p][q].p;
				U_import   [remote_id] = fluid_recv[p][q].U;
				dU_import  [remote_id] = fluid_recv[p][q].dU;
				Wrec_import[remote_id] = fluid_recv[p][q].Wrec;
			}
		}


		for (int i = 0; i < nimport; i++)
			if (site_import[i].is_active() || site_import[i].is_local_ngb())
			{
				const int remote_id = i;

				Wrec_import[remote_id].tlast = ptcl_import[remote_id].tlast; 
				Wrec_import[remote_id].vel   = ptcl_import[remote_id].vel;
				Wrec_import[remote_id].pos   = ptcl_import[remote_id].pos;
				Wrec_import[remote_id].bnd   = ptcl_import[remote_id].boundary;

				assert(site_import[remote_id].is_active() || site_import[remote_id].is_local_ngb());
				assert(ptcl_import[remote_id].boundary != -1289);

				assert(Wrec_import[remote_id].w[Fluid::DENS] > 0.0);
				assert(Wrec_import[remote_id].w[Fluid::ETHM] > 0.0);
			}

		if (!(nsend_cnt == nrecv_cnt))
			fprintf(stderr, " myproc= %d :: nsend_cnt= %d  nrecv_cnt= %d \n", myproc, nsend_cnt, nrecv_cnt);
		assert(nsend_cnt == nrecv_cnt);

		for (int i = 0; i < nimport; i++)
			if (site_import[i].is_active() || site_import[i].is_local_ngb())
				assert(ptcl_import[i].boundary != -1289);
			else
				assert(ptcl_import[i].boundary == -1289);
	}

	void system::get_site_list_batch(
			std::vector<             int  > &active_sites,
			std::vector< std::vector<int> > &active_list,
			const int batch)
	{
		assert(batch == 0);

		const int nleaves = import_tree.get_leaves();
		int np = 0;
		for (int leaf = 0; leaf < nleaves; leaf++)
		{
			std::vector<int> list;
			for (Octree::Body *bp = import_tree.leaf(leaf)->pfirst; bp != NULL; bp = bp->next)
			{

				np++;
				if (site_import[bp->id].is_active())
				{
					active_sites.push_back(bp->id);
					list.push_back(bp->id);
				}
			}

			assert(list.size() <= NLEAF_LOC);
			if (!list.empty())
				active_list.push_back(list);
		}

		assert(np == (int)site_import.size());
	}

	void system::collect_results()
	{
		std::vector<ParticleFluidStructLite> fluid_send[NMAXPROC];
		std::vector<ParticleFluidStructLite> fluid_recv[NMAXPROC];

		std::vector<int> ngb_send[NMAXPROC];
		std::vector<int> ngb_recv[NMAXPROC];

		const int nimport = site_import.size(); 

		/******* send particle data back ***/
		for (int i = 0; i < nimport; i++)
		{
			const SiteImport &s = site_import[i];
			const int proc = s.proc();
			const int id   = s.id();
			assert(proc < nproc);
			ptcl_import[i].local_id = id;
			ptcl_import[i].rmax     = s.rmax();
			if (s.is_active()) ptcl_import[i].  set_active();
			else               ptcl_import[i].unset_active();

			if (!(s.is_active() || s.is_local_ngb()))
				continue;

			if (proc == myproc) continue;

			fluid_send[proc].push_back(ParticleFluidStructLite(ptcl_import[i], U_import[i], dU_import[i]));


			const Cell &ci = cell_list[i];
			if (ci.Volume <= 0.0)
			{
				ngb_send[proc].push_back(0);
				continue;
			}

			const int nface = ci.faces().size();	
			ngb_send[proc].push_back(nface);
			for (int iface = 0; iface < nface; iface++)
			{
				const Face &face = face_list[ci.faces()[iface]];
				const int j = face.ngb<false>(i);
				assert(j >= 0);
				assert(j < nimport);
				ngb_send[proc].push_back(site_import[j].proc());
				ngb_send[proc].push_back(site_import[j].id());
			}
		}

		myMPI::all2all<false>(fluid_send, fluid_recv, myproc, nproc, 1, nsend2local_active_ngb, nrecv2local_active_ngb);
		myMPI::all2all<true >(ngb_send,   ngb_recv,   myproc, nproc, 1, nsend, nrecv);

		std::vector<int> id_recv;
		int nact = 0;
		for (int p = 0; p < nproc; p++)
			for (int q = 0; q < (const int)fluid_recv[p].size(); q++)
			{
				const Particle &pi = fluid_recv[p][q].p;
				assert(pi.local_id >= 0);
				assert(pi.local_id < (int)local_n);

				id_recv.push_back(pi.local_id);

				if (pi.is_active())
				{
					nact++;

					U_local[pi.local_id] = fluid_recv[p][q].U;
					divBi  [pi.local_id] = fluid_recv[p][q].dU.divB;
					dU_local[pi.local_id].U       = 0.0;
					dU_local[pi.local_id].divB    = 0.0;
					dU_local[pi.local_id].gradPsi = 0.0;

					ptcl_local[pi.local_id] = pi;
					correct(pi.local_id);
				}
				else if (!ptcl_local[pi.local_id].is_active() && !(ptcl_local[pi.local_id].tlast == t_global))
				{
					if (pi.tend - pi.tlast < ptcl_local[pi.local_id].tend - ptcl_local[pi.local_id].tlast)
					{
						const int rung = pi.rung < 0 ? -1-pi.rung : pi.rung;
						scheduler.move<true>(ptcl_local[pi.local_id].rung, rung, pi.local_id);
						ptcl_local[pi.local_id].rung = rung;
						ptcl_local[pi.local_id].tend = pi.tend;
						if (pi.rung < 0)
							Wrec_local[pi.local_id] = Fluid_rec(Wrec_local[pi.local_id].w);
					}

					if (!pi.is_boundary())
					{
						for (int k = 0; k < Fluid::NFLUID; k++)
							dU_local[pi.local_id].U[k]  += fluid_recv[p][q].dU.U[k];
						dU_local[pi.local_id].divB    += fluid_recv[p][q].dU.divB;
						dU_local[pi.local_id].gradPsi += fluid_recv[p][q].dU.gradPsi;
					}
				}
			}


		int irecv = 0;
		for (int p = 0; p < nproc; p++)
		{
			int q = 0;
			const int nq = ngb_recv[p].size();
			while (q < nq)
			{
				const int id = id_recv[irecv++];
				const int nj = ngb_recv[p][q++];
				if (nj == 0) continue;
				assert(nj != 0);
				assert(id >= 0);
				assert(id < (int)local_n);
				ngb_list[id].clear();
				for (int j = 0; j < nj; j++, q+=2)
					ngb_list[id].push_back(Ngb(ngb_recv[p][q], ngb_recv[p][q+1]));
				assert(q <= nq);
			}
		}

		for (int i = 0; i < nimport; i++)
		{
			const SiteImport &s = site_import[i];
			const int proc = s.proc();
			const int id   = s.id();
			assert(proc < nproc);
			ptcl_import[i].local_id = id;
			ptcl_import[i].rmax     = s.rmax();
			if (s.is_active()) ptcl_import[i].  set_active();
			else               ptcl_import[i].unset_active();

			if (!(s.is_active() || s.is_local_ngb()))
				continue;

			if (proc != myproc) continue;

			const int remote_id = i;
			const Particle &pi = ptcl_import[remote_id];

			if (pi.is_active())
			{
				nact++;

				U_local [pi.local_id] =  U_import[remote_id];
				divBi   [pi.local_id] = dU_import[remote_id].divB;
				dU_local[pi.local_id].U       = 0.0;
				dU_local[pi.local_id].divB    = 0.0;
				dU_local[pi.local_id].gradPsi = 0.0;

				ptcl_local[pi.local_id] = pi;
				correct(pi.local_id);

				const Cell &ci = cell_list[i];
				assert(ci.Volume > 0);
				const int nface = ci.faces().size();	
				ngb_list[pi.local_id].clear();
				for (int iface = 0; iface < nface; iface++)
				{
					const Face &face = face_list[ci.faces()[iface]];
					const int j = face.ngb<false>(i);
					assert(j >= 0);
					assert(j < nimport);
					ngb_list[pi.local_id].push_back(Ngb(site_import[j].proc(), site_import[j].id()));
				}
			}
			else if (!ptcl_local[pi.local_id].is_active() && !(ptcl_local[pi.local_id].tlast == t_global))
			{
				if (pi.tend - pi.tlast < ptcl_local[pi.local_id].tend - ptcl_local[pi.local_id].tlast)
				{
					const int rung = pi.rung < 0 ? -1-pi.rung : pi.rung;
					scheduler.move<true>(ptcl_local[pi.local_id].rung, rung, pi.local_id);
					ptcl_local[pi.local_id].rung = rung;
					ptcl_local[pi.local_id].tend = pi.tend;
					if (pi.rung < 0)
						Wrec_local[pi.local_id] = Fluid_rec(Wrec_local[pi.local_id].w);
				}

				if (!pi.is_boundary())
				{
					for (int k = 0; k < Fluid::NFLUID; k++)
						dU_local[pi.local_id].U[k]  += dU_import[remote_id].U[k];
					dU_local[pi.local_id].divB    += dU_import[remote_id].divB;
					dU_local[pi.local_id].gradPsi += dU_import[remote_id].gradPsi;
				}
			}
		}
		assert(nact == (int)active_ptcl.size());
	};

	void system::distribute_fluid_update_data()
	{
		std::vector<Fluid_st> fluid_send[NMAXPROC];
		std::vector<Fluid_st> fluid_recv[NMAXPROC];

		const int nimport = site_import.size();
		Wst_import.resize(nimport);

		const int ni = send2remote_fluid_data_idx.size();
		const int nj = send2remote_fluid_data_idx[0];
		int i = 1;
		for (int j = i; j < i+nj; j++)
		{
			const int  local_id = send2remote_fluid_data_idx[j];
			const int remote_id = recv2remote_fluid_data_idx[j];
			Fluid_st Wst;
			Wst.w   = U_local[local_id].to_primitive(ptcl_local[local_id].volume);
			Wst.bnd = ptcl_local[local_id].boundary;
			Wst.pos = ptcl_local[local_id].pos;
			Wst.vel = ptcl_local[local_id].vel;
			Wst.dt  = t_global - ptcl_local[local_id].tlast;
			Wst_import[remote_id] = Wst;
		}
		i += nj;

		int p = 0;
		while(i < ni)
		{
			if (p == myproc) p++;
			const int nj = send2remote_fluid_data_idx[i++];
			for (int j = i; j < i+nj; j++)
			{
				const int local_id = send2remote_fluid_data_idx[j];
				Fluid_st Wst;
				Wst.w   = U_local[local_id].to_primitive(ptcl_local[local_id].volume);
				Wst.bnd = ptcl_local[local_id].boundary;
				Wst.pos = ptcl_local[local_id].pos;
				Wst.vel = ptcl_local[local_id].vel;
				Wst.dt  = t_global - ptcl_local[local_id].tlast;
				Wst.bnd = local_id;
				fluid_send[p].push_back(Wst);
			}
			i += nj;
			p++;
			assert (i <= ni);
			assert (p <= nproc);
		}

		myMPI::all2all<false>(fluid_send, fluid_recv, myproc, nproc, 1, nsend2remote_active_ngb, nrecv2remote_active_ngb);

		i = recv2remote_fluid_data_idx[0] + 1;
		for (int p = 0; p < nproc; p++)
		{
			if (p == myproc) continue;
			const int nq = fluid_recv[p].size();
			assert(nq == recv2remote_fluid_data_idx[i++]);
			for (int q = 0; q < nq; q++)
			{
				const int remote_id = recv2remote_fluid_data_idx[i + q];
				assert(remote_id >= 0);
				assert(remote_id <  nimport);
				Wst_import[remote_id] = fluid_recv[p][q];

				assert(site_import[remote_id].id()   == Wst_import[remote_id].bnd);
				assert(site_import[remote_id].proc() == p);

				Wst_import[remote_id].bnd = ptcl_import[remote_id].boundary;
			}
			i += nq;
			assert(i <= (const int)recv2remote_fluid_data_idx.size());
		}
	}

	void system::collect_reconstruction()
	{
		std::vector<Fluid_rec> fluid_send[NMAXPROC];
		std::vector<Fluid_rec> fluid_recv[NMAXPROC];

		const int nimport = site_import.size(); 

		/******* send particle data back ***/
		int nsum = 0;
		for (int i = 0; i < nimport; i++)
		{
			const SiteImport &s = site_import[i];
			const int proc = s.proc();
			const int id   = s.id();
			assert(proc < nproc);
#if 0
			if (!s.is_active()) continue;
#else
			if (!(s.is_active() || s.is_local_ngb())) continue;
#endif
			
			Wrec_import[i].bnd = s.is_active() ? id : -1-id;

			if (proc == myproc)
			{
				if (!ptcl_local[id].is_active()) continue;
				nsum++;
				assert(id >= 0);
				assert(id < (int)local_n);
				assert(ptcl_local[id].is_active());
				Wrec_local[id] = Wrec_import[i];
				if (Wrec_local[id].tlast < 0.0)
				{
					ptcl_local[id].rung = std::max(ptcl_local[id].rung, scheduler.get_rung(-Wrec_local[id].tlast));
					ptcl_local[id].tend = t_global + scheduler.get_dt(ptcl_local[id].rung);
				}
			}
			else
				fluid_send[proc].push_back(Wrec_import[i]);
		}

#if 0
		myMPI::all2all<false>(fluid_send, fluid_recv, myproc, nproc, 1, nsend, nrecv);
#else
		myMPI::all2all<false>(fluid_send, fluid_recv, myproc, nproc, 1, nsend2local_active_ngb, nrecv2local_active_ngb);
#endif

		for (int p = 0; p < nproc; p++)
			for (int q = 0; q < (const int)fluid_recv[p].size(); q++)
			{
				const Fluid_rec &Wrec = fluid_recv[p][q];
				if (Wrec.bnd < 0) continue;
				nsum++;
				const int id = Wrec.bnd;
				assert(id >= 0);
				assert(id < (int)local_n);
				assert(ptcl_local[id].is_active());
				Wrec_local[id] = Wrec;
				if (Wrec.tlast < 0.0)
				{
					ptcl_local[id].rung = std::max(ptcl_local[id].rung, scheduler.get_rung(-Wrec_local[id].tlast));
					ptcl_local[id].tend = t_global + scheduler.get_dt(ptcl_local[id].rung);
				}
			}
		assert(nsum == (int)active_ptcl.size());
	}


}

