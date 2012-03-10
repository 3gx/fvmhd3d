#include "fvmhd3d.h"
#include "distributenew.h"

namespace fvmhd3d 
{

#define NMAXSAMPLE 2000000
	

	void system::determine_sampling_freq() 
	{
		const int nsample = std::min(NMAXSAMPLE, int((0.1*global_n)/nproc));
		sample_freq = (global_n + nsample - 1)/nsample;
	}

	void system::collect_sample_data(
			std::vector<vec3> &sample_pos,
			const std::vector<vec3> &ptcl_pos)
	{
		sample_pos.reserve(NMAXSAMPLE + local_n);
		sample_pos.clear();

		for (int i = 0; i < (int)local_n; i += sample_freq) 
			if (!ptcl_local[i].is_remove())
				sample_pos.push_back(ptcl_pos[i]);

		MPI_Status status;
		int nsample = sample_pos.size();

    assert(sizeof(vec3) == 3*sizeof(double));
		if (myproc != 0) 
		{
			MPI_Send(&nsample, 1, MPI_INT, 0, myproc * 2, MPI_COMM_WORLD);
			MPI_Send(&sample_pos[0], nsample*3, MPI_DOUBLE, 0, myproc*2 + 1, MPI_COMM_WORLD);
		} 
		else 
			for (int p = 1; p < nproc; p++) 
			{
				int nreceive;
				MPI_Recv(&nreceive, 1, MPI_INT, p, p*2, MPI_COMM_WORLD, &status);
				sample_pos.resize(nsample + nreceive);
				MPI_Recv(&sample_pos[nsample], 3*nreceive, MPI_DOUBLE, p, p*2 + 1, MPI_COMM_WORLD, &status);
				nsample += nreceive;
			}
	}

	void system::distribute_data(const bool FLUID, const bool GRADS, const bool NGB)
	{
		distribute_data_flag = true;

		ptcl_local.resize(local_n);
		std::vector<vec3> ptcl_pos(local_n);

		// compute integer coordinates for each position
		//
		for (int i = 0; i < (int)local_n; i++) 
		{
			ptcl_local[i].local_id = i;
			ptcl_local[i].orig_pos = periodic(ptcl_local[i].orig_pos);
			ptcl_pos[i] = ptcl_local[i].orig_pos;
		}


		// determine domain decomposition
		//
		std::vector<vec3> sample_pos;
		determine_sampling_freq();
		collect_sample_data(sample_pos, ptcl_pos);
    DistributeNew<real, vec3, boundary> distribute(nproc, global_domain);
		if (myproc == 0)
    {
      //			distribute_glb.determine_division(sample_pos);
#if 1
      distribute.determine_division(sample_pos, nproc*32);
#else
      distribute.determine_division(sample_pos, nproc*8);
#endif
    }

    myMPI::Bcast(distribute.tiles, 0, nproc);
    myMPI::Bcast(distribute.procs, 0, nproc);

    compute_proc_domain(distribute.tiles, distribute.procs);

    if (FLUID && GRADS)  
      for (int i = 0; i < (int)local_n; i++)
        Wrec_local[i].pos.x = divBi[i];

    int iloc = 0;

    std::vector<Particle> ptcl_send[NMAXPROC];
    std::vector<Particle> ptcl_recv[NMAXPROC];
    std::vector<ParticleFluidStruct> fluid_send[NMAXPROC];
    std::vector<ParticleFluidStruct> fluid_recv[NMAXPROC];
    std::vector<ParticleFluidStructLite> fluidlite_send[NMAXPROC];
    std::vector<ParticleFluidStructLite> fluidlite_recv[NMAXPROC];

#if 0
    std::vector<int> ngb_send[NMAXPROC];
    std::vector<int> ngb_recv[NMAXPROC];
#endif
    
    std::vector<int> remote_tiles;
    int nremove = 0;
    for (int i = 0; i < (int)local_n; i++)
    {
      remote_tiles.clear();	
      proc_tree.root.walk_boundary(boundary(ptcl_pos[i]), remote_tiles, global_domain_size);
      assert(remote_tiles.size() > 0);
      const int proc  = proc_procs[remote_tiles[0]];
      
      assert(proc >= 0);
      assert(proc < nproc);

      if (proc == myproc && !ptcl_local[i].is_remove())
      {
        std::swap(ptcl_local[i], ptcl_local[iloc]);
        std::swap(ptcl_pos  [i], ptcl_pos  [iloc]);
        if (FLUID)
        {
          std::swap(   U_local[i],    U_local[iloc]);
          std::swap(  dU_local[i],   dU_local[iloc]);
          if (GRADS)
            std::swap(Wrec_local[i], Wrec_local[iloc]);
        }
        iloc++;
      }
      else if (!ptcl_local[i].is_remove())
      {
        if (FLUID && GRADS)
          fluid_send[proc].push_back(ParticleFluidStruct(ptcl_local[i], U_local[i], dU_local[i], Wrec_local[i]));
        else if (FLUID)
          fluidlite_send[proc].push_back(ParticleFluidStructLite(ptcl_local[i], U_local[i], dU_local[i]));
        else
          ptcl_send[proc].push_back(ptcl_local[i]);
      }
      else
        nremove++;
    }


#if 0
    if (FLUID && GRADS)	myMPI::all2all(fluid_send,     fluid_recv,     myproc, nproc, mpi_debug_flag);
    else if (FLUID    )	myMPI::all2all(fluidlite_send, fluidlite_recv, myproc, nproc, mpi_debug_flag);
    else								myMPI::all2all(ptcl_send,      ptcl_recv,      myproc, nproc, mpi_debug_flag);
#else
		{
			static int nsend[NMAXPROC], nrecv[NMAXPROC];
			if (FLUID && GRADS)	myMPI::all2all<true>(fluid_send,     fluid_recv,     myproc, nproc, 1, nsend, nrecv);
			else if (FLUID    )	myMPI::all2all<true>(fluidlite_send, fluidlite_recv, myproc, nproc, 1, nsend, nrecv);
			else								myMPI::all2all<true>(ptcl_send,      ptcl_recv,      myproc, nproc, 1, nsend, nrecv);
		}
#endif

		int nrecv = 0;
		if (FLUID && GRADS)
			for (int p = 0; p < nproc; p++)
				nrecv += fluid_recv[p].size();
		else if (FLUID)
			for (int p = 0; p < nproc; p++)
				nrecv += fluidlite_recv[p].size();
		else
			for (int p = 0; p < nproc; p++)
				nrecv += ptcl_recv[p].size();

		{
			const int nloc = iloc + nrecv;

			ptcl_local.resize(nloc);    fit_vec(ptcl_local);
			U_local   .resize(nloc);    fit_vec(U_local);
			dU_local  .resize(nloc);    fit_vec(dU_local);
			Wrec_local.resize(nloc);    fit_vec(Wrec_local);
			divBi    .resize(nloc);     fit_vec(divBi);
			Wextra_local.resize(nloc);  fit_vec(Wextra_local);
		}

		for (int p = 0; p < nproc; p++) 
			for (size_t q = 0; q < (FLUID ? (GRADS ? fluid_recv[p].size() : fluidlite_recv[p].size()) : ptcl_recv[p].size()); q++) 
			{
				assert(p != myproc);

				if (FLUID && GRADS)
				{
					ptcl_local[iloc] = fluid_recv[p][q].p;
					U_local   [iloc] = fluid_recv[p][q].U;
					dU_local  [iloc] = fluid_recv[p][q].dU;
					Wrec_local[iloc] = fluid_recv[p][q].Wrec;
				}
				else if (FLUID)
				{
					ptcl_local[iloc] = fluidlite_recv[p][q].p;
					U_local   [iloc] = fluidlite_recv[p][q].U;
					dU_local  [iloc] = fluidlite_recv[p][q].dU;
				}
				else	
					ptcl_local[iloc] = ptcl_recv[p][q];

				assert(!ptcl_local[iloc].is_remove());
				iloc++;
			}

		local_n  = iloc;

		assert(iloc = (int)ptcl_local.size());

		if (FLUID && GRADS)  
			for (int i = 0; i < (int)local_n; i++)
			{
				divBi[i] = Wrec_local[i].pos.x;
				Wrec_local[i].pos.x = ptcl_local[i].pos.x;
			}


		unsigned long long nglob, nloc = local_n;
		unsigned long long nvirt_glob;
		virtual_n = nremove;
		MPI_Allreduce(&nloc, &nglob, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&virtual_n, &nvirt_glob, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

		unsigned long long local_n_min, local_n_max, local_n_mean;
		MPI_Allreduce(&local_n, &local_n_min,  1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
		MPI_Allreduce(&local_n, &local_n_max,  1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
		MPI_Allreduce(&local_n, &local_n_mean, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

		if (myproc == 0)
		{
			fprintf(stderr, "local_n= [min: %llu max: %llu ; mean: %llu ]  global_n= %llu nglob= %llu  remove_n_glob= %llu \n",
					local_n_min, local_n_max, local_n_mean/nproc, global_n, nglob, nvirt_glob);
		}

		MPI_Barrier(MPI_COMM_WORLD);
		assert(nglob == global_n - nvirt_glob);
		global_n = nglob;
		virtual_n = 0;

		sort_local_data();

		// build local tree
		// 
		global_domain_size = global_domain.hsize() * 2.0;
		local_tree.clear();
		local_tree.set_domain(
				boundary(
					global_domain.centre() - global_domain.hsize()*1.5,
					global_domain.centre() + global_domain.hsize()*1.5));

		std::vector<Octree::Particle> tree_ptcl(local_n);
		for (int i = 0; i < (int)local_n; i++)
		{
			assert(!ptcl_local[i].is_remove());
			tree_ptcl[i] = Octree::Particle(ptcl_local[i].orig_pos, i);
		}
		local_tree.insert(&tree_ptcl[0], local_n, 0, local_n);

		local_tree.get_leaves();
		local_tree.root.inner_boundary();

		if (NGB)
		{
			if (myproc == 0)
				fprintf(stderr, "---buidling mesh---\n");
			const double t10 = mytimer::get_wtime();
			clear_mesh(false);

			const double t15 = mytimer::get_wtime();
			build_mesh_global();
			double dt_mesh = mytimer::get_wtime() - t15;

			double volume_loc = 0.0;
			{
				std::vector<TREAL> v(local_n);
				for (int i = 0; i < (int)local_n; i++)
				{
					v[i] = cell_local[i].Volume;
					ptcl_local[i].volume_new = v[i];
					ptcl_local[i].local_id = i;
				}
				std::sort(v.begin(), v.end());  // sort volumes from low to high, to avoid roundoff errors
				for (int i = 0; i < (int)local_n; i++)
					volume_loc += v[i];
			}
			extract_ngb_from_mesh();
			clear_mesh(true);
			double dt = mytimer::get_wtime() - t10;

			double volume_glob = 0.0;	
			double dt_max = 0.0;
			double dt_mesh_max = 0.0;
			MPI_Allreduce(&volume_loc, &volume_glob,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&dt, &dt_max,  1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
			MPI_Allreduce(&dt_mesh, &dt_mesh_max,  1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
			const double volume_exact = global_domain_size.x*global_domain_size.y*global_domain_size.z;
			if (myproc == 0)
			{
				fprintf(stderr, " distribute::build_mesh:[ %g (all %g )  sec ::  %g cells/s/proc/thread ]\n",
						dt_mesh_max, dt_max,
						global_n/nproc/dt_mesh_max);
				fprintf(stderr, "   computed_volume= %g  exact_volume= %g diff= %g [ %g ] \n",
						volume_glob, volume_exact, 
						volume_glob - volume_exact,	(volume_glob - volume_exact)/volume_exact);
			}
		}
	}

	void system::compute_proc_domain(const std::vector<boundary> &tiles, const std::vector<int> &procs)
	{
		const int ntiles = tiles.size();
		assert(tiles.size() == procs.size());

		proc_tiles = tiles;
		proc_procs = procs;

		proc_tree.clear();
		const int node_n = (int)(ntiles * 100.0/NLEAF_GLB);	
		proc_tree.set_domain(
				boundary(
					global_domain.centre() - global_domain.hsize()*1.5,
					global_domain.centre() + global_domain.hsize()*1.5));
		proc_tree.insert(&proc_tiles[0], ntiles, 0, node_n);
		proc_tree.root.inner_boundary();
	}


}
