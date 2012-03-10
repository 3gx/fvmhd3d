#include "fvmhd3d.h"

// #define _SLOW_
//
namespace fvmhd3d {

#if 0
	void system::exchange_ptcl() 
	{
		assert(false);
		import_n = 0;
		if (no_export_return_flag) return;

		/// communicate ptcl data
		std::vector<Particle > ptcl_send[NMAXPROC];
		std::vector<Particle > ptcl_recv[NMAXPROC];
		for (int p = 0; p < nproc; p++) 
			ptcl_send[p].reserve(128);

		for (size_t i = 0; i < export_site_list.size(); i++) 
		{
			const std::pair<int, int> p = export_site_list[i];
			const int proc = p.first;
			const int jidx = p.second;

			assert(proc != myproc);
			assert(jidx >= 0);
			assert(jidx < local_n);

			ptcl_send[proc].push_back(ptcl[jidx]);
		}

#ifdef _SLOW_
		myMPI::all2all(ptcl_send, ptcl_recv,  myproc, nproc, mpi_debug_flag);
#else
		myMPI::all2all(ptcl_send, ptcl_recv, nsite_export, 1, myproc, nproc, mpi_debug_flag);
#endif
		
		int nrecv = 0;
		for (int p = 0; p < nproc; p++)
			nrecv += ptcl_recv[p].size();

		import_n = nrecv;

		ptcl.resize(local_n + import_n);

		int iloc = local_n;
		for (int p = 0; p < nproc; p++) 
			for (size_t q = 0; q < ptcl_recv[p].size(); q++) 
				ptcl[iloc++] = ptcl_recv[p][q];
	}
	
	void system::exchange_Wstate() 
	{
		if (no_export_return_flag) return;

		std::vector<Fluid_st> Fluid_send[NMAXPROC];
		std::vector<Fluid_st> Fluid_recv[NMAXPROC];
		for (int p = 0; p < nproc; p++) 
			Fluid_send[p].reserve(128);

		const int nexport = export_site_list.size();
		for (int i = 0; i < nexport; i++) 
		{
			const std::pair<int, int> p = export_site_list[i];
			const int proc = p.first;
			const int jidx = p.second;

			assert(jidx < local_n);
			Fluid_send[proc].push_back(Wst[jidx]);
		}

#ifdef _SLOW_
		myMPI::all2all(Fluid_send, Fluid_recv, myproc, nproc, mpi_debug_flag);
#else
		myMPI::all2all(Fluid_send, Fluid_recv, nsite_export, 1, myproc, nproc, mpi_debug_flag);
#endif

		int nrecv = 0;
		for (int p = 0; p < nproc; p++)
			nrecv += Fluid_recv[p].size();

		assert(nrecv == import_n);
		Wst.resize(local_n + import_n);

		int iloc = local_n;
		for (int p = 0; p < nproc; p++) 
			for (size_t q = 0; q < Fluid_recv[p].size(); q++)
				Wst[iloc++] = Fluid_recv[p][q];

		assert(iloc == local_n + nrecv);
	}

	void system::exchange_primitive() 
	{
		if (no_export_return_flag) return;

		std::vector<Fluid> Fluid_send[NMAXPROC];
		std::vector<Fluid> Fluid_recv[NMAXPROC];
		for (int p = 0; p < nproc; p++) 
			Fluid_send[p].reserve(128);

		const int nexport = export_site_list.size();
		for (int i = 0; i < nexport; i++) 
		{
			const std::pair<int, int> p = export_site_list[i];
			const int proc = p.first;
			const int jidx = p.second;

			assert(jidx < local_n);
			Fluid_send[proc].push_back(Wst [jidx].w);
			Fluid_send[proc].push_back(Wrec[jidx].t);
		}

#ifdef _SLOW_
		myMPI::all2all(Fluid_send, Fluid_recv, myproc, nproc, mpi_debug_flag);
#else
		myMPI::all2all(Fluid_send, Fluid_recv, nsite_export, 2, myproc, nproc, mpi_debug_flag);
#endif

		int nrecv = 0;
		for (int p = 0; p < nproc; p++)
			nrecv += Fluid_recv[p].size();

		assert(nrecv/2 == import_n);
		Wrec.resize(local_n + import_n);

		int iloc = local_n;
		for (int p = 0; p < nproc; p++) 
			for (size_t q = 0; q < Fluid_recv[p].size(); q += 2)
			{
				Wst [iloc].w = Fluid_recv[p][q + 0];
				Wrec[iloc].t = Fluid_recv[p][q + 1];
				iloc++;
			}

		assert(iloc == local_n + nrecv/2);
		assert(nrecv%2 == 0);
	}
	
	

	void system::exchange_reconstruction() 
	{
		if (no_export_return_flag) return;

		/// communicate gradients
		std::vector<Fluid_rec> Fluid_send[NMAXPROC];
		std::vector<Fluid_rec> Fluid_recv[NMAXPROC];
		for (int p = 0; p < nproc; p++) 
			Fluid_send[p].reserve(128);
		
		const int nexport = export_site_list.size();
		for (int i = 0; i < nexport; i++) 
		{
			const std::pair<int, int> p = export_site_list[i];
			const int proc = p.first;
			const int jidx = p.second;
			assert(jidx < local_n);

			Fluid_send[proc].push_back(Wrec[jidx]);
		}

#if 1
		myMPI::all2all(Fluid_send, Fluid_recv, myproc, nproc, mpi_debug_flag);
#else
		myMPI::all2all(Fluid_send, Fluid_recv, nsite_export, 1, myproc, nproc, mpi_debug_flag);
#endif

		int nrecv = 0;
		for (int p = 0; p < nproc; p++)
			nrecv += Fluid_recv[p].size();

		assert(local_n + nrecv == (int)ptcl.size());

		Wrec.resize(local_n + import_n);

		int iloc = local_n;
		for (int p = 0; p < nproc; p++) 
			for (size_t q = 0; q < Fluid_recv[p].size(); q++)
				Wrec[iloc++] = Fluid_recv[p][q];

		assert(iloc == local_n + nrecv);
		assert(nrecv == import_n);
	}


	void system::exchange_pvel() 
	{
		if (no_export_return_flag) return;

		/// communicate pvel data
		std::vector<double   > dble_send[NMAXPROC];
		std::vector<double   > dble_recv[NMAXPROC];
		for (int p = 0; p < nproc; p++) 
		{
			dble_send[p].clear();
			dble_recv[p].clear();
		}

		for (size_t i = 0; i < export_site_list.size(); i++) 
		{
			const std::pair<int, int> p = export_site_list[i];
			const int proc = p.first;
			const int jidx = p.second;
			assert(jidx < local_n);

			dble_send[proc].push_back(ptcl[jidx].vx);
			dble_send[proc].push_back(ptcl[jidx].vy);
			dble_send[proc].push_back(ptcl[jidx].vz);
		}

#ifdef _SLOW_
		myMPI::all2all(dble_send, dble_recv, myproc, nproc, mpi_debug_flag);
#else
		myMPI::all2all(dble_send, dble_recv, nsite_export, 3, myproc, nproc, mpi_debug_flag);
#endif

		int nrecv = 0;
		for (int p = 0; p < nproc; p++)
			nrecv += dble_recv[p].size();

		assert(nrecv % 3 == 0);
		assert(local_n + nrecv/3 == (int)ptcl.size());

		int iloc = local_n;
		for (int p = 0; p < nproc; p++) 
			for (size_t i = 0; i < dble_recv[p].size(); i += 3) 
			{
				ptcl[iloc].vx = dble_recv[p][i + 0];
				ptcl[iloc].vy = dble_recv[p][i + 1];
				ptcl[iloc].vz = dble_recv[p][i + 2];
				iloc++;
			}

		assert(iloc == local_n + nrecv/3);
		assert(nrecv%3 == 0);

		assert(nrecv/3 == import_n);
	}
#endif

}
