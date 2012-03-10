#include "fvmhd3d.h"

#if 1
#define PROFILE(exec, what) {{double t = mytimer::get_wtime(); exec; tprofile[what] += (mytimer::get_wtime() - t);}}
#else
#define PROFILE(exec, what) {{MPI_Barrier(MPI_COMM_WORLD); double t = mytimer::get_wtime(); exec; tprofile[what] += (mytimer::get_wtime() - t); MPI_Barrier(MPI_COMM_WORLD);}}
#endif

namespace fvmhd3d 
{
	
	void system::get_active_ptcl(const bool all)
	{
		active_ptcl.clear();

		if (!all)
		{	
			dt_global = scheduler.pull_active_list<false>(active_ptcl);
			t_prev    = t_global;
			t_global  = scheduler.get_tsys();
      
			const int nactive_ptcl = active_ptcl.size();
			active_n = nactive_ptcl;
      assert(nactive_ptcl >= 0);
      if (nactive_ptcl > (int)local_n)
      {
        fprintf(stderr, "myproc= %d :: nactive_ptcl= %d \n", myproc, nactive_ptcl);
      }
      assert(nactive_ptcl <= (int)local_n);

			for (int iptcl = 0; iptcl < nactive_ptcl; iptcl++)
			{
				const int i = active_ptcl[iptcl];
				ptcl_local[i].set_active();
        if (!(ptcl_local[i].tend == t_global))
        {
          fprintf(stderr, " nactive= %d\n", nactive_ptcl);
          fprintf(stderr, " myproc= %d: i= %d  tbeg= %g   tend= %g rung= %d  ( %g ) t_global= %g %g\n",
              myproc, i, ptcl_local[i].tlast,
              ptcl_local[i].tend, ptcl_local[i].rung, 
              ptcl_local[i].tlast + scheduler.get_dt(ptcl_local[i].rung),
              dt_global, t_global);
        }
        assert(ptcl_local[i].tend == t_global);
			}
		} 
		else
		{
      for (int i = 0; i < (int)local_n; i++)
      {
        ptcl_local[i].set_active();
        assert(ptcl_local[i].is_active());
        active_ptcl.push_back(i);
      }
    }
  }

  void system::get_active_faces()
  {
    const int nimport = site_import.size();
    face_active_list.clear();
    const int nactive = site_active_list.size();
    for (int isite = 0; isite < nactive; isite++)
    {
      const int i = site_active_list[isite];
      assert(i < nimport);
      if (ptcl_import[i].is_boundary()) continue;
      assert(site_import[i].is_active());
      const Cell &ci = cell_list[i];

      const int nface = ci.faces().size();
      for (int iface = 0; iface < nface; iface++)
      {
        const int i= ci.faces()[iface];
        assert(i >= 0);
        assert(i < (const int)face_list.size());
        Face &face = face_list[ci.faces()[iface]];
        if (face.s1 < 0) continue;
        face_active_list.push_back(&face);
        face.s1 = -1-face.s1;
      }
    }
    nface_active = face_active_list.size();

    const int nactive_face = face_active_list.size();
    for (int iface = 0; iface < nactive_face; iface++)
    {
      Face &face = *face_active_list[iface];
      assert(face.s1 < 0);
      face.s1 = -1-face.s1;
    }

  }

  void system::fluid_update(const bool flag)
  {
    if (flag)
    {
      PROFILE(slope_limiter(),                GRADIENTS);
      PROFILE(compute_update(),               UPDATE);
      PROFILE(compute_pvel(),                 MISC);
#if 0
      PROFILE(compute_timesteps(true),        TIMESTEP);
#else
      PROFILE(compute_timesteps(false),       TIMESTEP);
#endif
    }
    else
    {
      PROFILE(compute_update0(),              UPDATE);
      PROFILE(compute_reconstruction(),       RECONSTRCT);
    }
  } 

  void system::iterate_step() 
  {
    iteration++;

    PROFILE(get_active_ptcl(false), ACTIVEP);

    {   
      unsigned long long nactive_ptcl = active_ptcl.size();
      unsigned long long nvirtual_glb, boundary_glb;
      MPI_Allreduce(&nactive_ptcl, &nactive_glb, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&virtual_n, &nvirtual_glb, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&boundary_n, &boundary_glb, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
      assert(nvirtual_glb == 0);
      all_active = (nactive_glb == global_n); // - nvirtual_glb - boundary_glb);
      assert(nactive_glb > 0);
#if 0
      if (!(nactive_glb <= global_n - nvirtual_glb - boundary_glb))
      {
        fprintf(stderr, " proc= %d: nactive_glb= %d  global_n= %d  nvirt= %d  bnd= %d\n",
            (int)myproc, (int)nactive_glb, (int)global_n, (int)nvirtual_glb, (int)boundary_glb);
      }
      MPI_Barrier(MPI_COMM_WORLD);
      assert(nactive_glb <= global_n - nvirtual_glb - boundary_glb);
#endif
    }

    PROFILE(distribute_work(), DISTW);

    std::vector< std::vector<int> > active_list;
    site_active_list.clear();
    site_with_ngb_active_list.clear();

    PROFILE(get_site_list_batch(site_active_list, active_list, 0), BATCH);

		const double tmesh_beg = mytimer::get_wtime();
		PROFILE(build_mesh_active(site_active_list, active_list, site_with_ngb_active_list), MESH);

		if (all_active)
		{
			real vol = 0.0;
			for (int ip = 0; ip < (const int)site_active_list.size(); ip++)
			{
				const int i = site_active_list[ip];
				vol += cell_list[i].Volume;
			}

			double volume_glob;
			MPI_Allreduce(&vol, &volume_glob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			if (myproc == 0)
			{
				const double volume_exact = global_domain_size.x*global_domain_size.y*global_domain_size.z;
				fprintf(stderr, " done in %g sec ::  computed_volume= %g  exact_volume= %g diff= %g [ %g ] \n",
						mytimer::get_wtime() - tmesh_beg,
						volume_glob, volume_exact, 
						volume_glob - volume_exact,	(volume_glob - volume_exact)/volume_exact);
			}
		}

		PROFILE(distribute_particle_fluid_data(), MPICOM);

		PROFILE(get_active_faces(), ACTIVEF);

#if 0
		fprintf(stderr, "myproc= %d:: n= %d n_ngb= %d %d :: faces= %d faces_ngb= %d\n", myproc,
				(int)site_active_list.size(), (int)site_with_ngb_active_list.size(), (int)site_import.size(),
				(int)face_list.size(), (int)face_list.size() - nface_active);
#endif

		PROFILE(fluid_update(true), FLUID);

		PROFILE(collect_results(), COLLECT);

		PROFILE(distribute_fluid_update_data(), PVEL);

		PROFILE(fluid_update(false), FLUID);


		PROFILE(collect_reconstruction(), PREDICT);




#if 1
		{
			const double t0 = mytimer::get_wtime();
			const int nactive_loc = active_ptcl.size();
			double dt_min = HUGE;
			for (int ip = 0; ip < nactive_loc; ip++)
			{
				const int i = active_ptcl[ip];
				assert(ptcl_local[i].is_active());
				scheduler.push_particle(i, (int)ptcl_local[i].rung);
				assert(fmod(t_global, scheduler.get_dt(ptcl_local[i].rung)) == 0.0);
				ptcl_local[i].unset_active();
				dt_min = std::min(dt_min, scheduler.get_dt(ptcl_local[i].rung));
				assert(!ptcl_local[i].is_active());
			}
#if 0
			if (dt_min < HUGE) fprintf(stderr, " iterate:: myproc= %d  dt_min= %g \n", myproc, dt_min);
#endif
			tprofile[MISC] = mytimer::get_wtime() - t0;
		}
#else		
		{
			scheduler.flush_list();
			for (int i = 0; i < (int)local_n; i++)
				scheduler.push_particle(i, ptcl_local[i].rung);
		}
#endif



#if 0
		cell_list.clear();
		face_list.clear();
		ptcl_import.clear();
		U_import.clear();
		dU_import.clear();
		Wrec_import.clear();
		Wst_import.clear();
		Wextre_import.clear();
#else
		clear_vec(cell_list);
		clear_vec(face_list);
		clear_vec(ptcl_import);
		clear_vec(U_import);
		clear_vec(dU_import);
		clear_vec(Wrec_import);
		clear_vec(Wst_import);
		clear_vec(Wextra_import);
#endif
	}



	void system::iterate() 
	{
		distribute_data_flag = false;

		for (int i = 0; i < NPROFILING_TYPE; i++)
			tprofile[i] = 0;

		myMPI::data_inout = 0;
		myMPI::data_inout_time = 0;


		PROFILE(iterate_step(), ITERATE);

		{
			const float nfac = 3.3;
			static unsigned long long nactive_cum = (unsigned long long)(1.01*nfac*global_n);
			static int all_steps = 32;
			unsigned long long nactive_ptcl = active_ptcl.size();

			MPI_Allreduce(&nactive_ptcl, &nactive_glb, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
			nactive_cum += nactive_glb;

			if (all_active) all_steps++;

			const double t0 = mytimer::get_wtime();
			if (nactive_cum > (unsigned long long)(nfac*global_n) || problem_force_distribute())
//			if ((nactive_cum > (unsigned long long)(nfac*global_n) && all_active) || problem_force_distribute())
			{
				const unsigned long long n0 = nactive_cum;
				all_steps = 0;
				nactive_cum = 0;

				if (myproc == 0)
					fprintf(stderr, " ---- DISTRIBUTE  nactive_cum= %g %g  global_n= %llu  nactive_glb= %llu---- \n", 
							1.0*n0/global_n, 1.0*nactive_cum/global_n, global_n, nactive_glb);
				distribute_data(true, true, true);

				scheduler.flush_list();
				boundary_n = 0;
				for (int i = 0; i < (int)local_n; i++)
				{
					scheduler.push_particle(i, (int)ptcl_local[i].rung);
					if (ptcl_local[i].is_boundary())
						boundary_n++;
				}

				unsigned long long boundary_glb;
				MPI_Allreduce(&boundary_n, &boundary_glb, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

				if (myproc == 0)
					fprintf(stderr, "boundary_glb= %lld\n", boundary_glb);


#if 0
				if (all_active) 
					for (int i = 0; i < (int)local_n; i++)
					{
						if (ptcl_local[i].boundary == Particle::NO_BOUNDARY)
						{
							const real f = 1.0e-10;
							assert(ptcl_local[i].tlast == t_global);
#if 0
							if (!(std::abs(ptcl_local[i].volume - ptcl_local[i].volume_new) < f*ptcl_local[i].volume))
							{
								fprintf(stderr, "proc= %d : v= %g  vnew= %g dv =%g\n",
										myproc, ptcl_local[i].volume, ptcl_local[i].volume_new, 
										(ptcl_local[i].volume_new - ptcl_local[i].volume)/ptcl_local[i].volume_new);
							}
							assert(std::abs(ptcl_local[i].volume - ptcl_local[i].volume_new) < f*ptcl_local[i].volume);
							Ulocal[i] = Ulocal[i].to_primitive(ptcl_local[i].volume).to_conservative(ptcl_local[i].volume_new);
#endif
						}
						ptcl_local[i].volume = ptcl_local[i].volume_new;
					}
#endif
			}
			tprofile[DISTRIBUTE] = mytimer::get_wtime() - t0;
		}
	}

	void system::dump_profile_info() 
	{
		MPI_Barrier(MPI_COMM_WORLD);

		unsigned long data_inout_glob;
		MPI_Allreduce(&myMPI::data_inout, &data_inout_glob, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
		double data_inout_time_glob;
		MPI_Allreduce(&myMPI::data_inout_time, &data_inout_time_glob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);


		if (myproc == 0) 
		{
			fprintf(stderr, " PROC |");
			for (int i = 0; i < NPROFILING_TYPE; i++) 
				switch(i) 
				{
					case ACTIVEP:        fprintf(stderr, " ACTIVEP |"); break;
					case ACTIVEF:        fprintf(stderr, " ACTIVEF |"); break;
					case ITERATE:        fprintf(stderr, " ITERATE |"); break;
					case PREDICT:        fprintf(stderr, " PREDICT |"); break;
					case DISTRIBUTE:     fprintf(stderr, " DISTR   |"); break;
					case COLLECT:        fprintf(stderr, " COLLECT |"); break;
					case DISTW:          fprintf(stderr, " DISTW   |"); break;
					case BATCH:          fprintf(stderr, " BATCH   |"); break;
					case MESH:           fprintf(stderr, " MESH    |"); break;
					case FLUID:          fprintf(stderr, " FLUID   |"); break;
					case UPDATE:         fprintf(stderr, " UPDATE  |"); break;
					case GRADIENTS:      fprintf(stderr, " GRADS   |"); break;
					case RECONSTRCT:     fprintf(stderr, " RECNSTR |"); break;
					case PVEL:           fprintf(stderr, " PVELc   |"); break;
					case TIMESTEP:       fprintf(stderr, " TIMESTEP|"); break;
					case MISC:           fprintf(stderr, " MISC    |"); break;
					case MPICOM:         fprintf(stderr, " MPI     |"); break;
					default: break;
				};
			fprintf(stderr, "\n");

			fprintf(stderr, " %.3d: ",myproc);
			for (int i = 0; i < NPROFILING_TYPE; i++) 
				fprintf(stderr, "%9.3g ", tprofile[i]);
			fprintf(stderr, " %9.3g  [ %6.1f MB/s]  %8d %8d %8d \n", 
					data_inout_time_glob,
					(data_inout_time_glob > 0.0) ? (data_inout_glob/data_inout_time_glob/1e6) : 0.0,
					(int)active_n, (int)(import_n - active_n), (int)local_n);
		}

		MPI_Barrier(MPI_COMM_WORLD);

		tprofile[NPROFILING_TYPE    ] = myMPI::data_inout_time;
		if (myMPI::data_inout_time > 0)
			tprofile[NPROFILING_TYPE + 1] = myMPI::data_inout/myMPI::data_inout_time/1e6;
		else
			tprofile[NPROFILING_TYPE + 1] = 0.0;
		tprofile[NPROFILING_TYPE+2  ] = active_n;
		tprofile[NPROFILING_TYPE+3  ] = import_n - active_n;
		tprofile[NPROFILING_TYPE+4  ] = local_n;


		double tprof_loc[NPROFILING_TYPE + 5];
		for (int p = 1; p < nproc; p++) 
		{
			if (myproc == p) 
				MPI_Send(tprofile, NPROFILING_TYPE + 5, MPI_DOUBLE, 0, 123456789 + p, MPI_COMM_WORLD);

			if (myproc == 0) 
			{
				MPI_Status status;
				MPI_Recv(tprof_loc, NPROFILING_TYPE + 5, MPI_DOUBLE, p, 123456789 + p, MPI_COMM_WORLD, &status);
				fprintf(stderr, " %.3d: ", p);
				for (int i = 0; i < NPROFILING_TYPE; i++) 
					fprintf(stderr, "%9.3g ", tprof_loc[i]);
				fprintf(stderr, " %9.3g  [ %6.1f MB/s]  %8d %8d %8d \n", 
						tprof_loc[NPROFILING_TYPE    ],
						tprof_loc[NPROFILING_TYPE + 1],
						(int)tprof_loc[NPROFILING_TYPE + 2],
						(int)tprof_loc[NPROFILING_TYPE + 3],
						(int)tprof_loc[NPROFILING_TYPE + 4]
						);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}

		MPI_Barrier(MPI_COMM_WORLD);

		{

			Scheduler sc(scheduler);
			sc.flush_list();

			for (int i = 0; i < (int)local_n; i++)
			{
				if (ptcl_local[i].boundary == 0)
					sc.push_particle(i, (int)ptcl_local[i].rung);
				else
					sc.push_particle(i, 30);
			}

			int nrung[NMAXPROC][Scheduler::RUNGMAX + 1];
			for (int i = 0; i < Scheduler::RUNGMAX; i++)
			{
				//        nrung[myproc][i] = sc.list[i].size();
				nrung[myproc][i] = scheduler.list[i].size();
			}
			nrung[myproc][Scheduler::RUNGMAX] = active_ptcl.size();


			for (int p = 1; p < nproc; p++)
			{
				if (myproc == p)
				{
					MPI_Send(nrung[p], Scheduler::RUNGMAX+1, MPI_INT, 0, 234567891 + p, MPI_COMM_WORLD);
				}
				if (myproc == 0)
				{
					MPI_Status stat;
					MPI_Recv(nrung[p], Scheduler::RUNGMAX+1, MPI_INT, p, 234567891 + p, MPI_COMM_WORLD, &stat);
				}
			}

			MPI_Barrier(MPI_COMM_WORLD);

			if (myproc == 0)
			{
				fprintf(stderr, " r\\p: ");
				for (int p = 0; p < nproc; p++)
					fprintf(stderr, " %4d   ", p);
				fprintf(stderr, "\n");
				for (int i = 0; i < Scheduler::RUNGMAX; i++)
				{
					fprintf(stderr, " %2d : ", i);
					for (int p = 0; p < nproc; p++)
						if (nrung[p][i] > 0)
							fprintf(stderr, "%4.1e ", (double)nrung[p][i]);
						else
							fprintf(stderr, "      0 ");
					fprintf(stderr, " \n");
				}
				fprintf(stderr, "actv: ");
				const int i = Scheduler::RUNGMAX;
				for (int p = 0; p < nproc; p++)
					if (nrung[p][i] > 0)
						fprintf(stderr, "%4.1e ", (double)nrung[p][i]);
					else
						fprintf(stderr, "      0 ");
				fprintf(stderr, " \n");

			}

		}

		MPI_Barrier(MPI_COMM_WORLD);
	}  

	const Energy system::get_energy() const
	{
		return Energy(local_n, ptcl_local, U_local);//, gpot, gpot_body);
	}

}
