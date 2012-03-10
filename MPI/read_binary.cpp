#include "fvmhd3d.h"

namespace fvmhd3d {

	inline void myfread(void *f, const int sz, const int n, FILE *fin)
	{
		const size_t len = fread(f, sz, n, fin);
		assert((int)len == n);
	}

	void system::read_binary(const char *filename, const int n_files) 
	{
#if 1
		assert(n_files == 1);

		vec3 rmin, rmax;
		if (myproc == 0) 
		{

			FILE *fin; 
			if (!(fin = fopen(filename, "r"))) 
			{
				std::cerr << "Cannot open file " << filename << std::endl;
				exit(-1);
			}

			std::cerr << "proc= " << myproc << " read snapshot: " << filename << std::endl;

			int ival;
			float fval;

#define fload(x) { myfread(&fval, sizeof(float), 1, fin); x = fval; }
#define iload(x) { myfread(&ival, sizeof(int),   1, fin); x = ival;}

			float ftmp;
			int itmp, np0, npx, npy, npz;
			iload(itmp); // 20*4
			assert(itmp == 20*4);
			iload(itmp); // myid
			iload(np0);
			iload(npx);

			union 
			{
				unsigned long long uint_long;
				unsigned int       uint[2];
			} data;
			iload(data.uint[0]);
			iload(data.uint[1]);
			scheduler.tsysU = data.uint_long;

			float courant_No;
			int nglob, nloc, ndim;
			iload(nglob);
			iload(nloc);
			iload(ndim);  
			assert(ndim == 3);
			fload(t_global);
			fload(dt_global);
			iload(iteration);
			fload(courant_No);
			fload(gamma_gas);

			int periodic_on;
			iload(periodic_on);
			assert(periodic_on == -1);
			
			fload(rmin.x);
			fload(rmin.y);
			fload(rmin.z);
			fload(rmax.x);
			fload(rmax.y);
			fload(rmax.z);
			iload(itmp);     // 20*4
			assert(itmp == 20*4);

			ptcl_local.resize(nglob);
			U_local.resize(nglob);
			dU_local.resize(nglob);

			fprintf(stderr, "np =%d   nglob= %d \n", np0, nglob);

			int pc = 0;
			for (int pr = 0; pr < np0; pr++) 
			{
				fprintf(stderr, " p= %d out of %d; nloc= %d\n", pr, np0, nloc);
				for (int i = 0; i < nloc; i++) 
				{
					Particle p;
          p.tend  = t_global;
          p.rung  = 0.0;
          p.new_dt = 0.0;
          p.local_id = i;
					Fluid W(0.0);

					iload(ival);    
					assert(ival == 26*4);
					iload(ival); p.idx = ival;
				
					fload(p.pos.x); 
					fload(p.pos.y);
					fload(p.pos.z);

					p.pos = periodic(p.pos);

					assert(rmin.x <= p.pos.x);
					assert(rmax.x >= p.pos.x);
					assert(rmin.y <= p.pos.y);
					assert(rmax.y >= p.pos.y);
					assert(rmin.z <= p.pos.z);
					assert(rmax.z >= p.pos.z);
          p.orig_pos = p.pos;
          p.pot = 0;


					fload(p.vel.x);
					fload(p.vel.y);
					fload(p.vel.z);
          p.orig_vel = p.vel;
					fload(W[Fluid::DENS]);
					fload(W[Fluid::ETHM]);
					fload(ftmp); // compute_pressure(m.dens, m.ethm));
					fload(p.rmax);     //dump(    (sqr(m.B.x  ) + sqr(m.B.y  ) + sqr(m.B.z  ))*0.5f);
					iload(p.boundary); //				fload(ftmp); //dump(sqrt(sqr(m.vel.x) + sqr(m.vel.y) + sqr(m.vel.z))); 
					fload(W[Fluid::VELX]);
					fload(W[Fluid::VELY]);
					fload(W[Fluid::VELZ]);
					fload(W[Fluid::BX]);
					fload(W[Fluid::BY]);
					fload(W[Fluid::BZ]);
					float h;
					fload(h);
					fload(p.volume);
          p.volume_new = p.volume;
					fload(W[Fluid::PSI]);
					fload(ftmp); //L*divB_i[i]);
					fload(W[Fluid::ENTR]);
					fload(ftmp); // Jx
					fload(ftmp); // Jy
					fload(ftmp); // Jz
					iload(ival); 
					assert(ival == 26*4);

					p.tlast = t_global;

					ptcl_local[pc] = p;
					U_local    [pc] = W;
          dU_local   [pc] = 0.0;
          dU_local   [pc] = 0.0;
					pc++;
				}

				fprintf(stderr, "p= %d  np0= %d size= %d %d\n",
						pr, np0, (int)U_local.size(), (int)ptcl_local.size());
				if (!(pr < np0-1)) break;
				iload(itmp); // 20*4
				assert(itmp == 20*4);
				iload(itmp); // myid
				iload(np0);
				iload(npx);
				iload(npy);
				iload(npz);

				int nglob1;
				iload(nglob1);
				if (nglob != nglob1) {
					fprintf(stderr, "np; npx, npy, npz = %d; %d %d %d \n", 
							np0, npx, npy, npz);
					fprintf(stderr, "nglob= %d  nglob1= %d\n", nglob, nglob1);
				}
				assert(nglob == nglob1);
				iload(nloc);
				iload(ndim);
				fload(t_global);
				fload(dt_global);
				iload(iteration);
				fload(courant_No);
				fload(gamma_gas);

				iload(periodic_on);

				fload(rmin.x);
				fload(rmin.y);
				fload(rmin.z);
				fload(rmax.x);
				fload(rmax.y);
				fload(rmax.z);
				iload(itmp);     // 20*4
			}
			assert(pc == nglob);
			assert(nglob == (int)U_local.size());
			fclose(fin);


			local_n = U_local.size();
		}

		global_n = U_local.size();		      
		local_n = global_n;

		MPI_Bcast(&global_n,   1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
		MPI_Bcast(&iteration,  1, MPI_INT, 0, MPI_COMM_WORLD);

		double dt_glob = dt_global;
		double  t_glob =  t_global;
		MPI_Bcast(& t_glob,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&dt_glob,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&scheduler.tsysU, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
		dt_global = dt_glob;
		t_global  =  t_glob;
    scheduler.set_tsys(t_global);
		assert(t_global == scheduler.get_tsys());


//		scheduler.tsysU = (unsigned long long)(t_global / scheduler.dt_tick);
		scheduler.min_rung = 0;

		dt_global = 0.0f;
		
    distribute_data(true, false, true);
	
#if 1	
		fit_vec(ptcl_local);
		fit_vec(U_local);
		fit_vec(dU_local);
		fit_vec(Wrec_local);
#endif

		all_active = true;
		
		MPI_Barrier(MPI_COMM_WORLD);

		for (int i = 0; i < (int)local_n; i++)
		{
			ptcl_local[i].tlast = t_global;
//      ptcl_local[i].volume = cell_local[i].Volume;

      Wrec_local[i] = Fluid_rec(U_local[i]);
      U_local[i] = U_local[i].to_conservative(ptcl_local[i].volume);
      dU_local[i] = 0.0;
		}
    
    MPI_Barrier(MPI_COMM_WORLD);
    if (myproc == 0)
      fprintf(stderr , " pvel ... \n");

    get_active_ptcl(true);

    MPI_Barrier(MPI_COMM_WORLD);
    if (myproc == 0)
      fprintf(stderr , " pvel ... \n");


    cell_list.swap(cell_local);
    ptcl_import.swap(ptcl_local);
    U_import.swap(U_local);
    site_active_list.swap(active_ptcl);

    compute_pvel();
    compute_timesteps(true);

    cell_list.swap(cell_local);
    ptcl_import.swap(ptcl_local);
    U_import.swap(U_local);
    site_active_list.swap(active_ptcl);

    for (int i = 0; i < (int)local_n; i++)
    {
      ptcl_local[i].rung += 1;
      ptcl_local[i].tend  = ptcl_local[i].tlast + scheduler.get_dt(ptcl_local[i].rung);
      ptcl_local[i].orig_vel = ptcl_local[i].vel;
      ptcl_local[i].unset_active();
    }
    all_active = true;
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

    clear_mesh(true);

    MPI_Barrier(MPI_COMM_WORLD);
    if (myproc == 0) fprintf(stderr, " proc= %d: complete read_binary \n", myproc);
#endif
  }

}
