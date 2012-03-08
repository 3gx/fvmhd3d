#include "fvmhd3d.h"

namespace fvmhd3d
{

	void System::dump_binary(const CkVec<char> &fnVec, const int global_n, CkCallback& cb)
	{

#if 0
		return;
#endif

    char filename[256];
    sprintf(filename, "%s_part%.6d", &fnVec[0], thisIndex);

		FILE *fout; 
		if (!(fout = fopen(&filename[0], "w"))) 
    {
      CkPrintf("Cannot open file %s \n", &filename[0]);
      CkExit();
		}

		int ival;
		float fval;

#define fdump(x) { fval = x; fwrite(&fval, sizeof(float), 1, fout); }
#define idump(x) { ival = x; fwrite(&ival, sizeof(int),   1, fout); }

		idump(20*4);
		idump(thisIndex);
		idump(numElements);
		idump(numElements);
		union {
			unsigned long long uint_long;
			unsigned int       uint[2];
		} data;
		data.uint_long = scheduler.tsysU;
		idump(data.uint[0]);
		idump(data.uint[1]);


		idump(global_n);
		idump(local_n);
		idump(3);
		fdump(t_global);
		fdump(dt_global);
		idump(iteration);
		fdump(courant_no);
		fdump(gamma_gas);

		//const int periodic_on = (bc == BC_PERIODIC) ? 1 : 0;
		//  idump(periodic_on);
		idump(-1);

		fdump(global_domain.get_rmin().x);
		fdump(global_domain.get_rmin().y);
		fdump(global_domain.get_rmin().z);
		fdump(global_domain.get_rmax().x);
		fdump(global_domain.get_rmax().y);
		fdump(global_domain.get_rmax().z);
		idump(20*4);

		for (int i = 0; i < (int)local_n; i++) {
			const MeshPoint &pi = mesh_pnts[i];
			assert(pi.Volume > 0.0);
			const Fluid mi = U_list[i].to_primitive(pi.Volume);

			const real L = std::pow(pi.Volume*3.0/4.0/M_PI, 1.0/3.0);

			idump(26*4);
			ival = pi.idx; idump(ival);
			vec3 pos = periodic(pi.pos_orig);
			fdump(pos.x);
			fdump(pos.y);
			fdump(pos.z);
			fdump(pi.vel.x);
			fdump(pi.vel.y);
			fdump(pi.vel.z);
			fdump(mi[Fluid::DENS]);
			fdump(mi[Fluid::ETHM]);
			fdump(Problem_compute_pressure(mi));
			fdump(-1.0); //pi.rmax);      //   fdump(    (sqr(mi[Fluid::BX  ]) + sqr(mi[Fluid::BY  ]) + sqr(mi[Fluid::BZ  ]))*0.5f);
			idump(pi.boundary);  //		fdump(sqrt(sqr(mi[Fluid::VELX]) + sqr(mi[Fluid::VELY]) + sqr(mi[Fluid::VELZ]))*0.5f);
			fdump(mi[Fluid::VELX]);
			fdump(mi[Fluid::VELY]);
			fdump(mi[Fluid::VELZ]);
			fdump(mi[Fluid::BX]);
			fdump(mi[Fluid::BY]);
			fdump(mi[Fluid::BZ]);
			fdump(scheduler.get_dt(mesh_pnts[i].rung));
			fdump(pi.Volume);
			fdump(mi[Fluid::PSI]);
			fdump(L*divBi[i]);
			fdump(mi[Fluid::ENTR]);
			fdump(-1.0); //Wextra_local[i].J.x);
			fdump(-1.0); //Wextra_local[i].J.y);
			fdump(-1.0); //Wextra_local[i].J.z);
			idump(26*4);

		}

		fclose(fout);

    contribute(cb);
	}

	inline void myfread(void *f, const int sz, const int n, FILE *fin)
	{
		const size_t len = fread(f, sz, n, fin);
		assert((int)len == n);
	}

	void System::read_binary(const CkVec<char> &filename)
	{
    assert(false);
#if 0
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

#if 0
		global_n = U_local.size();		      
		local_n = global_n;
#endif

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
