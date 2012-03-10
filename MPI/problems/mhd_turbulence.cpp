#include "venom3d.h"

namespace venom3d {

	const char fin_desc[] = "mhd256.inp";
	const char fin_data[] = "mhd256.dat";

	const real cs2 = 1.0;

	real system::compute_pressure(const real dens, const real ethm) const
	{
		return cs2 * dens;
	}
	real system::compute_entropy_from_ethm(const real dens, const real ethm) const
	{
		return 1.0;
	}
	real system::compute_ethm_from_entropy(const real dens, const real entr) const
	{
		assert(false);
		return -1.0;
	}
	real system::compute_pressure_gradient(
			const real  dens, const real  ethm, const real pres,
			const real ddens, const real dethm) const
	{
		const real dpdethm = 0.0;
		const real dpddens = cs2;

		return dpddens * ddens + dpdethm * dethm;
	}

	void system::set_geometry(const bool init) 
	{
		const double dt_max = 1.0/512;
		scheduler = Scheduler(dt_max);

		int nx, ny, nz;
		double lx, ly, lz;
		if (myproc == 0)
		{
			std::ifstream fin(fin_desc);
			int nvar;
			fin >> nx >> ny >> nz;
			fin >> nvar;

			fin >> lx >> ly >> lz;

			fprintf(stderr, " nx= %d  ny= %d  nz= %d \n", nx, ny, nz);
			fprintf(stderr, " nvar= %d \n", nvar);
			fprintf(stderr, " lx= %g  ly= %g  lz= %g \n", lx, ly, lz);
		}

		MPI_Bcast(&nx,  1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&ny,  1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&nz,  1, MPI_INT, 0, MPI_COMM_WORLD);

		MPI_Bcast(&lx,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&ly,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&lz,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


		t_end   = 0.2;

		n_restart = 2;
		dt_restart = dt_max;

		dt_dump = 0.01;

		di_log = 100;

		global_n = local_n = 0;

//		eulerian = true;

		const vec3 rmin(0.0);
		const vec3 rmax(lx, ly, lz);
		global_domain = boundary(rmin, rmax);
		global_domain_size = global_domain.hsize() * 2.0;

		const vec3 Len3 = global_domain.hsize() * 2.0;
		pfloat<0>::set_scale(Len3.x);
		pfloat<1>::set_scale(Len3.y);
		pfloat<2>::set_scale(Len3.z);

		if (myproc == 0) {

			ptcl.clear();
			ptcl.reserve(128);

			const dvec3 dr = dvec3(Len3.x/nx, Len3.y/ny, Len3.z/nz);
			const real rmax = dr.abs() * 2.0;

			fprintf(stderr, "dr= %g %g %g \n", dr.x, dr.y, dr.z);

			for (int k = 0; k < nz; k++) 
				for (int j = 0; j < ny; j++) 
					for (int i = 0; i < nx; i++) {
						dvec3 pos = global_domain.rmin + dvec3(i*dr.x, j*dr.y, k*dr.z) + 0.5*dr;
						const int ijk = (k*ny + j)*nx + i;

						dvec3 vel(0.0, 0.0, 0.0);
						Particle p(pos.x, pos.y, pos.z, vel.x, vel.y, vel.z, ijk);
						p.rmax = rmax;
						p.unset_derefine();
						ptcl.push_back(p);
					}

			local_n  = ptcl.size();
			global_n = local_n;

			FILE *fin = fopen(fin_data, "r");

			U.resize(local_n);
			const int var_list[7] = {
				Fluid::DENS,
				Fluid::VELX,
				Fluid::VELY,
				Fluid::VELZ,
				Fluid::BX,
				Fluid::BY,
				Fluid::BZ};
		
			std::vector<float> data(local_n);
			for (int var = 0; var < 7; var++)
			{
				fprintf(stderr, " read var %d out of %d\n", var+1, 7);
				const size_t nread = fread(&data[0], sizeof(float), local_n, fin);
				assert((int)nread == local_n);
				for (int i = 0; i < local_n; i++)
					U[i][var_list[var]] = data[i];
			}
			for (int i = 0; i < local_n; i++)
			{
				assert(U[i][Fluid::DENS] > 0.0);
				U[i][Fluid::ETHM] = cs2 * U[i][Fluid::DENS];
			}
			fclose(fin);

			fprintf(stderr, "  *** proc= %d : local_n= %d  global_n= %d \n", myproc, local_n, global_n);
		} // myproc == 0

		MPI_Bcast(&global_n,  1, MPI_INT, 0, MPI_COMM_WORLD);

		fprintf(stderr, " proc= %d  distrubite \n", myproc);
		MPI_Barrier(MPI_COMM_WORLD);

		Distribute::int3 nt(1, 1, 1);
		switch(nproc) {
			case 1: break;
			case 2: nt.x = 2; nt.y = 1; nt.z = 1; break;
			case 4: nt.x = 2; nt.y = 2; nt.z = 1; break;
			case 6: nt.x = 3; nt.y = 2; nt.z = 1; break;
			case 8: nt.x = 2; nt.y = 2; nt.z = 2; break;
			case 16: nt.x = 4; nt.y = 2; nt.z = 2; break;
			case 32: nt.x = 4; nt.y = 4; nt.z = 2; break;
			case 64: nt.x = 4; nt.y = 4; nt.z = 4; break;
			case 128: nt.x = 8; nt.y = 4; nt.z = 4; break;
			case 256: nt.x = 8; nt.y = 8; nt.z = 4; break;
			case 512: nt.x = 8; nt.y = 8; nt.z = 8; break;
			default: assert(false);
		}

		const Distribute::int3 nt_glb(nt);
		const pBoundary pglobal_domain(pfloat3(0.0), pfloat3(Len3));
		distribute_glb.set(nproc, nt, pglobal_domain);

		for (int k = 0; k < 5; k++)
			distribute_data(true, false);
		
		const int nloc_reserve = (int)(2.0*global_n/nproc);
		fit_reserve_vec(ptcl,      nloc_reserve);
		fit_reserve_vec(ptcl_ppos, nloc_reserve);
		fit_reserve_vec(U,         nloc_reserve);
		fit_reserve_vec(dU,        nloc_reserve);
		fit_reserve_vec(Wgrad,     nloc_reserve);
		fit_reserve_vec(gradPsi,   nloc_reserve);
		fit_reserve_vec(cells,     nloc_reserve);

		MPI_Barrier(MPI_COMM_WORLD);

		fprintf(stderr, " *** proc= %d : local_n= %d  global_n= %d \n", myproc, local_n, global_n);
		fprintf(stderr, " proc= %d  building_mesh \n", myproc);

		MPI_Barrier(MPI_COMM_WORLD);



		const double t10 = mytimer::get_wtime();
		clear_mesh();
		int nattempt = build_mesh(true);
		double dt10 = mytimer::get_wtime() - t10;

		double volume_loc = 0.0;
		{
			std::vector<TREAL> v(local_n);
			for (int i = 0; i < local_n; i++)
				v[i] = cells[i].Volume;
			std::sort(v.begin(), v.end());  // sort volumes from low to high, to avoid roundoff errors
			for (int i = 0; i < local_n; i++)
				volume_loc += v[i];
		}


		double dt10max;
		MPI_Allreduce(&dt10, &dt10max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		double volume_glob = 0.0;	
		int    nattempt_max, nattempt_min;
		MPI_Allreduce(&volume_loc, &volume_glob,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&nattempt,   &nattempt_max, 1, MPI_INT,    MPI_MAX, MPI_COMM_WORLD);
		MPI_Allreduce(&nattempt,   &nattempt_min, 1, MPI_INT,    MPI_MIN, MPI_COMM_WORLD);

		const double volume_exact = global_domain_size.x*global_domain_size.y*global_domain_size.z;
		if (myproc == 0)
		{
			fprintf(stderr, "first call build_mesh:[ %g  sec ::  %g cells/s/proc/thread ]\n",
					dt10max,
					global_n/nproc/dt10max);
			fprintf(stderr, "   computed_volume= %g  exact_volume= %g diff= %g [ %g ]  nattempt= %d %d \n",
					volume_glob, volume_exact, 
					volume_glob - volume_exact,	(volume_glob - volume_exact)/volume_exact,
					nattempt_min, nattempt_max);
		}

		exchange_ptcl();

	}

	void system::set_problem(const bool init) 
	{
		if (myproc == 0)
			fprintf(stderr, " ********* Setting up MHD Turbulence ************* \n");

		const int reserve_n = (int)(1.25*local_n);
		U.reserve(reserve_n);
		dU.reserve(reserve_n);
		Wgrad.reserve(reserve_n);

		U.resize(local_n);
		dU.resize(local_n);
		Wgrad.resize(local_n);


		gamma_gas = 1.0;
		courant_no = 0.4;
		
		for (int i = 0; i < local_n; i++) 
		{
			assert(U[i][Fluid::DENS] > 0.0);
			U[i][Fluid::PSI ] = 0.0;
			
			for (int k = 0 ; k < Fluid::NSCALARS; k++)
				U[i].scal(k) = 1.0;
			
			dU[i] = Fluid(0.0);
			Wgrad[i] = 0.0;
			for (int k = 0; k < Fluid::NFLUID; k++)
				Wgrad[i].m[k] = U[i][k];
			U[i] = U[i].to_conservative(cells[i].Volume);

			ptcl[i].Volume = cells[i].Volume;
		}
		entropy_scalar = -1;
		isoeos_flag = true;

		MPI_Barrier(MPI_COMM_WORLD);
		if (myproc == 0)
			fprintf(stderr , " pvel ... \n");

		get_active_ptcl(true);

		MPI_Barrier(MPI_COMM_WORLD);
		if (myproc == 0)
			fprintf(stderr , " primitives ... \n");

		//		exchange_fluid();
		exchange_primitive_and_wdot();



		MPI_Barrier(MPI_COMM_WORLD);
		if (myproc == 0)
			fprintf(stderr , " gradients ... \n");

		compute_pvel();
		exchange_pvel();

		MPI_Barrier(MPI_COMM_WORLD);
		if (myproc == 0)
			fprintf(stderr , " tgradients ... \n");
		compute_tgradient();

		if (myproc == 0)
			fprintf(stderr , " timestep... \n");
		compute_timesteps(true);
		for (int i = 0; i < local_n; i++)
			ptcl[i].rung[0] += 3;

		all_active = true;
		scheduler.flush_list();
		for (int i = 0; i < local_n; i++)
			scheduler.push_particle(i, (int)ptcl[i].rung[0]);

		MPI_Barrier(MPI_COMM_WORLD);
		if (!eulerian)
			clear_mesh();

		if (myproc == 0) fprintf(stderr, " proc= %d: complete problem setup \n", myproc);
		MPI_Barrier(MPI_COMM_WORLD);


	}

};
