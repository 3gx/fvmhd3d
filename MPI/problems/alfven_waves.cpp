#include "fvmhd3d.h"

#if 1
namespace fvmhd3d
{


  bool system::compute_pvel1()
	{
		return false;
	}
	
	void system::predict(const int i)
	{
//		ptcl_local[i].vel = 0.0;
		ptcl_local[i].acc = 0.0;
		const real dt = t_global - ptcl_local[i].tlast;
		ptcl_local[i].pos_pred = ptcl_local[i].pos +
		 	(ptcl_local[i].vel  + ptcl_local[i].acc*(dt*0.5))*dt;
	}

	bool system::compute_update_prob(Fluid &U, const int i)
	{
		return false;
	}
  const vec3 system::gpot_acc(const int i, real &gpot) const
	{
		gpot = 0.0;
		return vec3(0.0);
	}
	real system::compute_ethm_update(const Fluid &W, const int i) const
	{
		return W[Fluid::ETHM];
	}
	real system::compute_pressure(const Fluid& W) const
	{
		return gamma_gas > 1.0 ? (gamma_gas - 1.0) * W[Fluid::ETHM] : W[Fluid::ETHM];
	}
	real system::compute_entropy_from_ethm(const Fluid& W) const
	{
		assert(gamma_gas > 1.0);
		return (gamma_gas - 1.0) * W[Fluid::ETHM]/std::pow(W[Fluid::DENS], gamma_gas);
	}
	real system::compute_ethm_from_entropy(const Fluid& W) const
	{
		assert(gamma_gas > 1.0);
		return W[Fluid::ENTR] * std::pow(W[Fluid::DENS], gamma_gas)/(gamma_gas - 1.0);
	}
	real system::compute_pressure_gradient(const Fluid &W, const Fluid &dW) const
	{
		const real dpdethm = (gamma_gas > 1.0) ? (gamma_gas - 1.0) : 0.0;
		const real dpddens = (gamma_gas > 1.0) ?              0.0  : 0.0;

		return dpddens * dW[Fluid::DENS] + dpdethm * dW[Fluid::ETHM];
	}

	bool system::refine(const int i) 
	{
		return false;
	}

	bool system::derefine(const int i)
	{
		return false;
	}


	void system::set_geometry(const bool init) 
	{
		const double dt_max = 1.0/64;
		scheduler = Scheduler(dt_max);

		t_end   = 5.0;
		
		n_restart = 1;
		dt_restart = dt_max * 2;

		dt_dump = dt_max / 64;
		dt_dump = dt_max * 4;

		di_log = 100;

		global_n = local_n = 0;
		int nx = 16;
		int ny = 16;
		int nz = 16;

		nx = ny = nz = 32;

		nx = ny = nz = 64;

		nx = ny = 256; nz = 16;

    nx = 32;
    ny = 16;
    nz = 16;

		dt_dump = dt_max;
		dt_restart = 1e10;

		const double Lx = 1.0;
		const vec3 rmin(0.0);
		const vec3 rmax(Lx, (Lx/nx)*ny, (Lx/nx)*nz);
		global_domain = boundary(rmin, rmax);
		global_domain_size = global_domain.hsize() * 2.0;

		const vec3 Len3 = global_domain.hsize() * 2.0;
		pfloat<0>::set_scale(Len3.x);
		pfloat<1>::set_scale(Len3.y);
		pfloat<2>::set_scale(Len3.z);
		
    Distribute::int3 nt(1, 1, 1);
		switch(nproc) 
    {
			case 1: break;
			case 2: nt.x = 2; nt.y = 1; nt.z = 1; break;
			case 4: nt.x = 2; nt.y = 2; nt.z = 1; break;
			case 6: nt.x = 3; nt.y = 2; nt.z = 1; break;
			case 8: nt.x = 4; nt.y = 2; nt.z = 1; break;
			case 10: nt.x = 5; nt.y = 2; nt.z = 1; break;
			case 12: nt.x = 4; nt.y = 3; nt.z = 1; break;
			case 14: nt.x = 7; nt.y = 2; nt.z = 1; break;
			case 16: nt.x = 4; nt.y = 4; nt.z = 1; break;
			case 32: nt.x = 8; nt.y = 4; nt.z = 1; break;
			case 64: nt.x = 8; nt.y = 8; nt.z = 1; break;
			case 128: nt.x = 8; nt.y = 8; nt.z = 2; break;
			default: assert(false);
		}

		const Distribute::int3 nt_glb(nt);
		const pBoundary pglobal_domain(pfloat3(0.0), pfloat3(Len3));
		distribute_glb.set(nproc, nt, pglobal_domain);

    assert(init);

		if (myproc == 0) 
    {

			ptcl_local.clear();
			ptcl_local.reserve(128);

			const dvec3 dr = dvec3(Len3.x/nx, Len3.y/ny, Len3.z/nz);
			const real rmax = dr.abs() * 1.0;

			fprintf(stderr, "dr= %g %g %g \n", dr.x, dr.y, dr.z);
			fprintf(stderr, "rmin= %g %g %g \n", 
					global_domain.rmin.x,
					global_domain.rmin.y,
					global_domain.rmin.z);
			fprintf(stderr, "rmax= %g %g %g \n", 
					global_domain.rmax.x,
					global_domain.rmax.y,
					global_domain.rmax.z);

			for (int k = 0; k < nz; k++) {
				for (int j = 0; j < ny; j++) {
					for (int i = 0; i < nx; i++) {
						dvec3 pos = global_domain.rmin + dvec3(i*dr.x, j*dr.y, k*dr.z) + 0.5*dr;
						const int ijk = (k*ny + j)*nx + i;
#if 1
						{
							const real f = 1.0e-6;
							pos += vec3(drand48()*dr.x*f, drand48()*dr.y*f, drand48()*dr.z*f);
						}
#endif


#if 1
						pos = global_domain.rmin + dvec3(
								drand48()*Len3.x,
								drand48()*Len3.y,
								drand48()*Len3.z);
#else
#define _UNIFORM_MESH_
#endif
						dvec3 vel(0.0, 0.0, 0.0);
						Particle p;
						p.set_pos(pos);
						p.vel = vel;
						p.boundary = 0;
						p.idx = ijk;
						p.rmax = rmax;
						ptcl_local.push_back(p);
					}
				}
			}
			local_n  = ptcl_local.size();
			global_n = local_n;

			fprintf(stderr, "  *** proc= %d : local_n= %llu  global_n= %llu \n", myproc, local_n, global_n);
		} // myproc == 0

		MPI_Bcast(&global_n,  1, MPI_INT, 0, MPI_COMM_WORLD);
		
    MPI_Barrier(MPI_COMM_WORLD);
		if (myproc == 0)
			fprintf(stderr, " ***  Distrubiting data \n");

		all_active = true;

		for (int k = 0; k < 5; k++)
			distribute_data(false,false);

#if 0
		std::vector< std::pair<int, TREAL> > rmax_list;
		local_tree.root.get_rmax(rmax_list);
		assert((int)rmax_list.size() == local_n);
		for (int i = 0; i < local_n; i++)
			ptcl[rmax_list[i].first].rmax = rmax_list[i].second;
#endif

		MPI_Barrier(MPI_COMM_WORLD);
		fprintf(stderr, " *** proc= %d : local_n= %llu  global_n= %llu \n", myproc, local_n, global_n);
		fprintf(stderr, " proc= %d  relax \n", myproc);

#ifndef _UNIFORM_MESH_
		relax_mesh(3);
#endif
		fprintf(stderr, " ---- done --- \n");
		{
			distribute_data(false, false);
			const double t10 = mytimer::get_wtime();
			clear_mesh(false);
			int nattempt = build_mesh_global();
			double dt10 = mytimer::get_wtime() - t10;

			double volume_loc = 0.0;
			{
				std::vector<TREAL> v(local_n);
				for (int i = 0; i < (int)local_n; i++)
					v[i] = cell_local[i].Volume;
				std::sort(v.begin(), v.end());  // sort volumes from low to high, to avoid roundoff errors
				for (int i = 0; i < (int)local_n; i++)
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
		}


		extract_ngb_from_mesh();
  }

	void system::set_problem(const bool init) 
	{
		if (myproc == 0)
			fprintf(stderr, " ********* Setting up acoustic wave problem ************* \n");


    assert(init);
		
    gamma_gas = 5.0/3.0;
		courant_no = 0.4;

		W_st.resize(local_n);

		double dt_min = HUGE;
		for (int i = 0; i < (int)local_n; i++) 
    {
			const Particle &pi = ptcl_local[i];

			real x = pi.pos.x;

			real d, p, vx, vy, vz, bx, by, bz;

			vx = vy = vz = 0.0;
			bx = by = bz = 0.0;

			d = 1.0;
			p = 0.1;

			vx = 0.0;
			bx = 1.0;

			const real kwave = 2.0;

			vy = by = 1.0e-4*sin(kwave * 2*M_PI*x);
			real scal = 1.0;


			Fluid m;

			m[Fluid::DENS] = d ;
			m[Fluid::ETHM] = p/(gamma_gas - 1.0);
			m[Fluid::VELX] = vx;
			m[Fluid::VELY] = vy;
			m[Fluid::VELZ] = vz;
			m[Fluid::BX  ] = bx;
			m[Fluid::BY  ] = by;
			m[Fluid::BZ  ] = bz;
			m[Fluid::PSI ] = 0.0;
			m[Fluid::ENTR] = compute_entropy_from_ethm(m);
			for (int k = 0 ; k < Fluid::NSCALARS; k++)
				m.scal(k) = scal;
			
      W_st[i].w  = m;
			Ulocal [i] = m.to_conservative(cell_local[i].Volume);
			dUlocal[i] = 0.0;
			ptcl_local[i].volume = cell_local[i].Volume;
			real gpot;
			ptcl_local[i].acc = gpot_acc(i, gpot);

			const double L      = std::pow(cell_local[i].Volume, 1.0/3);
			const double cs_est = std::sqrt((p*gamma_gas + (bx*bx+by*by+bz*bz))/d);
			const double v_est  = std::sqrt(vx*vx + vy*vy + vz*vz);
			const double dt_est = 0.1 * courant_no * L/(cs_est + v_est);

			ptcl_local[i].tlast = 0.0;
			ptcl_local[i].rung = scheduler.get_rung(dt_est);

			dt_min = std::min(dt_min, dt_est);
		}

		double dt_min_glob;
		MPI_Allreduce(&dt_min, &dt_min_glob, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);
		if (myproc == 0)
			fprintf(stderr , " pvel ... \n");

		get_active_ptcl(true);

		MPI_Barrier(MPI_COMM_WORLD);
		if (myproc == 0)
			fprintf(stderr , " pvel ... \n");


		cell_list.swap(cell_local);
		ptcl_import.swap(ptcl_local);
		site_active_list.swap(active_ptcl);


		compute_pvel();
		compute_timesteps(true);

		cell_list.swap(cell_local);
		ptcl_import.swap(ptcl_local);
		site_active_list.swap(active_ptcl);

		for (int i = 0; i < (int)local_n; i++)
			ptcl_local[i].rung += 4;
		all_active = true;
		scheduler.flush_list();
		for (int i = 0; i < (int)local_n; i++)
			scheduler.push_particle(i, (int)ptcl_local[i].rung);

		clear_mesh(true);

		MPI_Barrier(MPI_COMM_WORLD);
		if (myproc == 0) fprintf(stderr, " proc= %d: complete problem setup \n", myproc);
	}
}
#endif

