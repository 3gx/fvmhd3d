#include <sys/stat.h>
#include "fvmhd3d.h"

#ifndef __MACOSX_
#define __LINUX__
#endif

#ifdef __MACOSX__
#include <Accelerate/Accelerate.h>
#include <xmmintrin.h>
inline void fpe_catch() {
  _mm_setcsr( _MM_MASK_MASK &~
              (_MM_MASK_OVERFLOW|_MM_MASK_INVALID|_MM_MASK_DIV_ZERO) );
}
#elif defined __LINUX__
#include <fenv.h>
void fpe_catch(void) {
  /* Enable some exceptions. At startup all exceptions are masked. */
  feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}
#else
crap
void fpe_catch(void) {}
#endif

int gwait = 0;
// #define _SEQ_WRITE_
#define _CREATE_FOLDER_

int main(int argc, char * argv[])
{

	MPI_Init(&argc, &argv);

#ifdef _FPESIG_ENABLE_
  fpe_catch();
#endif


	int nproc, myproc;
	MPI_Comm_size(MPI_COMM_WORLD, & nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &myproc);
	
	srand48(12345);
  
  if (myproc == 1) {
    gwait = 0;
    while (gwait == 1) {
      std::cerr << " waiting ... zzzz \n";
      sleep(2);
    };
  }
		

#if 0
	{
		fvmhd3d::system s(myproc, nproc);

		s.set_geometry(true);

#if 0
		s.build_mesh(true);
		s.build_mesh(true);
#endif
		
		int nremove = 100;
		const double t0 = mytimer::get_wtime();
		int success = 0;
		for (int i = 0; i < nremove; i++)
		{
			int ix = s.local_n/(nremove - 1) * i;
			const fvmhd3d::Site si(fvmhd3d::vec3(s.ptcl[ix].x, s.ptcl[ix].y, s.ptcl[ix].z), ix, s.ptcl[ix].rmax);
			std::vector< std::pair<int, std::pair<fvmhd3d::real, fvmhd3d::real> > > volume;
			const bool success_flag =  s.remove_site(si, volume);
			if (!success_flag) continue;
			success++;
			fvmhd3d::real v = 0;
			for (int  i1 = 0; i1 < (int)volume.size(); i1++)
			{
				const double dv = volume[i1].second.second - volume[i1].second.first;
				assert(dv > -1.0e-10*volume[i1].second.first);
				const int j = volume[i1].first;
				assert(j != ix);
				if (j < s.local_n)
					s.cells[j].Volume = volume[i1].second.second;
				v += dv;
			}
			if (!(std::abs(s.cells[ix].Volume - v)/s.cells[ix].Volume < 1.0e-13))
			{
				fprintf(stderr, "i= %d [%d] ix= %d vi= %g  vex= %g  dv= %g [ %g ]  %g %g %g\n",
						i, nremove, ix,
						s.cells[ix].Volume, v,  
						(s.cells[ix].Volume - v),
						(s.cells[ix].Volume - v)/s.cells[ix].Volume,
						si.pos.x, si.pos.y, si.pos.z);
			}
			assert(std::abs(s.cells[ix].Volume - v)/s.cells[ix].Volume < 1.0e-10);
		}
		fprintf(stderr, " nremove= %d  success= %d ::  removing done in %g sec \n",
				nremove, success, mytimer::get_wtime() - t0);


	}
#endif

#if 1
	{
		char path[256] = "ZZZZZZ";
		if (argc > 1) {
			sprintf(path, "%s", argv[1]);
		} else {
			if (myproc == 0) {
				fprintf(stderr, " ./main data_path \n");
				exit(-1);
			}
		}

		fvmhd3d::system s(myproc, nproc);

		if (argc > 2)
		{
			if (myproc == 0) 
				fprintf(stderr, " ... Reading snapshot ... \n");
			s.set_geometry(false);
			s.set_problem(false);
			s.read_binary((const char*)argv[2], 1);
		}
		else
		{
			s.set_geometry(true);
			s.set_problem(true);
		}
		
		fvmhd3d::Energy E0(s.get_energy());
    {
      const double volume_exact = s.global_domain_size.x*s.global_domain_size.y*s.global_domain_size.z;
      E0.data[fvmhd3d::Energy::VOLUME] = volume_exact;
    }
    bool first_dE_flag = true;
    if (myproc == 0)
    {
      E0.print_energy  (stderr, " E0: ");
      E0.print_momentum(stderr, " M0: ");
    }

    double t_dump    = s.t_global + (s.t_global > 0 ? s.dt_dump : 0.0);
    t_dump = (int)(t_dump/s.dt_dump)*s.dt_dump;
    double t_restart = s.t_global + s.dt_restart;
    int    i_log     = s.iteration;

    const int niter_max = 5;
    int niter = 0;
    const double t_start = mytimer::get_wtime();


#ifdef _SEQ_WRITE_
    const char lockfn[256] = "write-lock.fvmhd3d";
#endif

    while(niter < niter_max)
    {
      //			niter++;

      bool out = false;
      //
      // ***** writing dump file
      //
      if (s.t_global >= t_dump)
      { 
        if (s.all_active)
        {
          out = true;
          char fn[256];
          const int iout = (t_dump + 0.001*s.dt_dump)/s.dt_dump;
#ifdef _CREATE_FOLDER_
					char filepath[256];
          sprintf(filepath, "%s/iter%.6d", path, iout);
					if (myproc == 0)
						mkdir(filepath, S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IWOTH);
					MPI_Barrier(MPI_COMM_WORLD);

					sprintf(fn, "%s/iter%.6d_proc%.3d.snap", filepath ,iout, myproc);
#else
					sprintf(fn, "%s/iter%.6d_proc%.3d.snap", path, iout, myproc);
#endif
					if (myproc == 0) 
					{
						fprintf(stderr, " *** dumping snapshot %s @ t= %g\n",
								fn, s.t_global);
					}
#ifdef _SEQ_WRITE_
					FILE *fd;
					while ((fd = fopen(lockfn, "r")))
					{
						if (myproc == 0)
							fprintf(stderr, " ... write-lock ... waiting ... zzz\n");
						fclose(fd);
						sleep(1.0);
					}
					MPI_Barrier(MPI_COMM_WORLD);
					if (myproc == 0)
					{
						fd = fopen(lockfn, "w");
						fclose(fd);
						fprintf(stderr, " write-lock removed ... writing ... \n");
					}
					MPI_Barrier(MPI_COMM_WORLD);
					for (int p = 0; p < nproc; p++)
					{
						MPI_Barrier(MPI_COMM_WORLD);
						if (myproc == p)
						{
							fd = fopen(lockfn, "w");
							fclose(fd);
							fprintf(stderr, " proc= %d dumps data into %s @ t = %g ... \n ",
									myproc, fn, s.t_global);
							s.dump_binary(fn);
							fd = fopen(lockfn, "w");
							fclose(fd);
							fprintf(stderr, " proc= %d ... done writing \n", myproc);
						}
					}
					MPI_Barrier(MPI_COMM_WORLD);
					if (myproc == 0)
						remove(lockfn);
#else
					s.dump_binary(fn);
					if (myproc == 0)
						fprintf(stderr, " proc= %d ... done writing \n", myproc);
#endif
					t_dump += s.dt_dump;
				}
			}

#if 0
			MPI_Barrier(MPI_COMM_WORLD);
#endif

			//
			// ***** writing restart file
			//
			if (s.t_global >= t_restart) 
			{
				if (!out && s.all_active)
				{
					char fn[256];
					const int irestart = (t_restart + 0.001*s.dt_restart)/s.dt_restart;
#ifdef _CREATE_FOLDER_
					char filepath[256];
          sprintf(filepath, "%s/restart%.3d", path, irestart%s.n_restart);
					if (myproc == 0)
						mkdir(filepath, S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IWOTH);
					MPI_Barrier(MPI_COMM_WORLD);
					sprintf(fn, "%s/restart%.3d_proc%.3d.snap", filepath ,irestart%s.n_restart, myproc);
#else
					sprintf(fn, "%s/restart%.3d_proc%.3d.snap", path, irestart%s.n_restart, myproc);
#endif

					if (myproc == 0) 
						fprintf(stderr, " *** dumping restart file %s @ t= %g\n",
								fn, s.t_global);
					if (myproc == 0) 
					{
						fprintf(stderr, " *** dumping snapshot %s @ t= %g\n",
								fn, s.t_global);
					}
#ifdef _SEQ_WRITE_
					FILE *fd;
					while ((fd = fopen(lockfn, "r")))
					{
						if (myproc == 0)
							fprintf(stderr, " ... write-lock ... waiting ... zzz\n");
						fclose(fd);
						sleep(1.0);
					}
					MPI_Barrier(MPI_COMM_WORLD);
					if (myproc == 0)
					{
						fd = fopen(lockfn, "w");
						fclose(fd);
						fprintf(stderr, " write-lock removed ... writing ... \n");
					}
					MPI_Barrier(MPI_COMM_WORLD);
					for (int p = 0; p < nproc; p++)
					{
						MPI_Barrier(MPI_COMM_WORLD);
						if (myproc == p)
						{
							fd = fopen(lockfn, "w");
							fclose(fd);
							fprintf(stderr, " proc= %d dumps data into %s @ t = %g ... \n ",
									myproc, fn, s.t_global);
							s.dump_binary(fn);
							fd = fopen(lockfn, "w");
							fclose(fd);
							fprintf(stderr, " proc= %d ... done writing \n", myproc);
						}
					}
					MPI_Barrier(MPI_COMM_WORLD);
					if (myproc == 0)
						remove(lockfn);
#else
					s.dump_binary(fn);
					if (myproc == 0)
						fprintf(stderr, " proc= %d ... done writing \n", myproc);
#endif
				}
				if (s.all_active)
					t_restart += s.dt_restart;
			}

#if 0
			MPI_Barrier(MPI_COMM_WORLD);
#endif

#if 0
#define _PROFILE_EVERY_
			if (myproc == 0)
				fprintf(stderr, "--------------------------------------   new iteration   ------------------------------------------------- \n");
#endif
			MPI_Barrier(MPI_COMM_WORLD);
			const double tbeg = mytimer::get_wtime();
			s.iterate();
			double dt_iter  = mytimer::get_wtime() - tbeg;

			const int n_active = s.nactive_glb;
			double dt_iter_max;
			MPI_Reduce(&dt_iter, &dt_iter_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			unsigned long long nvirtual_glb, boundary_glb;
			MPI_Reduce(&s.virtual_n, &nvirtual_glb,  1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&s.boundary_n, &boundary_glb, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

			if (myproc == 0) 
			{
				fprintf(stderr, " iter= %d : t= %g dt= %g [ %g ] n_act= %d [ %g ] [%lld %lld] running_time= %g h  [ %g sec :: %g cells/s/proc/threads ( %g // %g ) ] \n",
						s.iteration, s.t_global, s.dt_global, s.courant_no,
						n_active, n_active*1.0/s.global_n,
						nvirtual_glb, boundary_glb,
						(mytimer::get_wtime() - t_start)/3600.0, 
						dt_iter_max, 
						s.global_n/nproc/dt_iter, 
						n_active/dt_iter, n_active/dt_iter/nproc);
			}

			if (s.all_active)
			{
#if 1
				const fvmhd3d::Energy E1(s.get_energy());
				if (false && first_dE_flag)
				{
					E0 = s.get_energy();
					E0.print_energy  (stderr, " E0: ");
					E0.print_momentum(stderr, " M0: ");
					first_dE_flag = false;
				}
				const fvmhd3d::Energy dE = (E1 - E0)/E0.abs();
				if (myproc == 0)
				{
					E1.print_energy(stderr, " E1: ");
					dE.print_energy(stderr, " dE: ");
					E1.print_mass(stderr, " M1: ");
					dE.print_mass(stderr, " dM: "); 
					E1.print_momentum(stderr, "Mom1: ");
					dE.print_momentum(stderr, "dMom: ");
				}
#endif
			}

#if 1
			if (s.iteration >= i_log || s.all_active)
#endif
			{
				s.dump_profile_info();
				if (s.iteration >= i_log)
					i_log += s.di_log;
			}

			if (s.t_global >= s.t_end && s.all_active) break;
		}

		if (s.t_global >= t_dump) 
		{ 
			char fn[256];
			const int iout = (t_dump + 0.001*s.dt_dump)/s.dt_dump;
#ifdef _CREATE_FOLDER_
      char filepath[256];
      sprintf(filepath, "%s/iter%.6d", path, iout);
      if (myproc == 0)
        mkdir(filepath, S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IWOTH);
      MPI_Barrier(MPI_COMM_WORLD);

      sprintf(fn, "%s/iter%.6d_proc%.3d.snap", filepath ,iout, myproc);
#else
      sprintf(fn, "%s/iter%.6d_proc%.3d.snap", path, iout, myproc);
#endif
      if (myproc == 0) 
        fprintf(stderr, " *** dumping snapshot %s @ t= %g\n",
            fn, s.t_global);
      s.dump_binary(fn);
      t_dump += s.dt_dump;
    }

  }
#endif



  MPI_Finalize();
  fprintf(stderr, "end-of-program\n");
  return 0;
}
