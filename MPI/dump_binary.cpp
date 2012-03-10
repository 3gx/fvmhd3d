#include "fvmhd3d.h"

namespace fvmhd3d {

	void system::dump_binary(const char *filename) 
	{

#if 0
		return;
#endif

		FILE *fout; 
		if (!(fout = fopen(filename, "w"))) {
			std::cerr << "Cannot open file " << filename << std::endl;
			exit(-1);
		}

		int ival;
		float fval;

#define fdump(x) { fval = x; fwrite(&fval, sizeof(float), 1, fout); }
#define idump(x) { ival = x; fwrite(&ival, sizeof(int),   1, fout); }

		idump(20*4);
		idump(myproc);
		idump(nproc);
		idump(nproc);
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
			const Particle &pi = ptcl_local[i];
			assert(pi.volume > 0.0);
			const Fluid mi = U_local[i].to_primitive(pi.volume);

			const real L = std::pow(pi.volume*3.0/4.0/M_PI, 1.0/3.0);

			idump(26*4);
			ival = pi.idx; idump(ival);
			vec3 pos = pi.orig_pos;
			fdump(pos.x);
			fdump(pos.y);
			fdump(pos.z);
			fdump(pi.orig_vel.x);
			fdump(pi.orig_vel.y);
			fdump(pi.orig_vel.z);
			fdump(mi[Fluid::DENS]);
			fdump(mi[Fluid::ETHM]);
			fdump(compute_pressure(mi));
			fdump(pi.rmax);      //   fdump(    (sqr(mi[Fluid::BX  ]) + sqr(mi[Fluid::BY  ]) + sqr(mi[Fluid::BZ  ]))*0.5f);
			idump(pi.boundary);  //		fdump(sqrt(sqr(mi[Fluid::VELX]) + sqr(mi[Fluid::VELY]) + sqr(mi[Fluid::VELZ]))*0.5f);
			fdump(mi[Fluid::VELX]);
			fdump(mi[Fluid::VELY]);
			fdump(mi[Fluid::VELZ]);
			fdump(mi[Fluid::BX]);
			fdump(mi[Fluid::BY]);
			fdump(mi[Fluid::BZ]);
			fdump(scheduler.get_dt(ptcl_local[i].rung));
			fdump(pi.volume);
			fdump(mi[Fluid::PSI]);
			fdump(L*divBi[i]);
			fdump(mi[Fluid::ENTR]);
			fdump(Wextra_local[i].J.x);
			fdump(Wextra_local[i].J.y);
			fdump(Wextra_local[i].J.z);
			idump(26*4);

		}

		fclose(fout);

	}

}
