#ifndef __ENERGY_H__
#define __ENERGY_H__

struct Energy {
	enum {
		MASS  = 0,
		ETHM, EPOT , ETOT,
		EKINX, EKINY, EKINZ, EKIN,
		EMAGX, EMAGY, EMAGZ, EMAG,
		MOMX, MOMY, MOMZ, MOM,
		ANMX, ANMY, ANMZ, ANM,
		VOLUME,
		NVAR
	};
	double data[NVAR];
	Energy(const double _data[NVAR]) {
		for (int i = 0; i < NVAR; i++)
			data[i] = _data[i];
	}
	Energy(
			const int n,
			const std::vector<Particle> &ptcl,
			const std::vector<Fluid>    &U)
//			const std::vector<real>     &gpot,
//			const std::vector<real>     &gpot_body)
	{

		double ldata[NVAR];
		for (int i = 0; i < NVAR; i++)
			ldata[i] = 0;

		for (int i = 0; i < n; i++) {
			const Particle &pi = ptcl[i];
#if 0
			if (pi.is_virtual()) 
				continue;
			if (pi.bnd() > 0)
  			continue;
#endif

			const double     V = ptcl[i].volume;
			const Fluid     Wi = U[i].to_primitive(V);

			const double M = Wi[Fluid::DENS] * V;
			ldata[MASS ] += U[i][Fluid::MASS];	
			ldata[VOLUME] += V;
			ldata[EKINX] += M * sqr(Wi[Fluid::VELX]) * 0.5;
			ldata[EKINY] += M * sqr(Wi[Fluid::VELY]) * 0.5;
			ldata[EKINZ] += M * sqr(Wi[Fluid::VELZ]) * 0.5;
			ldata[EMAGX] += V * sqr(Wi[Fluid::BX  ]) * 0.5;
			ldata[EMAGY] += V * sqr(Wi[Fluid::BY  ]) * 0.5;
			ldata[EMAGZ] += V * sqr(Wi[Fluid::BZ  ]) * 0.5;
			ldata[ETHM ] += Wi[Fluid::ETHM] * V;
			ldata[EPOT ] += 0.0; //M * (0.5 * gpot[i] + gpot_body[i]);

			ldata[MOMX] += M * Wi[Fluid::VELX];
			ldata[MOMY] += M * Wi[Fluid::VELY];
			ldata[MOMZ] += M * Wi[Fluid::VELZ];

			ldata[ANMX] += M * (Wi[Fluid::VELY]*pi.pos.z - Wi[Fluid::VELZ]*pi.pos.y);
			ldata[ANMY] += M * (Wi[Fluid::VELZ]*pi.pos.x - Wi[Fluid::VELX]*pi.pos.z);
			ldata[ANMZ] += M * (Wi[Fluid::VELX]*pi.pos.y - Wi[Fluid::VELY]*pi.pos.z);
		}

		ldata[EKIN] = ldata[EKINX] + ldata[EKINY] + ldata[EKINZ];
		ldata[EMAG] = ldata[EMAGX] + ldata[EMAGY] + ldata[EMAGZ];
		ldata[ETOT] = ldata[ETHM] + ldata[EKIN] + ldata[EMAG] + ldata[EPOT];
		ldata[MOM ] = ldata[MOMX] + ldata[MOMY] + ldata[MOMZ];
		ldata[ANM ] = ldata[ANMX] + ldata[ANMY] + ldata[ANMZ];

		MPI_Allreduce(ldata, data, NVAR, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	}
	const double operator[] (const int i) const {return data[i];}

	const Energy abs() const 
	{
		double data1[NVAR];
		for (int i = 0; i < NVAR; i++) {
			data1[i] = std::abs(data[i]);
		}	
		return Energy(data1);
	}
	const Energy operator-(const Energy &rhs) const {
		double data1[NVAR];
		for (int i = 0; i < NVAR; i++) {
			data1[i] = data[i] - rhs[i];
		}	
		return Energy(data1);
	}
	const Energy operator/(const Energy &rhs) const {
		double data1[NVAR];
		for (int i = 0; i < NVAR; i++) {
			data1[i] = (rhs[i] != 0) ? data[i] / rhs[i] : data[i];
		}	
		return Energy(data1);
	}


	void print_energy(FILE *fout, const char *prefix= "##") const {
		fprintf(fout, "%s : Ekin= %g  Ethm= %g  Emag= %g  Epot= %g  :: Etot= %g \n",
				prefix, data[EKIN], data[ETHM], data[EMAG], data[EPOT], data[ETOT]);
	}
	void print_mass(FILE *fout, const char *prefix= "##") const {
		fprintf(fout, "%s : Mass= %g  Volume= %g \n", prefix, data[MASS], data[VOLUME]);
	}
	void print_momentum(FILE *fout, const char *prefix= "##") const {
		fprintf(fout, "%s : mom= %g %g %g; %g \n", prefix, 
				data[MOMX], data[MOMY], data[MOMZ], data[MOM]);
	}
};


#endif //  __ENERGY_H__

