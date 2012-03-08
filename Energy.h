#ifndef __ENERGY_H__
#define __ENERGY_H__

struct Energy 
{
	enum 
  {
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
  void pup(PUP::er &p)
  {
    PUParray(p, data, NVAR);
  }
  Energy() {};
	Energy(
			const std::vector<MeshPoint> &mesh_pnts,
			const std::vector<Fluid>     &   U_list)
//			const std::vector<real>     &gpot,
//			const std::vector<real>     &gpot_body)
	{

    const int n = mesh_pnts.size();
		
    for (int i = 0; i < NVAR; i++)
			data[i] = 0;

		for (int i = 0; i < n; i++) 
    {
			const MeshPoint &pi = mesh_pnts[i];
      const Fluid     &Ui =    U_list[i];
#if 0
			if (pi.is_virtual()) 
				continue;
			if (pi.bnd() > 0)
  			continue;
#endif

			const double     V = pi.Volume;
			const Fluid     Wi = Ui.to_primitive(V);

			const double M = Wi[Fluid::DENS] * V;
			data[MASS ] += Ui [Fluid::MASS];	
			data[VOLUME] += V;
			data[EKINX] += M * sqr(Wi[Fluid::VELX]) * 0.5;
			data[EKINY] += M * sqr(Wi[Fluid::VELY]) * 0.5;
			data[EKINZ] += M * sqr(Wi[Fluid::VELZ]) * 0.5;
			data[EMAGX] += V * sqr(Wi[Fluid::BX  ]) * 0.5;
			data[EMAGY] += V * sqr(Wi[Fluid::BY  ]) * 0.5;
			data[EMAGZ] += V * sqr(Wi[Fluid::BZ  ]) * 0.5;
			data[ETHM ] += Wi[Fluid::ETHM] * V;
			data[EPOT ] += 0.0; //M * (0.5 * gpot[i] + gpot_body[i]);

			data[MOMX] += M * Wi[Fluid::VELX];
			data[MOMY] += M * Wi[Fluid::VELY];
			data[MOMZ] += M * Wi[Fluid::VELZ];

			data[ANMX] += M * (Wi[Fluid::VELY]*pi.pos.z - Wi[Fluid::VELZ]*pi.pos.y);
			data[ANMY] += M * (Wi[Fluid::VELZ]*pi.pos.x - Wi[Fluid::VELX]*pi.pos.z);
			data[ANMZ] += M * (Wi[Fluid::VELX]*pi.pos.y - Wi[Fluid::VELY]*pi.pos.z);
		}

		data[EKIN] = data[EKINX] + data[EKINY] + data[EKINZ];
		data[EMAG] = data[EMAGX] + data[EMAGY] + data[EMAGZ];
		data[ETOT] = data[ ETHM] + data[ EKIN] + data[ EMAG] + data[EPOT];
		data[MOM ] = data[ MOMX] + data[ MOMY] + data[ MOMZ];
		data[ANM ] = data[ ANMX] + data[ ANMY] + data[ ANMZ];
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

  void print_energy(const char *prefix= "##") const 
  {
    CkPrintf("%s : Ekin= %g  Ethm= %g  Emag= %g  Epot= %g  :: Etot= %g \n",
        prefix, data[EKIN], data[ETHM], data[EMAG], data[EPOT], data[ETOT]);
  }
  void print_mass(const char *prefix= "##") const 
  {
    CkPrintf("%s : Mass= %g  Volume= %g \n", prefix, data[MASS], data[VOLUME]);
  }
  void print_momentum(const char *prefix= "##") const 
  {
    CkPrintf("%s : mom= %g %g %g; %g \n", prefix, 
        data[MOMX], data[MOMY], data[MOMZ], data[MOM]);
  }
};


#endif //  __ENERGY_H__

