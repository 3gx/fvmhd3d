#include "fvmhd3d.h"

namespace myMPI {

	unsigned long long data_inout;
	double             data_inout_time;
		
	static MPI_Datatype MPI_PFLOAT3   = 0;
	static MPI_Datatype MPI_PARTICLE  = 0;
	static MPI_Datatype MPI_PBOUNDARY = 0;
	static MPI_Datatype MPI_BOUNDARY  = 0;
	static MPI_Datatype MPI_PTCLFLUID = 0;
	static MPI_Datatype MPI_PTCLFLUIDLITE = 0;
	static MPI_Datatype MPI_FLUID     = 0;
	static MPI_Datatype MPI_FLUIDST   = 0;
	static MPI_Datatype MPI_FLUIDREC  = 0;
	
	template <> MPI_Datatype datatype<int   >() {return MPI_INT;    }
	template <> MPI_Datatype datatype<float >() {return MPI_FLOAT;  }
	template <> MPI_Datatype datatype<double>() {return MPI_DOUBLE; }

	template <> MPI_Datatype datatype<pBoundary>() 
	{
		if (MPI_PBOUNDARY) return MPI_PBOUNDARY;
		else {
			int ss = sizeof(pBoundary) / sizeof(float);
			assert(0 == sizeof(pBoundary) % sizeof(float));
			MPI_Type_contiguous(ss, MPI_FLOAT, &MPI_PBOUNDARY);
			MPI_Type_commit(&MPI_PBOUNDARY);
			return MPI_PBOUNDARY;
		}
	}
	
	template <> MPI_Datatype datatype<fvmhd3d::boundary>() 
	{
		if (MPI_BOUNDARY) return MPI_BOUNDARY;
		else {
			int ss = sizeof(fvmhd3d::boundary) / sizeof(float);
			assert(0 == sizeof(fvmhd3d::boundary) % sizeof(float));
			MPI_Type_contiguous(ss, MPI_FLOAT, &MPI_BOUNDARY);
			MPI_Type_commit(&MPI_BOUNDARY);
			return MPI_BOUNDARY;
		}
	}

	template <> MPI_Datatype datatype<fvmhd3d::Particle> () 
	{
		if (MPI_PARTICLE) return MPI_PARTICLE;
		else {
			int ss = sizeof(fvmhd3d::Particle) / sizeof(float);
			assert(0 == sizeof(fvmhd3d::Particle) % sizeof(float));
			MPI_Type_contiguous(ss, MPI_FLOAT, &MPI_PARTICLE);
			MPI_Type_commit(&MPI_PARTICLE);
			return MPI_PARTICLE;
		}
	}

	template <> MPI_Datatype datatype<fvmhd3d::ParticleFluidStruct> () 
	{
		if (MPI_PTCLFLUID) return MPI_PTCLFLUID;
		else {
			int ss = sizeof(fvmhd3d::ParticleFluidStruct) / sizeof(float);
			assert(0 == sizeof(fvmhd3d::ParticleFluidStruct) % sizeof(float));
			MPI_Type_contiguous(ss, MPI_FLOAT, &MPI_PTCLFLUID);
			MPI_Type_commit(&MPI_PTCLFLUID);
			return MPI_PTCLFLUID;
		}
	}
	
  template <> MPI_Datatype datatype<fvmhd3d::ParticleFluidStructLite> () 
	{
		if (MPI_PTCLFLUIDLITE) return MPI_PTCLFLUIDLITE;
		else {
			int ss = sizeof(fvmhd3d::ParticleFluidStructLite) / sizeof(float);
			assert(0 == sizeof(fvmhd3d::ParticleFluidStructLite) % sizeof(float));
			MPI_Type_contiguous(ss, MPI_FLOAT, &MPI_PTCLFLUIDLITE);
			MPI_Type_commit(&MPI_PTCLFLUIDLITE);
			return MPI_PTCLFLUIDLITE;
		}
	}
  
  template <> MPI_Datatype datatype<fvmhd3d::Fluid> () 
	{
		if (MPI_FLUID) return MPI_FLUID;
		else {
			int ss = sizeof(fvmhd3d::Fluid) / sizeof(float);
			assert(0 == sizeof(fvmhd3d::Fluid) % sizeof(float));
			MPI_Type_contiguous(ss, MPI_FLOAT, &MPI_FLUID);
			MPI_Type_commit(&MPI_FLUID);
			return MPI_FLUID;
		}
	}
	
  template <> MPI_Datatype datatype<fvmhd3d::Fluid_st> () 
	{
		if (MPI_FLUIDST) return MPI_FLUIDST;
		else {
			int ss = sizeof(fvmhd3d::Fluid_st) / sizeof(float);
			assert(0 == sizeof(fvmhd3d::Fluid_st) % sizeof(float));
			MPI_Type_contiguous(ss, MPI_FLOAT, &MPI_FLUIDST);
			MPI_Type_commit(&MPI_FLUIDST);
			return MPI_FLUIDST;
		}
	}
  
  template <> MPI_Datatype datatype<fvmhd3d::Fluid_rec> () 
	{
		if (MPI_FLUIDREC) return MPI_FLUIDREC;
		else {
			int ss = sizeof(fvmhd3d::Fluid_rec) / sizeof(float);
			assert(0 == sizeof(fvmhd3d::Fluid_rec) % sizeof(float));
			MPI_Type_contiguous(ss, MPI_FLOAT, &MPI_FLUIDREC);
			MPI_Type_commit(&MPI_FLUIDREC);
			return MPI_FLUIDREC;
		}
	}

	void free_type()
	{
		if (MPI_PFLOAT3  ) MPI_Type_free(&MPI_PFLOAT3);
		if (MPI_PARTICLE ) MPI_Type_free(&MPI_PARTICLE);
		if (MPI_PBOUNDARY) MPI_Type_free(&MPI_PBOUNDARY);
		if (MPI_BOUNDARY ) MPI_Type_free(&MPI_BOUNDARY);
		if (MPI_PTCLFLUID) MPI_Type_free(&MPI_PTCLFLUID);
		if (MPI_PTCLFLUIDLITE) MPI_Type_free(&MPI_PTCLFLUIDLITE);
		if (MPI_FLUID) MPI_Type_free(&MPI_FLUID);
		if (MPI_FLUIDST) MPI_Type_free(&MPI_FLUIDST);
		if (MPI_FLUIDREC) MPI_Type_free(&MPI_FLUIDREC);
	}

}
