#include "fvmhd3d.h"

namespace fvmhd3d {

	void system::compute_pvel() 
	{
		if (compute_pvel1()) return;

		const int nactive_site = site_active_list.size();
		for (int isite = 0; isite < nactive_site; isite++)
		{
			const int i = site_active_list[isite];

			Particle &p = ptcl_import[i];

			const Fluid W = U_import[i].to_primitive(ptcl_import[i].volume);

			p.vel = W.get_vel();

      const real B2   = W.get_B().norm2();
      const real pres = compute_pressure(W);
      const real cs2  = (gamma_gas * pres + B2)/W[Fluid::DENS];
      const real vel2 = W.get_vel().norm2();

      const real vabs = std::sqrt(cs2 + vel2);

      const vec3 centroid = cell_list[i].centroid - p.pos;
      const real d = centroid.abs();
      if (d == 0.0) continue;

      const real eta = 0.25f;
      const real ki  = 1.0f;

      const real f1  = 0.9f;
      const real f2  = 1.1f;

      const real R   = std::pow(cell_list[i].Volume * (3.0/(4.0*M_PI)), 1.0/3.0);
      const real fac = d/(eta*R);

      real f;
      if      (fac < f1) f = 0.0;
      else if (fac < f2) f = (d - f1*eta*R)/((f2 - f1)*eta*R);
      else               f = 1.0;

      real tau = d / vabs;

      f *= ki/tau;
      p.vel += centroid*f;

      if (p.is_boundary()) p.vel = 0.0;
    }
  }


}
