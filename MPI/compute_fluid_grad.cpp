#include "fvmhd3d.h"

namespace fvmhd3d {

	void system::compute_fluid_grad(const bool do_ngb)
	{
#if 0
		const int nimport = Uimport.size(); 
		assert(nimport == (int)site_import.size());
		Wextra_import.resize(nimport);

    const int nactive_site = site_active_list.size();
		const int nactive_site_with_ngb = nactive_site + (do_ngb ? site_with_ngb_active_list.size() : 0);
		
    for (int isite = 0; isite < nactive_site_with_ngb; isite++)
		{
			const int i = isite < nactive_site ? 
				site_active_list[isite] : site_with_ngb_active_list[isite - nactive_site];

			const Cell     &ci     = cell_list[i];
			const Fluid_st &Wst_i  = W_st  [i];
			const Fluid    &Wi     = Wst_i.w;
			const vec3     &ipos   = Wst_i.pos;

			vec3 Ji(0.0);     // current,  J = Del x B

			const vec3 Bi(Wi[Fluid::BX], Wi[Fluid::BY], Wi[Fluid::BZ]);
			
			const int nface = ci.faces().size();	
			for (int iface = 0; iface < nface; iface++)
			{
				const Face &face = face_list[ci.faces()[iface]];

				vec3 dri = face.centroid - ipos;
        if      (dri.x >  0.5*global_domain_size.x) dri.x -= global_domain_size.x;
        else if (dri.x < -0.5*global_domain_size.x) dri.x += global_domain_size.x;
        if      (dri.y >  0.5*global_domain_size.y) dri.y -= global_domain_size.y;
        else if (dri.y < -0.5*global_domain_size.y) dri.y += global_domain_size.y;
        if      (dri.z >  0.5*global_domain_size.z) dri.z -= global_domain_size.z;
        else if (dri.z < -0.5*global_domain_size.z) dri.z += global_domain_size.z;
        assert(std::abs(dri.x) < 0.5*global_domain_size.x);
        assert(std::abs(dri.y) < 0.5*global_domain_size.y);
        assert(std::abs(dri.z) < 0.5*global_domain_size.z);
        
        const vec3 centroid = dri;
				const real area     = face.area();
				const vec3 normal   = face.n * ((centroid * face.n < 0.0) ? (-1.0/area) : (1.0/area));

				const real dsh  = centroid * normal;
				assert(dsh > 0.0);

				const real ids  = (dsh != 0.0) ? 0.5/dsh : 0.0;
				const real idsA = area * ids;

				const vec3 drh = normal * (dsh * idsA);
				const vec3 fij = (centroid - drh) * idsA;
				const int j = site_map(face.ngb<false>(i));
				assert(j >= 0);

				const Fluid &Wj = W_st[j].w;
        const vec3 Bj(Wj[Fluid::BX], Wj[Fluid::BY], Wj[Fluid::BZ]);
        Ji += drh.cross(Bj+Bi) + fij.cross(Bj-Bi);
			}

			const real invV = 1.0/cell_list[i].Volume;
			Ji *= ptcl_import[i].etaOhm * invV;
      
      if (do_ngb) W_st[i].J = Ji;
      else     
      {
        const real dt = t_global - std::abs(W_st[i].tlast);
//        assert(dt > 0.0);
//        W_rec[i].J = (Ji - W_st[i].J) * (1.0/dt);
      }

		}
#endif	
	}

}
