#include "fvmhd3d.h"

#include "peano.h"


namespace fvmhd3d {

	template< class T >
		void reorder(const std::vector<int> &order, std::vector<T> &v)
		{
			const std::vector<T> v_orig(v);
			const int order_size = order.size();
			for (int i = 0; i < order_size; i++)
				v[i] = v_orig[order[i]];
		}

	void system::sort_local_data() {
		if (ptcl_local.empty()) return;
		static std::vector<int> order;
		order.resize(local_n);
		order.clear();

		{
			static std::vector<peano_hilbert::peano_struct> keys;
			keys.resize(local_n);
			keys.clear();

			const vec3 domain_len(global_domain.hsize() * 2.0);
			const float size = std::max(domain_len.x,
					std::max(domain_len.y, domain_len.z));

			const int domain_fac = 1.0f / size * (((peano_hilbert::peanokey)1) << (BITS_PER_DIMENSION));

			const double xmin = global_domain.get_rmin().x;
			const double ymin = global_domain.get_rmin().y;
			const double zmin = global_domain.get_rmin().z;

			for (int i = 0; i < (int)local_n; i++) {
				const int x = (int)((ptcl_local[i].pos.x - xmin) * domain_fac);
				const int y = (int)((ptcl_local[i].pos.y - ymin) * domain_fac);
				const int z = (int)((ptcl_local[i].pos.z - zmin) * domain_fac);
				keys[i].idx = i;
				keys[i].key = peano_hilbert::peano_hilbert_key(x, y, z, BITS_PER_DIMENSION);
			}

			std::sort(keys.begin(), keys.end(), peano_hilbert::cmp_peanokey_index());

			for (int i = 0; i < (int)local_n; i++)
				order[i] = keys[i].idx;
		}

		reorder (order, ptcl_local);
		reorder (order,    U_local);
		reorder (order,   dU_local);
		reorder (order, Wrec_local);
	}

}
