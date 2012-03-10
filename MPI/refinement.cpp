#include "fvmhd3d.h"

#define SMALLDIFF (1.0e-10)

namespace fvmhd3d
{

	void system::do_derefinement()
	{
#if 0
		const int nactive_ptcl = active_ptcl.size();
		int nderefine = 0;
		int nsuccess  = 0;
			
		std::vector<Fluid> Uj       (2*import_n, 0.0);
		std::vector<real>  gradPsi_j(4*import_n, 0.0);
		std::vector<int >  derefined;

		for (int iptcl = 0; iptcl < nactive_ptcl; iptcl++)
		{
			const int i = active_ptcl[iptcl];
			if (!ptcl[i].is_derefine()) continue;
			nderefine++;

			std::vector< std::pair<int, std::pair<fvmhd3d::real, fvmhd3d::real> > > volume;
			const bool success_flag = remove_site(Site(ptcl[i].pos(), i, ptcl[i].rmax), volume);
#if 0
			fprintf(stderr, " deref= %d   i= %d  success= %s \n",
					nderefine, i, success_flag ? "true" : "false");
#endif

			if (!success_flag) // || nsuccess > 0)
			{
				ptcl[i].unset_derefine();
				continue;
			}
			nsuccess++;
			derefined.push_back(i);

			real Vi = 0.0;
			for (int jptcl = 0; jptcl < (const int)volume.size(); jptcl++)
			{
				const real vold = volume[jptcl].second.first;
				const real vnew = volume[jptcl].second.second;
				const real dv   = vnew - vold;
				Vi += dv;
			}
#if 0
			if (!(std::abs(Vi - cells[i].Volume) < SMALLDIFF * Vi))
			{
				fprintf(stderr, " Vi= %g   cells[i].Volume= %g   diff= %g  %d %d \n",
						Vi, cells[i].Volume, Vi - cells[i].Volume, local_n, global_n);
			}
			assert(std::abs(Vi - cells[i].Volume) < SMALLDIFF * Vi);
#endif
			assert(Vi > 0.0);
			const real invVi = 1.0/Vi;

			for (int jptcl = 0; jptcl < (const int)volume.size(); jptcl++)
			{
				const int  j    = site_map(volume[jptcl].first);
				const real vold = volume[jptcl].second.first;
				const real vnew = volume[jptcl].second.second;
				const real dv   = vnew - vold;
				assert(vold > 0.0);
				if (!(dv > -SMALLDIFF*vold))
				{
					fprintf(stderr , "i= %d j= %d [ %d ] :: vnew= %g   vold= %g   dv= %g \n",
							i, j, jptcl,
							vnew, vold, dv);
				}
				assert(dv > -SMALLDIFF*vold);
				const real dv_over_Vi = dv * invVi;
				if (j >= local_n)
				{
					const int j1 = j - local_n;
					assert(false);
					for (int k = 0; k < Fluid::NFLUID; k++)
					{
						Uj[(j1 << 1) + 0][k] +=  U[i][k] * dv_over_Vi;
						Uj[(j1 << 1) + 1][k] += dU[i][k] * dv_over_Vi;
					}

					for (int k = 0; k < 4; k++)
						gradPsi_j[(j1 << 2) + k] += gradPsi[(i << 2) + k] * dv_over_Vi;
				}
				else
				{
					for (int k = 0; k < Fluid::NFLUID; k++)
					{
						U [j][k] +=  U[i][k] * dv_over_Vi;
						dU[j][k] += dU[i][k] * dv_over_Vi;
					}

					for (int k = 0; k < 4; k++)
						gradPsi[(j << 2) + k] += gradPsi[(i << 2) + k] * dv_over_Vi;
				}
			}
		}

#if 0
		std::vector<Fluid> Fluid_send[NMAXPROC];
		std::vector<Fluid> Fluid_recv[NMAXPROC];
		std::vector<real> gradPsi_send[NMAXPROC];
		std::vector<real> gradPsi_recv[NMAXPROC];

		for (int i = 0; i < (const int)return_site_list.size(); i++)
		{
			const std::pair<int, int> p = return_site_list[i];
			const int proc = p.first;
			const int jidx = p.second;

			assert(jidx >= local_n);
			const int j2 = (jidx - local_n) << 1;
			Fluid_send[proc].push_back(Uj[j2 + 0]);
			Fluid_send[proc].push_back(Uj[j2 + 1]);

			const int j4 = (jidx-local_n) << 2;
			gradPsi_send[proc].push_back(gradPsi_j[j4 + 0]);
			gradPsi_send[proc].push_back(gradPsi_j[j4 + 1]);
			gradPsi_send[proc].push_back(gradPsi_j[j4 + 2]);
			gradPsi_send[proc].push_back(gradPsi_j[j4 + 3]);
		}

		myMPI::all2all(Fluid_send, Fluid_recv, myproc, nproc, mpi_debug_flag);
		myMPI::all2all(gradPsi_send, gradPsi_recv, myproc, nproc, mpi_debug_flag);

		int iloc = 0;
		for (int p = 0; p < nproc; p++) 
			for (size_t q = 0; q < Fluid_recv[p].size(); q += 2)
			{
				assert(p != myproc);
				const int idx = return_site_map[iloc++];
				assert(idx < local_n);
				if (!ptcl[idx].is_active())
				{
					for (int k = 0; k < Fluid::NFLUID; k++)
					{
						U [idx][k] += Fluid_recv[p][q + 0][k];
						dU[idx][k] += Fluid_recv[p][q + 1][k];
					}
					gradPsi[(idx<<2) + 0] += gradPsi_recv[p][(q<<1) + 0];
					gradPsi[(idx<<2) + 1] += gradPsi_recv[p][(q<<1) + 1];
					gradPsi[(idx<<2) + 2] += gradPsi_recv[p][(q<<1) + 2];
					gradPsi[(idx<<2) + 3] += gradPsi_recv[p][(q<<1) + 3];
				}
			}

		assert(iloc == (int)export_site_list.size());
#endif

		fprintf(stderr, " nderefine= %d   nsuccess= %d  [ %g ]\n",
				nderefine, nsuccess, nderefine > 0 ? (real)nsuccess/nderefine : 0.0);
		assert(nsuccess == (int)derefined.size());

		for (int i = 0; i < nsuccess; i++)
		{
			ptcl[derefined[i]].set_virtual();
//			scheduler.remove<true>((int)ptcl[derefined[i]].rung[0], derefined[i]);
			assert(ptcl[derefined[i]].is_virtual());
			virtual_n++;
		}

#if 0
		if (derefined.size() == 0) return;

		fprintf(stderr , "local_n_old= %d ", local_n);

		int idx = 0;
		while (idx++ < local_n)
		{
			if (!ptcl[idx-1].is_derefine()) continue;

			idx--;

			assert(ptcl[idx].is_active());
			local_n--;
			std::swap(ptcl [idx], ptcl [local_n]);
			std::swap(ptcl_ppos[idx], ptcl_ppos[local_n]);
			std::swap( U   [idx],  U   [local_n]);
			std::swap(dU   [idx], dU   [local_n]);
			std::swap(Wgrad[idx], Wgrad[local_n]);
			std::swap(gradPsi[(idx<<2) + 0], gradPsi[(local_n<<2)+0]);
			std::swap(gradPsi[(idx<<2) + 1], gradPsi[(local_n<<2)+1]);
			std::swap(gradPsi[(idx<<2) + 2], gradPsi[(local_n<<2)+2]);
			std::swap(gradPsi[(idx<<2) + 3], gradPsi[(local_n<<2)+3]);
		}
		fprintf(stderr , "local_n_new= %d ", local_n);
			
		const int n_glob0 = global_n;
		MPI_Allreduce(&local_n, &global_n, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		if (myproc == 0)
			fprintf(stderr, " myproc= %d :  deref %d cells :: global_n_old= %d  global_n= %d \n",
					myproc, n_glob0 - global_n, n_glob0,  global_n);
#endif
#endif
	}

	void system::do_refinement()
	{
#if 0
		const int nactive_ptcl = active_ptcl.size();
		int nrefine = 0;
		int nsuccess  = 0;

		std::vector<vec3> new_sites_pos;
		std::vector<real> new_sites_rmax;
		for (int iptcl = 0; iptcl < nactive_ptcl; iptcl++)
		{
			const int i = active_ptcl[iptcl];
			if (!ptcl[i].is_refine()) continue;
			nrefine++;

			const Cell &ci = cells[i];
			const int nface = ci.faces.size();	
			const size_t size0 = new_sites_pos.size();
			for (int iface = 0; iface < nface; iface++)
			{
				Face &face = faces[ci.faces[iface]];
				if (face.s1 < 0) continue;
				face.s1 = int_map(face.s1);
				const int j = site_map(face.ngb<false>(i));
				if (j >= local_n)
				{
					new_sites_pos.resize(size0);
					new_sites_rmax.resize(size0);
					break;
				}
				else
				{
					new_sites_pos.push_back(face.centroid);
					new_sites_rmax.push_back(std::max(ptcl[i].rmax, ptcl[j].rmax));
				}
			}
			if (size0 != new_sites_pos.size())
				nsuccess++;
		}

		std::vector<Fluid> Uj       (2*import_n, 0.0);
		std::vector<real>  gradPsi_j(4*import_n, 0.0);

		const int n_new_sites = new_sites_pos.size();
		for (int i = 0; i < n_new_sites; i++)
		{
			std::vector< std::pair<int, std::pair<fvmhd3d::real, fvmhd3d::real> > > volume;
			const bool success_flag = insert_site(Site(new_sites_pos[i], -1-i, new_sites_rmax[i]), volume);
			assert(success_flag);

			for (int jsite = 0; jsite < (const int)volume.size(); jsite++)
			{
			}

		}
#endif
	}


}
