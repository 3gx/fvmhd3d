#include "venom3d.h"
#include "sol_kep/kep_drift.h"

/** ENABLE FOR CLOCKCYCLE TIMINGS ***/
// #define TICKPROFILE

#ifdef TICKPROFILE
#define RDTSCL0(x) rdtscl(x)
#else
#define RDTSCL0(x)
#endif

namespace venom3d {

	void system::compute_update() 
	{

		for (int i = local_n; i < local_n + import_n; i++)
			dU[i] = Fluid(0.0);

		unsigned long long ninteraction = 0;

		unsigned long t_outer_beg;	RDTSCL0(&t_outer_beg);

		unsigned long delta_outer = 0;
		unsigned long delta_inner = 0;
		unsigned long delta_read_simd = 0;
		unsigned long delta_flux = 0;

		const vreal tcur(t_global);

		const int nt = active_tiles.size();
		for (int at = 0; at < nt; at++) 
		{
			const std::vector<int> &ptcl_list = active_ptcl[at];
			const int np = ptcl_list.size();
			for (int idx = 0; idx < np; idx++)
			{
				const int i = ptcl_list[idx];
				
				const real tau = t_global - ptcl[i].tlast;
				assert(tau > 0.0);
				const vreal v_tau(tau);

				const Cell                 &ci = cells[i];
				const Fluid_grad_simd  gradWi(Wgrad[i]  );
				const Particle_simd        pi(ptcl [i]  );
				const vreal                ti(pi.tlast);

#if 0
				const simd::vfloat gaccxi(gaccx[i]);
				const simd::vfloat gaccyi(gaccy[i]);
				const simd::vfloat gacczi(gaccz[i]);
#endif

				Fluid_simd v_dU    (0.0);
				vreal v_dEtot(0.0);
				vreal      v_divB(0.0);

				vreal v_gradPsi_x(0.0);
				vreal v_gradPsi_y(0.0);
				vreal v_gradPsi_z(0.0);

//				simd::vfloat v_dVdt(0.0f);

				const int jbeg = ci.jbeg;
				const int jend = ci.jend;
				assert((jend - jbeg) % simd::vlen == 0);

				asm("#HOTSPOT1-BEG");

				simd::vfloat virel2max(0.0f);
				unsigned long t_inner_beg; RDTSCL0(&t_inner_beg);
				for (int j = jbeg; j < jend; j += simd::vlen) 
				{ 
					ninteraction++;
					asm("#LDST1-BEG");
					unsigned long t_read_simd_beg; RDTSCL0(&t_read_simd_beg);
					const Face_simd  fij(&faces[j]);
					const Fluid_grad_simd gradWj(&Wgrad[0], fij.idx);
					const Particle_simd   pj(&ptcl[0], fij.idx);
					const vreal           tj(pj.tlast);

					unsigned long t_read_simd_end; RDTSCL0(&t_read_simd_end);
					delta_read_simd += t_read_simd_end - t_read_simd_beg;

					asm("#LDST1-END");
					vreal psi_ij(0.0), Bn_ij(0.0);
					Fluid_simd fluxij;
					asm("#COMPUTEFLUX-BEG");
					unsigned long t_flux_beg; RDTSCL0(&t_flux_beg);
					simd::vfloat dti = (tcur - ti).to_vfloat();
					simd::vfloat dtj = (tcur - tj).to_vfloat();
					simd::vfloat dt  = simd::min(dti, dtj);
					dti += dt * (-0.5f);
					dtj += dt * (-0.5f);


					const simd::vfloat dsl = simd::vfloat(2.0f)*(fij.cx*fij.nx + fij.cy*fij.ny + fij.cz*fij.nz);
					const simd::vfloat dx  = fij.nx * dsl;
					const simd::vfloat dy  = fij.ny * dsl;
					const simd::vfloat dz  = fij.nz * dsl;
					//spflop 9
					const simd::vfloat dx_i = fij.cx;
					const simd::vfloat dy_i = fij.cy;
					const simd::vfloat dz_i = fij.cz;

					const simd::vfloat dx_j = fij.cx - dx;
					const simd::vfloat dy_j = fij.cy - dy;
					const simd::vfloat dz_j = fij.cz - dz;

					Fluid_simd Wi, Wj;
					for (int k = 0; k < Fluid::NFLUID; k++) 
					{
						const simd::vfloat fWi = gradWi.m[k].to_vfloat();
						const simd::vfloat fWj = gradWj.m[k].to_vfloat();
						const simd::vfloat dWi  = (gradWi.x[k]*dx_i + gradWi.y[k]*dy_i + gradWi.z[k]*dz_i);
						const simd::vfloat dWj  = (gradWj.x[k]*dx_j + gradWj.y[k]*dy_j + gradWj.z[k]*dz_j);
						const simd::vfloat adWi = simd::abs(dWi);
						const simd::vfloat adWj = simd::abs(dWj);
						simd::vfloat dWij = simd::abs(fWi - fWj);
						dWij *= (adWi + adWj <= dWij).select(1.0f, 0.5f);
						simd::vfloat fi = (adWi <= dWij).select(1.0f, 0.0f);
						simd::vfloat fj = (adWj <= dWij).select(1.0f, 0.0f);
					
						if (k == Fluid::DENS || k == Fluid::ETHM)	
						{
							fi = (adWi <= simd::abs(fWi) * 0.5f).select(fi, 0.0f);
							fj = (adWj <= simd::abs(fWj) * 0.5f).select(fj, 0.0f);
						}
            Wi[k] = gradWi.m[k] + vreal(fi*(dWi + gradWi.t[k]*dti));
            Wj[k] = gradWj.m[k] + vreal(fj*(dWj + gradWj.t[k]*dtj));

#if 0
						assert((Wi[Fluid::DENS] <= vreal(0.0)).select(-1.0, 0.0).sum() >= 0.0);
						assert((Wi[Fluid::ETHM] <= vreal(0.0)).select(-1.0, 0.0).sum() >= 0.0);
						assert((Wj[Fluid::DENS] <= vreal(0.0)).select(-1.0, 0.0).sum() >= 0.0);
						assert((Wj[Fluid::ETHM] <= vreal(0.0)).select(-1.0, 0.0).sum() >= 0.0);

#else
						const vreal qmin = simd::min(gradWi.m[k], gradWj.m[k]);
						const vreal qmax = simd::max(gradWi.m[k], gradWj.m[k]);
						Wi[k] = simd::min(qmax, simd::max(qmin, Wi[k]));
						Wj[k] = simd::min(qmax, simd::max(qmin, Wj[k]));
#endif
					}

#if 1
					extra_ptcl_update(Wi, gradWi.m, pi, i,       dti);
					extra_ptcl_update(Wj, gradWj.m, pj, fij.idx, dtj);
#endif

					const vreal ivx = pi.vx;
					const vreal ivy = pi.vy;
					const vreal ivz = pi.vz;
					const vreal jvx = pj.vx;
					const vreal jvy = pj.vy;
					const vreal jvz = pj.vz;

					//dpflop 24 + 3*20 + 1*20
					const simd::vfloat vijx = (jvx + ivx).to_vfloat();
					const simd::vfloat vijy = (jvy + ivy).to_vfloat();
					const simd::vfloat vijz = (jvz + ivz).to_vfloat();
					const simd::vfloat  dvx = (jvx - ivx).to_vfloat();
					const simd::vfloat  dvy = (jvy - ivy).to_vfloat();
					const simd::vfloat  dvz = (jvz - ivz).to_vfloat();
					//dpflop 6
					const simd::vfloat dsh = fij.cx*fij.nx + fij.cy*fij.ny + fij.cz*fij.nz;

					const simd::vfloat dxh = fij.nx * dsh;
					const simd::vfloat dyh = fij.ny * dsh;
					const simd::vfloat dzh = fij.nz * dsh;

					const simd::vfloat fx  = fij.cx - dxh;
					const simd::vfloat fy  = fij.cy - dyh;
					const simd::vfloat fz  = fij.cz - dzh;
					//spflop 7
					const simd::vfloat dvdij = dvx*fx + dvy*fy + dvz*fz;
					const simd::vfloat  ids2 = -simd::vfloat(0.5f)*dvdij * simd::rcp_safe(sqr(dsh));
					const vreal  wx = simd::vfloat(0.5f)*vijx + ids2 * dxh;
					const vreal  wy = simd::vfloat(0.5f)*vijy + ids2 * dyh;
					const vreal  wz = simd::vfloat(0.5f)*vijz + ids2 * dzh;

					compute_flux(wx, wy, wz, Wi, Wj, fij, psi_ij, Bn_ij, fluxij);
					unsigned long t_flux_end; RDTSCL0(&t_flux_end);
					delta_flux += t_flux_end - t_flux_beg;
					asm("#COMPUTEFLUX-END");

					const vreal  area( fij.area     );
					const vreal tarea(-fij.area * dt);
					Fluid_simd dUj(0.0);
					fluxij[Fluid::DIVB] = -Bn_ij;

					for (int k = 0; k < Fluid::NFLUID; k++) 
					{
						const vreal flux = fluxij[k] * tarea;
						v_dU [k] += flux;
						dUj  [k] -= flux;
					}


#ifdef __CONS_GRAV__
#if 0
					{
						const int *ii = fij.idx;
						const simd::vfloat gaccxj(gaccx[ii[0]], gaccx[ii[1]], gaccx[ii[2]], gaccx[ii[3]]);
						const simd::vfloat gaccyj(gaccy[ii[0]], gaccy[ii[1]], gaccy[ii[2]], gaccy[ii[3]]);
						const simd::vfloat gacczj(gaccz[ii[0]], gaccz[ii[1]], gaccz[ii[2]], gaccz[ii[3]]);

						const simd::vfloat gx = (gaccxi + gaccxj) * 0.5f;
						const simd::vfloat gy = (gaccyi + gaccyj) * 0.5f;
						const simd::vfloat gz = (gacczi + gacczj) * 0.5f;

#if 1
						const vreal gravij = -v_dU[Fluid::MASS] *
							(gx * fij.cx + gy * fij.cy + gz * fij.cz);
#else
						const vreal gravij = -v_dU[Fluid::MASS] *
							((gx * dx + gy * dy + gz * dz) * 0.5f);
#endif

						v_dU[Fluid::ENER] += gravij;
						dUj [Fluid::ENER] -= gravij;

					}
#endif
#endif


					v_divB      -= Bn_ij * tarea;

					psi_ij      *= area;
					v_gradPsi_x += psi_ij * fij.nx;
					v_gradPsi_y += psi_ij * fij.ny;
					v_gradPsi_z += psi_ij * fij.nz;


					asm("#WRITE_DATA1-BEG");
					const vreal divBj = Bn_ij * tarea;
					for (int ch = 0; ch < simd::vlen; ch++)
					{
						const int j = fij.idx[ch];
						if (ptcl[j].is_active()) continue;
						for (int k = 0; k < Fluid::NFLUID; k++) {
							dU[j][k] += dUj[k][ch];
						}
					}
					asm("#WRITE_DATA1-END");

#if 0
					const simd::vfloat dsh = dsl * 0.5f; //fij.cx*fij.nx + fij.cy*fij.ny + fij.cz*fij.nz;
					const simd::vfloat ids = simd::rcp_safe(dsh * simd::vfloat(2.0f));
					const simd::vfloat idsA = fij.area * ids;

					const simd::vfloat dxh = dx * 0.5f; //dsh * fij.nx;
					const simd::vfloat dyh = dy * 0.5f; //dsh * fij.ny;
					const simd::vfloat dzh = dz * 0.5f; //dsh * fij.nz;

					const simd::vfloat fx = fij.cx - dxh;
					const simd::vfloat fy = fij.cy - dyh;
					const simd::vfloat fz = fij.cz - dzh;

					const simd::vfloat vmx = (pj.vx + pi.vx).to_vfloat();
					const simd::vfloat vmy = (pj.vy + pi.vy).to_vfloat();
					const simd::vfloat vmz = (pj.vz + pi.vz).to_vfloat();
					const simd::vfloat dvx = (pj.vx - pi.vx).to_vfloat();
					const simd::vfloat dvy = (pj.vy - pi.vy).to_vfloat();
					const simd::vfloat dvz = (pj.vz - pi.vz).to_vfloat();

					const simd::vfloat vrelx = dvx;
					const simd::vfloat vrely = dvy; // - viy;
					const simd::vfloat vrelz = dvz; //- viz;
					const simd::vfloat vrel2 = sqr(vrelx) + sqr(vrely) + sqr(vrelz);

					virel2max = simd::max(virel2max, vrel2);


					v_dVdt -= idsA * ( 
							(fx *  dvx +  fy * dvy +  fz * dvz) -
							(dxh * vmx + dyh * vmy + dzh * vmz));
#endif
#if 1
#if 0
					virel2max = simd::max(virel2max, 
							(sqr(pj.vx - pi.vx) + sqr(pj.vy - pi.vy) + sqr(pj.vz - pi.vz)).to_vfloat());
#else
					virel2max = simd::max(virel2max, 
							(sqr(gradWj.m[Fluid::VELX] - gradWi.m[Fluid::VELX]) +
							 sqr(gradWj.m[Fluid::VELY] - gradWi.m[Fluid::VELY]) +
							 sqr(gradWj.m[Fluid::VELZ] - gradWi.m[Fluid::VELZ])).to_vfloat());
#endif
#endif
				}
				unsigned long t_inner_end; RDTSCL0(&t_inner_end);


				delta_inner += t_inner_end - t_inner_beg;

				//				const float dVdt = v_dVdt.sum(); 

				asm("#HOTSPOT1-END");
				Fluid dU(v_dU.sum());
				for (int k = 0; k < Fluid::NFLUID; k++)
				{
					dU[k]	+= this->dU[i][k];
					this->dU[i][k] = 0.0;
				}

				real divB = dU[Fluid::DIVB];

				real  gradPsi_x = v_gradPsi_x.sum() * tau;
				real  gradPsi_y = v_gradPsi_y.sum() * tau;
				real  gradPsi_z = v_gradPsi_z.sum() * tau;

				/*** predict privimite by half timestep ***/


#if 1 // divClean

#if 1
				Fluid Wh(Wgrad[i].m);
#if 1
				for (int k = 0; k < Fluid::NFLUID; k++)
					Wh[k] += Wgrad[i].t[k] * (tau * 0.5);
#endif

				dU[Fluid::WBX] -= Wh[Fluid::VELX] * divB;
				dU[Fluid::WBY] -= Wh[Fluid::VELY] * divB;
				dU[Fluid::WBZ] -= Wh[Fluid::VELZ] * divB;
#endif

#if 1       // subtracts force due to monopoles, but breaks conservation ... divBclean
				const real vB = 
					Wh[Fluid::VELX]*Wh[Fluid::BX] + 
					Wh[Fluid::VELY]*Wh[Fluid::BY] + 
					Wh[Fluid::VELZ]*Wh[Fluid::BZ];

				dU[Fluid::MOMX] -= Wh[Fluid::BX] * divB;
				dU[Fluid::MOMY] -= Wh[Fluid::BY] * divB;
				dU[Fluid::MOMZ] -= Wh[Fluid::BZ] * divB;
				dU[Fluid::ENER] -= vB     * divB;
#endif

#if 1     // Hyperbolic (due to Dender) divergence cleaning
				const real pres = compute_pressure(Wh[Fluid::DENS], Wh[Fluid::ETHM]);
				const real B2   = (sqr(Wh[Fluid::BX]) + sqr(Wh[Fluid::BY]) + sqr(Wh[Fluid::BZ]));
				const real dcs2 = (gamma_gas * pres + B2);

				dU[Fluid::WBX ] -= gradPsi_x;
				dU[Fluid::WBY ] -= gradPsi_y;
				dU[Fluid::WBZ ] -= gradPsi_z;
				dU[Fluid::MPSI] -= dcs2*divB;
				dU[Fluid::ENER] -= Wh[Fluid::BX]*gradPsi_x + Wh[Fluid::BY]*gradPsi_y + Wh[Fluid::BZ]*gradPsi_z;
#endif

#endif // divClean


				// compute conservative state at previous timestep
				Fluid Ui(U[i]);
				const real oldVol = cells[i].oldVol;

#ifdef __ENER_UB__
				Ui[Fluid::ENER] += (
						sqr(Ui[Fluid::MOMX]) +
						sqr(Ui[Fluid::MOMY]) +
						sqr(Ui[Fluid::MOMZ]) ) * 0.5 / Ui[Fluid::MASS];
#ifdef __ENER_U__
				Ui[Fluid::ENER] += ( 
						sqr(Ui[Fluid::WBX]) +
						sqr(Ui[Fluid::WBY]) +
						sqr(Ui[Fluid::WBZ]) ) * 0.5 / oldVol;
#endif
#endif

				const real m0 = Ui[Fluid::MASS];
				for (int k = 0; k < Fluid::NFLUID; k++)
					Ui[k] += dU[k];

				const real m1 = Ui[Fluid::MASS] > 0 ? Ui[Fluid::MASS] : m0;
				const real mh = (m0 + m1) * 0.5 * tau;


#if 0    // divClean

#if 1
				Fluid Wh(Wgrad[i].m);
				{
					const real m1o = Ui[Fluid::MASS];	
					Ui[Fluid::MASS] = m1;
					Fluid W1(Ui.to_primitive(ptcl[i].Volume));
					Ui[Fluid::MASS] = m1o;
					for (int k = 0; k < Fluid::NFLUID; k++)
						Wh[k] = (W1[k] + Wh[k]) * 0.5;
				}

				Ui[Fluid::WBX] -= Wh[Fluid::VELX] * divB;
				Ui[Fluid::WBY] -= Wh[Fluid::VELY] * divB;
				Ui[Fluid::WBZ] -= Wh[Fluid::VELZ] * divB;
#endif

#if 1       // subtracts force due to monopoles, but breaks conservation ... divBclean
				const real vB = 
					Wh[Fluid::VELX]*Wh[Fluid::BX] + 
					Wh[Fluid::VELY]*Wh[Fluid::BY] + 
					Wh[Fluid::VELZ]*Wh[Fluid::BZ];

				Ui[Fluid::MOMX] -= Wh[Fluid::BX] * divB;
				Ui[Fluid::MOMY] -= Wh[Fluid::BY] * divB;
				Ui[Fluid::MOMZ] -= Wh[Fluid::BZ] * divB;
				Ui[Fluid::ENER] -= vB     * divB;
#endif

#if 1     // Hyperbolic (due to Dender) divergence cleaning
				const real pres = compute_pressure(Wh[Fluid::DENS], Wh[Fluid::ETHM]);
				const real B2   = (sqr(Wh[Fluid::BX]) + sqr(Wh[Fluid::BY]) + sqr(Wh[Fluid::BZ]));
				const real dcs2 = (gamma_gas * pres + B2);

				Ui[Fluid::WBX ] -= gradPsi_x;
				Ui[Fluid::WBY ] -= gradPsi_y;
				Ui[Fluid::WBZ ] -= gradPsi_z;
				Ui[Fluid::MPSI] -= dcs2*divB;
				Ui[Fluid::ENER] -= Wh[Fluid::BX]*gradPsi_x + Wh[Fluid::BY]*gradPsi_y + Wh[Fluid::BZ]*gradPsi_z;
#endif
#endif    // divClean




#ifdef __ENER_UB__
				Ui[Fluid::ENER] -= (
						sqr(Ui[Fluid::MOMX]) +
						sqr(Ui[Fluid::MOMY]) +
						sqr(Ui[Fluid::MOMZ]) ) * 0.5 / m1;
#ifdef __ENER_U__
				Ui[Fluid::ENER] -= ( 
						sqr(Ui[Fluid::WBX]) +
						sqr(Ui[Fluid::WBY]) +
						sqr(Ui[Fluid::WBZ]) ) * 0.5 / ptcl[i].Volume;
#endif
#endif

#ifndef __CONS_GRAV__
#if 1
				const vec3 acc1(gaccx[i], gaccy[i], gaccz[i]);
				const vec3 acc0(pi.ax[0], pi.ay[0], pi.az[0]);
				const vec3 acc = (acc0 + acc1) * 0.5;
				Ui[Fluid::MOMX] += acc.x * mh;
				Ui[Fluid::MOMY] += acc.y * mh;
				Ui[Fluid::MOMZ] += acc.z * mh;
#endif
#endif

				U[i] = Ui;

				cells[i].divB = divB/tau/pi.Volume[0];
				ptcl[i].vrel2max = virel2max.max();
				//				ptcl[i].tlast    = t_global;
			}
		}

		unsigned long t_outer_end;	RDTSCL0(&t_outer_end);

#ifdef TICKPROFILE
		fprintf(stderr, "delta_outer      = %20ld \n", t_outer_end - t_outer_beg);
		fprintf(stderr, "delta_inner      = %20ld \n", delta_inner);
		fprintf(stderr, "delta_flux       = %20ld \n", delta_flux);
		fprintf(stderr, "delta_read_simd  = %20ld \n", delta_read_simd);
#endif

		unsigned long long ninter_glob;
		MPI_Allreduce(&ninteraction, &ninter_glob, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
		if (myproc == 0) 
			if (iteration >= i_log || all_active)
				fprintf(stderr, "ninteraction = %lld\n", ninter_glob);


		int dens_failed = 0, ethm_failed = 0, entr_failed = 0;
		for (int at = 0; at < nt; at++) 
		{
			const std::vector<int> &ptcl_list = active_ptcl[at];
			const int np = ptcl_list.size();
			for (int idx = 0; idx < np; idx++)
			{
				const int i = ptcl_list[idx];

				const Particle &pi =  ptcl[i];
				const Fluid W0 = Wgrad[i].m;
				assert(W0[Fluid::DENS] > 0.0);
				assert(W0[Fluid::ETHM] > 0.0);

				const double cr  = 0.025;
				U[i][Fluid::MPSI] *= std::exp(-courant_no*cr);

				Fluid W1 = U[i].to_primitive(pi.Volume);
				if (!(W1[Fluid::DENS] > 0.0))
				{
					if (dens_failed < 10)
					{
						const Fluid Ucopy(U[i]);
						const real  Vcopy(ptcl[i].Volume);
						dump_ptcl_info(i, " dens1 failed : ");
						ptcl[i].Volume = cells[i].oldVol;
						U[i] = Wgrad[i].m.to_conservative(ptcl[i].Volume);
						dump_ptcl_info(i, " dens0 failed : ");
						U[i] = Ucopy;
						ptcl[i].Volume = Vcopy;
					}
				}

#if 0
				if (W1[Fluid::DENS] <= 0.0) 
				{
					dens_failed++;
					W1 = W0;
					//					W1[Fluid::DENS] = W0[Fluid::DENS];
				}
				U[i] = W1.to_conservative(pi.Volume);
#endif

				set_boundary_ptcl(i);
				W1 = U[i].to_primitive(pi.Volume);

				assert(W1[Fluid::DENS] > 0.0);
				if (W1[Fluid::ETHM] <= 0.0) 
//				if (1)
				{
					ethm_failed++;
					if (entropy_scalar < 0)
					{
						W1[Fluid::ETHM] = W0[Fluid::ETHM]/W0[Fluid::DENS] * W1[Fluid::DENS];
					}
					else
					{
						real entropy = W1.scal(entropy_scalar);
						if (entropy <= 0.0) 
						{
							entropy = W0.scal(entropy_scalar);
							entr_failed++;
						}
						W1[Fluid::ETHM] = compute_ethm(W1[Fluid::DENS], entropy);
					}
				} 
				else if (entropy_scalar >= 0)
				{
					W1.scal(entropy_scalar) = compute_entropy(W1[Fluid::DENS], W1[Fluid::ETHM]);
				}

				if (!(W1[Fluid::ETHM] > 0.0))
				{
					if (entr_failed < 10)
						fprintf(stderr , "W1[Fluid::ETHM]= %g  entropy= %d %g  dens= %g\n",
								W1[Fluid::ETHM],
								entropy_scalar,
								W1.scal(entropy_scalar),
								W1[Fluid::DENS]);

				}
				assert(W1[Fluid::ETHM] > 0.0);
				U    [i]   = W1.to_conservative(pi.Volume);
				Wgrad[i].m = W1; 

				const real tau = t_global - ptcl[i].tlast;
				assert(tau > 0.0);
				Fluid dW;
				for (int k = 0; k < Fluid::NFLUID; k++)
				{
					dW[k] = (W1[k] - W0[k]);
				}

				assert(W0[Fluid::ETHM] > 0.0);
				if (!((dW[Fluid::ETHM] != W1[Fluid::ETHM]) && (dW[Fluid::DENS] != W1[Fluid::DENS])))
				{
					fprintf(stderr, " i= %d :: pos= %g DENS W0= %g  W1= %g  dW= %g\n",
							i,
							vec3(ptcl[i].x, ptcl[i].y, ptcl[i].z).abs(),
							W0[Fluid::DENS],
							W1[Fluid::DENS],
							dW[Fluid::DENS]);
					fprintf(stderr, " i= %d :: ETHM W0= %g  W1= %g  dW= %g\n",
							i,
							W0[Fluid::ETHM],
							W1[Fluid::ETHM],
							dW[Fluid::ETHM]);
				}
				assert(dW[Fluid::ETHM] != W1[Fluid::ETHM]);
				assert(dW[Fluid::DENS] != W1[Fluid::DENS]);
				extra_ptcl_update(W1, dW, i, tau);	

				U    [i]   = W1.to_conservative(pi.Volume);
				Wgrad[i].m = W1; 

			}
		}

		int dens_failed_glob, ethm_failed_glob, entr_failed_glob;
		MPI_Allreduce(&dens_failed, &dens_failed_glob, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&ethm_failed, &ethm_failed_glob, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&entr_failed, &entr_failed_glob, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		if (myproc == 0) {
			if (dens_failed_glob > 0 || ethm_failed_glob > 0) {
				if (iteration >= i_log || all_active || dens_failed_glob > 0 || entr_failed_glob > 0)
				{
					if (entropy_scalar >= 0)
						fprintf(ferr,  " >>>  failed dens= %d [ %g %c ]  ethm= %d [ %g %c ] entr= %d [ %g %c ]\n",
								dens_failed_glob, 100.0*dens_failed_glob/global_n, '%',
								ethm_failed_glob, 100.0*ethm_failed_glob/global_n, '%',
								entr_failed_glob, 100.0*entr_failed_glob/global_n, '%');
					else
						fprintf(ferr,  " >>>  failed dens= %d [ %g %c ]  ethm= %d [ %g %c ] \n",
								dens_failed_glob, 100.0*dens_failed_glob/global_n, '%',
								ethm_failed_glob, 100.0*ethm_failed_glob/global_n, '%');
				}
			}
		}


		// commit changed to imported particles as well ...	

		std::vector<Fluid> Fluid_send[NMAXPROC];
		std::vector<Fluid> Fluid_recv[NMAXPROC];
#if 0
		for (int p = 0; p < nproc; p++) 
		{
			Fluid_send[p].clear();
			Fluid_recv[p].clear();
		}
#endif

		for (size_t i = 0; i < return_ptcl_list.size(); i++) 
		{
			const std::pair<int, int> p = return_ptcl_list[i];
			const int proc = p.first;
			const int jidx = p.second;

			assert(jidx >= local_n);
			Fluid_send[proc].push_back(dU[jidx]);
		}

		myMPI::all2all(Fluid_send, Fluid_recv, nrecv_return, 1, myproc, nproc, mpi_debug_flag);

		int nrecv = 0;
		for (int p = 0; p < nproc; p++)
			nrecv += Fluid_recv[p].size();

		int iloc = 0;
		for (int p = 0; p < nproc; p++) 
			for (size_t i = 0; i < Fluid_recv[p].size(); i++)
			{
				const int idx = jidx_map[iloc++];
				assert(idx < local_n);
				if (!ptcl[idx].is_active())
					for (int k = 0; k < Fluid::NFLUID; k++)
						dU[idx][k] += Fluid_recv[p][i][k];
			}


		assert(iloc == (int)export_ptcl_list.size());

	}

	inline void system::compute_flux(
			const vreal      &wx, const vreal &wy, const vreal &wz,
			const Fluid_simd &Wi,
			const Fluid_simd &Wj,
			const Face_simd  &fij,
			vreal &psi_ij,
			vreal &Bn_ij,
			Fluid_simd &flux)
	{
		asm("#COMPUTEFLUX-CALL");
		// dpflop 24
		// spflop 120
		const simd::vfloat  ex = fij.nx;
		const simd::vfloat  ey = fij.ny;
		const simd::vfloat  ez = fij.nz;
		const simd::vfloat  ds = simd::sqrt(ex*ex + ey*ey);
		const simd::vfloat ids = simd::rcp_safe(ds);
		//spflop 3 + 1*20 + 1*10
		const simd::vfloat cosph = (ids == simd::vfloat(0.0f)).select(simd::vfloat(1.0f), ex * ids);
		const simd::vfloat sinph = (ids == simd::vfloat(0.0f)).select(simd::vfloat(0.0f), ey * ids);
		const simd::vfloat costh = ez;
		const simd::vfloat sinth = ds;
		//spflop 4
		const vreal Axx =  cosph*sinth;
		const vreal Axy =  sinth*sinph;
		const vreal Axz =  costh;
		const vreal Ayx = -sinph;
		const vreal Ayy =  cosph;
		const vreal Ayz =  vreal(0.0);
		const vreal Azx = -costh*cosph;
		const vreal Azy = -costh*sinph;
		const vreal Azz =  sinth;
		//spflop 7
		const vreal dens_L =     Wi[Fluid::DENS];
		const vreal ethm_L =     Wi[Fluid::ETHM];
		const vreal velx_L = Axx*Wi[Fluid::VELX] + Axy*Wi[Fluid::VELY] + Axz*Wi[Fluid::VELZ];
		const vreal vely_L = Ayx*Wi[Fluid::VELX] + Ayy*Wi[Fluid::VELY] + Ayz*Wi[Fluid::VELZ];
		const vreal velz_L = Azx*Wi[Fluid::VELX] + Azy*Wi[Fluid::VELY] + Azz*Wi[Fluid::VELZ];
		const vreal Bx_L   = Axx*Wi[Fluid::BX]   + Axy*Wi[Fluid::BY]   + Axz*Wi[Fluid::BZ];
		const vreal By_L   = Ayx*Wi[Fluid::BX]   + Ayy*Wi[Fluid::BY]   + Ayz*Wi[Fluid::BZ];
		const vreal Bz_L   = Azx*Wi[Fluid::BX]   + Azy*Wi[Fluid::BY]   + Azz*Wi[Fluid::BZ];
		const vreal psiL   =     Wi[Fluid::PSI];
		//dpflop 30
		const vreal dens_R =     Wj[Fluid::DENS];
		const vreal ethm_R =     Wj[Fluid::ETHM];
		const vreal velx_R = Axx*Wj[Fluid::VELX] + Axy*Wj[Fluid::VELY] + Axz*Wj[Fluid::VELZ];
		const vreal vely_R = Ayx*Wj[Fluid::VELX] + Ayy*Wj[Fluid::VELY] + Ayz*Wj[Fluid::VELZ];
		const vreal velz_R = Azx*Wj[Fluid::VELX] + Azy*Wj[Fluid::VELY] + Azz*Wj[Fluid::VELZ];
		const vreal Bx_R   = Axx*Wj[Fluid::BX]   + Axy*Wj[Fluid::BY]   + Axz*Wj[Fluid::BZ];
		const vreal By_R   = Ayx*Wj[Fluid::BX]   + Ayy*Wj[Fluid::BY]   + Ayz*Wj[Fluid::BZ];
		const vreal Bz_R   = Azx*Wj[Fluid::BX]   + Azy*Wj[Fluid::BY]   + Azz*Wj[Fluid::BZ];
		const vreal psiR   =     Wj[Fluid::PSI];
		//dpflop 30

		const vreal pres_L = compute_pressure(dens_L, ethm_L);
		const vreal pres_R = compute_pressure(dens_R, ethm_R);
		//dpflop  4
		const vreal cfl2 = (vreal(gamma_gas)*pres_L + sqr(Bx_L) + sqr(By_L) + sqr(Bz_L))/dens_L;
		const vreal cfr2 = (vreal(gamma_gas)*pres_R + sqr(Bx_R) + sqr(By_R) + sqr(Bz_R))/dens_R;
		const vreal ch   = simd::sqrt(simd::max(cfl2, cfr2));
		Bn_ij  = vreal(0.5)*((Bx_L + Bx_R) - (psiR - psiL)/ch);
		psi_ij = vreal(0.5)*((psiL + psiR) - (Bx_R - Bx_L)*ch);
		//spflop 15 + 1*10

		const vreal wn_ij = Axx*wx + Axy*wy + Axz*wz;
		//dpflop  5
		Fluid_simd F;
		asm("#CALL_RP-BEG");
		riemann_solver(
				F, 
				Bn_ij, wn_ij,
				dens_L, pres_L, ethm_L, velx_L, vely_L, velz_L, By_L, Bz_L,
				dens_R, pres_R, ethm_R, velx_R, vely_R, velz_R, By_R, Bz_R);
		asm("#CALL_RP-END");


		const vreal iAxx =  cosph*sinth;
		const vreal iAxy = -sinph;
		const vreal iAxz = -costh*cosph;
		const vreal iAyx =  sinth*sinph;
		const vreal iAyy =  cosph;
		const vreal iAyz = -costh*sinph;
		const vreal iAzx =  costh;
		const vreal iAzy =  vreal(0.0);
		const vreal iAzz =  sinth;
		//sflop 7
		flux[Fluid::MASS] =      F[Fluid::MASS];
		flux[Fluid::ENER] =      F[Fluid::ENER];
		flux[Fluid::MOMX] = iAxx*F[Fluid::MOMX] + iAxy*F[Fluid::MOMY] + iAxz*F[Fluid::MOMZ];
		flux[Fluid::MOMY] = iAyx*F[Fluid::MOMX] + iAyy*F[Fluid::MOMY] + iAyz*F[Fluid::MOMZ];
		flux[Fluid::MOMZ] = iAzx*F[Fluid::MOMX] + iAzy*F[Fluid::MOMY] + iAzz*F[Fluid::MOMZ];
		flux[Fluid::WBX ] = iAxx*F[Fluid::WBX ] + iAxy*F[Fluid::WBY ] + iAxz*F[Fluid::WBZ ];
		flux[Fluid::WBY ] = iAyx*F[Fluid::WBX ] + iAyy*F[Fluid::WBY ] + iAyz*F[Fluid::WBZ ];
		flux[Fluid::WBZ ] = iAzx*F[Fluid::WBX ] + iAzy*F[Fluid::WBY ] + iAzz*F[Fluid::WBZ ];
		flux[Fluid::MPSI] =      F[Fluid::MASS] * (F[Fluid::MASS] > vreal(0.0)).select(psiL, psiR);
		//dpflop 32

		for (int k = 0; k < Fluid::NSCALARS; k++) 
			flux.scal(k) = F[Fluid::MASS] * (F[Fluid::MASS] > vreal(0.0)).select(Wi.scal(k), Wj.scal(k));
		//dpflop 6
		asm("#COMPUTEFLUX-EXIT");
	}

	inline void system::riemann_solver(
			Fluid_simd &flux,
			const vreal Bx,
			const vreal w,
			const vreal dens_L, const vreal pres_L, const vreal ethm_L,
			const vreal velx_L, const vreal vely_L, const vreal velz_L, const vreal By_L, const vreal Bz_L,
			const vreal dens_R, const vreal pres_R, const vreal ethm_R,
			const vreal velx_R, const vreal vely_R, const vreal velz_R, const vreal By_R, const vreal Bz_R
			) const
	{
		asm("#RP_BEG");

#if 1
		{


			const vreal signBx = Bx * simd::abs(simd::rcp_safe(Bx));

			const vreal momx_L = dens_L*velx_L;
			const vreal momy_L = dens_L*vely_L;
			const vreal momz_L = dens_L*velz_L;

			const vreal momx_R = dens_R*velx_R;
			const vreal momy_R = dens_R*vely_R;
			const vreal momz_R = dens_R*velz_R;
			//flop 6
			const vreal chalf = vreal(0.5);

			const vreal B2_L   = sqr(Bx)     + sqr(By_L)   + sqr(Bz_L);
			const vreal v2_L   = sqr(velx_L) + sqr(vely_L) + sqr(velz_L);
			const vreal etot_L = ethm_L      + chalf*(dens_L*v2_L + B2_L);
			//flop 14
			const vreal B2_R   = sqr(Bx)     + sqr(By_R)   + sqr(Bz_R);
			const vreal v2_R   = sqr(velx_R) + sqr(vely_R) + sqr(velz_R);
			const vreal etot_R = ethm_R      + chalf*(dens_R*v2_R + B2_R);
			//flop 14
			const vreal gpl  = vreal(gamma_gas) * pres_L;
			const vreal gpr  = vreal(gamma_gas) * pres_R;
			const vreal gpbl = gpl + B2_L;
			const vreal gpbr = gpr + B2_R;
			//flop 4
			//flop//flop//flop//flop//flop

			const vreal cfl2  = gpbl * simd::rcp(dens_L);
			const vreal cfr2  = gpbr * simd::rcp(dens_R);
			const vreal cfmax = simd::sqrt(simd::max(cfl2, cfr2));
			//flop 6 + 1*10 + 2*10
			const vreal S_L = simd::min(velx_L, velx_R) - cfmax;
			const vreal S_R = simd::max(velx_L, velx_R) + cfmax;
			//flop 4
			//flop//flop//flop//flop//flop

			const vreal pT_L = pres_L + chalf * B2_L;
			const vreal pT_R = pres_R + chalf * B2_R;
			//flop 6
			const vreal iSM = simd::rcp((S_R - velx_R)*dens_R - (S_L - velx_L)*dens_L);
			const vreal S_M  =   iSM * ((S_R - velx_R)*momx_R - (S_L - velx_L)*momx_L - pT_R + pT_L);
			//flop 13 + 1*10
			const vreal ipTs = simd::rcp((S_R - velx_R)*dens_R - (S_L - velx_L)*dens_L);
			const vreal pT_s = ipTs * ((S_R - velx_R)*dens_R*pT_L - (S_L - velx_L)*dens_L*pT_R +
					dens_L*dens_R*(S_R - velx_R)*(S_L - velx_L)*(velx_R - velx_L));
			//flop 22 + 1*10
			const vreal velx_L_s  = S_M;
			const vreal velx_L_ss = S_M;
			const vreal velx_R_s  = S_M;
			const vreal velx_R_ss = S_M;
			const vreal B2x       = Bx*Bx;
			//flop 24
			const vreal iSLmSM   = simd::rcp(S_L - S_M);
			const vreal iSRmSM   = simd::rcp(S_R - S_M);
			const vreal dens_L_s = dens_L * (S_L - velx_L) * iSLmSM;
			const vreal dens_R_s = dens_R * (S_R - velx_R) * iSRmSM;
			const vreal divL     = dens_L * (S_L - velx_L)*(S_L - S_M) - B2x;
			const vreal divR     = dens_R * (S_R - velx_R)*(S_R - S_M) - B2x;
			const vreal idivL    = simd::rcp_safe(divL);
			const vreal idivR    = simd::rcp_safe(divR);
			const vreal vely_L_s = vely_L - Bx*By_L*(S_M - velx_L) * idivL;
			const vreal velz_L_s = velz_L - Bx*Bz_L*(S_M - velx_L) * idivL;
			const vreal   By_L_s = By_L * (dens_L*sqr(S_L - velx_L) - B2x) * idivL;
			const vreal   Bz_L_s = Bz_L * (dens_L*sqr(S_L - velx_L) - B2x) * idivL;
			//flop 42 + 2*10 + 2*20
			const vreal vely_R_s = vely_R - Bx*By_R*(S_M - velx_R) * idivR;
			const vreal velz_R_s = velz_R - Bx*Bz_R*(S_M - velx_R) * idivR;
			const vreal   By_R_s = By_R * (dens_R*sqr(S_R - velx_R) - B2x) * idivR;
			const vreal   Bz_R_s = Bz_R * (dens_R*sqr(S_R - velx_R) - B2x) * idivR;
			//flop 22
			const vreal   vB_L   = velx_L  *Bx + vely_L  *By_L   + velz_L  *Bz_L;
			const vreal   vB_L_s = velx_L_s*Bx + vely_L_s*By_L_s + velz_L_s*Bz_L_s;
			const vreal etot_L_s = ((S_L - velx_L)*etot_L - pT_L*velx_L + pT_s*S_M + Bx*(vB_L - vB_L_s)) * iSLmSM;
			//flop 20
			const vreal   vB_R   = velx_R  *Bx + vely_R  *By_R   + velz_R  *Bz_R;
			const vreal   vB_R_s = velx_R_s*Bx + vely_R_s*By_R_s + velz_R_s*Bz_R_s;
			const vreal etot_R_s = ((S_R - velx_R)*etot_R - pT_R*velx_R + pT_s*S_M + Bx*(vB_R - vB_R_s)) * iSRmSM;
			//flop 20
			const vreal dens_L_ss = dens_L_s;
			const vreal dens_R_ss = dens_R_s;
			const vreal sDens_L_s = simd::sqrt(dens_L_s);
			const vreal sDens_R_s = simd::sqrt(dens_R_s);
			//flop 2*10
			const vreal    S_L_s  = S_M - simd::abs(Bx)*simd::rcp(sDens_L_s);
			const vreal    S_R_s  = S_M + simd::abs(Bx)*simd::rcp(sDens_R_s);
			//flop 6 + 2*10
			const vreal idsqroot  = simd::rcp(sDens_L_s + sDens_R_s);
			const vreal  vely_ss = idsqroot*(sDens_L_s*vely_L_s + sDens_R_s*vely_R_s + (By_R_s - By_L_s)*signBx);
			const vreal  velz_ss = idsqroot*(sDens_L_s*velz_L_s + sDens_R_s*velz_R_s + (Bz_R_s - Bz_L_s)*signBx);
			//flop 15 + 1*10
			const vreal By_ss = idsqroot*(sDens_L_s*By_R_s + sDens_R_s*By_L_s + sDens_L_s*sDens_R_s*(vely_R_s - vely_L_s)*signBx);
			const vreal Bz_ss = idsqroot*(sDens_L_s*Bz_R_s + sDens_R_s*Bz_L_s + sDens_L_s*sDens_R_s*(velz_R_s - velz_L_s)*signBx);
			//flop 18
			const vreal vely_L_ss = vely_ss;
			const vreal velz_L_ss = velz_ss;
			const vreal   By_L_ss = By_ss;
			const vreal   Bz_L_ss = Bz_ss;

			const vreal vely_R_ss = vely_ss;
			const vreal velz_R_ss = velz_ss;
			const vreal   By_R_ss = By_ss;
			const vreal   Bz_R_ss = Bz_ss;

			const vreal vB_L_ss   = velx_L_ss*Bx + vely_L_ss*By_L_ss + velz_L_ss*Bz_L_ss;
			const vreal etot_L_ss = etot_L_s - sDens_L_s*(vB_L_s - vB_L_ss)*signBx;
			//flop 8
			const vreal vB_R_ss   = velx_R_ss*Bx + vely_R_ss*By_R_ss + velz_R_ss*Bz_R_ss;
			const vreal etot_R_ss = etot_R_s + sDens_R_s*(vB_R_s - vB_R_ss)*signBx;
			//flop 8
			const vreal Fdens_L = dens_L*velx_L;
			const vreal Fmomx_L = momx_L*velx_L + pT_L - B2x;
			const vreal Fmomy_L = momy_L*velx_L        - Bx*By_L;
			const vreal Fmomz_L = momz_L*velx_L        - Bx*Bz_L;
			const vreal Fetot_L = etot_L*velx_L + pT_L*velx_L - Bx*vB_L; 
			//flop 15
			const vreal Fdens_R = dens_R*velx_R;
			const vreal Fmomx_R = momx_R*velx_R + pT_R - B2x;
			const vreal Fmomy_R = momy_R*velx_R        - Bx*By_R;
			const vreal Fmomz_R = momz_R*velx_R        - Bx*Bz_R;
			const vreal Fetot_R = etot_R*velx_R + pT_R*velx_R - Bx*vB_R;
			//flop 15
			const vreal momx_L_s  = dens_L_s *velx_L_s;
			const vreal momy_L_s  = dens_L_s *vely_L_s;
			const vreal momz_L_s  = dens_L_s *velz_L_s;

			const vreal momx_L_ss = dens_L_ss*velx_L_ss;
			const vreal momy_L_ss = dens_L_ss*vely_L_ss;
			const vreal momz_L_ss = dens_L_ss*velz_L_ss;

			const vreal momx_R_s  = dens_R_s *velx_R_s;
			const vreal momy_R_s  = dens_R_s *vely_R_s;
			const vreal momz_R_s  = dens_R_s *velz_R_s;

			const vreal momx_R_ss = dens_R_ss*velx_R_ss;
			const vreal momy_R_ss = dens_R_ss*vely_R_ss;
			const vreal momz_R_ss = dens_R_ss*velz_R_ss;

			const vreal Fby_L  = By_L*velx_L - Bx * vely_L;
			const vreal Fbz_L  = Bz_L*velx_L - Bx * velz_L;

			const vreal Fby_R  = By_R*velx_R - Bx * vely_R;
			const vreal Fbz_R  = Bz_R*velx_R - Bx * velz_R;
			//flop 24
			const vreal w_lt_SL  = w <= S_L;
			const vreal w_le_SLs = w <= S_L_s;
			const vreal w_le_SM  = w <= S_M;
			const vreal w_le_SRs = w <= S_R_s;
			const vreal w_le_SR  = w <= S_R;
			//flop 5
			const vreal fdens = (w_le_SM).select(Fdens_L, Fdens_R);
			const vreal fetot = (w_le_SM).select(Fetot_L, Fetot_R);
			const vreal fmomx = (w_le_SM).select(Fmomx_L, Fmomx_R);
			const vreal fmomy = (w_le_SM).select(Fmomy_L, Fmomy_R);
			const vreal fmomz = (w_le_SM).select(Fmomz_L, Fmomz_R);
			const vreal fby   = (w_le_SM).select(Fby_L,   Fby_R);
			const vreal fbz   = (w_le_SM).select(Fbz_L,   Fbz_R);
			//flop 7
			const vreal cnull(0.0);
			const vreal a = 
				(w_le_SLs).select(
						cnull, (w_le_SM).select(
							S_L_s - w, (w_le_SRs).select(
								S_R_s - w, 
								cnull)));
			//flop 5
			const vreal b  = 
				(w_lt_SL).select(
						cnull, (w_le_SLs).select(
							S_L - w, (w_le_SM).select(
								- S_L_s + S_L, (w_le_SRs).select(
									-	S_R_s + S_R, (w_le_SR).select(
										S_R - w, 
										cnull)))));
			//flop 11
			const vreal c = 
				(w_lt_SL).select(
						-w, (w_le_SLs).select(
							-S_L, (w_le_SM).select(
								-S_L, (w_le_SRs).select(
									-S_R, (w_le_SR).select(
										-S_R, 
										-w)))));
			//flop 11
			const vreal dens    = (w_le_SM).select(dens_L,    dens_R   );
			const vreal dens_s  = (w_le_SM).select(dens_L_s,  dens_R_s );
			const vreal dens_ss = (w_le_SM).select(dens_L_ss, dens_R_ss);

			const vreal etot    = (w_le_SM).select(etot_L,    etot_R   );
			const vreal etot_s  = (w_le_SM).select(etot_L_s,  etot_R_s );
			const vreal etot_ss = (w_le_SM).select(etot_L_ss, etot_R_ss);

			const vreal momx    = (w_le_SM).select(momx_L,    momx_R   );
			const vreal momx_s  = (w_le_SM).select(momx_L_s,  momx_R_s );
			const vreal momx_ss = (w_le_SM).select(momx_L_ss, momx_R_ss);

			const vreal momy    = (w_le_SM).select(momy_L,    momy_R   );
			const vreal momy_s  = (w_le_SM).select(momy_L_s,  momy_R_s );
			const vreal momy_ss = (w_le_SM).select(momy_L_ss, momy_R_ss);

			const vreal momz    = (w_le_SM).select(momz_L,    momz_R   );
			const vreal momz_s  = (w_le_SM).select(momz_L_s,  momz_R_s );
			const vreal momz_ss = (w_le_SM).select(momz_L_ss, momz_R_ss);

			const vreal by    = (w_le_SM).select(By_L,    By_R   );
			const vreal by_s  = (w_le_SM).select(By_L_s,  By_R_s );
			const vreal by_ss = (w_le_SM).select(By_L_ss, By_R_ss);

			const vreal bz    = (w_le_SM).select(Bz_L,    Bz_R   );
			const vreal bz_s  = (w_le_SM).select(Bz_L_s,  Bz_R_s );
			const vreal bz_ss = (w_le_SM).select(Bz_L_ss, Bz_R_ss);
			//flop 21
			flux[Fluid::MASS] = fdens + a * dens_ss + b * dens_s + c * dens;
			flux[Fluid::ENER] = fetot + a * etot_ss + b * etot_s + c * etot;
			flux[Fluid::MOMX] = fmomx + a * momx_ss + b * momx_s + c * momx;
			flux[Fluid::MOMY] = fmomy + a * momy_ss + b * momy_s + c * momy;
			flux[Fluid::MOMZ] = fmomz + a * momz_ss + b * momz_s + c * momz;
			flux[Fluid::WBY ] = fby   + a * by_ss   + b * by_s   + c * by;
			flux[Fluid::WBZ ] = fbz   + a * bz_ss   + b * bz_s   + c * bz;
			flux[Fluid::WBX ] = -w * Bx;

#if 0	
			{
				const vreal w_lt_SL  = cnull <= S_L;
				const vreal w_le_SLs = cnull <= S_L_s;
				const vreal w_le_SM  = cnull <= S_M;
				const vreal w_le_SRs = cnull <= S_R_s;
				const vreal w_le_SR  = cnull <= S_R;

				const vreal fdens = (w_le_SM).select(Fdens_L, Fdens_R);

				const vreal dens    = (w_le_SM).select(dens_L,    dens_R   );
				const vreal dens_s  = (w_le_SM).select(dens_L_s,  dens_R_s );
				const vreal dens_ss = (w_le_SM).select(dens_L_ss, dens_R_ss);

				const vreal a = 
					(w_le_SLs).select(
							cnull, (w_le_SM).select(
								S_L_s , (w_le_SRs).select(
									S_R_s , 
									cnull)));

				const vreal b  = 
					(w_lt_SL).select(
							cnull, (w_le_SLs).select(
								S_L , (w_le_SM).select(
									- S_L_s + S_L, (w_le_SRs).select(
										-	S_R_s + S_R, (w_le_SR).select(
											S_R , 
											cnull)))));

				const vreal c = 
					(w_lt_SL).select(
							cnull, (w_le_SLs).select(
								-S_L, (w_le_SM).select(
									-S_L, (w_le_SRs).select(
										-S_R, (w_le_SR).select(
											-S_R, 
											cnull)))));

				mflux = fdens + a * dens_ss + b * dens_s + c * dens;
			}
#endif

		}
#else // HLLE
		{
			const vreal signBx = Bx * simd::abs(simd::rcp_safe(Bx));

			const vreal momx_L = dens_L*velx_L;
			const vreal momy_L = dens_L*vely_L;
			const vreal momz_L = dens_L*velz_L;

			const vreal momx_R = dens_R*velx_R;
			const vreal momy_R = dens_R*vely_R;
			const vreal momz_R = dens_R*velz_R;

			const vreal chalf = vreal(0.5);

			const vreal B2_L   = sqr(Bx)     + sqr(By_L)   + sqr(Bz_L);
			const vreal v2_L   = sqr(velx_L) + sqr(vely_L) + sqr(velz_L);
			const vreal etot_L = ethm_L      + chalf*(dens_L*v2_L + B2_L);
			//flop 14
			const vreal B2_R   = sqr(Bx)     + sqr(By_R)   + sqr(Bz_R);
			const vreal v2_R   = sqr(velx_R) + sqr(vely_R) + sqr(velz_R);
			const vreal etot_R = ethm_R      + chalf*(dens_R*v2_R + B2_R);
			
			
			const vreal pT_L = pres_L + chalf * B2_L;
			const vreal pT_R = pres_R + chalf * B2_R;

			const vreal gpl  = vreal(gamma_gas) * pres_L;
			const vreal gpr  = vreal(gamma_gas) * pres_R;
			const vreal gpbl = gpl + B2_L;
			const vreal gpbr = gpr + B2_R;

			const vreal cfl2  = gpbl * simd::rcp(dens_L);
			const vreal cfr2  = gpbr * simd::rcp(dens_R);
			const vreal cfmax = simd::sqrt(simd::max(cfl2, cfr2)) * vreal(1.2);
			
			//flop 6 + 1*10 + 2*10
			const vreal S_L = simd::min(velx_L, velx_R) - cfmax;
			const vreal S_R = simd::max(velx_L, velx_R) + cfmax;
			
			vreal vB_L = velx_L*Bx + vely_L*By_L   + velz_L*Bz_L;
			vreal vB_R = velx_R*Bx + vely_R*By_R   + velz_R*Bz_R;

			vreal Fdens_L = dens_L*velx_L;
			vreal Fmomx_L = momx_L*velx_L + pT_L - Bx*Bx;
			vreal Fmomy_L = momy_L*velx_L        - Bx*By_L; 
			vreal Fmomz_L = momz_L*velx_L        - Bx*Bz_L;
			vreal Fetot_L = etot_L*velx_L + pT_L*velx_L - Bx*vB_L;
			vreal Fby_L   = By_L  *velx_L - Bx*vely_L;
			vreal Fbz_L   = Bz_L  *velx_L - Bx*velz_L;

			vreal Fdens_R = dens_R*velx_R;
			vreal Fmomx_R = momx_R*velx_R + pT_R - Bx*Bx;
			vreal Fmomy_R = momy_R*velx_R        - Bx*By_R; 
			vreal Fmomz_R = momz_R*velx_R        - Bx*Bz_R; 
			vreal Fetot_R = etot_R*velx_R + pT_R*velx_R - Bx*vB_R;
			vreal Fby_R   = By_R  *velx_R        - Bx*vely_R;
			vreal Fbz_R   = Bz_R  *velx_R        - Bx*velz_R;


			vreal U_dens = (S_R*dens_R - S_L*dens_L + Fdens_L - Fdens_R)/(S_R - S_L);
			vreal U_momx = (S_R*momx_R - S_L*momx_L + Fmomx_L - Fmomx_R)/(S_R - S_L);
			vreal U_momy = (S_R*momy_R - S_L*momy_L + Fmomy_L - Fmomy_R)/(S_R - S_L);
			vreal U_momz = (S_R*momz_R - S_L*momz_L + Fmomz_L - Fmomz_R)/(S_R - S_L);
			vreal U_etot = (S_R*etot_R - S_L*etot_L + Fetot_L - Fetot_R)/(S_R - S_L);

			flux[Fluid::MASS] = (S_R*Fdens_L - S_L*Fdens_R + S_L*S_R*(dens_R - dens_L))/(S_R - S_L);
			flux[Fluid::MOMX] = (S_R*Fmomx_L - S_L*Fmomx_R + S_L*S_R*(momx_R - momx_L))/(S_R - S_L);
			flux[Fluid::MOMY] = (S_R*Fmomy_L - S_L*Fmomy_R + S_L*S_R*(momy_R - momy_L))/(S_R - S_L);
			flux[Fluid::MOMZ] = (S_R*Fmomz_L - S_L*Fmomz_R + S_L*S_R*(momz_R - momz_L))/(S_R - S_L);
			flux[Fluid::ENER] = (S_R*Fetot_L - S_L*Fetot_R + S_L*S_R*(etot_R - etot_L))/(S_R - S_L);

			vreal U_by = (S_R*By_R - S_L*By_L + Fby_L - Fby_R)/(S_R - S_L);
		  vreal U_bz = (S_R*Bz_R - S_L*Bz_L + Fbz_L - Fbz_R)/(S_R - S_L);

			flux[Fluid::WBY] = (S_R*Fby_L - S_L*Fby_R + S_L*S_R*(By_R - By_L))/(S_R - S_L);
			flux[Fluid::WBZ] = (S_R*Fbz_L - S_L*Fbz_R + S_L*S_R*(Bz_R - Bz_L))/(S_R - S_L);

			flux[Fluid::MASS] -= w*U_dens;
			flux[Fluid::MOMX] -= w*U_momx;
			flux[Fluid::MOMY] -= w*U_momy;
			flux[Fluid::MOMZ] -= w*U_momz;
			flux[Fluid::ENER] -= w*U_etot;
			flux[Fluid::WBY]  -= w*U_by;
			flux[Fluid::WBZ]  -= w*U_bz;
			
			flux[Fluid::WBX]   = -w*Bx;

		}
#endif

		//flop 45
		asm("#RP_END");
	}

};
