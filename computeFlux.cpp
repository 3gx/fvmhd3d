#include "fvmhd3d.h"

namespace fvmhd3d
{
  
  void System::computeFlux(
      const vec3  &wij,
      const Fluid &Wi,
      const Fluid &Wj,
      const int    jbnd,
      const vec3  &normal,
      real &psi_ij,
      real &Bn_ij,
      Fluid &flux)
  {
    asm("#COMPUTEFLUX-CALL");
    // dpflop 24
    // spflop 120
    const real  ex = normal.x;
    const real  ey = normal.y;
    const real  ez = normal.z;
    const real  ds = std::sqrt(ex*ex + ey*ey);
    const real ids = (ds != 0.0) ? 1.0/ds : 0.0;
    //spflop 3 + 1*20 + 1*10
    const real cosph = (ids == 0.0) ? 1.0 : ex * ids;
    const real sinph = (ids == 0.0) ? 0.0 : ey * ids;
    const real costh = ez;
    const real sinth = ds;
    //spflop 4
    const real Axx =  cosph*sinth;
    const real Axy =  sinth*sinph;
    const real Axz =  costh;
    const real Ayx = -sinph;
    const real Ayy =  cosph;
    const real Ayz =  0.0;
    const real Azx = -costh*cosph;
    const real Azy = -costh*sinph;
    const real Azz =  sinth;
    //spflop 7
    real dens_L =     Wi[Fluid::DENS];
    real ethm_L =     Wi[Fluid::ETHM];
    real velx_L = Axx*Wi[Fluid::VELX] + Axy*Wi[Fluid::VELY] + Axz*Wi[Fluid::VELZ];
    real vely_L = Ayx*Wi[Fluid::VELX] + Ayy*Wi[Fluid::VELY] + Ayz*Wi[Fluid::VELZ];
    real velz_L = Azx*Wi[Fluid::VELX] + Azy*Wi[Fluid::VELY] + Azz*Wi[Fluid::VELZ];
    real Bx_L   = Axx*Wi[Fluid::BX]   + Axy*Wi[Fluid::BY]   + Axz*Wi[Fluid::BZ];
    real By_L   = Ayx*Wi[Fluid::BX]   + Ayy*Wi[Fluid::BY]   + Ayz*Wi[Fluid::BZ];
    real Bz_L   = Azx*Wi[Fluid::BX]   + Azy*Wi[Fluid::BY]   + Azz*Wi[Fluid::BZ];
    real psiL   =     Wi[Fluid::PSI];
    real entrL  =     Wi[Fluid::ENTR];
    //dpflop 30
    real dens_R =     Wj[Fluid::DENS];
    real ethm_R =     Wj[Fluid::ETHM];
    real velx_R = Axx*Wj[Fluid::VELX] + Axy*Wj[Fluid::VELY] + Axz*Wj[Fluid::VELZ];
    real vely_R = Ayx*Wj[Fluid::VELX] + Ayy*Wj[Fluid::VELY] + Ayz*Wj[Fluid::VELZ];
    real velz_R = Azx*Wj[Fluid::VELX] + Azy*Wj[Fluid::VELY] + Azz*Wj[Fluid::VELZ];
    real Bx_R   = Axx*Wj[Fluid::BX]   + Axy*Wj[Fluid::BY]   + Axz*Wj[Fluid::BZ];
    real By_R   = Ayx*Wj[Fluid::BX]   + Ayy*Wj[Fluid::BY]   + Ayz*Wj[Fluid::BZ];
    real Bz_R   = Azx*Wj[Fluid::BX]   + Azy*Wj[Fluid::BY]   + Azz*Wj[Fluid::BZ];
    real psiR   =     Wj[Fluid::PSI];
    real entrR  =     Wj[Fluid::ENTR];
    //dpflop 30

    real pres_L = Problem_compute_pressure(Wi);
    real pres_R = Problem_compute_pressure(Wj);

    const real wn_ij = vec3(Axx, Axy, Axz)*wij;
    switch(jbnd)
    {
      case MeshPoint::REFLECTING:
        velx_R = -velx_R;
        break;
      case MeshPoint::DIOD:
        velx_R = std::max(velx_L - 1.0*wn_ij, 0.0) + 1.0*wn_ij;
        vely_R = vely_L;
        velz_R = velz_L;
        dens_R = dens_L;
        pres_R = pres_L;
        Bx_R   = Bx_L;
        By_R   = By_L;
        Bz_R   = Bz_L;
        psiR   = psiL;
        entrR  = entrL;
        break;
      case MeshPoint::INFLOW:
        break;
			case MeshPoint::OUTFLOW:
        dens_R = dens_L;
        pres_R = pres_L;
        Bx_R   = Bx_L;
        By_R   = By_L;
        Bz_R   = Bz_L;
        psiR   = psiL;
        entrR  = entrL;
				break;
      default:
        break;
    };


    //dpflop  4
    const real cfl2 = (gamma_gas*pres_L + sqr(Bx_L) + sqr(By_L) + sqr(Bz_L))/dens_L;
    const real cfr2 = (gamma_gas*pres_R + sqr(Bx_R) + sqr(By_R) + sqr(Bz_R))/dens_R;
    const real chL = std::sqrt(cfl2 + 0.0*sqr(velx_L - wn_ij));
    const real chR = std::sqrt(cfr2 + 0.0*sqr(velx_R - wn_ij));

    Bn_ij  = (Bx_L*chL + Bx_R*chR +         (psiL - psiR))/(chL + chR);
    psi_ij = (psiL*chR + psiR*chL + chR*chL*(Bx_L - Bx_R))/(chL + chR);
    //spflop 15 + 1*10


    //dpflop  5
    Fluid F;
    asm("#CALL_RP-BEG");
    riemannSolver(
        F, 
        Bn_ij, wn_ij,
        dens_L, pres_L, ethm_L, velx_L, vely_L, velz_L, By_L, Bz_L,
        dens_R, pres_R, ethm_R, velx_R, vely_R, velz_R, By_R, Bz_R);
    asm("#CALL_RP-END");

    const real iAxx =  cosph*sinth;
    const real iAxy = -sinph;
    const real iAxz = -costh*cosph;
    const real iAyx =  sinth*sinph;
    const real iAyy =  cosph;
    const real iAyz = -costh*sinph;
    const real iAzx =  costh;
    const real iAzy =  real(0.0);
    const real iAzz =  sinth;
    //sflop 7
    flux[Fluid::MASS] =      F[Fluid::MASS];
    flux[Fluid::ENER] =      F[Fluid::ENER];
    flux[Fluid::MOMX] = iAxx*F[Fluid::MOMX] + iAxy*F[Fluid::MOMY] + iAxz*F[Fluid::MOMZ];
    flux[Fluid::MOMY] = iAyx*F[Fluid::MOMX] + iAyy*F[Fluid::MOMY] + iAyz*F[Fluid::MOMZ];
    flux[Fluid::MOMZ] = iAzx*F[Fluid::MOMX] + iAzy*F[Fluid::MOMY] + iAzz*F[Fluid::MOMZ];
    flux[Fluid::WBX ] = iAxx*F[Fluid::WBX ] + iAxy*F[Fluid::WBY ] + iAxz*F[Fluid::WBZ ];
    flux[Fluid::WBY ] = iAyx*F[Fluid::WBX ] + iAyy*F[Fluid::WBY ] + iAyz*F[Fluid::WBZ ];
    flux[Fluid::WBZ ] = iAzx*F[Fluid::WBX ] + iAzy*F[Fluid::WBY ] + iAzz*F[Fluid::WBZ ];
    flux[Fluid::MPSI] =      F[Fluid::MASS] * (F[Fluid::MASS] > 0.0 ? psiL  : psiR);
    //    flux[Fluid::MPSI] = F[Fluid::MASS] * psi_ij;
    flux[Fluid::MENTR] =     F[Fluid::MASS] * (F[Fluid::MASS] > 0.0 ? entrL : entrR);
    //dpflop 32

    for (int k = 0; k < Fluid::NSCALARS; k++) 
      flux.scal(k) = F[Fluid::MASS] * (F[Fluid::MASS] > 0.0 ? Wi.scal(k) : Wj.scal(k));
    //dpflop 6
    asm("#COMPUTEFLUX-EXIT");
  }


  void System::riemannSolver(
      Fluid &flux,
      const real Bx,
      const real w,
      const real dens_L, const real pres_L, const real ethm_L,
      const real velx_L, const real vely_L, const real velz_L, const real By_L, const real Bz_L,
      const real dens_R, const real pres_R, const real ethm_R,
      const real velx_R, const real vely_R, const real velz_R, const real By_R, const real Bz_R
      ) const
  {
#if 1
    riemannSolver_HLLD(
        flux, Bx, w, 
        dens_L,  pres_L,  ethm_L,
        velx_L,  vely_L,  velz_L,  By_L,  Bz_L,
        dens_R,  pres_R,  ethm_R,
        velx_R,  vely_R,  velz_R,  By_R,  Bz_R);
#else
    riemannSolver_HLLE(
        flux, Bx, w, 
        dens_L,  pres_L,  ethm_L,
        velx_L,  vely_L,  velz_L,  By_L,  Bz_L,
        dens_R,  pres_R,  ethm_R,
        velx_R,  vely_R,  velz_R,  By_R,  Bz_R);
#endif
  }

  void System::riemannSolver_HLLD(
      Fluid &flux,
      const real Bx,
      const real w,
      const real dens_L, const real pres_L, const real ethm_L,
      const real velx_L, const real vely_L, const real velz_L, const real By_L, const real Bz_L,
      const real dens_R, const real pres_R, const real ethm_R,
      const real velx_R, const real vely_R, const real velz_R, const real By_R, const real Bz_R
      ) const
  {
    asm("#RP_BEG");

    {
      const real signBx = Bx == 0.0 ? 0.0 : (Bx > 0.0 ? +1.0 : -1.0);

      const real momx_L = dens_L*velx_L;
      const real momy_L = dens_L*vely_L;
      const real momz_L = dens_L*velz_L;

      const real momx_R = dens_R*velx_R;
      const real momy_R = dens_R*vely_R;
      const real momz_R = dens_R*velz_R;
      //flop 6
      const real chalf = 0.5;

      const real B2_L   = sqr(Bx)     + sqr(By_L)   + sqr(Bz_L);
      const real v2_L   = sqr(velx_L) + sqr(vely_L) + sqr(velz_L);
      const real etot_L = ethm_L      + chalf*(dens_L*v2_L + B2_L);
      //flop 14
      const real B2_R   = sqr(Bx)     + sqr(By_R)   + sqr(Bz_R);
      const real v2_R   = sqr(velx_R) + sqr(vely_R) + sqr(velz_R);
      const real etot_R = ethm_R      + chalf*(dens_R*v2_R + B2_R);
      //flop 14
      const real gpl  = gamma_gas * pres_L;
      const real gpr  = gamma_gas * pres_R;
      const real gpbl = gpl + B2_L;
      const real gpbr = gpr + B2_R;
      //flop 4
      //flop//flop//flop//flop//flop

      const real cfl2  = gpbl / dens_L;
      const real cfr2  = gpbr / dens_R;
#if 1
      const real cfmax = std::sqrt(std::max(cfl2, cfr2));
#else
      const real cfmax = std::max(std::sqrt(std::max(cfl2, cfr2)),
          std::max(std::abs(velx_L), std::abs(velx_R)));
#endif

      //flop 6 + 1*10 + 2*10
      const real S_L = std::min(velx_L, velx_R) - cfmax;
      const real S_R = std::max(velx_L, velx_R) + cfmax;

      //flop 4
      //flop//flop//flop//flop//flop

      const real pT_L = pres_L + chalf * B2_L;
      const real pT_R = pres_R + chalf * B2_R;
      //flop 6
      const real iSM = 1.0/((S_R - velx_R)*dens_R - (S_L - velx_L)*dens_L);
      const real S_M  =   iSM * ((S_R - velx_R)*momx_R - (S_L - velx_L)*momx_L - pT_R + pT_L);

      //flop 13 + 1*10
      const real ipTs = 1.0/((S_R - velx_R)*dens_R - (S_L - velx_L)*dens_L);
      const real pT_s = ipTs * ((S_R - velx_R)*dens_R*pT_L - (S_L - velx_L)*dens_L*pT_R +
          dens_L*dens_R*(S_R - velx_R)*(S_L - velx_L)*(velx_R - velx_L));
      //flop 22 + 1*10
      const real velx_L_s  = S_M;
      const real velx_L_ss = S_M;
      const real velx_R_s  = S_M;
      const real velx_R_ss = S_M;
      const real B2x       = Bx*Bx;
      //flop 24
      const real iSLmSM   = 1.0/(S_L - S_M);
      const real iSRmSM   = 1.0/(S_R - S_M);
      const real dens_L_s = dens_L * (S_L - velx_L) * iSLmSM;
      const real dens_R_s = dens_R * (S_R - velx_R) * iSRmSM;
      const real divL     = dens_L * (S_L - velx_L)*(S_L - S_M) - B2x;
      const real divR     = dens_R * (S_R - velx_R)*(S_R - S_M) - B2x;
      const real idivL    = (divL != 0.0) ? 1.0/divL : 0.0;
      const real idivR    = (divR != 0.0) ? 1.0/divR : 0.0;
      const real vely_L_s = vely_L - Bx*By_L*(S_M - velx_L) * idivL;
      const real velz_L_s = velz_L - Bx*Bz_L*(S_M - velx_L) * idivL;
      const real   By_L_s = By_L * (dens_L*sqr(S_L - velx_L) - B2x) * idivL;
      const real   Bz_L_s = Bz_L * (dens_L*sqr(S_L - velx_L) - B2x) * idivL;
      //flop 42 + 2*10 + 2*20
      const real vely_R_s = vely_R - Bx*By_R*(S_M - velx_R) * idivR;
      const real velz_R_s = velz_R - Bx*Bz_R*(S_M - velx_R) * idivR;
      const real   By_R_s = By_R * (dens_R*sqr(S_R - velx_R) - B2x) * idivR;
      const real   Bz_R_s = Bz_R * (dens_R*sqr(S_R - velx_R) - B2x) * idivR;
      //flop 22
      const real   vB_L   = velx_L  *Bx + vely_L  *By_L   + velz_L  *Bz_L;
      const real   vB_L_s = velx_L_s*Bx + vely_L_s*By_L_s + velz_L_s*Bz_L_s;
      const real etot_L_s = ((S_L - velx_L)*etot_L - pT_L*velx_L + pT_s*S_M + Bx*(vB_L - vB_L_s)) * iSLmSM;
      //flop 20
      const real   vB_R   = velx_R  *Bx + vely_R  *By_R   + velz_R  *Bz_R;
      const real   vB_R_s = velx_R_s*Bx + vely_R_s*By_R_s + velz_R_s*Bz_R_s;
      const real etot_R_s = ((S_R - velx_R)*etot_R - pT_R*velx_R + pT_s*S_M + Bx*(vB_R - vB_R_s)) * iSRmSM;
      //flop 20
      const real dens_L_ss = dens_L_s;
      const real dens_R_ss = dens_R_s;
      const real sDens_L_s = std::sqrt(dens_L_s);
      const real sDens_R_s = std::sqrt(dens_R_s);
      //flop 2*10
      const real    S_L_s  = S_M - std::abs(Bx/sDens_L_s);
      const real    S_R_s  = S_M + std::abs(Bx/sDens_R_s);
      //flop 6 + 2*10
      const real idsqroot  = 1.0/(sDens_L_s + sDens_R_s);
      const real  vely_ss = idsqroot*(sDens_L_s*vely_L_s + sDens_R_s*vely_R_s + (By_R_s - By_L_s)*signBx);
      const real  velz_ss = idsqroot*(sDens_L_s*velz_L_s + sDens_R_s*velz_R_s + (Bz_R_s - Bz_L_s)*signBx);
      //flop 15 + 1*10
      const real By_ss = idsqroot*(sDens_L_s*By_R_s + sDens_R_s*By_L_s + sDens_L_s*sDens_R_s*(vely_R_s - vely_L_s)*signBx);
      const real Bz_ss = idsqroot*(sDens_L_s*Bz_R_s + sDens_R_s*Bz_L_s + sDens_L_s*sDens_R_s*(velz_R_s - velz_L_s)*signBx);
      //flop 18
      const real vely_L_ss = vely_ss;
      const real velz_L_ss = velz_ss;
      const real   By_L_ss = By_ss;
      const real   Bz_L_ss = Bz_ss;

      const real vely_R_ss = vely_ss;
      const real velz_R_ss = velz_ss;
      const real   By_R_ss = By_ss;
      const real   Bz_R_ss = Bz_ss;

      const real vB_L_ss   = velx_L_ss*Bx + vely_L_ss*By_L_ss + velz_L_ss*Bz_L_ss;
      const real etot_L_ss = etot_L_s - sDens_L_s*(vB_L_s - vB_L_ss)*signBx;
      //flop 8
      const real vB_R_ss   = velx_R_ss*Bx + vely_R_ss*By_R_ss + velz_R_ss*Bz_R_ss;
      const real etot_R_ss = etot_R_s + sDens_R_s*(vB_R_s - vB_R_ss)*signBx;
      //flop 8
      const real Fdens_L = dens_L*velx_L;
      const real Fmomx_L = momx_L*velx_L + pT_L - B2x;
      const real Fmomy_L = momy_L*velx_L        - Bx*By_L;
      const real Fmomz_L = momz_L*velx_L        - Bx*Bz_L;
      const real Fetot_L = etot_L*velx_L + pT_L*velx_L - Bx*vB_L; 
      //flop 15
      const real Fdens_R = dens_R*velx_R;
      const real Fmomx_R = momx_R*velx_R + pT_R - B2x;
      const real Fmomy_R = momy_R*velx_R        - Bx*By_R;
      const real Fmomz_R = momz_R*velx_R        - Bx*Bz_R;
      const real Fetot_R = etot_R*velx_R + pT_R*velx_R - Bx*vB_R;
      //flop 15
      const real momx_L_s  = dens_L_s *velx_L_s;
      const real momy_L_s  = dens_L_s *vely_L_s;
      const real momz_L_s  = dens_L_s *velz_L_s;

      const real momx_L_ss = dens_L_ss*velx_L_ss;
      const real momy_L_ss = dens_L_ss*vely_L_ss;
      const real momz_L_ss = dens_L_ss*velz_L_ss;

      const real momx_R_s  = dens_R_s *velx_R_s;
      const real momy_R_s  = dens_R_s *vely_R_s;
      const real momz_R_s  = dens_R_s *velz_R_s;

      const real momx_R_ss = dens_R_ss*velx_R_ss;
      const real momy_R_ss = dens_R_ss*vely_R_ss;
      const real momz_R_ss = dens_R_ss*velz_R_ss;

      const real Fby_L  = By_L*velx_L - Bx * vely_L;
      const real Fbz_L  = Bz_L*velx_L - Bx * velz_L;

      const real Fby_R  = By_R*velx_R - Bx * vely_R;
      const real Fbz_R  = Bz_R*velx_R - Bx * velz_R;
      //flop 24
      const bool w_lt_SL  = w <= S_L;
      const bool w_le_SLs = w <= S_L_s;
      const bool w_le_SM  = w <= S_M;
      const bool w_le_SRs = w <= S_R_s;
      const bool w_le_SR  = w <= S_R;
      //flop 5
      const real fdens = (w_le_SM) ? Fdens_L : Fdens_R;
      const real fetot = (w_le_SM) ? Fetot_L : Fetot_R;
      const real fmomx = (w_le_SM) ? Fmomx_L : Fmomx_R;
      const real fmomy = (w_le_SM) ? Fmomy_L : Fmomy_R;
      const real fmomz = (w_le_SM) ? Fmomz_L : Fmomz_R;
      const real fby   = (w_le_SM) ? Fby_L   :   Fby_R;
      const real fbz   = (w_le_SM) ? Fbz_L   :   Fbz_R;
      //flop 7
      const real cnull(0.0);

      const real a = 
        (w_le_SLs ? 	cnull : 
         (w_le_SM ? S_L_s - w : 
          (w_le_SRs ? S_R_s - w :	cnull)
         ));

      //flop 5
      const real b  = 
        (w_lt_SL ? cnull :
         (w_le_SLs ? +S_L - w :
          (w_le_SM  ? -S_L_s + S_L :
           (w_le_SRs ? -S_R_s + S_R :
            (w_le_SR  ? +S_R - w : cnull)
           ))));

      //flop 11
      const real c = 
        (w_lt_SL ? -w : 
         (w_le_SLs ? -S_L : 
          (w_le_SM  ? -S_L : 
           (w_le_SRs ? -S_R : 
            (w_le_SR  ?	-S_R : -w)
           ))));

      //flop 11
      const real dens    = (w_le_SM ? dens_L :    dens_R   );
      const real dens_s  = (w_le_SM ? dens_L_s :  dens_R_s );
      const real dens_ss = (w_le_SM ? dens_L_ss : dens_R_ss);

      const real etot    = (w_le_SM ? etot_L :    etot_R   );
      const real etot_s  = (w_le_SM ? etot_L_s :  etot_R_s );
      const real etot_ss = (w_le_SM ? etot_L_ss : etot_R_ss);

      const real momx    = (w_le_SM ? momx_L :    momx_R   );
      const real momx_s  = (w_le_SM ? momx_L_s :  momx_R_s );
      const real momx_ss = (w_le_SM ? momx_L_ss : momx_R_ss);

      const real momy    = (w_le_SM ? momy_L :    momy_R   );
      const real momy_s  = (w_le_SM ? momy_L_s :  momy_R_s );
      const real momy_ss = (w_le_SM ? momy_L_ss : momy_R_ss);

      const real momz    = (w_le_SM ? momz_L :    momz_R   );
      const real momz_s  = (w_le_SM ? momz_L_s :  momz_R_s );
      const real momz_ss = (w_le_SM ? momz_L_ss : momz_R_ss);

      const real by    = (w_le_SM ? By_L :    By_R   );
      const real by_s  = (w_le_SM ? By_L_s :  By_R_s );
      const real by_ss = (w_le_SM ? By_L_ss : By_R_ss);

      const real bz    = (w_le_SM ? Bz_L :    Bz_R   );
      const real bz_s  = (w_le_SM ? Bz_L_s :  Bz_R_s );
      const real bz_ss = (w_le_SM ? Bz_L_ss : Bz_R_ss);

      //flop 21
      flux[Fluid::MASS] = fdens + a * dens_ss + b * dens_s + c * dens;
      flux[Fluid::ENER] = fetot + a * etot_ss + b * etot_s + c * etot;
      flux[Fluid::MOMX] = fmomx + a * momx_ss + b * momx_s + c * momx;
      flux[Fluid::MOMY] = fmomy + a * momy_ss + b * momy_s + c * momy;
      flux[Fluid::MOMZ] = fmomz + a * momz_ss + b * momz_s + c * momz;
      flux[Fluid::WBY ] = fby   + a * by_ss   + b * by_s   + c * by;
      flux[Fluid::WBZ ] = fbz   + a * bz_ss   + b * bz_s   + c * bz;
      flux[Fluid::WBX ] = -w * Bx;

    }

    //flop 45
    asm("#RP_END");
  }

  void System::riemannSolver_HLLE(
      Fluid &flux,
      const real Bx,
      const real w,
      const real dens_L, const real pres_L, const real ethm_L,
      const real velx_L, const real vely_L, const real velz_L, const real By_L, const real Bz_L,
      const real dens_R, const real pres_R, const real ethm_R,
      const real velx_R, const real vely_R, const real velz_R, const real By_R, const real Bz_R
      ) const
  {
    asm("#RP_BEG");

    {

      const real momx_L = dens_L*velx_L;
      const real momy_L = dens_L*vely_L;
      const real momz_L = dens_L*velz_L;

      const real momx_R = dens_R*velx_R;
      const real momy_R = dens_R*vely_R;
      const real momz_R = dens_R*velz_R;

      const real chalf = real(0.5);

      const real B2_L   = sqr(Bx)     + sqr(By_L)   + sqr(Bz_L);
      const real v2_L   = sqr(velx_L) + sqr(vely_L) + sqr(velz_L);
      const real etot_L = ethm_L      + chalf*(dens_L*v2_L + B2_L);
      //flop 14
      const real B2_R   = sqr(Bx)     + sqr(By_R)   + sqr(Bz_R);
      const real v2_R   = sqr(velx_R) + sqr(vely_R) + sqr(velz_R);
      const real etot_R = ethm_R      + chalf*(dens_R*v2_R + B2_R);

      const real pT_L = pres_L + chalf * B2_L;
      const real pT_R = pres_R + chalf * B2_R;

#if 0
      const real gpl  = real(gamma_gas) * pres_L;
      const real gpr  = real(gamma_gas) * pres_R;
      const real gpbl = gpl + B2_L;
      const real gpbr = gpr + B2_R;

      const real cfl2  = gpbl/ dens_L;
      const real cfr2  = gpbr/ dens_R;
      const real cfmax = std::max(std::sqrt(std::max(cfl2, cfr2)),
          std::max(std::abs(velx_L), std::abs(velx_R)));
#else
      const real cf2   = (std::max(pres_L, pres_R) + std::max(B2_L, B2_R)*0.5)/std::min(dens_L, dens_R);
      const real cfmax = std::sqrt(cf2);
#endif


      //flop 6 + 1*10 + 2*10
#if 0
      const real S_L = std::min(velx_L, velx_R) - cfmax;
      const real S_R = std::max(velx_L, velx_R) + cfmax;
#else
      const real S_R = std::max(std::abs(velx_L), std::abs(velx_R)) + cfmax;
      const real S_L = -S_R;
#endif

      real vB_L = velx_L*Bx + vely_L*By_L   + velz_L*Bz_L;
      real vB_R = velx_R*Bx + vely_R*By_R   + velz_R*Bz_R;

      real Fdens_L = dens_L*velx_L;
      real Fmomx_L = momx_L*velx_L + pT_L - Bx*Bx;
      real Fmomy_L = momy_L*velx_L        - Bx*By_L; 
      real Fmomz_L = momz_L*velx_L        - Bx*Bz_L;
      real Fetot_L = etot_L*velx_L + pT_L*velx_L - Bx*vB_L;
      real Fby_L   = By_L  *velx_L - Bx*vely_L;
      real Fbz_L   = Bz_L  *velx_L - Bx*velz_L;

      real Fdens_R = dens_R*velx_R;
      real Fmomx_R = momx_R*velx_R + pT_R - Bx*Bx;
      real Fmomy_R = momy_R*velx_R        - Bx*By_R; 
      real Fmomz_R = momz_R*velx_R        - Bx*Bz_R; 
      real Fetot_R = etot_R*velx_R + pT_R*velx_R - Bx*vB_R;
      real Fby_R   = By_R  *velx_R        - Bx*vely_R;
      real Fbz_R   = Bz_R  *velx_R        - Bx*velz_R;


      real U_dens = (S_R*dens_R - S_L*dens_L + Fdens_L - Fdens_R)/(S_R - S_L);
      real U_momx = (S_R*momx_R - S_L*momx_L + Fmomx_L - Fmomx_R)/(S_R - S_L);
      real U_momy = (S_R*momy_R - S_L*momy_L + Fmomy_L - Fmomy_R)/(S_R - S_L);
      real U_momz = (S_R*momz_R - S_L*momz_L + Fmomz_L - Fmomz_R)/(S_R - S_L);
      real U_etot = (S_R*etot_R - S_L*etot_L + Fetot_L - Fetot_R)/(S_R - S_L);

      flux[Fluid::MASS] = (S_R*Fdens_L - S_L*Fdens_R + S_L*S_R*(dens_R - dens_L))/(S_R - S_L);
      flux[Fluid::MOMX] = (S_R*Fmomx_L - S_L*Fmomx_R + S_L*S_R*(momx_R - momx_L))/(S_R - S_L);
      flux[Fluid::MOMY] = (S_R*Fmomy_L - S_L*Fmomy_R + S_L*S_R*(momy_R - momy_L))/(S_R - S_L);
      flux[Fluid::MOMZ] = (S_R*Fmomz_L - S_L*Fmomz_R + S_L*S_R*(momz_R - momz_L))/(S_R - S_L);
      flux[Fluid::ENER] = (S_R*Fetot_L - S_L*Fetot_R + S_L*S_R*(etot_R - etot_L))/(S_R - S_L);

      real U_by = (S_R*By_R - S_L*By_L + Fby_L - Fby_R)/(S_R - S_L);
      real U_bz = (S_R*Bz_R - S_L*Bz_L + Fbz_L - Fbz_R)/(S_R - S_L);

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

    //flop 45
    asm("#RP_END");
  }
}
