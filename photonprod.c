/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : photonprod.c                                   */
/*                                                                           */
/* Created:       2016/01/30 (JLe)                                           */
/* Last modified: 2020/06/05 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Produce photons from neutron collisions                      */
/*                                                                           */
/* Comments: - Probably does not work with implicit capture                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PhotonProd:"

/*****************************************************************************/

void PhotonProd(long part, long mat, long rea0, double x, double y, double z,
                double u0, double v0, double w0, double E0, double wgt,
                double t, long id)
{
  long nuc, mode, loc0, loc1, ptr, ntot, n, m, rea, new, mt, loc;
  double elaxs, totxs, inlxs, prodxs, yld, xs, mu, E, f, u, v, w, p, EGP, EGD;
  double wmin, wgt0;

  /* Check mode */

  if ((mode = (long)RDB[DATA_PHOTON_PRODUCTION]) == NO)
    return;

  /* Check particle type */

  if ((long)RDB[part + PARTICLE_TYPE] != PARTICLE_TYPE_NEUTRON)
    return;

  /* Skip inactive */

  if ((long)RDB[DATA_STOP_AFTER_PLOT] != STOP_AFTER_PLOT_TRACKS)
    if (RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP])
      return;

  /* Check reaction pointer */

  CheckPointer(FUNCTION_NAME, "(rea0)", DATA_ARRAY, rea0);

  /* Pointer to nuclide (parent pointer for S(a,b) reactions) */

  nuc = -1;
  if ((nuc = (long)RDB[rea0 + REACTION_PTR_PARENT_NUCLIDE]) < VALID_PTR)
    nuc = (long)RDB[rea0 + REACTION_PTR_NUCLIDE];

  /* Check */

  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

  /* Pointer to total */

  if ((loc0 = (long)RDB[nuc + NUCLIDE_PTR_PHOTON_PROD]) < VALID_PTR)
    return;

  /* Get nuclide total and elastic cross sections */

  ptr = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  totxs = MicroXS(ptr, E0, id);

  /* Check total xs */

  if (totxs == 0.0)
    return;

  /* Get nuclide photon production xs */

  ptr = (long)RDB[nuc + NUCLIDE_PTR_PHOTPRODXS];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  prodxs = MicroXS(ptr, E0, id);

  /* Check */

  if (prodxs == 0.0)
    return;

  /* Check ures mode */

  if ((long)RDB[nuc + NUCLIDE_URES_SAMPLING] == YES)
    {
      /* Loop over reactions */

      loc1 = NextItem(loc0);
      while (loc1 > VALID_PTR)
        {
          /* Check pointer to ures data */

          if ((rea = (long)RDB[loc1 + PHOTON_PROD_PTR_URES_REA]) > VALID_PTR)
            {
              /* Get factor */

              f = UresFactor(rea, E0, id);
              CheckValue(FUNCTION_NAME, "f", "", f, 0.0, 1E+6);

              /* Get cross section */

              ptr = (long)RDB[loc1 + PHOTON_PROD_PTR_PRODXS];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              xs = MicroXS(ptr, E0, id);

              /* Adjust */

              prodxs = prodxs + (f - 1.0)*xs;
            }

          /* Next */

          loc1 = NextItem(loc1);
        }
    }

  /* Check if production cross section is very low (limited numerical  */
  /* precision, this happens with ures probability table sampling, and */
  /* it may cause reaction sampling to fail later on) */

  /*
  if (prodxs < 1E-6)
    return;
  */

  if (prodxs < ZERO)
    return;

  /* Avoid compiler warning */

  ntot = -1;

  /* Russian roulette parameters */

  wmin = RDB[DATA_PHOTON_IMPL_WMIN];
  p = -1.0;

  /* Check mode */

  if (mode == PHOTON_PROD_ANA)
    {
      /* Analog mode, get elastic cross section */

      ptr = (long)RDB[nuc + NUCLIDE_PTR_ELAXS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      elaxs = MicroXS(ptr, E0, id);

      /* Reject elastic reactions */

      if (RandF(id) < elaxs/totxs)
        return;

      /* Calculate inelastic cross section */

      inlxs = totxs - elaxs;
      CheckValue(FUNCTION_NAME, "inlxs", "", inlxs, ZERO, INFTY);

      /* Calculate yield */

      yld = prodxs/inlxs;

      /* Sample number of emitted photons */

      ntot = (long)yld;

      if (RandF(id) < yld - (double)ntot)
        ntot++;

      /* Check multiplicity */

      if (ntot == 0)
        return;

      /* Set weight */

      wgt0 = wgt;
    }
  else if (mode == PHOTON_PROD_IMP)
    {
      /* Implicit mode, calculate weight */

      wgt = wgt*prodxs/totxs;

      /* Calculate number of photons */

      if ((ntot = (long)(wgt/wmin) + 1) > (long)RDB[DATA_PHOTON_IMPL_NMAX])
        ntot = (long)RDB[DATA_PHOTON_IMPL_NMAX];

      /* Calculate weight of single particle */

      wgt0 = wgt/((double)ntot);

      /* Calculate survival probability */

      p = wgt0/wmin;
    }

  /* Loop over photons */

  for (n = 0; n < ntot; n++)
    {
      /* Check mode and weight */

      if ((mode == PHOTON_PROD_IMP) && (wgt0 < wmin))
        {
          /* Russian roulette */

          if (RandF(id) < p)
            {
              /* Increase weight */

              wgt = wmin;
            }
          else
            {
              /* Skip particle */

              continue;
            }

        }
      else
        {
          /* Use original weight */

          wgt = wgt0;
        }

      /* Re-sampling loop (12024.03c default-kirjastossa feilaa) */

      for (m = 0; m < 10; m++)
        {
          /* Sample fraction of total production xs */

          xs = RandF(id)*prodxs;

          /* Loop over reactions */

          loc1 = NextItem(loc0);
          while (loc1 > VALID_PTR)
            {
              /* Check ures data */

              if ((rea = (long)RDB[loc1 + PHOTON_PROD_PTR_URES_REA]) >
                  VALID_PTR)
                {
                  /* Get factor */

                  f = UresFactor(rea, E0, id);
                  CheckValue(FUNCTION_NAME, "f", "", f, 0.0, 1E+6);
                }
              else
                f = 1.0;

              /* Subtract production xs */

              ptr = (long)RDB[loc1 + PHOTON_PROD_PTR_PRODXS];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              xs = xs - f*MicroXS(ptr, E0, id);

              /* Check */

              if (xs <= 0.0)
                break;

              /* Next reaction */

              loc1 = NextItem(loc1);
            }

          /* Check sample */

          if (loc1 > VALID_PTR)
            break;
        }

      /* Score sampled */

      ptr = (long)RDB[RES_PHOTON_SAMPLING_FAIL];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(1.0, 1.0, ptr, id, 1);

      /* Check count */

      if (m == 10)
        {
          /* Score failed */

          AddBuf1D(1.0, 1.0, ptr, id, 0);

#ifdef DEBUG

          Warn(FUNCTION_NAME, "Failed to sample reaction (E = %E, xs = %E %E)",
              E0, prodxs, totxs);

#endif

          /* Cycle loop */

          continue;
        }

      /* Sample direction (pointer is null if isotropic) */

      ptr = (long)RDB[loc1 + PHOTON_PROD_PTR_ANG];
      mu = SampleMu(-1, ptr, E0, NULL, NULL, id);

      /* Rotate direction cosines */

      u = u0;
      v = v0;
      w = w0;

      AziRot(mu, &u, &v, &w, id);

      /* Sample energy */

      ptr = (long)RDB[loc1 + PHOTON_PROD_PTR_ERG];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      SampleENDFLaw(-1, ptr, E0, &E, NULL, id);
      CheckValue(FUNCTION_NAME, "E", "", E, 0.0, INFTY);

      /* Check failure */

      if (E == 0.0)
        {
          /* Score sampled */

          ptr = (long)RDB[RES_PHOTON_SAMPLING_FAIL];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(1.0, 1.0, ptr, id, 0);

          /* Cycle loop */

          continue;
        }

      /* Store analog production rate */

      ptr = (long)RDB[RES_TOT_PHOTON_PRODRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(1.0, wgt, ptr, id, 1);

      /***********************************************************************/

      /***** Put photon in que ***********************************************/

      /* Get particle from stack */

      new = FromStack(PARTICLE_TYPE_GAMMA, id);

      /* Put variables */

      WDB[new + PARTICLE_X] = x;
      WDB[new + PARTICLE_Y] = y;
      WDB[new + PARTICLE_Z] = z;

      WDB[new + PARTICLE_U] = u;
      WDB[new + PARTICLE_V] = v;
      WDB[new + PARTICLE_W] = w;

      WDB[new + PARTICLE_E] = E;
      WDB[new + PARTICLE_WGT] = wgt;
      WDB[new + PARTICLE_T0] = t;
      WDB[new + PARTICLE_T] = t;
      WDB[new + PARTICLE_TD] = 0.0;
      WDB[new + PARTICLE_TT] = 0.0;
      WDB[new + PARTICLE_COL_IDX] = 0.0;

      /* Put material */

      WDB[new + PARTICLE_PTR_MAT] = (double)mat;

      /* Copy data from parent neutron */

      WDB[new + PARTICLE_RNG_IDX] = RDB[part + PARTICLE_RNG_IDX];
      WDB[new + PARTICLE_HISTORY_IDX] = RDB[part + PARTICLE_HISTORY_IDX];
      WDB[new + PARTICLE_DET_FLAGS] = RDB[part + PARTICLE_DET_FLAGS];
      WDB[new + PARTICLE_MPI_ID] = RDB[part + PARTICLE_MPI_ID];

      /* Set photon emission type */

      WDB[new + PARTICLE_PHOTON_TYPE] = (double)PHOTON_TYPE_NREA;

      /* Set energy deposition correction factor */

      if (((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT) &&
          ((long)RDB[DATA_EDEP_MODE] == EDEP_MODE_NEUTRON_PHOTON))
        {
          mt = (long)RDB[loc1 + PHOTON_PROD_MT];

          /* Separate fission gammas from other gammas */

          if (mt < 18001 || (mt > 22000 && mt < 38001) || mt > 39000)
            {
              if ((long)RDB[DATA_EDEP_KEFF_CORR])
                WDB[new + PARTICLE_EDEP_RENORM] = RDB[DATA_CYCLE_KEFF];
              else
                WDB[new + PARTICLE_EDEP_RENORM] = 1.0;
            }
          else
            {
              if ((long)RDB[DATA_EDEP_DELAYED] == YES &&
                  (long)RDB[DATA_EDEP_LOCAL_EGD] == NO &&
                  ((loc = (long)RDB[nuc + NUCLIDE_PTR_FISSE_DATA]) > VALID_PTR))
                {
                  /* Calculate prompt and delayed fission gamma energies */

                  EGP = FissEComp(loc, 0, FISSE_COMP_EGP, E0, id);
                  EGD = FissEComp(loc, 0, FISSE_COMP_EGD, E0, id);

                  /* Calculate correction factor */

                  WDB[new + PARTICLE_EDEP_RENORM] = (EGP + EGD)/EGP;
                }
              else
                WDB[new + PARTICLE_EDEP_RENORM] = 1.0;
            }
        }

      /* Check energy and weight cut-off */

      if (E > RDB[mat + MATERIAL_PHOTON_ECUT])
        {
          /* Score balance */

          ptr = (long)RDB[RES_G_BALA_SRC];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf(1.0, 1.0, ptr, id, -1, BALA_G_SRC_NREA, 0);
          AddBuf(wgt, 1.0, ptr, id, -1, BALA_G_SRC_NREA, 1);
          AddBuf(wgt*E, 1.0, ptr, id, -1, BALA_G_SRC_NREA, 2);

          /* Apply weight window */

          if (WeightWindow(-1, new, PARTICLE_TYPE_GAMMA, x, y, z, u, v, w,
                           E, &wgt, t, WWMESH_COLL, id) != TRACK_END_WCUT)
            {
              /* Put weight (may have been changed by russian roulette) */

              WDB[new + PARTICLE_WGT] = wgt;

              /* Check if source is written in file */

              if ((ptr = (long)RDB[DATA_PTR_NGAMMA_SRC_DET]) > VALID_PTR)
                {
                  /* Write source file */

                  WriteSourceFile(ptr, x, y, z, u, v, w, E, wgt, t, -1.0, id);

                  /* Put particle in que or stack */

                  if ((long)RDB[DATA_NGAMMA_SRC_SIM] == YES)
                    {
                      /* Score response matrix source term */

                      ScoreRMXSrc(new, x, y, z, E, wgt, id);

                      /* Photon is transported, put in que */

                      ToQue(new, id);
                    }
                  else
                    {
                      /* Not transported put to stack */

                      ToStack(new, id);

                      /* Score cut-off */

                      ptr = (long)RDB[RES_TOT_PHOTON_CUTRATE];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                      AddBuf1D(1.0, wgt, ptr, id, 0);

                      /* Score particle balance */

                      ptr = (long)RDB[RES_G_BALA_LOSS];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                      AddBuf(1.0, 1.0, ptr, id, -1, BALA_G_LOSS_CUT, 0);
                      AddBuf(wgt, 1.0, ptr, id, -1, BALA_G_LOSS_CUT, 1);
                    }
                }
              else
                {
                  /* Score response matrix source term */

                  ScoreRMXSrc(new, x, y, z, E, wgt, id);

                  /* Put particle in que */

                  ToQue(new, id);
                }
            }
        }
      else
        {
          /* Below cut-off, score photon heat detector */

          ScorePhotonHeat(new, mat, E, x, y, z, E, t, wgt, id);

          /* Score pulse-height detector */

          PulseDet(new, mat, E, x, y, z, wgt, id);

          /* Put photon back in stack */

          ToStack(new, id);
        }

      /***********************************************************************/
    }
}

/*****************************************************************************/
