/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : collision.c                                    */
/*                                                                           */
/* Created:       2011/02/28 (JLe)                                           */
/* Last modified: 2019/10/21 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Samples reaction, updates coordinates and scores estimators  */
/*                                                                           */
/* Comments: - Weights: wgt0 initial particle weight                         */
/*                      wgt1 particle weight after implicit capture          */
/*                      wgt2 particle weight after reactions                 */
/*                                                                           */
/*           - wgt1 = wgt0 in analog capture.                                */
/*                                                                           */
/*           - Implicit fission and (n,xn) may change weight, otherwise      */
/*             wgt2 = wgt1.                                                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "Collision:"

/*****************************************************************************/

long Collision(long mat, long part, double x, double y, double z, double *u,
               double *v, double *w, double *E, double *wgt, double t, long id)
{
  long type, rea, ptr, mt, scatt, icapt;
  double totxs, absxs, E0, u0, v0, w0, mu, wgt0, wgt1, wgt2, dE;

  /* Check input parameters */

  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);
  CheckValue(FUNCTION_NAME, "x", "", x, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "y", "", y, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "z", "", z, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "cos", "", *u**u+*v**v+*w**w - 1.0, -1E-5, 1E-5);
  CheckValue(FUNCTION_NAME, "E", "", *E, ZERO, INFTY);
  CheckValue(FUNCTION_NAME, "t", "", t, ZERO, INFTY);
  CheckValue(FUNCTION_NAME, "wgt", "", *wgt, ZERO, INFTY);

  /* Get particle type */

  type = (long)RDB[part + PARTICLE_TYPE];

  /* Get pointer to total xs */

  if (type == PARTICLE_TYPE_NEUTRON)
    ptr = (long)RDB[mat + MATERIAL_PTR_TOTXS];
  else
    ptr = (long)RDB[mat + MATERIAL_PTR_TOTPHOTXS];

  /* Get implicit capture flag */

  icapt = (long)RDB[DATA_OPT_IMPL_CAPT];

  /* Get initial weight and reset others */

  wgt0 = *wgt;
  wgt1 = -1.0;
  wgt2 = -1.0;

  /* Remember values before collision */

  E0 = *E;
  u0 = *u;
  v0 = *v;
  w0 = *w;

  /* Reset change in change in particle energy */

  dE = E0;

  /* Weight reduction by implicit capture */

  if ((icapt == YES) && (type == PARTICLE_TYPE_NEUTRON))
    {
      /* Get total xs */

      totxs = TotXS(mat, type, *E, id);

      /* Get material total absorption xs (may be zero for He) */

      if ((ptr = (long)RDB[mat + MATERIAL_PTR_ABSXS]) > VALID_PTR)
        absxs = MacroXS(ptr, *E, id);
      else
        absxs = 0.0;

      /* Score capture reaction */

      if (absxs > 0.0)
        ScoreCapture(mat, -1, *E, wgt0*absxs/totxs, id);

      /* Calculate weight reduction */

      wgt1 = wgt0*(1.0 - absxs/totxs);
    }
  else
    wgt1 = wgt0;

  /* Sample reaction */

  if ((rea = SampleReaction(mat, type, *E, wgt1, id)) < VALID_PTR)
    {
      /* Sample rejected, set final weight */

      *wgt = wgt1;

      /* Score efficiency */

      ptr = (long)RDB[RES_REA_SAMPLING_EFF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(1.0, 1.0, ptr, id, 4 - type);

      /* Return virtual */

      return TRACK_END_VIRT;
    }
  else
    {
      /* Score efficiency */

      ptr = (long)RDB[RES_REA_SAMPLING_EFF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(1.0, 1.0, ptr, id, 2 - type);
    }

  /* Update collision index */

  WDB[part + PARTICLE_COL_IDX] = RDB[part + PARTICLE_COL_IDX] + 1.0;

  /* Produce photons */

  PhotonProd(part, mat, rea, x, y, z, *u, *v, *w, *E, *wgt, t, id);

  /* Get reaction mt */

  mt = (long)RDB[rea + REACTION_MT];

  /* Tarkistetaan tää 9.8.2016 tehdyn muutoksen jäljiltä (JLe) */

  if (mt == 4)
    Die(FUNCTION_NAME, "WTF?");

  /* Reset scattering flag */

  scatt = NO;

  /* Check particle type */

  if (type == PARTICLE_TYPE_NEUTRON)
    {
      /***********************************************************************/

      /***** Neutron reactions ***********************************************/

      if ((long)RDB[DATA_SENS_MODE] != SENS_MODE_NONE)
        if (SensCollision(mat, rea, part, *E, id) == COLL_MISS)
          return TRACK_END_VIRT;

      /* Check reaction type */

      if (mt == 2)
        {
          /* Check sampling */

          if ((long)RDB[DATA_NPHYS_SAMPLE_SCATT] == NO)
            return TRACK_END_SCAT;

          /* Elastic scattering */

          ElasticScattering(mat, rea, part, E, u, v, w, id);

          /* Set scattering flag and calculate change in energy */

          scatt = YES;
          dE = dE - *E;

          /* Weight is preserved */

          wgt2 = wgt1;
        }
      else if ((mt == 1002) || (mt == 1004))
        {
          /* Check sampling */

          if ((long)RDB[DATA_NPHYS_SAMPLE_SCATT] == NO)
            return TRACK_END_SCAT;

          /* Scattering by S(a,b) laws */

          SabScattering(rea, E, u, v, w, id);

          /* Set scattering flag and calculate change in energy */

          scatt = YES;
          dE = dE - *E;

          /* Weight is preserved */

          wgt2 = wgt1;
        }
      else if ((mt == 2002) || (mt == 2004))
        {
          /* S(a,b) scattering with on-the-fly interpolation */

         /* Check sampling */

         if ((long)RDB[DATA_NPHYS_SAMPLE_SCATT] == NO)
           return TRACK_END_SCAT;

         /* Scattering by S(a,b) laws */

         OTFSabScattering(rea, E, u, v, w, id);

          /* Set scattering flag and calculate change in energy */

          scatt = YES;
          dE = dE - *E;

         /* Weight is preserved */

         wgt2 = wgt1;
        }
      else if (RDB[rea + REACTION_TY] == 0.0)
        {
          /* Capture */

          if (icapt == YES)
            {
              /* Not possible in implicit mode */

              Die(FUNCTION_NAME, "Capture reaction in implicit mode");
            }
          else if ((long)RDB[DATA_NPHYS_SAMPLE_CAPT] == NO)
            {
              /* Not sampled, return scattering */

              return TRACK_END_SCAT;
            }
          else
            {

              /* Score capture reaction */

              ScoreCapture(mat, rea, *E, wgt1, id);

              /* Put particle back in stack */

              ToStack(part, id);

              /* Exit subroutine */

              return TRACK_END_CAPT;
            }
        }
      else if (fabs(RDB[rea + REACTION_TY]) > 100.0)
        {
          /* Complex reaction */

          ComplexRea(rea, part, E, x, y, z, u, v, w, wgt1, &wgt2, t, &dE, id);

          /* Check weight */

          if (wgt2 > 0.0)
            {
              /* Set scattering */

              scatt = YES;
            }
          else
            {
              /* Score capture reaction */

              ScoreCapture(mat, rea, *E, wgt1, id);

              /* Put particle back in stack */

              ToStack(part, id);

              /* Exit subroutine */

              return TRACK_END_CAPT;
            }
        }
      else if (((mt > 17) && (mt < 22)) || (mt == 38))
        {
          /* Check sampling */

          if ((long)RDB[DATA_NPHYS_SAMPLE_FISS] == NO)
            {
              /* Tätä muutettiin 18.7.2013 / 2.1.15 sillai että   */
              /* readacefile.c:ssä fission TY laitetaan nollaksi, */
              /* eli koko listaa ei pitäisi luoda. */

              /* Check if capture is sampled */

              if ((long)RDB[DATA_NPHYS_SAMPLE_CAPT] == NO)
                return TRACK_END_SCAT;
              else
                {
                  /* Put particle back in stack */

                  ToStack(part, id);

                  /* Exit subroutine */

                  return TRACK_END_CAPT;
                }
            }

          /* Sample fission */

          Fission(mat, rea, part, E, t, x, y, z, u, v, w, wgt1, &wgt2, id);

          /* Put particle back in stack */

          ToStack(part, id);

          /* Exit subroutine */

          return TRACK_END_FISS;
        }
      else if ((mt > 50) && (mt < 91))
        {
          /* Check sampling */

          if ((long)RDB[DATA_NPHYS_SAMPLE_SCATT] == NO)
            return TRACK_END_SCAT;

          /* Inelastic level scattering */

          LevelScattering(rea, E, u, v, w, id);

          /* Set scattering flag and calculate change in energy */

          scatt = YES;
          dE = dE - *E;

          /* Weight is preserved */

          wgt2 = wgt1;
        }
      else if (RDB[rea + REACTION_WGT_F] > 1.0)
        {
          /* Check sampling */

          if ((long)RDB[DATA_NPHYS_SAMPLE_SCATT] == NO)
            return TRACK_END_SCAT;

          /* Multiplying scattering reaction */

          Nxn(rea, part, E, x, y, z, u, v, w, wgt1, &wgt2, t, &dE, id);

          /* Set scattering flag */

          scatt = YES;
        }
      else if (mt < 100)
        {
          /* Check sampling */

          if ((long)RDB[DATA_NPHYS_SAMPLE_SCATT] == NO)
            return TRACK_END_SCAT;

          /* Continuum single neutron reactions */

          InelasticScattering(rea, E, u, v, w, id);

          /* Set scattering flag and calculate change in energy */

          scatt = YES;
          dE = dE - *E;

          /* Weight is preserved */

          wgt2 = wgt1;
        }
      else
        {
          /* Unknown reaction mode */

          Die(FUNCTION_NAME, "Unknown reaction mode %ld sampled", mt);
        }

      /***********************************************************************/
    }
  else
    {
      /***** Photon reactions ************************************************/

      if (mt == 504)
        {
          /* Incoherent (Compton) scattering */

          ComptonScattering(mat, rea, part, E, x, y, z, u, v, w, *wgt, t, id);

          /* Set scattering flag and calculate change in energy */

          scatt = YES;
          dE = dE - *E;

          /* Switch fluorescence index flag off */

          if ((long)RDB[part + PARTICLE_PHOTON_TYPE] == PHOTON_TYPE_FLUOR)
            WDB[part + PARTICLE_PE_AR_IDX] = -1;
        }
      else if (mt == 502)
        {
          /* Coherent (Rayleigh) scattering */

          RayleighScattering(rea, *E, u, v, w, id);

          /* Set scattering flag and calculate change in energy */

          scatt = YES;
          dE = dE - *E;
        }
      else if (mt == 522)
        {
          /* Score capture rate */

          ptr = (long)RDB[RES_PHOTOELE_CAPT_RATE];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(1.0, *wgt, ptr, id, 0);

          /* Score particle balance */

          ptr = (long)RDB[RES_G_BALA_LOSS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf(1.0, 1.0, ptr, id, -1, BALA_G_LOSS_CAPT, 0);
          AddBuf(*wgt, 1.0, ptr, id, -1, BALA_G_LOSS_CAPT, 1);

          /* Photoelectric effect */

          Photoelectric(mat, rea, part, *E, x, y, z, *u, *v, *w, *wgt, t, id);

          /* Incident photon is killed */

          return TRACK_END_CAPT;
        }
      else if (mt == 516)
        {
          /* Score capture rate */

          ptr = (long)RDB[RES_PAIRPROD_CAPT_RATE];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(1.0, *wgt, ptr, id, 0);

          /* Score particle balance */

          ptr = (long)RDB[RES_G_BALA_LOSS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf(1.0, 1.0, ptr, id, -1, BALA_G_LOSS_CAPT, 0);
          AddBuf(*wgt, 1.0, ptr, id, -1, BALA_G_LOSS_CAPT, 1);

          /* Pair production */

          PairProduction(mat, rea, part, *E, x, y, z, *u, *v, *w, *wgt, t, id);

          /* Incident photon is killed */

          return TRACK_END_CAPT;
        }
      else
        Die(FUNCTION_NAME, "Invalid reaction mode %ld sampled", mt);

      /* Set weight */

      wgt2 = wgt1;

      /***********************************************************************/
    }

  /* Check final weight */

  if (wgt2 < 0.0)
    Die(FUNCTION_NAME, "Error in weight");

  /* Apply weight window */

  if (WeightWindow(-1, part, type, x, y, z, *u, *v, *w, *E, &wgt2, t,
                   WWMESH_COLL, id) == TRACK_END_WCUT)
    {
      /* Exit subroutine */

      return TRACK_END_WCUT;
    }

  /* Russian roulette */

  if (wgt2 < RDB[DATA_OPT_ROULETTE_W0])
    {
      if (RandF(id) < RDB[DATA_OPT_ROULETTE_P0])
        wgt2 = wgt2/RDB[DATA_OPT_ROULETTE_P0];
      else
        {
          /* Put particle back in stack */

          ToStack(part, id);

          /* Score cut-off */

          if (type == PARTICLE_TYPE_NEUTRON)
            ptr = (long)RDB[RES_TOT_NEUTRON_CUTRATE];
          else
            ptr = (long)RDB[RES_TOT_PHOTON_CUTRATE];

          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(1.0, wgt2, ptr, id, 0);

          /* Exit subroutine */

          return TRACK_END_WCUT;
        }
    }

  if (EneImportance(part, x, y, z, *u, *v, *w, E0, *E, &wgt2, t, id)
      == TRACK_END_WCUT)
    return TRACK_END_WCUT;

  /* Set final weight */

  *wgt = wgt2;

  /* Check energy, weight and cosines */

  CheckValue(FUNCTION_NAME, "E", "", *E, ZERO, INFTY);
  CheckValue(FUNCTION_NAME, "wgt", "", *wgt, ZERO, INFTY);
  CheckValue(FUNCTION_NAME, "r", "", *u**u + *v**v + *w**w - 1.0, -2E-5, 2E-5);

  /* Check with boundaries */

  if (type == PARTICLE_TYPE_NEUTRON)
    {
      /* Check cut-off (NOTE: this is a non-physical operation that */
      /* does not preserve the energy deposition, etc.) */

      if (*E < RDB[DATA_NEUTRON_ECUT])
        {
          /* Put particle back in stack */

          ToStack(part, id);

          /* Score cut-off */

          ptr = (long)RDB[RES_TOT_NEUTRON_CUTRATE];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(1.0, *wgt, ptr, id, 0);

          /* Score particle balance */

          ptr = (long)RDB[RES_N_BALA_LOSS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf(1.0, 1.0, ptr, id, -1, BALA_N_LOSS_CUT, 0);
          AddBuf(*wgt, 1.0, ptr, id, -1, BALA_N_LOSS_CUT, 1);

          /* Exit subroutine */

          return TRACK_END_ECUT;
        }

      /* Adjust neutron energy */

      if (*E < 1.0000001*RDB[DATA_NEUTRON_EMIN])
        *E = 1.000001*RDB[DATA_NEUTRON_EMIN];
      else if (*E > 0.999999*RDB[DATA_NEUTRON_EMAX])
        *E = 0.999999*RDB[DATA_NEUTRON_EMAX];

      /* Check scattering flag */

      if (scatt == YES)
        {
          /* Calculate scattering cosine */

          mu = *u*u0 + *v*v0 + *w*w0;

          /* Score scattering */

          ScoreScattering(mat, rea, mu, E0, *E, wgt1, wgt2, id);

          /* Score recoil detector */

          RecoilDet(mat, dE, x, y, z, u0, v0, w0, E0, t, wgt1, id);

          /* Score mesh */

          ScoreMesh(part, mat, 0.0, dE, x, y, z, E0, t, wgt1, 1.0, id);
        }

      /* Set time of thermalization */

      if ((RDB[part + PARTICLE_TT] == 0.0) && (*E < 0.625E-6))
        WDB[part + PARTICLE_TT] = t;
    }
  else
    {
      /* Do energy cut-off for photons and check the upper limit */

      if ((*E < RDB[DATA_PHOTON_EMIN]) ||
          (*E < RDB[mat + MATERIAL_PHOTON_ECUT]))
        {
          /* Score pulse-height detector */

          PulseDet(part, mat, *E, x, y, z, wgt2, id);

          /* Score photon heat detector */

          ScorePhotonHeat(part, mat, *E, x, y, z, E0, t, wgt2, id);

          /* Put particle back in stack */

          ToStack(part, id);

          /* Score cut-off */

          ptr = (long)RDB[RES_TOT_PHOTON_CUTRATE];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(1.0, *wgt, ptr, id, 0);

          /* Particle balance */

          ptr = (long)RDB[RES_G_BALA_LOSS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf(1.0, 1.0, ptr, id, -1, BALA_G_LOSS_CUT, 0);
          AddBuf(*wgt, 1.0, ptr, id, -1, BALA_G_LOSS_CUT, 1);

          /* Exit subroutine */

          return TRACK_END_ECUT;
        }
      else if ((*E > 0.999999*RDB[DATA_PHOTON_EMAX]) &&
               (*E <= RDB[DATA_PHOTON_EMAX]))
        {
          /* Adjust the energy */

          *E = 0.999999*RDB[DATA_PHOTON_EMAX];
        }
      else if (*E > RDB[DATA_PHOTON_EMAX])
        Die(FUNCTION_NAME, "Photon energy %E MeV above maximum %E MeV",
            *E, RDB[DATA_PHOTON_EMAX]);
    }

  /* Check that reaction was scattering */

  if (scatt == NO)
    Die(FUNCTION_NAME, "not a scattering reaction");

  /* Exit subroutine */

  return TRACK_END_SCAT;

  /***************************************************************************/
}

/*****************************************************************************/
