/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : score.c                                        */
/*                                                                           */
/* Created:       2011/03/03 (JLe)                                           */
/* Last modified: 2019/04/03 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Scores general parameters needed for every calculation       */
/*                                                                           */
/* Comments: - NOTE: "turhat" MacroXS() -kutsut vie tehokkuuden optimointi-  */
/*             moodi kakkoselta. Ei pitäisi kutsua täällä ellei välttämättä  */
/*             tarvita.                                                      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "Score:"

/*****************************************************************************/

void Score(long mat, long part, double flx, double x, double y, double z,
           double u, double v, double w, double E, double wgt, double t,
           double spd, double g, long id)
{
  long rea, ptr, prg, i, dng, n, norm;
  double val, fiss, tot, fissE, nsf, capt, ela, sprod, lambda, heat;

  /* Check input parameters */

  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);
  CheckValue(FUNCTION_NAME, "x", "", x, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "y", "", y, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "z", "", z, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "cos", "", u*u+v*v+w*w - 1.0, -1E-5, 1E-5);
  CheckValue(FUNCTION_NAME, "flx", "", flx, ZERO, INFTY);
  CheckValue(FUNCTION_NAME, "E", "", E, ZERO, INFTY);
  CheckValue(FUNCTION_NAME, "wgt", "", wgt, ZERO, INFTY);
  CheckValue(FUNCTION_NAME, "t", "", t, ZERO, INFTY);
  CheckValue(FUNCTION_NAME, "spd", "", spd, ZERO, INFTY);

  /***************************************************************************/

  /***** Photon transport ****************************************************/

  /* Get particle type */

  if ((long)RDB[part + PARTICLE_TYPE] == PARTICLE_TYPE_GAMMA)
    {
      /* Flux */

      ptr = (long)RDB[RES_TOT_PHOTON_FLUX];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(flx, wgt, ptr, id, 0);

      /* Avoid compiler warning */

      heat = 0.0;
      tot = 0.0;

      /* Check material pointer */

      if (mat > VALID_PTR)
        {
          /* Total reaction rate */

          if ((rea = (long)RDB[mat + MATERIAL_PTR_TOTPHOTXS]) > VALID_PTR)
            if ((tot = PhotonMacroXS(rea, E, id)*flx*g) >= 0.0)
              {
                ptr = (long)RDB[RES_TOT_PHOTON_RR];
                CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                AddBuf1D(tot, wgt, ptr, id, 0);
              }

          /* Photon heating reaction rate */

          if ((rea = (long)RDB[mat + MATERIAL_PTR_HEATPHOTXS]) > VALID_PTR)
            if ((heat = PhotonMacroXS(rea, E, id)*flx*g) >= 0.0)
              {
                ptr = (long)RDB[RES_TOT_PHOTON_HEATRATE];
                CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                AddBuf1D(heat, wgt, ptr, id, 0);
              }

          /* Pointer to normalization */

          if ((norm = (long)RDB[mat + MATERIAL_PTR_NORM]) > VALID_PTR)
            {
              /* Flux */

              ptr = (long)RDB[norm + NORM_PTR_PHOTON_FLUX];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf1D(flx, wgt, ptr, id, 0);

              /* Heat rate */

              ptr = (long)RDB[norm + NORM_PTR_PHOTON_HEATRATE];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf1D(heat, wgt, ptr, id, 0);
            }
        }

      /* Score collision detector */

      ColDet(part, mat, flx, x, y, z, u, v, w, E, t, wgt, g, id);

      /* Score mesh */

      ScoreMesh(part, mat, flx, 0.0, x, y, z, E, t, wgt, g, id);

      /* Score MFP for adaptive RMX */

      ScoreRMXMFP(part, mat, g, flx, tot, wgt, id);

      /* Exit subroutine */

      return;
    }

  /***************************************************************************/

  /***** Get commonly used values *********************************************/

  /* Check material pointer and calculation of implicit reaction rates */

  if ((mat > VALID_PTR) && ((long)RDB[DATA_OPTI_IMPLICIT_RR] == YES))
    {
      /* Total reaction rate */

      rea = (long)RDB[mat + MATERIAL_PTR_TOTXS];
      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
      tot = MacroXS(rea, E, id)*flx*g;

      /* Capture rate (materials consisting of He-4 have no absorption xs) */

      if ((rea = (long)RDB[mat + MATERIAL_PTR_ABSXS]) > VALID_PTR)
        capt = MacroXS(rea, E, id)*flx*g;
      else
        capt = 0.0;

      /* Elastic rate */

      rea = (long)RDB[mat + MATERIAL_PTR_ELAXS];
      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
      ela = MacroXS(rea, E, id)*flx*g;

      /* Inelastic scattering production rate */

      if ((rea = (long)RDB[mat + MATERIAL_PTR_INLPXS]) > VALID_PTR)
        sprod = MacroXS(rea, E, id)*flx*g;
      else
        sprod = 0.0;

      /* Pointer to fission cross section */

      if ((rea = (long)RDB[mat + MATERIAL_PTR_FISSXS]) > VALID_PTR)
        {
          /* Fission rate */

          fiss = MacroXS(rea, E, id)*flx*g;

          /* If fission is not sampled, it is included in capture */

          if ((long)RDB[DATA_NPHYS_SAMPLE_FISS] == NO)
            capt = capt - fiss;

          /* Fission neutron production rate */

          rea = (long)RDB[mat + MATERIAL_PTR_NSF];
          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
          nsf = MacroXS(rea, E, id)*flx*g;

          /* Fission energy production rate */

          rea = (long)RDB[mat + MATERIAL_PTR_FISSE];
          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
          fissE = MacroXS(rea, E, id)*flx*g;
        }
      else
        {
          fiss = 0.0;
          nsf = 0.0;
          fissE = 0.0;
        }
    }
  else
    {
      /* Reset values */

      tot = 0.0;
      capt = 0.0;
      fiss = 0.0;
      fissE = 0.0;
      nsf = 0.0;
      ela = 0.0;
      sprod = 0.0;
    }

  /***************************************************************************/

  /***** Score integral reaction rates ***************************************/

  /* Score flux */

  ptr = (long)RDB[RES_TOT_NEUTRON_FLUX];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddBuf1D(flx, wgt, ptr, id, 0);

  /* Check material pointer */

  if (mat > VALID_PTR)
    {
      /* Delayed neutron precursor production */

      PrecDet(part, mat, flx, g, wgt, x, y, z, t, E, id);

      /* Additional bins for burnable and non-burnable materials */

      if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
        i = 1;
      else
        i = 2;

      /* Flux in burnable and non-burnable materials */

      AddBuf1D(flx, wgt, ptr, id, i);

      /* Check implicit mode */

      if ((long)RDB[DATA_OPTI_IMPLICIT_RR] == YES)
        {
          /* Total reaction rate */

          ptr = (long)RDB[RES_TOT_NEUTRON_RR];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(tot, wgt, ptr, id, 0);

          /* Capture rate */

          ptr = (long)RDB[RES_TOT_CAPTRATE];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          AddBuf1D(capt, wgt, ptr, id, 0);
          AddBuf1D(capt, wgt, ptr, id, i);

          /* Elastic scattering rate */

          ptr = (long)RDB[RES_TOT_ELARATE];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(ela, wgt, ptr, id, 0);

          /* Inelastic scattering production rate */

          ptr = (long)RDB[RES_TOT_INLPRODRATE];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(sprod, wgt, ptr, id, 0);

          /* Check fission cross section */

          if (fiss > 0.0)
            {
              /* Fission rate */

              ptr = (long)RDB[RES_TOT_FISSRATE];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              AddBuf1D(fiss, wgt, ptr, id, 0);
              AddBuf1D(fiss, wgt, ptr, id, i);

              /*  Fission neutron production rate */

              ptr = (long)RDB[RES_TOT_NSF];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              AddBuf1D(nsf, wgt, ptr, id, 0);
              AddBuf1D(nsf, wgt, ptr, id, i);

              /*  Fission energy */

              ptr = (long)RDB[RES_TOT_FISSE];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              AddBuf1D(fissE, wgt, ptr, id, 0);
              AddBuf1D(fissE, wgt, ptr, id, i);

              /* Average energy of neutron causing fission */

              ptr = (long)RDB[RES_IMP_AFGE];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              AddBuf1D(E*fiss, wgt, ptr, id, 0);

              /* Average lethargy of neutron causing fission */

              ptr = (long)RDB[RES_IMP_ALF];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              AddBuf1D(fiss*log(RDB[DATA_NEUTRON_EMAX]/E), wgt, ptr, id, 0);

              /* Delayed neutrons */

              if ((dng = (long)RDB[part + PARTICLE_DN_GROUP]) > 0)
                {
                  /* Score fission rate */

                  ptr = (long)RDB[RES_ADJ_MEULEKAMP_BETA_EFF];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                  AddBuf1D(fiss, wgt, ptr, id, 0);
                  AddBuf1D(fiss, wgt, ptr, id, dng);

                  /* Get decay constant */

                  lambda = RDB[part + PARTICLE_DN_LAMBDA];

                  /* Score lambda */

                  ptr = (long)RDB[RES_ADJ_MEULEKAMP_LAMBDA];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                  AddBuf1D(fiss*lambda, wgt, ptr, id, 0);
                  AddBuf1D(fiss*lambda, wgt, ptr, id, dng);
                }

#ifdef OLD_IFP

              /* Beta-eff using IFP */

              if ((prg = (long)RDB[part + PARTICLE_PTR_FISS_PROG]) > VALID_PTR)
                {
                  /* Loop over progenies */

                  for (n = 0; n < (long)RDB[DATA_IFP_CHAIN_LENGTH]; n++)
                    {
                      /* Check pointer */

                      CheckPointer(FUNCTION_NAME, "(prg)", DATA_ARRAY, prg);

                      /* Get delayed neutron group */

                      if ((dng = (long)RDB[prg + FISS_PROG_DN_GROUP]) > 0)
                        {
                          /* Score beta-eff */

                          ptr = (long)RDB[RES_ADJ_IFP_IMP_BETA_EFF];
                          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                          AddBuf(fiss, wgt, ptr, id, -1, 0, n);
                          AddBuf(fiss, wgt, ptr, id, -1, dng, n);

                          /* Get decay constant */

                          lambda = RDB[prg + FISS_PROG_LAMBDA];

                          /* Score lambda */

                          ptr = (long)RDB[RES_ADJ_IFP_IMP_LAMBDA];
                          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                          AddBuf(fiss*lambda, wgt, ptr, id, -1, 0, n);
                          AddBuf(fiss*lambda, wgt, ptr, id, -1, dng, n);
                        }

                      /* Pointer to next */

                      prg = NextItem(prg);
                    }
                }

#else

              /* Beta-eff using IFP */

              n = 0;

              /* Loop over events */

              prg = (long)RDB[part + PARTICLE_PTR_EVENTS];
              while (prg > VALID_PTR)
                {
                  /* Check if fission */

                  if ((long)RDB[prg + EVENT_TYPE] == EVENT_TYPE_FISS)
                    {
                      /* Get delayed neutron group */

                      if ((dng = (long)RDB[prg + EVENT_DN_GROUP]) > 0)
                        {
                          /* Score beta-eff */

                          ptr = (long)RDB[RES_ADJ_IFP_IMP_BETA_EFF];
                          CheckPointer(FUNCTION_NAME, "(ptr)",DATA_ARRAY, ptr);

                          AddBuf(fiss, wgt, ptr, id, -1, 0, n);
                          AddBuf(fiss, wgt, ptr, id, -1, dng, n);

                          /* Get decay constant */

                          lambda = RDB[prg + EVENT_LAMBDA];

                          /* Score lambda */

                          ptr = (long)RDB[RES_ADJ_IFP_IMP_LAMBDA];
                          CheckPointer(FUNCTION_NAME, "(ptr)",DATA_ARRAY, ptr);

                          AddBuf(fiss*lambda, wgt, ptr, id, -1, 0, n);
                          AddBuf(fiss*lambda, wgt, ptr, id, -1, dng, n);
                        }

                      /* Update count */

                      n++;
                    }

                  /* Check chain length */

                  if (n == (long)RDB[DATA_IFP_CHAIN_LENGTH])
                    break;

                  /* Next event */

                  prg = NextItem(prg);

                }

#endif

            }

          /* Pointer to normalization (tää korvaa myöhemmin ton binningin) */

          if ((norm = (long)RDB[mat + MATERIAL_PTR_NORM]) > VALID_PTR)
            {
              /* Check fission rate */

              if (fiss > 0.0)
                {
                  /* Fission rate */

                  ptr = (long)RDB[norm + NORM_PTR_FISSRATE];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                  AddBuf1D(fiss, wgt, ptr, id, 0);

                  /* Fission neutron production rate */

                  ptr = (long)RDB[norm + NORM_PTR_NSF];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                  AddBuf1D(nsf, wgt, ptr, id, 0);

                  /* Fission energy */

                  ptr = (long)RDB[norm + NORM_PTR_FISSE];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                  AddBuf1D(fissE, wgt, ptr, id, 0);
                }

              /* Flux */

              ptr = (long)RDB[norm + NORM_PTR_NEUTRON_FLUX];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf1D(flx, wgt, ptr, id, 0);
            }
        }
    }

  /* 1/v */

  ptr = (long)RDB[RES_TOT_RECIPVEL];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddBuf1D(flx/spd, wgt, ptr, id, 0);

  /***************************************************************************/

  /***** Score other stuff ***************************************************/

  /* Score MFP for adaptive RMX */

  ScoreRMXMFP(part, mat, g, flx, tot, wgt, id);

  /*
  if (tot < ZERO)
    return;
  */

  /* Reaction rates for uniform fission source method */

  ScoreUFS(flx, mat, wgt, x, y, z, E, g, id);

  /* Poisons */

  ScorePoison(flx, mat, E, wgt, g, id);

  /* Score cross sections for on-the-fly burnup routine */

  ScoreOTFBurn(flx, mat, E, wgt, id);

  /* Iterable nuclides */

  ScoreIterNucs(flx, mat, E, wgt, g, id);

  /* Source rate for source convergence acceleration */

  if (nsf > 0.0)
    ScoreSCAResp(part, wgt*nsf, x, y, z, id);

  /* Score neutron heating for normalization */

  heat = 0.0;

  if (((long)RDB[DATA_EDEP_MODE] > EDEP_MODE_MT458) && (mat > VALID_PTR))
    {
      rea = (long)RDB[mat + MATERIAL_PTR_HEATTXS];
      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
      heat = MacroXS(rea, E, id)*flx*g;

      /* Additional bins for burnable and non-burnable materials */

      if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
        i = 1;
      else
        i = 2;

      /* Correction for energy deposition */

      if (((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT) &&
          (long)RDB[DATA_EDEP_KEFF_CORR])
        heat = heat*RDB[DATA_CYCLE_KEFF];

      ptr = (long)RDB[RES_TOT_FISSE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      AddBuf1D(heat, wgt, ptr, id, 0);
      AddBuf1D(heat, wgt, ptr, id, i);

      /* Normalization to energy deposition in single material */

      if ((norm = (long)RDB[mat + MATERIAL_PTR_NORM]) > VALID_PTR)
        {
          ptr = (long)RDB[norm + NORM_PTR_FISSE];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(heat, wgt, ptr, id, 0);
        }
    }

  /* Score interface flux and power */

  if (((long)RDB[DATA_PTR_IFC0] > VALID_PTR) && (mat > VALID_PTR))
    {
      /* Score flux for multi-physics interface */

      ScoreInterfaceFlux(wgt, flx, E, x, y, z, 0.0, id);

      /* Get fission energy */

      rea = (long)RDB[mat + MATERIAL_PTR_FISSE];

      if ((rea > VALID_PTR) && (fissE == 0.0))
        fissE = MacroXS(rea, E, id)*flx*g;

      /* Score power for multi-physics interface */

      if ((long)RDB[DATA_EDEP_MODE] > EDEP_MODE_MT458)
        ScoreInterfacePower(heat + fissE, wgt, x, y, z, 0.0, part, id);
      else if (rea > VALID_PTR)
        ScoreInterfacePower(fissE, wgt, x, y, z, 0.0, part, id);
    }

  /* Check if active cycle */

  if (RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP])
    return;

  /* Score transmutation cross sections */
  /* Multiply with density factor if material volume has grown due to */
  /* thermal expansion */

  ScoreTransmuXS(g*flx, mat, E, wgt, id);

  /* Score activation detectors */

  ScoreActiDet(g*flx, mat, part, x, y, z, E, wgt, t, id);

  /* Check for corrector step and B1 mode. Removed for dev (AIs) */

#ifndef STAB_BURN

  if (((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP) &&
      ((long)RDB[DATA_B1_BURNUP_CORR] == NO))
    return;

#endif

  /* Score ICM stuff (vaatii implicit rr) */

  ScoreICMCol(flx, capt, fiss, nsf, fissE, sprod, mat, part, x, y, z, E, wgt,
              g, id);

  /* Score group constant data (must be called before the check for */
  /* corrector calculation to get cross sections for B1). */

  ScoreGC(flx, tot, capt, fiss, fissE, nsf, mat, part, E, spd, wgt, u, v, w,
          id);

  /* Check for corrector step. Removed for dev. (AIs) */

#ifndef STAB_BURN

  if ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP)
    return;

#endif

  /* Score collision detector */

  ColDet(part, mat, flx, x, y, z, u, v, w, E, t, wgt, g, id);

  /* Check material pointer */

  if (mat > VALID_PTR)
    {
      /* Score pin-power distribution */

      if (fissE > 0.0)
        ScorePinPower(mat, part, E, fissE, wgt, id);

      /* Score mesh */

      ScoreMesh(part, mat, flx, 0.0, x, y, z, E, t, wgt, g, id);

      /* Get fission energy */

      if (((long)RDB[DATA_CORE_PDE_DEPTH] > 0.0) ||
          ((long)RDB[DATA_PTR_PB0] > VALID_PTR))
        if ((rea = (long)RDB[mat + MATERIAL_PTR_FISSE]) > VALID_PTR)
          if ((fissE = MacroXS(rea, E, id)*flx*g) > 0.0)
            {
              /* Score full-core power distribution */

              ScoreCPD(fissE, wgt, z, id);

              /* Score PB power distribution */

              ScorePB(fissE, wgt, id);
            }

      /* Photon production rate */

      if ((rea = (long)RDB[mat + MATERIAL_PTR_PHOTPXS]) > VALID_PTR)
        if ((val = MacroXS(rea, E, id)*flx*g) > 0.0)
          {
            ptr = (long)RDB[RES_TOT_PHOTON_PRODRATE];
            CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
            AddBuf1D(val, wgt, ptr, id, 0);
          }
    }

  /***************************************************************************/
}

/*****************************************************************************/
