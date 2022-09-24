/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : precdet.c                                      */
/*                                                                           */
/* Created:       2015/05/08 (VVa)                                           */
/* Last modified: 2018/09/28 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Scores implicit precursor production detectors               */
/*                                                                           */
/* Comments:   -There might be many things here that can be done faster      */
/*             -Maybe name this to ScorePrecProd at some point?              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrecDet:"

/*****************************************************************************/

void PrecDet(long part, long mat, double flx, double g, double wgt0, double x0,
             double y0, double z0, double t0, double E0, long id)
{
  long loc0, loc1, rea, nuc, ptr, rls, lst, n, np, i, j, ibin, tbin, gbin;
  long stp, erg, ng, new, idx0, idx1;
  double dnu, rr, f, P, P0, P1, Pemit, Ptal, val, adens, dt, lambda, nprec;
  double macroxs, Nemit, Wemit, Wtresh, prob, temit, Etal, Etresh, Dt;
  double x, y, z, t, u, v, w, E, mu, wgt2, norm, Wtal, impP, tresh;
#ifdef OLD_IFP
  long prg, prev;
#endif

  /* Get pointer to precursor detector */

  if ((loc0 = (long)RDB[DATA_PTR_PREC_DET]) < VALID_PTR)
    return;

  /* Proportion of events to skip */

  impP = 1;

  if (RandF(id) > impP)
    return;
  else
    wgt0 = wgt0/impP;

  /*
  printf("Material %s\n", GetText(mat + MATERIAL_PTR_NAME));
  */

  /* Get normalization */

  if (RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
    {
      /* In criticality source mode values will be */
      /* normalized only in collectprecdet.c */

      norm = 1.0;
    }
  else
    {
      /* In dynamic mode the values are */
      /* multiplied by normalization already in precdet.c */

      norm = RDB[DATA_NORM_COEF_N];
    }

  /* Get treshold for storing additional precursors in transients */

  tresh = RDB[DATA_PREC_STORE_TRESH];

  /* Get pointer to the precursor population stat object */

  stp = (long)RDB[loc0 + PRECDET_PTR_STAT];
  CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);

  /**********************************************/
  /* Get spatial bin and time bin for the tally */

  /* Get pointer to spatial mesh */

  ptr = (long)RDB[loc0 + PRECDET_PTR_MESH];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get mesh bin */

  ibin = MeshIndex(ptr, x0, y0, z0, -1.0);

  /* Get time bin */
  /* +1 since the scores are added to EOI values */

  tbin = (long)RDB[DATA_DYN_TB] + 1;

  /* Check if material is fissile */

  if ((rls = (long)RDB[mat + MATERIAL_PTR_FISS_REA_LIST]) < VALID_PTR)
    return;

  /* Calculate interval length */

  Dt = RDB[DATA_TIME_CUT_TMAX] - RDB[DATA_TIME_CUT_TMIN];

  /*
  printf("time bin %ld, mesh bin %ld\n", tbin, ibin);
  */
  rls = (long)RDB[rls + RLS_PTR_REA0];
  CheckPointer(FUNCTION_NAME, "(rls)", DATA_ARRAY, rls);

  /* Loop over fission list and score contributions from different reactions */
  /* Loop over reactions */

  while (rls > VALID_PTR)
    {

      /* Get pointer to reaction */

      rea = (long)RDB[rls + RLS_DATA_PTR_REA];
      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

      /* Get pointer to nuclide */

      nuc = (long)RDB[rls + RLS_DATA_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Get delayed neutron nubar for reaction */
      /* so we don't have to think about beta   */

      if ((ptr = (long)RDB[rea + REACTION_PTR_DNUBAR]) > VALID_PTR)
        {
          dnu = Nubar(ptr, E0, id);
        }
      else
        {
          /* Next reaction */

          rls = NextItem(rls);

          /* Cycle loop */

          continue;
        }

      /* Scale nubar with external source k-eff */

      dnu = dnu*RDB[DATA_EXT_SRC_NUBAR_F];

      /* Get incident delayed neutron group and decay constant */

      ng = (long)RDB[part + PARTICLE_DN_GROUP];
      lambda = RDB[part + PARTICLE_DN_LAMBDA];

      /* Get atomic density of nuclide in material    */
      /* TODO: There might be a faster way to do this */

      ptr = (long)RDB[mat + MATERIAL_PTR_COMP];

      /* Loop over composition to find correct nuclide */

      while (ptr > VALID_PTR)
        {
          /* Check if match */

          if (nuc == (long)RDB[ptr + COMPOSITION_PTR_NUCLIDE])
            break;

          /* No match, next nuclide */

          ptr = NextItem(ptr);
        }

      /* We should have found the correct nuclide from composition */

      if (ptr < VALID_PTR)
        Die(FUNCTION_NAME, "Could not find nuclide from mat. comp.");

      /* Get atomic density from composition */

      adens = RDB[ptr + COMPOSITION_ADENS];

      /* Calculate macroscopic cross section */

      macroxs = MicroXS(rea, E0, id)*adens;

      /* Get total number of reactions caused by this interaction */

      rr = macroxs*flx*g;

      /* Number of precursor produced in total is dnu*rr */

      nprec = dnu*rr;

      /************************************************/
      /* Move on to tally contributions to each group */
      /************************************************/

      /* Loop over groups to tally contributions */

      /* Pointer to precursor group data */

      lst = (long)RDB[rea + REACTION_PTR_PREC_LIST];
      CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);

      /* Loop over groups */

      n = 0;

      while ((loc1 = ListPtr(lst, n++)) > VALID_PTR)
        {
          /* Get decay constant */

          lambda = RDB[loc1 + PREC_LAMBDA];

          /*****************/
          /* Get group bin */
          /*****************/

          /* Get number of stored precursor groups */

          ng = (long)RDB[loc0 + PRECDET_NG];

          /* Get pointer to stored lambdas */

          ptr = (long)RDB[loc0 + PRECDET_PTR_LAM_ARRAY];

          /* Find correct lambda */

          for (gbin = 0; gbin < ng; gbin++)
            if (RDB[ptr + gbin] == lambda)
              break;

          /* Check that group was found */

          if (gbin == ng)
            Die(FUNCTION_NAME, "Could not find precursor group to score");

          /*********************************************/
          /* Calculate proportion to add to this group */
          /*********************************************/

          /* Pointer to energy grid */

          erg = (long)RDB[loc1 + PREC_PTR_EGRID];
          CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

          /* Get interpolation factor for energy grid */

          if ((f = GridFactor(erg, E0, id)) < 0.0)
            break;

          /* Get number of points in energy grid */

          np = (long)RDB[erg + ENERGY_GRID_NE];

          /* Separate integer and decimal parts of interpolation factor */

          i = (long)f;
          f = f - (double)i;

          /* Check boundaries */

          if ((i < 0) || (i > np - 2))
            break;
          else
            {
              /* Pointer to probabilities */

              ptr = (long)RDB[loc1 + PREC_PTR_PTS];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Get tabulated points */

              P0 = RDB[ptr + i];
              P1 = RDB[ptr + i + 1];

              /* Interpolate between energy points */

              P = f*(P1 - P0) + P0;
            }

          /* Contribution to group n is P*nprec */

          /**************************************************************/
          /* Emit the part of the weight to be emitted on this interval */
          /**************************************************************/

          /* Calculate time to end of interval */

          dt = RDB[DATA_TIME_CUT_TMAX] - t0;

          /* Avoid compiler warning */

          Ptal = -1;

          /* Calculate proportion of value that is left at the end of */
          /* the time interval*/

          if (RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
            Ptal = 1.0;
          else if ((RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DELDYN) ||
                   (RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN))
            Ptal = exp(-lambda*dt);
          else
            Die(FUNCTION_NAME, "Unknown simulation mode");

          /* Store rest of the value to the tally */

          val = Ptal*P*nprec;

          /* Tally value multiplied by normalization */
          /* to get physical neutrons */

          AddBuf(val*norm, wgt0, stp, id, -1, tbin, gbin, ibin);

          /* Create new precursor to list if we pass russian roulette */

          if ((RDB[DATA_SIMULATION_MODE] != SIMULATION_MODE_CRIT) &&
              (RDB[DATA_PRECURSOR_TRANSPORT_MODE] == PREC_MODE_POINT))
            {
              /* Calculate weight to be tallied */

              Wtal = val*wgt0;

              /* Calculate emission of the to-be-tallied population     */
              /* This is how many neutrons the to-be-tallied population */
              /* would emit on an interval that is as long as this one  */
              /* The previously banked precursors have been normalized  */
              /* to emit an equal number of delayed neutrons during the */
              /* interval */

              Etal = Wtal*norm*(1-exp(-lambda*Dt));

              /* Set weight that we want to emit */
              /* Currently set to the incoming neutron weight */
              /* TODO: Test other thresholds */
              /* such a weight that wgt0 will be emitted before */
              /* end of simulation */

              /* Set activity treshold */
              /* Calculated first in normalizeprecdet.c */
              /* And later in precursorpopcontrol.c */

              Etresh = RDB[loc0 + PRECDET_AVE_EMIT];

              /* Calculate the probability to store */

              if (tresh > 0)
                prob = Etal/(tresh*Etresh);
              else if (tresh < 0)
                prob = Wtal/(-tresh*wgt0);
              else
                prob = 0;

              /* Precursors that should be stored with larger than unity probability    */
              /* are stored with 100% probability at their larger than AVE_EMIT weight  */
              /* weights will be normalized before next interval in precursorpopcontrol */

              if (prob > 1.0)
                prob = 1.0;

              /* Increase the weight to threshold if not already over */

              if (RandF(id) < prob)
                {
                  /* If we survive russian roulette create new precursor to bank */

                  new = FromStack(PARTICLE_TYPE_PRECURSOR, id);

                  /* Put stuff */

#ifdef DNPRINT
                  /*
                  fprintf(outp, "Emission from this would have been %E\n", Etal);
                  fprintf(outp, "Emission from this will be %E\n", Etal/prob);
                  */
#endif

                  WDB[new + PARTICLE_X] = x0;
                  WDB[new + PARTICLE_Y] = y0;
                  WDB[new + PARTICLE_Z] = z0;

                  WDB[new + PARTICLE_PTR_MAT] = (double)mat;

                  /* Global DN group (not nuclide wise) starts from 0 */

                  WDB[new + PARTICLE_DN_GROUP] = (double)gbin;
                  WDB[new + PARTICLE_DN_LAMBDA] = lambda;

                  /* Precursor _WGT is physical neutrons */

                  WDB[new + PARTICLE_WGT] = Wtal/prob*norm;
                  WDB[new + PARTICLE_T0] = t0;
                  WDB[new + PARTICLE_T] = RDB[DATA_TIME_CUT_TMAX];

                  /* Update generation index */

                  WDB[new + PARTICLE_GEN_IDX] = RDB[part + PARTICLE_GEN_IDX] + 1.0;

                  /* Put to bank */
                  /* These will get mixed up with prompt neutons */
                  /* They will be separated in: getbankedprecursors.c */

                  ToBank(new, id);
                }
            }

          /* In criticality source mode we are not emitting anything */

          if (RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
            {
              /* Write point and precursor information to file */
              /* TODO: Should we divide by lambda here?                */
              /* Now we divide by lambda just when storing the weight  */
              /* meaning that the weights stored are stable weights    */
              /* val should actually correspond to activity (in steady */
              /* state) as the incoming neutrons are 1/s in CRIT-mode  */
              /* Seems to give physical results this way^              */

              ptr = (long)RDB[loc0 + PRECDET_PTR_PREC_DET];

              /* Get pointer to detector */

              /* Set weight that we want to tally             */
              /* Currently set to the incoming neutron weight */
              /* without biasing correction                   */

              Wtresh = wgt0*impP;

              /* Precursor production rate is val*wgt0           */
              /* In steady state it's also the emission rate     */
              /* for delayed neutrons                            */

              /* This gives good group distribution for short    */
              /* transients / time intervals, but overrepresents */
              /* short lived precursors in the source for long   */
              /* transients */

              Wemit = val*wgt0;

              if (Wemit < Wtresh)
                {
                  /* Calculate the probability to store this */

                  prob = Wemit/Wtresh;

                  /* Increase the weight to emit */

                  Wemit = Wemit/prob;

                  /* Sample random number from unit interval */
                  /* Store with weight Wtresh if random number < prob */
                  /* Otherwise do not store */

                  if (RandF(id) < prob)
                    WriteSourceFile(ptr, x0, y0, z0, 1.0, 0.0, 0.0, E0,
                                    Wemit/lambda, (double)gbin, -1.0, id);
                }
              else
                Die(FUNCTION_NAME, "Too large precursor production");

              /* Move to tally production for next group */

              continue;
            }

          /* Get fission matrix indexes */

          /* Score emitted neutron */
          /* Scoring of fissions probably will be messed up */
          /* Due to implicit treatment of delayed neutrons  */
          /* Should we score here? */

          idx0 = (long)RDB[part + PARTICLE_FMTX_IDX];
          idx1 = FissMtxIndex(mat, id);

          ScoreFission(mat, rea, dnu, dnu, lambda, ng, E0, 0.0, wgt0, 0.0,
                       idx0, idx1, id);

          /****************************************/
          /* Move on to emit the remaining weight */

          /* Calculate the proportion to emit now */

          Pemit = 1.0 - Ptal;

          /* Calculate the number to emit now */

          Nemit = Pemit*P*nprec;

          /* Calculate the weight to emit now */

          Wemit = Nemit*wgt0;

#ifdef mmmaaa
          /* Get generation index */

          i = (long)RDB[part + PARTICLE_GEN_IDX] + 1;

          /* Score new source weight */

          ptr = (long)RDB[RES_NEW_SRC_WGT];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(P*nprec*wgt0, 1.0, ptr, id, 0);

          if (i < MAX_EXT_K_GEN + 1)
            AddBuf1D(P*nprec*wgt0, 1.0, ptr, id, i);

#endif
          /*
#ifdef ajflkb
          */
          /* Get generation index */

          i = (long)RDB[part + PARTICLE_GEN_IDX] + 1;

          /* Score new source weight */

          ptr = (long)RDB[RES_NEW_SRC_WGT];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(Wemit, 1.0, ptr, id, 0);

          if (i < MAX_EXT_K_GEN + 1)
            AddBuf1D(Wemit, 1.0, ptr, id, i);
          /*
#endif
          */
          /* Use Russian Roulette */

          /* Set weight that we want to emit */
          /* Currently set to the incoming neutron weight */

          Wtresh = wgt0*impP;

          if (Wemit < Wtresh)
            {
              /* Calculate the probability to emit */

              prob = Wemit/Wtresh;

              /* Sample random number from unit interval */
              /* Emit with weight Wtresh if random number < prob */
              /* Otherwise do not emit */

              if (RandF(id) < prob)
                {

                  /* Handle emission */

                  /* Sample emission time */
                  /* (time from this moment to emission) */

                  temit = -1/lambda*log(1 - (1 - exp(-lambda*dt))*RandF(id));

                  /* Calculate absolute time of emission */

                  t = t0 + temit;

                  /* Emission coordinates are these interaction coordinates */

                  x = x0;
                  y = y0;
                  z = z0;

                  /* Sample emission energy */
                  /* Pointer to energy grid */

                  erg = (long)RDB[loc1 + PREC_PTR_ERG];
                  CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

                  /* Reset mu */

                  mu = 0.1;

                  /* Sample energy */
                  /* Using incoming neutron energy */
                  /* Sampling is probably independent from incoming neutron */
                  /* energy, but it has to be checked at some point */

                  SampleENDFLaw(rea, erg, E0, &E, &mu, id);

                  /* Check that direction cosine was not changed */
                  /* If it was changed, the delayed neutron emission reaction */
                  /* actually has some kind of directional data, which would  */
                  /* be weird */

                  if (mu != 0.1)
                    Die(FUNCTION_NAME, "mu was changed by SampleENDFLaw");

                  /* Sample emission direction cosines isotropically */
                  /* The direction for the delayed neutron should not */
                  /* depend on the direction of the neutron that caused */
                  /* the fission (right?) */

                  IsotropicDirection(&u, &v, &w, id);

                  /* Set outgoing neutron weight */

                  wgt2 = Wtresh;

                  /* Create the neutron to be emitted and set the parameters */

                  /* Duplicate incident neutron */

                  new = DuplicateParticle(part, id);

                  /* Put parameters */

                  WDB[new + PARTICLE_X] = x;
                  WDB[new + PARTICLE_Y] = y;
                  WDB[new + PARTICLE_Z] = z;

                  WDB[new + PARTICLE_U] = u;
                  WDB[new + PARTICLE_V] = v;
                  WDB[new + PARTICLE_W] = w;

                  WDB[new + PARTICLE_E] = E;
                  WDB[new + PARTICLE_WGT] = wgt2;
                  WDB[new + PARTICLE_PTR_MAT] = (double)mat;
                  WDB[new + PARTICLE_DN_GROUP] = RDB[loc1 + PREC_IDX];
                  WDB[new + PARTICLE_DN_LAMBDA] = lambda;

                  /* Store absolute emission time */
                  /* PARTICLE_T0 is the birth time of the neutron*/
                  /* PARTICLE_T  is the current time of the neutron */

                  WDB[new + PARTICLE_T0] = t;
                  WDB[new + PARTICLE_T] = t;

                  if(t > RDB[DATA_TIME_CUT_TMAX])
                    Die(FUNCTION_NAME, "Delayed neutron emitted on next interval?!");

                  /* Delayed neutron emission time */

                  WDB[new + PARTICLE_TD] = temit;

                  /* Put fission matrix index */
                  /*
                  WDB[new + PARTICLE_FMTX_IDX] = (double)idx1;
                  */
                  /* Reset thermalization time */

                  WDB[new + PARTICLE_TT] = 0.0;

                  /* Update generation index */

                  WDB[new + PARTICLE_GEN_IDX] = RDB[new + PARTICLE_GEN_IDX] + 1.0;

                  /* Score time-dependent source rate detectors */

                  ScoreTimeSource(part, mat, MT_SECONDARY_DN_SOURCE,
                                  x, y, z, E, t, 1, wgt2, id);

                  /* Update weight for ICM */
                  /* I don't know if this should be done */

                  ptr = (long)RDB[DATA_PTR_CYCLE_EIG_KEFF];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                  j = (long)RDB[new + PARTICLE_EIG_IDX];
                  WDB[new + PARTICLE_ICM_WGT] =
                    RDB[new + PARTICLE_ICM_WGT]*RDB[ptr + j];

                  /* Put lifetimes and delayed neutron groups to progenies */
                  /* I don't know if this should be done either */

#ifdef OLD_IFP
                  if ((prg = (long)RDB[new + PARTICLE_PTR_FISS_PROG]) > VALID_PTR)
                    {
                      /* Loop from last to first */

                      prg = LastItem(prg);
                      while (prg > VALID_PTR)
                        {
                          /* Get pointer to previous */

                          if ((prev = PrevItem(prg)) > VALID_PTR)
                            {
                              /* Copy data from previous */

                              WDB[prg + FISS_PROG_DN_GROUP] =
                                WDB[prev + FISS_PROG_DN_GROUP];
                              WDB[prg + FISS_PROG_LIFETIME] =
                                WDB[prev + FISS_PROG_LIFETIME];
                              WDB[prg + FISS_PROG_LAMBDA] =
                                WDB[prev + FISS_PROG_LAMBDA];
                            }
                          else
                            {
                              /* Put new */

                              WDB[prg + FISS_PROG_DN_GROUP] =
                                WDB[part + PARTICLE_DN_GROUP];
                              WDB[prg + FISS_PROG_LIFETIME] =
                                t0 - RDB[part + PARTICLE_T0];
                              WDB[prg + FISS_PROG_LAMBDA] =
                                WDB[part + PARTICLE_DN_LAMBDA];
                            }

                          /* Previous progeny */

                          prg = prev;
                        }
                    }
#endif

                  /* Bank or que neutron */

                  if (RDB[new + PARTICLE_T] < RDB[DATA_TIME_CUT_TMIN])
                    Die(FUNCTION_NAME, "Error in time");
                  else if (RDB[new + PARTICLE_T] >= RDB[DATA_TIME_CUT_TMAX])
                    Die(FUNCTION_NAME, "Delayed neutron emitted on next interval?!");
                  else
                    ToQue(new, id);

                  /* Score collision detector */
                  /*
                  ColDet(new, mat, 1.0, 0.0, x, y, z, u, v, w, E, t, wgt2, -2.0, id);
                  */
                  /* Score mesh plotter */
                  /*
                  ScoreMesh(new, mat, 0.0, -2.0, x, y, z, E, t, wgt2, 1.0, id);
                  */

                }

            }
          else
            Die(FUNCTION_NAME, "Too large precursor emission");

        }

      rls = NextItem(rls);
    }

}
