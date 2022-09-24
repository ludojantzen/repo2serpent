/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : fission.c                                      */
/*                                                                           */
/* Created:       2011/03/07 (JLe)                                           */
/* Last modified: 2018/09/28 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Handles fission                                              */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "Fission:"

/*****************************************************************************/

void Fission(long mat, long rea, long part, double *E0, double t0, double x,
             double y, double z, double *u, double *v, double *w, double wgt1,
             double *wgt2, long id)
{
  long ptr, nmax, n, i, ieig, new, mode, ng, idx0, idx1, gcu, ncol;
  long imat, izai, ienePart, ieneNewRej, gidx;
  long ieneNew, ptr1, mt, sens;
  double tnu, dnu, beta, E, f, mu, u0, v0, w0, lambda, td;
  double Erej, Eacc, muReject, t;

#ifdef OLD_IFP

  long prg, prev;

#endif

  /* Reset some indices */

  ieneNewRej = 0;
  ieneNew = 0;
  Erej = 0.0;
  muReject = 0.0;
  Eacc = 0.0;

  /* Check pointers */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Check coordinates, energy and weight */

  CheckValue(FUNCTION_NAME, "x", "", x, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "y", "", y, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "z", "", z, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "E0", "", *E0, ZERO, INFTY);
  CheckValue(FUNCTION_NAME, "wgt1", "", wgt1, ZERO, INFTY);

  /* Avoid compiler warning */

  tnu = -1.0;

  /* Get total nubar */

  if ((ptr = (long)RDB[rea + REACTION_PTR_TNUBAR]) > VALID_PTR)
    tnu = Nubar(ptr, *E0, id);
  else
    Die(FUNCTION_NAME, "No prompt nubar data");

  /* Get delayed nubar */

  if ((ptr = (long)RDB[rea + REACTION_PTR_DNUBAR]) > VALID_PTR)
    dnu = Nubar(ptr, *E0, id);
  else
    dnu = 0.0;

  /* Get eigenvalue index */

  ieig = (long)RDB[part + PARTICLE_EIG_IDX];
  CheckValue(FUNCTION_NAME, "(ieig)", "", ieig, 0, 9999);

  /* Check delayed neutron flag */

  if ((long)RDB[DATA_USE_DELNU] != YES)
    {
      /* Subtract delayed nubar from total */

      tnu = tnu - dnu;

      /* Set delayed nubar to zero */

      dnu = 0.0;
    }

  /* Calculate delayed neutron fraction */

  beta = dnu/tnu;
  CheckValue(FUNCTION_NAME, "beta", "", beta, 0.0, 3E-2);

  /* Original direction cosines */

  u0 = *u;
  v0 = *v;
  w0 = *w;

  /* Get simulation mode */

  mode = (long)RDB[DATA_SIMULATION_MODE];

  /* Get incident delayed neutron group and decay constant */

  ng = (long)RDB[part + PARTICLE_DN_GROUP];
  lambda = RDB[part + PARTICLE_DN_LAMBDA];

  /* Get fission matrix indexes */

  idx0 = (long)RDB[part + PARTICLE_FMTX_IDX];
  idx1 = FissMtxIndex(mat, id);

  /* Score incident neutron */

  ScoreFission(mat, rea, tnu, dnu, lambda, ng, *E0, 0.0, wgt1, 0.0,
               idx0, idx1, id);

  /* Sample number of emitted neutrons and set weight */

  if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
    {
      /* Criticality source simulation, get UFS factor */

      f = UFSFactor(x, y, z);
      CheckValue(FUNCTION_NAME, "f", "", f, ZERO, INFTY);

      /* Check if weight windows are used */

      if ((long)RDB[DATA_USE_WEIGHT_WINDOWS] == YES)
       {
         /* Get factor (assume no energy dependence) */

         f = WWImportance(PARTICLE_TYPE_NEUTRON, x, y, z, *u, *v, *w, 1.0,
                          WWMESH_SRC);

         /* Check zero */

         if (f == 0.0)
           f = 1.0;
       }

      /* Sample number of fission neutrons and set weight */

      nmax = SampleNu(wgt1*tnu*f, ieig, id);
      *wgt2 = 1.0/f;
    }
  else if ((long)RDB[DATA_OPT_IMPL_FISS] == YES)
    {
      /* External source simulation with implicit fission, */
      /* get given nubar */

      f = RDB[DATA_OPT_IMPL_FISS_NUBAR];
      CheckValue(FUNCTION_NAME, "f", "", f, ZERO, INFTY);

      /* Sample number of fission neutrons and set weight */

      nmax = SampleNu(f, ieig, id);
      *wgt2 = tnu*wgt1/f;
    }
  else
    {
      /* External source simulation with analog fission,  */
      /* sample number of fission neutrons and set weight */

      nmax = SampleNu(tnu, ieig, id);
      *wgt2 = wgt1;
    }

  /* Check for group constant generation */

  if ((long)RDB[DATA_OPTI_GC_CALC] == YES)
    {
      /* Get collision number */

      ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
      CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
      ncol = (long)GetPrivateData(ptr, id);

      /* Get universe for group constant generation */

      gcu = (long)TestValuePair(DATA_GCU_PTR_UNI, (double)ncol, id);
    }
  else
    gcu = -1;

  /* Check if events are recorded */

  if (((long)RDB[DATA_EVENT_RECORD_FLAGS] & RECORD_EVENT_IFP) ||
      ((long)RDB[DATA_EVENT_RECORD_FLAGS] & RECORD_EVENT_SENS))
    {
      /****************************************************************/
      /* Record fission for old (copied for new in duplicateparticle) */
      /****************************************************************/

      /* New event from bank */

      ptr = EventFromBank(part, id);

      /* Put type */

      WDB[ptr + EVENT_TYPE] = (double)EVENT_TYPE_FISS;

      /* Put time */

      WDB[ptr + EVENT_T] = t0;

      /* Put lifetime and delayed neutron group */

      WDB[ptr + EVENT_DN_GROUP] = RDB[part + PARTICLE_DN_GROUP];
      WDB[ptr + EVENT_LIFETIME] = t0 - RDB[part + PARTICLE_T0];
      WDB[ptr + EVENT_LAMBDA] = RDB[part + PARTICLE_DN_LAMBDA];

      /* Put number of accepted collisions during neutron life */
      /* (can this be inferred from EVENTS?) */

      WDB[ptr + EVENT_ACCEPTED_COLLS] = RDB[part + PARTICLE_NA_COLL];

      /* Put material */

      WDB[ptr + EVENT_PTR_MAT] = (double)mat;

      if ((sens = (long)RDB[DATA_PTR_SENS0]) > VALID_PTR)
        {
          /* Get indices */

          FindSensIndices(mat, rea, *E0, &imat, &izai, &ienePart);

          if ((izai >= 0) && (imat >= 0))
            {
              if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_XS)
                {

                  /* Get pointer to index list */

                  ptr1 = (long)RDB[sens + SENS_PTR_PERT_INDICES];
                  CheckPointer(FUNCTION_NAME, "(ptr1)", DATA_ARRAY, ptr1);

                  /* Calculate label */

                  gidx = CompressSensLabel(imat, izai, (long)RDB[ptr1 + FISS_IDX], ienePart, COLL_HIT);

                  /* Put sensitivity label */

                  WDB[ptr + EVENT_LABEL] = (double)gidx;

                  /* Put sensitivity value */

                  WDB[ptr + EVENT_VAL] = 1.0;

                  /* Create a new event (if scoring direct) */

                  if ((long)RDB[DATA_SENS_SCORE_TYPE] == SENS_SCORE_TYPE_DIRECT)
                    StoreSensEvent(part, gidx, *E0, 1.0, id);
                }
              else if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_XSMT)
                {

                  mt = (long)RDB[rea + REACTION_MT];

                  /* Get pointer to index list */

                  ptr1 = (long)RDB[sens + SENS_PTR_MT_INDICES];
                  CheckPointer(FUNCTION_NAME, "(ptr1)", DATA_ARRAY, ptr1);

                  /* Calculate label */

                  gidx = CompressSensLabel(imat, izai, (long)RDB[ptr1 + mt], ienePart, COLL_HIT);

                  /* Put sensitivity label */

                  WDB[ptr + EVENT_LABEL] = (double)gidx;

                  /* Put sensitivity value */

                  WDB[ptr + EVENT_VAL] = 1.0;

                  /* Create a new event (if scoring direct) */

                  if ((long)RDB[DATA_SENS_SCORE_TYPE] == SENS_SCORE_TYPE_DIRECT)
                    StoreSensEvent(part, gidx, *E0, 1.0, id);
                }
            }
        }
      /* Put energy */

      WDB[ptr + EVENT_E] = *E0;
    }

  /* Loop over source neutrons */

  for (n = 0; n < nmax; n++)
    {
      /* Reset mu */

      mu = 0.1;

      /* Reset precursor group and decay constant */

      ng = 0;
      lambda = -1.0;

      /* Sample between prompt and delayed neutron */

      if (RandF(id) > beta)
        {
          /* Prompt neutron, sample energy */

          SampleENDFLaw(rea, -1, *E0, &E, &mu, id);

          /*********** Fission spectrum sensitivities *****************/

          if ((sens = (long)RDB[DATA_PTR_SENS0]) > VALID_PTR)
            {
              /* chi sensitivity */

              if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_CHI)
                {
                  /* Find energy grid index for the accepted fission */
                  /* neutron energy */

                  ieneNew = FindSensEIndex(E);
                  Eacc = E;

                  /* internally constraining chi sensitivity */

                  if (ieneNew > 0)
                    {
                      /* Reset mu */

                      muReject = 0.1;

                      /* sample E (and mu if part of sampling law) */

                      SampleENDFLaw(rea, -1, *E0, &Erej, &muReject, id);

                      /* Find energy index for the rejected fission */
                      /* neutron energy */

                      ieneNewRej = FindSensEIndex(Erej);
                    }
                }
            }

          /***************************/

        }
      else
        {
          /* Delayed neutron, sample precursor group */

          if ((ptr = SamplePrecursorGroup(rea, *E0, id)) > VALID_PTR)
            {
              /* Get precursor group */

              ng = (long)RDB[ptr + PREC_IDX];
              CheckValue(FUNCTION_NAME, "ng", "", ng, 1, 8);

              /* Get decay constant */

              lambda = RDB[ptr + PREC_LAMBDA];
              CheckValue(FUNCTION_NAME, "lambda", "", lambda, ZERO, INFTY);

              /* Get pointer to energy distribution */

              ptr = (long)RDB[ptr + PREC_PTR_ERG];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Sample energy */

              SampleENDFLaw(rea, ptr, *E0, &E, &mu, id);

              /*********** Fission spectrum sensitivity *****************/

              if ((sens = (long)RDB[DATA_PTR_SENS0]) > VALID_PTR)
                {
                  /* Get number of energy groups */

                  if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_CHI)
                    {
                      /* Find energy grid index for the accepted fission */
                      /* neutron energy */

                      ieneNew = FindSensEIndex(E);
                      Eacc = E;

                      /* internally constraining chi sensitivity */

                      if (ieneNew > 0)
                        {
                          /* Reset mu */

                          muReject = 0.1;

                          /* sample E (and mu if part of sampling law) */

                          SampleENDFLaw(rea, ptr, *E0, &Erej, &muReject, id);

                          /* Find energy index for the rejected fission */
                          /* neutron energy */

                          ieneNewRej = FindSensEIndex(Erej);
                        }
                    }
                }

              /***************************/

            }
          else
            {
              /* Sampling failed, sample prompt neutron energy */

              SampleENDFLaw(rea, -1, *E0, &E, &mu, id);

              /*********** Fission spectrum sensitivity *****************/

              if ((sens = (long)RDB[DATA_PTR_SENS0]) > VALID_PTR)
                {
                  /* Get number of energy groups */

                  if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_CHI)
                    {
                      /* Find energy grid index for the accepted fission */
                      /* neutron energy */

                      ieneNew = FindSensEIndex(E);
                      Eacc = E;

                      /* internally constraining chi sensitivity */

                      if (ieneNew > 0)
                        {
                          /* Reset mu */

                          muReject = 0.1;

                          /* sample E (and mu if part of sampling law) */

                          SampleENDFLaw(rea, -1, *E0, &Erej, &muReject, id);

                          /* Find energy index for the rejected fission */
                          /* neutron energy */

                          ieneNewRej = FindSensEIndex(Erej);
                        }
                    }
                }
              /***************************/

            }
        }

      /* Adjust minimum and maximum energy */

      if (E < 1.0000001*RDB[DATA_NEUTRON_EMIN])
        E = 1.000001*RDB[DATA_NEUTRON_EMIN];
      else if (E > 0.999999*RDB[DATA_NEUTRON_EMAX])
        E = 0.999999*RDB[DATA_NEUTRON_EMAX];

      /* Check if mu is sampled from distribution or sample isotropic */

      if (mu == 0.1)
        IsotropicDirection(u, v, w, id);
      else
        {
          /* Get incident direction cosines */

          *u = u0;
          *v = v0;
          *w = w0;

          /* Sanity check for mu and direction vectors (for NAN's etc.) */

          CheckValue(FUNCTION_NAME, "mu", "", mu, -1.01, 1.01);
          CheckValue(FUNCTION_NAME, "u", "", *u, -1.01, 1.01);
          CheckValue(FUNCTION_NAME, "v", "", *v, -1.01, 1.01);
          CheckValue(FUNCTION_NAME, "w", "", *w, -1.01, 1.01);

          /* Rotate */

          AziRot(mu, u, v, w, id);
        }

      /* Sample delayed neutron emission time */

      if (lambda > 0.0)
        td = -log(RandF(id))/lambda;
      else
        td = 0.0;

      /* Score emitted neutron */

      ScoreFission(mat, rea, 0.0, 0.0, lambda, ng, 0.0, E, wgt1, *wgt2,
                   idx0, idx1, id);

      /* Duplicate incident neutron */

      new = DuplicateParticle(part, id);

      /* Put variables */

      WDB[new + PARTICLE_X] = x;
      WDB[new + PARTICLE_Y] = y;
      WDB[new + PARTICLE_Z] = z;

      WDB[new + PARTICLE_U] = *u;
      WDB[new + PARTICLE_V] = *v;
      WDB[new + PARTICLE_W] = *w;

      WDB[new + PARTICLE_E] = E;
      WDB[new + PARTICLE_WGT] = *wgt2;
      WDB[new + PARTICLE_PTR_MAT] = (double)mat;
      WDB[new + PARTICLE_DN_GROUP] = (double)ng;
      WDB[new + PARTICLE_DN_LAMBDA] = lambda;
      WDB[new + PARTICLE_PTR_GCU] = (double)gcu;
      WDB[new + PARTICLE_NA_COLL] = 0.0;

      /* Put next generation event block if needed */

      if (((long)RDB[DATA_EVENT_RECORD_FLAGS] & RECORD_EVENT_SENS)
          && ((long)RDB[DATA_SENS_SCORE_TYPE] == SENS_SCORE_TYPE_DIRECT))
          (void)EBlockFromBank(new, id);

      /* Update weight for ICM */

      ptr = (long)RDB[DATA_PTR_CYCLE_EIG_KEFF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      WDB[new + PARTICLE_ICM_WGT] =
        RDB[new + PARTICLE_ICM_WGT]*RDB[ptr + ieig];

      /**********************************************/
      /* Add events for nubar and chi perturbations */
      /**********************************************/

      if (((sens = (long)RDB[DATA_PTR_SENS0]) > VALID_PTR) &&
          (izai >= 0) && (imat >= 0))
        {
          /* Get pointer to index list */

          ptr1 = (long)RDB[sens + SENS_PTR_PERT_INDICES];

          /* Score sensitivity to nubar */

          if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_NUBAR)
            {
              if (ng > 0)
                {
                  /* Calculate label */

                  gidx = CompressSensLabel(imat, izai,
                                           (long)RDB[ptr1 + NUBAR_DEL_IDX], ienePart, COLL_HIT);

                  /* Create event */

                  StoreSensEvent(new, gidx, *E0, 1.0, id);

                  /* Score custom perturbations for delayed nubar*/

                  ScorePerturbations(new, mat, rea, (long)RDB[ptr1 + NUBAR_DEL_IDX], *E0, (double)COLL_HIT, 1.0, id);
                }
              else
                {
                  /* Calculate label */

                  gidx = CompressSensLabel(imat, izai, (long)RDB[ptr1 + NUBAR_PRO_IDX], ienePart, COLL_HIT);

                  /* Create event */

                  StoreSensEvent(new, gidx, *E0, 1.0, id);

                  /* Score custom perturbations for prompt nubar*/

                  ScorePerturbations(new, mat, rea, (long)RDB[ptr1 + NUBAR_PRO_IDX], *E0, (double)COLL_HIT, 1.0, id);
                }

              /* Score custom perturbations for total nubar*/

              ScorePerturbations(new, mat, rea, (long)RDB[ptr1 + NUBAR_TOT_IDX], *E0, (double)COLL_HIT, 1.0, id);
            }

          /* score sensitivity to chi (internally constrained) */

          if (((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_CHI)
              && (ieneNew > 0) && (ieneNewRej > 0))
            {

              if (ng > 0)
                {
                  /* Calculate label for accepted */

                  gidx = CompressSensLabel(imat, izai, (long)RDB[ptr1 + CHI_DEL_IDX], ieneNew, COLL_HIT);

                  /* Create accepted event */

                  StoreSensEvent(new, gidx, Eacc, 1.0, id);

                  /* Calculate label for rejected */

                  gidx = CompressSensLabel(imat, izai, (long)RDB[ptr1 + CHI_DEL_IDX], ieneNewRej, COLL_MISS);

                  /* Create rejected event */

                  StoreSensEvent(new, gidx, Erej, 1.0, id);

                  /* Score custom perturbations for delayed chi */

                  ScorePerturbations(new, mat, rea, (long)RDB[ptr1 + CHI_DEL_IDX], Eacc, (double)COLL_HIT, 1.0, id);
                  ScorePerturbations(new, mat, rea, (long)RDB[ptr1 + CHI_DEL_IDX], Erej, (double)COLL_MISS, 1.0, id);
                }
              else
                {
                  /* Calculate label for accepted */

                  gidx = CompressSensLabel(imat, izai, (long)RDB[ptr1 + CHI_PRO_IDX], ieneNew, COLL_HIT);

                  /* Create accepted event */

                  StoreSensEvent(new, gidx, Eacc, 1.0, id);

                  /* Calculate label for rejected */

                  gidx = CompressSensLabel(imat, izai, (long)RDB[ptr1 + CHI_PRO_IDX], ieneNewRej, COLL_MISS);

                  /* Create rejected event */

                  StoreSensEvent(new, gidx, Erej, 1.0, id);

                  /* Score custom perturbations for prompt chi */

                  ScorePerturbations(new, mat, rea, (long)RDB[ptr1 + CHI_PRO_IDX], Eacc, (double)COLL_HIT, 1.0, id);
                  ScorePerturbations(new, mat, rea, (long)RDB[ptr1 + CHI_PRO_IDX], Erej, (double)COLL_MISS, 1.0, id);
                }

              /* Score custom perturbations for total chi */

              ScorePerturbations(new, mat, rea, (long)RDB[ptr1 + CHI_TOT_IDX], Eacc, (double)COLL_HIT, 1.0, id);
              ScorePerturbations(new, mat, rea, (long)RDB[ptr1 + CHI_TOT_IDX], Erej, (double)COLL_MISS, 1.0, id);
            }
        }

      /******************************************************/

#ifdef OLD_IFP

      /* Put lifetimes and delayed neutron groups to progenies */

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
                  WDB[prg + FISS_PROG_NA_COLL] =
                    WDB[prev + FISS_PROG_NA_COLL];
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
                  WDB[prg + FISS_PROG_NA_COLL] =
                    WDB[part + PARTICLE_NA_COLL];
                }

              /* Previous progeny */

              prg = prev;
            }
        }

#endif

      /* Put time */

      if (mode == SIMULATION_MODE_CRIT)
        {
          WDB[new + PARTICLE_T0] = 0.0;
          WDB[new + PARTICLE_T] = 0.0;
        }
      else
        {
          WDB[new + PARTICLE_T0] = t0 + td;
          WDB[new + PARTICLE_T] = t0 + td;
        }

      /* Delayed neutron emission time */

      WDB[new + PARTICLE_TD] = td;

      /* Put fission matrix index */

      WDB[new + PARTICLE_FMTX_IDX] = (double)idx1;

      /* Reset thermalization time */

      WDB[new + PARTICLE_TT] = 0.0;

      /* Update generation index */

      WDB[new + PARTICLE_GEN_IDX] = RDB[new + PARTICLE_GEN_IDX] + 1.0;

      /* Reset collision index */

      WDB[new + PARTICLE_COL_IDX] = 0.0;

      /* Check mode and store particle */

      if (mode == SIMULATION_MODE_CRIT)
        {
          /* Reset detector flags */

          WDB[new + PARTICLE_DET_FLAGS] = 0.0;

          /* Score response matrix source term */

          ScoreRMXSrc(new, x, y, z, E, *wgt2, id);

          /* Check Wielandt method */

          if ((long)RDB[DATA_WIELANDT_MODE] != WIELANDT_MODE_NONE)
            {
              /* Bank or que neutron */

              if (RandF(id) < RDB[DATA_WIELANDT_P])
                ToBank(new, id);
              else
                ToQue(new, id);
            }
          else
            ToBank(new, id);

          /* Write source file */

          if ((RDB[DATA_CYCLE_IDX] > RDB[DATA_CRIT_SKIP]) &&
              ((ptr = (long)RDB[DATA_PTR_CRIT_SRC_DET]) > VALID_PTR))
            WriteSourceFile(ptr, x, y, z, *u, *v, *w, E, 1.0, t0 + td, -1.0,
                            id);
        }
      else
        {
          /* Bank or que neutron */

          if (RDB[new + PARTICLE_T] < RDB[DATA_TIME_CUT_TMIN])
            Die(FUNCTION_NAME, "Error in time");
          else if (RDB[new + PARTICLE_T] >= RDB[DATA_TIME_CUT_TMAX])
            ToBank(new, id);
          else
            ToQue(new, id);

          /* Score secondary prompt source for time-dependent simulations */

          if (ng == 0)
            {
              /* Get particle emission time */

              t = RDB[new + PARTICLE_T];

              /* Score time-dependent source rate detectors */

              ScoreTimeSource(part, mat, MT_SECONDARY_PROMPT_SOURCE,
                              x, y, z, E, t, 1, *wgt2, id);
            }

          /* Score particle balance (in criticality source mode this is */
          /* done in normalizecritsource.c) */

          ptr = (long)RDB[RES_N_BALA_SRC];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf(1.0, 1.0, ptr, id, -1, BALA_N_SRC_FISS, 0);
          AddBuf(*wgt2, 1.0, ptr, id, -1, BALA_N_SRC_FISS, 1);
        }

      /* Get generation index */

      i = (long)RDB[new + PARTICLE_GEN_IDX];

      /* Score new source weight */

      ptr = (long)RDB[RES_NEW_SRC_WGT];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(*wgt2, 1.0, ptr, id, 0);

      if (i < MAX_EXT_K_GEN + 1)
        AddBuf1D(*wgt2, 1.0, ptr, id, i);
    }
}

/*****************************************************************************/
