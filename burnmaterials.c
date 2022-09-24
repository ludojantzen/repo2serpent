/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : burnmaterials.c                                */
/*                                                                           */
/* Created:       2011/05/22 (JLe)                                           */
/* Last modified: 2019/11/08 (ARi)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Performs burnup calculation for materials                    */
/*                                                                           */
/* Comments: - OpenMP parallelization revised 4.6.2012 (2.1.6)               */
/*           - Added separate subroutine to be used with corrector iteration */
/*             for convergence criterion calculation 4.11.2014 (2.1.22)      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "BurnMaterials:"

/* Use local function to simplify OpenMP implementation */

void BurnMaterials0(long, long, long, long, long);

void BurnMaterialsMSR(long, long, long, long, long);

void BurnMaterialsCI(long, long, long, long, long);

void PrintDepMatrixMSR(long, struct ccsMatrix *, double, double *,
                       double *, long);

/*****************************************************************************/

void BurnMaterials(long dep, long step)
{
  long mat, nss, type, mode;

  /***************************************************************************/

  /***** Get parameters ******************************************************/

  /* Get the number of substeps. Constant extrapolation is forced to one */
  /* substep (no point using more) */

  if ( (long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP )
    {
      if ((long)RDB[DATA_BURN_FIT_TYPE] == PRED_TYPE_CONSTANT)
        nss = 1;
      else
        nss = (long)RDB[DATA_BURN_PRED_NSS];
    }
  else
    {
      if ((long)RDB[DATA_BURN_FIT_TYPE] == CORR_TYPE_CONSTANT)
        nss = 1;
      else
        nss = (long)RDB[DATA_BURN_CORR_NSS];
    }

  /* Get step type */

  type = (long)RDB[dep + DEP_HIS_STEP_TYPE];

  /* Set burnup mode (override if decay step ) */

  if ((type == DEP_STEP_DEC_STEP) || (type == DEP_STEP_DEC_TOT))
    mode = BUMODE_TTA;
  else
    mode = (long)RDB[DATA_BURN_BUMODE];

  /***************************************************************************/

  /***** Print output ********************************************************/

  fprintf(outp, "Running burnup calculation:\n\n");

  fprintf(outp, " - Step %ld / %ld ", (long)RDB[DATA_BURN_STEP] + 1,
          (long)RDB[DATA_BURN_TOT_STEPS] + 1);

  if (((long)RDB[DATA_BURN_CORR_TYPE] == CORR_TYPE_NONE) ||
      (type == DEP_STEP_DEC_TOT) || (type == DEP_STEP_DEC_STEP))
    fprintf(outp, "\n");
  else if ((long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP)
    fprintf(outp, "(predictor)\n");
  else
    fprintf(outp, "(corrector %ld/%ld)\n",
            (long)WDB[DATA_BURN_CI_I]+1, (long)WDB[DATA_BURN_CI_MAXI]);

  if ((type != DEP_STEP_DEC_TOT) && (type != DEP_STEP_DEC_STEP))
    {
      fprintf(outp, " - Algorithm: ");

      if ((long)RDB[DATA_BURN_PRED_TYPE] == PRED_TYPE_CONSTANT)
        fprintf(outp, "CE");
      else
        fprintf(outp, "LE");

      if ((long)RDB[DATA_BURN_CORR_TYPE] == CORR_TYPE_CONSTANT)
        fprintf(outp, "/CE\n");
      else if ((long)RDB[DATA_BURN_CORR_TYPE] == CORR_TYPE_LINEAR)
        fprintf(outp, "/LI\n");
      else if ((long)RDB[DATA_BURN_CORR_TYPE] == CORR_TYPE_QUADRATIC)
        fprintf(outp, "/QI\n");
      else
        fprintf(outp, "\n");
    }

  fprintf(outp, " - Time interval: %s\n",
          TimeIntervalStr(RDB[DATA_BURN_TIME_INTERVAL]));
  fprintf(outp, " - Burnup interval: %1.2f MWd/kgU\n",
          RDB[DATA_BURN_BURNUP_INTERVAL]);

  fprintf(outp, " - Cumulative burn time after step: %s\n",
          TimeIntervalStr(RDB[DATA_BURN_CUM_BURNTIME] +
                          RDB[DATA_BURN_TIME_INTERVAL]));
  fprintf(outp, " - Cumulative burnup after step: %1.2f MWd/kgU\n",
          RDB[DATA_BURN_CUM_BURNUP] + RDB[DATA_BURN_BURNUP_INTERVAL]);

  if (nss > 1)
    fprintf(outp, " - Interval divided into %ld substeps\n", nss);

  if ((type == DEP_STEP_DEC_STEP) || (type == DEP_STEP_DEC_TOT))
    fprintf(outp,
            " - Decay step, transmutation cross sections not calculated\n");
  else if (((type == DEP_STEP_ACT_STEP) || (type == DEP_STEP_ACT_TOT)) &&
           (step > 0))
    fprintf(outp,
            " - Activation step, transmutation cross sections not updated\n");
  else if ((long)RDB[DATA_BU_SPECTRUM_COLLAPSE] == YES)
    fprintf(outp, " - Transmutation cross sections from spectrum\n");
  else
    fprintf(outp, " - Transmutation cross sections from direct tallies\n");

  if (mode == BUMODE_TTA)
    fprintf(outp, " - Bateman equations solved by linear chains method\n");
  else if (mode == BUMODE_CRAM)
    fprintf(outp, " - Bateman equations solved by CRAM\n");

  /* Reduce private results */

  ReducePrivateRes();

  /* If using corrector iteration, calculate initial extrapolations for      */
  /* nuclide field using non-averaged flux & xs                              */

  if (RDB[DATA_BURN_SIE] == YES)
    {
      /* Set flag indicating that we are in CI branch of the program */
      /* (checked in storetransmuxs.c) */

      WDB[DATA_BURN_CI_FLAG] = (double)YES;

      /* Reset thread numbers */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Check thread number */

          WDB[mat + MATERIAL_OMP_ID] = -1.0;

          /* Next material */

          mat = NextItem(mat);
        }

      if ((long)RDB[DATA_N_BURN_MATERIALS] > 1)
        fprintf(outp, "\nBurning %ld materials (CI):\n\n",
                (long)RDB[DATA_N_BURN_MATERIALS]);
      else
        fprintf(outp, "\nBurning material (CI):\n\n");

      /* Write momentary (iteration based) fluxes to file */

      WriteCIMomFluxes(step);

      /* Print */

      PrintProgress(0, 0);

      /* Start parallel timer */

      StartTimer(TIMER_OMP_PARA);

#ifdef OPEN_MP
#pragma omp parallel private (mat)
#endif
      {
        /* Loop over materials */

        mat = (long)RDB[DATA_PTR_M0];
        while (mat > VALID_PTR)
          {
            /* Check burn flag and test parallel id's */

            if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
              if (MyParallelMat(mat, YES) == YES)
                {
                  /* Burn */

                  BurnMaterialsCI(mat, step, nss, type, mode);

                  /* Print */

                  PrintProgress(mat, 2);
                }

            /* Next material */

            mat = NextItem(mat);
          }
      }

      /* Stop parallel timer */

      StopTimer(TIMER_OMP_PARA);

      /* Print */

      PrintProgress(0, 100);


      WDB[DATA_BURN_CI_FLAG] = (double)NO;
    }

  /***************************************************************************/

  /***** Main loop ***********************************************************/

  /* Reset thread numbers */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check thread number */

      WDB[mat + MATERIAL_OMP_ID] = -1.0;

      /* Next material */

      mat = NextItem(mat);
    }

  if ((long)RDB[DATA_N_BURN_MATERIALS] > 1)
    fprintf(outp, "\nBurning %ld materials:\n\n",
            (long)RDB[DATA_N_BURN_MATERIALS]);
  else
    fprintf(outp, "\nBurning material:\n\n");

  /* Print */

  PrintProgress(0, 0);

  /* Start parallel timer */

  StartTimer(TIMER_OMP_PARA);

#ifdef OPEN_MP
#pragma omp parallel private (mat)
#endif
  {
    /* Loop over materials */

    mat = (long)RDB[DATA_PTR_M0];
    while (mat > VALID_PTR)
      {
        /* Check burn flag and inflow, and test parallel id's */
        /* NOTE: No inflo in conventional burnup calculation. */

        if (((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT) &&
            ((long)RDB[mat + MATERIAL_PTR_INFLOW] < VALID_PTR))
          if (MyParallelMat(mat, YES) == YES)
            {
              /* Check flow index and burn */

              if ((long)RDB[mat + MATERIAL_FLOW_IDX] == 0)
                {
                  /* Not involved in continuous reprocessing */

                  BurnMaterials0(mat, step, nss, type, mode);

                  /* Print */

                  PrintProgress(mat, 2);
                }
              else if ((long)RDB[mat + MATERIAL_FLOW_IDX] == 1)
                {
                  /* First material in chain */

                  BurnMaterialsMSR(mat, dep, step, type, mode);

                  /* Print */

                  PrintProgress(mat, 2);
                }
            }

        /* Next material */

        mat = NextItem(mat);
      }
  }

  /* Stop parallel timer */

  StopTimer(TIMER_OMP_PARA);

  /* Print */

  PrintProgress(0, 100);

  /***************************************************************************/
}

/*****************************************************************************/

/*****************************************************************************/

void BurnMaterials0(long mat, long step, long nss, long type, long mode)
{
  long iso, ptr, lst, i, sz, ss, id;
  double t, t1, t2, tot, *N, *N0;
  struct ccsMatrix *A;

  /* This is used to move EOS compositions to MaterialBurnup.
     The approach is a bit hacky but not that messy. (AIs) */

  double *Neos;

  /* Check divisor type */

  if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT)
    Die(FUNCTION_NAME, "Divided parent material");

  /* Check MPI id */

  if (((long)RDB[mat + MATERIAL_MPI_ID] > -1) &&
      ((long)RDB[mat + MATERIAL_MPI_ID] != mpiid))
    Die(FUNCTION_NAME, "Error in MPI id");

  /* Avoid compiler warning */

  A = NULL;
  N = NULL;
  N0 = NULL;
  Neos = NULL;

  /* Get OpenMP id */

  id = OMP_THREAD_NUM;

  /* Get time step length */

  t = RDB[DATA_BURN_TIME_INTERVAL];
  CheckValue(FUNCTION_NAME, "t", "", t, ZERO, INFTY);

  /* Check type  */

  if ((type != DEP_STEP_DEC_STEP) && (type != DEP_STEP_DEC_TOT) &&
      (type != DEP_STEP_ACT_STEP) && (type != DEP_STEP_ACT_TOT))
    {
      /* Calculate transmutation cross sections */

      CalculateTransmuXS(mat, id);

      /* Store the xs from above */

      StoreTransmuXS(mat, step, type, id);
    }

  /* Size of composition vector */

  ptr = (long)RDB[mat + MATERIAL_PTR_COMP];
  sz = ListSize(ptr);

  /* Allocate memory for composition vectors */

  N0 = (double *)Mem(MEM_ALLOC, sz, sizeof(double));
  Neos = (double *)Mem(MEM_ALLOC, sz, sizeof(double));
  /*

  N0 = WorkArray(DATA_PTR_WORK_COMP1, PRIVA_ARRAY, sz, id);
  */
  /* Copy composition to N0. On predictor also store it for corrector */

  i = 0;
  lst = (long)RDB[mat + MATERIAL_PTR_COMP];

  if ((long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP)
    {
      /* Predictor step, copy composition to N0 and BOS for corrector */

      while ((iso = ListPtr(lst, i)) > VALID_PTR)
        {
          WDB[iso + COMPOSITION_ADENS_BOS] =
            RDB[iso + COMPOSITION_ADENS];
          N0[i++] = RDB[iso + COMPOSITION_ADENS_BOS];
        }
    }
  else
    {
      /* Corrector step */

      while ((iso = ListPtr(lst, i)) > VALID_PTR)
        {
          /* this should get the predicted EOS atomic densities */

          Neos[i] = RDB[iso + COMPOSITION_ADENS];

          /* copy composition from BOS to N0 */

          N0[i++] = RDB[iso + COMPOSITION_ADENS_BOS];
        }
    }

  /* Check size */

  if (i != sz)
    Die(FUNCTION_NAME, "Mismatch in size");

  /***************************************************************************/

  /***** Loop over substeps **************************************************/

  for (ss = 0; ss < nss; ss++)
    {
      /* substep beginning and end times */

      t1 = t*((double)ss)/((double)nss);
      t2 = t*((double)ss + 1.0)/((double)nss);

      /* Check for OTF mode */

      if ((long)RDB[DATA_OTF_BURN_MODE] == YES)
        {
          /* Check for sub-steps */

          if (nss > 1)
            Die(FUNCTION_NAME, "OTF mode does not work with substeps");

          /* Special treatment for predictor step in OTF calculation */

          if ((long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP)
            t2 = 0.5*t2;
        }

      /* calculate weighted xs and flux (if-lauseke lisätty */
      /* 16.8.2013 / 2.1.16)*/

      if ((type != DEP_STEP_DEC_STEP) && (type != DEP_STEP_DEC_TOT))
        AverageTransmuXS(mat, t1, t2, id);

      /* Create burnup matrix */

      A = MakeBurnMatrix(mat, id);

      /* Check size (sz == i tarkistettiin jo aikaisemmin) */

      if (sz != A->n)
        Die(FUNCTION_NAME, "Mismatch in size");

      /* Print depletion matrix */

      PrintDepMatrix(mat, A, t2 - t1, N0, N0, id);

      /* Calculate material-wise burnup (JLe 21.5.2015 / 2.1.24: Tonne  */
      /* jää jotain vanhoja vuoarvoja aikaisemmilta askeleilta jotka    */
      /* antaa nollasta poikkeavan tehon jos tätä kutsutaan decay-      */
      /* stepistä. Predictor-corrector laskussa if-lauseke estää niiden */
      /* lukemisen, mutta ilman pcc:tä ne aiheuttaa kummallisia arvoja  */
      /* materiaalikohtaisiin palamiin. Pitäisi selvittää että miten    */
      /* noi nollaamattomat arvot vaikuttaa muissa tilanteissa, esim.   */
      /* jäähtymisajan jälkeen ajetaan taas teholla.) */

      if ((type != DEP_STEP_DEC_STEP) && (type != DEP_STEP_DEC_TOT))
        MaterialBurnup(mat, N0, Neos, t1, t2, ss, id);

      /* Start burnup equation timers */

      ResetTimer(TIMER_BATEMAN);
      StartTimer(TIMER_BATEMAN);
      StartTimer(TIMER_BATEMAN_TOTAL);

      /* Override burnup mode if flux is very low */

      if (RDB[mat + MATERIAL_BURN_FLUX_SSA] < 1E-6)
        mode = BUMODE_TTA;

      /* Solve depletion equations */

      if (mode == BUMODE_TTA)
        N = TTA(A, N0, t2 - t1);
      else if (mode == BUMODE_CRAM)
        N = MatrixExponential(A, N0, t2 - t1);
      else
        Die(FUNCTION_NAME, "Invalid burnup mode");

      /* Stop timers */

      StopTimer(TIMER_BATEMAN);
      StopTimer(TIMER_BATEMAN_TOTAL);

      /* Print depletion matrix with final composition */

      PrintDepMatrix(mat, A, t2 - t1, N0, N, id);

      /* Check negative */

      for (i = 0; i < sz; i++)
         if (N[i] < 0.0)
          {
            if (N[i] < -1E-12)
              Warn(FUNCTION_NAME, "N[%ld] = %E\n", i, N[i]);
            N[i] = 0.0;
          }

      /* Free burnup matrix and N0 (which is no longer needed)*/

      ccsMatrixFree(A);
      Mem(MEM_FREE, N0);

      /* the final composition becomes initial for the next substep */

      N0 = N;
    }

  /* free the EOS composition array */

  Mem(MEM_FREE, Neos);

  /***************************************************************************/

  /***** Put final composition ***********************************************/

  /* Reset index and total atomic density */

  i = 0;
  tot = 0.0;

  /* Update composition */

  iso = (long)RDB[mat + MATERIAL_PTR_COMP];
  while (iso > VALID_PTR)
    {
      /* Put atomic density */
      /* with SIE, only update atomic densities on corrector */
      /* (or first predictor)                                */
      if((RDB[DATA_BURN_SIE] == (double)NO) ||
         (!((RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP) &&
            (RDB[DATA_BURN_STEP] != 0.0))))
        WDB[iso + COMPOSITION_ADENS] = N[i++];

      /* Add to total */

      tot = tot + RDB[iso + COMPOSITION_ADENS];

      /* Next nuclide */

      iso = NextItem(iso);
    }

  /* Put total atomic density */
  /* with SIE, only update atomic densities on corrector (or first predictor) */

  if((RDB[DATA_BURN_SIE] == (double)NO) ||
     (!((RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP) &&
        (RDB[DATA_BURN_STEP] != 0.0))))
    WDB[mat + MATERIAL_ADENS] = tot;

  /* Free composition vectors (at this point N = N0) */

  Mem(MEM_FREE, N);

}

/*****************************************************************************/

void BurnMaterialsCI(long mat, long step, long nss, long type, long mode)
{
  long iso, ptr, lst, i, sz, ss, id;
  double t, t1, t2, *N, *N0, n;
  struct ccsMatrix *A;

  /* Don't burn on predictor step unless it is the very first step */

  if (((long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP) && ((long)RDB[DATA_BURN_STEP] != 0))
    return;

  /* Check divisor type */

  if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT)
    Die(FUNCTION_NAME, "Divided parent material");

  /* Avoid compiler warning */

  A = NULL;
  N = NULL;
  N0 = NULL;

  /* Get OpenMP id */

  id = OMP_THREAD_NUM;

  /* Get time step length */

  t = RDB[DATA_BURN_TIME_INTERVAL];
  CheckValue(FUNCTION_NAME, "t", "", t, ZERO, INFTY);

  /* Check if decay step (if-lauseketta muutettu 30.7.2013 / 2.1.15 (JLe)) */

  if ((type != DEP_STEP_DEC_STEP) && (type != DEP_STEP_DEC_TOT))
    {
      /* Calculate transmutation cross sections */

      CalculateTransmuXS(mat, id);

      /* Store the xs from above */

      StoreTransmuXS(mat, step, type, id);
    }

  /* Size of composition vector */

  ptr = (long)RDB[mat + MATERIAL_PTR_COMP];
  sz = ListSize(ptr);

  /* Allocate memory for composition vectors */

  N0 = (double *)Mem(MEM_ALLOC, sz, sizeof(double));
  /*

  N0 = WorkArray(DATA_PTR_WORK_COMP1, PRIVA_ARRAY, sz, id);
  */
  /* Copy composition to N0. On predictor also store it for corrector */

  i = 0;
  lst = (long)RDB[mat + MATERIAL_PTR_COMP];

  if ((long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP)
    {
      /* Predictor step, returning */
      if((long)RDB[DATA_BURN_STEP] != 0)
        return;

      /* Predictor step, copy composition to N0 and for corrector */

      while ((iso = ListPtr(lst, i)) > VALID_PTR)
          N0[i++] = RDB[iso + COMPOSITION_ADENS];

    }
  else
    {
      /* Corrector step, copy composition from BOS to N0 */

      while ((iso = ListPtr(lst, i)) > VALID_PTR)
        N0[i++] = RDB[iso + COMPOSITION_ADENS_BOS];
    }

  /* Check size */

  if (i != sz)
    Die(FUNCTION_NAME, "Mismatch in size");

  /***************************************************************************/

  /***** Loop over substeps **************************************************/

  for (ss = 0; ss < nss; ss++)
    {
      /* substep beginning and end times */

      t1 = t*((double)ss)/((double)nss);
      t2 = t*((double)ss + 1.0)/((double)nss);

      /* calculate weighted xs and flux */

      AverageTransmuXS(mat, t1, t2, id);

      /* Create burnup matrix */

      A = MakeBurnMatrix(mat, id);

      /* Check size (sz == i tarkistettiin jo aikaisemmin) */

      if (sz != A->n)
        Die(FUNCTION_NAME, "Mismatch in size");

      /* Print depletion matrix */

      PrintDepMatrix(mat, A, t2 - t1, N0, N0, id);

      /* Calculate material-wise burnup */
      /* Not needed here */
      /*MaterialBurnup(mat, N0, t2 - t1, id); */

      /* Start burnup equation timers */

      ResetTimer(TIMER_BATEMAN);
      StartTimer(TIMER_BATEMAN);
      StartTimer(TIMER_BATEMAN_TOTAL);

      /* Override burnup mode if flux is very low */

      if (RDB[mat + MATERIAL_BURN_FLUX_SSA] < 1E-6)
        mode = BUMODE_TTA;

      /* Solve depletion equations */

      if (mode == BUMODE_TTA)
        N = TTA(A, N0, t2 - t1);
      else if (mode == BUMODE_CRAM)
        N = MatrixExponential(A, N0, t2 - t1);
      else
        Die(FUNCTION_NAME, "Invalid burnup mode");

      /* Stop timers */

      StopTimer(TIMER_BATEMAN);
      StopTimer(TIMER_BATEMAN_TOTAL);

      /* Print depletion matrix with final composition */

      PrintDepMatrix(mat, A, t2 - t1, N0, N, id);

      /* Check negative */

      for (i = 0; i < sz; i++)
        if (N[i] < 0.0)
          {
            if (N[i] < -1E-12)
              Warn(FUNCTION_NAME, "N[%ld] = %E\n", i, N[i]);
            N[i] = 0.0;
          }

      /* Free burnup matrix and N0 (which is no longer needed)*/

      ccsMatrixFree(A);

      /* Free N0 matrix */

      Mem(MEM_FREE, N0);

      /* the final composition becomes initial for the next substep */

      N0 = N;
    }

  /***************************************************************************/

  UpdateCIStop(mat, N, id);

  i = 0;
  lst = (long)RDB[mat + MATERIAL_PTR_COMP];
  n=RDB[DATA_BURN_CI_I]+1.0;

  /* Store averaged atomic density                   */
  /* Used for calculating stopping criterion         */
  /* Could also be used for nuclide field relaxation */

  /* NOTE: Tuo poistettiin sieltä rakenteesta kun sitä ei enää käytetä */
  /*       (JLe / 2.1.31 / 2.3.2018) */

  /*
  while ((iso = ListPtr(lst, i)) > VALID_PTR)
    {
      WDB[iso + COMPOSITION_ADENS_AVE] = (n-1)/n*RDB[iso + COMPOSITION_ADENS_AVE] + 1/n*N[i++];
    }
  */
  /* Return without putting new composition */

  Mem(MEM_FREE, N);

}

/*****************************************************************************/

void BurnMaterialsMSR(long mat, long dep, long step, long type, long mode)
{
  long iso, ptr, i, j, sz, nm, id, mat0;
  double t, t1, t2, tot, *N, *N0, **R;
  struct ccsMatrix *A;

  /* Check divisor type */

  if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT)
    Die(FUNCTION_NAME, "Divided parent material");

  /* Avoid compiler warning */

  A = NULL;
  N = NULL;
  N0 = NULL;
  R = NULL;

  /* Get OpenMP id */

  id = OMP_THREAD_NUM;

  /* Get time step length */

  t = RDB[DATA_BURN_TIME_INTERVAL];
  CheckValue(FUNCTION_NAME, "t", "", t, ZERO, INFTY);

  /* Size of composition vector (NOTE: All materials in the chain must have */
  /* the exact same nuclides */

  ptr = (long)RDB[mat + MATERIAL_PTR_COMP];
  sz = ListSize(ptr);

  /* Multiply by number of materials in chain */

  nm = (long)RDB[mat + MATERIAL_FLOW_N];
  sz = sz*nm;

  /* Allocate memory for composition vector */

  N0 = (double *)Mem(MEM_ALLOC, sz + 1, sizeof(double));

  /* Check that predictor-corrector calculation is not in use */

  if ((long)RDB[DATA_BURN_STEP_PC] != PREDICTOR_STEP)
    Die(FUNCTION_NAME, "Doesn't work with predictor-corrector");

  /* Avoid compiler warning */

  i = -1;

  /* Copy composition to N0, loop over materials */

  mat0 = (long)RDB[DATA_PTR_M0];
  while (mat0 > VALID_PTR)
    {
      /* Check that material is linked to first in chain */

      if ((long)RDB[mat0 + MATERIAL_FLOW_PTR_FIRST] != mat)
        {
          /* Pointer to next */

          mat0 = NextItem(mat0);

          /* Cycle loop */

          continue;
        }

      /* Get material index */

      j = (long)RDB[mat0 + MATERIAL_FLOW_IDX] - 1;
      CheckValue(FUNCTION_NAME, "j", "", j, 0, 20);

      /* Reset array index */

      i = 0;

      /* Loop over composition */

      iso = (long)RDB[mat0 + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
        {
          /* Check index */

          if (i*nm + j > sz + 1)
            Die(FUNCTION_NAME, "Mismatch in size 1");

          /* Copy composition to N0 */

          if (N0[i*nm + j] != 0.0)
            Die(FUNCTION_NAME, "Not zero");
          else
            N0[i*nm + j] = RDB[iso + COMPOSITION_ADENS]*
              RDB[mat0 + MATERIAL_VOLUME];

          /* Update index */

          i++;

          /* Next */

          iso = NextItem(iso);
        }

      /* Pointer to next material */

      mat0 = NextItem(mat0);
    }

  /* Check indexing and set dummy concentration for constant flow rate */

  if (i*nm != sz)
    Die(FUNCTION_NAME, "Indexing error");
  else
    N0[sz] = 1.0;

  /* Beginning and end times */

  t1 = 0.0;
  t2 = t;

  /* Calculate weighted xs and flux */

  mat0 = (long)RDB[DATA_PTR_M0];
  while (mat0 > VALID_PTR)
    {
      /* Check that material is linked to first in chain */

      if ((long)RDB[mat0 + MATERIAL_FLOW_PTR_FIRST] != mat)
        {
          /* Pointer to next */

          mat0 = NextItem(mat0);

          /* Cycle loop */

          continue;
        }

      if ((type != DEP_STEP_DEC_STEP) && (type != DEP_STEP_DEC_TOT))
        {
          /* Calculate transmutation cross sections */

          CalculateTransmuXS(mat0, id);

          /* Store the xs from above */

          StoreTransmuXS(mat0, step, type, id);
        }

      /* Calculate average */

      AverageTransmuXS(mat0, t1, t2, id);

      /* Put data in table */

      R = MSRReaList(R, mat0, id);

      /* Next material */

      mat0 = NextItem(mat0);
    }

  /* Create burnup matrix */

  A = MakeBurnMatrixMSR(mat, dep, R, id);

  /* Check size */

  if (sz + 1 != A->n)
    Die(FUNCTION_NAME, "Mismatch in size 2");

  /* Print depletion matrix */

  PrintDepMatrixMSR(mat, A, t2 - t1, N0, N0, sz + 1);

  /* Calculate material-wise burnup */
  /*
  MaterialBurnup(mat, N0, Neos, t1, t2, ss, id);
  */
  /* Start burnup equation timers */

  ResetTimer(TIMER_BATEMAN);
  StartTimer(TIMER_BATEMAN);
  StartTimer(TIMER_BATEMAN_TOTAL);

  /* Override burnup mode if flux is very low */
  /*
  if (RDB[mat + MATERIAL_BURN_FLUX_SSA] < 1E-6)
    mode = BUMODE_TTA;
  */
  /* Solve depletion equations */

  if (mode == BUMODE_TTA)
    N = TTA(A, N0, t2 - t1);
  else if (mode == BUMODE_CRAM)
    N = MatrixExponential(A, N0, t2 - t1);
  else
    Die(FUNCTION_NAME, "Invalid burnup mode");

  /* Reset dummy */

  N[sz] = 0.0;

  /* Stop timers */

  StopTimer(TIMER_BATEMAN);
  StopTimer(TIMER_BATEMAN_TOTAL);

  /* Print depletion matrix with final composition */

  PrintDepMatrixMSR(mat, A, t2 - t1, N0, N, sz + 1);

  /* Check negative */

  for (i = 0; i < sz; i++)
    if (N[i] < 0.0)
      {
        /* Negative values may result from continuous flows */

        if (N[i] < -1E-12)
          {
            /* Check step type */

            if ((type == DEP_STEP_BU_STEP) || (type == DEP_STEP_BU_TOT))
              Error(dep,
                    "Numerical error, try smaller step size (max %1.1f MWd/kgU)",
                    RDB[DATA_BURN_BURNUP_INTERVAL]*N0[i]/(N0[i] - N[i]));
            else
              Error(dep,
                    "Numerical error, try smaller step size (max %1.1f days)",
                    RDB[DATA_BURN_TIME_INTERVAL]*N0[i]/(N0[i] - N[i])
                    /(60.0*60.0*24.0));
          }

        N[i] = 0.0;
      }

  /* Free burnup matrix and N0 (which is no longer needed)*/

  ccsMatrixFree(A);
  Mem(MEM_FREE, N0);

  /***************************************************************************/

  /***** Put final composition ***********************************************/

  /* Loop over materials */

  mat0 = (long)RDB[DATA_PTR_M0];
  while (mat0 > VALID_PTR)
    {
      /* Check that material is linked to first in chain */

      if ((long)RDB[mat0 + MATERIAL_FLOW_PTR_FIRST] != mat)
        {
          /* Pointer to next */

          mat0 = NextItem(mat0);

          /* Cycle loop */

          continue;
        }

      /* Get material index */

      j = (long)RDB[mat0 + MATERIAL_FLOW_IDX] - 1;
      CheckValue(FUNCTION_NAME, "j", "", j, 0, 20);

      /* Reset array index */

      i = 0;

      /* Loop over composition */

      iso = (long)RDB[mat0 + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
        {
          /* Check index */

          if (i*nm + j > sz + 1)
            Die(FUNCTION_NAME, "Mismatch in size 3");

          /* Copy composition to structure */

          WDB[iso + COMPOSITION_ADENS] =
            N[i*nm + j]/RDB[mat0 + MATERIAL_VOLUME];

          /* Update index */

          i++;

          /* Next */

          iso = NextItem(iso);
        }

      /* Pointer to next material */

      mat0 = NextItem(mat0);
    }

  /* Calculate totals */

  mat0 = (long)RDB[DATA_PTR_M0];
  while (mat0 > VALID_PTR)
    {
      /* Check that material is linked to first in chain */

      if ((long)RDB[mat0 + MATERIAL_FLOW_PTR_FIRST] != mat)
        {
          /* Pointer to next */

          mat0 = NextItem(mat0);

          /* Cycle loop */

          continue;
        }

      /* Reset total */

      tot = 0.0;

      /* Loop over composition */

      iso = (long)RDB[mat0 + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
        {
          /* Add to total */

          tot = tot + RDB[iso + COMPOSITION_ADENS];

          /* Next */

          iso = NextItem(iso);
        }

      /* Put total atomic density */

      WDB[mat0 + MATERIAL_ADENS] = tot;

      /* Pointer to next material */

      mat0 = NextItem(mat0);
    }

  /* Free R */

  for (i = 0; i < 4; i++)
    Mem(MEM_FREE, R[i]);

  Mem(MEM_FREE, R);

  /* Free composition vectors (at this point N = N0) */

  Mem(MEM_FREE, N);
}

/*****************************************************************************/

void PrintDepMatrixMSR(long mat, struct ccsMatrix *A, double t, double *N0,
                       double *N1, long sz)
{
  long i, j;
  char tmpstr[MAX_STR];
  FILE *fp;

  /* Check print flag */

  if ((long)RDB[DATA_BURN_PRINT_DEPMTX] == NO)
    return;

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* File name */

  sprintf(tmpstr,"%s_depmtx_%s%ld.m", GetText(DATA_PTR_INPUT_FNAME),
          GetText(mat + MATERIAL_PTR_NAME), (long)RDB[DATA_BURN_STEP]);

  /* Open file for writing */

  if ((fp = fopen(tmpstr, "w")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file for writing");

  /* Check that matrix is square */

  if (A->n != A->m)
    Die(FUNCTION_NAME, "Matrix is not square");

  /* Print time step */

  /*  fprintf(fp, "t = %E;\n", t);*/
  fprintf(fp, "t = %22.15E;\n", t);

  /* Print size */

  /*
  fprintf(fp, "%ld %ld\n", A->n, A->m);
  */

  /* Print flux */

  fprintf(fp, "flx = %22.15E;\n", Truncate(RDB[mat + MATERIAL_BURN_FLUX_SSA], 6));

  /* Loop over composition and print */

  fprintf(fp, "N0 = [\n");
  for (i = 0; i < sz; i++)
    fprintf(fp, "%E\n", N0[i]);
  fprintf(fp,"];\n");

  fprintf(fp, "N1 = [\n");
  for (i = 0; i < sz; i++)
    fprintf(fp, "%E\n", N1[i]);
  fprintf(fp,"];\n");

  /* Print matrix size and number of non-zero elements */
  /*
  fprintf(fp, "%ld %ld\n", A->n, A->nnz);
  */
  /* Print matrix */

  fprintf(fp, "A = zeros(%ld, %ld);\n", A->n, A->m);

  for (i = 0; i < A->m; i++)
    for (j = A->colptr[i]; j < A->colptr[i + 1]; j++)
      /*      fprintf(fp, "A(%4ld, %4ld) = %17.10E;\n", A->rowind[j] + 1, i + 1, */
      fprintf(fp, "A(%4ld, %4ld) = %22.15E;\n", A->rowind[j] + 1, i + 1,
              A->values[j].re);


  /* Close file */

  fclose(fp);
}

/*****************************************************************************/
