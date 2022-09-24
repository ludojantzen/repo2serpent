/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : solvermx.c                                     */
/*                                                                           */
/* Created:       2016/04/28 (JLe)                                           */
/* Last modified: 2019/10/22 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Solves forward and adjoint currents and importances using    */
/*              the response matrix method.                                  */
/*                                                                           */
/* Comments: - Noita GVR ja MULTI -moodien ratkaisuja säädettiin             */
/*             8.-11.9.2018. Backuppeja löytyy myös kansioista 327-329.      */
/*             Uusin versio kertoo importancet suoraan MC-laskun virtoihin.  */
/*             Tuolla on myös sellainen versio joka laskee forward-laskun    */
/*             kertoo ne siellä.                                             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SolveRMX:"

/*****************************************************************************/

void SolveRMX(long rmx, long refe)
{
  long loc0, ptr, det, loc1, loc2, sz, i, j, k, it, imax, nmax, compl;
  long n, m, ng;
  double flim, max, sum, chk, val, *R, *Jin, *Jout, *I, bala, Rtot, frac;
  double dir, wgt;
  const double *s, *r, *S, *rs, *alpha;

  /* Check that mesh exists */

  if ((loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA]) < VALID_PTR)
    Die(FUNCTION_NAME, "no mesh data");

  /* Get number of mesh cells */

  sz = ListSize(loc0);
  CheckValue(FUNCTION_NAME, "sz", "", sz, 2, 1000000000);

  /* Get number of energy groups */

  ng = (long)RDB[rmx + RMX_NG];
  CheckValue(FUNCTION_NAME, "ng", "", ng, 1, 100000);

  /***************************************************************************/

  /***** Iteration loop ******************************************************/

  /* Set convergence limit */

  flim = RDB[rmx + RMX_CONV];

  /* Set maximum number of iterations and reset completion */

  imax = (long)RDB[rmx + RMX_N_ITER];
  CheckValue(FUNCTION_NAME, "imax", "", imax, -1, 1000000);

  /* Print (refe is set to -1 to supress printing this) */

  if ((refe == -1) && ((long)RDB[DATA_RMTX_SOLVE_FORWARD] == YES))
    refe = NO;
  else if (refe != -1)
    {
      if (imax > 0)
        fprintf(outp, "Number of iterations : %ld\n", imax);

      if (flim > 0.0)
        fprintf(outp, "Convergence limit : %1.2E\n", flim);
    }

  /***************************************************************************/

  /***** Forward solution ****************************************************/

  if ((long)RDB[DATA_RMTX_SOLVE_FORWARD] == YES)
    {
      /* Check that test mode is on */

      if ((long)RDB[DATA_RMX_TEST_MODE] == NO)
        Die(FUNCTION_NAME, "Forward solution without test mode");

      /* Start timer */

      ResetTimer(TIMER_MISC);
      StartTimer(TIMER_MISC);

      /* Reset completed flag */

      compl = 0;

      /* Print */

      fprintf(outp, "\nIterating forward solution:\n\n");

      /***********************************************************************/

      /***** Source to forward currents **************************************/

      /* Reset current balance */

      bala = 0.0;

      /* Loop over mesh */

      loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      while (loc0 > VALID_PTR)
        {
          /* Get matrix size */

          nmax = (long)RDB[loc0 + RMX_CELL_MTX_SIZE];
          CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

          /* Check reference solution mode */

          if (refe == YES)
            {
              /* Set current balance */

              bala = 1.0;

              /* Get outward current vector */

              ptr = (long)RDB[loc0 + RMX_CELL_WRK_OUT_CURR];
              CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
              Jout = &WDB[ptr];

              /* Reset */

              memset(Jout, 0.0, nmax*ng*sizeof(double));
            }
          else
            {
              /* Get source coefficient vector */

              ptr = (long)RDB[loc0 + RMX_CELL_COEF_FWD_SRCC];
              CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
              s = &RDB[ptr];

              /* Get total source rate */

              ptr = (long)RDB[loc0 + RMX_CELL_MC_SRCC_TOT];
              CheckPointer(FUNCTION_NAME, "ptr", RES2_ARRAY, ptr);
              S = &RES2[ptr];

              /* Get outward current vector */

              ptr = (long)RDB[loc0 + RMX_CELL_WRK_OUT_CURR];
              CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
              Jout = &WDB[ptr];

              /* Reset */

              memset(Jout, 0.0, nmax*ng*sizeof(double));

              /* Copy values */

              for (n = 0; n < ng; n++)
                for (i = 0; i < nmax; i++)
                  for (m = 0; m < ng; m++)
                    if (S[m] > 0.0)
                      {
                        /* Put source particles to outward currents */

                        Jout[i*ng + n] = Jout[i*ng + n] +
                          s[i*ng*ng + m*ng + n]*S[m];

                        /* Add to current balance */

                        bala = bala + s[i*ng*ng + m*ng + n]*S[m];
                      }

              /* Get inward current vector */

              ptr = (long)RDB[loc0 + RMX_CELL_WRK_IN_CURR];
              CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
              Jin = &WDB[ptr];

              /* Reset */

              memset(Jin, 0.0, nmax*ng*sizeof(double));
            }

          /* Reset solution vector */

          ptr = (long)RDB[loc0 + RMX_CELL_FWD_SOL_IN_CURR];
          CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
          memset(&WDB[ptr], 0.0, nmax*ng*sizeof(double));

          /* Next */

          loc0 = NextItem(loc0);
        }

      /***********************************************************************/

      /***** Main loop *******************************************************/

      for (it = 0; it != imax; it++)
        {
          /*******************************************************************/

          /***** Distribute forward current to neighbours ********************/

          /* Skip this for cycle in reference mode */

          if ((refe == NO) || (it > 0))
            {
#ifdef OPEN_MP
#pragma omp parallel private(loc0, loc1, loc2, ptr, i, j, k, n, m, nmax, Jin, Jout)
#endif
              {

#ifdef OPEN_MP
#pragma omp for
#endif
                /* Loop over mesh */

                for (k = 0; k < sz; k++)
                  {
                    /* Pointer to data */

                    loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
                    CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

                    loc0 = ListPtr(loc0, k);
                    CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

                    /* Check empty flag */

                    if ((long)RDB[loc0 + RMX_CELL_EMPTY] == YES)
                      continue;

                    /* Get matrix size */

                    nmax = (long)RDB[loc0 + RMX_CELL_MTX_SIZE];
                    CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

                    /* Get outward current vector */

                    ptr = (long)RDB[loc0 + RMX_CELL_WRK_OUT_CURR];
                    CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
                    Jout = &WDB[ptr];

                    /* Loop over neighbours */

                    loc2 = (long)RDB[loc0 + RMX_CELL_PTR_BOUNDS];
                    CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

                    for (i = 0; i < nmax; i++)
                      {
                        /* Pointer to cell */

                        loc1 = (long)RDB[loc2 + RMX_CELL_BOUND_PTR_CELL];
                        CheckPointer(FUNCTION_NAME, "loc1", DATA_ARRAY, loc1);

                        /* Check pointers */

                        if (loc0 == loc1)
                          Die(FUNCTION_NAME, "Pointer to self");

                        /* Boundary index */

                        j = (long)RDB[loc2 + RMX_CELL_BOUND_FWD_IDX];

                        /* Get inward current vector */

                        ptr = (long)RDB[loc1 + RMX_CELL_WRK_IN_CURR];
                        CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
                        Jin = &WDB[ptr];

                        /* Put data */

                        for (n = 0; n < ng; n++)
                          {
                            /* Add to current */

                            Jin[j*ng + n] += Jout[i*ng + n];
                          }

                        /* Next */

                        loc2 = NextItem(loc2);
                      }
                  }
              }
            }

          /*******************************************************************/

          /***** Move forward currents from inward to outward buffer *********/

#ifdef OPEN_MP
#pragma omp parallel private(loc0, loc1, ptr, i, j, k, n, m, nmax, Jin, Jout, alpha)
#endif
          {

#ifdef OPEN_MP
#pragma omp for
#endif
            /* Loop over mesh */

            for (k = 0; k < sz; k++)
              {
                /* Pointer to data */

                loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
                CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

                loc0 = ListPtr(loc0, k);
                CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

                /* Check empty flag */

                if ((long)RDB[loc0 + RMX_CELL_EMPTY] == YES)
                  continue;

                /* Get matrix size */

                nmax = (long)RDB[loc0 + RMX_CELL_MTX_SIZE];
                CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

                /* Get inward current vector */

                ptr = (long)RDB[loc0 + RMX_CELL_WRK_IN_CURR];
                CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
                Jin = &WDB[ptr];

                /* Add inward current to solution */

                ptr = (long)RDB[loc0 + RMX_CELL_FWD_SOL_IN_CURR];
                CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

                for (n = 0; n < ng; n++)
                  for (i = 0; i < nmax; i++)
                    WDB[ptr + i*ng + n] = RDB[ptr + i*ng + n] + Jin[i*ng + n];

                /* Get outward current vector */

                ptr = (long)RDB[loc0 + RMX_CELL_WRK_OUT_CURR];
                CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
                Jout = &WDB[ptr];

                /* Reset outward current */

                memset(Jout, 0.0, nmax*ng*sizeof(double));

                /* Get coefficient matrix */

                ptr = (long)RDB[loc0 + RMX_CELL_COEF_FWD_TRANS_MTX];
                CheckPointer(FUNCTION_NAME, "(ptr4)", DATA_ARRAY, ptr);
                alpha = &RDB[ptr];

                /* Loop over matrix */

                for (n = 0; n < ng; n++)
                  for (i = 0; i < nmax; i++)
                    if (Jin[i*ng + n] != 0.0)
                      for (m = 0; m < ng; m++)
                        for (j = 0; j < nmax; j++)
                          Jout[j*ng + m] = Jout[j*ng + m] +
                            alpha[i*nmax*ng*ng + n*nmax*ng + j*ng + m]*
                            Jin[i*ng + n];

                /* Reset inward current */

                memset(Jin, 0.0, nmax*ng*sizeof(double));
              }
          }

          /*******************************************************************/

          /***** Check convergence *******************************************/

          /* Reset sum of outward currents and current fraction */

          sum = 0.0;
          frac = INFTY;

          /* Loop over mesh */

          loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
          while (loc0 > VALID_PTR)
            {
              /* Get matrix size */

              nmax = (long)RDB[loc0 + RMX_CELL_MTX_SIZE];
              CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

              /* Get outward current vector */

              ptr = (long)RDB[loc0 + RMX_CELL_WRK_OUT_CURR];
              CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
              Jout = &WDB[ptr];

              /* Add to sum */

              for (n = 0; n < ng; n++)
                for (i = 0; i < nmax; i++)
                  sum = sum + Jout[i*ng + n];

              /* Next cell */

              loc0 = NextItem(loc0);
            }

          /* Calculate remaining fraction */

          if (bala == 0.0)
            Die(FUNCTION_NAME, "Error in current balance");
          else
            frac = sum/bala;

          /* Print progress */

          if ((it > 0) && (imax < 0))
            {
              /* Only convergence limit given */

              if ((k = (long)log10(frac)) != compl)
                {
                  /* Print */

                  fprintf(outp, "%4ld iterations completed: frac = %1.2E\n",
                          it, frac);
                  compl = k;
                }
            }
          else if (it > 0)
            {
              /* Maximum number of iterations given */

              k = (long)(100.0*((double)it + 1.0)/((double)imax));
              if ((frac < flim) || (!(k % 10) && (k != compl)))
                {
                  /* Print */

                  fprintf(outp,
                          "%4ld%% of iterations completed: frac = %1.2E\n",
                          k, frac);
                  compl = k;
                }
            }

          /* Check convergence */

          if ((it > 10) && (frac < flim))
            break;

          /*******************************************************************/
        }

      /***********************************************************************/

      /***** Calculate forward solution **************************************/

      /* Reset maximum difference */

      max = 0.0;

      /* Loop over mesh */

      loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      while (loc0 > VALID_PTR)
        {
          /* Get matrix size */

          nmax = (long)RDB[loc0 + RMX_CELL_MTX_SIZE];
          CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

          /* Get total source rate */

          ptr = (long)RDB[loc0 + RMX_CELL_MC_SRCC_TOT];
          CheckPointer(FUNCTION_NAME, "ptr", RES2_ARRAY, ptr);
          S = &RES2[ptr];

          /* Loop over detectors */

          det = (long)RDB[loc0 + RMX_CELL_PTR_DET];
          while (det > VALID_PTR)
            {
              /* Get direct source-to-response coefficient */

              ptr = (long)RDB[det + RMX_DET_COEF_DIR_SRCC_RES];
              CheckPointer(FUNCTION_NAME, "(ptr5)", DATA_ARRAY, ptr);
              rs = &RDB[ptr];

              /* Put direct contribution */

              WDB[det + RMX_DET_FWD_CHK] = 0.0;
              for (n = 0; n < ng; n++)
                WDB[det + RMX_DET_FWD_CHK] =
                  RDB[det + RMX_DET_FWD_CHK] + rs[n]*S[n];

              /* Get response vector */

              ptr = (long)RDB[det + RMX_DET_COEF_FWD_RES];
              CheckPointer(FUNCTION_NAME, "(ptr5)", DATA_ARRAY, ptr);
              r = &RDB[ptr];

              /* Get inward current solution */

              ptr = (long)RDB[loc0 + RMX_CELL_FWD_SOL_IN_CURR];
              CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
              Jin = &WDB[ptr];

              /* Add contribution */

              for (n = 0; n < ng; n++)
                for (i = 0; i < nmax; i++)
                  WDB[det + RMX_DET_FWD_CHK] = RDB[det + RMX_DET_FWD_CHK]
                    + r[i*ng + n]*Jin[i*ng + n];

              /* Reset total response */

              Rtot = 0.0;

              /* Add to total response */

              ptr = (long)RDB[det + RMX_DET_MC_RES_TOT];
              CheckPointer(FUNCTION_NAME, "(ptr6)", RES2_ARRAY, ptr);

              for (n = 0; n < ng; n++)
                Rtot = Rtot + RES2[ptr + n];

              /* Compare differences to maximum */

              if (Rtot != 0.0)
                {
                  /* Calculate difference */

                  val = RDB[det + RMX_DET_FWD_CHK]/Rtot - 1.0;

                  /* Compare absolute values */

                  if (fabs(val) > fabs(max))
                    max = val;
                }

              /* Next detector */

              det = NextItem(det);
            }

          /* Next */

          loc0 = NextItem(loc0);
        }

      /***********************************************************************/

      /***** Simulation completed ********************************************/

      /* Print */

      fprintf(outp, "\nCompleted in %s. ",
              TimeIntervalStr(TimerVal(TIMER_MISC)));

      /* Check number of iterations */

      if (it == imax)
        fprintf(outp,
                "Maximum number of iterations reached (no convergence).\n");
      else
        fprintf(outp, "Convergence reached after %ld iterations.\n", it);

      /* Print maximum differences */

      fprintf(outp, "Maximum relative difference to local response: %1.2E.\n",
              max);

      /* Stop timer */

      StopTimer(TIMER_MISC);

      /***********************************************************************/
    }

  /***************************************************************************/

  /***** Adjoint solution ****************************************************/

  if ((long)RDB[DATA_RMTX_SOLVE_ADJOINT] == YES)
    {
      /* Start timer */

      ResetTimer(TIMER_MISC);
      StartTimer(TIMER_MISC);

      /* Reset completed flag */

      compl = 0;

      /* Print */

      if (refe == -1)
        fprintf(outp, "Iterating adjoint solution:\n\n");
      else
        fprintf(outp, "\nIterating adjoint solution:\n\n");

      /***********************************************************************/

      /***** Response to adjoint currents ************************************/

      /* Reset current balance */

      bala = 0.0;

      /* Loop over mesh */

      loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      while (loc0 > VALID_PTR)
        {
          /* Get matrix size */

          nmax = (long)RDB[loc0 + RMX_CELL_MTX_SIZE];
          CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

          /* Get total source rate */

          ptr = (long)RDB[loc0 + RMX_CELL_MC_SRCC_TOT];
          CheckPointer(FUNCTION_NAME, "ptr", RES2_ARRAY, ptr);
          S = &RES2[ptr];

          /* Get inward current vector */

          ptr = (long)RDB[loc0 + RMX_CELL_WRK_IN_CURR];
          CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
          Jin = &WDB[ptr];

          /* Reset */

          memset(Jin, 0.0, nmax*ng*sizeof(double));

          /* Loop over detectors */

          det = (long)RDB[loc0 + RMX_CELL_PTR_DET];
          while (det > VALID_PTR)
            {
              /* Get weighting factor */

              wgt = RDB[det + RMX_DET_WGT];
              CheckValue(FUNCTION_NAME, "wgt", "", wgt, 0.0, INFTY);

              /* Get direct source-to-response coefficient */

              ptr = (long)RDB[det + RMX_DET_COEF_DIR_SRCC_RES];
              CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
              rs = &RDB[ptr];

              /* Get response coefficient vector */

              ptr = (long)RDB[det + RMX_DET_COEF_ADJ_RES];
              CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
              r = &RDB[ptr];

              /* Get total response */

              ptr = (long)RDB[det + RMX_DET_MC_RES_TOT];
              CheckPointer(FUNCTION_NAME, "(ptr7)", RES2_ARRAY, ptr);
              R = &RES2[ptr];

              /* Add to balance (toi että lähteen kontribuutiota ei  */
              /* miinusteta kuin alla supistaa pois sen erotustermin */
              /* ja tulos on kontribuutio responssien summaan).      */

              for (n = 0; n < ng; n++)
                for (i = 0; i < nmax; i++)
                  if (r[i*ng + n] > 0.0)
                    bala = bala + R[n]*wgt/r[i*ng + n];

              /* Copy values (source contribution not included) */

              for (n = 0; n < ng; n++)
                for (i = 0; i < nmax; i++)
                  if (r[i*ng + n] > 0.0)
                    Jin[i*ng + n] = Jin[i*ng + n] +
                      (R[n] - S[n]*rs[n])*wgt/r[i*ng + n];

              /* Next detector */

              det = NextItem(det);
            }

          /* Get outward current vector */

          ptr = (long)RDB[loc0 + RMX_CELL_WRK_OUT_CURR];
          CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
          Jout = &WDB[ptr];

          /* Reset */

          memset(Jout, 0.0, nmax*ng*sizeof(double));

          /* Reset solution vectors */

          ptr = (long)RDB[loc0 + RMX_CELL_ADJ_SOL_IN_CURR];
          CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
          memset(&WDB[ptr], 0.0, nmax*ng*sizeof(double));

          ptr = (long)RDB[loc0 + RMX_CELL_ADJ_SOL_OUT_CURR];
          CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
          memset(&WDB[ptr], 0.0, nmax*ng*sizeof(double));

          /* Next */

          loc0 = NextItem(loc0);
        }

      /***********************************************************************/

      /* Main loop */

      for (it = 0; it != imax; it++)
        {
          /*******************************************************************/

          /***** Distribute adjoint current to neighbours ********************/

#ifdef OPEN_MP
#pragma omp parallel private(loc0, loc1, loc2, ptr, i, j, k, n, m, nmax, Jin, Jout)
#endif
          {

#ifdef OPEN_MP
#pragma omp for
#endif
            /* Loop over mesh */

            for (k = 0; k < sz; k++)
              {
                /* Pointer to data */

                loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
                CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

                loc0 = ListPtr(loc0, k);
                CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

                /* Check empty flag */

                if ((long)RDB[loc0 + RMX_CELL_EMPTY] == YES)
                  continue;

                /* Get matrix size */

                nmax = (long)RDB[loc0 + RMX_CELL_MTX_SIZE];
                CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

                /* Get inward current vector */

                ptr = (long)RDB[loc0 + RMX_CELL_WRK_IN_CURR];
                CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
                Jin = &WDB[ptr];

                /* Add inward current to solution */

                ptr = (long)RDB[loc0 + RMX_CELL_ADJ_SOL_IN_CURR];
                CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

                for (n = 0; n < ng; n++)
                  for (i = 0; i < nmax; i++)
                    if (Jin[i*ng + n] > 0.0)
                      WDB[ptr + i*ng + n] = RDB[ptr + i*ng + n]
                        + Jin[i*ng + n];

                /* Loop over neighbours */

                loc2 = (long)RDB[loc0 + RMX_CELL_PTR_BOUNDS];
                CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

                for (i = 0; i < nmax; i++)
                  {
                    /* Pointer to cell */

                    loc1 = (long)RDB[loc2 + RMX_CELL_BOUND_PTR_CELL];
                    CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

                    /* Check pointers */

                    if (loc0 == loc1)
                      Die(FUNCTION_NAME, "Pointer to self");

                    /* Boundary index */

                    j = (long)RDB[loc2 + RMX_CELL_BOUND_ADJ_IDX];

                    /* Get outward current vector */

                    ptr = (long)RDB[loc1 + RMX_CELL_WRK_OUT_CURR];
                    CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
                    Jout = &WDB[ptr];

                    /* Put data */

                    for (m = 0; m < ng; m++)
                      if (Jin[i*ng + m] > 0.0)
                        Jout[j*ng + m] += Jin[i*ng + m];

                    /* Next */

                    loc2 = NextItem(loc2);
                  }

                /* Reset inward current */

                memset(Jin, 0.0, nmax*ng*sizeof(double));
              }
          }

          /*******************************************************************/

          /***** Move adjoint currents from outward to inward buffer *********/

#ifdef OPEN_MP
#pragma omp parallel private(loc0, loc1, ptr, i, j, k, n, m, nmax, Jin, Jout, alpha)
#endif
          {

#ifdef OPEN_MP
#pragma omp for
#endif
            /* Loop over mesh */

            for (k = 0; k < sz; k++)
              {
                /* Pointer to data */

                loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
                CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

                loc0 = ListPtr(loc0, k);
                CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

                /* Get matrix size */

                nmax = (long)RDB[loc0 + RMX_CELL_MTX_SIZE];
                CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

                /* Get outward current vector */

                ptr = (long)RDB[loc0 + RMX_CELL_WRK_OUT_CURR];
                CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
                Jout = &WDB[ptr];

                /* Loop over groups */

                for (n = 0; n < ng; n++)
                  {
                    /* Loop over bounds */

                    for (i = 0; i < nmax; i++)
                      if (Jout[i*ng + n] > 0.0)
                        break;

                    /* Check index */

                    if (i < nmax)
                      break;
                  }

                /* Check index */

                if (n == ng)
                  continue;

                /* Set empty flag */

                WDB[loc0 + RMX_CELL_EMPTY] = (double)YES;

                /* Add outward current to solution */

                ptr = (long)RDB[loc0 + RMX_CELL_ADJ_SOL_OUT_CURR];
                CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

                for (n = 0; n < ng; n++)
                  for (i = 0; i < nmax; i++)
                    if (Jout[i*ng + n] > 0.0)
                      {
                        /* Add current */

                        WDB[ptr + i*ng + n] = RDB[ptr + i*ng + n]
                          + Jout[i*ng + n];

                        /* Reset empty flag */

                        WDB[loc0 + RMX_CELL_EMPTY] = (double)NO;
                      }

                /* Get inward current vector */

                ptr = (long)RDB[loc0 + RMX_CELL_WRK_IN_CURR];
                CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
                Jin = &WDB[ptr];

                /* Get coefficient matrix */

                ptr = (long)RDB[loc0 + RMX_CELL_COEF_ADJ_TRANS_MTX];
                CheckPointer(FUNCTION_NAME, "(ptr8)", DATA_ARRAY, ptr);
                alpha = &RDB[ptr];

                /* Loop over matrix */

                for (m = 0; m < ng; m++)
                  for (j = 0; j < nmax; j++)
                    if (Jout[j*ng + m] != 0.0)
                      for (n = 0; n < ng; n++)
                        for (i = 0; i < nmax; i++)
                          Jin[i*ng + n] = Jin[i*ng + n] +
                            alpha[i*nmax*ng*ng + n*nmax*ng + j*ng + m]*
                            Jout[j*ng + m];

                /* Reset outward current */

                memset(Jout, 0.0, nmax*ng*sizeof(double));
              }
          }

          /*******************************************************************/

          /***** Check convergence *******************************************/

          /* Reset sum of outward currents and current fraction */

          sum = 0.0;
          frac = INFTY;

          /* Loop over mesh */

          loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
          while (loc0 > VALID_PTR)
            {
              /* Get matrix size */

              nmax = (long)RDB[loc0 + RMX_CELL_MTX_SIZE];
              CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

              /* Get inward current vector */

              ptr = (long)RDB[loc0 + RMX_CELL_WRK_IN_CURR];
              CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
              Jin = &WDB[ptr];

              /* Calculate total inward current */

              chk = 0.0;
              for (n = 0; n < ng; n++)
                for (i = 0; i < nmax; i++)
                  chk = chk + Jin[i*ng + n];

               /* Add to sum */

              sum = sum + chk;

              /* Next cell */

              loc0 = NextItem(loc0);
            }

          /* Calculate remaining fraction */

          if (bala == 0.0)
            Die(FUNCTION_NAME, "Error in current balance");
          else
            frac = sum/bala;

          /* Print progress */

          if ((it > 0) && (imax < 0))
            {
              /* Only convergence limit given */

              if ((k = (long)log10(frac)) != compl)
                {
                  /* Print */

                    fprintf(outp, "%4ld iterations completed: frac = %1.2E\n",
                            it, frac);
                  compl = k;
                }
            }
          else if (it > 0)
            {
              /* Maximum number of iterations given */

              k = (long)(100.0*((double)it + 1.0)/((double)imax));
              if ((frac < flim) || (!(k % 10) && (k != compl)))
                {
                  /* Print */

                  fprintf(outp,
                          "%4ld%% of iterations completed: frac = %1.2E\n",
                          k, frac);
                  compl = k;
                }
            }

          /* Check convergence */

          if ((it > 10) && (frac < flim))
            break;

          /*******************************************************************/
        }

      /***********************************************************************/

      /***** Calculate importances *******************************************/

      /* Loop over mesh and clear results */

      loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      while (loc0 > VALID_PTR)
        {
          /* Source importance */

          ptr = (long)RDB[loc0 + RMX_CELL_IMP_SRC];
          CheckPointer(FUNCTION_NAME, "(ptr 1)", DATA_ARRAY, ptr);

          for (n = 0; n < ng; n++)
            WDB[ptr + n] = 0.0;

          /* Current importance */

          ptr = (long)RDB[loc0 + RMX_CELL_IMP_CURR];
          CheckPointer(FUNCTION_NAME, "(ptr 2)", DATA_ARRAY, ptr);

          for (n = 0; n < ng; n++)
            WDB[ptr + n] = 0.0;

          /* Next */

          loc0 = NextItem(loc0);
        }

      /* Loop over mesh and calculate total response */

      Rtot = 0.0;

      loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      while (loc0 > VALID_PTR)
        {
          /* Get matrix size */

          nmax = (long)RDB[loc0 + RMX_CELL_MTX_SIZE];
          CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

          /* Loop over detectors */

          det = (long)RDB[loc0 + RMX_CELL_PTR_DET];
          while (det > VALID_PTR)
            {
              /* Get weighting factor */

              wgt = RDB[det + RMX_DET_WGT];
              CheckValue(FUNCTION_NAME, "wgt", "", wgt, 0.0, INFTY);

              /* Add to total response */

              ptr = (long)RDB[det + RMX_DET_MC_RES_TOT];
              CheckPointer(FUNCTION_NAME, "(ptr 3)", RES2_ARRAY, ptr);

              for (n = 0; n < ng; n++)
                Rtot = Rtot + RES2[ptr + n]*wgt;

              /* Next detector */

              det = NextItem(det);
            }

          /* Next */

          loc0 = NextItem(loc0);
        }

      /* Reset sum of contributions */

      chk = 0.0;

      /* Loop over mesh */

      loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      while (loc0 > VALID_PTR)
        {
          /* Get matrix size */

          nmax = (long)RDB[loc0 + RMX_CELL_MTX_SIZE];
          CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

          /* Get total inward net currents */

          ptr = (long)RDB[loc0 + RMX_CELL_MC_CURR_NET_IN];
          CheckPointer(FUNCTION_NAME, "ptr", RES2_ARRAY, ptr);
          Jin = &RES2[ptr];

          /* Pointer to solution */

          ptr = (long)RDB[loc0 + RMX_CELL_ADJ_SOL_IN_CURR];
          CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
          I = &WDB[ptr];

          /* Divide importance by current */

          for (n = 0; n < ng; n++)
            for (i = 0; i < nmax; i++)
              if (Jin[i*ng + n] > 0.0)
                I[i*ng + n] = I[i*ng + n]/Jin[i*ng + n]/bala*Rtot;

          /* Calculate average */

          ptr = (long)RDB[loc0 + RMX_CELL_IMP_CURR];
          CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

          for (n = 0; n < ng; n++)
            {
              /* Calculate sum of currents */

              sum = 0;
              for (i = 0; i < nmax; i++)
                sum = sum + Jin[i*ng + n];

              /* Calculate weighted average */

              if (sum > 0.0)
                for (i = 0; i < nmax; i++)
                  WDB[ptr + n] = RDB[ptr + n] + I[i*ng + n]*Jin[i*ng + n]/sum;
            }

          /* Get total source */

          ptr = (long)RDB[loc0 + RMX_CELL_MC_SRCC_TOT];
          CheckPointer(FUNCTION_NAME, "ptr", RES2_ARRAY, ptr);
          S = &RES2[ptr];

          /* Loop over energies */

          for (m = 0; m < ng; m++)
            if (S[m] > 0.0)
              {
                /* Reset direct contribution */

                dir = 0.0;

                /* Loop over detectors */

                det = (long)RDB[loc0 + RMX_CELL_PTR_DET];
                while (det > VALID_PTR)
                  {
                    /* Get weighting factor */

                    wgt = RDB[det + RMX_DET_WGT];
                    CheckValue(FUNCTION_NAME, "wgt", "", wgt, 0.0, INFTY);

                    /* Add to direct contribution */

                    ptr = (long)RDB[det + RMX_DET_COEF_DIR_SRCC_RES];
                    CheckPointer(FUNCTION_NAME, "(ptr 4)", DATA_ARRAY, ptr);
                    dir = dir + RDB[ptr + m]*wgt;

                    /* Next detector */

                    det = NextItem(det);
                  }

                /* Put direct contribution to source importance */

                ptr = (long)RDB[loc0 + RMX_CELL_IMP_SRC];
                CheckPointer(FUNCTION_NAME, "(ptr 5)", DATA_ARRAY, ptr);

                if (dir > 0.0)
                  WDB[ptr + m] = dir;
                else
                  WDB[ptr + m] = 0.0;

                /* Get source coefficient vector */

                ptr = (long)RDB[loc0 + RMX_CELL_COEF_ADJ_SRCC];
                CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
                s = &RDB[ptr];

                /* Get outward current solution */

                ptr = (long)RDB[loc0 + RMX_CELL_ADJ_SOL_OUT_CURR];
                CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
                Jout = &WDB[ptr];

                /* Add contribution to source importance */

                ptr = (long)RDB[loc0 + RMX_CELL_IMP_SRC];
                CheckPointer(FUNCTION_NAME, "(ptr 6)", DATA_ARRAY, ptr);

                for (n = 0; n < ng; n++)
                  for (i = 0; i < nmax; i++)
                    WDB[ptr + m] = RDB[ptr + m]
                      + Jout[i*ng + n]*s[i*ng*ng + m*ng + n]/S[m]/bala*Rtot;

                /* Add to sum for checking */

                chk = chk + S[m]*RDB[ptr + m];
              }

          /* Compare current importances to previous */

          ptr = (long)RDB[loc0 + RMX_CELL_IMP_CURR];
          CheckPointer(FUNCTION_NAME, "(ptr 9)", DATA_ARRAY, ptr);
          I = &WDB[ptr];

          ptr = (long)RDB[loc0 + RMX_CELL_IMP_CURR_KEEP];
          CheckPointer(FUNCTION_NAME, "(ptr 10)", DATA_ARRAY, ptr);

          for (n = 0; n < ng; n++)
            {
              /* Check zero */

              if ((I[n] == 0.0) && (RDB[ptr + n] > 0.0))
                I[n] = RDB[ptr + n];
              else
                WDB[ptr + n] = I[n];
            }

          /* Next */

          loc0 = NextItem(loc0);
        }

      /***********************************************************************/

      /***** Filter results **************************************************/

      /* Loop over mesh */

      loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      while (loc0 > VALID_PTR)
        {
          /* Skip if not GVR mode */

          if ((long)RDB[rmx + RMX_MODE] != RMX_MODE_GVR)
            break;

          /* Pointer to average importances */

          ptr = (long)RDB[loc0 + RMX_CELL_IMP_CURR];
          CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
          I = &WDB[ptr];

          /* Loop over energy groups */

          for (n = 0; n < ng; n++)
            if (I[n] > 0.0)
              {
                /* Loop over neighbours */

                loc1 = (long)RDB[loc0 + RMX_CELL_PTR_BOUNDS];
                CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc1);

                while (loc1 > VALID_PTR)
                  {
                    /* Pointer to mesh cell */

                    loc2 = (long)RDB[loc1 + RMX_CELL_BOUND_PTR_CELL];
                    CheckPointer(FUNCTION_NAME, "loc2", DATA_ARRAY, loc2);

                    /* Pointer to average importances */

                    ptr = (long)RDB[loc2 + RMX_CELL_IMP_CURR];
                    CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

                    /* Compare */

                    if ((RDB[ptr + n] > 0.0) && (RDB[ptr + n]/I[n] < 1E-6))
                      {
                        /* Reset current importances */

                        ptr = (long)RDB[loc2 + RMX_CELL_IMP_CURR];
                        CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

                        for (m = 0; m < ng; m++)
                          WDB[ptr + n] = 0.0;

                        /* Reset current importances */

                        ptr = (long)RDB[loc2 + RMX_CELL_IMP_SRC];
                        CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

                        for (m = 0; m < ng; m++)
                          WDB[ptr + n] = 0.0;

                        /* Reset partial importances */

                        ptr = (long)RDB[loc2 + RMX_CELL_ADJ_SOL_IN_CURR];
                        CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

                        for (m = 0; m < ng*nmax; m++)
                          WDB[ptr + n] = 0.0;
                      }

                    /* Next */

                    loc1 = NextItem(loc1);
                  }
              }

          /* Next */

          loc0 = NextItem(loc0);
        }

      /***********************************************************************/

      /***** Calculation completed *******************************************/

      /* Print */

      fprintf(outp, "\nCompleted in %s. ",
              TimeIntervalStr(TimerVal(TIMER_MISC)));

      /* Check number of iterations */

      if (it == imax)
        fprintf(outp,
                "Maximum number of iterations reached (no convergence).\n");
      else
        fprintf(outp, "Convergence reached after %ld iterations.\n", it);

      /* Print maximum differences */

      fprintf(outp, "Relative difference to total response: %1.2E.\n",
              chk/Rtot - 1.0);

      /* Stop timer */

      StopTimer(TIMER_MISC);

      /***********************************************************************/
    }

  /***************************************************************************/

  /* Newline */

  fprintf(outp, "\n");

  /***************************************************************************/
}

/*****************************************************************************/
