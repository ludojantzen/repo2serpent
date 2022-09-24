/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : importancesolver.c                             */
/*                                                                           */
/* Created:       2016/04/21 (JLe)                                           */
/* Last modified: 2020/03/04 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Calculates importance maps for weight windows                */
/*                                                                           */
/* Comments: - Batchien lukumäärän pitää olla jaollinen 50:llä noissa testi- */
/*             moodeissa, tai muuten detektorioutputti ei vastaa täysin      */
/*             tämän laskemaa (exit-käskyt lopettaa ajon ennen viimeistä     */
/*             printtausta).                                                 */
/*                                                                           */
/*           - Detektorien päällekkäisyys on mahdollinen ongelma (lisää      */
/*             tarkistus?)                                                   */
/*                                                                           */
/*           - Jos ton MC-datan halutaan säilyvän useamman iteraation yli,   */
/*             niin toi normeeraus pitää asettaa ykköseksi, tms.             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ImportanceSolver:"

void RefTest(long);

void PrintRMX(long);

void WriteWWMfile(long);

void PrintCurr(long);

/*****************************************************************************/

void ImportanceSolver()
{
  long rmx, msh, i, j, n, m, ptr, nmax, nr, nx, ny, nz, det, ng;
  double *J, *Js, *Jin, *Jout, *Ri, *Rs, *R, *S, chk, sum, norm, Rtot, div;
  double *F, *Fs;

  /* Check if response matrix calculation is on */

  if ((long)RDB[DATA_RMTX_CALC] == NO)
    return;

  /* Check if mfp iteration */

  if ((long)RDB[DATA_RMTX_MFP_CALC] == YES)
    return;

  /* Check corrector step */

  if (((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP) &&
      ((long)RDB[DATA_BURN_SIE] == NO))
    return;

  /* Start timer */

  StartTimer(TIMER_RMX);

  /* Reduce RES2 array (just in case) */

  ReducePrivateRes();

  /***************************************************************************/

  /***** Process data for response matrix calculation ************************/

  /* Reset number of cells with responses */

  nr = 0;

  /* Pointer to response matrix solver */

  rmx = (long)RDB[DATA_PTR_RMX0];
  while (rmx > VALID_PTR)
    {
      /* Pointer to mesh */

      if ((long)RDB[rmx + RMX_PTR_MESH] < VALID_PTR)
        {
          /* Pointer to next */

          rmx = NextItem(rmx);

          /* Cycle loop */

          continue;
        }

      /* Print */

      if ((long)RDB[rmx + RMX_MODE] == RMX_MODE_WWG)
        fprintf(outp, "Calculating importances...\n");
      else if ((long)RDB[DATA_RMX_CONVG_ACC] == NO)
        fprintf(outp, "Running response matrix solver...\n\n");

      /***********************************************************************/

      /***** Get normalization coefficient ***********************************/

      /* Get number of energy groups */

      ng = (long)RDB[rmx + RMX_NG];
      CheckValue(FUNCTION_NAME, "ng", "", ng, 1, 100000);

      /* Divider */

      if ((long)RDB[DATA_RMX_CONVG_ACC] == YES)
        div = 1.0;
      else if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
        div = RDB[DATA_CRIT_CYCLES]/RDB[DATA_BATCH_INTERVAL];
      else
        div = RDB[DATA_SRC_BATCHES]/RDB[DATA_BATCH_INTERVAL];

      /* Check value */

      CheckValue(FUNCTION_NAME, "div", "", div, ZERO, INFTY);

      /* Normalization coefficient */

      ptr = (long)RDB[RES_NORM_COEF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* NOTE: tohon matriisiin pitää liittää jotenki partikkelin */
      /* tyyppi. */

      if ((long)RDB[DATA_RMX_CONVG_ACC] == YES)
        norm = 1.0;
      else if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES)
        norm = Mean(ptr, 1)/div;
      else
        norm = Mean(ptr, 0)/div;

      /* Check value */

      CheckValue(FUNCTION_NAME, "norm", "", norm, ZERO, INFTY);

      /***********************************************************************/

      /***** WWG mode ********************************************************/

      if ((long)RDB[rmx + RMX_MODE] == RMX_MODE_WWG)
        {
          /* Loop over cells */

          msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
          while (msh > VALID_PTR)
            {
              /* Get matrix size */

              nmax = (long)RDB[msh + RMX_CELL_MTX_SIZE];
              CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

              /* Calculate source importance */

              for (n = 0; n < ng; n++)
                {
                  /* Total source */

                  ptr = (long)RDB[msh + RMX_CELL_MC_WWG_SRC0];
                  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
                  div = RES2[ptr + n];

                  /* Contributing source */

                  ptr = (long)RDB[msh + RMX_CELL_MC_WWG_SRC1];
                  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
                  chk = RES2[ptr + n];

                  /* Importance */

                  ptr = (long)RDB[msh + RMX_CELL_IMP_SRC];
                  CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

                  if (div > 0.0)
                    WDB[ptr + n] = chk/div;
                  else if (chk > 0.0)
                    Die(FUNCTION_NAME, "Error in source terms");
                }

              /* Calculate segment-wise current importance */

              for (n = 0; n < nmax*ng; n++)
                {
                  /* Total current */

                  ptr = (long)RDB[msh + RMX_CELL_MC_WWG_CURR0];
                  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
                  div = RES2[ptr + n];

                  /* Contributing current */

                  ptr = (long)RDB[msh + RMX_CELL_MC_WWG_CURR1];
                  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
                  chk = RES2[ptr + n];

                  /* Importance */

                  ptr = (long)RDB[msh + RMX_CELL_ADJ_SOL_IN_CURR];
                  CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

                  if (div > 0.0)
                    WDB[ptr + n] = chk/div;
                  else if (chk > 0.0)
                    Die(FUNCTION_NAME, "Error in current terms");
                }

              /* Calculate total current importance */

              for (n = 0; n < ng; n++)
                {
                  /* Total current */

                  ptr = (long)RDB[msh + RMX_CELL_MC_WWG_CURR0];
                  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);

                  div = 0.0;
                  for (i = 0; i < nmax; i++)
                    div = div + RES2[ptr + i*ng + n];

                  /* Contributing current */

                  ptr = (long)RDB[msh + RMX_CELL_MC_WWG_CURR1];
                  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);

                  chk = 0.0;
                  for (i = 0; i < nmax; i++)
                    chk = chk + RES2[ptr + i*ng + n];

                  /* Importance */

                  ptr = (long)RDB[msh + RMX_CELL_IMP_CURR];
                   CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

                  if (div > 0.0)
                    WDB[ptr + n] = chk/div;
                  else if (chk > 0.0)
                    Die(FUNCTION_NAME, "Error in current terms");
                }

              /* Next cell */

              msh = NextItem(msh);
            }

          /* OK */

          fprintf(outp, "OK.\n\n");

          /* Normalize importances */

          NormalizeImp(rmx);

          /* Check for iteration mode */

          if ((long)RDB[DATA_RUN_VR_ITER] == YES)
            break;
          else
            WriteWWMesh(rmx);

          /* Pointer to next */

          rmx = NextItem(rmx);

          /* Cycle loop */

          continue;
        }

      /***********************************************************************/

      /***** Calculate response in GVR mode if not provided ******************/

      /* Tää on vaihtoehtoinen yritelmä TLE:lle */

      if (1 == 2)
        {
          /* Pointer to detector */

          ptr = (long)RDB[rmx + RMX_PTR_DET];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Check if detector is not linked */

          if (((long)RDB[rmx + RMX_MODE] == RMX_MODE_GVR) &&
              ((long)RDB[ptr + RMX_DET_PTR_DET] < VALID_PTR))
            {
              /* Loop over cells */

              msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
              while (msh > VALID_PTR)
                {
                  /* Get matrix size */

                  nmax = (long)RDB[msh + RMX_CELL_MTX_SIZE];
                  CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

                  /* Get inward current and source rate */

                  ptr = (long)RDB[msh + RMX_CELL_MC_CURR_NET_IN];
                  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
                  Jin = &RES2[ptr];

                  ptr = (long)RDB[msh + RMX_CELL_MC_SRCC_TOT];
                  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
                  S = &RES2[ptr];

                  /* Get pointer to detector */

                  det = (long)RDB[msh + RMX_CELL_PTR_DET];
                  CheckPointer(FUNCTION_NAME, "(det)", DATA_ARRAY, det);

                  /* Divide by cell volume */

                  if ((div = RDB[msh + RMX_CELL_RVOL]) < ZERO)
                    Die(FUNCTION_NAME, "Error in volume");

                  /* Get result vectors */

                  ptr = (long)RDB[det + RMX_DET_MC_RES_IN_CURR];
                  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
                  Ri = &RES2[ptr];

                  ptr = (long)RDB[det + RMX_DET_MC_RES_SRCC];
                  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
                  Rs = &RES2[ptr];

                  ptr = (long)RDB[det + RMX_DET_MC_RES_TOT];
                  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
                  R = &RES2[ptr];

                  /* Loop over energy groups */

                  for (n = 0; n < ng; n++)
                    {
                      /* Add source contribution */

                      Rs[n] = Rs[n] + S[n]/div;
                      R[n] = R[n] + S[n]/div;

                      /* Add current contributions */

                      for (i = 0; i < nmax; i++)
                        {
                          Ri[i*ng + n] = Ri[i*ng + n] + Jin[i*ng + n]/div;
                          R[n] = R[n] + Jin[i*ng + n]/div;
                        }
                    }

                  /* Put number of scores */

                  ptr = (long)RDB[det + RMX_DET_MC_RES_SCORE_N];
                  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
                  RES2[ptr] = 1E6;

                  /* Next cell */

                  msh = NextItem(msh);
                }
            }
        }

      /***********************************************************************/

      /***** Remove GVR scores with bad statistics ***************************/

      /* Check mode */

      if ((long)RDB[rmx + RMX_MODE] == RMX_MODE_GVR)
        {
          /* Reset sum and count */

          sum = 0.0;
          n = 0;

          /* Loop over cells and count number of scores */

          msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
          while (msh > VALID_PTR)
            {
              /* Get pointer to detector */

              det = (long)RDB[msh + RMX_CELL_PTR_DET];
              CheckPointer(FUNCTION_NAME, "(det)", DATA_ARRAY, det);

              /* Add to sum */

              ptr = (long)RDB[det + RMX_DET_MC_RES_SCORE_N];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              sum = sum + RES2[ptr];

              /* Add to count */

              n++;

              /* Next cell */

              msh = NextItem(msh);
            }

          /* Get average */

          if (n > 0)
            chk = sum/((double)n);

          /* Loop over cells */

          msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
          while (msh > VALID_PTR)
            {
              /* Get matrix size */

              nmax = (long)RDB[msh + RMX_CELL_MTX_SIZE];
              CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

              /* Get pointer to detector */

              det = (long)RDB[msh + RMX_CELL_PTR_DET];
              CheckPointer(FUNCTION_NAME, "(det)", DATA_ARRAY, det);

              /* Number of scores */

              ptr = (long)RDB[det + RMX_DET_MC_RES_SCORE_N];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);

              /* Check (toi yllä oleva on turha jos käytetään */
              /* fiksattu lukua). */

              if (RES2[ptr] < 10)
                {
                  /* Get result vectors */

                  ptr = (long)RDB[det + RMX_DET_MC_RES_IN_CURR];
                  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
                  Ri = &RES2[ptr];

                  ptr = (long)RDB[det + RMX_DET_MC_RES_SRCC];
                  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
                  Rs = &RES2[ptr];

                  ptr = (long)RDB[det + RMX_DET_MC_RES_TOT];
                  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
                  R = &RES2[ptr];

                  /* Loop over energy groups */

                  for (n = 0; n < ng; n++)
                    {
                      /* Reset total and source contributions */

                      R[n] = 0.0;
                      Rs[n] = 0.0;

                      /* Loop over boundaries and reset current */
                      /* contributions */

                      for (i = 0; i < nmax; i++)
                        Ri[i*ng + n] = 0.0;
                    }
                }

              /* Next cell */

              msh = NextItem(msh);
            }
        }

      /***********************************************************************/

      /***** Put local results to global structures **************************/

      /* Loop over RMX detectors and reset */

      ptr = (long)RDB[rmx + RMX_PTR_DET];
      while (ptr > VALID_PTR)
        {
          /* Reset temporary array */

          WDB[ptr + RMX_DET_MC_RES_TOT] = 0.0;

          /* Next */

          ptr = NextItem(ptr);
        }

      /* Loop over mesh */

      msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      while (msh > VALID_PTR)
        {
          /* Loop over detectors */

          det = (long)RDB[msh + RMX_CELL_PTR_DET];
          while (det > VALID_PTR)
            {
              /* Initialize weights */

              if ((long)RDB[rmx + RMX_MODE] != RMX_MODE_SINGLE)
                WDB[det + RMX_DET_WGT] = 1.0;

              /* Get response */

              ptr = (long)RDB[det + RMX_DET_MC_RES_TOT];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              R = &RES2[ptr];

              /* Reset previous solution */

              WDB[det + RMX_DET_FWD_PREV] = 0.0;

              /* Pointer to RMX detector */

              ptr = (long)RDB[det + RMX_DET_PTR_DET0];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Loop over energy groups and add to total response */

              for (n = 0; n < ng; n++)
                WDB[ptr + RMX_DET_MC_RES_TOT] =
                  RDB[ptr + RMX_DET_MC_RES_TOT] + R[n];

              /* Next detector */

              det = NextItem(det);
            }

          /* Next mesh cell */

          msh = NextItem(msh);
        }

      /***********************************************************************/

      /***** Calculate sum of responses and outward currents *****************/

      /* Reset sums */

      Rtot = 0.0;
      sum = 0.0;

      /* Loop over mesh cells */

      msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      while (msh > VALID_PTR)
        {
          /* Loop over detectors */

          det = (long)RDB[msh + RMX_CELL_PTR_DET];
          while (det > VALID_PTR)
            {
              /* Get total response in cell */

              ptr = (long)RDB[det + RMX_DET_MC_RES_TOT];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              R = &RES2[ptr];

              /* Add to total */

              for (n = 0; n < ng; n++)
                Rtot = Rtot + R[n];

              /* Get matrix size */

              nmax = (long)RDB[msh + RMX_CELL_MTX_SIZE];
              CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

              /* Get outward current */

              ptr = (long)RDB[msh + RMX_CELL_MC_CURR_NET_OUT];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              Jout = &RES2[ptr];

              /* Add to sum */

              for (n = 0; n < ng*nmax; n++)
                sum = sum + Jout[n];

              /* Next detector */

              det = NextItem(det);
            }

          /* Next cell */

          msh = NextItem(msh);
        }

      /* Check */

      if (Rtot == 0.0)
        {
          /* Print error */

          fprintf(outp, "Calculation failed: no detector contributions.\n\n");

          /* pointer to next */

          rmx = NextItem(rmx);

          /* Cycle loop */

          continue;
        }
      else if (sum == 0.0)
        {
          /* Print error */

          fprintf(outp, "Calculation failed: no boundary currents.\n\n");

          /* pointer to next */

          rmx = NextItem(rmx);

          /* Cycle loop */

          continue;
        }

      /***********************************************************************/

      /***** Calculate coupling coefficients *********************************/

      /* Loop over mesh cells */

      msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      while (msh > VALID_PTR)
        {
          /* Get matrix size */

          nmax = (long)RDB[msh + RMX_CELL_MTX_SIZE];
          CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

          /*******************************************************************/

          /***** Get pointers used in the calculation ************************/

          ptr = (long)RDB[msh + RMX_CELL_MC_CURR_NET_IN];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          Jin = &RES2[ptr];

          ptr = (long)RDB[msh + RMX_CELL_MC_CURR_NET_OUT];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          Jout = &RES2[ptr];

          ptr = (long)RDB[msh + RMX_CELL_MC_CURR_THROUGH];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          J = &RES2[ptr];

          ptr = (long)RDB[msh + RMX_CELL_MC_CURR_SRCC_OUT];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          Js = &RES2[ptr];

          ptr = (long)RDB[msh + RMX_CELL_MC_SRCC_TOT];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          S = &RES2[ptr];

          /*******************************************************************/

          /***** Sanity checks ***********************************************/

          /* Check outward currents */

          for (n = 0; n < ng; n++)
            for (i = 0; i < nmax; i++)
              {
                /* Get current coming from source */

                chk = 0.0;
                for (m = 0; m < ng; m++)
                  chk = chk + Js[i*ng*ng + m*ng + n];

                /* Add currents passing through */

                for (m = 0; m < ng; m++)
                  for (j = 0; j < nmax; j++)
                    chk = chk + J[j*nmax*ng*ng + m*nmax*ng + i*ng + n];

                /* Check */

                if (Jout[i*ng + n] > 0.0)
                  if (fabs(chk/Jout[i*ng + n] - 1.0) > 1E-6)
                    Warn(FUNCTION_NAME, "Mismatch in outward current");
              }

          /* Loop over detectors */

          det = (long)RDB[msh + RMX_CELL_PTR_DET];
          while (det > VALID_PTR)
            {
              /* Get Monte Carlo results */

              ptr = (long)RDB[det + RMX_DET_MC_RES_IN_CURR];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              Ri = &RES2[ptr];

              ptr = (long)RDB[det + RMX_DET_MC_RES_SRCC];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              Rs = &RES2[ptr];

              ptr = (long)RDB[det + RMX_DET_MC_RES_TOT];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              R = &RES2[ptr];

              /* Check response */

              for (n = 0; n < ng; n++)
                if (R[n] > 0.0)
                  {
                    /* Add to count */

                    nr++;

                    /* Get direct contribution */

                    chk = Rs[n];

                    /* Add current contributions */

                    for (i = 0; i < nmax; i++)
                      chk = chk + Ri[i*ng + n];

                    /* Check */

                    if (fabs(chk/R[n] - 1.0) > 1E-6)
                      Die(FUNCTION_NAME, "Mismatch in response %E %E %E", chk,
                          R[n], chk/R[n] - 1.0);
                  }

              /* Next detector */

              det = NextItem(det);
            }

          /* Check for non-zero currents */

          for (n = 0; n < nmax*ng; n++)
            if ((Jin[n] > 0.0) || (Jout[n] > 0.0))
              break;

          /* Check index */

          if (n < nmax*ng)
            {
              /* Reset empty flag */

              WDB[msh + RMX_CELL_EMPTY] = (double)NO;
            }
          else
            {
              /* Set empty flag */

              WDB[msh + RMX_CELL_EMPTY] = (double)YES;

              /* Pointer to next mesh cell */

              msh = NextItem(msh);

              /* Cycle loop */

              continue;
            }

          /*******************************************************************/

          /***** Calculate coefficient vectors and matrices ******************/

          /* Forward source coefficient */

          if (((long)RDB[DATA_RMX_TEST_MODE] != NO) ||
              ((long)RDB[DATA_RMX_CONVG_ACC] == YES))
            {
              ptr = (long)RDB[msh + RMX_CELL_COEF_FWD_SRCC];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Loop over boundaries */

              for (n = 0; n < ng; n++)
                for (i = 0; i < nmax; i++)
                  for (m = 0; m < ng; m++)
                    if (S[m] > 0.0)
                      WDB[ptr + i*ng*ng + m*ng + n] =
                        Js[i*ng*ng + m*ng + n]/S[m];
            }

          /* Adjoint source coefficient */

          ptr = (long)RDB[msh + RMX_CELL_COEF_ADJ_SRCC];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Loop over boundaries */

          for (n = 0; n < ng; n++)
            for (i = 0; i < nmax; i++)
              for (m = 0; m < ng; m++)
                if (Jout[i*ng + n] > 0.0)
                  WDB[ptr + i*ng*ng + m*ng + n] =
                    Js[i*ng*ng + m*ng + n]/Jout[i*ng + n];

          /* Loop over detectors */

          det = (long)RDB[msh + RMX_CELL_PTR_DET];
          while (det > VALID_PTR)
            {
              /* Get Monte Carlo results */

              ptr = (long)RDB[det + RMX_DET_MC_RES_IN_CURR];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              Ri = &RES2[ptr];

              ptr = (long)RDB[det + RMX_DET_MC_RES_SRCC];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              Rs = &RES2[ptr];

              ptr = (long)RDB[det + RMX_DET_MC_RES_TOT];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              R = &RES2[ptr];

              /* Forward response coefficient */

              ptr = (long)RDB[det + RMX_DET_COEF_FWD_RES];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Loop over boundaries */

              for (n = 0; n < ng; n++)
                for (i = 0; i < nmax; i++)
                  if (Jin[i*ng + n] > 0.0)
                    WDB[ptr + i*ng + n] = Ri[i*ng + n]/Jin[i*ng + n];

              /* Adjoint response coefficient */

              ptr = (long)RDB[det + RMX_DET_COEF_ADJ_RES];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Loop over boundaries */

              for (n = 0; n < ng; n++)
                for (i = 0; i < nmax; i++)
                  if (Ri[i*ng + n] > 0.0)
                    WDB[ptr + i*ng + n] = (R[n] - Rs[n])/Ri[i*ng + n]/Rtot;

              /* Direct source-to-response coefficient */

              ptr = (long)RDB[det + RMX_DET_COEF_DIR_SRCC_RES];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              for (n = 0; n < ng; n++)
                {
                  /* Check */

                  if (S[n] > 0.0)
                    WDB[ptr + n] = Rs[n]/S[n];
                  else if (Rs[n] != 0.0)
                    Die(FUNCTION_NAME, "Error in direct coefficient");
                }

              /* Next detector */

              det = NextItem(det);
            }

          /* Forward current transfer matrix */

          if (((long)RDB[DATA_RMX_TEST_MODE] != NO) ||
              ((long)RDB[DATA_RMX_CONVG_ACC] == YES))
            {
              /* Get pointer */

              ptr = (long)RDB[msh + RMX_CELL_COEF_FWD_TRANS_MTX];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Loop over boundaries */

              for (m = 0; m < ng; m++)
                for (j = 0; j < nmax; j++)
                  for (n = 0; n < ng; n++)
                    for (i = 0; i < nmax; i++)
                      if (Jin[i*ng + n] > 0.0)
                        WDB[ptr + i*nmax*ng*ng + n*nmax*ng + j*ng + m] =
                          J[i*nmax*ng*ng + n*nmax*ng + j*ng + m]/
                          Jin[i*ng + n];
            }

          /* Adjoint current transfer matrix */

          ptr = (long)RDB[msh + RMX_CELL_COEF_ADJ_TRANS_MTX];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Loop over boundaries */

          for (m = 0; m < ng; m++)
            for (j = 0; j < nmax; j++)
              if (Jout[j*ng + m] > 0.0)
                for (n = 0; n < ng; n++)
                  for (i = 0; i < nmax; i++)
                    WDB[ptr + i*nmax*ng*ng + n*nmax*ng + j*ng + m] =
                      J[i*nmax*ng*ng + n*nmax*ng + j*ng + m]/Jout[j*ng + m];

          /* Loop over outward boundaries and check currents */

          for (m = 0; m < ng; m++)
            for (j = 0; j < nmax; j++)
              {
                /* Calculate sum of matrix coefficients */

                chk = 0.0;
                if (Jout[j*ng + m] > 0.0)
                  for (n = 0; n < ng; n++)
                    for (i = 0; i < nmax; i++)
                      chk = chk +
                        J[i*nmax*ng*ng + n*nmax*ng + j*ng + m]/Jout[j*ng + m];

                /* Calculate total current from source to face j and group m */

                sum = 0.0;
                for (n = 0; n < ng; n++)
                  sum = sum + Js[j*ng*ng + n*ng + m];

                /* Check current from source to boundary */

                if (sum > 0.0)
                  {
                    /* Check sum */

                    if (fabs(1.0 - chk) < ZERO)
                      Warn(FUNCTION_NAME, "Error in sum (S > 0) %E", 1.0 - chk);
                  }
                else if ((chk > 0.0) && (fabs(1.0 - chk) > 1E-6))
                  Warn(FUNCTION_NAME, "Error in sum (S = 0) %E", chk);
              }

          /*******************************************************************/

          /***** Normalize sources and responses *****************************/

          /* Normalize total source */

          ptr = (long)RDB[msh + RMX_CELL_MC_SRCC_TOT];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);

          for (n = 0; n < ng; n++)
            RES2[ptr + n] = S[n]*norm;

          /* Normalize net currents */

          ptr = (long)RDB[msh + RMX_CELL_MC_CURR_NET_IN];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);

          for (n = 0; n < ng; n++)
            for (i = 0; i < nmax; i++)
              RES2[ptr + i*ng + n] = Jin[i*ng + n]*norm;

          ptr = (long)RDB[msh + RMX_CELL_MC_CURR_NET_OUT];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);

          for (n = 0; n < ng; n++)
            for (i = 0; i < nmax; i++)
              RES2[ptr + i*ng + n] = Jout[i*ng + n]*norm;

          ptr = (long)RDB[msh + RMX_CELL_MC_CURR_SRCC_OUT];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);

          for (n = 0; n < ng; n++)
            for (i = 0; i < nmax; i++)
              for (m = 0; m < ng; m++)
                RES2[ptr + i*ng*ng + m*ng + n] = Js[i*ng*ng + m*ng + n]*norm;

          ptr = (long)RDB[msh + RMX_CELL_MC_CURR_THROUGH];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);

          for (m = 0; m < ng; m++)
            for (j = 0; j < nmax; j++)
              for (n = 0; n < ng; n++)
                for (i = 0; i < nmax; i++)
                  RES2[ptr + i*nmax*ng*ng + n*nmax*ng + j*ng + m] =
                    J[i*nmax*ng*ng + n*nmax*ng + j*ng + m]*norm;

          /* Loop over detectors */

          det = (long)RDB[msh + RMX_CELL_PTR_DET];
          while (det > VALID_PTR)
            {
              /* Normalize total response */

              ptr = (long)RDB[det + RMX_DET_MC_RES_TOT];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);

              /* Normalize to physical results */

              for (n = 0; n < ng; n++)
                RES2[ptr + n] = RES2[ptr + n]*norm;

              /* Response current */

              ptr = (long)RDB[det + RMX_DET_MC_RES_IN_CURR];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);

              for (n = 0; n < ng; n++)
                for (i = 0; i < nmax; i++)
                  RES2[ptr + i*ng + n] = RES2[ptr + i*ng + n]*norm;

              /* Next detector */

              det = NextItem(det);
            }

          /* Initialize importances to 1 (this is needed for initial   */
          /* weight factors in MULTI and GVR modes). The solutions are */
          /* set to zero before running the adjoint solver. */

          ptr = (long)RDB[msh + RMX_CELL_ADJ_SOL_IN_CURR];
          CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

          for (n = 0; n < nmax*ng; n++)
            WDB[ptr + n] = 1.0;

          ptr = (long)RDB[msh + RMX_CELL_IMP_SRC];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          for (n = 0; n < ng; n++)
            WDB[ptr + n] = 1.0;

          /*******************************************************************/

          /***** Local form factors ******************************************/

          /* Check source convergence acceleration */

          if ((long)RDB[DATA_RMX_CONVG_ACC] == YES)
            {
              /* Check number of energy groups */

              if (ng > 1)
                Die(FUNCTION_NAME, "Multiple energy groups");

              /* Local mesh size */

              nx = (long)RDB[rmx + RMX_CONVG_SUBMESH_NX];
              ny = (long)RDB[rmx + RMX_CONVG_SUBMESH_NY];
              nz = (long)RDB[rmx + RMX_CONVG_SUBMESH_NZ];

              /* Pointer to detector */

              det = (long)RDB[msh + RMX_CELL_PTR_DET];
              CheckPointer(FUNCTION_NAME, "(det)", DATA_ARRAY, det);

              /* Get response coefficients for checking */

              ptr = (long)RDB[det + RMX_DET_COEF_DIR_SRCC_RES];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              Rs = &WDB[ptr];

              ptr = (long)RDB[det + RMX_DET_COEF_FWD_RES];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              Ri = &WDB[ptr];

              /* Pointer to results */

              ptr = (long)RDB[msh + RMX_CELL_MC_SRC_FF];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              F = &RES2[ptr];

              ptr = (long)RDB[msh + RMX_CELL_MC_SURF_FF];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              Fs = &RES2[ptr];

              /* Source form factors */

              if (*S > 0.0)
                {
                  /* Pointer to source coefficients */

                  ptr = (long)RDB[msh + RMX_CELL_COEF_SRC_FF];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                  /* Reset sum */

                  sum = 0.0;

                  /* Loop over sub-regions */

                  for (n = 0; n < nx*ny*nz; n++)
                    {
                      WDB[ptr + n] = F[n]/(*S);
                      sum = sum + RDB[ptr + n];
                    }

                  /* Check sum */

                  if (*Rs > 0.0)
                    if (fabs(sum/(*Rs) - 1.0) > 1E-6)
                      Die(FUNCTION_NAME, "Error in source coefficients");
                }

              /* Pointer to surface coefficients */

              ptr = (long)RDB[msh + RMX_CELL_COEF_SURF_FF];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Loop over surfaces */

              for (i = 0; i < nmax; i++)
                if (Jin[i] > 0.0)
                  {
                    /* Reset sum */

                    sum = 0.0;

                    /* Loop over sub-regions */

                    for (n = 0; n < nx*ny*nz; n++)
                      {
                        WDB[ptr + i*nx*ny*nz + n] = Fs[i*nx*ny*nz + n]/Jin[i];
                        sum = sum + RDB[ptr + i*nx*ny*nz + n];
                      }

                    /* Check sum */

                    if (Ri[i] > 0.0)
                      if (fabs(sum/(Ri[i]) - 1.0) > 1E-6)
                        Die(FUNCTION_NAME, "Error in source coefficients");
                  }
            }

          /*******************************************************************/

          /* Next */

          msh = NextItem(msh);
        }

      /* Check response count */

      if (nr == 0)
        Die(FUNCTION_NAME, "No detector cells with contributions");

      /***********************************************************************/

      /***** Special modes for debugging, etc. *******************************/

      /* This is for calculating reference solutions */

      if ((long)RDB[DATA_RMX_TEST_MODE] == 1)
        RefTest(rmx);

      /* This is for printing the solutions */

      else if ((long)RDB[DATA_RMX_TEST_MODE] == 2)
        PrintRMX(rmx);

      /* This is for printing the currents */

      else if (1 == 2)
        PrintCurr(rmx);

      /***********************************************************************/

      /***** Run response matrix solver **************************************/

      /* Check mode */

      if ((long)RDB[DATA_RMX_CONVG_ACC] == YES)
        {
          /* Set solution flags */

          WDB[DATA_RMTX_SOLVE_FORWARD] = NO;
          WDB[DATA_RMTX_SOLVE_ADJOINT] = YES;

          /* Criticality source mode solver */

          SolveRMXCrit(rmx);

          /* Stop timer */

          StopTimer(TIMER_RMX);

          /* Exit subroutine */

          return;
        }
      else if ((long)RDB[rmx + RMX_MODE] == RMX_MODE_SINGLE)
        {
          /* Set solution flags */

          WDB[DATA_RMTX_SOLVE_FORWARD] = NO;
          WDB[DATA_RMTX_SOLVE_ADJOINT] = YES;

          /* External source solver */

          SolveRMX(rmx, NO);
        }
      else
        {
          /* Iterations for multi-response and GVR mode */

          fprintf(outp, "Begin outer iterations...\n\n");

          /* Remember original convergence criteria */

          WDB[DATA_RMX_CONV0] = RDB[rmx + RMX_CONV];
          WDB[DATA_RMX_N_ITER0] = RDB[rmx + RMX_N_ITER];

          /* Iteration loops */

          for (n = 0; n < (long)RDB[DATA_RMX_MULT_MAX_ITER] + 1; n++)
            {
              /* Update weights */

              chk = UpdateRMXWgt(rmx, n);

              /* Check convergence */

              if (n > 0)
                {
                  fprintf(outp,
                          "Convergence test for cycle %ld: %4.1f%% pass\n\n",
                          n, chk*100.0);

                  /* Check */

                  if (chk > RDB[DATA_RMX_MULT_CONVG_PASS])
                    break;
                  else if (n ==(long)RDB[DATA_RMX_MULT_MAX_ITER])
                    break;
                }

              /* Set solution flags */

              WDB[DATA_RMTX_SOLVE_FORWARD] = NO;
              WDB[DATA_RMTX_SOLVE_ADJOINT] = YES;

              /* Solve adjoint */

              SolveRMX(rmx, NO);
            }

          /* Iterations complete */

          if (n < (long)RDB[DATA_RMX_MULT_MAX_ITER])
              fprintf(outp, "Outer iterations completed.\n\n");
          else
            fprintf(outp, "Maximum number of outer iterations reached.\n\n");
        }

      /***********************************************************************/

      /* Normalize importances */

      NormalizeImp(rmx);

      /* Check for iteration mode */

      if ((long)RDB[DATA_RUN_VR_ITER] == YES)
        break;
      else
        WriteWWMesh(rmx);

      /* Check if written in ASCII file */

      if ((long)RDB[DATA_RMX_TEST_MODE] == 3)
        WriteWWMfile(rmx);

      /* Next */

      rmx = NextItem(rmx);
    }

  /* Stop timer */

  StopTimer(TIMER_RMX);
}

/*****************************************************************************/

/***** Function for calculating reference results ****************************/

void RefTest(long rmx)
{
  long msh, i, j, ptr, det, nmax, ng, n;
  double sum, *Jin, *S, div;
  char outfile[MAX_STR];
  FILE *fp;

  /* Open file */

  sprintf(outfile, "%s_ref.m", GetText(DATA_PTR_INPUT_FNAME));
  fp = fopen(outfile, "w");

  /* Switch forward mode off and adjoint mode on */

  WDB[DATA_RMTX_SOLVE_FORWARD] = NO;
  WDB[DATA_RMTX_SOLVE_ADJOINT] = YES;

  /* Solve response matrix */

  SolveRMX(rmx, NO);

  /* Allow memory allocation */

  Mem(MEM_ALLOW);

  /* Get number of energy groups */

  ng = (long)RDB[rmx + RMX_NG];
  CheckValue(FUNCTION_NAME, "ng", "", ng, 1, 100000);

  fprintf(fp, "ng = %ld;\n", ng);

  /* Loop over mesh and print source importances  */

  fprintf(fp, "src_imp = [\n");

  msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (msh > VALID_PTR)
    {
      /* Pointer to data */

      ptr = (long)RDB[msh + RMX_CELL_IMP_SRC];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Print */

      for (n = 0; n < ng; n++)
        fprintf(fp, "%1.5E ", RDB[ptr + n]);

      fprintf(fp, "\n");

      /* Pointer to next */

      msh = NextItem(msh);
    }

  fprintf(fp, "];\n");

  /* Loop over mesh and print importances of inward currents */

  fprintf(fp, "curr_imp = [\n");

  msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (msh > VALID_PTR)
    {
      ptr = (long)RDB[msh + RMX_CELL_IMP_CURR];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Print */

      for (n = 0; n < ng; n++)
        fprintf(fp, "%1.5E ", RDB[ptr + n]);

      fprintf(fp, "\n");

      /* Pointer to next */

      msh = NextItem(msh);
    }

  fprintf(fp, "];\n");

  /* Switch forward mode on and adjoint mode off */

  WDB[DATA_RMTX_SOLVE_FORWARD] = YES;
  WDB[DATA_RMTX_SOLVE_ADJOINT] = NO;

  /* Loop over mesh and reset all sources  */

  msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (msh > VALID_PTR)
    {
      /* Pointer to source vector */

      ptr = (long)RDB[msh + RMX_CELL_MC_SRCC_TOT];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      S = &RES2[ptr];

      /* Get matrix size */

      nmax = (long)RDB[msh + RMX_CELL_MTX_SIZE];
      CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

      /* Allocate memory for storing source vectors */

      ptr = ReallocMem(DATA_ARRAY, ng);
      WDB[msh + RMX_CELL_TMP_SRC] = (double)ptr;

      /* Copy values and reset */

      for (n = 0; n < ng; n++)
        {
          WDB[ptr + n] = S[n];
          S[n] = 0.0;
        }

      /* Pointer to next */

      msh = NextItem(msh);
    }

  /* Loop over energy groups */

  for (n = 0; n < ng; n++)
    {
      /* Source importances */

      fprintf(fp, "ref_src_imp(:,%ld) = [\n", n + 1);

      /* Loop over mesh */

      msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      while (msh > VALID_PTR)
        {
          /* Pointer to source vector */

          ptr = (long)RDB[msh + RMX_CELL_TMP_SRC];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          S = &WDB[ptr];

          /* Check source size */

          if (S[n] == 0.0)
            {
              /* No source */

              fprintf(fp, "%1.5E\n", 0.0);
            }
          else
            {
              /* Put unit source */

              ptr = (long)RDB[msh + RMX_CELL_MC_SRCC_TOT];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              RES2[ptr + n] = 1.0;

              /* Solve response matrix */

              SolveRMX(rmx, NO);

              /* Reset sum */

              sum = 0.0;

              /* Loop over mesh */

              ptr = (long)RDB[rmx + RMX_PTR_MESH_DATA];
              while (ptr > VALID_PTR)
                {
                  /* Loop over detectors */

                  det = (long)RDB[ptr + RMX_CELL_PTR_DET];
                  while (det > VALID_PTR)
                    {
                      /* Add calculated response to sum */

                      sum = sum + RDB[det + RMX_DET_FWD_CHK];

                      /* Next detector */

                      det = NextItem(det);
                    }

                  /* Next mesh cell */

                  ptr = NextItem(ptr);
                }

              /* Print */

              fprintf(fp, "%1.5E\n", sum);

              /* Reset source */

              ptr = (long)RDB[msh + RMX_CELL_MC_SRCC_TOT];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              RES2[ptr + n] = 0.0;
            }

          /* Pointer to next */

          msh = NextItem(msh);
        }

      fprintf(fp, "];\n");
    }

  /* Loop over mesh and reset all sources and inward currents  */

  msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (msh > VALID_PTR)
    {
      /* Reset source term */

      ptr = (long)RDB[msh + RMX_CELL_MC_SRCC_TOT];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      memset(&RES2[ptr], 0.0, ng*sizeof(double));

      /* Get matrix size */

      nmax = (long)RDB[msh + RMX_CELL_MTX_SIZE];
      CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

      /* Reset inward current vector */

      ptr = (long)RDB[msh + RMX_CELL_WRK_IN_CURR];
      CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
      memset(&WDB[ptr], 0.0, nmax*ng*sizeof(double));

      /* Reset outward current vector */

      ptr = (long)RDB[msh + RMX_CELL_WRK_OUT_CURR];
      CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
      memset(&WDB[ptr], 0.0, nmax*ng*sizeof(double));

      /* Pointer to next */

      msh = NextItem(msh);
    }

  /* Loop over energy groups */

  for (n = 0; n < ng; n++)
    {
      /* Reference current importances */

      fprintf(fp, "ref_curr_imp(:,%ld) = [\n", n + 1);

      /* Loop over mesh */

      msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      while (msh > VALID_PTR)
        {
          /* Get matrix size */

          nmax = (long)RDB[msh + RMX_CELL_MTX_SIZE];
          CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

          /* Get total inward net currents */

          ptr = (long)RDB[msh + RMX_CELL_MC_CURR_NET_IN];
          CheckPointer(FUNCTION_NAME, "ptr", RES2_ARRAY, ptr);
          Jin = &RES2[ptr];

          /* Calculate divider */

          div = 0.0;
          for (j = 0; j < nmax; j++)
            div = div + Jin[j*ng + n];

          /* Reset sum and divider */

          sum = 0.0;

          /* Loop over boundaries */

          for (i = 0; i < nmax; i++)
            {
              /* Put unit source on boundary */

              ptr = (long)RDB[msh + RMX_CELL_WRK_IN_CURR];
              CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
              WDB[ptr + i*ng + n] = 1.0;

              /* Solve response matrix */

              SolveRMX(rmx, YES);

              /* Reset source */

              WDB[ptr + i*ng + n] = 0.0;

              /* Loop over mesh */

              ptr = (long)RDB[rmx + RMX_PTR_MESH_DATA];
              while (ptr > VALID_PTR)
                {
                  /* Loop over detectors */

                  det = (long)RDB[ptr + RMX_CELL_PTR_DET];
                  while (det > VALID_PTR)
                    {
                      /* Add calculated response to sum */

                      if (div > 0.0)
                        sum = sum + RDB[det + RMX_DET_FWD_CHK]
                          *Jin[i*ng + n]/div;

                      /* Next detector */

                      det = NextItem(det);
                    }

                  /* Next mesh cell */

                  ptr = NextItem(ptr);
                }
            }

          /* Print */

          fprintf(fp, "%1.5E\n", sum);

          /* Pointer to next */

          msh = NextItem(msh);
        }

      fprintf(fp, "];\n");
    }

  /* Close file */

  fclose(fp);

  /* Exit */

  exit(0);
}

/*****************************************************************************/

/***** Function for printing results in file *********************************/

void PrintRMX(long rmx)
{
  long msh, nmax, ptr, det, i, n, ng;
  double R, *S;
  char outfile[MAX_STR];
  FILE *fp;

  /* Switch both modes on */

  WDB[DATA_RMTX_SOLVE_FORWARD] = YES;
  WDB[DATA_RMTX_SOLVE_ADJOINT] = YES;

  /* Solve response matrix */

  SolveRMX(rmx, NO);

  /* Get number of energy groups */

  ng = (long)RDB[rmx + RMX_NG];
  CheckValue(FUNCTION_NAME, "ng", "", ng, 1, 100000);

  /* Open file */

  sprintf(outfile, "%s_rmx.m", GetText(DATA_PTR_INPUT_FNAME));
  fp = fopen(outfile, "w");

  /* Print number of energy groups */

  fprintf(fp, "ng = %ld;\n\n", ng);

  fprintf(fp, "Jin = [\n");

  /* Loop over mesh */

  msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (msh > VALID_PTR)
    {
      /* Get number of boundaries */

      nmax = (long)RDB[msh + RMX_CELL_MTX_SIZE];
      CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

      /* Pointer to solution */

      ptr = (long)RDB[msh + RMX_CELL_FWD_SOL_IN_CURR];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Loop over energy groups */

      for (i = 0; i < nmax; i++)
        {
          /* Print current */

          for (n = 0; n < ng; n++)
            fprintf(fp, "%1.5E ", RDB[ptr + i*ng + n]);

          /* Newline */

          fprintf(fp, "\n");
        }

      /* Next */

      msh = NextItem(msh);
    }

  fprintf(fp, "];\n\n");

  fprintf(fp, "src = [\n");

  /* Loop over mesh */

  msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (msh > VALID_PTR)
    {
      /* Get number of boundaries */

      nmax = (long)RDB[msh + RMX_CELL_MTX_SIZE];
      CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

      /* Get source */

      ptr = (long)RDB[msh + RMX_CELL_MC_SRCC_TOT];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      S = &RES2[ptr];

      /* Print */

      for (n = 0; n < ng; n++)
        fprintf(fp, "%1.5E ", S[n]);

      /* Newline */

      fprintf(fp, "\n");

      /* Next */

      msh = NextItem(msh);
    }

  fprintf(fp, "];\n");

  fprintf(fp, "fwd = [\n");

  /* Loop over mesh */

  msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (msh > VALID_PTR)
    {
      /* Loop over detectors */

      det = (long)RDB[msh + RMX_CELL_PTR_DET];
      while (det > VALID_PTR)
        {
          /* Get response */

          ptr = (long)RDB[det + RMX_DET_MC_RES_TOT];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);

          R = 0.0;
          for (n = 0; n < ng; n++)
            R = R + RES2[ptr + n];

          /* Print */

          fprintf(fp, " %1.5E %1.5E", RDB[det + RMX_DET_FWD_CHK], R);

          /* Next detector */

          det = NextItem(det);
        }

      /* Newline */

      fprintf(fp, "\n");

      /* Next */

      msh = NextItem(msh);
    }

  fprintf(fp, "];\n\n");

  fprintf(fp, "imp_src = [\n");

  /* Loop over mesh */

  msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (msh > VALID_PTR)
    {
      /* Source importance */

      ptr = (long)RDB[msh + RMX_CELL_IMP_SRC];
      CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

      /* Print */

      for (n = 0; n < ng; n++)
        fprintf(fp, "%1.5E ", RDB[ptr + n]);

      fprintf(fp, "\n");

      /* Next */

      msh = NextItem(msh);
    }

  fprintf(fp, "];\n");

  fprintf(fp, "imp_curr = [\n");

  /* Loop over mesh */

  msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (msh > VALID_PTR)
    {
      /* Source importance */

      ptr = (long)RDB[msh + RMX_CELL_IMP_CURR];
      CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

      /* Print */

      for (n = 0; n < ng; n++)
        fprintf(fp, "%1.5E ", RDB[ptr + n]);

      fprintf(fp, "\n");

      /* Next */

      msh = NextItem(msh);
    }

  fprintf(fp, "];\n");

  /* Close file */

  fclose(fp);

  /* Exit */

  exit(0);
}

/*****************************************************************************/

/***** Function for printing results in Matlab-format file *******************/

void WriteWWMfile(long rmx)
{
  long msh, ptr, n, ng;
  char outfile[MAX_STR];
  FILE *fp;

  /* Get number of energy groups */

  ng = (long)RDB[rmx + RMX_NG];
  CheckValue(FUNCTION_NAME, "ng", "", ng, 1, 100000);

  /* Open file */

  sprintf(outfile, "%s_wwd%ld.m", GetText(DATA_PTR_INPUT_FNAME),
          (long)RDB[DATA_BURN_STEP]);
  fp = fopen(outfile, "w");

  /* Print number of energy groups */

  fprintf(fp, "ng = %ld;\n\n", ng);

  fprintf(fp, "imp_src = [\n");

  /* Loop over mesh */

  msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (msh > VALID_PTR)
    {
      /* Source importance */

      ptr = (long)RDB[msh + RMX_CELL_IMP_SRC];
      CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

      /* Print */

      for (n = 0; n < ng; n++)
        fprintf(fp, "%1.5E ", RDB[ptr + n]);

      fprintf(fp, "\n");

      /* Next */

      msh = NextItem(msh);
    }

  fprintf(fp, "];\n");

  fprintf(fp, "imp_curr = [\n");

  /* Loop over mesh */

  msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (msh > VALID_PTR)
    {
      /* Source importance */

      ptr = (long)RDB[msh + RMX_CELL_IMP_CURR];
      CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

      /* Print */

      for (n = 0; n < ng; n++)
        fprintf(fp, "%1.5E ", RDB[ptr + n]);

      fprintf(fp, "\n");

      /* Next */

      msh = NextItem(msh);
    }

  fprintf(fp, "];\n");

  /* Close file */

  fclose(fp);
}

/*****************************************************************************/

/***** Function for printing currents in file ********************************/

void PrintCurr(long rmx)
{
  long msh, nmax, ptr, n, ng;
  char outfile[MAX_STR];
  FILE *fp;

  /* Get number of energy groups */

  ng = (long)RDB[rmx + RMX_NG];
  CheckValue(FUNCTION_NAME, "ng", "", ng, 1, 100000);

  /* Open file */

  sprintf(outfile, "%s_curr.txt", GetText(DATA_PTR_INPUT_FNAME));
  fp = fopen(outfile, "w");

  /* Print number of energy groups */

  fprintf(fp, "ng: %ld\n", ng);

  /* Print inward currents */

  fprintf(fp, "\nJin:\n");

  /* Loop over mesh */

  msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (msh > VALID_PTR)
    {
      /* Get number of boundaries */

      nmax = (long)RDB[msh + RMX_CELL_MTX_SIZE];
      CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

      /* Get inward current */

      ptr = (long)RDB[msh + RMX_CELL_MC_CURR_NET_IN];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);

      /* Print vector */

      for (n = 0; n < ng*nmax; n++)
        fprintf(fp, "%1.5E ", RES2[ptr + n]);

      /* Newline */

      fprintf(fp, "\n");

      /* Next */

      msh = NextItem(msh);
    }

  /* Print outward currents */

  fprintf(fp, "\nJout:\n");

  /* Loop over mesh */

  msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (msh > VALID_PTR)
    {
      /* Get number of boundaries */

      nmax = (long)RDB[msh + RMX_CELL_MTX_SIZE];
      CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

      /* Get inward current */

      ptr = (long)RDB[msh + RMX_CELL_MC_CURR_NET_OUT];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);

      /* Print vector */

      for (n = 0; n < ng*nmax; n++)
        fprintf(fp, "%1.5E ", RES2[ptr + n]);

      /* Newline */

      fprintf(fp, "\n");

      /* Next */

      msh = NextItem(msh);
    }

  /* Close file */

  fclose(fp);

  /* Exit */

  exit(0);
}

/*****************************************************************************/
