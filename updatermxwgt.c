/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : updatermxwgt.c                                 */
/*                                                                           */
/* Created:       2018/09/11 (JLe)                                           */
/* Last modified: 2018/10/29 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Updates response weights with importance-weighted forward    */
/*              solutions.                                                   */
/*                                                                           */
/* Comments: - Used for RMX_MULTI and RMX_GVR modes.                         */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "UpdateRMXWgt:"

/*****************************************************************************/

double UpdateRMXWgt(long rmx, long cyc)
{
  long loc0, det, ng, nmax, ptr, n, mode;
  double *S, *Jin, sum, wgt, max, chk;
  const double *Is, *Iin, *rs, *r;
  char fname[MAX_STR];
  FILE *fp;

  /* Get mode */

  mode = (long)RDB[rmx + RMX_MODE];

  /* Get number of energy groups */

  ng = (long)RDB[rmx + RMX_NG];
  CheckValue(FUNCTION_NAME, "ng", "", ng, 1, 100000);

  /***************************************************************************/

  /***** Calculate importance-weighted results *******************************/

  /* Loop over mesh */

  loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (loc0 > VALID_PTR)
    {
      /* Get source */

      ptr = (long)RDB[loc0 + RMX_CELL_MC_SRCC_TOT];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      S = &RES2[ptr];

      /* Get inward current */

      ptr = (long)RDB[loc0 + RMX_CELL_MC_CURR_NET_IN];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      Jin = &RES2[ptr];

      /* Get matrix size */

      nmax = (long)RDB[loc0 + RMX_CELL_MTX_SIZE];
      CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

      /* Get source importance */

      ptr = (long)RDB[loc0 + RMX_CELL_IMP_SRC];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      Is = &RDB[ptr];

      /* Get inward current importance */

      ptr = (long)RDB[loc0 + RMX_CELL_ADJ_SOL_IN_CURR];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      Iin = &RDB[ptr];

      /* Loop over detectors */

      det = (long)RDB[loc0 + RMX_CELL_PTR_DET];
      while (det > VALID_PTR)
        {
          /* Get direct source-to-response coefficient */

          ptr = (long)RDB[det + RMX_DET_COEF_DIR_SRCC_RES];
          CheckPointer(FUNCTION_NAME, "(ptr5)", DATA_ARRAY, ptr);
          rs = &RDB[ptr];

          /* Get response vector */

          ptr = (long)RDB[det + RMX_DET_COEF_FWD_RES];
          CheckPointer(FUNCTION_NAME, "(ptr5)", DATA_ARRAY, ptr);
          r = &RDB[ptr];

          /* Put direct contribution from source */

          WDB[det + RMX_DET_FWD_CHK] = 0.0;
          for (n = 0; n < ng; n++)
            WDB[det + RMX_DET_FWD_CHK] =
              RDB[det + RMX_DET_FWD_CHK] + rs[n]*S[n]*Is[n];

          /* Add contribution from inward currents */

          for (n = 0; n < ng*nmax; n++)
            WDB[det + RMX_DET_FWD_CHK] = RDB[det + RMX_DET_FWD_CHK]
              + r[n]*Jin[n]*Iin[n];

          /* Next detector */

          det = NextItem(det);
        }

      /* Next mesh cell */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Print weights (for debugging) ***************************************/

  if (1 == 2)
    {
      /* Open file */

      sprintf(fname, "%s_wgt.m", GetText(DATA_PTR_INPUT_FNAME));

      if ((cyc == 0) && ((long)RDB[DATA_VR_ITER_IDX] == 0))
        fp = fopen(fname, "w");
      else
        fp = fopen(fname, "a");

      fprintf(fp, "w%ld(:,%ld) = [\n", (long)RDB[DATA_VR_ITER_IDX] + 1,
              cyc + 1);

      /* Loop over mesh */

      loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      while (loc0 > VALID_PTR)
        {
          /* Loop over detectors */

          det = (long)RDB[loc0 + RMX_CELL_PTR_DET];
          while (det > VALID_PTR)
            {
              /* Print weight */

              fprintf(fp, "%E\n", RDB[det + RMX_DET_WGT]);

              /* Next detector */

              det = NextItem(det);
            }

          /* Next mesh cell */

          loc0 = NextItem(loc0);
        }

      fprintf(fp, "];\n");

      fprintf(fp, "f%ld(:,%ld) = [\n", (long)RDB[DATA_VR_ITER_IDX] + 1,
              cyc + 1);

      /* Loop over mesh */

      loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      while (loc0 > VALID_PTR)
        {
          /* Loop over detectors */

          det = (long)RDB[loc0 + RMX_CELL_PTR_DET];
          while (det > VALID_PTR)
            {
              /* Print weight */

              fprintf(fp, "%E\n", RDB[det + RMX_DET_FWD_CHK]);

              /* Next detector */

              det = NextItem(det);
            }

          /* Next mesh cell */

          loc0 = NextItem(loc0);
        }

      fprintf(fp, "];\n");

      fclose(fp);
    }

  /***************************************************************************/

  /***** Update weights ******************************************************/

  /* Loop over detectors and reset */

  ptr = (long)RDB[rmx + RMX_PTR_DET];
  while (ptr > VALID_PTR)
    {
      /* Reset temporary variable */

      WDB[ptr + RMX_DET_MC_RES_TOT] = 0.0;

      /* Next */

      ptr = NextItem(ptr);
    }

  /* Loop over mesh and get total responses */

  loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (loc0 > VALID_PTR)
    {
      /* Loop over detectors */

      det = (long)RDB[loc0 + RMX_CELL_PTR_DET];
      while (det > VALID_PTR)
        {
          /* Pointer to RMX detector */

          ptr = (long)RDB[det + RMX_DET_PTR_DET0];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Add to total */

          WDB[ptr + RMX_DET_MC_RES_TOT] =
            RDB[ptr + RMX_DET_MC_RES_TOT] + RDB[det + RMX_DET_FWD_CHK];

          /* Next */

          det = NextItem(det);
        }

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Reset sum and number of non-zero values */

  sum = 0.0;
  n = 0;

  /* Loop over mesh and adjust weights */

  loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (loc0 > VALID_PTR)
    {
      /* Check mode */

      if (mode == RMX_MODE_MULTI)
        {
          /* Loop over detectors */

          det = (long)RDB[loc0 + RMX_CELL_PTR_DET];
          while (det > VALID_PTR)
            {
              /* Pointer to RMX detector */

              ptr = (long)RDB[det + RMX_DET_PTR_DET0];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Check nonzero and adjust */

              if (RDB[ptr + RMX_DET_MC_RES_TOT] > 0.0)
                {
                  /* Calculate adjustment */

                  wgt = 1.0/RDB[ptr + RMX_DET_MC_RES_TOT];
                  CheckValue(FUNCTION_NAME, "(wgt)", "", wgt, ZERO, INFTY);

                  /* Adjust */

                  WDB[det + RMX_DET_WGT] = RDB[det + RMX_DET_WGT]*wgt;

                  /* Add to counts */

                  sum = sum + RDB[det + RMX_DET_WGT];
                  n++;
                }

              /* Next */

              det = NextItem(det);
            }
        }
      else if (mode == RMX_MODE_GVR)
        {
          /* Pointer to detector */

          det = (long)RDB[loc0 + RMX_CELL_PTR_DET];
          CheckPointer(FUNCTION_NAME, "(det)", DATA_ARRAY, det);

          /* Check multiple */

          if (NextItem(det) > VALID_PTR)
            Die(FUNCTION_NAME, "Multiple responses in GVR mode");

          /* Check nonzero and adjust */

          if (RDB[det + RMX_DET_FWD_CHK] > 0.0)
            {
              /* Calculate adjustment */

              wgt = 1.0/RDB[det + RMX_DET_FWD_CHK];
              CheckValue(FUNCTION_NAME, "(wgt)", "", wgt, ZERO, INFTY);

              /* Adjust */

              WDB[det + RMX_DET_WGT] = RDB[det + RMX_DET_WGT]*wgt;

              /* Add to counts */

              sum = sum + RDB[det + RMX_DET_WGT];
              n++;
            }
        }
      else
        Die(FUNCTION_NAME, "Error in mode");

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Check */

  if ((sum == 0.0) || (n == 0))
    Die(FUNCTION_NAME, "Error in sum");

  /* Loop over mesh and renormalize */

  loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (loc0 > VALID_PTR)
    {
      /* Check mode */

      if (mode == RMX_MODE_MULTI)
        {
          /* Loop over detectors */

          det = (long)RDB[loc0 + RMX_CELL_PTR_DET];
          while (det > VALID_PTR)
            {
              /* Pointer to RMX detector */

              ptr = (long)RDB[det + RMX_DET_PTR_DET0];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Re-normalize weight */

              if (RDB[ptr + RMX_DET_MC_RES_TOT] > 0.0)
                {
                  wgt = ((double)n)*RDB[det + RMX_DET_WGT]/sum;
                  CheckValue(FUNCTION_NAME, "wgt", "", wgt, 0.0, INFTY);

                  /* Put value */

                  WDB[det + RMX_DET_WGT] = wgt;
                }

              /* Next detector */

              det = NextItem(det);
            }
        }
      else if (mode == RMX_MODE_GVR)
        {
          /* Pointer to detector */

          det = (long)RDB[loc0 + RMX_CELL_PTR_DET];
          CheckPointer(FUNCTION_NAME, "(det)", DATA_ARRAY, det);

          /* Re-normalize weight */

          if (RDB[det + RMX_DET_FWD_CHK] > 0.0)
            {
              wgt = ((double)n)*RDB[det + RMX_DET_WGT]/sum;
              CheckValue(FUNCTION_NAME, "wgt", "", wgt, 0.0, INFTY);

              /* Put value */

              WDB[det + RMX_DET_WGT] = wgt;
            }
        }

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Check maximum difference in forward solution ************************/

  /* Reset maximum */

  max = -1.0;

  n = 0;
  nmax = 0;

  /* Loop over mesh */

  loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (loc0 > VALID_PTR)
    {
      /* Loop over detectors */

      det = (long)RDB[loc0 + RMX_CELL_PTR_DET];
      while (det > VALID_PTR)
        {
          /* Check cycle index */

          if ((cyc > 0) &&(RDB[det + RMX_DET_FWD_PREV] != 0.0))
            {
              /* Calculate relative difference to previous */

              chk = fabs(RDB[det + RMX_DET_FWD_CHK]/
                         RDB[det + RMX_DET_FWD_PREV] - 1.0);

              /* Compare to maximum */

              if (chk > max)
                max = chk;

              /* Compare to limit */

              if (chk < RDB[DATA_RMX_MULT_CONVG_LIM])
                n++;

              /* Add total count */

              nmax++;
            }

          /* Remember previous */

          WDB[det + RMX_DET_FWD_PREV] = RDB[det + RMX_DET_FWD_CHK];

          /* Next detector */

          det = NextItem(det);
        }

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Calculate pass fraction */

  if (nmax > 0)
    chk = ((double)n)/((double)nmax);
  else
    chk = 0.0;

  /* Put convergence criteria */

  if (chk < 1E-2)
    {
      WDB[rmx + RMX_CONV] = 1E-2;
      WDB[rmx + RMX_N_ITER] = 1000;
    }
  else if (chk < 5E-2)
    {
      WDB[rmx + RMX_CONV] = 1E-3;
      WDB[rmx + RMX_N_ITER] = 1000;
    }
  else if (chk < 10E-2)
    {
      WDB[rmx + RMX_CONV] = 1E-4;
      WDB[rmx + RMX_N_ITER] = 1000;
    }
  else if (chk < 50E-2)
    {
      WDB[rmx + RMX_CONV] = 1E-5;
      WDB[rmx + RMX_N_ITER] = 1000;
    }
  else
    {
      WDB[rmx + RMX_CONV] = WDB[DATA_RMX_CONV0];
      WDB[rmx + RMX_N_ITER] = WDB[DATA_RMX_N_ITER0];
    }

  /*
  if (cyc == 0)
    {
      WDB[rmx + RMX_CONV] = 1E-2;
      WDB[rmx + RMX_N_ITER] = 1000;
    }
  else if (cyc == 1)
    {
      WDB[rmx + RMX_CONV] = 1E-3;
      WDB[rmx + RMX_N_ITER] = 1000;
    }
  else if (cyc == 2)
    {
      WDB[rmx + RMX_CONV] = 1E-4;
      WDB[rmx + RMX_N_ITER] = 1000;
    }
  else if (cyc == 3)
    {
      WDB[rmx + RMX_CONV] = 1E-5;
      WDB[rmx + RMX_N_ITER] = 1000;
    }
  else
    {
      WDB[rmx + RMX_CONV] = WDB[DATA_RMX_CONV0];
      WDB[rmx + RMX_N_ITER] = WDB[DATA_RMX_N_ITER0];
    }
  */

  /* Check limiting */

  if (RDB[rmx + RMX_CONV] < RDB[DATA_RMX_CONV0])
    WDB[rmx + RMX_CONV] = RDB[DATA_RMX_CONV0];

  if (RDB[rmx + RMX_N_ITER] > RDB[DATA_RMX_N_ITER0])
    WDB[rmx + RMX_N_ITER] = RDB[DATA_RMX_N_ITER0];

  /* Return value */

  return chk;

  /***************************************************************************/
}

/*****************************************************************************/
