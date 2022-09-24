/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : detectoroutput.c                               */
/*                                                                           */
/* Created:       2011/03/03 (JLe)                                           */
/* Last modified: 2020/04/10 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Prints detector output                                       */
/*                                                                           */
/* Comments: - Time binejä ei pidä sallia classic lookin kanssa              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DetectorOutput:"

/*****************************************************************************/

void DetectorOutput()
{
  long det0, erg, ptr, n0, n1, n2, n, m, i, j, loc0, tmp;
  long ebins0, ubins0, cbins0, mbins0, lbins0, rbins0, zbins0, ybins0, xbins0;
  long tbins0, eb0, ub0, cb0, mb0, lb0, rb0, zb0, yb0, xb0, tb0, idx0;
  double fetValue, fetRelUnc;
  double min0, max0, min1, max1, min2, max2, x, y, x0, y0, pitch;
  FILE *fp;
  char tmpstr[MAX_STR];

  /* Check detector definitions */

  if ((long)RDB[DATA_PTR_DET0] < VALID_PTR)
    return;

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Check if in active cycles */

  if (RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP])
    return;

#ifdef STAB_BURN

  /* Open file for writing. Modified to save pred. and corr. separately (AIs)*/

  if ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP)
    sprintf(tmpstr, "%s_det%ldc%d.m", GetText(DATA_PTR_INPUT_FNAME),
            (long)RDB[DATA_BURN_STEP], (int)RDB[DATA_BURN_CI_I]);
  else
    sprintf(tmpstr, "%s_det%ldp.m", GetText(DATA_PTR_INPUT_FNAME),
            (long)RDB[DATA_BURN_STEP]);

#else

  /* Check corrector step */

  if (((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP) &&
      ((long)RDB[DATA_BURN_SIE] == NO))
    return;

  /* Check special indexes */

  if ((long)RDB[DATA_RUN_VR_ITER] == YES)
    sprintf(tmpstr, "%s_det%ld.m", GetText(DATA_PTR_INPUT_FNAME),
            (long)RDB[DATA_VR_ITER_IDX]);
  else if ((long)RDB[DATA_COEF_CALC_IDX] < 0)
    {
      if (RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN)
        sprintf(tmpstr, "%s_det%ld.m", GetText(DATA_PTR_INPUT_FNAME),
                (long)RDB[DATA_DYN_TB]);
      else
        sprintf(tmpstr, "%s_det%ld.m", GetText(DATA_PTR_INPUT_FNAME),
                (long)RDB[DATA_BURN_STEP]);
    }
  else
    sprintf(tmpstr, "%s_det%ldb%ld.m", GetText(DATA_PTR_INPUT_FNAME),
            (long)RDB[DATA_COEF_CALC_BU_IDX],
            (long)RDB[DATA_COEF_CALC_IDX]);

#endif

  if ((fp = fopen(tmpstr, "w")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file for writing");

  /* Loop over detectors */

  det0 = (long)RDB[DATA_PTR_DET0];
  while(det0 > VALID_PTR)
    {
      /* Get number of bins */

      ebins0 = (long)RDB[det0 + DET_N_EBINS];
      ubins0 = (long)RDB[det0 + DET_N_UBINS];
      cbins0 = (long)RDB[det0 + DET_N_CBINS];
      mbins0 = (long)RDB[det0 + DET_N_MBINS];
      lbins0 = (long)RDB[det0 + DET_N_LBINS];
      rbins0 = (long)RDB[det0 + DET_N_RBINS];
      tbins0 = (long)RDB[det0 + DET_N_TBINS];

      /* Mesh bins */

      if ((ptr = (long)RDB[det0 + DET_PTR_MESH]) > VALID_PTR)
        {
          xbins0 = (long)RDB[ptr + MESH_N0];
          ybins0 = (long)RDB[ptr + MESH_N1];
          zbins0 = (long)RDB[ptr + MESH_N2];
        }
      else
        {
          xbins0 = 1;
          ybins0 = 1;
          zbins0 = 1;
        }

      /* FET bins */

      if ((ptr = (long)RDB[det0 + DET_FET_PTR_PARAMS]) > VALID_PTR)
        {
          if ((long)RDB[ptr + FET_PARAM_NCOEF0] > 0)
            xbins0 = (long)RDB[ptr + FET_PARAM_NCOEF0];
          if ((long)RDB[ptr + FET_PARAM_NCOEF1] > 0)
            ybins0 = (long)RDB[ptr + FET_PARAM_NCOEF1];
          if ((long)RDB[ptr + FET_PARAM_NCOEF2] > 0)
            zbins0 = (long)RDB[ptr + FET_PARAM_NCOEF2];
        }

      /* Check mode */

      if (tbins0 == 1)
        {
          /*******************************************************************/

          /***** Serpent 1 type output (no time bins) ************************/

          fprintf(fp, "\n");

          fprintf(fp, "DET%s = [\n", GetText(det0 + DET_PTR_NAME));

          /* Time bins are not used */

          tb0 = 0;

          /* Pointer to statistics */

          ptr = (long)RDB[det0 + DET_PTR_STAT];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Loop over bins */

          n0 = 0;

          for (eb0 = 0; eb0 < ebins0; eb0++)
          for (ub0 = 0; ub0 < ubins0; ub0++)
          for (cb0 = 0; cb0 < cbins0; cb0++)
          for (mb0 = 0; mb0 < mbins0; mb0++)
          for (lb0 = 0; lb0 < lbins0; lb0++)
          for (rb0 = 0; rb0 < rbins0; rb0++)
          for (zb0 = 0; zb0 < zbins0; zb0++)
          for (yb0 = 0; yb0 < ybins0; yb0++)
          for (xb0 = 0; xb0 < xbins0; xb0++)
            {
              /* Print indexes */

              fprintf(fp, "%5ld ", n0 + 1);
              fprintf(fp, "%4ld ", eb0 + 1);
              fprintf(fp, "%4ld ", ub0 + 1);
              fprintf(fp, "%4ld ", cb0 + 1);
              fprintf(fp, "%4ld ", mb0 + 1);
              fprintf(fp, "%4ld ", lb0 + 1);
              fprintf(fp, "%4ld ", rb0 + 1);
              fprintf(fp, "%4ld ", zb0 + 1);
              fprintf(fp, "%4ld ", yb0 + 1);
              fprintf(fp, "%4ld ", xb0 + 1);

              if ((tmp = (long)RDB[det0 + DET_FET_PTR_PARAMS]) > VALID_PTR)
                {
                  /* FET detector, get the corresponding bin index */

                  idx0 = FETIdx(&RDB[tmp], zb0, yb0, xb0);

                  FETFinalize(&RDB[tmp], (long)RDB[det0 + DET_PTR_STAT], idx0,
                              0, &fetValue, &fetRelUnc);

                  /* Print coefficient and relative uncertainty */

                  fprintf(fp, "%12.5E %12.5E ", fetValue, fetRelUnc);
                }
              else
                {
                  /* Get index */

                  idx0 = DetIdx(det0, eb0, ub0, cb0, mb0, lb0, zb0, yb0,
                                xb0, tb0);

                  /* Print mean */

                  fprintf(fp, "%12.5E ", Mean(ptr, idx0, rb0));

                  /* Print relative statistical error */

                  fprintf(fp, "%7.5f ", RelErr(ptr, idx0, rb0));
                }

              /* Print newline */

              fprintf(fp, "\n");

              /* Update index */

              n0++;
            }

          fprintf(fp, "];\n\n");

          /*******************************************************************/
        }
      else
        {
          /*******************************************************************/

          /***** Serpent 1 type output (with time bins) **********************/

          fprintf(fp, "\n");

          fprintf(fp, "DET%s = [\n", GetText(det0 + DET_PTR_NAME));

          /* Pointer to statistics */

          ptr = (long)RDB[det0 + DET_PTR_STAT];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Loop over bins */

          n0 = 0;

          for (tb0 = 0; tb0 < tbins0; tb0++)
          for (eb0 = 0; eb0 < ebins0; eb0++)
          for (ub0 = 0; ub0 < ubins0; ub0++)
          for (cb0 = 0; cb0 < cbins0; cb0++)
          for (mb0 = 0; mb0 < mbins0; mb0++)
          for (lb0 = 0; lb0 < lbins0; lb0++)
          for (rb0 = 0; rb0 < rbins0; rb0++)
          for (zb0 = 0; zb0 < zbins0; zb0++)
          for (yb0 = 0; yb0 < ybins0; yb0++)
          for (xb0 = 0; xb0 < xbins0; xb0++)
            {
              /* Print indexes */

              fprintf(fp, "%5ld ", n0 + 1);
              fprintf(fp, "%4ld ", tb0 + 1);
              fprintf(fp, "%4ld ", eb0 + 1);
              fprintf(fp, "%4ld ", ub0 + 1);
              fprintf(fp, "%4ld ", cb0 + 1);
              fprintf(fp, "%4ld ", mb0 + 1);
              fprintf(fp, "%4ld ", lb0 + 1);
              fprintf(fp, "%4ld ", rb0 + 1);
              fprintf(fp, "%4ld ", zb0 + 1);
              fprintf(fp, "%4ld ", yb0 + 1);
              fprintf(fp, "%4ld ", xb0 + 1);

              /* Get index */

              idx0 = DetIdx(det0, eb0, ub0, cb0, mb0, lb0, zb0, yb0, xb0, tb0);

              /* Print mean */

              fprintf(fp, "%12.5E ", Mean(ptr, idx0, rb0));

              /* Print relative statistical error */

              fprintf(fp, "%7.5f ", RelErr(ptr, idx0, rb0));

              /* Print newline */

              fprintf(fp, "\n");

              /* Update index */

              n0++;
            }

          fprintf(fp, "];\n\n");

          /*******************************************************************/
        }

      /***********************************************************************/

      /***** Print energy intervals ******************************************/

      if ((erg = (long)RDB[det0 + DET_PTR_EGRID]) > VALID_PTR)
        {
          /* Pointer to values */

          ptr = (long)RDB[erg + ENERGY_GRID_PTR_DATA];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Loop over energy bins and print values */

          fprintf(fp, "\nDET%sE = [\n", GetText(det0 + DET_PTR_NAME));

          for (n0 = 0; n0 < (long)RDB[erg + ENERGY_GRID_NE] - 1; n0++)
            fprintf(fp, "%12.5E %12.5E %12.5E\n", RDB[ptr + n0],
                    RDB[ptr + n0 + 1],
                    (RDB[ptr + n0] + RDB[ptr + n0 + 1])/2.0);

          fprintf(fp, "];\n");
        }

      /***********************************************************************/

      /***** Print time intervals ********************************************/

      if ((ptr = (long)RDB[det0 + DET_PTR_TME]) > VALID_PTR)
        {
          /* Loop over time bins and print values */

          fprintf(fp, "\nDET%sT = [\n", GetText(det0 + DET_PTR_NAME));

          for (n0 = 0; n0 < (long)RDB[det0 + DET_N_TBINS]; n0++)
            fprintf(fp, "%12.5E %12.5E %12.5E\n", RDB[ptr + n0],
                    RDB[ptr + n0 + 1],
                    (RDB[ptr + n0] + RDB[ptr + n0 + 1])/2.0);

          fprintf(fp, "];\n");
        }

      /***********************************************************************/

      /***** Print mesh intervals ********************************************/

      /* Pointer to mesh */

      if ((ptr = (long)RDB[det0 + DET_PTR_MESH]) > VALID_PTR)
        {
          /* Get sizes */

          n0 = (long)RDB[ptr + MESH_N0];
          n1 = (long)RDB[ptr + MESH_N1];
          n2 = (long)RDB[ptr + MESH_N2];

          /* Get dimensions */

          min0 = RDB[ptr + MESH_MIN0];
          max0 = RDB[ptr + MESH_MAX0];
          min1 = RDB[ptr + MESH_MIN1];
          max1 = RDB[ptr + MESH_MAX1];
          min2 = RDB[ptr + MESH_MIN2];
          max2 = RDB[ptr + MESH_MAX2];

          /* Check type */

          if ((long)RDB[ptr + MESH_TYPE] == MESH_TYPE_CARTESIAN)
            {
              /* x-direction */

              if (n0 > 0)
                {
                  fprintf(fp, "\nDET%sX = [\n", GetText(det0 + DET_PTR_NAME));

                  for (n = 0; n < n0; n++)
                    fprintf(fp, "%12.5E %12.5E %12.5E\n",
                            ((double)n)/((double)n0)*(max0 - min0) + min0,
                            ((double)n + 1.0)/((double)n0)*(max0 - min0) + min0,
                            ((double)n + 0.5)/((double)n0)*(max0 - min0) + min0);

                  fprintf(fp, "];\n");
                }

              /* y-direction */

              if (n1 > 0)
                {
                  fprintf(fp, "\nDET%sY = [\n", GetText(det0 + DET_PTR_NAME));

                  for (n = 0; n < n1; n++)
                    fprintf(fp, "%12.5E %12.5E %12.5E\n",
                            ((double)n)/((double)n1)*(max1 - min1) + min1,
                            ((double)n + 1.0)/((double)n1)*(max1 - min1) + min1,
                            ((double)n + 0.5)/((double)n1)*(max1 - min1) + min1);

                  fprintf(fp, "];\n");
                }

              /* z-direction */

              if (n2 > 0)
                {
                  fprintf(fp, "\nDET%sZ = [\n", GetText(det0 + DET_PTR_NAME));

                  for (n = 0; n < n2; n++)
                    fprintf(fp, "%12.5E %12.5E %12.5E\n",
                            ((double)n)/((double)n2)*(max2 - min2) + min2,
                            ((double)n + 1.0)/((double)n2)*(max2 - min2) + min2,
                            ((double)n + 0.5)/((double)n2)*(max2 - min2) + min2);

                  fprintf(fp, "];\n");
                }
            }
          else if (((long)RDB[ptr + MESH_TYPE] == MESH_TYPE_HEXX) ||
                   ((long)RDB[ptr + MESH_TYPE] == MESH_TYPE_HEXY))
            {
              /* Print cell center coordinates */

              if ((n0 > 0) && (n1 > 0))
                {
                  /* Get pitch */

                  pitch = RDB[ptr + MESH_MAX0];

                  /* Get center coordinates */

                  x0 = RDB[ptr + MESH_MIN0] + (1 - (n0 % 2))*0.5*pitch;
                  y0 = RDB[ptr + MESH_MIN1] + (1 - (n1 % 2))*0.5*pitch;

                  /* Print */

                  fprintf(fp, "\nDET%sCOORD = [\n",
                          GetText(det0 + DET_PTR_NAME));

                  /* Avoid compiler warning */

                  x = 0.0;
                  y = 0.0;

                  /* Loop over lattice */

                  j = -(long)((double)n1/2.0);
                  for (n = 0; n < n1; n++)
                    {
                      i = -(long)((double)n0/2.0);
                      for (m = 0; m < n0; m++)
                        {
                          if ((long)RDB[ptr + MESH_TYPE] == MESH_TYPE_HEXX)
                            {
                              x = x0 + (i + COS60*j)*pitch;
                              y = y0 + j*SIN60*pitch;
                            }
                          else if ((long)RDB[ptr + MESH_TYPE] == MESH_TYPE_HEXY)
                            {
                              x = x0 + j*SIN60*pitch;
                              y = y0 + (i + COS60*j)*pitch;
                            }

                          i++;

                          fprintf(fp, "%E %E\n", x, y);
                        }
                      j++;
                    }

                  fprintf(fp, "];\n");
                }

              /* z-direction */

              if (n2 > 0)
                {
                  fprintf(fp, "\nDET%sZ = [\n", GetText(det0 + DET_PTR_NAME));

                  for (n = 0; n < n2; n++)
                    fprintf(fp, "%12.5E %12.5E %12.5E\n",
                            ((double)n)/((double)n2)*(max2 - min2) + min2,
                            ((double)n + 1.0)/((double)n2)*(max2 - min2) + min2,
                            ((double)n + 0.5)/((double)n2)*(max2 - min2) + min2);

                  fprintf(fp, "];\n");
                }
            }
          else if ((long)RDB[ptr + MESH_TYPE] == MESH_TYPE_CYLINDRICAL)
            {
              /* r-direction */

              if (n0 > 0)
                {
                  fprintf(fp, "\nDET%sR = [\n", GetText(det0 + DET_PTR_NAME));

                  for (n = 0; n < n0; n++)
                    fprintf(fp, "%12.5E %12.5E %12.5E\n",
                            ((double)n)/((double)n0)*(max0 - min0) + min0,
                            ((double)n + 1.0)/((double)n0)*(max0 - min0) + min0,
                            ((double)n + 0.5)/((double)n0)*(max0 - min0) + min0);

                  fprintf(fp, "];\n");
                }

              /* phi-direction */

              if (n1 > 0)
                {
                  fprintf(fp, "\nDET%sPHI = [\n", GetText(det0 + DET_PTR_NAME));

                  min1 = min1*180.0/PI;
                  max1 = max1*180.0/PI;

                  for (n = 0; n < n1; n++)
                    fprintf(fp, "%12.5E %12.5E %12.5E\n",
                            ((double)n)/((double)n1)*(max1 - min1) + min1,
                            ((double)n + 1.0)/((double)n1)*(max1 - min1) + min1,
                            ((double)n + 0.5)/((double)n1)*(max1 - min1) + min1);

                  fprintf(fp, "];\n");
                }

              /* z-direction */

              if (n2 > 0)
                {
                  fprintf(fp, "\nDET%sZ = [\n", GetText(det0 + DET_PTR_NAME));

                  for (n = n2 - 1; n > -1; n--)
                    fprintf(fp, "%12.5E %12.5E %12.5E\n",
                            ((double)n)/((double)n2)*(max2 - min2) + min2,
                            ((double)n + 1.0)/((double)n2)*(max2 - min2) + min2,
                            ((double)n + 0.5)/((double)n2)*(max2 - min2) + min2);

                  fprintf(fp, "];\n");
                }
            }
          else if ((long)RDB[ptr + MESH_TYPE] == MESH_TYPE_SPHERICAL)
            {
              /* r-direction */

              if (n0 > 0)
                {
                  fprintf(fp, "\nDET%sR = [\n", GetText(det0 + DET_PTR_NAME));

                  for (n = 0; n < n0; n++)
                    fprintf(fp, "%12.5E %12.5E %12.5E\n",
                            ((double)n)/((double)n0)*(max0 - min0) + min0,
                            ((double)n + 1.0)/((double)n0)*(max0 - min0) + min0,
                            ((double)n + 0.5)/((double)n0)*(max0 - min0) + min0);

                  fprintf(fp, "];\n");
                }

              /* phi-direction */

              if (n1 > 0)
                {
                  fprintf(fp, "\nDET%sPHI = [\n", GetText(det0 + DET_PTR_NAME));

                  min1 = min1*180.0/PI;
                  max1 = max1*180.0/PI;

                  for (n = 0; n < n1; n++)
                    fprintf(fp, "%12.5E %12.5E %12.5E\n",
                            ((double)n)/((double)n1)*(max1 - min1) + min1,
                            ((double)n + 1.0)/((double)n1)*(max1 - min1) + min1,
                            ((double)n + 0.5)/((double)n1)*(max1 - min1) + min1);

                  fprintf(fp, "];\n");
                }

              /* theta-direction */

              if (n2 > 0)
                {
                  fprintf(fp, "\nDET%sTHETA = [\n", GetText(det0 + DET_PTR_NAME));

                  min2 = min2*180.0/PI;
                  max2 = max2*180.0/PI;

                  for (n = 0; n < n2; n++)
                    fprintf(fp, "%12.5E %12.5E %12.5E\n",
                            ((double)n)/((double)n2)*(max2 - min2) + min2,
                            ((double)n + 1.0)/((double)n2)*(max2 - min2) + min2,
                            ((double)n + 0.5)/((double)n2)*(max2 - min2) + min2);

                  fprintf(fp, "];\n");
                }
            }
          else if ((long)RDB[ptr + MESH_TYPE] == MESH_TYPE_ICYL)
            {
              /* Radial direction */

              if (n0 > 0)
                {
                  fprintf(fp, "\nDET%sR = [\n", GetText(det0 + DET_PTR_NAME));

                  loc0 = (long)RDB[ptr + MESH_ORTHO_PTR_XLIM];
                  CheckPointer(FUNCTION_NAME, "loc0", DATA_ARRAY, loc0);

                  for (n = 0; n < n0; n++)
                    fprintf(fp, "%12.5E %12.5E %12.5E\n",
                            RDB[loc0 + n], RDB[loc0 + n + 1],
                            0.5*(RDB[loc0 + n] + RDB[loc0 + n + 1]));

                  fprintf(fp, "];\n");
                }

              /* Angular direction */

              if (n1 > 0)
                {
                  fprintf(fp, "\nDET%sPHI = [\n",
                          GetText(det0 + DET_PTR_NAME));

                  loc0 = (long)RDB[ptr + MESH_ORTHO_PTR_YLIM];
                  CheckPointer(FUNCTION_NAME, "loc0", DATA_ARRAY, loc0);

                  for (n = 0; n < n1; n++)
                    fprintf(fp, "%12.5E %12.5E %12.5E\n",
                            RDB[loc0 + n]*180.0/PI,
                            RDB[loc0 + n + 1]*180.0/PI,
                            0.5*(RDB[loc0 + n] + RDB[loc0 + n + 1])*180.0/PI);

                  fprintf(fp, "];\n");
                }

              /* z-direction */

              if (n2 > 0)
                {
                  fprintf(fp, "\nDET%sZ = [\n", GetText(det0 + DET_PTR_NAME));

                  loc0 = (long)RDB[ptr + MESH_ORTHO_PTR_ZLIM];
                  CheckPointer(FUNCTION_NAME, "loc0", DATA_ARRAY, loc0);

                  for (n = 0; n < n2; n++)
                    fprintf(fp, "%12.5E %12.5E %12.5E\n",
                            RDB[loc0 + n], RDB[loc0 + n + 1],
                            0.5*(RDB[loc0 + n] + RDB[loc0 + n + 1]));

                  fprintf(fp, "];\n");
                }
            }
          else if ((long)RDB[ptr + MESH_TYPE] == MESH_TYPE_ORTHOGONAL)
            {
              /* x-direction */

              if (n0 > 0)
                {
                  fprintf(fp, "\nDET%sX = [\n", GetText(det0 + DET_PTR_NAME));

                  loc0 = (long)RDB[ptr + MESH_ORTHO_PTR_XLIM];
                  CheckPointer(FUNCTION_NAME, "loc0", DATA_ARRAY, loc0);

                  for (n = 0; n < n0; n++)
                    fprintf(fp, "%12.5E %12.5E %12.5E\n",
                            RDB[loc0 + n], RDB[loc0 + n + 1],
                            0.5*(RDB[loc0 + n] + RDB[loc0 + n + 1]));

                  fprintf(fp, "];\n");
                }

              /* y-direction */

              if (n1 > 0)
                {
                  fprintf(fp, "\nDET%sY = [\n", GetText(det0 + DET_PTR_NAME));

                  loc0 = (long)RDB[ptr + MESH_ORTHO_PTR_YLIM];
                  CheckPointer(FUNCTION_NAME, "loc0", DATA_ARRAY, loc0);

                  for (n = 0; n < n1; n++)
                    fprintf(fp, "%12.5E %12.5E %12.5E\n",
                            RDB[loc0 + n], RDB[loc0 + n + 1],
                            0.5*(RDB[loc0 + n] + RDB[loc0 + n + 1]));

                  fprintf(fp, "];\n");
                }

              /* z-direction */

              if (n2 > 0)
                {
                  fprintf(fp, "\nDET%sZ = [\n", GetText(det0 + DET_PTR_NAME));

                  loc0 = (long)RDB[ptr + MESH_ORTHO_PTR_ZLIM];
                  CheckPointer(FUNCTION_NAME, "loc0", DATA_ARRAY, loc0);

                  for (n = 0; n < n2; n++)
                    fprintf(fp, "%12.5E %12.5E %12.5E\n",
                            RDB[loc0 + n], RDB[loc0 + n + 1],
                            0.5*(RDB[loc0 + n] + RDB[loc0 + n + 1]));

                  fprintf(fp, "];\n");
                }
            }
          else
            Die(FUNCTION_NAME, "Invalid mesh type");
        }


      /***********************************************************************/

      /***** Print FET degrees ***********************************************/

      /* Pointer to FET */

      if ((ptr = (long)RDB[det0 + DET_FET_PTR_PARAMS]) > VALID_PTR)
        {
          switch (n = (long)RDB[ptr + FET_PARAM_TYPE])
            {
            case FET_TYPE_CARTESIAN:
              {
                fprintf(fp, "\nDET%sFETCART = [\n",
                        GetText(det0 + DET_PTR_NAME));
                fprintf(fp, "%12.5E %12.5E %4ld\n",
                        RDB[ptr + FET_CART_MIN_X],
                        RDB[ptr + FET_CART_MAX_X],
                        (long)RDB[ptr + FET_CART_NCOEF_X]);
                fprintf(fp, "%12.5E %12.5E %4ld\n",
                        RDB[ptr + FET_CART_MIN_Y],
                        RDB[ptr + FET_CART_MAX_Y],
                        (long)RDB[ptr + FET_CART_NCOEF_Y]);
                fprintf(fp, "%12.5E %12.5E %4ld\n",
                        RDB[ptr + FET_CART_MIN_Z],
                        RDB[ptr + FET_CART_MAX_Z],
                        (long)RDB[ptr + FET_CART_NCOEF_Z]);
                fprintf(fp, "];\n");
                break;
              }

            case FET_TYPE_CYLINDRICAL:
              {
                /* Zernike (circular) configuration */

                fprintf(fp, "\nDET%sFETCYL_Z = [\n",
                        GetText(det0 + DET_PTR_NAME));
                fprintf(fp, "%12.5E %4ld\n",
                        RDB[ptr + FET_CYL_MAX_R],
                        (long)RDB[ptr + FET_CYL_NCOEF_R]);
                fprintf(fp, "];\n");

                /* Legendre (axial) configuration */

                fprintf(fp, "\nDET%sFETCYL_L = [\n",
                        GetText(det0 + DET_PTR_NAME));
                fprintf(fp, "%12.5E %12.5E %4ld\n",
                        RDB[ptr + FET_CYL_MIN_A],
                        RDB[ptr + FET_CYL_MAX_A],
                        (long)RDB[ptr + FET_CYL_NCOEF_A]);
                fprintf(fp, "];\n");

                /* Orientation and center */

                fprintf(fp, "\nDET%sFETCYL_O = [\n",
                        GetText(det0 + DET_PTR_NAME));
                fprintf(fp, "%4ld %12.5E %12.5E\n",
                        (long)RDB[ptr + FET_CYL_ORIENTATION_A],
                        RDB[ptr + FET_CYL_CENTER_A0],
                        RDB[ptr + FET_CYL_CENTER_A1]);
                fprintf(fp, "];\n");
                break;
              }

#ifdef DEBUG
            default:
              {
                Die(FUNCTION_NAME, "Unsupported FET type %ld", n);
              }
#endif /* DEBUG */
            }
        }

      /***********************************************************************/

      /* Next detector */

      det0 = NextItem(det0);
    }

  /* Close file */

  fclose(fp);

  /* Write iteration output (10.3.2020 / 2.1.32 / JLe) */

  UserIter(3);
}

/*****************************************************************************/
