/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printinterfaceoutput.c                         */
/*                                                                           */
/* Created:       2012/02/15 (JLe)                                           */
/* Last modified: 2018/11/14 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Prints output for multi-physics interface                    */
/*                                                                           */
/* Comments:    -2014/03/13 Added angular and time dependence for fuel perf. */
/*               interface                                                   */
/*              -TODO: relaxation for fast flux                              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrintInterfaceOutput:"

/*****************************************************************************/

void PrintInterfaceOutput()
{
  long loc0, loc1, loc2, nz, nr, n, ptr, i, j, k, uni, mfile, na, idx;
  long N0, N1, N2, m, n0, n1, n2;
  double zmin, zmax, zl, vol, value, relUnc;
  const double *params;
  char outfile[MAX_STR], object[MAX_STR], location[MAX_STR];
  FILE *fp;

  /* Check that interfaces are defined */

  if ((long)RDB[DATA_PTR_IFC0] < VALID_PTR)
    return;

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Check corrector step */
  /* Print interfaceoutput on corrector if running a */
  /* coupled calculation */

  if (((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP) &&
      (RDB[DATA_RUN_CC] == (double)NO))
    return;

  /* Loop over interfaces */

  loc0 = (long)RDB[DATA_PTR_IFC0];
  while (loc0 > VALID_PTR)
    {
      /* Check if printed */

      if ((long)RDB[loc0 + IFC_CALC_OUTPUT] == NO)
        {
          /* Pointer to next */

          loc0 = NextItem(loc0);

          /* Cycle loop */

          continue;
        }

      /***********************************************************************/

      /***** Parse output file name ******************************************/

      /* Get file name */

      sprintf(outfile, "%s", GetText(loc0 + IFC_PTR_OUTPUT_FNAME));

      /* Check for Matlab m-file format */

      mfile = NO;

      if (strlen(outfile) > 2)
        if ((outfile[strlen(outfile) - 1] == 'm') &&
            (outfile[strlen(outfile) - 2] == '.'))
          {
            /* Adjust name */

            outfile[strlen(outfile) - 2] = '\0';

            /* Set type flag */

            mfile = YES;
          }

      /* Check if burnup mode */

      if ((long)RDB[DATA_BURN_PTR_DEP] > VALID_PTR)
        sprintf(&outfile[strlen(outfile)], "%ld",
                (long)RDB[DATA_BURN_STEP]+(long)RDB[DATA_BURN_STEP_PC]);

      /* Tää indeksi lisättiin UGM 2016 plottauksia varten */
      /*
      if (RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN)
        sprintf(&outfile[strlen(outfile)], "%ld",
                (long)RDB[DATA_DYN_TB]);
      */

      /* Check for matlab format */

      if (mfile == YES)
        {
          /* Adjust file name */

          sprintf(&outfile[strlen(outfile)], ".m");
        }

      /***********************************************************************/

      /***** Write data ******************************************************/

      /* Open file for writing (m*/

      if ((fp = fopen(outfile, "w")) == NULL)
        Error(loc0, "Unable to open output file for writing");

      /* Print variable name in Matlab mode */

      if (mfile == YES)
        fprintf(fp, "IFC%ld = [\n", (long)RDB[loc0 + IFC_IDX]);

      /* Check type */

      if (((long)RDB[loc0 + IFC_TYPE] == IFC_TYPE_FUEP) ||
          ((long)RDB[loc0 + IFC_TYPE] == IFC_TYPE_FPIP))
        {
          /***************************************************************/

          /***** Interface to fuel perfomance codes **********************/

          /* Get pointer */

          loc1 = (long)RDB[loc0 + IFC_PTR_FUEP];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          /* Loop over pins */

          while (loc1 > VALID_PTR)
            {
              /* Pointer to universe */

              uni = (long)RDB[loc1 + IFC_FUEP_PTR_UNI];
              CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

              /* Get dimensions */

              ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_LIM];

              nz = (long)RDB[ptr + FUEP_NZ];
              na = (long)RDB[ptr + FUEP_NA];
              nr = (long)RDB[ptr + FUEP_NR];

              /* Loop over axial and radial zones */
              for (i = 0; i < nz; i++)
                for (j = 0; j < na; j++)
                  for (k = 0; k < nr; k++)
                    {
                      /* Print universe name and indexes */

                      fprintf(fp, "1 %6s %4ld %4ld %4ld ",
                              GetText(uni + UNIVERSE_PTR_NAME), i + 1,
                              j + 1, k + 1);

                      /* Print z-coordinates */

                      ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_Z];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                      fprintf(fp, "%12.5E %12.5E ", RDB[ptr + i],
                              RDB[ptr + i + 1]);

                      /* Print angular limits */

                      ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_PHI];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                      fprintf(fp, "%12.5E %12.5E ", RDB[ptr + j]*360.0/2/PI,
                              RDB[ptr + j + 1]*360.0/2/PI);

                      /* Print radii */

                      ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_R2];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                      fprintf(fp, "%12.5E %12.5E ", sqrt(RDB[ptr + k]),
                              sqrt(RDB[ptr + k + 1]));


                      /* Print result */

                      if (RDB[DATA_RUN_CC] == YES)
                        {
                          /* In coupled calculation print relaxed power */

                          /* Calculate bin index */

                          idx = i + j*nz + k*nz*na;

                          /* Get pointer to relaxed power */

                          ptr = (long)RDB[loc1 + IFC_FUEP_PTR_POWER_REL];
                          CheckPointer(FUNCTION_NAME, "(Pptr1)", DATA_ARRAY, ptr);

                          /* Print relaxed power */

                          fprintf(fp, "%12.5E %7.5f\n", RDB[ptr + idx],
                                  -1.0);
                        }
                      else
                        {
                          /* Otherwise print normal power */

                          ptr = (long)RDB[loc1 + IFC_FUEP_PTR_POWER];
                          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                          fprintf(fp, "%12.5E %7.5f\n", Mean(ptr, i, j, k),
                                  RelErr(ptr, i, j, k));
                        }
                    }

              /* Next pin */

              loc1 = NextItem(loc1);
            }

          /* Get pointer */

          loc1 = (long)RDB[loc0 + IFC_PTR_FUEP];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          /* Loop over pins */

          while (loc1 > VALID_PTR)
            {
              /* Pointer to universe */

              uni = (long)RDB[loc1 + IFC_FUEP_PTR_UNI];
              CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

              /* Fast flux */
              /* Get dimensions */

              ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_FLIM];

              nz = (long)RDB[ptr + FUEP_NZ];
              na = (long)RDB[ptr + FUEP_NA];
              nr = (long)RDB[ptr + FUEP_NR];

              /* Loop over axial and radial zones */
              for (i = 0; i < nz; i++)
                for (j = 0; j < na; j++)
                  for (k = 0; k < nr; k++)
                    {
                      /* Print universe name and indexes */

                      fprintf(fp, "2 %6s %4ld %4ld %4ld ",
                              GetText(uni + UNIVERSE_PTR_NAME), i + 1,
                              j + 1, k + 1);

                      /* Print z-coordinates */

                      ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_FZ];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                      fprintf(fp, "%12.5E %12.5E ", RDB[ptr + i],
                              RDB[ptr + i + 1]);

                      /* Print angular limits */

                      ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_FPHI];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                      fprintf(fp, "%12.5E %12.5E ", RDB[ptr + j]*360.0/2/PI,
                              RDB[ptr + j + 1]*360.0/2/PI);

                      /* Print radii */

                      ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_FR2];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                      fprintf(fp, "%12.5E %12.5E ", sqrt(RDB[ptr + k]),
                              sqrt(RDB[ptr + k + 1]));

                      /* Print result */

                      ptr = (long)RDB[loc1 + IFC_FUEP_PTR_FLUX];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                      fprintf(fp, "%12.5E %7.5f\n", Mean(ptr, i, j, k),
                              RelErr(ptr, i, j, k));
                    }

              /* Next pin */

              loc1 = NextItem(loc1);
            }

          /***************************************************************/
        }
      else if ((long)RDB[loc0 + IFC_TYPE] == IFC_TYPE_TET_MESH)
        {
          /***** Unstructured tetrahedral mesh *******************************/

          /*******************************************************************/

          /* Get pointer to statistics */

          if(RDB[DATA_RUN_CC] == YES)
            {
              /* In coupled calculation get pointer to the relaxed power */

              ptr = (long)RDB[loc0 + IFC_PTR_STAT_REL];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
            }
          else
            {
              ptr = (long)RDB[loc0 + IFC_PTR_STAT];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
            }

          /* Get pointer to cells */

          loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH_PARENTS];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          /* Get pointer to volumes */

          loc2 = (long)RDB[loc0 + IFC_PTR_STAT_VOL];
          CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

          /* Check Matlab format */

          if (mfile == NO)
            {
              /* Get output file path */

              sprintf(location, "%s", GetText(loc0 + IFC_PTR_OUTPUT_FNAME));

              /* Separate object and location */

              i = strlen(location) - 1;
              while (i > 0)
                {
                  if (location[i] == '/')
                    break;

                  i--;
                }

              /* Copy string */

              if (i > 0)
                sprintf(object, "%s", &location[i + 1]);
              else
                sprintf(object, "%s", &location[i]);

              location[i] = '\0';

              /* Print some header data */

              fprintf(fp, "FoamFile\n{\n");
              fprintf(fp, "    version     2.0;\n");
              fprintf(fp, "    format      ascii;\n");
              fprintf(fp, "    class       volScalarField;\n");
              fprintf(fp, "    location    \"%s\";\n", location);
              fprintf(fp, "    object      %s;\n}\n", object);
              fprintf(fp, "\ndimensions      [1 -1 -3 0 0 0 0];\n");
              fprintf(fp, "\ninternalField   nonuniform List<scalar>\n");
              fprintf(fp, "%ld\n(\n", (long)RDB[loc0 + IFC_STAT_NREG]);
            }

          /* Loop over parents */

          while (loc1 > VALID_PTR)
            {

              /* Check stat index */

              if ((i = (long)RDB[loc1 + IFC_TET_PRNT_STAT_IDX]) > -1)
                {
                  /* Get volume */

                  vol = RDB[loc2 + i];

                  /* Check Matlab mode */

                  if (mfile == YES)
                    {
                      if(RDB[DATA_RUN_CC] == NO )
                        fprintf(fp, "%6ld %6ld %12.5E %12.5E %7.5f\n",
                                (long)RDB[loc1 + IFC_TET_PRNT_IDX] + 1,
                                i + 1, vol, Mean(ptr, i), RelErr(ptr, i));

                      else
                        fprintf(fp, "%6ld %6ld %12.5E %12.5E %7.5f\n",
                                (long)RDB[loc1 + IFC_TET_PRNT_IDX] + 1,
                                i + 1, vol, RDB[ptr + i], 0.0);
                    }
                  else if (vol > 0.0)
                    {
                      /* Print in W/m3 */


                      if(RDB[DATA_RUN_CC] == NO )
                        fprintf(fp, "%12.5E\n", 1E+6*Mean(ptr, i)/vol);
                      else
                        fprintf(fp, "%12.5E\n", 1E+6*RDB[ptr + i]/vol);

                    }
                  else if (Mean(ptr, i) > 0.0)
                    Die(FUNCTION_NAME, "Zero volume");
                  else
                    fprintf(fp, "%12.5E\n", 0.0);
                }

              /* Next cell */

              loc1 = NextItem(loc1);
            }

          /* Print remaining OpenFOAM stuff */

          if (mfile == NO)
            {
              fprintf(fp, ")\n;\n\n");

              /* Check pointer to batches */

              if ((loc1 = (long)RDB[loc0 + IFC_PTR_OF_PATCHES]) > VALID_PTR)
                {
                  fprintf(fp, "boundaryField\n{\n");

                  /* Loop over batches */

                  while (loc1 > VALID_PTR)
                    {
                      fprintf(fp, "    %s\n    {\n",
                              GetText(loc1 + IFC_OF_PATCH_PTR_NAME));
                      fprintf(fp, "        type            %s;\n",
                              GetText(loc1 + IFC_OF_PATCH_PTR_TYPE));
                      fprintf(fp, "        value           uniform 0;\n");
                      fprintf(fp, "    }\n");

                      /* Pointer to next */

                      loc1 = NextItem(loc1);
                    }

                  fprintf(fp, "}\n\n");
                }
            }

          /*******************************************************************/
        }
      else if (RDB[loc0 + IFC_OUTPUT_TYPE] == (double)IFC_OUTPUT_SAME_MESH)
        {
          /*******************************************************************/

          /***** Results tallied on the interface regular mesh ***************/

          if(RDB[DATA_RUN_CC] == YES)
            {
              /* In coupled calculation get pointer to the relaxed power */

              ptr = (long)RDB[loc0 + IFC_PTR_STAT_REL];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
            }
          else
            {
              /* Otherwise get direct pointer to statistics */

              ptr = (long)RDB[loc0 + IFC_PTR_STAT];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
            }

          /* Get number of bins */

          n = (long)RDB[loc0 + IFC_STAT_NREG];

          /* Loop over output regions */

          for (i = 0; i < n; i++)
            {
              /* Print results */

              if(RDB[DATA_RUN_CC] == (double)NO)
                fprintf(fp, "%12.5E %7.5f ", Mean(ptr, i),
                        RelErr(ptr, i));
              else
                fprintf(fp, "%12.5E 0.0 ", RDB[ptr + i]);

              /* Newline */

              fprintf(fp, "\n");
            }

        }
      else if (((long)RDB[loc0 + IFC_TYPE] == IFC_TYPE_FET_DENSITY)
          || ((long)RDB[loc0 + IFC_TYPE] == IFC_TYPE_FET_TEMP))
        {
          if(RDB[DATA_RUN_CC] == YES)
            {
              /* In coupled calculation get pointer to the relaxed power */

              ptr = (long)RDB[loc0 + IFC_PTR_STAT_REL];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
            }
          else
            {
              /* Otherwise get direct pointer to statistics */

              ptr = (long)RDB[loc0 + IFC_PTR_STAT];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
            }

          /* Output FETs */

          loc1 = (long)RDB[loc0 + IFC_FET_OUTPUT_PARAMS_PTR];
          params = &RDB[loc1];

          /* Get the basic parameters */

          N0 = (long)params[FET_PARAM_NCOEF0];
          N1 = (long)params[FET_PARAM_NCOEF1];
          N2 = (long)params[FET_PARAM_NCOEF2];

          switch ((long)params[FET_PARAM_TYPE])
            {
            case FET_TYPE_CARTESIAN:
              {
                /* Print a descriptive header */

                fprintf(fp, "%-12s %12s %12s %8s %12s %12s %8s %12s %12s %8s %12s\n",
                        "Type",
                        "X_Min",
                        "X_Max",
                        "X_Order",
                        "Y_Min",
                        "Y_Max",
                        "Y_Order",
                        "Z_Min",
                        "Z_Max",
                        "Z_Order",
                        "#_Regions");
                /*           %-12s %12s %12s %8s %12s %12s %8s %12s %12s %8s %12s */
                fprintf(fp, "%-12s %12.5E %12.5E %8ld %12.5E %12.5E %8ld %12.5E %12.5E %8ld %12ld\nn",
                        "Cartesian",
                        params[FET_CART_MIN_X],
                        params[FET_CART_MAX_X],
                        (long)params[FET_CART_ORDER_X],
                        params[FET_CART_MIN_Y],
                        params[FET_CART_MAX_Y],
                        (long)params[FET_CART_ORDER_Y],
                        params[FET_CART_MIN_Z],
                        params[FET_CART_MAX_Z],
                        (long)params[FET_CART_ORDER_Z],
                        (long)RDB[loc0 + IFC_STAT_NREG]);
                /*           %-12s %12.5E %12.5E %8ld %12.5E %12.5E %8ld %12.5E %12.5E %8ld %12ld */

                /* Print header for data */

                fprintf(fp, "%-8s %-12s %8s %8s %8s %12s   %s\n",
                        "Region",
                        "Linear_Index",
                        "Z_Order",
                        "Y_Order",
                        "X_Order",
                        "Coefficient",
                        "Relative_Uncertainty");
                /*           %-8s %-12s %8s %8s %8s %12s   %s*/

                loc1 = (long)RDB[loc0 + IFC_PTR_OUT];
                while (loc1 > VALID_PTR)
                  {
                    /* Process this region */

                    if ((loc2 = (long)RDB[loc1 + IFC_OUT_PTR_SCORE]) > VALID_PTR)
                      {
                        i = (long)RDB[loc2 + IFC_SCORE_STAT_IDX];

                        for (n0 = 0; n0 < N0; ++n0)
                          for (n1 = 0; n1 < N1; ++n1)
                            for (n2 = 0; n2 < N2; ++n2)
                              {
                                idx = FETIdx(params, n2, n1, n0);
                                FETFinalize(params, ptr, idx, i, &value, &relUnc);
                                fprintf(fp, "%-8ld %-12ld %8ld %8ld %8ld %12.5E   %-12.5E\n",
                                        i,
                                        idx,
                                        n2,
                                        n1,
                                        n0,
                                        value,
                                        relUnc);
                                /*           %-8ld %-12ld %8ld %8ld %8ld %12.5E   %12.5E */
                              }
                      }

                    /* Get the next region */

                    loc1 = NextItem(loc1);
                  }

                break;
              }

            case FET_TYPE_CYLINDRICAL:
              {
                /* Print a descriptive header */

                fprintf(fp, "%-12s %12s %12s %8s %12s %12s %8s %12s\n",
                        "Type",
                        "R_Min",
                        "R_Max",
                        "R_Order",
                        "A_Min",
                        "A_Max",
                        "A_Order",
                        "#_Regions");
                /*           %-12s %12s %12s %8s %12s %12s %8s %12s */
                fprintf(fp, "%-12s %12.5E %12.5E %8ld %12.5E %12.5E %8ld %12ld\n\n",
                        "Cylindrical",
                        0.0,
                        params[FET_CYL_MAX_R],
                        (long)params[FET_CYL_ORDER_R],
                        params[FET_CYL_MIN_A],
                        params[FET_CYL_MAX_A],
                        (long)params[FET_CYL_ORDER_A],
                        (long)RDB[loc0 + IFC_STAT_NREG]);
                /*           %-12s %12.5E %12.5E %8ld %12.5E %12.5E %8ld %12ld */

                /* Print header for data */

                fprintf(fp, "%-8s %-12s %8s %8s %8s %12s   %s\n",
                        "Region",
                        "Linear_Index",
                        "R_Order",
                        "R_Rank",
                        "A_Order",
                        "Coefficient",
                        "Relative_Uncertainty");
                /*           %-8s %-12s %8s %8s %8s %12s   %s*/

                loc1 = (long)RDB[loc0 + IFC_PTR_OUT];
                while (loc1 > VALID_PTR)
                  {
                    /* Process this region */

                    if ((loc2 = (long)RDB[loc1 + IFC_OUT_PTR_SCORE]) > VALID_PTR)
                      {
                        i = (long)RDB[loc2 + IFC_SCORE_STAT_IDX];

                        for (n0 = 0, n = 0, m = 0; n0 < N0; ++n0)
                          {
                            for (n1 = 0; n1 < N1; ++n1)
                              {
                                idx = FETIdx(params, 0, n1, n0);
                                FETFinalize(params, ptr, idx, i, &value, &relUnc);
                                fprintf(fp, "%-8ld %-12ld %8ld %8ld %8ld %12.5E   %-12.5E\n",
                                        i,
                                        idx,
                                        n,
                                        m,
                                        n1,
                                        value,
                                        relUnc);
                                /*           %-8ld %-12ld %8ld %8ld %8ld %12.5E   %12.5E */
                              }

                            /* Track the Zernike rank and order rather than computing it  */
                            /* from n0 à la:                                              */
                            /*   n = ceil((sqrt(9 + 8 * n0) - 3) / 2)                     */
                            /*   m = 2 * n0 - n * (n + 2)                                 */

                            if (m == n)
                              m = -(++n);
                            else
                              m += 2;
                          }
                      }

                    /* Get the next region */

                    loc1 = NextItem(loc1);
                  }

                break;
              }

#ifdef DEBUG
            default:
              {
                Die(FUNCTION_NAME, "Unsupported FET type %ld", (long)params[FET_PARAM_TYPE]);
              }
#endif /* DEBUG */
            }
        }
      else
        {
          /*******************************************************************/

          /***** Other types *************************************************/

          if(RDB[DATA_RUN_CC] == YES)
            {
              /* In coupled calculation get pointer to the relaxed power */

              ptr = (long)RDB[loc0 + IFC_PTR_STAT_REL];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
            }
          else
            {
              /* Otherwise get direct pointer to statistics */

              ptr = (long)RDB[loc0 + IFC_PTR_STAT];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
            }

          /* get parameters */

          nz = (long)RDB[loc0 + IFC_NZ];
          zmin = RDB[loc0 + IFC_ZMIN];
          zmax = RDB[loc0 + IFC_ZMAX];
          nr = (long)RDB[loc0 + IFC_NR];
          n = (long)RDB[loc0 + IFC_STAT_NREG];

          /* Calculate length */

          zl = (zmax - zmin)/((double)nz);

          /* Loop over output regions */

          loc1 = (long)RDB[loc0 + IFC_PTR_OUT];
          while (loc1 > VALID_PTR)
            {
              /* Pointer to statistics */

              if ((loc2 = (long)RDB[loc1 + IFC_OUT_PTR_SCORE]) > VALID_PTR)
                {
                  /* Get stat index */

                  i = (long)RDB[loc2 + IFC_SCORE_STAT_IDX];
                  CheckValue(FUNCTION_NAME, "i", "", i, 0, n);

                  /* Loop over axial and radial segments */

                  for (j = 0; j < nz; j++)
                    for (k = 0; k < nr; k++)
                      {

                        /* Calculate index for relaxed tally */

                        idx = i + j*n + k*n*nz;

                        /* Print bottom coordinates */

                        fprintf(fp, "%12.5E %12.5E %12.5E ",
                                RDB[loc1 + IFC_OUT_X0],
                                RDB[loc1 + IFC_OUT_Y0], zmin + (double)j*zl);

                        /* Print top coordinates */

                        fprintf(fp, "%12.5E %12.5E %12.5E ",
                                RDB[loc1 + IFC_OUT_X0],
                                RDB[loc1 + IFC_OUT_Y0],
                                zmin + (double)(j + 1)*zl);

                        /* Print radii */

                        fprintf(fp, "%12.5E ", sqrt((double)k/((double)nr))*
                                RDB[loc1 + IFC_OUT_R]);
                        fprintf(fp, "%12.5E ",
                                sqrt((double)(k + 1)/((double)nr))*
                                RDB[loc1 + IFC_OUT_R]);

                        /* Print results */

                        if(RDB[DATA_RUN_CC] == (double)NO)
                          fprintf(fp, "%12.5E %7.5f ", Mean(ptr, i, j, k),
                                  RelErr(ptr, i, j, k));
                        else
                          fprintf(fp, "%12.5E 0.0 ", RDB[ptr + idx]);

                        /* Newline */

                        fprintf(fp, "\n");
                      }
                }

              /* Next region */

              loc1 = NextItem(loc1);
            }

          /*******************************************************************/
        }

      /* Print delimiter in Matlab mode */

      if (mfile == YES)
        fprintf(fp, "];\n");

      /* Close file */

      fclose(fp);

      /***********************************************************************/

      /* Print errors for OpenFOAM interfaces */

      if ((long)RDB[loc0 + IFC_TYPE] == IFC_TYPE_TET_MESH)
        {
          /* Print relative errors */

          if ((((long)RDB[loc0 + IFC_PRINT_ERROR] == IFC_PRINT_ERROR_REL) ||
               ((long)RDB[loc0 + IFC_PRINT_ERROR] == IFC_PRINT_ERROR_BOTH))
              && (mfile == NO))
            {

              /********************************************************************/

              /***** Parse output file name ***************************************/

              /* Get file name */

              sprintf(outfile, "%s_relerr", GetText(loc0 + IFC_PTR_OUTPUT_FNAME));

              /* Check if burnup mode */

              if ((long)RDB[DATA_BURN_PTR_DEP] > VALID_PTR)
                sprintf(&outfile[strlen(outfile)], "%ld",
                        (long)RDB[DATA_BURN_STEP]);

              /***********************************************************************/

              /***** Write data ******************************************************/

              /* Open file for writing (m*/

              if ((fp = fopen(outfile, "w")) == NULL)
                Error(loc0, "Unable to open output file for writing");

              /* Get pointer to statistics (latest iteration) */

              ptr = (long)RDB[loc0 + IFC_PTR_STAT];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Get pointer to cells */

              loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH_PARENTS];
              CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

              /* Get pointer to volumes */

              loc2 = (long)RDB[loc0 + IFC_PTR_STAT_VOL];
              CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

              /* Get output file path */

              sprintf(location, "%s_relerr", GetText(loc0 + IFC_PTR_OUTPUT_FNAME));

              /* Separate object and location */

              i = strlen(location) - 1;
              while (i > 0)
                {
                  if (location[i] == '/')
                    break;

                  i--;
                }

              /* Copy string */

              if (i > 0)
                sprintf(object, "%s", &location[i + 1]);
              else
                sprintf(object, "%s", &location[i]);

              location[i] = '\0';

              /* Print some header data */

              fprintf(fp, "FoamFile\n{\n");
              fprintf(fp, "    version     2.0;\n");
              fprintf(fp, "    format      ascii;\n");
              fprintf(fp, "    class       volScalarField;\n");
              fprintf(fp, "    location    \"%s\";\n", location);
              fprintf(fp, "    object      %s;\n}\n", object);
              fprintf(fp, "\ndimensions      [0 0 0 0 0 0 0];\n");
              fprintf(fp, "\ninternalField   nonuniform List<scalar>\n");
              fprintf(fp, "%ld\n(\n", (long)RDB[loc0 + IFC_STAT_NREG]);

              /* Loop over parents */

              while (loc1 > VALID_PTR)
                {
                  /* Check stat index */

                  if ((i = (long)RDB[loc1 + IFC_TET_PRNT_STAT_IDX]) > -1)
                    {
                      /* Get volume */

                      vol = RDB[loc2 + i];

                      if (vol > 0.0)
                        {
                          /* Relative error shouldn't be divided by volume */

                          fprintf(fp, "%12.5E\n", RelErr(ptr, i));

                        }
                      else if (Mean(ptr, i) > 0.0)
                        {
                          /* Zero volume / unable to calculate cell volume */

                          /* Print in W */

                          fprintf(fp, "%12.5E\n", RelErr(ptr, i));
                        }
                      else
                        fprintf(fp, "%12.5E\n", 0.0);
                    }

                  /* Next cell */

                  loc1 = NextItem(loc1);
                }

              /* Print remaining OpenFOAM stuff */

              fprintf(fp, ")\n;\n\n");

              /* Check pointer to batches */

              if ((loc1 = (long)RDB[loc0 + IFC_PTR_OF_PATCHES]) > VALID_PTR)
                {
                  fprintf(fp, "boundaryField\n{\n");

                  /* Loop over batches */

                  while (loc1 > VALID_PTR)
                    {
                      fprintf(fp, "    %s\n    {\n",
                              GetText(loc1 + IFC_OF_PATCH_PTR_NAME));
                      fprintf(fp, "        type            %s;\n",
                              GetText(loc1 + IFC_OF_PATCH_PTR_TYPE));
                      fprintf(fp, "        value           uniform 0;\n");
                      fprintf(fp, "    }\n");

                      /* Pointer to next */

                      loc1 = NextItem(loc1);
                    }

                  fprintf(fp, "}\n\n");
                }


              /*******************************************************************/

              /* Close file */

              fclose(fp);

            }

          /* Print absolute errors */

          if ((((long)RDB[loc0 + IFC_PRINT_ERROR] == IFC_PRINT_ERROR_ABS) ||
               ((long)RDB[loc0 + IFC_PRINT_ERROR] == IFC_PRINT_ERROR_BOTH))
              && (mfile == NO))
            {

              /********************************************************************/

              /***** Parse output file name ***************************************/

              /* Get file name */

              sprintf(outfile, "%s_abserr", GetText(loc0 + IFC_PTR_OUTPUT_FNAME));

              /* Check if burnup mode */

              if ((long)RDB[DATA_BURN_PTR_DEP] > VALID_PTR)
                sprintf(&outfile[strlen(outfile)], "%ld",
                        (long)RDB[DATA_BURN_STEP]);

              /***********************************************************************/

              /***** Write data ******************************************************/

              /* Open file for writing (m*/

              if ((fp = fopen(outfile, "w")) == NULL)
                Error(loc0, "Unable to open output file for writing");

              /* Get pointer to statistics (latest iteration) */

              ptr = (long)RDB[loc0 + IFC_PTR_STAT];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Get pointer to cells */

              loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH_PARENTS];
              CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

              /* Get pointer to volumes */

              loc2 = (long)RDB[loc0 + IFC_PTR_STAT_VOL];
              CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

              /* Get output file path */

              sprintf(location, "%s_abserr", GetText(loc0 + IFC_PTR_OUTPUT_FNAME));

              /* Separate object and location */

              i = strlen(location) - 1;
              while (i > 0)
                {
                  if (location[i] == '/')
                    break;

                  i--;
                }

              /* Copy string */

              if (i > 0)
                sprintf(object, "%s", &location[i + 1]);
              else
                sprintf(object, "%s", &location[i]);

              location[i] = '\0';

              /* Print some header data */

              fprintf(fp, "FoamFile\n{\n");
              fprintf(fp, "    version     2.0;\n");
              fprintf(fp, "    format      ascii;\n");
              fprintf(fp, "    class       volScalarField;\n");
              fprintf(fp, "    location    \"%s\";\n", location);
              fprintf(fp, "    object      %s;\n}\n", object);
              fprintf(fp, "\ndimensions      [1 2 -3 0 0 0 0];\n");
              fprintf(fp, "\ninternalField   nonuniform List<scalar>\n");
              fprintf(fp, "%ld\n(\n", (long)RDB[loc0 + IFC_STAT_NREG]);

              /* Loop over parents */

              while (loc1 > VALID_PTR)
                {
                  /* Check stat index */

                  if ((i = (long)RDB[loc1 + IFC_TET_PRNT_STAT_IDX]) > -1)
                    {
                      /* Get volume */

                      vol = RDB[loc2 + i];

                      if (vol > 0.0)
                        {
                          /* Print in W/m3 */

                          fprintf(fp, "%12.5E\n", Mean(ptr,i)*RelErr(ptr, i)/vol);

                        }
                      else if (Mean(ptr, i) > 0.0)
                        {
                          /* Zero volume / unable to calculate cell volume */

                          /* Print in W */

                          fprintf(fp, "%12.5E\n", Mean(ptr,i)*RelErr(ptr, i));
                        }
                      else
                        fprintf(fp, "%12.5E\n", 0.0);
                    }

                  /* Next cell */

                  loc1 = NextItem(loc1);
                }

              /* Print remaining OpenFOAM stuff */

              fprintf(fp, ")\n;\n\n");

              /* Check pointer to batches */

              if ((loc1 = (long)RDB[loc0 + IFC_PTR_OF_PATCHES]) > VALID_PTR)
                {
                  fprintf(fp, "boundaryField\n{\n");

                  /* Loop over batches */

                  while (loc1 > VALID_PTR)
                    {
                      fprintf(fp, "    %s\n    {\n",
                              GetText(loc1 + IFC_OF_PATCH_PTR_NAME));
                      fprintf(fp, "        type            %s;\n",
                              GetText(loc1 + IFC_OF_PATCH_PTR_TYPE));
                      fprintf(fp, "        value           uniform 0;\n");
                      fprintf(fp, "    }\n");

                      /* Pointer to next */

                      loc1 = NextItem(loc1);
                    }

                  fprintf(fp, "}\n\n");
                }


              /*******************************************************************/

              /* Close file */

              fclose(fp);

            }
        }

      /* Next interface */

      loc0 = NextItem(loc0);
    }
}

/*****************************************************************************/
