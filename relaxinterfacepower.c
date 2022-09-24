/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : relaxinterfacepower.c                          */
/*                                                                           */
/* Created:       2014/07/07 (VVa)                                           */
/* Last modified: 2018/08/09 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Relaxes power solution for multi-physics interface           */
/*                                                                           */
/* Comments: -Add relaxation for interface flux                              */
/*           -Add convergence output for all interface types                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "RelaxInterfacePower:"

/*****************************************************************************/

void RelaxInterfacePower()
{
  long loc0, loc1, ptr, ptr1, ptr2, n, nz, nr, i, j, k, na, type, idx;
  long nst, uni;
  double alpha, sumpow, pow1, pow2, D, mult;
  double momnew, momold, relnew, relold;
  double maxeps, rel_maxeps, maxdiff, rel_maxdiff;
  double L2abs, rel_L2abs, L2rel, rel_L2rel;
  char tmpstr[MAX_STR];
  FILE *fout;

  /* Check that interfaces are defined */

  if ((loc0 = (long)RDB[DATA_PTR_IFC0]) < VALID_PTR)
    return;

  fprintf(outp, "Relaxing interface powers...\n");

  /* Get relaxation alpha from memory */

  alpha = RDB[DATA_SOL_REL_ALPHA];
  /*
  fprintf(outp, "Relaxing powers, SOL_REL_ITER %ld alpha = %f\n", (long)RDB[DATA_SOL_REL_ITER], alpha);
  */
  /* Get relaxation factor (< 1.0 underrelaxes) */

  /* Check if no relaxation is wanted */

  if ((D = RDB[DATA_SOL_REL_FACT]) == 0.0)
    {
      D = 1;
      alpha = 1;
    }

  /* Loop over interfaces */

  while (loc0 > VALID_PTR)
    {

      /* Reset power counters */

      sumpow = 0.0;
      pow1 = 0.0;
      pow2 = 0.0;

      /* Reset maximum of convergence criterion */

      maxeps  = 0.0;
      maxdiff = 0.0;
      L2abs   = 0.0;
      L2rel   = 0.0;

      rel_maxeps  = 0.0;
      rel_maxdiff = 0.0;
      rel_L2abs   = 0.0;
      rel_L2rel   = 0.0;

      /* Cycle loop if no output is requested */

      if (RDB[loc0 + IFC_CALC_OUTPUT] == (double)NO)
        {
          /* Next interface */

          loc0 = NextItem(loc0);

          /* Cycle loop*/

          continue;
        }

      /* Get interface type */

      type = (long)RDB[loc0 + IFC_TYPE];

      if ((type == IFC_TYPE_FUEP) || (type == IFC_TYPE_FPIP))
        {
          /* Fuel behavior interface */

          /* Get pointer to first rod */

          loc1 = (long)RDB[loc0 + IFC_PTR_FUEP];

          /* Loop over pins to relax power */

          while (loc1 > VALID_PTR)
            {
              /* Pointer to universe */

              uni = (long)RDB[loc1 + IFC_FUEP_PTR_UNI];
              CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

              /* Pointer to nest */

              nst = (long)RDB[uni + UNIVERSE_PTR_NEST];
              CheckPointer(FUNCTION_NAME, "(nst)", DATA_ARRAY, nst);

              /* Get pointer to output limits */

              ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_LIM];
              CheckPointer(FUNCTION_NAME, "(limptr)", DATA_ARRAY, ptr);

              /* get number of bins*/

              nz = (long)RDB[ptr + FUEP_NZ];

              na = (long)RDB[ptr + FUEP_NA];

              nr = (long)RDB[ptr + FUEP_NR];

              /* Get number of nests from nest count or FINIX definition */

              mult = RDB[nst + NEST_COUNT];

              if ((ptr = (long)RDB[loc1 + IFC_FUEP_PTR_FINIX]) > VALID_PTR)
                if (RDB[ptr + FINIX_N_RODS] > 0.0)
                  mult = RDB[ptr + FINIX_N_RODS];

              /* Get pointer to power tally */

              ptr = (long)RDB[loc1 + IFC_FUEP_PTR_POWER];
              CheckPointer(FUNCTION_NAME, "(Pptr)", DATA_ARRAY, ptr);

              /* Get pointer to relaxed power */

              ptr1 = (long)RDB[loc1 + IFC_FUEP_PTR_POWER_REL];
              CheckPointer(FUNCTION_NAME, "(Pptr1)", DATA_ARRAY, ptr1);

              /* Get pointer to power at previous iteration */

              ptr2 = (long)RDB[loc1 + IFC_FUEP_PTR_POWER_PREV];
              CheckPointer(FUNCTION_NAME, "(Pptr2)", DATA_ARRAY, ptr2);

              for (i = 0; i < nz; i++)
                for (j = 0; j < na; j++)
                  for (k = 0; k < nr; k++)
                    {

                      /* Calculate index for relaxed tally */

                      idx = i + j*nz + k*nz*na;

                      /* Add to sum of relaxed power */

                      pow2 = pow2 + RDB[ptr1 + idx]*mult;

                      /* Get new and old momentary powers */

                      momnew = Mean(ptr, i, j, k);
                      momold = RDB[ptr2 + idx];

                      /* Get previous relaxed power */

                      relold = RDB[ptr1 + idx];

                      /* Calculate new relaxed power */

                      if (WDB[DATA_SOL_REL_ITER] == (double)0)
                        relnew = momnew;
                      else
                        relnew = relold - alpha*D*(relold - momnew);

                      /********************/
                      /* Store new values */
                      /********************/

                      /* Store new relaxed power */

                      WDB[ptr1 + idx] = relnew;

                      /* Store new momentary power to be used as previous */
                      /* momentary power */

                      WDB[ptr2 + idx] = momnew;

                      /**********************************************************/
                      /* Convergence criterions based on momentary distribution */
                      /**********************************************************/

                      CalcConvCriteria(momnew, momold,
                                       &maxdiff, &maxeps, &L2abs, &L2rel);

                      /********************************************************/
                      /* Convergence criterions based on relaxed distribution */
                      /********************************************************/

                      CalcConvCriteria(relnew, relold,
                                       &rel_maxdiff, &rel_maxeps,
                                       &rel_L2abs, &rel_L2rel);

                      /*************************************/
                      /* Add to sums of power or debugging */
                      /*************************************/

                      /* Add to sum of power */

                      sumpow = sumpow + RDB[ptr1 + idx]*mult;

                      /* Add to sum of momentary power */

                      pow1 = pow1 + momnew*mult;

                    }

              /* Next pin */

              loc1 = NextItem(loc1);
            }
        }

      else if ((RDB[loc0 + IFC_OUTPUT_TYPE] == (double)IFC_OUTPUT_SAME_MESH) ||
               ((long)RDB[loc0 + IFC_TYPE] == IFC_TYPE_TET_MESH))
        {
          /* Simple 1D array */

          /* Get pointer to power tally */

          ptr = (long)RDB[loc0 + IFC_PTR_STAT];
          CheckPointer(FUNCTION_NAME, "(Pptr)", DATA_ARRAY, ptr);

          /* Get pointer to relaxed power */

          ptr1 = (long)RDB[loc0 + IFC_PTR_STAT_REL];
          CheckPointer(FUNCTION_NAME, "(Pptr1)", DATA_ARRAY, ptr1);

          /* Get pointer to previous momentary power */

          ptr2 = (long)RDB[loc0 + IFC_PTR_STAT_PREV];
          CheckPointer(FUNCTION_NAME, "(Pptr1)", DATA_ARRAY, ptr1);

          /* Get number of regions */

          n = (long)RDB[loc0 + IFC_STAT_NREG];

          for (i = 0; i < n; i++)
            {

              /* Get new momentary power */

              momnew = Mean(ptr, i);

              /* Get old momentary power */

              momold = RDB[ptr2 + i];

              /* Get old relaxed power */

              relold = RDB[ptr1 + i];

              /* Calculate new relaxed power */
              /* If this is the first iteration, then use plain value */
              /* Otherwise relax                                      */

              if(WDB[DATA_SOL_REL_ITER] == (double)0)
                relnew = momnew;
              else
                relnew = relold - alpha*D*(relold - momnew);

              /* Store new relaxed power */

              WDB[ptr1 + i] = relnew;

              /* Store new momentary power to be used as previous */
              /* momentary power */

              WDB[ptr2 + i] = momnew;

              /**********************************************************/
              /* Convergence criterions based on momentary distribution */
              /**********************************************************/

              CalcConvCriteria(momnew, momold,
                               &maxdiff, &maxeps, &L2abs, &L2rel);

              /********************************************************/
              /* Convergence criterions based on relaxed distribution */
              /********************************************************/

              CalcConvCriteria(relnew, relold,
                               &rel_maxdiff, &rel_maxeps,
                               &rel_L2abs, &rel_L2rel);

              /* Add to sum of previous relaxed power */

              pow2 = pow2 + relold;

              /* Add to sum of relaxed power */

              sumpow = sumpow + relnew;

              /* Add to sum of momentary power */

              pow1 = pow1 + momnew;

            }

        }
      else
        {

          /* Get pointer to power tally */

          ptr = (long)RDB[loc0 + IFC_PTR_STAT];
          CheckPointer(FUNCTION_NAME, "(Pptr)", DATA_ARRAY, ptr);

          /* Get pointer to relaxed power */

          ptr1 = (long)RDB[loc0 + IFC_PTR_STAT_REL];
          CheckPointer(FUNCTION_NAME, "(Pptr1)", DATA_ARRAY, ptr1);

          /* Get pointer to previous momentary power */

          ptr2 = (long)RDB[loc0 + IFC_PTR_STAT_PREV];
          CheckPointer(FUNCTION_NAME, "(Pptr1)", DATA_ARRAY, ptr1);

          /* Get number of regions */

          n = (long)RDB[loc0 + IFC_STAT_NREG];

          /* Get number of regions */

          nz = (long)RDB[loc0 + IFC_NZ];

          /* Get number of regions */

          nr = (long)RDB[loc0 + IFC_NR];

          for (i = 0; i < n; i++)
            {

              for (j = 0; j < nz; j++)
                {

                  for (k = 0; k < nr; k++)
                    {

                      /* Calculate index for relaxed tally */

                      idx = i + j*n + k*n*nz;

                      /* Get new momentary power */

                      momnew = Mean(ptr, i, j, k);

                      /* Get old momentary power */

                      momold = RDB[ptr2 + idx];

                      /* Get old relaxed power */

                      relold = RDB[ptr1 + idx];

                      /* Calculate new relaxed power */
                      /* If this is the first iteration, then use plain value */
                      /* Otherwise relax                                      */

                      if(WDB[DATA_SOL_REL_ITER] == (double)0)
                        relnew = momnew;
                      else
                        relnew = relold - alpha*D*(relold - momnew);

                      /* Store new relaxed power */

                      WDB[ptr1 + idx] = relnew;

                      /* Store new momentary power to be used as previous */
                      /* momentary power */

                      WDB[ptr2 + idx] = momnew;

                      /**********************************************************/
                      /* Convergence criterions based on momentary distribution */
                      /**********************************************************/

                      CalcConvCriteria(momnew, momold,
                                       &maxdiff, &maxeps, &L2abs, &L2rel);

                      /********************************************************/
                      /* Convergence criterions based on relaxed distribution */
                      /********************************************************/

                      CalcConvCriteria(relnew, relold,
                                       &rel_maxdiff, &rel_maxeps,
                                       &rel_L2abs, &rel_L2rel);

                      /* Add to sum of previous relaxed power */

                      pow2 = pow2 + relold;

                      /* Add to sum of relaxed power */

                      sumpow = sumpow + relnew;

                      /* Add to sum of momentary power */

                      pow1 = pow1 + momnew;

                    }

                }

            }

        }

      /* Open output file for convergence */

      if (WDB[DATA_SOL_REL_ITER] == (double)0)
        {
          /* On first iteration open a new file */

          sprintf(tmpstr, "%s_Pconv%ld.m", GetText(loc0 + IFC_PTR_INPUT_FNAME),
                  (long)RDB[DATA_BURN_STEP]);

          fout = fopen(tmpstr, "w");

          /* Reset idx */

          fprintf(fout,"\nidx = 1;\n\n");
        }
      else
        {
          /* On subsequent iterations append to old file */

          sprintf(tmpstr, "%s_Pconv%ld.m", GetText(loc0 + IFC_PTR_INPUT_FNAME),
                  (long)RDB[DATA_BURN_STEP]);

          fout = fopen(tmpstr, "a");

          /* Increment idx */

          fprintf(fout,"\nidx = idx + 1;\n\n");
        }

      /* Write out convergence criteria of momentary power distribution */

      fprintf(fout, "P_eps(idx) = %E;\n", maxeps);
      fprintf(fout, "P_delta(idx) = %E;\n", maxdiff);
      fprintf(fout, "P_L2_of_absolute(idx) = %E;\n", sqrt(L2abs));
      fprintf(fout, "P_L2_of_relative(idx) = %E;\n", sqrt(L2rel));

      /* Write out convergence criteria of relaxed power distribution */

      fprintf(fout, "relaxedP_eps(idx) = %E;\n", rel_maxeps);
      fprintf(fout, "relaxedP_delta(idx) = %E;\n", rel_maxdiff);
      fprintf(fout, "relaxedP_L2_of_absolute(idx) = %E;\n", sqrt(rel_L2abs));
      fprintf(fout, "relaxedP_L2_of_relative(idx) = %E;\n", sqrt(rel_L2rel));

      /* Close output file */

      fclose(fout);

      /* Get timestep beginning time */
      /*
      t0 = RDB[DATA_TIME_CUT_TMIN];
      */
      /* Get timestep end time */
      /*
      t1 = RDB[DATA_TIME_CUT_TMAX];
      */
      /*
      fprintf(outp, "\nPower = %E (was: %E mom: %E) alpha %E D %E\n", sumpow/(t1-t0), pow2/(t1-t0), pow1/(t1-t0), alpha, D);
      */

      fprintf(outp, "\nPower = %E (was: %E mom: %E) alpha %E D %E\n", sumpow, pow2, pow1, alpha, D);


      /* Next interface */

      loc0 = NextItem(loc0);
    }

  fprintf(outp, "OK.\n\n");

  return;

  /***************************************************************************/
}

/*****************************************************************************/
