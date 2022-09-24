/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : useriter.c                                     */
/*                                                                           */
/* Created:       2020/03/10 (JLe)                                           */
/* Last modified: 2020/03/12 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: User-defined iteration routines                              */
/*                                                                           */
/* Comments: - Toimii OK jos sauvat ei mene ääripäiden yli.                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "UserIter:"

/*****************************************************************************/

void UserIter(long part)
{
  long np, loc0, ptr, n, i, nt, uni, tra, det, mod;
  double zmin, zmax, keff, prof[1000], tot;
  const double *vol;
  char tmpstr[MAX_STR];
  FILE *fp;

  /* Check mode */

  if ((long)RDB[DATA_ITER_MODE] != ITER_MODE_USER)
    return;

  /* Get number of parameters */

  np = (long)RDB[DATA_ITER_USR_N_PARAM];
  CheckValue(FUNCTION_NAME, "(np)", "", np, 0, 1000000);

  /* Get pointer */

  if (np > 0)
    {
      /* Get pointer */

      loc0 = (long)RDB[DATA_ITER_USR_PTR_PARAM];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);
    }
  else
    loc0 = -1;

  /* Check part */

  if (part == 1)
    {
      /***********************************************************************/

      /***** Preprocessing stage *********************************************/

      /* Get mode */

      mod = AtoF(GetText(loc0++), FUNCTION_NAME, NULL, -1);

      /* Check pointer */

      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      /* Minimum and maximum dimensions */

      zmin = AtoF(GetText(loc0), FUNCTION_NAME, NULL, -1);
      WDB[loc0++] = zmin;

      zmax = AtoF(GetText(loc0), FUNCTION_NAME, NULL, -1);
      WDB[loc0++] = zmax;

      /* Read number of iterated transformations */

      nt = AtoI(GetText(loc0), FUNCTION_NAME, NULL, -1);
      WDB[loc0++] = (double)nt;

      /* Loop over transformations */

      for (n = 0; n < nt; n++)
        {
          /* Loop over universes */

          uni = (long)RDB[DATA_PTR_U0];
          while (uni > VALID_PTR)
            {
              /* Compare names */

              if (CompareStr(uni + UNIVERSE_PTR_NAME, loc0))
                break;

              /* next */

              uni = NextItem(uni);
            }

          /* Avoid compiler warning */

          tra = -1;

          /* Check pointer */

          if (uni < VALID_PTR)
            Error(0, "Universe %s is not defined", GetText(loc0));
          else if ((tra = (long)RDB[uni + UNIVERSE_PTR_TRANS]) < VALID_PTR)
            Error(0, "Universe %s is not associated with a transformation",
                  GetText(loc0));
          else
            WDB[loc0++] = tra;

          /* Put initial height */

          WDB[tra + TRANS_Z0] = 0.5*(zmin + zmax);
        }

      /* Allocate memory for stats */

      ptr = NewStat("pos", 1, 2*nt);
      AllocStatHistory(ptr);
      WDB[DATA_ITER_USR_PTR_DATA] = (double)ptr;

      /* Allocate memory for profile */

      ptr = ReallocMem(DATA_ARRAY, nt);
      WDB[DATA_ITER_USR_PTR_PROFILE] = (double)ptr;

      if (mod == 1)
        {
          /*******************************************************************/

          /***** Use assembly-averaged powers ********************************/

          /* Find detector */

          det = (long)RDB[DATA_PTR_DET0];
          while (det > VALID_PTR)
            {
              /* Compare names */

              if (CompareStr(det + DET_PTR_NAME, loc0))
                break;

              /* next */

              det = NextItem(det);
            }

          /* Check pointer */

          if (det < VALID_PTR)
            Error(0, "Detector %s is not defined", GetText(loc0));
          else
            WDB[loc0++] = det;

          /* Check number of bins */

          if ((long)RDB[det + DET_N_TOT_BINS] != nt)
            Error(0, "Mismatch in number of detector bins");

          /* Convert volumes from text to values */

          for (n = 0; n < nt; n++)
            WDB[loc0 + n] = AtoF(GetText(loc0 + n), FUNCTION_NAME, NULL, -1);

          /*******************************************************************/
        }
      else if (mod == 2)
        {
          /*******************************************************************/

          /***** Use maximum pin-powers **************************************/

          /* Loop over detectors */

          for (n = 0; n < nt; n++)
            {
              /* Find detector */

              det = (long)RDB[DATA_PTR_DET0];
              while (det > VALID_PTR)
                {
                  /* Compare names */

                  if (CompareStr(det + DET_PTR_NAME, loc0))
                    break;

                  /* next */

                  det = NextItem(det);
                }

              /* Check pointer */

              if (det < VALID_PTR)
                Error(0, "Detector %s is not defined", GetText(loc0));
              else
                WDB[loc0++] = det;
            }

          /*******************************************************************/
        }
      else
        Die(FUNCTION_NAME, "Invalid mode");

       /***********************************************************************/
    }
  else if (part == 3)
    {
      /***********************************************************************/

      /***** Post-processing stage *******************************************/

      /* Skip mode */

      loc0++;

      /* Check pointer */

      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      /* Minimum and maximum dimensions */

      zmin = RDB[loc0++];
      zmax = RDB[loc0++];

      /* Read number of iterated transformation */

      nt = (long)RDB[loc0++];

      /* Pointer to data */

      ptr = (long)RDB[DATA_ITER_USR_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Set file name */

      sprintf(tmpstr, "%s_rods%ld.m", GetText(DATA_PTR_INPUT_FNAME),
              (long)RDB[DATA_BURN_STEP]);

      /* Open file */

      if ((fp = fopen(tmpstr, "w")) == NULL)
        Die(FUNCTION_NAME, "Unable to open file for writing");

      /* Print data */

      fprintf(fp, "rods = [\n");

      for (n = 0; n < nt; n++)
        fprintf(fp, "%1.5E %1.5f %1.5E %1.5f\n", Mean(ptr, n), RelErr(ptr, n),
                Mean(ptr, nt + n), RelErr(ptr, nt + n));

      fprintf(fp, "];\n");

      /* Close file */

      fclose(fp);

      /***********************************************************************/
    }
  else if (part == 2)
    {
      /***********************************************************************/

      /***** Iteration stage *************************************************/

      /* Get mode */

      mod = AtoF(GetText(loc0++), FUNCTION_NAME, NULL, -1);

      /* Check pointer */

      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      /* Minimum and maximum dimensions */

      zmin = RDB[loc0++];
      zmax = RDB[loc0++];

      /* Get cycle k-eff */

      keff = RDB[DATA_CYCLE_KEFF];

      /* Read number of iterated transformations */

      nt = (long)RDB[loc0++];

      /* Check cycle */

      if (RDB[DATA_CYCLE_IDX] == RDB[DATA_CRIT_SKIP]
          + 2.0*RDB[DATA_ITER_NCYC] - 1.0)
        {
          /* Pointer to stored positions */

          ptr = (long)RDB[DATA_ITER_USR_PTR_DATA];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Calculate profile */

          tot = 0;
          for (n = 0; n < nt; n++)
            {
              prof[n] = Mean(ptr,n);
              tot = tot + prof[n];
            }

          /* Pointer to profile */

          ptr = (long)RDB[DATA_ITER_USR_PTR_PROFILE];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          for (n = 0; n < nt; n++)
            {
              /* Check */

              if (tot == 0.0)
                prof[n] = 1.0;
              else if (prof[n] > 0.0)
                prof[n] = tot/prof[n]/((double)nt);
              else
                prof[n] = INFTY;
            }

         /* Loop over transformations */

          for (n = 0; n < nt; n++)
            {
              /* Get pointer to transformation */

              tra = (long)RDB[loc0 + n];
              CheckPointer(FUNCTION_NAME, "(tra)", DATA_ARRAY, tra);

              /* Update */

              WDB[tra + TRANS_Z0] = RDB[tra + TRANS_Z0]*prof[n]/keff;
                WDB[tra + TRANS_Z0] = zmax;
            }
        }

      /* Check mode */

      if (mod == 1)
        {
          /*******************************************************************/

          /***** Use assembly-wise powers ************************************/

          /* Check cycle */

          if (RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP]
              + 2.0*RDB[DATA_ITER_NCYC])
            {
              /* Obtain profile from detectors */

              det = (long)RDB[loc0 + nt];
              CheckPointer(FUNCTION_NAME, "(det)", DATA_ARRAY, det);

              /* Pointer to volumes */

              vol = &RDB[loc0 + nt + 1];

              /* Pointer to stat */

              ptr = (long)RDB[det + DET_PTR_STAT];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Read values */

              tot = 0.0;
              for (n = 0; n < nt; n++)
                {
                  /* Get power */

                  if (vol[n] > 0.0)
                    prof[n] = BufVal(ptr, n, 0)/vol[n];
                  else
                    prof[n] = 0.0;

                  tot = tot + prof[n];
                }

              /* Normalize */

              for (n = 0; n < nt; n++)
                {
                  /* Check */

                  if (tot == 0.0)
                    prof[n] = 1.0;
                  else if (prof[n] > 0.0)
                    prof[n] = tot/prof[n]/((double)nt);
                  else
                    prof[n] = INFTY;
                }
            }
          else
            {
              /* Use fixed profile */

              for (n = 0; n < nt; n++)
                prof[n] = 1.0;
            }

          /* Loop over transformations */

          for (n = 0; n < nt; n++)
            {
              /* Get pointer to transformation */

              tra = (long)RDB[loc0++];
              CheckPointer(FUNCTION_NAME, "(tra)", DATA_ARRAY, tra);

              /* Update */

              WDB[tra + TRANS_Z0] = RDB[tra + TRANS_Z0]*prof[n]/keff;

              /* Check */
              /*
              if (RDB[tra + TRANS_Z0] < zmin)
                WDB[tra + TRANS_Z0] = zmin;
              else if (RDB[tra + TRANS_Z0] > zmax)
                WDB[tra + TRANS_Z0] = zmax;
              */
              /* Store data */

              ptr = (long)RDB[DATA_ITER_USR_PTR_DATA];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddStat(RDB[tra + TRANS_Z0], ptr, n);
              AddStat(prof[n], ptr, n + nt);
            }

          /*******************************************************************/
        }
      else if (mod == 2)
        {
          /*******************************************************************/

          /***** Use maximum pin-powers **************************************/

          /* Check cycle */

          if (RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP]
              + 2.0*RDB[DATA_ITER_NCYC])
            {
              /* Reset total */

              tot = 0.0;

              /* Obtain profile from detectors */

              for (n = 0; n < nt; n++)
                {
                  /* Pointer to detector */

                  det = (long)RDB[loc0 + nt + n];
                  CheckPointer(FUNCTION_NAME, "(det)", DATA_ARRAY, det);

                  /* Pointer to stat */

                  ptr = (long)RDB[det + DET_PTR_STAT];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                  /* Reset maximum */

                  prof[n] = 0.0;

                  /* Loop over values */

                  for (i = 0; i < (long)RDB[det + DET_N_TOT_BINS]; i++)
                    if (BufVal(ptr, i, 0) > prof[n])
                      prof[n] = BufVal(ptr, i, 0);

                  /* Add to total */

                  tot = tot + prof[n];
                }

              /* Normalize */

              for (n = 0; n < nt; n++)
                {
                  /* Check */

                  if (tot == 0.0)
                    prof[n] = 1.0;
                  else if (prof[n] > 0.0)
                    prof[n] = tot/prof[n]/((double)nt);
                  else
                    prof[n] = INFTY;
                }
            }
          else
            {
              /* Use fixed profile */

              for (n = 0; n < nt; n++)
                prof[n] = 1.0;
            }

          /* Loop over transformations */

          for (n = 0; n < nt; n++)
            {
              /* Get pointer to transformation */

              tra = (long)RDB[loc0++];
              CheckPointer(FUNCTION_NAME, "(tra)", DATA_ARRAY, tra);

              /* Update */

              WDB[tra + TRANS_Z0] = RDB[tra + TRANS_Z0]*prof[n]/keff;

              /* Check */
              /*
              if (RDB[tra + TRANS_Z0] < zmin)
                WDB[tra + TRANS_Z0] = zmin;
              else if (RDB[tra + TRANS_Z0] > zmax)
                WDB[tra + TRANS_Z0] = zmax;
              */
              /* Store data */

              ptr = (long)RDB[DATA_ITER_USR_PTR_DATA];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddStat(RDB[tra + TRANS_Z0], ptr, n);
              AddStat(prof[n], ptr, n + nt);
            }

          /*******************************************************************/
        }
      else
        Die(FUNCTION_NAME, "Invalid mode");
    }
  else
    Die(FUNCTION_NAME, "Invalid part");

  /***************************************************************************/
}

/*****************************************************************************/
