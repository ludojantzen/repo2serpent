/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readoftemperatures.c                           */
/*                                                                           */
/* Created:       2017/08/03 (VVa)                                           */
/* Last modified: 2017/11/29 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Reads OpenFOAM multi-physics interface temperature file      */
/*                                                                           */
/* Comments:   -                                                             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadOFTemperatures:"

/*****************************************************************************/

void ReadOFTemperatures(long ifc, long update)
{
  char *line, tmpstr[MAX_STR], fname[MAX_STR];
  long n, nc, dim[7], i, tmplist, cgns, type, disttype;
  double mul, T, Tmax, Tmin, T0, old, new;
  double maxeps, maxdiff, L2abs, L2rel, value;
  FILE *fp;

   /* Get temperature file filename */

  sprintf(fname, "%s", GetText(ifc + IFC_PTR_OF_TFILE));

  /* Open temperature file for reading */

  if ((fp = fopen(fname, "r")) == NULL)
    Error(ifc, "Temperature file \"%s\" does not exist", fname);

  /* Read header data */

  ReadOFHeader(fp, &n, &i, (long *)dim, &disttype, &value);

  /* Reset multiplier */

  mul = 1.0;

  /* Get number of cells */

  nc = (long)RDB[ifc + IFC_NC_PARENTS];

  /* Check units */

  if (dim[3] == 1)
    {
      /* K */

      mul = 1.0;
    }
  else if (dim[3] == 0)
    {
      /* Header might be missing altogether (K) */

      mul = 1.0;

    }
  else
    Die(FUNCTION_NAME, "Undefined dimensions %ld", dim[3]);

  /* Check type */
  /*
    if (n != OF_FILE_TEMP)
    Note(ifc, "Possibly invalid format for temperature file");
  */
  /* Check size */

  if ((i != nc) && (disttype != OF_INTERNAL_FIELD_UNIFORM))
    Error(ifc, "Invalid number of entries in temperature file");

  /* Create temperature list or read from memory */

  if (update == NO)
    {
      /* Allocate memory for list */

      tmplist = ReallocMem(DATA_ARRAY, nc);

      /* Store pointer to list */

      WDB[ifc + IFC_PTR_TMP_LIST] = (double)tmplist;
    }
  else
    {
      /* Get pointer to list */

      tmplist = (long)RDB[ifc + IFC_PTR_TMP_LIST];
    }

  /* Reset maximum of convergence criterion */

  maxeps  = 0.0;
  maxdiff = 0.0;
  L2abs   = 0.0;
  L2rel   = 0.0;

  /* Reset maximum and minimum temperature */

  Tmax = -1.0;
  Tmin = INFTY;

  /* Get temperature type flag */

  type = (long)RDB[ifc + IFC_TEMPERATURE_TYPE];

  /* Get nominal temperature */

  T0 = RDB[ifc + IFC_NOMINAL_TEMPERATURE];

  /* Loop over parent cells */

  cgns = (long)RDB[ifc + IFC_PTR_TET_MSH_PARENTS];

  while (cgns > VALID_PTR)
    {
      /* Get tet index */

      i = (long)RDB[cgns + IFC_TET_PRNT_IDX];

      /* Read temperature or use uniform value */

      if (disttype == OF_INTERNAL_FIELD_UNIFORM)
        T = value;
      else
        {
          line = ReadOFData(fp, OF_FILE_TEMP);

          if (sscanf(line, "%lf", &T) == EOF)
            Die(FUNCTION_NAME, "Not enough entries in temp. file");
        }

      /* Check type flag */

      if (type == 2)
        T = T + T0;
      else if (type != 1)
        Die(FUNCTION_NAME, "Invalid temperature type flag");

      /* Compare to maximum and minimum */

      if (T > Tmax)
        Tmax = T;
      if ((T > 0.0) && (T < Tmin))
        Tmin = T;

      /**********************************************************/
      /* Convergence criterions based on momentary distribution */
      /**********************************************************/

      if (update)
        {
          old = RDB[tmplist + i];
          new = T;

          CalcConvCriteria(new, old,
                           &maxdiff, &maxeps, &L2abs, &L2rel);
        }

      /* Put value */

      WDB[tmplist + i] = T;

      /* Next cell */

      cgns = NextItem(cgns);
    }

  /* Close file */

  fclose(fp);

  /* Put maximum and minimum temperature */

  if (update == NO)
    {
      WDB[ifc + IFC_MAX_TEMP] = Tmax;
      WDB[ifc + IFC_MIN_TEMP] = Tmin;

      /* Set TMS on */

      WDB[DATA_TMS_MODE] = (double)TMS_MODE_CE;
    }

  /* Print convergence criteria if updating */

  if (update == YES)
    {
      /* Open output file for convergence */

      if (WDB[DATA_SOL_REL_ITER] == (double)0)
        {
          /* On first iteration open a new file */

          sprintf(tmpstr, "%s_Tconv%ld.m", GetText(ifc + IFC_PTR_INPUT_FNAME),
                  (long)RDB[DATA_BURN_STEP]);

          fp = fopen(tmpstr, "w");

          /* Reset idx */

          fprintf(fp,"\nidx = 1;\n\n");
        }
      else
        {
          /* On subsequent iterations append to old file */

          sprintf(tmpstr, "%s_Tconv%ld.m", GetText(ifc + IFC_PTR_INPUT_FNAME),
                  (long)RDB[DATA_BURN_STEP]);

          fp = fopen(tmpstr, "a");

          /* Increment idx */

          fprintf(fp,"\nidx = idx + 1;\n\n");
        }

      /* Write out convergence criteria of momentary temperature distribution */

      fprintf(fp, "T_eps(idx) = %E;\n", maxeps);
      fprintf(fp, "T_delta(idx) = %E;\n", maxdiff);
      fprintf(fp, "T_L2_of_absolute(idx) = %E;\n", sqrt(L2abs));
      fprintf(fp, "T_L2_of_relative(idx) = %E;\n", sqrt(L2rel));

      /* Close output file */

      fclose(fp);
    }
}
