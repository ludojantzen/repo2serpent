/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readofdensities.c                              */
/*                                                                           */
/* Created:       2018/01/26 (VVa)                                           */
/* Last modified: 2018/08/13 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Reads OpenFOAM multi-physics interface density file.         */
/*                                                                           */
/* Comments:   -                                                             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadOFDensities:"

/*****************************************************************************/

void ReadOFDensities(long ifc, long update)
{
  char *line, tmpstr[MAX_STR], fname[MAX_STR];
  long n, nc, dim[7], i, dflist, cgns, disttype, mat, cell;
  double mul, d, dmax, dmin, d0, old, new, value;
  double maxeps, maxdiff, L2abs, L2rel, matdens;
  FILE *fp;

  /* Get density file filename */

  sprintf(fname, "%s", GetText(ifc + IFC_PTR_OF_RFILE));

  /* Open density file for reading */

  if ((fp = fopen(fname, "r")) == NULL)
    Error(ifc, "Density file \"%s\" does not exist", fname);

  /* Read header data */

  ReadOFHeader(fp, &n, &i, (long *)dim, &disttype, &value);

  /* Get number of cells */

  nc = (long)RDB[ifc + IFC_NC_PARENTS];

  /* Reset multiplier */

  mul = 1.0;

  /* Get nominal density */

  d0 = RDB[ifc + IFC_NOMINAL_DENSITY];

  /* Check units */

  if ((dim[0] == 0) && (dim[1] == 0))
    {
      /* Relative values, multiply by nominal density */

      mul = d0;
    }
  else if ((dim[0] == 1) && (dim[1] == -3))
    {
      /* kg/m3, convert to g/cm3 */

      mul = -0.001;
    }
  else
    Die(FUNCTION_NAME, "Undefined dimensions %ld/%ld", dim[0], dim[1]);

  /* Check type */
  /*
    if (n != OF_FILE_DENSITY)
    Note(ifc, "Possibly invalid format for density file");
  */
  /* Check size */

  if ((i != nc) && (disttype != OF_INTERNAL_FIELD_UNIFORM))
    Error(ifc, "Invalid number of entries in density file");

  /* Create density factor list or get pointer from memory */

  if (update == NO)
    {
      /* Allocate memory for list */

      dflist = ReallocMem(DATA_ARRAY, nc);

      /* Store pointer to list */

      WDB[ifc + IFC_PTR_DF_LIST] = (double)dflist;
    }
  else
    {
      /* Get pointer from memory */

      dflist = (long)RDB[ifc + IFC_PTR_DF_LIST];
    }

  /* Reset maximum of convergence criterion */

  maxeps  = 0.0;
  maxdiff = 0.0;
  L2abs   = 0.0;
  L2rel   = 0.0;

  /* Reset maximum and minimum density */

  dmax = 0.0;
  dmin = INFTY;

  /* Loop over parent cells */

  cgns = (long)RDB[ifc + IFC_PTR_TET_MSH_PARENTS];

  while (cgns > VALID_PTR)
    {
      /* Read density or use uniform value */

      if (disttype == OF_INTERNAL_FIELD_UNIFORM)
        d = value;
      else
        {

          line = ReadOFData(fp, OF_FILE_DENSITY);

          if (sscanf(line, "%lf", &d) == EOF)
            Die(FUNCTION_NAME, "Not enough entries in density file");
        }

      /* Convert to g/cm3 */

      d = d*mul;

      /* Compare to maximum and minimum */

      if (fabs(d) > fabs(dmax))
        dmax = d;
      if (fabs(d) < fabs(dmin))
        dmin = d;

      /* Get tet index */

      i = (long)RDB[cgns + IFC_TET_PRNT_IDX];

      /**********************************************************/
      /* Convergence criterions based on momentary distribution */
      /**********************************************************/

      if (update)
        {
          /* Get old density factor */

          old = RDB[dflist + i];

          /* Get geometry cell */

          cell = (long)RDB[cgns + IFC_TET_PRNT_PTR_CELL];
          CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

          /* Get material */

          mat = (long)RDB[cell + CELL_PTR_MAT];
          CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

          /* Check whether atomic or mass densities are given in interface */

          if (RDB[ifc + IFC_MAX_DENSITY] > 0)
            {
              /* Atomic densities given */

              matdens = RDB[mat + MATERIAL_ADENS];
            }
          else
            {
              /* Mass densities given */

              matdens = -RDB[mat + MATERIAL_MDENS];
            }

          /* Calculate old density */

          old = matdens*old;

          /* Put new density into a separate variable */

          new = d;

          /* Call the function to update convergence criteria */

          CalcConvCriteria(new, old,
                           &maxdiff, &maxeps, &L2abs, &L2rel);
        }

      /* Put value */

      WDB[dflist + i] = d;

      /* Next cell */

      cgns = NextItem(cgns);
    }

  /* Close file */

  fclose(fp);

  /* Put maximum  density (minimum not needed for majorant) */

  if (update == NO)
    WDB[ifc + IFC_MAX_DENSITY] = dmax;

  /* Print convergence criteria if updating */

  if (update == YES)
    {
      /* Open output file for convergence */

      if (WDB[DATA_SOL_REL_ITER] == (double)0)
        {
          /* On first iteration open a new file */

          sprintf(tmpstr, "%s_rhoconv%ld.m", GetText(ifc + IFC_PTR_INPUT_FNAME),
                  (long)RDB[DATA_BURN_STEP]);

          fp = fopen(tmpstr, "w");

          /* Reset idx */

          fprintf(fp,"\nidx = 1;\n\n");
        }
      else
        {
          /* On subsequent iterations append to old file */

          sprintf(tmpstr, "%s_rhoconv%ld.m", GetText(ifc + IFC_PTR_INPUT_FNAME),
                  (long)RDB[DATA_BURN_STEP]);

          fp = fopen(tmpstr, "a");

          /* Increment idx */

          fprintf(fp,"\nidx = idx + 1;\n\n");
        }

      /* Write out convergence criteria of momentary density distribution */

      fprintf(fp, "rho_eps(idx) = %E;\n", maxeps);
      fprintf(fp, "rho_delta(idx) = %E;\n", maxdiff);
      fprintf(fp, "rho_L2_of_absolute(idx) = %E;\n", sqrt(L2abs));
      fprintf(fp, "rho_L2_of_relative(idx) = %E;\n", sqrt(L2rel));

      /* Close output file */

      fclose(fp);
    }
}
