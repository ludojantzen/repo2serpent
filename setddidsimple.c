/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : setddidsimple.c                                */
/*                                                                           */
/* Created:       2018/03/12 (JLe)                                           */
/* Last modified: 2019/04/04 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Sets MPI id's for decomposed materials when domain           */
/*              decomposition is in use, using the material indexes.         */
/*                                                                           */
/* Comments: - Used with domain decomposition.                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SetDDIDSimple:"

/*****************************************************************************/

void SetDDIDSimple()
{
  long mat, nmat, nmatrest, n, id, mode, sum, *nr, subn;
  double x, y, f, rmin, rmax, r, phi0, phi;

  /* Check if domain decomposition is in use */

  if ((long)RDB[DATA_DD_DECOMPOSE] == NO)
    return;

  /* Count number of decomposed materials */

  nmat = 0;

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check if divided and add counter */

      if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_NEW)
        nmat++;

      /* Reset ID */

       WDB[mat + MATERIAL_MPI_ID] = -1.0;

      /* Next material */

      mat = NextItem(mat);
    }

  /* Print */

  fprintf(outp, "Decomposing %ld divided materials into %d domains:\n\n",
          nmat, mpitasks);

  /* Number of sub-segments for calculation */

  subn = 100000;

  /* Allocate memory for temporary array */

  nr = (long *)Mem(MEM_ALLOC, subn, sizeof(long));

  /* Check mode */

  if ((mode = (long)RDB[DATA_DD_MODE]) == DD_MODE_SIMPLE)
    {
      /***********************************************************************/

      /***** Simple mode based on indexing ***********************************/

      /* Calculate number of materials per task */

      nmatrest = nmat % mpitasks;
      nmat = nmat / mpitasks + 1;

      if (!nmatrest)
        {
          nmatrest = mpitasks;
          nmat--;
        }

      /* Reset index and count */

      id = 0;
      n = 0;

      /* Loop over materials */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Check if divided */

          if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_NEW)
            {
              /* Set id */

              WDB[mat + MATERIAL_MPI_ID] = (double)id;

              /* Add count */

              if (++n == nmat)
                {
                  /* Update id */

                  if (++id == nmatrest)
                    nmat--;

                  /* Reset count */

                  n = 0;
                }
            }

          /* Next material */

          mat = NextItem(mat);
        }

      /***********************************************************************/
    }
  else if ((mode == DD_MODE_SECTOR) ||
           ((mode == DD_MODE_AUTO) && (mpitasks < 5)))
    {
      /***********************************************************************/

      /***** Angular division only *******************************************/

      /* Reset counts */

      for (n = 0; n < subn; n++)
        nr[n] = 0;

      /* Loop over materials */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Check if divided and add counter */

          if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_NEW)
            {
              /* Get coordinates */

              x = RDB[mat + MATERIAL_DD_X0] - RDB[DATA_DD_ORIG_X0];
              y = RDB[mat + MATERIAL_DD_Y0] - RDB[DATA_DD_ORIG_Y0];

              /* Get angular sector */

              f = (PolarAngle(x, y) + RDB[DATA_DD_SECT0]/180.0*PI)/(2.0*PI);

              /* Truncate */

              f = f - (double)((long)f);

              /* Index */

              n = (long)(((double)subn - 1.0)*f);

              /* Add to count */

              nr[n]++;
            }

          /* Next material */

          mat = NextItem(mat);
        }

      /* Reset */

      id = 0;
      sum = 0;

      /* Loop over angular segments */

      for (n = 0; n < subn; n++)
        {
          /* Calculate cumulative number materials in segment */

          sum = sum + nr[n];

          /* Put id */

          nr[n] = id;

          /* Compare to number per domain */

          if (sum > (long)((double)nmat)/((double)mpitasks) - 1)
            {
              /* Reset sum */

              sum = 0;

              /* Update id */

              if (id < mpitasks - 1)
                id++;
            }
        }

      /* Loop over materials */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Check if divided and add counter */

          if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_NEW)
            {
              /* Get coordinates */

              x = RDB[mat + MATERIAL_DD_X0] - RDB[DATA_DD_ORIG_X0];
              y = RDB[mat + MATERIAL_DD_Y0] - RDB[DATA_DD_ORIG_Y0];

              /* Get angular sector */

              f = (PolarAngle(x, y) + RDB[DATA_DD_SECT0]/180.0*PI)/(2.0*PI);

              /* Truncate */

              f = f - (double)((long)f);

              /* Index */

              n = (long)(((double)subn - 1.0)*f);

              /* Check and put index */

              if ((nr[n] < 0) || (nr[n] > mpitasks - 1))
                Die(FUNCTION_NAME, "indexing error: %ld %ld", nr[n], n);
              else
                WDB[mat + MATERIAL_MPI_ID] = (double)nr[n];
            }

          /* Next material */

          mat = NextItem(mat);
        }

      /***********************************************************************/
    }
  else if (mode == DD_MODE_AUTO)
    {
      /***********************************************************************/

      /***** Radial / angular division ***************************************/

      /* Sector width */

      phi0 = 2.0*PI/((double)mpitasks - 1);

      /* Loop over materials and get minimum and maximum radius */

      rmin = INFTY;
      rmax = -INFTY;

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Check if divided and add counter */

          if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_NEW)
            {
              /* Get coordinates */

              x = RDB[mat + MATERIAL_DD_X0] - RDB[DATA_DD_ORIG_X0];
              y = RDB[mat + MATERIAL_DD_Y0] - RDB[DATA_DD_ORIG_Y0];

              /* Calculate radius */

              r = sqrt(x*x + y*y);

              /* Calculate translated angle */

              phi = PolarAngle(x, y);
              f = phi/phi0;
              phi = phi - floor(f)*phi0 - 0.5*phi0;

              /* Re-calculate x-coordinate */

              x = r*cos(phi);

              /* Compare to limits */

              if (x < rmin)
                rmin = x;

              if (x > rmax)
                rmax = x;
            }

          /* Next material */

          mat = NextItem(mat);
        }

      /* Check */

      if (rmin > rmax)
        Error(0, "Unable to perform radial division");

      /* Reset counts */

      for (n = 0; n < subn; n++)
        nr[n] = 0;

      /* Count number of particles in radial segments */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Check if divided and add counter */

          if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_NEW)
            {
              /* Get coordinates */

              x = RDB[mat + MATERIAL_DD_X0] - RDB[DATA_DD_ORIG_X0];
              y = RDB[mat + MATERIAL_DD_Y0] - RDB[DATA_DD_ORIG_Y0];

              /* Calculate radius */

              r = sqrt(x*x + y*y);

              /* Calculate translated angle */

              phi = PolarAngle(x, y);
              f = phi/phi0;
              phi = phi - floor(f)*phi0 - 0.5*phi0;

              /* Re-calculate x-coordinate */

              x = r*cos(phi);

              /* Calculate factor */

              f = (x - rmin)/(rmax - rmin);
              CheckValue(FUNCTION_NAME, "ff", "", f, 0.0, 1.0);

              /* Index */

              n = (long)(((double)subn - 1.0)*f);

              /* Add to count */

              nr[n]++;
            }

          /* Next material */

          mat = NextItem(mat);
        }

      /* Reset */

      id = 0;
      sum = 0;

      /* Loop over radial segments */

      for (n = 0; n < subn; n++)
        {
          /* Calculate cumulative number materials in segment */

          sum = sum + nr[n];

          /* Put id */

          nr[n] = id;

          /* Compare to number per domain */

          if (sum > (long)((double)nmat)/((double)mpitasks) - 1)
            {
              /* Reset sum */

              sum = 0;

              /* Update id */

              if (id < mpitasks - 1)
                id++;
            }
        }

      /* Reset remaining number of materials */

      nmatrest = nmat;

      /* Put indexes */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Check if divided and add counter */

          if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_NEW)
            {
              /* Get coordinates */

              x = RDB[mat + MATERIAL_DD_X0] - RDB[DATA_DD_ORIG_X0];
              y = RDB[mat + MATERIAL_DD_Y0] - RDB[DATA_DD_ORIG_Y0];

              /* Calculate radius */

              r = sqrt(x*x + y*y);

              /* Calculate translated angle */

              phi = PolarAngle(x, y);
              f = phi/phi0;
              phi = phi - floor(f)*phi0 - 0.5*phi0;

              /* Re-calculate x-coordinate */

              x = r*cos(phi);

              /* Calculate factor */

              f = (x - rmin)/(rmax - rmin);
              CheckValue(FUNCTION_NAME, "ff", "", f, 0.0, 1.0);

              /* Index */

              n = (long)(((double)subn - 1.0)*f);

              /* Put ID (only first zone) */

              if (nr[n] == 0)
                {
                  /* Check */

                  if ((nr[n] < 0) || (nr[n] > mpitasks - 1))
                    Die(FUNCTION_NAME, "indexing error: %ld %ld", nr[n], n);
                  else
                    WDB[mat + MATERIAL_MPI_ID] = (double)nr[n];

                  /* Subtract from remainng */

                  nmatrest--;
                }
            }

          /* Next material */

          mat = NextItem(mat);
        }

      /* Check remaining */

      if (nmatrest <= 0)
        Die(FUNCTION_NAME, "Error in DD");

      /* Reset counts */

      for (n = 0; n < subn; n++)
        nr[n] = 0;

      /* Loop over materials */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Check if divided and not ID'd and add counter */

          if (((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_NEW) &&
              ((long)RDB[mat + MATERIAL_MPI_ID] < 0.0))
            {
              /* Get coordinates */

              x = RDB[mat + MATERIAL_DD_X0] - RDB[DATA_DD_ORIG_X0];
              y = RDB[mat + MATERIAL_DD_Y0] - RDB[DATA_DD_ORIG_Y0];

              /* Get angular sector */

              f = (PolarAngle(x, y) + RDB[DATA_DD_SECT0]/180.0*PI)/(2.0*PI);

              /* Truncate */

              f = f - (double)((long)f);

              /* Index */

              n = (long)(((double)subn - 1.0)*f);

              /* Add to count */

              nr[n]++;
            }

          /* Next material */

          mat = NextItem(mat);
        }

      /* Reset */

      id = 0;
      sum = 0;

      /* Loop over angular segments */

      for (n = 0; n < subn; n++)
        {
          /* Calculate cumulative number materials in segment */

          sum = sum + nr[n];

          /* Put id */

          nr[n] = id + 1;

          /* Compare to number per domain */

          if (sum > (long)((double)nmatrest)/((double)(mpitasks - 1)) - 1)
            {
              /* Reset sum */

              sum = 0;

              /* Update id */

              if (id < mpitasks - 2)
                id++;
            }
        }

      /* Loop over materials */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Check if divided and add counter */

          if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_NEW)
            {
              /* Get coordinates */

              x = RDB[mat + MATERIAL_DD_X0] - RDB[DATA_DD_ORIG_X0];
              y = RDB[mat + MATERIAL_DD_Y0] - RDB[DATA_DD_ORIG_Y0];

              /* Get angular sector */

              f = (PolarAngle(x, y) + RDB[DATA_DD_SECT0]/180.0*PI)/(2.0*PI);

              /* Truncate */

              f = f - (double)((long)f);

              /* Index */

              n = (long)(((double)subn - 1.0)*f);

              /* Check and put index */

              if ((nr[n] < 0) || (nr[n] > mpitasks - 1))
                Die(FUNCTION_NAME, "indexing error: %ld %ld", nr[n], n);
              else if ((long)RDB[mat + MATERIAL_MPI_ID] < 0.0)
                WDB[mat + MATERIAL_MPI_ID] = (double)nr[n];
            }

          /* Next material */

          mat = NextItem(mat);
        }

      /***********************************************************************/
    }
  else
    Die(FUNCTION_NAME, "Invalid DD mode");

  /***************************************************************************/

  /***** Check number of materials per domain ********************************/

  /* Reset sum */

  sum = 0;

  /* Loop over domains */

  for (id = 0; id < mpitasks; id++)
    {
      /* Reset counter */

      n = 0;

      /* Loop over materials */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Check domain */

          if ((long)RDB[mat + MATERIAL_MPI_ID] == id)
            n++;

          /* Next material */

          mat = NextItem(mat);
        }

      /* Check */

      if (n == 0)
        Error(0, "Decomposition failed: no materials assigned to task %ld",
              id + 1);

      /* Add to sum */

      sum = sum + n;

      /* Print */

      fprintf(outp, "Domain %ld: %ld materials (%1.1f%%)\n", id + 1, n,
              100.0*((double)n/((double)nmat)));
    }

  /* Free allocated memory */

  Mem(MEM_FREE, nr);

  /* Check undivided */

  if (sum < nmat)
    Note(0, "%ld materials not decomposed", nmat - sum);
  else if ((sum > nmat) && (mode != DD_MODE_SIMPLE))
    Die(FUNCTION_NAME, "Not possible");

  /* Newline */

  fprintf(outp, "\n");

  /***************************************************************************/
}

/*****************************************************************************/
