/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : calculatetetcenter.c                           */
/*                                                                           */
/* Created:       2015/02/19 (VVa)                                           */
/* Last modified: 2018/01/25 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Calculates cell centerpoints for all tet-mesh cells in list  */
/*              and face centerpoints                                        */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CalculateTetCenter:"

/*****************************************************************************/

void CalculateTetCenter(long cgns, long surflist, long cellpts, long facepts)
{
  long loc1, ptr, surf, pt;
  long nf, np, i, j, k, n;
  double p1[3], p2[3], nc;

  /* Get pointer to first cgns cell */

  cgns = FirstItem(cgns);
  CheckPointer(FUNCTION_NAME, "cgns", DATA_ARRAY, cgns);

  /* Reset index */

  i = 0;

  /* Loop over cells to create centerpoints */

  while (cgns > VALID_PTR)
    {

      /* Get number of faces */

      nf = (long)RDB[cgns + IFC_TET_MSH_NF];

      /* Get pointer to face list */

      loc1 = (long)WDB[cgns + IFC_TET_MSH_PTR_FACES];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Loop over faces to calculate and store centerpoints */

      p2[0] = 0.0;
      p2[1] = 0.0;
      p2[2] = 0.0;

      /* Reset number of averaged points */

      nc = 1.0;

      for (j = 0; j < nf; j++)
        {

          /* Get index of face */

          n = (long)RDB[loc1 + j];

          /* Get pointer to surface */

          surf = (long)RDB[surflist + n];
          CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

          /* Get pointer to surface parameters */

          ptr = (long)RDB[surf + UMSH_SURF_PTR_POINTS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Get number of points on the face */

          np = (long)RDB[surf + UMSH_SURF_N_POINTS];

          /* Reset face centerpoint */

          p1[0] = 0.0;
          p1[1] = 0.0;
          p1[2] = 0.0;

          /* Calculate face centerpoint */

          for (k = 0; k < np; k++)
            {

              /* Get pointer to beginning of point */

              pt = (long)RDB[ptr + k];

              /* Loop over xyz and average for cell centerpoint */

              p2[0] = p2[0]*(nc-1.0)/nc + RDB[pt + 0]/nc;
              p2[1] = p2[1]*(nc-1.0)/nc + RDB[pt + 1]/nc;
              p2[2] = p2[2]*(nc-1.0)/nc + RDB[pt + 2]/nc;

              /* Loop over xyz and average for face centerpoint */

              p1[0] += RDB[pt + 0]/(double)np;
              p1[1] += RDB[pt + 1]/(double)np;
              p1[2] += RDB[pt + 2]/(double)np;

              /* Increment number of points stored for cell centerpoint */

              nc++;

            }

          /* Put face centerpoint if requested */

          if (facepts > VALID_PTR)
            {
              WDB[facepts + n*3 + 0] = p1[0];
              WDB[facepts + n*3 + 1] = p1[1];
              WDB[facepts + n*3 + 2] = p1[2];

              /* printf("Face point %E %E %E\n", p1[0], p1[1], p1[2]);*/
            }

        }

      /* Put cell centerpoint if requested */

      if (cellpts > VALID_PTR)
        {
          WDB[cellpts + i*3 + 0] = p2[0];
          WDB[cellpts + i*3 + 1] = p2[1];
          WDB[cellpts + i*3 + 2] = p2[2];
        }
      /*   printf("Cell point %E %E %E\n", p2[0], p2[1], p2[2]);*/
      /* Next cell */

      cgns = NextItem(cgns);
      i++;
    }

}

/*****************************************************************************/
