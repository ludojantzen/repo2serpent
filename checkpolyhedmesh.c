/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : checkpolyhedmesh.c                             */
/*                                                                           */
/* Created:       2015/01/12 (VVa)                                           */
/* Last modified: 2018/01/25 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Does some sanity checks to a divided polyhedral mesh         */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CheckPolyhedMesh:"

/*****************************************************************************/

void CheckPolyhedMesh(long ifc)
{
  long cgns, ownrs, nbrs, loc1, surf, surflist;
  long i, j, k, n, nsurf;
  long nc, nfo, nfn, np, fail, failnbr, failown, failownside, failnbrside;
  long ptr, pt;
  double x,y,z, ncalc;
  long prnt;

  /* Number of faces per cell and number of points per face is already */
  /* checked in fixpolyhedmesh.c                                       */

  /* Check neighbour relations */

  /* Get owner list */

  ownrs = (long)WDB[ifc + IFC_PTR_OWNR_LIST_PARENTS];
  CheckPointer(FUNCTION_NAME, "(ownrs)", DATA_ARRAY, ownrs);

  /* Get neighbour list */

  nbrs =  (long)WDB[ifc + IFC_PTR_NBR_LIST_PARENTS];
  CheckPointer(FUNCTION_NAME, "(nbrs)", DATA_ARRAY, nbrs);

  /* Get pointer to surface list */

  surflist = (long)RDB[ifc + IFC_PTR_SURF_LIST_PARENTS];
  CheckPointer(FUNCTION_NAME, "(surflist)", DATA_ARRAY, surflist);

  /* Get number of surfaces in mesh */

  nsurf = (long)RDB[ifc + IFC_NF_PARENTS];

  /* Loop over surfaces */

  failown = 0;
  failownside = 0;
  failnbr = 0;
  failnbrside = 0;

  for (n = 0; n < nsurf; n++)
    {
      /* Get pointer to n'th surface */

      surf = (long)RDB[surflist + n];

      /**********************/
      /* Check neighbourity */
      /**********************/

      /* Get pointer to owner cell */

      cgns = (long)RDB[ownrs + n];
      CheckPointer(FUNCTION_NAME, "(cgns)", DATA_ARRAY, cgns);

      /* Get pointer to face list */

      loc1 = (long)WDB[cgns + IFC_TET_MSH_PTR_FACES];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Get number of faces */

      nfo = (long)RDB[cgns + IFC_TET_MSH_NF];

      /* Loop over face list to find this surfaces index */

      for (i = 0; i < nfo; i++)
        if ((long)RDB[loc1 + i] == n)
          break;

      /* Check if this surface was not in owners face list */

      if (i == nfo)
        {
#ifdef DEBUG
          Warn(FUNCTION_NAME, "Could not find face from owner");
#endif

          /* Count to failed */

          failown++;

          /* Next surface */

          continue;
        }

      /* Get pointer to side list */

      loc1 = (long)WDB[cgns + IFC_TET_MSH_PTR_SIDES];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Check that side for owner is -1 */

      if ((long)RDB[loc1 + i] != -1)
        {
#ifdef DEBUG
          Warn(FUNCTION_NAME, "Side for owner is not -1");
#endif

          /* Count to failed */

          failownside++;

          /* Next surface */

          continue;
        }
      /* Get pointer to neighbour cell */

      cgns = (long)RDB[nbrs + n];

      /* If neighbour cell is -1 there is no neighbour */

      if(cgns == -1)
        {

          /* Handle next surface */

          continue;
        }

      CheckPointer(FUNCTION_NAME, "(cgns)", DATA_ARRAY, cgns);

      /* Get number of faces (neighbor) */

      nfn = (long)RDB[cgns + IFC_TET_MSH_NF];

      /* Get pointer to face list */

      loc1 = (long)WDB[cgns + IFC_TET_MSH_PTR_FACES];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Loop over face list to find this surfaces index */

      for (i = 0; i < nfn; i++)
        if ((long)RDB[loc1 + i] == n)
          break;

      if (i == nfn)
        {

#ifdef DEBUG
          fprintf(outp, "Surface index %ld not owned by owner and neighbour\n",
                  n);

          /* Get pointer to owner cell */

          cgns = (long)RDB[ownrs + n];
          CheckPointer(FUNCTION_NAME, "(cgns)", DATA_ARRAY, cgns);

          /* Get pointer to face list */

          loc1 = (long)WDB[cgns + IFC_TET_MSH_PTR_FACES];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          fprintf(outp, "Surfaces for owner:\n");
          for (i = 0; i < nfo; i++)
            fprintf(outp, "%ld\n", (long)RDB[loc1 + i]);


          cgns = (long)RDB[nbrs + n];
          CheckPointer(FUNCTION_NAME, "(cgns)", DATA_ARRAY, cgns);

          /* Get pointer to face list */

          loc1 = (long)WDB[cgns + IFC_TET_MSH_PTR_FACES];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          fprintf(outp, "Surfaces for neighbor:\n");
          for (i = 0; i < nfn; i++)
            fprintf(outp, "%ld\n", (long)RDB[loc1 + i]);

#endif

          /* Count to failed */

          failnbr++;

          /* Cycle loop */

          continue;

        }

      /* Get pointer to side list */

      loc1 = (long)WDB[cgns + IFC_TET_MSH_PTR_SIDES];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Check that side for neighbour is 1 */

      if ((long)RDB[loc1 + i] != 1)
        {
#ifdef DEBUG
          Warn(FUNCTION_NAME, "Side for neighbour is not 1");
#endif

          /* Count to failed */

          failnbrside++;

          /* Next surface */

          continue;
        }
    }

  /* Loop over cells */

  nc = 0;

  fail = 0;

  cgns = (long)RDB[ifc + IFC_PTR_TET_MSH];

  while (cgns > VALID_PTR)
    {
      /* Get pointer to parent or self */

      prnt = cgns;

      /* Get number of faces */

      nfo = (long)RDB[cgns + IFC_TET_MSH_NF];

      /* Get pointer to face list */

      loc1 = (long)WDB[cgns + IFC_TET_MSH_PTR_FACES];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Loop over faces to check centerpoints */

      x = 0.0;
      y = 0.0;
      z = 0.0;

      ncalc = 1.0;

      for (j = 0; j < nfo; j++)
        {

          /* Get index of face */

          n = (long)RDB[loc1 + j];

          /* Get pointer to surface */

          surf = (long)RDB[surflist + n];
          CheckPointer(FUNCTION_NAME, "(surf 1)", DATA_ARRAY, surf);

          /* Get number of points on the face */

          np = (long)RDB[surf + UMSH_SURF_N_POINTS];

          /* Get pointer to surface parameters */

          ptr = (long)RDB[surf + UMSH_SURF_PTR_POINTS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Calculate face centerpoint */

          for (k = 0; k < np; k++)
            {

              /* Get pointer to beginning of point */

              pt = (long)RDB[ptr + k];

              /* Loop over xyz and average for cell centerpoint */

              x = x*(ncalc - 1)/ncalc + RDB[pt + 0]/ncalc;
              y = y*(ncalc - 1)/ncalc + RDB[pt + 1]/ncalc;
              z = z*(ncalc - 1)/ncalc + RDB[pt + 2]/ncalc;

              ncalc++;
            }

        }

      if (!InUMSHCell(ifc, cgns, x, y, z, YES))
        {

#ifdef DEBUG

          /* Get number of faces */

          nfo = (long)RDB[cgns + IFC_TET_MSH_NF];

          /* Get pointer to face list */

          loc1 = (long)WDB[cgns + IFC_TET_MSH_PTR_FACES];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          /* Loop over faces and print points */

          ncalc = 1.0;

          for (j = 0; j < nfo; j++)
            {

              /* Get index of face */

              n = (long)RDB[loc1 + j];

              fprintf(outp, "F%ld = [ ", j+1);

              /* Get pointer to surface */

              surf = ListPtr(surflist, n);
              CheckPointer(FUNCTION_NAME, "(surf 2)", DATA_ARRAY, surf);

              /* Get number of points on the face */

              np = (long)RDB[surf + SURFACE_N_PARAMS];

              /* Get pointer to surface parameters */

              ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Print out points */

              for (k = 0; k < np; k++)
                {

                  /* Get pointer to beginning of point */

                  pt = (long)RDB[ptr + k];

                  fprintf(outp, "%E %E %E; ",
                          RDB[pt + 0],
                          RDB[pt + 1],
                          RDB[pt + 2]);

                }

              fprintf(outp, "]\n");

            }

          fprintf(outp, "figure();\n");
          fprintf(outp, "hold on;\n");

          /* Get pointer to side list */

          ptr = (long)WDB[cgns + IFC_TET_MSH_PTR_SIDES];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          for (j = 0; j < nfo; j++)
            {

              /* Get index of face */

              n = (long)RDB[loc1 + j];

              /* Get pointer to surface */

              surf = ListPtr(surflist, n);
              CheckPointer(FUNCTION_NAME, "(surf 3)", DATA_ARRAY, surf);

              /* Get number of points on the face */

              np = (long)RDB[surf + SURFACE_N_PARAMS];

              fprintf(outp,
                      "fill3(F%ld(:,1),F%ld(:,2),F%ld(:,3),ones(%ld,1)*%f);\n",
                      j+1, j+1, j+1, np, (1-j*1.0/(2.0*(double)nfo))*RDB[ptr + j]);
            }

          fprintf(outp, "hold off;\n");

          fprintf(outp, "title('Cell %ld')\n",
                  (long)RDB[prnt + IFC_TET_PRNT_IDX]);

          fprintf(outp, "Barycenter (%E %E %E)\n", x, y, z);
          Warn(FUNCTION_NAME, "Barycenter not inside cell %ld", (long)RDB[prnt + IFC_TET_PRNT_IDX]);
#endif
          fail++;
        }

      /* Increment number of cells tested */

      nc++;

      /* Next cell */

      cgns = NextItem(cgns);
    }

  fprintf(outp, " - Checked neighbority for %ld surfaces:\n",
          nsurf);
  fprintf(outp, "   %ld fails due to face not in owner cells face list.\n", failown);
  fprintf(outp, "   %ld fails due to side of face not being -1 for owner cell.\n", failownside);
  fprintf(outp, "   %ld fails due to face not in neighbor cells face list.\n", failnbr);
  fprintf(outp, "   %ld fails due to side of face not being +1 for neighbor cell.\n", failnbrside);

  fprintf(outp, " - Checked that barycenter is inside cell for %ld cells, %ld failed\n", nc, fail);

  /* Die if there were errors */

  if (failnbr + fail > 0)
    Die(FUNCTION_NAME, "There were errors in the mesh");

  fprintf(outp, "\n");
}

/*****************************************************************************/
