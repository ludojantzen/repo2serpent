/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : finalizeumsh.c                                 */
/*                                                                           */
/* Created:       2017/08/03 (VVa)                                           */
/* Last modified: 2019/10/11 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Checks type of unstructured mesh and splits it if necessary. */
/*                                                                           */
/* Comments:   -                                                             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FinalizeUMSH:"

/*****************************************************************************/

void FinalizeUMSH(long ifc)
{
  long nhexc, ntetc, npyramc, nprismc, npolyc, sslist, nf, nc;
  long nfacepts[3], ptr, cgns, nd, i, n, np, surf;

  /*******************************************************************/

  /***** Error check and finalize the mesh ***************************/

  /* Reset number of different cell types */

  nhexc = 0;
  ntetc = 0;
  nprismc = 0;
  npolyc = 0;
  npyramc = 0;

  /* Get pointer to surface list */

  sslist = (long)RDB[ifc + IFC_PTR_SURF_LIST_PARENTS];

  /* Loop over cells */

  cgns = (long)RDB[ifc + IFC_PTR_TET_MSH];

  while (cgns > VALID_PTR)
    {

      /* Check boundaries */

      if (RDB[cgns + IFC_TET_PRNT_XMIN] < RDB[ifc + IFC_MESH_XMIN])
        Die(FUNCTION_NAME, "Error in boundaries");
      if (RDB[cgns + IFC_TET_PRNT_XMAX] > RDB[ifc + IFC_MESH_XMAX])
        Die(FUNCTION_NAME, "Error in boundaries");
      if (RDB[cgns + IFC_TET_PRNT_YMIN] < RDB[ifc + IFC_MESH_YMIN])
        Die(FUNCTION_NAME, "Error in boundaries");
      if (RDB[cgns + IFC_TET_PRNT_YMAX] > RDB[ifc + IFC_MESH_YMAX])
        Die(FUNCTION_NAME, "Error in boundaries");
      if (RDB[cgns + IFC_TET_PRNT_ZMIN] < RDB[ifc + IFC_MESH_ZMIN])
        Die(FUNCTION_NAME, "Error in boundaries");
      if (RDB[cgns + IFC_TET_PRNT_ZMAX] > RDB[ifc + IFC_MESH_ZMAX])
        Die(FUNCTION_NAME, "Error in boundaries");

      /* Get number of faces */

      nd = (long)RDB[cgns + IFC_TET_MSH_NF];

      /* Reset number of different face types */

      /* 3 point faces */
      nfacepts[0] = 0;
      /* 4 point faces */
      nfacepts[1] = 0;
      /* other faces */
      nfacepts[2] = 0;

      /* Get pointers to face and list */

      ptr = (long)WDB[cgns + IFC_TET_MSH_PTR_FACES];

      /* Loop over faces */

      for (i = 0; i < nd ; i++)
        {

          /* Get surface index */

          n = (long)RDB[ptr + i];

          surf = (long)RDB[sslist + n];

          np = (long)RDB[surf + UMSH_SURF_N_POINTS];

          /* Store number of points */

          if (np == 3)
            nfacepts[0]++;
          else if (np == 4)
            nfacepts[1]++;
          else
            nfacepts[2]++;

        }

      /* Check if polyhedral cell */

      /* Less than four or more than six faces */
      if ((nd < 4) || (nd > 6))
        npolyc++;
      /* Contains faces that have less than three or more than */
      /* four points */
      else if (nfacepts[2] > 0)
        npolyc++;
      /* Four sides with three points in each */
      else if ((nd == 4) && (nfacepts[0] == 4))
        ntetc++;
      /* Five sides with three points in 4, four points in 1 */
      else if ((nd == 5) && ((nfacepts[0] == 4) && (nfacepts[1] == 1)))
        npyramc++;
      /* Five sides with three points in 2, four points in 3 */
      else if ((nd == 5) && ((nfacepts[0] == 2) && (nfacepts[1] == 3)))
        nprismc++;
      /* Six sides with four points in each */
      else if ((nd == 6) && (nfacepts[1] == 6))
        nhexc++;
      else
        npolyc++;

      /* Next */

      cgns = NextItem(cgns);
    }

  /* Read batches */

  ReadOFPatches(ifc);

  /* Print out number of different cells */
  /*
  fprintf(outp, "Composition of mesh:\n\n");
  */
  fprintf(outp, " - Tetrahedrons:      %ld\n", ntetc);
  fprintf(outp, " - Pyramids:          %ld\n", npyramc);
  fprintf(outp, " - Prisms:            %ld\n", nprismc);
  fprintf(outp, " - Hexahedrons:       %ld\n", nhexc);
  fprintf(outp, " - Other polyhedrons: %ld\n", npolyc);

  /*****************************************************/

  /* Check initial mesh */

  fprintf(outp, "\nChecking initial mesh:\n\n");

  CheckPolyhedMesh(ifc);

  /* Convert mesh to a tetmesh */

  /* Reset tet-list for children */

  WDB[ifc + IFC_PTR_TET_MSH] = (double)NULLPTR;

   /* Reset point list for children */

  WDB[ifc + IFC_PTR_POINT_LIST] = (double)NULLPTR;

  /* Reset total number of cells */

  WDB[ifc + IFC_NC] = 0;

  if (npolyc == 0)
    {
      /******* Calculate number of new faces ********************/
      /* Each hex cell will have 6*2 external faces and 6 internal faces        */
      /* Each prism cell will have 3*2+2 external faces and 2(?) internal faces */
      /* Each pyramid cell will have 6 external faces and 1 internal face       */
      /* Each tet cell will have 4 external faces                               */
      /* NB; Some hex cells only have 4 internal faces                          */

      nf = nhexc*(6*2 + 6) + nprismc*(3*2 + 2 + 2) + npyramc*(6+1) +ntetc*4;

      /******* Calculate number of new cells *****/
      /* Each hex cell will have 5 or 6 subcells */
      /* Each prism cell will have 3 subcells    */
      /* Each pyramid cell will have 2 subcells  */
      /* Each tet cell will have 1 subcells      */

      nc = nhexc*6 + nprismc*3 + npyramc*2 + ntetc*1;

      fprintf(outp, "Dividing simple mesh to tetrahedrons...\n");
      FixHexMesh(ifc, nf, nc);
    }
  else
    {

      fprintf(outp, "\nDividing polyhedral mesh to tetrahedrons.\n");
      fprintf(outp, "Total mem use %f MB\n", RDB[DATA_TOTAL_BYTES]/MEGA);

      FixPolyhedMesh(ifc);

      fprintf(outp, "Divided mem use %f MB\n", RDB[DATA_TOTAL_BYTES]/MEGA);
    }
  /*
  fprintf(outp, "Checking divided mesh.\n");
  */
  /* New subroutine must be written for checking tet meshes */

  /*
    CheckTetMesh(ifc);
  */

  /* Switch type to regular tet mesh */

  WDB[ifc + IFC_TYPE] = (double)IFC_TYPE_TET_MESH;
}
