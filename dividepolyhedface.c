/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : dividepolyhedface.c                            */
/*                                                                           */
/* Created:       2014/02/24 (VVa)                                           */
/* Last modified: 2017/11/29 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Creates subcells on a polyhedral cell-face by adding points  */
/*              to face centerpoints as well as to the centerpoint of        */
/*              face centerpoints ("cell centerpoint")                       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DividePolyhedFace:"

/*****************************************************************************/

void DividePolyhedFace(long ifc, long cidx, long fidx, long side, long np,
                       long *newcells)
{
  long p0, p1, p2, p3;
  long sparams;
  long surflist, cgns;
  long i, j, ptr;
  long *perimeter;
  long surf, tet;
  long tet1, nbridx;
  long p1a, p2a, p3a, p0b, p2b, p3b;
  double u, u0, u1, v, v0, v1, w, w0, w1;

  /* Get pointer to parent list */

  cgns = (long)RDB[ifc + IFC_PTR_TET_MSH_PARENTS];
  CheckPointer(FUNCTION_NAME, "(cgnslist)", DATA_ARRAY, cgns);

  /* Get pointer to cell to be divided */

  cgns = ListPtr(cgns, cidx);
  CheckPointer(FUNCTION_NAME, "(cgns)", DATA_ARRAY, cgns);

  /* Get face surface */
  /* Get pointer to interface surface list */

  surflist = (long)RDB[ifc + IFC_PTR_SURF_LIST_PARENTS];
  CheckPointer(FUNCTION_NAME, "(surflist)", DATA_ARRAY, surflist);

  /* Get pointer to surface of this face */

  surf = (long)RDB[surflist + fidx];
  CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

  /* Get pointer of surface parameters (points) of old face */

  sparams = (long)RDB[surf + UMSH_SURF_PTR_POINTS];
  CheckPointer(FUNCTION_NAME, "(sparams)", DATA_ARRAY, sparams);

  /* Get old neighbour index */

  if(side == -1)
    {
      /* This cell owns this surface */

      ptr = (long)WDB[ifc + IFC_PTR_NBR_LIST_PARENTS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Pointer to neighbour (may be -1 if outside) */

      ptr = (long)RDB[ptr + fidx];

      /* Neighbour's index */

      if (ptr > VALID_PTR)
        {
          nbridx = (long)RDB[ptr + IFC_TET_PRNT_IDX];
        }
      else if (ptr < 0)
        {
          nbridx = -1;
        }
      else
        {
          nbridx = -1;
          Die(FUNCTION_NAME, "Index error");
        }
    }
  else
    {
      /* This cell is this surfaces neighbour */

      ptr = (long)WDB[ifc + IFC_PTR_OWNR_LIST_PARENTS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Pointer to owner of the surface (should exist) */

      ptr = (long)RDB[ptr + fidx];
      CheckPointer(FUNCTION_NAME, "(ptr other)", DATA_ARRAY, ptr);

      /* Neighbour's index */

      nbridx = (long)RDB[ptr + IFC_TET_PRNT_IDX];
    }

  /* Create np tets */

  for (i = 0; i < np; i++)
    {
      /* Create new tet cell */

      tet = ReallocMem(DATA_ARRAY, TET_BLOCK_SIZE);
      CheckPointer(FUNCTION_NAME, "(tet)", DATA_ARRAY, tet);

      /* Put pointer to parent */

      WDB[tet + TET_PTR_PARENT] = (double)cgns;

      /* Store the created cell to newcells */

      newcells[i] = tet;
    }

  /***************************************************************/
  /* Get the face perimeter (pointers to points) to an array     */
  /* in such a manner that cross product of two subsequent edges */
  /* will point inside of the parent cell, e.g. same ordering as */
  /* the bottom face of this hexahedron                          */
  /*   7-----6                                                   */
  /*  /|    /|                                                   */
  /* 4-----5 |                                                   */
  /* | 3---|-2                                                   */
  /* |/    |/                                                    */
  /* 0-----1                                                     */
  /***************************************************************/

  perimeter = (long *)Mem(MEM_ALLOC, np, sizeof(long));

  if (side == -1)
    {
      /* This cell owns this surface */

      j = 0;
      for (i = np-1; i >= 0; i--)
        perimeter[j++] = (long)RDB[sparams + i];
    }
  else
    {
      /* This cell does not own this surface */

      for (i = 0; i < np; i++)
        perimeter[i] = (long)RDB[sparams + i];
    }

  /* Get pointer to face centerpoint list */

  ptr = (long)RDB[ifc + IFC_PTR_FACE_CP_LIST_PARENTS];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Create face centerpoint pointer */

  p2 = ptr + fidx*3;

  /* Get pointer to cell centerpoint list */

  ptr = (long)RDB[ifc + IFC_PTR_CELL_CP_LIST_PARENTS];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Create cell centerpoint pointer */

  p3 = ptr + cidx*3;

  /* Loop over face perimeter */

  for(i = 0; i < np; i++)
    {

      /* Get i:th new tet */

      tet = newcells[i];

      /* Get next two points from surface perimeter */

      p0 = perimeter[i];

      if(i < np - 1)
        p1 = perimeter[i + 1];
      else
        p1 = perimeter[0];

      /* Test that tetrahedron will be like */

      /*         3                      3            */
      /*        /|\                    /|\           */
      /*       / | \                  / | \          */
      /*  OK: 0 -|- 2        Not OK: 0 -|- 1         */
      /*       \ | /                  \ | /          */
      /*        \|/                    \|/           */
      /*         1                      2            */
      /*                                             */

      u0 = RDB[p1 + 0] - RDB[p0 + 0];
      v0 = RDB[p1 + 1] - RDB[p0 + 1];
      w0 = RDB[p1 + 2] - RDB[p0 + 2];

      u1 = RDB[p2 + 0] - RDB[p1 + 0];
      v1 = RDB[p2 + 1] - RDB[p1 + 1];
      w1 = RDB[p2 + 2] - RDB[p1 + 2];

      /* From p0 to p3 */

      u = RDB[p3 + 0] - RDB[p0 + 0];
      v = RDB[p3 + 1] - RDB[p0 + 1];
      w = RDB[p3 + 2] - RDB[p0 + 2];

      /* Check scalar triple product */

      if (u*(v0*w1-w0*v1) + v*(w0*u1-u0*w1) + w*(u0*v1-v0*u1) <= 0.0)
        Die(FUNCTION_NAME, "Not OK %f", u*(v0*w1-w0*v1) + v*(w0*u1-u0*w1) + w*(u0*v1-v0*u1));

      /* Put points to Tet */

      WDB[tet + TET_POINTS + 0] = (double)p0;
      WDB[tet + TET_POINTS + 1] = (double)p1;
      WDB[tet + TET_POINTS + 2] = (double)p2;
      WDB[tet + TET_POINTS + 3] = (double)p3;

      /***********************************************************/
      /* Numbering of tet faces:                                 */
      /* First face  (0,2,1) is out of cell (out of face).       */
      /* Second face (1,2,3) is forward on old face perimeter    */
      /* Third face  (0,3,2) is backward on old face perimeter   */
      /* Fourth face (0,1,3) is inside cell but not on this face */
      /***********************************************************/

      /* Put first face neighbour as the neighbour parent cell */

      WDB[tet + TET_NEIGHBOURS + 0] = (double)nbridx;

    }

  /* Set internal neighbours for this face (loop) */

  for (i = 0; i < np; i++)
    {

      tet = newcells[i];

      /* Get next cell */

      if(i < np - 1)
        tet1 = newcells[i+1];
      else
        tet1 = newcells[0];

      /* tet1 is the second neighbour of tet */
      /* tet  is the third neighbour of tet1 */

      WDB[tet  + TET_NEIGHBOURS + 1] = (double)tet1;
      WDB[tet1 + TET_NEIGHBOURS + 2] = (double)tet;

      /* Check that the faces match just to be sure */

      p1a = (long)RDB[tet + TET_POINTS + 1];
      p2a = (long)RDB[tet + TET_POINTS + 2];
      p3a = (long)RDB[tet + TET_POINTS + 3];

      p0b = (long)RDB[tet1 + TET_POINTS + 0];
      p2b = (long)RDB[tet1 + TET_POINTS + 2];
      p3b = (long)RDB[tet1 + TET_POINTS + 3];

      if (!((p1a == p0b) && (p2a == p2b) && (p3a == p3b)))
        Die(FUNCTION_NAME, "Forward and backward faces don't match %ld %ld.", p1a, p0b);
    }

  Mem(MEM_FREE, perimeter);
}

/*****************************************************************************/
