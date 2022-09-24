/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : findtetcell.c                                  */
/*                                                                           */
/* Created:       2012/09/11 (JLe)                                           */
/* Last modified: 2019/04/03 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Finds tetrahedral mesh cell                                  */
/*                                                                           */
/* Comments: -NOTE: Return value is pointer to tet cell, not geometry cell.  */
/*            (tuonne neighbour cell -rakenteeseen tallennetaan tässä        */
/*             moodissa tet-cellin pointteri, muulloin geometriacellin       */
/*             pointteri.)                                                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FindTetCell:"

#define FAST_MODE

/*****************************************************************************/

long FindTetCell(long ifc, double x, double y, double z, long id)
{
#ifdef FAST_MODE
  long msh, lst, loc0, ptr, next, prev, i, pt0, pt1, pt2;
  double bb[6];
  long *face[4];
  long face1[3] = {0,2,1};
  long face2[3] = {1,2,3};
  long face3[3] = {0,3,2};
  long face4[3] = {0,1,3};
  long nc;
  double dx, dy, dz, facex, facey, facez, r2, min;

  /* Store point list for easy access (faceX are constant lists) */

  face[0] = face1;
  face[1] = face2;
  face[2] = face3;
  face[3] = face4;
#else
  long tetlist, tet, ntet, i, prnt;
#endif

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ifc)", DATA_ARRAY, ifc);

#ifdef FAST_MODE

  /***************************************************************************/

  /***** Try previous cell ***************************************************/

  /* Check previous cell */

  ptr = (long)RDB[ifc + IFC_PTR_PREV_CELL];
  CheckPointer(FUNCTION_NAME, "(ptr1)", PRIVA_ARRAY, ptr);

  if ((loc0 = (long)GetPrivateData(ptr, id)) > VALID_PTR)
    {

      /* Test cell */

      if (InTetCell(loc0, x, y, z, YES) == YES)
        {
          /* Return pointer */

          return loc0;
        }
    }

  /***************************************************************************/

  /***** Try neighbour search ************************************************/

  /* Check if previous cell was given (NOTE: Tää on nyt disabloitu kun  */
  /* se ei toimi kovin hyvin) */

  prev = -1;
  next = -1;

  if (loc0 > VALID_PTR)
    {
      /* Reset minimum distance */

      min = INFTY;

      /* Reset counter */

      nc = 0;

      /* Loop */

      while (1 != 2)
        {
          /* Check pointer */

          CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

          /* Test cell */

          if (InTetCell(loc0, x, y, z, YES) == YES)
            {
              /* Store pointer */

              ptr = (long)RDB[ifc + IFC_PTR_PREV_CELL];
              CheckPointer(FUNCTION_NAME, "(ptr2)", PRIVA_ARRAY, ptr);
              PutPrivateData(ptr, loc0, id);

              /* Return pointer */

              return loc0;
            }

          /* Loop over faces and see which face has its barycenter closest to */
          /* target coordinates */

          for (i = 0; i < 4; i++)
            {
              /* Get pointers to points */

              pt0 = (long)RDB[loc0 + TET_POINTS + face[i][0]];
              pt1 = (long)RDB[loc0 + TET_POINTS + face[i][1]];
              pt2 = (long)RDB[loc0 + TET_POINTS + face[i][2]];

              /* Calculate barycenter of face */

              facex = (RDB[pt0 + 0] + RDB[pt1 + 0] + RDB[pt2 + 0])/3.0;
              facey = (RDB[pt0 + 1] + RDB[pt1 + 1] + RDB[pt2 + 1])/3.0;
              facez = (RDB[pt0 + 2] + RDB[pt1 + 2] + RDB[pt2 + 2])/3.0;

              /* Calculate distance components */

              dx = (x - facex);
              dy = (y - facey);
              dz = (z - facez);

              /* Compare square distance to minimum */

              if ((r2 = dx*dx + dy*dy + dz*dz) == 0.0)
                Die(FUNCTION_NAME, "Error in vector lenght");
              else if (r2 < min)
                {
                  /* Get neighbour through current face */

                  ptr = (long)RDB[loc0 + TET_NEIGHBOURS + i];

                  /* Update pointer and minimum */

                  next = ptr;
                  min = r2;
                }

            }

          /* Check if minimum was found */

          if (next < VALID_PTR)
            {
              /* No neighbour cell (point may be outside), break loop */

              break;
            }
          else if (next == prev)
            {
              /* Would return to previous cell leading to infinite oscillation */

              break;
            }
          else
            {
              /* Store previous tet */

              prev = loc0;

              /* Go to next tet */

              loc0 = next;
            }

          /* Check maximum number of iterations */

          if (nc++ > 10)
            break;
        }
    }

  /***************************************************************************/

  /***** Try list search *****************************************************/

  /* Get pointer to search mesh */

  msh = (long)RDB[ifc + IFC_PTR_SEARCH_MESH_LIST];
  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

  /* Get pointer */

  if ((lst = MeshPtr(msh, x, y, z)) > VALID_PTR)
    lst = (long)RDB[lst];
  else
    return NULLPTR;

  /* Check pointer */

  if (lst < VALID_PTR)
    return NULLPTR;

  /* Loop over content */

  while (lst > VALID_PTR)
    {
      /* Pointer to tet */

      loc0 = (long)RDB[lst + SEARCH_MESH_CELL_CONTENT];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      /* Calculate bounding box */

      CalculateTetBoundingBox(loc0, bb);

      /* Check with bounding box */

      if ((x < bb[0]) ||
          (x > bb[1]) ||
          (y < bb[2]) ||
          (y > bb[3]) ||
          (z < bb[4]) ||
          (z > bb[5]))
        {
          /* Cannot be inside, pointer to next */

          lst = NextItem(lst);

          /* Cycle loop */

          continue;
        }

      /* Test cell */

      if (InTetCell(loc0, x, y, z, YES) == YES)
        {
          /* Increment count (collisions in this search mesh cell) only      */
          /* during transport and only by id 0 so that no priva array needs  */
          /* to be allocated (would increase the memory size of search mesh) */

          if (((long)RDB[DATA_PLOTTER_MODE] == NO) && (id == 0))
            WDB[lst + SEARCH_MESH_PTR_CELL_COUNT] = RDB[lst + SEARCH_MESH_PTR_CELL_COUNT] + 1.0;

          /* Store pointer */

          ptr = (long)RDB[ifc + IFC_PTR_PREV_CELL];
          CheckPointer(FUNCTION_NAME, "(ptr4)", PRIVA_ARRAY, ptr);
          PutPrivateData(ptr, loc0, id);

          /* Return pointer */

          return loc0;
        }

      /* Next */

      lst = NextItem(lst);

      /***********************************************************************/
    }

  /***************************************************************************/

#else

  /***************************************************************************/

  /***** Safe mode (test all cells) ******************************************/

  /* Loop over all tet cells */

  tetlist = (long)RDB[ifc + IFC_PTR_TET_LIST];
  CheckPointer(FUNCTION_NAME, "(tetlist)", DATA_ARRAY, tetlist);

  ntet = (long)RDB[ifc + IFC_NC];

  for (i = 0; i < ntet; i++)
    {
      /* Get pointer to tet */

      tet = (long)RDB[tetlist + i];
      CheckPointer(FUNCTION_NAME, "(tet)", DATA_ARRAY, tet);

      /* Get pointer to parent */

      prnt = (long)RDB[tet + TET_PTR_PARENT];
      CheckPointer(FUNCTION_NAME, "(prnt)", DATA_ARRAY, prnt);

      /* Check with bounding box */

      if ((x < RDB[prnt + IFC_TET_PRNT_XMIN]) ||
          (x > RDB[prnt + IFC_TET_PRNT_XMAX]) ||
          (y < RDB[prnt + IFC_TET_PRNT_YMIN]) ||
          (y > RDB[prnt + IFC_TET_PRNT_YMAX]) ||
          (z < RDB[prnt + IFC_TET_PRNT_ZMIN]) ||
          (z > RDB[prnt + IFC_TET_PRNT_ZMAX]))
        {
          /* Cannot be inside, cycle loop */

          continue;
        }

      /* Test */

      if (InTetCell(tet, x, y, z, YES) == YES)
        return tet;

      /* Next tet cell */
    }

  /***************************************************************************/

#endif

  /* Not in any, return null */

  return NULLPTR;
}

/*****************************************************************************/
