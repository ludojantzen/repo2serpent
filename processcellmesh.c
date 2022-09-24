/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processcellmesh.c                              */
/*                                                                           */
/* Created:       2013/10/15 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Creates super-imposed cell search meshes                     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessCellMesh:"

/*****************************************************************************/

void ProcessCellMesh()
{
  long loc0, loc1, uni, type, n0, n1, n2, n, msh, lst, new, ptr, sz;
  double params[6];

  /* Check pointer */

  if ((long)RDB[DATA_PTR_CSM0] < VALID_PTR)
    return;

  fprintf(outp, "Creating super-imposed cell search meshes...\n");

  /* Loop over structrues */

  loc0 = (long)RDB[DATA_PTR_CSM0];
  while (loc0 > VALID_PTR)
    {
      /* Find universe */
      
      uni = (long)RDB[DATA_PTR_U0];
      while (uni > VALID_PTR)
        {
          /* Compare name */

          if (CompareStr(uni + UNIVERSE_PTR_NAME, loc0 + CELL_MESH_PTR_UNI))
            break;
          
          /* Next universe */

          uni = NextItem(uni);
        }

      /* Check */

      if (uni < VALID_PTR)
        Error(0, "Universe %s in cell search mesh does not exist",
              GetText(loc0 + CELL_MESH_PTR_UNI));

      /* Put pointers */

      WDB[loc0 + CELL_MESH_PTR_UNI] = (double)uni;
      WDB[uni + UNIVERSE_PTR_CELL_MESH] = (double)loc0;

      /* Get data */

      type = (long)RDB[loc0 + CELL_MESH_TYPE];
      n0 = (long)RDB[loc0 + CELL_MESH_N0];
      n1 = (long)RDB[loc0 + CELL_MESH_N1];
      n2 = (long)RDB[loc0 + CELL_MESH_N2];
      params[0] = RDB[loc0 + CELL_MESH_MIN0];
      params[1] = RDB[loc0 + CELL_MESH_MAX0];
      params[2] = RDB[loc0 + CELL_MESH_MIN1]*PI/180.0;
      params[3] = RDB[loc0 + CELL_MESH_MAX1]*PI/180.0;
      params[4] = RDB[loc0 + CELL_MESH_MIN2];
      params[5] = RDB[loc0 + CELL_MESH_MAX2];

      /* Create mesh */

      msh = CreateMesh(type, MESH_CONTENT_PTR, -1, n0, n1, n2, params, -1);
      WDB[loc0 + CELL_MESH_PTR_MESH] = (double)msh;

      /* Pre-allocate memory */

      lst = (long)RDB[uni + UNIVERSE_PTR_CELL_LIST];
      CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);
      sz = ListSize(lst);

      PreallocMem(sz*n0*n1*n2*CELL_LIST_BLOCK_SIZE, DATA_ARRAY);
      PreallocMem(sz*n0*n1*n2, PRIVA_ARRAY);
      
      /* Get pointer to content */

      loc1 = (long)RDB[msh + MESH_PTR_PTR];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Loop over mesh */

      for (n = 0; n < n0*n1*n2; n++)
        {
          /* Loop over universe cell list */
        
          lst = (long)RDB[uni + UNIVERSE_PTR_CELL_LIST];
          while (lst > VALID_PTR)
            {
              /* Allocate memory */

              new = NewItem(loc1, CELL_LIST_BLOCK_SIZE);
              
              /* Copy data */

              memcpy(&WDB[new + LIST_DATA_SIZE], &RDB[lst + LIST_DATA_SIZE],
                     (CELL_LIST_BLOCK_SIZE - LIST_DATA_SIZE)*sizeof(double));
              
              /* Allocate memory for counter */

              ptr = AllocPrivateData(1, PRIVA_ARRAY);
              WDB[new + CELL_LIST_PTR_COUNT] = (double)ptr;

              /* Next */

              lst = NextItem(lst);
            }          

          /* Close list */

          ptr = (long)RDB[loc1];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          CloseList(ptr);

          /* Next */

          loc1++;
        }

      /* Pointer to next */

      loc0 = NextItem(loc0);
    }

  fprintf(outp, "OK.\n\n");
}

/*****************************************************************************/
