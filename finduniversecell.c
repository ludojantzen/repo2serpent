/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : finduniversecell.c                             */
/*                                                                           */
/* Created:       2010/10/13 (JLe)                                           */
/* Last modified: 2018/06/26 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Find cell in universe                                        */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FindUniverseCell:"

/*****************************************************************************/

long FindUniverseCell(long uni, double x, double y, double z, long *ridx, 
                      long id)
{
  long cell, lst, ptr, msh, found, loc0, loc1, loc2, n, nc, opt, surf;
  

  /* Check universe pointer */

  CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

  /* Get search list option */

  ptr = (long)RDB[DATA_CELL_SEARCH_LIST];
  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
  opt = (long)GetPrivateData(ptr, id);

  /* Reset cell pointer */

  cell = -1;

  /***************************************************************************/

  /***** Check for surface crossing from ST **********************************/

  if (opt == CELL_SEARCH_LIST_SURF)
    {
      Die(FUNCTION_NAME, "This should not be used");

      /* Get pointer to nearest surface */

      ptr = (long)RDB[uni + UNIVERSE_PTR_NEAREST_SURF];
      CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);        

      /* Check surface pointer */

      if ((surf = (long)GetPrivateData(ptr, id)) > VALID_PTR)
        {
          /* Get pointer to search list */

          ptr = (long)RDB[surf + SURFACE_PTR_CELL_LIST];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          
          /* Loop over list */
          
          n = 0;
          while ((loc0 = ListPtr(ptr, n++)) > VALID_PTR)
            {
              /* Pointer to cell */
              
              loc1 = (long)RDB[loc0 + CELL_LIST_PTR_CELL];
              CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
              
              /* Check universe */
              
              if ((long)RDB[loc1 + CELL_PTR_UNI] != uni)
                continue;
              
              /* Test cell */
              
              if (InCell(loc1, x, y, z, NO, id) == YES)
                {
                  /* Put region index */
                  
                  *ridx = (long)RDB[loc0 + CELL_LIST_REG_IDX];
                  
                  /* Pointer to search list count */
                  
                  loc2 = (long)RDB[loc0 + CELL_LIST_PTR_COUNT];
                  CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);
                  
                  /* Add counter */
                  
                  AddPrivateData(loc2, 1, id);
                  
                  /* Put previous pointer */
                  
                  ptr = (long)RDB[uni + UNIVERSE_PTR_PREV_REG];
                  PutPrivateData(ptr, loc0, id);
                  
                  /* Return cell pointer */
                  
                  return loc1;
                }
            }
          
          /* Cell not found */

          Warn(FUNCTION_NAME, "Cell not found in surface neighbours");
        }
      
      /* Override option */

      opt = CELL_SEARCH_LIST_CELL;
    }

  /***************************************************************************/

  /***** Check for previous point ********************************************/

  /* Check plotter mode */
  
  if (((long)RDB[DATA_PLOTTER_MODE] == NO) || 
      ((long)RDB[DATA_QUICK_PLOT_MODE] == YES))
    {
      /* Check if source point search mode */

      if (opt == CELL_SEARCH_LIST_SRC)
        ptr = (long)RDB[uni + UNIVERSE_PTR_SRC_REG];
      else
        ptr = (long)RDB[uni + UNIVERSE_PTR_PREV_REG];

      /* Check previous */

      if ((lst = (long)GetPrivateData(ptr, id)) > VALID_PTR)
        {
          /* Get cell pointer */
          
          cell = (long)RDB[lst + CELL_LIST_PTR_CELL];
          CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);      

          /* Test cell */
          
          if (InCell(cell, x, y, z, NO, id) == YES)
            {
              /* Put region index */
              
              *ridx = (long)RDB[lst + CELL_LIST_REG_IDX];

              /* Check plotter mode */

              if ((long)RDB[DATA_PLOTTER_MODE] == NO)
                {
                  /* Pointer to search list count */

                  loc2 = (long)RDB[lst + CELL_LIST_PTR_COUNT];
                  CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

                  /* Add counter */
                  
                  AddPrivateData(loc2, 1, id);
                }

              /* Return cell pointer */
              
              return cell;
            }
        }
    }

  /***************************************************************************/

  /***** Adaptive cell search list *******************************************/

  /* Check if mesh exists */

  if ((msh = (long)RDB[uni + UNIVERSE_PTR_SEARCH_MESH]) > VALID_PTR)
    {
      /* Get mesh pointer */

      if ((loc0 = MeshPtr(msh, x, y, z)) > VALID_PTR)
        {
          /* Pointer to pointer */
      
          loc0 = (long)RDB[loc0];
          CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

          /* Loop over candidates */

          nc = (long)RDB[loc0 + CELL_SEARCH_MESH_N];
          for (n = 0; n < nc; n++)      
            {
              /* Get pointer to cell */
              
              cell = (long)RDB[loc0 + CELL_SEARCH_MESH_PTR_C1 + n];
              CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);
              
              /* Test cell */
              
              if (InCell(cell, x, y, z, NO, id) == YES)
                {
                  /* Put region index */
                 
                  /* HUOM!!!!!! */                  
 
                  *ridx = 0;
                  
                  /* Return cell pointer */
                
                  return cell;
                }
            }
          if ((long)RDB[DATA_PLOTTER_MODE] == YES)
            return GEOM_ERROR_NO_CELL;
        }

      /* Reset cell pointer */

      cell = -1;
    }

  /***************************************************************************/

  /***** Select search list **************************************************/

  /* Default to universe search lists */

  if (opt == CELL_SEARCH_LIST_SRC)
    ptr = (long)RDB[uni + UNIVERSE_PTR_SRC_CELL_LIST];
  else
    ptr = (long)RDB[uni + UNIVERSE_PTR_CELL_LIST];

  /* Check for search mesh */
  
  if ((msh = (long)RDB[uni + UNIVERSE_PTR_CELL_MESH]) > VALID_PTR)
    {
      /* Pointer to mesh */

      msh = (long)RDB[msh + CELL_MESH_PTR_MESH];
      CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

      /* Get pointer from search mesh */

      if ((loc0 = MeshPtr(msh, x, y, z)) > VALID_PTR)
        ptr = (long)RDB[loc0];
    }
  else if ((cell > VALID_PTR) && (opt == CELL_SEARCH_LIST_CELL))
    {
      /* Get pointer from cell search list of previous cell */

      if ((loc0 = (long)RDB[cell + CELL_PTR_SEARCH_LIST]) > VALID_PTR)
        ptr = loc0;
    }

  /* Check pointer */

  if (ptr < VALID_PTR)
    Die(FUNCTION_NAME, "Pointer error");

  /***************************************************************************/
  
  /***** Find cell in list ***************************************************/
      
  /* Reset flag and cell pointer */

  found = NO;
  cell = GEOM_ERROR_NO_CELL;
  lst = -1;

  /* Loop over list */
  
  n = 0;
  while ((loc0 = ListPtr(ptr, n++)) > VALID_PTR)
    {
      /* Pointer to cell */

      loc1 = (long)RDB[loc0 + CELL_LIST_PTR_CELL];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Test cell */

      if (InCell(loc1, x, y, z, NO, id) == YES)
        {
          /* Put region index */

          *ridx = (long)RDB[loc0 + CELL_LIST_REG_IDX];

          /* Check plotter mode */

          if ((long)RDB[DATA_PLOTTER_MODE] == NO)
            {   
              /* Pointer to search list count */

              loc2 = (long)RDB[loc0 + CELL_LIST_PTR_COUNT];
              CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

              /* Add counter */

              AddPrivateData(loc2, 1, id);

              /* Put previous pointer */
                  
              if (opt == CELL_SEARCH_LIST_SRC)
                ptr = (long)RDB[uni + UNIVERSE_PTR_SRC_REG];
              else
                ptr = (long)RDB[uni + UNIVERSE_PTR_PREV_REG];

              PutPrivateData(ptr, loc0, id);
              
              /* Return cell pointer */

              return loc1;
            }
          else if (found == YES)
            {
              /* Duplicate definition of geometry region */

              return GEOM_ERROR_MULTIPLE_CELLS;
            }
          else
            {
              /* Put pointer */

              cell = loc1;
              lst = loc0;

              /* Put found flag */

              found = YES;

              /* Check quick plotter mode */

              if ((long)RDB[DATA_QUICK_PLOT_MODE] == YES)
                break;
            }
        }
    }

  /* Put previous pointer */
  
  if (opt == CELL_SEARCH_LIST_SRC)
    ptr = (long)RDB[uni + UNIVERSE_PTR_SRC_REG];
  else
    ptr = (long)RDB[uni + UNIVERSE_PTR_PREV_REG];

  PutPrivateData(ptr, lst, id);

  /* Return cell pointer */

  return cell;
}

/*****************************************************************************/
