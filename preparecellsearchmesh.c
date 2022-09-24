/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : preparecellsearchmesh.c                        */
/*                                                                           */
/* Created:       2018/06/20 (JLe)                                           */
/* Last modified: 2018/06/26 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Prepares an adaptive cell search mesh for cell universes     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrepareCellSearchMesh:"

void SplitCell(long, long, long);

/*****************************************************************************/

void PrepareCellSearchMesh()
{
  long uni, loc0, loc1, loc2, ptr, n, i, j, k, l, idx, msh, szp, sz, last;
  long it, np, nc, cell, id, depth, split, sz0;
  unsigned long seed;
  double xmin, xmax, ymin, ymax, zmin, zmax, lims[6], x, y, z;

  /* Get number of points */
  return;
  fprintf(outp, "Preparing adaptive cell search lists:\n\n");

  /* Expand PRIVA, BUF and RES2 arrays for OpenMP parallel calculation */
  /* (this is needed to enable call to cell test routine) */

  ExpandPrivateArrays();

  /* Put mesh depth, size, split criterion and number of sampled points */

  depth = 20;
  sz0 = 5;
  sz = 2;
  split = MAX_CELL_SEARCH_MESH_SZ;
  np = 1000;

  /* Loop over universes */

  uni = (long)RDB[DATA_PTR_U0];
  while (uni > VALID_PTR)
    {
      /* Check type */

      if ((long)RDB[uni + UNIVERSE_TYPE] != UNIVERSE_TYPE_CELL)
        {
          /* Next universe */

          uni = NextItem(uni);

          /* Cycle loop */

          continue;
        }

      /* Check size of cell list */

      ptr = (long)RDB[uni + UNIVERSE_PTR_CELL_LIST];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      
      if (ListSize(ptr) <= split)
        {
          /* Next universe */

          uni = NextItem(uni);

          /* Cycle loop */

          continue;
        }

      /***********************************************************************/

      /***** Prepare adaptive search mesh ************************************/

      /* Get boundaries */
      
      xmin = RDB[uni + UNIVERSE_MINX];
      xmax = RDB[uni + UNIVERSE_MAXX];
      ymin = RDB[uni + UNIVERSE_MINY];
      ymax = RDB[uni + UNIVERSE_MAXY];
      zmin = RDB[uni + UNIVERSE_MINZ];
      zmax = RDB[uni + UNIVERSE_MAXZ];

      /* Compare to geometry boundaries */
      /*
      if (xmin < RDB[DATA_GEOM_MINX])
        xmin = RDB[DATA_GEOM_MINX];
      if (xmax < RDB[DATA_GEOM_MAXX])
        xmax = RDB[DATA_GEOM_MAXX];

      if (ymin < RDB[DATA_GEOM_MINY])
        ymin = RDB[DATA_GEOM_MINY];
      if (ymax < RDB[DATA_GEOM_MAXY])
        ymax = RDB[DATA_GEOM_MAXY];

      if (zmin < RDB[DATA_GEOM_MINZ])
        zmin = RDB[DATA_GEOM_MINZ];
      if (zmax < RDB[DATA_GEOM_MAXZ])
        zmax = RDB[DATA_GEOM_MAXZ];
      */
      /* Allocate memory for size vector */

      szp = ReallocMem(DATA_ARRAY, depth + 1);
      
      /* Put values */

      WDB[szp] = (double)sz0;
      for (n = 1; n < depth; n++)
        WDB[szp + n] = (double)sz;
      
      WDB[szp + n] = -1.0;      

      /* Put boundaries */

      lims[0] = xmin;
      lims[1] = xmax;
      lims[2] = ymin;
      lims[3] = ymax;
      lims[4] = zmin;
      lims[5] = zmax;

      /* Create search mesh */

      msh = CreateMesh(MESH_TYPE_ADAPTIVE, MESH_CONTENT_PTR,
                       MESH_CONTENT_PTR, 1, 0, 0, lims, szp);
      WDB[uni + UNIVERSE_PTR_SEARCH_MESH] = (double)msh;

      /* Init data */

      for (k = 0; k < sz0; k++)
        for (j = 0; j < sz0; j++)
          for (i = 0; i < sz0; i++)
            {
              /* Create search item */

              loc1 = NewItem(uni + UNIVERSE_PTR_SEARCH_MESH_DATA, 
                             CELL_SEARCH_MESH_SIZE);

              /* Get mesh pointer */
              
              loc0 = ReadMeshPtr(msh, i, j, k);
              CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);
              
              /* Put pointer */

              WDB[loc0] = (double)loc1;

              /* Calculate index */

              idx = i + j*sz0 + k*sz0*sz0;              
              
              /* Put mesh pointer and cell index */

              WDB[loc1 + CELL_SEARCH_MESH_PTR_MSH] = (double)msh;
              WDB[loc1 + CELL_SEARCH_MESH_IDX] = (double)idx;
            }
        
      /***********************************************************************/

      /***** Adapt mesh by sampling random points ****************************/

      /* Loop over iterations */

      for (it = 0; it < depth - 1; it++)
        {
          /******************************************************************/
      
          /***** Sample random points in mesh *******************************/

          /* Pointer to mesh */

          msh = (long)RDB[uni + UNIVERSE_PTR_SEARCH_MESH];
          CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

          /* Number of points */


          id = 0;

          /* Loop over mesh */

#ifdef OPEN_MP
#pragma omp parallel for private(id, idx, xmin, xmax, ymin, ymax, zmin, zmax, l, seed, x, y, z, loc0, loc1, ptr, nc, n, cell)
#endif

          for (idx = 0; idx < sz0*sz0*sz0; idx++)
            {
              /* Get boundaries */
              
              MeshCellBounds(msh, idx, &xmin, &xmax, &ymin, &ymax, 
                             &zmin, &zmax);
              
              /* Get Open MP thread id */

              id = OMP_THREAD_NUM;

              /* Loop over random points */
              
              for (l = 0; l < np; l++)
                {
                  /* Init random number sequence */
                  
                  seed = ReInitRNG(it*sz0*sz0*sz0*np + idx*np + l); 
                  SEED[id*RNG_SZ] = seed;
                  
                  /* Sample point */
                  
                  x = RandF(id)*(xmax - xmin) + xmin;
                  y = RandF(id)*(ymax - ymin) + ymin;
                  
                  if ((long)RDB[uni + UNIVERSE_DIM] == 3)
                    z = RandF(id)*(zmax - zmin) + zmin;
                  else
                    z = 0.0;
                  
                  /* Get mesh pointer */
                  
                  loc0 = MeshPtr(msh, x, y, z);
                  CheckPointer(FUNCTION_NAME, "(loc02)", DATA_ARRAY, loc0);
                  
                  loc0 = (long)RDB[loc0];
                  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);
                  
                  /* Get number of previously found cells */
                  
                  nc = (long)RDB[loc0 + CELL_SEARCH_MESH_N];

                  /* Loop over previous points */
                  
                  for (n = 0; n < nc; n++)
                    {
                      /* Get pointer to cell */
                      
                      cell = (long)RDB[loc0 + CELL_SEARCH_MESH_PTR_C1 + n];
                      CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);
                      
                      /* Test cell */
                      
                      if (InCell(cell, x, y, z, NO, id) == YES)
                        break;
                    }
                  
                  /* Check if found */
                  
                  if (n < nc)
                    continue;
                  
                  /* Loop over all cells in universe */
                  
                  ptr = (long)RDB[uni + UNIVERSE_PTR_CELL_LIST];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                  
                  n = 0;
                  while ((loc1 = ListPtr(ptr, n++)) > VALID_PTR)
                    {
                      /* Pointer to cell */
                      
                      cell = (long)RDB[loc1 + CELL_LIST_PTR_CELL];
                      CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);
          
                      /* Test cell */
                          
                      if (InCell(cell, x, y, z, NO, id) == YES)
                        {
                          /* Check count */
                          
                          if (nc < split)
                            {
                              /* Add pointer in list */
                              
                              WDB[loc0 + CELL_SEARCH_MESH_PTR_C1 + nc] = 
                                (double)cell;
                              
                              /* Add count */
                              
                              WDB[loc0 + CELL_SEARCH_MESH_N] = 
                                (double)(++nc);
                            }
                          
                          /* Break loop */
                          
                          break;
                        }
                    }
                }
            }
          
          /******************************************************************/

          /***** Check splits ***********************************************/

          /* Get pointer to mesh data */

          loc1 = (long)RDB[uni + UNIVERSE_PTR_SEARCH_MESH_DATA];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          /* Pointer to last */

          if ((last = LastItem(loc1)) < VALID_PTR)
            break;

          /* Loop over data */ 
          
          while (loc1 > VALID_PTR)
            {
              /* Check splitting */
              
              if ((long)RDB[loc1 + CELL_SEARCH_MESH_N] == split)
                {
                  /* Get mesh pointer and index */

                  msh = (long)RDB[loc1 + CELL_SEARCH_MESH_PTR_MSH];
                  idx = (long)RDB[loc1 + CELL_SEARCH_MESH_IDX];
                  
                  /* Get pointer to pointer */
                  
                  loc2 = (long)RDB[msh + MESH_PTR_PTR] + idx;
                  CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);
                  
                  /* Split cell */

                  SplitCell(uni, msh, idx);
                }
              
              /* Check if last */

              if (loc1 == last)
                break;
                
              /* Next cell */
                  
              loc1 = NextItem(loc1);
            }

          /* Put new list pointer */
          
          WDB[uni + UNIVERSE_PTR_SEARCH_MESH_DATA] = (double)NextItem(last);
 
          /******************************************************************/
        }

      ptr = (long)RDB[uni + UNIVERSE_PTR_CELL_LIST];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      fprintf(outp, " - Universe %s -- OK (%ld cells, %ld levels)\n", 
              GetText(uni + UNIVERSE_PTR_NAME), ListSize(ptr), it + 1);

      /***********************************************************************/
    
      /* Next universe */

      uni = NextItem(uni);
    }

  /* Exit subroutine */

  fprintf(outp, "\nOK.\n\n");

  /***************************************************************************/
}

/*****************************************************************************/

/***** Split cell ************************************************************/

void SplitCell(long uni, long msh, long idx)
{
  long szp, sz, loc0, loc1, loc2, i, j, k;
  double lims[6], x0, y0, dx, dy;

  /* Get pointer to location of splitted cell */
                  
  loc2 = (long)RDB[msh + MESH_PTR_PTR] + idx;
  CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

  /* Get boundaries */
  
  MeshCellBounds(msh, idx, &lims[0], &lims[1], &lims[2], &lims[3], &lims[4],
                 &lims[5]);

  /* Get size pointer */
  
  szp = (long)RDB[msh + MESH_ADA_PTR_SZ];
  CheckPointer(FUNCTION_NAME, "(szp)", DATA_ARRAY, szp);
  
  /* Get size */

  sz = (long)RDB[szp];
  CheckValue(FUNCTION_NAME, "sz", "", sz, 1, 10000000);

  x0 = lims[0];
  y0 = lims[2];
  dx = (lims[1] - lims[0])/sz;
  dy = (lims[3] - lims[2])/sz;
  
  /*
  for (i = 0; i < sz; i++)
  for (j = 0; j < sz; j++)
    printf("rectangle(\"position\", [%E,%E,%E,%E]);\n",
           x0 + i*dx, y0 + j*dy, dx, dy);
  */

  /* Create new search mesh */
              
  msh = CreateMesh(MESH_TYPE_ADAPTIVE, MESH_CONTENT_PTR, MESH_CONTENT_PTR, 
                   1, 0, 0, lims, szp);

  /* Init data */

  for (k = 0; k < sz; k++)
    for (j = 0; j < sz; j++)
      for (i = 0; i < sz; i++)
        {
          /* Create search item */

          loc1 = NewItem(uni + UNIVERSE_PTR_SEARCH_MESH_DATA, 
                         CELL_SEARCH_MESH_SIZE);

          /* Get mesh pointer */
              
          loc0 = ReadMeshPtr(msh, i, j, k);
          CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);
              
          /* Put pointer */

          WDB[loc0] = (double)loc1;
          
          /* Calculate index */
          
          idx = i + j*sz + k*sz*sz;              
              
          /* Put mesh pointer and cell index */
          
          WDB[loc1 + CELL_SEARCH_MESH_PTR_MSH] = (double)msh;
          WDB[loc1 + CELL_SEARCH_MESH_IDX] = (double)idx;
        }

  /* Put pointer */

  WDB[loc2] = -msh;
}

/*****************************************************************************/

