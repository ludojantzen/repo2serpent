/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : wwinsrc.c                                      */
/*                                                                           */
/* Created:       2017/04/02 (JLe)                                           */
/* Last modified: 2018/11/08 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Sample source from weight window mesh                        */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "WWinSrc:"

/*****************************************************************************/

void WWinSrc(long src, long type, double *x, double *y, double *z, double E, 
             double *wgt, long id)
{
  long wwd, loc0, ptr, msh, idx, mtype, n, erg, ne;
  double xmin, xmax, ymin, ymax, zmin, zmax, Itot, Imax, I;

  /* Pointer to weight window structure */
  
  wwd = (long)RDB[DATA_PTR_WWD0];
  CheckPointer(FUNCTION_NAME, "(wwd)", DATA_ARRAY, wwd);

  /* Loop over weight window structures */

  while (wwd > VALID_PTR)
    {
      /* Check particle type */

      if (((long)RDB[wwd + WWD_PARTICLE_TYPE] > 0) &&
          ((long)RDB[wwd + WWD_PARTICLE_TYPE] != type))
        {
          /* Pointer to next */

          wwd = NextItem(wwd);

          /* Cycle loop */

          continue;
        }

      /* Check source biasing */

      if ((long)RDB[wwd + WWD_SRC_BIAS] == NO)
        return;

      Die(FUNCTION_NAME, "Tähän tarvii sampling listan");

      /* Check starting weight */

      if (*wgt != 1.0)
        Error(wwd, "Non-analog source sampling cannot be used with biasing");

      /* Number of energy groups */
              
      ne = (long)RDB[wwd + WWD_NE];
      CheckValue(FUNCTION_NAME, "ne", "", ne, 1, 10000);
      
      /* Reset energy group */
      
      n = 0;
      
      /* Find energy group if defined */
      
      if ((erg = (long)RDB[wwd + WWD_PTR_ERG]) > VALID_PTR)
        if ((n = SearchArray(&RDB[erg], E, ne + 1)) < 0)
          Error(wwd, "Particle energy %1.2E beyond energy boundaries", E);

      /* Get pointer to importance vector */

      ptr = (long)RDB[wwd + WWD_PTR_SUM_SRC_IMP];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      
      /* Get total importance */

      Imax = RDB[ptr + n];
      CheckValue(FUNCTION_NAME, "Imax", "", Imax, ZERO, INFTY);

      /* Avoid compiler warning */

      I = 1.0;

      /* Sample fraction of total importance */

      Itot = Imax*RandF(id);

      /* Loop over data */

      loc0 = (long)RDB[wwd + WWD_PTR_MESH_DATA];
      while (loc0 > VALID_PTR)
        {
          /* Pointer to source importance */

          ptr = (long)RDB[loc0 + RMX_CELL_IMP_SRC];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Get importance */

          I = RDB[ptr + n];
          CheckValue(FUNCTION_NAME, "I", "", I, 0.0, INFTY);

          /* Subtract importance */

          if ((Itot = Itot - I) < 0.0)
            break;
          
          /* Next */

          loc0 = NextItem(loc0);
        }

      /* Pointer to mesh */

      msh = (long)RDB[wwd + WWD_PTR_MESH];
      CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

      /* Get mesh type */

      mtype = (long)RDB[msh + MESH_TYPE];

      /* Check importance */

      CheckValue(FUNCTION_NAME, "I", "", I, ZERO, Imax);

      /* Set weight */

      *wgt = 1.0/I;
      CheckValue(FUNCTION_NAME, "*wgt", "", *wgt, ZERO, INFTY);

      /* Get mesh index */

      idx = (long)RDB[loc0 + RMX_CELL_MESH_IDX];
      
      /* Get boundaries */

      MeshCellBounds(msh, idx, &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);
      
      /* Check type */

      if (mtype == MESH_TYPE_CARTESIAN)
        {
          /* Sample point */
          
          *x = RandF(id)*(xmax - xmin) + xmin;
          *y = RandF(id)*(ymax - ymin) + ymin;
          *z = RandF(id)*(zmax - zmin) + zmin;
        }
      else
        Die(FUNCTION_NAME, "Mesh type not supported");
      
#ifdef DEBUG
      
      /* Check that index matches */
      
      if (idx != MeshIndex(msh, *x, *y, *z, -1.0))
        Die(FUNCTION_NAME, "Mismatch in mesh index");

#endif

      /* Break loop */

      break;
    }
}

/*****************************************************************************/
