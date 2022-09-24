/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : meshcelvol.c                                   */
/*                                                                           */
/* Created:       2014/12/12 (JLe)                                           */
/* Last modified: 2018/01/26 (JLe)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Returns volume of mesh cell at (x,y,z)                       */
/*                                                                           */
/* Comments: - Laskee sylinterihilalle koko rinkulan tilavuuden              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MeshCellVol:"

/*****************************************************************************/

double MeshCellVol(long msh, double x, double y, double z)
{
  long n0, n1, n2, i;
  double min0, max0, min1, max1, min2, max2, d0, d1, d2, V, r, r0, r1;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

  /* Get sizes */
  
  n0 = (long)RDB[msh + MESH_N0];
  n1 = (long)RDB[msh + MESH_N1];
  n2 = (long)RDB[msh + MESH_N2];

  CheckValue(FUNCTION_NAME, "n0", "", n0, 1, 10000000);
  CheckValue(FUNCTION_NAME, "n1", "", n1, 1, 10000000);
  CheckValue(FUNCTION_NAME, "n2", "", n2, 1, 10000000);

  /* Get dimensions */

  min0 = RDB[msh + MESH_MIN0];
  max0 = RDB[msh + MESH_MAX0];
  min1 = RDB[msh + MESH_MIN1];
  max1 = RDB[msh + MESH_MAX1];
  min2 = RDB[msh + MESH_MIN2];
  max2 = RDB[msh + MESH_MAX2];

  /* Check type */

  if ((long)RDB[msh + MESH_TYPE] == MESH_TYPE_CARTESIAN)
    {
      /* Check that point is inside */
      Die(FUNCTION_NAME, "Poista tää");
      if ((x < min0) || (x >= max0) || (y < min1) || (y >= max1) ||
          (z < min2) || (z >= max2))
        return 0.0;

      /* Calculate dimensions */
      
      d0 = (max0 - min0)/((double)n0);
      d1 = (max1 - min1)/((double)n1);
      d2 = (max2 - min2)/((double)n2);
      
      /* Calculate volume */
      
      V = d0*d1*d2;
      CheckValue(FUNCTION_NAME, "V", "", V, 0.0, INFTY);

      /* Return value */
      
      return V;
    }
  else if ((long)RDB[msh + MESH_TYPE] != MESH_TYPE_CYLINDRICAL)
    Die(FUNCTION_NAME, "Type not supported");

  /* Cylindrical mesh, calculate radius */

  r = sqrt(x*x + y*y);

  /* Get local normalised co-ordinates */
  
  if (max0 - min0 > 0.0)
    r = (r - min0)/(max0 - min0);
  else
    r = 0.0;

  /* Check that point is inside region */
      
  if ((r < 0.0) || (r >= 1.0) || (z < min2) || (z >= max2))
    return 0.0;

  /* Calculate index */
      
  i = (long)(r*n0);

  /* Get radii */

  r0 = ((double)i)/((double)n0)*(max0 - min0) + min0;
  r1 = ((double)i + 1.0)/((double)n0)*(max0 - min0) + min0;

  /* Calculate height */

  d2 = (max2 - min2)/((double)n2);

  /* Calculate volume */
      
  V = PI*(r1*r1 - r0*r0)*d2; 
  CheckValue(FUNCTION_NAME, "V", "", V, 0.0, INFTY);

  /* Return value */
      
  return V;
}

/*****************************************************************************/
