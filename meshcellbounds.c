/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : meshcellbounds.c                               */
/*                                                                           */
/* Created:       2017/04/02 (JLe)                                           */
/* Last modified: 2018/06/30 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Calculates the boundaries of mesh cell                       */
/*                                                                           */
/* Comments: - First level boundaries are returned for adaptive mesh. This   */
/*             is used at least in preparesearchmesh.c.                      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MeshCellBounds:"

/*****************************************************************************/

void MeshCellBounds(long msh, long idx, double *min0, double *max0, 
                    double *min1, double *max1, double *min2, double *max2)
{
  long n0, n1, n2, i, j, k, type;
  double mind0, maxd0, mind1, maxd1, mind2, maxd2, p0, p1, p2, x0, y0;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

  /* Get dimensions */

  mind0 = RDB[msh + MESH_MIN0];
  maxd0 = RDB[msh + MESH_MAX0];
  mind1 = RDB[msh + MESH_MIN1];
  maxd1 = RDB[msh + MESH_MAX1];
  mind2 = RDB[msh + MESH_MIN2];
  maxd2 = RDB[msh + MESH_MAX2];

  /* Get sizes */
  
  n0 = (long)RDB[msh + MESH_N0];
  n1 = (long)RDB[msh + MESH_N1];
  n2 = (long)RDB[msh + MESH_N2];

  /* Calculate indexes */

  k = (long)((double)idx/((double)(n0*n1)));
  idx = idx - k*n0*n1;

  j = (long)((double)idx/((double)n0));
  i = idx - j*n0;

  /* Check values */
  
  CheckValue(FUNCTION_NAME, "i", "", i, 0, n0 - 1);
  CheckValue(FUNCTION_NAME, "j", "", j, 0, n1 - 1);
  CheckValue(FUNCTION_NAME, "k", "", k, 0, n2 - 1);
  
  /* Calculate mesh pitch */
      
  p0 = (maxd0 - mind0)/((double)n0);
  p1 = (maxd1 - mind1)/((double)n1);
  p2 = (maxd2 - mind2)/((double)n2);

  /* Get mesh type */

  type = (long)RDB[msh + MESH_TYPE];

  /* Check type */

  if ((type == MESH_TYPE_CARTESIAN) || (type == MESH_TYPE_ADAPTIVE))
    {
      /***********************************************************************/

      /***** Cartesian mesh **************************************************/

      /* Check pitches */
      
      CheckValue(FUNCTION_NAME, "p0", "", p0, ZERO, INFTY);
      CheckValue(FUNCTION_NAME, "p1", "", p1, ZERO, INFTY);
      CheckValue(FUNCTION_NAME, "p2", "", p2, ZERO, INFTY);
      
      /* Get boundaries */

      *min0 = mind0 + (double)i*p0;
      *max0 = *min0 + p0;

      *min1 = mind1 + (double)j*p1;
      *max1 = *min1 + p1;

      *min2 = mind2 + (double)k*p2;
      *max2 = *min2 + p2;

      /***********************************************************************/
    }
  else if ((type == MESH_TYPE_HEXX) || (type == MESH_TYPE_HEXY))
    {
      /***********************************************************************/

      /***** Hexagonal mesh **************************************************/

      /* Center coordinates */

      x0 = mind0;
      y0 = mind1;

      /* Get hexagonal pitch */

      p0 = maxd0;

      /* Check values */
      
      CheckValue(FUNCTION_NAME, "p0", "", p0, ZERO, INFTY);
      CheckValue(FUNCTION_NAME, "p2", "", p2, ZERO, INFTY);

      /* If even number of cells, shift origin by pitch/2. This results */
      /* from the fact that the local origin is not in the centre of a  */
      /* lattice element. */
      
      x0 = x0 + (1 - (n0 % 2))*0.5*p0;
      y0 = y0 + (1 - (n1 % 2))*0.5*p0;

      /* Shift indexes */

      i = i - (long)((double)n0*0.5);
      j = j - (long)((double)n1*0.5);
  
      /* Check type */
      
      if (type == MESH_TYPE_HEXX)
        {
          /* Calculate origin */

          x0 = x0 + (i + COS60*j)*p0; 
          y0 = y0 + j*SIN60*p0; 

          /* Calculate boundaries */

          *min0 = x0 - 0.5*p0;
          *max0 = x0 + 0.5*p0;
          *min1 = y0 - 0.5*p0/COS30;
          *max1 = y0 + 0.5*p0/COS30;


        }
      else if (type == MESH_TYPE_HEXY)
        {
          /* Calculate origin */

          x0 = x0 + j*SIN60*p0; 
          y0 = y0 + (i + COS60*j)*p0; 

          /* Calculate boundaries */

          *min0 = x0 - 0.5*p0/COS30;
          *max0 = x0 + 0.5*p0/COS30;
          *min1 = y0 - 0.5*p0;
          *max1 = y0 + 0.5*p0;
        }
      else
        Die(FUNCTION_NAME, "Error in hex mesh type");

      /* Get z-boundary */

      *min2 = mind2 + (double)k*p2;
      *max2 = *min2 + p2;

      /***********************************************************************/
    }
  else if (type == MESH_TYPE_CYLINDRICAL)
    {
      /***********************************************************************/

      /***** Cylindrical mesh ************************************************/

      /* Check pitches */
      
      CheckValue(FUNCTION_NAME, "p0", "", p0, ZERO, INFTY);
      CheckValue(FUNCTION_NAME, "p1", "", p1, ZERO, 2*PI);
      CheckValue(FUNCTION_NAME, "p2", "", p2, ZERO, INFTY);

      /* Get boundaries */

      *min0 = mind0 + (double)i*p0;
      *max0 = *min0 + p0;

      *min1 = mind1 + (double)j*p1;
      *max1 = *min1 + p1;

      *min2 = mind2 + (double)k*p2;
      *max2 = *min2 + p2;

      /***********************************************************************/
    }
  else
    Die(FUNCTION_NAME, "Mesh type %ld not supported", type);
}

/*****************************************************************************/
