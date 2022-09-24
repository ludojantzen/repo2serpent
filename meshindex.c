/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : meshindex.c                                    */
/*                                                                           */
/* Created:       2011/05/13 (JLe)                                           */
/* Last modified: 2017/12/08 (JLe)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Returns index to mesh structure                              */
/*                                                                           */
/* Comments: - Tolerance parameter is used to test if point is too close to  */
/*             mesh boundaries (not used for hex lattices)                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MeshIndex:"

/*****************************************************************************/

long MeshIndex(long msh, double x, double y, double z, double tol)
{
  long n0, n1, n2, i, j, k, idx, ptr;
  double min0, max0, min1, max1, min2, max2, r, phi, theta;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

  /* Get sizes */
  
  n0 = (long)RDB[msh + MESH_N0];
  n1 = (long)RDB[msh + MESH_N1];
  n2 = (long)RDB[msh + MESH_N2];

  /* Reset index */

  idx = -1;

  /* Get dimensions */

  min0 = RDB[msh + MESH_MIN0];
  max0 = RDB[msh + MESH_MAX0];
  min1 = RDB[msh + MESH_MIN1];
  max1 = RDB[msh + MESH_MAX1];
  min2 = RDB[msh + MESH_MIN2];
  max2 = RDB[msh + MESH_MAX2];
  
  /* Check type */

  if (((long)RDB[msh + MESH_TYPE] == MESH_TYPE_CARTESIAN) ||
      ((long)RDB[msh + MESH_TYPE] == MESH_TYPE_ADAPTIVE))
    {
      /***********************************************************************/
      
      /***** Cartesian mesh **************************************************/

      /* Get local normalised co-ordinates */
      
      if (max0 - min0 > 0.0)
        x = (x - min0)/(max0 - min0);
      else
        x = 0.0;
      
      if (max1 - min1 > 0.0)
        y = (y - min1)/(max1 - min1);
      else
        y = 0.0;
      
      if (max2 - min2 > 0.0)
        z = (z - min2)/(max2 - min2);
      else
        z = 0.0;

      /* Check that point is inside region */
      
      if ((x < 0.0) || (x >= 1.0) ||
          (y < 0.0) || (y >= 1.0) ||
          (z < 0.0) || (z >= 1.0))
        return -1;
      
      /* Calculate indexes */
      
      i = (long)(x*n0);
      j = (long)(y*n1);
      k = (long)(z*n2);

      /* Check tolerance */

      if (tol > 0.0)
        if ((fabs(x*n0 - (double)i) < tol) || 
            (fabs(y*n1 - (double)j) < tol) ||
            (fabs(z*n2 - (double)k) < tol))
          return -2;

      /* Calculate index */

      idx = i + j*n0 + k*n0*n1;
      
      /***********************************************************************/
    }
  else if ((long)RDB[msh + MESH_TYPE] == MESH_TYPE_CYLINDRICAL)
    {
      /***********************************************************************/

      /***** Cylindrical mesh ************************************************/

      /* Calculate radius */

      r = sqrt(x*x + y*y);

      /* Calculate angle */

      phi = PolarAngle(x, y);

      /* Check phi */

      CheckValue(FUNCTION_NAME, "phi", "", phi, 0.0, 2.0*PI);

      /* Get local normalised co-ordinates */
      
      if (max0 - min0 > 0.0)
        r = (r - min0)/(max0 - min0);
      else
        r = 0.0;
      
      if (max1 - min1 > 0.0)
        phi = (phi - min1)/(max1 - min1);
      else
        phi = 0.0;
      
      /* (invert z-co-ordinate) */
      
      if (max2 - min2 > 0.0)
        z = 1.0 - (z - min2)/(max2 - min2);
      else
        z = 0.0;
      
      /* Check that point is inside region */
      
      if ((r < 0.0) || (r >= 1.0) ||
          (phi < 0.0) || (phi >= 1.0) ||
          (z < 0.0) || (z >= 1.0))
        return -1;
      
      /* Calculate indexes */
      
      i = (long)(r*n0);
      j = (long)(phi*n1);
      k = (long)(z*n2);

      /* Check tolerance */

      if (tol > 0.0)
        if ((fabs(r*n0 - (double)i) < tol) || 
            (fabs(phi*n1 - (double)j) < tol) ||
            (fabs(z*n2 - (double)k) < tol))
          return -2;

      /* Calculate index */

      idx = i + j*n0 + k*n0*n1;

      /***********************************************************************/
    }
  else if ((long)RDB[msh + MESH_TYPE] == MESH_TYPE_SPHERICAL)
    {
      /***********************************************************************/

      /***** Spherical mesh **************************************************/

      /* Calculate radius */

      r = sqrt(x*x + y*y + z*z);

      /* Calculate angles */

      phi = PolarAngle(x, y);

      if (r == 0.0)
        theta = 0.0;
      else
        theta = acos(z/r);

      /* Set the maximum theta to be slightly below the defined value of PI.
       * (otherwise theta=PI returns -1)  */

      if ((theta >= PI) && (theta - PI < 1e-14))
        theta = 0.9999999999*PI;

      /* Check */

      CheckValue(FUNCTION_NAME, "phi", "", phi, 0.0, 2.0*PI);
      CheckValue(FUNCTION_NAME, "theta", "", theta, 0.0, PI);

      /* Get local normalised co-ordinates */
      
      if (max0 - min0 > 0.0)
        r = (r - min0)/(max0 - min0);
      else
        r = 0.0;
      
      if (max1 - min1 > 0.0)
        phi = (phi - min1)/(max1 - min1);
      else
        phi = 0.0;
      
      if (max2 - min2 > 0.0)
        theta = (theta - min2)/(max2 - min2);
      else
        theta = z;
      
      /* Check that point is inside region */
      
      if ((r < 0.0) || (r >= 1.0) ||
          (phi < 0.0) || (phi >= 1.0) ||
          (theta < 0.0) || (theta >= 1.0))
        return -1;
      
      /* Calculate indexes */
      
      i = (long)(r*n0);
      j = (long)(phi*n1);
      k = (long)(theta*n2);

      /* Check tolerance */

      if (tol > 0.0)
        if ((fabs(r*n0 - (double)i) < tol) || 
            (fabs(phi*n1 - (double)j) < tol) ||
            (fabs(theta*n2 - (double)k) < tol))
          return -2;

      /* Calculate index */

      idx = i + j*n0 + k*n0*n1;

      /***********************************************************************/
    }
  else if (((long)RDB[msh + MESH_TYPE] == MESH_TYPE_HEXX) ||
           ((long)RDB[msh + MESH_TYPE] == MESH_TYPE_HEXY))
    {
      /***********************************************************************/

      /***** Hexagonal mesh **************************************************/

      /* Variables: min0 = x0, min1 = y0, max0 = pitch */
      /* Coordinate transformation relative to origin  */

      x = x - min0;
      y = y - min1;

      /* Adjust if even number of cells */

      x = x - (1 - (n0 % 2))*0.5*max0;
      y = y - (1 - (n1 % 2))*0.5*max0;

      /* Adjust axial coordinate */
      
      if (max2 - min2 > 0.0)
        z = (z - min2)/(max2 - min2);
      else
        z = 0.0;

      /* Check */

      if ((z < 0.0) || (z >= 1.0))
        return -1;

      /* Get hex indexes */

      if ((long)RDB[msh + MESH_TYPE] == MESH_TYPE_HEXX)
        GetLatticeIndexes(max0, max0, 1.0, x, y, 0.0, &i, &j, &k, LAT_TYPE_HX);
      else
        GetLatticeIndexes(max0, max0, 1.0, x, y, 0.0, &i, &j, &k, LAT_TYPE_HY);

      /* Calculate indexes */
      
      i = i + (long)(((double)n0 - 1.0)/2.0);
      j = j + (long)(((double)n1 - 1.0)/2.0);
      k = (long)(z*n2);

      /* Check */

      if ((i < 0) || (i > n0 - 1))
        return -1;
      else if ((j < 0) || (j > n1 - 1))
        return -1;
      else if ((k < 0) || (k > n2 - 1))
        Die(FUNCTION_NAME, "Indexing error");

      /* Calculate index */

      idx = i + j*n0 + k*n0*n1;
      
      /***********************************************************************/
    }
  else if ((long)RDB[msh + MESH_TYPE] == MESH_TYPE_ORTHOGONAL)
    {
      /***********************************************************************/

      /***** Unevenly-spaced orthogonal mesh *********************************/

      /* Search x-index */
      
      if (max0 - min0 > 0.0)
        {
          /* Pointer to array */

          ptr = (long)RDB[msh + MESH_ORTHO_PTR_XLIM];
          
          /* Search array */

          i = SearchArray(&RDB[ptr], x, n0 + 1);
        }
      else
        i = 0;

      /* Search y-index */
      
      if (max1 - min1 > 0.0)
        {
          /* Pointer to array */

          ptr = (long)RDB[msh + MESH_ORTHO_PTR_YLIM];
          
          /* Search array */

          j = SearchArray(&RDB[ptr], y, n1 + 1);
        }
      else
        j = 0;
      
      /* Search z-index */
      
      if (max2 - min2 > 0.0)
        {
          /* Pointer to array */

          ptr = (long)RDB[msh + MESH_ORTHO_PTR_ZLIM];
          
          /* Search array */

          k = SearchArray(&RDB[ptr], z, n2 + 1);
        }
      else
        k = 0;

      /* Check that point is inside region */
      
      if ((i < 0) || (j < 0) || (k < 0))
        return -1;

      /* Calculate index */

      idx = i + j*n0 + k*n0*n1;
      
      /***********************************************************************/
    }
  else if ((long)RDB[msh + MESH_TYPE] == MESH_TYPE_ICYL)
    {
      /***********************************************************************/

      /***** Unevenly-spaced cylindrical mesh ********************************/

      /* Calculate radius */

      r = sqrt(x*x + y*y);

      /* Calculate angle */

      phi = PolarAngle(x, y);

      /* Check phi */

      CheckValue(FUNCTION_NAME, "phi", "", phi, 0.0, 2.0*PI);

      /* Search radial index */
      
      if (max0 - min0 > 0.0)
        {
          /* Pointer to array */

          ptr = (long)RDB[msh + MESH_ORTHO_PTR_XLIM];
          
          /* Search array */

          i = SearchArray(&RDB[ptr], r, n0 + 1);
        }
      else
        i = 0;

      /* Search angular index */
      
      if (max1 - min1 > 0.0)
        {
          /* Pointer to array */

          ptr = (long)RDB[msh + MESH_ORTHO_PTR_YLIM];
          
          /* Search array */

          j = SearchArray(&RDB[ptr], phi, n1 + 1);
        }
      else
        j = 0;
      
      /* Search axial index */
      
      if (max2 - min2 > 0.0)
        {
          /* Pointer to array */

          ptr = (long)RDB[msh + MESH_ORTHO_PTR_ZLIM];
          
          /* Search array */

          k = SearchArray(&RDB[ptr], z, n2 + 1);
        }
      else
        k = 0;

      /* Check that point is inside region */
      
      if ((i < 0) || (j < 0) || (k < 0))
        return -1;

      /* Calculate index */

      idx = i + j*n0 + k*n0*n1;
      
      /***********************************************************************/
    }
  else
    Die(FUNCTION_NAME, "Invalid mesh type");

  /* Return index */

  return idx;
}

/*****************************************************************************/
