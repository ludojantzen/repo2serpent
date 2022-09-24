/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processsurfaces.c                              */
/*                                                                           */
/* Created:       2016/10/24 (JLe)                                           */
/* Last modified: 2016/11/14 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Additional processing for surfaces                           */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessSurfaces:"

/*****************************************************************************/

void ProcessSurfaces()
{
  double r0[3], v1[3], v2[3], v3[3], t0[3], t1[3], t2[3], p1[3], p2[3], p3[3];
  double d, rad, l, A, B, C, D;
  long surf, ptr, n, m, i, type, np;

  /* Loop over surfaces */

  surf = (long)RDB[DATA_PTR_S0];
  while (surf > VALID_PTR)
    {
      /* Get type and number of parameters */

      type = (long)RDB[surf + SURFACE_TYPE];
      np = (long)RDB[surf + SURFACE_N_PARAMS];

      /* Check type */

      if ((type == SURF_PLANE) && (np == 9))
        {
          /*******************************************************************/

          /***** Plane defined by three points *******************************/

          /* Pointer to parameter list */

          ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
          CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

          /* Two vectors on plane */

          v1[0] = RDB[ptr + 3] - RDB[ptr];
          v1[1] = RDB[ptr + 4] - RDB[ptr + 1];
          v1[2] = RDB[ptr + 5] - RDB[ptr + 2];

          v2[0] = RDB[ptr + 6] - RDB[ptr + 3];
          v2[1] = RDB[ptr + 7] - RDB[ptr + 4];
          v2[2] = RDB[ptr + 8] - RDB[ptr + 5];

          /* Normal vector */

          VectorProd(v3, v1, v2);

          /* First three constants are given by direction cosines */

          A = v3[0];
          B = v3[1];
          C = v3[2];
          
          /* Fourth constant is given by putting one of the points */
          /* on surface */

          D = A*RDB[ptr] + B*RDB[ptr + 1] + C*RDB[ptr + 2];

          /* Change number of parameters */

          WDB[surf + SURFACE_N_PARAMS] = 4.0;

          /* Put parameters */

          WDB[ptr] = A;
          WDB[ptr + 1] = B;
          WDB[ptr + 2] = C;
          WDB[ptr + 3] = D;

          /*******************************************************************/
        }
      else if ((type == SURF_MPLANE) && (np == 9))
        {
          /*******************************************************************/

          /***** MCNP-style Plane defined by three points ********************/

          /* Pointer to parameter list */

          ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
          CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

          /* Two vectors on plane */

          v1[0] = RDB[ptr + 3] - RDB[ptr];
          v1[1] = RDB[ptr + 4] - RDB[ptr + 1];
          v1[2] = RDB[ptr + 5] - RDB[ptr + 2];

          v2[0] = RDB[ptr + 6] - RDB[ptr + 3];
          v2[1] = RDB[ptr + 7] - RDB[ptr + 4];
          v2[2] = RDB[ptr + 8] - RDB[ptr + 5];

          /* Normal vector */

          VectorProd(v3, v1, v2);

          /* First three constants are given by direction cosines */

          A = v3[0];
          B = v3[1];
          C = v3[2];
          
          /* Fourth constant is given by putting one of the points */
          /* on surface */

          D = A*RDB[ptr] + B*RDB[ptr + 1] + C*RDB[ptr + 2];

          /* Reset invert flag */

          i = NO;

          /* Check with respect to origin */

          if (D == 0.0)
            {
              /* Plane passes through origin, special case */

              if (C < 0.0)
                i = YES;
              else if (C > 0.0)
                i = NO;
              else if (B < 0.0)
                i = YES;
              else if (B > 0.0)
                i = NO;
              else if (A < 0.0)
                i = YES;
              else if (A > 0.0)
                i = NO;
              else
                Error(surf, "Error in surface definition");
            }
          else if (D < 0.0)
            {
              /* Origin is on other side */

              i = YES;
            }

          /* Check if inverted */

          if (i == YES)
            {
              /* Normal vector (reverse order of product) */

              VectorProd(v3, v2, v1);
              
              /* First three constants are given by direction cosines */
              
              A = v3[0];
              B = v3[1];
              C = v3[2];
              
              /* Fourth constant is given by putting one of the points */
              /* on surface */
              
              D = A*RDB[ptr] + B*RDB[ptr + 1] + C*RDB[ptr + 2];
            }
          
          /* Change type and number of parameters */

          WDB[surf + SURFACE_TYPE] = (double)SURF_PLANE;
          WDB[surf + SURFACE_N_PARAMS] = 4.0;

          /* Put parameters */

          WDB[ptr] = A;
          WDB[ptr + 1] = B;
          WDB[ptr + 2] = C;
          WDB[ptr + 3] = D;

          /*******************************************************************/
        }
      else if (type == SURF_BOX)
        {
          /*******************************************************************/

          /***** MCNP-type box ***********************************************/

          /* Pointer to parameter list */

          ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
          CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

          /* Read data */

          r0[0] = RDB[ptr++];
          r0[1] = RDB[ptr++];
          r0[2] = RDB[ptr++];

          v1[0] = RDB[ptr++];
          v1[1] = RDB[ptr++];
          v1[2] = RDB[ptr++];

          v2[0] = RDB[ptr++];
          v2[1] = RDB[ptr++];
          v2[2] = RDB[ptr++];

          v3[0] = RDB[ptr++];
          v3[1] = RDB[ptr++];
          v3[2] = RDB[ptr++];

          /* Check lengths */

          d = ScalarProd(v1,v1);
          if (fabs(d) < 1E-6)
            Error(surf, "Vector 1 has zero length");

          d = ScalarProd(v2,v2);
          if (fabs(d) < 1E-6)
            Error(surf, "Vector 2 has zero length");

          d = ScalarProd(v3,v3);
          if (fabs(d) < 1E-6)
            Error(surf, "Vector 3 has zero length");  

          /* Avoid compiler warning */

          t0[0] = 0.0;
          t0[1] = 0.0;
          t0[2] = 0.0;

          t1[0] = 0.0;
          t1[1] = 0.0;
          t1[2] = 0.0;

          t2[0] = 0.0;
          t2[1] = 0.0;
          t2[2] = 0.0;

          /* Allocate memory for new data */

          ptr = ReallocMem(DATA_ARRAY, 54);
          WDB[surf + SURFACE_PTR_PARAMS] = (double)ptr;

          /* Put number of parameters */

          WDB[surf + SURFACE_N_PARAMS] = 54.0;

          /* Loop over two corner points */

          for (m = 0; m < 2; m++)
            {
              /* Loop over three planes */

              for (n = 0; n < 3; n++)
                {
                  /* Loop over coordinates */
                  
                  for (i = 0; i < 3; i++)
                    {
                      /* Select three points on surface */
                      
                      if (n == 0)
                        {
                          p1[i] = r0[i] + v1[i];
                          p2[i] = r0[i];
                          p3[i] = r0[i] + v2[i];
                        }
                      else if (n == 1)
                        {
                          p1[i] = r0[i] + v1[i];
                          p2[i] = r0[i];
                          p3[i] = r0[i] + v3[i];
                        }
                      else
                        {
                          p1[i] = r0[i] + v2[i];
                          p2[i] = r0[i];
                          p3[i] = r0[i] + v3[i];
                        }

                      /* Calculate temporary vectors for testing */
                      
                      t1[i] = p2[i] - p3[i];
                      t2[i] = p1[i] - p2[i];
                    }

                  /* Check lengths */
                      
                  d = ScalarProd(t1, t1);
                  if (fabs(d) < 1E-6)
                    Die(FUNCTION_NAME, "Vector t1 has zero lenght");
                  
                  d = ScalarProd(t2, t2);
                  if (fabs(d) < 1E-6)
                    Die(FUNCTION_NAME, "Vector t2 has zero lenght");
                  
                  /* Calculate surface normal */
                      
                  VectorProd(t0, t1, t2);

                  /* Check lenght */

                  d = ScalarProd(t0, t0);
                  if (fabs(d) < 1E-6)
                    Die(FUNCTION_NAME, "Vectors 1 and 2 in same direction");

                  /* Check direction with third vector */
                  
                  if (n == 0)
                    d = ScalarProd(t0, v3);
                  else if (n == 1)
                    d = ScalarProd(t0, v2);
                  else
                    d = ScalarProd(t0, v1);

                  /* Check scalar product */

                  if (fabs(d) < 1E-6)
                    Die(FUNCTION_NAME, "WTF");
                  else if (d < 0.0)
                    {
                      /* Store points in order 1, 2, 3 */

                      WDB[ptr++] = p1[0];
                      WDB[ptr++] = p1[1];
                      WDB[ptr++] = p1[2];

                      WDB[ptr++] = p2[0];
                      WDB[ptr++] = p2[1];
                      WDB[ptr++] = p2[2];

                      WDB[ptr++] = p3[0];
                      WDB[ptr++] = p3[1];
                      WDB[ptr++] = p3[2];
                    }
                  else
                    {
                      /* Store points in order 3, 2, 1 */

                      WDB[ptr++] = p3[0];
                      WDB[ptr++] = p3[1];
                      WDB[ptr++] = p3[2];

                      WDB[ptr++] = p2[0];
                      WDB[ptr++] = p2[1];
                      WDB[ptr++] = p2[2];

                      WDB[ptr++] = p1[0];
                      WDB[ptr++] = p1[1];
                      WDB[ptr++] = p1[2];
                    }
                }
              
              /* Update vectors */

              for (i = 0; i < 3; i++)
                {
                  /* Corner position */

                  r0[i] = r0[i] + v1[i] + v2[i]+ v3[i];

                  /* Invert directions */

                  v1[i] = -v1[i];
                  v2[i] = -v2[i];
                  v3[i] = -v3[i];
                }
            }

          /*******************************************************************/
        }

      else if (type == SURF_RCC)
        {
          /*******************************************************************/

          /***** MCNP-type truncated cylinder ********************************/

          /* Pointer to parameter list */

          ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
          CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

          /* Get position vector */

          r0[0] = RDB[ptr++];
          r0[1] = RDB[ptr++];
          r0[2] = RDB[ptr++];

          /* Get height vector */

          v1[0] = RDB[ptr++];
          v1[1] = RDB[ptr++];
          v1[2] = RDB[ptr++];

          /* Get radius */

          if ((rad = RDB[ptr++]) < ZERO)
            Error(surf, "Error in cylinder radius");

          /* Calculate length */

          if ((l = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2])) < ZERO)
            Error(surf, "Error in cylinder height");
          
          /* Normalize direction vector */

          v1[0] = v1[0]/l;
          v1[1] = v1[1]/l;
          v1[2] = v1[2]/l;

          /* Allocate memory for new data */

          ptr = ReallocMem(DATA_ARRAY, 9);
          WDB[surf + SURFACE_PTR_PARAMS] = (double)ptr;

          /* Put number of parameters */

          WDB[surf + SURFACE_N_PARAMS] = 9.0;

          /* Put data */

          WDB[ptr++] = r0[0];
          WDB[ptr++] = r0[1];
          WDB[ptr++] = r0[2];
          WDB[ptr++] = v1[0];
          WDB[ptr++] = v1[1];
          WDB[ptr++] = v1[2];
          WDB[ptr++] = rad;
          WDB[ptr++] = v1[0]*r0[0] + v1[1]*r0[1] + v1[2]*r0[2];
          WDB[ptr++] = v1[0]*(r0[0] + v1[0]*l) + v1[1]*(r0[1] + v1[1]*l)
            + v1[2]*(r0[2] + v1[2]*l);

          /*******************************************************************/
        }

      /* Next surface */

      surf = NextItem(surf);
    }
}

/*****************************************************************************/
