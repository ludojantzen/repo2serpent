/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : Torusdis.c                                     */
/*                                                                           */
/* Created:       2015/07/20 (JLe)                                           */
/* Last modified: 2015/10/07 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Calculates distance to elliptical torus.                     */
/*                                                                           */
/* Comments: - Centered at (0,0,0) on xy-plane, toroidal radius R, poloidal  */
/*             radii r1 (R) and r2 (z)                                       */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "TorusDis:"

/* Additional functions */

double Iterate(double min, double max, const double *coe);
double PolyVal(const double *coe, const long n, const double x);

/*****************************************************************************/

double TorusDis(double x, double y, double z, double u, double v, double w, 
                double R, double r1, double r2) 
{
  double A, B, C, D, dis, rad, x0, y0, z0;
  double a, b, c, d, e, a0, a1, a2, Q, S, T, t, theta, max;
  double coe[5], lims[5], zmin, zmax, Rmax;
  long n, i;
 
  /* Nää debuggaukseen: */

  x0 = x;
  y0 = y;
  z0 = z;

  /* Check radii */

  CheckValue(FUNCTION_NAME, "R", "", R, ZERO, INFTY);
  CheckValue(FUNCTION_NAME, "r1", "", r1, ZERO, INFTY);
  CheckValue(FUNCTION_NAME, "r2", "", r2, ZERO, INFTY);

  /* Calculate bounding planes and radius */

  zmin = -r1 - 1E-3;
  zmax = r1 + 1E-3;
  Rmax = R + r2 + 1E-3;

  /* Reset distance */

  dis = 0.0;

  /***************************************************************************/

  /***** Trivial cases *******************************************************/

  /* Check axial planes */
  
  if (z > zmax)
    {
      /* Above, check line-of-sight */

      if (w >= 0.0)
        return INFTY;

      /* Calculate distance */

      d = ZDis(z, w, zmax);
      CheckValue(FUNCTION_NAME, "d1", "", d, ZERO, 1E+12);
          
      /* Move over boundary */

      x = x + (d + EXTRAP_L)*u;
      y = y + (d + EXTRAP_L)*v;
      z = z + (d + EXTRAP_L)*w;
      
      /* Add to total */
      
      dis = dis + d + EXTRAP_L;
    }
  else if (z < zmin)
    {
      /* Below, check line-of-sight */

      if (w <= 0.0)
        return INFTY;

      /* Calculate distance */

      d = ZDis(z, w, zmin);
      CheckValue(FUNCTION_NAME, "d2", "", d, ZERO, 1E+12);
      
      /* Move over boundary */
      
      x = x + (d + EXTRAP_L)*u;
      y = y + (d + EXTRAP_L)*v;
      z = z + (d + EXTRAP_L)*w;
      
      /* Add to total */
      
      dis = dis + d + EXTRAP_L;
    }

  /* Check outer radius */

  if ((rad = sqrt(x*x + y*y)) > Rmax)
    {
      /* Outside, Calculate distance */

      d = CylDis(x, y, u, v, Rmax);
      CheckValue(FUNCTION_NAME, "d3", "", d, ZERO, INFTY);

      /* Check line-of-sight */
          
      if (d == INFTY)
        return INFTY;

      /* Compare to distances to axial planes */

      if (w > 0.0)
        { 
          /* Upper */

          if (d - ZDis(z, w, zmax) > -1E-6)
            return INFTY;
        }
      else if (w < 0.0)
        { 
          /* Lower */

          if (d - ZDis(z, w, zmin) > -1E-6)
            return INFTY;
        }

      /* Move over boundary */
      
      x = x + (d + EXTRAP_L)*u;
      y = y + (d + EXTRAP_L)*v;
      z = z + (d + EXTRAP_L)*w;

      /* Add to total */
          
      dis = dis + d + EXTRAP_L;
    }

  /* Sanity checks */

  if ((z - zmax > 1E-6) || (z - zmin < -1E-6))
    Die(FUNCTION_NAME, "Error in axial position: %E %E %E : %E %E", zmin, z,
        zmax, z - zmin, z - zmax);
  else if (sqrt(x*x + y*y) - Rmax > 1E-6)
    Die(FUNCTION_NAME, "Error in radial position");

  /* Calculate distance to nearest boundary */

  max = INFTY;

  if (((d = ZDis(z, w, zmin)) > 0.0) && (d < max))
    max = d;

  if (((d = ZDis(z, w, zmax)) > 0.0) && (d < max))
    max = d;

  if (((d = CylDis(x, y, u, v, Rmax)) > 0.0) && (d < max))
    max = d;

  if (max == INFTY)
    return INFTY;

  /* Check value */

  if (max > 2.0*sqrt(Rmax*Rmax + zmax*zmax))
    printf("%1.14E %1.14E %1.14E : %1.14E %1.14E %1.14E\n", x0, y0, z0,
           u, v, w);

  CheckValue(FUNCTION_NAME, "max", "", max, ZERO, 
             2.0*sqrt(Rmax*Rmax + zmax*zmax));

  /***************************************************************************/

  /***** Calculate polynomial coefficients ***********************************/
  
  /* Check type and calculate intermediate constants */

  if (r1 != r2)
    {
      /* Elliptical torus */

      r1 = r1*r1;
      r2 = r2*r2;
      A = w*w/r1 + (1.0 - w*w)/r2;
      B = 2.0*(u*x + v*y)/r2 + (2.0*w*z)/r1;
      C = z*z/r1 + (R*R + x*x + y*y)/r2 - 1.0;
      D = 4.0*R*R/(r2*r2);
    }
  else
    {
      /* Ring torus */
      
      A = 1.0;
      B = 2.0*(u*x + v*y + w*z);
      C = R*R - r1*r1 + x*x + y*y + z*z;
      D = 4.0*R*R;
    }
  
  /* Polynomial coefficients : at^4 + bt^3 + ct^2 + dt + e */
  
  a = A*A;
  b = 2.0*A*B;
  c = B*B + 2.0*A*C - D*(1.0 - w*w);
  d = 2.0*(B*C - D*(u*x + v*y));
  e = C*C - D*(x*x + y*y);

  /* Check values */

  CheckValue(FUNCTION_NAME, "a", "", a, ZERO, INFTY);
  CheckValue(FUNCTION_NAME, "b", "", b, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "c", "", c, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "d", "", d, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "e", "", e, -INFTY, INFTY);

  /* No roots if all coefficients are positive (a is always > 0, so */
  /* no need to test if all are negative) */

  if (b > 0.0)
    if (c > 0.0)
      if (d > 0.0)
        if (e > 0.0)
          return INFTY;

  /***************************************************************************/

  /***** Find local minima and maxima ****************************************/

  /* Polynomial coefficients for derivative: t^3 + a2t^2 + a1t + a0 */

  a2 = 3.0*b/(4.0*a);
  a1 = 2.0*c/(4.0*a);
  a0 = d/(4.0*a);

  /* Coefficients for calculating local minima and maxima */

  Q = (3.0*a1 - a2*a2)/9.0;
  R = (9.0*a2*a1 - 27.0*a0 - 2.0*a2*a2*a2)/54.0;
  D = Q*Q*Q + R*R;

  /* Reset count */

  n = 0.0;

  /* Add zero */

  lims[n++] = 0.0;

  /* Check discriminant */

  if (D > 0.0)
    {
      /* One root, calculate additional coefficents */

      S = cbrt(R + sqrt(D));
      T = cbrt(R - sqrt(D));
      
      /* Calculate root */

      t = -a2/3.0 + S + T;
      CheckValue(FUNCTION_NAME, "t", "", t, -INFTY, INFTY);

      if ((t > 0.0) && (t < max))
        lims[n++] = t;
    }
  else if (D == 0.0)
    {
      /* Two roots, calculate additional coefficients */
  
      S = cbrt(R);
      A = -a2/3.0;

      /* Calculate roots */

      t = A + 2.0*S;
      CheckValue(FUNCTION_NAME, "t", "", t, -INFTY, INFTY);

      if ((t > 0.0) && (t < max))
        lims[n++] = t;

      t = A - S;
      CheckValue(FUNCTION_NAME, "t", "", t, -INFTY, INFTY);

      if ((t > 0.0) && (t < max))
        lims[n++] = t;
    }
  else
    {
      /* Three roots, check coefficients */
 
      CheckValue(FUNCTION_NAME, "Q", "", Q, -INFTY, -ZERO);
      CheckValue(FUNCTION_NAME, "acos", "", R/sqrt(-Q*Q*Q), -1.0, 1.0);

      /* Calculate additional coefficients */

      theta = acos(R/sqrt(-Q*Q*Q));

      A = 2.0*sqrt(-Q);
      B = a2/3.0;

      /* Calculate roots */

      t = A*cos(theta/3.0) - B;
      CheckValue(FUNCTION_NAME, "t", "", t, -INFTY, INFTY);

      if ((t > 0.0) && (t < max))
        lims[n++] = t;
        
      t = A*cos((theta + 2.0*PI)/3.0) - B;
      CheckValue(FUNCTION_NAME, "t", "", t, -INFTY, INFTY);

      if ((t > 0.0) && (t < max))
        lims[n++] = t;
      
      t = A*cos((theta + 4.0*PI)/3.0) - B;
      CheckValue(FUNCTION_NAME, "t", "", t, -INFTY, INFTY);

      if ((t > 0.0) && (t < max))
        lims[n++] = t;
    }

  /* Add maximum and sort array */

  lims[n++] = max;
  SortArray(lims, n);

  /* Check size */

  CheckValue(FUNCTION_NAME, "n", "", n, 2, 5);

  /***************************************************************************/

  /***** Find smallest positive distance *************************************/

  /* Put coefficients in vector */

  coe[0] = a;
  coe[1] = b;
  coe[2] = c;
  coe[3] = d;
  coe[4] = e;

  /* Loop over intervals and iterate */

  t = -1.0;
  for (i = 1; i < n; i++)
    if ((t = Iterate(lims[i - 1], lims[i], coe)) > 0.0)
      break;
  
  /* Check value and add to total distance */

  if (t < 0.0)
    return INFTY;
  else
    dis = dis + t;

  /* Return distance */

  return dis;

  /***************************************************************************/
}

/*****************************************************************************/

/***** Evaluate N:th order polynomial ****************************************/

double PolyVal(const double *coe, const long N, const double x)
{
  long n;
  double y, z;

  /* Reset value */

  y = 0.0;
  z = 1.0;

  /* Loop over terms */

  for (n = 0; n < N; n++)
    {     
      y = y + z*coe[N - n];
      z = z*x;
    }

  /* Highest order term */

  y = y + z*coe[0];
           
  /* Return value */

  return y;
}

/*****************************************************************************/

/***** Iterate root **********************************************************/

/* Bisection search (NOTE: No major speed-up with Newton) */

double Iterate(double min, double max, const double *coe)
{
  double x, f0, f1, f;

  /* Initial values */

  f0 = PolyVal(coe, 4, min);
  f1 = PolyVal(coe, 4, max);

  /* Check */

  if (f0*f1 > 0.0)
    return -1.0;

  /* Iterate */
  
  while (1 != 2)
    {
      /* Make a guess */

      x = 0.5*(min + max);
      f = PolyVal(coe, 4, x);
   
      /* Swap limiting value */

      if (f*f0 > 0.0)
        {
          min = x;
          f0 = f;
        }
      else if (f*f1 > 0.0)
        {
          max = x;
          f1 = f;
        }
      else if (fabs(f) < 1E-9)
        break;
      else
        Die(FUNCTION_NAME, "WTF? %E %E", f*f0, f*f1);

      /* Check convergence */

      if (fabs(max - min) < 1E-9)
        break;
    }

  /* Return distance */

  return x;
}

/*****************************************************************************/
