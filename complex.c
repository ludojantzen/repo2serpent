
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : complex.c                                      */
/*                                                                           */
/* Created:       2014/06/07 (MPu)                                           */
/* Last modified: 2016/09/30 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Elementary functions for complex numbers                     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "Complex:"

/*****************************************************************************/

complex c_add(complex x, complex y)
{
  complex z;

  /* Check values */

  CheckValue(FUNCTION_NAME, "x.re (c_add)", "", x.re, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "x.im (c_add)", "", x.im, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "y.re (c_add)", "", y.re, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "y.im (c_add)", "", y.im, -INFTY, INFTY);
  
  /* Calculate */
  
  z.re = x.re + y.re;
  z.im = x.im + y.im; 
  
  return z;
}

complex c_sub(complex x, complex y)
{
  complex z;
  
  /* Check values */
  
  CheckValue(FUNCTION_NAME, "x.re (c_sub)", "", x.re, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "x.im (c_sub)", "", x.im, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "y.re (c_sub)", "", y.re, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "y.im (c_sub)", "", y.im, -INFTY, INFTY);
  
  /* Calculate */
  
  z.re = x.re - y.re;
  z.im = x.im - y.im;
  
  return z;
}

complex c_mul(complex x, complex y)
{
  complex z;
  
  /* Check values */
  
  CheckValue(FUNCTION_NAME, "x.re (c_mul)", "", x.re, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "x.im (c_mul)", "", x.im, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "y.re (c_mul)", "", y.re, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "y.im (c_mul)", "", y.im, -INFTY, INFTY);
  
  /* Calculate */
  
  z.re = x.re*y.re - x.im*y.im; 
  z.im = x.re*y.im + x.im*y.re; 
  
  return z;
}

double c_norm(complex x)
{
  double norm; 
  
  /* Check values */
  
  CheckValue(FUNCTION_NAME, "x.re (c_norm)", "", x.re, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "x.im (c_norm)", "", x.im, -INFTY, INFTY);
  
  /* Calculate */
  
  norm = x.re*x.re + x.im*x.im;
  norm = sqrt(norm); 
  
  return norm;
}

complex c_div(complex x, complex y)
{
  complex z;
  double n;
  
  /* Check values */
  
  CheckValue(FUNCTION_NAME, "x.re (c_div)", "", x.re, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "x.im (c_div)", "", x.im, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "y.re (c_div)", "", y.re, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "y.im (c_div)", "", y.im, -INFTY, INFTY);
  
  /* Calculate */
  
  n = c_norm(y);
  n = n*n;

  z.re = (x.re*y.re + x.im*y.im)/n;
  z.im = (x.im*y.re - x.re*y.im)/n; 

  /* Check division by zero */

  if (n == 0.0)
    Die(FUNCTION_NAME, "Division by zero (c_div)");
  
  return z;
}

complex c_con(complex x)
{
  complex z;
  
  /* Check values */
  
  CheckValue(FUNCTION_NAME, "x.re (c_conv)", "", x.re, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "x.im (c_conv)", "", x.im, -INFTY, INFTY);
  
  /* Calculate */
  
  z.re = x.re;
  z.im = -x.im; 
  
  return z;
}

complex c_sqrt(complex z)
{
  double x,y,r,theta; 
  double eps; 
  complex w;
  
  /* Check values */
  
  CheckValue(FUNCTION_NAME, "z.re (c_sqrt)", "", z.re, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "z.im (c_sqrt)", "", z.im, -INFTY, INFTY);

  /* Define epsilon */

  eps = 1E-30; 
  
  /* Get real and imaginary part */
  
  x = z.re; 
  y = z.im; 

  /* Avoid compiler warning */

  theta = 0.0;

  /* Calculate r */

  r = sqrt(x*x + y*y); 

  /* Check special cases */

  if (fabs(x) < eps)
    {
      /* Only imaginary part */

      if (y > eps)
        {         
          /* Pi/2 */

          theta = 0.5*PI;
        }
      else if (y < -eps)
        {
          /* -Pi/2 */

          theta = -0.5*PI;
        }
      else if (y > -eps && y < eps)
        { 
          /* 0.0 */
          
          r = 0.0;
          theta = 0.0;
        }
      else
        Die(FUNCTION_NAME, "WTF?");
    }
  else
    {
      /* Real part exists */
  
      if (x > eps)
        {
          /* Right half-plane of C */
          
          theta = atan(y/x);
        }
      else if (x < -eps)
        { 
          if (y > eps)
            { 
              /* second quarter, Pi/2 < theta < Pi  */

              theta = -atan(fabs(y/x)) + PI; 
            }
          else if (y < -eps)
            { 
              /* third quarter, -Pi < theta < -Pi/2 */

              theta = atan(fabs(y/x)) - PI; 
            }
          else if (-eps < y && y < eps)
            { 
              /* y = 0, x < 0 => theta = Pi */
              
              theta = PI; 
            }
          else
            Die(FUNCTION_NAME, "WTF?");
        }
      else
        Die(FUNCTION_NAME, "WTF?");
    }

  /* Calculate real and imaginary parts */

  w.re = sqrt(r)*cos(0.5*theta); 
  w.im = sqrt(r)*sin(0.5*theta); 

  /* Check values */
  
  CheckValue(FUNCTION_NAME, "w.re (c_sqrt)", "", w.re, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "w.im (c_sqrt)", "", w.im, -INFTY, INFTY);

  /* Return value */

  return w;
}

complex c_exp(complex z)
{
  double x,y; 
  complex w;

  /* Check values */
  
  CheckValue(FUNCTION_NAME, "z.re (c_exp)", "", z.re, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "z.im (c_exp)", "", z.im, -INFTY, INFTY);
  
  /* Calculate */
  
  x = z.re; 
  y = z.im; 
  w.re = exp(x)*cos(y); 
  w.im = exp(x)*sin(y); 
  
  return w;
}

/*****************************************************************************/
