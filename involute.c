/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : involute.c                                     */
/*                                                                           */
/* Created:       2019/11/07 (JLe)                                           */
/* Last modified: 2019/11/07 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Geometry routines for involute surfaces.                     */
/*                                                                           */
/* Comments: - Based on usersurf.c subroutine written at TUM.                */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "Involute:"

/* Additional functions */

double modulus(double, double);
double newtons_method(double, double, int, double, double, double, double,
                      double, double);
char sign(double);
double f(double, double, double, double, double, double, double);
double df(double, double, double, double, double, double, double);
double x(double, double, double);
double y(double, double, double);
void Schnittstellen_finden(double, double, double, double, double, double,
                           double, double, double, int, int *, double *,
                           double *, double *, double *);

/*****************************************************************************/

void Involute(long part, long np, const double *params, long *in,
              double *dmin, double x, double y, double z, double u, double v,
              double w)
{
  double d, n, l;
  int i_s;
  double x_s1, y_s1, x_s2, y_s2;
  double r0, theta1, theta2, r1, r2, r_point, thetapoint, phi1, phi2;
  int n_max = 1000;
  double e = .00000001;
  int b;

  /* Check parameters pointer */

  if ((np > 0) && (params == NULL))
    Die(FUNCTION_NAME, "Pointer error");

  /* Check call type */

  if (part == 1)
    {
      /***********************************************************************/

      /***** Test if point is inside surface *********************************/

      /* Assign parameters for better readability */

      r0 = params[0];
      theta1 = params[1];
      theta2 = params[2];
      r1 = params[3];
      r2 = params[4];

      r_point = sqrt(x*x + y*y);

      if ((r_point < r1) || (r_point > r2))
        *in = NO;
      else
        {
          if (((x < 0.0) && (y < -r0)) || ((x > 0.0) && (y < r0)))
            b = -1;
          else
            b = 1;

          thetapoint = b*acos((x*sqrt(x*x + y*y - r0*r0) + r0*y)/
                              (x*x + y*y)) - sqrt(x*x + y*y - r0*r0)/r0;

          if (thetapoint < -M_PI)
            thetapoint = thetapoint + 2.0*M_PI;

          /* See whether the point lies between the two involutes or not */

          if ((theta1 <= thetapoint) && (thetapoint < theta2))
            {
              /* Inside */

              *in = YES;
            }
          else
            {
              /* Outside */

              *in = NO;
            }
        }

      /***********************************************************************/
    }
  else if (part == 2)
    {
      /***********************************************************************/

      /***** Calculate minimum distance to surface ***************************/

      /* Reset the minimum distance */

      *dmin = INFTY;

      /* Assign parameters for better readability */

      r0 = params[0];
      theta1 = params[1];
      theta2 = params[2];
      r1 = params[3];
      r2 = params[4];

      phi1 = sqrt((r1*r1)/(r0*r0) - 1.0);
      phi2 = sqrt((r2*r2)/(r0*r0) - 1.0);

      /* länge Bewegungsvektor (u, v, w) */

      l=sqrt(u*u + v*v + w*w);

      /* Distance to Evolute_1 (theta1) */

      Schnittstellen_finden(r0, theta1, x, y, u, v, phi1, phi2, e, n_max,
                            &i_s, &x_s1, &y_s1, &x_s2, &y_s2);

      if (i_s >= 1)
        {
          if (u != 0)
            {
              n = (x_s1 - x)/u;
              d = n*l;

              if ((d > 0) && (d < *dmin))
                *dmin = d;
            }
          else
            {
              if (v != 0)
                {
                  n = (y_s1 - y)/v;
                  d = n*l;

                  if ((d > 0) && (d < *dmin))
                    *dmin = d;
                }
            }

          if (i_s == 2)
            {
              if (u != 0)
                {
                  n = (x_s2 - x)/u;
                  d = n*l;

                  if ((d > 0) && (d < *dmin))
                    *dmin = d;
                }
              else
                {
                  if (v != 0)
                    {
                      n = (y_s2 - y)/v;
                      d = n*l;

                      if ((d > 0) && (d < *dmin))
                        *dmin = d;
                    }

                }
            }
        }

      /* Distance to Evolute_2 (theta2) */

      Schnittstellen_finden(r0, theta2, x, y, u, v, phi1, phi2, e, n_max,
                            &i_s, &x_s1, &y_s1, &x_s2, &y_s2);

      if (i_s >= 1)
        {
          if (u!= 0)
            {
              n = (x_s1 - x)/u;
              d = n*l;

              if ((d > 0) && (d < *dmin))
                *dmin = d;
            }
          else
            {
              if (v!=0)
                {
                  n = (y_s1 - y)/v;
                  d = n*l;
                  if ((d > 0) && (d < *dmin))
                    *dmin = d;
                }
            }
          if (i_s == 2)
            {
              if (u != 0)
                {
                  n = (x_s2 - x)/u;
                  d = n*l;

                  if ((d > 0) && (d < *dmin))
                    *dmin = d;
                }
              else
                {
                  if (v != 0)
                    {
                      n = (y_s2 - y)/v;
                      d = n*l;

                      if ((d > 0) && (d < *dmin))
                        *dmin = d;
                    }
                }
            }
        }

      if ((u != 0) || (v != 0))
        {
          /* Distance to inner Cylinder (r1)*/

	  n = (-u*x - v*y + sqrt(r1*r1*(u*u + v*v) - (v*x - u*y)*(v*x - u*y)))
            /(u*u + v*v);
          d = n*l;

          if ((d > 0) && (d < *dmin))
            *dmin = d;

          n = (-u*x - v*y - sqrt(r1*r1*(u*u + v*v) - (v*x - u*y)*(v*x - u*y)))
            /(u*u+v*v);
          d=n*l;

          if ((d > 0) && (d < *dmin))
            *dmin = d;

          if (*dmin == INFTY)
            {
              /* Distance to outer Cylinder (r2)*/

              n = (-u*x - v*y + sqrt(r2*r2*(u*u + v*v) -
                                     (v*x - u*y)*(v*x - u*y)))/(u*u + v*v);
              d=n*l;

              if ((d > 0) && (d < *dmin))
                *dmin = d;
            }
        }

      /***********************************************************************/
    }
  else
    Die(FUNCTION_NAME, "Invalid part");

  /***************************************************************************/
}

/*****************************************************************************/

double modulus(double a, double b)
{
  int result = a / b;
  return a - (floor(result) * b);
}

/*****************************************************************************/

double newtons_method(double x_start, double e, int n, double r0, double theta,
                      double x0, double y0, double dx, double dy)
{
  double delta = 1000;
  int n_zaehler=0;
  double x = x_start;
  double x1;

  while (delta > e)
    {
      if (n_zaehler > n)
        {
          return 0;
          break;
        }

      x1 = x;
      x = x - f(x, r0, theta, x0, y0, dx, dy)/df(x, r0, theta, x0, y0, dx, dy);
      delta = fabs(x1 - x);
      n_zaehler++;
    }

  return x;
}

/*****************************************************************************/

char sign(double x)
{
  if (x > 0)
    return 1;
  else if(x < 0)
    return -1;
  else
    return 0;
}

/*****************************************************************************/

double f(double phi, double r0, double theta, double x0, double y0,
         double dx, double dy)
{
  return y0*dx + (r0*(phi*cos(phi + theta) - sin(phi + theta)) - x0)*dy -
    r0*(cos(phi + theta) + phi*sin(phi + theta))*dx;
}

/*****************************************************************************/

double df(double phi, double r0, double theta, double x0, double y0,
          double dx, double dy)
{
  return (-1.0)*dx*phi*r0*cos(phi + theta) - dy*phi*r0*sin(phi + theta);
}

/*****************************************************************************/

double x(double phi, double r0, double theta)
{
  return r0*(phi*cos(phi + theta) - sin(phi + theta));
}

/*****************************************************************************/

double y(double phi, double r0, double theta)
{
  return r0*(cos(phi + theta) + phi*sin(phi + theta));
}

/*****************************************************************************/

void Schnittstellen_finden(double r0, double theta, double x0, double y0,
                           double dx, double dy, double phi1, double phi2,
                           double e, int n, int *i_s, double *x_s1,
                           double *y_s1, double *x_s2, double *y_s2)
{
  double phi_s, phi_s1, phi_s2;
  int l = 0, l_z = 0;
  int i_z; /* gesuchte Nullstellen */

  if (sign(f(phi1, r0, theta, x0, y0, dx, dy)) !=
      sign(f(phi2, r0, theta, x0, y0, dx, dy)))

    /* eine NS im Intervall phi1, phi2 */
    {
      i_z = 1;
      phi_s = newtons_method(phi2, e, n, r0, theta, x0, y0, dx, dy);
      if ((phi_s > phi1)  && (phi_s < phi2))
        {
          *i_s = 1;
          *x_s1 = x(phi_s, r0, theta);
          *y_s1 = y(phi_s, r0, theta);
          *x_s2 = 0;
          *y_s2 = 0;
        }
      else
        {
          l = 2;
          /* Unterteilung des intervals phi1 - phi2 in l Teile (als Startwert) */
          while (i_z == 1)
            {
              l_z = 1;
              while (l_z < l)
                {
                  phi_s = newtons_method(phi2-((phi2-phi1)*l_z/l), e, n, r0,
                                         theta, x0, y0, dx, dy);
                  if ((phi_s > phi1) && (phi_s < phi2))
                    {
                      i_z = 0;
                      break;
                    }
                  l_z++;
                }
              l++;
            }

          *i_s = 1;
          *x_s1 = x(phi_s, r0, theta);
          *y_s1 = y(phi_s, r0, theta);
          *x_s2 = 0;
          *y_s2 = 0;
        }
    }
  else /* null oder zwei NS im Intervall phi1, phi2 */
    {
      phi_s1 = newtons_method(phi2, e, n, r0, theta, x0, y0, dx, dy);
      if ((phi_s1 > phi1)  && (phi_s1 < phi2)) /* dann zwei Schnittpunkte */
        {
          i_z = 1; /* ein fehlender Schnittpunkt */

          if (sign(f(phi_s1+e, r0, theta, x0, y0, dx, dy)) ==
              sign(f(phi_s1-e, r0, theta, x0, y0, dx, dy)))
            /* dann beide NS im Rundungsbereich */

            {
              phi_s2 = newtons_method(phi_s1-e, e, n, r0, theta, x0, y0, dx, dy);
            }
          else /* zweite NS außerhalb e */
            {
              phi_s2 = newtons_method(phi1, e, n, r0, theta, x0, y0, dx, dy);

              if ((phi_s2 > phi1) && ((phi_s2 + e) < phi_s1))
                {
                  i_z = 0; /* alle NS gefunden */
                }
              else
                {
                  l = 2;
                  while (i_z == 1)
                    {
                      l_z = 1;
                      while (l_z < l)
                        {
                          phi_s2 =
                            newtons_method(phi_s1 - ((phi_s1 - phi1)*l_z/l),
                                           e, n, r0, theta, x0, y0, dx, dy);

                          if ((phi_s2 > phi1)  && ((phi_s2 + e) < phi_s1))
                            {
                              i_z = 0;
                              break;
                            }
                          l_z++;
                        }
                      l++;
                    }
                }
            }

          *i_s = 2;
          *x_s1 = x(phi_s1, r0, theta);
          *y_s1 = y(phi_s1, r0, theta);
          *x_s2 = x(phi_s2, r0, theta);
          *y_s2 = y(phi_s2, r0, theta);
        }
      else /* keine NS in Intervall phi1, phi2 */
        {
          *i_s = 0;
          *x_s1 = 0;
          *y_s1 = 0;
          *x_s2 = 0;
          *y_s2 = 0;
        }
    }
}

/*****************************************************************************/
