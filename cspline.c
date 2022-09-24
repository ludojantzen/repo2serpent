/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : cspline.c                                      */
/*                                                                           */
/* Created:       2014/11/11 (TKa)                                           */
/* Last modified: 2017/03/09 (TKa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Functions for cubic spline interpolation                     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"


/* TODO: Better function descriptions */

/* Local function definitions */
static double CSplineF(double **coeff, double x, long idx);


/*****************************************************************************/

void CSplineConstruct(const double *x0, const double *f0, long n0, double d2f0,
                     double d2fend, double **coeff) {
  /* Construct cubic splines (uses Thomas algorithm for solving tridiagonal
   *  system)
   * TODO: better description,  A d2f = D
   * */
  static char * const FUNCTION_NAME = "CSplineConstruct:";
  long i, usestatic;
  double *d2f, *h, *d, *H, *D, *u, *v;
  double Hhu, d2f11, d2f12, d2f21, d2f22, h2i, h6i, hp6;
  double fd2f[CSPLINE_FAST_SZ];
  double fh[CSPLINE_FAST_SZ];
  double fd[CSPLINE_FAST_SZ];
  double fH[CSPLINE_FAST_SZ];
  double fD[CSPLINE_FAST_SZ];
  double fu[CSPLINE_FAST_SZ];
  double fv[CSPLINE_FAST_SZ];


  /* Check the array size */
  if (n0 < 3)
    Die(FUNCTION_NAME, "CSplineConstruct: Number of elements too small: %ld", n0);

  /* Check the array is sorted */
  for (i = 1; i < n0; i++)
    if (x0[i-1] > x0[i])
      Die(FUNCTION_NAME, "CSplineConstruct: Input array not sorted");

  usestatic = (n0 < CSPLINE_FAST_SZ);

  if (usestatic) {
    /* Use static arrays for speed up*/
    d2f = fd2f;
    h = fh;
    d = fd;
    H = fH;
    D = fD;
    u = fu;
    v = fv;
  }
  else {
    /* Allocate memory */
    d2f = (double *)Mem(MEM_ALLOC, n0, sizeof(double));
    h = (double *)Mem(MEM_ALLOC, n0, sizeof(double));
    d = (double *)Mem(MEM_ALLOC, n0, sizeof(double));
    H = (double *)Mem(MEM_ALLOC, n0, sizeof(double));
    D = (double *)Mem(MEM_ALLOC, n0, sizeof(double));
    u = (double *)Mem(MEM_ALLOC, n0, sizeof(double));
    v = (double *)Mem(MEM_ALLOC, n0, sizeof(double));
  }

  /* Calculate the elements of the tridiagonal matrix (h and H) and vector d */
  for (i = 0; i < n0-1; i++) {
    h[i] = x0[i+1] - x0[i];
    d[i] = (f0[i+1] - f0[i])/h[i];
  }

  for (i = 1; i < n0-1; i++)
    H[i] = 2.0*(h[i-1] + h[i]);

  /* Set the second derivative of the first and last points */
  d2f[0] = d2f0;
  d2f[n0-1] = d2fend;

  /* More stuff */
  D[1] = 6.0*(d[1] - d[0]) - h[0]*d2f[0];
  D[n0-2] = 6.0*(d[n0-2] - d[n0-3]) - h[n0-2]*d2f[n0-1];

  for (i = 2; i < n0-2; i++)
    D[i] = 6.0*(d[i] - d[i-1]);

  /* Forward substitution */
  u[1] = h[1]/H[1];
  v[1] = D[1]/H[1];
  for (i = 2; i < n0-1; i++) {
    if ((Hhu = H[i] - h[i-1]*u[i-1]) == 0.) {
      /* This should not happen! */
      Die(FUNCTION_NAME, "CSplineConstruct: H[i] - h[i-1]*u[i-1] is zero");
    }
    u[i] = h[i]/Hhu;
    v[i] = (D[i] - h[i-1]*v[i-1])/Hhu;
  }

  /* Back substitution */
  d2f[n0-2] = v[n0-2];
  for (i = n0-3; i > 0; i--) {
    d2f[i] = v[i] - u[i]*d2f[i+1];
  }

  /* Calculate spline coefficients */
  for (i = 0; i < n0-1; i++) {
    d2f11 = d2f[i]*x0[i+1];
    d2f12 = d2f11*x0[i+1];
    d2f21 = d2f[i+1]*x0[i];
    d2f22 = d2f21*x0[i];
    h2i = 1.0/(2.0*h[i]);
    h6i = h2i/3.0;
    hp6 = h[i]/6.0;

    coeff[i][0] = (d2f12*x0[i+1] - d2f22*x0[i] + 6.0*(f0[i]*x0[i+1]
                  - f0[i+1]*x0[i]))*h6i + hp6*(d2f21 - d2f11);
    coeff[i][1] = (d2f22 - d2f12 + 2.0*(f0[i+1] - f0[i]))
                  *h2i + hp6*(d2f[i] - d2f[i+1]);
    coeff[i][2] = (d2f11 - d2f21)*h2i;
    coeff[i][3] = (d2f[i+1] - d2f[i])*h6i;
  }


  if (!usestatic) {
    /* Free memory */
    Mem(MEM_FREE, d2f);
    Mem(MEM_FREE, h);
    Mem(MEM_FREE, d);
    Mem(MEM_FREE, H);
    Mem(MEM_FREE, D);
    Mem(MEM_FREE, u);
    Mem(MEM_FREE, v);
  }

}

/*****************************************************************************/


/*****************************************************************************/

void CSplineInterpolate(double **coeff, const double *x0, long n0,
                        const double *x, double *f, long n) {
  /* Interpolation using cubic splines.
   *
   * Arguments:
   * - coeff :  spline coefficients (to be used with x0)
   * - x0    :  sample points
   * - n0    :  size of x0
   * - x     :  query points
   * - f     :  query values (calculated)
   * - n     :  size of x and f
   *
   * NOTE:
   * - Argument 'x' must be in ascending order
   * */
  static char * const FUNCTION_NAME = "CSplineInterpolate:";
  long i, j;

  /* Check the array size */
  if (n0 < 3)
    Die(FUNCTION_NAME, "CSplineInterpolate: Number of elements too small: %ld", n0);

  j = 0;

  /* Loop over x0 values */
  for (i = 0; i < n; i++) {

    /* Check boundaries */
    if ((x[i] < x0[0]) || (x[i] > x0[n0-1]))
      Die(FUNCTION_NAME, "CSplineInterpolate: x0=%E out of bounds", x[i]);

    /* Find correct interval */
    while (j < n0) {
      if (x[i] <= x0[j+1])
        break;
      if ((i > 0) && (x[i] < x[i-1]))
        Die(FUNCTION_NAME, "CSplineInterpolate: x is not sorted");
      j++;
    }

    /* Interpolate */
    f[i] = coeff[j][0] + x[i]*(coeff[j][1] + x[i]*(coeff[j][2]
            + x[i]*coeff[j][3]));
  }

}

/*****************************************************************************/


/*****************************************************************************/

void CSplineInterpolate0(const double *x0, const double *f0, long n0,
                         double d2f0, double d2fend, const double *x, double *f,
                         long n, long mode) {
  /* Wrapper for interpolation when spline coefficients etc. are not needed.
   * NOTE:
   * - Parameter 'mode' defines the interpolation mode the same way as in
   *   function ENDFInterp().
   * */
  static char * const FUNCTION_NAME = "CSplineInterpolate0:";
  long i;
  double *lx0, *lf0, *lx;
  double **coeff;

  /* Check the array size */
  if (n0 < 3)
    Die(FUNCTION_NAME, "CSplineInterpolate0: Number of elements too small: %ld", n0);

  /* Allocate memory for the coefficients */
  coeff = (double **)Mem(MEM_ALLOC, n0-1, sizeof(double*));
  for (i = 0; i < n0-1; i++)
    coeff[i] = (double *)Mem(MEM_ALLOC, 4, sizeof(double));


  /* Select interpolation mode */

  switch (mode) {
    case 0:
    case 2:
      {
        /* Interpolation in x0 and f0 */

        /* Construct splines */
        CSplineConstruct(x0, f0, n0, d2f0, d2fend, coeff);

        /* Interpolate */
        CSplineInterpolate(coeff, x0, n0, x, f, n);

        break;
      }
    case 3:
      {
        /* Interpolation in log(x0) and f0 */

        /* Allocate memory */
        lx0 = (double *)Mem(MEM_ALLOC, n0, sizeof(double));
        for (i = 0; i < n0; i++)
          lx0[i] = log(x0[i]);
        lx = (double *)Mem(MEM_ALLOC, n, sizeof(double));
        for (i = 0; i < n; i++)
          lx[i] = log(x[i]);

        /* Construct splines */
        CSplineConstruct(lx0, f0, n0, d2f0, d2fend, coeff);

        /* Interpolate */
        CSplineInterpolate(coeff, lx0, n0, lx, f, n);

        /* Free memory */
        Mem(MEM_FREE, lx0);
        Mem(MEM_FREE, lx);

        break;
      }
    case 4:
      {
        /* Interpolation in x0 and log(f0) */

        /* Allocate memory */
        lf0 = (double *)Mem(MEM_ALLOC, n0, sizeof(double));
        for (i = 0; i < n0; i++)
          lf0[i] = log(f0[i]);

        /* Construct splines */
        CSplineConstruct(x0, lf0, n0, d2f0, d2fend, coeff);

        /* Interpolate */
        CSplineInterpolate(coeff, x0, n0, x, f, n);

        /* log(f) -> f */
        for (i = 0; i < n; i++)
          f[i] = exp(f[i]);

        /* Free memory */
        Mem(MEM_FREE, lf0);

        break;
      }
    case 5:
      {
        /* Interpolation in log(x0) and log(f0) */

        /* Allocate memory */
        lx0 = (double *)Mem(MEM_ALLOC, n0, sizeof(double));
        lf0 = (double *)Mem(MEM_ALLOC, n0, sizeof(double));
        for (i = 0; i < n0; i++) {
          lx0[i] = log(x0[i]);
          lf0[i] = log(f0[i]);
        }
        lx = (double *)Mem(MEM_ALLOC, n, sizeof(double));
        for (i = 0; i < n; i++)
          lx[i] = log(x[i]);

        /* Construct splines */
        CSplineConstruct(lx0, lf0, n0, d2f0, d2fend, coeff);

        /* Interpolate */
        CSplineInterpolate(coeff, lx0, n0, lx, f, n);

        /* log(f) -> f */
        for (i = 0; i < n; i++)
          f[i] = exp(f[i]);

        /* Free memory */
        Mem(MEM_FREE, lx0);
        Mem(MEM_FREE, lf0);
        Mem(MEM_FREE, lx);

        break;
      }
    default:
      Die(FUNCTION_NAME, "CSplineInterpolate0: Invalid interpolation mode");
  }

  /* Free memory */
  for (i = 0; i < n0-1; i++)
    Mem(MEM_FREE, coeff[i]);
  Mem(MEM_FREE, coeff);

}

/*****************************************************************************/


/*****************************************************************************/

double CSplineF(double **coeff, double x, long idx) {
  /* Returns the spline integral.
   * */
  return x*(coeff[idx][0] + x*(coeff[idx][1]/2.0 + x*(coeff[idx][2]/3.0
         + x*coeff[idx][3]/4.0)));
}

/*****************************************************************************/


/*****************************************************************************/

double CSplineIntegrate(double **coeff, const double *x0, double *F, long n0,
                     double xa, double xb) {
  /* Integration using cubic splines
   * */
  static char * const FUNCTION_NAME = "CSplineIntegrate:";
  long ia, ib;
  double Fa, Fb, Fab;

  /* Check the array size */
  if (n0 < 3)
    Die(FUNCTION_NAME, "CSplineIntegrate: Number of elements too small: %ld", n0);

  /* Find interval for xa */
  if ((ia = SearchArray(x0, xa, n0)) == -1) {
    if (xa == x0[n0-1])
      ia = n0 - 2;
    else
      Die(FUNCTION_NAME, "CSplineIntegrate: xa %.5E not found\n", xa);
  }

  /* Find interval for xb */
  if ((ib = SearchArray(x0, xb, n0)) == -1) {
    if (xb == x0[n0-1])
      ib = n0 - 2;
    else
      Die(FUNCTION_NAME, "CSplineIntegrate: xb %.5E not found\n", xb);
  }

  /* Calculate integral */
  Fa = CSplineF(coeff, xa, ia);
  Fb = CSplineF(coeff, xb, ib);

  if (ia == ib) {
    Fab = Fb - Fa;
  }
  else {
    Fa = CSplineF(coeff, x0[ia+1], ia) - Fa;
    Fb = Fb - CSplineF(coeff, x0[ib], ib);
    Fab = F[ib-1] - F[ia] + Fa + Fb;
  }

  return Fab;
}

/*****************************************************************************/


/*****************************************************************************/

double CSplineIntegrate0(const double *x0, double *f0, long n0, double d2f0,
                         double d2fend, double xa, double xb) {
  /* Wrapper for integration when spline coefficients etc. are not needed
   * */
  static char * const FUNCTION_NAME = "CSplineIntegrate0:";
  long i;
  double Fab;
  double *F0;
  double **coeff;

  /* Check the array size */
  if (n0 < 3)
    Die(FUNCTION_NAME, "Number of elements too small: %ld", n0);

  /* Allocate memory for the coefficients */
  coeff = (double **)Mem(MEM_ALLOC, n0-1, sizeof(double*));
  for (i = 0; i < n0-1; i++)
    coeff[i] = (double *)Mem(MEM_ALLOC, 4, sizeof(double));

  /* Allocate memory for integral array */
  F0 = (double *)Mem(MEM_ALLOC, n0, sizeof(double));

  /* Construct splines */
  CSplineConstruct(x0, f0, n0, d2f0, d2fend, coeff);

  /* Integrate */
  Fab = CSplineIntegrate(coeff, x0, F0, n0, xa, xb);

  /* Free memory */
  for (i = 0; i < n0-1; i++)
    Mem(MEM_FREE, coeff[i]);
  Mem(MEM_FREE, coeff);
  Mem(MEM_FREE, F0);

  return Fab;
}

/*****************************************************************************/


/*****************************************************************************/

void CSplineCumIntegral(double **coeff, const double *x0, long n0,
                        const double *x, double *F, long n) {
  /* Calculates the cumulative integral of the spline at given points.
   * The first value F[0] is set to zero.
   * */
  static char * const FUNCTION_NAME = "CSplineCumIntegral:";
  long i, j;
  double Fa, Fb;

  /* Check the array size */
  if (n0 < 3)
    Die(FUNCTION_NAME, "CSplineCumIntegral: Number of elements too"
                       " small: %ld", n0);

  /* Calculate the integral */
  F[0] = 0.;
  j = 0;

  /* Loop over x-values */
  for (i = 1; i < n; i++) {

    /* Find correct interval */
    while (j < n0) {
      if (x[i] <= x0[j+1])
        break;
      if ((i > 0) && (x[i] < x[i-1]))
        Die(FUNCTION_NAME, "CSplineCumIntegral: x is not sorted");
      j++;
    }

    Fa = CSplineF(coeff, x[i-1], j);
    Fb = CSplineF(coeff, x[i], j);
    F[i] = F[i-1] + Fb - Fa;
  }
}

/*****************************************************************************/


/*****************************************************************************/

void CSplineCumIntegral0(const double *x, const double *f, double *F, long n,
                         double d2f0, double d2fend, long mode,
                         const double *lxflag, double **coeff) {
  /* Wrapper for integration when spline coefficients etc. are not needed
   * */
  static char * const FUNCTION_NAME = "CSplineCumIntegral0:";
  long i, j;
  double Fa, Fb, a1, b1, c1, a0, b0, c0, dF;
  double *lx;

  /* Check the array size */
  if (n < 3)
    Die(FUNCTION_NAME, "CSplineCumIntegral0: Number of elements too"
                       " small: %ld", n);

  /* Select interpolation mode */
  switch (mode) {
    case 0:
    case 2:
      {
        /* Interpolation in x0 and f0 */

        /* Construct splines */
        CSplineConstruct(x, f, n, d2f0, d2fend, coeff);

        /* Calculate the integral */
        F[0] = 0.;
        j = 0;

        /* Loop over x-values */
        for (i = 1; i < n; i++) {
          Fa = CSplineF(coeff, x[i-1], i-1);
          Fb = CSplineF(coeff, x[i], i-1);
          F[i] = F[i-1] + Fb - Fa;
        }

        break;
      }
    case 3:
      {
        /* Interpolation in log(x0) and f0 */

        /* Allocate memory and create log arrays if needed */
        if (lxflag) {
          lx = (double *)lxflag;
        }
        else {
          lx = (double *)Mem(MEM_ALLOC, n, sizeof(double));
          for (i = 0; i < n; i++)
            lx[i] = log(x[i]);
        }

        /* Construct splines */
        CSplineConstruct(lx, f, n, d2f0, d2fend, coeff);

        F[0] = 0.;
        j = 0;
        i = 0;
        a1 = (lx[i] - 1.0)*x[i];
        b1 = lx[i]*(a1 + x[i]) - 2.0*a1;
        c1 = (b1 + 2.0*a1)*(lx[i] - 3.0) + 6.0*a1;

        /* Loop over x-values */
        for (i = 1; i < n; i++) {
          j = i-1;
          a0 = a1;
          b0 = b1;
          c0 = c1;
          a1 = (lx[i] - 1.0)*x[i];
          b1 = lx[i]*(a1 + x[i]) - 2.0*a1;
          c1 = (b1 + 2.0*a1)*(lx[i] - 3.0) + 6.0*a1;

          dF = coeff[j][0]*(x[i] - x[i-1])
               + coeff[j][1]*(a1 - a0)
               + coeff[j][2]*(b1 - b0)
               + coeff[j][3]*(c1 - c0);

          F[i] = F[i-1] + dF;

        }

        /* Free memory if needed */
        if (!lxflag)
          Mem(MEM_FREE, lx);

        break;
      }
    case 4:
      {
        /* Interpolation in x0 and log(f0) */
        Die(FUNCTION_NAME, "CSplineInterpolate0: log-lin interpolation"
                           " not implemented");
        break;
      }
    case 5:
      {
        /* Interpolation in log(x0) and log(f0) */
        Die(FUNCTION_NAME, "CSplineInterpolate0: log-log interpolation not"
                           " supported");
        break;
      }
    default:
      Die(FUNCTION_NAME, "CSplineInterpolate0: Invalid interpolation mode");
  }

}

/*****************************************************************************/
