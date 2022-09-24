/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : Boundaryconditions.c                           */
/*                                                                           */
/* Created:       2010/10/10 (JLe)                                           */
/* Last modified: 2019/08/24 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Reflective and periodic boundary conditions by coordinate    */
/*              transformations                                              */
/*                                                                           */
/* Comments: - From Serpent 1.1.12                                           */
/*                                                                           */
/*           - Luodaan erillinen aliohjelma reflection.c niille mieli-       */
/*             valtaisille heijastuksille?                                   */
/*                                                                           */
/*           - Toi cell pointteri on NULL lattice ja pbed -tyyppisille       */
/*             universumeille, mutta niiden ei pitäisi koskaan olla          */
/*             uloimpana, eli tarkistus OK.                                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "BoundaryConditions:"

/*****************************************************************************/

long BoundaryConditions(long *cell, double *x0, double *y0, double *z0,
                        double *u0, double *v0, double *w0, double *xt0,
                        double *yt0, double *zt0, double *wgt, long id)
{
  double x, y, z, u, v, w, px, py, pz, xc, yc, zc, x1, y1, a1, a2, a3, f;
  double xt, yt, zt;
  long surf, param, type, i, j, k, bc0, bc1, bc2, bc3;

  /* Check cell pointer */

  CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, *cell);

  /* Check cell type */

  if ((long)RDB[*cell + CELL_TYPE] != CELL_TYPE_OUTSIDE)
    return NO;

  /* Check black boundary */

  if ((bc0 = (long)RDB[DATA_GEOM_BC0]) == BC_BLACK)
    return -1;

  /* Pointer to surface */

  if ((surf = (long)RDB[*cell + CELL_PTR_BC_SURF]) < VALID_PTR)
    return -1;

  /* Get surface type */

  type = (long)RDB[surf + SURFACE_TYPE];

  /* Get pointer to surface parameters */

  param = (long)RDB[surf + SURFACE_PTR_PARAMS];
  CheckPointer(FUNCTION_NAME, "(param)", DATA_ARRAY, param);

  /* Avoid compiler warning */

  px = 0.0;
  py = 0.0;
  pz = 0.0;

  i = 0;
  j = 0;
  k = 0;

  /* Get partial boundary conditions */

  bc1 = (long)RDB[DATA_GEOM_BC1];
  bc2 = (long)RDB[DATA_GEOM_BC2];
  bc3 = (long)RDB[DATA_GEOM_BC3];

  /* Get coordinates and direction cosines */

  x = *x0;
  y = *y0;
  z = *z0;
  xt = *xt0;
  yt = *yt0;
  zt = *zt0;
  u = *u0;
  v = *v0;
  w = *w0;

  /* Check type */

  if ((type = (long)RDB[surf + SURFACE_TYPE]) == SURF_SQC)
    {
      /***********************************************************************/

      /***** Square lattice **************************************************/

      /* Lattice pitch */

      px = 2.0*RDB[param + 2];
      py = 2.0*RDB[param + 2];
      pz = INFTY;

      /* Centered co-ordinates */

      x = x - RDB[param];
      y = y - RDB[param + 1];

      xt = xt - RDB[param];
      yt = yt - RDB[param + 1];

      /* Get indexes */

      GetLatticeIndexes(px, py, pz, x, y, z, &i, &j, &k, LAT_TYPE_S);

      /* Check leakage from partial boundary conditions */

      if ((bc1 == BC_BLACK) && (i != 0))
        return -1;
      else if ((bc2 == BC_BLACK) && (j != 0))
        return -1;

      /* Calculate new position */

      x = x - i*px;
      y = y - j*py;

      /* Update source coordinates */

      xt = xt - i*px;
      yt = yt - j*py;

      /* Handle Reflection */

      if ((i % 2) && (bc1 == BC_REFLECTIVE))
        {
          /* Odd number of x-surface crossings. Swap vectors. */

          x = -x;
          xt = -xt;
          u = -u;
        }
      if ((j % 2) && (bc2 == BC_REFLECTIVE))
        {
          /* Odd number of y-surface crossings. Swap vectors. */

          y = -y;
          yt = -yt;
          v = -v;
        }

      x = x + RDB[param];
      y = y + RDB[param + 1];

      xt = xt + RDB[param];
      yt = yt + RDB[param + 1];

      /***********************************************************************/
    }
  else if (type == SURF_CUBE)
    {
      /***********************************************************************/

      /***** Cubical 3D-lattice **********************************************/

      /* Lattice pitch */

      px = 2.0*RDB[param + 3];
      py = 2.0*RDB[param + 3];
      pz = 2.0*RDB[param + 3];

      /* Centered co-ordinates */

      x = x - RDB[param];
      y = y - RDB[param + 1];
      z = z - RDB[param + 2];

      xt = xt - RDB[param];
      yt = yt - RDB[param + 1];
      zt = zt - RDB[param + 2];

      /* Get indexes */

      GetLatticeIndexes(px, py, pz, x, y, z, &i, &j, &k, LAT_TYPE_S);

      /* Check leakage from partial boundary conditions */

      if ((bc1 == BC_BLACK) && (i != 0))
        return -1;
      else if ((bc2 == BC_BLACK) && (j != 0))
        return -1;
      else if ((bc3 == BC_BLACK) && (k != 0))
        return -1;

      /* Calculate new position */

      x = x - i*px;
      y = y - j*py;
      z = z - k*pz;

      /* Update source coordinates */

      xt = xt - i*px;
      yt = yt - j*py;
      zt = zt - k*pz;

      /* Handle Reflection */

      if ((i % 2) && (bc1 == BC_REFLECTIVE))
        {
          /* Odd number of x-surface crossings. Swap vectors. */

          x = -x;
          xt = -xt;
          u = -u;
        }
      if ((j % 2) && (bc2 == BC_REFLECTIVE))
        {
          /* Odd number of y-surface crossings. Swap vectors. */

          y = -y;
          yt = -yt;
          v = -v;
        }
      if ((k % 2) && (bc3 == BC_REFLECTIVE))
        {
          /* Odd number of z-surface crossings. Swap vectors. */

          z = -z;
          zt = -zt;
          w = -w;
        }

      x = x + RDB[param];
      y = y + RDB[param + 1];
      z = z + RDB[param + 2];

      xt = xt + RDB[param];
      yt = yt + RDB[param + 1];
      zt = zt + RDB[param + 2];

      /***********************************************************************/
    }
  else if (type == SURF_CUBOID)
    {
      /***********************************************************************/

      /***** Cuboidal 3D-lattice *********************************************/

      /* Get pitches */

      px = (RDB[param + 1] - RDB[param]);
      py = (RDB[param + 3] - RDB[param + 2]);
      pz = (RDB[param + 5] - RDB[param + 4]);

      /* Center co-ordinates */

      xc = (RDB[param] + RDB[param + 1])/2.0;
      yc = (RDB[param + 2] + RDB[param + 3])/2.0;
      zc = (RDB[param + 4] + RDB[param + 5])/2.0;

      /* Transfer co-ordinates */

      x = x - xc;
      y = y - yc;
      z = z - zc;

      xt = xt - xc;
      yt = yt - yc;
      zt = zt - zc;

      /* Get indexes */

      GetLatticeIndexes(px, py, pz, x, y, z, &i, &j, &k, LAT_TYPE_S);

      /* Check leakage from partial boundary conditions */

      if ((bc1 == BC_BLACK) && (i != 0))
        return -1;
      else if ((bc2 == BC_BLACK) && (j != 0))
        return -1;
      else if ((bc3 == BC_BLACK) && (k != 0))
        return -1;

      /* Calculate new position */

      x = x - i*px;
      y = y - j*py;
      z = z - k*pz;

      /* Update source coordinates */

      xt = xt - i*px;
      yt = yt - j*py;
      zt = zt - k*pz;

      /* Handle Reflection */

      if ((i % 2) && (bc1 == BC_REFLECTIVE))
        {
          /* Odd number of x-surface crossings. Swap vectors. */

          x = -x;
          xt = -xt;
          u = -u;
        }
      if ((j % 2) && (bc2 == BC_REFLECTIVE))
        {
          /* Odd number of y-surface crossings. Swap vectors. */

          y = -y;
          yt = -yt;
          v = -v;
        }
      if ((k % 2) && (bc3 == BC_REFLECTIVE))
        {
          /* Odd number of z-surface crossings. Swap vectors. */

          z = -z;
          zt = -zt;
          w = -w;
        }

      /* Transfer co-ordinates */

      x = x + xc;
      y = y + yc;
      z = z + zc;

      xt = xt + xc;
      yt = yt + yc;
      zt = zt + zc;

      /***********************************************************************/
    }
  else if (type == SURF_RECT)
    {
      /***********************************************************************/

      /***** Rectangular 2D-lattice ******************************************/

      /* Get pitches */

      px = (RDB[param + 1] - RDB[param]);
      py = (RDB[param + 3] - RDB[param + 2]);
      pz = INFTY;

      /* Center co-ordinates */

      xc = (RDB[param] + RDB[param + 1])/2.0;
      yc = (RDB[param + 2] + RDB[param + 3])/2.0;

      /* Transfer co-ordinates */

      x = x - xc;
      y = y - yc;

      xt = xt - xc;
      yt = yt - yc;

      /* Get indexes */

      GetLatticeIndexes(px, py, pz, x, y, z, &i, &j, &k, LAT_TYPE_S);

      /* Check leakage from partial boundary conditions */

      if ((bc1 == BC_BLACK) && (i != 0))
        return -1;
      else if ((bc2 == BC_BLACK) && (j != 0))
        return -1;

      /* Calculate new position */

      x = x - i*px;
      y = y - j*py;

      /* Update source coordinates */

      xt = xt - i*px;
      yt = yt - j*py;

      /* Handle Reflection */

      if ((i % 2) && (bc1 == BC_REFLECTIVE))
        {
          /* Odd number of x-surface crossings. Swap vectors. */

          x = -x;
          xt = -xt;
          u = -u;
        }
      if ((j % 2) && (bc2 == BC_REFLECTIVE))
        {
          /* Odd number of y-surface crossings. Swap vectors. */

          y = -y;
          yt = -yt;
          v = -v;
        }

      /* Transfer co-ordinates */

      x = x + xc;
      y = y + yc;

      xt = xt + xc;
      yt = yt + yc;

      /***********************************************************************/
    }
  else if ((type == SURF_HEXYC) || (type == SURF_HEXYPRISM))
    {
      /***********************************************************************/

      /***** Hexagonal Y-type lattice ****************************************/

      /* Check that boundary conditions match */

      if (bc1 != bc2)
        Error(0, "Radial hexagonal boundary conditions must match");

      /* Center co-ordinates */

      xc = RDB[param];
      yc = RDB[param + 1];
      zc = 0.0;

      /* Check */

      if ((xc != 0.0) || (yc != 0.0))
        Error(0, "Only centred hexagonal repeated boundaries allowed");

      /* Lattice pitch */

      px = 2.0*RDB[param + 2];
      py = 2.0*RDB[param + 2];

      /* Z-pitch */

      if (type == SURF_HEXYPRISM)
        {
          pz = (RDB[param + 4] - RDB[param + 3]);
          zc = (RDB[param + 3] + RDB[param + 4])/2.0;
        }
      else
        pz = INFTY;

      /* Transfer */

      x = x - xc;
      y = y - yc;
      z = z - zc;

      xt = xt - xc;
      yt = yt - yc;
      zt = zt - zc;

      /* Get indexes */

      GetLatticeIndexes(px, py, pz, x, y, z, &i, &j, &k, LAT_TYPE_HY);

      /* Axial leakage */

      if ((bc3 == BC_BLACK) && (k != 0))
        return -1;

      /* Calculate new position */

      y = y - (i*py + COS60*j*py);
      x = x - SIN60*j*px;
      z = z - k*pz;

      /* Update source coordinates */

      yt = yt - (i*py + COS60*j*py);
      xt = xt - SIN60*j*px;
      zt = zt - k*pz;

      /* Axial reflection */

      if ((bc3 == BC_REFLECTIVE) && (k % 2))
        {
          z = -z;
          zt = -zt;
          w = -w;
        }

      /* Radial reflection */

      if ((bc1 == BC_REFLECTIVE) && (bc2 == BC_REFLECTIVE))
        {
          if (!(j % 2) && (i % 2))
            {
              y = -y;
              yt = -yt;
              v = -v;
            }
          else if (!(i % 2) && (j % 2))
            {
              y1 = x;
              x1 = y;

              y = x1*COS60 - y1*SIN60;
              x = -x1*SIN60 - y1*COS60;

              y1 = xt;
              x1 = yt;

              yt = x1*COS60 - y1*SIN60;
              xt = -x1*SIN60 - y1*COS60;

              y1 = u;
              x1 = v;

              v = x1*COS60 - y1*SIN60;
              u = -x1*SIN60 - y1*COS60;
            }
          else if ((i % 2) && (j % 2))
            {
              y1 = x;
              x1 = y;

              y = x1*COS60 + y1*SIN60;
              x = x1*SIN60 - y1*COS60;

              y1 = xt;
              x1 = yt;

              yt = x1*COS60 + y1*SIN60;
              xt = x1*SIN60 - y1*COS60;

              y1 = u;
              x1 = v;

              v = x1*COS60 + y1*SIN60;
              u = x1*SIN60 - y1*COS60;
            }
        }

      /* Transfer */

      x = x + xc;
      y = y + yc;
      z = z + zc;

      xt = xt + xc;
      yt = yt + yc;
      zt = zt + zc;

      /***********************************************************************/
    }
  else if ((type == SURF_HEXXC) || (type == SURF_HEXXPRISM))
    {
      /***********************************************************************/

      /***** Hexagonal X-type lattice ****************************************/

      /* Check that boundary conditions match */

      if (bc1 != bc2)
        Error(0, "Radial hexagonal boundary conditions must match");

      /* Center co-ordinates */

      xc = RDB[param];
      yc = RDB[param + 1];
      zc = 0.0;

      /* Check */

      if ((xc != 0.0) || (yc != 0.0))
        Error(0, "Only centred hexagonal repeated boundaries allowed");

      /* Lattice pitch */

      px = 2.0*RDB[param + 2];
      py = 2.0*RDB[param + 2];

      /* Z-pitch */

      if (type == SURF_HEXXPRISM)
        {
          pz = (RDB[param + 4] - RDB[param + 3]);
          zc = (RDB[param + 3] + RDB[param + 4])/2.0;
        }
      else
        pz = INFTY;

      /* Transfer */

      x = x - xc;
      y = y - yc;
      z = z - zc;

      xt = xt - xc;
      yt = yt - yc;
      zt = zt - zc;

      /* Get indexes */

      GetLatticeIndexes(px, py, pz, x, y, z, &i, &j, &k, LAT_TYPE_HX);

      /* Axial leakage */

      if ((bc3 == BC_BLACK) && (k != 0))
        return -1;

      /* Calculate new position */

      x = x - (i*px + COS60*j*px);
      y = y - SIN60*j*py;
      z = z - k*pz;

      /* Update source coordinates */

      xt = xt - (i*px + COS60*j*px);
      yt = yt - SIN60*j*py;
      zt = zt - k*pz;

      /* Axial reflection */

      if ((bc3 == BC_REFLECTIVE) && (k % 2))
        {
          z = -z;
          zt = -zt;
          w = -w;
        }

      /* Radial reflection */

      if ((bc1 == BC_REFLECTIVE) && (bc2 == BC_REFLECTIVE))
        {
          if (!(j % 2) && (i % 2))
            {
              x = -x;
              xt = -xt;
              u = -u;
            }
          else if (!(i % 2) && (j % 2))
            {
              x1 = x;
              y1 = y;

              x = x1*COS60 - y1*SIN60;
              y = -x1*SIN60 - y1*COS60;

              x1 = xt;
              y1 = yt;

              xt = x1*COS60 - y1*SIN60;
              yt = -x1*SIN60 - y1*COS60;

              x1 = u;
              y1 = v;

              u = x1*COS60 - y1*SIN60;
              v = -x1*SIN60 - y1*COS60;
            }
          else if ((i % 2) && (j % 2))
            {
              x1 = x;
              y1 = y;

              x = x1*COS60 + y1*SIN60;
              y = x1*SIN60 - y1*COS60;

              x1 = xt;
              y1 = yt;

              xt = x1*COS60 + y1*SIN60;
              yt = x1*SIN60 - y1*COS60;

              x1 = u;
              y1 = v;

              u = x1*COS60 + y1*SIN60;
              v = x1*SIN60 - y1*COS60;
            }
        }

      /* Transfer */

      x = x + xc;
      y = y + yc;
      z = z + zc;

      xt = xt + xc;
      yt = yt + yc;
      zt = zt + zc;

      /***********************************************************************/
    }
  else
    {
      /* Invalid surface type */

      Die(FUNCTION_NAME, "Invalid surface type");
    }

  /* Check that lattice indexes have been changed */

  if ((i == 0) && (j == 0) && (k == 0))
    {
      /* Check if universe symmetries are defined */

      if ((long)RDB[DATA_PTR_SYM0] < VALID_PTR)
        Die(FUNCTION_NAME, "No change in lattice indexes (cell %s at %E %E %E",
            GetText(*cell + CELL_PTR_NAME), *x0, *y0, *z0);
      else
        Error(0, "Universe symmetry not compatible with boundary conditions");
    }

  /* Check albedo and adjust weight (weight is set to -1 in plotter mode) */

  a1 = RDB[DATA_GEOM_ALBEDO1];
  a2 = RDB[DATA_GEOM_ALBEDO2];
  a3 = RDB[DATA_GEOM_ALBEDO3];

  if (((a1 != 1.0) || (a2 != 1.0) || (a3 != 1.0)) && (*wgt > 0.0))
    {
      /* Indexes allowed to change by +/- 1 only (particle must */
      /* be stopped at the outer boundary) */

      CheckValue(FUNCTION_NAME, "i", "", i, -1, 1);
      CheckValue(FUNCTION_NAME, "j", "", j, -1, 1);
      CheckValue(FUNCTION_NAME, "k", "", k, -1, 1);

      /* Avoid compiler warning */

      f = -1.0;

      /* Calculate changes in weight */

      if (i != 0)
        f = a1;
      else if (j != 0)
        f = a2;
      else if (k != 0)
        f = a3;
      else
        Die(FUNCTION_NAME, "Impossible");

      /* Check value */

      CheckValue(FUNCTION_NAME, "f", "", f, 0.01, 100.0);

      /* Adjust weight */

      *wgt = *wgt*f;
    }

  /* Potential numerical problems */

  if (((i != 0) && (fabs(px/2.0 - x) < 1E-15)) ||
      ((j != 0) && (fabs(py/2.0 - y) < 1E-15)) ||
      ((k != 0) && (fabs(pz/2.0 - z) < 1E-15)))
    Die(FUNCTION_NAME, "Infinite loop? (position too close to boudary)");

  /* Find location (ei testata koska muuten kosahtaa plotterimoodissa) */

  *cell = WhereAmI(x, y, z, u, v, w, id);

  /* Put coordinates and direction cosines */

  *x0 = x;
  *y0 = y;
  *z0 = z;

  *xt0 = xt;
  *yt0 = yt;
  *zt0 = zt;

  *u0 = u;
  *v0 = v;
  *w0 = w;

  /* Exit subroutine */

  return YES;
}

/*****************************************************************************/
