/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : coordtrans.c                                   */
/*                                                                           */
/* Created:       2010/10/07 (JLe)                                           */
/* Last modified: 2019/12/04 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Transfers and rotates coordinates and direction cosines      */
/*                                                                           */
/* Comments: - Converted from rotatecoord.c in Serpent 1.1.14                */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CoordTrans:"

/*****************************************************************************/

void CoordTrans(long ptr, double *x, double *y, double *z, double *u,
                double *v, double *w, long id)
{
  long lvl, n, order, priv;
  double x0, y0, z0, u0, v0, w0, t, t2, vx, vy, vz, dx, dy, dz, r2, bu;
  double ut, vt, wt, cosang, sinang, ang, vel, active;

  /* Check coordinates and direction cosines (HUOM: pretrans.c kutsuu */
  /* tätä rutiinia siten että suuntakosinien paikalla on pintojen eri */
  /* parametreja --> ei voi käyttää normaalia rajatarkastusta. */

  CheckValue(FUNCTION_NAME, "x", "", *x, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "y", "", *y, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "z", "", *z, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "u", "", *u, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "v", "", *v, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "w", "", *w, -INFTY, INFTY);

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Loop over transformations */

  while (ptr > VALID_PTR)
    {
      /* Check burnup transformation (JLe 23.7.2020 / 2.1.32) */

      if (((bu = RDB[ptr + TRANS_BURNUP]) > -INFTY) && (bu < 0.0))
        {
          /* Compare to time */

          if (fabs(-RDB[DATA_BURN_CUM_BURNTIME]/bu/86400.0 - 1.0) > 1E-4)
            {
              /* Skip transformation */

              ptr = NextItem(ptr);

              /* Cycle loop */

              continue;
            }
        }
      else if (bu >= 0.0)
        {
          /* Compare to burnup */

          if (fabs(RDB[DATA_BURN_CUM_BURNUP]/bu - 1.0) > 1E-4)
              {
              /* Skip transformation */

              ptr = NextItem(ptr);

              /* Cycle loop */

              continue;
            }
        }

      /* Check level pointer */

      if ((lvl = (long)RDB[ptr + TRANS_PTR_LVL]) > VALID_PTR)
        {
          /* Get coordinates and direction cosines */

          *x = GetPrivateData(lvl + LVL_PRIV_X, id);
          *y = GetPrivateData(lvl + LVL_PRIV_Y, id);
          *z = GetPrivateData(lvl + LVL_PRIV_Z, id);
          *u = GetPrivateData(lvl + LVL_PRIV_U, id);
          *v = GetPrivateData(lvl + LVL_PRIV_V, id);
          *w = GetPrivateData(lvl + LVL_PRIV_W, id);
        }

      /* Get current time */

      if (RDB[DATA_TDEP_TRANS_TYPE] == (double)TDEP_TRANS_TYPE_INTERVAL)
        t = RDB[DATA_TIME_CUT_TMIN];
      else if (RDB[DATA_TDEP_TRANS_TYPE] == (double)TDEP_TRANS_TYPE_PARTICLE)
        {
          priv = (long)RDB[DATA_PTR_PARTICLE_INITIAL_TIME];
          t = GetPrivateData(priv, id);
        }
      else
        {
          t = 0;
          Die(FUNCTION_NAME, "Unknown transformation time type (set transtime) %ld.",
              (long)RDB[DATA_TDEP_TRANS_TYPE]);
        }

      /* Reset time after end of transformation */

      t2 = 0.0;

      /* Reset active flag */

      active = YES;

      /* Check if time limits are associated with the transformation */

      if (RDB[ptr + TRANS_T_TYPE] != (double)TRANS_TLIM_NONE)
        {
          /* Some time limits are associated with this transformation */

          /* Check initial time */

          if (t < RDB[ptr + TRANS_T0])
            active = NO;
          else if ((t > RDB[ptr + TRANS_T1]) &&
                   (RDB[ptr + TRANS_T_TYPE] == (double)TRANS_TLIM_RESET))
            active = NO;
          else if (t > RDB[ptr + TRANS_T1])
            {
              /* Calculate time after transformation end */

              t2 = t - RDB[ptr + TRANS_T1];

              /* Calculate effect time of transformation */

              t = RDB[ptr + TRANS_T1] - RDB[ptr + TRANS_T0];
            }
          else
            t = t - RDB[ptr + TRANS_T0];

          /* Use post-transformation time only with release */
          /* type transformations */

          if (RDB[ptr + TRANS_T_TYPE] != (double)TRANS_TLIM_RELEASE)
            t2 = 0.0;
        }

      /* Get order */

      order = (long)RDB[ptr + TRANS_ORDER];
      CheckValue(FUNCTION_NAME, "order", "", order, 0, 1);

      /* Loop over order */

      for (n = 0; n < 2; n++)
        {
          /* Check order (0 = translations first, 1 = rotations first) */

          if (n != order)
            {
              /* Check rotations */

              if (RDB[ptr + TRANS_ROT_ANG] != 0.0)
                {
                  /* Rotation wrt general axis */

                  /* Copy values to temporary variables */

                  x0 = *x;
                  y0 = *y;
                  z0 = *z;
                  u0 = *u;
                  v0 = *v;
                  w0 = *w;

                  /* Get unit vector of rotation axis */

                  ut = RDB[ptr + TRANS_ROT_AX_U];
                  vt = RDB[ptr + TRANS_ROT_AX_V];
                  wt = RDB[ptr + TRANS_ROT_AX_W];

                  /* Get rotation angle */

                  if (active == NO)
                    {
                      /* Do not rotate */

                      ang = 0.0;
                    }
                  else if (RDB[ptr + TRANS_MOVE_TYPE] == TRANS_MOVE_NONE)
                    {
                      /* Do transformation */

                      ang = RDB[ptr + TRANS_ROT_ANG];
                    }
                  else if (RDB[ptr + TRANS_MOVE_TYPE] == TRANS_MOVE_VEL)
                    {

                      /* Calculate current angle */

                      ang = RDB[ptr + TRANS_ROT_ANG]*t;
                    }
                  else if (RDB[ptr + TRANS_MOVE_TYPE] == TRANS_MOVE_ACC)
                    {
                      /* Calculate angular velocity collected during */
                      /* acceleration */

                      vel = RDB[ptr + TRANS_ROT_ANG]*t;

                      /* Calculate current angle */

                      ang = 0.5*RDB[ptr + TRANS_ROT_ANG]*t*t + vel*t2;
                    }
                  else
                    ang = 0.0;

                  /* Invert angle since we are transforming neutron */
                  /* and not the geometry */

                  ang = -ang;

                  /* Calculate cosine and sine of current angle */

                  cosang = cos(ang);
                  sinang = sin(ang);

                  /* Create explicit rotation matrix */

                  WDB[ptr + TRANS_RX1] = ut*ut + (vt*vt + wt*wt)*cosang;
                  WDB[ptr + TRANS_RX2] = ut*vt*(1-cosang) - wt*sinang;
                  WDB[ptr + TRANS_RX3] = ut*wt*(1-cosang) + vt*sinang;
                  WDB[ptr + TRANS_RX4] = ut*vt*(1-cosang) + wt*sinang;
                  WDB[ptr + TRANS_RX5] = vt*vt + (ut*ut + wt*wt)*cosang;
                  WDB[ptr + TRANS_RX6] = vt*wt*(1-cosang) - ut*sinang;
                  WDB[ptr + TRANS_RX7] = ut*wt*(1-cosang) - vt*sinang;
                  WDB[ptr + TRANS_RX8] = vt*wt*(1-cosang) + ut*sinang;
                  WDB[ptr + TRANS_RX9] = wt*wt + (ut*ut + vt*vt)*cosang;

                  /* Do translation */

                  x0 = x0 - RDB[ptr + TRANS_X0];
                  y0 = y0 - RDB[ptr + TRANS_Y0];
                  z0 = z0 - RDB[ptr + TRANS_Z0];

                  /* Do matrix multiplication (rotation) */

                  *x = x0*RDB[ptr + TRANS_RX1] + y0*RDB[ptr + TRANS_RX2]
                    + z0*RDB[ptr + TRANS_RX3];
                  *y = x0*RDB[ptr + TRANS_RX4] + y0*RDB[ptr + TRANS_RX5]
                    + z0*RDB[ptr + TRANS_RX6];
                  *z = x0*RDB[ptr + TRANS_RX7] + y0*RDB[ptr + TRANS_RX8]
                    + z0*RDB[ptr + TRANS_RX9];
                  *u = u0*RDB[ptr + TRANS_RX1] + v0*RDB[ptr + TRANS_RX2]
                    + w0*RDB[ptr + TRANS_RX3];
                  *v = u0*RDB[ptr + TRANS_RX4] + v0*RDB[ptr + TRANS_RX5]
                    + w0*RDB[ptr + TRANS_RX6];
                  *w = u0*RDB[ptr + TRANS_RX7] + v0*RDB[ptr + TRANS_RX8]
                    + w0*RDB[ptr + TRANS_RX9];

                  /* Inverse translation */

                  *x = *x + RDB[ptr + TRANS_X0];
                  *y = *y + RDB[ptr + TRANS_Y0];
                  *z = *z + RDB[ptr + TRANS_Z0];

                }
              else if ((long)RDB[ptr + TRANS_ROT] == YES)
                {
                  /* Copy values to temporary variables */

                  x0 = *x;
                  y0 = *y;
                  z0 = *z;
                  u0 = *u;
                  v0 = *v;
                  w0 = *w;

                  /* Do matrix multiplication */

                  *x = x0*RDB[ptr + TRANS_RX1] + y0*RDB[ptr + TRANS_RX2]
                    + z0*RDB[ptr + TRANS_RX3];
                  *y = x0*RDB[ptr + TRANS_RX4] + y0*RDB[ptr + TRANS_RX5]
                    + z0*RDB[ptr + TRANS_RX6];
                  *z = x0*RDB[ptr + TRANS_RX7] + y0*RDB[ptr + TRANS_RX8]
                    + z0*RDB[ptr + TRANS_RX9];
                  *u = u0*RDB[ptr + TRANS_RX1] + v0*RDB[ptr + TRANS_RX2]
                    + w0*RDB[ptr + TRANS_RX3];
                  *v = u0*RDB[ptr + TRANS_RX4] + v0*RDB[ptr + TRANS_RX5]
                    + w0*RDB[ptr + TRANS_RX6];
                  *w = u0*RDB[ptr + TRANS_RX7] + v0*RDB[ptr + TRANS_RX8]
                    + w0*RDB[ptr + TRANS_RX9];

                  /* Check explicit definition */

                  if ((long)RDB[ptr + TRANS_EXPLI] == NO)
                    {
                      /* Copy values to temporary variables */

                      x0 = *x;
                      y0 = *y;
                      z0 = *z;
                      u0 = *u;
                      v0 = *v;
                      w0 = *w;

                      /* Do matrix multiplication */

                      *x = x0*RDB[ptr + TRANS_RY1] + y0*RDB[ptr + TRANS_RY2]
                        + z0*RDB[ptr + TRANS_RY3];
                      *y = x0*RDB[ptr + TRANS_RY4] + y0*RDB[ptr + TRANS_RY5]
                        + z0*RDB[ptr + TRANS_RY6];
                      *z = x0*RDB[ptr + TRANS_RY7] + y0*RDB[ptr + TRANS_RY8]
                        + z0*RDB[ptr + TRANS_RY9];
                      *u = u0*RDB[ptr + TRANS_RY1] + v0*RDB[ptr + TRANS_RY2]
                        + w0*RDB[ptr + TRANS_RY3];
                      *v = u0*RDB[ptr + TRANS_RY4] + v0*RDB[ptr + TRANS_RY5]
                        + w0*RDB[ptr + TRANS_RY6];
                      *w = u0*RDB[ptr + TRANS_RY7] + v0*RDB[ptr + TRANS_RY8]
                        + w0*RDB[ptr + TRANS_RY9];

                      /* Copy values to temporary variables */

                      x0 = *x;
                      y0 = *y;
                      z0 = *z;
                      u0 = *u;
                      v0 = *v;
                      w0 = *w;

                      /* Do matrix multiplication */

                      *x = x0*RDB[ptr + TRANS_RZ1] + y0*RDB[ptr + TRANS_RZ2]
                        + z0*RDB[ptr + TRANS_RZ3];
                      *y = x0*RDB[ptr + TRANS_RZ4] + y0*RDB[ptr + TRANS_RZ5]
                        + z0*RDB[ptr + TRANS_RZ6];
                      *z = x0*RDB[ptr + TRANS_RZ7] + y0*RDB[ptr + TRANS_RZ8]
                        + z0*RDB[ptr + TRANS_RZ9];
                      *u = u0*RDB[ptr + TRANS_RZ1] + v0*RDB[ptr + TRANS_RZ2]
                        + w0*RDB[ptr + TRANS_RZ3];
                      *v = u0*RDB[ptr + TRANS_RZ4] + v0*RDB[ptr + TRANS_RZ5]
                        + w0*RDB[ptr + TRANS_RZ6];
                      *w = u0*RDB[ptr + TRANS_RZ7] + v0*RDB[ptr + TRANS_RZ8]
                        + w0*RDB[ptr + TRANS_RZ9];
                    }
                }
            }
          else
            {
              /* Transfer coordinates */

              if ((active == YES) && (RDB[ptr + TRANS_ROT_ANG] == 0))
                {
                  if (RDB[ptr + TRANS_MOVE_TYPE] == TRANS_MOVE_NONE)
                    {
                      /* No time limits applied to this transformation */

                      *x = *x - RDB[ptr + TRANS_X0];
                      *y = *y - RDB[ptr + TRANS_Y0];
                      *z = *z - RDB[ptr + TRANS_Z0];
                    }
                  else if (RDB[ptr + TRANS_MOVE_TYPE] == TRANS_MOVE_VEL)
                    {

                      /* Do transformation with initial velocity */

                      *x = *x - RDB[ptr + TRANS_X0]*t;
                      *y = *y - RDB[ptr + TRANS_Y0]*t;
                      *z = *z - RDB[ptr + TRANS_Z0]*t;
                    }
                  else if (RDB[ptr + TRANS_MOVE_TYPE] == TRANS_MOVE_ACC)
                    {
                      /* Calculate velocity after end of */
                      /* acceleration */

                      vx = RDB[ptr + TRANS_X0]*t;
                      vy = RDB[ptr + TRANS_Y0]*t;
                      vz = RDB[ptr + TRANS_Z0]*t;

                      /* Do transformation with current velocity */

                      dx = vx*t2;
                      dy = vy*t2;
                      dz = vz*t2;

                      /* Add acceleration transformation during */
                      /* active time */

                      dx += 0.5*RDB[ptr + TRANS_X0]*t*t;
                      dy += 0.5*RDB[ptr + TRANS_Y0]*t*t;
                      dz += 0.5*RDB[ptr + TRANS_Z0]*t*t;

                      /* Do transformation */

                      *x = *x - dx;
                      *y = *y - dy;
                      *z = *z - dz;

                    }
                }
            }
        }

      /* Check coordinates and direction cosines (HUOM: pretrans.c kutsuu */
      /* tätä rutiinia siten että suuntakosinien paikalla on pintojen eri */
      /* parametreja --> ei voi käyttää normaalia rajatarkastusta. */

      CheckValue(FUNCTION_NAME, "x", "", *x, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "y", "", *y, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "z", "", *z, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "u", "", *u, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "v", "", *v, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "w", "", *w, -INFTY, INFTY);

      /* Re-normalize direction cosines */

      r2 = sqrt(*u**u + *v**v + *w**w);

      *u = *u/r2;
      *v = *v/r2;
      *w = *w/r2;

      /* Break loop if lattice transformation */

      if ((long)RDB[ptr + TRANS_TYPE] == TRANS_TYPE_LAT)
        break;

      /* Next */

      ptr = NextItem(ptr);
    }
}

/*****************************************************************************/
