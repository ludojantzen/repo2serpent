/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : pretrans.c                                     */
/*                                                                           */
/* Created:       2017/02/01 (JLe)                                           */
/* Last modified: 2019/12/21 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Convert certain transformed surface types to new surfaces    */
/*              to speed-up the calculation.                                 */
/*                                                                           */
/* Comments: - Toimiikohan tää ihan oikein branch-laskennan kanssa?          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PreTrans:"

/*****************************************************************************/

void PreTrans()
{
  long surf, tra, type, ptr, i, rot;
  double *parm, x, y, z, u, v, w, l, r, T[12];

  /* Loop over surfaces */

  surf = (long)RDB[DATA_PTR_S0];
  while (surf > VALID_PTR)
    {
      /* Pointer to transformation */

      if (((tra = (long)RDB[surf + SURFACE_PTR_TRANS]) < VALID_PTR) ||
          ((long)RDB[tra + TRANS_EXPLI] == NO) ||
          (NextItem(tra) > VALID_PTR) ||
          (RDB[tra + TRANS_MOVE_TYPE] != TRANS_MOVE_NONE))
        {
          /* Next surface */

          surf = NextItem(surf);

          /* Cycle loop */

          continue;
        }

      /* Get type */

      type = (long)RDB[surf + SURFACE_TYPE];

      /* Check */

      if ((type != SURF_SPH) &&
          (type != SURF_CYLX) &&
          (type != SURF_CYLY) &&
          (type != SURF_CYLZ) &&
          (type != SURF_CYLV) &&
          (type != SURF_PX) &&
          (type != SURF_PY) &&
          (type != SURF_PZ) &&
          (type != SURF_PLANE) &&
          (type != SURF_RCC))
        {
          /* Next surface */

          surf = NextItem(surf);

          /* Cycle loop */

          continue;
        }

      /* Read data */

      for (i = 0; i < 12; i++)
        T[i] = RDB[tra + TRANS_X0 + i];

      /* Set rotation flag */

      if ((T[3] == 1.0) && (T[4] == 0.0) && (T[5] == 0.0) &&
          (T[6] == 0.0) && (T[7] == 1.0) && (T[8] == 0.0) &&
          (T[9] == 0.0) && (T[10] == 0.0) && (T[11] == 1.0))
        rot = NO;
      else
        rot = YES;

      /* Invert position vector */

      WDB[tra + TRANS_X0] = -T[0];
      WDB[tra + TRANS_Y0] = -T[1];
      WDB[tra + TRANS_Z0] = -T[2];

      /* Transpose matrix */

      WDB[tra + TRANS_RX1] = T[3];
      WDB[tra + TRANS_RX2] = T[6];
      WDB[tra + TRANS_RX3] = T[9];
      WDB[tra + TRANS_RX4] = T[4];
      WDB[tra + TRANS_RX5] = T[7];
      WDB[tra + TRANS_RX6] = T[10];
      WDB[tra + TRANS_RX7] = T[5];
      WDB[tra + TRANS_RX8] = T[8];
      WDB[tra + TRANS_RX9] = T[11];

      /* Switch order */

      if ((long)RDB[tra + TRANS_ORDER] == TRANS_ORDER_ROT_TRANS)
        WDB[tra + TRANS_ORDER] = TRANS_ORDER_TRANS_ROT;
      else
        WDB[tra + TRANS_ORDER] = TRANS_ORDER_ROT_TRANS;

      /* Reset pointer */

      WDB[surf + SURFACE_PTR_TRANS] = -1.0;

      /* Get pointer to surface parameters */

      if ((ptr = (long)RDB[surf + SURFACE_PTR_PARAMS]) > VALID_PTR)
        parm = &WDB[ptr];
      else
        parm = NULL;

      /* Check type */

      switch (type)
        {
          /*******************************************************************/

          /***** Sphere ******************************************************/

        case SURF_SPH:
          {
            /* Set parameters */

            x = parm[0];
            y = parm[1];
            z = parm[2];
            u = 1.0;
            v = 0.0;
            w = 0.0;

            /* Apply transformation */

            CoordTrans(tra, &x, &y, &z, &u, &v, &w, 0);

            /* Put parameters */

            parm[0] = x;
            parm[1] = y;
            parm[2] = z;

            /* Break case */

            break;
          }

          /*******************************************************************/

          /***** Infinite cylinders ******************************************/

        case SURF_CYLX:
          {
            /* Check rotation */

            if (rot == NO)
              {
                /* Transfer position */

                parm[0] = parm[0] - T[1];
                parm[1] = parm[1] - T[2];
              }
            else
              {
                /* Check number of parameters */

                if ((long)RDB[surf + SURFACE_N_PARAMS] == 3)
                  {
                    /* Set parameters for general equation (CYLV) */

                    x = 0.0;
                    y = parm[0];
                    z = parm[1];
                    u = 1.0;
                    v = 0.0;
                    w = 0.0;
                    r = parm[2];

                    /* Apply transformation */

                    CoordTrans(tra, &x, &y, &z, &u, &v, &w, 0);

                    /* Change type */

                    WDB[surf + SURFACE_TYPE] = (double)SURF_CYLV;

                    /* Put number of parameters */

                    WDB[surf + SURFACE_N_PARAMS] = 7.0;

                    /* Allocate memory for parameters */

                    ptr = ReallocMem(DATA_ARRAY, 7);
                    WDB[surf + SURFACE_PTR_PARAMS] = (double)ptr;

                    /* Put parameters */

                    WDB[ptr++] = x;
                    WDB[ptr++] = y;
                    WDB[ptr++] = z;
                    WDB[ptr++] = u;
                    WDB[ptr++] = v;
                    WDB[ptr++] = w;
                    WDB[ptr++] = r;
                  }
                else
                  {
                    /* Set parameters for general equation (RCC) */

                    x = -parm[4];
                    y = parm[0];
                    z = parm[1];
                    u = 1.0;
                    v = 0.0;
                    w = 0.0;
                    r = parm[2];

                    /* Apply transformation */

                    CoordTrans(tra, &x, &y, &z, &u, &v, &w, 0);

                    /* Change type */

                    WDB[surf + SURFACE_TYPE] = (double)SURF_RCC;

                    /* Put number of parameters */

                    WDB[surf + SURFACE_N_PARAMS] = 9.0;

                    /* Allocate memory for parameters */

                    ptr = ReallocMem(DATA_ARRAY, 9);
                    WDB[surf + SURFACE_PTR_PARAMS] = (double)ptr;

                    /* Put parameters */

                    WDB[ptr++] = x;
                    WDB[ptr++] = y;
                    WDB[ptr++] = z;
                    WDB[ptr++] = u;
                    WDB[ptr++] = v;
                    WDB[ptr++] = w;
                    WDB[ptr++] = r;

                    /* Get length */

                    l = parm[4] - parm[3];

                    /* Re-calculate position of truncation planes */

                    WDB[ptr] = u*x + v*y + w*z;
                    WDB[ptr + 1] = WDB[ptr] + l;
                  }
              }

            /* Break case */

            break;
          }
        case SURF_CYLY:
          {
            /* Check rotation */

            if (rot == NO)
              {
                /* Transfer position */

                parm[0] = parm[0] - T[0];
                parm[1] = parm[1] - T[2];
              }
            else
              {
                /* Check number of parameters */

                if ((long)RDB[surf + SURFACE_N_PARAMS] == 3)
                  {
                    /* Set parameters for general equation (CYLV) */

                    x = parm[0];
                    y = 0.0;
                    z = parm[1];
                    u = 0.0;
                    v = 1.0;
                    w = 0.0;
                    r = parm[2];

                    /* Apply transformation */

                    CoordTrans(tra, &x, &y, &z, &u, &v, &w, 0);

                    /* Change type */

                    WDB[surf + SURFACE_TYPE] = (double)SURF_CYLV;

                    /* Put number of parameters */

                    WDB[surf + SURFACE_N_PARAMS] = 7.0;

                    /* Allocate memory for parameters */

                    ptr = ReallocMem(DATA_ARRAY, 7);
                    WDB[surf + SURFACE_PTR_PARAMS] = (double)ptr;

                    /* Put parameters */

                    WDB[ptr++] = x;
                    WDB[ptr++] = y;
                    WDB[ptr++] = z;
                    WDB[ptr++] = u;
                    WDB[ptr++] = v;
                    WDB[ptr++] = w;
                    WDB[ptr++] = r;
                  }
                else
                  {
                    /* Set parameters for general equation (RCC) */

                    x = parm[0];
                    y = -parm[4];
                    z = parm[1];
                    u = 0.0;
                    v = 1.0;
                    w = 0.0;
                    r = parm[2];

                    /* Apply transformation */

                    CoordTrans(tra, &x, &y, &z, &u, &v, &w, 0);

                    /* Change type */

                    WDB[surf + SURFACE_TYPE] = (double)SURF_RCC;

                    /* Put number of parameters */

                    WDB[surf + SURFACE_N_PARAMS] = 9.0;

                    /* Allocate memory for parameters */

                    ptr = ReallocMem(DATA_ARRAY, 9);
                    WDB[surf + SURFACE_PTR_PARAMS] = (double)ptr;

                    /* Put parameters */

                    WDB[ptr++] = x;
                    WDB[ptr++] = y;
                    WDB[ptr++] = z;
                    WDB[ptr++] = u;
                    WDB[ptr++] = v;
                    WDB[ptr++] = w;
                    WDB[ptr++] = r;

                    /* Get length */

                    l = parm[4] - parm[3];

                    /* Re-calculate position of truncation planes */

                    WDB[ptr] = u*x + v*y + w*z;
                    WDB[ptr + 1] = WDB[ptr] + l;
                  }
              }

            /* Break case */

            break;
          }
        case SURF_CYLZ:
          {
            /* Check rotation */

            if (rot == NO)
              {
                /* Transfer position */

                parm[0] = parm[0] - T[0];
                parm[1] = parm[1] - T[1];
              }
            else
              {
                /* Check number of parameters */

                if ((long)RDB[surf + SURFACE_N_PARAMS] == 3)
                  {
                    /* Set parameters for general equation (CYLV) */

                    x = parm[0];
                    y = parm[1];
                    z = 0.0;
                    u = 0.0;
                    v = 0.0;
                    w = 1.0;
                    r = parm[2];

                    /* Apply transformation */

                    CoordTrans(tra, &x, &y, &z, &u, &v, &w, 0);

                    /* Change type */

                    WDB[surf + SURFACE_TYPE] = (double)SURF_CYLV;

                    /* Put number of parameters */

                    WDB[surf + SURFACE_N_PARAMS] = 7.0;

                    /* Allocate memory for parameters */

                    ptr = ReallocMem(DATA_ARRAY, 7);
                    WDB[surf + SURFACE_PTR_PARAMS] = (double)ptr;

                    /* Put parameters */

                    WDB[ptr++] = x;
                    WDB[ptr++] = y;
                    WDB[ptr++] = z;
                    WDB[ptr++] = u;
                    WDB[ptr++] = v;
                    WDB[ptr++] = w;
                    WDB[ptr++] = r;
                  }
                else
                  {
                    /* Set parameters for general equation (RCC) */

                    x = parm[0];
                    y = parm[1];
                    z = -parm[4];
                    u = 0.0;
                    v = 0.0;
                    w = 1.0;
                    r = parm[2];

                    /* Apply transformation */

                    CoordTrans(tra, &x, &y, &z, &u, &v, &w, 0);

                    /* Change type */

                    WDB[surf + SURFACE_TYPE] = (double)SURF_RCC;

                    /* Put number of parameters */

                    WDB[surf + SURFACE_N_PARAMS] = 9.0;

                    /* Allocate memory for parameters */

                    ptr = ReallocMem(DATA_ARRAY, 9);
                    WDB[surf + SURFACE_PTR_PARAMS] = (double)ptr;

                    /* Put parameters */

                    WDB[ptr++] = x;
                    WDB[ptr++] = y;
                    WDB[ptr++] = z;
                    WDB[ptr++] = u;
                    WDB[ptr++] = v;
                    WDB[ptr++] = w;
                    WDB[ptr++] = r;

                    /* Get length */

                    l = parm[4] - parm[3];

                    /* Re-calculate position of truncation planes */

                    WDB[ptr] = u*x + v*y + w*z;
                    WDB[ptr + 1] = WDB[ptr] + l;
                  }
              }

            /* Break case */

            break;
          }
        case SURF_CYLV:
          {
            /* Set parameters */

            x = parm[0];
            y = parm[1];
            z = parm[2];
            u = parm[3];
            v = parm[4];
            w = parm[5];
            r = parm[6];

            /* Apply transformation */

            CoordTrans(tra, &x, &y, &z, &u, &v, &w, 0);

            /* Put parameters */

            parm[0] = x;
            parm[1] = y;
            parm[2] = z;
            parm[3] = u;
            parm[4] = v;
            parm[5] = w;
            parm[6] = r;

            /* Break case */

            break;
          }

          /*******************************************************************/

          /***** Planes ******************************************************/

        case SURF_PX:
          {
            /* Check rotation */

            if (rot == NO)
              parm[0] = parm[0] + T[0];
            else
              {
                /* Set parameters for general equation */

                x = parm[0];
                y = 0.0;
                z = 0.0;
                u = 1.0;
                v = 0.0;
                w = 0.0;

                /* Apply transformation */

                CoordTrans(tra, &x, &y, &z, &u, &v, &w, 0);

                /* Change type */

                WDB[surf + SURFACE_TYPE] = (double)SURF_PLANE;

                /* Put number of parameters */

                WDB[surf + SURFACE_N_PARAMS] = 4.0;

                /* Allocate memory for parameters */

                ptr = ReallocMem(DATA_ARRAY, 4);
                WDB[surf + SURFACE_PTR_PARAMS] = (double)ptr;

                /* Put parameters */

                WDB[ptr++] = u;
                WDB[ptr++] = v;
                WDB[ptr++] = w;
                WDB[ptr++] = u*x + v*y + w*z;
              }

            /* Break case */

            break;
          }

        case SURF_PY:
          {
            /* Check rotation */

            if (rot == NO)
              parm[0] = parm[0] + T[1];
            else
              {
                /* Set parameters for general equation */

                x = 0.0;
                y = parm[0];
                z = 0.0;
                u = 0.0;
                v = 1.0;
                w = 0.0;

                /* Apply transformation */

                CoordTrans(tra, &x, &y, &z, &u, &v, &w, 0);

                /* Change type */

                WDB[surf + SURFACE_TYPE] = (double)SURF_PLANE;

                /* Put number of parameters */

                WDB[surf + SURFACE_N_PARAMS] = 4.0;

                /* Allocate memory for parameters */

                ptr = ReallocMem(DATA_ARRAY, 4);
                WDB[surf + SURFACE_PTR_PARAMS] = (double)ptr;

                /* Put parameters */

                WDB[ptr++] = u;
                WDB[ptr++] = v;
                WDB[ptr++] = w;
                WDB[ptr++] = u*x + v*y + w*z;
              }

            /* Break case */

            break;
          }

        case SURF_PZ:
          {
            /* Check rotation */

            if (rot == NO)
              parm[0] = parm[0] + T[2];
            else
              {
                /* Set parameters for general equation */

                x = 0.0;
                y = 0.0;
                z = parm[0];
                u = 0.0;
                v = 0.0;
                w = 1.0;

                /* Apply transformation */

                CoordTrans(tra, &x, &y, &z, &u, &v, &w, 0);

                /* Change type */

                WDB[surf + SURFACE_TYPE] = (double)SURF_PLANE;

                /* Put number of parameters */

                WDB[surf + SURFACE_N_PARAMS] = 4.0;

                /* Allocate memory for parameters */

                ptr = ReallocMem(DATA_ARRAY, 4);
                WDB[surf + SURFACE_PTR_PARAMS] = (double)ptr;

                /* Put parameters */

                WDB[ptr++] = u;
                WDB[ptr++] = v;
                WDB[ptr++] = w;
                WDB[ptr++] = u*x + v*y + w*z;
              }

            /* Break case */

            break;
          }

        case SURF_PLANE:
          {
            /* Check number of parameters */

            if ((long)RDB[surf + SURFACE_N_PARAMS] != 4)
              Die(FUNCTION_NAME, "Unable to proceed");

            /* Set parameters */

            u = parm[0];
            v = parm[1];
            w = parm[2];

            if (u != 0.0)
              {
                x = parm[3]/u;
                y = 0.0;
                z = 0.0;
              }
            else if (v != 0.0)
              {
                x = 0.0;
                y = parm[3]/v;
                z = 0.0;
              }
            else if (w != 0.0)
              {
                x = 0.0;
                y = 0.0;
                z = parm[3]/w;
              }
            else
              Die(FUNCTION_NAME, "Error in definition");

            /* Apply transformation */

            CoordTrans(tra, &x, &y, &z, &u, &v, &w, 0);

            /* Put parameters */

            parm[0] = u;
            parm[1] = v;
            parm[2] = w;
            parm[3] = u*x + v*y + w*z;

            /* Break case */

            break;
          }

          /*******************************************************************/

          /***** MCNP-type truncated cylinde *********************************/

        case SURF_RCC:
          {
            /* Get position and direction */

            x = parm[0];
            y = parm[1];
            z = parm[2];
            u = parm[3];
            v = parm[4];
            w = parm[5];

            /* Get length */

            l = parm[8] - parm[7];

            /* Apply transformation */

            CoordTrans(tra, &x, &y, &z, &u, &v, &w, 0);

            /* Put coordinates */

            parm[0] = x;
            parm[1] = y;
            parm[2] = z;
            parm[3] = u;
            parm[4] = v;
            parm[5] = w;

            /* Re-calculate position of truncation planes */

            parm[7] = u*x + v*y + w*z;
            parm[8] = parm[7] + l;

            /* Break case */

            break;
          }

          /*******************************************************************/

        default:
          {
            Die(FUNCTION_NAME, "Invalid type");
          }
        }

      /* Next surface */

      surf = NextItem(surf);
    }

}

/******************************************************************************/
