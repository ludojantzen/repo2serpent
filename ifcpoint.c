/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : ifcpoint.c                                     */
/*                                                                           */
/* Created:       2013/03/17 (JLe)                                           */
/* Last modified: 2019/03/29 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Retrieves material density factor and temperature at a point */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "IFCPoint:"

/*****************************************************************************/

void IFCPoint(long mat, double *f0, double *T0, double t, long id)
{
  long loc0, loc1, msh, type, dim, ptr, lst, np, uni, ncol;
  long ang, nr, rad, i, tmplist0, tmplist, dflist, prnt, idx;
  double f, T, wgt, r, r2, d2, d, ex, w, x, y, z, dx, dy, dz;
  double t0, t1, fBOI, TBOI;
  double Temp1, Temp0, phi, phi2;

  /* Check if interfaces are defined */

  if ((long)RDB[DATA_PTR_IFC0] < VALID_PTR)
    return;

  /* Check material pointer */

  if (mat < VALID_PTR)
    return;

  /* Check pointer to material-wise interface */

  if ((loc0 = (long)RDB[mat + MATERIAL_PTR_IFC]) > VALID_PTR)
    {
      /* Get pointer to root universe */

      uni = (long)RDB[DATA_PTR_ROOT_UNIVERSE];
      CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

      /* Get interface type */

      type = (long)RDB[loc0 + IFC_TYPE];
    }
  else
    {
      /* Get collision universe */

      ptr = (long)RDB[DATA_PTR_COLLISION_UNI];
      uni = (long)GetPrivateData(ptr, id);

      CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

      /* Check pointer to universe-wise interface */

      if ((loc0 = (long)RDB[uni + UNIVERSE_PTR_IFC_FUEP]) < VALID_PTR)
        return;

      /* Check whether the point lies in a nest region with material */
      /* without TMS treatment */

      if (RDB[mat + MATERIAL_TMS_MODE] == TMS_MODE_NONE)
        return;

      /* Put type */

      type = (long)RDB[loc0 + IFC_FUEP_TYPE];
    }

  /* Get coordinates */

  ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_X];
  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
  x = GetPrivateData(ptr, id);

  ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_Y];
  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
  y = GetPrivateData(ptr, id);

  ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_Z];
  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
  z = GetPrivateData(ptr, id);

  /* Get time from collision universe as it is not put to root universe */

  /* Get collision universe */

  ptr = (long)RDB[DATA_PTR_COLLISION_UNI];
  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
  ptr = (long)GetPrivateData(ptr, id);

  /* Get time */

  if (t < 0.0)
    {
      ptr = (long)RDB[ptr + UNIVERSE_PTR_PRIVA_T];
      CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
      t = GetPrivateData(ptr, id);
    }

  /* Reset density factor and temperature */

  f = -1.0;
  T = -1.0;

  /* Check type */

  if (type == IFC_TYPE_PT_AVG)
    {
      /***********************************************************************/

      /***** Average of point-wise values ************************************/

      /* Get pointer to search mesh */

      msh = (long)RDB[loc0 + IFC_PTR_SEARCH_MESH_LIST];
      CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

      /* Get dimensions */

      dim = (long)RDB[loc0 + IFC_DIM];

      /* Reset mean density factor, temperature and weight */

      f = 0.0;
      T = 0.0;
      wgt = 0.0;

      /* Get exclusion radius and square of radius */

      r = RDB[loc0 + IFC_EXCL_RAD];
      r2 = r*r;

      /* Get exponent */

      ex = RDB[loc0 + IFC_EXP];

      /* Get pointer to search mesh */

      if ((lst = MeshPtr(msh, x, y, z)) > VALID_PTR)
        lst = (long)RDB[lst];

      /* Loop over content */

      while (lst > VALID_PTR)
        {
          /* Pointer to point */

          loc1 = (long)RDB[lst + SEARCH_MESH_CELL_CONTENT];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          /* Get parameters */

          dx = x - RDB[loc1 + IFC_PT_X];
          dy = y - RDB[loc1 + IFC_PT_Y];
          dz = z - RDB[loc1 + IFC_PT_Z];

          /* Avoid compiler warning */

          d2 = -1.0;

          /* Calculate square distance */

          if (dim == 3)
            d2 = dx*dx + dy*dy + dz*dz;
          else if (dim == 2)
            d2 = dx*dx + dy*dy;
          else if (dim == 1)
            d2 = dz*dz;
          else
            Die(FUNCTION_NAME, "Invalid dimension");

          /* Compare to exclusion radius */

          if (d2 < r2)
            {
              /* Calculate distance */

              d = sqrt(d2);

              /* Calculate weight factor */

              w = pow(d, ex);
              CheckValue(FUNCTION_NAME, "w", "", w, 0.0, INFTY);

              /* Invert */

              if (w < 1E-3)
                w = 1E+3;
              else
                w = 1.0/w;

              /* Add to values */

              f = f + RDB[loc1 + IFC_PT_DF]*w;
              T = T + RDB[loc1 + IFC_PT_TMP]*w;
              wgt = wgt + w;
            }

          /* Next */

          lst = NextItem(lst);
        }

      /* Calculate mean */

      if (wgt != 0.0)
        {
          f = f/wgt;
          T = T/wgt;
        }

      /* Avoid round-off errors */

      if (T > RDB[loc0 + IFC_MAX_TEMP])
        T = RDB[loc0 + IFC_MAX_TEMP];
      else if (T < RDB[loc0 + IFC_MIN_TEMP])
        T = RDB[loc0 + IFC_MIN_TEMP];

      /***********************************************************************/
    }
  else if ((type == IFC_TYPE_REG_MESH) || (type == IFC_TYPE_REG_MESH_MULTIMAT)
           || (type == IFC_TYPE_REG_MESH_MULTILVL))
    {
      /***********************************************************************/

      /***** Regular mesh based distribution *********************************/

      /* Get pointer to mesh */

      msh = (long)RDB[loc0 + IFC_PTR_SEARCH_MESH_LIST];
      CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

      /* Use dedicated subroutine to find index */

      idx = FindRegMeshIFCIndex(loc0, x, y, z, id);

      /* Get temperature and density if we found a mesh cell */

      if (idx >= 0)
        {
          /* Get pointer to temperature list */

          tmplist = (long)RDB[loc0 + IFC_PTR_TMP_LIST];

          /* Get pointer to density list */

          dflist = (long)RDB[loc0 + IFC_PTR_DF_LIST];

          /* Get values */

          f = RDB[dflist + idx];
          T = RDB[tmplist + idx];
          /*
          fprintf(outp, "Returning from %ld: %f %f\n", i, f, T);
          */
          /* Check for time interpolation */

          if (((tmplist = (long)RDB[loc0 + IFC_PTR_TMP_LIST_BOI]) > VALID_PTR)
              && ((dflist = (long)RDB[loc0 + IFC_PTR_DF_LIST_BOI]) > VALID_PTR))
            {

              /* Get time-interval beginning time */

              t0 = RDB[DATA_TIME_CUT_TMIN];

              /* Get time-interval end time */

              t1 = RDB[DATA_TIME_CUT_TMAX];

              /* Get BOI density factor */

              fBOI = RDB[dflist + idx];

              /* Interpolate density factor */

              f = fBOI + (f - fBOI)*(t - t0)/(t1 - t0);

              /* Get BOI temperature */

              TBOI = RDB[tmplist + idx];

              /* Interpolate temperature */

              T = TBOI + (T - TBOI)*(t - t0)/(t1 - t0);

#ifdef DEBUG
             /* Check that the interpolated temperature is inside TMS limits */

              if ((T < RDB[mat + MATERIAL_TMS_TMIN]) ||
                  (T > RDB[mat + MATERIAL_TMS_TMAX]))
                {
                  /* Print some kind of a warning */

                  fprintf(errp, "Time-interpolated temperature out of TMS-limits for material %s:\n",
                          GetText(mat + MATERIAL_PTR_NAME));
                  fprintf(errp, "Limits in Kelvin: [%E, %E]\n",
                          RDB[mat + MATERIAL_TMS_TMIN], RDB[mat + MATERIAL_TMS_TMAX]);
                  fprintf(errp, "BOI    : time %E s, temperature %E K, density factor %E\n",
                          t0, TBOI, fBOI);

                  tmplist = (long)RDB[loc0 + IFC_PTR_TMP_LIST];
                  dflist = (long)RDB[loc0 + IFC_PTR_DF_LIST];
                  fprintf(errp, "EOI    : time %E s, temperature %E K, density factor %E\n",
                          t1, RDB[tmplist + idx], RDB[dflist + idx]);
                  fprintf(errp, "Current: time %E s, temperature %E K, density factor %E\n",
                          t, T, f);
                }
#endif
            }
        }

      /***********************************************************************/
    }
  else if (type == IFC_TYPE_FUNC)
    {
      /***********************************************************************/

      /***** User-defined functional dependence ******************************/

      /* Get number of parameters */

      np = (long)RDB[loc0 + IFC_FUNC_NP];

      /* Get pointer to parameters */

      ptr = (long)RDB[loc0 + IFC_FUNC_PTR_PARAM];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get density factor and temperature */

      UserIFC(mat, &f, &T, x, y, z, t, np, &RDB[ptr]);

      /***********************************************************************/
    }
  else if ((type == IFC_TYPE_FET_DENSITY) || (type == IFC_TYPE_FET_TEMP))
    {
      /***********************************************************************/

      /***** Functional expansions *******************************************/

      /* Get pointer to parameters */

      ptr = (long)RDB[loc0 + IFC_FET_INPUT_PARAMS_PTR];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get pointer to coefficient list */

      lst = (long)RDB[loc0 + IFC_FET_PTR_INPUT_COEFS_ARRAY];
      CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);

      /* Get density factor or temperature, depending on 'type' */

      FEIFC(&RDB[ptr], &RDB[lst], mat, &f, &T, x, y, z, type, id);

      /***********************************************************************/
    }
  else if (type == IFC_TYPE_TET_MESH)
    {
      /***********************************************************************/

      /***** Unstructured tetrahedral mesh ***********************************/

      /* If density file is not given, the density factor will be 1.0 */
      /* even if the tet-cell is not found */

      if (RDB[loc0 + IFC_PTR_OF_RFILE] < VALID_PTR)
        f = 1.0;

      /* Find region */

      if ((loc1 = FindTetCell(loc0, x, y, z, id)) > VALID_PTR)
        {
          /* Get collision number */

          ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
          ncol = (long)GetPrivateData(ptr, id);

          /* Put collision flag */

          StoreValuePair(loc0 + IFC_PTR_PREV_COL_CELL, (double)ncol,
                         (double)loc1, id);

          /* Get parent pointer */

          prnt = (long)RDB[loc1 + TET_PTR_PARENT];
          CheckPointer(FUNCTION_NAME, "(prnt)", DATA_ARRAY, prnt);

          /* This is the index of the parent cell  */

          i = (long)RDB[prnt + IFC_TET_PRNT_IDX];

          /* Get pointer to temperature list */

          tmplist = (long)RDB[loc0 + IFC_PTR_TMP_LIST];

          /* Get pointer to density list */

          dflist = (long)RDB[loc0 + IFC_PTR_DF_LIST];

          /* Get values */

          f = RDB[dflist + i];
          T = RDB[tmplist + i];

          /* Check for time interpolation */

          if (((tmplist = (long)RDB[loc0 + IFC_PTR_TMP_LIST_BOI]) > VALID_PTR)
              && ((dflist = (long)RDB[loc0 + IFC_PTR_DF_LIST_BOI]) > VALID_PTR))
            {

              /* Get time-interval beginning time */

              t0 = RDB[DATA_TIME_CUT_TMIN];

              /* Get time-interval end time */

              t1 = RDB[DATA_TIME_CUT_TMAX];

              /* Get BOI density factor */

              fBOI = RDB[dflist + i];

              /* Interpolate density factor */

              f = fBOI + (f - fBOI)*(t - t0)/(t1 - t0);

              /* Get BOI temperature factor */

              TBOI = RDB[tmplist + i];

              /* Interpolate density factor */

              T = TBOI + (T - TBOI)*(t - t0)/(t1 - t0);

#ifdef DEBUG
             /* Check that the interpolated temperature is inside TMS limits */

              if ((T < RDB[mat + MATERIAL_TMS_TMIN]) ||
                  (T > RDB[mat + MATERIAL_TMS_TMAX]))
                {
                  /* Print some kind of a warning */

                  fprintf(errp, "Time-interpolated temperature out of TMS-limits for material %s:\n",
                          GetText(mat + MATERIAL_PTR_NAME));
                  fprintf(errp, "Limits in Kelvin: [%E, %E]\n",
                          RDB[mat + MATERIAL_TMS_TMIN], RDB[mat + MATERIAL_TMS_TMAX]);
                  fprintf(errp, "BOI    : time %E s, temperature %E K, density factor %E\n",
                          t0, TBOI, fBOI);

                  tmplist = (long)RDB[loc0 + IFC_PTR_TMP_LIST];
                  dflist = (long)RDB[loc0 + IFC_PTR_DF_LIST];
                  fprintf(errp, "EOI    : time %E s, temperature %E K, density factor %E\n",
                          t1, RDB[tmplist + i], RDB[dflist + i]);
                  fprintf(errp, "Current: time %E s, temperature %E K, density factor %E\n",
                          t, T, f);
                }
#endif
            }

        }

      /***********************************************************************/
    }
  else if ((type == IFC_TYPE_FUEP) || (type == IFC_TYPE_FPIP))
    {
      /***********************************************************************/

      /***** Interface to fuel performance codes *****************************/

      /* Reset density factor and temperature */

      *f0 = -1.0;
      *T0 = -1.0;

      loc1 = (long)RDB[loc0 + IFC_FUEP_PTR_AX];

      /* Loop over axial zones */

      while (loc1 > VALID_PTR)
        {
          /* Compare coordinates */

          if ((z >= RDB[loc1 + IFC_FUEP_AX_ZMIN]) &&
              (z < RDB[loc1 + IFC_FUEP_AX_ZMAX]))
            break;

          /* Next */

          loc1 = NextItem(loc1);
        }

      /* Check pointer */

      if (loc1 < VALID_PTR)
        {
          /*
          if(RDB[mat + MATERIAL_USE_IFC] == YES)
            Warn(FUNCTION_NAME,"Ax not found");
          */
          return;
        }

      /* Find angular zone */

      ang = (long)RDB[loc1 + IFC_FUEP_AX_PTR_ANG];

      /* Get polar angle */

      phi = PolarAngle(x,y);

      while (ang > VALID_PTR)
        {

          /* Rotate if needed */
          if(phi > 2.0*PI+RDB[ang + IFC_FUEP_ANG_AMIN])
            phi2 = phi - 2.0*PI;
          else
            phi2 = phi;

          /* Compare coordinates */

          if ((phi2 >= RDB[ang + IFC_FUEP_ANG_AMIN]) &&
              (phi2 < RDB[ang + IFC_FUEP_ANG_AMAX]))
            break;

          /* Next */

          ang = NextItem(ang);
        }

      /* Check if found */

      if (ang < VALID_PTR)
        {
          if(RDB[mat + MATERIAL_USE_IFC] == YES)
            Warn(FUNCTION_NAME,"Ang not found");
          return;
        }

      nr = (long)RDB[ang + IFC_FUEP_ANG_N_RAD];

      /* Calculate square radius */

      r2 = x*x + y*y;

      /* Find radial zone, use BOI geometry */

      rad = (long)RDB[ang + IFC_FUEP_ANG_PTR_HOT_R2_BOI];

      CheckPointer(FUNCTION_NAME, "(rad)", DATA_ARRAY, rad);

      i = SearchArray(&RDB[rad], r2, nr);

      if (i < 0)
        {
          if(RDB[mat + MATERIAL_USE_IFC] == YES)
            {

              Warn(FUNCTION_NAME, "Rad not found %E", sqrt(r2));
              for(i=0; i < nr; i++)
                {
                  fprintf(outp, "%E ", sqrt(RDB[rad + i]));
                }
              fprintf(outp, "\n nr = %ld \n", nr);
            }
          return;
        }

      /**************************************/
      /* Set density factor and temperature */
      /**************************************/

      /* Use BOI density factor since BOI geometry is also used */

      dflist = (long)RDB[ang + IFC_FUEP_ANG_PTR_DF_BOI];
      CheckPointer(FUNCTION_NAME, "(dflist)", DATA_ARRAY, dflist);

      /* Get density factor */

      f = RDB[dflist + i + 1];

      /* Get pointer to EOI temperature list */

      tmplist = (long)RDB[ang + IFC_FUEP_ANG_PTR_TMP];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get EOI temperature */

      T = RDB[tmplist + i + 1];

      /* Check if BOI temperature list exists */

      if ((tmplist0 = (long)RDB[ang + IFC_FUEP_ANG_PTR_TMP_BOI]) == tmplist)

        {
          /* Is equal to the EOI list, do not interpolate */

          if ((type == IFC_TYPE_FPIP) && (i >= 0))
            {

              T = RDB[tmplist + i] +
                (T - RDB[tmplist + i])*
                (sqrt(r2) - sqrt(RDB[rad + i]))/
                (sqrt(RDB[rad + i + 1]) - sqrt(RDB[rad + i]));

            }

        }
      else
        {
          /* BOI temperature list exists */

          /* Get time-interval beginning time */

          t0 = RDB[DATA_TIME_CUT_TMIN];

          /* Get time-interval end time */

          t1 = RDB[DATA_TIME_CUT_TMAX];

          /* Get BOI temperature */

          TBOI = RDB[tmplist0 + i + 1];

          /* Interpolate between timesteps */

          T = TBOI + (T - TBOI)*
            (t - t0)/(t1 - t0);

          /* Interpolate with FPIP */

          if((type == IFC_TYPE_FPIP) && (i >= 0))
            {

              /* Interpolate outer temperature wrt time */

              Temp1 = RDB[tmplist0 + i + 1] +
                (RDB[tmplist + i + 1] - RDB[tmplist0 + i + 1])*
                (t - t0)/(t1 - t0);

              /* Interpolate inner temperature wrt time */

              Temp0 = RDB[tmplist0 + i] +
                (RDB[tmplist + i] - RDB[tmplist0 + i])*
                (t - t0)/(t1 - t0);

              /* Interpolate the both wrt radius */

              T = Temp0 +
                (Temp1 - Temp0)*
                (sqrt(r2) - sqrt(RDB[rad + i]))/
                (sqrt(RDB[rad + i + 1]) - sqrt(RDB[rad + i]));

            }
        }

      /***********************************************************************/
    }
  else
    Die(FUNCTION_NAME, "Invalid interface type");

  /* Put density factor */

  *f0 = f;

  /* Check temperature */

  if ((T < ZERO) && (RDB[mat + MATERIAL_TMS_MODE] != TMS_MODE_NONE))
    {
      /* Use TMS minimum (this is used if the point */
      /* lies outside the distribution) */

      T = RDB[mat + MATERIAL_TMS_TMIN];
      CheckValue(FUNCTION_NAME, "T", "", T, ZERO, 1E+12);

    }

  /* Put temperature */

  *T0 = T;
}

/*****************************************************************************/
