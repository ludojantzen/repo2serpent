/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : normalizeimp.c                                 */
/*                                                                           */
/* Created:       2019/03/10 (JLe)                                           */
/* Last modified: 2019/10/05 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Normalizes importances in mesh.                              */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "NormalizeImp:"

/*****************************************************************************/

void NormalizeImp(long rmx)
{
  long msh, nr, n, nmax, ptr, ng, wwd;
  double norm, chk, min, max, tot, sum, f, wmax, wmin;
  const double *S;

  fprintf(outp, "Normalizing importances...\n");

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "rmx", DATA_ARRAY, rmx);

  /* Find wwd structure associated with rmx */

  wwd = (long)RDB[DATA_PTR_WWD0];
  while (wwd > VALID_PTR)
    {
      /* Check pointer */

      if ((long)RDB[wwd + WWD_PTR_RMX] == rmx)
        break;

      /* Next */

      wwd = NextItem(wwd);
    }

  /* Get limits */

  wmin = -1.0;
  wmax = -1.0;

  if (wwd > VALID_PTR)
    {
      wmin = RDB[wwd + WWD_LIM_MIN];
      wmax = RDB[wwd + WWD_LIM_MAX];
    }

  /* Check */

  if (wmin < 0.0)
    wmin = 0.0;
  if (wmax < 0.0)
    wmax = INFTY;

  /* Get number of energy groups */

  ng = (long)RDB[rmx + RMX_NG];
  CheckValue(FUNCTION_NAME, "ng", "", ng, 1, 100000);

  /* Reset normalization factor */

  norm = 0.0;
  nr = 0;

  /* Loop over mesh and calculate maximum source importance */

  msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (msh > VALID_PTR)
    {
      /* Pointer to importance */

      ptr = (long)RDB[msh + RMX_CELL_IMP_SRC];
      CheckPointer(FUNCTION_NAME, "(ptr1)", DATA_ARRAY, ptr);

      /* Compare importances */

      for (n = 0; n < ng; n++)
        if (RDB[ptr + n] > norm)
          norm = RDB[ptr + n];

      /* Count nonzero (multiple energy groups as one) */

      for (n = 0; n < ng; n++)
        if (RDB[ptr + n] > 0.0)
          {
            /* Increase count */

            nr++;

            /* Break loop */

            break;
          }

      /* Next */

      msh = NextItem(msh);
    }

  /* Check number of non-zero cells (NOTE: tota rajaa pitää vielä kattoa) */

  if (nr > 20)
    {
      /* Reset total */

      tot = 0.0;

      /* Loop over mesh */

      msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      while (msh > VALID_PTR)
        {
          /* Get source vector */

          ptr = (long)RDB[msh + RMX_CELL_MC_SRCC_TOT];
          CheckPointer(FUNCTION_NAME, "(ptr2)", RES2_ARRAY, ptr);
          S = &RES2[ptr];

          /* Calculate source weighted importance */

          ptr = (long)RDB[msh + RMX_CELL_IMP_SRC];
          CheckPointer(FUNCTION_NAME, "(ptr3)", DATA_ARRAY, ptr);

          sum = 0.0;
          for (n = 0; n < ng; n++)
            sum = sum + RDB[ptr + n]*S[n];

          /* Put value (use empty indicator for temporary space) */

          WDB[msh + RMX_CELL_EMPTY] = sum;
          tot = tot + sum;

          /* Next */

          msh = NextItem(msh);
        }

      /* Check total */

      if (tot == 0.0)
        Die(FUNCTION_NAME, "Total is zero");

      /* Sort array based on importance */

      msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);
      SortList(msh, RMX_CELL_EMPTY, SORT_MODE_ASCEND);

      /* Reset cumulative */

      chk = 0.0;

      /* Loop over mesh */

      msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      while (msh > VALID_PTR)
        {
          /* Add to cumulative */

          chk = chk + RDB[msh + RMX_CELL_EMPTY];

          /* Check fraction (NOTE: tota pitää vielä kattoa) */

          if ((chk/tot > 0.9) && (norm > 0.0))
            break;

          /* Next */

          msh = NextItem(msh);
        }

      /* Check pointer */

      if (msh < VALID_PTR)
        Die(FUNCTION_NAME, "Error in search");

      /* Get maximum group-wise importance */

      ptr = (long)RDB[msh + RMX_CELL_IMP_SRC];
      CheckPointer(FUNCTION_NAME, "(ptr4)", DATA_ARRAY, ptr);

      norm = 0.0;
      for (n = 0; n < ng; n++)
        if (RDB[ptr + n] > norm)
          norm = RDB[ptr + n];

      /* Check */

      if (norm == 0.0)
        Die(FUNCTION_NAME, "Error in normalization");

      /* Loop over mesh */

      msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      while (msh > VALID_PTR)
        {
          /* Reset empty criterion */

          WDB[msh + RMX_CELL_EMPTY] = 0.0;

          /* Next */

          msh = NextItem(msh);
        }

      /* Sort array based on index */

      msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);
      SortList(msh, RMX_CELL_IDX, SORT_MODE_ASCEND);
    }

  fprintf(outp, "OK.\n\n");

  /* Check value */

  CheckValue(FUNCTION_NAME, "norm", "", norm, ZERO, INFTY);

  /* Reset minimum and maximum */

  min = INFTY;
  max = -INFTY;

  /* Normalize distribution */

  msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (msh > VALID_PTR)
    {
      /* Divide importances by normalization factor */

      ptr = (long)RDB[msh + RMX_CELL_IMP_SRC];
      CheckPointer(FUNCTION_NAME, "(ptr5)", DATA_ARRAY, ptr);

      for (n = 0; n < ng; n++)
        {
          /* Normalize */

          f = RDB[ptr + n]/norm;
          CheckValue(FUNCTION_NAME, "f", "", f, 0.0, INFTY);

          /* Compare to limits */

          if (f > wmax)
            WDB[ptr + n] = wmax;
          else if (f < wmin)
            WDB[ptr + n] = wmin;
          else
            WDB[ptr + n] = f;
        }

      ptr = (long)RDB[msh + RMX_CELL_IMP_CURR];
      CheckPointer(FUNCTION_NAME, "(ptr6)", DATA_ARRAY, ptr);

      for (n = 0; n < ng; n++)
        {
          /* Normalize */

          f = RDB[ptr + n]/norm;
          CheckValue(FUNCTION_NAME, "f", "", f, 0.0, INFTY);

          /* Compare to limits */

          if (f > wmax)
            WDB[ptr + n] = wmax;
          else if (f < wmin)
            WDB[ptr + n] = wmin;
          else
            WDB[ptr + n] = f;
        }

      /* Get matrix size */

      nmax = (long)RDB[msh + RMX_CELL_MTX_SIZE];
      CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

      /* Normalize partials */

      ptr = (long)RDB[msh + RMX_CELL_ADJ_SOL_IN_CURR];
      CheckPointer(FUNCTION_NAME, "(ptr7)", DATA_ARRAY, ptr);

      for (n = 0; n < nmax*ng; n++)
        {
          /* Normalize */

          f = RDB[ptr + n]/norm;
          CheckValue(FUNCTION_NAME, "f", "", f, 0.0, INFTY);

          /* Compare to limits */

          if (f > wmax)
            WDB[ptr + n] = wmax;
          else if (f < wmin)
            WDB[ptr + n] = wmin;
          else
            WDB[ptr + n] = f;
        }

      /* Compare source importances to minimum and maximum */

      ptr = (long)RDB[msh + RMX_CELL_IMP_SRC];
      CheckPointer(FUNCTION_NAME, "(ptr8)", DATA_ARRAY, ptr);

      for (n = 0; n < ng; n++)
        {
          if ((RDB[ptr + n] > 0.0) && (RDB[ptr + n] < min))
            min = RDB[ptr + n];

          if (RDB[ptr + n] > max)
            max = RDB[ptr + n];
        }

      /* Next */

      msh = NextItem(msh);
    }

  if (min == wmin)
    fprintf(outp, "Minimum source importance: %1.5E (lim)\n", min);
  else
    fprintf(outp, "Minimum source importance: %1.5E\n", min);

  if (max == wmax)
    fprintf(outp, "Maximum source importance: %1.5E (lim)\n\n", max);
  else
    fprintf(outp, "Maximum source importance: %1.5E\n\n", max);

  /* Reset minimum and maximum */

  min = INFTY;
  max = -INFTY;

  /* Normalize distribution */

  msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (msh > VALID_PTR)
    {
      /* Compare current importances to minimum and maximum */

      ptr = (long)RDB[msh + RMX_CELL_IMP_CURR];
      CheckPointer(FUNCTION_NAME, "(ptr9)", DATA_ARRAY, ptr);

      for (n = 0; n < ng; n++)
        {
          if ((RDB[ptr + n] > 0.0) && (RDB[ptr + n] < min))
            min = RDB[ptr + n];

          if (RDB[ptr + n] > max)
            max = RDB[ptr + n];
        }

      /* Next */

      msh = NextItem(msh);
    }

  if (min == wmin)
    fprintf(outp, "Minimum cell importance: %1.5E (lim)\n", min);
  else
    fprintf(outp, "Minimum cell importance: %1.5E\n", min);

  if (max == wmax)
    fprintf(outp, "Maximum cell importance: %1.5E (lim)\n\n", max);
  else
    fprintf(outp, "Maximum cell importance: %1.5E\n\n", max);

  /* Estimate source sampling efficiency */

  tot = 0.0;
  sum = 0.0;

  /* Loop over mesh */

  msh = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (msh > VALID_PTR)
    {
      /* Get source vector */

      ptr = (long)RDB[msh + RMX_CELL_MC_SRCC_TOT];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      S = &RES2[ptr];

      /* Get source importance importance */

      ptr = (long)RDB[msh + RMX_CELL_IMP_SRC];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Add to source rate and importance-weighted source */

      for (n = 0; n < ng; n++)
        {
          sum = sum + RDB[ptr + n]*S[n];
          tot = tot + S[n];
        }

      /* Next */

      msh = NextItem(msh);
    }

  /* Check */

  if (tot == 0.0)
    Die(FUNCTION_NAME, "Zero source rate");
  else if (sum/tot < 0.02)
    Note(rmx, "Low source sampling efficiency");
}

/*****************************************************************************/
