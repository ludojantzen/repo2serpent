/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : contribsplit.c                                 */
/*                                                                           */
/* Created:       2019/04/17 (JLe)                                           */
/* Last modified: 2019/04/17 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Sets split flags to adaptive RMX mesh cells based on high    */
/*              contribution,                                                */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ContribSplit:"

/*****************************************************************************/

void ContribSplit(long rmx)
{
  long loc0, loc1, loc2, ng, nmax, n;
  double f, avg, count, val;

  /* Loop over data and reset flags */

  loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (loc0 > VALID_PTR)
    {
      /* Reset common flag */

      WDB[loc0 + RMX_CELL_SPLIT_FLAG] = (double)NO;

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Get number of energy groups */

  ng = (long)RDB[rmx + RMX_NG];
  CheckValue(FUNCTION_NAME, "ng", "", ng, 1, 100000);

  /* Reset sum and number of nonzero elements */

  avg = 0.0;
  count = 0.0;

  /* Loop over mesh */

  loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (loc0 > VALID_PTR)
    {
      /* Reset value */

      val = 0.0;

      /* Get matrix size */

      nmax = (long)RDB[loc0 + RMX_CELL_MTX_SIZE];
      CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

      /* Pointer to forward net currents */

      loc1 = (long)RDB[loc0 + RMX_CELL_MC_CURR_NET_IN];
      CheckPointer(FUNCTION_NAME, "(loc1)", RES2_ARRAY, loc1);

      /* Pointer to current importances  */

      loc2 = (long)RDB[loc0 + RMX_CELL_ADJ_SOL_IN_CURR];
      CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

      /* Add current components */

      for (n = 0; n < nmax*ng; n++)
        if ((f = RES2[loc1 + n]*RDB[loc2 + n]) > 0.0)
          val = val + f;

      /* Pointer to total forward source  */

      loc1 = (long)RDB[loc0 + RMX_CELL_MC_SRCC_TOT];
      CheckPointer(FUNCTION_NAME, "(loc1)", RES2_ARRAY, loc1);

      /* Pointer to source importance */

      loc2 = (long)RDB[loc0 + RMX_CELL_IMP_SRC];
      CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

      /* Add source components */

      for (n = 0; n < nmax*ng; n++)
        if ((f = RES2[loc1 + n]*RDB[loc2 + n]) > 0.0)
          val = val + f;

      /* Add to mean */

      if (val > 0.0)
        {
          avg = avg + val;
          count = count + 1.0;
        }

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Calaulate average */

  if (count == 0.0)
    return;
  else
    avg = avg/count;

  /* Loop over mesh and set split flags */

  loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (loc0 > VALID_PTR)
    {
      /* Reset value */

      val = 0.0;

      /* Get matrix size */

      nmax = (long)RDB[loc0 + RMX_CELL_MTX_SIZE];
      CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

      /* Pointer to forward net currents */

      loc1 = (long)RDB[loc0 + RMX_CELL_MC_CURR_NET_IN];
      CheckPointer(FUNCTION_NAME, "(loc1)", RES2_ARRAY, loc1);

      /* Pointer to current importances  */

      loc2 = (long)RDB[loc0 + RMX_CELL_ADJ_SOL_IN_CURR];
      CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

      /* Add current components */

      for (n = 0; n < nmax*ng; n++)
        if ((f = RES2[loc1 + n]*RDB[loc2 + n]) > 0.0)
          val = val + f;

      /* Pointer to total forward source  */

      loc1 = (long)RDB[loc0 + RMX_CELL_MC_SRCC_TOT];
      CheckPointer(FUNCTION_NAME, "(loc1)", RES2_ARRAY, loc1);

      /* Pointer to source importance */

      loc2 = (long)RDB[loc0 + RMX_CELL_IMP_SRC];
      CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

      /* Add source components */

      for (n = 0; n < nmax*ng; n++)
        if ((f = RES2[loc1 + n]*RDB[loc2 + n]) > 0.0)
          val = val + f;

      /* Compare to mean */

      if (val/avg > 1.0)
        {
          /* Set split flag */

          printf("%E %E\n", avg, val);
          WDB[loc0 + RMX_CELL_SPLIT_FLAG] = (double)YES;
        }

      /* Next */

      loc0 = NextItem(loc0);
    }
}

/*****************************************************************************/
