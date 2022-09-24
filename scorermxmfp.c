/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scorermxmfp.c                                 */
/*                                                                           */
/* Created:       2018/07/21 (JLe)                                           */
/* Last modified: 2018/10/09 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Scores flux and total xs for response matrix solver          */
/*                                                                           */
/* Comments: - Used with adaptive mesh                                       */
/*           - Käytetään tätä nyt noiden tiheyksien proubaamiseen            */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreRMXMFP:"

/*****************************************************************************/

void ScoreRMXMFP(long part, long mat, double g, double flx, double tot,
                 double wgt, long id)
{
  long rmx, loc0, loc1, itp, ptr, i, splitmax;
  long rhmax, rho;

  /* Check if mfp iteration */

  if ((long)RDB[DATA_RMTX_MFP_CALC] == NO)
    return;

  /* Check material pointer */

  if (mat < VALID_PTR)
    return;

  /* Pointer to weight window structure */

  rmx = (long)RDB[DATA_PTR_RMX0];
  CheckPointer(FUNCTION_NAME, "(rmx)", DATA_ARRAY, rmx);

  /* Get pointer to cell */

  loc0 = (long)RDB[part + PARTICLE_ICM_PTR_ICM];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Maximum number of splits */

  splitmax = (long)RDB[DATA_MAX_RMX_SPLIT_FLAGS];
  CheckValue(FUNCTION_NAME, "splitmax", "", splitmax, 1, 100);

  /* Loop over split flags */

  ptr = (long)RDB[loc0 + RMX_CELL_PTR_SPLIT];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  for (i = 0; i < splitmax; i++)
    if ((long)RDB[ptr++] == YES)
      return;

  /* Get pointer to iteration data */

  itp = (long)RDB[rmx + RMX_PTR_ITER];
  CheckPointer(FUNCTION_NAME, "(itp)", DATA_ARRAY, itp);

  /* Pointer to split vector */

  ptr = (long)RDB[loc0 + RMX_CELL_PTR_SPLIT];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Loop over split criteria */

  loc1 = (long)RDB[itp + WWD_ITER_SPLIT_PTR_DENS];
  while (loc1 > VALID_PTR)
    {
      /* Get maximum density */

      rhmax = RDB[loc1 + WWD_ITER_SPLIT_DENS_RHO];

      /* Get density */

      if (rhmax < 0.0)
        rho = RDB[mat + MATERIAL_MDENS];
      else
        rho = RDB[mat + MATERIAL_ADENS];

      /* Check material density */

      if (rho >= fabs(rhmax))
        {
          /* Set split flag */

#ifdef OPEN_MP
#pragma omp critical
#endif
          {
            WDB[ptr] = (double)YES;
          }

          /* Break loop */

          break;
        }

      /* Update pointer */

      ptr++;

      /* Next */

      loc1 = NextItem(loc1);
    }

  /* Tästä eteenpäin vanhaa ifp:hen perustuvaa rutiinia */

#ifdef mmmmmmmmmmmmmmmmmmmmmmm

  long rmx, loc0, ptr;
  double rho, rhmax;

  /* Check if response matrix calculation is on */

  if ((long)RDB[DATA_RMTX_CALC] == NO)
    return;

  /* Check if mfp iteration */

  if ((long)RDB[DATA_RMTX_MFP_CALC] == NO)
    return;

  /* Check if source convergence acceleration is on */

  if ((long)RDB[DATA_RMX_CONVG_ACC] == YES)
    return;

  /* Check if active cycle */

  if (RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP])
    return;

  /* Check particle pointer */

  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Pointer to weight window structure */

  rmx = (long)RDB[DATA_PTR_RMX0];
  CheckPointer(FUNCTION_NAME, "(rmx)", DATA_ARRAY, rmx);

  /* Pointer to iteration */

  ptr = (long)RDB[rmx + RMX_PTR_ITER];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get limiting density */

  rhmax = RDB[ptr + WWD_ITER_SPLIT_MFP_DENS];

  /* Get material density */

  if (mat < VALID_PTR)
    rho = 0.0;
  else if (rhmax < 0.0)
    rho = g*RDB[mat + MATERIAL_MDENS];
  else
    rho = g*RDB[mat + MATERIAL_ADENS];

  /* Get pointer to cell */

  loc0 = (long)RDB[part + PARTICLE_ICM_PTR_ICM];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Check density score separate results for thick and thin materials */

  if (rho > fabs(rhmax))
    {
      /* Score flux */

      ptr = (long)RDB[loc0 + RMX_CELL_MC_MFP_FLX0];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      AddPrivateRes(ptr, wgt*flx, id);

      /* Score total reaction rate */

      ptr = (long)RDB[loc0 + RMX_CELL_MC_MFP_TOTRR0];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      AddPrivateRes(ptr, wgt*tot, id);
    }
  else
    {
      /* Score flux */

      ptr = (long)RDB[loc0 + RMX_CELL_MC_MFP_FLX1];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      AddPrivateRes(ptr, wgt*flx, id);

      /* Score total reaction rate */

      ptr = (long)RDB[loc0 + RMX_CELL_MC_MFP_TOTRR1];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      AddPrivateRes(ptr, wgt*tot, id);
    }

#endif
}

/*****************************************************************************/
