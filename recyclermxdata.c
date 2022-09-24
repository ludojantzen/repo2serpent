/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : recyclermxdata.c                               */
/*                                                                           */
/* Created:       2018/10/10 (JLe)                                           */
/* Last modified: 2019/04/16 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Recycles previously allocated RMX data blocks.               */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "RecycleRMXData:"

/*****************************************************************************/

void RecycleRMXData(long rmx, long loc1, long ng, long nmax, long file)
{
  long loc0, ptr, det, det0, mode, splitmax;

  /* Get pointer to dump */

  if ((loc0 = (long)RDB[DATA_RMX_ADA_PTR_CELL_DUMP]) < VALID_PTR)
    Die(FUNCTION_NAME, "Pointer error");

  /* Maximum number of splits */

  splitmax = (long)RDB[DATA_MAX_RMX_SPLIT_FLAGS];
  CheckValue(FUNCTION_NAME, "splitmax", "", splitmax, 0, 100);

  /* Get mode (ei määritetty convergence acceleration -moodissa) */

  mode = (long)RDB[rmx + RMX_MODE];
  CheckValue(FUNCTION_NAME, "mode", "", mode, 0, 4);

  /* Sort list */

  SortList(loc0, RMX_CELL_MTX_SIZE, SORT_MODE_ASCEND);

  /* Find cell with matching size */

  loc0 = (long)RDB[DATA_RMX_ADA_PTR_CELL_DUMP];
  while (loc0 > VALID_PTR)
    {
      /* Compare size */

      if ((long)RDB[loc0 + RMX_CELL_MTX_SIZE] == nmax)
        break;

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Check pointer */

  if (loc0 < VALID_PTR)
    {
      /***********************************************************************/

      /***** Matching size not found, create new *****************************/

      /* Check mode */

      if (mode == RMX_MODE_WWG)
        {
          /* Allocate memory for Monte Carlo results (wwg mode) */

          ptr = AllocPrivateData(nmax*ng, RES2_ARRAY);
          WDB[loc1 + RMX_CELL_MC_WWG_CURR0] = (double)ptr;

          ptr = AllocPrivateData(nmax*ng, RES2_ARRAY);
          WDB[loc1 + RMX_CELL_MC_WWG_CURR1] = (double)ptr;

          ptr = AllocPrivateData(ng, RES2_ARRAY);
          WDB[loc1 + RMX_CELL_MC_WWG_SRC0] = (double)ptr;

          ptr = AllocPrivateData(ng, RES2_ARRAY);
          WDB[loc1 + RMX_CELL_MC_WWG_SRC1] = (double)ptr;

          /* Check if partial solution was read from file */

          if (file == NO)
            {
              ptr = ReallocMem(DATA_ARRAY, nmax*ng);
              WDB[loc1 + RMX_CELL_ADJ_SOL_IN_CURR] = (double)ptr;
            }

          ptr = AllocPrivateData(ng, RES2_ARRAY);
          WDB[loc1 + RMX_CELL_MC_SRCC_TOT] = (double)ptr;
        }
      else
        {
          /* Allocate memory for Monte Carlo results (rmx modes) */

          ptr = AllocPrivateData(nmax*ng, RES2_ARRAY);
          WDB[loc1 + RMX_CELL_MC_CURR_NET_IN] = (double)ptr;

          ptr = AllocPrivateData(nmax*ng, RES2_ARRAY);
          WDB[loc1 + RMX_CELL_MC_CURR_NET_OUT] = (double)ptr;

          ptr = AllocPrivateData(nmax*nmax*ng*ng, RES2_ARRAY);
          WDB[loc1 + RMX_CELL_MC_CURR_THROUGH] = (double)ptr;

          ptr = AllocPrivateData(nmax*ng*ng, RES2_ARRAY);
          WDB[loc1 + RMX_CELL_MC_CURR_SRCC_OUT] = (double)ptr;

          ptr = AllocPrivateData(ng, RES2_ARRAY);
          WDB[loc1 + RMX_CELL_MC_SRCC_TOT] = (double)ptr;

          ptr = AllocPrivateData(1, RES2_ARRAY);
          WDB[loc1 + RMX_CELL_MC_MFP_TOTRR0] = (double)ptr;

          ptr = AllocPrivateData(1, RES2_ARRAY);
          WDB[loc1 + RMX_CELL_MC_MFP_FLX0] = (double)ptr;

          ptr = AllocPrivateData(1, RES2_ARRAY);
          WDB[loc1 + RMX_CELL_MC_MFP_TOTRR1] = (double)ptr;

          ptr = AllocPrivateData(1, RES2_ARRAY);
          WDB[loc1 + RMX_CELL_MC_MFP_FLX1] = (double)ptr;

          /* Allocate memory for calculated coefficients */

          if (((long)RDB[DATA_RMX_TEST_MODE] != NO) ||
              ((long)RDB[DATA_RMX_CONVG_ACC] == YES))
            {
              ptr = ReallocMem(DATA_ARRAY, nmax*ng*ng);
              WDB[loc1 + RMX_CELL_COEF_FWD_SRCC] = (double)ptr;

              ptr = ReallocMem(DATA_ARRAY, nmax*nmax*ng*ng);
              WDB[loc1 + RMX_CELL_COEF_FWD_TRANS_MTX] = (double)ptr;
            }

          ptr = ReallocMem(DATA_ARRAY, nmax*ng*ng);
          WDB[loc1 + RMX_CELL_COEF_ADJ_SRCC] = (double)ptr;

          ptr = ReallocMem(DATA_ARRAY, nmax*nmax*ng*ng);
          WDB[loc1 + RMX_CELL_COEF_ADJ_TRANS_MTX] = (double)ptr;

          /* Vectors used for solution */

          ptr = ReallocMem(DATA_ARRAY, nmax*ng);
          WDB[loc1 + RMX_CELL_WRK_IN_CURR] = (double)ptr;

          ptr = ReallocMem(DATA_ARRAY, nmax*ng);
          WDB[loc1 + RMX_CELL_WRK_OUT_CURR] = (double)ptr;

          ptr = ReallocMem(DATA_ARRAY, nmax*ng);
          WDB[loc1 + RMX_CELL_FWD_SOL_IN_CURR] = (double)ptr;

          /* Check if partial solution was read from file */

          if (file == NO)
            {
              ptr = ReallocMem(DATA_ARRAY, nmax*ng);
              WDB[loc1 + RMX_CELL_ADJ_SOL_IN_CURR] = (double)ptr;
            }

          ptr = ReallocMem(DATA_ARRAY, nmax*ng);
          WDB[loc1 + RMX_CELL_ADJ_SOL_OUT_CURR] = (double)ptr;

          /* Loop over detectors */

          det = (long)RDB[loc1 + RMX_CELL_PTR_DET];
          while (det > VALID_PTR)
            {
              /* Allocate memory for Monte Carlo results */

              ptr = AllocPrivateData(nmax*ng, RES2_ARRAY);
              WDB[det + RMX_DET_MC_RES_IN_CURR] = (double)ptr;

              ptr = AllocPrivateData(ng, RES2_ARRAY);
              WDB[det + RMX_DET_MC_RES_SRCC] = (double)ptr;

              ptr = AllocPrivateData(ng, RES2_ARRAY);
              WDB[det + RMX_DET_MC_RES_TOT] = (double)ptr;

              ptr = AllocPrivateData(1, RES2_ARRAY);
              WDB[det + RMX_DET_MC_RES_SCORE_N] = (double)ptr;

              /* Allocate memory for calculated coefficients */

              ptr = ReallocMem(DATA_ARRAY, nmax*ng);
              WDB[det + RMX_DET_COEF_FWD_RES] = (double)ptr;

              ptr = ReallocMem(DATA_ARRAY, nmax*ng);
              WDB[det + RMX_DET_COEF_ADJ_RES] = (double)ptr;

              ptr = ReallocMem(DATA_ARRAY, ng);
              WDB[det + RMX_DET_COEF_DIR_SRCC_RES] = (double)ptr;

              /* Next detector */

              det = NextItem(det);
            }
        }

      /* Split flags */

      if (splitmax > 0)
        {
          ptr = ReallocMem(DATA_ARRAY, splitmax);
          WDB[loc1 + RMX_CELL_PTR_SPLIT] = (double)ptr;
        }

      /***********************************************************************/
    }
  else
    {
      /***********************************************************************/

      /***** Recycle memory blocks *******************************************/

     /* Check mode */

      if (mode == RMX_MODE_WWG)
        {
          /* Allocate memory for Monte Carlo results (wwg mode) */

          ptr = (long)RDB[loc0 + RMX_CELL_MC_WWG_CURR0];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmax*ng*sizeof(double));
          WDB[loc1 + RMX_CELL_MC_WWG_CURR0] = (double)ptr;

          ptr = (long)RDB[loc0 + RMX_CELL_MC_WWG_CURR1];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmax*ng*sizeof(double));
          WDB[loc1 + RMX_CELL_MC_WWG_CURR1] = (double)ptr;

          ptr = (long)RDB[loc0 + RMX_CELL_MC_WWG_SRC0];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, ng*sizeof(double));
          WDB[loc1 + RMX_CELL_MC_WWG_SRC0] = (double)ptr;

          ptr = (long)RDB[loc0 + RMX_CELL_MC_WWG_SRC1];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, ng*sizeof(double));
          WDB[loc1 + RMX_CELL_MC_WWG_SRC1] = (double)ptr;

          ptr = (long)RDB[loc0 + RMX_CELL_ADJ_SOL_IN_CURR];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          memset(&WDB[ptr], 0.0, nmax*ng*sizeof(double));
          WDB[loc1 + RMX_CELL_ADJ_SOL_IN_CURR] = (double)ptr;

          ptr = (long)RDB[loc0 + RMX_CELL_MC_SRCC_TOT];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, ng*sizeof(double));
          WDB[loc1 + RMX_CELL_MC_SRCC_TOT] = (double)ptr;
        }
      else
        {
          /* Allocate memory for Monte Carlo results (rmx modes) */

          ptr = (long)RDB[loc0 + RMX_CELL_MC_CURR_NET_IN];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmax*ng*sizeof(double));
          WDB[loc1 + RMX_CELL_MC_CURR_NET_IN] = (double)ptr;

          ptr = (long)RDB[loc0 + RMX_CELL_MC_CURR_NET_OUT];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmax*ng*sizeof(double));
          WDB[loc1 + RMX_CELL_MC_CURR_NET_OUT] = (double)ptr;

          ptr = (long)RDB[loc0 + RMX_CELL_MC_CURR_THROUGH];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmax*nmax*ng*ng*sizeof(double));
          WDB[loc1 + RMX_CELL_MC_CURR_THROUGH] = (double)ptr;

          ptr = (long)RDB[loc0 + RMX_CELL_MC_CURR_SRCC_OUT];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, nmax*ng*ng*sizeof(double));
          WDB[loc1 + RMX_CELL_MC_CURR_SRCC_OUT] = (double)ptr;

          ptr = (long)RDB[loc0 + RMX_CELL_MC_SRCC_TOT];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, ng*sizeof(double));
          WDB[loc1 + RMX_CELL_MC_SRCC_TOT] = (double)ptr;

          ptr = (long)RDB[loc0 + RMX_CELL_MC_MFP_TOTRR0];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, sizeof(double));
          WDB[loc1 + RMX_CELL_MC_MFP_TOTRR0] = (double)ptr;

          ptr = (long)RDB[loc0 + RMX_CELL_MC_MFP_FLX0];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, sizeof(double));
          WDB[loc1 + RMX_CELL_MC_MFP_FLX0] = (double)ptr;

          ptr = (long)RDB[loc0 + RMX_CELL_MC_MFP_TOTRR1];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, sizeof(double));
          WDB[loc1 + RMX_CELL_MC_MFP_TOTRR1] = (double)ptr;

          ptr = (long)RDB[loc0 + RMX_CELL_MC_MFP_FLX1];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          memset(&RES2[ptr], 0.0, sizeof(double));
          WDB[loc1 + RMX_CELL_MC_MFP_FLX1] = (double)ptr;

          /* Calculated coefficients */

          if (((long)RDB[DATA_RMX_TEST_MODE] != NO) ||
              ((long)RDB[DATA_RMX_CONVG_ACC] == YES))
            {
              ptr = (long)RDB[loc0 + RMX_CELL_COEF_FWD_SRCC];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              memset(&WDB[ptr], 0.0, nmax*ng*ng*sizeof(double));
              WDB[loc1 + RMX_CELL_COEF_FWD_SRCC] = (double)ptr;

              ptr = (long)RDB[loc0 + RMX_CELL_COEF_FWD_TRANS_MTX];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              memset(&WDB[ptr], 0.0, nmax*nmax*ng*ng*sizeof(double));
              WDB[loc1 + RMX_CELL_COEF_FWD_TRANS_MTX] = (double)ptr;
            }

          ptr = (long)RDB[loc0 + RMX_CELL_COEF_ADJ_SRCC];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          memset(&WDB[ptr], 0.0, nmax*ng*ng*sizeof(double));
          WDB[loc1 + RMX_CELL_COEF_ADJ_SRCC] = (double)ptr;

          ptr = (long)RDB[loc0 + RMX_CELL_COEF_ADJ_TRANS_MTX];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          memset(&WDB[ptr], 0.0, nmax*nmax*ng*ng*sizeof(double));
          WDB[loc1 + RMX_CELL_COEF_ADJ_TRANS_MTX] = (double)ptr;

          /* Vectors used for solution */

          ptr = (long)RDB[loc0 + RMX_CELL_WRK_IN_CURR];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          memset(&WDB[ptr], 0.0, nmax*ng*sizeof(double));
          WDB[loc1 + RMX_CELL_WRK_IN_CURR] = (double)ptr;

          ptr = (long)RDB[loc0 + RMX_CELL_WRK_OUT_CURR];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          memset(&WDB[ptr], 0.0, nmax*ng*sizeof(double));
          WDB[loc1 + RMX_CELL_WRK_OUT_CURR] = (double)ptr;

          ptr = (long)RDB[loc0 + RMX_CELL_FWD_SOL_IN_CURR];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          memset(&WDB[ptr], 0.0, nmax*ng*sizeof(double));
          WDB[loc1 + RMX_CELL_FWD_SOL_IN_CURR] = (double)ptr;

          ptr = (long)RDB[loc0 + RMX_CELL_ADJ_SOL_IN_CURR];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          memset(&WDB[ptr], 0.0, nmax*ng*sizeof(double));
          WDB[loc1 + RMX_CELL_ADJ_SOL_IN_CURR] = (double)ptr;

          ptr = (long)RDB[loc0 + RMX_CELL_ADJ_SOL_OUT_CURR];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          memset(&WDB[ptr], 0.0, nmax*ng*sizeof(double));
          WDB[loc1 + RMX_CELL_ADJ_SOL_OUT_CURR] = (double)ptr;

          /* Pointer to original detector */

          det0 = (long)RDB[loc0 + RMX_CELL_PTR_DET];
          CheckPointer(FUNCTION_NAME, "(det0)", DATA_ARRAY, det0);

          /* Loop over detectors */

          det = (long)RDB[loc1 + RMX_CELL_PTR_DET];
          while (det > VALID_PTR)
            {
              /* Monte Carlo results */

              ptr = (long)RDB[det0 + RMX_DET_MC_RES_IN_CURR];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, nmax*ng*sizeof(double));
              WDB[det + RMX_DET_MC_RES_IN_CURR] = (double)ptr;

              ptr = (long)RDB[det0 + RMX_DET_MC_RES_SRCC];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, ng*sizeof(double));
              WDB[det + RMX_DET_MC_RES_SRCC] = (double)ptr;

              ptr = (long)RDB[det0 + RMX_DET_MC_RES_TOT];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, ng*sizeof(double));
              WDB[det + RMX_DET_MC_RES_TOT] = (double)ptr;

              ptr = (long)RDB[det0 + RMX_DET_MC_RES_SCORE_N];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              memset(&RES2[ptr], 0.0, 1*sizeof(double));
              WDB[det + RMX_DET_MC_RES_SCORE_N] = (double)ptr;

              /* Calculated coefficients */

              ptr = (long)RDB[det0 + RMX_DET_COEF_FWD_RES];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              memset(&WDB[ptr], 0.0, nmax*ng*sizeof(double));
              WDB[det + RMX_DET_COEF_FWD_RES] = (double)ptr;

              ptr = (long)RDB[det0 + RMX_DET_COEF_ADJ_RES];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              memset(&WDB[ptr], 0.0, nmax*ng*sizeof(double));
              WDB[det + RMX_DET_COEF_ADJ_RES] = (double)ptr;

              ptr = (long)RDB[det0 + RMX_DET_COEF_DIR_SRCC_RES];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              memset(&WDB[ptr], 0.0, ng*sizeof(double));
              WDB[det + RMX_DET_COEF_DIR_SRCC_RES] = (double)ptr;

              /* Next detector */

              det = NextItem(det);
            }
        }

      /* Split flags */

      if (splitmax > 0)
        {
          ptr = (long)RDB[loc0 + RMX_CELL_PTR_SPLIT];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          memset(&WDB[ptr], 0.0, splitmax*sizeof(double));
          WDB[loc1 + RMX_CELL_PTR_SPLIT] = (double)ptr;
        }

      /* Remove item */

      RemoveItem(loc0);

      /***********************************************************************/
    }
}

/*****************************************************************************/
