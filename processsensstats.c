/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processsensstats.c                             */
/*                                                                           */
/* Created:       2017/05/02 (VVa)                                           */
/* Last modified: 2019/01/18 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Allocates memory for sensitivity calculation statistics.     */
/*                                                                           */
/* Comments: -Needs to be executed after processsenszais.c which is executed */
/*            after processnuclides.c                                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessSensStats:"

/*****************************************************************************/

void ProcessSensStats()
{
  long sens, nmat, nzai, nrea, nene, nxyz, nmu, maxi, np, ptr, resp, type, nblock;
  long loc1, idx, i;
  char tmpstr[MAX_STR];

  /* Get pointer to sensitivity block or return -1*/

  if ((sens = (long)RDB[DATA_PTR_SENS0]) < VALID_PTR)
    return;

  /* Calculate number of scores for Sens */

  nmat = (long)RDB[sens + SENS_N_MAT];
  nzai = (long)RDB[sens + SENS_N_ZAI];
  nrea = (long)RDB[sens + SENS_N_PERT] + 1;
  nene = (long)RDB[sens + SENS_N_ENE] + 1;
  nmu  = (long)RDB[sens + SENS_N_MU];
  nxyz = (long)RDB[sens + SENS_N_SPAT];

  maxi = (1 + (nmat*nzai*nrea*nene*nmu*nxyz));

  /***********************************/
  /* Process basic result estimators */
  /***********************************/

  if ((np = (long)RDB[DATA_SENS_LAST_GEN]) > 0)
    {
      /* keff sensitivity */

      ptr = NewStat("ADJ_PERT_KEFF_SENS", 2, np, maxi);
      WDB[RES_ADJ_PERT_KEFF_SENS] = (double)ptr;

      /* total normalization buffer */

      ptr = NewStat("RES_SENS_BUFF", 2, np, maxi);
      WDB[RES_SENS_BUFF] = (double)ptr;

      /* delayed neutron fraction sensitivity */

      ptr = NewStat("ADJ_PERT_BEFF_SENS", 2, np, maxi);
      WDB[RES_ADJ_PERT_BEFF_SENS] = (double)ptr;

      /* Group-wise delayed neutron fraction sensitivity */

      for (i = 0; i < 8; i++)
        {
          switch (i)
            {
            case 0:
              loc1 = RES_ADJ_PERT_BEFF_G1_SENS;
              break;
            case 1:
              loc1 = RES_ADJ_PERT_BEFF_G2_SENS;
              break;
            case 2:
              loc1 = RES_ADJ_PERT_BEFF_G3_SENS;
              break;
            case 3:
              loc1 = RES_ADJ_PERT_BEFF_G4_SENS;
              break;
            case 4:
              loc1 = RES_ADJ_PERT_BEFF_G5_SENS;
              break;
            case 5:
              loc1 = RES_ADJ_PERT_BEFF_G6_SENS;
              break;
            case 6:
              loc1 = RES_ADJ_PERT_BEFF_G7_SENS;
              break;
            case 7:
              loc1 = RES_ADJ_PERT_BEFF_G8_SENS;
              break;
            default:
              Die(FUNCTION_NAME, "WTF!");
            }

          sprintf(tmpstr, "ADJ_PERT_BEFF_GROUP_%ld_SENS", i+1);
          ptr = NewStat(tmpstr, 2, np, maxi);
          WDB[loc1] = (double)ptr;
        }

      /* delayed neutron decay constant sensitivity */

      ptr = NewStat("ADJ_PERT_LAMBDA_SENS", 2, np, maxi);
      WDB[RES_ADJ_PERT_LAMBDA_SENS] = (double)ptr;

      /* Group-wise delayed neutron fraction sensitivity */

      for (i = 0; i < 8; i++)
        {
          switch (i)
            {
            case 0:
              loc1 = RES_ADJ_PERT_LAMBDA_G1_SENS;
              break;
            case 1:
              loc1 = RES_ADJ_PERT_LAMBDA_G2_SENS;
              break;
            case 2:
              loc1 = RES_ADJ_PERT_LAMBDA_G3_SENS;
              break;
            case 3:
              loc1 = RES_ADJ_PERT_LAMBDA_G4_SENS;
              break;
            case 4:
              loc1 = RES_ADJ_PERT_LAMBDA_G5_SENS;
              break;
            case 5:
              loc1 = RES_ADJ_PERT_LAMBDA_G6_SENS;
              break;
            case 6:
              loc1 = RES_ADJ_PERT_LAMBDA_G7_SENS;
              break;
            case 7:
              loc1 = RES_ADJ_PERT_LAMBDA_G8_SENS;
              break;
            default:
              Die(FUNCTION_NAME, "WTF!");
            }

          sprintf(tmpstr, "ADJ_PERT_LAMBDA_GROUP_%ld_SENS", i+1);
          ptr = NewStat(tmpstr, 2, np, maxi);
          WDB[loc1] = (double)ptr;
        }

      /* prompt generation time sensitivity */

      ptr = NewStat("ADJ_PERT_LEFF_SENS", 2, np, maxi);
      WDB[RES_ADJ_PERT_LEFF_SENS] = (double)ptr;

      /* void reactivity coefficient sensitivity */

      ptr = NewStat("ADJ_PERT_VOID_SENS", 2, np, maxi);
      WDB[RES_ADJ_PERT_VOID_SENS] = (double)ptr;

      /* Store last global stat */

      WDB[DATA_LAST_GLOBAL_STAT] = (double)ptr;
    }

  /* Check pointer to the first covariance block */

  if ((ptr = (long)RDB[DATA_PTR_COVBLOCK0]) > VALID_PTR)
    {
      /* Count the covariance blocks */

      nblock = 0;

      while (ptr > VALID_PTR)
        {
          /* Increment number of covariance blocks */

          nblock++;

          /* Next block */

          ptr = NextItem(ptr);
        }

      /* Store the number of covariance blocks */

      WDB[sens + SENS_N_COV_BLOCK] = (double)nblock;

      /* Allocate memory for uncertainty results if needed */

      resp = (long)RDB[sens + SENS_PTR_RESP0];

      while (resp > VALID_PTR)
        {
          /* Print statistic name */

          sprintf(tmpstr, "%s_COV_DATA_UNCERTAINTY", GetText(resp + SENS_RESP_PTR_NAME));

          /* Allocate memory for uncertainty result and store pointer */

          ptr = NewStat(tmpstr, 2, np, nblock+1);
          WDB[resp + SENS_RESP_PTR_UNC_STAT] = (double)ptr;

          /* Print statistic name */

          sprintf(tmpstr, "%s_COV_DATA_VARIANCE", GetText(resp + SENS_RESP_PTR_NAME));

          /* Allocate memory for variance result and store pointer */

          ptr = NewStat(tmpstr, 2, np, nblock+1);
          WDB[resp + SENS_RESP_PTR_VAR_STAT] = (double)ptr;

          /* Next response */

          resp = NextItem(resp);
        }
    }

  /* Link pointers for basic response blocks */

  resp = (long)RDB[sens + SENS_PTR_RESP0];

  /* Loop over all responses and link based on type */

  while (resp > VALID_PTR)
    {
      /* Get response type */

      type = (long)RDB[resp + SENS_RESP_TYPE];

      /* Link pointers based on type */

      if (type == SENS_RESP_TYPE_KEFF)
        {
          /* K-eff will not be divided by anything */

          WDB[resp + SENS_RESP_PTR_STAT] = RDB[RES_ADJ_PERT_KEFF_SENS];
        }
      else if (type == SENS_RESP_TYPE_BEFF)
        {
          /* Check the bin index (0 = total) */

          idx = (long)RDB[resp + SENS_RESP_DET_BIN_IDX];

          /* Link the correct statistic */

          switch (idx)
            {
            case 0:
              WDB[resp + SENS_RESP_PTR_STAT] = RDB[RES_ADJ_PERT_BEFF_SENS];
              break;
            case 1:
              WDB[resp + SENS_RESP_PTR_STAT] = RDB[RES_ADJ_PERT_BEFF_G1_SENS];
              break;
            case 2:
              WDB[resp + SENS_RESP_PTR_STAT] = RDB[RES_ADJ_PERT_BEFF_G2_SENS];
              break;
            case 3:
              WDB[resp + SENS_RESP_PTR_STAT] = RDB[RES_ADJ_PERT_BEFF_G3_SENS];
              break;
            case 4:
              WDB[resp + SENS_RESP_PTR_STAT] = RDB[RES_ADJ_PERT_BEFF_G4_SENS];
              break;
            case 5:
              WDB[resp + SENS_RESP_PTR_STAT] = RDB[RES_ADJ_PERT_BEFF_G5_SENS];
              break;
            case 6:
              WDB[resp + SENS_RESP_PTR_STAT] = RDB[RES_ADJ_PERT_BEFF_G6_SENS];
              break;
            case 7:
              WDB[resp + SENS_RESP_PTR_STAT] = RDB[RES_ADJ_PERT_BEFF_G7_SENS];
              break;
            case 8:
              WDB[resp + SENS_RESP_PTR_STAT] = RDB[RES_ADJ_PERT_BEFF_G8_SENS];
              break;
            default:
              Die(FUNCTION_NAME, "WTF!");
            }

          /* beta-eff will be divided by buff */

          WDB[resp + SENS_RESP_PTR_STAT_DIVIDER] = RDB[RES_SENS_BUFF];
        }
      else if (type == SENS_RESP_TYPE_LAMBDA)
        {
          /* Check the bin index (0 = total) */

          idx = (long)RDB[resp + SENS_RESP_DET_BIN_IDX];

          /* Link the correct statistic */

          switch (idx)
            {
            case 0:
              WDB[resp + SENS_RESP_PTR_STAT] = RDB[RES_ADJ_PERT_LAMBDA_SENS];
              break;
            case 1:
              WDB[resp + SENS_RESP_PTR_STAT] = RDB[RES_ADJ_PERT_LAMBDA_G1_SENS];
              break;
            case 2:
              WDB[resp + SENS_RESP_PTR_STAT] = RDB[RES_ADJ_PERT_LAMBDA_G2_SENS];
              break;
            case 3:
              WDB[resp + SENS_RESP_PTR_STAT] = RDB[RES_ADJ_PERT_LAMBDA_G3_SENS];
              break;
            case 4:
              WDB[resp + SENS_RESP_PTR_STAT] = RDB[RES_ADJ_PERT_LAMBDA_G4_SENS];
              break;
            case 5:
              WDB[resp + SENS_RESP_PTR_STAT] = RDB[RES_ADJ_PERT_LAMBDA_G5_SENS];
              break;
            case 6:
              WDB[resp + SENS_RESP_PTR_STAT] = RDB[RES_ADJ_PERT_LAMBDA_G6_SENS];
              break;
            case 7:
              WDB[resp + SENS_RESP_PTR_STAT] = RDB[RES_ADJ_PERT_LAMBDA_G7_SENS];
              break;
            case 8:
              WDB[resp + SENS_RESP_PTR_STAT] = RDB[RES_ADJ_PERT_LAMBDA_G8_SENS];
              break;
            default:
              Die(FUNCTION_NAME, "WTF!");
            }

          /* beta-eff will be divided by buff */

          WDB[resp + SENS_RESP_PTR_STAT_DIVIDER] = RDB[RES_SENS_BUFF];
        }
      else if (type == SENS_RESP_TYPE_LEFF)
        {
          /* ell-eff will be divided by buff */

          WDB[resp + SENS_RESP_PTR_STAT] = RDB[RES_ADJ_PERT_LEFF_SENS];
          WDB[resp + SENS_RESP_PTR_STAT_DIVIDER] = RDB[RES_SENS_BUFF];
        }
      else if (type == SENS_RESP_TYPE_VOID)
        {
          /* void sensitivity will be divided by buff */

          WDB[resp + SENS_RESP_PTR_STAT] = RDB[RES_ADJ_PERT_VOID_SENS];
          WDB[resp + SENS_RESP_PTR_STAT_DIVIDER] = RDB[RES_SENS_BUFF];
        }

      /* Next response */

      resp = NextItem(resp);
    }

  /*************************************************************/
  /* Process stats of detectors linked as sensitivity response */
  /*************************************************************/
}
