/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : getrequiredcovmatrices.c                       */
/*                                                                           */
/* Created:       2018/06/20 (VVa)                                           */
/* Last modified: 2018/06/21 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Removes additional covariance matrices from list, e.g.       */
/*              total cross section when partials are available. Returns an  */
/*              array of the unique ZAI-mts present in the calculation       */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "GetRequiredCovMatrices:"

/*****************************************************************************/

void GetRequiredCovMatrices(long *ZAImtArray)
{
  long cmtx, next, sens, remove, zai, mt, mtarr, maxmt, i;

  /* Get pointer to the first covariance matrix or return */

  if ((cmtx = (long)RDB[DATA_PTR_COVMTX0]) < VALID_PTR)
    return;

  /* Get pointer to the sensitivity data block or return */

  if ((sens = (long)RDB[DATA_PTR_SENS0]) < VALID_PTR)
    return;

  /************************************************************************/
  /* Currently we remove every MT not included in sensitivity calculation */
  /************************************************************************/

  /* Get maximum MT flagged for sensitivty calculations */

  maxmt = (long)RDB[sens + SENS_MAX_MT];

  /* Get pointer to MT array for sensitivity calculations */

  mtarr = (long)RDB[sens + SENS_PTR_MT_INDICES];

  /************************************************************/
  /* Loop over covariance matrices and remove not needed ones */
  /************************************************************/

  next = cmtx;
  while (next > VALID_PTR)
    {
      /* Get pointer to current covariance matrix */

      cmtx = next;

      /* Test ZAI and MT against the ones flagged for sensitivity calculation */

      remove = 0;
      for (i = 0; i < 2; i++)
        {
          /* Get ZAI and mt*/

          if (i == 0)
            {
              zai = (long)RDB[cmtx + COVMTX_ZAI1];
              mt = (long)RDB[cmtx + COVMTX_MT1];
            }
          else
            {
              zai = (long)RDB[cmtx + COVMTX_ZAI2];
              mt = (long)RDB[cmtx + COVMTX_MT2];
            }

          /* Check if ZAI is not flagged for sensitivity calculation */

          if (FindSensZAIIndex(zai) < 0)
            {
              remove = 1;
              break;
            }

          /* Check if MT is larger than the maximum MT flagged for */
          /* sensitivity calculation */

          if (mt > maxmt)
            {
              remove = 1;

              /* Can still be nubar 452-456 or chi 1018 */

              /* Check against nubar */

              if (((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_NUBAR)
                  && ((451 < mt) && (mt < 457)))
                remove = 0;

              /* Check against chi */

              if (((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_CHI)
                  && (mt == 1018))
                remove = 0;

              /* Don't test further if this should be removed */

              if (remove == 1)
                break;
            }

          /* MT is not larger than maximum but might not be flagged     */
          /* MT 1 sensitivity is actually calculated as sum from others */
          /* and those that were not included during post processing    */

          if ((RDB[mtarr + mt] == 0) && (mt != 1))
            {
              remove = 1;
              break;
            }
        }

      /* Get pointer to the next covariance matrix */

      next = NextItem(cmtx);

      /* Remove current from the covariance matrix list if needed */

      if (remove == 1)
        {
          /* Remove pointer from list (pointer is lost now) */

          RemoveItem(cmtx);

          /* Print out something */

          fprintf(outp, "  Removed covariance matrix for (%ld mt %ld) -- (%ld mt %ld) "
                  "due to some ZAI or MT not being flagged for sensitivity "
                  "calculation.\n", (long)RDB[cmtx + COVMTX_ZAI1], (long)RDB[cmtx + COVMTX_MT1],
                  (long)RDB[cmtx + COVMTX_ZAI2], (long)RDB[cmtx + COVMTX_MT2]);

        }
    }

  /***************************************************************/
  /* Loop over covariance matrices and remove sum reaction modes */
  /***************************************************************/

  /* Get pointer to the first covariance matrix or return */

  if ((cmtx = (long)RDB[DATA_PTR_COVMTX0]) < VALID_PTR)
    return;

  next = cmtx;
  while (next > VALID_PTR)
    {
      /* Get pointer to current covariance matrix */

      cmtx = next;

      /* Test ZAI and MT against the ones flagged for sensitivity calculation */

      remove = 0;
      for (i = 0; i < 2; i++)
        {
          /* Get ZAI and mt*/

          if (i == 0)
            {
              zai = (long)RDB[cmtx + COVMTX_ZAI1];
              mt = (long)RDB[cmtx + COVMTX_MT1];
            }
          else
            {
              zai = (long)RDB[cmtx + COVMTX_ZAI2];
              mt = (long)RDB[cmtx + COVMTX_MT2];
            }

          /* Check if MT is 1 (total) */

          if (mt == 1)
            remove = 1;

          /* Check if MT is 452 (total nubar) */

          if (mt == 452)
            remove = 1;
        }

      /* Get pointer to the next covariance matrix */

      next = NextItem(cmtx);

      /* Remove current from the covariance matrix list if needed */

      if (remove == 1)
        {
          /* Remove pointer from list (pointer is lost now) */

          RemoveItem(cmtx);

          /* Print out something */

          fprintf(outp, "  Removed covariance matrix for (%ld mt %ld) -- (%ld mt %ld) "
                  "due to MT being for total XS or total nubar\n", (long)RDB[cmtx + COVMTX_ZAI1], (long)RDB[cmtx + COVMTX_MT1],
                  (long)RDB[cmtx + COVMTX_ZAI2], (long)RDB[cmtx + COVMTX_MT2]);

        }
    }
}
