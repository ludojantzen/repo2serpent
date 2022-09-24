/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : combinefissionyields.c                         */
/*                                                                           */
/* Created:       2010/09/12 (JLe)                                           */
/* Last modified: 2019/09/05 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Combines fission yields into a single list of nuclides       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CombineFissionYields:"

/*****************************************************************************/

void CombineFissionYields()
{
  long nuc, rea, yld, loc0, ZAI, Z, A, nfp, n, tot, ace;
  double *zai, *all;

  /* Check burnup calculation mode */

  if ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == NO)
    return;

  /* Allocate memory for vectors */

  zai = (double *)Mem(MEM_ALLOC, MAX_FP_NUCLIDES, sizeof(double));
  all = NULL;

  /* Reset count */

  tot = 0;

  /* Loop over nuclides */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR)
    {
      if ((((long)RDB[nuc + NUCLIDE_PTR_NFY_DATA] > VALID_PTR) ||
          ((long)RDB[nuc + NUCLIDE_PTR_SFY_DATA] > VALID_PTR)) &&
          ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_NEW_DAUGHTERS))
        {
          /* Loop over reactions */

          rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
          while (rea > VALID_PTR)
            {
              /* Get pointer to yield */

              if ((yld = (long)RDB[rea + REACTION_PTR_FISSY]) > VALID_PTR)
                {
                  /* Get pointer to distribution */

                  if ((yld = (long)RDB[yld + FISSION_YIELD_PTR_DISTR]) <
                      VALID_PTR)
                    Die(FUNCTION_NAME, "No pointer to distribution");

                  /* Loop over yield */

                  nfp = 0;
                  while ((loc0 = ListPtr(yld, nfp)) > VALID_PTR)
                    {
                      /* Get ZAI */

                      ZAI = (long)RDB[loc0 + FY_TGT_ZAI];

                      /* Separate Z and A for check */

                      Z = (long)((double)ZAI/10000.0);
                      A = (long)((double)ZAI/10.0 - (double)Z*1000.0);

                      /* Check values */

                      if ((Z < 1) || (Z > 80))
                        Die(FUNCTION_NAME, "Error in Z");

                      if ((A < 1) || (A > 210))
                        Die(FUNCTION_NAME, "Error in A");

                      /* Add to vector */

                      zai[nfp++] = (double)ZAI;
                    }

                  /* Sort array */

                  SortArray(zai, nfp);

                  /* Add nuclides to total vector */

                  all = AddPts(all, &tot, zai, nfp);
                }

              /* Next reaction */

              rea = NextItem(rea);
            }
        }

      /* Next nuclide */

      nuc = NextItem(nuc);
    }

  /* Reset count */

  nfp = 0;

  /* Check if fission product yields are read */

  if ((long)RDB[DATA_PTR_NFYDATA_FNAME_LIST] > VALID_PTR)
    {
      /* Add nuclides from ACE data */

      ace = (long)RDB[DATA_PTR_ACE0];
      while (ace > VALID_PTR)
        {
          /* Get ZAI */

          ZAI = (long)ACE[ace + ACE_ZAI];

          /* Separate Z and A for check */

          Z = (long)((double)ZAI/10000.0);
          A = (long)((double)ZAI/10.0 - (double)Z*1000.0);

          /* Check mass chain to see that nuclide is in one of the two peaks. */

          if ((((A > 79) && (A < 111)) || ((A > 129) && (A < 151))) &&
              ((long)ACE[ace + ACE_TYPE] == NUCLIDE_TYPE_TRANSPORT))
            {
              /* Check if nuclide was already added */

              for (n = 0; n < nfp; n++)
                if (ACE[ace + ACE_ZAI] == zai[n])
                  break;

              /* Check count */

              if (n == nfp)
                zai[nfp++] = ACE[ace + ACE_ZAI];

              /* Check count */

              if (nfp == MAX_FP_NUCLIDES)
                Die(FUNCTION_NAME, "Array overflow");
            }

          /* Next nuclide */

          ace = (long)ACE[ace + ACE_PTR_NEXT];
        }

      /* Sort array */

      SortArray(zai, nfp);

      /* Add nuclides to total vector */

      all = AddPts(all, &tot, zai, nfp);
    }

  /* Check count */

  if (tot > 0)
    {
      /* Check duplicates and order */

      for (n = 1; n < tot; n++)
        if (all[n] <= all[n - 1])
          Die(FUNCTION_NAME, "Error in list");

      /* Allocate memory for list */

      loc0 = ReallocMem(DATA_ARRAY, tot + 1);
      WDB[DATA_PTR_FP_ZAI_LIST] = (double)loc0;

      /* Read data */

      for (n = 0; n < tot; n++)
        WDB[loc0++] = all[n];

      /* Null terminator */

      WDB[loc0] = -1.0;

      /* Set size */

      WDB[DATA_TOT_FP_NUCLIDES] = (double)tot;

      /* Free total vector */

      Mem(MEM_FREE, all);
    }

  /* Free temporary vector */

  Mem(MEM_FREE, zai);
}

/*****************************************************************************/
