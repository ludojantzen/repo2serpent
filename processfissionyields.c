/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processfissionyields.c                         */
/*                                                                           */
/* Created:       2010/09/12 (JLe)                                           */
/* Last modified: 2020/06/11 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Adds fission yield data to nuclides and perform yield        */
/*              cut-off.                                                     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessFissionYields:"

/*****************************************************************************/

void ProcessFissionYields(long nuc)
{
  long rea, mt, ZAI, N, yld, loc0, loc1, n, i, ptr, nfp;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

  /* Check that nuclide has fission channel or decay mode */

  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
  while (rea > 0)
    {
      /* Get mt */

      mt = (long)RDB[rea + REACTION_MT];

      /* Check */

      if ((mt == 18) || (mt == 19) || (mt == 20) || (mt == 21) || (mt == 38) ||
          (mt - 10000 == 6))
        break;

      /* Next */

      rea = NextItem(rea);
    }

  /* Check if found */

  if (rea < 0)
    return;

  /***************************************************************************/

  /**** Copy data ************************************************************/

  /* Check type */

  if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_TRANSPORT)
    if ((long)RDB[DATA_PTR_ACE_NFY_DATA] > 0)
      {
        /*********************************************************************/

        /***** NFY data ******************************************************/

        /* Match ZAI */

        ZAI = (long)RDB[nuc + NUCLIDE_ZAI];

        /* Loop over nuclides */

        yld = (long)RDB[DATA_PTR_ACE_NFY_DATA];
        while (yld > 0)
          {
            /* Compare */

            if ((long)ACE[yld + FISSION_YIELD_PARENT_ZAI] == ZAI)
              break;

            /* Next */

            yld = (long)ACE[yld + FISSION_YIELD_PTR_NEXT];
          }

        /* Check if found */

        if (yld < 0)
          {
            /* Loop over nuclides */

            yld = (long)RDB[DATA_PTR_ACE_NFY_DATA];
            while (yld > 0)
              {
                /* Match ZA (isomeric states) */

                if ((long)(ACE[yld + FISSION_YIELD_PARENT_ZAI]/10.0)
                    == (long)RDB[nuc + NUCLIDE_ZA])
                  {
                    ZAI = (long)ACE[yld + FISSION_YIELD_PARENT_ZAI];
                    break;
                  }
                /* Next */

                yld = (long)ACE[yld + FISSION_YIELD_PTR_NEXT];
              }
          }

        /* Check if found */

        if (yld < 0)
          {
            /* Loop over nuclides */

            yld = (long)RDB[DATA_PTR_ACE_NFY_DATA];
            while (yld > 0)
              {
                /* Match neutron number (same thing done in Serpent 1) */

                N = (long)(ACE[yld + FISSION_YIELD_PARENT_ZAI]/10000.0);
                N = (long)ACE[yld + FISSION_YIELD_PARENT_ZAI] - 10000*N;
                N = (long)((double)N/10.0);
                N = N - (long)(ACE[yld + FISSION_YIELD_PARENT_ZAI]/10000.0);

                /* Compare */

                if (N == (long)RDB[nuc + NUCLIDE_A] -
                    (long)RDB[nuc + NUCLIDE_Z])
                  {
                    ZAI = (long)ACE[yld + FISSION_YIELD_PARENT_ZAI];
                    break;
                  }

                /* Next */

                yld = (long)ACE[yld + FISSION_YIELD_PTR_NEXT];
              }
          }

        /* Check if found */

        if (yld < 0)
          {
            /* Loop over nuclides */

            yld = (long)RDB[DATA_PTR_ACE_NFY_DATA];
            while (yld > 0)
              {
                /* Match U-235, Pu-239 or U-233 data */

                if (((long)ACE[yld + FISSION_YIELD_PARENT_ZAI] == 922350) ||
                    ((long)ACE[yld + FISSION_YIELD_PARENT_ZAI] == 942390) ||
                    ((long)ACE[yld + FISSION_YIELD_PARENT_ZAI] == 922330))
                  {
                    ZAI = (long)ACE[yld + FISSION_YIELD_PARENT_ZAI];
                    break;
                  }

                /* Next */

                yld = (long)ACE[yld + FISSION_YIELD_PTR_NEXT];
              }
          }

        /* Check if found */

        if (yld < 0)
          {
            /* Match any */

            yld = (long)RDB[DATA_PTR_ACE_NFY_DATA];

            if(yld > 0)
              ZAI = (long)ACE[yld + FISSION_YIELD_PARENT_ZAI];
          }

        /* Check pointer */

        if (yld < 1)
          Die(FUNCTION_NAME, "NFY data not found");

        /* Check Pointer to distribution */

        if ((long)ACE[yld + FISSION_YIELD_PTR_DISTR] < 1)
          Die(FUNCTION_NAME, "No NFY distribution");

        /* Loop over distributions */

        while (yld > 0)
          {
            /* Break if next nuclide */

            if (ZAI != (long)ACE[yld + FISSION_YIELD_PARENT_ZAI])
              break;

            /* Allocate memory */

            loc0 = NewItem(nuc + NUCLIDE_PTR_NFY_DATA,
                           FISSION_YIELD_BLOCK_SIZE);

            /* Set data flag */

            SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_NFY_DATA);

            /* Add number of interpolation energies */

            WDB[nuc + NUCLIDE_NFY_NE] = RDB[nuc + NUCLIDE_NFY_NE] + 1.0;

            /* Copy data */

            for (i = LIST_DATA_SIZE; i < FISSION_YIELD_BLOCK_SIZE - 1; i++)
              WDB[loc0 + i] = ACE[yld + i];

            /* Preserve original nfp */

            WDB[loc0 + FISSION_YIELD_ORIG_NFP] = ACE[yld + FISSION_YIELD_NFP];

            /* Put null to distribution pointer (this is needed for */
            /* calling NewItem() */

            WDB[loc0 + FISSION_YIELD_PTR_DISTR] = -1;

            /* Reset count */

            nfp = 0;

            /* Loop over distribution */

            n = (long)ACE[yld + FISSION_YIELD_PTR_DISTR];
            while ((ptr = (long)ACE[n++]) > 0)
              {
                /* Check that independent yield >= 0.0 and compare */
                /* cumulative yield to cut-off */

                if ((ACE[ptr + FY_INDEPENDENT_FRAC] > 0.0) &&
                    (ACE[ptr + FY_CUMULATIVE_FRAC] >=
                     RDB[DATA_DEP_FP_YIELD_CUTOFF]))
                  {
                    /* Allocate memory */

                    loc1 = NewItem(loc0 + FISSION_YIELD_PTR_DISTR,
                                   FY_BLOCK_SIZE);

                    /* Copy data */

                    for (i = LIST_DATA_SIZE; i < FY_BLOCK_SIZE; i++)
                      WDB[loc1 + i] = ACE[ptr + i];

                    /* Add counter */

                    nfp++;
                  }
              }

            /* Set count */

            WDB[loc0 + FISSION_YIELD_NFP] = (double)nfp;

            /* Pointer to distribution */

            loc1 = (long)RDB[loc0 + FISSION_YIELD_PTR_DISTR];

            /* Close list or remove distribution */

            if (loc1 > VALID_PTR)
              CloseList(loc1);
            else
              RemoveItem(loc0);

            /* Next */

            yld = (long)ACE[yld + FISSION_YIELD_PTR_NEXT];
          }

        /*********************************************************************/
      }

  /* Check if spontaneous yields are given */

  if ((long)RDB[DATA_PTR_ACE_SFY_DATA] > 0)
    {
      /***********************************************************************/

      /***** SFY data ********************************************************/

      /* Match ZAI */

      ZAI = (long)RDB[nuc + NUCLIDE_ZAI];

      /* Loop over nuclides */

      yld = (long)RDB[DATA_PTR_ACE_SFY_DATA];
      while (yld > 0)
        {
          /* Compare */

          if ((long)ACE[yld + FISSION_YIELD_PARENT_ZAI] == ZAI)
            break;

          /* Next */

          yld = (long)ACE[yld + FISSION_YIELD_PTR_NEXT];
        }

      /* Check pointer */

      if (yld > 0)
        {
          /* Allocate memory */

          loc0 = NewItem(nuc + NUCLIDE_PTR_SFY_DATA,
                         FISSION_YIELD_BLOCK_SIZE);

          /* Set data flag */

          SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_SFY_DATA);

          /* Copy data */

          for (i = LIST_DATA_SIZE; i < FISSION_YIELD_BLOCK_SIZE - 1; i++)
            WDB[loc0 + i] = ACE[yld + i];

          /* Preserve original nfp */

          WDB[loc0 + FISSION_YIELD_ORIG_NFP] = ACE[yld + FISSION_YIELD_NFP];

          /* Put null to distribution pointer (this is needed for */
          /* calling NewItem() */

          WDB[loc0 + FISSION_YIELD_PTR_DISTR] = -1;

          /* Reset count */

          nfp = 0;

          /* Loop over distribution */

          n = (long)ACE[yld + FISSION_YIELD_PTR_DISTR];
          while ((ptr = (long)ACE[n++]) >= 0)
            {
              /* Check that independent yield > 0.0 and compare */
              /* cumulative yield to cut-off */

              if ((ACE[ptr + FY_INDEPENDENT_FRAC] > 0.0) &&
                  (ACE[ptr + FY_CUMULATIVE_FRAC] >=
                   RDB[DATA_DEP_FP_YIELD_CUTOFF]))
                {
                  /* Allocate memory */

                  loc1 = NewItem(loc0 + FISSION_YIELD_PTR_DISTR,
                                 FY_BLOCK_SIZE);

                  /* Copy data */

                  for (i = LIST_DATA_SIZE; i < FY_BLOCK_SIZE; i++)
                    WDB[loc1 + i] = ACE[ptr + i];

                  /* Add counter */

                  nfp++;
                }
            }

          /* Set count */

          WDB[loc0 + FISSION_YIELD_NFP] = (double)nfp;

          /* Pointer to distribution */

          loc1 = (long)RDB[loc0 + FISSION_YIELD_PTR_DISTR];

          /* Close list or remove distribution */

          if (loc1 > VALID_PTR)
            CloseList(loc1);
          else
            RemoveItem(loc0);
        }

      /***********************************************************************/
    }

  /***************************************************************************/

  /***** Copy reaction pointers **********************************************/

  /* Loop over reactions */

  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
  while (rea > 0)
    {
      /* Get mt */

      mt = (long)RDB[rea + REACTION_MT];

      /* Check */

      if ((mt == 18) || (mt == 19) || (mt == 20) || (mt == 21) || (mt == 38))
        {
          /* Neutron-induced fission */

          if ((yld = (long)RDB[nuc + NUCLIDE_PTR_NFY_DATA]) > 0)
            WDB[rea + REACTION_PTR_FISSY] = (double)yld;
        }
      else if (mt - 10000 == 6)
        {
          /* Spontaneous fission */

          if ((yld = (long)RDB[nuc + NUCLIDE_PTR_SFY_DATA]) > 0)
            WDB[rea + REACTION_PTR_FISSY] = (double)yld;
        }

      /* Next */

      rea = NextItem(rea);
    }

  /***************************************************************************/
}

/*****************************************************************************/
