/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processtimebins.c                              */
/*                                                                           */
/* Created:       2012/09/22 (JLe)                                           */
/* Last modified: 2019/03/26 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Processes time binnings                                      */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessTimeBins:"

/*****************************************************************************/

void ProcessTimeBins()
{
  long tme, nt, n, ptr;
  double *t, tmin, tmax;

  /***************************************************************************/

  /***** Process time binnings ***********************************************/

  /* Loop over binnings */

  tme = (long)RDB[DATA_PTR_TME0];
  while (tme > VALID_PTR)
    {
      /* Avoid compiler warning */

      t = NULL;

      /* Check type */

      switch((long)RDB[tme + TME_TYPE])
        {
        default:
          {
            /* Invalid type */

            Error(tme, "Invalid time binning type %ld",
                  (long)RDB[tme + TME_TYPE]);
          }
        case TB_TYPE_ARB:
          {
            /* Number of times */

            if ((nt = (long)RDB[tme + TME_NB] + 1) < 2)
              Error(tme, "Not enough values for time binning");

            /* Pointer to data */

            ptr = (long)RDB[tme + TME_PTR_BINS];
            CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

            /* Check that values are in ascending order */

            for (n = 0; n < nt - 1; n++)
              if (RDB[ptr + n] >= RDB[ptr + n + 1])
                Error(tme, "Time bin intervals be in ascending order");

            /* Put minimum and maximum times */

            WDB[tme + TME_TMIN] = RDB[ptr];
            WDB[tme + TME_TMAX] = RDB[ptr + nt - 1];

            /* Break case */

            break;
          }
        case TB_TYPE_UNI_T:
        case TB_TYPE_UNI_LOGT:
          {
            /* Number of times */

            nt = (long)RDB[tme + TME_NB] + 1;

            /* Get minimum and maximum tmergy */

            tmin = RDB[tme + TME_TMIN];
            tmax = RDB[tme + TME_TMAX];

            /* Check type and make array */

            if ((long)RDB[tme + TME_TYPE] == TB_TYPE_UNI_T)
              t = MakeArray(tmin, tmax, nt, 1);
            else if ((long)RDB[tme + TME_TYPE] == TB_TYPE_UNI_LOGT)
              {
                /* Check for zero (tää on aika turha) */

                if (tmin == 0.0)
                  {
                    t = MakeArray(ZERO, tmax, nt, 2);
                    t[0] = 0.0;
                  }
                else
                  t = MakeArray(tmin, tmax, nt, 2);
              }
            else
              Die(FUNCTION_NAME, "Invalid array type");

            /* Allocate memory for values */

            ptr = ReallocMem(DATA_ARRAY, nt);
            WDB[tme + TME_PTR_BINS] = (double)ptr;

            /* Loop over values */

            for (n = 0; n < nt; n++)
              WDB[ptr++] = t[n];

            /* Free temporary array */

            Mem(MEM_FREE, t);

            /* Break case */

            break;
          }
        }

      /* Next binning */

      tme = NextItem(tme);
    }

  /***************************************************************************/

  /***** RIA simulation ******************************************************/

  if ((ptr = (long)RDB[DATA_PTR_RIA0]) > VALID_PTR)
    {
      /* Check pointer to time binning */

      if ((long)RDB[ptr + RIA_PTR_TME] < VALID_PTR)
        Error(ptr, "No time binning defined");

      /* Find bins */

      tme = (long)RDB[DATA_PTR_TME0];
      while (tme > VALID_PTR)
        {
          /* Compare names */

          if (CompareStr(ptr + RIA_PTR_TME, tme + TME_PTR_NAME))
            {
              /* Put pointer */

              WDB[ptr + RIA_PTR_TME] = (double)tme;

              /* Break loop */

              break;
            }

          /* Next binning */

          tme = NextItem(tme);
        }

      /* Check if found */

      if (tme < VALID_PTR)
        Error(ptr, "Time binning %s is not defined",
              GetText(ptr + RIA_PTR_TME));

      /* Get pointer to bins */

      ptr = (long)RDB[tme + TME_PTR_BINS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get number of bins */

      nt = (long)RDB[tme + TME_NB];
      CheckValue(FUNCTION_NAME, "nt", "", nt, 1, 1000000000);

      /* Put values */

      WDB[DATA_DYN_PTR_TIME_BINS] = (double)tme;
      WDB[DATA_DYN_TMIN] = RDB[ptr];
      WDB[DATA_DYN_TMAX] = RDB[ptr + nt];
      WDB[DATA_DYN_NB] = (double)nt;
    }

  /***************************************************************************/

  /***** Time-dependent simulation *******************************************/

  /* Check simulation mode (jos RIA-moodi, alustus on tehty */
  /* jo tossa ylempänä) */

  if (((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_SRC) &&
      ((long)RDB[DATA_PTR_RIA0]) < VALID_PTR)
    {
      /* Check if interval is defined */

      if (RDB[DATA_DYN_DT] > 0.0)
        {
          /* NOTE: Tässä viritetään dynaamisen simulaation aikavälit    */
          /* sellaiseksi että binejä enemmän kuin yksi, jolloin ajo     */
          /* menee transportcycle.c:ssä oikeaan silmukkaan. Aikavälejä  */
          /* ei kuitenkaan käytetä, vaan pysäytys tehdään DATA_DYN_DT:n */
          /* moninkertoina. */

          /* Use 10000 bins */

          nt = 10000;

          /* Create new item */

          tme = NewItem(DATA_PTR_TME0, TME_BLOCK_SIZE);

          /* Put type and number of bins */

          WDB[tme + TME_TYPE] = (double)TB_TYPE_ARB;
          WDB[tme + TME_NB] = (double)nt;

          /* Allocate memory for bins */

          ptr = ReallocMem(DATA_ARRAY, nt + 1);
          WDB[tme + TME_PTR_BINS] = (double)ptr;

          /* Put values */

          WDB[ptr++] = RDB[DATA_TIME_CUT_TMIN];
          WDB[ptr] = RDB[DATA_TIME_CUT_TMAX];

          /* Put pointer */

          WDB[DATA_DYN_PTR_TIME_BINS] = (double)tme;

          /* Put values */

          WDB[DATA_DYN_PTR_TIME_BINS] = (double)tme;
          WDB[DATA_DYN_TMIN] = RDB[DATA_TIME_CUT_TMIN];
          WDB[DATA_DYN_TMAX] = RDB[DATA_TIME_CUT_TMAX];
          WDB[DATA_DYN_NB] = 1.0;
          WDB[DATA_DYN_NB] = (double)nt;
        }

      /* Check if time binning is defined */

      else if ((long)RDB[DATA_DYN_PTR_TIME_BINS] < VALID_PTR)
        {
          /* Create new item */

          tme = NewItem(DATA_PTR_TME0, TME_BLOCK_SIZE);

          /* Put type and number of bins */

          WDB[tme + TME_TYPE] = (double)TB_TYPE_ARB;
          WDB[tme + TME_NB] = 1.0;

          /* Allocate memory for bins */

          ptr = ReallocMem(DATA_ARRAY, 2);
          WDB[tme + TME_PTR_BINS] = (double)ptr;

          /* Put values */

          WDB[ptr++] = RDB[DATA_TIME_CUT_TMIN];
          WDB[ptr] = RDB[DATA_TIME_CUT_TMAX];

          /* Put pointer */

          WDB[DATA_DYN_PTR_TIME_BINS] = (double)tme;

          /* Put values */

          WDB[DATA_DYN_TMIN] = RDB[DATA_TIME_CUT_TMIN];
          WDB[DATA_DYN_TMAX] = RDB[DATA_TIME_CUT_TMAX];
          WDB[DATA_DYN_NB] = 1.0;
         }
      else
        {
          /* Find bins */

          tme = (long)RDB[DATA_PTR_TME0];
          while (tme > VALID_PTR)
            {
              /* Compare names */

              if (CompareStr((long)DATA_DYN_PTR_TIME_BINS, tme + TME_PTR_NAME))
                {
                  /* Put pointer */

                  WDB[DATA_DYN_PTR_TIME_BINS] = (double)tme;

                  /* Break loop */

                  break;
                }

              /* Next binning */

              tme = NextItem(tme);
            }

          /* Check if found */

          if (tme < VALID_PTR)
            Error(0, "Time binning %s entered in the nps card is not defined",
                  GetText(DATA_DYN_PTR_TIME_BINS));

          /* Get pointer to bins */

          ptr = (long)RDB[tme + TME_PTR_BINS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Get number of bins */

          nt = (long)RDB[tme + TME_NB];
          CheckValue(FUNCTION_NAME, "nt", "", nt, 1, 1000000000);

          /* Put values */

          WDB[DATA_DYN_TMIN] = RDB[ptr];
          WDB[DATA_DYN_TMAX] = RDB[ptr + nt];
          WDB[DATA_DYN_NB] = (double)nt;
        }
    }

  /***************************************************************************/
}

/*****************************************************************************/
