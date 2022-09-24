/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : usertransmucut.c                               */
/*                                                                           */
/* Created:       2018/03/05 (JLe)                                           */
/* Last modified: 2018/03/05 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Performs transmutation reaction cut-off based on user-       */
/*              defined list.                                                */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "UserTransmuCut:"

/*****************************************************************************/

long UserTransmuCut(long nuc, long rea)
{
  long loc0, mt, mt0, n, m, ZAI, tst;
  char *str;

  /* Pointer to list */

  if ((loc0 = (long)RDB[DATA_DEP_TRANSMU_REA_LIST]) < VALID_PTR)
    return NO;

  /* Get mt */

  mt = (long)RDB[rea + REACTION_MT];
  mt0 = (long)RDB[rea + REACTION_BRANCH_MT];

  /* Check decay reactions */

  if (((mt > 10000) && (mt < 20000)) || ((mt0 > 10000) && (mt0 < 20000)))
    return NO;

  /* Loop over list and compare name */

  while ((long)RDB[loc0] > -1)
    {
      /* Pointer to text */

      str = &ASCII[(long)RDB[loc0]];

      /* Remove dashes */

      n = 0;
      m = 0;

      while (str[n] != '\0')
        {
          if (str[n] != '-')
            str[m++] = str[n];

          n++;
        }
      
      str[m] = '\0';

      /* Get ZAI */

      if ((ZAI = IsotoZAI(str)) < 0)
        ZAI = atoi(str);

      /* Update pointer */

      loc0++;

      /* Compare ZAI */

      if ((ZAI == (long)RDB[nuc + NUCLIDE_ZAI]) || (!strcasecmp(str, "all")))
        {
          /* Check all */

          if (!strcmp(GetText(loc0), "all"))
            break;
          else
            tst = atol(GetText(loc0));

          /* Compare MT */

          if (tst == mt)
            break;
          else if (tst == mt0)
            break;
          else if ((tst == 103) && ((mt > 599) && (mt < 650)))
            break;
          else if ((tst == 104) && ((mt > 649) && (mt < 700)))
            break;
          else if ((tst == 105) && ((mt > 699) && (mt < 750)))
            break;
          else if ((tst == 106) && ((mt > 749) && (mt < 800)))
            break;
          else if ((tst == 107) && ((mt > 799) && (mt < 850)))
            break;
          else if ((tst == 16) && ((mt > 874) && (mt < 892)))
            break;
          else if ((tst == 18) && ((mt == 19) || (mt == 20) ||
                                   (mt == 21) || (mt == 38)))
            break;
        }

      /* Next entry */

      loc0++;
    }
  
  /* Check */

  if ((long)RDB[loc0] < 0)
    return YES;

  /* Accept reaction */

  return NO;
}

/*****************************************************************************/
