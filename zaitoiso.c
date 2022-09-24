/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : zaitoiso.c                                     */
/*                                                                           */
/* Created:       2010/11/19 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Converts ZAI to string variable                              */
/*                                                                           */
/* Comments: - From Serpent 1.1.14                                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "element_data.h"

#define FUNCTION_NAME "ZAtoIso:"

/*****************************************************************************/

char *ZAItoIso(long ZAI, long type)
{                
  static char iso[MAX_STR];
  long Z, A, I;

  /* Check elemental */

  if ((ZAI > 0) && (ZAI < 112))
    {
      sprintf(iso, "%s", element_names[ZAI]);
      return iso;
    }

  /* Check unknown */

  else if (ZAI < 10000)
    {
      sprintf(iso, "unknown");
      return iso;
    }

  /* Get I */

  Z = (long)((double)ZAI/10.0);
  I = ZAI - 10*Z;
  
  /* Get Z and A */

  Z = (long)((double)ZAI/10000.0);
  A = (long)((double)ZAI/10.0) - 1000*Z;

  if ((Z < 0) || (A < 0) || (Z > NUMBER_OF_ELEMENTS))
    sprintf(iso, "<No such isotope>");
  else if (I == 0)
    {
      if (A == 0)
        {
          if (type == 1)
            sprintf(iso, "%s-nat", element_symbols[Z]);
          else if (type == 3)
            sprintf(iso, "%snat", element_symbols[Z]);
          else
            sprintf(iso, "natural %s", element_names[Z]);
        }
      else
        {
          if (type == 1)
            sprintf(iso, "%s-%ld", element_symbols[Z], A);
          else if (type == 3)
            sprintf(iso, "%s%ld", element_symbols[Z], A);
          else
            {
              /* Hydrogen isotopes are special */

              if ((Z == 1) && (A == 1))
                sprintf(iso, "hydrogen");
              else if ((Z == 1) && (A == 2))
                sprintf(iso, "deuterium");
              else if ((Z == 1) && (A == 3))
                sprintf(iso, "tritium");
              else
                sprintf(iso, "%s %ld", element_names[Z], A);
            }
        }
    }
  else
    {
      if (type == 1)
        sprintf(iso, "%s-%ldm", element_symbols[Z], A);
      else if (type == 3)
        sprintf(iso, "%s%ldm", element_symbols[Z], A);
      else
              sprintf(iso, "%s %ldm", element_names[Z], A);
    }

  /* Return name */

  return iso;
}

/*****************************************************************************/
