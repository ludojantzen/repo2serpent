/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : alpha.c                                        */
/*                                                                           */
/* Created:       2012/11/01 (JLe)                                           */
/* Last modified: 2014/08/14 (JLe)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/* Description: Handles weight adjustment in alpha-eigenvalue mode           */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "Alpha:"

/*****************************************************************************/

void Alpha(double E, double majorant, double *wgt)
{
  double alpha, xs;

  /* Get eigenvalue */

  if ((alpha = RDB[DATA_ALPHA_EIG]) == 0.0)
    return;
  Die(FUNCTION_NAME, "Tää ei pelaa ollenkaan ST:n kanssa, ja surface tallyt antaa myös ihan vääriä tuloksia kun tota painoa muutetaan miten sattuu");
  /* Get alpha cross section */
  
  if ((xs = AlphaXS(E)) != 0.0)
    {
      /* Calculate weight reduction */
      
      if (alpha > 0.0)
        *wgt = *wgt*(1.0 - xs/majorant);
      else
        *wgt = *wgt*(1.0 + xs/majorant);
    }
  else
    Die(FUNCTION_NAME, "Error");
}

/*****************************************************************************/
