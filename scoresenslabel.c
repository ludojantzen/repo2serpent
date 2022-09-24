/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scoresenslabel.c                               */
/*                                                                           */
/* Created:       2017/04/05 (VVa)                                           */
/* Last modified: 2018/06/20 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Scores sensitivity event based on label                      */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreSensLabel:"

/*****************************************************************************/

void ScoreSensLabel(long res, long gen, long idx,
                    double val, double wgt, long id)
{
  /* Add to correct bin */

  AddBuf(val, wgt, res, id, -1, gen, idx);
}
