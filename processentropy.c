/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processentropy.c                               */
/*                                                                           */
/* Created:       2011/05/18 (JLe)                                           */
/* Last modified: 2012/12/28 (JLe)                                           */
/* Version:       2.1.12                                                     */
/*                                                                           */
/* Description: Processes meshes used for fission source entropy calculation */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessEntropy:"

/*****************************************************************************/

void ProcessEntropy()
{
  long ptr;

  /* Check if entropies are calculated */

  if ((long)RDB[DATA_OPTI_ENTROPY_CALC] == NO)
    return;

  /* Set boundaries if not already set */

  if (RDB[DATA_ENTROPY_XMIN] == -INFTY)
    WDB[DATA_ENTROPY_XMIN] = RDB[DATA_GEOM_MINX];
  
  if (RDB[DATA_ENTROPY_XMAX] == INFTY)
    WDB[DATA_ENTROPY_XMAX] = RDB[DATA_GEOM_MAXX];
  
  if (RDB[DATA_ENTROPY_YMIN] == -INFTY)
    WDB[DATA_ENTROPY_YMIN] = RDB[DATA_GEOM_MINY];
  
  if (RDB[DATA_ENTROPY_YMAX] == INFTY)
    WDB[DATA_ENTROPY_YMAX] = RDB[DATA_GEOM_MAXY];
  
  if (RDB[DATA_ENTROPY_ZMIN] == -INFTY)
    WDB[DATA_ENTROPY_ZMIN] = RDB[DATA_GEOM_MINZ];
  
  if (RDB[DATA_ENTROPY_ZMAX] == INFTY)
    WDB[DATA_ENTROPY_ZMAX] = RDB[DATA_GEOM_MAXZ];

  /* Check mode */

  if ((long)RDB[DATA_SIMULATION_MODE] != SIMULATION_MODE_CRIT)
    return;

  /* Allocate memory for statistics */
  
  ptr = NewStat("ENTR_SPT", 1, 4);
  AllocStatHistory(ptr);
  WDB[DATA_ENTROPY_PTR_SPT_STAT] = (double)ptr;

  ptr = NewStat("ENTR_SWG", 1, 4);
  AllocStatHistory(ptr);
  WDB[DATA_ENTROPY_PTR_SWG_STAT] = (double)ptr;
}

/*****************************************************************************/
