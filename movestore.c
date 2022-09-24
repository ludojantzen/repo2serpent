/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : movestore.c                                    */
/*                                                                           */
/* Created:       2015/01/19 (VVa)                                           */
/* Last modified: 2016/09/16 (VVa)                                           */
/* Version:       2.1.27                                                     */
/*                                                                           */
/* Description: Moves dynamic mode neutrons from EOI store to BOI store      */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MoveStore:"

/*****************************************************************************/

void MoveStore()
{
  long boiflag, eoiflag;

  /* Stores are only used in coupled transient calculations */

  if (RDB[DATA_SIMULATION_MODE] != SIMULATION_MODE_DYN)
    return;

  /*********************/
  /* File based stores */
  /*********************/

  /* Check that EOI store and BOI store are different */

  boiflag = (long)RDB[DATA_BOI_STORE_NAME];
  eoiflag = (long)RDB[DATA_EOI_STORE_NAME];

  if ((boiflag + eoiflag != 1) ||
      (boiflag*eoiflag != 0))
    Die(FUNCTION_NAME, "Something wrong with store names %ld and %ld\n",
        boiflag, eoiflag);

  /***********************************************/
  /* Move file based stores by changing the flag */
  /***********************************************/

 if (eoiflag == 0)
    {
      WDB[DATA_BOI_STORE_NAME] = (double)0;
      WDB[DATA_EOI_STORE_NAME] = (double)1;
    }
  else
    {
      WDB[DATA_BOI_STORE_NAME] = (double)1;
      WDB[DATA_EOI_STORE_NAME] = (double)0;
    }

}
