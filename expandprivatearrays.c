/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : expandprivatearrays.c                          */
/*                                                                           */
/* Created:       2014/02/25 (JLe)                                           */
/* Last modified: 2014/02/25 (JLe)                                           */
/* Version:       2.1.18                                                     */
/*                                                                           */
/* Description: Expands PRIVA, BUF and RES2 arrays for OpenMP                */
/*              parallelization.                                             */
/*                                                                           */
/* Comments: - Must be called before data is accessed                        */
/*                                                                           */
/*           - Jos hyvin käy niin tän voi korvata Mem() -funktion optiolla   */
/*             joka sallii noiden blokkien käytön (riippuu siitä hävittääkö  */
/*             allocprivatedata nyt jotain tärkeää kun sisältöä ei enää      */
/*             kopioida uuteen muistialueeseen).                             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ExpandPrivateArrays:"

/*****************************************************************************/

void ExpandPrivateArrays()
{
  /* Set flag */

  WDB[DATA_PRIVA_MEM_READY] = (double)YES;
}

/*****************************************************************************/
