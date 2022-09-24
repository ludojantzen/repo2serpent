/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : removestorefiles.c                             */
/*                                                                           */
/* Created:       2016/06/08 (VVa)                                           */
/* Last modified: 2016/06/08 (VVa)                                           */
/* Version:       2.1.27                                                     */
/*                                                                           */
/* Description: Removes files used to store particles between different      */
/*              batches in DYN simulation mode                               */
/* Comments:   -                                                             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "RemoveStoreFiles:"

/*****************************************************************************/

void RemoveStoreFiles()
{
  char tmpstr[MAX_STR];

  /* Remove live neutron stores */

  sprintf(tmpstr, "%s.storeA%d", GetText(DATA_PTR_INPUT_FNAME), mpiid);
  remove(tmpstr);

  sprintf(tmpstr, "%s.storeB%d", GetText(DATA_PTR_INPUT_FNAME), mpiid);
  remove(tmpstr);

  /* Remove mesh-based precursor stores */

  sprintf(tmpstr, "%s.storepmeshA%d", GetText(DATA_PTR_INPUT_FNAME), mpiid);
  remove(tmpstr);

  sprintf(tmpstr, "%s.storepmeshB%d", GetText(DATA_PTR_INPUT_FNAME), mpiid);
  remove(tmpstr);

  /* Remov point-wise precursor stores */

  sprintf(tmpstr, "%s.storeprecA%d", GetText(DATA_PTR_INPUT_FNAME), mpiid);
  remove(tmpstr);

  sprintf(tmpstr, "%s.storeprecB%d", GetText(DATA_PTR_INPUT_FNAME), mpiid);
  remove(tmpstr);

  /***************************************************************************/
}

/*****************************************************************************/
