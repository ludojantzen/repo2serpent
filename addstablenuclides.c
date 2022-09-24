/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : addstablenuclides.c                            */
/*                                                                           */
/* Created:       2010/09/26 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Adds nuclides defined only by mass                           */
/*                                                                           */
/* Comments: - Muuta aliohjelman nimi                                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"
#include "nuclide_masses.h"

#define FUNCTION_NAME "AddStableNuclides:"

/*****************************************************************************/

void AddStableNuclides()
{
  long ace, ptr, n;
  char tmpstr[MAX_STR];

  /***************************************************************************/

  /***** Add nuclides to ACE data ********************************************/

  /* Loop over array */

  n = 0;
  while (nuc[n][0] > 0.0)
    {
      /* Allocate memory for data */

      ace = ReallocMem(ACE_ARRAY, ACE_BLOCK_SIZE);
      
      /* Set null pointer to next */
      
      ACE[ace + ACE_PTR_NEXT] = NULLPTR;
      
      /* Check if previous exists (ei voi käyttää VALID_PTR) */
      
      if ((ptr = (long)RDB[DATA_PTR_ACE0]) < 1)
        {
          /* First definition, set pointer */
          
          WDB[DATA_PTR_ACE0] = (double)ace;
        }
      else
        {
          /* Find last block  (tohon ei VALID_PTR) */
          
          while ((long)ACE[ptr + ACE_PTR_NEXT] > 0)
            ptr = (long)ACE[ptr + ACE_PTR_NEXT];
            
          /* Set pointer to new */
          
          ACE[ptr + ACE_PTR_NEXT] = (double)ace;
        }

      /* Set type */
      
      ACE[ace + ACE_TYPE] = (double)NUCLIDE_TYPE_DECAY;

      /* Put values */

      ACE[ace + ACE_ZAI] = nuc[n][0]*10.0;
      ACE[ace + ACE_ZA] = nuc[n][0];
      ACE[ace + ACE_I] = 0.0;
      ACE[ace + ACE_AWR] = nuc[n][1]/M_NEUTRON;

      /* Put name and alias */

      sprintf(tmpstr, "%ld", (long)ACE[ace + ACE_ZAI]);
      ACE[ace + ACE_PTR_NAME] = (double)PutText(tmpstr);
      ACE[ace + ACE_PTR_ALIAS] = ACE[ace + ACE_PTR_NAME];

      /* Reset pointers */

      ACE[ace + ACE_PTR_NXS] = -1.0;
      ACE[ace + ACE_PTR_JXS] = -1.0;
      ACE[ace + ACE_PTR_XSS] = -1.0;
      ACE[ace + ACE_PTR_FILE] = -1.0;
      ACE[ace + ACE_LAMBDA] = -1.0;
      ACE[ace + ACE_DECAY_E] = -1.0;
      ACE[ace + ACE_PTR_DECAY_LIST] = -1.0;

      /* Temperature */

      ACE[ace + ACE_TEMP] = -1.0;

      /* Library id */

      ACE[ace + ACE_PTR_LIB_ID] = (double)PutText("str");
      
      /* Next */

      n++;
    }

  /***************************************************************************/

  /***** Add elements to ACE data ********************************************/

  /* Stabiileja alkuaineita tarvitaan vain silloin kun käytetään       */
  /* hajoamislähdettä ja materiaalikoostumus luetaan aikaisemmasta     */
  /* neutronilaskusta jonka alkukoostumuksessa on luonnonmateriaaleja. */
  /* Todella spesiaalikeissi, mutta tuottaa muuten kryptisen errorin.  */

  if (((long)RDB[DATA_USE_DECAY_SRC] == NO) ||
      ((long)RDB[DATA_READ_RESTART_FILE] == NO))
    return;

  /* Loop over array */

  n = 0;
  while (atom[n][0] > 0.0)
    {
      /* Allocate memory for data */

      ace = ReallocMem(ACE_ARRAY, ACE_BLOCK_SIZE);
      
      /* Set null pointer to next */
      
      ACE[ace + ACE_PTR_NEXT] = NULLPTR;
      
      /* Check if previous exists (ei voi käyttää VALID_PTR) */
      
      if ((ptr = (long)RDB[DATA_PTR_ACE0]) < 1)
        {
          /* First definition, set pointer */
          
          WDB[DATA_PTR_ACE0] = (double)ace;
        }
      else
        {
          /* Find last block  (tohon ei VALID_PTR) */
          
          while ((long)ACE[ptr + ACE_PTR_NEXT] > 0)
            ptr = (long)ACE[ptr + ACE_PTR_NEXT];
            
          /* Set pointer to new */
          
          ACE[ptr + ACE_PTR_NEXT] = (double)ace;
        }

      /* Set type */
      
      ACE[ace + ACE_TYPE] = (double)NUCLIDE_TYPE_DECAY;

      /* Put values */

      ACE[ace + ACE_ZAI] = atom[n][0]*10000.0;
      ACE[ace + ACE_ZA] = atom[n][0]*1000.0;
      ACE[ace + ACE_I] = 0.0;
      ACE[ace + ACE_AWR] = atom[n][1]/M_NEUTRON;

      /* Put name and alias */

      sprintf(tmpstr, "%ld", (long)ACE[ace + ACE_ZAI]);
      ACE[ace + ACE_PTR_NAME] = (double)PutText(tmpstr);
      ACE[ace + ACE_PTR_ALIAS] = ACE[ace + ACE_PTR_NAME];

      /* Reset pointers */

      ACE[ace + ACE_PTR_NXS] = -1.0;
      ACE[ace + ACE_PTR_JXS] = -1.0;
      ACE[ace + ACE_PTR_XSS] = -1.0;
      ACE[ace + ACE_PTR_FILE] = -1.0;
      ACE[ace + ACE_LAMBDA] = -1.0;
      ACE[ace + ACE_DECAY_E] = -1.0;
      ACE[ace + ACE_PTR_DECAY_LIST] = -1.0;

      /* Temperature */

      ACE[ace + ACE_TEMP] = -1.0;

      /* Library id */

      ACE[ace + ACE_PTR_LIB_ID] = (double)PutText("str");
      
      /* Next */

      n++;
    }

  /***************************************************************************/
}

/*****************************************************************************/
