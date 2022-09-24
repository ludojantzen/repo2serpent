/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : reactioncutoff.c                               */
/*                                                                           */
/* Created:       2010/12/17 (JLe)                                           */
/* Last modified: 2012/12/28 (JLe)                                           */
/* Version:       2.1.12                                                     */
/*                                                                           */
/* Description: Removes unused reactions from nuclide lists                  */
/*                                                                           */
/* Comments: - Used to simplify the structure of processxsdata.c             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReactionCutoff:"

/*****************************************************************************/

void ReactionCutoff()
{
  long extras, decay, nuc, ptr, rea, type, mt;

  /* Flags */

  extras = (long)RDB[DATA_OPTI_INCLUDE_SPECIALS];
  decay = YES;

  /* Tää ei toimi gammadatan kanssa */

  return;

  /* Loop over nuclides */

  nuc = (long)RDB[DATA_PTR_NUC0];

  while (nuc > VALID_PTR)
    {
      /***********************************************************************/

      /***** Decay nuclides **************************************************/

      /* Check type */
      
      if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DECAY)
        {
          /* Check if data is used */

          if (decay == NO)
            {
              /* Pointer to next next nuclide */
              
              ptr = nuc;
              nuc = NextItem(nuc);
              
              /* Remove nuclide from list */
              
              RemoveItem(ptr);
            }

          /* Next nuclide */
          
          continue;
        }

      /***********************************************************************/

      /***** Other types *****************************************************/

      /* Pointer to reaction data */

      rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
                  
      /* Loop over reactions */

      while (rea > VALID_PTR)
        {
          /* Energy cut-off (NOTE: Tossa pitää olla toi yhtäläisyysmerkki). */

          if (RDB[rea + REACTION_EMIN] >= RDB[DATA_NEUTRON_EMAX])
            {
              /* Pointer to next reaction */

              ptr = rea;                  
              rea = NextItem(rea);
              
              /* Remove reaction from list */
              
              RemoveItem(ptr);

              /* Next reaction */
          
              continue;
            }

          /* Get type */

          type = (long)RDB[rea + REACTION_TYPE];

          /* Check type */
          
          if (type == REACTION_TYPE_DECAY)
            {
              /***** Decay reaction ******************************************/
              
              /* Remember reaction pointer */
              
              ptr = rea;
              
              /* Pointer to next reaction */
                  
              rea = NextItem(rea);
              
              /* Check if data is used */

              if (decay == NO)
                {
                  /* Remove reaction from list */
              
                  RemoveItem(ptr);
                }

              /* Next reaction */
          
              continue;

              /***************************************************************/
            }
          else if (type == REACTION_TYPE_SPECIAL)
            {
              /***** Special non-transport reactions *************************/

              /* Remember reaction pointer */
              
              ptr = rea;
              
              /* Pointer to next reaction */
                  
              rea = NextItem(rea);
              
              /* Check if data is used */

              if ((extras == NO) && ((long)RDB[nuc + NUCLIDE_TYPE] !=
                                     NUCLIDE_TYPE_DOSIMETRY))
                {
                  /* Remove reaction from list */
              
                  RemoveItem(ptr);
                }

              /* Next reaction */
          
              continue;

              /***************************************************************/
            }

          /***** Transport reactions *****************************************/
          
          /* Get mt */

          mt = (long)RDB[rea + REACTION_MT];

          /* Check */

          if ((mt != 2) && (mt != 1002) && (mt != 1004) &&
              ((mt < 16) || (mt > 199)))
            Die(FUNCTION_NAME, "Invalid reaction mt %ld (%s)",
                mt, GetText(nuc + NUCLIDE_PTR_NAME));
  
          /* Next reaction */

          rea = NextItem(rea);

          /*******************************************************************/
        }
      
      /***********************************************************************/
      
      /* Next nuclide */

      nuc = NextItem(nuc);
    }
}

/*****************************************************************************/
