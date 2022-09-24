/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : reactionccount.c                               */
/*                                                                           */
/* Created:       2010/12/17 (JLe)                                           */
/* Last modified: 2012/12/28 (JLe)                                           */
/* Version:       2.1.12                                                     */
/*                                                                           */
/* Description: Counts number of nuclides and reaction modes                 */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReactionCount:"

/*****************************************************************************/

void ReactionCount()
{
  long nuc, rea, type;

  /* Reset nuclide counters */

  WDB[DATA_N_TRANSPORT_NUCLIDES] = 0.0;
  WDB[DATA_N_DOSIMETRY_NUCLIDES] = 0.0;
  WDB[DATA_N_DECAY_NUCLIDES] = 0.0;
  WDB[DATA_N_PHOTON_NUCLIDES] = 0.0;
  WDB[DATA_N_TOT_NUCLIDES] = 0.0;
     
  /* Reset reaction counters */
      
  WDB[DATA_N_TRANSPORT_REA] = 0.0;
  WDB[DATA_N_SPECIAL_REA] = 0.0;
  WDB[DATA_N_TRANSMUTATION_REA] = 0.0;
  WDB[DATA_N_TRANSPORT_BRANCH] = 0.0;
  WDB[DATA_N_DECAY_REA] = 0.0;
  WDB[DATA_N_DECAY_BRANCH] = 0.0;

  /* Loop over nuclides */

  nuc = (long)RDB[DATA_PTR_NUC0];

  while (nuc > VALID_PTR)
    {
      /* Reset counters */
      
      WDB[nuc + NUCLIDE_N_TRANSPORT_REA] = 0.0; 
      WDB[nuc + NUCLIDE_N_SPECIAL_REA] = 0.0;
      WDB[nuc + NUCLIDE_N_DECAY_REA] = 0.0;
      WDB[nuc + NUCLIDE_N_TRANSPORT_BRANCH] = 0.0; 
      WDB[nuc + NUCLIDE_N_DECAY_BRANCH] = 0.0;
      WDB[nuc + NUCLIDE_N_TRANSMUTATION_REA] = 0.0; 

      /* Check type */
      
      if (((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_TRANSPORT) ||
          ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_SAB))
        WDB[DATA_N_TRANSPORT_NUCLIDES] = 
          RDB[DATA_N_TRANSPORT_NUCLIDES] + 1.0;
      else if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DOSIMETRY)
        WDB[DATA_N_DOSIMETRY_NUCLIDES] = 
          RDB[DATA_N_DOSIMETRY_NUCLIDES] + 1.0;
      else if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DECAY)
        WDB[DATA_N_DECAY_NUCLIDES] = RDB[DATA_N_DECAY_NUCLIDES] + 1.0;
      else if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_PHOTON)
        WDB[DATA_N_PHOTON_NUCLIDES] = RDB[DATA_N_PHOTON_NUCLIDES] + 1.0;
      else
        Die(FUNCTION_NAME, "Invalid nuclide type");
                  
      /* Add to total */

      WDB[DATA_N_TOT_NUCLIDES] = RDB[DATA_N_TOT_NUCLIDES] + 1.0;

      /* Loop over reactions */

      rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
      while (rea > VALID_PTR)
        {
          /* Get type */

          type = (long)RDB[rea + REACTION_TYPE];

          /* Check type */
          
          if (type == REACTION_TYPE_DECAY)
            WDB[nuc + NUCLIDE_N_DECAY_REA] = 
              RDB[nuc + NUCLIDE_N_DECAY_REA] + 1.0;
          else if (type == REACTION_TYPE_SPECIAL)
            WDB[nuc + NUCLIDE_N_SPECIAL_REA] = 
              RDB[nuc + NUCLIDE_N_SPECIAL_REA] + 1.0;
          else if (type == REACTION_TYPE_DEC_BRANCH)
            WDB[nuc + NUCLIDE_N_DECAY_BRANCH] = 
              RDB[nuc + NUCLIDE_N_DECAY_BRANCH] + 1.0;
          else if (type == REACTION_TYPE_TRA_BRANCH)
            WDB[nuc + NUCLIDE_N_TRANSPORT_BRANCH] = 
              RDB[nuc + NUCLIDE_N_TRANSPORT_BRANCH] + 1.0;
          else
            {
              WDB[nuc + NUCLIDE_N_TRANSPORT_REA] = 
                RDB[nuc + NUCLIDE_N_TRANSPORT_REA] + 1.0;
              
              /* Check transmutation reaction */
              
              if ((long)RDB[rea + REACTION_TGT_ZAI] != 0.0)
                WDB[nuc + NUCLIDE_N_TRANSMUTATION_REA] = 
                  RDB[nuc + NUCLIDE_N_TRANSMUTATION_REA] + 1.0;
            }

          /* Next reaction */

          rea = NextItem(rea);
        }
      
      /* Add reaction */

      WDB[DATA_N_TRANSPORT_REA] = RDB[DATA_N_TRANSPORT_REA] +
        RDB[nuc + NUCLIDE_N_TRANSPORT_REA];
      WDB[DATA_N_TRANSMUTATION_REA] = RDB[DATA_N_TRANSMUTATION_REA] +
        RDB[nuc + NUCLIDE_N_TRANSMUTATION_REA];
      WDB[DATA_N_SPECIAL_REA] = RDB[DATA_N_SPECIAL_REA] +
        RDB[nuc + NUCLIDE_N_SPECIAL_REA];
      WDB[DATA_N_TRANSPORT_BRANCH] = RDB[DATA_N_TRANSPORT_BRANCH] +
        RDB[nuc + NUCLIDE_N_TRANSPORT_BRANCH];
      WDB[DATA_N_DECAY_REA] = RDB[DATA_N_DECAY_REA] +
        RDB[nuc + NUCLIDE_N_DECAY_REA];
      WDB[DATA_N_DECAY_BRANCH] = RDB[DATA_N_DECAY_BRANCH] +
        RDB[nuc + NUCLIDE_N_DECAY_BRANCH];
      
      /* Next nuclide */

      nuc = NextItem(nuc);
    }
}

/*****************************************************************************/
