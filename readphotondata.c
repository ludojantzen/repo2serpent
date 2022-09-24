/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readphotondata.c                               */
/*                                                                           */
/* Created:       2011/04/07 (JLe)                                           */
/* Last modified: 2017/07/06 (JLe)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Adds photon nuclides to materials defined with neutron cross */
/*              section composition                                          */
/*                                                                           */
/* Comments: - Nimi on huono, eikä kuvaa ollenkaan sitä mitä tää aliohjelma  */
/*             tekee.                                                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadPhotonData:"

/*****************************************************************************/

void ReadPhotonData()
{
  long elem[110], nuc, ptr, n;

  /* Check photon mode */

  if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == NO)
    return;

  /* Reset element array */

  for (n = 0; n < 110; n++)
    elem[n] = 0;

  /* Reset counter */

  n = 0;

  /* Loop over nuclides */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR)
    {
      /* Set element index */

      if ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_PHOTON)
        {
          /* Mark nuclide */

          elem[(long)RDB[nuc + NUCLIDE_Z]]++;

          /* Add counter */
          
          n++;
        }

      /* Next nuclide */

      nuc = NextItem(nuc);
    }

  /* Check count */

  if (n == 0)
    {
      /* Check if decay source is used */

      if ((long)RDB[DATA_USE_DECAY_SRC] == YES)
        Error(0, "Radioactive decay source requires isotopic material data");
      
      /* Exit subroutine */
      
      return;
    } 
  else
    fprintf(outp, "Adding photon interaction data...\n\n");

  /* Loop over element array */

  for (n = 0; n < 110; n++)
    if (elem[n] > 0)
      {
        /* Only materials Z < 100 are processed for now */

        if (n > 99)
          Note(0, "Photon physics models not available for %s (Z = %ld)",
               ZAItoIso(n, 1), n);
        
        /* Find nuclide in data */

        else if ((nuc = AddNuclide(NULL, n*10000, NULL, 0.0, 
                                   NUCLIDE_TYPE_PHOTON, NO)) > VALID_PTR)
          {
            /* Set used- and initial-flags */

            SetOption(nuc + NUCLIDE_OPTIONS, OPT_USED);
            SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_INITIAL);

            /* Set pointer */

            elem[n] = nuc;
          }
        else if (nuc == 0)
          Note(0, "Photon interaction data not available for %s (Z = %ld)",
               ZAItoIso(n, 1), n);
      }

  /* Link data */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR)
    {
      /* Check type */

      if ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_PHOTON)
        {
          /* Link photon data */

          if ((ptr = elem[(long)RDB[nuc + NUCLIDE_Z]]) > VALID_PTR)
            WDB[nuc + NUCLIDE_PTR_PHOTON_DATA] = (double)ptr;
        }

      /* Next nuclide */

      nuc = NextItem(nuc);
    }

  /* Print newline */

  fprintf(outp, "\n");
  
  /* Check if no neutron transport mode */

  if ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES)
    return;

  /* Remove neutron data */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR)
    { 
      /* Check type and reset used-flag */
      
      if ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_PHOTON)
        ResetOption(nuc + NUCLIDE_OPTIONS, OPT_USED);

     /* Next nuclide */

      nuc = NextItem(nuc);
    }

  /* Remove unused */

  nuc = (long)RDB[DATA_PTR_NUC0];
  RemoveFlaggedItems(nuc, NUCLIDE_OPTIONS, OPT_USED, NO);
}

/*****************************************************************************/
