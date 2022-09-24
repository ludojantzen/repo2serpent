/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : normalizecompositions.c                        */
/*                                                                           */
/* Created:       2011/11/07 (JLe)                                           */
/* Last modified: 2015/06/06 (JLe)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Normalizes material compositions before transport cycle      */
/*                                                                           */
/* Comments: - Nuklidit joilla on S(a,b) dataa eivät välttämättä kuulu       */
/*             mihinkään materiaaliin (kopioidaan addsabdata.c:ssä). Jos     */
/*             tällaisella nuklidilla on URES-dataa, se aiheuttaa errorin    */
/*             normalizecompositions.c:ssä jos noita minimi ja maksimi       */
/*             konsentraatioita ei aseteta erikseen. (käsittelyä muutettiin  */
/*             myöhemmin versiossa 2.1.14, en tiedä onko tuo enää ongelma).  */
/*                                                                           */
/*           - Tää ei enää tee sitä mitä aliohjelman nimi antaa ymmärtää.    */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "NormalizeCompositions:"

/*****************************************************************************/

void NormalizeCompositions()
{
  long mat, nuc, iso, ptr;
  double f;

  /* Reset nuclide minimum and maximum fractions */
  
  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR)
    {
      /* Reset values */

      WDB[nuc + NUCLIDE_MIN_AFRAC] = INFTY;
      WDB[nuc + NUCLIDE_MAX_AFRAC] = -INFTY;
      
      /* Next nuclide */

      nuc = NextItem(nuc);
    }

  /* Normalize material compositions */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Calculate fractions (tätä ei pidä eikä saa kutsua?) */
      /*
      IsotopeFractions(mat);
      */
      /* Next material */

      mat = NextItem(mat);
    }

  /* Find minimum and maximum fractions */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check zero density and continuous reprocessing */

      if ((RDB[mat + MATERIAL_ADENS] < ZERO) &&
          ((long)RDB[mat + MATERIAL_FLOW_PTR_FIRST] > VALID_PTR))
        {
          /* Pointer to next */
          
          mat = NextItem(mat);
          
          /* Cycle loop */
          
          continue;
        }

      /* Check divided */
      
      if ((long)RDB[mat + MATERIAL_DIV_TYPE] != MAT_DIV_TYPE_PARENT)
        {
          /* Loop over composition */

          iso = (long)RDB[mat + MATERIAL_PTR_COMP];
          while (iso > VALID_PTR)
            {
              /* Pointer to nuclide data */
              
              nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);
              
              /* Calculate atomic fraction */
              
              f = RDB[iso + COMPOSITION_ADENS]/RDB[mat + MATERIAL_ADENS];
              CheckValue(FUNCTION_NAME, "f", "", f, 0.0, 1.0);
              
              /* Compare to minimum */
              
              if (f < RDB[nuc + NUCLIDE_MIN_AFRAC])
                WDB[nuc + NUCLIDE_MIN_AFRAC] = f;
              
              /* Compare to minimum */
              
              if (f > RDB[nuc + NUCLIDE_MAX_AFRAC])
                WDB[nuc + NUCLIDE_MAX_AFRAC] = f;
              
              /* Nuclides linked with S(a,b) data */
              
              if ((ptr = (long)RDB[nuc + NUCLIDE_PTR_SAB_NUC]) > VALID_PTR)
                {
                  /* Compare to minimum */
                  
                  if (f < RDB[ptr + NUCLIDE_MIN_AFRAC])
                    WDB[ptr + NUCLIDE_MIN_AFRAC] = f;
                  
                  /* Compare to minimum */
                  
                  if (f > RDB[ptr + NUCLIDE_MAX_AFRAC])
                    WDB[ptr + NUCLIDE_MAX_AFRAC] = f;
                }
              
              /* Next nuclide */
              
              iso = NextItem(iso);
            }
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /* Check maximum and minimum fractions */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR)
    {
      /* Check if not set */
          
      if (RDB[nuc + NUCLIDE_MIN_AFRAC] == INFTY)
        WDB[nuc + NUCLIDE_MIN_AFRAC] = 0.0;
      if (RDB[nuc + NUCLIDE_MAX_AFRAC] == -INFTY)
        WDB[nuc + NUCLIDE_MAX_AFRAC] = 0.0;
      
      /* Check minimum (account for numerical round-off error) */
      
      if ((RDB[nuc + NUCLIDE_MIN_AFRAC] - 1.0 > 1E-12) || 
          (RDB[nuc + NUCLIDE_MIN_AFRAC] < 0.0))
        Die(FUNCTION_NAME, "nuclide %s min afrac = %E\n",
            GetText(nuc + NUCLIDE_PTR_NAME), RDB[nuc + NUCLIDE_MIN_AFRAC]);
      
      /* Check maximum (account for numerical round-off error) */
      
      if ((RDB[nuc + NUCLIDE_MAX_AFRAC] - 1.0 > 1E-12) || 
          (RDB[nuc + NUCLIDE_MAX_AFRAC] < 0.0))
        Die(FUNCTION_NAME, "nuclide %s max afrac = %E\n",
            GetText(nuc + NUCLIDE_PTR_NAME), RDB[nuc + NUCLIDE_MAX_AFRAC]);
      
      /* Check against each other */
      
      if (RDB[nuc + NUCLIDE_MIN_AFRAC] > RDB[nuc + NUCLIDE_MAX_AFRAC])
        Die(FUNCTION_NAME, "nuclide %s min / max afrac = %E / %E\n",
            GetText(nuc + NUCLIDE_PTR_NAME), RDB[nuc + NUCLIDE_MIN_AFRAC],
            RDB[nuc + NUCLIDE_MAX_AFRAC]);
      
      /* Next nuclide */

      nuc = NextItem(nuc);
    }
}

/*****************************************************************************/
