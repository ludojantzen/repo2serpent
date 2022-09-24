/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : isotopefractions.c                             */
/*                                                                           */
/* Created:       2010/12/28 (JLe)                                           */
/* Last modified: 2015/04/02 (JLe)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Calculates isotope fractions and material densities          */
/*                                                                           */
/* Comments: - Kaikkia kombinaatioita ei oo vielä testattu                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "IsotopeFractions:"

/*****************************************************************************/

void IsotopeFractions(long mat)
{
  long iso, nuc;
  double sum, adens, mass;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "mat", DATA_ARRAY, mat);

  /* Materials that were produced by division do not yet have composition */
  /* lists, so they must be skipped. */

  if ((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] > VALID_PTR)
    return;

  /* Check pointer to composition */

  if ((long)RDB[mat + MATERIAL_PTR_COMP] < VALID_PTR)
    {
      /* Check pointer to mixture */

      if ((long)RDB[mat + MATERIAL_PTR_MIX] > VALID_PTR)
        return;
      else
        Error(mat, "Material %s has no composition", 
              GetText(mat + MATERIAL_PTR_NAME));
    }

  /***************************************************************************/

  /***** Normalize composition ***********************************************/

  /* Loop until composition is atomic (once if atomic fractions or densities */
  /* given, twice if mass fractions or densities given). */

  do
    {
      /* Calculate sum */

      sum = 0.0;

      /* Loop over composition */

      iso = (long)RDB[mat + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
        {
          /* Pointer to nuclide */
          
          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Set fissile material flag (this is set for the first time */
          /* in processnuclide.c) */

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_FISSILE)
            SetOption(mat + MATERIAL_OPTIONS, OPT_FISSILE_MAT);
          
          /* Get atomic density */

          adens = RDB[iso + COMPOSITION_ADENS];
          CheckValue(FUNCTION_NAME, "adens", "", adens, -INFTY, INFTY);

          /* Check type and add to sum */
          
          if ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_SAB)
            sum = sum + adens;
          
          /* Next isotope */
          
          iso = NextItem(iso);
        }
      
      /* Check sum */
      
      if (sum == 0.0)
        {
          /* Check reprocessing mode */
       
          if ((long)RDB[mat + MATERIAL_FLOW_PTR_FIRST] < VALID_PTR)
            Error(mat, "Zero total density in material %s ", 
                  GetText(mat + MATERIAL_PTR_NAME));
        }
      else
        {
          /* Normalise fractions */
          
          iso = (long)RDB[mat + MATERIAL_PTR_COMP];
          while (iso > VALID_PTR)
            {
              /* Pointer to nuclide */
              
              nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);
              
              /* Divide by sum */
              
              if ((WDB[iso + COMPOSITION_ADENS] = 
                   RDB[iso + COMPOSITION_ADENS]/sum) < 0.0)
                Error(mat, "Mixed atomic and mass fractions in material %s",
                      GetText(mat + MATERIAL_PTR_NAME));
              
              /* Convert to atomic fraction if mass fraction given */
              
              if (sum < 0.0)
                {
                  /* Atomic weight is zero for lost nuclide */
                  
                  if (RDB[nuc + NUCLIDE_AW] > 0.0)
                    WDB[iso + COMPOSITION_ADENS] = 
                      RDB[iso + COMPOSITION_ADENS]/RDB[nuc + NUCLIDE_AW];
                  else
                    WDB[iso + COMPOSITION_ADENS] = 0.0;                
                }
              
              /* Next isotope */
              
              iso = NextItem(iso);
            }
        }
    }
  while (sum < 0.0);

  /***************************************************************************/

  /***** Calculate material densities ****************************************/

  sum = 0.0;

  /* Loop over composition */
  
  iso = (long)RDB[mat + MATERIAL_PTR_COMP];
  while (iso > VALID_PTR)
    {
      /* Pointer to nuclide */
      
      nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);
    
      /* Check type and add to sum */
      
      if ((long)WDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_SAB)
        sum = sum + RDB[iso + COMPOSITION_ADENS]*RDB[nuc + NUCLIDE_AW];
      
      /* Next isotope */
      
      iso = NextItem(iso);
    }

  /* Check given density */

  if (RDB[mat + MATERIAL_ADENS] > 0.0)
    {
      /* Put mass density */

      WDB[mat + MATERIAL_MDENS] = sum*RDB[mat + MATERIAL_ADENS]/N_AVOGADRO;
    }
  else if (RDB[mat + MATERIAL_ADENS] < 0.0)
    {
      /* Put mass density */

      WDB[mat + MATERIAL_MDENS] = -RDB[mat + MATERIAL_ADENS];

      /* Put atomic density */

      WDB[mat + MATERIAL_ADENS] = N_AVOGADRO*RDB[mat + MATERIAL_MDENS]/sum;

    }
  else if ((long)RDB[mat + MATERIAL_FLOW_PTR_FIRST] < VALID_PTR)
    Error(mat, "Zero density in material %s", 
          GetText(mat + MATERIAL_PTR_NAME));

  /***************************************************************************/

  /***** Convert atomic fractions to densities *******************************/

  /* Loop over composition */

  iso = (long)RDB[mat + MATERIAL_PTR_COMP];
  while (iso > VALID_PTR)
    {
      /* Multiply by atomic density */
      
      WDB[iso + COMPOSITION_ADENS] = 
        RDB[iso + COMPOSITION_ADENS]*RDB[mat + MATERIAL_ADENS];

      /* Next isotope */
          
      iso = NextItem(iso);
    }

  /***************************************************************************/

  /***** Convert given mass to given volume **********************************/

  if ((mass = RDB[mat + MATERIAL_MASS_GIVEN]) > 0.0)
    {
      /* Check that volume is not given */

      if (RDB[mat + MATERIAL_VOLUME_GIVEN] > 0.0)
        Error(mat, "Both mass and volume given in input");

      /* Convert */

      if (RDB[mat + MATERIAL_MDENS] > 0.0)
        {
          WDB[mat + MATERIAL_VOLUME_GIVEN] = mass/RDB[mat + MATERIAL_MDENS];
          WDB[mat + MATERIAL_MASS_GIVEN] = 0.0;
        }
      else
        Die(FUNCTION_NAME, "Zero mass density");
    }

  /***************************************************************************/

  /***** Calculate initial fissile mass density ******************************/

  /* Loop over composition */

  iso = (long)RDB[mat + MATERIAL_PTR_COMP];
  while (iso > VALID_PTR)
    {
      /* Get pointer to nuclide */

      nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Check fissile flag and add to mass */

      if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_FISSILE)
        WDB[mat + MATERIAL_INI_FISS_MDENS] = RDB[mat + MATERIAL_INI_FISS_MDENS]
          + RDB[iso + COMPOSITION_ADENS]*RDB[nuc + NUCLIDE_AW]/N_AVOGADRO;

      /* Next isotope */
          
      iso = NextItem(iso);
    }

  /***************************************************************************/
}

/*****************************************************************************/
