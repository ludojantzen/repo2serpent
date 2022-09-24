/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processmixtures.c                              */
/*                                                                           */
/* Created:       2011/01/20 (JLe)                                           */
/* Last modified: 2015/11/02 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Converts mixtures into materials                             */
/*                                                                           */
/* Comments: - Burnable materials should not be allowed?                     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessMixtures:"

/*****************************************************************************/

void ProcessMixture(long mat, long recu)
{
  long mix, ptr, iso;
  double mem, sum, vol, mass, dens, f;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

  /* Get pointer to mixture */

  if ((mix = (long)RDB[mat + MATERIAL_PTR_MIX]) < VALID_PTR)
    return;

  /* Avoid processing divided mixtures (17.1.2014 / 2.1.17) */

  if ((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] > VALID_PTR)
    return;

  /* Check burn flag */

  if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
    Die(FUNCTION_NAME, "Burn flag in mixture %s", 
        GetText(mat + MATERIAL_PTR_NAME));

  /* Call recursively if mixture contains another mixture */

  while (mix > VALID_PTR)
    {
      /* Pointer to material */

      ptr = (long)RDB[mix + MIXTURE_PTR_MAT];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Compare pointer */

      if ((ptr == mat) || (recu > 10000))
        Error(mat, "Mixture diluted with itself (this is not homeopathy)");

      /* Check burn flag */
      
      if ((long)RDB[ptr + MATERIAL_OPTIONS] & OPT_BURN_MAT)
        Error(mat, "Burnable materials cannot be used in mixtures");

      /* Call recursively */
      
      ProcessMixture(ptr, recu + 1);

      /* Next */

      mix = NextItem(mix);
    }

  /***************************************************************************/

  /***** Create composition **************************************************/

  /* Check if composition exists */

  if ((long)RDB[mat + MATERIAL_PTR_COMP] < VALID_PTR)
    {
      /* Get memory size */

      mem = RDB[DATA_TOTAL_BYTES];

      /* Reset volume and mass */

      vol = 0.0;
      mass = 0.0;

      /* Loop over mixture and allocate memory for composition */

      mix = (long)RDB[mat + MATERIAL_PTR_MIX];
      while (mix > VALID_PTR)
        {
          /* Loop over composition */

          ptr = (long)RDB[mix + MIXTURE_PTR_MAT];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          ptr = (long)RDB[ptr + MATERIAL_PTR_COMP];
          while (ptr > VALID_PTR)
            {
              /* Allocate memory */

              iso = NewItem(mat + MATERIAL_PTR_COMP, COMPOSITION_BLOCK_SIZE);

              /* Next nuclide */

              ptr = NextItem(ptr);
            }

          /* Pointer to material */

          ptr = (long)RDB[mix + MIXTURE_PTR_MAT];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          
          /* Check fraction and add to sum */

          if (vol*RDB[mix + MIXTURE_VFRAC] < 0.0)
            Error(mat, "Mixed volume and mass fractions");
          else if (RDB[mix + MIXTURE_VFRAC] > 0.0)
            {
              vol = vol + RDB[mix + MIXTURE_VFRAC];
              mass = mass + 
                RDB[mix + MIXTURE_VFRAC]*RDB[ptr + MATERIAL_MDENS];
            }
          else
            {
              vol = vol + RDB[mix + MIXTURE_VFRAC]/RDB[ptr + MATERIAL_MDENS];
              mass = mass - RDB[mix + MIXTURE_VFRAC];
            }

          /* Next material */

          mix = NextItem(mix);
        }

      /* Check sum */

      if (vol == 0.0)
        Error(mat, "No nonzero fractions");

      /* Normalize fractions */

      mix = (long)RDB[mat + MATERIAL_PTR_MIX];
      while (mix > VALID_PTR)
        {
          /* Pointer to material */

          ptr = (long)RDB[mix + MIXTURE_PTR_MAT];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          if (vol > 0)
            {
              WDB[mix + MIXTURE_MFRAC] = 
                RDB[mix + MIXTURE_VFRAC]*RDB[ptr + MATERIAL_MDENS]/mass;
              WDB[mix + MIXTURE_VFRAC] = RDB[mix + MIXTURE_VFRAC]/vol;
            }
          else
            {
              WDB[mix + MIXTURE_MFRAC] = -RDB[mix + MIXTURE_VFRAC]/mass;
              WDB[mix + MIXTURE_VFRAC] = 
                RDB[mix + MIXTURE_VFRAC]/RDB[ptr + MATERIAL_MDENS]/vol;
            }

          /* Next material */

          mix = NextItem(mix);
        }

      /* Update memory size */
  
      WDB[mat + MATERIAL_MEM_SIZE] = RDB[mat + MATERIAL_MEM_SIZE] 
        + RDB[DATA_TOTAL_BYTES] - mem;
    }

  /***************************************************************************/

  /***** Calculate material density ******************************************/

  /* Reset sum */

  sum = 0.0;

  /* Pointer to composition */

  iso = (long)RDB[mat + MATERIAL_PTR_COMP];

  /* Loop over mixture */

  mix = (long)RDB[mat + MATERIAL_PTR_MIX];
  while (mix > VALID_PTR)
    {
      /* Check composition pointer */

      CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

      /* Get pointer to material */

      ptr = (long)RDB[mix + MIXTURE_PTR_MAT];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Check fraction */

      if (RDB[mix + MIXTURE_VFRAC] < 0.0)
        Die(FUNCTION_NAME, "Negative fraction");
      
      /* Add atomic density to sum */

      sum = sum + RDB[mix + MIXTURE_VFRAC]*RDB[ptr + MATERIAL_ADENS];

      /* Loop over composition */
      
      ptr = (long)RDB[ptr + MATERIAL_PTR_COMP];
      while (ptr > VALID_PTR)
        {
          /* Copy nuclide pointer */
          
          WDB[iso + COMPOSITION_PTR_NUCLIDE] = 
            RDB[ptr + COMPOSITION_PTR_NUCLIDE];
          
          /* Set atomic density */
          
          WDB[iso + COMPOSITION_ADENS] = 
            RDB[mix + MIXTURE_VFRAC]*RDB[ptr + COMPOSITION_ADENS];
          
          /* Next */
          
          iso = NextItem(iso);
          ptr = NextItem(ptr);
        }
      
      /* Next material */

      mix = NextItem(mix);
    }
  
  /* Save original density (needed if overridden by IFC) */

  dens = RDB[mat + MATERIAL_ADENS];

  /* Put density */

  WDB[mat + MATERIAL_ADENS] = sum;

  /***************************************************************************/

  /***** Combine nuclides ****************************************************/

  /* Loop over nuclides */

  iso = (long)RDB[mat + MATERIAL_PTR_COMP];
  while(iso > VALID_PTR)
    {
      /* Loop over remaining */

      ptr = NextItem(iso);
      while (ptr > VALID_PTR)
        {
          /* Compare pointers */

          if (RDB[iso + COMPOSITION_PTR_NUCLIDE] == 
              RDB[ptr + COMPOSITION_PTR_NUCLIDE])
            {
              /* Add to density */

              WDB[iso + COMPOSITION_ADENS] = RDB[iso + COMPOSITION_ADENS] 
                + RDB[ptr + COMPOSITION_ADENS];

              /* Reset */

              WDB[ptr + COMPOSITION_ADENS] = 0.0;
            }

          /* Next nuclide */

          ptr = NextItem(ptr);
        }

      /* Next nuclide */

      iso = NextItem(iso);
    }

  /* Calculate isotope fractions */

  IsotopeFractions(mat);

  /***************************************************************************/
  
  /***** Put TMS data ********************************************************/

  if (recu == 0)
    {
      /* Loop over constituent materials */

      mix = (long)RDB[mat + MATERIAL_PTR_MIX];
      while (mix > VALID_PTR)
        {
          /* Pointer to material */

          ptr = (long)RDB[mix + MIXTURE_PTR_MAT];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Set TMS flag and copy parameters */
                
          WDB[mat + MATERIAL_TMS_MODE] = RDB[ptr + MATERIAL_TMS_MODE];
          WDB[mat + MATERIAL_TMS_TMIN] = RDB[ptr + MATERIAL_TMS_TMIN];
          WDB[mat + MATERIAL_TMS_TMAX] = RDB[ptr + MATERIAL_TMS_TMAX];

          /* Next */

          mix = NextItem(mix);
        }

      /* Calculate density if overridden */

      if (dens != 0.0)
        {
          /* Calculate adjustment fraction */

          if (dens > 0.0)
            f = dens/RDB[mat + MATERIAL_ADENS];
          else
            f = -dens/RDB[mat + MATERIAL_MDENS];

          /* Check value */

          CheckValue(FUNCTION_NAME, "f", "", f, ZERO, 100.0);

          /* Loop over composition */
          
          iso = (long)RDB[mat + MATERIAL_PTR_COMP];
          while (iso > VALID_PTR)
            {
              /* Adjust density */

              WDB[iso + COMPOSITION_ADENS] = f*RDB[iso + COMPOSITION_ADENS];

              /* Next isotope */
              
              iso = NextItem(iso);
            }
          
          /* Adjust material densities */

          WDB[mat + MATERIAL_MDENS] = f*RDB[mat + MATERIAL_MDENS];
          WDB[mat + MATERIAL_ADENS] = f*RDB[mat + MATERIAL_ADENS];
        }
    }

  /***************************************************************************/
  
  /***** Check poison data ***************************************************/

  if (recu == 0)
    {
      /* Loop over constituent materials */

      mix = (long)RDB[mat + MATERIAL_PTR_MIX];
      while (mix > VALID_PTR)
        {
          /* Pointer to material */

          ptr = (long)RDB[mix + MIXTURE_PTR_MAT];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          
          /* Check flags */
          
          if ((long)RDB[ptr + MATERIAL_XENON_EQUIL_CALC] == YES)
            Error(mat, 
                  "Fissile mixtures not allowed in equilibrium xenon mode");

          if ((long)RDB[ptr + MATERIAL_SAMARIUM_EQUIL_CALC] == YES)
            Error(mat, 
                  "Fissile mixtures not allowed in equilibrium samarium mode");

          /* Next */

          mix = NextItem(mix);
        }
    }

  /***************************************************************************/
}

/*****************************************************************************/
