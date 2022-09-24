/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : putcompositions.c                              */
/*                                                                           */
/* Created:       2012/05/11 (JLe)                                           */
/* Last modified: 2018/07/17 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Generates independent material compositions for materials    */
/*              divided into Serpent 1 -type regions and depletion zones.    */
/*                                                                           */
/* Comments: - NOTE: tässä ei pidä testat burn-flagia koska toi jälkimmäinen */
/*             jako tehdään myös transport-moodissa.                         */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PutCompositions:"

/*****************************************************************************/

void PutCompositions()
{
  long mat, mat0, loc0, loc1, mem;

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check if material was produced by division */

      if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) < VALID_PTR)
        {
          /* No, get pointer to next */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Check if domain decomposition is in use */

      if (((long)RDB[DATA_DD_DECOMPOSE] == YES) && 
          ((long)RDB[mat + MATERIAL_MPI_ID] > -1) &&
          ((long)RDB[mat + MATERIAL_MPI_ID] != mpiid))
        {
          /* Pointer to next */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Calculate bytes */

      CalculateBytes();

      /* Get memory size */
      
      mem = (long)RDB[DATA_TOTAL_BYTES];

      /* Check pointer to composition */

      if ((long)RDB[mat + MATERIAL_PTR_COMP] > VALID_PTR)
        Die(FUNCTION_NAME, "Divided material already has composition");

      /* Get pointer to original composition */

      loc0 = (long)RDB[mat0 + MATERIAL_PTR_COMP];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      /* Loop over composition and copy */

      while (loc0 > VALID_PTR)
        {
          /* Create new item */

          loc1 = NewItem(mat + MATERIAL_PTR_COMP, COMPOSITION_BLOCK_SIZE);
          
          /* Copy data */
          
          memcpy(&WDB[loc1 + LIST_DATA_SIZE], &RDB[loc0 + LIST_DATA_SIZE], 
                 (COMPOSITION_BLOCK_SIZE - LIST_DATA_SIZE)*sizeof(double));
          
          /* Next */
          
          loc0 = NextItem(loc0);
        }

      /* Original compositions converted to elemental */

      loc0 = (long)RDB[mat0 + MATERIAL_PTR_ORIG_NUC_COMP];
      while (loc0 > VALID_PTR)
        {
          /* Create new item */

          loc1 = NewItem(mat + MATERIAL_PTR_ORIG_NUC_COMP, 
                         COMPOSITION_BLOCK_SIZE);
          
          /* Copy data */
          
          memcpy(&WDB[loc1 + LIST_DATA_SIZE], &RDB[loc0 + LIST_DATA_SIZE], 
                 (COMPOSITION_BLOCK_SIZE - LIST_DATA_SIZE)*sizeof(double));
          
          /* Next */
          
          loc0 = NextItem(loc0);
        }

      /* Put additional options */

      SetOption(mat + MATERIAL_OPTIONS, (long)RDB[mat0 + MATERIAL_OPTIONS]);
      SetOption(mat + MATERIAL_OPTIONS, OPT_PHYSICAL_MAT);

      /* Copy physical parameters */

      WDB[mat + MATERIAL_ADENS] = RDB[mat0 + MATERIAL_ADENS];
      WDB[mat + MATERIAL_MDENS] = RDB[mat0 + MATERIAL_MDENS];
      WDB[mat + MATERIAL_INI_FISS_MDENS] = RDB[mat0 + MATERIAL_INI_FISS_MDENS];

      /* Copy poison flags */

      WDB[mat + MATERIAL_XENON_EQUIL_CALC] = 
        RDB[mat0 + MATERIAL_XENON_EQUIL_CALC];
      WDB[mat + MATERIAL_SAMARIUM_EQUIL_CALC] = 
        RDB[mat0 + MATERIAL_SAMARIUM_EQUIL_CALC];

      /* Copy TMS limits */

      WDB[mat + MATERIAL_TMS_TMIN] = RDB[mat0 + MATERIAL_TMS_TMIN];
      WDB[mat + MATERIAL_TMS_TMAX] = RDB[mat0 + MATERIAL_TMS_TMAX];

      /* Copy TMS flag */

      WDB[mat + MATERIAL_TMS_MODE] = RDB[mat0 + MATERIAL_TMS_MODE];

      /* Calculate bytes */

      CalculateBytes();

      /* Update memory size */
      
      WDB[mat + MATERIAL_MEM_SIZE] = RDB[mat + MATERIAL_MEM_SIZE] + 
        RDB[DATA_TOTAL_BYTES] - (double)mem;

      /* Next material */

      mat = NextItem(mat);
    }

  /* Find I-135, Xe-135, Pm-149 and Sm-149 for iteration and poison */
  /* calculation */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
       /* Loop over composition and find Xe-135 */

      loc0 = (long)RDB[mat + MATERIAL_PTR_COMP];
      while (loc0 > VALID_PTR)
        {
          /* Pointer to nuclide */

          loc1 = (long)RDB[loc0 + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          /* Check ZAI */

          if ((long)RDB[loc1 + NUCLIDE_ZAI] == 531350)
            WDB[mat + MATERIAL_PTR_I135_ISO] = (double)loc0;
          else if ((long)RDB[loc1 + NUCLIDE_ZAI] == 541350)
            WDB[mat + MATERIAL_PTR_XE135_ISO] = (double)loc0;
          else if (((long)RDB[loc1 + NUCLIDE_ZAI] == 541351) &&
                   ((long)RDB[DATA_OPTI_POISON_CALC_XE135M] == YES))
            WDB[mat + MATERIAL_PTR_XE135M_ISO] = (double)loc0;
          else if ((long)RDB[loc1 + NUCLIDE_ZAI] == 611470)
            WDB[mat + MATERIAL_PTR_PM147_ISO] = (double)loc0;
          else if ((long)RDB[loc1 + NUCLIDE_ZAI] == 611480)
            WDB[mat + MATERIAL_PTR_PM148_ISO] = (double)loc0;
          else if ((long)RDB[loc1 + NUCLIDE_ZAI] == 611481)
            WDB[mat + MATERIAL_PTR_PM148M_ISO] = (double)loc0;
          else if ((long)RDB[loc1 + NUCLIDE_ZAI] == 611490)
            WDB[mat + MATERIAL_PTR_PM149_ISO] = (double)loc0;
          else if ((long)RDB[loc1 + NUCLIDE_ZAI] == 621490)
            WDB[mat + MATERIAL_PTR_SM149_ISO] = (double)loc0;

          /* Next */

          loc0 = NextItem(loc0);
        }

      /* Check pointers */

      if ((long)RDB[mat + MATERIAL_XENON_EQUIL_CALC] == YES)
        {
          if ((long)RDB[mat + MATERIAL_PTR_I135_ISO] < VALID_PTR)
            Die(FUNCTION_NAME, "I-135 not found");
          if ((long)RDB[mat + MATERIAL_PTR_XE135_ISO] < VALID_PTR)
            Die(FUNCTION_NAME, "Xe-135 not found");
        }

      if ((long)RDB[mat + MATERIAL_SAMARIUM_EQUIL_CALC] == YES)
        {
          if ((long)RDB[mat + MATERIAL_PTR_PM149_ISO] < VALID_PTR)
            Die(FUNCTION_NAME, "Pm-149 not found");
          if ((long)RDB[mat + MATERIAL_PTR_SM149_ISO] < VALID_PTR)
            Die(FUNCTION_NAME, "Sm-149 not found");
        }
  
      /* Next material */

      mat = NextItem(mat);
    }
}

/*****************************************************************************/
