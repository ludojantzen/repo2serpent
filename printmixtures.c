/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printmixtures.c                                */
/*                                                                           */
/* Created:       2016/12/10 (JLe)                                           */
/* Last modified: 2017/11/29 (JLe)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Prints decomposed mixtures into material cards               */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrintMixtures:"

/*****************************************************************************/

void PrintMixtures()
{
  long mat, mix, ptr, iso, nuc, nm;
  char outfile[MAX_STR];
  FILE *fp;

  /* Check mode */

  if ((long)RDB[DATA_DECOMPOSE_MIXTURES] == NO)
    return;

  /* Count mixtures */

  nm = 0;

  /* Loop over materials */
  
  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Get pointer to mixture */

      if ((long)RDB[mat + MATERIAL_PTR_MIX] > VALID_PTR)
        nm++;
      
      /* Next material */
      
      mat = NextItem(mat);
    }

  /* Check count */

  if (nm == 0)
    {
      /* No mixtures */

      fprintf(outp, "No mixtures.\n\n");

      /* Terminate run */

      exit(-1);
    }

 /* Set output file */
  
  sprintf(outfile, "%s.mix", GetText(DATA_PTR_INPUT_FNAME));

  /* Print */

  fprintf(outp, "Writing compositions of %ld mixtures in \"%s\"...\n", nm,
          outfile);

  /* Open file for writing */

  if ((fp = fopen(outfile, "w")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file for writing");

  /* Loop over materials */
  
  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Get pointer to mixture */

      if ((mix = (long)RDB[mat + MATERIAL_PTR_MIX]) < VALID_PTR)
        {
          /* Next material */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Print */

      fprintf(fp, "\n%% Material %s is a mixture of %ld components:\n\n",
              GetText(mat + MATERIAL_PTR_NAME), ListSize(mix));
      
      fprintf(fp, "%% -----------------------------------------\n");
      fprintf(fp, "%% Material             v. frac      m. frac\n"); 
      fprintf(fp, "%% -----------------------------------------\n");

      /* Loop over mixture */
      
      while (mix > VALID_PTR)
        {
          /* Pointer to material */

          ptr = (long)RDB[mix + MIXTURE_PTR_MAT];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          fprintf(fp, "%% %-15s %12.5E %12.5E\n", 
                  GetText(ptr + MATERIAL_PTR_NAME),
                  RDB[mix + MIXTURE_VFRAC], 
                  RDB[mix + MIXTURE_MFRAC]);
          
          /* Next material */
          
          mix = NextItem(mix);
        }

      fprintf(fp, "%% -----------------------------------------\n\n");
      
      /* Print */

      fprintf(fp, "mat %s  %1.5E", GetText(mat + MATERIAL_PTR_NAME),
              RDB[mat + MATERIAL_ADENS]);

      /* Loop over composition and print S(a,b) entries */
      
      iso = (long)RDB[mat + MATERIAL_PTR_COMP];  
      while (iso > VALID_PTR)
        {
          /* Pointer to nuclide */
          
          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          if ((ptr = (long)RDB[nuc + NUCLIDE_PTR_ORIG_THERM]) > VALID_PTR)
            fprintf(fp, " moder %s %ld", GetText(ptr + THERM_PTR_ALIAS), 
                    (long)RDB[ptr + THERM_ZA]);

          /* Next isotope */
          
          iso = NextItem(iso);
        }

       /* Newline */

      fprintf(fp, "\n\n");

      /* Loop over composition */
      
      iso = (long)RDB[mat + MATERIAL_PTR_COMP];  
      while (iso > VALID_PTR)
        {
          /* Pointer to nuclide */
          
          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Check if S(a,b) data */
          
          if ((ptr = (long)RDB[nuc + NUCLIDE_SAB_PTR_FREE]) > VALID_PTR)
            nuc = ptr;
          
          /* Print density */
          
          if (RDB[iso + COMPOSITION_ADENS] > 0.0)
            fprintf(fp, "%10s %12.5E\n", GetText(nuc + NUCLIDE_PTR_NAME),
                    RDB[iso + COMPOSITION_ADENS]);
                     
          /* Next isotope */
          
          iso = NextItem(iso);
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /* Close file */

  fclose(fp);

  /* Exit OK */

  fprintf(outp, "OK.\n\n");

  /* Terminate run */

  exit(-1);
}

/*****************************************************************************/
