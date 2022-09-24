/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printgammaspectra.c                            */
/*                                                                           */
/* Created:       2017/03/02 (JLe)                                           */
/* Last modified: 2017/03/12 (JLe)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Prints gamma spectra for different materials                 */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrintGammaSpectra:"

/*****************************************************************************/

void PrintGammaSpectra()
{
  long mat, nuc, src, ptr, rad, erg, n, ne;
  char outfile[MAX_STR];
  double cum, *spec1, *spec2;
  FILE *fp;

  fprintf(outp, "Printing photon spectra in file...\n");

  /* File name */

  sprintf(outfile, "%s_gsrc.m", GetText(DATA_PTR_INPUT_FNAME));

  /* Open file */

  if ((fp = fopen(outfile, "w")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file for writing");

  /* Print columns */

  fprintf(fp, "\n%% Column 1: Nuclide ZAI\n");
  fprintf(fp, "%% Column 2: Specific intensity (photons per decay)\n");
  fprintf(fp, "%% Column 3: Total emission rate (photons/s)\n");
  fprintf(fp, "%% Column 4: Cumulative fraction of material total\n");
  fprintf(fp, "%% Column 5: Emission line energy\n");
  fprintf(fp, "%% Column 6: Relative intensity (photons per decay)\n");
  fprintf(fp, "%% Column 7: Cumulative fraction of nuclide total\n");

  /* Newline */

  fprintf(fp, "\n");

  /* Check pointer to energy group structure */

  if ((erg = (long)RDB[DATA_GSPEC_PTR_EGRID]) > VALID_PTR)
    {
      /* Number of energy groups */

      ne = (long)RDB[erg + ENERGY_GRID_NE];
      
      /* Allocate memory for spectra */

      spec1 = Mem(MEM_ALLOC, ne - 1, sizeof(double));
      spec2 = Mem(MEM_ALLOC, ne - 1, sizeof(double));
    }
  else
    {
      /* Avoid compiler warning */

      ne = 0;
      spec1 = NULL;
      spec2 = NULL;
    }

  /* Loop over materials */
      
  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check decay source */

      if ((src = (long)RDB[mat + MATERIAL_PTR_PHOTON_DECAY_SRC]) < VALID_PTR)
        {
          /* Next material */
          
          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }
      
      /* Check if spectra is calculated and reset spectra */
      
      if (erg > VALID_PTR)
        for (n = 0; n < ne; n++)
          {
            spec1[n] = 0.0;
            spec2[n] = 0.0;
          }

      /* Print material name */

      fprintf(fp, "mat_%s = [\n", GetText(mat + MATERIAL_PTR_NAME));

      /* Loop over source */

      while (src > VALID_PTR)
        {
          /* Pointer to nuclide */

          nuc = (long)RDB[src + SRC_DECCAY_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Loop over radiations to find gamma emission */

          rad = (long)RDB[nuc + NUCLIDE_PTR_RADIATIONS];
          while (rad > VALID_PTR)
            {
              /* Check type */

              if ((long)RDB[rad + NUCLIDE_RAD_TYPE] == PARTICLE_TYPE_GAMMA)
                break;

              /* Next */

              rad = NextItem(rad);
            }

          /* Check pointer */

          if (rad < VALID_PTR)
            Die(FUNCTION_NAME, "No gamma radiation");

          /* Pointer to line spectra */

          ptr = (long)RDB[src + SRC_DECCAY_PTR_SPEC];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          fprintf(fp, "\n%% --- %s :\n\n", 
                  ZAItoIso((long)RDB[nuc + NUCLIDE_ZAI], 1));
                                
          /* Reset sum */

          cum = 0.0;
          
          /* Loop over spectra */

          while (ptr > VALID_PTR)
            {
              /* Check spectrum and particle type */

              if (((long)RDB[ptr + DECAY_SPEC_PARTICLE] != PARTICLE_TYPE_GAMMA)
                  || ((long)RDB[ptr + DECAY_SPEC_TYPE] != DECAY_SPEC_LINE))
                {
                  /* Next */

                  ptr = NextItem(ptr);

                  /* Cycle loop */

                  continue;
                }

              /* Print ZAI, specific and total intensity and fraction */
              
              fprintf(fp, "%6ld %11.5E %11.5E %11.5E ", 
                      (long)RDB[nuc + NUCLIDE_ZAI],
                      RDB[rad + NUCLIDE_RAD_SPEC_I], 
                      RDB[src + SRC_DECCAY_I], 
                      RDB[src + SRC_DECCAY_CUM_P]);

              /* Add to cumulative */

              cum = cum + RDB[ptr + DECAY_SPEC_RI];

              /* Print energy, relative intensity and total */

              fprintf(fp, "%11.5E %11.5E %11.5E", 
                      RDB[ptr + DECAY_SPEC_LINE_E],
                      RDB[ptr + DECAY_SPEC_RI],
                      cum/RDB[rad + NUCLIDE_RAD_SPEC_I]);

              /* Newline */

              fprintf(fp, "\n");

              /* Check if spectra is calculated and get index */
              
              if (erg > VALID_PTR)
                if ((n = GridSearch(erg, RDB[ptr + DECAY_SPEC_LINE_E])) > -1)
                  {
                    /* Check */
                    
                    CheckValue(FUNCTION_NAME, "n", "", n, 0, ne - 2);

                    /* Store */

                    spec1[n] = spec1[n] + RDB[ptr + DECAY_SPEC_RI]*
                      RDB[src + SRC_DECCAY_I]/RDB[rad + NUCLIDE_RAD_SPEC_I];
                    spec2[n] = spec2[n] + RDB[ptr + DECAY_SPEC_RI]*
                      RDB[ptr + DECAY_SPEC_LINE_E]*RDB[src + SRC_DECCAY_I]/
                      RDB[rad + NUCLIDE_RAD_SPEC_I];
                  }
              
              /* Next line */

              ptr = NextItem(ptr);
            }

          /* Next */

          src = NextItem(src);
        }

      fprintf(fp, "];\n\n");
      
      /* Print volume and total */

      fprintf(fp, "mat_%s_vol = %11.5E;\n", GetText(mat + MATERIAL_PTR_NAME),
              RDB[mat + MATERIAL_VOLUME]);

      fprintf(fp, "mat_%s_tot = %11.5E;\n\n", GetText(mat + MATERIAL_PTR_NAME),
              RDB[mat + MATERIAL_PHOTON_DEC_SRC_RATE]);

      /* Check if spectra is calculated and print spectra */
      
      if (erg > VALID_PTR)
        {

          fprintf(fp, "mat_%s_gspec = [\n", GetText(mat + MATERIAL_PTR_NAME));

          for (n = 0; n < ne - 1; n++)
            fprintf(fp, "%E %E\n", spec1[n], spec2[n]);
          
          fprintf(fp, "];\n\n");                  
        }

      /* Next material */

      mat = NextItem(mat);
    }      

  /* Print total */

  fprintf(fp, "tot = %11.5E;\n\n", RDB[DATA_TOT_PHOTON_DEC_SRC_RATE]);

  /* Check energy group structure */
  
  if (erg > VALID_PTR)
    {
      /* Pointer to data */

      ptr = (long)RDB[erg + ENERGY_GRID_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Loop over data and print */
      
      fprintf(fp, "\nEg = [\n");

      for (n = 0; n < ne; n++)
        fprintf(fp, "%E\n", RDB[ptr + n]);

      fprintf(fp, "];\n");

      /* Free memory */

      Mem(MEM_FREE, spec1);
      Mem(MEM_FREE, spec2);
    }

  /* Close file */

  fclose(fp);

  fprintf(outp, "OK.\n\n");
}
 
/*****************************************************************************/
