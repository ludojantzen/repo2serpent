/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printcompositions.c                            */
/*                                                                           */
/* Created:       2010/11/21 (JLe)                                           */
/* Last modified: 2015/08/27 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Prints material compositions after each burnup step          */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrintCompostions:"

/*****************************************************************************/

void PrintCompositions(long step)
{
  long mat, iso, nuc;
  FILE *fp;
  char tmpstr[MAX_STR];

  /* Check mode */

  if (((long)RDB[DATA_BURN_PRINT_COMP] == NO) || (mpiid != 0))
    return;

  /* Check step type */

  if ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP)
    return;

  /* File name */
  
  sprintf(tmpstr,"%s.bumat%ld", GetText(DATA_PTR_INPUT_FNAME), step);
  
  /* Open file */
  
  if ((fp = fopen(tmpstr, "w")) == NULL)      
    Die(FUNCTION_NAME, "Unable to open file for writing");
  
  /* Print comment */
  
  fprintf(fp,"\n%% Material compositions (%1.2f MWd/kgU / %1.2f days)\n\n",
          RDB[DATA_BURN_CUM_BURNUP], RDB[DATA_BURN_CUM_BURNTIME]/24.0/3600.0);
  
  if (RDB[DATA_BURN_PRINT_COMP_LIM] < 1.0)
    fprintf(fp, "%% Decay isotopes with afrac >= %1.5E included\n\n",
            RDB[DATA_BURN_PRINT_COMP_LIM]);
  else
    fprintf(fp, "%% Decay isotopes not included\n\n");

  fprinf(fp,"\n%% Initial fissile mass: %1.20E kg\n\n",
            RDB[DATA_INI_BURN_FMASS]);

  /* Loop over materials */
  
  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check burn flag */
      
      if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
        {
          /* Print material name and density */
          
          fprintf(fp, "mat  %s  %1.14E", GetText(mat + MATERIAL_PTR_NAME),
                  RDB[mat + MATERIAL_ADENS]);
          
          /* Print default id */

          if ((long)RDB[mat + MATERIAL_DEFAULT_PTR_LIB_ID] > VALID_PTR)
            fprintf(fp, " fix %s %1.5E", 
                    GetText(mat + MATERIAL_DEFAULT_PTR_LIB_ID),
                    RDB[mat + MATERIAL_DEFAULT_TMP]);

          /* Print volume */

          fprintf(fp, " vol %1.5E", RDB[mat + MATERIAL_VOLUME]);

          /* Newline */

          fprintf(fp, "\n");

          /* Loop over isotopes */
          
          iso = (long)RDB[mat + MATERIAL_PTR_COMP];
          while (iso > VALID_PTR)
            {
              /* Pointer to nuclide data */
              
              nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
              
              sprintf(tmpstr, "%s", GetText(nuc + NUCLIDE_PTR_NAME));

              if (tmpstr[strlen(tmpstr) - 1] == 's')
                tmpstr[strlen(tmpstr) - 1] = 'c';

              /* Check limit and type and print */
              
              if (((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_TRANSPORT) ||
                  (((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DECAY) &&
                    (RDB[iso + COMPOSITION_ADENS]/RDB[mat + MATERIAL_ADENS] >= 
                     RDB[DATA_BURN_PRINT_COMP_LIM])))
                fprintf(fp, "%20s  %1.14E\n", tmpstr,
                        RDB[iso + COMPOSITION_ADENS]);
              
              /* Next isotope */
              
              iso = NextItem(iso);
            }
        }          
      
      /* Next material */
      
      mat = NextItem(mat);
    }
  
  /* Close file */
  
  fclose(fp);
}

/*****************************************************************************/
