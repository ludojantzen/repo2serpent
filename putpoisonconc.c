/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : putpoisonconc.c                                */
/*                                                                           */
/* Created:       2012/12/05 (JLe)                                           */
/* Last modified: 2012/12/05 (JLe)                                           */
/* Version:       2.1.10                                                     */
/*                                                                           */
/* Description: Puts poison concenctrations from equilibrium iteration into  */
/*              compositions                                                 */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PutPoisonConc:"

/*****************************************************************************/

void PutPoisonConc()
{
  long mat, mat0, iso, ptr;
  double adens;

  /* Check burnup flag */

  if ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == NO)
    return;

  /* Reset parent concenctrations */
  
  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check iterations */
      
      if ((((long)RDB[mat + MATERIAL_XENON_EQUIL_CALC] == NO) &&
           ((long)RDB[mat + MATERIAL_SAMARIUM_EQUIL_CALC] == NO)) ||
          ((long)RDB[mat + MATERIAL_DIV_TYPE] != MAT_DIV_TYPE_PARENT))
        {
          /* Pointer to next */
          
          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Reset compositions */

      if ((long)RDB[mat + MATERIAL_XENON_EQUIL_CALC] == YES)
        {  
          if ((iso = (long)RDB[mat + MATERIAL_PTR_I135_ISO]) > VALID_PTR)
            WDB[iso + COMPOSITION_ADENS] = 0.0;
          if ((iso = (long)RDB[mat + MATERIAL_PTR_XE135_ISO]) > VALID_PTR)
            WDB[iso + COMPOSITION_ADENS] = 0.0;
        }

      if ((long)RDB[mat + MATERIAL_SAMARIUM_EQUIL_CALC] == YES)
        {  
          if ((iso = (long)RDB[mat + MATERIAL_PTR_PM149_ISO]) > VALID_PTR)
            WDB[iso + COMPOSITION_ADENS] = 0.0;
          if ((iso = (long)RDB[mat + MATERIAL_PTR_SM149_ISO]) > VALID_PTR)
            WDB[iso + COMPOSITION_ADENS] = 0.0;
        }

      /* Next material */
      
      mat = NextItem(mat);
    }

  /* Put poisons */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check iterations */
      
      if ((((long)RDB[mat + MATERIAL_XENON_EQUIL_CALC] == NO) &&
           ((long)RDB[mat + MATERIAL_SAMARIUM_EQUIL_CALC] == NO)) ||
          ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT))
        {
          /* Pointer to next */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Pointer to parent */

      mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT];

      /* Xenon iteration */
          
      if ((long)RDB[mat + MATERIAL_XENON_EQUIL_CALC] == YES)
        {
          /* Get I-135 concentration */

          ptr = (long)RDB[mat + MATERIAL_PTR_I135_CONC];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          adens = Mean(ptr, 0);

          /* Put value */

          iso = (long)RDB[mat + MATERIAL_PTR_I135_ISO];
          CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);
          WDB[iso + COMPOSITION_ADENS] = adens;

          /* Put parent */

          if (mat0 > VALID_PTR)
            {
              iso = (long)RDB[mat0 + MATERIAL_PTR_I135_ISO];
              CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);
              WDB[iso + COMPOSITION_ADENS] = RDB[iso + COMPOSITION_ADENS] +
                adens*RDB[mat + MATERIAL_VOLUME]/RDB[mat0 + MATERIAL_VOLUME];
            }

          /* Get Xe-135 concentration */

          ptr = (long)RDB[mat + MATERIAL_PTR_XE135_CONC];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          adens = Mean(ptr, 0);

          /* Put value */

          iso = (long)RDB[mat + MATERIAL_PTR_XE135_ISO];
          CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);
          WDB[iso + COMPOSITION_ADENS] = adens;

          /* Put parent */

          if (mat0 > VALID_PTR)
            {
              iso = (long)RDB[mat0 + MATERIAL_PTR_XE135_ISO];
              CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);
              WDB[iso + COMPOSITION_ADENS] = RDB[iso + COMPOSITION_ADENS] +
                adens*RDB[mat + MATERIAL_VOLUME]/RDB[mat0 + MATERIAL_VOLUME];
            }
        }

      /* Samarium iteration */
      
      if ((long)RDB[mat + MATERIAL_SAMARIUM_EQUIL_CALC] == YES)
        {
          /* Get Pm-149 concentration */

          ptr = (long)RDB[mat + MATERIAL_PTR_PM149_CONC];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          adens = Mean(ptr, 0);

          /* Put value */

          iso = (long)RDB[mat + MATERIAL_PTR_PM149_ISO];
          CheckPointer(FUNCTION_NAME, "(iso1)", DATA_ARRAY, iso);
          WDB[iso + COMPOSITION_ADENS] = adens;

          /* Put parent */

          if (mat0 > VALID_PTR)
            {
              iso = (long)RDB[mat0 + MATERIAL_PTR_PM149_ISO];
              CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);
              WDB[iso + COMPOSITION_ADENS] = RDB[iso + COMPOSITION_ADENS] +
                adens*RDB[mat + MATERIAL_VOLUME]/RDB[mat0 + MATERIAL_VOLUME];
            }

          /* Get Sm-149 concentration */

          ptr = (long)RDB[mat + MATERIAL_PTR_SM149_CONC];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          adens = Mean(ptr, 0);

          /* Put value */

          iso = (long)RDB[mat + MATERIAL_PTR_SM149_ISO];
          CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);
          WDB[iso + COMPOSITION_ADENS] = adens;

          /* Put parent */

          if (mat0 > VALID_PTR)
            {
              iso = (long)RDB[mat0 + MATERIAL_PTR_SM149_ISO];
              CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);
              WDB[iso + COMPOSITION_ADENS] = RDB[iso + COMPOSITION_ADENS] +
                adens*RDB[mat + MATERIAL_VOLUME]/RDB[mat0 + MATERIAL_VOLUME];
            }
        }

      /* Next material */
      
      mat = NextItem(mat);
    }


  /***************************************************************************/
}

/*****************************************************************************/
