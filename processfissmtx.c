/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processfissmtx.c                               */
/*                                                                           */
/* Created:       2012/09/02 (JLe)                                           */
/* Last modified: 2013/03/20 (JLe)                                           */
/* Version:       2.1.13                                                     */
/*                                                                           */
/* Description: Allocates memory for fission matrixes, etc.                  */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessFissMtx:"

/*****************************************************************************/

void ProcessFissMtx()
{
  long loc0, mat, mat0, uni, uni0, ptr, n, type;

  /* Loop over materials and reset indexes */
          
  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Reset index */
      
      WDB[mat + MATERIAL_FMTX_IDX] = -1.0;
      
      /* Next material */
      
      mat = NextItem(mat);
    }

  /* Loop over universes and reset indexes */
          
  uni = (long)RDB[DATA_PTR_U0];
  while (uni > VALID_PTR)
    {
      /* Reset index */
      
      WDB[uni + UNIVERSE_FMTX_IDX] = -1.0;
      
      /* Next material */
      
      uni = NextItem(uni);
    }

  /* Get pointer */

  if ((loc0 = (long)RDB[DATA_PTR_FMTX]) < VALID_PTR)
    return;
  
  /* Reset index */

  n = 0;
  
  /* Check type */
  
  if ((type = (long)RDB[DATA_FMTX_TYPE]) == FISSION_MATRIX_TYPE_MAT)
    {          
      /************************************************************************/
          
      /***** Material-based ***************************************************/
          
      /* Check if material list is given */
      
      if ((mat0 = (long)RDB[loc0 + FMTX_PTR_MAT]) > VALID_PTR)
        {
          while((long)RDB[mat0] > VALID_PTR)
            {
              /* Find match */
              
              mat = (long)RDB[DATA_PTR_M0];
              while (mat > VALID_PTR)
                {
                  /* Compare names */
                  
                  if (CompareStr(mat0, mat + MATERIAL_PTR_NAME))
                    {
                      /* Set index */

                      WDB[mat + MATERIAL_FMTX_IDX] = (double)(n++);
                      
                      /* Break loop */
                      
                      break;
                    }
                  
                  /* Next material */
                  
                  mat = NextItem(mat);
                }
              
              /* Check match */
              
              if (mat < VALID_PTR)
                Error(0, "Fission matrix material %s not defined",
                      GetText(mat0));
              
              /* Next */
              
              mat0++;
            }
        }
      else
        {
          /* Loop over all materials and set indexes and pointers */
          
          mat = (long)RDB[DATA_PTR_M0];
          while (mat > VALID_PTR)
            {
              /* Check fissile flag */
              
              if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_FISSILE_MAT)
                {
                  /* Set index */

                  WDB[mat + MATERIAL_FMTX_IDX] = (double)(n++);
                }
              
              /* Next material */
              
              mat = NextItem(mat);
            }
        }
    }
  else if (type == FISSION_MATRIX_TYPE_UNI)
    {          
      /************************************************************************/

      /***** Universe-based ***************************************************/

      /* Check if universe list is given */
      
      if ((uni0 = (long)RDB[loc0 + FMTX_PTR_UNI]) > VALID_PTR)
        {
          while((long)RDB[uni0] > VALID_PTR)
            {
              /* Find match */
              
              uni = (long)RDB[DATA_PTR_U0];
              while (uni > VALID_PTR)
                {
                  /* Compare names */
                  
                  if (CompareStr(uni0, uni + UNIVERSE_PTR_NAME))
                    {
                      /* Set index */

                      WDB[uni + UNIVERSE_FMTX_IDX] = (double)(n++);
                      
                      /* Break loop */
                      
                      break;
                    }
                  
                  /* Next universe */
                  
                  uni = NextItem(uni);
                }
              
              /* Check match */
              
              if (uni < VALID_PTR)
                Error(0, "Fission matrix universe %s not defined",
                      GetText(uni0));
              
              /* Next */
              
              uni0++;
            }
        }
      else
        {
          /* Loop over all universes and set indexes and pointers */
          
          uni = (long)RDB[DATA_PTR_U0];
          while (uni > VALID_PTR)
            {
              /* Set index */
              
              WDB[uni + UNIVERSE_FMTX_IDX] = (double)(n++);
              
              /* Next universe */
              
              uni = NextItem(uni);
            }
        }
      
      /************************************************************************/
    }
  else if (type == FISSION_MATRIX_TYPE_LVL)
    {          
      /************************************************************************/
          
      /***** Level-based ******************************************************/
          
      /* Loop over all universes and set indexes and pointers */
      
      uni = (long)RDB[DATA_PTR_U0];
      while (uni > VALID_PTR)
        {
          /* Check level */
          
          if ((long)RDB[loc0 + FMTX_LVL] == (long)RDB[uni + UNIVERSE_LEVEL])
            {
              /* Set index */
              
              WDB[uni + UNIVERSE_FMTX_IDX] = (double)(n++);
            }
          
          /* Next universe */
          
          uni = NextItem(uni);
        }
      
      /* Check count */
      
      if (n == 0)
        Error(0, "No universes in fission matrix level %ld",
              (long)RDB[loc0 + FMTX_LVL]);
      
      /* Set type (just to be safe) */
      
      WDB[DATA_FMTX_TYPE] = (double)FISSION_MATRIX_TYPE_UNI;
      
      /************************************************************************/
    }
  else if (type == FISSION_MATRIX_TYPE_XYZ)
    {          
      /************************************************************************/
          
      /***** Mesh based *******************************************************/

      /* Pointer to mesh */

      ptr = (long)RDB[loc0 + FMTX_PTR_MESH];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get mesh size */

      n = (long)(RDB[ptr + MESH_N0]*RDB[ptr + MESH_N1]*RDB[ptr + MESH_N2]);
      
      /************************************************************************/
    }
  else
    Die(FUNCTION_NAME, "Invalid matrix type");
  
  /* Allocate memory for stats */
  
  ptr = NewStat("FMTX", 3, 3, n, n); 
  WDB[loc0 + FMTX_PTR_MTX] = (double)ptr;
  
  ptr = NewStat("FMTX_SRC", 2, 3, n); 
  WDB[loc0 + FMTX_PTR_SRC] = (double)ptr;
  
  /* Put size */
  
  WDB[loc0 + FMTX_SIZE] = (double)n;
}

/*****************************************************************************/
