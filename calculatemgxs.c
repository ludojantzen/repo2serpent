/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : calculatemgxs.c                                */
/*                                                                           */
/* Created:       2011/01/06 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Calculates coarse multi-group majorant cross sections for    */
/*              nuclides                                                     */
/*                                                                           */
/* Comments: - Rutiinia on osittain siivottu 4.11, lähinnä poistettu toinen  */
/*             mg-moodi                                                      */
/*                                                                           */
/*           - NOTE: Rajapisteiden tarkistuksesta korjattiin pieni häikkä    */
/*                   5.11.2011. Piste lasketaan joko tavallisesta tai        */
/*                   ures-majorantti -vaikutusalasta. Alueen viimeisen       */
/*                   tai ensimmäisen pisteen kanssa voi tulla ongelmia,      */
/*                   jotka ehkä korjaantuu sillä että tarkistus tehdään      */
/*                   molemmille eikä joko-tai.                               */
/*                                                                           */
/*           - 7.2.2013 / 2.1.13 eteenpäin lämpötilamajorantti lasketaan     */
/*             nuklidin CE-vaikutusalalle, josta puolestaan lasketaan        */
/*             MG-majorantti.                                                */
/*                                                                           */
/*           - Ures-majorantin laskemista muutettiin 15.2.2016 / 2.1.25      */
/*             siten että nuklidin oma gridi vaihdettiin unionized gridiin.  */
/*             Ei pitäisi vaikuttaa tämän aliohjelman toimintaan sillä       */
/*             optimointimoodissa 2 gridi on jo valmiiksi union gridi.       */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CalculateMGXS:"

/*****************************************************************************/

void CalculateMGXS()
{
  long nuc, rea, loc0, loc1, loc3, ptr, np, i0, ne, n, m0, mp, i;
  double E, *E0, *xs0, val, f, xs1, xs2;
  
  /* Check if data is to be calculated */
  
  if ((long)RDB[DATA_OPTI_MG_MODE] == NO)
    return;

  fprintf(outp, "\nCalculating coarse multi-group cross sections...\n");

  /***************************************************************************/

  /***** Common multi-group energy grid **************************************/

  /* Get number of groups */

  np = (long)RDB[DATA_COARSE_MG_NE];
  CheckValue(FUNCTION_NAME, "np", "", np, 10, 50000);

  /* Pointer to grid */

  ptr = (long)RDB[DATA_COARSE_MG_PTR_GRID];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Pointer to data */

  ptr = (long)RDB[ptr + ENERGY_GRID_PTR_DATA];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Allocate memory */

  E0 = Mem(MEM_ALLOC, np, sizeof(double));

  /* Copy data */

  memcpy(E0, &RDB[ptr], np*sizeof(double));

  /***************************************************************************/

  /***** Calculate cross sections ********************************************/

  /* Allocate memory for temporary data */
  
  xs0 = Mem(MEM_ALLOC, np, sizeof(double));

  /* Loop over nuclides */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR)
    {
      /* Get pointer to total cross section */

      if ((rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS]) < VALID_PTR)
        {
          /* No cross sections, get pointer to next nuclide */

          nuc = NextItem(nuc);

          /* Cycle loop */

          continue;
        }

      /* Get pointer to temperature majorant if given */

      if ((long)RDB[rea + REACTION_PTR_TMP_MAJORANT] > VALID_PTR)
        rea = (long)RDB[rea + REACTION_PTR_TMP_MAJORANT];

      /* Reset temporary array */

      memset(xs0, 0.0, np*sizeof(double));

      /* Array size and index to first point */
      
      ne = (long)RDB[rea + REACTION_XS_NE];
      i0 = (long)RDB[rea + REACTION_XS_I0];
      
      /* Pointer to energy array */
      
      loc0 = (long)RDB[rea + REACTION_PTR_EGRID];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);
      
      /* Check that nuclide uses unionized grid */

      if (loc0 != (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID])
        Die(FUNCTION_NAME, "Mismatch in grid");

      loc0 = (long)RDB[loc0 + ENERGY_GRID_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);
      
      /* Pointer to xs data */
      
      loc1 = (long)RDB[rea + REACTION_PTR_XS];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
      
      /* Check pointer to ures majorant data */
      
      if ((loc3 = (long)RDB[rea + REACTION_PTR_URES_MAX]) > VALID_PTR)
        {
          /* Get first point and number of points */
          
          m0 = (long)RDB[rea + REACTION_URES_MAX_N0];
          mp = (long)RDB[rea + REACTION_URES_MAX_NP];
        }
      else
        {
          /* Reset values */
          
          m0 = 1000000000;
          mp = -1;
        }

      /* Reset index */
      
      n = 0;
      
      /***********************************************************************/
      
      /***** Calculate MG majorants for CE cross sections ********************/

      /* Loop over continuous grid */
      
      for (i = 0; i < ne; i++)
        {
          /* Get energy */
          
          E = RDB[loc0 + i0 + i];
          
          /* Find interval */
          
          while ((E0[n + 1] <= E) && (n < np - 2))
            n++;
          
#ifdef DEBUG
          
          /* Check index */

          ptr = (long)RDB[DATA_COARSE_MG_PTR_GRID];
          if (n != GridSearch(ptr, E))
            Die(FUNCTION_NAME, "Error in coarse grid index search");

#endif
          
          /* Get value (ZERO limit is needed for checking */
          /* the first point) */
          
          if ((val = RDB[loc1 + i]) < ZERO)
            val = ZERO;
          
          /* Adjust ures */
          
          if ((i >= m0) && (i < m0 + mp))
            {
              /* Check pointer */

              CheckPointer(FUNCTION_NAME, "(loc3)", DATA_ARRAY, loc3);

              /* Compare */

              if (RDB[loc3 + i - m0] >= val)
                val = RDB[loc3 + i - m0];     
            }

          /* Compare values */
          
          if (val > xs0[n])
            xs0[n] = val;
        }

      /***********************************************************************/

      /***** Interval boundaries *********************************************/

      /* Allow memory access to avoid error in GridFactor() */
          
      WDB[DATA_PRIVA_MEM_READY] = YES;
      
      /* Loop over coarse grid to get boundary values */
      
      for (n = 0; n < np; n++)
        {
          /* Get grid factor for lower limit */

          ptr = (long)RDB[rea + REACTION_PTR_EGRID];
          f = GridFactor(ptr, E0[n], 0);

          /* Separate integer and decimal parts of interpolation factor */
          
          i = (long)f;
          f = f - (double)i;
          
          /* Get relative index */
          
          i = i - i0;
          
          /* Check boundaries */
          
          if ((i < 0) || (i > ne - 1))
            val = 0.0;
          else
            {      
              /* Get tabulated cross sections */
              
              if ((i + i0 < m0) || (i + i0 >= m0 + mp - 1))              
                {
                  /* Not in ures region, get normal cross sections */
                  
                  xs1 = RDB[loc1 + i];
                  xs2 = RDB[loc1 + i + 1];
                  
                  /* Interpolate */
                  
                  if (i == ne - 1)
                    val = (1.0 - f)*xs1;
                  else
                    val = f*(xs2 - xs1) + xs1;
                }
              else
                {
                  /* Get ures majorant cross sections */
                  
                  xs1 = RDB[loc3 + i0 - m0 + i];
                  xs2 = RDB[loc3 + i0 - m0 + i + 1];
                  
                  /* Interpolate */
                  
                  if (i == mp - 1)
                    val = (1.0 - f)*xs1;
                  else
                    val = f*(xs2 - xs1) + xs1;
                }
            }
                  
          /* Compare */
          
          if (val > xs0[n])
            xs0[n] = val;
          
          /* Check if last point */
          
          if (n == np - 1)
            break;
          
          /* Get grid factor for upper limit */
          
          ptr = (long)RDB[rea + REACTION_PTR_EGRID];
          f = GridFactor(ptr, E0[n + 1], 0);
          
          /* Separate integer and decimal parts of interpolation factor */
          
          i = (long)f;
          f = f - (double)i;
          
          /* Get relative index */
          
          i = i - i0;
          
          /* Check boundaries */
          
          if ((i < 0) || (i > ne - 1))
            val = 0.0;
          else
            {      
              /* Get tabulated cross sections */
              
              if ((i + i0 < m0) || (i + i0 >= m0 + mp - 1))              
                {
                  /* Not in ures region, get normal cross sections */
                  
                  xs1 = RDB[loc1 + i];
                  xs2 = RDB[loc1 + i + 1];
                  
                  /* Interpolate */
                  
                  if (i == ne - 1)
                    val = (1.0 - f)*xs1;
                  else
                    val = f*(xs2 - xs1) + xs1;
                }
              else
                {
                  /* Get ures majorant cross sections */
                  
                  xs1 = RDB[loc3 + i0 - m0 + i];
                  xs2 = RDB[loc3 + i0 - m0 + i + 1];
                  
                  /* Interpolate */
                  
                  if (i == mp - 1)
                    val = (1.0 - f)*xs1;
                  else
                    val = f*(xs2 - xs1) + xs1;
                }
            }
          
          /* Compare */
          
          if (val > xs0[n])
            xs0[n] = val;
        }

      /* Deny memory access */

      WDB[DATA_PRIVA_MEM_READY] = NO;
      
      /***********************************************************************/
      
      /***** Store data ******************************************************/
      
      /* Allocate memory for data */
      
      ptr = ReallocMem(DATA_ARRAY, np);
      
      /* Copy data */
      
      memcpy(&WDB[ptr], xs0, np*sizeof(double));
      
      /* Get pointer to total xs (may be majorant now) */

      rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

      /* Put pointer */
      
      WDB[rea + REACTION_PTR_MGXS] = (double)ptr;
      
      /***********************************************************************/
      
      /* Next nuclide */
      
      nuc = NextItem(nuc);
    }
  
  /***************************************************************************/

  fprintf(outp, "OK.\n");

  /* Free temporary arrays */
  
  Mem(MEM_FREE, E0);
  Mem(MEM_FREE, xs0);
}

/*****************************************************************************/
