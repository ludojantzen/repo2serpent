/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : adjustsabdata.c                                */
/*                                                                           */
/* Created:       2011/01/23 (JLe)                                           */
/* Last modified: 2016/02/05 (TVi)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: - Adjusts cross sections of nuclides with S(a,b) data        */
/*                                                                           */
/* Comments: - Ton nuklidin uus nimi pitää miettiä (muuta c s:ksi?)          */
/*           - Jos 1002 ja 1004 -reaktioiden Emaxit eroavat, välienergioilla */
/*             voimassaolevan reaktion vaikutusala vähennetään vapaan ytimen */
/*             elastisesta sironnasta.                                       */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AdjustSabData:"

/*****************************************************************************/

void AdjustSabData(long nuc)
{
  long rea, pte, ptr, ne, n, i0, act, ptea, nn, i0a;
  double Emax, Emaxlow, f, xs;

  /* Check S(a,b) flag */

  if (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_SAB_DATA))
    return;  

  /* Skip for SAB nuclides (frac set to zero -> problems)      */
  /* Nailla nuklideilla ei ole myoskaan elastista sirontaa,    */
  /* vaan ainoastaan S(a,b) -reaktiot joten mitaan ei tarvitse */ 
  /* tehda. (nämä siis .xxt -nuklideita) */

  if((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_SAB)
    return;
  
  /* Reset maximum energies */
  /* Emaxlow = maximum energy of the S(a,b) reaction with lower
     maximum energy (maximum energies of reactions differ for some 
     nuclides like H in ZrH) */

  /* HUOM! ZrH-kirjastoilla saattaa näkyä pieni ero MCNP-spektriin 
     9.3 eV tienoilla. Ero syntyy, jos S(a,b) -vaikusualojen summa ei 
     ole täsmälleen sama kuin elastisen sironnan vaikutusala rajaenergioilla, 
     kuten saattaa olla esim. käytettäessä makxsf-interpolointia. MCNP 
     ilmeisesti huolehtii rajasta jollakin toisella tavalla kuin Serpent. */
     
  Emax = 0.0;
  Emaxlow = 0.0;

  /* Loop over reactions */

  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
  while (rea > VALID_PTR)
    {

      /* Check mt */

      if (((long)RDB[rea + REACTION_MT] == 1002) ||
          ((long)RDB[rea + REACTION_MT] == 1004))
        {

          /* Check if EMAX is different for mt = 1002 and mt = 1004 */
          /* and set both Emax and Emaxlow according to the limits */
          
          if(Emax > 0.0 && RDB[rea + REACTION_EMAX] != Emax){

            if(RDB[rea + REACTION_EMAX] > Emax){
              Emaxlow = Emax; 
              Emax = RDB[rea + REACTION_EMAX];
            }
            else{
              Emaxlow = RDB[rea + REACTION_EMAX];
              Emax = Emax;
            }
          }

          /* Compare maximum energy */

          else if(RDB[rea + REACTION_EMAX] > Emax){            
              Emax = RDB[rea + REACTION_EMAX];
          }

          /* Number of points */
          
          ne = (long)RDB[rea + REACTION_XS_NE];

          /* Pointer to data */

          ptr = (long)RDB[rea + REACTION_PTR_XS];

          /* Get factor */

          f = RDB[rea + REACTION_SAB_FRAC];

          /* Check */

          CheckValue(FUNCTION_NAME, "f", "", f, 0.0, 1.0);
          
          /* Loop over points and adjust */

          for (n = 0; n < ne; n++)
            WDB[ptr + n] = RDB[ptr + n]*f;
        }

      /* Next reaction */

      rea = NextItem(rea);
    }

  /* In case of OTF S(a,b) interpolation, do not set elastic xs */
  /* to zero. --> skip following */
  
  if( (long)RDB[nuc + NUCLIDE_PTR_SAB] > VALID_PTR)
    return;

  /* Pointer to elastic */  
  
  rea = (long)RDB[nuc + NUCLIDE_PTR_ELAXS];
  CheckPointer(FUNCTION_NAME, "(ela)", DATA_ARRAY, rea);
  
  /* Get first energy point and number of points */
  
  i0 = (long)RDB[rea + REACTION_XS_I0];
  ne = (long)RDB[rea + REACTION_XS_NE];
  
  /* Pointer to energy data */
  
  pte = (long)RDB[rea + REACTION_PTR_EGRID];
  pte = (long)RDB[pte + ENERGY_GRID_PTR_DATA];
  
  /* Pointer to data */
    
  ptr = (long)RDB[rea + REACTION_PTR_XS];  
  
  /* If emax is the same for all S(a,b) reactions, 
     simply reset elastic scattering cross section */

  if(Emaxlow == 0.0){
    
    for (n = 0; n < ne; n++)
      if (RDB[pte + i0 + n] < Emax)
        WDB[ptr + n] = 0.0;  
  }

  /* Otherwise, subtract the active reaction xs from 
     elastic scattering between Emaxlow and Emax.  */

  else{
   
    /* Loop over reactions */

    rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
    while (rea > VALID_PTR)
      {
        /* Check mt */
        
        if (((long)RDB[rea + REACTION_MT] == 1002) ||
            ((long)RDB[rea + REACTION_MT] == 1004))
          {
            
            /* Recognize active reaction */

            if (RDB[rea + REACTION_EMAX] == Emax){
              
              /* Set elastic cross section to zero below Emaxlow */

              for (n = 0; n < ne; n++){

                  if (RDB[pte + i0 + n] >= Emaxlow)
                    break; 

                  WDB[ptr + n] = 0.0; 

              }
              
              /* Pointers to active reaction xs and egrid */

              act = (long)RDB[rea + REACTION_PTR_XS]; 
              ptea = (long)RDB[rea + REACTION_PTR_EGRID];
              ptea = (long)RDB[ptea + ENERGY_GRID_PTR_DATA];
              
              /* i0 */

              i0a = (long)RDB[rea + REACTION_XS_I0];

              /* Reset index for energy grid of active reaction */

              nn = 0;

              for(; n < ne; n++){
                
                /* Stop at higher upper limit */
                
                if(RDB[pte + i0 + n] >= Emax)
                  break;                

                /* Find energy from grid of active reaction */
                
                while(RDB[ptea + i0a + nn + 1] <= RDB[pte + i0 + n])
                  nn++;

                /* Interpolation factor */
                
                f = (RDB[pte + i0 + n] - RDB[ptea + i0a + nn])/
                  (RDB[ptea + i0a + nn + 1] - RDB[ptea + i0a + nn]);
               
                /* Interpolate active xs */
        
                xs=(1-f)*RDB[act + nn] + f*RDB[act + nn + 1]; 

                /* Subtract active xs from elastic scattering */
                /* ...unless the active xs exceeds elastic scattering
                   which may happen */ 

                if( xs > RDB[ptr + n] )
                  WDB[ptr + n] = 0.0;
                else
                  WDB[ptr + n] = RDB[ptr + n] - xs;

                CheckValue(FUNCTION_NAME, "elaxs", "", RDB[ptr + n], 0.0, INFTY);

              }
              
            }

          }

        /* Next reaction */

        rea = NextItem(rea);

      }
  }       
}

/*****************************************************************************/
