/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processdecaydata.c                             */
/*                                                                           */
/* Created:       2010/09/11 (JLe)                                           */
/* Last modified: 2018/11/02 (RTu)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Processes data read from ENDF format decay library and       */
/*              and adds reactions in NUCLIDE block                          */
/*                                                                           */
/* Comments: - Käyttää stability cut-offia suoraan lambdaan. Vanhassa        */
/*             serpentissä cut-off on exp(-lamda*t). Korvaa tää uudella      */
/*             hlfcut-parametrilla joka ottaa argumentiksi puoliintumis-     */
/*             ajan (vuosissa tai vapaissa yksiköissä)                       */
/*                                                                           */
/*           - Jatkuvien energiaspektrien log-interpolaatiot konvertoidaan   */
/*             nyt lineaarisiksi koska niiden sämpläystä ei ole vielä        */
/*             implementoitu.                                                */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"
#include "radiotox.h"

#define FUNCTION_NAME "ProcessDecayData:"

/*****************************************************************************/

void ProcessDecayData(long nuc)
{
  long ace, ptr, rea, dec, spec, n, loc0, loc1, loc2, loc3, ne, i, type, part;
  long ne0, m, np, rad;
  double T, norm, *E, *EE, Emin, Emax, chk;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "nuc", DATA_ARRAY, nuc);

  /* Get pointer to ace data */

  if ((ace = (long)RDB[nuc + NUCLIDE_PTR_DECAY_ACE]) < 1)
    Die(FUNCTION_NAME, "Nuclide %s has no decay data", 
        GetText(nuc + NUCLIDE_PTR_NAME));

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "ace", ACE_ARRAY, ace);

  /***************************************************************************/

  /***** Copy and process nuclide data ***************************************/

  /* Copy decay data */

  WDB[nuc + NUCLIDE_LAMBDA] = ACE[ace + ACE_LAMBDA];
  WDB[nuc + NUCLIDE_DECAY_E] = ACE[ace + ACE_DECAY_E];
  WDB[nuc + NUCLIDE_SF_BR] = ACE[ace + ACE_SF_BR];
  WDB[nuc + NUCLIDE_SF_NUBAR] = ACE[ace + ACE_SF_NUBAR];

  /* Set delayed neutron precursor flag */

  if ((long)ACE[ace + ACE_DELNU_PREC] == YES)
    SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_DELNU_PREC);

  /* Get specific radiotoxicities */

  n = 0;
  while ((long)toxdat[n][0] > 0)
    {
      /* Compare */
      
      if ((long)toxdat[n][0] == (long)ACE[ace + ACE_ZAI])
        break;
      else
        n++;
    }

  /* Check */
  
  if ((long)toxdat[n][0] > 0)
    {
      WDB[nuc + NUCLIDE_SPEC_ING_TOX] = toxdat[n][1];
      WDB[nuc + NUCLIDE_SPEC_INH_TOX] = toxdat[n][2];
    }
  else
    {
      WDB[nuc + NUCLIDE_SPEC_ING_TOX] = 0.0;
      WDB[nuc + NUCLIDE_SPEC_INH_TOX] = 0.0;
    }

  /***************************************************************************/

  /***** Copy data to reaction structures ************************************/

  /* Check that nuclide and ACE ZAI match */

  if (RDB[nuc + NUCLIDE_ZAI] != ACE[ace + ACE_ZAI])
    Die(FUNCTION_NAME, "Mismatch in ZAI");

  /* Copy remaining parameters for decay type nuclide */

  if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DECAY)
    {
      WDB[nuc + NUCLIDE_AWR] = ACE[ace + ACE_AWR];
      WDB[nuc + NUCLIDE_AW] = M_NEUTRON*ACE[ace + ACE_AWR];
    }
  else
    {
      /* Compare ACE and ENDF atomic weight ratios */
    }
  
  /* Reset reaction pointer */

  rea = -1;

  /* Loop over decay modes */

  if ((ptr = (long)ACE[ace + ACE_PTR_DECAY_LIST]) > 0)
    while ((dec = (long)ACE[ptr++]) > 0)
      {
        /* Get half-life for mode */

        T = log(2.0)/(ACE[ace + ACE_LAMBDA]*ACE[dec + DECAY_BR] + ZERO);
        
        /* Compare to cut-off */

        if (T < RDB[DATA_DEP_HALF_LIFE_CUTOFF])
          {
            /* Allocate memory for reaction data */

            rea = NewItem(nuc + NUCLIDE_PTR_REA, REACTION_BLOCK_SIZE);

            /* Put nuclide pointer */

            WDB[rea + REACTION_PTR_NUCLIDE] = (double)nuc;

            /* Types (use MT to store first value) */
            
            WDB[rea + REACTION_MT] = ACE[dec + DECAY_RTYP1];
            WDB[rea + REACTION_RTYP2] = ACE[dec + DECAY_RTYP2];
            WDB[rea + REACTION_RTYP3] = ACE[dec + DECAY_RTYP3];
            WDB[rea + REACTION_RTYP4] = ACE[dec + DECAY_RTYP4];
            WDB[rea + REACTION_RTYP5] = ACE[dec + DECAY_RTYP5];
            
            /* Target state (ground or isomeric) */
            
            WDB[rea + REACTION_RFS] = ACE[dec + DECAY_RFS];
            
            /* Q-value */
            
            WDB[rea + REACTION_Q] = ACE[dec + DECAY_Q];
            
            /* Branching */
            
            WDB[rea + REACTION_BR] = ACE[dec + DECAY_BR];
            
            /* Set type */

            WDB[rea + REACTION_TYPE] = (double)REACTION_TYPE_DECAY;
          }
      }

  /* If pointer is < 0, all decay modes are above cut-off */
  
  if (rea < 0)
    return;

  /***************************************************************************/

  /***** Process radiation for decay source **********************************/

  /* Loop over radiation spectra */

  if ((ptr = (long)ACE[ace + ACE_PTR_RAD_SPEC]) > 0)
    while ((spec = (long)ACE[ptr++]) > 0)
      {
        /* Types: 0 -- Gamma rays                                         */
        /*        1 -- Beta minus                                         */
        /*        2 -- Beta plus / EC                                     */
        /*        4 -- Alpha particles                                    */ 
        /*        5 -- Neutrons                                           */
        /*        6 -- Spontaneous fission fragments                      */
        /*        7 -- Protons                                            */
        /*        8 -- Discrete electrons                                 */
        /*        9 -- X-rays and annihilation radiation (photons not     */
        /*             arising as transitions between nuclear states)     */

        /* Get spectrum type */

        type = (long)ACE[spec + RAD_SPEC_TYPE];

        /* Include only certain radiation types for now. */

        if ((type == 0) || (type == 9) || (type == 5) || (type == 1) ||
            (type == 8) || (type == 4))
          {
            /* Avoid compiler warning */

            part = -1;
            
            /* Set particle type */
            
            if ((type == 0) || (type == 9))
              part = PARTICLE_TYPE_GAMMA;
            else if (type == 5)
              part = PARTICLE_TYPE_NEUTRON;
            else if (type == 4)
              part = PARTICLE_TYPE_ALPHA;
            else if ((type == 1) || (type == 8))
              part = PARTICLE_TYPE_ELECTRON;
            else
              Die(FUNCTION_NAME, "Invalid particle type");

            /* Find existing radiation */
            
            rad = (long)RDB[nuc + NUCLIDE_PTR_RADIATIONS];
            while (rad > VALID_PTR)
              {
                /* Check type */
                
                if ((long)RDB[rad + NUCLIDE_RAD_TYPE] == part)
                  break;
                
                /* Next */
                
                rad = NextItem(rad);
            }
            
            /* Check pointer */
            
            if (rad < VALID_PTR)
              {
                /* Create new */

                rad = NewItem(nuc + NUCLIDE_PTR_RADIATIONS, 
                              NUCLIDE_RAD_BLOCK_SIZE);

                /* Put type */
                
                WDB[rad + NUCLIDE_RAD_TYPE] = (double)part;
              }
            
            /* Check number of lines */

            if ((ne = (long)ACE[spec + RAD_SPEC_DISC_NE]) > 0)
              {
                /*************************************************************/

                /***** Discrete lines ****************************************/

                /* Pointer to lines */

                loc0 = (long)ACE[spec + RAD_SPEC_PTR_DISC_E];
                CheckPointer(FUNCTION_NAME, "(loc0)", ACE_ARRAY, loc0);
                
                loc1 = (long)ACE[spec + RAD_SPEC_PTR_DISC_RI];
                CheckPointer(FUNCTION_NAME, "(loc1)", ACE_ARRAY, loc1);
                
                /* Loop over spectral lines */
                
                for (n = 0; n < ne; n++)
                  {
                    /* Add line */
                    
                    loc2 = NewItem(rad + NUCLIDE_RAD_PTR_SPEC, 
                                   DECAY_SPEC_BLOCK_SIZE);

                    /* Set particle type */

                    WDB[loc2 + DECAY_SPEC_PARTICLE] = (double)part; 

                    /* Set emission type */

                    WDB[loc2 + DECAY_SPEC_TYPE] = (double)DECAY_SPEC_LINE;

                    /* Store data */

                    WDB[loc2 + DECAY_SPEC_LINE_E] = ACE[loc0++];
                    WDB[loc2 + DECAY_SPEC_RI] = 
                      ACE[loc1++]*ACE[spec + RAD_SPEC_DISC_NORM];
                    
                    /* Add to specific intensity */
                    
                    WDB[rad + NUCLIDE_RAD_SPEC_I] = 
                      RDB[rad + NUCLIDE_RAD_SPEC_I] +
                      RDB[loc2 + DECAY_SPEC_RI]; 
                  }          
                
                /*************************************************************/
              }
                        
            /* Check  number of points in continuous spectrum */
            
            if ((ne = (long)ACE[spec + RAD_SPEC_CONT_NE]) > 0)
              {
                /*************************************************************/

                /***** Continuous spectrum ***********************************/

                /* Add radiation */

                loc0 = NewItem(rad + NUCLIDE_RAD_PTR_SPEC, 
                               DECAY_SPEC_BLOCK_SIZE);

                /* Set particle type */
                
                WDB[loc0 + DECAY_SPEC_PARTICLE] = (double)part; 

                /* Set emission type */

                WDB[loc0 + DECAY_SPEC_TYPE] = (double)DECAY_SPEC_CONT;

                /* Check number of points */
                
                CheckValue(FUNCTION_NAME, "ne", "", ne, 2, 10000);

                /* Get interpolation scheme */
                
                i = (long)ACE[spec + RAD_SPEC_CONT_INTT];
                CheckValue(FUNCTION_NAME, "ne", "", i, 1, 5);
                
                /* Check interpolation type */

                if ((i == 1) || (i == 2) || (i == 4))
                  {
                    /* Histogram, lin-lin or log-lin, Copy energy grid */
                    
                    loc1 = (long)ACE[spec + RAD_SPEC_PTR_CONT_E];
                    CheckPointer(FUNCTION_NAME, "(loc1)", ACE_ARRAY, loc1);
                    
                    loc3 = ReallocMem(DATA_ARRAY, ne);
                    WDB[loc0 + DECAY_SPEC_CONT_PTR_E] = (double)loc3;
                    
                    for (n = 0; n < ne; n++)
                      {
                        /* Check for coincident points */
                        
                        if ((n > 0) && (ACE[loc1 + n] <= ACE[loc1 + n - 1])) 
                          Die(FUNCTION_NAME, "Coincident points");
                        else
                          WDB[loc3++] = ACE[loc1 + n];
                      }
                    
                    /* Copy PDF */
                    
                    loc2 = (long)ACE[spec + RAD_SPEC_PTR_CONT_PDF];
                    CheckPointer(FUNCTION_NAME, "(loc2)", ACE_ARRAY, loc2);
                    
                    loc3 = ReallocMem(DATA_ARRAY, ne);
                    WDB[loc0 + DECAY_SPEC_CONT_PTR_PDF] = (double)loc3;
                    
                    for (n = 0; n < ne; n++)
                      WDB[loc3++] = ACE[loc2 + n];

                    /* Reset check */

                    chk = -1.0;
                  }
                else
                  {
                    /* Other types, print error to note */

                    Die(FUNCTION_NAME, "Unsupported interpolation type");

                    /* Pointer to energy grid */

                    loc1 = (long)ACE[spec + RAD_SPEC_PTR_CONT_E];
                    CheckPointer(FUNCTION_NAME, "(loc1)", ACE_ARRAY, loc1);
                    
                    /* Reset pointer and number of points */

                    E = NULL;
                    ne0 = 0;

                    /* Loop over intervals */

                    for (m = 1; m < ne; m++)
                      {
                        /* Minimum and maximum energy */

                        Emin = ACE[loc1 + m - 1];
                        Emax = ACE[loc1 + m];
                        
                        /* Set number of points (try different methods) */
                        
                        np = (long)(100.0*log(Emax/Emin) + 50.0);
                        np = 100;

                        /* Make a temporary energy array */

                        EE = MakeArray(Emin, Emax, np, 2);

                        /* Merge with original */

                        E = AddPts(E, &ne0, EE, np);

                        /* Free temporary array */

                        Mem(MEM_FREE, EE);
                      }

                    /* Check number of points */

                    if (ne0 > 10000)
                      Die(FUNCTION_NAME, "ne = %ld ne0 = %ld", ne, ne0);

                    /* Copy energy grid */

                    loc3 = ReallocMem(DATA_ARRAY, ne0);
                    WDB[loc0 + DECAY_SPEC_CONT_PTR_E] = (double)loc3;
                    
                    for (n = 0; n < ne0; n++)
                      {
                        /* Check for coincident points */
                        
                        if ((n > 0) && (E[n] <= E[n - 1])) 
                          Die(FUNCTION_NAME, "Coincident points");
                        else
                          WDB[loc3++] = E[n];
                      }

                    /* Pointer to original PDF */

                    loc2 = (long)ACE[spec + RAD_SPEC_PTR_CONT_PDF];
                    CheckPointer(FUNCTION_NAME, "(loc2)", ACE_ARRAY, loc2);

                    /* Allocate memory */
                    
                    loc3 = ReallocMem(DATA_ARRAY, ne0);
                    WDB[loc0 + DECAY_SPEC_CONT_PTR_PDF] = (double)loc3;

                    /* Reconstruct cross section */

                    InterpolateData(E, &WDB[loc3], ne0, &ACE[loc1], &ACE[loc2],
                                    ne, 6, NULL, NULL, NO);

                    /* Free allocated memory */

                    Mem(MEM_FREE, E);

                    /* Calculate integral for checking */

                    chk = TrapzReal(&ACE[loc1], &ACE[loc2], ne, NULL, i);
                    CheckValue(FUNCTION_NAME, "chm", "", chk, ZERO, INFTY);

                    /* Change interpolation to lin-lin and set number of */
                    /* energy points */

                    i = 2;
                    ne = ne0;
                  }
                
                /* Get pointer to energy array */

                loc1 = (long)RDB[loc0 + DECAY_SPEC_CONT_PTR_E];
                CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

                /* Get pointer to PDF */

                loc2 = (long)RDB[loc0 + DECAY_SPEC_CONT_PTR_PDF];
                CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

                /* Allocate memory for CDF */
                
                loc3 = ReallocMem(DATA_ARRAY, ne);
                WDB[loc0 + DECAY_SPEC_CONT_PTR_CDF] = (double)loc3;
                
                /* Integrate */
                
                norm = TrapzReal(&RDB[loc1], &RDB[loc2], ne, &WDB[loc3], i);
                CheckValue(FUNCTION_NAME, "norm", "", norm, ZERO, INFTY);

                /* Check */

                if ((chk > 0.0) && (fabs(chk/norm - 1.0) > 1E-5))
                  Die(FUNCTION_NAME, "Error in checksum %E (%s)", 
                      chk/norm - 1.0, GetText(nuc + NUCLIDE_PTR_NAME));
                  
                /* Normalize CDF */

                for (n = 0; n < ne; n++)
                  WDB[loc3 + n] = WDB[loc3 + n]/norm;

                /* Normalize pdf */

                loc3 = (long)RDB[loc0 + DECAY_SPEC_CONT_PTR_PDF];
                CheckPointer(FUNCTION_NAME, "(loc3)", DATA_ARRAY, loc3);

                for (n = 0; n < ne; n++)
                  WDB[loc3 + n] = WDB[loc3 + n]/norm;

                /* Store other variables */

                WDB[loc0 + DECAY_SPEC_CONT_NE] = (double)ne;
                WDB[loc0 + DECAY_SPEC_CONT_INTT] = (double)i;
                WDB[loc0 + DECAY_SPEC_RI] = 
                  ACE[spec + RAD_SPEC_CONT_NORM];

                /* Add to specific intensity (NOTE: kun tätä ei kerrota */
                /* tolla norm:illa, niin SF nubar menee suht. oikein. */
                
                WDB[rad + NUCLIDE_RAD_SPEC_I] = RDB[rad + NUCLIDE_RAD_SPEC_I] +
                  RDB[loc0 + DECAY_SPEC_RI]; 

                /* Check if spontaneous fission. NOTE: JEF-2.2 -datassa */
                /* noi nubarit on nollia koska sitä blokkia ei lueta    */
                /* readdecaydata.c:ssä. Toi jakolasku antaa silti ihan  */
                /* oikean näköisiä tuloksia (JLe / 22.7.2017 / 2.1.30). */
                
                if ((part == PARTICLE_TYPE_NEUTRON) && 
                    (RDB[nuc + NUCLIDE_SF_BR] > 0.0))
                  {
                    /* Check nubar */
                    
                    if (RDB[nuc + NUCLIDE_SF_NUBAR] == 0.0)
                      {
                        /* Not defined, calculate */
                        
                        WDB[nuc + NUCLIDE_SF_NUBAR] = 
                          RDB[loc0 + DECAY_SPEC_RI]/RDB[nuc + NUCLIDE_SF_BR];
                      }
                    else if (fabs(RDB[loc0 + DECAY_SPEC_RI]/
                                  RDB[nuc + NUCLIDE_SF_BR] 
                                  - RDB[nuc + NUCLIDE_SF_NUBAR]) > 1E-3)
                      Die(FUNCTION_NAME, "Mismatch in SF nubar");
                  }
                
                /*************************************************************/
              }
          }
      }

  /* Loop over radiations */

  rad = (long)RDB[nuc + NUCLIDE_PTR_RADIATIONS];
  while (rad > VALID_PTR)
    {
      /* Get pointer to spectra */

      if ((ptr = (long)RDB[rad + NUCLIDE_RAD_PTR_SPEC]) < VALID_PTR)
        {
          /* Copy pointer */

          ptr = rad;

          /* Pointer to Next */

          rad = NextItem(rad);

          /* Remove radiation */

          RemoveItem(ptr);
        }
      else
        {
          /* Close list */

          CloseList(ptr);
          
          /* Sort by intensity */
          
          SortList(ptr, DECAY_SPEC_RI, SORT_MODE_DESCEND);
          
          /* Next */
          
          rad = NextItem(rad);
        }
    }

  /***************************************************************************/
}

/*****************************************************************************/
