/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processxsdata.c                                */
/*                                                                           */
/* Created:       2010/12/13 (JLe)                                           */
/* Last modified: 2018/11/02 (RTu)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Processes cross sections and ENDF reaction laws              */
/*                                                                           */
/* Comments: - Unionize ei voi olla sama optio kuin gridin generoinnissa     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessXSData:"

/*****************************************************************************/

void ProcessXSData()
{
  long nuc, ace, rea, ptr, type, NES, L0, L1, I0, n, ptp;
  long recon, ne, erg, i0, mt, pte, dum1, dum2;
  const double *XSS;
  double *xs, *xs0, *E0, *tmp, Emin, Emax, mem;

  /* Check decay only mode */

  if ((long)RDB[DATA_BURN_DECAY_CALC] == YES)
    return;

  fprintf(outp, "Processing cross sections and ENDF reaction laws...\n\n");

  /* Do reaction cut-offs */
  
  ReactionCutoff();

  /* Pre-allocate memory for data */

  AllocMicroXS();

  /* Reconstruction option */

  recon = (long)RDB[DATA_OPTI_RECONSTRUCT_MICROXS];

  /* Loop over nuclides */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR)
    {
      /* Skip nuclide if no ACE data */

      if ((long)RDB[nuc + NUCLIDE_PTR_ACE] < VALID_PTR)
        {
          /* Pointer to next nuclide */
          
          nuc = NextItem(nuc);

          /* Loop */

          continue;          
        }

      /* Calculate allocated memory in bytes */

      CalculateBytes();

      /* Get memory size */

      mem = RDB[DATA_TOTAL_BYTES];

      /***********************************************************************/

      /***** Energy grid *****************************************************/

      /* Check nuclide type and unionization */
      
      if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_PHOTON)
        {
          /* Get pointer to unionized grid */
          
          if ((erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_PGRID]) < VALID_PTR)
            Die(FUNCTION_NAME, "No unionized photon grid");
          
          /* Link to nuclide-wise grid */
          
          WDB[nuc + NUCLIDE_PTR_EGRID] = (double)erg;
        }
      else if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DOSIMETRY)
        {
          /* Loop over reactions */
          
          rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
          while (rea > VALID_PTR)
            {
              /* Check type */
                
              if ((long)RDB[rea + REACTION_TYPE] != REACTION_TYPE_SPECIAL)
                Die (FUNCTION_NAME, "Invalid reaction type: %s %ld %ld\n",
                     GetText(nuc + NUCLIDE_PTR_NAME), 
                     (long)RDB[rea + REACTION_MT], 
                     (long)RDB[rea + REACTION_TYPE]); 
                     
              /* Number of energy points and pointer to array */
          
              NES = (long)RDB[rea + REACTION_XS_NE];
              L0 = (long)RDB[rea + REACTION_PTR_EGRID];

              /* Pointer to ACE data */
              
              ace = (long)RDB[nuc + NUCLIDE_PTR_ACE];

              /* Get pointer to XSS array */
          
              XSS = &ACE[(long)ACE[ace + ACE_PTR_XSS]];

              /* Reset grid and number of points */
          
              E0 = NULL;
              ne = 0;
          
              /* Add points to array */
          
              E0 = AddPts(E0, &ne, &XSS[L0], NES);
              
              /* Add extra points just below and above boundaries (en tiedä */
              /* onko tällä mitään väliä) */
          
              tmp = (double *)Mem(MEM_ALLOC, 2, sizeof(double));
          
              tmp[0] = 0.999999999*RDB[rea + REACTION_EMIN];
              tmp[1] = 1.000000001*RDB[rea + REACTION_EMAX];
              
              E0 = AddPts(E0, &ne, tmp, 2);
              
              Mem(MEM_FREE, tmp);
              
              /* Generate grid */
          
              ptr = MakeEnergyGrid(ne, 0, 0, -1, E0, EG_INTERP_MODE_LIN);
              WDB[rea + REACTION_PTR_EGRID] = (double)ptr;
              
              /* Next reaction */

              rea = NextItem(rea);
            }
        }
      else if (recon == YES)
        {
          /* Get pointer to unionized grid */
          
          if ((erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID]) < VALID_PTR)
            Die(FUNCTION_NAME, "No unionized neutron grid");
          
          /* Link to nuclide-wise grid */
          
          WDB[nuc + NUCLIDE_PTR_EGRID] = (double)erg;
        }
      else 
        {
          /* Reset grid and number of points */
          
          E0 = NULL;
          ne = 0;          

          /* Nuclide energy grid not (yet) defined for S(a,b) nuclides */

          if ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_SAB)
            {            
              /* Number of energy points and pointer to array */
            
              NES = (long)RDB[nuc + NUCLIDE_EGRID_NE];
              L0 = (long)RDB[nuc + NUCLIDE_PTR_EGRID];
              
              /* Pointer to ACE data */
              
              ace = (long)RDB[nuc + NUCLIDE_PTR_ACE];
              
              /* Get pointer to XSS array */
              
              XSS = &ACE[(long)ACE[ace + ACE_PTR_XSS]];
              
              /* Reset grid and number of points */
              
              E0 = NULL;
              ne = 0;
              
              /* Add points to array */
            
              E0 = AddPts(E0, &ne, &XSS[L0], NES);
              
              /* Add extra points just below and above boundaries (en tiedä */
              /* onko tällä mitään väliä) */
            
              tmp = (double *)Mem(MEM_ALLOC, 2, sizeof(double));
              
              tmp[0] = 0.999999999*RDB[nuc + NUCLIDE_EMIN];
              tmp[1] = 1.000000001*RDB[nuc + NUCLIDE_EMAX];
              
              E0 = AddPts(E0, &ne, tmp, 2);
              
              Mem(MEM_FREE, tmp);
            }

          /* Loop over reactions */
          
          rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
          while (rea > VALID_PTR)
            {
              /* Get mt */
              
              mt = (long)RDB[rea + REACTION_MT];
              
              /* Check for S(a,b) data */
              
              if ((mt == 1002) || (mt == 1004))
                {
                  /* Get pointer to ACE data */
                  
                  ptr = (long)RDB[rea + REACTION_PTR_NUCLIDE];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                  
                  ace = (long)RDB[ptr + NUCLIDE_PTR_ACE];
                  
                  /* Get pointer to XSS array */
                  
                  XSS = &ACE[(long)ACE[ace + ACE_PTR_XSS]];
                  
                  /* Number of energy points */
                  
                  NES = (long)RDB[rea + REACTION_XS_NE];
                  
                  /* Pointer to energy grid */
                  
                  L0 = (long)RDB[rea + REACTION_PTR_EGRID];
                  
                  /* Add points to array */
                  
                  E0 = AddPts(E0, &ne, &XSS[L0], NES);
                  
                  /* Add extra points for xs/E format data */

                  if ((long)RDB[rea + REACTION_ITP] == 4)
                    {
                      /* Get minimum and maximum energy */
                      
                      Emin = XSS[L0];
                      Emax = RDB[rea + REACTION_SAB_EMAX];
                      
                      /* Make an array of extra points */
                      
                      tmp = MakeArray(Emin, Emax, 1000, 2);
                      
                      /* Add points */
                      
                      E0 = AddPts(E0, &ne, tmp, 1000);
                      
                      /* Free memory */
                      
                      Mem(MEM_FREE, tmp);
                    }
                  
                  /* Add extra points just below and above S(a,b) boundaries */
                  /* (korjattu 3.4.2014 / 2.1.20 ja 27.4 / 2.1.21) */
                  
                  tmp = (double *)Mem(MEM_ALLOC, 4, sizeof(double));

                  tmp[0] = 0.999999999*RDB[rea + REACTION_EMIN];
                  tmp[1] = 1.000000001*RDB[rea + REACTION_EMIN];
                  tmp[2] = 0.999999999*RDB[rea + REACTION_EMAX];
                  tmp[3] = 1.000000001*RDB[rea + REACTION_EMAX];

                  E0 = AddPts(E0, &ne, tmp, 4);

                  Mem(MEM_FREE, tmp);
                }
            
              /* Next reaction */

              rea = NextItem(rea);
            }          

          /* Generate nuclide-wise grid */
          
          ptr = MakeEnergyGrid(ne, 0, 0, -1, E0, EG_INTERP_MODE_LIN);
          WDB[nuc + NUCLIDE_PTR_EGRID] = (double)ptr;

          /* Free temporary array */

          Mem(MEM_FREE, E0);
        }
      
      /***********************************************************************/

      /***** Process interaction data ****************************************/
      
      /* Check type */

      if ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_DOSIMETRY)
        {
          /* Get pointer to energy grid */
          
          erg = (long)RDB[nuc + NUCLIDE_PTR_EGRID];
          CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);
          
          /* Get pointer to grid array */
          
          pte = (long)RDB[erg + ENERGY_GRID_PTR_DATA];
          CheckPointer(FUNCTION_NAME, "(pte)", DATA_ARRAY, pte);
          
          /* Get number of grid points */
          
          ne = (long)RDB[erg + ENERGY_GRID_NE];
          
          /* Allocate memory for temporary cross section array */
          
          xs = (double *)Mem(MEM_ALLOC, ne, sizeof(double));
          
          /* Get pointer to reaction list */
          
          rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
          
          /* Sort reaction list to get branch reactions at the end */
          
          SortList(rea, REACTION_TYPE, SORT_MODE_ASCEND);
        }
      else
        {
          /* Avoid compiler warning */
          
          xs = NULL;
        }

      /* Check data type */

      switch((long)RDB[nuc + NUCLIDE_TYPE])
        {
        case NUCLIDE_TYPE_SAB:
        case NUCLIDE_TYPE_TRANSPORT:
        case NUCLIDE_TYPE_DECAY:
          {                
            /*****************************************************************/

            /***** Transport data ********************************************/

            /* Check transmu-flag if decay data */

            if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DECAY)
              if (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & 
                    NUCLIDE_FLAG_TRANSMU_DATA))
                Die(FUNCTION_NAME, "Decay nuclide without transmu flag");
                
            /* Pointer to reaction data */
            
            rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
            CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
            
            /* Loop over reactions */
            
            while (rea > VALID_PTR)
              {
                /* Energy cut-off (yht�l�isyysmerkki pit�� olla) */

                if (RDB[rea + REACTION_EMIN] >= RDB[DATA_NEUTRON_EMAX])
                  {
                    /* Copy pointer */

                    ptr = rea;

                    /* Pointer to next reaction */
                    
                    rea = NextItem(rea);

                    /* Remove channel */

                    RemoveItem(ptr);
                    
                    /* Loop */
                    
                    continue;
                  }
                
                /* Get type and mt */
                
                type = (long)RDB[rea + REACTION_TYPE];
                mt = (long)RDB[rea + REACTION_MT];

                /* Skip everything but partial and special reactions */
                
                if ((type != REACTION_TYPE_PARTIAL) &&
                    (type != REACTION_TYPE_SPECIAL))
                  {
                    /* Pointer to next reaction */
                    
                    rea = NextItem(rea);
                    
                    /* Loop */
                    
                    continue;
                  }

                /* Skip mt 4 if not pre-calculated (ei voi k�ytt�� VALID_PTR */
                /* koska toi osoittaa ACE-blokkiin) */

                if ((mt == 4) && (RDB[rea + REACTION_PTR_EGRID] < 0.0))
                  {
                    /* Pointer to next reaction */

                    rea = NextItem(rea);
                    
                    /* Loop */
                    
                    continue;
                  }

                /* Get pointer to ACE data */

                ptr = (long)RDB[rea + REACTION_PTR_NUCLIDE];
                CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                ace = (long)RDB[ptr + NUCLIDE_PTR_ACE];

                /* Get pointer to XSS array */
      
                XSS = &ACE[(long)ACE[ace + ACE_PTR_XSS]];
                
                /* Number of energy points and index to first point */
                
                NES = (long)RDB[rea + REACTION_XS_NE];
                I0 = (long)RDB[rea + REACTION_XS_I0];
                
                /* Pointer to energy grid and data */
                
                L0 = (long)RDB[rea + REACTION_PTR_EGRID];
                L1 = (long)RDB[rea + REACTION_PTR_XS];

                /* Check that data begins with a zero */

                if ((I0 > 0) && (XSS[L1] != 0.0) && (mt < 3000))
                  Warn(FUNCTION_NAME, "No leading zero (%s E0 = %E mt = %ld)",
                      GetText(nuc + NUCLIDE_PTR_NAME), XSS[L0], 
                      (long)RDB[rea + REACTION_MT]);

                /* Check interpolation mode */
                
                if ((long)RDB[rea + REACTION_ITP] == 4)
                  {
                    /* S(a,b) data with xs/E form, allocate memory */

                    E0 = (double *)Mem(MEM_ALLOC, NES + 1, sizeof(double));
                    xs0 = (double *)Mem(MEM_ALLOC, NES + 1, sizeof(double));
                    
                    /* Copy arrays */

                    memcpy(E0, &XSS[L0 + I0], NES*sizeof(double));
                    memcpy(xs0, &XSS[L1], NES*sizeof(double));

                    /* Put extra point */

                    E0[NES] = RDB[rea + REACTION_SAB_EMAX];
                    xs0[NES] = 0.0;

                    /* Interpolate */

                    n = InterpolateData(&WDB[pte], xs, ne, E0, xs0, NES + 1, 
                                        4, &dum1, &dum2, NO);

                    /* Free temporary arrays */

                    Mem(MEM_FREE, E0);
                    Mem(MEM_FREE, xs0);
                  }
                else
                  {
                    /* Linear interpolation */

                    if (mt == 301)
                      n = InterpolateData(&WDB[pte], xs, ne, &XSS[L0 + I0], 
                                          &XSS[L1], NES, 0, &dum1, &dum2, YES);
                    else
                      n = InterpolateData(&WDB[pte], xs, ne, &XSS[L0 + I0], 
                                          &XSS[L1], NES, 0, &dum1, &dum2, NO);              
                  }
                
                /* Check negative points */
                
                if (n > 0)
                  Warn(FUNCTION_NAME, 
                       "%ld negative xs points in data (%s mt %ld)", 
                       n, GetText(nuc + NUCLIDE_PTR_NAME), mt);

                /* Negative KERMA points can not be discarded (energy balance) */ 

                if (mt == 301)
                  {
                    i0 = 0;
                    n = ne;
                  }
                else
                  {
                    /* Find first non-zero point */

                    for (i0 = 0; i0 < ne; i0++)
                      if (xs[i0] > 0.0)
                        break;

                    /* T�� tarkistaa vaan sen ett� toi InterpolateData():n */
                    /* laskema indeksi menee oikein. Jostain syyst� se ei  */
                    /* toimi jos i0 == ne. */

                    if ((i0 != ne) && (i0 != dum1))
                      Die(FUNCTION_NAME, "Mismatch %ld %ld %ld", i0, dum1, ne);

                    /* One step back if not at the beginning */
                
                    if (i0 > 0)
                      i0--;

                    /********* T�� EI V�LTT�M�TT� TOIMI ********************/

                    /* Find last point */

                    n = ne - i0 - 1;
                    while (1 != 2)
                      {
                        /* Check if nonzero */
                    
                        if (xs[i0 + n] > 0.0)
                          break;

                        /* Next */

                        if (--n < 0)
                          break;
                      }

                    /* One step forward if not at the end */
                    /*
                      if (n < ne - i0 - 1)
                    */
                
                    /* NOTE: ton if-lauseen mukana olo pudottaa aina (?)   */
                    /* viimeisen pisteen pois, jolloin vikan ja toka vikan */
                    /* v�liss� interpoloidaan kohti nollaa. Yl�rajan yli   */
                    /* s�mpl�tyt neutronit korjataan tolle alueelle -->    */
                    /* aika paha v��ristym� vuohon. (30.1.2012) */

                    n++;
                  }
                /*******************************************************/
                
                /* Discard data if less than two (non-zero) points */
                
                if (n < 2)
                  {
                    /* Check if removed channel has other pointers */

                    if (rea == (long)RDB[nuc + NUCLIDE_PTR_TOTXS])
                      Die(FUNCTION_NAME, "Total xs removed");
                    else if (rea == (long)RDB[nuc + NUCLIDE_PTR_ELAXS])
                      Die(FUNCTION_NAME, "Elastic xs removed");
                    else if (rea == (long)RDB[nuc + NUCLIDE_PTR_FISSXS])
                      Die(FUNCTION_NAME, "Fission xs removed");
                    else if (rea == (long)RDB[nuc + NUCLIDE_PTR_HEATPRODXS])
                      {
                        /* Reset pointer */

                        WDB[nuc + NUCLIDE_PTR_HEATPRODXS] = -1.0;
                        
                        /* Print warning */

                        Note(0, "%s has no non-zero heat production xs",
                             GetText(nuc + NUCLIDE_PTR_NAME));
                      }
                    else if (rea == (long)RDB[nuc + NUCLIDE_PTR_PHOTPRODXS])
                      {
                        /* Reset pointer */

                        WDB[nuc + NUCLIDE_PTR_PHOTPRODXS] = -1.0;

                        /* Print warning */
                        
                        Note(0, "%s has no non-zero photon production xs",
                             GetText(nuc + NUCLIDE_PTR_NAME));
                      }

                    /* Check that reaction is not first in list (may */
                    /* cause problems with S(a,b) data) */

                    if (PrevItem(rea) < VALID_PTR)
                      Die(FUNCTION_NAME, "Reaction is first in list");

                    /* Remember reaction pointer */
                    
                    ptr = rea;
                    
                    /* Pointer to next reaction */
                    
                    rea = NextItem(rea);
                    
                    /* Remove reaction from list */
                    
                    RemoveItem(ptr);
                    
                    /* Mark reaction removed by setting mt to -1 */
                    
                    WDB[ptr + REACTION_MT] = -1.0;
                    
                    /* Cycle loop */
                    
                    continue;
                  }

                /* Allocate memory */
                
                ptr = ReallocMem(DATA_ARRAY, n);
                
                /* Copy data */
                
                /* ...tai nollaa vaikutusala placeholder-reaktioille */

                if(mt == 2002 || mt == 2004)
                  memset(&WDB[ptr], 0.0, n*sizeof(double));
                else
                  memcpy(&WDB[ptr], &xs[i0], n*sizeof(double));
                
                /* Put pointer */
                
                WDB[rea + REACTION_PTR_XS] = (double)ptr;
                
                /* Update first energy point */
                
                WDB[rea + REACTION_XS_I0] = (double)i0;
                WDB[rea + REACTION_XS_NE] = (double)n;
                
                /* Update minimum and maximum energy */
                
                WDB[rea + REACTION_EMIN] = RDB[pte + i0];
                WDB[rea + REACTION_EMAX] = RDB[pte + i0 + n - 1];
                
                /* Link nuclide-wise energy grid to reactiond data */
                
                WDB[rea + REACTION_PTR_EGRID] = (double)erg;

                /* Add reaction in partial list */

                if ((type == REACTION_TYPE_PARTIAL) && 
                    (((long)RDB[DATA_OPT_IMPL_CAPT] == NO) || 
                     ((long)RDB[rea + REACTION_TY] != 0)))
                  {
                    /* Allocate memory for data */
                    
                    ptr = NewItem(nuc + NUCLIDE_PTR_SAMPLE_REA_LIST, 
                                  RLS_DATA_BLOCK_SIZE);
                    
                    /* Allocate memory for reaction counter */
                    
                    ptp = AllocPrivateData(1, PRIVA_ARRAY);
                    WDB[ptr + RLS_DATA_PTR_COUNT] = (double)ptp;

                    /* Put data */

                    WDB[ptr + RLS_DATA_PTR_NUCLIDE] = (double)nuc;
                    WDB[ptr + RLS_DATA_PTR_REA] = (double)rea;
                    WDB[ptr + RLS_DATA_EMIN] = RDB[rea + REACTION_EMIN];
                    WDB[ptr + RLS_DATA_EMAX] = RDB[rea + REACTION_EMAX];
                    WDB[ptr + RLS_DATA_COMP_IDX] = -1.0;
                  }
                
                /* Link elastic, fission and capture cross sections */

                if (mt == 2)
                  WDB[nuc + NUCLIDE_PTR_ELAXS] = (double)rea;
                else if ((mt == 18) || (mt == 19))
                  WDB[nuc + NUCLIDE_PTR_FISSXS] = (double)rea;
                else if ((mt == 102) && ((long)RDB[rea + REACTION_RFS] == 0))
                  WDB[nuc + NUCLIDE_PTR_NGAMMAXS] = (double)rea;

                /* Check type */

                if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_TRANSPORT ||
                    (long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_SAB)
                  {
                    /* Process angular distribution */
                    
                    ProcessMuDistributions(rea);
                    
                    /* Process fission nubar data */
                    
                    ProcessNubarData(rea);
                    
                    /* Process energy distributions */
                    
                    ProcessEDistributions(rea, -1);
                  }

                /* Allocate memory for previous value */

                AllocValuePair(rea + REACTION_PTR_PREV_XS);

                /* Allocate memory for analog reaction rate tally */

                if ((((long)RDB[DATA_ANA_RR_NCALC] == ARR_MODE_ALL) ||
                     (((long)RDB[DATA_ANA_RR_NCALC] == ARR_MODE_BALA) &&
                      (fabs(RDB[rea + REACTION_TY]) != 1.0))) &&
                    ((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_PARTIAL))
                  {
                    ptr = NewStat("ANALOG_REACTION_RATE", 1, 1);
                    WDB[rea + REACTION_PTR_ANA_RATE] = (double)ptr;
                  }

                /* Next reaction */
                
                rea = NextItem(rea);
              }
            
            /* Break case */
            
            break;
            
            /*****************************************************************/
          }
        case NUCLIDE_TYPE_PHOTON:
          {                
            /***** Photon interaction data ***********************************/
                
            /* Pointer to reaction data */
            
            rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
            CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
            
            /* Loop over reactions */
            
            while (rea > VALID_PTR)
              {
                /* Get type and mt */
                
                type = (long)RDB[rea + REACTION_TYPE];
                mt = (long)RDB[rea + REACTION_MT];

                /* Get pointer to ACE data */

                ptr = (long)RDB[rea + REACTION_PTR_NUCLIDE];
                CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                ace = (long)RDB[ptr + NUCLIDE_PTR_ACE];

                /* Get pointer to XSS array */
      
                XSS = &ACE[(long)ACE[ace + ACE_PTR_XSS]];
                
                /* Number of energy points */
                
                NES = (long)RDB[rea + REACTION_XS_NE];
                
                /* Pointer to energy grid and data */
                
                L0 = (long)RDB[rea + REACTION_PTR_EGRID];
                L1 = (long)RDB[rea + REACTION_PTR_XS];

                /* Interpolation (heating numbers are not log) */

                if (mt == 301)
                  InterpolateData(&RDB[pte], xs, ne, &XSS[L0], 
                                  &XSS[L1], NES, 2, NULL, NULL, YES);
                else
                  InterpolateData(&RDB[pte], xs, ne, &XSS[L0], 
                                  &XSS[L1], NES, 1, NULL, NULL, NO);

                /* Special treatment for pair production */

                if (mt == 516)
                  {
                    /* Reset minimum energy */

                    Emin = -1.0;

                    /* Find threshold energy */

                    for (n = 0; n < NES; n++)
                      if(XSS[L1 + n] != 0.0)
                        break;

                    /* Get minimum energy */

                    if (n < NES)
                      Emin = exp(XSS[L0 + n]);
                    else
                      Die(FUNCTION_NAME, "Error");

                    /* Reset reconstructed cross section */

                    for (n = 0; n < ne; n++)
                      if (RDB[pte + n] < Emin)
                        xs[n] = 0.0;
                  }

                /* Allocate memory */
                
                ptr = ReallocMem(DATA_ARRAY, ne);
                
                /* Copy data */
                
                memcpy(&WDB[ptr], xs, ne*sizeof(double));
                
                /* Put pointer */
                
                WDB[rea + REACTION_PTR_XS] = (double)ptr;
                
                /* Update first energy point */
                
                WDB[rea + REACTION_XS_I0] = 0.0;
                WDB[rea + REACTION_XS_NE] = (double)ne;
                
                /* Update minimum and maximum energy */
                
                WDB[rea + REACTION_EMIN] = RDB[pte];
                WDB[rea + REACTION_EMAX] = RDB[pte + ne - 1];
                
                /* Link nuclide-wise energy grid to reactiond data */
                
                WDB[rea + REACTION_PTR_EGRID] = (double)erg;

                /* Add reaction in partial list */

                ptr = NewItem(nuc + NUCLIDE_PTR_SAMPLE_REA_LIST, 
                              RLS_DATA_BLOCK_SIZE);
                
                /* Allocate memory for reaction counter */
                
                ptp = AllocPrivateData(1, PRIVA_ARRAY);
                WDB[ptr + RLS_DATA_PTR_COUNT] = (double)ptp;

                /* Put data */

                WDB[ptr + RLS_DATA_PTR_NUCLIDE] = (double)nuc;
                WDB[ptr + RLS_DATA_PTR_REA] = (double)rea;
                WDB[ptr + RLS_DATA_EMIN] = RDB[rea + REACTION_EMIN];
                WDB[ptr + RLS_DATA_EMAX] = RDB[rea + REACTION_EMAX];
                WDB[ptr + RLS_DATA_COMP_IDX] = -1.0;

                /* Link direct pointers */

                if (mt == 504)
                  WDB[nuc + NUCLIDE_PTR_PHOTON_INCOHEXS] = (double)rea;
                else if (mt == 502)
                  WDB[nuc + NUCLIDE_PTR_PHOTON_COHEXS] = (double)rea;
                else if (mt == 522)
                  WDB[nuc + NUCLIDE_PTR_PHOTON_PHOTOELXS] = (double)rea;
                else if (mt == 516)
                  WDB[nuc + NUCLIDE_PTR_PHOTON_PAIRPRODXS] = (double)rea;
                else if (mt == 301)
                  WDB[nuc + NUCLIDE_PTR_PHOTON_HEATPRODXS] = (double)rea;
                else
                  Die(FUNCTION_NAME, "Invalid mt %ld", mt);

                /* Allocate memory for previous value */

                AllocValuePair(rea + REACTION_PTR_PREV_XS);

                /* Allocate memory for analog reaction rate tally */

                if ((long)RDB[DATA_ANA_RR_PCALC] != ARR_MODE_NONE) 
                  {
                    ptr = NewStat("ANALOG_REACTION_RATE", 1, 1);
                    WDB[rea + REACTION_PTR_ANA_RATE] = (double)ptr;
                  }

                /* Process reaction data */

                if (mt != 301)
                  ProcessPhotonRea(rea);

                /* Next reaction */
                
                rea = NextItem(rea);
              }
            
            /* Break case */
            
            break;
            
            /*****************************************************************/
          }
        case NUCLIDE_TYPE_DOSIMETRY:
          {                
            /***** Dosimetry data (used with detectors only) *****************/ 

            /* Loop over reactions */
          
            rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
            while (rea > VALID_PTR)
              {
                /* Number of energy points and pointer to array */
                
                NES = (long)RDB[rea + REACTION_XS_NE];
                L0 = (long)RDB[rea + REACTION_PTR_EGRID];

                /* Pointer to ACE data */
              
                ace = (long)RDB[nuc + NUCLIDE_PTR_ACE];
                
                /* Get pointer to XSS array */
                
                XSS = &ACE[(long)ACE[ace + ACE_PTR_XSS]];

                /* Allocate memory (ylim��r�iset pisteet yl�- ja alap��ss�) */
                
                ptr = ReallocMem(DATA_ARRAY, NES + 2);

                /* Pointer to reaction data */
                                
                L1 = (long)RDB[rea + REACTION_PTR_XS];

                /* Copy data */
                
                memcpy(&WDB[ptr + 1], &XSS[L1], NES*sizeof(double));
                
                /* Put zeros */

                WDB[ptr] = 0.0;
                WDB[ptr + NES + 1] = 0.0;

                /* Put pointer */
                
                WDB[rea + REACTION_PTR_XS] = (double)ptr;
                
                /* Update first energy point */
                
                WDB[rea + REACTION_XS_I0] = 0.0;
                WDB[rea + REACTION_XS_NE] = (double)(NES + 2);

                /* Allocate memory for previous value */

                AllocValuePair(rea + REACTION_PTR_PREV_XS);

                /* Next reaction */
                
                rea = NextItem(rea);
              }

            /* Break case */
            
            break;
            
            /*****************************************************************/
          }
        }

      /* Free temporary array */

      if (xs != NULL)
        Mem(MEM_FREE, xs);

      /***********************************************************************/

      /***** Remove branch reactions for which the parent is removed *********/

      /* Loop over reactions */

      rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
      while (rea > VALID_PTR)
        {
          /* Check pointer to parent and parent mt (this was set to -1) */

          if ((ptr = (long)RDB[rea + REACTION_PTR_BRANCH_PARENT]) > VALID_PTR)
            if ((long)RDB[ptr + REACTION_MT] == -1)
              {
                /* Copy pointer */

                ptr = rea;
                
                /* Pointer to next reaction */
                
                rea = NextItem(rea);
                
                /* Remove channel */
                
                RemoveItem(ptr);
                
                /* Loop */
                
                continue;
              }

          /* Next reaction */

          rea = NextItem(rea);
        }      

      /***********************************************************************/

      /***** Process distributions, etc. *************************************/

      /* Check type */
      
      if ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_DOSIMETRY)
        {
          /* Close partial reaction list */

          ptr = (long)RDB[nuc + NUCLIDE_PTR_SAMPLE_REA_LIST];
          CloseList(ptr);
          
          /* Adjust cross sections for S(a,b) data */

          AdjustSabData(nuc);
          
          /* Calculate total cross section */
          
          SumTotXS(nuc);
          
          /* Process total inelastic xs */

          ProcessInelaXS(nuc);

          /* Check type */
          
          if ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_DECAY)
            {  
              /* Process unresolved resonance probability table data */
              
              ProcessUresData(nuc);
              
              /* Calculate ures majorants */
              
              CalculateUresMajorants(nuc);
            }
        }

      /* Process photon production */

      if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_TRANSPORT)
        ProcessPhotonProd(nuc);

      /* Set search indexes for double-indexing */

      ProcessDIX(nuc);
      
      /* Update memory size */
  
      WDB[nuc + NUCLIDE_MEMSIZE] = RDB[nuc + NUCLIDE_MEMSIZE] 
        + RDB[DATA_TOTAL_BYTES] - mem;

      /* Update reaction count */

      ReactionCount();

      /* Print data */
  
      PrintNuclideData(nuc, NO);

      /* Calculate allocated memory in bytes */

      CalculateBytes();

      /* Next nuclide */

      nuc = NextItem(nuc);

      /***********************************************************************/
    }
  
  /* Process DBRC and TMS data */

  ProcessTmpData();

  /* Process coarse multi-group majorants */

  CalculateMGXS();
  
  /* Print summary */
  
  PrintNuclideData(-1, 0);

  /* Check that neutron and gamma data exists */

  if (((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES) &&
      ((long)RDB[DATA_N_TRANSPORT_NUCLIDES] < 1))
    Error(0, "No neutron transport data in neutron transport problem");

  if (((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES) &&
      ((long)RDB[DATA_N_PHOTON_NUCLIDES] < 1))
    Error(0, "No photon transport data in photon transport problem");

  /* Free ACE array */

  Mem(MEM_FREE, ACE);

  /* Reset data size */

  WDB[DATA_ALLOC_ACE_SIZE] = 0.0;
  WDB[DATA_REAL_ACE_SIZE] = 0.0;

  /* Calculate allocated memory in bytes */

  CalculateBytes();

  /* Set null pointer */

  ACE = NULL;
}

/*****************************************************************************/
