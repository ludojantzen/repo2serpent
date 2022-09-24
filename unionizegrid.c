/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : unionizegrid.c                                 */
/*                                                                           */
/* Created:       2010/12/10 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Generates unionized energy grid                              */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "UnionizeGrid:"

#define XTRA_GAMMA_PTS 5000

/*****************************************************************************/

void UnionizeGrid()
{
  long nuc, rea, ace, ptr, n, NXS[16], JXS[32], NES, L0, na, ni, np, nr, tot;
  long ng;
  const double *XSS;
  double *Eg, *Ena, *Eni, *tmp, tol, Emin, Emax, mem;

  /* Check decay only mode */

  if ((long)RDB[DATA_BURN_DECAY_CALC] == YES)
    return;

  /* Check unionization */

  if ((long)RDB[DATA_OPTI_UNIONIZE_GRID] == NO)
    {
      fprintf(outp, "Energy grid unionization not performed. Minimum and\n");
      fprintf(outp, "maximum energy set to %1.2E and %1.1f MeV.\n\n",
              RDB[DATA_NEUTRON_EMIN], RDB[DATA_NEUTRON_EMAX]);
            
      /* Exit subroutine */

      return;
    }

  fprintf(outp, "Generating unionize energy grids...\n");
  fprintf(outp, "\nAdding points:\n\n");

  /* Get memory size */

  mem = RDB[DATA_TOTAL_BYTES];

  /* Reset pointers */

  Ena = NULL;
  Eni = NULL;
  Eg = NULL;

  /***************************************************************************/
  
  /***** Process nuclide-wise grids and add points to neutron grid **********/

  /* Reset number of all and important energy points */

  na = 0;
  ni = 0;

  /* Loop over nuclides */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR)
    {
      /* Check type and get pointer to ACE data */

      if (((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_PHOTON) &&
          ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_DECAY) &&
          ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_DOSIMETRY) &&
          ((ace = (long)RDB[nuc + NUCLIDE_PTR_ACE]) > VALID_PTR))
        {
          /* Reset number of added points */
          
          tot = -1;

          /* Read data to NXS array */
          
          ptr = (long)ACE[ace + ACE_PTR_NXS];
          
          for (n = 0; n < 16; n++)
            NXS[n] = (long)ACE[ptr++];
          
          /* Read data to JXS array */
          
          ptr = (long)ACE[ace + ACE_PTR_JXS];
          
          for (n = 0; n < 32; n++)
            JXS[n] = (long)ACE[ptr++];
          
          /* Get pointer to XSS array */
          
          XSS = &ACE[(long)ACE[ace + ACE_PTR_XSS]];

          /* Check type */

          switch((long)RDB[nuc + NUCLIDE_TYPE])
            {
            case NUCLIDE_TYPE_SAB:
            case NUCLIDE_TYPE_TRANSPORT:
              {
                /*************************************************************/

                /***** Transport data ****************************************/

                /* Reset number of added points */

                tot = na;

                /* SAB-nuklideilla on ainoastaan reaktioiden energiagridit */
                /* --> skipataan nämä */

                if ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_SAB)
                  {
                    /* Number of energy points */
                  
                    NES = (long)RDB[nuc + NUCLIDE_EGRID_NE];
                    
                    /* Pointer to energies */
                    
                    L0 = (long)RDB[nuc + NUCLIDE_PTR_EGRID];
                    
                    /* Add points to main array */
                    
                    Ena = AddPts(Ena, &na, &XSS[L0], NES);
                    
                    /* Get important points */
                    
                    tmp = GetImportantPts(ace, NXS, JXS, &n);
                    
                    /* Add points to important array */
                    
                    Eni = AddPts(Eni, &ni, tmp, n);
                    
                    /* Free temporary array */
                    
                    Mem(MEM_FREE, tmp);
                    
                    /* Add extra points just below and above boundaries  */
                    
                    tmp = (double *)Mem(MEM_ALLOC, 2, sizeof(double));
                    
                    tmp[0] = 0.999999999*RDB[nuc + NUCLIDE_EMIN];
                    tmp[1] = 1.000000001*RDB[nuc + NUCLIDE_EMAX];
                    
                    Ena = AddPts(Ena, &na, tmp, 2);
                    Eni = AddPts(Eni, &ni, tmp, 2);
                    
                    Mem(MEM_FREE, tmp);
                  }
                
                /* Check ures flag */
                
                if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] &
                    NUCLIDE_FLAG_URES_USED)
                  {
                    /* Pointer to ures energy grid */
                
                    if ((L0 = JXS[22] - 1) < 0)
                      Die(FUNCTION_NAME, "No pointer to ures data");

                    /* Number of energy points in ures grid */
              
                    NES = (long)XSS[L0];

                    /* Pointer to values */

                    L0 = L0 + 6;

                    /* Allocate memory for temporary array */

                    tmp = (double *)Mem(MEM_ALLOC, NES, sizeof(double));

                    /* Read absolute value to array (NOTE: PURR sets     */
                    /* negative energies for materials consisting of     */
                    /* multiple isotopes with overlapping ures regions.) */

                    for (n = 0; n < NES; n++)
                      tmp[n] = fabs(XSS[L0 + n]);

                    /* Add points to main and important array */

                    Ena = AddPts(Ena, &na, tmp, NES);
                    Eni = AddPts(Eni, &ni, tmp, NES);
 
                    /* Free memory */

                    Mem(MEM_FREE, tmp);
                  }

                /* Check S(a,b) flag */

                if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & 
                    NUCLIDE_FLAG_SAB_DATA)
                  {
                    /* S(a,b) data for standard treatment or interpolation */
                    /* in pre-processing phase */
                    
                    /* Pointer to reaction data */

                    rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
                    CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

                    /* Avoid compiler warning */

                    Emax = -1.0;

                    /* Loop over reactions */

                    while (rea > VALID_PTR)
                      {
                        /* Check mt */

                        if ((long)RDB[rea + REACTION_MT] == 1004)
                          {
                            /* Inelastic channel, set reaction counter */

                            nr = 0;
                
                            /* Reset minimum and maximum energy */
                            
                            Emin = -1.0;
                            Emax = -1.0;
                          }
                        else if ((long)RDB[rea + REACTION_MT] == 1002)
                          {
                            /* Eelastic channel, set reaction counter */
                            
                            nr = 1;
                          }
                        else
                          {
                            /* Not an S(a,b) reaction, skip it */
                            
                            rea = NextItem(rea);
                            continue;
                          }
                        
                        /* Get pointer to nuclide */

                        ptr = (long)RDB[rea + REACTION_PTR_NUCLIDE];
                        CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                        /* Pointer to ACE data */

                        ace = (long)RDB[ptr + NUCLIDE_PTR_ACE];

                        /* Read data to NXS array */
          
                        ptr = (long)ACE[ace + ACE_PTR_NXS];
                        
                        for (n = 0; n < 16; n++)
                          NXS[n] = (long)ACE[ptr++];
                        
                        /* Read data to JXS array */
                        
                        ptr = (long)ACE[ace + ACE_PTR_JXS];
                        
                        for (n = 0; n < 32; n++)
                          JXS[n] = (long)ACE[ptr++];
                        
                        /* Get pointer to XSS array */
                        
                        XSS = &ACE[(long)ACE[ace + ACE_PTR_XSS]];

                        /* Get number of grid points */

                        NES = (long)RDB[rea + REACTION_XS_NE];

                        /* Get pointer to grid */
                        
                        L0 = (long)RDB[rea + REACTION_PTR_EGRID];

                        /* Add points to main and important array */
                    
                        Ena = AddPts(Ena, &na, &XSS[L0], NES);
                        Eni = AddPts(Eni, &ni, &XSS[L0], NES);

                        /* Add extra points for xs/E format data */

                        if ((long)RDB[rea + REACTION_ITP] == 4)
                          {
                            /* Get minimum and maximum energy */

                            Emin = XSS[L0];
                            Emax = RDB[rea + REACTION_SAB_EMAX];

                            /* Make an array of extra points */
                            
                            tmp = MakeArray(Emin, Emax, 1000, 2);

                            /* Add points */

                            Ena = AddPts(Ena, &na, tmp, 1000);
                            Eni = AddPts(Eni, &ni, tmp, 1000);
                        
                            /* Free memory */

                            Mem(MEM_FREE, tmp);
                          }

                        /* Add extra points just below and above boundaries  */

                        tmp = (double *)Mem(MEM_ALLOC, 4, sizeof(double));

                        tmp[0] = 0.999999999*RDB[rea + REACTION_EMIN];
                        tmp[1] = 1.000000001*RDB[rea + REACTION_EMIN];
                        tmp[2] = 0.999999999*RDB[rea + REACTION_EMAX];
                        tmp[3] = 1.000000001*RDB[rea + REACTION_EMAX];
                        
                        Ena = AddPts(Ena, &na, tmp, 4);
                        Eni = AddPts(Eni, &ni, tmp, 4);
                        
                        Mem(MEM_FREE, tmp);

                        /* Update counter */

                        nr++;

                        /* Next reaction */
                    
                        rea = NextItem(rea);
                      }
                  }

                /* Update number of added points */
                
                tot = na - tot;
                
                /* Break case */
                
                break;
                
                /*************************************************************/
              }
            }
          
          fprintf(outp, "%10s -- Points added in neutron grid: %ld\n", 
                  GetText(nuc + NUCLIDE_PTR_NAME), tot);
        }
      
      /* Next nuclide */

      nuc = NextItem(nuc);
    }

  /***************************************************************************/

  /***** Process nuclide-wise grids and add points to photon grid ************/

  /* Reset number of energy points */

  ng = 0;

  /* Loop over nuclides */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR)
    {
      /* Check type and get pointer to ACE data */

      if (((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_PHOTON) &&
          ((ace = (long)RDB[nuc + NUCLIDE_PTR_ACE]) > VALID_PTR))
        {
          /* Get pointer to XSS array */
          
          XSS = &ACE[(long)ACE[ace + ACE_PTR_XSS]];

          /* Reset number of added points */
          
          tot = ng;

          /* Number of energy points */
              
          NES = (long)RDB[nuc + NUCLIDE_EGRID_NE];

          /* Pointer to energies */
                
          L0 = (long)RDB[nuc + NUCLIDE_PTR_EGRID];

          /* Add points to array */
                
          Eg = AddPts(Eg, &ng, &XSS[L0], NES);

          /* Update number of added points */
                
          tot = ng - tot;
          
          /* Print */

          fprintf(outp, "%10s -- Points added in photon grid: %ld\n", 
                  GetText(nuc + NUCLIDE_PTR_NAME), tot);
        }
      
      /* Next nuclide */

      nuc = NextItem(nuc);
    }

  /* Check if points were added and add bacground to enable linear */
  /* interpolation */

  if (ng > 0)
    {
      /* Make an array of exponentially distributed points */
  
      tmp = MakeArray(RDB[DATA_PHOTON_EMIN], RDB[DATA_PHOTON_EMAX],
                      XTRA_GAMMA_PTS, 2);
                  
      /* Convert to log */
  
      for (n = 0; n < XTRA_GAMMA_PTS; n++)
        tmp[n] = log(tmp[n]);
      
      /* Add to grid */
      
      Eg = AddPts(Eg, &ng, tmp, XTRA_GAMMA_PTS);

      /* Free memory */
      
      Mem(MEM_FREE, tmp);
      
      /* Convert the entire grid to linear */

      for (n = 0; n < ng; n++)
        Eg[n] = exp(Eg[n]);
    }

  /* Newline */

  fprintf(outp, "\n");

  /***************************************************************************/

  /***** Add microgroup energy grid ******************************************/

  /* Check if grid is defined */

  if ((ptr = (long)RDB[DATA_MICRO_PTR_EGRID]) > VALID_PTR)
    {
      /* Number of points and pointer to grid */
      
      np = (long)RDB[ptr + ENERGY_GRID_NE];

      ptr = (long)RDB[ptr + ENERGY_GRID_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Add points */

      Ena = AddPts(Ena, &na, &RDB[ptr], np);
      Eni = AddPts(Eni, &ni, &RDB[ptr], np);
    }

 /***************************************************************************/

  /***** Generate coarse multi-group grid ************************************/
 
  /* Check mode */
  
  if ((long)RDB[DATA_OPTI_MG_MODE] == YES)
    { 
      /* Get number of groups */

      np = (long)RDB[DATA_COARSE_MG_NE];
      CheckValue(FUNCTION_NAME, "np", "", np, 10, 50000);
      
      /* Get minimum and maximum values */
      
      Emin = 0.999999*RDB[DATA_NEUTRON_XS_EMIN];
      Emax = 1.000001*RDB[DATA_NEUTRON_XS_EMAX];
      
      /* Generate array with uniform lethargy intervals */
      
      tmp = MakeArray(Emin, Emax, np, 2);
      
      /* Make grid */
      
      ptr = MakeEnergyGrid(np, 0, 0, -1, tmp, EG_INTERP_MODE_LOG);
      
      /* Put pointer */
      
      WDB[DATA_COARSE_MG_PTR_GRID] = (double)ptr;

     /* Add points */

      Ena = AddPts(Ena, &na, tmp, np);
      Eni = AddPts(Eni, &ni, tmp, np);

      /* Free array */

      Mem(MEM_FREE, tmp);
     }

  /***************************************************************************/

  /***** Perform grid thinning ***********************************************/

  /* Check for no grid */

  if (na + ng < 1)
    Die(FUNCTION_NAME, "No grid generated");

  /* Check neutron transport mode */

  if ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES)
    {
      /* Check number of grid points */
      
      CheckValue(FUNCTION_NAME, "na", "", na, 2, 10000000);
        
      /* Get thinning tolerance */

      tol = RDB[DATA_ERG_TOL];

      /* Print */
      
      if (tol > 0.0)
        {
          fprintf(outp, "Performing grid thinning for neutron grid:\n\n");
          fprintf(outp, " - %ld total grid points in unionized grid, with\n", 
                  na);
          fprintf(outp, "   %ld local minima and maxima.\n\n", ni);
        }
      else
        fprintf(outp, "Generating unionized energy grid:\n\n");
  
      /* Store initial size */
  
      WDB[DATA_ERG_INITIAL_PTS] = (double)na; 
  
      /* Perform thinning */

      if (tol > 0.0)
        Ena = ThinGrid(Ena, &na, tol);
      
      /* Store number of important points */
      
      WDB[DATA_ERG_IMPORTANT_PTS] = (double)ni; 
      
      /* Add important points */
      
      if (tol > 0.0)
        Ena = AddPts(Ena, &na, Eni, ni);
      
      /* Re-thin grid with zero tolerance */
      
      if (tol > ZERO)
        Ena = ThinGrid(Ena, &na, ZERO);  
      
      /* Allocate memory for final data */
      
      tmp = (double *)Mem(MEM_ALLOC, na + 2, sizeof(double));

      /* Reset index */
      
      np = 0;

      /* Get minimum and maximum values */

      Emin = RDB[DATA_NEUTRON_EMIN];
      Emax = RDB[DATA_NEUTRON_EMAX];
  
      /* Put lower limit */
      
      tmp[np++] = Emin;
      
      /* Copy values */
      
      for (n = 0; n < na; n++)
        if ((Ena[n] > Emin) && (Ena[n] < Emax))
          tmp[np++] = Ena[n];
      
      /* Put upper limit */
      
      tmp[np++] = Emax;
      
      /* Create unionized neutron grid */
      
      ptr = MakeEnergyGrid(np, 0, 0, -1, tmp, EG_INTERP_MODE_LIN);
      
      /* Put pointer */
      
      WDB[DATA_ERG_PTR_UNIONIZED_NGRID] = (double)ptr;
      
      /* Get final number of grid points */
      
      np = (long)RDB[ptr + ENERGY_GRID_NE];
      
      /* Check ascending order and co-incident energy points */
      
      ptr = (long)RDB[ptr + ENERGY_GRID_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      for (n = 1; n < np; n++)
        if (RDB[ptr + n] <= RDB[ptr + n - 1])
          Die(FUNCTION_NAME, "Error in unionized neutron grid construction");
      
      /* Print */
      
      if (tol > 0.0)
        fprintf(outp, " - Thinning performed with %1.2E fractional tolerance\n",
                tol);
      else
        fprintf(outp, " - Unionization performed without grid thinning\n");
      
      fprintf(outp, "   between %1.2E and %1.1f MeV.\n\n", Emin, Emax);     
      
      fprintf(outp, " - Final neutron grid size: %ld points.\n\n", np);

      /* Free memory */

      Mem(MEM_FREE, tmp);    
    }
  else if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES)
    {
      /* Only photon grid */

      fprintf(outp, "Only photon grid generated:\n\n");
    }

  /* Check photon tranport mode */

  if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES)
    {
      /* Check if photon grid was created */

      if (ng < 2)
        Die(FUNCTION_NAME, "Photon grid not created");

      /* Thinning with zero tolerance */
      
      Eg = ThinGrid(Eg, &ng, ZERO);  

      /* Allocate memory for final data */
      
      tmp = (double *)Mem(MEM_ALLOC, ng + 2, sizeof(double));

      /* Reset index */
      
      np = 0;

      /* Get minimum and maximum values */

      Emin = RDB[DATA_PHOTON_EMIN];
      Emax = RDB[DATA_PHOTON_EMAX];
  
      /* Put lower limit */
      
      tmp[np++] = Emin;
      
      /* Copy values */
      
      for (n = 0; n < ng; n++)
        if ((Eg[n] > Emin) && (Eg[n] < Emax))
          tmp[np++] = Eg[n];

      /* Put upper limit */
      
      tmp[np++] = Emax;


      /* Make grid */

      ptr = MakeEnergyGrid(np, 0, 0, -1, tmp, EG_INTERP_MODE_LOG);
      
      /* Put pointer */
  
      WDB[DATA_ERG_PTR_UNIONIZED_PGRID] = (double)ptr;

      /* Copy number of points */

      ng = np;

      /* Check ascending order and co-incident energy points */
      
      ptr = (long)RDB[ptr + ENERGY_GRID_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      
      for (n = 1; n < ng; n++)
        if (RDB[ptr + n] <= RDB[ptr + n - 1])
          Die(FUNCTION_NAME, "Error in unionized photon grid construction");

      /* Free memory */

      Mem(MEM_FREE, tmp);
    }      
  
  /* Free arrays */

  if (Ena != NULL)
    Mem(MEM_FREE, Ena);

  if (Eni != NULL)
    Mem(MEM_FREE, Eni);

  if (Eg != NULL)
    Mem(MEM_FREE, Eg);

  if (ng > 0)
    fprintf(outp, " - Final photon grid size: %ld points.\n\n", ng);

  if ((mem = (RDB[DATA_TOTAL_BYTES] - mem)/KILO) < KILO)
    fprintf(outp, " - %1.2f kb of memory allocated for grid data\n\n", 
            mem);
  else
    fprintf(outp, " - %1.2f Mb of memory allocated for grid data\n\n", 
            mem/KILO);

  /***************************************************************************/

  fprintf(outp, "OK.\n\n");
}

/*****************************************************************************/
