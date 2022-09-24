/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : tmpmajorants.c                                 */
/*                                                                           */
/* Created:       2013/02/06 (JLe)                                           */
/* Last modified: 2017/05/04 (JLe)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Calculates / updates majorant cross sections needed for      */
/*              DBRC, TMS and TMS with S(a,b) interpolation.                 */
/*                                                                           */
/* Comments:   - S(a,b) interpoloinnin kanssa majorantti lasketaan           */
/*               vähentämällä normaalista kokonaisvaikutusalaan              */
/*               majorantista elastisen sironnan kontribuutio ja lisäämällä  */
/*               sitten S(a,b)-vaikutusalojen yli laskettu majorantti        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "TmpMajorants:"

/*****************************************************************************/

void TmpMajorants()
{
  long nuc, rea, rea2, ptr, i0, ne, n, i, sab, iso, loc, erg;
  const double *E0, *E1, *xs0, *ela;
  double *xs1, *maj, *min, *max, T0, T1, awr, f, maj0, *tempxs, elaxs, *sabtot, *mino;
  double qparam, sab_Emax, Emaxlow, smax;
  long sidx;
#ifdef TRADMAJO
  double ar;
#else 
  long iter; 
  double v,x, x_old, aux, kT, denom;
#endif
  
  /* Avoid complier warning */ 
  
  mino = NULL; 
    
  fprintf(outp, "\nAdjusting majorant temperatures:\n\n");
  
  /* Loop over nuclides */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR)
    {
      /* Get reaction pointer and temperatures */

      if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DBRC)
        {
          rea = (long)RDB[nuc + NUCLIDE_PTR_ELAXS];
          T0 = RDB[nuc + NUCLIDE_XS_TEMP];
          T1 = RDB[nuc + NUCLIDE_DBRC_MAX_TEMP];

          /* Check TMS flag */

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS)
            Die(FUNCTION_NAME, "Both DBRC and TMS flags set");
        }
      else if (((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS) &&
               ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_DECAY))
        {
          rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
          T0 = RDB[nuc + NUCLIDE_TMS_MIN_TEMP];
          T1 = RDB[nuc + NUCLIDE_TMS_MAX_TEMP];

          /* Check DBRC flag */

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DBRC)
            Die(FUNCTION_NAME, "Both DBRC and TMS flags set");
        }
      else
        T1 = RDB[nuc + NUCLIDE_MAJORANT_TEMP];
      
      /* Check if majorant is already at correct temperature */
      
      if (T1 == RDB[nuc + NUCLIDE_MAJORANT_TEMP]) 
        {
          /* Pointer to next */
          
          nuc = NextItem(nuc);
          
          /* Cycle loop */
          
          continue;
        }
      
      /* Check for zero temperature difference */
      /* This happens if Doppler-preprocessor is able to remove all temperature
         differences, i.e. nuclide minimum temperature = nuclide maximum
         temperature */

      if (fabs(T1 - T0) < ZERO)
        {
          /* Print warning */
          
          if ((long)RDB[nuc + NUCLIDE_PTR_SAB] < VALID_PTR)
            Note(0, 
                 "Nuclide %s is adjusted to a single temperature using TMS.\n Consider using Doppler preprocessor instead.", GetText(nuc + NUCLIDE_PTR_NAME));

          /* Increase majorant temperature artificially by 5 Kelvin to
             force TMS treatment. (temporary solution) */

          WDB[nuc + NUCLIDE_TMS_MAX_TEMP] = T1 + 5.0;
          T1 = RDB[nuc + NUCLIDE_TMS_MAX_TEMP];
        } 

      if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DBRC)
        {
          fprintf(outp, "Nuclide %10s DBRC majorant to %1.0fK...\n", 
                  GetText(nuc + NUCLIDE_PTR_NAME), T1);
          qparam=RDB[DATA_QPARAM_DBRC];
        }
      else
        {
          fprintf(outp, "Nuclide %10s TMS majorant to %1.0fK...\n", 
                  GetText(nuc + NUCLIDE_PTR_NAME), T1);
          qparam=RDB[DATA_QPARAM_TMS];
        }

      /* Check lower temperature */

      if (T0 != RDB[nuc + NUCLIDE_XS_TEMP])
        Die(FUNCTION_NAME, "Incorrect XS temperature");

      /***********************************************************************/

      /***** Get cross sections and energy grid in arrays ********************/

      /* Check pointer */

      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

      /* Get first energy point and number of points */

      i0 = (long)RDB[rea + REACTION_XS_I0];
      ne = (long)RDB[rea + REACTION_XS_NE];

      /* Get cross section array */

      ptr = (long)RDB[rea + REACTION_PTR_XS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      xs0 = &RDB[ptr];

      /* Get pointer to energy grid */

      ptr = (long)RDB[rea + REACTION_PTR_EGRID];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get energy array */

      ptr = (long)RDB[ptr + ENERGY_GRID_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      E0 = &RDB[ptr + i0];

      /* Get majorant reaction */

      ptr = (long)RDB[rea + REACTION_PTR_TMP_MAJORANT];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Pointer to cross section array */

      ptr = (long)RDB[ptr + REACTION_PTR_XS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      xs1 = &WDB[ptr];

      /* Get work arrays for limits amd temporary majorant xs */
    
      min = WorkArray(DATA_PTR_WORK_GRID1, DATA_ARRAY, ne, 0);
      max = WorkArray(DATA_PTR_WORK_GRID2, DATA_ARRAY, ne, 0);
      tempxs = WorkArray(DATA_PTR_WORK_GRID3, DATA_ARRAY, ne, 0);

      /***********************************************************************/

      /***** Equal temperatures **********************************************/

      if (T0 == T1)
        {
          /* Copy data */

          memcpy(xs1, xs0, ne*sizeof(double));

          /* Update majorant temperature */

          WDB[nuc + NUCLIDE_MAJORANT_TEMP] = T1;

          /* Pointer to next nuclide */
          
          nuc = NextItem(nuc);
          
          /* Cycle loop */
          
          continue;
        }

      /***********************************************************************/
      
      /***** Generate limits *************************************************/

      /* Get atomic weight ratio */
      
      awr = RDB[nuc + NUCLIDE_AWR];

#ifdef TRADMAJO

      /************************************************************************/
      
      /***** Traditional Majorant (Dagan, Becker, Rothenstein et al.) *********/

      /* Check that QPARAM is within acceptable limits */

      CheckValue(FUNCTION_NAME, "QPARAM", "", qparam, 2.5, 8.0);

      /* AWR/kT */
      
      ar = awr/(T1 - T0)/KELVIN;          

      /* Loop over energy grid */
      
      for (n = 0; n < ne; n++)
        {
          /* Calculate factor */
          
          f = qparam/sqrt(ar*E0[n]);

          /* Calculate minimum energy */

          if (E0[n] < qparam*qparam/ar)
            min[n] = E0[0];
          else
            min[n] = E0[n]*(1.0 - f)*(1.0 - f);
              
          /* Calculate maximum energy */
          
          max[n] = E0[n]*(1.0 + f)*(1.0 + f);

          /* Check grid boundaries */
          
          if (min[n] < E0[0])
            min[n] = E0[0];
          if (max[n] > E0[ne - 1])
            max[n] = E0[ne - 1];
          
          /* Check */
          
          CheckValue(FUNCTION_NAME, "E", "", E0[n], min[n], max[n]);
        }
      
      /************************************************************************/

#else          

      /************************************************************************/
      
      /***** Revisited limits *************************************************/

      /* Check that QPARAM is within acceptable limits */

      CheckValue(FUNCTION_NAME, "QPARAM", "", qparam, 1E-10, 3.2E-4);

      kT = (T1 - T0)*KELVIN;

      /* Loop over energies */
      
#ifdef OPEN_MP
#pragma omp parallel for private(v, denom, x, iter, x_old, aux)
#endif          
      for (n = 0; n < ne; n++)
        {
          /* Lasketaan rajoja vastaavat suhteelliset nopeudet (kertaa gamma) */
          /* ja niista edelleen vastaavat energiat */
          
          /*** Init some variables */
          
          /* neutron velocity * gamma */
          
          v = sqrt(E0[n]*awr/kT);
          
          /* The same as "PotCorr" or g-factor, but written for neutron speed */
          /* instead of energy */

          denom = ( (1 + 1/(2*v*v)) * erf(v) + exp(-v*v)/(SQRTPI*v) ) 
            * (2*SQRTPI*v*v)*qparam/2;
              
          /************************************/
          /*** Solve lower boundary (E_min) ***/
          /************************************/
          
          /* Initial guess for solution */
          
          x = 1.00;             
                    
          /* Solve x using Newton's method */
          /* ... joka on toivottavasti riittavan hyva tahan tarkoitukseen */
          /* ja toivottavasti suppenee niin kuin pitaakin */ 

          iter = 0;
              
          do
            {
              x_old = x;
                
              /* Nominator : the function */
              
              x = exp(-(v+x_old)*(v+x_old))*(x_old-v)-exp(-(v-x_old)*(v-x_old))*
                (v+x_old)+2*exp(-v*v)*v+( SQRTPI*v*v +SQRTPI/2 ) * 
                ( erfc(v-x_old) + erfc(v+x_old) -2*erfc(v))-denom;
              
              /* Denominator: the derivative */
              
              if( fabs( (aux = 2*x_old*x_old*(exp(-(v-x_old)*(v-x_old)) 
                                              - exp(-(v+x_old)*(v+x_old))))) < 1E-50 || 
                  aux > INFTY){

                /* Nuo eksponentit rajahtaa kasiin jos alkuarvaus on kovin    */
                /* huono, niin kuin tapahtuu isoilla energioilla joilla gridi */
                /* on harva-> jos nain kay, niin arvataan uudestaan */ 
                
                  x = v;
                
                  continue;
              }                                
              
              /* Divide nominator by denominator */
              
              x /= aux;
              
              /* Subtract from previous value */
              
              x = x_old - x;
              
              iter++;
              
              if (iter > 1E5)
                {
                  Warn(FUNCTION_NAME, "Maximum number of iterations exceeded (Emin loop), Emin = %E, E=%E", x*x*kT/awr, E0[n]);
                  break; 
                }
              
              /* Loop until converged (relative difference < 1E-9) */
              
            } 
          while (fabs(x/x_old-1) > 1E-9);

          /* 1E-8 naytti olevan riittavasti, mutta pannaan varmuuden vuoksi */
          /* 1E-9 */
          
          /* To energy units */
          
          min[n] = x*x*kT/awr;

          /************************************/
          /*** Solve upper boundary (E_max) ***/
          /************************************/
                
          /* Initial guess */

          x = 1.0;
          
          iter = 0;
          
          do
            {
              x_old = x;

              /* Nominator : the function */
              
              x = exp(-(v+x_old)*(v+x_old))*(v-x_old)
                + exp(-(v-x_old)*(v-x_old))*(v+x_old)+( SQRTPI*v*v +SQRTPI/2 ) 
                * ( erf(v-x_old) + erf(v+x_old) )-denom;
              
              /* Denominator: the derivative */
              
              if( fabs( (aux= 2*x_old*x_old*( exp(-(v+x_old)*(v+x_old)) 
                                              - exp(-(v-x_old)*(v-x_old))))) < 1E-50 ||
                  aux > INFTY )
                {

                  /* Nuo eksponentit rajahtaa kasiin jos alkuarvaus on kovin */
                  /* huono -> jos nain kay, niin arvataan uudestaan */

                  x=v;
                  
                  continue;                 
                }
              
              /* Divide nominator by denominator (Newton)*/
              
              x /= aux;
              
              /* Subtract from previous value (Newton) */
              
              x = x_old - x;
              
              iter++;
              
              if(iter > 1E5)
                {
                  Warn(FUNCTION_NAME, "Maximum number of iterations exceeded (Emax loop), Emax = %E, E=%E", x*x*kT/awr, E0[n]);
                  break; 
                }
              
              /* Loop until converged */
            
            } 
          while(fabs(x/x_old-1)>1E-9 );
              
          max[n] = x*x*kT/awr;
          
          /*** A few important checks */
          
          /* Nain voi kayda raskaille ytimille */
          
          if(min[n] < E0[0])
            min[n] = E0[0];
          
          if(min[n] > E0[ne -1])
            Die(FUNCTION_NAME, "Emin exceeds the highest point of the energy grid. (min[n] = %E, E = %E, dT = %.1f K)\n", min[n], E0[n], 
                kT/KELVIN);
          
          /* Nain tapahtuu usein */

          if(max[n] > E0[ne -1])
            max[n] = E0[ne -1];

          /* Check */
          /* Vaikka noi tarkistukset nayttaa hassuilta, niissa on jarkea: */
          /* - minimiraja on aina isompi tai yhtasuuri kuin energiagridin */
          /*   alin piste ja lisaksi pienempi kuin maksimi                */
          /* - maksimienergia on aina suurempi kuin neutronin energia ja  */
          /*   pienempi kuin gridin suurin piste                          */
          /* HUOM! Kun rajat lasketaan talla tavalla, E_min ei            */
          /* valttamatta ole pienempi kuin neutronin energia matalilla    */
          /* energioilla */
          
          CheckValue(FUNCTION_NAME, "E_min", "", min[n], E0[0], max[n]);
          CheckValue(FUNCTION_NAME, "E_max", "", max[n], E0[n], E0[ne-1]);
        }

#endif 
      
      /***********************************************************************/
      
      /***** Calculate indexes from limit energies ***************************/
      
      /* Calculate indexes for lower limit */
      
      n = 0;
      for (i = 0; i < ne; i++)
        {
          while(E0[n] <= min[i])
            n++;
          
          /* Check minimum, maximum and interval */

          if (n == 0)
            min[i] = 0.0;
          else if (n > ne - 1)
            Die(FUNCTION_NAME, "Error (n > ne -1) (1) (min: %E)", min[i] );
          else
            {
              /* Check interval */
              
              CheckValue(FUNCTION_NAME, "min", "", min[i], E0[n - 1], E0[n]);
              
              /* Put index */
              
              min[i] = (double)(n - 1) + (min[i]-E0[n-1])/(E0[n]-E0[n-1]);
            }
        }

      /* Calculate indexes for upper */

      n = 0;
      for (i = 0; i < ne; i++)
        {
          while(E0[n] < max[i])
            n++;
          
          /* Check minimum, maximum and interval */

          if (n == 0)
            Die(FUNCTION_NAME, "Error (n == 0) (1)");
          else if (n > ne - 1)
            Die(FUNCTION_NAME, "Error (n > ne -1) (2)");
          else
            {
              /* Check interval */
              
              CheckValue(FUNCTION_NAME, "min", "", max[i], E0[n - 1], E0[n]);

              /* Put index */

              max[i] = (double)(n - 1) + (max[i]-E0[n-1])/(E0[n]-E0[n-1]);
            }
        }

      /* If the upper energy limits (Emax) of S(a,b) reactions of a nuclide 
         differ, it has to be taken into account in majorant generation */
      
      /* Find out Emaxlow */

      /* Initialize */ 
      Emaxlow = 0.0; 
      
      if( (sab = (long)RDB[nuc + NUCLIDE_PTR_SAB]) > VALID_PTR){
                  
        while(sab > VALID_PTR){
          
          /* Pointer to S(a,b) nuclide */
          
          iso = (long)RDB[sab + SAB_PTR_ISO];
          
          /* Loop over reactions */

          rea = (long)RDB[iso + NUCLIDE_PTR_REA];

          while( rea > VALID_PTR){
            
            if(( (long)RDB[rea + REACTION_MT] == 1002 || 
                 (long)RDB[rea + REACTION_MT] == 1004 ) &&
               RDB[rea + REACTION_EMAX] != RDB[nuc + NUCLIDE_SAB_EMAX]){
             
              if(RDB[rea + REACTION_EMAX] < RDB[nuc + NUCLIDE_SAB_EMAX]){
                
                /* Check that Emaxlow is the same for all temperatures */
                if(Emaxlow != 0.0 && Emaxlow != RDB[rea + REACTION_EMAX])
                  Die(FUNCTION_NAME, "Upper E-limits of S(a,b) reactions must be same for different temperatures");                
                else
                  Emaxlow = RDB[rea + REACTION_EMAX];
              }
              else
                Die(FUNCTION_NAME, "EMAX of S(a,b) reaction above NUCLIDE_SAB_EMAX");
            }                          
            rea = NextItem(rea);
          }                                      
          sab = NextItem(sab); 
        }        
      }

      /* Set Emaxlow of nuclide */

      WDB[nuc + NUCLIDE_SAB_EMAXLOW] = Emaxlow;

      /* In OTF S(a,b) mode, get pointer to elastic scattering data 
         for subtraction */

      if (RDB[nuc + NUCLIDE_PTR_SAB] > VALID_PTR) 
        {
          rea2 = (long)RDB[nuc + NUCLIDE_PTR_ELAXS];
          ptr = (long)RDB[rea2 + REACTION_PTR_XS];
          
          CheckPointer(FUNCTION_NAME, "(ptr2)", DATA_ARRAY, ptr);
          ela = &RDB[ptr]; 
        }
      else
        {
          /* This pointer was not set until 4.5.2017 / 2.1.29 (JLe) */

          ela = NULL;        
        }

      /* OpenMP parallel block */      
      /* Loop over energy grid */

#ifdef OPEN_MP
#pragma omp parallel for private(n, i, f, maj0) 
#endif          
        for (n = 0; n < ne; n++)
          {
            maj0 = -1.0;
            
            /* Loop over interval */    

            for (i = (long)min[n] + 1; i <= (long)max[n]; i++)
              {
                if (E0[n] < RDB[nuc + NUCLIDE_SAB_EMAX])
                  {
                    /* Between Emaxlow and Emax, do not subtract elastic */
                    /* scattering from majorant */

                    /* Tästä pitää vähentää vielä aktiivisen reaktion  */
                    /* minimivaikutusala. Tämä tehdään vasta alempana, */
                    /* S(a,b) -majorantin lisäämisen yhteydessä */            
                    
                    if ((Emaxlow != 0.0) && (E0[n] > Emaxlow))
                      { 
                        if (xs0[i] > maj0)
                          maj0 = xs0[i];
                      }
                    else if (ela == NULL)
                      {
                        /* Lisättiin 4.5.2017 / 2.1.29 (Jle) */

                        Die(FUNCTION_NAME, "Null pointer (ela)");
                      }
                    else
                      { 
                        if ((xs0[i] - ela[i]) > maj0)
                          maj0 = xs0[i] - ela[i];
                      }

                  }
                else if (xs0[i] > maj0)
                  maj0 = xs0[i];
              }
            
            /* Check boundary energies by interpolating */
            
            /* Emin */

            i = (long)min[n];
            f = min[n]-(double)i;
            
            /* Get elastic xs for subtraction in the case of OTF S(a,b) */

            if(E0[n] < RDB[nuc + NUCLIDE_SAB_EMAX])
              {
                if (ela == NULL)
                  Die(FUNCTION_NAME, "Null pointer (ela)");
                else
                  elaxs = ela[i] + f*(ela[i + 1] - ela[i]);
              }
            else 
              elaxs = 0.0;

            if (xs0[i] + f*(xs0[i + 1] - xs0[i]) - elaxs > maj0)
              maj0 = xs0[i] + f*(xs0[i + 1] - xs0[i]) - elaxs;
            
            /* Emax */

            i = (long)max[n];
            f = max[n]-(double)i;

            /* Get elastic xs for subtraction in the case of OTF S(a,b) */

            if(E0[n] < RDB[nuc + NUCLIDE_SAB_EMAX])
              {
                if (ela == NULL)
                  Die(FUNCTION_NAME, "Null pointer (ela)");
                else
                  elaxs = ela[i] + f*(ela[i + 1] - ela[i]);
              }
            else 
              elaxs = 0.0;

            if (xs0[i] + f*(xs0[i + 1] - xs0[i]) - elaxs > maj0)
              maj0 = xs0[i] + f*(xs0[i + 1] - xs0[i]) - elaxs;           

            /* T0:n miinustaminen hoidetaan vasta PotCorrissa */

            if (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DBRC))
              maj0 = maj0*PotCorr(nuc, E0[n], T1*KELVIN);
 
            /* Check value */
            
            CheckValue(FUNCTION_NAME, "maj", "", maj0, 0.0, INFTY); 
            
            /* Put value to temporary array*/
            
            tempxs[n] = maj0;
          }
        
        /*********************************************************************/
        
        /********* Generate S(a,b) majorants for OTF *************************/
        
        /* On-the-fly S(a,b) -nuklideille:                                  */
        /* TOTXS pitää sisällään vapaan nuklidin reaktiot ml. vapaan        */
        /* elastisen sironnan. Edellä poistettiin vapaan nuklidin elastisen */ 
        /* osuus tmp-majorantista. Seuraavassa lasketaan S(a,b)-alueella    */
        /* lämpötilamajorantti nuklidiin liitetyistä S(a,b) -nuklideista ja */
        /* lisätään tämä majorantti vanhaan TMP_MAJORANTttiin */
                
        /* Check if OTF sab data exists */
        
        if ((sab = (long)RDB[nuc + NUCLIDE_PTR_SAB]) > VALID_PTR)
          {
            /* Alloc & init temporary array for majorant and minorant */

            maj = (double *)Mem(MEM_ALLOC, ne, sizeof(double));          
            memset(maj, 0.0, ne*sizeof(double));
            
            if(Emaxlow > 0.0)
              {
                mino = (double *)Mem(MEM_ALLOC, ne, sizeof(double));
                memset(mino, 0.0, ne*sizeof(double));
              }
            
            sab_Emax = 0.0;
            
            /* Calculate temperature majorant by looping over S(a,b) data */
            
            while (sab > VALID_PTR)
              {
                iso = (long)RDB[sab + SAB_PTR_ISO];
                loc = (long)RDB[iso + NUCLIDE_PTR_TOTXS];
                CheckPointer(FUNCTION_NAME, "(totxs1)", DATA_ARRAY, loc);
                
                erg = (long)RDB[loc + REACTION_PTR_EGRID];
                E1 = &RDB[(long)RDB[erg + ENERGY_GRID_PTR_DATA]];

                /* Upper limit for S(a,b) treatment */

                sab_Emax = RDB[iso + NUCLIDE_EMAX];

                /* Pointer to S(a,b) xs data */
                
                ptr = (long)RDB[loc + REACTION_PTR_XS];
                sabtot=&WDB[ptr]; 
                
                /* Check if energy grids differ (opti modes 1 and 3) */
                /* If the grids differ, the majorant is determined   */
                /* as the maximum S(a,b) total xs in the vicinity of */
                /* majorant energy grid point E0[i] */

                if ((ne != (long)RDB[loc + REACTION_XS_NE]) ||
                    (i0 != (long)RDB[loc + REACTION_XS_I0]))
                  {
                    for (i = 0; i < ne; i++)
                      {
                        if (E0[i] > sab_Emax)
                          break;

                        /* Find maximum cross section around E0[i] */
             
                        sidx = 0;

                        while (E1[sidx + 1] < E0[i - 1])
                          sidx++;

                        smax = 0.0; 

                        while (E1[sidx] < E0[i + 1])
                          {
                            if (smax < sabtot[sidx]) 
                              smax = sabtot[sidx];
                            sidx++;
                          }

                        /* If the maximum is higher than majorant, */
                        /* update majorant */

                        if (maj[i] <= smax)
                          maj[i] = smax; 
                      }                    
                  }

                /* If S(a,b) reactions and majorant have the same energy   */
                /* grid, the majorant can be calculated simply by checking */
                /* each energy grid point individually */

                else
                  {
                    /* Find maximum total (sum of S(a,b) reactions) cross */
                    /* section */
              
                    for (i = 0; i < ne; i++)
                      {                            
                        if (E0[i] > sab_Emax)
                          break;
                
                        if (maj[i] <= sabtot[i])
                          maj[i] = sabtot[i]; 
                
                        /* Find also minimum for subtraction */ 
                        
                        if (Emaxlow != 0.0)
                          { 
                            if (mino[i] == 0.0)
                              mino[i] = sabtot[i]; 
                            else if (mino[i] > sabtot[i])
                              mino[i] = sabtot[i];
                          }                  
                      }
                  }
                
                sab = NextItem(sab);
              }
            
            /* Add contribution to majorant, calculated for */
            /* total - free_elastic */
            
            for (i = 0; i < ne; i++)
              {            
                /* Subtract (minimum) active reaction from majorant */

                if ((Emaxlow != 0.0) && (E0[i] > Emaxlow))
                  tempxs[i] -= mino[i];

                if (E0[i] > sab_Emax)
                  {
                    /* If energy grids differ, add maj also to the first */
                    /* xs point above boundary */
                 
                    if ((ne != (long)RDB[loc + REACTION_XS_NE]) ||
                        (i0 != (long)RDB[loc + REACTION_XS_I0]))
                      tempxs[i] += maj[i - 1];
                    
                    break;
                  }
            
                tempxs[i] += maj[i];
              }
          
            Mem(MEM_FREE, maj);
          }
        
        /* Koska majoranttivaikutusaloista kaivataan jatkossa kullekin  */
        /* energiavalille ainoastaan kahden valia reunustavan energia-  */   
        /* gridipisteen maksimiarvoa, voidaan nama maksimit laskea jo   */    
        /* etukateen. Talla tavalla materiaalimajorantit tulevat oikein */ 
        /* lasketuksi ja plussana majoranttia muistista haettaessa      */
        /* tarvitsee noutaa vain yksi arvo ilman interpolointia.        */
        /* Majorantti on siis histogrammityyppia. */ 
        
#ifdef OPEN_MP
#pragma omp parallel for private(n)
#endif          
        for (n = 0; n < ne - 1; n++)
          {
            if (tempxs[n] > tempxs[n + 1])
              xs1[n] = tempxs[n];
            else
              xs1[n] = tempxs[n + 1];
          }
        
        /* Last energy grid point */
        
        xs1[ne - 1] = tempxs[ne - 1];
        
        /* Update majorant temperature */
        
        WDB[nuc + NUCLIDE_MAJORANT_TEMP] = T1;
        
        /*********************************************************************/        
        
        /* Next nuclide */
        
        nuc = NextItem(nuc);
    }
}

/*****************************************************************************/
