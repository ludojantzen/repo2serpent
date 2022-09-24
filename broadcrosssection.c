/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : broadcrosssection.c                            */
/*                                                                           */
/* Created:       2011/07/28 (JLe)                                           */
/* Last modified: 2013/01/05 (JLe)                                           */
/* Version:       2.1.13                                                     */
/*                                                                           */
/* Description: Perfrorms integreation for DopplerBroad()                    */
/*                                                                           */
/* Comments: - From Serpent 1.1.15                                           */
/*           - Original routine developed by Tuomas Viitanen                 */
/*           - OpenMP -rinnakkaistus lis‰tty ja muuttujien nimi‰ lyhennetty  */
/*             4.2.2013 / 2.1.13                                             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifdef INTELCC
#include <mathimf.h>
#else
#include <tgmath.h>
#endif

#define FUNCTION_NAME "BroadCrossSection:"

/*****************************************************************************/

void BroadCrossSection(long L0, long L1, long NES, double Ti, double Tf, 
                       double awr, double *XSS, long *lims, long idx, long sz)
{
  long i, i0, mlow, intidx;
  double alpha, aux;
  register double x1, x2, v;
  double *intg, *xs;
  
  /* Varataan muistia */
  
  intg = (double *)Mem(MEM_ALLOC, NES, sizeof(double));
  xs = (double *)Mem(MEM_ALLOC, NES, sizeof(double));
  
  /* Varmuuden vuoksi tarkistetaan ett‰ ei jaeta nollalla ja tuoteta */
  /* siten h‰m‰‰vi‰ tuloksia. T‰t‰ ei pit‰isi ikin‰ tapahtua. En‰‰.  */
  /* Lis‰tty 27.4.2009 (TVi) */
  
  if (Tf == Ti)
    Die(FUNCTION_NAME, "Cross section already in desired temperature");
  
#ifdef DEBUG
  
  /* Tarkistetaan negatiiviset, raja on sama kuin interpolatedata.c:ss‰ */
  /* (JLE 31.8.2011) */
  
  for (i = 0; i < NES; i++)
    if (XSS[L1 + i] < -1.0)
      Warn(FUNCTION_NAME, "Negative cross sections before broadening %E",
           XSS[L1 + i]);
  
#endif
  
  /* Alffa = m/2k(delta T) Kannattaa huomioida, ett‰ SPD_C on      */
  /* senttimetriyksikˆiss‰, mink‰ vuoksi kymppitonnilla kertominen */
  /* osoittajassa */
  
  alpha = (awr*NEUTRON_E0*1E4)/(SPD_C*SPD_C*2.0*KELVIN*(Tf - Ti));
  
  /* OpenMP parallel block */
  
#ifdef OPEN_MP
#pragma omp parallel private(i, i0, v, intidx, x1, x2, aux, mlow)
#endif
  {
    /* Kynnyseneriallisten reaktioiden tapauksessa ensimm‰inen piste pit‰‰  */
    /* jattaa leventamatta eli kaytannossa nollaksi (muuten processxsdata.c */
    /* herjaa) (TVi 2011-12-01) (toi luupin alkuindeksi pit‰‰ kirjoittaa    */
    /* noin ett‰ OpenMP silmukka hyv‰ksyy sen (JLe)) */
    
    if (idx > 0)
      {
        i0 = 1;
        xs[0] = XSS[L1 + 0];
      }
    else
      i0 = 0;
    
    /* Begin parallel loop */
    
#ifdef OPEN_MP
#pragma omp for 
#endif          
    
    for(i = i0; i < NES; i++)
      {  
        /* Hypp‰‰ pisteen yli (j‰t‰ levent‰m‰tt‰), jos integrointi menee */ 
        /* yli energiagridin (TVi 2.9.2011) */
        
        if (lims[i*3 + 1] >= sz - 1)
          {
            xs[i] = XSS[L1 + i];
            continue;
          }
        
        /* v = neutronin nopeus */
        
        v = sqrt(2*XSS[L0 + i]*SPD_C*SPD_C/(1E4*NEUTRON_E0));
        
        /*********************************************************************/
        
        /***** Varsinainen integrointi alkaapi t‰st‰ *************************/
      
        /* Integrointirajat on valmiiksi annettu lims -taulukossa */
    
        intg[i] = 0.0;
        
        /* In case of a threshold reaction, the indexes must be shifted by */
        /* idx (TVi 1.9.2011) */
        
        intidx = lims[i*3] - idx;
        
        if (intidx < 0)
          intidx = 0;
        
        /* Ekstrapoloidaan tarvittaessa nollasta (0, ei ZERO) ekaan */
        /* energiapisteeseen, vaikutusalalla 1/v -approksimaatio    */
        /* j‰lkimm‰inen ehto poistaa ekstrapolaation k‰ytˆst‰       */
        /* threshold-reaktioilla. T‰m‰n kanssa ei tosin havaittu    */
        /* ongelmia, joten ehdon voisi yht‰ hyvin poistaa...        */
        
        if ((lims[i*3] == 0) && (XSS[L0] < 1E-8)) 
          {
            x2 = sqrt(2*XSS[L0]*SPD_C*SPD_C/(1E4*NEUTRON_E0));
            
            intg[i] += XSS[L1]*x2/(2*alpha)*(exp(-alpha*pow(v,2))
                     - exp(-alpha*pow(x2 - v,2))+ sqrt(alpha*PI)*v
                      *(erf(sqrt(alpha)*(x2 - v)) + erf(sqrt(alpha)*v)));
          }

        /* Energiapistev‰lin plow->phigh yli integrointi */
        /* Added threshold correction (TVi 1.9.2011) */
        
        while (intidx < lims[i*3 + 1] - idx)
          {
            x1 = sqrt(2*XSS[L0 + intidx]*SPD_C*SPD_C/(1E4*NEUTRON_E0));
            x2 = sqrt(2*XSS[L0 + intidx + 1]*SPD_C*SPD_C/(1E4*NEUTRON_E0));

            /* Check coincident energy points (10.6.2011 / 1.1.15 JLe) */
      
            if (x1 == x2)
              {
                intidx++;
                continue;
              }
            
            /* Integraalin analyyttinen ratkaisu */
            
            aux = sqrt(alpha*PI)*(v*(3 + 2*alpha*pow(v,2))*(XSS[L1 + intidx] 
                - XSS[L1 + intidx+1]) + (1 + 2*alpha*pow(v,2))
                 *(-x2*XSS[L1 + intidx] + x1*XSS[L1 + intidx+1]));
            
            /* Voi pojat t‰st‰ oli lysti‰ etsi‰ virheit‰ */
      
            intg[i] += 1/(4*pow(alpha,2)*(x2 - x1))*(-2*exp(-alpha
                      *pow(v - x1,2))*((XSS[L1 + intidx] 
                     - XSS[L1 + intidx + 1])*(1 + alpha*pow(v,2)) 
                     + XSS[L1 + intidx]*alpha*(v + x1)*(x1 - x2))
                     - aux*erf(sqrt(alpha)*(v - x1)) 
                     + 2*exp(-alpha*pow(v - x2,2))*((XSS[L1 + intidx]
                     - XSS[L1 + intidx + 1])*(1 + alpha*pow(v,2)) + 
                       alpha*XSS[L1 + intidx+1]*(v + x2)*(x1 - x2))
                     + aux*erf(sqrt(alpha)*(v - x2)));
            intidx++;
          }

        /* Yll‰ oleva integrointi oli +e^(-alffa(V-v)^2) -termille, seuraava */
        /* on -e^(-alffa(V+v)^2)-termille */
      
        /* T‰m‰ integroidaan nollasta (tai thresholdista = indeksi 0) */
        /* intlim-v :hen */
    
        mlow = 0;
        
        /* Ekstrapoloidaan tarvittaessa nollasta (0, ei ZERO) ekaan    */
        /* energiapisteeseen, vaikutusalalla 1/v -approksimaatio       */    
        /* Ehto poistaa ekstrapoloinnin k‰ytˆst‰ threshold-reaktioilla */
        
        if((XSS[L0] < 1E-8) && (lims[i*3 + 2]>0))
          {
            x2 = sqrt(2*XSS[L0]*SPD_C*SPD_C/(1E4*NEUTRON_E0));

            intg[i] -= XSS[L1]*x2/(2*alpha)*(exp(-alpha*pow(v,2))
                     - exp(-alpha*pow(x2 + v,2))- sqrt(alpha*PI)*v
                      *erf(sqrt(alpha)*(x2+v)) + sqrt(alpha*PI)*v
                      *erf(sqrt(alpha)*v));    
          }
   
        /* Energiapistev‰lin mlow->mhigh yli integrointi*/
        
        while(mlow < lims[i*3 + 2] - idx)
          {
            x1 = sqrt(2*XSS[L0 + mlow]*SPD_C*SPD_C/(1E4*NEUTRON_E0));
            x2 = sqrt(2*XSS[L0 + mlow + 1]*SPD_C*SPD_C/(1E4*NEUTRON_E0));
 
            /* Check coincident energy points (10.6.2011 / 1.1.15 JLe) */
            
            if (x1 == x2)
              {
                mlow++;
                continue;
              }
            
            /* Analyyttinen integraalin ratkaisu */
        
            aux = sqrt(alpha*PI)*(v*(-3 - 2*alpha*pow(v,2))*(XSS[L1 + mlow]
                - XSS[L1 + mlow + 1]) + (1 + 2*alpha*pow(v,2))
                 *(-x2*XSS[L1 + mlow] + x1*XSS[L1 + mlow + 1]));
      
            intg[i] -= 1/(4*pow(alpha,2)*(x2 - x1))
                     *(-2*exp(-alpha*pow(v + x1,2))*((XSS[L1 + mlow]
                    - XSS[L1 + mlow + 1])*(1 + alpha*pow(v,2)) 
                    + XSS[L1 + mlow]*alpha*(x1 - v)*(x1 - x2))
                     -aux*erf(sqrt(alpha)*(-v - x1)) 
                    + 2*exp(-alpha*pow(v + x2,2))*((XSS[L1 + mlow]
                    - XSS[L1 + mlow + 1])*(1 + alpha*pow(v,2)) + alpha
                     *XSS[L1 + mlow + 1]*(x2 - v)*(x1 - x2)) + aux
                     *erf(sqrt(alpha)*(-v - x2)));
            mlow++;
          }

        /* Ensin v‰liaikaiseen s‰ilˆˆn.  */

        xs[i] = (1/pow(v,2))*sqrt(alpha/PI)*intg[i];
      }
  }

  /* Siirret‰‰n vaikutusalat */

  memmove(&XSS[L1], xs, NES*sizeof(double));
  
#ifdef DEBUG

  /* Tarkistetaan negatiiviset, raja on sama kuin interpolatedata.c:ss‰ */
  /* (JLE 31.8.2011) */
  
  for (i = 0; i < NES; i++)
    if (XSS[L1 + i] < -1.0)
      Warn(FUNCTION_NAME, "Negative cross sections after broadening");
  
#endif

  Mem(MEM_FREE, intg);
  Mem(MEM_FREE, xs);
}
/*****************************************************************************/
