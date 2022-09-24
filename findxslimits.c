/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : findxslimits.c                                 */
/*                                                                           */
/* Created:       2011/07/28 (JLe)                                           */
/* Last modified: 2011/12/16 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Laskee limits-vektoriin integrointirajat ja huolehtii        */ 
/*              tarvittaessa muistin varauksesta. Muistia varataan           */
/*              oletusarvoisesti STD_LIMITS_SIZE verran, mutta funktio osaa  */
/*              hoitaa myös tilanteen, jossa vaikutusalan energiapisteiden   */
/*              määrä on suurempi kuin STD_LIMITS_SIZE. (esmes nopeissa      */
/*              kirjastoissa täysin mahdollista, jopa todennäköistä)         */
/*                                                                           */
/*              Limits-vektori sisältää rajat:                               */
/*              i*3+0 = plus-termin integroinnin alaraja                     */
/*              i*3+1 = - '' - yläraja                                       */
/*              i*3+2 = miinus-termin yläraja                                */
/*                                                                           */
/* Comments: - From Serpent 1.1.15                                           */
/*           - Original routine developed by Tuomas Viitanen                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "FindXSLimits:"

/* Tämä vakio kertoo limits-vektorin "standardikoon", eli kuinka monelle
   vaikutusalapisteelle varataan oletusarvoisesti muistia. Homma toimii
   vaikka arvo olisi liian pieni, mutta hitusen hitaammin. Vastaavasti myös
   turhan suuret arvot hidastavat ohjelmaa ja lisäävät muistin tarvetta. */

#define STD_LIMITS_SIZE 70000

/*****************************************************************************/

long *FindXSLimits(long L0, long NES, double T, double Td, double awr, 
		   double *XSS, long *limits)
{
  long i;
  double alpha, intlim, v;

  /* Varataan muistia, jos limits ei ole alustettu tai se on liian pieni */

  if ((limits == NULL) && (NES < STD_LIMITS_SIZE))
    limits = (long *)Mem(MEM_ALLOC, STD_LIMITS_SIZE*3, sizeof(long));
  else if ((NES > STD_LIMITS_SIZE) && (limits != NULL))
    limits = (long *)Mem(MEM_REALLOC, limits, NES*3*sizeof(long));
  else if ((NES > STD_LIMITS_SIZE) && (limits == NULL))
    limits=(long *)Mem(MEM_ALLOC, NES*3, sizeof(long));
  
  alpha=(awr*NEUTRON_E0*1E4)/(SPD_C*SPD_C*2*KELVIN*(Td-T));  

  /* Integrointialueeksi valitaan n=4 "keskihajontaa", jonka havaittiin
   olevan riittävän suuri lukema. (>=3.5 on about ok, NJOYlla n=4) */

  intlim=4/sqrt(alpha);

  limits[0]=0;
  
  /* Etsitään sopivat integrointirajat jokaiselle energiapisteelle */
  
  for(i=0; i<NES; i++){
    
    v=sqrt(2*XSS[L0+i]*SPD_C*SPD_C/(1E4*NEUTRON_E0));

    /* Tämä if-lause säästää paljon iterointiaikaa, kun etsintää ei tarvitse
       aloittaa nollasta. Alarajat kasvavat monotonisesti. */
    if(i!=0)
      limits[i*3]=limits[(i-1)*3];
        
    /*e^-alpha(V-v)^2 -termin integrointi v-4/sqrt(alpha) :sta v+4/sqrt(alpha):n */
    
    while(sqrt(2*XSS[L0+limits[i*3]+1]*SPD_C*SPD_C/(1E4*NEUTRON_E0))<v-intlim){
      limits[i*3]++;
    }

    limits[i*3+1]=limits[i*3];

    while(sqrt(2*XSS[L0+limits[i*3+1]]*SPD_C*SPD_C/(1E4*NEUTRON_E0))<v+intlim && 
	  limits[i*3+1]<NES){
      limits[i*3+1]++;
    }
 
    /*e^-alpha(V+v)^2 -termin integrointi nollasta 4/sqrt(alpha)-v :hen */

    limits[i*3+2]=0;

    while(sqrt(2*XSS[L0+limits[i*3+2]]*SPD_C*SPD_C/(1E4*NEUTRON_E0))<intlim-v){
	limits[i*3+2]++;
    }
  }  
  return limits;
}

/*****************************************************************************/
