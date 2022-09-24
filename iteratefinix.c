/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : iteratefinix.c                                 */
/*                                                                           */
/* Created:       2013/04/11 (VVa)                                           */
/* Last modified: 2017/02/01 (VVa)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Iterates Finix for temperature feedback                      */
/*                                                                           */
/* Comments: - Tää ei välttämättä vielä toimi reproducible MPI -moodissa     */
/*           - Erotettu iteratetfb:stä                                       */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifdef FINIX

#define FUNCTION_NAME "IterateFinix:"

#define MAX_ITER 1000

/*****************************************************************************/

void IterateFinix()
{
  long fib, fpe;

  if((fib = (long)RDB[DATA_PTR_FIN0]) < VALID_PTR)
    return;

  fprintf(outp, "Iterating fuel behavior module FINIX\n");


  /* Loop over pins */

  while (fib > VALID_PTR)
    {

      /* Get pointer to interface pin block */

      fpe = (long)RDB[fib + FINIX_PTR_FUEP];

      /* Update power densities */

      UpdateFinixPower(fib, fpe);

      /* Start timer */
  
      StartTimer(TIMER_FINIX);

      /* Run Finix */

      RunFinix(fib, fpe);

      /* Stop timer */

      StopTimer(TIMER_FINIX);

      /* Next feedback */
      
      fib = NextItem(fib);
    }

  /* Update interface structures */

  UpdateFinixIFC();

  /* Print output TODO: Maybe do not print at every iteration? */
  /*
  PrintFinix();
  */

  WriteFinixIFC();

  /* Print momentary power distributions TODO: Maybe do not print at every iteration? */

  PrintInterfaceOutput();

  /* Save current power and temperature distributions */

  /*
  if((((long)WDB[DATA_CYCLE_IDX] - (long)RDB[DATA_CRIT_SKIP] + 1)*8 % 128) == 0)
    {
      sprintf(str,"cp ./Finout.m ./Finout%ld.m",(long)WDB[DATA_SOL_REL_ITER]);

      if(system(str) != 0)
        Die(FUNCTION_NAME,"Could not copy serpent input");

      sprintf(str,"cp ./Finixifcout.m ./Finixifcout%ld.m",(long)WDB[DATA_SOL_REL_ITER]);

      if(system(str) != 0)
        Die(FUNCTION_NAME,"Could not copy serpent input");
    }
  */

  fprintf(outp, "OK.\n\n");

}

#endif 

/*****************************************************************************/
