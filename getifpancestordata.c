/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : getifpancestordata.c                           */
/*                                                                           */
/* Created:       2018/06/12 (VVa)                                           */
/* Last modified: 2018/09/17 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Retrieves some data from IFP ancestors.                      */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "GetIFPAncestorData:"

/*****************************************************************************/

void GetIFPAncestorData(long part, double *lifetime, long *dngroup,
                        double *lambda, double *nvoidcoll)
{
  long loc1, gen;

#ifdef OLD_IFP

      /* Get delayed neutron group, lifetime and net number of accepted */
      /* collisions from predecessors */

      if ((loc1 = (long)RDB[part + PARTICLE_PTR_FISS_PROG]) > VALID_PTR)
        {
          /* Get first ancestor  */

          /* gen starts from 1 because if the loop runs once, the progeny */
          /* with index 1 is obtained using NextItem etc. Latest progeny  */
          /* has index 0 */

          for (gen = 1; gen < (long)RDB[DATA_IFP_CHAIN_LENGTH] - 1; gen++)
            {
              if (loc1 < VALID_PTR)
                break;

              loc1 = NextItem(loc1);
            }

          if (loc1 > VALID_PTR)
            {
              /* Get data from first ancestor */

              *lifetime = (double)RDB[loc1 + FISS_PROG_LIFETIME];
              *lambda = (double)RDB[loc1 + FISS_PROG_LAMBDA];
              *dngroup = (long)RDB[loc1 + FISS_PROG_DN_GROUP];
              *nvoidcoll = (double)RDB[loc1 + FISS_PROG_NA_COLL];
            }
        }
      else
        Die(FUNCTION_NAME, "Particle has no fission progeny?");

#else

      /* Get delayed neutron group, lifetime and net number of accepted */
      /* collisions from predecessors */

      if ((loc1 = (long)RDB[part + PARTICLE_PTR_EVENTS]) > VALID_PTR)
        {
          /* Reset current generation */

          gen = -1;

          /* Loop over event history */

          while (loc1 > VALID_PTR)
            {
              /* Check if event was a fission */

              if ((long)RDB[loc1 + EVENT_TYPE] == EVENT_TYPE_FISS)
                {
                  /* Increment generation due to fission */

                  gen++;

                  /* Check index */

                  if (gen == (long)RDB[DATA_IFP_CHAIN_LENGTH] - 1)
                    {
                      /* Get data from the earliest ancestor */

                      /* Get delayed neutron group */

                      if ((ng = (long)RDB[loc1 + EVENT_DN_GROUP]) > 0)
                        *dngroup = ng;

                      /* Get neutron lifetime */

                      if ((t = (double)RDB[loc1 + EVENT_LIFETIME]) > 0)
                        *lifetime = t;

                      /* Get decay constant */

                      if ((t = (double)RDB[loc1 + EVENT_LAMBDA]) > 0)
                        *lambda = t;

                      /* Get number of collisions in void fraction */
                      /* material */

                      *nvoidcoll = (long)RDB[loc1 + EVENT_ACCEPTED_COLLS];

                      break;
                    }

                }
              /* Pointer to next event */

              loc1 = NextItem(loc1);
            }
        }

#endif

}
