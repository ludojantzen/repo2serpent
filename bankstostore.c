/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : bankstostore.c                                 */
/*                                                                           */
/* Created:       2015/01/19 (VVa)                                           */
/* Last modified: 2018/01/25 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Stores dynamic mode neutrons in the end of time interval     */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "BanksToStore:"

/*****************************************************************************/

void BanksToStore()
{
  long ptr, id, nb, tt, tot, flag;
  char tmpstr[MAX_STR];
  FILE *fp;

  /* Stores are only used in coupled transient calculations */

  if (RDB[DATA_SIMULATION_MODE] != SIMULATION_MODE_DYN)
    return;

#ifdef DNPRINT
  fprintf(outp, "bankstostore.c -->\n");
#endif

  /* Loop over threads */

  /*printf("Storing neutrons to 3rd bank\n");*/

  nb = (long)RDB[DATA_CYCLE_IDX];

  /* Get total number of threads */

  tt = (long)RDB[DATA_OMP_MAX_THREADS];

  /* Calculate total number of particles in banks */

  tot = 0;

  for (id = 0; id < tt; id++)
    {
      /* Add bank size */

      ptr = (long)RDB[OMPPtr(DATA_PART_PTR_BANK, id)];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      tot = tot + ListSize(ptr) - 1;
    }

  /* Create store file name*/

  flag = (long)RDB[DATA_EOI_STORE_NAME];
  CheckValue(FUNCTION_NAME, "flag", "", flag, 0, 1);

  if (flag == 0)
    sprintf(tmpstr, "%s.storeA%d", GetText(DATA_PTR_INPUT_FNAME), mpiid);
  else
    sprintf(tmpstr, "%s.storeB%d", GetText(DATA_PTR_INPUT_FNAME), mpiid);

  /* Open file for writing */

  if (nb == 0)
    fp = fopen(tmpstr, "w");
  else
    fp = fopen(tmpstr, "a");

  /* Write batch number */

  fwrite(&nb, sizeof(long), 1, fp);

  /* Write number of particles that follows */

  fwrite(&tot, sizeof(long), 1, fp);

#ifdef DNPRINT
  fprintf(outp, "<-- Storing %ld particles for batch %ld to %s.\n", tot, nb, tmpstr);
#endif

  /* Loop over threads to write particles to file */

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
   {

     while ((ptr = FromBank(id)) > VALID_PTR)
       {
         /* Check type */

         if ((long)RDB[ptr + PARTICLE_TYPE] != PARTICLE_TYPE_NEUTRON)
           Die(FUNCTION_NAME, "Invalid particle type");

         /* Write particle to file */

         fwrite(&RDB[ptr + LIST_DATA_SIZE], sizeof(double), (PARTICLE_BLOCK_SIZE - LIST_DATA_SIZE), fp);

         /* Put original to stack */

         ToStack(ptr,id);

         /* Count particle*/

         tot--;
       }

   }

  /* Close the store file */

  fclose(fp);

  /* Check that we stored as many particles as we were thinking to */

  if (tot != 0)
    Die(FUNCTION_NAME, "Stored %ld particles too few", tot);

#ifdef DNPRINT
  fprintf(outp, "<-- bankstostore.c\n");
#endif

}
