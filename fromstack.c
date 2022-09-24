/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : fromstack.c                                    */
/*                                                                           */
/* Created:       2011/03/09 (JLe)                                           */
/* Last modified: 2019/04/03 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Retrieves neutron / photon from stack                        */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FromStack:"

/*****************************************************************************/

long FromStack(long type, long id)
{
  long ptr, loc0, sz, n;

#ifdef OLD_IFP

  long prg;

#endif

  /* Check id */

  if ((id < 0) || (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1))
    Die(FUNCTION_NAME, "Error in thread id");

  /* Avoid compiler warning */

  ptr = -1;

  /* Check type and get pointer to stack */

  if (type == PARTICLE_TYPE_NEUTRON)
    ptr = (long)RDB[OMPPtr(DATA_PART_PTR_NSTACK, id)];
  else if (type == PARTICLE_TYPE_GAMMA)
    ptr = (long)RDB[OMPPtr(DATA_PART_PTR_GSTACK, id)];
  else if (type == PARTICLE_TYPE_PRECURSOR)
    ptr = (long)RDB[OMPPtr(DATA_PART_PTR_PSTACK, id)];
  else
    Die(FUNCTION_NAME, "Invalid particle type");

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get size */

  sz = ListSize(ptr) - 2;

  /* Avoid compiler warning */

  loc0 = -1;

  /* Compare to minimum */

  if (type == PARTICLE_TYPE_NEUTRON)
    loc0 = (long)RDB[DATA_PART_PTR_MIN_NSTACK];
  else if (type == PARTICLE_TYPE_GAMMA)
    loc0 = (long)RDB[DATA_PART_PTR_MIN_GSTACK];
  else if (type == PARTICLE_TYPE_PRECURSOR)
    loc0 = (long)RDB[DATA_PART_PTR_MIN_PSTACK];
  else
    Die(FUNCTION_NAME, "Invalid particle type");

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(loc0)", PRIVA_ARRAY, loc0);

  /* Compare to minimum */

  if ((n = (long)GetPrivateData(loc0, id)) == 0)
    PutPrivateData(loc0, sz, id);
  if (sz < n)
    PutPrivateData(loc0, sz, id);

  /* Get pointer to last item */

  ptr = LastItem(ptr);
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Check type */

  if ((long)RDB[ptr + PARTICLE_TYPE] == PARTICLE_TYPE_DUMMY)
    {
      if (type == PARTICLE_TYPE_NEUTRON)
        Error(0, "Insufficient neutron buffer size, increase value of\nparameter \"nbuf\" (currently set to %1.1f)", RDB[DATA_PART_NBUF_FACTOR]);
      else if (type == PARTICLE_TYPE_GAMMA)
          Error(0, "Insufficient photon buffer size, increase value of\nparameter \"gbuf\" (currently set to %1.1f)", RDB[DATA_PART_GBUF_FACTOR]);
      else
          Error(0, "Insufficient precursor buffer size, increase value of\nparameter \"pbuf\" (currently set to %1.1f)", RDB[DATA_PART_PBUF_FACTOR]);
    }

  /* Remove particle from stack */

  RemoveItem(ptr);

  /* Remember pointers to fission progeny and history data */

#ifdef OLD_IFP

  prg = (long)RDB[ptr + PARTICLE_PTR_FISS_PROG];

#endif

  /* Wipe data */

  memset(&WDB[ptr + LIST_DATA_SIZE], 0.0,
         (PARTICLE_BLOCK_SIZE - LIST_DATA_SIZE)*sizeof(double));

  /* Put type */

  WDB[ptr + PARTICLE_TYPE] = (double)type;

  /* Put pointers */

  WDB[ptr + PARTICLE_PTR_EVENTS] = NULLPTR;
  WDB[ptr + PARTICLE_PTR_SENS_EBLOCK] = NULLPTR;

#ifdef OLD_IFP

  WDB[ptr + PARTICLE_PTR_FISS_PROG] = (double)prg;

#endif

  /* Loop over progeny data */

#ifdef OLD_IFP

  while (prg > VALID_PTR)
    {
      /* Reset data */

      WDB[prg + FISS_PROG_DN_GROUP] = 0.0;
      WDB[prg + FISS_PROG_LIFETIME] = 0.0;
      WDB[prg + FISS_PROG_LAMBDA] = 0.0;
      WDB[prg + FISS_PROG_NA_COLL] = 0.0;

      /* Pointer to next */

      prg = NextItem(prg);
    }

#endif

  /* Reset ICM data */

  WDB[ptr + PARTICLE_ICM_PTR_ICM] = -1.0;
  WDB[ptr + PARTICLE_ICM_IDX] = -1.0;
  WDB[ptr + PARTICLE_ICM_MUA] = -1.0;
  WDB[ptr + PARTICLE_ICM_MUS] = -1.0;
  WDB[ptr + PARTICLE_ICM_G] = -1.0;
  WDB[ptr + PARTICLE_ICM_WGT] = -1.0;

  /* Reset albedo stuff */

  WDB[ptr + PARTICLE_ALB_PTR_GCU] = -1.0;
  WDB[ptr + PARTICLE_ALB_SURF_IDX] = -1.0;
  WDB[ptr + PARTICLE_ALB_G] = -1.0;

  /* Return pointer */

  return ptr;
}

/*****************************************************************************/
