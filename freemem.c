/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : freemem.c                                      */
/*                                                                           */
/* Created:       2010/09/15 (JLe)                                           */
/* Last modified: 2018/03/14 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Frees memory allocated to data blocks                        */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FreeMem:"

/*****************************************************************************/

void FreeMem()
{
  /* Free FINIX data */

#ifdef FINIX

  FreeFinix();

#endif

  /* Free data arrays */

  if (ACE != NULL)
    Mem(MEM_FREE, ACE);

  if (WDB != NULL)
    Mem(MEM_FREE, WDB);

  if (RES1 != NULL)
    Mem(MEM_FREE, RES1);

  if (RES2 != NULL)
    Mem(MEM_FREE, RES2);

  if (RES3 != NULL)
    Mem(MEM_FREE, RES3);

  if (BUF != NULL)
    Mem(MEM_FREE, BUF);

  if (PRIVA != NULL)
    Mem(MEM_FREE, PRIVA);

  if (ASCII != NULL)
    Mem(MEM_FREE, ASCII);

  if (SEED != NULL)
    Mem(MEM_FREE, SEED);

  if (SEED0 != NULL)
    Mem(MEM_FREE, SEED0);

  if (mpiid > 0)
    fclose(outp);
}

/*****************************************************************************/
