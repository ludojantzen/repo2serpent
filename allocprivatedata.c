/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : allocprivatedata.c                             */
/*                                                                           */
/* Created:       2011/11/10 (JLe)                                           */
/* Last modified: 2018/03/14 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Allocates and re-arranges memory for OpenMP private data     */
/*              and scoring buffers.                                         */
/*                                                                           */
/* Comments: - NOTE: memory sizes stored in DATA array are per private       */
/*                   segment.                                                */
/*                                                                           */
/*           - Jos BUF-blokin koko ylittää cachen koon niin BufVal, ym.      */
/*             funktiot jotka summaa sen yli hidastuu ihan järjettömästi.    */
/*                                                                           */
/*           - Tätä muutettiin 25.2.2014 (2.1.18) siten että muistialueen    */
/*             koon kasvaessa olemassa olevaa dataa ei kopioida uuteen.      */
/*             Oletus on se, että PRIVA, BUF ja RES2 -blokkien data on       */
/*             aina luonteeltaan väliaikaista. Silloin kun dataa käytetään,  */
/*             sitä ei ole tarvetta varata lisää, ja silloin kun varataan    */
/*             sisältöä ei enää tarvita muuhun. Muistinvaraus väärässä       */
/*             paikassa on estetty DATA_ALLOW_MEM_OP -flagillä. Sen lisäksi  */
/*             tässä funktiossa varatut blokit pitää vapauttaa käyttöön      */
/*             ExpandPrivateArrays() -funktiolla, joka ei tällä hetkellä     */
/*             tee muuta kuin muuttaa DATA_PRIVA_MEM_READY -flagin tilaa.    */
/*             Jos sisällön säilyttäminen muistinvarauksen jälkeen ei        */
/*             tosiaan ole tarpeen, tuon funktion voinee korvata uudella     */
/*             Mem()-funktion allow -optiolla.                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AllocPrivateData:"

/*****************************************************************************/

long AllocPrivateData(long sz, long type)
{
  long memsize, realsize, loc0, block, nseg;
  double *dat;

  /* Check allow flag */

  if (RDB != NULL)
    if ((long)RDB[DATA_ALLOW_MEM_OP] == NO)
      Die(FUNCTION_NAME, "Memory allocation not allowed");

  /* Check OpenMP thread number */

  if (OMP_THREAD_NUM > 0)
    Die(FUNCTION_NAME, "Not allowed in threads");
  
  /* Number of segments */
  
  if ((type == BUF_ARRAY) && ((long)RDB[DATA_OPTI_SHARED_BUF] == YES))
    nseg = 1;
  else if ((type == RES2_ARRAY) && ((long)RDB[DATA_OPTI_SHARED_RES2] == YES))
    nseg = 1;
  else if (type == RES3_ARRAY)
    nseg = 1;
  else
    nseg = (long)RDB[DATA_OMP_MAX_THREADS];
  
  /* Avoid compiler warning */
  
  memsize = -1;
  realsize = -1;
  dat = NULL;
  
  /* Get memory size and pointer */
  
  if (type == PRIVA_ARRAY)
    {
      memsize = (long)RDB[DATA_ALLOC_PRIVA_SIZE];
      realsize = (long)RDB[DATA_REAL_PRIVA_SIZE];

      dat = PRIVA;
    }
  else if (type == BUF_ARRAY)
    {
      memsize = (long)RDB[DATA_ALLOC_BUF_SIZE];
      realsize = (long)RDB[DATA_REAL_BUF_SIZE];

      dat = BUF;
    }
  else if (type == RES2_ARRAY)
    {
      memsize = (long)RDB[DATA_ALLOC_RES2_SIZE];
      realsize = (long)RDB[DATA_REAL_RES2_SIZE];

      dat = RES2;
    }
  else if (type == RES3_ARRAY)
    {
      memsize = (long)RDB[DATA_ALLOC_RES3_SIZE];
      realsize = (long)RDB[DATA_REAL_RES3_SIZE];

      dat = RES3;
    }
  else
    Die(FUNCTION_NAME, "Invalid array type");

  /* Set block size */

  block = (long)(5*(long)MEGA/sizeof(double));

  /* Get pointer to free memory */
  
  loc0 = memsize;
  
  /* Check if data can fit into existing segment */
  
  if (memsize + sz < realsize)
    {
      /***********************************************************************/
      
      /***** Sufficient segment size *****************************************/
      
      /* Update data size */
      
      memsize = memsize + sz;
      
      /***********************************************************************/
    }
  else
    {
      /***********************************************************************/

      /***** Allocate more memory ********************************************/

      /* Calculate new size */
      
      if (sz > block)
        realsize = realsize + sz + block;
      else
        realsize = realsize + block;
      
      /* Free array */
      
      if (dat != NULL)
        Mem(MEM_FREE, dat);

      /* Allocate more memory and clear data */

      dat = (double *)Mem(MEM_ALLOC, nseg*realsize, sizeof(double));

      /* Update size */
      
      memsize = memsize + sz;
         
      /***********************************************************************/
    }
  
  /* Put new sizes and pointer */
  
  if (type == PRIVA_ARRAY)
    {
      WDB[DATA_ALLOC_PRIVA_SIZE] = (double)memsize;
      WDB[DATA_REAL_PRIVA_SIZE] = (double)realsize;
      PRIVA = dat;
    }
  else if (type == BUF_ARRAY)
    {
      WDB[DATA_ALLOC_BUF_SIZE] = (double)memsize;
      WDB[DATA_REAL_BUF_SIZE] = (double)realsize;
      BUF = dat;
    }
  else if (type == RES2_ARRAY)
    {
      WDB[DATA_ALLOC_RES2_SIZE] = (double)memsize;
      WDB[DATA_REAL_RES2_SIZE] = (double)realsize;
      RES2 = dat;
    }
  else
    {
      WDB[DATA_ALLOC_RES3_SIZE] = (double)memsize;
      WDB[DATA_REAL_RES3_SIZE] = (double)realsize;
      RES3 = dat;
    }
  
  /* Check pointer */
  
  if (memsize > realsize)
    Die(FUNCTION_NAME, "Real data size exceeded");
  else if (loc0 > memsize)
    Die(FUNCTION_NAME, "Pointer beyond memory size");
  
  /* Calculate data size */
  
  CalculateBytes();

 /* Deny access until blocks are expanded (typerä termi, korvaa) */

  WDB[DATA_PRIVA_MEM_READY] = (double)NO;

  /* Return pointer to allocated memory */
  
  return loc0;
}

/*****************************************************************************/
