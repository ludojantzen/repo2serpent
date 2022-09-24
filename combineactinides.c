/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : combineactinides.c                             */
/*                                                                           */
/* Created:       2013/01/24 (JLe)                                           */
/* Last modified: 2018/04/18 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Combines actinides produced by neutron reactions into a      */
/*              single list of nuclides                                      */
/*                                                                           */
/* Comments: - NOTE: tässä oletetaan että aktinideilla ei ole peräkkäisiä    */
/*                   hajoamismoodeja                                         */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CombineActinides:"

/* Recursive subroutine */

void CombineActinides0(long, double *, long *);

/*****************************************************************************/

void CombineActinides()
{
  long mat, iso, nuc, ZAI, sz, ptr, n;
  double *zai;

  /* Check burnup calculation mode */

  if ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == NO)
    return;

  /* Calculate array size */

  sz = 0;
  ptr = (long)RDB[DATA_PTR_ACE0];
  while (ptr > 0)
    {
      /* Add count */

      sz++;
      
      /* next */
      
      ptr = (long)ACE[ptr + ACE_PTR_NEXT];
    }

  /* Allocate memory for ZAI vector */

  zai = (double *)Mem(MEM_ALLOC, sz, sizeof(double));
  sz = 0;

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check burn-flag, division and fixed flag */
          
      if (((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT) &&
          ((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] < VALID_PTR))
        {
          /* Loop over composition */
              
          iso = (long)RDB[mat + MATERIAL_PTR_COMP];
          while (iso > VALID_PTR)
            {
              /* Pointer to nuclide */

              nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Get ZAI */

              ZAI = (long)RDB[nuc + NUCLIDE_ZAI];

              /* Find daughters */

              if ((ZAI > 890000) && ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & 
                                     NUCLIDE_FLAG_NEW_DAUGHTERS))
                CombineActinides0(ZAI, zai, &sz);

              /* Next nuclide in composition */
      
              iso = NextItem(iso);
            }
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /* Check count */

  if (sz > 0)
    {
      /* Sort array */

      SortArray(zai, sz);

      /* Check duplicates and order */

      for (n = 1; n < sz; n++)
        if (zai[n] <= zai[n - 1])
          Die(FUNCTION_NAME, "Error in list");

      /* Allocate memory for list */
      
      ptr = ReallocMem(DATA_ARRAY, sz + 1);
      WDB[DATA_PTR_AC_ZAI_LIST] = (double)ptr;

      /* Read data */
      
      for (n = 0; n < sz; n++)
        WDB[ptr++] = zai[n];
      
      /* Null terminator */
      
      WDB[ptr] = -1.0;
    }

  /* Free memory */

  Mem(MEM_FREE, zai);
}

/*****************************************************************************/

/***** Recursive subroutine **************************************************/

void CombineActinides0(long ZAI, double *zai, long *sz)
{
  long ace, n, ZA, Z, ptr, dec, rtyp, rfs, new;
  double T;

  /* Find existing */

  for (n = 0; n < *sz; n++)
    if (zai[n] == (double)ZAI)
      return;

  /* Separate Z */
  
  Z = (long)((double)ZAI/10000.0);

  /* Check limits */

  if ((Z < (long)RDB[DATA_BU_ACT_MIN_Z]) || (Z > (long)RDB[DATA_BU_ACT_MAX_Z]))
    return;

  /* Add to list */

  zai[*sz] = (double)ZAI;
  *sz = *sz + 1;

  /* Find transport data */
              
  ace = (long)RDB[DATA_PTR_ACE0];          
  while (ace > 0)
    {
      /* Compare ZAI */
 
      if (ZAI == (long)ACE[ace + ACE_ZAI])
        break;
      
      /* This assumes that decay data follows transport data */

      if ((long)ACE[ace + ACE_TYPE] == NUCLIDE_TYPE_DECAY)
        {
          /* Reset pointer */

          ace = -1;

          /* Break loop */

          break;
        }

      /* next */
      
      ace = (long)ACE[ace + ACE_PTR_NEXT];
    }

  /* Check pointer */
  
  if (ace > VALID_PTR)
    {
      /* Get ZA */

      ZA = (long)ACE[ace + ACE_ZA];

      /* (n,gamma) reaction */

      CombineActinides0(10*ZA + 10, zai, sz);
      CombineActinides0(10*ZA + 11, zai, sz);
      
      /* (n,2n) reaction */

      CombineActinides0(10*ZA - 10, zai, sz);
    }

  /* Find decay data */
              
  if ((ace = (long)RDB[DATA_PTR_DECAY_ACE0]) < 1)
    ace = (long)RDB[DATA_PTR_ACE0];
    
  while (ace > 0)
    {
      /* Compare ZAI */

      if ((ZAI == (long)ACE[ace + ACE_ZAI]) && 
          ((long)ACE[ace + ACE_TYPE] == NUCLIDE_TYPE_DECAY))
        break;
      
      /* next */
      
      ace = (long)ACE[ace + ACE_PTR_NEXT];
    }

  /* Check pointer */

  if (ace > VALID_PTR)
    {
      /* Get ZA */

      ZA = (long)ACE[ace + ACE_ZA];

      /* Loop over decay modes */

      if ((ptr = (long)ACE[ace + ACE_PTR_DECAY_LIST]) > 0)
        while ((dec = (long)ACE[ptr++]) > 0)
          {
            /* Get half-life for mode */
            
            T = log(2.0)/(ACE[ace + ACE_LAMBDA]*ACE[dec + DECAY_BR] + ZERO);
            
            /* Check limit (10 years) */

            if (T < 10*365*24*60*60)
              {
                /* Get type */

                rtyp = (long)ACE[dec + DECAY_RTYP1] - 10000;
                rfs = (long)ACE[dec + DECAY_RFS];

                /* Check second mode */

                if ((long)ACE[dec + DECAY_RTYP2] > 10000)
                  Note(0, "Secondary decay mode (%s) for %ld ignored", 
                       ReactionMT((long)ACE[dec + DECAY_RTYP2], NO), ZAI);

                /* Check */

                if (rtyp == 1)
                  {
                    /* Beta decay */

                    new = 10*(ZA + 1000) + rfs;
                  }
                else if (rtyp == 2)
                  {
                    /* Electron capture / positron emission */

                    new = 10*(ZA - 1000) + rfs;
                  }
                else if (rtyp == 3)
                  {
                    /* Isomeric transition */

                    new = 10*ZA;
                  }
                else if (rtyp == 4)
                  {
                    /* Alpha decay */

                    new = 10*(ZA - 2004) + rfs;
                  }
                else if (rtyp == 5)
                  {
                    /* Neutron emission */

                    new = 10*(ZA - 1) + rfs;
                  }
                else if (rtyp == 7)
                  {
                    /* Proton emission */

                    new = 10*(ZA - 1001) + rfs;
                  }
                else
                  new = -1;

                /* Check new and call recursively */

                if ((new > 0) && (labs(new - ZAI) > 20041))
                  Die(FUNCTION_NAME, "Invalid decay mode (%ld %ld %ld %ld)", 
                      ZA, new, rtyp, rfs);
                else if (new > 0)
                  CombineActinides0(new, zai, sz);
              }
          }
    }
}

/*****************************************************************************/
