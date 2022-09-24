/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : addbranching.c                                 */
/*                                                                           */
/* Created:       2010/09/11 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Creates branches to isomeric states, secondary products      */
/*              and fission product yields at different energies             */
/*                                                                           */
/* Comments: - Noiden fixed-arvojen lisäksi data pitää lukea käyttäjältä     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AddBranching:"

/*****************************************************************************/

void AddBranching(long nuc)
{
  long mt0[20], mt, mt1, rea, loc0, ptr, yld, n, ni;
  double br0[20], Q;

  /* Check DBRC flag */

  if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DBRC)
    return;

  /***************************************************************************/

  /***** Energy-dependent branching ratios ***********************************/

  /* Check type */

  if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_TRANSPORT)
    {
      /* Loop over data */

      loc0 = (long)RDB[DATA_PTR_BRA_LIST];
      while (loc0 > VALID_PTR)
        {
          /* Compare ZAI */

          if ((long)RDB[loc0 + BRA_LIST_ZAI] != (long)RDB[nuc + NUCLIDE_ZAI])
            {
              /* Next in list */

              loc0 = NextItem(loc0);

              /* Cycle loop */

              continue;
            }

          /* Find reaction */

          rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
          while (rea > VALID_PTR)
            {
              /* Compare mt */

              if (((long)RDB[loc0 + BRA_LIST_MT] == 
                   (long)RDB[rea + REACTION_MT]) && 
                  ((long)RDB[rea + REACTION_RFS] == 0))
                break;
              
              /* Next reaction */
              
              rea = NextItem(rea);
            }

          /* Check if found */

          if (rea < VALID_PTR)
            {
              /* Next in list */

              loc0 = NextItem(loc0);

              /* Cycle loop */

              continue;
            }
                    
          /* Set flag and type */

          SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_BRA_DATA);

          /* Put state */

          WDB[rea + REACTION_RFS] = 0.0;
          
          /* Put pointer to data */

          WDB[rea + REACTION_PTR_ISO_BRA] = (double)loc0;

          /* Create duplicate */
          
          ptr = DuplicateItem(rea);
          
          /* Put state */

          WDB[ptr + REACTION_RFS] = 1.0;

          /* Set pointer to parent reaction */
          
          WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
          
          /* Set type */
          
          WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
          
          /* Set branch mt */
          
          WDB[ptr + REACTION_BRANCH_MT] = RDB[ptr + REACTION_MT];

          /* Switch types if inelastic scattering from isomeric to ground */
          /* state */

          if (((long)RDB[rea + REACTION_MT] == 4) && 
              ((long)RDB[nuc + NUCLIDE_I] == 1))
            {
              WDB[rea + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
              WDB[rea + REACTION_PTR_BRANCH_PARENT] = (double)ptr;
              WDB[rea + REACTION_BRANCH_MT] = RDB[rea + REACTION_MT];
              WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_SPECIAL;
              WDB[ptr + REACTION_PTR_BRANCH_PARENT] = 0.0;
              WDB[ptr + REACTION_BRANCH_MT] = 0.0;
            }

          /* Next in list */

          loc0 = NextItem(loc0);
        }
    }

  /***************************************************************************/

  /***** Fixed branching ratios to isomeric states ***************************/

  /* Check type and if flag is already set */

  if (((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_TRANSPORT) &&
      (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_BRA_DATA)))
    {
      /* Reset count */
      
      ni = 0;

      /* Collect branching reactions */

      ptr = (long)RDB[DATA_PTR_FIX_BRA0];
      while (ptr > VALID_PTR)
        {
          /* Compare ZAI */

          if ((long)RDB[ptr + FIX_BRA_ZAI] == (long)RDB[nuc + NUCLIDE_ZAI])
            {
              /* Check count (leads to segfault) */

              if (ni > 19)
                Error(0, "Maximum number of fixed branching ratios exceeded");

              /* Set fraction and mt */

              br0[ni] = RDB[ptr + FIX_BRA_FRAC];
              mt0[ni++] = (long)RDB[ptr + FIX_BRA_MT];
            }
          
          /* Pointer to next */

          ptr = NextItem(ptr);
        }
      
      /* Loop over branch reactions */
      
      for (n = 0; n < ni; n++)
        {
          /* Loop over reaction channels to find correct mt */
          
          rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
          while (rea > VALID_PTR)
            {
              /* Check MT */
              
              if ((long)RDB[rea + REACTION_MT] == mt0[n])
                break;
              
              /* Next */
              
              rea = NextItem(rea);
            }
          
          /* Check if found */

          if (rea < VALID_PTR)
            continue;
                      
          /* Skip if data was already read */

          if ((long)RDB[rea + REACTION_PTR_ISO_BRA] > VALID_PTR)
            continue;

          /* Set flag and type */
          
          SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_BRA_DATA);

          /* Put state */

          WDB[rea + REACTION_RFS] = 0.0;
          
          /* Create structure */

          loc0 = NewItem(rea + REACTION_PTR_ISO_BRA, BRA_LIST_BLOCK_SIZE);
          
          /* Put fixed ratio */

          WDB[loc0 + BRA_LIST_FIX_FRAC] = br0[n];

          /* Create duplicate */
          
          ptr = DuplicateItem(rea);
          
          /* Put state */

          WDB[ptr + REACTION_RFS] = 1.0;

          /* Set pointer to parent reaction */
          
          WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
          
          /* Set type */
          
          WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
          
          /* Set branch mt */
          
          WDB[ptr + REACTION_BRANCH_MT] = RDB[ptr + REACTION_MT];
        }
    }

  /* Check burnup mode */

  if (((long)RDB[DATA_BURNUP_CALCULATION_MODE] == NO) &&
      ((long)RDB[DATA_XENON_EQUILIBRIUM_MODE] == -1) &&
      ((long)RDB[DATA_SAMARIUM_EQUILIBRIUM_MODE] == -1))
    return;
  
  /***************************************************************************/
  
  /***** Branching to secondary products *************************************/
  
  /* Loop over reactions */
  
  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
  while (rea > VALID_PTR)
    {
      /* Check type (noi isomeerisiin tiloihin haarautumiset pitää ottaa */
      /* huomioon molempien reaktioiden sekundäärituotteiden osalta nyt  */
      /* reaktiot jaettiin kokonaan kahteen osaan. */

      if ((((long)RDB[rea + REACTION_TYPE] != REACTION_TYPE_TRA_BRANCH) &&
           ((long)RDB[rea + REACTION_TYPE] != REACTION_TYPE_SPECIAL)) ||
          ((long)RDB[rea + REACTION_PTR_ISO_BRA] > VALID_PTR))
        {
          /* Loop over successive decay modes */
          
          for (n = 0; n < 5; n++)
            {
              /* Avoid compiler warning */
              
              mt = -1;
              
              /* Get mt */
              
              if (n == 0)
                mt = (long)RDB[rea + REACTION_MT];
              else if (n == 1)
                mt = (long)RDB[rea + REACTION_RTYP2];
              else if (n == 2)
                mt = (long)RDB[rea + REACTION_RTYP3];
              else if (n == 3)
                mt = (long)RDB[rea + REACTION_RTYP4];
              else if (n == 4)
                mt = (long)RDB[rea + REACTION_RTYP5];
              else
                Die(FUNCTION_NAME, "Overflow");
              
              /* Break if no secondary modes */
              
              if ((mt == 10000) || (mt == 0))
                break;
              
              /* Put mt for partial modes */
              
              if ((mt > 599) && (mt < 650))
                mt1 = 103;
              else if ((mt > 649) && (mt < 700))
                mt1 = 104;
              else if ((mt > 699) && (mt < 750))
                mt1 = 105;
              else if ((mt > 749) && (mt < 800))
                mt1 = 106;
              else if ((mt > 799) && (mt < 850))
                mt1 = 107;
              else
                mt1 = mt;
              
              /* Neutron reactions */
              
              switch (mt1)
                {
                case 22:
                  {
                    /*(n,na) */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 2004;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;

                    break;
                  }
                case 23:
                  {
                    /* (n, n3a) */
                
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 2004;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BR] = 3.0*RDB[rea + REACTION_BR];
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    break;
                  }
                case 24:
                  {
                    /* (n, 2na) */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 2004;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    break;
                  }
                case 25:
                  {
                    /* (n, 3na) */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 2004;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    break;
                  }
                case 28:
                  {
                    /* (n, np) */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 1001;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    break;
                  }
                case 29:
                  {
                    /* (n, n2a) */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 2004;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BR] = 2.0*RDB[rea + REACTION_BR];
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    break;
                  }
                case 30:
                  {
                    /* (n, 2n2a) */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 2004;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BR] = 2.0*RDB[rea + REACTION_BR];
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;

                    break;
                  }
                case 32:
                  {
                    /* (n, nd) */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 1002;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    break;
                  }
                case 33:
                  {
                    /* (n, nt) */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 1003;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    break;
                  }
                case 34:
                  {
                    /* (n, nHe-3) */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 2003;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    break;
                  }
                case 35:
                  {
                    /* (n, nd2a) */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 1002;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 2004;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BR] = 2.0*RDB[rea + REACTION_BR];
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    break;
                  }
                case 36:
                  {
                    /* (n, nt2a) */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 1003;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;                
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 2004;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BR] = 2.0*RDB[rea + REACTION_BR];
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    break;
                  }
                case 41:
                  {
                    /* (n,2np) */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 1001;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    break;
                  }
                case 42:
                  {
                    /* (n,3np) */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 1001;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    break;
                  }
                case 44:
                  {
                    /* (n,n2p) */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 1001;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BR] = 2.0*RDB[rea + REACTION_BR];
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    break;
                  }
                case 45:
                  {
                    /* (n, npa) */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 1001;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 2004;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    break;
                  }
                case 103:
                  {
                    /* (n,p) */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 1001;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    break;
                  }
                case 104:
                  {
                    /* (n,d) */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 1002;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    break;
                  }
                case 105:
                  {
                    /* (n, t) */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 1003;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    break;
                  }
                case 106:
                  {
                    /* (n, He-3) */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 2003;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    break;
                  }
                case 107:
                  {
                    /* (n, a) */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 2004;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    break;
                  }
                case 108:
                  {
                    /* (n, 2a) */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 2004;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BR] = 2.0*RDB[rea + REACTION_BR];
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    break;
                  }
                case 109:
                  {
                    /* (n, 3a) */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 2004;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BR] = 3.0*RDB[rea + REACTION_BR];
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    break;
                  }
                case 111:
                  {
                    /* (n,2p) */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 1001;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BR] = 2.0*RDB[rea + REACTION_BR];
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    break;
                  }
                case 112:
                  {
                    /* (n, pa) */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 1001;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;                
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 2004;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    break;
                  }
                case 113:
                  {
                    /* (n, t2a) - Special treatment for B-10 */
                    
                    if ((long)RDB[nuc + NUCLIDE_ZAI] != 50100)
                      {
                        ptr = DuplicateItem(rea);            
                        WDB[ptr + REACTION_MT] = 20000 + 1003;
                        WDB[ptr + REACTION_RFS] = 0.0;
                        WDB[ptr + REACTION_TYPE] = 
                          (double)REACTION_TYPE_TRA_BRANCH;
                        WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                        WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                        WDB[ptr + REACTION_RTYP2] = 0.0;
                        WDB[ptr + REACTION_RTYP3] = 0.0;
                        WDB[ptr + REACTION_RTYP4] = 0.0;
                        WDB[ptr + REACTION_RTYP5] = 0.0;
                        WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;

                      }
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 2004;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BR] = 2.0*RDB[rea + REACTION_BR];
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;

                    break;
                  }
                case 114:
                  {
                    /* (n, d2a) */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 1002;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;                
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 2004;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BR] = 2.0*RDB[rea + REACTION_BR];
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    break;
                  }
                case 115:
                  {
                    /* (n, pd) */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 1001;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 1002;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    break;
                  }
                case 116:
                  {
                    /* (n, pt) */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 1001;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 1003;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    break;
                  }
                case 117:
                  {
                    /* (n, da) */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 1002;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 2004;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    WDB[ptr + REACTION_PTR_ISO_BRA] = NULLPTR;
                    
                    break;
                  }
                }
              
              /* Decay reactions */
              
              switch (mt - 10000)
                {
                case 4:
                  {
                    /* Alpha decay */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 2004;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_DEC_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    
                    break;
                  }
                case 7: 
                  {
                    /* Proton emission */
                    
                    ptr = DuplicateItem(rea);            
                    WDB[ptr + REACTION_MT] = 20000 + 1001;
                    WDB[ptr + REACTION_RFS] = 0.0;
                    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_DEC_BRANCH;
                    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
                    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
                    WDB[ptr + REACTION_RTYP2] = 0.0;
                    WDB[ptr + REACTION_RTYP3] = 0.0;
                    WDB[ptr + REACTION_RTYP4] = 0.0;
                    WDB[ptr + REACTION_RTYP5] = 0.0;
                    
                    break;
                  }
                }
            }
        }
      
      /* Next */
      
      rea = NextItem(rea);
    }

  /***************************************************************************/

  /***** Fission yields ******************************************************/

  /* Create duplicates */

  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
  while (rea > VALID_PTR)
    {
      /* Check type (NOTE: spontaneous fissiond shouldn't have branches) */

      if (((long)RDB[rea + REACTION_TYPE] != REACTION_TYPE_TRA_BRANCH) && 
          ((long)RDB[rea + REACTION_TYPE] != REACTION_TYPE_DEC_BRANCH))
        {
          /* Pointer to first yield */
          
          if ((yld = (long)RDB[rea + REACTION_PTR_FISSY]) > VALID_PTR)
            {
              /* Put yield energies */

              WDB[rea + REACTION_FISSY_IE0] = -INFTY;
              WDB[rea + REACTION_FISSY_IE1] = RDB[yld + FISSION_YIELD_E];
              WDB[rea + REACTION_FISSY_IE2] = INFTY;

              /* Check if energy-dependence is included */

              if ((long)RDB[DATA_FISSY_ENE_DEP] == YES)
                yld = NextItem(yld);
              else
                yld = -1;

              /* Loop over remaining */

              while(yld > VALID_PTR)
                {
                  /* Check type */

                  if ((long)RDB[rea + REACTION_TYPE] != REACTION_TYPE_PARTIAL)
                    Die(FUNCTION_NAME, "Invalid reaction type");

                  /* Duplicate block */
                  
                  ptr = DuplicateItem(rea);            

                  /* Set yield pointer */
                  
                  WDB[ptr + REACTION_PTR_FISSY] = (double)yld;

                  /* Put yield energies */

                  WDB[ptr + REACTION_FISSY_IE0] = -INFTY;
                  WDB[ptr + REACTION_FISSY_IE1] = RDB[yld + FISSION_YIELD_E];
                  WDB[ptr + REACTION_FISSY_IE2] = INFTY;

                  /* Set type */
                  
                  WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
                  
                  /* Set pointer to parent */
                  
                  WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;

                  /* Set branch mt */

                  WDB[ptr + REACTION_BRANCH_MT] = RDB[rea + REACTION_MT];

                  /* Reset secondary decay modes */

                  WDB[ptr + REACTION_RTYP2] = 0.0;
                  WDB[ptr + REACTION_RTYP3] = 0.0;
                  WDB[ptr + REACTION_RTYP4] = 0.0;
                  WDB[ptr + REACTION_RTYP5] = 0.0;
                  
                  /* Next yield */
                  
                  yld = NextItem(yld);
                }
            }
        }
      
      /* Next */
      
      rea = NextItem(rea);
    }

  /* Put boundaries (tän voisi tehdä myös tossa edellisessä luupissa) */

  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
  while (rea > VALID_PTR)
    {
      /* Pointer to yield */
          
      if ((yld = (long)RDB[rea + REACTION_PTR_FISSY]) > VALID_PTR)
        {
          /* Check previous yield */

          if ((ptr = PrevItem(yld)) > VALID_PTR)
            WDB[rea + REACTION_FISSY_IE0] = RDB[ptr + FISSION_YIELD_E];
              
          /* Check next yield */

          if ((ptr = NextItem(yld)) > VALID_PTR)
            WDB[rea + REACTION_FISSY_IE2] = RDB[ptr + FISSION_YIELD_E];
        }
      
      /* Next */
      
      rea = NextItem(rea);
    }

  /* Check data */

  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
  while (rea > VALID_PTR)
    {
      /* Pointer to first yield */
      
      if ((yld = (long)RDB[rea + REACTION_PTR_FISSY]) > VALID_PTR)  
        {
          /* Check that energy matches yield energy */

          if (RDB[yld + FISSION_YIELD_E] != RDB[rea + REACTION_FISSY_IE1])
            Die(FUNCTION_NAME, "Mismatch in energy");
          
          /* Check boundaries */

          if (RDB[rea + REACTION_FISSY_IE0] >= RDB[rea + REACTION_FISSY_IE1])
            Die(FUNCTION_NAME, "Error in boundaries");

          if (RDB[rea + REACTION_FISSY_IE1] >= RDB[rea + REACTION_FISSY_IE2])
            Die(FUNCTION_NAME, "Error in boundaries");

          /* Check zeros (joissain jakaumissa on nollaenergioita) */
          /*
          if (RDB[rea + REACTION_FISSY_IE0] == 0.0)
            Die(FUNCTION_NAME, "Zero yield energy");

          if (RDB[rea + REACTION_FISSY_IE1] == 0.0)
            Die(FUNCTION_NAME, "Zero yield energy");

          if (RDB[rea + REACTION_FISSY_IE2] == 0.0)
            Die(FUNCTION_NAME, "Zero yield energy");
          */
        }

      /* Next */
      
      rea = NextItem(rea);
    }

  /***************************************************************************/

  /***** Allocate memory one-group transmutation reactions *******************/

  /* Transmutation and fission branching reactions */

  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
  while (rea > VALID_PTR)
    {
      /* Get mt */

      mt = (long)RDB[rea + REACTION_MT];

      /* Allocate memory transmutation reactions */

      if (((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_PARTIAL) &&
          (((mt > 15) && (mt < 50)) || ((mt > 101) && (mt < 200)) ||
           ((mt > 599) && (mt < 900))))
        AllocValuePair(rea + REACTION_PTR_TRANSMUXS);

      /* Allocate memory for isomeric branch reactions */

      if (((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_TRA_BRANCH) &&
          ((long)RDB[rea + REACTION_PTR_ISO_BRA] > VALID_PTR))
        AllocValuePair(rea + REACTION_PTR_TRANSMUXS);
      
      /* Allocate memory for fission branch reactions */

      if ((((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_TRA_BRANCH) &&
           ((long)RDB[rea + REACTION_PTR_FISSY] > VALID_PTR)))
        AllocValuePair(rea + REACTION_PTR_TRANSMUXS);
      
      /* Next */
      
      rea = NextItem(rea);
    }

  /* Link pointer for other branch reactions */

  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
  while (rea > VALID_PTR)
    {
      /* Check type and allocate memory (NOTE: pointteri PRIVA-blokkiin)*/

      if (((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_TRA_BRANCH) &&
          ((long)RDB[rea + REACTION_PTR_TRANSMUXS] < 1))
        {
          /* Pointer to parent reaction */

          ptr = (long)RDB[rea + REACTION_PTR_BRANCH_PARENT];
          
          /* Check pointer */

          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Link pointer */

          WDB[rea + REACTION_PTR_TRANSMUXS] = 
            RDB[ptr + REACTION_PTR_TRANSMUXS];
        }

      /* Next */
      
      rea = NextItem(rea);
    }

  /* Add reaction structure for total fission */
  
  if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_FISSILE)
    {
      /* Reset Q-value */

      Q = -1.0;

      /* Loop over reactions to get Q-value */
      
      rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
      while (rea > VALID_PTR)
        {
          /* Check mt */

          if (((long)RDB[rea + REACTION_MT] == 18) ||
              ((long)RDB[rea + REACTION_MT] == 19))
            {
              /* Get value */

              Q = RDB[rea + REACTION_Q];
              
              /* Break loop */

              break;
            }
          
          /* Next reaction */
          
          rea = NextItem(rea);
        }

      /* Check Q-value */

      if (Q < ZERO)
        Die(FUNCTION_NAME, "Nuclide %s has no fission channel",
            GetText(nuc + NUCLIDE_PTR_NAME));

      /* Add reaction */

      rea = NewItem(nuc + NUCLIDE_PTR_TOTFISS_REA, REACTION_BLOCK_SIZE);

      /* Set mt and Q-value */

      WDB[rea + REACTION_MT] = -INFTY;
      WDB[rea + REACTION_Q] = Q;

      /* Set nuclide pointer */

      WDB[rea + REACTION_PTR_NUCLIDE] = (double)nuc;

      /* Allocate memory for one-group transmutation xs */
      
      AllocValuePair(rea + REACTION_PTR_TRANSMUXS);
    }

  /***************************************************************************/
}

/*****************************************************************************/
