/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processinventory.c                             */
/*                                                                           */
/* Created:       2011/05/24 (JLe)                                           */
/* Last modified: 2018/06/08 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Puts inventory flags in nuclides, etc.                       */
/*                                                                           */
/* Comments: - Transport-dataa ei lueta inventory-listasta samaan tapaan     */
/*             kuin ykkösessä.                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessInventory:"

/* Local function prototype */

/*****************************************************************************/

void ProcessInventory()
{
  long nuc, loc0, ZAI, i, n, m, ptr, loc1, nf, Z, A, I;
  char name[MAX_STR], *str;

  /* Nuclides in pre-defined lists */

  long ZAI0[][150] = {

    /* Nuclides with significant health impact & high migration probability */
    /* in accident conditions                                               */
    /* Source: http://www.oecd-nea.org/rp/chernobyl/c02.html                */
    /* (Chernobyl: Assessment of Radiological and Health Impact)            */

    { 541330, 531310, 551340, 551370, 521320, 380890, 
      380900, 561400, 400950, 420990, 441030, 441060, 
      581410, 581440, 932390, 942380, 942390, 942400, 
      942410, 962420, -1 },

    /* All actinides for which cross section data are found in JEFF-3.1.1   */
    
    { 892250, 892260, 892270, 902270, 902280, 902290, 
      902300, 902320, 902330, 902340, 912310, 912320,
      912330, 922320, 922330, 922340, 922350, 922360,
      922370, 922380, 932350, 932360, 932370, 932380, 
      932390, 942360, 942370, 942380, 942390, 942400, 
      942410, 942420, 942430, 942440, 942460, 952410, 
      952420, 952421, 952430, 952440, 952441, 962400,
      962410, 962420, 962430, 962440, 962450, 962460, 
      962470, 962480, 962490, 962500, 972470, 972490, 
      972500, 982490, 982500, 982510, 982520, 992530,
      992540, 992550, 1002550, -1},

    /* Nuclides commonly considered in burnup credit criticality analyses   */
    /* for PWR fuels                                                        */
    /* (source: www.oecd-nea.org/science/wpncs/ADSNF/SOAR_final.pdf         */

    { 922340, 922350, 922360, 922380, 942380, 942390,
      942400, 942410, 942420, 932370, 952410, 952430,
      962430, 962440, 962450, 551330, 601430, 601350,
      621470, 621490, 621500, 621510, 621520, 631530,
      641550, 420950, 430990, 441010, 451030, 471090,
      481130, -1 },

    /* Burnup indicators (commonly measured from spent fuel)                */
    /* (source: www.oecd-nea.org/science/wpncs/ADSNF/SOAR_final.pdf         */

    { 551370, 571390, 601480, -1 },

    /* Inventory list used by the COSI6 code (other fission products are    */
    /* lumped)                                                              */
    
    { 922320, 922330, 922340, 922350, 922360, 922380, 
      942360, 942380, 942390, 942400, 942410, 942420, 
      952410, 952421, 952430, 932370, 932390, 962420, 
      962430, 962440, 962450, 962460, 962470, 962480, 
      902270, 902280, 902290, 902300, 902310, 902320, 
      912310, 912330, 551330, 551340, 551350, 551360, 
      551370, 531290, 461040, 461050, 461060, 461070, 
      461080, 461090, 461100, 461120, 340760, 340770, 
      340780, 340790, 340800, 340820, 501150, 501160, 
      501170, 501180, 501190, 501200, 501210, 501220, 
      501230, 501240, 501250, 501260, 430990, 400900, 
      400910, 400920, 400930, 400940, 400950, 400960, 
      400970, -1 },
    
    /* All lanthanides for which cross sections data are found in           */
    /* JEFF-3.1.1                                                           */
    
    { 571380, 571390, 571400, 581400, 581410, 581420, 
      581430, 581440, 591410, 591420, 591430, 601420, 
      601430, 601440, 601450, 601460, 601470, 601480, 
      601500, 611470, 611480, 611481, 611490, 611510, 
      621440, 621470, 621480, 621490, 621500, 621510, 
      621520, 621530, 621540, 631510, 631520, 631530, 
      631540, 631550, 631560, 631570, 641520, 641540, 
      641550, 641560, 641570, 641580, 641600, 651590, 
      651600, 661600, 661610, 661620, 661630, 661640, 
      671650, 681640, 681660, 681670, 661680, 661700, 
      711750, 711760, -1 },

    /* Relevant radionuclides in long-term waste analyses                   */
    /* Source: Nirex Report, "The Identification of Radionuclides           */
    /* Relevant to Long-Term Waste Management in the United Kingdom",       */
    /* Nirex Report no. N/105, (2004)                                       */

    { 10030, 40100, 60140, 170360, 180390, 180420,
      190400, 200410, 250530, 250540, 260550, 270600,
      280590, 280630, 300650, 340790, 360810, 360850,
      370870, 380900, 400930, 410910, 410920, 410931,
      410940, 420930, 430970, 430990, 441060, 461070,
      471081, 471101, 481090, 481131, 501191, 501211,
      501230, 501250, 501260, 521251, 521271, 531290,
      551340, 551350, 551370, 561330, 571370, 571380,
      581440, 611450, 611470, 621470, 621510, 631520,
      631540, 631550, 641530, 671630, 671661, 691700,
      691710, 711740, 711760, 721781, 721820, 781930,
      812040, 822050, 822100, 832080, 832101, 842100,
      882230, 882250, 882260, 882280, 892270, 902270,
      902280, 902290, 902300, 902320, 902340, 912310, 
      912330, 922320, 922330, 922340, 922350, 922360, 
      922380, 932370, 942360, 942380, 942390, 942400, 
      942410, 942420, 952410, 952421, 952430, 962420, 
      952430, 952440, 952450, 952460, 952480, 962490, 
      962500, 962510, 962520, -1 },
    
    /* All minor actinides for which xs data is found in JEFF-3.1.1         */
    /* (Here MA = actinides - thorium - uranium - plutonium)                */

    { 892250, 892260, 892270, 912310, 912320, 912330, 
      932350, 932360, 932370, 932380, 932390, 952410, 
      952420, 952421, 952430, 952440, 952441, 962400,
      962410, 962420, 962430, 962440, 962450, 962460, 
      962470, 962480, 962490, 962500, 972470, 972490, 
      972500, 982490, 982500, 982510, 982520, 992530,
      992540, 992550, 1002550, -1},

    /* List of actinide decay products taken from Serpent 1                 */

    { 802060, 812060, 812060, 812070, 812080, 812090, 
      812100, 822060, 822060, 822070, 822070, 822080, 
      822080, 822090, 822090, 822100, 822100, 822110, 
      822120, 822140, 832090, 832100, 832110, 832120,  
      832130, 832140, 832140, 842100, 842110, 842120, 
      842130, 842130, 842140, 842140, 842150, 842160, 
      842180, 852170, 852180, 862170, 862180, 862190, 
      862200, 862220, 872210, 872230, 882230, 882230, 
      882240, 882250, 882260, 882270, 882280, 892250,
      892270, 892280, -1},

    /* List terminator */

    {-1} 
  };
  
  /* Get pointer */

  if ((loc0 = (long)RDB[DATA_BURN_PTR_INVENTORY]) < VALID_PTR)
    return;

  /***************************************************************************/
  
  /***** Check for all-option ************************************************/

  /* Loop over list */

  loc0 = (long)RDB[DATA_BURN_PTR_INVENTORY];
  while (loc0 > VALID_PTR)
    {
      /* Check option */

      if (!strcasecmp(GetText(loc0 + INVENTORY_PTR_NAME), "all"))
        {
          fprintf(outp, "Adding all nuclides in inventory list...\n");
          
          /* Reset pointer */
          
          WDB[DATA_BURN_PTR_INVENTORY] = NULLPTR;
          
          /* Reset index */
          
          i = 0;
          
          /* Loop over nuclides */

          nuc = (long)RDB[DATA_PTR_NUC0];
          while (nuc > VALID_PTR)
            {
              /* Check depleted flag (added 6.11.2017 / 2.1.30 / JLe) */
              
              if (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DEP))
                {
                  /* Pointer to next */

                  nuc = NextItem(nuc);

                  /* Cycle loop */

                  continue;
                }

              /* Check existing */

              loc1 = (long)RDB[DATA_BURN_PTR_INVENTORY];
              while (loc1 > VALID_PTR)
                {
                  /* Compare ZAI */

                  if (!strcmp(GetText(loc1 + INVENTORY_PTR_NAME),
                              ZAItoIso((long)RDB[nuc + NUCLIDE_ZAI], 3)))
                    break;

                  /* Next */

                  loc1 = NextItem(loc1);
                }

              /* Check if found */

              if (loc1 < VALID_PTR)
                {
                  /* Add nuclide to inventory */
                  
                  loc0 = NewItem(DATA_BURN_PTR_INVENTORY, 
                                 INVENTORY_BLOCK_SIZE);
              
                  /* Put name */
                  
                  sprintf(name, "%s", ZAItoIso((long)RDB[nuc + NUCLIDE_ZAI], 
                                               3));
                  WDB[loc0 + INVENTORY_PTR_NAME] = (double)PutText(name);
                  
                  /* Put ZAI */
                  
                  WDB[loc0 + INVENTORY_ZAI] = RDB[nuc + NUCLIDE_ZAI];
              
                  /* Set index */
              
                  WDB[nuc + NUCLIDE_INVENTORY_IDX] = (double)i;
                  
                  /* Add counter */
                  
                  i++;
                }
              
              /* Next nuclide */
              
              nuc = NextItem(nuc);
            }
          
          fprintf(outp, "OK.\n\n");
          
          /* Put number of nuclides */
          
          WDB[DATA_BURN_INVENTORY_NUCLIDES] = (double)i;
          
          /* Exit subroutine */
          
          return;
        }

      /* Next item */
          
      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Check for fission products ******************************************/

  /* Loop over list */

  loc0 = (long)RDB[DATA_BURN_PTR_INVENTORY];
  while (loc0 > VALID_PTR)
    {
      /* Get name */
      
      sprintf(name, "%s", GetText(loc0 + INVENTORY_PTR_NAME));

      /* Check key word */

      if (!strcasecmp(name, "fp"))
        {
          /* Loop over nuclides */

          nuc = (long)RDB[DATA_PTR_NUC0];
          while (nuc > VALID_PTR)
            {
              /* Check flags */
              
              if (((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DEP) &&
                  ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_FP))
                {
                  /* Create new entry */
                    
                  loc1 = NewItem(DATA_BURN_PTR_INVENTORY, 
                                 INVENTORY_BLOCK_SIZE);

                  /* Put name */
                  
                  sprintf(name, "%s", ZAItoIso((long)RDB[nuc + NUCLIDE_ZAI], 
                                               3));
                  WDB[loc1 + INVENTORY_PTR_NAME] = (double)PutText(name);
                  
                  /* Put ZAI */
                  
                  WDB[loc1 + INVENTORY_ZAI] = RDB[nuc + NUCLIDE_ZAI];
                }

              /* Next nuclide */

              nuc = NextItem(nuc);
            }

          /* Skip entry */

          loc1 = loc0;
          loc0 = NextItem(loc0);

          /* Remove from list */

          RemoveItem(loc1);
        }
      else
        {
          /* Next item */
          
          loc0 = NextItem(loc0);
        }
    }
  
  /***************************************************************************/

  /***** Check for noble gases ***********************************************/

  /* Loop over list */

  loc0 = (long)RDB[DATA_BURN_PTR_INVENTORY];
  while (loc0 > VALID_PTR)
    {
      /* Get name */
      
      sprintf(name, "%s", GetText(loc0 + INVENTORY_PTR_NAME));

      /* Check key word */

      if (!strcasecmp(name, "ng"))
        {
          /* Loop over nuclides */

          nuc = (long)RDB[DATA_PTR_NUC0];
          while (nuc > VALID_PTR)
            {
              /* Check flags */
              
              if (((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DEP) &&
                  (((long)RDB[nuc + NUCLIDE_Z] == 10) ||
                   ((long)RDB[nuc + NUCLIDE_Z] == 18) ||
                   ((long)RDB[nuc + NUCLIDE_Z] == 36) ||
                   ((long)RDB[nuc + NUCLIDE_Z] == 54)))
                {
                  /* Create new entry */
                    
                  loc1 = NewItem(DATA_BURN_PTR_INVENTORY, 
                                 INVENTORY_BLOCK_SIZE);

                  /* Put name */
                  
                  sprintf(name, "%s", ZAItoIso((long)RDB[nuc + NUCLIDE_ZAI], 
                                               3));
                  WDB[loc1 + INVENTORY_PTR_NAME] = (double)PutText(name);
                  
                  /* Put ZAI */
                  
                  WDB[loc1 + INVENTORY_ZAI] = RDB[nuc + NUCLIDE_ZAI];
                }

              /* Next nuclide */

              nuc = NextItem(nuc);
            }

          /* Skip entry */

          loc1 = loc0;
          loc0 = NextItem(loc0);

          /* Remove from list */

          RemoveItem(loc1);
        }
      else
        {
          /* Next item */
          
          loc0 = NextItem(loc0);
        }
    }
  
  /***************************************************************************/
  
  /***** Add nuclides from pre-defined lists ********************************/
  
  /* Loop over list */
  
  loc0 = (long)RDB[DATA_BURN_PTR_INVENTORY];
  while (loc0 > VALID_PTR)
    {
      /* Get name */

      sprintf(name, "%s", GetText(loc0 + INVENTORY_PTR_NAME));

      /* Check key word */

      if (!strcasecmp(name, "accident"))
        m = 0;
      else if ((!strcasecmp(name, "actinides")) ||
               (!strcasecmp(name, "act")))
        m = 1;
      else if (!strcasecmp(name, "burnupcredit"))
        m = 2;
      else if (!strcasecmp(name, "burnupindicators"))
        m = 3;
      else if (!strcasecmp(name, "cosi6"))
        m = 4;
      else if (!strcasecmp(name, "lanthanides"))
        m = 5;
      else if (!strcasecmp(name, "longterm"))
        m = 6;
      else if (!strcasecmp(name, "minoractinides"))
        m = 7;
      else if (!strcasecmp(name, "dp"))
        m = 8;
      else
        m = -1;

      /* Check index */

      if (m > -1)
        {
          /* Loop over list */

          n = 0;
          while ((ZAI = ZAI0[m][n++]) > 0)
            {
              /* Create new entry */
                    
              loc1 = NewItem(DATA_BURN_PTR_INVENTORY, INVENTORY_BLOCK_SIZE);
              WDB[loc1 + INVENTORY_PTR_NAME] 
                = (double)PutText(ZAItoIso(ZAI, 1));
            }

          /* Skip entry */

          loc1 = loc0;
          loc0 = NextItem(loc0);

          /* Remove from list */

          RemoveItem(loc1);
        }
      else
        {
          /* Next item */
          
          loc0 = NextItem(loc0);
        }
    }
  
  /***************************************************************************/
  
  /***** Process list ********************************************************/

  fprintf(outp, "Processing nuclide inventory list...\n");

  /* Put ZAI's */

  loc0 = (long)RDB[DATA_BURN_PTR_INVENTORY];
  while (loc0 > VALID_PTR)
    {
      /* Pointer to text */

      str = &ASCII[(long)RDB[loc0 + INVENTORY_PTR_NAME]];

      /* Remove dashes */

      n = 0;
      m = 0;

      while (str[n] != '\0')
        {
          if (str[n] != '-')
            str[m++] = str[n];

          n++;
        }
      
      str[m] = '\0';

      /* Get name */

      sprintf(name, "%s", GetText(loc0 + INVENTORY_PTR_NAME));

      /* Get ZA */

     if ((ZAI = IsotoZAI(name)) < 0)
        ZAI = atoi(name);
      else if (ZAI < 1)
        Error(0, "Error in inventory list: \"%s\" is not a valid entry", name);

      /* Check value */

      if (ZAI < 1000)
        {
          /* Elemental */

          Z = ZAI;
          A = 0;
          I = 0;
        }
    else
        {
          /* Decompose ZAI */

          Z = (long)((double)ZAI/10000.0);
          A = (long)((double)ZAI/10.0) - 1000*Z;
          I = ZAI - 10000*Z - 10*A;
        }

      /* Check */

      if ((Z < 1) || (Z > 111) || (A < 0) || (A > 300) || (I < 0) || (I > 10))
        Error(0, "Invalid ZAI \"%s\" in inventory list (interpreted as %s)",
              name, ZAItoIso(ZAI,1));

      /* Put ZAI */

      WDB[loc0 + INVENTORY_ZAI] = (double)ZAI;

      /* Put name */

      WDB[loc0 + INVENTORY_PTR_NAME] = (double)PutText(ZAItoIso(ZAI, 3));

      /* Next item */
      
      loc0 = NextItem(loc0);
    }

  /* Remove duplicates */

  do
    {
      /* Reset count */

      i = 0;

      /* First loop */

      loc0 = (long)RDB[DATA_BURN_PTR_INVENTORY];
      while (loc0 > VALID_PTR)
        {
          /* Loop over remaining */
          
          loc1 = NextItem(loc0);
          while (loc1 > VALID_PTR)
            {
              /* Compare ZAI */
              
              if (RDB[loc0 + INVENTORY_ZAI] == RDB[loc1 + INVENTORY_ZAI])
                {
                  /* Remove second */

                  RemoveItem(loc1);

                  /* Add to count */

                  i++;

                  /* Break loop */

                  break;                  
                }

              /* Next item */
          
              loc1 = NextItem(loc1);
            }
          
          /* Next item */
      
          loc0 = NextItem(loc0);
        }
    }
  while (i > 0);

  /* Reset index and number of nuclides not found */

  i = 0;
  nf = 0;
  
  /* Loop over list */

  loc0 = (long)RDB[DATA_BURN_PTR_INVENTORY];
  while (loc0 > VALID_PTR)
    {
      /* Get ZAI */

      ZAI = (long)RDB[loc0 + INVENTORY_ZAI];

      /* Reset count */

      n = 0;

      /* Find nuclides */

      nuc = (long)RDB[DATA_PTR_NUC0];
      while (nuc > VALID_PTR)
        {
          /* Compare ZAI and Z */
          
          if (((long)RDB[nuc + NUCLIDE_ZAI] == ZAI) ||
              ((long)RDB[nuc + NUCLIDE_Z] == ZAI))
            {
              /* Set index */

              WDB[nuc + NUCLIDE_INVENTORY_IDX] = (double)i;

              /* Add counter */

              n++;
            }

          /* Next nuclide */

          nuc = NextItem(nuc);
        }

      /* Check count */

      if ((n == 0) && ((long)RDB[DATA_PTR_NUC0] > VALID_PTR))
        {
          /* Not found */

          if (nf++ == 0)
            {
              fprintf(outp, "\nThe following nuclides were not found in ");
              fprintf(outp, "compositions and\nwere removed from the list:\n\n");
            }

          if ((long)RDB[loc0 + INVENTORY_PTR_ENTRY] < VALID_PTR)
            fprintf(outp, " %ld (%s)\n", (long)RDB[loc0 + INVENTORY_ZAI],
                    ZAItoIso((long)RDB[loc0 + INVENTORY_ZAI],1));
          else if (atol(GetText(loc0 + INVENTORY_PTR_ENTRY)) == 
              (long)RDB[loc0 + INVENTORY_ZAI])
            fprintf(outp, " %s (%s)\n", GetText(loc0 + INVENTORY_PTR_ENTRY),
                    ZAItoIso((long)RDB[loc0 + INVENTORY_ZAI],1));
          else
            fprintf(outp, " %s\n", GetText(loc0 + INVENTORY_PTR_ENTRY));

          /* Copy pointer */

          ptr = loc0;

          /* Next item */

          loc0 = NextItem(loc0);

          /* Remove nuclide from list */

          RemoveItem(ptr);
        }
      else
        {
          /* Update index */

          i++;

          /* Next item */

          loc0 = NextItem(loc0);
        }
    }

  /* Put number of nuclides */

  WDB[DATA_BURN_INVENTORY_NUCLIDES] = (double)i;

  /***************************************************************************/

  /* Done */
  
  if (nf == 0)
    fprintf(outp, "OK.\n\n");
  else
    fprintf(outp, "\n");
}

/*****************************************************************************/
