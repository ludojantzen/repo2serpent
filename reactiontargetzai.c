/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : reactiontargetzai.c                            */
/*                                                                           */
/* Created:       2010/09/10 (JLe)                                           */
/* Last modified: 2019/10/23 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Calculates ZAI for daughter nuclide after decay or           */
/*              transmutation reaction. Return ZAI, -1 if neutron-induced    */
/*              fission, -2 if spontaneous fission or 0 if reaction doesn't  */
/*              change the nuclide (used also for natural elements).         */
/*                                                                           */
/* Comments: - Reaktiomoodeja ei oo tarkistettu                              */
/*           - Decay reation MT's are defined as MT = 10000 + RTYP           */
/*           - Secondaries are defined as MT = product ZA + 20000            */
/*           - Tää ei ota nyt huomioon sitä että nuklidi voi "hävitä"        */
/*             kokonaan, esim B10(n,t2a) tuottaa ZAI:n 0. Mutta vaikuttaako  */
/*             se?                                                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReactionTargetZAI:"

/*****************************************************************************/

long ReactionTargetZAI(long rea)
{
  long nuc, ZAI, ZA, Z, A, I, new, diff, mt, n;

  /* Exit if sum or special reaction */

  if (((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_SPECIAL) ||
      ((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_SUM))
    return 0;

  /* Get pointer to nuclide */

  nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

  /* Get ZAI */

  ZAI = (long)RDB[nuc + NUCLIDE_ZAI];

  /* Get ZA */

  ZA = (long)((double)ZAI/10.0);

  /* Separate Z and A */

  Z = (long)((double)ZA/1000.0);
  A = ZA - 1000*Z;

  /* Get I */

  I = ZAI - 10*ZA;

  /* Check natural (branching must be excluded) */

  if (A == 0)
    return 0;

  /* Loop over reactions */

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

      /* Change mt of (n,t2a) to (n,2a) for B-10 */

      if ((10000*Z + 10*A + I == 50100) && (mt == 113))
        mt = 108;

      /* Neutron reactions */

      switch (mt)
        {
        case 4:
          {
            /* Total inelastic scattering */

            break;
          }
        case 16:
          {
            /* (n,2n) */

            A = A - 1;

            break;
          }
        case 17:
          {
            /* (n,3n) */

            A = A - 2;

            break;
          }
        case 18:
        case 19:
        case 20:
        case 21:
        case 38:
          {
            /* (n,f) */

            return FISSION_YIELD_TYPE_NFY;
          }
        case 22:
          {
            /* (n,na) */

            A = A - 4;
            Z = Z - 2;

            break;
          }
        case 23:
          {
            /* (n,n3a) */

            A = A - 12;
            Z = Z - 6;

            break;
          }
        case 24:
          {
            /* (n,2na) */

            A = A - 5;
            Z = Z - 2;

            break;
          }
        case 25:
          {
            /* (n,3na) */

            A = A - 6;
            Z = Z - 2;

            break;
          }
        case 28:
          {
            /* (n,np) */

            A = A - 1;
            Z = Z - 1;

            break;
          }
        case 29:
          {
            /* (n,n2a) */

            A = A - 8;
            Z = Z - 4;

            break;
          }
        case 30:
          {
            /* (n,2n2a) */

            A = A - 9;
            Z = Z - 4;

            break;
          }
        case 32:
          {
            /* (n,nd) */

            A = A - 2;
            Z = Z - 1;

            break;
          }
        case 33:
          {
            /* (n,nt) */

            A = A - 3;
            Z = Z - 1;

            break;
          }
        case 34:
          {
            /* (n,nHe-3) */

            A = A - 3;
            Z = Z - 2;

            break;
          }
        case 35:
          {
            /* (n,nd2a) */

            A = A - 10;
            Z = Z - 5;

            break;
          }
        case 36:
          {
            /* (n,nt2a) */

            A = A - 11;
            Z = Z - 5;

            break;
          }
        case 37:
          {
            /* (n,4n) */

            A = A - 3;

            break;
          }
        case 41:
          {
            /* (n,2np) */

            A = A - 2;
            Z = Z - 1;

            break;
          }
        case 42:
          {
            /* (n,3np) */

            A = A - 3;
            Z = Z - 1;

            break;
          }
        case 44:
          {
            /* (n,n2p) */

            A = A - 2;
            Z = Z - 2;

            break;
          }
        case 45:
          {
            /* (n,npa) */

            A = A - 5;
            Z = Z - 3;

            break;
          }
        case 102:
          {
            /* (n,g) */

            A = A + 1;

            break;
          }
        case 103:
          {
            /* (n,p) */

            Z = Z - 1;

            break;
          }
        case 104:
          {
            /* (n,d) */

            A = A - 1;
            Z = Z - 1;

            break;
          }
        case 105:
          {
            /* (n,t) */

            A = A - 2;
            Z = Z - 1;

            break;
          }
        case 106:
          {
            /* (n,He-3) */

            A = A - 2;
            Z = Z - 2;

            break;
          }
        case 107:
          {
            /* (n,a) */

            A = A - 3;
            Z = Z - 2;

            break;
          }
        case 108:
          {
            /* (n,2a) */

            A = A - 7;
            Z = Z - 4;

            break;
          }
        case 109:
          {
            /* (n,3a) */

            A = A - 11;
            Z = Z - 6;

            break;
          }
        case 111:
          {
            /* (n,2p) */

            A = A - 1;
            Z = Z - 2;

            break;
          }
        case 112:
          {
            /* (n,pa) */

            A = A - 4;
            Z = Z - 3;

            break;
          }
        case 113:
          {
            /* (n,t2a) */

            A = A - 10;
            Z = Z - 5;

            break;
          }
        case 114:
          {
            /* (n,d2a) */

            A = A - 9;
            Z = Z - 5;

            break;
          }
        case 115:
          {
            /* (n,pd) */

            A = A - 2;
            Z = Z - 2;

            break;
          }
        case 116:
          {
            /* (n,pt) */

            A = A - 3;
            Z = Z - 2;

            break;
          }
        case 117:
          {
            /* (n,da) */

            A = A - 5;
            Z = Z - 3;

            break;
          }
        }
      /* Decay reactions */

      switch (mt - 10000)
        {
        case 1:
          {
            /* Beta decay */

            Z = Z + 1;

            break;
          }
        case 2:
          {
            /* Electron capture / positron emission */

            Z = Z - 1;

            break;
          }
        case 3:
          {
            /* Isomeric transition */

            if (I == 0)
              Die(FUNCTION_NAME, "Isomeric transition from ground state");
            else
              I = 0;

            break;
          }
        case 4:
          {
            /* Alpha decay */

            A = A - 4;
            Z = Z - 2;

            break;
          }
        case 5:
          {
            /* Neutron emission */

            A = A - 1;

            break;
          }
        case 6:
          {
            /* Spontaneous fission */

            return FISSION_YIELD_TYPE_SFY;
          }
        case 7:
          {
            /* Proton emission */

            A = A - 1;
            Z = Z - 1;

            break;
          }
        }

      /* Additional branching modes to account for secondaries */

      if (mt > 20000)
        return (10*(mt - 20000));

      /* Partial (n,p) reactions */

      if ((mt > 599) && (mt < 650))
        Z = Z - 1;

      /* Partial (n,d) reactions */

      if ((mt > 649) && (mt < 700))
        {
          A = A - 1;
          Z = Z - 1;
        }

      /* Partial (n,t) reactions */

      if ((mt > 699) && (mt < 750))
        {
          A = A - 2;
          Z = Z - 1;
        }

      /* Partial (n,He-3) reactions */

      if ((mt > 749) && (mt < 800))
        {
          A = A - 2;
          Z = Z - 2;
        }

      /* Partial (n,a) reactions */

      if ((mt > 799) && (mt < 850))
        {
          A = A - 3;
          Z = Z - 2;
        }

      /* Partial (n,2n) reactions */

      if (((mt > 874) && (mt < 892)))
        A = A - 1;
    }

  /* Compare to old (MT 4 is special) */

  if ((long)RDB[rea + REACTION_MT] == 4)
    {
      /* Check that incident and target states are different */

      if ((long)RDB[rea + REACTION_RFS] == (long)RDB[nuc + NUCLIDE_I])
        Die(FUNCTION_NAME, "MT 4 to same state");
    }
  else if (Z == 0)
    {
      /* Product has no protons */

      Die(FUNCTION_NAME, "product Z = 0 (%ld)", ZAI);
    }
  else if (10000*Z + 10*A + I == ZAI)
    return 0;

  /* Get new state */

  I = (long)RDB[rea + REACTION_RFS];

  /* Calculate new ZAI and get old */

  new = 10000*Z + 10*A + I;
  ZAI = (long)RDB[nuc + NUCLIDE_ZAI];

  /* Calculate difference between old and new ZA */

  diff = (long)((double)new/10.0) - (long)((double)ZAI/10.0);

  /* Check isomeric transition */

  if ((long)RDB[rea + REACTION_MT] != 4)
    if ((diff == 0) && (ZAI - new != 1) &&
        (ZAI - new != 2) && (ZAI - new != 3))
      Die(FUNCTION_NAME, "Error in isomeric transition %ld --> %ld", ZAI, new);

  /* Check other values */

  if ((diff == 1) ||     /* (n,gamma) */
      (diff == -1) ||    /* (n,2n) */
      (diff == -2) ||    /* (n,3n) */
      (diff == -3) ||    /* (n,4n) */
      (diff == 999) ||   /* Beta- + neutron emission */
      (diff == -1000) || /* (n,p) */
      (diff == 1000) ||  /* Beta- */
      (diff == -1001) || /* (n,d) */
      (diff == -1002) || /* (n,t) */
      (diff == -1004) || /* Beta- + alpha decay */
      (diff == -2003) || /* (n,alpha) */
      (diff == -2004) || /* (n,nalpha) */
      (diff == -2005) || /* (n,2nalpha) */
      (diff == -2006) || /* (n,3nalpha) */
      (diff == -6012) || /* (n,n3alpha) */
      (diff == -4007) || /* B-10 (n,t2alpha) */
      (diff == -2002) || /* (n,n2p) */
      (diff == -3005) || /* (n,npalpha) */
      (diff == -3004) || /* (n,palpha) */
      (diff == -1003) || /* (n,nt) */
      (diff == -2001) || /* (n,2p) */
      (diff == -4008) || /* (n,n2alpha) */
      (diff == -4009) || /* (n,2n2alpha) */
      (diff == -5009) || /* (n,d2alpha) */
      (diff == -5010) || /* (n,t2alpha) */
      (diff == -6011) || /* (n,3alpha) */
      (diff == 0) ||     /* IT */
      (diff == -2000) || /* EC + EC */
      (diff == 2000) ||  /* Beta- + Beta- */
      (diff == 998) ||   /* Beta- + 2x neutron emission */
      (diff == 998) ||   /* Beta- + 2x neutron emission */
      (diff == 997) ||   /* Beta- + 3x neutron emission */
      (diff == 996))     /* Beta- + 4x neutron emission */
    return new;
  else
    Warn(FUNCTION_NAME, "%ld %ld --> %ld (%ld %ld %ld %ld %ld) %s", diff, ZAI,
         new, (long)RDB[rea + REACTION_MT], (long)RDB[rea + REACTION_RTYP2],
         (long)RDB[rea + REACTION_RTYP3], (long)RDB[rea + REACTION_RTYP4],
         (long)RDB[rea + REACTION_RTYP5],
         ReactionMT((long)RDB[rea + REACTION_MT], NO));

  /* Return new ZAI */

  return new;
}

/*****************************************************************************/
