/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : reactionmt.c                                   */
/*                                                                           */
/* Created:       2010/09/12 (JLe)                                           */
/* Last modified: 2017/08/23 (JLe)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Returns the name of transmutation or decay reaction          */
/*                                                                           */
/* Comments: - MT for decay modes is RTYP*10000                              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "ReactionMT:"

/*****************************************************************************/

char *ReactionMT(long mt, long dosi)
{
  long n;
  static char returnval[MAX_STR];

  /* Reset */
  
  sprintf(returnval, "unknown");

  /* Check mt and possible dosimetry flag */

  if ((mt > 1000) && (dosi == YES))
    {
      n = (long)((double)mt/1000.0);
      mt = mt - n*1000;
    }
  else
    n = 0;

  /* Neutron reactions */

  switch (mt)
    {
    case 0:
      {
        sprintf(returnval, "majorant");
        break;
      }
    case 4:
      {
        sprintf(returnval, "total inelastic scattering");
        break;
      }
    case 5:
      {
        sprintf(returnval, "anything");
        break;
      }
    case 6:
    case 7:
    case 8:
    case 9:
      {
        sprintf(returnval, "Be-9 (n,2n) in ENDF version 5");
        break;
      }
    case 1:
    case -1:
      {
        sprintf(returnval, "total");
        break;
      }
    case -2:
    case 101:
      {
        sprintf(returnval, "total absorption");
        break;
      }
    case 2:
    case -3:
      {
        sprintf(returnval, "elastic scattering");
        break;
      }
    case 3:
      {
        sprintf(returnval, "nonelastic");
        break;
      }
    case -4:
      {
        sprintf(returnval, "average heating number");
        break;
      }
    case 11:
      {
        sprintf(returnval, "(n,2nd)");
        break;
      }
    case 16:
      {
        sprintf(returnval, "(n,2n)");
        break;
      }
    case 17:
      {
        sprintf(returnval, "(n,3n)");
        break;
      }
    case -6:
    case 18:
      {
        sprintf(returnval, "total fission");
        break;
      }
    case 19:
      {
        sprintf(returnval, "(n,f)");
        break;
      }
    case 20:
      {
        sprintf(returnval, "(n,nf)");
        break;
      }
    case 21:
      {
        sprintf(returnval, "(n,2nf)");
        break;
      }
    case 22:
      {
        sprintf(returnval, "(n,nalpha)");
        break;
      }
    case 23:
      {
        sprintf(returnval, "(n,n3alpha)");
        break;
      }
    case 24:
      {
        sprintf(returnval, "(n,2nalpha)");
        break;
      }
    case 25:
      {
        sprintf(returnval, "(n,3nalpha)");
        break;
      }
    case 27:
      {
        sprintf(returnval, "absorption");
        break;
      }
    case 28:
      {
        sprintf(returnval, "(n,np)");
        break;
      }
    case 29:
      {
        sprintf(returnval, "(n,n2alpha)");
        break;
      }
    case 30:
      {
        sprintf(returnval, "(n,2n2alpha)");
        break;
      }
    case 32:
      {
        sprintf(returnval, "(n,nd)");
        break;
      }
    case 33:
      {
        sprintf(returnval, "(n,nt)");
        break;
      }
    case 34:
      {
        sprintf(returnval, "(n,nHe-3)");
        break;
      }
    case 35:
      {
        sprintf(returnval, "(n,nd2alpha)");
        break;
      }
    case 36:
      {
        sprintf(returnval, "(n,nt2alpha)");
        break;
      }
    case 37:
      {
        sprintf(returnval, "(n,4n)");
        break;
      }
    case 38:
      {
        sprintf(returnval, "(n,3nf)");
        break;
      }
    case 41:
      {
        sprintf(returnval, "(n,2np)");
        break;
      }
    case 42:
      {
        sprintf(returnval, "(n,3np)");
        break;
      }
    case 44:
      {
        sprintf(returnval, "(n,n2p)");
        break;
      }
    case 45:
      {
        sprintf(returnval, "(n,npalpha)");
        break;
      }
    case 91:
      {
        sprintf(returnval, "inelastic scattering to continuum");
        break;
      }
    case 102:
      {
        sprintf(returnval, "(n,gamma)");
        break;
      }
    case 103:
      {
        sprintf(returnval, "(n,p)");
        break;
      }
    case 104:
      {
        sprintf(returnval, "(n,d)");
        break;
      }
    case 105:
      {
        sprintf(returnval, "(n,t)");
        break;
      }
    case 106:
      {
        sprintf(returnval, "(n,He-3)");
        break;
      }
    case 107:
      {
        sprintf(returnval, "(n,alpha)");
        break;
      }
    case 108:
      {
        sprintf(returnval, "(n,2alpha)");
        break;
      }
    case 109:
      {
        sprintf(returnval, "(n,3alpha)");
        break;
      }
    case 111:
      {
        sprintf(returnval, "(n,2p)");
        break;
      }
    case 112:
      {
        sprintf(returnval, "(n,palpha)");
        break;
      }
    case 113:
      {
        sprintf(returnval, "(n,t2alpha)");
        break;
      }
    case 114:
      {
        sprintf(returnval, "(n,d2alpha)");
        break;
      }
    case 115:
      {
        sprintf(returnval, "(n,pd)");
        break;
      }
    case 116:
      {
        sprintf(returnval, "(n,pt)");
        break;
      }
    case 117:
      {
        sprintf(returnval, "(n,dalpha)");
        break;
      }
    case 201:
      {
        sprintf(returnval, "total neutron production");
        break;
      }
    case 202:
      {
        sprintf(returnval, "total photon production");
        break;
      }
    case 203:
      {
        sprintf(returnval, "total proton production");
        break;
      }
    case 204:
      {
        sprintf(returnval, "total deuteron production");
        break;
      }
    case 205:
      {
        sprintf(returnval, "total triton production");
        break;
      }
    case 206:
      {
        sprintf(returnval, "total He-3 production");
        break;
      }
    case 207:
      {
        sprintf(returnval, "total alpha production");
        break;
      }
    case 301:
      {
        sprintf(returnval, "total heat production");
        break;
      }
    case 443:
      {
        sprintf(returnval, "kinematic KERMA");
        break;
      }
    case 444:
      {
        sprintf(returnval, "damage-energy production");
        break;
      }
    case 501:
      {
        sprintf(returnval, "total");
        break;
      }
    case 502:
      {
        sprintf(returnval, "Rayleigh scattering");
        break;
      }
    case 504:
      {
        sprintf(returnval, "Compton scattering");
        break;
      }
    case 516:
      {
        sprintf(returnval, "pair production");
        break;
      }
    case 522:
      {
        sprintf(returnval, "photoelectric effect");
        break;
      }
    case 600:
      {
        sprintf(returnval, "(n,p) to ground state");
        break;
      }
    case 649:
      {
        sprintf(returnval, "(n,p) to continuum");
        break;
      }
    case 650:
      {
        sprintf(returnval, "(n,d) to ground state");
        break;
      }
    case 699:
      {
        sprintf(returnval, "(n,d) to continuum");
        break;
      }
    case 700:
      {
        sprintf(returnval, "(n,t) to ground state");
        break;
      }
    case 749:
      {
        sprintf(returnval, "(n,t) to continuum");
        break;
      }
    case 750:
      {
        sprintf(returnval, "(n,He-3) to ground state");
        break;
      }
    case 799:
      {
        sprintf(returnval, "(n,He-3) to continuum");
        break;
      }
    case 800:
      {
        sprintf(returnval, "(n,alpha) to ground state");
        break;
      }
    case 849:
      {
        sprintf(returnval, "(n,alpha) to continuum");
        break;
      }
    case 875:
      {
        sprintf(returnval, "(n,2n) to ground state");
        break;
      }
    case 891:
      {
        sprintf(returnval, "(n,2n) to continuum");
        break;
      }
    case 1002:
      {
        sprintf(returnval, "S(a,b) elastic scattering");
        break;
      }
    case 1004:
      {
        sprintf(returnval, "S(a,b) inelastic scattering");
        break;
      }
    }

  if ((mt > 50) && (mt < 91))
    sprintf(returnval, "inelastic scatt. to %ld. excited state", mt - 50);

  if ((mt > 600) && (mt < 649))
    sprintf(returnval, "(n,p%ld)", mt - 600);

  if ((mt > 650) && (mt < 699))
    sprintf(returnval, "(n,d%ld)", mt - 650);

  if ((mt > 700) && (mt < 749))
    sprintf(returnval, "(n,t%ld)", mt - 700);

  if ((mt > 750) && (mt < 799))
    sprintf(returnval, "(n,He-3_%ld)", mt - 750);

  if ((mt > 800) && (mt < 849))
    sprintf(returnval, "(n,a%ld)", mt - 800);

  if ((mt > 875) && (mt < 891))
    sprintf(returnval, "(n,2n) to %ld. excited state", mt - 875);

  /* Check state number from dosimetry data */

  if (n == 10)
    sprintf(returnval, "%s to ground state", returnval);
  else if (n > 10)
    sprintf(returnval, "%s to %ld. excited state", returnval, n - 10);

  /* Decay reactions */

  switch (mt - 10000)
    {
    case 1:  
      {
        sprintf(returnval, "beta-");
        break;
      }
    case 3:
      {
        sprintf(returnval, "IT");
        break;
      }
    case 4:
      {
        sprintf(returnval, "alpha");
        break;
      }
    case 7: 
      {
        sprintf(returnval, "proton emission");
        break;
      }
    case 5:
      {
        sprintf(returnval, "neutron emission");
        break;
      }
    case 6:
      {
        sprintf(returnval, "spontaneous fission");
        break;
      }
    case 2:
      {
        sprintf(returnval, "EC/beta+");
        break;
      }
    }

  return returnval;
}

/*****************************************************************************/
