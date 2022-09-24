/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : defaultbradata.c                               */
/*                                                                           */
/* Created:       2016/09/16 (JLe)                                           */
/* Last modified: 2016/09/24 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: - sets default isomeric branching ratio data                 */
/*                                                                           */
/* Comments: - These values are calculated from energy-dependent branching   */
/*             ratios in the JEFF-3.1 activation file in PWR flux spectrum.  */
/*                                                                           */
/*           - Päivitä tää taulukko Wikiin jos sitä muutetaan                */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DefaultBraData:"

/*****************************************************************************/

void DefaultBraData()
{
  long loc0, loc1, ptr;

  /***************************************************************************/

  /***** Hard-coded values ***************************************************/

  /* These values are calculated from energy-dependent branching  */
  /* ratios in the JEFF-3.1 activation file in PWR flux spectrum. */

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 110230;
  WDB[loc0 + FIX_BRA_FRAC] = 0.2320;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 170370;
  WDB[loc0 + FIX_BRA_FRAC] = 0.8809;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 210450;
  WDB[loc0 + FIX_BRA_FRAC] = 0.5560;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 270590;
  WDB[loc0 + FIX_BRA_FRAC] = 0.4440;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 320720;
  WDB[loc0 + FIX_BRA_FRAC] = 0.5012;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 320740;
  WDB[loc0 + FIX_BRA_FRAC] = 0.6660;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 320760;
  WDB[loc0 + FIX_BRA_FRAC] = 0.4005;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 340760;
  WDB[loc0 + FIX_BRA_FRAC] = 0.7409;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 340780;
  WDB[loc0 + FIX_BRA_FRAC] = 0.1178;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 340800;
  WDB[loc0 + FIX_BRA_FRAC] = 0.8454;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 340820;
  WDB[loc0 + FIX_BRA_FRAC] = 0.1402;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 350790;
  WDB[loc0 + FIX_BRA_FRAC] = 0.7687;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 350810;
  WDB[loc0 + FIX_BRA_FRAC] = 0.0914;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 360780;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9704;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 360800;
  WDB[loc0 + FIX_BRA_FRAC] = 0.6031;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 360820;
  WDB[loc0 + FIX_BRA_FRAC] = 0.3330;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 360840;
  WDB[loc0 + FIX_BRA_FRAC] = 0.1839;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 370850;
  WDB[loc0 + FIX_BRA_FRAC] = 0.8791;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 380840;
  WDB[loc0 + FIX_BRA_FRAC] = 0.2530;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 380860;
  WDB[loc0 + FIX_BRA_FRAC] = 0.1988;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 390890;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9979;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 390900;
  WDB[loc0 + FIX_BRA_FRAC] = 0.7496;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 410930;
  WDB[loc0 + FIX_BRA_FRAC] = 0.3101;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 410940;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9610;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 420920;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9978;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 451030;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9240;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 451050;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9040;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 461060;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9527;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 461080;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9779;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 461100;
  WDB[loc0 + FIX_BRA_FRAC] = 0.8500;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 471070;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9898;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 471090;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9540;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 481100;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9945;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 481120;
  WDB[loc0 + FIX_BRA_FRAC] = 0.8685;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 481140;
  WDB[loc0 + FIX_BRA_FRAC] = 0.8812;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 481160;
  WDB[loc0 + FIX_BRA_FRAC] = 0.6660;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 491130;
  WDB[loc0 + FIX_BRA_FRAC] = 0.4191;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 501120;
  WDB[loc0 + FIX_BRA_FRAC] = 0.7253;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 501160;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9568;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 501180;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9794;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 501200;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9875;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 501220;
  WDB[loc0 + FIX_BRA_FRAC] = 0.0112;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 501240;
  WDB[loc0 + FIX_BRA_FRAC] = 0.0375;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 501260;
  WDB[loc0 + FIX_BRA_FRAC] = 0.3018;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 511210;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9369;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 521200;
  WDB[loc0 + FIX_BRA_FRAC] = 0.8871;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 521220;
  WDB[loc0 + FIX_BRA_FRAC] = 0.6448;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 521240;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9912;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 521260;
  WDB[loc0 + FIX_BRA_FRAC] = 0.8689;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 521280;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9245;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 521300;
  WDB[loc0 + FIX_BRA_FRAC] = 0.8559;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 521320;
  WDB[loc0 + FIX_BRA_FRAC] = 0.8517;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 531290;
  WDB[loc0 + FIX_BRA_FRAC] = 0.4130;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 531310;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9839;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 541240;
  WDB[loc0 + FIX_BRA_FRAC] = 0.8300;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 541260;
  WDB[loc0 + FIX_BRA_FRAC] = 0.8691;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 541280;
  WDB[loc0 + FIX_BRA_FRAC] = 0.8923;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 541300;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9164;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 541320;
  WDB[loc0 + FIX_BRA_FRAC] = 0.8867;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 541330;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9600;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 541340;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9853;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 551330;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9070;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 551340;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9960;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 551350;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9840;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 551370;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9021;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 561300;
  WDB[loc0 + FIX_BRA_FRAC] = 0.8871;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 561320;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9175;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 561340;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9263;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 561350;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9978;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 561360;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9731;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 581360;
  WDB[loc0 + FIX_BRA_FRAC] = 0.8662;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 581380;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9787;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 591410;
  WDB[loc0 + FIX_BRA_FRAC] = 0.6519;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 591430;
  WDB[loc0 + FIX_BRA_FRAC] = 0.3100;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 611470;
  WDB[loc0 + FIX_BRA_FRAC] = 0.5330;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 631530;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9840;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 661640;
  WDB[loc0 + FIX_BRA_FRAC] = 0.3700;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 671650;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9490;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 681660;
  WDB[loc0 + FIX_BRA_FRAC] = 0.2503;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 711750;
  WDB[loc0 + FIX_BRA_FRAC] = 0.3331;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 711760;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9990;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 721790;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9910;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 741820;
  WDB[loc0 + FIX_BRA_FRAC] = 0.8699;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 741840;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9983;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 751850;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9990;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 751870;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9729;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 791970;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9990;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 801960;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9660;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 801980;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9918;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 822060;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9783;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 832090;
  WDB[loc0 + FIX_BRA_FRAC] = 0.6791;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 912330;
  WDB[loc0 + FIX_BRA_FRAC] = 0.4871;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 922340;
  WDB[loc0 + FIX_BRA_FRAC] = 0.5000;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 932350;
  WDB[loc0 + FIX_BRA_FRAC] = 0.4000;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 932390;
  WDB[loc0 + FIX_BRA_FRAC] = 0.3573;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 942360;
  WDB[loc0 + FIX_BRA_FRAC] = 0.5001;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 952410;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9190;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 952430;
  WDB[loc0 + FIX_BRA_FRAC] = 0.0626;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 972470;
  WDB[loc0 + FIX_BRA_FRAC] = 0.4000;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 992530;
  WDB[loc0 + FIX_BRA_FRAC] = 0.0320;
  WDB[loc0 + FIX_BRA_MT] =  102;

  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);
  WDB[loc0 + FIX_BRA_ZAI] = 992550;
  WDB[loc0 + FIX_BRA_FRAC] = 0.9840;
  WDB[loc0 + FIX_BRA_MT] =  102;

  /***************************************************************************/

  /***** Remove duplicates ***************************************************/

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_FIX_BRA0];
  while (loc0 > VALID_PTR)
    {
      /*
      printf("|%s\n", ZAItoIso((long)RDB[loc0 + FIX_BRA_ZAI], 1));
      printf("|%ld\n", (long)RDB[loc0 + FIX_BRA_ZAI]);
      printf("|%ld\n", (long)RDB[loc0 + FIX_BRA_MT]);
      printf("|%1.5f\n", RDB[loc0 + FIX_BRA_FRAC]);
      printf("|%1.5f\n", 1.0 - RDB[loc0 + FIX_BRA_FRAC]);
      printf("|-\n");
      */
      /* Loop over remaining */

      loc1 = NextItem(loc0);
      while (loc1 > VALID_PTR)
        {
          /* Compare ZAI and mt */

          if (((long)RDB[loc0 + FIX_BRA_ZAI] == 
               (long)RDB[loc1 + FIX_BRA_ZAI]) && 
              ((long)RDB[loc0 + FIX_BRA_MT] == 
               (long)RDB[loc1 + FIX_BRA_MT]))
            {
              /* Get pointer */

              ptr = loc1;

              /* Pointer to next */
          
              loc1 = NextItem(loc1);

              /* Remove duplicate */

              RemoveItem(ptr);
            }
          else
            {
              /* Next */
          
              loc1 = NextItem(loc1);
            }
        }

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/
}

