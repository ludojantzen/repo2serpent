/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : stdcomp.c                                      */
/*                                                                           */
/* Created:       2012/05/19 (JLe)                                           */
/* Last modified: 2019/10/26 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Prints standard material compositions that can be            */
/*              copy-pasted in the input.                                    */
/*                                                                           */
/* Comments: - Material compositions from:                                   */
/*                                                                           */
/*             R. J. McConn, Jr. et al. "Compendium of Material Composition  */
/*             Data for Radiation Transport Modeling." Pacific Northwest     */
/*             National Laboratory, PIET-43741-TM-963, PNNL-15870 Rev. 1.    */
/*             March 4, 2011.                                                */
/*                                                                           */
/*             S. Penttilän VTT-raportti (nro?)                              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"
#include "natural_elements.h"

#define FUNCTION_NAME "StdComp:"

#define MAX_STD_COMP 360

/*****************************************************************************/

void StdComp(char *mat, char *id)
{
  long n, m, i, itot, A;
  double dens[MAX_STD_COMP], comp[MAX_STD_COMP][112], sum[112];
  char comment[MAX_STD_COMP][MAX_STR], src[MAX_STD_COMP][MAX_STR];
  char alias[MAX_STD_COMP][MAX_STR], name[MAX_STR];

  /***************************************************************************/

  /***** Put compositions ****************************************************/

  /* Reset arrays */

  for (i = 0; i < MAX_STD_COMP; i++)
    {
      /* Reset density */

      dens[i] = 0.0;

      /* Reset composition */

      for (n = 0; n < 112; n++)
        comp[i][n] = 0.0;
    }

  /* Reset index */

  i = 0;

  sprintf(alias[i], " ");
  sprintf(comment[i], "A-150 Tissue-Equivalent Plastic (A150TEP)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.127000E+00;
  comp[i][ 1] = -0.101327;
  comp[i][ 6] = -0.775501;
  comp[i][ 7] = -0.035057;
  comp[i][ 8] = -0.052316;
  comp[i][ 9] = -0.017422;
  comp[i][20] = -0.018378;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Acetone");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.899000E-01;
  comp[i][ 1] = -0.104122;
  comp[i][ 6] = -0.620405;
  comp[i][ 8] = -0.275473;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Acetylene");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.097000E-03;
  comp[i][ 1] = -0.077418;
  comp[i][ 6] = -0.922582;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Air (Dry, Near Sea Level)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.205000E-03;
  comp[i][ 6] = -0.000124;
  comp[i][ 7] = -0.755268;
  comp[i][ 8] = -0.231781;
  comp[i][18] = -0.012827;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Alanine");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.420000E+00;
  comp[i][ 1] = -0.079190;
  comp[i][ 6] = -0.404439;
  comp[i][ 7] = -0.157213;
  comp[i][ 8] = -0.359159;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Aluminum");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.698900E+00;
  comp[i][13] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Aluminum Oxide");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -3.970000E+00;
  comp[i][ 8] = -0.470749;
  comp[i][13] = -0.529251;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Aluminum, Alloy 2024-O");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.780000E+00;
  comp[i][12] = -0.015000;
  comp[i][13] = -0.927000;
  comp[i][14] = -0.002830;
  comp[i][22] = -0.000850;
  comp[i][24] = -0.000570;
  comp[i][25] = -0.006000;
  comp[i][26] = -0.002830;
  comp[i][29] = -0.043500;
  comp[i][30] = -0.001420;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Aluminum, Alloy 2090-T83");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.590000E+00;
  comp[i][ 3] = -0.022500;
  comp[i][12] = -0.001630;
  comp[i][13] = -0.944000;
  comp[i][14] = -0.000650;
  comp[i][22] = -0.000980;
  comp[i][24] = -0.000330;
  comp[i][25] = -0.000330;
  comp[i][26] = -0.000780;
  comp[i][29] = -0.027000;
  comp[i][30] = -0.000650;
  comp[i][40] = -0.001150;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Aluminum, Alloy 3003");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.730000E+00;
  comp[i][13] = -0.978500;
  comp[i][14] = -0.003320;
  comp[i][25] = -0.012500;
  comp[i][26] = -0.003880;
  comp[i][29] = -0.001250;
  comp[i][30] = -0.000550;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Aluminum, Alloy 4043-O");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.690000E+00;
  comp[i][ 4] = -0.000005;
  comp[i][12] = -0.000280;
  comp[i][13] = -0.939000;
  comp[i][14] = -0.052500;
  comp[i][22] = -0.001130;
  comp[i][25] = -0.000280;
  comp[i][26] = -0.004530;
  comp[i][29] = -0.001700;
  comp[i][30] = -0.000570;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Aluminum, Alloy 5086-O");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.660000E+00;
  comp[i][12] = -0.040000;
  comp[i][13] = -0.946500;
  comp[i][40] = -0.002140;
  comp[i][22] = -0.000800;
  comp[i][24] = -0.001500;
  comp[i][25] = -0.004500;
  comp[i][26] = -0.002680;
  comp[i][29] = -0.000540;
  comp[i][30] = -0.001340;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Aluminum, Alloy 6061-O");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.700000E+00;
  comp[i][12] = -0.010000;
  comp[i][13] = -0.972000;
  comp[i][14] = -0.006000;
  comp[i][22] = -0.000880;
  comp[i][24] = -0.001950;
  comp[i][25] = -0.000880;
  comp[i][26] = -0.004090;
  comp[i][29] = -0.002750;
  comp[i][30] = -0.001460;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Aluminum, Alloy 7075-O");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.810000E+00;
  comp[i][12] = -0.025000;
  comp[i][13] = -0.892500;
  comp[i][14] = -0.002340;
  comp[i][22] = -0.001170;
  comp[i][24] = -0.002300;
  comp[i][25] = -0.001760;
  comp[i][26] = -0.002930;
  comp[i][29] = -0.016000;
  comp[i][30] = -0.056000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Ammonia (Liquid at T= -79°C)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.710000E-01;
  comp[i][ 1] = -0.177547;
  comp[i][ 7] = -0.822453;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Anthracene");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.250000E+00;
  comp[i][ 1] = -0.056553;
  comp[i][ 6] = -0.943447;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Argon");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.662000E-03;
  comp[i][18] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Asphalt");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.300000E+00;
  comp[i][ 1] = -0.103725;
  comp[i][ 6] = -0.848050;
  comp[i][ 7] = -0.006050;
  comp[i][ 8] = -0.004050;
  comp[i][16] = -0.037700;
  comp[i][23] = -0.000393;
  comp[i][28] = -0.000034;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Asphalt Pavement");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.578400E+00;
  comp[i][ 1] = -0.007781;
  comp[i][ 6] = -0.076175;
  comp[i][ 7] = -0.000363;
  comp[i][ 8] = -0.459103;
  comp[i][11] = -0.011659;
  comp[i][12] = -0.021757;
  comp[i][13] = -0.051009;
  comp[i][14] = -0.231474;
  comp[i][16] = -0.002804;
  comp[i][19] = -0.017058;
  comp[i][20] = -0.084471;
  comp[i][22] = -0.003403;
  comp[i][23] = -0.000024;
  comp[i][25] = -0.000362;
  comp[i][26] = -0.031375;
  comp[i][28] = -0.000002;
  comp[i][82] = -0.001179;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Bakelite");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.250000E+00;
  comp[i][ 1] = -0.057444;
  comp[i][ 6] = -0.774589;
  comp[i][ 8] = -0.167968;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Barium Fluoride");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -4.890000E+00;
  comp[i][ 9] = -0.216720;
  comp[i][56] = -0.783280;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Barium Sulfate");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -4.500000E+00;
  comp[i][ 8] = -0.274212;
  comp[i][16] = -0.137368;
  comp[i][56] = -0.588420;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Benzene");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -8.765000E-01;
  comp[i][ 1] = -0.077418;
  comp[i][ 6] = -0.922582;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Beryllium");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.848000E+00;
  comp[i][ 4] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Beryllium Carbide");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.900000E+00;
  comp[i][ 4] = -0.600111;
  comp[i][ 6] = -0.399889;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Beryllium Oxide");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -3.010000E+00;
  comp[i][ 4] = -0.360320;
  comp[i][ 8] = -0.639680;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Bismuth");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -9.747000E+00;
  comp[i][83] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Bismuth Germanate (BGO)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.130000E+00;
  comp[i][ 8] = -0.154126;
  comp[i][32] = -0.174820;
  comp[i][83] = -0.671054;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Blood (ICRP)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.060000E+00;
  comp[i][ 1] = -0.101866;
  comp[i][ 6] = -0.100020;
  comp[i][ 7] = -0.029640;
  comp[i][ 8] = -0.759414;
  comp[i][11] = -0.001850;
  comp[i][12] = -0.000040;
  comp[i][14] = -0.000030;
  comp[i][15] = -0.000350;
  comp[i][16] = -0.001850;
  comp[i][17] = -0.002780;
  comp[i][19] = -0.001630;
  comp[i][20] = -0.000060;
  comp[i][26] = -0.000460;
  comp[i][30] = -0.000010;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Bone Equivalent Plastic, B-100");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.450000E+00;
  comp[i][ 1] = -0.065471;
  comp[i][ 6] = -0.536945;
  comp[i][ 7] = -0.021500;
  comp[i][ 8] = -0.032085;
  comp[i][ 9] = -0.167411;
  comp[i][20] = -0.176589;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Bone Equivalent Plastic, B-110");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.785000E+00;
  comp[i][ 1] = -0.035500;
  comp[i][ 6] = -0.367300;
  comp[i][ 7] = -0.039700;
  comp[i][ 8] = -0.045300;
  comp[i][ 9] = -0.249300;
  comp[i][20] = -0.262900;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Bone, Compact (ICRU)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.850000E+00;
  comp[i][ 1] = -0.063984;
  comp[i][ 6] = -0.278000;
  comp[i][ 7] = -0.027000;
  comp[i][ 8] = -0.410016;
  comp[i][12] = -0.002000;
  comp[i][15] = -0.070000;
  comp[i][16] = -0.002000;
  comp[i][20] = -0.147000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Bone, Cortical (ICRP)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.850000E+00;
  comp[i][ 1] = -0.047234;
  comp[i][ 6] = -0.144330;
  comp[i][ 7] = -0.041990;
  comp[i][ 8] = -0.446096;
  comp[i][12] = -0.002200;
  comp[i][15] = -0.104970;
  comp[i][16] = -0.003150;
  comp[i][20] = -0.209930;
  comp[i][30] = -0.000100;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Boral (65%% Al-35%% B4C)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.530000E+00;
  comp[i][ 5] = -0.274000;
  comp[i][ 6] = -0.076000;
  comp[i][13] = -0.650000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Boral (Aluminum 10%% Boron Alloy)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.600000E+00;
  comp[i][ 5] = -0.100000;
  comp[i][11] = -0.005000;
  comp[i][13] = -0.879000;
  comp[i][14] = -0.002500;
  comp[i][19] = -0.010000;
  comp[i][22] = -0.000500;
  comp[i][26] = -0.003000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Boral (Aluminum 5%% Boron Alloy)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.600000E+00;
  comp[i][ 5] = -0.050000;
  comp[i][11] = -0.005000;
  comp[i][13] = -0.929500;
  comp[i][14] = -0.002000;
  comp[i][19] = -0.010000;
  comp[i][22] = -0.000500;
  comp[i][26] = -0.003000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Borax");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.730000E+00;
  comp[i][ 1] = -0.052859;
  comp[i][ 5] = -0.113391;
  comp[i][ 8] = -0.713187;
  comp[i][11] = -0.120563;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Boric Acid");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.500000E+00;
  comp[i][ 1] = -0.048903;
  comp[i][ 5] = -0.174842;
  comp[i][ 8] = -0.776255;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Boron");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.370000E+00;
  comp[i][ 5] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Boron Carbide");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.520000E+00;
  comp[i][ 5] = -0.782610;
  comp[i][ 6] = -0.217390;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Boron Fluoride (B2F4)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -4.058000E-03;
  comp[i][ 5] = -0.221501;
  comp[i][ 9] = -0.778499;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Boron Fluoride (BF3)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.831000E-03;
  comp[i][ 5] = -0.159440;
  comp[i][ 9] = -0.840560;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Boron Oxide");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.812000E+00;
  comp[i][ 5] = -0.310551;
  comp[i][ 8] = -0.689449;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Brain (ICRP)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.030000E+00;
  comp[i][ 1] = -0.110667;
  comp[i][ 6] = -0.125420;
  comp[i][ 7] = -0.013280;
  comp[i][ 8] = -0.737723;
  comp[i][11] = -0.001840;
  comp[i][12] = -0.000150;
  comp[i][15] = -0.003540;
  comp[i][16] = -0.001770;
  comp[i][17] = -0.002360;
  comp[i][19] = -0.003100;
  comp[i][20] = -0.000090;
  comp[i][26] = -0.000050;
  comp[i][30] = -0.000010;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Brass (Typical Composition)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -8.070000E+00;
  comp[i][26] = -0.000868;
  comp[i][29] = -0.665381;
  comp[i][30] = -0.325697;
  comp[i][50] = -0.002672;
  comp[i][82] = -0.005377;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Brick, Common Silica");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.800000E+00;
  comp[i][ 8] = -0.525000;
  comp[i][13] = -0.005000;
  comp[i][14] = -0.449000;
  comp[i][20] = -0.014000;
  comp[i][26] = -0.007000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Brick, Fire");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.100000E+00;
  comp[i][ 8] = -0.497000;
  comp[i][12] = -0.006000;
  comp[i][13] = -0.212000;
  comp[i][14] = -0.252000;
  comp[i][20] = -0.007000;
  comp[i][22] = -0.012000;
  comp[i][26] = -0.014000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Brick, Kaolin (White)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.100000E+00;
  comp[i][ 8] = -0.500318;
  comp[i][12] = -0.001205;
  comp[i][13] = -0.240568;
  comp[i][14] = -0.242823;
  comp[i][20] = -0.000714;
  comp[i][22] = -0.010179;
  comp[i][26] = -0.004192;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Bronze (Typical Composition)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -8.400000E+00;
  comp[i][13] = -0.028528;
  comp[i][14] = -0.003339;
  comp[i][25] = -0.003555;
  comp[i][26] = -0.010208;
  comp[i][28] = -0.006718;
  comp[i][29] = -0.874157;
  comp[i][30] = -0.036037;
  comp[i][50] = -0.024503;
  comp[i][82] = -0.012957;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "C-552 Air-Equivalent Plastic");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.760000E+00;
  comp[i][ 1] = -0.024680;
  comp[i][ 6] = -0.501610;
  comp[i][ 8] = -0.004527;
  comp[i][ 9] = -0.465209;
  comp[i][14] = -0.003973;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Cadmium");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -8.650000E+00;
  comp[i][48] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Cadmium Nitrate Tetrahydrate");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.450000E+00;
  comp[i][ 1] = -0.026139;
  comp[i][ 7] = -0.090811;
  comp[i][ 8] = -0.518650;
  comp[i][48] = -0.364401;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Cadmium Telluride");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -6.200000E+00;
  comp[i][48] = -0.468355;
  comp[i][52] = -0.531645;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Cadmium Tungstate (CWO)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.900000E+00;
  comp[i][ 8] = -0.177644;
  comp[i][48] = -0.312027;
  comp[i][74] = -0.510329;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Calcium Carbonate");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.800000E+00;
  comp[i][ 6] = -0.120003;
  comp[i][ 8] = -0.479554;
  comp[i][20] = -0.400443;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Calcium Fluoride");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -3.180000E+00;
  comp[i][ 9] = -0.486659;
  comp[i][20] = -0.513341;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Calcium Oxide");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -3.300000E+00;
  comp[i][ 8] = -0.285299;
  comp[i][20] = -0.714701;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Calcium Sulfate");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.960000E+00;
  comp[i][ 8] = -0.470095;
  comp[i][16] = -0.235497;
  comp[i][20] = -0.294408;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Carbon Dioxide");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.842000E-03;
  comp[i][ 6] = -0.272912;
  comp[i][ 8] = -0.727088;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Carbon Tetrachloride");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.594000E+00;
  comp[i][ 6] = -0.078083;
  comp[i][17] = -0.921917;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Carbon, Activated");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -3.200000E-01;
  comp[i][ 5] = -0.000001;
  comp[i][ 6] = -0.999999;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Carbon, Amorphous");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.000000E+00;
  comp[i][ 5] = -0.000001;
  comp[i][ 6] = -0.999999;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Carbon, Graphite (Reactor Grade)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.700000E+00;
  comp[i][ 5] = -0.000001;
  comp[i][ 6] = -0.999999;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Cat Litter (Clumping)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.100000E+00;
  comp[i][ 1] = -0.040400;
  comp[i][ 8] = -0.641100;
  comp[i][11] = -0.008400;
  comp[i][13] = -0.098300;
  comp[i][14] = -0.204600;
  comp[i][20] = -0.007300;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Cat Litter (Non-clumping)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.100000E+00;
  comp[i][ 1] = -0.013732;
  comp[i][ 8] = -0.539919;
  comp[i][11] = -0.043271;
  comp[i][12] = -0.050466;
  comp[i][13] = -0.052132;
  comp[i][14] = -0.293185;
  comp[i][19] = -0.003765;
  comp[i][20] = -0.001341;
  comp[i][26] = -0.002188;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Cellulose Acetate");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.420000E+00;
  comp[i][ 1] = -0.062162;
  comp[i][ 6] = -0.444462;
  comp[i][ 8] = -0.493376;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Celotex");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.400000E-01;
  comp[i][ 1] = -0.062165;
  comp[i][ 6] = -0.444455;
  comp[i][ 8] = -0.493380;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Ceric Sulfate Dosimeter Solution");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.030000E+00;
  comp[i][ 1] = -0.107596;
  comp[i][ 7] = -0.000800;
  comp[i][ 8] = -0.874976;
  comp[i][16] = -0.014627;
  comp[i][58] = -0.002001;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Cerium Fluoride");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -6.160000E+00;
  comp[i][ 9] = -0.289153;
  comp[i][58] = -0.710847;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Cesium Iodide");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -4.510000E+00;
  comp[i][53] = -0.488451;
  comp[i][55] = -0.511549;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Chromium");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.180000E+00;
  comp[i][24] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Clay");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.200000E+00;
  comp[i][ 8] = -0.484345;
  comp[i][11] = -0.007608;
  comp[i][12] = -0.010691;
  comp[i][13] = -0.122125;
  comp[i][14] = -0.294194;
  comp[i][15] = -0.000113;
  comp[i][19] = -0.020427;
  comp[i][20] = -0.018957;
  comp[i][22] = -0.004668;
  comp[i][25] = -0.000064;
  comp[i][26] = -0.036804;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Coal, Anthracite");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -8.400000E-01;
  comp[i][ 1] = -0.024000;
  comp[i][ 6] = -0.937000;
  comp[i][ 7] = -0.009000;
  comp[i][ 8] = -0.024000;
  comp[i][16] = -0.006000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Coal, Bituminous");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.500000E-01;
  comp[i][ 1] = -0.056000;
  comp[i][ 6] = -0.845000;
  comp[i][ 7] = -0.016000;
  comp[i][ 8] = -0.070000;
  comp[i][16] = -0.013000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Coal, Lignite");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.500000E-01;
  comp[i][ 1] = -0.042000;
  comp[i][ 6] = -0.727000;
  comp[i][ 7] = -0.012000;
  comp[i][ 8] = -0.213000;
  comp[i][16] = -0.006000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Concrete, Barite (Type BA)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -3.350000E+00;
  comp[i][ 1] = -0.003585;
  comp[i][ 8] = -0.311622;
  comp[i][12] = -0.001195;
  comp[i][13] = -0.004183;
  comp[i][14] = -0.010457;
  comp[i][16] = -0.107858;
  comp[i][20] = -0.050194;
  comp[i][26] = -0.047505;
  comp[i][56] = -0.463400;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Concrete, Barytes-limonite");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -3.360000E+00;
  comp[i][ 1] = -0.010240;
  comp[i][ 8] = -0.378476;
  comp[i][11] = -0.000904;
  comp[i][12] = -0.002309;
  comp[i][13] = -0.005020;
  comp[i][14] = -0.013553;
  comp[i][16] = -0.076097;
  comp[i][20] = -0.053910;
  comp[i][25] = -0.001405;
  comp[i][26] = -0.137135;
  comp[i][56] = -0.320952;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Concrete, Boron Frits-baryte");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -3.100000E+00;
  comp[i][ 1] = -0.005626;
  comp[i][ 5] = -0.010449;
  comp[i][ 8] = -0.339596;
  comp[i][ 9] = -0.002311;
  comp[i][11] = -0.012157;
  comp[i][12] = -0.002311;
  comp[i][13] = -0.006430;
  comp[i][14] = -0.033256;
  comp[i][16] = -0.091932;
  comp[i][19] = -0.001005;
  comp[i][20] = -0.062896;
  comp[i][25] = -0.000201;
  comp[i][26] = -0.022003;
  comp[i][30] = -0.006631;
  comp[i][56] = -0.403195;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Concrete, Colemanite-baryte");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -3.200000E+00;
  comp[i][ 1] = -0.008564;
  comp[i][ 5] = -0.009874;
  comp[i][ 8] = -0.351537;
  comp[i][11] = -0.001108;
  comp[i][12] = -0.002217;
  comp[i][13] = -0.006146;
  comp[i][14] = -0.017733;
  comp[i][16] = -0.097028;
  comp[i][20] = -0.085239;
  comp[i][25] = -0.000101;
  comp[i][26] = -0.010378;
  comp[i][56] = -0.410076;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Concrete, Ferro-phosphorus");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -4.800000E+00;
  comp[i][ 1] = -0.005000;
  comp[i][ 8] = -0.104000;
  comp[i][12] = -0.002000;
  comp[i][13] = -0.004000;
  comp[i][14] = -0.034000;
  comp[i][15] = -0.197000;
  comp[i][20] = -0.042000;
  comp[i][26] = -0.612000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Concrete, Hanford Dry");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.180000E+00;
  comp[i][ 1] = -0.004000 ;
  comp[i][ 8] = -0.482102 ;
  comp[i][11] = -0.002168 ;
  comp[i][12] = -0.014094 ;
  comp[i][13] = -0.069387 ;
  comp[i][14] = -0.277549 ;
  comp[i][19] = -0.013010 ;
  comp[i][20] = -0.080229 ;
  comp[i][26] = -0.057461 ;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Concrete, Hanford Wet");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.350000E+00;
  comp[i][ 1] = -0.012309;
  comp[i][ 8] = -0.513359;
  comp[i][11] = -0.002001;
  comp[i][12] = -0.013009;
  comp[i][13] = -0.064045;
  comp[i][14] = -0.256179;
  comp[i][19] = -0.012008;
  comp[i][20] = -0.074052;
  comp[i][26] = -0.053037;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Concrete, Iron-limonite");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -4.400000E+00;
  comp[i][ 1] = -0.000500;
  comp[i][ 8] = -0.179910;
  comp[i][12] = -0.001999;
  comp[i][13] = -0.004998;
  comp[i][14] = -0.013993;
  comp[i][16] = -0.001000;
  comp[i][20] = -0.060970;
  comp[i][25] = -0.015992;
  comp[i][26] = -0.720640;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Concrete, Iron-Portland");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -5.900000E+00;
  comp[i][ 1] = -0.003321;
  comp[i][ 8] = -0.058563;
  comp[i][12] = -0.001308;
  comp[i][13] = -0.003321;
  comp[i][14] = -0.009157;
  comp[i][16] = -0.000503;
  comp[i][20] = -0.039847;
  comp[i][25] = -0.003522;
  comp[i][26] = -0.880459;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Concrete, Limonite and Steel");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -4.540000E+00;
  comp[i][ 1] = -0.006840 ;
  comp[i][ 8] = -0.156222 ;
  comp[i][12] = -0.001545 ;
  comp[i][13] = -0.006399 ;
  comp[i][14] = -0.014784 ;
  comp[i][19] = -0.000883 ;
  comp[i][20] = -0.057590 ;
  comp[i][23] = -0.000883 ;
  comp[i][26] = -0.754854 ;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Concrete, Los Alamos (MCNP)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.250000E+00;
  comp[i][ 1] = -0.004530;
  comp[i][ 8] = -0.512600;
  comp[i][11] = -0.015270;
  comp[i][13] = -0.035550;
  comp[i][14] = -0.360360;
  comp[i][20] = -0.057910;
  comp[i][26] = -0.013780;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Concrete, Luminite-colemanite-baryte");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -3.100000E+00;
  comp[i][ 1] = -0.010957;
  comp[i][ 5] = -0.008846;
  comp[i][ 8] = -0.371431;
  comp[i][11] = -0.001106;
  comp[i][12] = -0.001407;
  comp[i][13] = -0.017692;
  comp[i][14] = -0.009650;
  comp[i][16] = -0.091074;
  comp[i][20] = -0.055086;
  comp[i][22] = -0.012766;
  comp[i][25] = -0.001206;
  comp[i][26] = -0.030860;
  comp[i][56] = -0.387917;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Concrete, Luminite-Portland-colemanite-baryte");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -3.100000E+00;
  comp[i][ 1] = -0.011126 ;
  comp[i][ 5] = -0.010316 ;
  comp[i][ 8] = -0.374023 ;
  comp[i][11] = -0.001113 ;
  comp[i][12] = -0.002023 ;
  comp[i][13] = -0.013351 ;
  comp[i][14] = -0.015070 ;
  comp[i][16] = -0.090724 ;
  comp[i][20] = -0.077576 ;
  comp[i][22] = -0.000718 ;
  comp[i][25] = -0.000405 ;
  comp[i][26] = -0.018914 ;
  comp[i][56] = -0.384643 ;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Concrete, M-1");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -4.500000E+00;
  comp[i][ 1] = -0.008000;
  comp[i][ 5] = -0.009000;
  comp[i][ 8] = -0.107000;
  comp[i][12] = -0.043000;
  comp[i][17] = -0.021000;
  comp[i][25] = -0.003000;
  comp[i][20] = -0.011000;
  comp[i][26] = -0.798000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Concrete, Magnetite");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -3.530000E+00;
  comp[i][ 1] = -0.003113;
  comp[i][ 8] = -0.330504;
  comp[i][12] = -0.009338;
  comp[i][13] = -0.023486;
  comp[i][14] = -0.025750;
  comp[i][16] = -0.001415;
  comp[i][20] = -0.071024;
  comp[i][22] = -0.054329;
  comp[i][23] = -0.003113;
  comp[i][24] = -0.001698;
  comp[i][25] = -0.001981;
  comp[i][26] = -0.474250;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Concrete, Magnetite and Steel");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -4.640000E+00;
  comp[i][ 1] = -0.002374 ;
  comp[i][ 8] = -0.137678 ;
  comp[i][12] = -0.003669 ;
  comp[i][13] = -0.010358 ;
  comp[i][14] = -0.015753 ;
  comp[i][20] = -0.055675 ;
  comp[i][22] = -0.015969 ;
  comp[i][23] = -0.000647 ;
  comp[i][26] = -0.757877 ;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Concrete, Magnuson");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.147000E+00;
  comp[i][ 1] = -0.003319;
  comp[i][ 6] = -0.105320;
  comp[i][ 8] = -0.499428;
  comp[i][11] = -0.001411;
  comp[i][12] = -0.094200;
  comp[i][13] = -0.007859;
  comp[i][14] = -0.042101;
  comp[i][16] = -0.002483;
  comp[i][17] = -0.000523;
  comp[i][19] = -0.009445;
  comp[i][20] = -0.226317;
  comp[i][22] = -0.001488;
  comp[i][25] = -0.000512;
  comp[i][26] = -0.005595;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Concrete, MO");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -5.500000E+00;
  comp[i][ 1] = -0.005000;
  comp[i][ 8] = -0.060000;
  comp[i][12] = -0.037000;
  comp[i][25] = -0.004000;
  comp[i][17] = -0.013000;
  comp[i][26] = -0.881000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Concrete, Oak Ridge (ORNL)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.300000E+00;
  comp[i][ 1] = -0.006187;
  comp[i][ 6] = -0.175193;
  comp[i][ 8] = -0.410184;
  comp[i][11] = -0.000271;
  comp[i][12] = -0.032649;
  comp[i][13] = -0.010830;
  comp[i][14] = -0.034479;
  comp[i][19] = -0.001138;
  comp[i][20] = -0.321287;
  comp[i][26] = -0.007784;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Concrete, Ordinary (NBS 03)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.350000E+00;
  comp[i][ 1] = -0.008485;
  comp[i][ 6] = -0.050064;
  comp[i][ 8] = -0.473483;
  comp[i][12] = -0.024183;
  comp[i][13] = -0.036063;
  comp[i][14] = -0.145100;
  comp[i][16] = -0.002970;
  comp[i][19] = -0.001697;
  comp[i][20] = -0.246924;
  comp[i][26] = -0.011031;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Concrete, Ordinary (NBS 04)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.350000E+00;
  comp[i][ 1] = -0.005558;
  comp[i][ 8] = -0.498076;
  comp[i][11] = -0.017101;
  comp[i][12] = -0.002565;
  comp[i][13] = -0.045746;
  comp[i][14] = -0.315092;
  comp[i][16] = -0.001283;
  comp[i][19] = -0.019239;
  comp[i][20] = -0.082941;
  comp[i][26] = -0.012398;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Concrete, Ordinary (NIST)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.300000E+00;
  comp[i][ 1] = -0.022100;
  comp[i][ 6] = -0.002484;
  comp[i][ 8] = -0.574930;
  comp[i][11] = -0.015208;
  comp[i][12] = -0.001266;
  comp[i][13] = -0.019953;
  comp[i][14] = -0.304627;
  comp[i][19] = -0.010045;
  comp[i][20] = -0.042951;
  comp[i][26] = -0.006435;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Concrete, Portland");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.300000E+00;
  comp[i][ 1] = -0.010000 ;
  comp[i][ 6] = -0.001000 ;
  comp[i][ 8] = -0.529107 ;
  comp[i][11] = -0.016000 ;
  comp[i][12] = -0.002000 ;
  comp[i][13] = -0.033872 ;
  comp[i][14] = -0.337021 ;
  comp[i][19] = -0.013000 ;
  comp[i][20] = -0.044000 ;
  comp[i][26] = -0.014000 ;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Concrete, Regular");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.300000E+00;
  comp[i][ 1] = -0.010000;
  comp[i][ 8] = -0.532000;
  comp[i][11] = -0.029000;
  comp[i][13] = -0.034000;
  comp[i][14] = -0.337000;
  comp[i][20] = -0.044000;
  comp[i][26] = -0.014000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Concrete, Rocky Flats");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.320000E+00;
  comp[i][ 1] = -0.007500;
  comp[i][ 6] = -0.055200;
  comp[i][ 7] = -0.000200;
  comp[i][ 8] = -0.484900;
  comp[i][11] = -0.006300;
  comp[i][12] = -0.012500;
  comp[i][13] = -0.021700;
  comp[i][14] = -0.155000;
  comp[i][16] = -0.001900;
  comp[i][19] = -0.013700;
  comp[i][20] = -0.230000;
  comp[i][22] = -0.001000;
  comp[i][26] = -0.010100;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Concrete, Serpentine");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.100000E+00;
  comp[i][ 1] = -0.015909 ;
  comp[i][ 6] = -0.000909 ;
  comp[i][ 8] = -0.511818 ;
  comp[i][11] = -0.004091 ;
  comp[i][12] = -0.135000 ;
  comp[i][13] = -0.019091 ;
  comp[i][14] = -0.209091 ;
  comp[i][19] = -0.004091 ;
  comp[i][20] = -0.068182 ;
  comp[i][24] = -0.000909 ;
  comp[i][26] = -0.030909 ;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Copper");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -8.960000E+00;
  comp[i][29] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Diatomaceous Earth");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.200000E-01;
  comp[i][ 1] = -0.008956;
  comp[i][ 8] = -0.546579;
  comp[i][11] = -0.009896;
  comp[i][12] = -0.002774;
  comp[i][13] = -0.015581;
  comp[i][14] = -0.394761;
  comp[i][19] = -0.011074;
  comp[i][20] = -0.003945;
  comp[i][26] = -0.006434;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Earth, Typical Western U.S.");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.520000E+00;
  comp[i][ 1] = -0.023834;
  comp[i][ 8] = -0.598898;
  comp[i][13] = -0.080446;
  comp[i][14] = -0.296821;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Earth, U.S. Average");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.520000E+00;
  comp[i][ 8] = -0.513713 ;
  comp[i][11] = -0.006140 ;
  comp[i][12] = -0.013303 ;
  comp[i][13] = -0.068563 ;
  comp[i][14] = -0.271183 ;
  comp[i][19] = -0.014327 ;
  comp[i][20] = -0.051167 ;
  comp[i][22] = -0.004605 ;
  comp[i][25] = -0.000716 ;
  comp[i][26] = -0.056283 ;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Ethane");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.253000E-03;
  comp[i][ 1] = -0.201125;
  comp[i][ 6] = -0.798875;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Ethyl Acetate");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -9.010000E-01;
  comp[i][ 1] = -0.091522;
  comp[i][ 6] = -0.545290;
  comp[i][ 8] = -0.363189;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Ethyl Alcohol");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.893000E-01;
  comp[i][ 1] = -0.131269;
  comp[i][ 6] = -0.521438;
  comp[i][ 8] = -0.347294;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Ethylene");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.175000E-03;
  comp[i][ 1] = -0.143711;
  comp[i][ 6] = -0.856289;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Ethylene Glycol");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.114000E+00;
  comp[i][ 1] = -0.097436;
  comp[i][ 6] = -0.387018;
  comp[i][ 8] = -0.515546;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Explosive Compound, AN");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.720000E+00;
  comp[i][ 1] = -0.050370;
  comp[i][ 7] = -0.349978;
  comp[i][ 8] = -0.599652;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Explosive Compound, EGDN");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.490000E+00;
  comp[i][ 1] = -0.026514;
  comp[i][ 6] = -0.157970;
  comp[i][ 7] = -0.184222;
  comp[i][ 8] = -0.631294;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Explosive Compound, HMX");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.890000E+00;
  comp[i][ 1] = -0.027227;
  comp[i][ 6] = -0.162222;
  comp[i][ 7] = -0.378361;
  comp[i][ 8] = -0.432190;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Explosive Compound, NC");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.490000E+00;
  comp[i][ 1] = -0.029216;
  comp[i][ 6] = -0.271296;
  comp[i][ 7] = -0.121276;
  comp[i][ 8] = -0.578212;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Explosive Compound, NG");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.600000E+00;
  comp[i][ 1] = -0.022193;
  comp[i][ 6] = -0.158671;
  comp[i][ 7] = -0.185040;
  comp[i][ 8] = -0.634096;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Explosive Compound, PETN");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.770000E+00;
  comp[i][ 1] = -0.025506;
  comp[i][ 6] = -0.189961;
  comp[i][ 7] = -0.177223;
  comp[i][ 8] = -0.607310;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Explosive Compound, RDX");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.820000E+00;
  comp[i][ 1] = -0.027227;
  comp[i][ 6] = -0.162222;
  comp[i][ 7] = -0.378361;
  comp[i][ 8] = -0.432190;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Explosive Compound, TNT");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.650000E+00;
  comp[i][ 1] = -0.022189;
  comp[i][ 6] = -0.370160;
  comp[i][ 7] = -0.185004;
  comp[i][ 8] = -0.422648;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Eye Lens (ICRP)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.100000E+00;
  comp[i][ 1] = -0.099269;
  comp[i][ 6] = -0.193710;
  comp[i][ 7] = -0.053270;
  comp[i][ 8] = -0.653751;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Felt");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.850000E-01;
  comp[i][ 1] = -0.044200 ;
  comp[i][ 6] = -0.434600 ;
  comp[i][ 7] = -0.176500 ;
  comp[i][ 8] = -0.344700 ;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Ferric Oxide");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -5.200000E+00;
  comp[i][ 8] = -0.300567;
  comp[i][26] = -0.699433;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Ferrous Sulfate Dosimeter Solution");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.024000E+00;
  comp[i][ 1] = -0.108259;
  comp[i][ 7] = -0.000027;
  comp[i][ 8] = -0.878636;
  comp[i][11] = -0.000022;
  comp[i][16] = -0.012968;
  comp[i][17] = -0.000034;
  comp[i][26] = -0.000054;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Fertilizer (Muriate of Potash)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.270000E+00;
  comp[i][ 1] = -0.000050;
  comp[i][ 8] = -0.000718;
  comp[i][11] = -0.008487;
  comp[i][12] = -0.000206;
  comp[i][16] = -0.000159;
  comp[i][17] = -0.477922;
  comp[i][19] = -0.511852;
  comp[i][20] = -0.000276;
  comp[i][35] = -0.000330;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Fiberglass, Type C");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.490000E+00;
  comp[i][ 5] = -0.018579 ;
  comp[i][ 8] = -0.478631 ;
  comp[i][11] = -0.059171 ;
  comp[i][12] = -0.018037 ;
  comp[i][13] = -0.021107 ;
  comp[i][14] = -0.302924 ;
  comp[i][16] = -0.000399 ;
  comp[i][20] = -0.099757 ;
  comp[i][26] = -0.001395 ;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Fiberglass, Type E");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.565000E+00;
  comp[i][ 5] = -0.022803;
  comp[i][ 8] = -0.471950;
  comp[i][ 9] = -0.004895;
  comp[i][11] = -0.007262;
  comp[i][12] = -0.014759;
  comp[i][13] = -0.072536;
  comp[i][14] = -0.247102;
  comp[i][19] = -0.008127;
  comp[i][20] = -0.143428;
  comp[i][22] = -0.004400;
  comp[i][26] = -0.002739;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Fiberglass, Type R");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.550000E+00;
  comp[i][ 8] = -0.486722;
  comp[i][12] = -0.036182;
  comp[i][13] = -0.132313;
  comp[i][14] = -0.280461;
  comp[i][20] = -0.064322;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Freon-12");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.120000E+00;
  comp[i][ 6] = -0.099335;
  comp[i][ 9] = -0.314247;
  comp[i][17] = -0.586418;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Freon-12B2");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.800000E+00;
  comp[i][ 6] = -0.057245;
  comp[i][ 9] = -0.181096;
  comp[i][35] = -0.761659;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Freon-13");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -9.500000E-01;
  comp[i][ 6] = -0.114983;
  comp[i][ 9] = -0.545622;
  comp[i][17] = -0.339396;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Freon-13B1");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.500000E+00;
  comp[i][ 6] = -0.080659;
  comp[i][ 9] = -0.382749;
  comp[i][35] = -0.536592;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Freon-13I1");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.800000E+00;
  comp[i][ 6] = -0.061309;
  comp[i][ 9] = -0.290924;
  comp[i][53] = -0.647767;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Gadolinium");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.900400E+00;
  comp[i][64] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Gadolinium Oxysulfide");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.440000E+00;
  comp[i][ 8] = -0.084528;
  comp[i][16] = -0.084690;
  comp[i][64] = -0.830782;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Gadolinium Silicate (GSO)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -6.710000E+00;
  comp[i][ 8] = -0.189305;
  comp[i][14] = -0.066462;
  comp[i][64] = -0.744233;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Gafchromic Sensor (GS)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.300000E+00;
  comp[i][ 1] = -0.089700;
  comp[i][ 6] = -0.605800;
  comp[i][ 7] = -0.112200;
  comp[i][ 8] = -0.192300;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Gallium Arsenide");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -5.310000E+00;
  comp[i][31] = -0.482030;
  comp[i][33] = -0.517970;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Gasoline");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.210000E-01;
  comp[i][ 1] = -0.157000;
  comp[i][ 6] = -0.843000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Germanium, High Purity");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -5.323000E+00;
  comp[i][32] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Glass, Foam");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.280000E-01;
  comp[i][ 1] = -0.001000;
  comp[i][ 5] = -0.015000;
  comp[i][ 8] = -0.534000;
  comp[i][11] = -0.161000;
  comp[i][14] = -0.279000;
  comp[i][16] = -0.010000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Glass, Lead");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -6.220000E+00;
  comp[i][ 8] = -0.156453;
  comp[i][14] = -0.080866;
  comp[i][22] = -0.008092;
  comp[i][33] = -0.002651;
  comp[i][82] = -0.751938;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Glass, Plate");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.400000E+00;
  comp[i][ 8] = -0.459800;
  comp[i][11] = -0.096441;
  comp[i][14] = -0.336553;
  comp[i][20] = -0.107205;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Glycerol");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.261300E+00;
  comp[i][ 1] = -0.087554;
  comp[i][ 6] = -0.391262;
  comp[i][ 8] = -0.521185;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Gold");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.932000E+01;
  comp[i][79] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Gypsum (Plaster of Paris)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.320000E+00;
  comp[i][ 1] = -0.023416;
  comp[i][ 8] = -0.557572;
  comp[i][16] = -0.186215;
  comp[i][20] = -0.232797;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "He-3 Proportional Gas");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.250000E-04;
  comp[i][ 2] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Helium, Natural");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.660000E-04;
  comp[i][ 2] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Hydrogen");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -8.400000E-05;
  comp[i][ 1] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Incoloy-800");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.940000E+00;
  comp[i][ 6] = -0.000650;
  comp[i][13] = -0.003750;
  comp[i][14] = -0.006500;
  comp[i][16] = -0.000100;
  comp[i][22] = -0.003750;
  comp[i][24] = -0.210000;
  comp[i][25] = -0.009750;
  comp[i][26] = -0.435630;
  comp[i][28] = -0.325000;
  comp[i][29] = -0.004880;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Inconel-600");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -8.470000E+00;
  comp[i][ 6] = -0.000980;
  comp[i][14] = -0.003250;
  comp[i][16] = -0.000100;
  comp[i][24] = -0.155000;
  comp[i][25] = -0.006500;
  comp[i][26] = -0.080000;
  comp[i][28] = -0.750930;
  comp[i][29] = -0.003250;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Inconel-625");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.320000E+00;
  comp[i][ 1] = -0.007500;
  comp[i][ 6] = -0.055200;
  comp[i][ 7] = -0.000200;
  comp[i][ 8] = -0.484900;
  comp[i][11] = -0.006300;
  comp[i][12] = -0.012500;
  comp[i][13] = -0.021700;
  comp[i][14] = -0.155000;
  comp[i][16] = -0.001900;
  comp[i][19] = -0.013700;
  comp[i][20] = -0.230000;
  comp[i][22] = -0.001000;
  comp[i][26] = -0.010100;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Inconel-718");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -8.190000E+00;
  comp[i][ 5] = -0.000050;
  comp[i][ 6] = -0.000730;
  comp[i][13] = -0.005000;
  comp[i][14] = -0.003180;
  comp[i][15] = -0.000140;
  comp[i][16] = -0.000140;
  comp[i][22] = -0.009000;
  comp[i][24] = -0.190000;
  comp[i][25] = -0.003180;
  comp[i][26] = -0.170000;
  comp[i][28] = -0.525000;
  comp[i][27] = -0.009100;
  comp[i][29] = -0.002730;
  comp[i][41] = -0.051250;
  comp[i][42] = -0.030500;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Indium");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.310000E+00;
  comp[i][49] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Iron");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.874000E+00;
  comp[i][26] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Iron Boride (Fe2B)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.300000E+00;
  comp[i][ 5] = -0.088252;
  comp[i][26] = -0.911748;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Iron Boride (FeB)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.150000E+00;
  comp[i][ 5] = -0.162174;
  comp[i][26] = -0.837826;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Iron, Armco Ingot");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.866000E+00;
  comp[i][ 6] = -0.000120;
  comp[i][ 8] = -0.001100;
  comp[i][15] = -0.000050;
  comp[i][16] = -0.000250;
  comp[i][25] = -0.000170;
  comp[i][26] = -0.998310;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Iron, Cast (Gray)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.150000E+00;
  comp[i][ 6] = -0.034000;
  comp[i][14] = -0.026000;
  comp[i][15] = -0.003000;
  comp[i][16] = -0.001000;
  comp[i][25] = -0.006500;
  comp[i][26] = -0.929500;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Iron, Wrought (Byers No. 1)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.700000E+00;
  comp[i][ 6] = -0.000810 ;
  comp[i][14] = -0.001599 ;
  comp[i][15] = -0.000628 ;
  comp[i][16] = -0.000101 ;
  comp[i][25] = -0.000152 ;
  comp[i][26] = -0.996711 ;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Kaowool");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -9.600000E-02;
  comp[i][ 5] = -0.000248;
  comp[i][ 8] = -0.500064;
  comp[i][13] = -0.238163;
  comp[i][14] = -0.243627;
  comp[i][20] = -0.000715;
  comp[i][22] = -0.010189;
  comp[i][26] = -0.006994;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Kapton Polyimide Film");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.420000E+00;
  comp[i][ 1] = -0.026362;
  comp[i][ 6] = -0.691133;
  comp[i][ 7] = -0.073270;
  comp[i][ 8] = -0.209235;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Kennertium");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.680000E+01;
  comp[i][28] = -0.090000;
  comp[i][29] = -0.150000;
  comp[i][74] = -0.760000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Kernite");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.950000E+00;
  comp[i][ 1] = -0.029506;
  comp[i][ 5] = -0.158240;
  comp[i][ 8] = -0.644003;
  comp[i][11] = -0.168250;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Kerosene");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -8.190000E-01;
  comp[i][ 1] = -0.160000;
  comp[i][ 6] = -0.840000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Krypton");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -3.478000E-03;
  comp[i][36] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Kynar");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.790000E+00;
  comp[i][ 1] = -0.031481;
  comp[i][ 6] = -0.375135;
  comp[i][ 9] = -0.593384;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Lead");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.135000E+01;
  comp[i][82] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Lead Tungstate (PWO)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -8.240000E+00;
  comp[i][ 8] = -0.140642;
  comp[i][74] = -0.404011;
  comp[i][82] = -0.455347;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Lithium");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -5.340000E-01;
  comp[i][ 3] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Lithium Amide");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.178000E+00;
  comp[i][ 1] = -0.087783;
  comp[i][ 3] = -0.302262;
  comp[i][ 7] = -0.609955;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Lithium Fluoride");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.635000E+00;
  comp[i][ 3] = -0.267585;
  comp[i][ 9] = -0.732415;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Lithium Gadrium Borate (LGB)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -3.500000E+00;
  comp[i][ 3] = -0.098240;
  comp[i][ 5] = -0.081766;
  comp[i][ 8] = -0.391956;
  comp[i][64] = -0.428038;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Lithium Hydride");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -8.200000E-01;
  comp[i][ 1] = -0.126797;
  comp[i][ 3] = -0.873203;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Lithium Iodide (High Density)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -4.080000E+00;
  comp[i][ 3] = -0.051858;
  comp[i][53] = -0.948142;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Lithium Iodide (Low Density)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -3.494000E+00;
  comp[i][ 3] = -0.051858;
  comp[i][53] = -0.948142;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Lithium Oxide");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.013000E+00;
  comp[i][ 3] = -0.464570;
  comp[i][ 8] = -0.535430;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Lithium Tetraborate");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.440000E+00;
  comp[i][ 3] = -0.082085;
  comp[i][ 5] = -0.255680;
  comp[i][ 8] = -0.662235;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Lucite");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.190000E+00;
  comp[i][ 1] = -0.080538 ;
  comp[i][ 6] = -0.599848 ;
  comp[i][ 8] = -0.319614 ;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Lutetium Aluminum Garnet (LuAG)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -6.730000E+00;
  comp[i][ 8] = -0.225396;
  comp[i][13] = -0.158379;
  comp[i][71] = -0.616225;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Lutetium Orthoaluminate (LuAP)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.012000E-01;
  comp[i][ 8] = -0.192034;
  comp[i][13] = -0.107949;
  comp[i][71] = -0.700017;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Lutetium Oxyorthosilicate (LSO)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.400000E+00;
  comp[i][ 8] = -0.174660;
  comp[i][14] = -0.061320;
  comp[i][71] = -0.764021;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Lutetium Yttrium OxyorthoSilicate (LYSO)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.300000E+00;
  comp[i][ 8] = -0.125815;
  comp[i][14] = -0.044172;
  comp[i][39] = -0.279654;
  comp[i][71] = -0.550359;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Magnesium");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.740000E+00;
  comp[i][12] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Magnesium Oxide");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -3.580000E+00;
  comp[i][ 8] = -0.396964;
  comp[i][12] = -0.603036;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Magnesium Tetraborate");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.530000E+00;
  comp[i][ 5] = -0.240837 ;
  comp[i][ 8] = -0.623790 ;
  comp[i][12] = -0.135373 ;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Masonite");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.300000E+00;
  comp[i][ 1] = -0.062165;
  comp[i][ 6] = -0.444455;
  comp[i][ 8] = -0.493380;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Melamine");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.350000E+00;
  comp[i][ 1] = -0.046680;
  comp[i][ 6] = -0.397313;
  comp[i][ 7] = -0.556008;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Mercury");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.354600E+01;
  comp[i][80] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Mercury Iodide");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -6.360000E+00;
  comp[i][53] = -0.558560 ;
  comp[i][80] = -0.441440 ;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Methane");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -6.670000E-04;
  comp[i][ 1] = -0.251318;
  comp[i][ 6] = -0.748682;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Methanol");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.914000E-01;
  comp[i][ 1] = -0.125822;
  comp[i][ 6] = -0.374852;
  comp[i][ 8] = -0.499326;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Methylene Chloride");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.326600E+00;
  comp[i][ 1] = -0.023735;
  comp[i][ 6] = -0.141415;
  comp[i][17] = -0.834850;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Molybdenum");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.022000E+01;
  comp[i][42] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Monosodium Titanate, MST");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.000000E+00;
  comp[i][ 1] = -0.005047;
  comp[i][ 8] = -0.400528;
  comp[i][11] = -0.115105;
  comp[i][22] = -0.479320;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Muscle Equivalent-Liquid, with Sucrose");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.110000E+00;
  comp[i][ 1] = -0.098234;
  comp[i][ 6] = -0.156214;
  comp[i][ 7] = -0.035451;
  comp[i][ 8] = -0.710101;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Muscle Equivalent-Liquid, without Sucrose");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.070000E+00;
  comp[i][ 1] = -0.101969;
  comp[i][ 6] = -0.120058;
  comp[i][ 7] = -0.035451;
  comp[i][ 8] = -0.742522;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Muscle, Skeletal (ICRP)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.040000E+00;
  comp[i][ 1] = -0.100637;
  comp[i][ 6] = -0.107830;
  comp[i][ 7] = -0.027680;
  comp[i][ 8] = -0.754773;
  comp[i][11] = -0.000750;
  comp[i][12] = -0.000190;
  comp[i][15] = -0.001800;
  comp[i][16] = -0.002410;
  comp[i][17] = -0.000790;
  comp[i][19] = -0.003020;
  comp[i][20] = -0.000030;
  comp[i][26] = -0.000040;
  comp[i][30] = -0.000050;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Muscle, Striated (ICRU)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.040000E+00;
  comp[i][ 1] = -0.101997;
  comp[i][ 6] = -0.123000;
  comp[i][ 7] = -0.035000;
  comp[i][ 8] = -0.729003;
  comp[i][11] = -0.000800;
  comp[i][12] = -0.000200;
  comp[i][15] = -0.002000;
  comp[i][16] = -0.005000;
  comp[i][19] = -0.003000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Neon");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -8.390000E-04;
  comp[i][10] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Nickel");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -8.902000E+00;
  comp[i][28] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Niobium");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -8.570000E+00;
  comp[i][41] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Nitrogen");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.165000E-03;
  comp[i][ 7] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Nylon, Dupont ELVAmide 8062");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.080000E+00;
  comp[i][ 1] = -0.103509;
  comp[i][ 6] = -0.648416;
  comp[i][ 7] = -0.099536;
  comp[i][ 8] = -0.148539;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Nylon, Type 11 (Rilsan)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.425000E+00;
  comp[i][ 1] = -0.115476;
  comp[i][ 6] = -0.720819;
  comp[i][ 7] = -0.076417;
  comp[i][ 8] = -0.087289;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Nylon, Type 6 and Type 6/6");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.140000E+00;
  comp[i][ 1] = -0.097976;
  comp[i][ 6] = -0.636856;
  comp[i][ 7] = -0.123779;
  comp[i][ 8] = -0.141389;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Nylon, Type 6/10");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.140000E+00;
  comp[i][ 1] = -0.107062;
  comp[i][ 6] = -0.680449;
  comp[i][ 7] = -0.099189;
  comp[i][ 8] = -0.113300;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Oil, Crude (Heavy, Cold Lake, Canada)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -9.700000E-01;
  comp[i][ 1] = -0.104000 ;
  comp[i][ 6] = -0.837000 ;
  comp[i][ 7] = -0.004000 ;
  comp[i][ 8] = -0.011000 ;
  comp[i][16] = -0.044000 ;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Oil, Crude (Heavy, Mexican)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -9.750000E-01;
  comp[i][ 1] = -0.104039;
  comp[i][ 6] = -0.853733;
  comp[i][16] = -0.042228;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Oil, Crude (Heavy, Qayarah, Iraq)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -9.700000E-01;
  comp[i][ 1] = -0.102000;
  comp[i][ 6] = -0.807000;
  comp[i][ 7] = -0.007000;
  comp[i][16] = -0.084000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Oil, Crude (Light, Texas)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -8.750000E-01;
  comp[i][ 1] = -0.123246;
  comp[i][ 6] = -0.852204;
  comp[i][ 7] = -0.007014;
  comp[i][16] = -0.017535;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Oil, Fuel (California)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -9.550000E-01;
  comp[i][ 1] = -0.125878;
  comp[i][ 6] = -0.862308;
  comp[i][16] = -0.011814;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Oil, Hydraulic");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -8.710000E-01;
  comp[i][ 1] = -0.040495;
  comp[i][ 6] = -0.584904;
  comp[i][ 8] = -0.077915;
  comp[i][15] = -0.037709;
  comp[i][17] = -0.258977;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Oil, Lard");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -9.150000E-01;
  comp[i][ 1] = -0.117621;
  comp[i][ 6] = -0.778655;
  comp[i][ 8] = -0.103724;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Oxygen");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.332000E-03;
  comp[i][ 8] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "P-10 Gas");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.561000E-03;
  comp[i][ 1] = -0.010735;
  comp[i][ 6] = -0.031980;
  comp[i][18] = -0.957286;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "P-5 Gas");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.611000E-03;
  comp[i][ 1] = -0.005202;
  comp[i][ 6] = -0.015497;
  comp[i][18] = -0.979302;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Palladium");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.202000E+01;
  comp[i][46] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Photographic Emulsion, Gel in");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.291400E+00;
  comp[i][ 1] = -0.081180;
  comp[i][ 6] = -0.416060;
  comp[i][ 7] = -0.111240;
  comp[i][ 8] = -0.380640;
  comp[i][16] = -0.010880;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Photographic Emulsion, Kodak Type AA");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.200000E+00;
  comp[i][ 1] = -0.030500;
  comp[i][ 6] = -0.210700;
  comp[i][ 7] = -0.072100;
  comp[i][ 8] = -0.163200;
  comp[i][35] = -0.222800;
  comp[i][47] = -0.300700;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Photographic Emulsion, Standard Nuclear");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -3.815000E+00;
  comp[i][ 1] = -0.014100;
  comp[i][ 6] = -0.072261;
  comp[i][ 7] = -0.019320;
  comp[i][ 8] = -0.066101;
  comp[i][16] = -0.001890;
  comp[i][35] = -0.349104;
  comp[i][47] = -0.474105;
  comp[i][53] = -0.003120;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Platinum");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.145000E+01;
  comp[i][78] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Polycarbonate");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.200000E+00;
  comp[i][ 1] = -0.055491;
  comp[i][ 6] = -0.755751;
  comp[i][ 8] = -0.188758;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Polyethylene Terephthalate (PET)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.380000E+00;
  comp[i][ 1] = -0.041960;
  comp[i][ 6] = -0.625016;
  comp[i][ 8] = -0.333024;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Polyethylene, Borated");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.000000E+00;
  comp[i][ 1] = -0.125355;
  comp[i][ 5] = -0.100000;
  comp[i][ 6] = -0.774645;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Polyethylene, Non-borated");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -9.300000E-01;
  comp[i][ 1] = -0.143716;
  comp[i][ 6] = -0.856284;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Polyisocyanurate (PIR)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -4.820000E-02;
  comp[i][ 1] = -0.040277 ;
  comp[i][ 6] = -0.719916 ;
  comp[i][ 7] = -0.111941 ;
  comp[i][ 8] = -0.127866 ;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Polypropylene (PP)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -9.000000E-01;
  comp[i][ 1] = -0.143711;
  comp[i][ 6] = -0.856289;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Polystyrene (PS)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.060000E+00;
  comp[i][ 1] = -0.077421;
  comp[i][ 6] = -0.922579;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Polytetrafluoroethylene (PTFE)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.250000E+00;
  comp[i][ 6] = -0.240183;
  comp[i][ 9] = -0.759818;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Polyurethane Foam (PUR)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.100000E-02;
  comp[i][ 1] = -0.041000;
  comp[i][ 6] = -0.544000;
  comp[i][ 7] = -0.121000;
  comp[i][ 8] = -0.294000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Polyvinyl Acetate (PVA)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.190000E+00;
  comp[i][ 1] = -0.070245;
  comp[i][ 6] = -0.558066;
  comp[i][ 8] = -0.371689;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Polyvinyl Chloride (PVC)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.406000E+00;
  comp[i][ 1] = -0.048382;
  comp[i][ 6] = -0.384361;
  comp[i][17] = -0.567257;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Polyvinyl Toluene (PVT)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.032000E+00;
  comp[i][ 1] = -0.085000;
  comp[i][ 6] = -0.915000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Polyvinylidene Chloride (PVDC)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.700000E+00;
  comp[i][ 1] = -0.020793;
  comp[i][ 6] = -0.247793;
  comp[i][17] = -0.731413;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Potassium Aluminum Silicate");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.100000E+00;
  comp[i][ 8] = -0.459866 ;
  comp[i][13] = -0.096940 ;
  comp[i][14] = -0.302720 ;
  comp[i][19] = -0.140474 ;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Potassium Iodide");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -3.130000E+00;
  comp[i][19] = -0.235528;
  comp[i][53] = -0.764472;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Potassium Oxide");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.320000E+00;
  comp[i][ 8] = -0.169852;
  comp[i][19] = -0.830148;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Propane (Gas)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.879000E-03;
  comp[i][ 1] = -0.182855;
  comp[i][ 6] = -0.817145;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Propane (Liquid)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -4.300000E-01;
  comp[i][ 1] = -0.182855;
  comp[i][ 6] = -0.817145;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "P-terphenyl");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.230000E+00;
  comp[i][ 1] = -0.056553;
  comp[i][ 6] = -0.943447;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Radiochromic Dye Film, Nylon Base (RDF: NB)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.080000E+00;
  comp[i][ 1] = -0.101996;
  comp[i][ 6] = -0.654396;
  comp[i][ 7] = -0.098915;
  comp[i][ 8] = -0.144693;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Rock (Average of 5 Types)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.662000E+00;
  comp[i][ 1] = -0.001657;
  comp[i][ 6] = -0.026906;
  comp[i][ 8] = -0.488149;
  comp[i][11] = -0.012403;
  comp[i][12] = -0.023146;
  comp[i][13] = -0.054264;
  comp[i][14] = -0.246249;
  comp[i][16] = -0.000577;
  comp[i][19] = -0.018147;
  comp[i][20] = -0.089863;
  comp[i][22] = -0.003621;
  comp[i][25] = -0.000386;
  comp[i][26] = -0.033377;
  comp[i][82] = -0.001255;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Rock, Basalt");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -3.010000E+00;
  comp[i][ 8] = -0.441115 ;
  comp[i][11] = -0.021700 ;
  comp[i][12] = -0.041878 ;
  comp[i][13] = -0.083934 ;
  comp[i][14] = -0.232811 ;
  comp[i][19] = -0.008920 ;
  comp[i][20] = -0.068973 ;
  comp[i][22] = -0.011151 ;
  comp[i][25] = -0.001541 ;
  comp[i][26] = -0.085141 ;
  comp[i][82] = -0.002835 ;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Rock, Granite");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.690000E+00;
  comp[i][ 8] = -0.484170;
  comp[i][11] = -0.027328;
  comp[i][12] = -0.004274;
  comp[i][13] = -0.076188;
  comp[i][14] = -0.336169;
  comp[i][19] = -0.034144;
  comp[i][20] = -0.012985;
  comp[i][22] = -0.001795;
  comp[i][25] = -0.000387;
  comp[i][26] = -0.021555;
  comp[i][82] = -0.001004;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Rock, Limestone");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.610000E+00;
  comp[i][ 1] = -0.000899 ;
  comp[i][ 6] = -0.113782 ;
  comp[i][ 8] = -0.497802 ;
  comp[i][11] = -0.000373 ;
  comp[i][12] = -0.047860 ;
  comp[i][13] = -0.004254 ;
  comp[i][14] = -0.024419 ;
  comp[i][16] = -0.000201 ;
  comp[i][19] = -0.000334 ;
  comp[i][20] = -0.305865 ;
  comp[i][22] = -0.000361 ;
  comp[i][26] = -0.003513 ;
  comp[i][82] = -0.000337 ;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Rock, Sandstone");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.320000E+00;
  comp[i][ 1] = -0.001791;
  comp[i][ 6] = -0.013652;
  comp[i][ 8] = -0.519609;
  comp[i][11] = -0.002969;
  comp[i][12] = -0.007240;
  comp[i][13] = -0.025417;
  comp[i][14] = -0.366185;
  comp[i][16] = -0.000280;
  comp[i][19] = -0.011628;
  comp[i][20] = -0.039328;
  comp[i][22] = -0.001199;
  comp[i][26] = -0.010031;
  comp[i][82] = -0.000671;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Rock, Shale");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.680000E+00;
  comp[i][ 1] = -0.005597;
  comp[i][ 6] = -0.007098;
  comp[i][ 8] = -0.498049;
  comp[i][11] = -0.009647;
  comp[i][12] = -0.014477;
  comp[i][13] = -0.081529;
  comp[i][14] = -0.271661;
  comp[i][16] = -0.002404;
  comp[i][19] = -0.035707;
  comp[i][20] = -0.022162;
  comp[i][22] = -0.003597;
  comp[i][26] = -0.046646;
  comp[i][82] = -0.001425;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Rubber, Butyl");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -9.200000E-01;
  comp[i][ 1] = -0.143711;
  comp[i][ 6] = -0.856289;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Rubber, Natural");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -9.200000E-01;
  comp[i][ 1] = -0.118371;
  comp[i][ 6] = -0.881629;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Rubber, Neoprene");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.230000E+00;
  comp[i][ 1] = -0.056920;
  comp[i][ 6] = -0.542646;
  comp[i][17] = -0.400434;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Rubber, Silicon");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.018500E+00;
  comp[i][ 1] = -0.080716;
  comp[i][ 6] = -0.321164;
  comp[i][ 8] = -0.223545;
  comp[i][14] = -0.374575;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Salt Water (T = 0°C)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.209865E+00;
  comp[i][ 1] = -0.082491;
  comp[i][ 8] = -0.654709;
  comp[i][11] = -0.103378;
  comp[i][17] = -0.159422;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Salt Water (T = 20°C)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.022394E+00;
  comp[i][ 1] = -0.108114;
  comp[i][ 8] = -0.858069;
  comp[i][11] = -0.013302;
  comp[i][17] = -0.020514;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Sand");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.700000E+00;
  comp[i][ 1] = -0.007833;
  comp[i][ 6] = -0.003360;
  comp[i][ 8] = -0.536153;
  comp[i][11] = -0.017063;
  comp[i][13] = -0.034401;
  comp[i][14] = -0.365067;
  comp[i][19] = -0.011622;
  comp[i][20] = -0.011212;
  comp[i][26] = -0.013289;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Sea Water, Simple Artificial");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.023343E+00;
  comp[i][ 1] = -0.107974;
  comp[i][ 8] = -0.858765;
  comp[i][11] = -0.010785;
  comp[i][12] = -0.001284;
  comp[i][16] = -0.000906;
  comp[i][17] = -0.019472;
  comp[i][19] = -0.000399;
  comp[i][20] = -0.000415;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Sea Water, Standard");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.023343E+00;
  comp[i][ 1] = -0.107979 ;
  comp[i][ 5] = -0.000005 ;
  comp[i][ 8] = -0.858803 ;
  comp[i][ 9] = -0.000001 ;
  comp[i][11] = -0.010784 ;
  comp[i][12] = -0.001284 ;
  comp[i][16] = -0.000905 ;
  comp[i][17] = -0.019352 ;
  comp[i][19] = -0.000399 ;
  comp[i][20] = -0.000412 ;
  comp[i][35] = -0.000067 ;
  comp[i][38] = -0.000008 ;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Sepiolite");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.140000E+00;
  comp[i][ 1] = -0.021782;
  comp[i][ 8] = -0.568029;
  comp[i][12] = -0.150070;
  comp[i][14] = -0.260119;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Silicon");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.330000E+00;
  comp[i][14] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Silicon Carbide (Hexagonal)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -3.210000E+00;
  comp[i][ 6] = -0.299547;
  comp[i][14] = -0.700453;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Silicon Dioxide (Alpha-quartz)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.648000E+00;
  comp[i][ 8] = -0.532565;
  comp[i][14] = -0.467435;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Silicon Dioxide (Silica)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.320000E+00;
  comp[i][ 8] = -0.532565;
  comp[i][14] = -0.467435;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Silver");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.050000E+01;
  comp[i][47] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Skin (ICRP)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.100000E+00;
  comp[i][ 1] = -0.100588;
  comp[i][ 6] = -0.228250;
  comp[i][ 7] = -0.046420;
  comp[i][ 8] = -0.619002;
  comp[i][11] = -0.000070;
  comp[i][12] = -0.000060;
  comp[i][15] = -0.000330;
  comp[i][16] = -0.001590;
  comp[i][17] = -0.002670;
  comp[i][19] = -0.000850;
  comp[i][20] = -0.000150;
  comp[i][26] = -0.000010;
  comp[i][30] = -0.000010;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Sodium");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -9.710000E-01;
  comp[i][11] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Sodium Bismuth Tungstate (NBWO)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.570000E+00;
  comp[i][ 8] = -0.175903;
  comp[i][11] = -0.031595;
  comp[i][74] = -0.505301;
  comp[i][83] = -0.287201;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Sodium Chloride");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.170000E+00;
  comp[i][11] = -0.393372;
  comp[i][17] = -0.606628;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Sodium Iodide");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -3.667000E+00;
  comp[i][11] = -0.153373;
  comp[i][53] = -0.846627;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Sodium Nitrate");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.261000E+00;
  comp[i][ 7] = -0.164795;
  comp[i][ 8] = -0.564720;
  comp[i][11] = -0.270485;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Sodium Oxide");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.270000E+00;
  comp[i][ 8] = -0.258143;
  comp[i][11] = -0.741857;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Steel, Boron Stainless");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.870000E+00;
  comp[i][ 5] = -0.010000;
  comp[i][ 6] = -0.000396;
  comp[i][14] = -0.004950;
  comp[i][15] = -0.000228;
  comp[i][16] = -0.000149;
  comp[i][24] = -0.188100;
  comp[i][25] = -0.009900;
  comp[i][26] = -0.694713;
  comp[i][28] = -0.091575;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Steel, Carbon");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.820000E+00;
  comp[i][ 6] = -0.005000;
  comp[i][26] = -0.995000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Steel, HT9 Stainless");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.874000E+00;
  comp[i][ 6] = -0.002000 ;
  comp[i][14] = -0.004000 ;
  comp[i][15] = -0.000300 ;
  comp[i][16] = -0.000200 ;
  comp[i][23] = -0.003000 ;
  comp[i][24] = -0.115000 ;
  comp[i][25] = -0.006000 ;
  comp[i][26] = -0.849500 ;
  comp[i][28] = -0.005000 ;
  comp[i][42] = -0.010000 ;
  comp[i][74] = -0.005000 ;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Steel, Stainless 202");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.860000E+00;
  comp[i][ 6] = -0.000750;
  comp[i][ 7] = -0.001250;
  comp[i][14] = -0.005000;
  comp[i][15] = -0.000300;
  comp[i][16] = -0.000150;
  comp[i][24] = -0.180000;
  comp[i][25] = -0.087500;
  comp[i][26] = -0.675050;
  comp[i][28] = -0.050000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Steel, Stainless 302");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.860000E+00;
  comp[i][ 6] = -0.001400;
  comp[i][14] = -0.009300;
  comp[i][15] = -0.000420;
  comp[i][16] = -0.000280;
  comp[i][24] = -0.180000;
  comp[i][25] = -0.018600;
  comp[i][26] = -0.700000;
  comp[i][28] = -0.090000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Steel, Stainless 304");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -8.000000E+00;
  comp[i][ 6] = -0.000400;
  comp[i][14] = -0.005000;
  comp[i][15] = -0.000230;
  comp[i][16] = -0.000150;
  comp[i][24] = -0.190000;
  comp[i][25] = -0.010000;
  comp[i][26] = -0.701730;
  comp[i][28] = -0.092500;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Steel, Stainless 304L");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -8.000000E+00;
  comp[i][ 6] = -0.000150;
  comp[i][14] = -0.005000;
  comp[i][15] = -0.000230;
  comp[i][16] = -0.000150;
  comp[i][24] = -0.190000;
  comp[i][25] = -0.010000;
  comp[i][26] = -0.694480;
  comp[i][28] = -0.100000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Steel, Stainless 316");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -8.000000E+00;
  comp[i][ 6] = -0.000410 ;
  comp[i][14] = -0.005070 ;
  comp[i][15] = -0.000230 ;
  comp[i][16] = -0.000150 ;
  comp[i][24] = -0.170000 ;
  comp[i][25] = -0.010140 ;
  comp[i][26] = -0.669000 ;
  comp[i][28] = -0.120000 ;
  comp[i][42] = -0.025000 ;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Steel, Stainless 316L");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -8.000000E+00;
  comp[i][ 6] = -0.000300;
  comp[i][14] = -0.010000;
  comp[i][15] = -0.000450;
  comp[i][16] = -0.000300;
  comp[i][24] = -0.170000;
  comp[i][25] = -0.020000;
  comp[i][26] = -0.653950;
  comp[i][28] = -0.120000;
  comp[i][42] = -0.025000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Steel, Stainless 321");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -8.000000E+00;
  comp[i][ 6] = -0.000800;
  comp[i][14] = -0.010000;
  comp[i][15] = -0.000450;
  comp[i][16] = -0.000300;
  comp[i][22] = -0.001500;
  comp[i][24] = -0.180000;
  comp[i][25] = -0.020000;
  comp[i][26] = -0.676950;
  comp[i][28] = -0.110000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Steel, Stainless 347");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -8.000000E+00;
  comp[i][ 6] = -0.000800 ;
  comp[i][14] = -0.010000 ;
  comp[i][15] = -0.000450 ;
  comp[i][16] = -0.000300 ;
  comp[i][24] = -0.170000 ;
  comp[i][25] = -0.020000 ;
  comp[i][26] = -0.680450 ;
  comp[i][28] = -0.110000 ;
  comp[i][41] = -0.004000 ;
  comp[i][73] = -0.004000 ;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Steel, Stainless 409");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.800000E+00;
  comp[i][ 6] = -0.000790;
  comp[i][14] = -0.009830;
  comp[i][15] = -0.000440;
  comp[i][16] = -0.000440;
  comp[i][22] = -0.007370;
  comp[i][24] = -0.111300;
  comp[i][25] = -0.009830;
  comp[i][26] = -0.860000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Steel, Stainless 440");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.800000E+00;
  comp[i][ 6] = -0.006750;
  comp[i][14] = -0.006500;
  comp[i][15] = -0.000260;
  comp[i][16] = -0.000200;
  comp[i][24] = -0.170000;
  comp[i][25] = -0.006500;
  comp[i][26] = -0.795050;
  comp[i][42] = -0.004880;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Sterotex");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -8.620000E-01;
  comp[i][ 1] = -0.124370;
  comp[i][ 6] = -0.767948;
  comp[i][ 8] = -0.107682;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Stilbene (Trans-stilbene Isomer)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.220000E+00;
  comp[i][ 1] = -0.056553 ;
  comp[i][ 6] = -0.943447 ;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Sulphur");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.000000E+00;
  comp[i][16] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Tantalum");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.665400E+01;
  comp[i][73] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Thorium");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.172000E+01;
  comp[i][90] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Thorium Dioxide");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.000000E+01;
  comp[i][ 8] = -0.121191;
  comp[i][90] = -0.878809;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Tin");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.310000E+00;
  comp[i][50] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Tissue Equivalent, MS20");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.000000E+00;
  comp[i][ 1] = -0.081192;
  comp[i][ 6] = -0.583442;
  comp[i][ 7] = -0.017798;
  comp[i][ 8] = -0.186381;
  comp[i][12] = -0.130287;
  comp[i][17] = -0.000900;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Tissue Equivalent-Gas, Methane Based (TEG: MB)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.064000E-03;
  comp[i][ 1] = -0.101869;
  comp[i][ 6] = -0.456179;
  comp[i][ 7] = -0.035172;
  comp[i][ 8] = -0.406780;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Tissue Equivalent-Gas, Propane Based (TEG: PB)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.826000E-03;
  comp[i][ 1] = -0.102672;
  comp[i][ 6] = -0.568940;
  comp[i][ 7] = -0.035022;
  comp[i][ 8] = -0.293366;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Tissue, Adipose (ICRP)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -9.200000E-01;
  comp[i][ 1] = -0.119477;
  comp[i][ 6] = -0.637240;
  comp[i][ 7] = -0.007970;
  comp[i][ 8] = -0.232333;
  comp[i][11] = -0.000500;
  comp[i][12] = -0.000020;
  comp[i][15] = -0.000160;
  comp[i][16] = -0.000730;
  comp[i][17] = -0.001190;
  comp[i][19] = -0.000320;
  comp[i][20] = -0.000020;
  comp[i][26] = -0.000020;
  comp[i][30] = -0.000020;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Tissue, Breast");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.020000E+00;
  comp[i][ 1] = -0.106000;
  comp[i][ 6] = -0.332000;
  comp[i][ 7] = -0.030000;
  comp[i][ 8] = -0.527000;
  comp[i][11] = -0.001000;
  comp[i][15] = -0.001000;
  comp[i][16] = -0.002000;
  comp[i][17] = -0.001000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Tissue, Lung (ICRP)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.050000E+00;
  comp[i][ 1] = -0.101278;
  comp[i][ 6] = -0.102310;
  comp[i][ 7] = -0.028650;
  comp[i][ 8] = -0.757072;
  comp[i][11] = -0.001840;
  comp[i][12] = -0.000730;
  comp[i][15] = -0.000800;
  comp[i][16] = -0.002250;
  comp[i][17] = -0.002660;
  comp[i][19] = -0.001940;
  comp[i][20] = -0.000090;
  comp[i][26] = -0.000370;
  comp[i][30] = -0.000010;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Tissue, Ovary");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.050000E+00;
  comp[i][ 1] = -0.105000;
  comp[i][ 6] = -0.093000;
  comp[i][ 7] = -0.024000;
  comp[i][ 8] = -0.768000;
  comp[i][11] = -0.002000;
  comp[i][15] = -0.002000;
  comp[i][16] = -0.002000;
  comp[i][17] = -0.002000;
  comp[i][19] = -0.002000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Tissue, Soft (ICRP)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.000000E+00;
  comp[i][ 1] = -0.104472 ;
  comp[i][ 6] = -0.232190 ;
  comp[i][ 7] = -0.024880 ;
  comp[i][ 8] = -0.630238 ;
  comp[i][11] = -0.001130 ;
  comp[i][12] = -0.000130 ;
  comp[i][15] = -0.001330 ;
  comp[i][16] = -0.001990 ;
  comp[i][17] = -0.001340 ;
  comp[i][19] = -0.001990 ;
  comp[i][20] = -0.000230 ;
  comp[i][26] = -0.000050 ;
  comp[i][30] = -0.000030 ;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Tissue, Soft (ICRU Four Component)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.000000E+00;
  comp[i][ 1] = -0.101172;
  comp[i][ 6] = -0.111000;
  comp[i][ 7] = -0.026000;
  comp[i][ 8] = -0.761828;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Tissue, Testes (ICRP)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.040000E+00;
  comp[i][ 1] = -0.104166;
  comp[i][ 6] = -0.092270;
  comp[i][ 7] = -0.019940;
  comp[i][ 8] = -0.773884;
  comp[i][11] = -0.002260;
  comp[i][12] = -0.000110;
  comp[i][15] = -0.001250;
  comp[i][16] = -0.001460;
  comp[i][17] = -0.002440;
  comp[i][19] = -0.002080;
  comp[i][20] = -0.000100;
  comp[i][26] = -0.000020;
  comp[i][30] = -0.000020;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Tissue, Testis (ICRU)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.040000E+00;
  comp[i][ 1] = -0.106000;
  comp[i][ 6] = -0.099000;
  comp[i][ 7] = -0.020000;
  comp[i][ 8] = -0.766000;
  comp[i][11] = -0.002000;
  comp[i][15] = -0.001000;
  comp[i][16] = -0.002000;
  comp[i][17] = -0.002000;
  comp[i][19] = -0.002000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Titanium");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -4.540000E+00;
  comp[i][22] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Titanium Alloy, Grade 5");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -4.430000E+00;
  comp[i][ 1] = -0.000110;
  comp[i][ 6] = -0.000570;
  comp[i][ 7] = -0.000210;
  comp[i][ 8] = -0.001410;
  comp[i][13] = -0.061250;
  comp[i][22] = -0.893630;
  comp[i][23] = -0.040000;
  comp[i][26] = -0.002830;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Titanium Dioxide");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -4.260000E+00;
  comp[i][ 8] = -0.400592;
  comp[i][22] = -0.599408;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Titanium Hydride");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -3.750000E+00;
  comp[i][ 1] = -0.040412;
  comp[i][22] = -0.959588;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Toluene");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -8.669000E-01;
  comp[i][ 1] = -0.087510;
  comp[i][ 6] = -0.912490;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Tributyl Borate");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -8.640000E-01;
  comp[i][ 1] = -0.118245;
  comp[i][ 5] = -0.046973;
  comp[i][ 6] = -0.626231;
  comp[i][ 8] = -0.208550;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Tributyl Phosphate (TBP)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -9.724000E-01;
  comp[i][ 1] = -0.102189 ;
  comp[i][ 6] = -0.541197 ;
  comp[i][ 8] = -0.240309 ;
  comp[i][15] = -0.116305 ;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Tungsten");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.930000E+01;
  comp[i][74] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Vermiculite, Exfoliated");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -8.500000E-02;
  comp[i][ 1] = -0.011835;
  comp[i][ 8] = -0.496356;
  comp[i][12] = -0.133383;
  comp[i][13] = -0.063151;
  comp[i][14] = -0.189668;
  comp[i][19] = -0.021668;
  comp[i][20] = -0.016353;
  comp[i][22] = -0.009854;
  comp[i][26] = -0.057732;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Viton Fluoroelastomer");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.800000E+00;
  comp[i][ 1] = -0.009417;
  comp[i][ 6] = -0.280555;
  comp[i][ 9] = -0.710028;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Water, Heavy");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.105340E+00;
  comp[i][ 1] = -0.201133;
  comp[i][ 8] = -0.798867;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Water, Liquid");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -9.982070E-01;
  comp[i][ 1] = -0.111894;
  comp[i][ 8] = -0.888106;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Water, Vapor");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.560000E-04;
  comp[i][ 1] = -0.111894;
  comp[i][ 8] = -0.888106;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Wax, M3");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -1.050000E+00;
  comp[i][ 1] = -0.114318;
  comp[i][ 6] = -0.655823;
  comp[i][ 8] = -0.092183;
  comp[i][12] = -0.134792;
  comp[i][20] = -0.002883;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Wax, Mix D");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -9.900000E-01;
  comp[i][ 1] = -0.134040;
  comp[i][ 6] = -0.777960;
  comp[i][ 8] = -0.035020;
  comp[i][12] = -0.038594;
  comp[i][22] = -0.014386;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Wax, Paraffin");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -9.300000E-01;
  comp[i][ 1] = -0.148605;
  comp[i][ 6] = -0.851395;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Wood, Ash (black)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -0.55;
  comp[i][ 1] = -0.059642;
  comp[i][ 6] = -0.497018;
  comp[i][ 7] = -0.004970;
  comp[i][ 8] = -0.427435;
  comp[i][12] = -0.001988;
  comp[i][16] = -0.004970;
  comp[i][19] = -0.001988;
  comp[i][20] = -0.001988;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Wood, Ash (white)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -0.67;
  comp[i][ 1] = -0.059642;
  comp[i][ 6] = -0.497018;
  comp[i][ 7] = -0.004970;
  comp[i][ 8] = -0.427435;
  comp[i][12] = -0.001988;
  comp[i][16] = -0.004970;
  comp[i][19] = -0.001988;
  comp[i][20] = -0.001988;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Wood, Balsa");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -0.12;
  comp[i][ 1] = -0.059642;
  comp[i][ 6] = -0.497018;
  comp[i][ 7] = -0.004970;
  comp[i][ 8] = -0.427435;
  comp[i][12] = -0.001988;
  comp[i][16] = -0.004970;
  comp[i][19] = -0.001988;
  comp[i][20] = -0.001988;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Wood, Birch");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -0.71;
  comp[i][ 1] = -0.059642;
  comp[i][ 6] = -0.497018;
  comp[i][ 7] = -0.004970;
  comp[i][ 8] = -0.427435;
  comp[i][12] = -0.001988;
  comp[i][16] = -0.004970;
  comp[i][19] = -0.001988;
  comp[i][20] = -0.001988;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Wood, Ceder");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -0.35;
  comp[i][ 1] = -0.059642;
  comp[i][ 6] = -0.497018;
  comp[i][ 7] = -0.004970;
  comp[i][ 8] = -0.427435;
  comp[i][12] = -0.001988;
  comp[i][16] = -0.004970;
  comp[i][19] = -0.001988;
  comp[i][20] = -0.001988;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Wood, Cherry");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -0.43;
  comp[i][ 1] = -0.059642;
  comp[i][ 6] = -0.497018;
  comp[i][ 7] = -0.004970;
  comp[i][ 8] = -0.427435;
  comp[i][12] = -0.001988;
  comp[i][16] = -0.004970;
  comp[i][19] = -0.001988;
  comp[i][20] = -0.001988;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Wood, Fir (douglas)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -0.51;
  comp[i][ 1] = -0.059642;
  comp[i][ 6] = -0.497018;
  comp[i][ 7] = -0.004970;
  comp[i][ 8] = -0.427435;
  comp[i][12] = -0.001988;
  comp[i][16] = -0.004970;
  comp[i][19] = -0.001988;
  comp[i][20] = -0.001988;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Wood, Elm");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -0.56;
  comp[i][ 1] = -0.059642;
  comp[i][ 6] = -0.497018;
  comp[i][ 7] = -0.004970;
  comp[i][ 8] = -0.427435;
  comp[i][12] = -0.001988;
  comp[i][16] = -0.004970;
  comp[i][19] = -0.001988;
  comp[i][20] = -0.001988;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Wood, Hickory");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -0.77;
  comp[i][ 1] = -0.059642;
  comp[i][ 6] = -0.497018;
  comp[i][ 7] = -0.004970;
  comp[i][ 8] = -0.427435;
  comp[i][12] = -0.001988;
  comp[i][16] = -0.004970;
  comp[i][19] = -0.001988;
  comp[i][20] = -0.001988;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Wood, Mahogany");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -0.70;
  comp[i][ 1] = -0.059642;
  comp[i][ 6] = -0.497018;
  comp[i][ 7] = -0.004970;
  comp[i][ 8] = -0.427435;
  comp[i][12] = -0.001988;
  comp[i][16] = -0.004970;
  comp[i][19] = -0.001988;
  comp[i][20] = -0.001988;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Wood, Maple (sugar)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -0.68;
  comp[i][ 1] = -0.059642;
  comp[i][ 6] = -0.497018;
  comp[i][ 7] = -0.004970;
  comp[i][ 8] = -0.427435;
  comp[i][12] = -0.001988;
  comp[i][16] = -0.004970;
  comp[i][19] = -0.001988;
  comp[i][20] = -0.001988;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Wood, Maple (white)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -0.53;
  comp[i][ 1] = -0.059642;
  comp[i][ 6] = -0.497018;
  comp[i][ 7] = -0.004970;
  comp[i][ 8] = -0.427435;
  comp[i][12] = -0.001988;
  comp[i][16] = -0.004970;
  comp[i][19] = -0.001988;
  comp[i][20] = -0.001988;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Wood, Oak (black or red)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -0.67;
  comp[i][ 1] = -0.059642;
  comp[i][ 6] = -0.497018;
  comp[i][ 7] = -0.004970;
  comp[i][ 8] = -0.427435;
  comp[i][12] = -0.001988;
  comp[i][16] = -0.004970;
  comp[i][19] = -0.001988;
  comp[i][20] = -0.001988;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Wood, Oak (white)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -0.77;
  comp[i][ 1] = -0.059642;
  comp[i][ 6] = -0.497018;
  comp[i][ 7] = -0.004970;
  comp[i][ 8] = -0.427435;
  comp[i][12] = -0.001988;
  comp[i][16] = -0.004970;
  comp[i][19] = -0.001988;
  comp[i][20] = -0.001988;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Wood, Pine (white)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -0.43;
  comp[i][ 1] = -0.059642;
  comp[i][ 6] = -0.497018;
  comp[i][ 7] = -0.004970;
  comp[i][ 8] = -0.427435;
  comp[i][12] = -0.001988;
  comp[i][16] = -0.004970;
  comp[i][19] = -0.001988;
  comp[i][20] = -0.001988;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Wood, Pine (yellow)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -0.71;
  comp[i][ 1] = -0.059642;
  comp[i][ 6] = -0.497018;
  comp[i][ 7] = -0.004970;
  comp[i][ 8] = -0.427435;
  comp[i][12] = -0.001988;
  comp[i][16] = -0.004970;
  comp[i][19] = -0.001988;
  comp[i][20] = -0.001988;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Wood, Poplar");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -0.43;
  comp[i][ 1] = -0.059642;
  comp[i][ 6] = -0.497018;
  comp[i][ 7] = -0.004970;
  comp[i][ 8] = -0.427435;
  comp[i][12] = -0.001988;
  comp[i][16] = -0.004970;
  comp[i][19] = -0.001988;
  comp[i][20] = -0.001988;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Wood, Redwood");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -0.42;
  comp[i][ 1] = -0.059642;
  comp[i][ 6] = -0.497018;
  comp[i][ 7] = -0.004970;
  comp[i][ 8] = -0.427435;
  comp[i][12] = -0.001988;
  comp[i][16] = -0.004970;
  comp[i][19] = -0.001988;
  comp[i][20] = -0.001988;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Wood, Spruce");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -0.45;
  comp[i][ 1] = -0.059642;
  comp[i][ 6] = -0.497018;
  comp[i][ 7] = -0.004970;
  comp[i][ 8] = -0.427435;
  comp[i][12] = -0.001988;
  comp[i][16] = -0.004970;
  comp[i][19] = -0.001988;
  comp[i][20] = -0.001988;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Wood, Southern Pine");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -6.400000E-01;
  comp[i][ 1] = -0.059642;
  comp[i][ 6] = -0.497018;
  comp[i][ 7] = -0.004970;
  comp[i][ 8] = -0.427435;
  comp[i][12] = -0.001988;
  comp[i][16] = -0.004970;
  comp[i][19] = -0.001988;
  comp[i][20] = -0.001988;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Wood, Walnut");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -0.59;
  comp[i][ 1] = -0.059642;
  comp[i][ 6] = -0.497018;
  comp[i][ 7] = -0.004970;
  comp[i][ 8] = -0.427435;
  comp[i][12] = -0.001988;
  comp[i][16] = -0.004970;
  comp[i][19] = -0.001988;
  comp[i][20] = -0.001988;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Wood, Coarse sawdust");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -0.29;
  comp[i][ 1] = -0.059642;
  comp[i][ 6] = -0.497018;
  comp[i][ 7] = -0.004970;
  comp[i][ 8] = -0.427435;
  comp[i][12] = -0.001988;
  comp[i][16] = -0.004970;
  comp[i][19] = -0.001988;
  comp[i][20] = -0.001988;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Wood, Fine sawdust");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -0.40;
  comp[i][ 1] = -0.059642;
  comp[i][ 6] = -0.497018;
  comp[i][ 7] = -0.004970;
  comp[i][ 8] = -0.427435;
  comp[i][12] = -0.001988;
  comp[i][16] = -0.004970;
  comp[i][19] = -0.001988;
  comp[i][20] = -0.001988;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Xenon");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -5.485000E-03;
  comp[i][54] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Yttrium Aluminum Garnet (YAG)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -4.560000E+00;
  comp[i][ 8] = -0.323428;
  comp[i][13] = -0.227263;
  comp[i][39] = -0.449308;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Yttrium Aluminum Perovskite (YAP)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -5.370000E+00;
  comp[i][ 8] = -0.292876;
  comp[i][13] = -0.164636;
  comp[i][39] = -0.542487;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Yttrium OxyorthoSilicate (YSO)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -4.450000E+00;
  comp[i][ 8] = -0.279813;
  comp[i][14] = -0.098237;
  comp[i][39] = -0.621949;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Zeolite (Natrolite)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -2.250000E+00;
  comp[i][ 1] = -0.010604;
  comp[i][ 8] = -0.504947;
  comp[i][11] = -0.120928;
  comp[i][13] = -0.141925;
  comp[i][14] = -0.221597;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Zinc");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -7.133000E+00;
  comp[i][30] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Zinc Selenide");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -5.420000E+00;
  comp[i][30] = -0.453068;
  comp[i][34] = -0.546932;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Zinc Sulfide");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -4.090000E+00;
  comp[i][16] = -0.328960;
  comp[i][30] = -0.671040;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Zircaloy-2");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -6.560000E+00;
  comp[i][ 8] = -0.001197;
  comp[i][24] = -0.000997;
  comp[i][26] = -0.000997;
  comp[i][28] = -0.000499;
  comp[i][40] = -0.982348;
  comp[i][50] = -0.013962;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Zircaloy-4");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -6.560000E+00;
  comp[i][ 8] = -0.001196;
  comp[i][24] = -0.000997;
  comp[i][26] = -0.001994;
  comp[i][40] = -0.981858;
  comp[i][50] = -0.013955;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Zirconium");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -6.506000E+00;
  comp[i][40] = -1.000000;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Zirconium Hydride (Zr5H8)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -5.610000E+00;
  comp[i][ 1] = -0.017371;
  comp[i][40] = -0.982629;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "Zirconium Hydride (ZrH2)");
  sprintf(src[i], "PNNL-15870, Rev. 1");
  dens[i] = -5.610000E+00;
  comp[i][ 1] = -0.021620;
  comp[i][40] = -0.978380;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "316NG");
  sprintf(src[i], "VTT (S. Penttila)");
  dens[i] = -6.560000E+00;
  comp[i][ 6] = -1.40000E-04;
  comp[i][14] = -4.20000E-03;
  comp[i][25] = -8.00000E-03;
  comp[i][16] = -1.00000E-05;
  comp[i][15] = -3.20000E-04;
  comp[i][24] = -1.66000E-01;
  comp[i][28] = -1.13000E-01;
  comp[i][42] = -2.11000E-02;
  comp[i][29] = -2.30000E-03;
  comp[i][13] = -7.00000E-05;
  comp[i][74] = -2.00000E-04;
  comp[i][23] = -4.00000E-04;
  comp[i][22] = -4.00000E-05;
  comp[i][27] = -7.00000E-04;
  comp[i][41] = -1.00000E-04;
  comp[i][26] = -6.83420E-01;
  i++;

  sprintf(alias[i], " ");
  sprintf(comment[i], "800H");
  sprintf(src[i], "VTT (S. Penttila)");
  dens[i] = -6.560000E+00;
  comp[i][ 6] = -6.00000E-04;
  comp[i][14] = -3.60000E-03;
  comp[i][25] = -6.70000E-03;
  comp[i][16] = -3.00000E-05;
  comp[i][15] = -1.00000E-04;
  comp[i][24] = -2.05000E-01;
  comp[i][28] = -3.08000E-01;
  comp[i][42] = -1.30000E-03;
  comp[i][29] = -6.00000E-04;
  comp[i][13] = -2.60000E-03;
  comp[i][74] = -1.00000E-04;
  comp[i][23] = -4.00000E-04;
  comp[i][22] = -3.60000E-03;
  comp[i][27] = -9.00000E-04;
  comp[i][41] = -1.00000E-04;
  comp[i][26] = -4.66370E-01;
  i++;

  /* Get total and check */

  if ((itot = i) > MAX_STD_COMP)
    Die(FUNCTION_NAME, "Too many compositions, increase maximum");

  /* Check that compositions are consistent */

  for (i = 0; i < itot; i++)
    {
      /* Reset previous index */

      m = -1;

      /* Loop over composition */

      for (n = 0; n < 112; n++)
        {
          if (comp[i][n] != 0.0)
            {
              /* Check that both entries are either positive or negative */

              if (m != -1)
                if (comp[i][n]*comp[i][m] < 0.0)
                  Die(FUNCTION_NAME, "Error in composition \"%s\"",
                      comment[i]);

              /* Remember previous index */

              m = n;
            }
        }
    }

  /* Check listing */

  if (!strcasecmp(mat, "list"))
    {
      /* Print all */

      fprintf(outp, "\nList of available material compositions:\n\n");

      for (i = 0; i < itot; i++)
        fprintf(outp, "%3ld  \"%s\"\n", i + 1, comment[i]);

      fprintf(outp, "\n");

      /* Exit subroutine */

      return;
    }

  /* Match name */

  for (i = 0; i < itot; i++)
    if (!strcasecmp(mat, alias[i]))
      {
        sprintf(name, "%s", mat);
        break;
      }

  /* Match number */

  if (i == itot)
    for (i = 0; i < itot; i++)
      {
        sprintf(name, "%ld", i + 1);
        if (!strcasecmp(mat, name))
          {
            sprintf(name, "m%ld", i + 1);
            break;
          }
      }

  /* Check if found */

  if (i == itot)
    {
      /* Nope */

      printf("\nComposition %s is not defined\n\n", mat);

      /* Exit subroutine */

      return;
    }

  /***************************************************************************/

  /***** Convert atomic to mass fractions ************************************/

  /* Check type of composition */

  for (n = 0; n < 112; n++)
    if (comp[i][n] != 0.0)
      break;

  /* Reset sum */

  for (n = 0; n < 112; n++)
    sum[n] = 0.0;

  /* Calculate sums and products */

  m = 0;
  do
    {
      /* Get atomic number */

      n = (long)nat_frac[m][0];

      /* Multiply fraction and mass if mass densities given */

      if (comp[i][n] < 0.0)
        nat_frac[m][3] = nat_frac[m][3]*nat_frac[m][2];

      /* Add to total */

      sum[n] = sum[n] + nat_frac[m][3];
    }
  while ((long)nat_frac[++m][0] != -1);

  /* Divide by total */

  m = 0;
  do
    {
      /* Get atomic number */

      n = (long)nat_frac[m][0];

      /* Multiply fraction and mass */

      if (sum[n] > 0.0)
        nat_frac[m][3] = nat_frac[m][3]/sum[n];
      else
        Die(FUNCTION_NAME, "WTF?");
    }
  while ((long)nat_frac[++m][0] != -1);

  /***************************************************************************/

  /***** Print output ********************************************************/

  /* Print name and comments */

  fprintf(outp, "\n%% --- \"%s\" [%s]\n\n", comment[i], src[i]);
  fprintf(outp, "mat %s %1.5E\n\n", name, dens[i]);

  /* Print elemental composition */

  for (n = 0; n < 112; n++)
    if (comp[i][n] != 0.0)
      {
        if (id == NULL)
          fprintf(outp, "%5ld  %1.5E\n", 1000*n, comp[i][n]);
        else
          fprintf(outp, "%5ld.%s  %1.5E\n", 1000*n, id, comp[i][n]);
      }

  fprintf(outp, "\n");

  /* Print isotopic composition */

  for (n = 0; n < 112; n++)
    if (comp[i][n] != 0.0)
      {
        /* Loop over natural composition */

        m = 0;
        do {
          if ((long)nat_frac[m][0] == n)
            {
              /* Get A */

              A = (long)nat_frac[m][1];

              /* Ta-180 is actually isomeric */

              if ((n == 73) && (A == 180))
                A = A + 200;

              /* Check if is given */

              if (id == NULL)
                fprintf(outp, "%5ld  %1.5E\n", 1000*n + A,
                        comp[i][n]*nat_frac[m][3]);
              else
                fprintf(outp, "%5ld.%s  %1.5E\n", 1000*n + A,
                        id, comp[i][n]*nat_frac[m][3]);
            }
          m++;
        }
        while ((long)nat_frac[m][0] != -1);
      }

  fprintf(outp, "\n");

  /***************************************************************************/
}

/*****************************************************************************/
