#region Copyright
///<remarks>
/// <GRAMM Mesoscale Model>
/// Copyright (C) [2019]  [Dietmar Oettl, Markus Kuntner]
/// This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
/// the Free Software Foundation version 3 of the License
/// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
/// You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
///</remarks>
#endregion

using System;
using System.Threading.Tasks;
using System.Runtime.CompilerServices;

namespace GRAMM_2001
{
    partial class Program
    {
        ///<summary>
        ///Calculate the diffusion and advection terms for the implicit scheme (Patankar 1980, p52)
        ///</summary>
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        public static void TERMIVterms(int NI, int NJ, int NK)
        {

            Parallel.For(1, NI + 1, Program.pOptions, i =>
             {
                 double DS; double FS; double PS; double DB; double FB; double PB; double DE; double FE; double PE; double DW1; double FW1; double PW1;
                 double DS2; double FS2; double PS2I; double DT2; double FT2; double PT2; double DE2; double FE2; double PE2; double DN2; double FN2; double PN2;
                 double VISKM, VISKP, DDZP, DDZM, DWDXDUDZDZ, DWDYDVDZDZ;

                 double V2_L_M = 0, VISH = 0, VISH_I = 0, VISH_J = 0, VISH_IM = 0, VISH_JM = 0, VISH_T1 = 0, VISH_T2 = 0;
                 float ZSP_LL = 0, ZSP_L_M = 0, ZSP_L_P = 0;

                 double stretch = Program.STRETCH * 2, stretchRzp = 1 / Program.STRETCH * 2;
                 double VISH_R1 = 0, VISH_R2 = 0, VISH_KM = 0;
                 double AEREA_LL = 0, AEREAZX_LL = 0, AEREAZY_LL = 0, AREA_T1 = 0, AREA_T2 = 0, AREA_T3 = 0, AREA_T4 = 0;
                 double DDX_Rez = 0.0001, DDY_Rez = 0.0001;
                 double RHO = 0, RHO_KP = 0, RHO_L_KM1 = 0, RHOKM, RHOKP;

                 int NK_P = NK; int NJ_P = NJ; int NI_P = NI;

                 for (int j = 1; j <= NJ_P; j++)
                 {
                     float[] AE2_L = Program.AE2[i][j];
                     float[] AN2_L = Program.AN2[i][j];
                     float[] AIM_L = Program.AIM[i][j];
                     ReadOnlySpan <float> AREA_L = Program.AREAImm[i][j];
                     ReadOnlySpan <float> AREAX_L = Program.AREAXImm[i][j];
                     ReadOnlySpan <float> AREAZX_L = Program.AREAZXImm[i][j];
                     ReadOnlySpan <float> AREAZY_L = Program.AREAZYImm[i][j];
                     ReadOnlySpan <float> AREAY_L = Program.AREAYImm[i][j];
                     float[] AP0_L = Program.AP0[i][j];
                     float[] AS1_L = Program.AS1[i][j];
                     float[] AW1_L = Program.AW1[i][j];
                     float[] BIM_L = Program.BIM[i][j];
                     float[] CIM_L = Program.CIM[i][j];
                     double[] DIMU_L = Program.DIMU[i][j];
                     double[] DIMV_L = Program.DIMV[i][j];
                     double[] DIMW_L = Program.DIMW[i][j];
                     ReadOnlySpan <float> RHO_L = Program.RHOImm[i][j];
                     ReadOnlySpan <double> U_L = Program.U[i][j];
                     ReadOnlySpan <double> U1_L = Program.U1[i][j];
                     ReadOnlySpan <double> U2_L = Program.U2[i][j];
                     ReadOnlySpan <double> V_L = Program.V[i][j];
                     ReadOnlySpan <double> V1_L = Program.V1[i][j];
                     ReadOnlySpan <double> V2_L = Program.V2[i][j];
                     ReadOnlySpan <double> W1_L = Program.W1[i][j];
                     ReadOnlySpan <double> W2_L = Program.W2[i][j];
                     double[] VISH_L = Program.VISH[i][j];
                     double[] VISV_L = Program.VISV[i][j];
                     ReadOnlySpan <float> VOL_L = Program.VOLImm[i][j];
                     ReadOnlySpan <float> ZSP_L = Program.ZSPImm[i][j];

                     int m = 2;

                     if ((i > 1) && (j > 1) && (i < NI_P) && (j < NJ_P))
                     {
                         DDX_Rez = 1 / Program.DDXImm[i]; DDY_Rez = 1 / Program.DDYImm[j];
                     }

                     for (int kn = 1; kn <= 2 * NK_P; kn++)
                     {

                         RHOKM = 0; VISKM = 0; RHOKP = 0; VISKP = 0; DDZP = 0; DDZM = 0; DWDXDUDZDZ = 0; DWDYDVDZDZ = 0;

                         int k = 1 + kn >> 1;  // (int) (kn * 0.5 + 0.5)

                         //As array access is a slow process, certain often-used array-values are used as constant
                         if (kn % 2 != 0) // ODD numbers 1,3,5,..- new k value - set local constants new
                         {
                             ZSP_L_M = ZSP_L[k - 1];
                             ZSP_LL = ZSP_L[k];
                             V2_L_M = Program.V2[i][j - 1][k];
                             RHO_L_KM1 = RHO_L[k - 1];

                             VISH = VISH_L[k];
                             RHO = RHO_L[k];
                             RHO_KP = 0;
                             if (k < NK_P)
                             {
                                 RHO_KP = RHO_L[k + 1];
                                 ZSP_L_P = ZSP_L[k + 1];
                             }

                             VISH_I = 0; VISH_J = 0;
                             VISH_IM = 0; VISH_JM = 0;


                             if (i > 1)
                             {
                                 VISH_IM = Program.VISH[i - 1][j][k];
                             }

                             if (i < NI_P)
                             {
                                 VISH_I = Program.VISH[i + 1][j][k];
                             }

                             if (j > 1)
                             {
                                 VISH_JM = Program.VISH[i][j - 1][k];
                             }

                             if (j < NJ_P)
                             {
                                 VISH_J = Program.VISH[i][j + 1][k];
                             }

                             VISH_R1 = 1 / (VISH + VISH_J);
                             VISH_R2 = 1 / (VISH + VISH_I);
                             VISH_KM = VISH_L[k - 1];

                             AEREA_LL = AREA_L[k];
                             if (k < NK_P)
                             {
                                 AEREAZX_LL = AREAZX_L[k + 1];
                                 AEREAZY_LL = AREAZY_L[k + 1];
                             }
                             AREA_T1 = AREAX_L[k] + AREAZX_L[k];
                             AREA_T2 = AREAY_L[k] + AREAZY_L[k];
                             AREA_T3 = Program.AREAXImm[i + 1][j][k];
                             AREA_T4 = Program.AREAYImm[i][j + 1][k];
                             VISH_T1 = 1 / (VISH + VISH_IM);
                             VISH_T2 = 1 / (VISH + VISH_JM);

                             if ((i > 1) && (j > 1) && (i < NI_P) && (j < NJ_P) && (k < NK_P))
                             {


                                 RHOKP = 0.5 * (RHO_KP + RHO);
                                 VISKP = 0.5 * (VISV_L[k + 1] + VISV_L[k]);

                                 //Diffusion terms for the u-component
                                 DDZP = 1 / (ZSP_L_P - ZSP_LL);

                                 if (k > 1)
                                 {
                                     RHOKM = 0.5 * (RHO_L_KM1 + RHO);
                                     VISKM = 0.5 * (VISV_L[k - 1] + VISV_L[k]);
                                     DDZM = 1 / (ZSP_LL - ZSP_L_M);

                                     DWDXDUDZDZ = ((RHOKP * VISKP * (U_L[k + 1] - U_L[k]) * DDZP) * AREA_L[k + 1] -
                                                  (RHOKM * VISKM * (U_L[k] - U_L[k - 1]) * DDZM) * AEREA_LL);
                                 }
                                 else
                                 {
                                     DWDXDUDZDZ = (RHOKP * VISKP * (U_L[k + 1] - U_L[k]) * DDZP) * AREA_L[k + 1];
                                 }

                                 //Diffusion terms for the v-component
                                 if (k > 1)
                                 {
                                     DWDYDVDZDZ = ((RHOKP * VISKP * (V_L[k + 1] - V_L[k]) * DDZP) * AREA_L[k + 1] -
                                                  (RHOKM * VISKM * (V_L[k] - V_L[k - 1]) * DDZM) * AEREA_LL);
                                 }
                                 else
                                 {
                                     DWDYDVDZDZ = (RHOKP * VISKP * (V_L[k + 1] - V_L[k]) * DDZP) * AREA_L[k + 1];
                                 }
                             }
                         }

                         //Assemble TDMA coefficients
                         //Half-cell below
                         if ((m - 1) == 1)
                         {
                             m--;
                             //note that the diffusion terms at the cell-faces are computed as harmoNI_Pc mean rather than arithmetic mean according to PataNK_Par 1980, chap. 4.2.3
                             DE = 0;
                             DE = 2 * VISH * (AREA_T1 * DDX_Rez + AREA_T2 * DDY_Rez) * RHO;

                             DB = 0;
                             if (k > 1)
                             {
                                 DB = 4 * VISH * VISH_KM / (VISH + VISH_KM) *
                                (AREAZX_L[k] * DDX_Rez + AREAZY_L[k] * DDY_Rez) * RHO_L_KM1;
                             }

                             DW1 = 0;
                             if (i > 1)
                             {
                                 DW1 = 4 * VISH * VISH_IM * VISH_T1 *
                                AREAX_L[k] * DDX_Rez * 0.5 * (RHO + Program.RHOImm[i - 1][j][k]);
                             }

                             DS = 0;
                             if (j > 1)
                             {
                                 DS = 4 * VISH * VISH_JM * VISH_T2 *
                                AREAY_L[k] * DDY_Rez * 0.5 * (RHO + Program.RHOImm[i][j - 1][k]);
                             }

                             //Advection terms
                             FE = (0.5 * RHO * (U1_L[k] + U2_L[k]) * AREA_T1 +
                                   0.5 * RHO * (V1_L[k] + V2_L[k]) * AREA_T2 +
                                   0.5 * RHO * (W1_L[k] + W2_L[k]) * AEREA_LL);

                             FB = 0;
                             if (k > 1)
                             {
                                 FB = (0.5 * (W1_L[k] * RHO + W2_L[k - 1] * RHO_L_KM1) * AEREA_LL +
                                      0.5 * (U1_L[k] * RHO + U2_L[k - 1] * RHO_L_KM1) * AREAZX_L[k] +
                                      0.5 * (V1_L[k] * RHO + V2_L[k - 1] * RHO_L_KM1) * AREAZY_L[k]);
                             }

                             if (i > 1)
                             {
                                 FW1 = 0.5 * (U1_L[k] * RHO + Program.U2[i - 1][j][k] * Program.RHOImm[i - 1][j][k]) * AREAX_L[k];
                             }
                             else
                             {
                                 FW1 = U1_L[k] * RHO * AREAX_L[k];
                             }

                             if (j > 1)
                             {
                                 FS = 0.5 * (V1_L[k] * RHO + V2_L_M * Program.RHOImm[i][j - 1][k]) * AREAY_L[k];
                             }
                             else
                             {
                                 FS = V1_L[k] * RHO * AREAY_L[k];
                             }

                             //Peclet numbers
                             DE = Math.Max(DE, 0.0001);
                             DB = Math.Max(DB, 0.0001);
                             DW1 = Math.Max(DW1, 0.0001);
                             DS = Math.Max(DS, 0.0001);

                             if (k < NK_P)
                             {
                                 PE = Math.Abs(FE / DE);
                             }
                             else
                             {
                                 PE = 0;
                             }

                             if (k > 1)
                             {
                                 PB = Math.Abs(FB / DB);
                             }
                             else
                             {
                                 PB = 0;
                             }

                             if (i > 1)
                             {
                                 PW1 = Math.Abs(FW1 / DW1);
                             }
                             else
                             {
                                 PW1 = 0;
                             }

                             if (j > 1)
                             {
                                 PS = Math.Abs(FS / DS);
                             }
                             else
                             {
                                 PS = 0;
                             }

                             //calculate coefficients of source terms
                             /*double SMP = (FE - FW1 - FB - FS) * 0;
                               double CPI = Math.Max(0, SMP);
                               double SP = -CPI;
                               Program.SU[i][j][kn] = CPI;
                              */

                             //advection scheme "power-law" by PataNK_Par 1980, p90
                             BIM_L[kn] = (float)(DE * Math.Max(0, Pow5(1 - 0.1 * PE)) + Math.Max(-FE, 0));
                             CIM_L[kn] = (float)(DB * Math.Max(0, Pow5(1 - 0.1 * PB)) + Math.Max(FB, 0));
                             AW1_L[k] = (float)(DW1 * Math.Max(0, Pow5(1 - 0.1 * PW1)) + Math.Max(FW1, 0));
                             AS1_L[k] = (float)(DS * Math.Max(0, Pow5(1 - 0.1 * PS)) + Math.Max(FS, 0));

                             AP0_L[k] = (float)(VOL_L[k] * 0.5 / Program.DT * RHO);
                             AIM_L[kn] = BIM_L[kn] + CIM_L[kn] + AW1_L[k] + AS1_L[k] + AP0_L[k]; // -SP;

                             //additional diffusion terms due to the usage of an eddy-diffusivity model
                             //note that in the first grid cell, the effect of diffusivity is computed using the friction velocity
                             DIMU_L[kn] = 0;
                             DIMV_L[kn] = 0;
                             DIMW_L[kn] = 0;

                             if ((i > 1) && (j > 1) && (i < NI_P) && (j < NJ_P) && (k < NK_P) && (k > 1))
                             {
                                 double DUDXE = 2 * (U2_L[k] - U1_L[k]) * DDX_Rez;
                                 double DUDXW = 2 * (U1_L[k] - Program.U2[i - 1][j][k]) * DDX_Rez;
                                 double DUDXB = 2 * 0.5 * ((U2_L[k] - U1_L[k]) + (U2_L[k - 1] - U1_L[k - 1])) * DDX_Rez;
                                 double DUDYB = 2 * 0.5 * ((U2_L[k] - U1_L[k]) + (U2_L[k - 1] - U1_L[k - 1])) * DDY_Rez;

                                 double DVDXS = 2 * 0.5 * ((V2_L[k] - V1_L[k]) + (V2_L_M - Program.V1[i][j - 1][k])) * DDX_Rez;
                                 double DVDXE = 2 * (V2_L[k] - V1_L[k]) * DDX_Rez;
                                 double DVDXB = 2 * 0.5 * ((V2_L[k] - V1_L[k]) + (V2_L[k - 1] - V1_L[k - 1])) * DDX_Rez;

                                 double DWDXE = 2 * (W2_L[k] - W1_L[k]) * DDX_Rez;
                                 double DWDXB = 2 * 0.5 * ((W2_L[k] - W1_L[k]) + (W2_L[k - 1] - W1_L[k - 1])) * DDX_Rez;
                                 double VISB = 2 * VISH_L[k] * VISH_L[k - 1] / (VISH_L[k] + VISH_L[k - 1]);

                                 DIMU_L[kn] = (VISH * DUDXE * AREA_T1
                                              - 2 * VISH * VISH_IM * VISH_T1 * DUDXW
                                              - VISB * DUDXB * AREAZX_L[k]
                                              - 2 * VISH * VISH_JM * VISH_T2 * DVDXS * AREAY_L[k]
                                              + VISH * (DVDXE * AREA_T2 + DWDXE * AEREA_LL)
                                              - VISB * DWDXB * AEREA_LL
                                              - VISB * DVDXB * AREAZY_L[k]) * RHO
                                    + DWDXDUDZDZ;

                                 double DVDYS = 2 * 0.5 * ((V2_L[k] - V1_L[k]) + (V2_L_M - Program.V1[i][j - 1][k])) * DDY_Rez;
                                 double DVDYE = 2 * (V2_L[k] - V1_L[k]) * DDY_Rez;
                                 double DVDYB = 2 * 0.5 * ((V2_L[k] - V1_L[k]) + (V2_L[k - 1] - V1_L[k - 1])) * DDY_Rez;

                                 double DUDYW = 2 * 0.5 * ((U2_L[k] - U1_L[k]) + (Program.U2[i - 1][j][k] - Program.U1[i - 1][j][k])) * DDY_Rez;
                                 double DUDYE = 2 * (U2_L[k] - U1_L[k]) * DDY_Rez;

                                 double DWDYE = 2 * (W2_L[k] - W1_L[k]) * DDY_Rez;
                                 double DWDYB = 2 * 0.5 * ((W2_L[k] - W1_L[k]) + (W2_L[k - 1] - W1_L[k - 1])) * DDY_Rez;

                                 DIMV_L[kn] = (VISH * DVDYE * AREA_T2
                                              - 2 * VISH * VISH_JM * VISH_T2 * DVDYS * AREAY_L[k]
                                              - VISB * DVDYB * AREAZY_L[k]
                                              - 2 * VISH * VISH_IM * VISH_T1 * DUDYW * AREAX_L[k]
                                              + VISH * (DUDYE * AREA_T1 + DWDYE * AEREA_LL)
                                              - VISB * DWDYB * AEREA_LL
                                              - VISB * DUDYB * AREAZX_L[k]) * RHO
                                    + DWDYDVDZDZ;

                                 double f1 = 1 / (ZSP_L_P - ZSP_LL) * stretch;
                                 double f2 = 1 / (ZSP_LL - ZSP_L_M) * stretchRzp;

                                 double DWDZE = (W2_L[k] - W1_L[k]) * f1;
                                 double DWDZB = (W1_L[k] - W2_L[k - 1]) * f2;

                                 double DUDZE = (U2_L[k] - U1_L[k]) * f1;
                                 double DUDZW = 0.5 * ((U1_L[k] - Program.U2[i - 1][j][k]) + (U2_L[k] - U1_L[k])) * f1;
                                 double DUDZB = (U1_L[k] - U2_L[k - 1]) * f2;

                                 double DVDZE = (V2_L[k] - V1_L[k]) * f1;
                                 double DVDZS = 0.5 * ((V1_L[k] - V2_L_M) + (V2_L[k] - V1_L[k])) * f1;
                                 double DVDZB = (V1_L[k] - V2_L[k - 1]) * f2;

                                 DIMW_L[kn] = (VISH * DWDZE * AEREA_LL
                                              - VISB * DWDZB * AEREA_LL
                                              - 2 * VISH * VISH_IM * VISH_T1 * DUDZW * AREAX_L[k]
                                              + VISH * (DVDZE * AREA_T2
                                                        + DUDZE * AREA_T1)
                                              - VISB * (DVDZB * AREAZY_L[k] + DUDZB * AREAZX_L[k])
                                              - 2 * VISH_JM * VISH * VISH_T2 * DVDZS * AREAY_L[k]) * RHO;
                             }

                             if (k == 1)
                             {
                                 DIMU_L[kn] = DWDXDUDZDZ;
                                 DIMV_L[kn] = DWDYDVDZDZ;
                             }
                         }
                         //Half-cell above
                         else
                         {
                             m = 2;

                             DE2 = 0;
                             if (i < NI_P)
                             {
                                 DE2 = 4 * VISH * VISH_I * VISH_R2 *
                                AREA_T3 * DDX_Rez * 0.5 * (RHO + Program.RHOImm[i + 1][j][k]);
                             }

                             DS2 = 0;
                             DS2 = 2 * VISH * (AREA_T1 * DDX_Rez + AREA_T2 * DDY_Rez) * RHO;

                             DN2 = 0;
                             if (j < NJ_P)
                             {
                                 DN2 = 4 * VISH * VISH_J * VISH_R1 *
                                AREA_T4 * DDY_Rez * 0.5 * (RHO + Program.RHOImm[i][j + 1][k]);
                             }

                             DT2 = 0;
                             if (k < NK_P)
                             {
                                 DT2 = 4 * VISH * VISH_L[k + 1] / (VISH + VISH_L[k + 1]) *
                                (AEREAZX_LL * DDX_Rez + AEREAZY_LL * DDY_Rez) * RHO_L[k + 1];
                             }

                             //Advection terms
                             if (i < NI_P)
                             {
                                 FE2 = 0.5 * (U2_L[k] * RHO + Program.U1[i + 1][j][k] * Program.RHOImm[i + 1][j][k]) * AREA_T3;
                             }
                             else
                             {
                                 FE2 = U2_L[k] * RHO * AREA_T3;
                             }

                             FS2 = (0.5 * RHO * (U1_L[k] + U2_L[k]) * AREA_T1 +
                                   0.5 * RHO * (V1_L[k] + V2_L[k]) * AREA_T2 +
                                   0.5 * RHO * (W1_L[k] + W2_L[k]) * AEREA_LL);

                             if (j < NJ_P)
                             {
                                 FN2 = 0.5 * (V2_L[k] * RHO + Program.V1[i][j + 1][k] * Program.RHOImm[i][j + 1][k]) * AREA_T4;
                             }
                             else
                             {
                                 FN2 = V2_L[k] * RHO * AREA_T4;
                             }

                             FT2 = 0;
                             if (k < NK_P)
                             {
                                 FT2 = (0.5 * (W2_L[k] * RHO + W1_L[k + 1] * RHO_KP) * AEREA_LL +
                                       0.5 * (U2_L[k] * RHO + U1_L[k + 1] * RHO_KP) * AEREAZX_LL +
                                       0.5 * (V2_L[k] * RHO + V1_L[k + 1] * RHO_KP) * AEREAZY_LL);
                             }
                             else
                             {
                                 FT2 = W2_L[k] * RHO * AEREA_LL +
                                    U2_L[k] * RHO * AEREAZX_LL +
                                    V2_L[k] * RHO * AEREAZY_LL;
                             }

                             //Peclet numbers
                             DE2 = Math.Max(DE2, 0.0001);
                             DT2 = Math.Max(DT2, 0.0001);
                             DN2 = Math.Max(DN2, 0.0001);
                             DS2 = Math.Max(DS2, 0.0001);

                             if (i < NI_P)
                             {
                                 PE2 = Math.Abs(FE2 / DE2);
                             }
                             else
                             {
                                 PE2 = 0;
                             }

                             if (k < NK_P)
                             {
                                 PT2 = Math.Abs(FT2 / DT2);
                             }
                             else
                             {
                                 PT2 = 0;
                             }

                             if (j < NJ_P)
                             {
                                 PN2 = Math.Abs(FN2 / DN2);
                             }
                             else
                             {
                                 PN2 = 0;
                             }

                             if (k < NK_P)
                             {
                                 PS2I = Math.Abs(FS2 / DS2);
                             }
                             else
                             {
                                 PS2I = 0;
                             }

                             //calculate coefficients of source terms
                             /*
                         double SMP = FE2 + FN2 + FT2 - FS2;
                         double CPI = Math.Max(0, SMP);
                         double SP = -CPI;
                         Program.SU[i][j][kn] = CPI;
                              */

                             //advection scheme "power-law" by PataNK_Par 1980, p90
                             BIM_L[kn] = (float)(DT2 * Math.Max(0, Pow5(1 - 0.1 * PT2)) + Math.Max(-FT2, 0));
                             CIM_L[kn] = (float)(DS2 * Math.Max(0, Pow5(1 - 0.1 * PS2I)) + Math.Max(FS2, 0));
                             AE2_L[k] = (float)(DE2 * Math.Max(0, Pow5(1 - 0.1 * PE2)) + Math.Max(-FE2, 0));
                             AN2_L[k] = (float)(DN2 * Math.Max(0, Pow5(1 - 0.1 * PN2)) + Math.Max(-FN2, 0));

                             AIM_L[kn] = BIM_L[kn] + CIM_L[kn] + AE2_L[k] + AN2_L[k] + AP0_L[k]; // -SP;

                             //additional diffusion terms due to the usage of an eddy-diffusivity model
                             //note that in the first grid cell, the effect of diffusivity is computed using the friction velocity
                             DIMU_L[kn] = 0;
                             DIMV_L[kn] = 0;
                             DIMW_L[kn] = 0;

                             if ((i > 1) && (j > 1) && (i < NI_P) && (j < NJ_P) && (k < NK_P) && (k > 1))
                             {
                                 double DUDXE = 2 * (Program.U1[i + 1][j][k] - U2_L[k]) * DDX_Rez;
                                 double DUDXS = 2 * (U2_L[k] - U1_L[k]) * DDX_Rez;
                                 double DUDXT = 2 * 0.5 * ((U2_L[k] - U1_L[k]) + (U2_L[k + 1] - U1_L[k + 1])) * DDX_Rez;
                                 double DUDYT = 2 * 0.5 * ((U2_L[k] - U1_L[k]) + (U2_L[k + 1] - U1_L[k + 1])) * DDY_Rez;

                                 double DVDXN = 2 * 0.5 * ((V2_L[k] - V1_L[k]) + (Program.V2[i][j + 1][k] - Program.V1[i][j + 1][k])) * DDX_Rez;
                                 double DVDXS = 2 * (V2_L[k] - V1_L[k]) * DDX_Rez;
                                 double DVDXT = 2 * 0.5 * ((V2_L[k] - V1_L[k]) + (V2_L[k + 1] - V1_L[k + 1])) * DDX_Rez;

                                 double DWDXS = 2 * (W2_L[k] - W1_L[k]) * DDX_Rez;
                                 double DWDXT = 2 * 0.5 * ((W2_L[k] - W1_L[k]) + (W2_L[k + 1] - W1_L[k + 1])) * DDX_Rez;

                                 double VIST = 2 * VISH * VISH_L[k + 1] / (VISH + VISH_L[k + 1]);

                                 DIMU_L[kn] = (2 * VISH * VISH_I * VISH_R2 * DUDXE * AREA_T3
                                              - VISH * DUDXS * AREA_T1
                                              + VIST * DUDXT * (AREAX_L[k + 1] + AEREAZX_LL)
                                              + 2 * VISH * VISH_J * VISH_R1 * DVDXN * AREA_T4
                                              - VISH * (DVDXS * AREA_T2 + DWDXS * AEREA_LL)
                                              + VIST * (DVDXT * AEREAZY_LL + DWDXT * AEREA_LL)) * RHO
                                    + DWDXDUDZDZ;

                                 double DVDYN = 2 * 0.5 * ((V2_L[k] - V1_L[k]) + (Program.V2[i][j + 1][k] - Program.V1[i][j + 1][k])) * DDY_Rez;
                                 double DVDYW = 2 * (V2_L[k] - V1_L[k]) * DDY_Rez;
                                 double DVDYT = 2 * 0.5 * ((V2_L[k] - V1_L[k]) + (V2_L[k + 1] - V1_L[k + 1])) * DDY_Rez;

                                 double DUDYE = 2 * 0.5 * ((U2_L[k] - U1_L[k]) + (Program.U2[i + 1][j][k] - Program.U1[i + 1][j][k])) * DDY_Rez;
                                 double DUDYW = 2 * (U2_L[k] - U1_L[k]) * DDY_Rez;

                                 double DWDYW = 2 * (W2_L[k] - W1_L[k]) * DDY_Rez;
                                 double DWDYT = 2 * 0.5 * ((W2_L[k] - W1_L[k]) + (W2_L[k + 1] - W1_L[k + 1])) * DDY_Rez;

                                 DIMV_L[kn] = (2 * VISH * VISH_J * VISH_R1 * DVDYN * AREA_T4
                                              + VIST * DVDYT * AEREAZY_LL
                                              - VISH * DVDYW * AREA_T2
                                              - VISH * (DWDYW * AEREA_LL + DUDYW * AREA_T1)
                                              + VIST * (DWDYT * AEREA_LL + DUDYT * AEREAZX_LL)
                                              + 2 * VISH * VISH_I * VISH_R2 * DUDYE * AREA_T3) * RHO
                                    + DWDYDVDZDZ;

                                 double f1 = 1 / (ZSP_L_P - ZSP_LL) * stretch;

                                 double DWDZW = (W2_L[k] - W1_L[k]) * f1;
                                 double DWDZT = (W1_L[k] - W2_L[k + 1]) / (ZSP_LL - ZSP_L_P) * stretchRzp;

                                 double DUDZW = (U2_L[k] - U1_L[k]) * f1;
                                 double DUDZE = 0.5 * ((Program.U1[i + 1][j][k] - U2_L[k]) + (U2_L[k] - U1_L[k])) * f1;
                                 double DUDZT = (U1_L[k + 1] - U2_L[k]) * f1;

                                 double DVDZW = (V2_L[k] - V1_L[k]) * f1;
                                 double DVDZN = 0.5 * ((Program.V1[i][j + 1][k] - V2_L[k]) + (V2_L[k] - V1_L[k])) * f1;
                                 double DVDZT = (V1_L[k + 1] - V2_L[k]) * f1;

                                 DIMW_L[kn] = (VIST * DWDZT * AEREA_LL
                                              - VISH * DWDZW * AEREA_LL
                                              + 2 * VISH * VISH_J * VISH_R1 * DVDZN * AREA_T4
                                              - VISH * (DVDZW * AREA_T2
                                                        - DUDZW * AREA_T1)
                                              + VIST * (DVDZT * AEREAZY_LL + DUDZT * AEREAZX_LL)
                                              + 2 * VISH * VISH_I * VISH_R2 * DUDZE * AREA_T3) * RHO;

                             }

                             if (k == 1)
                             {
                                 DIMU_L[kn] = DWDXDUDZDZ;
                                 DIMV_L[kn] = DWDYDVDZDZ;
                             }
                         }
                     }
                 }
             });
        }
    }
}
