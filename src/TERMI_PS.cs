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
using System.Runtime.CompilerServices;
using System.Threading.Tasks;

namespace GRAMM_2001
{
    partial class Program
    {
        ///<summary>
        ///Calculate the diffusion and advection terms for the implicit scheme (Patankar 1980, p52)
        ///</summary>
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        public static void TERMIPSterms(int NI, int NJ, int NK)
        {
            Parallel.For(2, NI, Program.pOptions, i =>
            {
                double DEAST, DWEST, DSOUTH, DNORTH, DBOTTOM, DTOP, FEAST, FWEST, FSOUTH, FNORTH, FBOTTOM, FTOP;
                double PEAST, PWEST, PSOUTH, PNORTH, PBOTTOM, PTOP;
                double VISH1, AREA;
                double RHO_I1, RHO_IM1, RHO_J1, RHO_JM1;
                double DDX_Rez = 1 / Program.DDX[i];
                float RHO;

                for (int j = 2; j <= NJ - 1; j++)
                {
                    ReadOnlySpan <float> AREA_L = Program.AREA[i][j];
                    ReadOnlySpan <float> AREAX_L = Program.AREAXImm[i][j];
                    ReadOnlySpan <float> AREAY_L = Program.AREAY[i][j];
                    ReadOnlySpan <float> AREAZX_L = Program.AREAZX[i][j];
                    ReadOnlySpan <float> AREAZY_L = Program.AREAZY[i][j];
                    ReadOnlySpan <float> RHO_L = Program.RHO[i][j];
                    ReadOnlySpan <float> ZSP_L = Program.ZSP[i][j];
                    ReadOnlySpan <float> VOL_L = Program.VOL[i][j];
                    ReadOnlySpan <float> AREAXiP = Program.AREAXImm[i + 1][j];
                    ReadOnlySpan <float> AREAXjP = Program.AREAY[i][j + 1];
                    ReadOnlySpan <float> RhoiP = Program.RHO[i + 1][j];
                    ReadOnlySpan <float> RhoiM = Program.RHO[i - 1][j];
                    ReadOnlySpan <float> RhojP = Program.RHO[i][j + 1];
                    ReadOnlySpan <float> RhojM = Program.RHO[i][j - 1];

                    ReadOnlySpan <double> UN_L = Program.UN[i][j];
                    ReadOnlySpan <double> VN_L = Program.VN[i][j];
                    ReadOnlySpan <double> WN_L = Program.WN[i][j];
                    ReadOnlySpan <double> VISH_L = Program.VISH[i][j];
                    ReadOnlySpan <double> VISHiP = Program.VISH[i + 1][j];
                    ReadOnlySpan <double> VISHiM = Program.VISH[i - 1][j];
                    ReadOnlySpan <double> VISHjP = Program.VISH[i][j + 1];
                    ReadOnlySpan <double> VISHjM = Program.VISH[i][j - 1];
                    ReadOnlySpan <double> UNiP = Program.UN[i + 1][j];
                    ReadOnlySpan <double> UNiM = Program.UN[i - 1][j];
                    ReadOnlySpan <double> VNjP = Program.VN[i][j + 1];
                    ReadOnlySpan <double> VNjM = Program.VN[i][j - 1];

                    float[] AEAST_PS_L = Program.AEAST_PS[i][j];
                    float[] ANORTH_PS_L = Program.ANORTH_PS[i][j];
                    float[] ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                    float[] AWEST_PS_L = Program.AWEST_PS[i][j];
                    float[] AP0_PS_L = Program.AP0_PS[i][j];
                    double[] A_PS_L = Program.A_PS[i][j];
                    double[] B_PS_L = Program.B_PS[i][j];
                    double[] C_PS_L = Program.C_PS[i][j];
                    
                    double DDY_Rez = 1 / Program.DDY[j];

                    for (int k = 1; k <= NK; k++)
                    {
                        DEAST = DWEST = DSOUTH = DNORTH = DBOTTOM = DTOP = 0;
                        FEAST = FWEST = FSOUTH = FNORTH = FBOTTOM = FTOP = 0;
                        PEAST = PWEST = PSOUTH = PNORTH = PBOTTOM = PTOP = 0;

                        RHO = RHO_L[k];
                        VISH1 = VISH_L[k];
                        AREA = AREA_L[k];

                        RHO_I1 = RhoiP[k]; RHO_IM1 = RhoiM[k];
                        RHO_J1 = RhojP[k]; RHO_JM1 = RhojM[k];

                        //DIFFUSION TERMS (NOTE THAT VIS IS COMPUTED AS HARMONIC MEAN RATHER THAN ARITHMETIC MEAN AT BORDERS ACCORDING TO PATANKAR 1980, CHAP. 4.2.3)
                        DEAST = 0.5 * (VISH1 + VISHiP[k]) * AREAXiP[k] * DDX_Rez * 0.5 * (RHO + RHO_I1);
                        DWEST = 0.5 * (VISH1 + VISHiM[k]) * AREAX_L[k] * DDX_Rez * 0.5 * (RHO + RHO_IM1);
                        DNORTH = 0.5 * (VISH1 + VISHjP[k]) * AREAXjP[k] * DDY_Rez * 0.5 * (RHO + RHO_J1);
                        DSOUTH = 0.5 * (VISH1 + VISHjM[k]) * AREAY_L[k] * DDY_Rez * 0.5 * (RHO + RHO_JM1);
                        if (k > 1)
                        {
                            DBOTTOM = (0.5 * (VISH1 + VISH_L[k - 1]) * (AREA / (ZSP_L[k] - ZSP_L[k - 1]) +
                              AREAZX_L[k] * DDX_Rez + AREAZY_L[k] * DDY_Rez)) * 0.5 * (RHO + RHO_L[k - 1]);
                        }

                        if (k < NK)
                        {
                            DTOP = (0.5 * (VISH1 + VISH_L[k + 1]) * (AREA / (ZSP_L[k + 1] - ZSP_L[k]) +
                              AREAZX_L[k + 1] * DDX_Rez + AREAZY_L[k + 1] * DDY_Rez)) * 0.5 * (RHO + RHO_L[k + 1]);
                        }

                        //Advection terms
                        FEAST = 0.5 * (UNiP[k] + UN_L[k]) * AREAXiP[k] * 0.5 * (RHO + RHO_I1);
                        FWEST = 0.5 * (UNiM[k] + UN_L[k]) * AREAX_L[k] * 0.5 * (RHO + RHO_IM1);
                        FNORTH = 0.5 * (VNjP[k] + VN_L[k]) * AREAXjP[k] * 0.5 * (RHO + RHO_J1);
                        FSOUTH = 0.5 * (VNjM[k] + VN_L[k]) * AREAY_L[k] * 0.5 * (RHO + RHO_JM1);
                        if (k > 1)
                        {
                            FBOTTOM = (0.5 * (WN_L[k - 1] + WN_L[k]) * AREA
                                  + 0.5 * (UN_L[k - 1] + UN_L[k]) * AREAZX_L[k]
                                  + 0.5 * (VN_L[k - 1] + VN_L[k]) * AREAZY_L[k]) *
                                  0.5 * (RHO + RHO_L[k - 1]);
                        }

                        if (k < NK)
                        {
                            FTOP = (0.5 * (WN_L[k + 1] + WN_L[k]) * AREA
                                  + 0.5 * (UN_L[k + 1] + UN_L[k]) * AREAZX_L[k + 1]
                                  + 0.5 * (VN_L[k + 1] + VN_L[k]) * AREAZY_L[k + 1]) *
                                  0.5 * (RHO + RHO_L[k + 1]);
                        }

                        //Peclet numbers
                        DEAST = Math.Max(DEAST, 0.0001);
                        DWEST = Math.Max(DWEST, 0.0001);
                        DSOUTH = Math.Max(DSOUTH, 0.0001);
                        DNORTH = Math.Max(DNORTH, 0.0001);
                        DBOTTOM = Math.Max(DBOTTOM, 0.0001);
                        DTOP = Math.Max(DTOP, 0.0001);
                        PEAST = Math.Abs(FEAST / DEAST);
                        PWEST = Math.Abs(FWEST / DWEST);
                        PSOUTH = Math.Abs(FSOUTH / DSOUTH);
                        PNORTH = Math.Abs(FNORTH / DNORTH);
                        PBOTTOM = Math.Abs(FBOTTOM / DBOTTOM);
                        PTOP = Math.Abs(FTOP / DTOP);

                        //calculate coefficients of source terms
                        /*double SMP = (FEAST - FWEST + FNORTH - FSOUTH - FBOTTOM + FTOP) * 0;
                        double CPI = Math.Max(0, SMP);
                        double SP = -CPI;
                        Program.SU_PS[i][j][k] = CPI;
                         */

                        //advection scheme "power-law" by Patankar 1980, p90
                        B_PS_L[k] = DTOP * Math.Max(0, Pow5(1 - 0.1 * PTOP)) + Math.Max(-FTOP, 0);
                        C_PS_L[k] = DBOTTOM * Math.Max(0, Pow5(1 - 0.1 * PBOTTOM)) + Math.Max(FBOTTOM, 0);
                        AWEST_PS_L[k] = (float)(DWEST * Math.Max(0, Pow5(1 - 0.1 * PWEST)) + Math.Max(FWEST, 0));
                        ASOUTH_PS_L[k] = (float)(DSOUTH * Math.Max(0, Pow5(1 - 0.1 * PSOUTH)) + Math.Max(FSOUTH, 0));
                        AEAST_PS_L[k] = (float)(DEAST * Math.Max(0, Pow5(1 - 0.1 * PEAST)) + Math.Max(-FEAST, 0));
                        ANORTH_PS_L[k] = (float)(DNORTH * Math.Max(0, Pow5(1 - 0.1 * PNORTH)) + Math.Max(-FNORTH, 0));

                        AP0_PS_L[k] = (float)(VOL_L[k] / Program.DT * RHO);
                        A_PS_L[k] = B_PS_L[k] + C_PS_L[k] + AWEST_PS_L[k] + ASOUTH_PS_L[k]
                          + ANORTH_PS_L[k] + AEAST_PS_L[k] + AP0_PS_L[k];  // -SP;
                    }
                }
            });
        }
    }
}
