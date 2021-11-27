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

namespace GRAMM_2001
{
    partial class Program
    {
        /// <summary>
        /// Calculate the new soil temperature
        /// </summary>
        /// <param name="NI"></param>
        /// <param name="NJ"></param>
        /// <param name="NK"></param>
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        public static void Calctb_calculate(int NI, int NJ, int NK)
        {
            double EINTRAG1 = 0;
            double EINTRAG2 = 0;
            double EINTRAG3 = 0;
            double EINTRAG4 = 0;
            double EINTRAG5 = 0;

            /*
            double WQUmax = 0;
            double USTmax = 0;
            int imax = 1;
            int jmax = 1;
            */

            //Surface energy balance is not computed at the edges, where topography is smoothed - OETTL 23.Dec 13
            for (int i = 2; i <= NI - 1; i++)
            {
                double TB2_L;
                for (int j = 2; j <= NJ - 1; j++)
                {
                    TB2_L = Program.TB[i][j][2];
                    Program.TB[i][j][1] = TB2_L;
                    EINTRAG1 += Program.WQU[i][j] / Program.DDXImm[i] / Program.DDYImm[j];
                    EINTRAG2 += Program.RSOLG[i][j];
                    EINTRAG3 += (Program.RL[i][j] - Program.EPSG[i][j] * Program.SIGMA * Pow4(TB2_L));

                    /*
                    if (Math.Abs(Program.WQU[i][j] / Program.DDXImm[i] / Program.DDYImm[j]) > WQUmax)
                    {
                        WQUmax = Math.Abs(Program.WQU[i][j] / Program.DDXImm[i] / Program.DDYImm[j]);
                        imax = i;
                        jmax = j;
                        USTmax = Program.UST[i][j] * Program.USTV[i][j];
                    }
                    */
                }
            }
            //Console.WriteLine("max. sens. Wärmefluss: " + (Program.WQU[imax][jmax] / Program.DDXImm[imax] / Program.DDYImm[jmax]).ToString("0.0") + " i: " + imax.ToString() + " j: " + jmax.ToString() + " USTxUSTV: " + USTmax.ToString("0.000") + " XWQ :" + Program.XWQ[imax][jmax].ToString("0.000") + " TB: " + (Program.TB[imax][jmax][2]).ToString("0.0") + " TLUFT: " + (Program.TABS[imax][jmax][1]).ToString("0.0") + " OL: " + Program.OL[imax][jmax].ToString("0.0"));

            int NLAND = 0;
            for (int i = 2; i <= NI - 1; i++)
            {
                for (int j = 2; j <= NJ - 1; j++)
                {
                    NLAND++;
                    for (int k = 2; k <= Program.NZB - 1; k++)
                    {
                        double DTBDZ = 0;
                        double TERMTB = 0;
                        if (k > 2)
                        {
                            DTBDZ = (Program.ALAMBDA[i][j] * (Program.TB[i][j][k + 1] - Program.TB[i][j][k]) / Program.DWB[k + 1] -
                                   Program.ALAMBDA[i][j] * (Program.TB[i][j][k] - Program.TB[i][j][k - 1]) / Program.DWB[k]) / Program.DZB[k];
                            TERMTB = DTBDZ * Program.DT;
                        }
                        else if (k == 2)
                        {
                            DTBDZ = (Program.ALAMBDA[i][j] * (Program.TB[i][j][k + 1] - Program.TB[i][j][k]) / Program.DWB[k + 1] -
                                     Program.ALAMBDA[i][j] * (Program.TB[i][j][k] - Program.TB[i][j][k - 1]) / Program.DWB[k]) / Program.DZB[k];


                            //influence of snow cover, except in urban areas and water bodies
                            if (Program.SNOW[i][j] > 0.15 && Program.SNOW[i][j] > Program.Z0[i][j] && Program.ALAMBDA[i][j] != 4.0 && Program.FW[i][j] < 0.99)
                            {
                                DTBDZ = 0.0;
                            }

                            //saturation vapour pressure - Magnus formulae
                            double EGSAT = 611.2 * Math.Exp(17.269 * (Program.TB[i][j][2] - 273.16) / (Program.TB[i][j][2] - 30.04));
                            //maximum water vapour content possible -> assuming condensation above 80% relative humidity
                            double QGSAT = EGSAT / 461.5 / Program.TB[i][j][2];
                            /*
                            double FW1 = 1.0;
                            if (QUGSAT > (Program.QU[i][j][1] * 0.001))
                            {
                                FW1 = Program.FW[i][j];
                            }
                            else
                            {
                                FW1 = 1.0;
                            }
                            */
                            Program.QUG[i][j] = Program.FW[i][j] * QGSAT * 1000.0 + (1 - Program.FW[i][j]) * Program.QU[i][j][1];
                            double ALW1 = Program.ALW - 2300 * (Program.TB[i][j][2] - 273.15);
                            double WQU_L = Program.XWQ[i][j] * Program.RHOImm[i][j][1] * Program.UST[i][j] * ALW1 * (Program.QU[i][j][1] - Program.QUG[i][j]) * 0.001;

                            //amount of water evaporated from the soil
                            Program.QU[i][j][1] += -WQU_L / ALW1 * 1000.0 * Program.DT * Program.AREAImm[i][j][1] / Program.VOLImm[i][j][1] / Program.RHOImm[i][j][1];

                            WQU_L = Math.Min(WQU_L, 20);
                            WQU_L = Math.Max(WQU_L, -400);
                            EINTRAG4 += WQU_L;

                            TERMTB = (DTBDZ + (Program.WQU[i][j] / Program.DDXImm[i] / Program.DDYImm[j] + WQU_L + Program.RSOLG[i][j] +
                                    (Program.RL[i][j] - Program.EPSG[i][j] * Program.SIGMA * Pow4(Program.TB[i][j][2]))
                                    + Program.AWQ[i][j]) / Program.DZB[2]) * Program.DT;
                            EINTRAG5 += DTBDZ * Program.DZB[2];

                            /*
                            if (i == Program.inrec[27] && j == Program.jnrec[27])
                            {
                                Console.WriteLine(Program.TB[i][j][k].ToString("0.00") + "  " + (TERMTB / Program.CPBOD / Program.RHOB[i][j]).ToString("0.000"));
                                Console.WriteLine("Boden: "+ (DTBDZ* Program.DZB[k]).ToString("0.00") + " Sensibler: " + (Program.WQU[i][j] / Program.DDXImm[i] / Program.DDYImm[j]).ToString("0.000"));
                                Console.WriteLine(" Latent: "+WQU_L.ToString("0.00") + "  Sonne: " + Program.RSOLG[i][j].ToString("0.000"));
                                Console.WriteLine(" langwellige Strahlungsbilanz: "+((Program.RL[i][j] - Program.EPSG[i][j] * Program.SIGMA * Pow4(Program.TB[i][j][2]))).ToString("0.00"));
                            }
                            */

                        }

                        Program.TBN[i][j][k] = Program.TB[i][j][k] + 1.0 * (TERMTB / Program.CPBOD / Program.RHOB[i][j]); //Relaxation der Bodentemperaturänderung

                        //influence of snow cover, except in urban areas and water bodies
                        if (Program.SNOW[i][j] > 0.15 && Program.SNOW[i][j] > Program.Z0[i][j] && Program.ALAMBDA[i][j] != 4.0 && Program.FW[i][j] < 0.99)
                        {
                            Program.TBN[i][j][k] = Math.Min(Program.TBN[i][j][k], 273);
                        }
                    }
                }
            }

            EINTRAG1 /= (NI - 1) * (NJ - 1);
            EINTRAG2 /= (NI - 1) * (NJ - 1);
            EINTRAG3 /= (NI - 1) * (NJ - 1);
            EINTRAG4 /= (NI - 1) * (NJ - 1);
            EINTRAG5 /= (NI - 1) * (NJ - 1);

            if (Program.TerminalOut >= TerminalThreshold) // 
            {
                Console.WriteLine(" WQU-S[W/m**2] WQU-L[W/m**2] RSOLG[W/m**2] RTERR[W/m**2] HEATF[W/m**2] RELAXV RELAXT");
                Console.WriteLine(EINTRAG1.ToString("0.0").PadLeft(14) + EINTRAG4.ToString("0.0").PadLeft(14) + EINTRAG2.ToString("0.0").PadLeft(14)
                                  + EINTRAG3.ToString("0.0").PadLeft(14) + EINTRAG5.ToString("0.0").PadLeft(14) + Program.RELAXV.ToString("0.000".PadLeft(7)) + Program.RELAXT.ToString("0.000".PadLeft(7)));
                Console.WriteLine("-------------------------------------------------------------------------------------");
                Program.TerminalOut = 0;
            }

        }
    }
}
