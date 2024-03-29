﻿#region Copyright
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
using System.Collections.Concurrent;
using System.Runtime.CompilerServices;
using System.Threading.Tasks;

namespace GRAMM_2001
{
    partial class Program
    {
        /// <summary>
        /// Calculate the Richardson number and the stability class for each cell
        /// </summary>
        /// <param name="NI"></param>
        /// <param name="NJ"></param>
        /// <param name="NK"></param>
        public static void Prandtl_calculate(int NI, int NJ, int NK)
        {
            Parallel.ForEach(PartitionerI[Program.StripeCounter++ % PartitionerI.Count], range =>
            {
                int NJ_P = NJ;
                for (int i = range.Item1; i < range.Item2; i++)
                {
                    for (int j = 2; j <= NJ_P - 1; j++)
                    {
                        //if (Program.GLOBRAD[i][j] > 150) Program.OL[i][j] = -20; //ACHTUNG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        ref double UST = ref Program.UST[i][j];
                        ref float OL = ref Program.OL[i][j];
                        ref double GlobRad = ref Program.GLOBRAD[i][j];
                        ReadOnlySpan<float> ZSP_L = Program.ZSPImm[i][j].AsSpan();
                        ReadOnlySpan<double> U_L = Program.U[i][j];
                        ReadOnlySpan<double> V_L = Program.V[i][j];
                        ReadOnlySpan<double> QU_L = Program.QU[i][j];

                        int k = 1;
                        double ANUE = Program.VISEL / Program.RHO[i][j][k];
                        double Z1 = ZSP_L[k] - Program.AHImm[i][j];
                        double ZETA = Math.Max(Z1 - Program.Z0[i][j], 0.1) / OL;

                        double PHIM = 0;
                        double PHIH = 0;

                        //flux profile relationships for stable...
                        double PSIM = -4.7 * ZETA;
                        double PSIH = -4.7 * ZETA;
                        //and convective conditions
                        if (ZETA <= 0)
                        {
                            PHIM = Math.Pow(1 - 16 * ZETA, 0.25);
                            PHIH = Pow2(PHIM);
                            PSIM = Math.Log((1 + Math.Pow(PHIM, -2)) * 0.5F * Pow2(1 + Math.Pow(PHIM, -1)) * 0.5F) - 2 * Math.Atan(1 / PHIM) + Math.PI * 0.5;
                            PSIH = 2 * Math.Log(0.5F + 0.5F * Math.Pow(PHIH, -1));
                        }

                        //computation of friction velocity and temperature scale
                        double VP = Math.Sqrt(Pow2(U_L[k]) + Pow2(V_L[k]) + Pow2(Program.W[i][j][k]));

                        UST = Program.CK * VP / (Math.Log(Z1 / Program.Z0[i][j]) - PSIM);
                        UST = Math.Min(UST, 2.0);
                        UST = Math.Max(UST, 0.15);

                        Program.USTV[i][j] = UST / Math.Max(VP, 0.1);
                        Program.TST[i][j] = Program.CK / (Math.Log(Z1 / Program.Z0[i][j]) - PSIM * ZETA) * ((Program.T[i][j][k] + Program.TBZ1) * (1 + 0.00061 * QU_L[k])
                                                   - Program.TB[i][j][2] * Program.FAC[i][j][k] * (1 + 0.00061 * Program.QUG[i][j]));
                        Program.XWQ[i][j] = (float)(Program.CK / (Math.Log(Z1 / Program.Z0[i][j]) - PSIH));

                        //Obukhov length
                        OL = (float)((Program.T[i][j][k] + Program.TBZ1) * (1 + 0.00061 * QU_L[k]) * Pow2(Program.UST[i][j]) / (Program.CK * Program.GERD * Program.TST[i][j]));
                        if (OL <= 0)
                        {
                            OL = Math.Min(-5, OL);
                        }

                        if (OL > 0)
                        {
                            OL = Math.Max(5, OL);
                        }

                        if (OL > 1000)
                        {
                            OL = Math.Min(1000, OL);
                        }

                        if (OL < -1000)
                        {
                            OL = Math.Max(-1000, OL);
                        }

                        //stability class
                        int ischnitt = 1;
                        int ischnitt_plus1 = 1;
                        double Uoben = 0;
                        double Voben = 0;
                        double Uunten = 0;
                        double Vunten = 0;
                        double Umittel = 0;
                        double Vmittel = 0;
                        //compute wind speed at 10m above ground level
                        for (int k1 = 1; k1 <= NK; k1++)
                        {
                            if (ZSP_L[k1] - Program.AHImm[i][j] >= 10)
                            {
                                ischnitt = k1;
                                break;
                            }
                        }
                        //compute temperature gradient over a vertical distance of >=90m
                        for (int k1 = 1; k1 <= NK; k1++)
                        {
                            if (ZSP_L[k1] - Program.AHImm[i][j] >= 90)
                            {
                                ischnitt_plus1 = Math.Max(k1, 2);
                                break;
                            }
                        }

                        Uoben = U_L[ischnitt];
                        Voben = V_L[ischnitt];
                        if (ischnitt > 1)
                        {
                            Uunten = U_L[ischnitt - 1];
                            Vunten = V_L[ischnitt - 1];
                            Umittel = Uunten + (Uoben - Uunten) / (ZSP_L[ischnitt] - ZSP_L[ischnitt - 1]) *
                            (10 - ZSP_L[ischnitt - 1] + Program.AHImm[i][j]);
                            Vmittel = Vunten + (Voben - Vunten) / (ZSP_L[ischnitt] - ZSP_L[ischnitt - 1]) *
                            (10 - ZSP_L[ischnitt - 1] + Program.AHImm[i][j]);
                        }
                        else
                        {
                            Umittel = Uoben / (ZSP_L[ischnitt] - Program.AHImm[i][j]) * 10;
                            Vmittel = Voben / (ZSP_L[ischnitt] - Program.AHImm[i][j]) * 10;
                        }
                        double wg = Math.Sqrt(Umittel * Umittel + Vmittel * Vmittel);

                        //int ischnitt_plus1 = Math.Max(ischnitt, 2);
                        double dt = Program.TABS[i][j][ischnitt_plus1] - Program.TABS[i][j][1];

                        if (wg <= 2)
                        {
                            if ((GlobRad <= 20) && (dt >= 0))
                            {
                                Program.stabilityclass[i][j] = 7;
                            }

                            if ((GlobRad <= 20) && (dt < 0))
                            {
                                Program.stabilityclass[i][j] = 6;
                            }

                            if ((GlobRad > 20) && (GlobRad <= 175))
                            {
                                Program.stabilityclass[i][j] = 4;
                            }

                            if ((GlobRad > 175) && (GlobRad <= 675))
                            {
                                Program.stabilityclass[i][j] = 2;
                            }

                            if (GlobRad > 675)
                            {
                                Program.stabilityclass[i][j] = 1;
                            }
                        }
                        if ((wg > 2) && (wg <= 3))
                        {
                            if ((GlobRad <= 20) && (dt >= 0))
                            {
                                Program.stabilityclass[i][j] = 6;
                            }

                            if ((GlobRad <= 20) && (dt < 0))
                            {
                                Program.stabilityclass[i][j] = 5;
                            }

                            if ((GlobRad > 20) && (GlobRad <= 175))
                            {
                                Program.stabilityclass[i][j] = 4;
                            }

                            if ((GlobRad > 175) && (GlobRad <= 675))
                            {
                                Program.stabilityclass[i][j] = 3;
                            }

                            if ((GlobRad > 675) && (GlobRad <= 925))
                            {
                                Program.stabilityclass[i][j] = 2;
                            }

                            if (GlobRad > 925)
                            {
                                Program.stabilityclass[i][j] = 1;
                            }
                        }
                        if ((wg > 3) && (wg <= 5))
                        {
                            if (GlobRad <= 175)
                            {
                                Program.stabilityclass[i][j] = 4;
                            }

                            if ((GlobRad > 175) && (GlobRad <= 675))
                            {
                                Program.stabilityclass[i][j] = 3;
                            }

                            if (GlobRad > 675)
                            {
                                Program.stabilityclass[i][j] = 2;
                            }
                        }
                        else if (wg > 5)
                        {
                            if (GlobRad <= 925)
                            {
                                Program.stabilityclass[i][j] = 4;
                            }

                            if (GlobRad > 925)
                            {
                                Program.stabilityclass[i][j] = 3;
                            }
                        }
                    }
                }
            });

            //Richardson number
            Parallel.ForEach(PartitionerI[Program.StripeCounter++ % PartitionerI.Count], range =>
            {
                int NJ_P = NJ; int NK_P = NK;

                for (int i = range.Item1; i < range.Item2; i++)
                {
                    for (int j = 2; j <= NJ_P - 1; j++)
                    {
                        ReadOnlySpan<float> ZSP_L = Program.ZSPImm[i][j].AsSpan();
                        ReadOnlySpan<double> U_L = Program.U[i][j];
                        ReadOnlySpan<double> V_L = Program.V[i][j];
                        ReadOnlySpan<double> QU_L = Program.QU[i][j];
                        ReadOnlySpan<double> T_L = Program.T[i][j];

                        double DZZZ = 1 / (ZSP_L[1] - Program.AHImm[i][j]);
                        double DTDZ = ((Program.T[i][j][1] + Program.TBZ1) * (1 + 0.00061 * QU_L[1]) -
                                       Program.TB[i][j][2] * Program.FAC[i][j][1] * (1 + 0.00061 * Program.QUG[i][j])) * DZZZ;
                        double DUDZ = Pow2(U_L[1] * DZZZ);
                        double DVDZ = Pow2(V_L[1] * DZZZ);
                        double DGDZ = DUDZ + DVDZ;
                        double[] RITSCH_L = Program.RITSCH[i][j];

                        RITSCH_L[0] = Program.GERD / (T_L[1] + Program.TBZ1) / (1 + 0.00061 * QU_L[1]);
                        RITSCH_L[0] *= DTDZ;
                        RITSCH_L[0] /= (DGDZ + 0.0000000001);
                        if (RITSCH_L[0] > 0)
                        {
                            RITSCH_L[0] = Math.Min(0.2, RITSCH_L[0]);
                        }

                        for (int k = 1; k < RITSCH_L.Length; ++k)
                        {
                            DZZZ = 1 / (ZSP_L[k + 1] - ZSP_L[k]);
                            DTDZ = (T_L[k + 1] * (1 + 0.00061 * QU_L[k + 1]) -
                                   T_L[k] * (1 + 0.00061 * QU_L[k])) * DZZZ;
                            DUDZ = Pow2((U_L[k + 1] - U_L[k]) * DZZZ);
                            DVDZ = Pow2((V_L[k + 1] - V_L[k]) * DZZZ);
                            DGDZ = DUDZ + DVDZ;
                            RITSCH_L[k] = Program.GERD / (T_L[k] + Program.TBZ1);
                            RITSCH_L[k] *= DTDZ;
                            RITSCH_L[k] /= (DGDZ + 0.0000000001);
                            /*
                             if (i == Program.inrec[63] && j == Program.jnrec[63])
                             {
                                 Console.WriteLine(DTDZ.ToString("0.000") + " , " + k.ToString() + " , " + RITSCH_L[k].ToString("0.000") + " , " + DGDZ.ToString("0.000"));
                             }
                             */
                        }
                    }
                }
            });

            /*PRANDTLSCHICHTTHEORIE (SOMIESKI, DETERING, FLASSAK) 
* 
*       SOLARE UND TERRESTRISCHE STRAHLUNG
*
*       SOLARE STRAHLUNG                RSOLG(NI,NJ_P)              [ W/m**2 ]
*                                       EPSG(NI,NJ_P)               [       ]
*                                       RL(NI,NJ_P)                 [       ]
*       TEMPERATUR DER BODENZELLE       TB(NI,NJ_P,2)               [ K ]
*       WAERMELEITFAEHIGKEIT DES BODENS ALAMBDA(NI,NJ_P)            [ W/m/K ]
*       ANTROPGOGENE WAERMEQUELLEN      AWQ(NI,NJ_P)                [ W/ZELLE ]
*
*       TERRESTRISCHE STRAHLUNG         EPSG * (RL-SIGMA*TB**4)   [ W/m**2 ]
*       STEFAN-BOLTZMANN KONSTANTE      SIGMA = 5.6697 E-8
*
*       WAERMEUEBERGANG AUS DER LUFT    WQU(NI,NJ_P)                [ W/ZELLE ]
*       ( OHNE FEUCHTE !!!! )
*       WAERMEKAPAZITAET DES BODENS     CPBOD = 900               [ J/kg/K ]
*
*       GLEICHGEWICHT FUER DIE BODENTEMPERATUR TB(NI,NJ_P,2) IN DER OBERSTEN 
*       BODENZELLE - WAERMELEISTUNGEN IN [ W / m**2 ]
            */

            Parallel.For(2, NI, Program.pOptions, i =>
            {
                int NJ_P = NJ;
                double DYB = Math.Abs(Program.DZB[2] + Program.DZB[3]) * 0.5;
                bool meteopgt_all = ((Program.METEO == "Y") || (Program.METEO == "y")) && (Program.ISTAT == 0);
                //            for (int i = 2; i <= NI - 1; i++)
                //            {
                for (int j = 2; j <= NJ_P - 1; j++)
                {
                    int k = 1;
                    double DDX_DDY = Program.DDXImm[i] * Program.DDYImm[j];
                    float RHO_L = (float)(Program.RHO[i][j][k] * Program.UST[i][j]);
                    float XWQ_L = Program.XWQ[i][j];
                    // double QU_L = Program.QU[i][j][k];
                    // float PBZ_L = Program.PBZ[i][j][k];
                    // double TGI = Program.TB[i][j][2] + 0.5F; //assures that the following loop will be entered
                    // double f1 = DYB * Program.CPBOD * Program.RHOB[i][j] / Program.DT;
                    // double f2 = Program.ALAMBDA[i][j] / DYB;
                    double FAC_L = Program.FAC[i][j][k];
                    //double[] TB_L = Program.TB[i][j];
                    double T = Program.T[i][j][k];
                    double EPSG_L = Program.EPSG[i][j];
                    double RL_L = Program.RL[i][j];
                    double TBA2_L = Program.TBA[i][j][2];
                    ref double WQU = ref Program.WQU[i][j];

                    WQU = XWQ_L * RHO_L * Program.CPLUFT *
                                      (Program.T[i][j][k] + Program.TBZ1 - Program.TB[i][j][2] * FAC_L) *
                                      DDX_DDY;

                    WQU = Math.Min(WQU, 200 * DDX_DDY);
                    WQU = Math.Max(WQU, -500 * DDX_DDY);

                    if (((Program.METEO == "Y") || (Program.METEO == "y")) && (Program.ISTAT == 0) && (Program.AKLA == 4))
                    {
                        WQU = 0;
                    }

                    if (((Program.METEO == "Y") || (Program.METEO == "y")) && (Program.ISTAT == 0) && (Program.ICATAFORC == 1))
                    {
                        if (Program.AKLA == 5)
                        {
                            WQU = 15 * DDX_DDY;
                        }

                        if (Program.AKLA == 6)
                        {
                            WQU = 25 * DDX_DDY;
                        }

                        if (Program.AKLA == 7)
                        {
                            WQU = 30 * DDX_DDY;
                        }
                    }
                }
            });
        }

        static Func<double, int, double> MyPow = (double num, int exp) =>
        {
            double result = 1.0;
            while (exp > 0)
            {
                if (exp % 2 == 1)
                {
                    result *= num;
                }

                exp >>= 1;
                num *= num;
            };
            return result;
        };

        //[MethodImpl(MethodImplOptions.NoInlining)]
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        //[MethodImpl(MethodImplOptions.NoInlining)]
        public static double Pow5(double x1)
        {
            return (x1 * x1 * x1 * x1 * x1);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        //[MethodImpl(MethodImplOptions.NoInlining)]
        public static double Pow4(double x1)
        {
            return (x1 * x1 * x1 * x1);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        //[MethodImpl(MethodImplOptions.NoInlining)]
        public static double Pow3(double x1)
        {
            return (x1 * x1 * x1);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        // [MethodImpl(MethodImplOptions.NoInlining)]
        public static double Pow2(double x1)
        {
            return (x1 * x1);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        // [MethodImpl(MethodImplOptions.NoInlining)]
        public static float Pow2(float x1)
        {
            return (x1 * x1);
        }
    }
}
