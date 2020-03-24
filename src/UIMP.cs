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
using System.Collections.Concurrent;
using System.Runtime.CompilerServices;
using System.Threading.Tasks;

namespace GRAMM_2001
{
    partial class Program
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void UIMPcalculate(int NI, int NJ, int NK)
        {
            if (Program.ICU)
            {
                //computation of Coriolis-term
                Parallel.For(2, NI, Program.pOptions, i =>
                {
                    int NK_P = NK; int NJ_P = NJ;
                    for (int j = 2; j <= NJ_P - 1; j++)
                    {
                        float[] AP0_L = Program.AP0[i][j];
                        double[] U1_L = Program.U1[i][j];
                        double[] U2_L = Program.U2[i][j];
                        double[] DIMU_L = Program.DIMU[i][j];
                        float[] F1U_L = Program.F1U[i][j];
                        float[] F2U_L = Program.F2U[i][j];
                        float[] RHO_L = Program.RHO[i][j];
                        double[] V1N_L = Program.V1N[i][j];
                        float[] VOL_L = Program.VOL[i][j];
                        double[] VG_L = Program.VG[i][j];
                        double[] W1N_L = Program.W1N[i][j];
                        double[] TPDX_L = Program.TPDX[i][j];
                        double[] TPDX_M_L = Program.TPDX[i - 1][j];

                        double f1; double f2 = 0.0;
                        for (int kn = 1; kn <= 2 * (NK_P - 1); kn++)
                        {
                            int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                            f1 = RHO_L[k] * (Program.FN * (V1N_L[k] - VG_L[k]) - Program.FH * W1N_L[k]) * 0.5F * VOL_L[k];
                            f1 += DIMU_L[kn];
                            if (ICPSI)
                            {
                                f1 += TPDX_M_L[k];
                                f2 = TPDX_L[k];
                            }
                            F1U_L[kn] = (float)(AP0_L[k] * U1_L[k] + f1); // Temp array
                            F2U_L[kn] = (float)(AP0_L[k] * U2_L[k] + f1 + f2); // Temp array
                        }
                    }
                });


                int range_parallel = (int)(NI / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
                range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
                range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
                                                               //Iterative solution using an implicit scheme and the TDMA or Thomas-Algorithm
                                                               //Parallel.For(2, NI, Program.pOptions, i =>
                Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
                {
                    int NK_P = NK; int NJ_P = NJ;
                    double DIM;
                    double[] PIM = new double[2 * NK_P];
                    double[] QIM = new double[2 * NK_P];
                    double help;

                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        for (int j = 2; j <= NJ_P - 1; j++)
                        {
                            float relaxv = (float)(Program.RELAXV * Program.Relax_Border_factor[i][j]);
                            float[] AE2_L = Program.AE2[i][j];
                            float[] AN2_L = Program.AN2[i][j];
                            float[] AIM_L = Program.AIM[i][j];
                            float[] AREA_L = Program.AREA[i][j];
                            float[] AS1_L = Program.AS1[i][j];
                            float[] AW1_L = Program.AW1[i][j];
                            float[] BIM_L = Program.BIM[i][j];
                            float[] CIM_L = Program.CIM[i][j];
                            float[] F1U_L = Program.F1U[i][j];
                            float[] F2U_L = Program.F2U[i][j];
                            float[] RHO_L = Program.RHO[i][j];
                            double[] U1N_L = Program.U1N[i][j]; double[] U1Ni_L = Program.U1N[i + 1][j]; double[] U1NJ_P_L = Program.U1N[i][j + 1];
                            double[] U2N_L = Program.U2N[i][j]; double[] U2Ni_L = Program.U2N[i - 1][j]; double[] U2NJ_P_L = Program.U2N[i][j - 1];
                            float USTxUSTV = (float)(Program.UST[i][j] * Program.USTV[i][j]);

                            int m = 2;
                            for (int kn = 1; kn <= 2 * (NK_P - 1); kn++)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                                      //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = AW1_L[k] * U2Ni_L[k] + AS1_L[k] * U2NJ_P_L[k] + F1U_L[kn];
                                    if (k == 1) DIM -= RHO_L[k] * U1N_L[k] * USTxUSTV * AREA_L[k];

                                    //Recurrence formula
                                    if (k > 1)
                                    {
                                        help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                        PIM[kn] = BIM_L[kn] * help;
                                        QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    }
                                    else
                                    {
                                        help = 1 / AIM_L[kn];
                                        PIM[kn] = BIM_L[kn] * help;
                                        QIM[kn] = DIM * help;
                                    }
                                    m--;
                                }

                                //Coefficients for the upper half-cell
                                else
                                {
                                    DIM = AE2_L[k] * U1Ni_L[k] + AN2_L[k] * U1NJ_P_L[k] + F2U_L[kn];
                                    if (k == 1) DIM -= RHO_L[k] * U2N_L[k] * USTxUSTV * AREA_L[k];

                                    //Recurrence formula
                                    help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                    PIM[kn] = BIM_L[kn] * help;
                                    QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    m = 2;
                                }
                            }

                            //Obtain new U-components
                            m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    U1N_L[k] += (relaxv * (PIM[kn] * U2N_L[k] + QIM[kn] - U1N_L[k]));
                                    m--;
                                }
                                else
                                {
                                    U2N_L[k] += (relaxv * (PIM[kn] * U1N_L[k + 1] + QIM[kn] - U2N_L[k]));
                                    m = 2;
                                }

                            }
                        }
                    }
                });

                range_parallel = (int)(NI / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
                range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
                range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
                                                               //Iterative solution using an implicit scheme and the TDMA or Thomas-Algorithm
                                                               //Parallel.For(2, NI, Program.pOptions, ih =>
                Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
                {
                    int NK_P = NK; int NJ_P = NJ;
                    double DIM;
                    double[] PIM = new double[2 * NK_P];
                    double[] QIM = new double[2 * NK_P];
                    double help;
                    for (int ih = range.Item1; ih < range.Item2; ih++)
                    {
                        int i = NI + 1 - ih;
                        for (int j = NJ_P - 1; j >= 2; j--)
                        {
                            float relaxv = (float)(Program.RELAXV * Program.Relax_Border_factor[i][j]);
                            float[] AE2_L = Program.AE2[i][j];
                            float[] AN2_L = Program.AN2[i][j];
                            float[] AIM_L = Program.AIM[i][j];
                            float[] AREA_L = Program.AREA[i][j];
                            float[] AS1_L = Program.AS1[i][j];
                            float[] AW1_L = Program.AW1[i][j];
                            float[] BIM_L = Program.BIM[i][j];
                            float[] CIM_L = Program.CIM[i][j];
                            float[] F1U_L = Program.F1U[i][j];
                            float[] F2U_L = Program.F2U[i][j];
                            float[] RHO_L = Program.RHO[i][j];
                            double[] U1N_L = Program.U1N[i][j]; double[] U1Ni_L = Program.U1N[i + 1][j]; double[] U1NJ_P_L = Program.U1N[i][j + 1];
                            double[] U2N_L = Program.U2N[i][j]; double[] U2Ni_L = Program.U2N[i - 1][j]; double[] U2NJ_P_L = Program.U2N[i][j - 1];
                            float USTxUSTV = (float)(Program.UST[i][j] * Program.USTV[i][j]);

                            int m = 2;
                            for (int kn = 1; kn <= 2 * (NK_P - 1); kn++)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                                      //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = AW1_L[k] * U2Ni_L[k] + AS1_L[k] * U2NJ_P_L[k] + F1U_L[kn];
                                    if (k == 1) DIM -= RHO_L[k] * U1N_L[k] * USTxUSTV * AREA_L[k];

                                    //Recurrence formula
                                    if (k > 1)
                                    {
                                        help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                        PIM[kn] = BIM_L[kn] * help;
                                        QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    }
                                    else
                                    {
                                        help = 1 / AIM_L[kn];
                                        PIM[kn] = BIM_L[kn] * help;
                                        QIM[kn] = DIM * help;
                                    }
                                    m--;
                                }

                                //Coefficients for the upper half-cell
                                else
                                {
                                    DIM = AE2_L[k] * U1Ni_L[k] + AN2_L[k] * U1NJ_P_L[k] + F2U_L[kn];
                                    if (k == 1) DIM -= RHO_L[k] * U2N_L[k] * USTxUSTV * AREA_L[k];


                                    //Recurrence formula
                                    help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                    PIM[kn] = BIM_L[kn] * help;
                                    QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    m = 2;
                                }
                            }

                            //Obtain new U-components
                            m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    U1N_L[k] += (relaxv * (PIM[kn] * U2N_L[k] + QIM[kn] - U1N_L[k]));
                                    m--;
                                }
                                else
                                {
                                    U2N_L[k] += (relaxv * (PIM[kn] * U1N_L[k + 1] + QIM[kn] - U2N_L[k]));
                                    m = 2;
                                }

                            }
                        }
                    }
                });

                range_parallel = (int)(NJ / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
                range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
                range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
                                                               //Iterative solution using an implicit scheme and the TDMA or Thomas-Algorithm
                                                               //Parallel.For(2, NJ, Program.pOptions, jh =>
                Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
                {
                    int NK_P = NK; int NI_P = NI;
                    double DIM;
                    double[] PIM = new double[2 * NK_P];
                    double[] QIM = new double[2 * NK_P];
                    double help;

                    for (int jh = range.Item1; jh < range.Item2; jh++)
                    {
                        int j = NJ + 1 - jh;

                        for (int i = NI_P - 1; i >= 2; i--)
                        {
                            float relaxv = (float)(Program.RELAXV * Program.Relax_Border_factor[i][j]);
                            float[] AE2_L = Program.AE2[i][j];
                            float[] AN2_L = Program.AN2[i][j];
                            float[] AIM_L = Program.AIM[i][j];
                            float[] AREA_L = Program.AREA[i][j];
                            float[] AS1_L = Program.AS1[i][j];
                            float[] AW1_L = Program.AW1[i][j];
                            float[] BIM_L = Program.BIM[i][j];
                            float[] CIM_L = Program.CIM[i][j];
                            float[] F1U_L = Program.F1U[i][j];
                            float[] F2U_L = Program.F2U[i][j];
                            float[] RHO_L = Program.RHO[i][j];
                            double[] U1N_L = Program.U1N[i][j]; double[] U1Ni_L = Program.U1N[i + 1][j]; double[] U1NJ_P_L = Program.U1N[i][j + 1];
                            double[] U2N_L = Program.U2N[i][j]; double[] U2Ni_L = Program.U2N[i - 1][j]; double[] U2NJ_P_L = Program.U2N[i][j - 1];
                            float USTxUSTV = (float)(Program.UST[i][j] * Program.USTV[i][j]);

                            int m = 2;
                            for (int kn = 1; kn <= 2 * (NK_P - 1); kn++)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                                      //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = AW1_L[k] * U2Ni_L[k] + AS1_L[k] * U2NJ_P_L[k] + F1U_L[kn];
                                    if (k == 1) DIM -= RHO_L[k] * U1N_L[k] * USTxUSTV * AREA_L[k];


                                    //Recurrence formula
                                    if (k > 1)
                                    {
                                        help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                        PIM[kn] = BIM_L[kn] * help;
                                        QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    }
                                    else
                                    {
                                        help = 1 / AIM_L[kn];
                                        PIM[kn] = BIM_L[kn] * help;
                                        QIM[kn] = DIM * help;
                                    }
                                    m--;
                                }

                                //Coefficients for the upper half-cell
                                else
                                {
                                    DIM = AE2_L[k] * U1Ni_L[k] + AN2_L[k] * U1NJ_P_L[k] + F2U_L[kn];
                                    if (k == 1) DIM -= RHO_L[k] * U2N_L[k] * USTxUSTV * AREA_L[k];


                                    //Recurrence formula
                                    help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                    PIM[kn] = BIM_L[kn] * help;
                                    QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    m = 2;
                                }
                            }

                            //Obtain new U-components
                            m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    U1N_L[k] += (relaxv * (PIM[kn] * U2N_L[k] + QIM[kn] - U1N_L[k]));
                                    m--;
                                }
                                else
                                {
                                    U2N_L[k] += (relaxv * (PIM[kn] * U1N_L[k + 1] + QIM[kn] - U2N_L[k]));
                                    m = 2;
                                }

                            }
                        }
                    }
                });

                range_parallel = (int)(NJ / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
                range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
                range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
                                                               //Iterative solution using an implicit scheme and the TDMA or Thomas-Algorithm
                                                               //Parallel.For(2, NJ, Program.pOptions, j =>
                Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
                {
                    int NK_P = NK; int NI_P = NI;
                    double DIM;
                    double[] PIM = new double[2 * NK_P];
                    double[] QIM = new double[2 * NK_P];
                    double help;

                    for (int j = range.Item1; j < range.Item2; j++)
                    {
                        for (int i = 2; i <= NI_P - 1; i++)
                        {
                            float relaxv = (float)(Program.RELAXV * Program.Relax_Border_factor[i][j]);
                            float[] AE2_L = Program.AE2[i][j];
                            float[] AN2_L = Program.AN2[i][j];
                            float[] AIM_L = Program.AIM[i][j];
                            float[] AREA_L = Program.AREA[i][j];
                            float[] AS1_L = Program.AS1[i][j];
                            float[] AW1_L = Program.AW1[i][j];
                            float[] BIM_L = Program.BIM[i][j];
                            float[] CIM_L = Program.CIM[i][j];
                            float[] F1U_L = Program.F1U[i][j];
                            float[] F2U_L = Program.F2U[i][j];
                            float[] RHO_L = Program.RHO[i][j];
                            double[] U1N_L = Program.U1N[i][j]; double[] U1Ni_L = Program.U1N[i + 1][j]; double[] U1NJ_P_L = Program.U1N[i][j + 1];
                            double[] U2N_L = Program.U2N[i][j]; double[] U2Ni_L = Program.U2N[i - 1][j]; double[] U2NJ_P_L = Program.U2N[i][j - 1];
                            float USTxUSTV = (float)(Program.UST[i][j] * Program.USTV[i][j]);

                            int m = 2;
                            for (int kn = 1; kn <= 2 * (NK_P - 1); kn++)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                                      //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = AW1_L[k] * U2Ni_L[k] + AS1_L[k] * U2NJ_P_L[k] + F1U_L[kn];
                                    if (k == 1) DIM -= RHO_L[k] * U1N_L[k] * USTxUSTV * AREA_L[k];


                                    //Recurrence formula
                                    if (k > 1)
                                    {
                                        help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                        PIM[kn] = BIM_L[kn] * help;
                                        QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    }
                                    else
                                    {
                                        help = 1 / AIM_L[kn];
                                        PIM[kn] = BIM_L[kn] * help;
                                        QIM[kn] = DIM * help;
                                    }
                                    m--;
                                }

                                //Coefficients for the upper half-cell
                                else
                                {
                                    DIM = AE2_L[k] * U1Ni_L[k] + AN2_L[k] * U1NJ_P_L[k] + F2U_L[kn];
                                    if (k == 1) DIM -= RHO_L[k] * U2N_L[k] * USTxUSTV * AREA_L[k];


                                    //Recurrence formula
                                    help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                    PIM[kn] = BIM_L[kn] * help;
                                    QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    m = 2;
                                }
                            }

                            //Obtain new U-components
                            m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    U1N_L[k] += (relaxv * (PIM[kn] * U2N_L[k] + QIM[kn] - U1N_L[k]));
                                    m--;
                                }
                                else
                                {
                                    U2N_L[k] += (relaxv * (PIM[kn] * U1N_L[k + 1] + QIM[kn] - U2N_L[k]));
                                    m = 2;
                                }

                            }
                        }
                    }
                });

                range_parallel = (int)(NI / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
                range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
                range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
                                                               //Iterative solution using an implicit scheme and the TDMA or Thomas-Algorithm
                                                               //Parallel.For(2, NI, Program.pOptions, ih =>
                Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
                {
                    int NK_P = NK; int NJ_P = NJ;
                    double DIM;
                    double[] PIM = new double[2 * NK_P];
                    double[] QIM = new double[2 * NK_P];
                    double help;

                    for (int ih = range.Item1; ih < range.Item2; ih++)
                    {
                        int i = NI + 1 - ih;
                        for (int j = 2; j <= NJ_P - 1; j++)
                        {
                            float relaxv = (float)(Program.RELAXV * Program.Relax_Border_factor[i][j]);
                            float[] AE2_L = Program.AE2[i][j];
                            float[] AN2_L = Program.AN2[i][j];
                            float[] AIM_L = Program.AIM[i][j];
                            float[] AREA_L = Program.AREA[i][j];
                            float[] AS1_L = Program.AS1[i][j];
                            float[] AW1_L = Program.AW1[i][j];
                            float[] BIM_L = Program.BIM[i][j];
                            float[] CIM_L = Program.CIM[i][j];
                            float[] F1U_L = Program.F1U[i][j];
                            float[] F2U_L = Program.F2U[i][j];
                            float[] RHO_L = Program.RHO[i][j];
                            double[] U1N_L = Program.U1N[i][j]; double[] U1Ni_L = Program.U1N[i + 1][j]; double[] U1NJ_P_L = Program.U1N[i][j + 1];
                            double[] U2N_L = Program.U2N[i][j]; double[] U2Ni_L = Program.U2N[i - 1][j]; double[] U2NJ_P_L = Program.U2N[i][j - 1];
                            float USTxUSTV = (float)(Program.UST[i][j] * Program.USTV[i][j]);

                            int m = 2;
                            for (int kn = 1; kn <= 2 * (NK_P - 1); kn++)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                                      //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = AW1_L[k] * U2Ni_L[k] + AS1_L[k] * U2NJ_P_L[k] + F1U_L[kn];
                                    if (k == 1) DIM -= RHO_L[k] * U1N_L[k] * USTxUSTV * AREA_L[k];


                                    //Recurrence formula
                                    if (k > 1)
                                    {
                                        help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                        PIM[kn] = BIM_L[kn] * help;
                                        QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    }
                                    else
                                    {
                                        help = 1 / AIM_L[kn];
                                        PIM[kn] = BIM_L[kn] * help;
                                        QIM[kn] = DIM * help;
                                    }
                                    m--;
                                }

                                //Coefficients for the upper half-cell
                                else
                                {
                                    DIM = AE2_L[k] * U1Ni_L[k] + AN2_L[k] * U1NJ_P_L[k] + F2U_L[kn];
                                    if (k == 1) DIM -= RHO_L[k] * U2N_L[k] * USTxUSTV * AREA_L[k];


                                    //Recurrence formula
                                    help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                    PIM[kn] = BIM_L[kn] * help;
                                    QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    m = 2;
                                }
                            }

                            //Obtain new U-components
                            m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    U1N_L[k] += (relaxv * (PIM[kn] * U2N_L[k] + QIM[kn] - U1N_L[k]));
                                    m--;
                                }
                                else
                                {
                                    U2N_L[k] += (relaxv * (PIM[kn] * U1N_L[k + 1] + QIM[kn] - U2N_L[k]));
                                    m = 2;
                                }

                            }
                        }
                    }
                });

                range_parallel = (int)(NI / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
                range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
                range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
                                                               //Iterative solution using an implicit scheme and the TDMA or Thomas-Algorithm
                                                               //Parallel.For(2, NI, Program.pOptions, i =>
                Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
                {
                    int NK_P = NK; int NJ_P = NJ;
                    double DIM;
                    double[] PIM = new double[2 * NK_P];
                    double[] QIM = new double[2 * NK_P];
                    double help;
                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        for (int j = NJ_P - 1; j >= 2; j--)
                        {
                            float relaxv = (float)(Program.RELAXV * Program.Relax_Border_factor[i][j]);
                            float[] AE2_L = Program.AE2[i][j];
                            float[] AN2_L = Program.AN2[i][j];
                            float[] AIM_L = Program.AIM[i][j];
                            float[] AREA_L = Program.AREA[i][j];
                            float[] AS1_L = Program.AS1[i][j];
                            float[] AW1_L = Program.AW1[i][j];
                            float[] BIM_L = Program.BIM[i][j];
                            float[] CIM_L = Program.CIM[i][j];
                            float[] F1U_L = Program.F1U[i][j];
                            float[] F2U_L = Program.F2U[i][j];
                            float[] RHO_L = Program.RHO[i][j];
                            double[] U1N_L = Program.U1N[i][j]; double[] U1Ni_L = Program.U1N[i + 1][j]; double[] U1NJ_P_L = Program.U1N[i][j + 1];
                            double[] U2N_L = Program.U2N[i][j]; double[] U2Ni_L = Program.U2N[i - 1][j]; double[] U2NJ_P_L = Program.U2N[i][j - 1];
                            float USTxUSTV = (float)(Program.UST[i][j] * Program.USTV[i][j]);

                            int m = 2;
                            for (int kn = 1; kn <= 2 * (NK_P - 1); kn++)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                                      //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = AW1_L[k] * U2Ni_L[k] + AS1_L[k] * U2NJ_P_L[k] + F1U_L[kn];
                                    if (k == 1) DIM -= RHO_L[k] * U1N_L[k] * USTxUSTV * AREA_L[k];


                                    //Recurrence formula
                                    if (k > 1)
                                    {
                                        help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                        PIM[kn] = BIM_L[kn] * help;
                                        QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    }
                                    else
                                    {
                                        help = 1 / AIM_L[kn];
                                        PIM[kn] = BIM_L[kn] * help;
                                        QIM[kn] = DIM * help;
                                    }
                                    m--;
                                }

                                //Coefficients for the upper half-cell
                                else
                                {
                                    DIM = AE2_L[k] * U1Ni_L[k] + AN2_L[k] * U1NJ_P_L[k] + F2U_L[kn];
                                    if (k == 1) DIM -= RHO_L[k] * U2N_L[k] * USTxUSTV * AREA_L[k];


                                    //Recurrence formula
                                    help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                    PIM[kn] = BIM_L[kn] * help;
                                    QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    m = 2;
                                }
                            }

                            //Obtain new U-components
                            m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    U1N_L[k] += (relaxv * (PIM[kn] * U2N_L[k] + QIM[kn] - U1N_L[k]));
                                    m--;
                                }
                                else
                                {
                                    U2N_L[k] += (relaxv * (PIM[kn] * U1N_L[k + 1] + QIM[kn] - U2N_L[k]));
                                    m = 2;
                                }

                            }
                        }
                    }
                });

                range_parallel = (int)(NJ / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
                range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
                range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
                                                               //Iterative solution using an implicit scheme and the TDMA or Thomas-Algorithm
                                                               //Parallel.For(2, NJ, Program.pOptions, jh =>
                Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
                {
                    int NK_P = NK; int NI_P = NI;
                    double DIM;
                    double[] PIM = new double[2 * NK_P];
                    double[] QIM = new double[2 * NK_P];
                    double help;

                    for (int jh = range.Item1; jh < range.Item2; jh++)
                    {
                        int j = NJ + 1 - jh;
                        for (int i = 2; i <= NI_P - 1; i++)
                        {
                            float relaxv = (float)(Program.RELAXV * Program.Relax_Border_factor[i][j]);
                            float[] AE2_L = Program.AE2[i][j];
                            float[] AN2_L = Program.AN2[i][j];
                            float[] AIM_L = Program.AIM[i][j];
                            float[] AREA_L = Program.AREA[i][j];
                            float[] AS1_L = Program.AS1[i][j];
                            float[] AW1_L = Program.AW1[i][j];
                            float[] BIM_L = Program.BIM[i][j];
                            float[] CIM_L = Program.CIM[i][j];
                            float[] F1U_L = Program.F1U[i][j];
                            float[] F2U_L = Program.F2U[i][j];
                            float[] RHO_L = Program.RHO[i][j];
                            double[] U1N_L = Program.U1N[i][j]; double[] U1Ni_L = Program.U1N[i + 1][j]; double[] U1NJ_P_L = Program.U1N[i][j + 1];
                            double[] U2N_L = Program.U2N[i][j]; double[] U2Ni_L = Program.U2N[i - 1][j]; double[] U2NJ_P_L = Program.U2N[i][j - 1];
                            float USTxUSTV = (float)(Program.UST[i][j] * Program.USTV[i][j]);

                            int m = 2;
                            for (int kn = 1; kn <= 2 * (NK_P - 1); kn++)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                                      //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = AW1_L[k] * U2Ni_L[k] + AS1_L[k] * U2NJ_P_L[k] + F1U_L[kn];
                                    if (k == 1) DIM -= RHO_L[k] * U1N_L[k] * USTxUSTV * AREA_L[k];


                                    //Recurrence formula
                                    if (k > 1)
                                    {
                                        help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                        PIM[kn] = BIM_L[kn] * help;
                                        QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    }
                                    else
                                    {
                                        help = 1 / AIM_L[kn];
                                        PIM[kn] = BIM_L[kn] * help;
                                        QIM[kn] = DIM * help;
                                    }
                                    m--;
                                }

                                //Coefficients for the upper half-cell
                                else
                                {
                                    DIM = AE2_L[k] * U1Ni_L[k] + AN2_L[k] * U1NJ_P_L[k] + F2U_L[kn];
                                    if (k == 1) DIM -= RHO_L[k] * U2N_L[k] * USTxUSTV * AREA_L[k];


                                    //Recurrence formula
                                    help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                    PIM[kn] = BIM_L[kn] * help;
                                    QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    m = 2;
                                }
                            }

                            //Obtain new U-components
                            m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    U1N_L[k] += (relaxv * (PIM[kn] * U2N_L[k] + QIM[kn] - U1N_L[k]));
                                    m--;
                                }
                                else
                                {
                                    U2N_L[k] += (relaxv * (PIM[kn] * U1N_L[k + 1] + QIM[kn] - U2N_L[k]));
                                    m = 2;
                                }

                            }
                        }
                    }
                });

                range_parallel = (int)(NJ / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
                range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
                range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
                                                               //Iterative solution using an implicit scheme and the TDMA or Thomas-Algorithm
                                                               //Parallel.For(2, NJ, Program.pOptions, j =>
                Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
                {
                    int NK_P = NK; int NI_P = NI;
                    double DIM;
                    double[] PIM = new double[2 * NK_P];
                    double[] QIM = new double[2 * NK_P];
                    double help;

                    for (int j = range.Item1; j < range.Item2; j++)
                    {
                        for (int i = 2; i <= NI_P - 1; i++)
                        {
                            float relaxv = (float)(Program.RELAXV * Program.Relax_Border_factor[i][j]);
                            float[] AE2_L = Program.AE2[i][j];
                            float[] AN2_L = Program.AN2[i][j];
                            float[] AIM_L = Program.AIM[i][j];
                            float[] AREA_L = Program.AREA[i][j];
                            float[] AS1_L = Program.AS1[i][j];
                            float[] AW1_L = Program.AW1[i][j];
                            float[] BIM_L = Program.BIM[i][j];
                            float[] CIM_L = Program.CIM[i][j];
                            float[] F1U_L = Program.F1U[i][j];
                            float[] F2U_L = Program.F2U[i][j];
                            float[] RHO_L = Program.RHO[i][j];
                            double[] U1N_L = Program.U1N[i][j]; double[] U1Ni_L = Program.U1N[i + 1][j]; double[] U1NJ_P_L = Program.U1N[i][j + 1];
                            double[] U2N_L = Program.U2N[i][j]; double[] U2Ni_L = Program.U2N[i - 1][j]; double[] U2NJ_P_L = Program.U2N[i][j - 1];
                            float USTxUSTV = (float)(Program.UST[i][j] * Program.USTV[i][j]);

                            int m = 2;
                            for (int kn = 1; kn <= 2 * (NK_P - 1); kn++)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                                      //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = AW1_L[k] * U2Ni_L[k] + AS1_L[k] * U2NJ_P_L[k] + F1U_L[kn];
                                    if (k == 1) DIM -= RHO_L[k] * U1N_L[k] * USTxUSTV * AREA_L[k];


                                    //Recurrence formula
                                    if (k > 1)
                                    {
                                        help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                        PIM[kn] = BIM_L[kn] * help;
                                        QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    }
                                    else
                                    {
                                        help = 1 / AIM_L[kn];
                                        PIM[kn] = BIM_L[kn] * help;
                                        QIM[kn] = DIM * help;
                                    }
                                    m--;
                                }

                                //Coefficients for the upper half-cell
                                else
                                {
                                    DIM = AE2_L[k] * U1Ni_L[k] + AN2_L[k] * U1NJ_P_L[k] + F2U_L[kn];
                                    if (k == 1) DIM -= RHO_L[k] * U2N_L[k] * USTxUSTV * AREA_L[k];


                                    //Recurrence formula
                                    help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                    PIM[kn] = BIM_L[kn] * help;
                                    QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    m = 2;
                                }
                            }

                            //Obtain new U-components
                            m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    U1N_L[k] += (relaxv * (PIM[kn] * U2N_L[k] + QIM[kn] - U1N_L[k]));
                                    m--;
                                }
                                else
                                {
                                    U2N_L[k] += (relaxv * (PIM[kn] * U1N_L[k + 1] + QIM[kn] - U2N_L[k]));
                                    m = 2;
                                }

                            }
                        }
                    }
                });
            }
        }
    }
}
