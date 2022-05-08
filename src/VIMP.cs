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
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        public static void VIMPcalculate(int NI, int NJ, int NK)
        {
            if (Program.ICV)
            {
                //computation of Coriolis-term
                Parallel.For(2, NI, Program.pOptions, i =>
                {
                    int NK_P = NK; int NJ_P = NJ;
                    for (int j = 2; j <= NJ_P - 1; j++)
                    {
                        float[] AP0_L = Program.AP0[i][j];
                        double[] V1_L = Program.V1[i][j];
                        double[] V2_L = Program.V2[i][j];
                        double[] DIMV_L = Program.DIMV[i][j];
                        float[] F1V_L = Program.F1V[i][j];
                        float[] F2V_L = Program.F2V[i][j];
                        ReadOnlySpan<float> RHO_L = Program.RHO[i][j];
                        double[] U1N_L = Program.U1N[i][j];
                        double[] UG_L = Program.UG[i][j];
                        ReadOnlySpan<float> VOL_L = Program.VOLImm[i][j].AsSpan();
                        double[] TPDY_L = Program.TPDY[i][j];
                        double[] TPDY_M_L = Program.TPDY[i][j - 1];

                        double f1; double f2 = 0.0;

                        for (int kn = 1; kn <= 2 * (NK_P - 1); kn++)
                        {
                            int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                            f1 = -RHO_L[k] * Program.FN * (U1N_L[k] - UG_L[k]) * 0.5F * VOL_L[k];
                            f1 += DIMV_L[kn];
                            if (ICPSI)
                            {
                                f1 += TPDY_M_L[k];
                                f2 = TPDY_L[k];
                            }
                            F1V_L[kn] = (float)(AP0_L[k] * V1_L[k] + f1); // Temp array
                            F2V_L[kn] = (float)(AP0_L[k] * V2_L[k] + f1 + f2); // Temp array
                        }
                    }
                });

                int range_parallel = NI / Program.pOptions.MaxDegreeOfParallelism - (StripeCounter % 6);
                range_parallel = Math.Max(Program.StripeWidth - (StripeCounter % 6), range_parallel); // min. Program.StripeWidth cells per processor
                StripeCounter++;
                range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
                                                               //Iterative solution using an implicit scheme and the TDMA or Thomas-Algorithm
                                                               //Parallel.For(2, NI, Program.pOptions, i =>
                Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
                {
                    int NK_P = NK; int NJ_P = NJ;
                    double DIM;
                    Span<double> PIM = stackalloc double[2 * NK_P];
                    Span<double> QIM = stackalloc double[2 * NK_P];
                    double help;
                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        for (int j = 2; j <= NJ_P - 1; j++)
                        {
                            float relaxv = (float)(Program.RELAXV * Program.Relax_Border_factor[i][j]);
                            float[] AE2_L = Program.AE2[i][j];
                            float[] AIM_L = Program.AIM[i][j];
                            float[] AN2_L = Program.AN2[i][j];
                            float[] AS1_L = Program.AS1[i][j];
                            ReadOnlySpan<float> AREA_L = Program.AREAImm[i][j].AsSpan();
                            float[] AW1_L = Program.AW1[i][j];
                            float[] BIM_L = Program.BIM[i][j];
                            float[] CIM_L = Program.CIM[i][j];
                            float[] F1V_L = Program.F1V[i][j];
                            float[] F2V_L = Program.F2V[i][j];
                            ReadOnlySpan<float> RHO_L = Program.RHO[i][j];
                            double[] V1N_L = Program.V1N[i][j]; double[] V1Ni_L = Program.V1N[i + 1][j]; double[] V1NJ_P_L = Program.V1N[i][j + 1];
                            double[] V2N_L = Program.V2N[i][j]; double[] V2Ni_L = Program.V2N[i - 1][j]; double[] V2NJ_P_L = Program.V2N[i][j - 1];

                            float USTxUSTV = (float)(Program.UST[i][j] * Program.USTV[i][j]);

                            int m = 2;
                            for (int kn = 1; kn <= 2 * (NK_P - 1); kn++)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                                      //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = AW1_L[k] * V2Ni_L[k] + AS1_L[k] * V2NJ_P_L[k] + F1V_L[kn];

                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * V1N_L[k] * USTxUSTV * AREA_L[k];
                                    }


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
                                    DIM = AE2_L[k] * V1Ni_L[k] + AN2_L[k] * V1NJ_P_L[k] + F2V_L[kn];

                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * V2N_L[k] * USTxUSTV * AREA_L[k];
                                    }


                                    //Recurrence formula
                                    help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                    PIM[kn] = BIM_L[kn] * help;
                                    QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    m = 2;
                                }
                            }

                            //Obtain new V-components
                            m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    V1N_L[k] += (relaxv * (PIM[kn] * V2N_L[k] + QIM[kn] - V1N_L[k]));
                                    m--;
                                }
                                else
                                {
                                    V2N_L[k] += (relaxv * (PIM[kn] * V1N_L[k + 1] + QIM[kn] - V2N_L[k]));
                                    m = 2;
                                }

                            }
                        }
                    }
                });

                range_parallel = NI / Program.pOptions.MaxDegreeOfParallelism - (StripeCounter % 6);
                range_parallel = Math.Max(Program.StripeWidth - (StripeCounter % 6), range_parallel); // min. Program.StripeWidth cells per processor
                StripeCounter++;
                range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
                                                               //Iterative solution using an implicit scheme and the TDMA or Thomas-Algorithm
                                                               //Parallel.For(2, NI, Program.pOptions, ih =>
                Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
                {
                    int NK_P = NK; int NJ_P = NJ;
                    double DIM;
                    Span<double> PIM = stackalloc double[2 * NK_P];
                    Span<double> QIM = stackalloc double[2 * NK_P];
                    double help;

                    for (int ih = range.Item1; ih < range.Item2; ih++)
                    {
                        int i = NI + 1 - ih;
                        for (int j = NJ_P - 1; j >= 2; j--)
                        {
                            float relaxv = (float)(Program.RELAXV * Program.Relax_Border_factor[i][j]);
                            float[] AE2_L = Program.AE2[i][j];
                            float[] AIM_L = Program.AIM[i][j];
                            float[] AN2_L = Program.AN2[i][j];
                            float[] AS1_L = Program.AS1[i][j];
                            ReadOnlySpan<float> AREA_L = Program.AREAImm[i][j].AsSpan();
                            float[] AW1_L = Program.AW1[i][j];
                            float[] BIM_L = Program.BIM[i][j];
                            float[] CIM_L = Program.CIM[i][j];
                            float[] F1V_L = Program.F1V[i][j];
                            float[] F2V_L = Program.F2V[i][j];
                            ReadOnlySpan<float> RHO_L = Program.RHO[i][j];
                            double[] V1N_L = Program.V1N[i][j]; double[] V1Ni_L = Program.V1N[i + 1][j]; double[] V1NJ_P_L = Program.V1N[i][j + 1];
                            double[] V2N_L = Program.V2N[i][j]; double[] V2Ni_L = Program.V2N[i - 1][j]; double[] V2NJ_P_L = Program.V2N[i][j - 1];

                            float USTxUSTV = (float)(Program.UST[i][j] * Program.USTV[i][j]);

                            int m = 2;
                            for (int kn = 1; kn <= 2 * (NK_P - 1); kn++)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                                      //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = AW1_L[k] * V2Ni_L[k] + AS1_L[k] * V2NJ_P_L[k] + F1V_L[kn];

                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * V1N_L[k] * USTxUSTV * AREA_L[k];
                                    }

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
                                    DIM = AE2_L[k] * V1Ni_L[k] + AN2_L[k] * V1NJ_P_L[k] + F2V_L[kn];

                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * V2N_L[k] * USTxUSTV * AREA_L[k];
                                    }



                                    //Recurrence formula
                                    help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                    PIM[kn] = BIM_L[kn] * help;
                                    QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    m = 2;
                                }
                            }

                            //Obtain new V-components
                            m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    V1N_L[k] += (relaxv * (PIM[kn] * V2N_L[k] + QIM[kn] - V1N_L[k]));
                                    m--;
                                }
                                else
                                {
                                    V2N_L[k] += (relaxv * (PIM[kn] * V1N_L[k + 1] + QIM[kn] - V2N_L[k]));
                                    m = 2;
                                }

                            }
                        }
                    }
                });

                range_parallel = NJ / Program.pOptions.MaxDegreeOfParallelism - (StripeCounter % 6);
                range_parallel = Math.Max(Program.StripeWidth - (StripeCounter % 6), range_parallel); // min. Program.StripeWidth cells per processor
                StripeCounter++;
                range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
                                                               //Iterative solution using an implicit scheme and the TDMA or Thomas-Algorithm
                                                               //Parallel.For(2, NJ, Program.pOptions, jh =>
                Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
                {
                    int NK_P = NK; int NI_P = NI;
                    double DIM;
                    Span<double> PIM = stackalloc double[2 * NK_P];
                    Span<double> QIM = stackalloc double[2 * NK_P];
                    double help;

                    for (int jh = range.Item1; jh < range.Item2; jh++)
                    {
                        int j = NJ + 1 - jh;
                        for (int i = NI_P - 1; i >= 2; i--)
                        {
                            float relaxv = (float)(Program.RELAXV * Program.Relax_Border_factor[i][j]);
                            float[] AE2_L = Program.AE2[i][j];
                            float[] AIM_L = Program.AIM[i][j];
                            float[] AN2_L = Program.AN2[i][j];
                            float[] AS1_L = Program.AS1[i][j];
                            ReadOnlySpan<float> AREA_L = Program.AREAImm[i][j].AsSpan();
                            float[] AW1_L = Program.AW1[i][j];
                            float[] BIM_L = Program.BIM[i][j];
                            float[] CIM_L = Program.CIM[i][j];
                            float[] F1V_L = Program.F1V[i][j];
                            float[] F2V_L = Program.F2V[i][j];
                            ReadOnlySpan<float> RHO_L = Program.RHO[i][j];
                            double[] V1N_L = Program.V1N[i][j]; double[] V1Ni_L = Program.V1N[i + 1][j]; double[] V1NJ_P_L = Program.V1N[i][j + 1];
                            double[] V2N_L = Program.V2N[i][j]; double[] V2Ni_L = Program.V2N[i - 1][j]; double[] V2NJ_P_L = Program.V2N[i][j - 1];

                            float USTxUSTV = (float)(Program.UST[i][j] * Program.USTV[i][j]);

                            int m = 2;
                            for (int kn = 1; kn <= 2 * (NK_P - 1); kn++)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                                      //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = AW1_L[k] * V2Ni_L[k] + AS1_L[k] * V2NJ_P_L[k] + F1V_L[kn];

                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * V1N_L[k] * USTxUSTV * AREA_L[k];
                                    }

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
                                    DIM = AE2_L[k] * V1Ni_L[k] + AN2_L[k] * V1NJ_P_L[k] + F2V_L[kn];

                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * V2N_L[k] * USTxUSTV * AREA_L[k];
                                    }


                                    //Recurrence formula
                                    help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                    PIM[kn] = BIM_L[kn] * help;
                                    QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    m = 2;
                                }
                            }

                            //Obtain new V-components
                            m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    V1N_L[k] += (relaxv * (PIM[kn] * V2N_L[k] + QIM[kn] - V1N_L[k]));
                                    m--;
                                }
                                else
                                {
                                    V2N_L[k] += (relaxv * (PIM[kn] * V1N_L[k + 1] + QIM[kn] - V2N_L[k]));
                                    m = 2;
                                }

                            }
                        }
                    }
                });

                range_parallel = NJ / Program.pOptions.MaxDegreeOfParallelism - (StripeCounter % 6);
                range_parallel = Math.Max(Program.StripeWidth - (StripeCounter % 6), range_parallel); // min. Program.StripeWidth cells per processor
                StripeCounter++;
                range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
                                                               //Iterative solution using an implicit scheme and the TDMA or Thomas-Algorithm
                                                               //Parallel.For(2, NJ, Program.pOptions, j =>
                Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
                {
                    int NK_P = NK; int NI_P = NI;
                    double DIM;
                    Span<double> PIM = stackalloc double[2 * NK_P];
                    Span<double> QIM = stackalloc double[2 * NK_P];
                    double help;

                    for (int j = range.Item1; j < range.Item2; j++)
                    {
                        for (int i = 2; i <= NI_P - 1; i++)
                        {
                            float relaxv = (float)(Program.RELAXV * Program.Relax_Border_factor[i][j]);
                            float[] AE2_L = Program.AE2[i][j];
                            float[] AIM_L = Program.AIM[i][j];
                            float[] AN2_L = Program.AN2[i][j];
                            float[] AS1_L = Program.AS1[i][j];
                            ReadOnlySpan<float> AREA_L = Program.AREAImm[i][j].AsSpan();
                            float[] AW1_L = Program.AW1[i][j];
                            float[] BIM_L = Program.BIM[i][j];
                            float[] CIM_L = Program.CIM[i][j];
                            float[] F1V_L = Program.F1V[i][j];
                            float[] F2V_L = Program.F2V[i][j];
                            ReadOnlySpan<float> RHO_L = Program.RHO[i][j];
                            double[] V1N_L = Program.V1N[i][j]; double[] V1Ni_L = Program.V1N[i + 1][j]; double[] V1NJ_P_L = Program.V1N[i][j + 1];
                            double[] V2N_L = Program.V2N[i][j]; double[] V2Ni_L = Program.V2N[i - 1][j]; double[] V2NJ_P_L = Program.V2N[i][j - 1];

                            float USTxUSTV = (float)(Program.UST[i][j] * Program.USTV[i][j]);

                            int m = 2;
                            for (int kn = 1; kn <= 2 * (NK_P - 1); kn++)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                                      //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = AW1_L[k] * V2Ni_L[k] + AS1_L[k] * V2NJ_P_L[k] + F1V_L[kn];

                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * V1N_L[k] * USTxUSTV * AREA_L[k];
                                    }

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
                                    DIM = AE2_L[k] * V1Ni_L[k] + AN2_L[k] * V1NJ_P_L[k] + F2V_L[kn];

                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * V2N_L[k] * USTxUSTV * AREA_L[k];
                                    }


                                    //Recurrence formula
                                    help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                    PIM[kn] = BIM_L[kn] * help;
                                    QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    m = 2;
                                }
                            }

                            //Obtain new V-components
                            m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    V1N_L[k] += (relaxv * (PIM[kn] * V2N_L[k] + QIM[kn] - V1N_L[k]));
                                    m--;
                                }
                                else
                                {
                                    V2N_L[k] += (relaxv * (PIM[kn] * V1N_L[k + 1] + QIM[kn] - V2N_L[k]));
                                    m = 2;
                                }

                            }
                        }
                    }
                });

                range_parallel = NI / Program.pOptions.MaxDegreeOfParallelism - (StripeCounter % 6);
                range_parallel = Math.Max(Program.StripeWidth - (StripeCounter % 6), range_parallel); // min. Program.StripeWidth cells per processor
                StripeCounter++;
                range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
                                                               //Iterative solution using an implicit scheme and the TDMA or Thomas-Algorithm
                                                               //Parallel.For(2, NI, Program.pOptions, ih =>
                Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
                {
                    int NK_P = NK; int NJ_P = NJ;
                    double DIM;
                    Span<double> PIM = stackalloc double[2 * NK_P];
                    Span<double> QIM = stackalloc double[2 * NK_P];
                    double help;

                    for (int ih = range.Item1; ih < range.Item2; ih++)
                    {
                        int i = NI + 1 - ih;
                        for (int j = 2; j <= NJ_P - 1; j++)
                        {
                            float relaxv = (float)(Program.RELAXV * Program.Relax_Border_factor[i][j]);
                            float[] AE2_L = Program.AE2[i][j];
                            float[] AIM_L = Program.AIM[i][j];
                            float[] AN2_L = Program.AN2[i][j];
                            float[] AS1_L = Program.AS1[i][j];
                            ReadOnlySpan<float> AREA_L = Program.AREAImm[i][j].AsSpan();
                            float[] AW1_L = Program.AW1[i][j];
                            float[] BIM_L = Program.BIM[i][j];
                            float[] CIM_L = Program.CIM[i][j];
                            float[] F1V_L = Program.F1V[i][j];
                            float[] F2V_L = Program.F2V[i][j];
                            ReadOnlySpan<float> RHO_L = Program.RHO[i][j];
                            double[] V1N_L = Program.V1N[i][j]; double[] V1Ni_L = Program.V1N[i + 1][j]; double[] V1NJ_P_L = Program.V1N[i][j + 1];
                            double[] V2N_L = Program.V2N[i][j]; double[] V2Ni_L = Program.V2N[i - 1][j]; double[] V2NJ_P_L = Program.V2N[i][j - 1];

                            float USTxUSTV = (float)(Program.UST[i][j] * Program.USTV[i][j]);

                            int m = 2;
                            for (int kn = 1; kn <= 2 * (NK_P - 1); kn++)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                                      //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = AW1_L[k] * V2Ni_L[k] + AS1_L[k] * V2NJ_P_L[k] + F1V_L[kn];

                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * V1N_L[k] * USTxUSTV * AREA_L[k];
                                    }

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
                                    DIM = AE2_L[k] * V1Ni_L[k] + AN2_L[k] * V1NJ_P_L[k] + F2V_L[kn];

                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * V2N_L[k] * USTxUSTV * AREA_L[k];
                                    }


                                    //Recurrence formula
                                    help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                    PIM[kn] = BIM_L[kn] * help;
                                    QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    m = 2;
                                }
                            }

                            //Obtain new V-components
                            m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    V1N_L[k] += (relaxv * (PIM[kn] * V2N_L[k] + QIM[kn] - V1N_L[k]));
                                    m--;
                                }
                                else
                                {
                                    V2N_L[k] += (relaxv * (PIM[kn] * V1N_L[k + 1] + QIM[kn] - V2N_L[k]));
                                    m = 2;
                                }

                            }
                        }
                    }
                });

                range_parallel = NI / Program.pOptions.MaxDegreeOfParallelism - (StripeCounter % 6);
                range_parallel = Math.Max(Program.StripeWidth - (StripeCounter % 6), range_parallel); // min. Program.StripeWidth cells per processor
                StripeCounter++;
                range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
                                                               //Iterative solution using an implicit scheme and the TDMA or Thomas-Algorithm
                                                               //Parallel.For(2, NI, Program.pOptions, i =>
                Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
                {
                    int NK_P = NK; int NJ_P = NJ;
                    double DIM;
                    Span<double> PIM = stackalloc double[2 * NK_P];
                    Span<double> QIM = stackalloc double[2 * NK_P];
                    double help;
                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        for (int j = NJ_P - 1; j >= 2; j--)
                        {
                            float relaxv = (float)(Program.RELAXV * Program.Relax_Border_factor[i][j]);
                            float[] AE2_L = Program.AE2[i][j];
                            float[] AIM_L = Program.AIM[i][j];
                            float[] AN2_L = Program.AN2[i][j];
                            float[] AS1_L = Program.AS1[i][j];
                            ReadOnlySpan<float> AREA_L = Program.AREAImm[i][j].AsSpan();
                            float[] AW1_L = Program.AW1[i][j];
                            float[] BIM_L = Program.BIM[i][j];
                            float[] CIM_L = Program.CIM[i][j];
                            float[] F1V_L = Program.F1V[i][j];
                            float[] F2V_L = Program.F2V[i][j];
                            ReadOnlySpan<float> RHO_L = Program.RHO[i][j];
                            double[] V1N_L = Program.V1N[i][j]; double[] V1Ni_L = Program.V1N[i + 1][j]; double[] V1NJ_P_L = Program.V1N[i][j + 1];
                            double[] V2N_L = Program.V2N[i][j]; double[] V2Ni_L = Program.V2N[i - 1][j]; double[] V2NJ_P_L = Program.V2N[i][j - 1];

                            float USTxUSTV = (float)(Program.UST[i][j] * Program.USTV[i][j]);

                            int m = 2;
                            for (int kn = 1; kn <= 2 * (NK_P - 1); kn++)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                                      //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = AW1_L[k] * V2Ni_L[k] + AS1_L[k] * V2NJ_P_L[k] + F1V_L[kn];

                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * V1N_L[k] * USTxUSTV * AREA_L[k];
                                    }

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
                                    DIM = AE2_L[k] * V1Ni_L[k] + AN2_L[k] * V1NJ_P_L[k] + F2V_L[kn];

                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * V2N_L[k] * USTxUSTV * AREA_L[k];
                                    }


                                    //Recurrence formula
                                    help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                    PIM[kn] = BIM_L[kn] * help;
                                    QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    m = 2;
                                }
                            }

                            //Obtain new V-components
                            m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    V1N_L[k] += (relaxv * (PIM[kn] * V2N_L[k] + QIM[kn] - V1N_L[k]));
                                    m--;
                                }
                                else
                                {
                                    V2N_L[k] += (relaxv * (PIM[kn] * V1N_L[k + 1] + QIM[kn] - V2N_L[k]));
                                    m = 2;
                                }

                            }
                        }
                    }
                });

                range_parallel = NJ / Program.pOptions.MaxDegreeOfParallelism - (StripeCounter % 6);
                range_parallel = Math.Max(Program.StripeWidth - (StripeCounter % 6), range_parallel); // min. Program.StripeWidth cells per processor
                StripeCounter++;
                range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
                                                               //Iterative solution using an implicit scheme and the TDMA or Thomas-Algorithm
                                                               //Parallel.For(2, NJ, Program.pOptions, jh =>
                Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
                {
                    int NK_P = NK; int NI_P = NI;
                    double DIM;
                    Span<double> PIM = stackalloc double[2 * NK_P];
                    Span<double> QIM = stackalloc double[2 * NK_P];
                    double help;

                    for (int jh = range.Item1; jh < range.Item2; jh++)
                    {
                        int j = NJ + 1 - jh;
                        for (int i = 2; i <= NI_P - 1; i++)
                        {
                            float relaxv = (float)(Program.RELAXV * Program.Relax_Border_factor[i][j]);
                            float[] AE2_L = Program.AE2[i][j];
                            float[] AIM_L = Program.AIM[i][j];
                            float[] AN2_L = Program.AN2[i][j];
                            float[] AS1_L = Program.AS1[i][j];
                            ReadOnlySpan<float> AREA_L = Program.AREAImm[i][j].AsSpan();
                            float[] AW1_L = Program.AW1[i][j];
                            float[] BIM_L = Program.BIM[i][j];
                            float[] CIM_L = Program.CIM[i][j];
                            float[] F1V_L = Program.F1V[i][j];
                            float[] F2V_L = Program.F2V[i][j];
                            ReadOnlySpan<float> RHO_L = Program.RHO[i][j];
                            double[] V1N_L = Program.V1N[i][j]; double[] V1Ni_L = Program.V1N[i + 1][j]; double[] V1NJ_P_L = Program.V1N[i][j + 1];
                            double[] V2N_L = Program.V2N[i][j]; double[] V2Ni_L = Program.V2N[i - 1][j]; double[] V2NJ_P_L = Program.V2N[i][j - 1];

                            float USTxUSTV = (float)(Program.UST[i][j] * Program.USTV[i][j]);

                            int m = 2;
                            for (int kn = 1; kn <= 2 * (NK_P - 1); kn++)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                                      //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = AW1_L[k] * V2Ni_L[k] + AS1_L[k] * V2NJ_P_L[k] + F1V_L[kn];

                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * V1N_L[k] * USTxUSTV * AREA_L[k];
                                    }

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
                                    DIM = AE2_L[k] * V1Ni_L[k] + AN2_L[k] * V1NJ_P_L[k] + F2V_L[kn];

                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * V2N_L[k] * USTxUSTV * AREA_L[k];
                                    }


                                    //Recurrence formula
                                    help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                    PIM[kn] = BIM_L[kn] * help;
                                    QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    m = 2;
                                }
                            }

                            //Obtain new V-components
                            m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    V1N_L[k] += (relaxv * (PIM[kn] * V2N_L[k] + QIM[kn] - V1N_L[k]));
                                    m--;
                                }
                                else
                                {
                                    V2N_L[k] += (relaxv * (PIM[kn] * V1N_L[k + 1] + QIM[kn] - V2N_L[k]));
                                    m = 2;
                                }

                            }
                        }
                    }
                });

                range_parallel = NJ / Program.pOptions.MaxDegreeOfParallelism - (StripeCounter % 6);
                range_parallel = Math.Max(Program.StripeWidth - (StripeCounter % 6), range_parallel); // min. Program.StripeWidth cells per processor
                StripeCounter++;
                range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
                                                               //Iterative solution using an implicit scheme and the TDMA or Thomas-Algorithm
                                                               //Parallel.For(2, NJ, Program.pOptions, j =>
                Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
                {
                    int NK_P = NK; int NI_P = NI;
                    double DIM;
                    Span<double> PIM = stackalloc double[2 * NK_P];
                    Span<double> QIM = stackalloc double[2 * NK_P];
                    double help;
                    for (int j = range.Item1; j < range.Item2; j++)
                    {
                        for (int i = 2; i <= NI_P - 1; i++)
                        {
                            float relaxv = (float)(Program.RELAXV * Program.Relax_Border_factor[i][j]);
                            float[] AE2_L = Program.AE2[i][j];
                            float[] AIM_L = Program.AIM[i][j];
                            float[] AN2_L = Program.AN2[i][j];
                            float[] AS1_L = Program.AS1[i][j];
                            ReadOnlySpan<float> AREA_L = Program.AREAImm[i][j].AsSpan();
                            float[] AW1_L = Program.AW1[i][j];
                            float[] BIM_L = Program.BIM[i][j];
                            float[] CIM_L = Program.CIM[i][j];
                            float[] F1V_L = Program.F1V[i][j];
                            float[] F2V_L = Program.F2V[i][j];
                            ReadOnlySpan<float> RHO_L = Program.RHO[i][j];
                            double[] V1N_L = Program.V1N[i][j]; double[] V1Ni_L = Program.V1N[i + 1][j]; double[] V1NJ_P_L = Program.V1N[i][j + 1];
                            double[] V2N_L = Program.V2N[i][j]; double[] V2Ni_L = Program.V2N[i - 1][j]; double[] V2NJ_P_L = Program.V2N[i][j - 1];

                            float USTxUSTV = (float)(Program.UST[i][j] * Program.USTV[i][j]);

                            int m = 2;
                            for (int kn = 1; kn <= 2 * (NK_P - 1); kn++)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                                      //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = AW1_L[k] * V2Ni_L[k] + AS1_L[k] * V2NJ_P_L[k] + F1V_L[kn];

                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * V1N_L[k] * USTxUSTV * AREA_L[k];
                                    }

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
                                    DIM = AE2_L[k] * V1Ni_L[k] + AN2_L[k] * V1NJ_P_L[k] + F2V_L[kn];

                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * V2N_L[k] * USTxUSTV * AREA_L[k];
                                    }

                                    //Recurrence formula
                                    help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                    PIM[kn] = BIM_L[kn] * help;
                                    QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    m = 2;
                                }
                            }

                            //Obtain new V-components
                            m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    V1N_L[k] += (relaxv * (PIM[kn] * V2N_L[k] + QIM[kn] - V1N_L[k]));
                                    m--;
                                }
                                else
                                {
                                    V2N_L[k] += (relaxv * (PIM[kn] * V1N_L[k + 1] + QIM[kn] - V2N_L[k]));
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
