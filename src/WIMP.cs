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
        public static void WVELcalculate(int NI, int NJ, int NK)
        {
            if (Program.ICW)
            {
                //computation of Coriolis-term  
                Parallel.For(2, NI, Program.pOptions, i =>
                {
                    int NK_P = NK; int NJ_P = NJ;
                    for (int j = 2; j <= NJ_P - 1; j++)
                    {
                        float[] AP0_L = Program.AP0[i][j];
                        double[] W2_L = Program.W2[i][j];
                        double[] DIMW_L = Program.DIMW[i][j];
                        double[] W1_L = Program.W1[i][j];
                        float[] F1W_L = Program.F1W[i][j];
                        float[] F2W_L = Program.F2W[i][j];
                        double[] QBZ_L = Program.QBZ[i][j];
                        double[] QUN_L = Program.QUN[i][j];
                        ReadOnlySpan<float> RHO_L = Program.RHO[i][j];
                        float[] RHOBZ_L = Program.RHOBZ[i][j];
                        double[] TBZ_L = Program.TBZ[i][j];
                        double[] TN_L = Program.TN[i][j];
                        double[] U1N_L = Program.U1N[i][j];
                        ReadOnlySpan<float> VOL_L = Program.VOLImm[i][j].AsSpan();
                        double f1, f2;

                        for (int kn = 1; kn <= 2 * (NK_P - 1); kn++)
                        {
                            int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)

                            f1 = RHO_L[k] * Program.FH * U1N_L[k] * 0.5F * VOL_L[k];
                            f2 = 0;
                            if (Program.ICPSI == false)
                            {
                                f2 = 0.5F * VOL_L[k] * ((TN_L[k] + Program.TBZ1) * (1 + 0.00061 * QUN_L[k]) - (TBZ_L[k] * (1 + 0.00061 * QBZ_L[k]))) /
                                (TBZ_L[k] * (1 + 0.00061 * QBZ_L[k])) * RHOBZ_L[k] * Program.GERD;
                            }


                            //no bouyancy, where terrain is smoothed at the border
                            /*
                            if ((i < nr_cell_smooth) || (i > NI - nr_cell_smooth)||(j < nr_cell_smooth) || (j > NJ - nr_cell_smooth))
                            {
                                f2 = 0;
                            }
                            */

                            f1 += f2;
                            f1 += DIMW_L[kn];
                            F1W_L[kn] = (float)(AP0_L[k] * W1_L[k] + f1);
                            F2W_L[kn] = (float)(AP0_L[k] * W2_L[k] + f1);
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
                            float[] AN2_L = Program.AN2[i][j];
                            float[] AIM_L = Program.AIM[i][j];
                            ReadOnlySpan<float> AREAZX_L = Program.AREAZXImm[i][j].AsSpan();
                            ReadOnlySpan<float> AREAZY_L = Program.AREAZYImm[i][j].AsSpan();
                            float[] AS1_L = Program.AS1[i][j];
                            float[] AW1_L = Program.AW1[i][j];
                            float[] BIM_L = Program.BIM[i][j];
                            float[] CIM_L = Program.CIM[i][j];
                            float[] F1W_L = Program.F1W[i][j];
                            float[] F2W_L = Program.F2W[i][j];
                            ReadOnlySpan<float> RHO_L = Program.RHO[i][j];
                            double[] W1N_L = Program.W1N[i][j]; double[] W1Ni_L = Program.W1N[i + 1][j]; double[] W1NJ_P_L = Program.W1N[i][j + 1];
                            double[] W2N_L = Program.W2N[i][j]; double[] W2Ni_L = Program.W2N[i - 1][j]; double[] W2NJ_P_L = Program.W2N[i][j - 1];
                            float USTxUSTV = (float)(Program.UST[i][j] * Program.USTV[i][j]);

                            int m = 2;
                            for (int kn = 1; kn <= 2 * (NK_P - 1); kn++)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                                      //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = AW1_L[k] * W2Ni_L[k] + AS1_L[k] * W2NJ_P_L[k] + F1W_L[kn];
                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * W1N_L[k] * USTxUSTV * MathF.Sqrt(Pow2(AREAZX_L[k]) + Pow2(AREAZY_L[k]));
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
                                    DIM = AE2_L[k] * W1Ni_L[k] + AN2_L[k] * W1NJ_P_L[k] + F2W_L[kn];
                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * W2N_L[k] * USTxUSTV * MathF.Sqrt(Pow2(AREAZX_L[k]) + Pow2(AREAZY_L[k]));
                                    }

                                    //Recurrence formula
                                    help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                    PIM[kn] = BIM_L[kn] * help;
                                    QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    m = 2;
                                }
                            }

                            //Obtain new W-components
                            m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    W1N_L[k] += (relaxv * (PIM[kn] * W2N_L[k] + QIM[kn] - W1N_L[k]));
                                    m--;
                                }
                                else
                                {
                                    W2N_L[k] += (relaxv * (PIM[kn] * W1N_L[k + 1] + QIM[kn] - W2N_L[k]));
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
                // Parallel.For(2, NI, Program.pOptions, ih =>
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
                            float[] AN2_L = Program.AN2[i][j];
                            float[] AIM_L = Program.AIM[i][j];
                            ReadOnlySpan<float> AREAZX_L = Program.AREAZXImm[i][j].AsSpan();
                            ReadOnlySpan<float> AREAZY_L = Program.AREAZYImm[i][j].AsSpan();
                            float[] AS1_L = Program.AS1[i][j];
                            float[] AW1_L = Program.AW1[i][j];
                            float[] BIM_L = Program.BIM[i][j];
                            float[] CIM_L = Program.CIM[i][j];
                            float[] F1W_L = Program.F1W[i][j];
                            float[] F2W_L = Program.F2W[i][j];
                            ReadOnlySpan<float> RHO_L = Program.RHO[i][j];
                            double[] W1N_L = Program.W1N[i][j]; double[] W1Ni_L = Program.W1N[i + 1][j]; double[] W1NJ_P_L = Program.W1N[i][j + 1];
                            double[] W2N_L = Program.W2N[i][j]; double[] W2Ni_L = Program.W2N[i - 1][j]; double[] W2NJ_P_L = Program.W2N[i][j - 1];
                            float USTxUSTV = (float)(Program.UST[i][j] * Program.USTV[i][j]);

                            int m = 2;
                            for (int kn = 1; kn <= 2 * (NK_P - 1); kn++)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                                      //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = AW1_L[k] * W2Ni_L[k] + AS1_L[k] * W2NJ_P_L[k] + F1W_L[kn];
                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * W1N_L[k] * USTxUSTV * MathF.Sqrt(Pow2(AREAZX_L[k]) + Pow2(AREAZY_L[k]));
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
                                    DIM = AE2_L[k] * W1Ni_L[k] + AN2_L[k] * W1NJ_P_L[k] + F2W_L[kn];
                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * W2N_L[k] * USTxUSTV * MathF.Sqrt(Pow2(AREAZX_L[k]) + Pow2(AREAZY_L[k]));
                                    }

                                    //Recurrence formula
                                    help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                    PIM[kn] = BIM_L[kn] * help;
                                    QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    m = 2;
                                }
                            }

                            //Obtain new W-components
                            m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    W1N_L[k] += (relaxv * (PIM[kn] * W2N_L[k] + QIM[kn] - W1N_L[k]));
                                    m--;
                                }
                                else
                                {
                                    W2N_L[k] += (relaxv * (PIM[kn] * W1N_L[k + 1] + QIM[kn] - W2N_L[k]));
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
                            float[] AN2_L = Program.AN2[i][j];
                            float[] AIM_L = Program.AIM[i][j];
                            ReadOnlySpan<float> AREAZX_L = Program.AREAZXImm[i][j].AsSpan();
                            ReadOnlySpan<float> AREAZY_L = Program.AREAZYImm[i][j].AsSpan();
                            float[] AS1_L = Program.AS1[i][j];
                            float[] AW1_L = Program.AW1[i][j];
                            float[] BIM_L = Program.BIM[i][j];
                            float[] CIM_L = Program.CIM[i][j];
                            float[] F1W_L = Program.F1W[i][j];
                            float[] F2W_L = Program.F2W[i][j];
                            ReadOnlySpan<float> RHO_L = Program.RHO[i][j];
                            double[] W1N_L = Program.W1N[i][j]; double[] W1Ni_L = Program.W1N[i + 1][j]; double[] W1NJ_P_L = Program.W1N[i][j + 1];
                            double[] W2N_L = Program.W2N[i][j]; double[] W2Ni_L = Program.W2N[i - 1][j]; double[] W2NJ_P_L = Program.W2N[i][j - 1];
                            float USTxUSTV = (float)(Program.UST[i][j] * Program.USTV[i][j]);

                            int m = 2;
                            for (int kn = 1; kn <= 2 * (NK_P - 1); kn++)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                                      //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = AW1_L[k] * W2Ni_L[k] + AS1_L[k] * W2NJ_P_L[k] + F1W_L[kn];
                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * W1N_L[k] * USTxUSTV * MathF.Sqrt(Pow2(AREAZX_L[k]) + Pow2(AREAZY_L[k]));
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
                                    DIM = AE2_L[k] * W1Ni_L[k] + AN2_L[k] * W1NJ_P_L[k] + F2W_L[kn];
                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * W2N_L[k] * USTxUSTV * MathF.Sqrt(Pow2(AREAZX_L[k]) + Pow2(AREAZY_L[k]));
                                    }

                                    //Recurrence formula
                                    help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                    PIM[kn] = BIM_L[kn] * help;
                                    QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    m = 2;
                                }
                            }

                            //Obtain new W-components
                            m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    W1N_L[k] += (relaxv * (PIM[kn] * W2N_L[k] + QIM[kn] - W1N_L[k]));
                                    m--;
                                }
                                else
                                {
                                    W2N_L[k] += (relaxv * (PIM[kn] * W1N_L[k + 1] + QIM[kn] - W2N_L[k]));
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
                                                               // Parallel.For(2, NJ, Program.pOptions, j =>
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
                            float[] AN2_L = Program.AN2[i][j];
                            float[] AIM_L = Program.AIM[i][j];
                            ReadOnlySpan<float> AREAZX_L = Program.AREAZXImm[i][j].AsSpan();
                            ReadOnlySpan<float> AREAZY_L = Program.AREAZYImm[i][j].AsSpan();
                            float[] AS1_L = Program.AS1[i][j];
                            float[] AW1_L = Program.AW1[i][j];
                            float[] BIM_L = Program.BIM[i][j];
                            float[] CIM_L = Program.CIM[i][j];
                            float[] F1W_L = Program.F1W[i][j];
                            float[] F2W_L = Program.F2W[i][j];
                            ReadOnlySpan<float> RHO_L = Program.RHO[i][j];
                            double[] W1N_L = Program.W1N[i][j]; double[] W1Ni_L = Program.W1N[i + 1][j]; double[] W1NJ_P_L = Program.W1N[i][j + 1];
                            double[] W2N_L = Program.W2N[i][j]; double[] W2Ni_L = Program.W2N[i - 1][j]; double[] W2NJ_P_L = Program.W2N[i][j - 1];
                            float USTxUSTV = (float)(Program.UST[i][j] * Program.USTV[i][j]);

                            int m = 2;
                            for (int kn = 1; kn <= 2 * (NK_P - 1); kn++)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                                      //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = AW1_L[k] * W2Ni_L[k] + AS1_L[k] * W2NJ_P_L[k] + F1W_L[kn];
                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * W1N_L[k] * USTxUSTV * MathF.Sqrt(Pow2(AREAZX_L[k]) + Pow2(AREAZY_L[k]));
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
                                    DIM = AE2_L[k] * W1Ni_L[k] + AN2_L[k] * W1NJ_P_L[k] + F2W_L[kn];
                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * W2N_L[k] * USTxUSTV * MathF.Sqrt(Pow2(AREAZX_L[k]) + Pow2(AREAZY_L[k]));
                                    }

                                    //Recurrence formula
                                    help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                    PIM[kn] = BIM_L[kn] * help;
                                    QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    m = 2;
                                }
                            }

                            //Obtain new W-components
                            m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    W1N_L[k] += (relaxv * (PIM[kn] * W2N_L[k] + QIM[kn] - W1N_L[k]));
                                    m--;
                                }
                                else
                                {
                                    W2N_L[k] += (relaxv * (PIM[kn] * W1N_L[k + 1] + QIM[kn] - W2N_L[k]));
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
                            float[] AN2_L = Program.AN2[i][j];
                            float[] AIM_L = Program.AIM[i][j];
                            ReadOnlySpan<float> AREAZX_L = Program.AREAZXImm[i][j].AsSpan();
                            ReadOnlySpan<float> AREAZY_L = Program.AREAZYImm[i][j].AsSpan();
                            float[] AS1_L = Program.AS1[i][j];
                            float[] AW1_L = Program.AW1[i][j];
                            float[] BIM_L = Program.BIM[i][j];
                            float[] CIM_L = Program.CIM[i][j];
                            float[] F1W_L = Program.F1W[i][j];
                            float[] F2W_L = Program.F2W[i][j];
                            ReadOnlySpan<float> RHO_L = Program.RHO[i][j];
                            double[] W1N_L = Program.W1N[i][j]; double[] W1Ni_L = Program.W1N[i + 1][j]; double[] W1NJ_P_L = Program.W1N[i][j + 1];
                            double[] W2N_L = Program.W2N[i][j]; double[] W2Ni_L = Program.W2N[i - 1][j]; double[] W2NJ_P_L = Program.W2N[i][j - 1];
                            float USTxUSTV = (float)(Program.UST[i][j] * Program.USTV[i][j]);

                            int m = 2;
                            for (int kn = 1; kn <= 2 * (NK_P - 1); kn++)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                                      //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = AW1_L[k] * W2Ni_L[k] + AS1_L[k] * W2NJ_P_L[k] + F1W_L[kn];
                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * W1N_L[k] * USTxUSTV * MathF.Sqrt(Pow2(AREAZX_L[k]) + Pow2(AREAZY_L[k]));
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
                                    DIM = AE2_L[k] * W1Ni_L[k] + AN2_L[k] * W1NJ_P_L[k] + F2W_L[kn];
                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * W2N_L[k] * USTxUSTV * MathF.Sqrt(Pow2(AREAZX_L[k]) + Pow2(AREAZY_L[k]));
                                    }

                                    //Recurrence formula
                                    help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                    PIM[kn] = BIM_L[kn] * help;
                                    QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    m = 2;
                                }
                            }

                            //Obtain new W-components
                            m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    W1N_L[k] += (relaxv * (PIM[kn] * W2N_L[k] + QIM[kn] - W1N_L[k]));
                                    m--;
                                }
                                else
                                {
                                    W2N_L[k] += (relaxv * (PIM[kn] * W1N_L[k + 1] + QIM[kn] - W2N_L[k]));
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
                Parallel.For(2, NI, Program.pOptions, i =>
                {
                    int NK_P = NK; int NJ_P = NJ;
                    double DIM;
                    Span<double> PIM = stackalloc double[2 * NK_P];
                    Span<double> QIM = stackalloc double[2 * NK_P];
                    double help;

                    for (int j = NJ_P - 1; j >= 2; j--)
                    {
                        float relaxv = (float)(Program.RELAXV * Program.Relax_Border_factor[i][j]);
                        float[] AE2_L = Program.AE2[i][j];
                        float[] AN2_L = Program.AN2[i][j];
                        float[] AIM_L = Program.AIM[i][j];
                        ReadOnlySpan<float> AREAZX_L = Program.AREAZXImm[i][j].AsSpan();
                        ReadOnlySpan<float> AREAZY_L = Program.AREAZYImm[i][j].AsSpan();
                        float[] AS1_L = Program.AS1[i][j];
                        float[] AW1_L = Program.AW1[i][j];
                        float[] BIM_L = Program.BIM[i][j];
                        float[] CIM_L = Program.CIM[i][j];
                        float[] F1W_L = Program.F1W[i][j];
                        float[] F2W_L = Program.F2W[i][j];
                        ReadOnlySpan<float> RHO_L = Program.RHO[i][j];
                        double[] W1N_L = Program.W1N[i][j]; double[] W1Ni_L = Program.W1N[i + 1][j]; double[] W1NJ_P_L = Program.W1N[i][j + 1];
                        double[] W2N_L = Program.W2N[i][j]; double[] W2Ni_L = Program.W2N[i - 1][j]; double[] W2NJ_P_L = Program.W2N[i][j - 1];
                        float USTxUSTV = (float)(Program.UST[i][j] * Program.USTV[i][j]);

                        int m = 2;
                        for (int kn = 1; kn <= 2 * (NK_P - 1); kn++)
                        {
                            int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                            //Coefficients for the lower half-cell
                            if (m == 2)
                            {
                                DIM = AW1_L[k] * W2Ni_L[k] + AS1_L[k] * W2NJ_P_L[k] + F1W_L[kn];
                                if (k == 1)
                                {
                                    DIM -= RHO_L[k] * W1N_L[k] * USTxUSTV * MathF.Sqrt(Pow2(AREAZX_L[k]) + Pow2(AREAZY_L[k]));
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
                                DIM = AE2_L[k] * W1Ni_L[k] + AN2_L[k] * W1NJ_P_L[k] + F2W_L[kn];
                                if (k == 1)
                                {
                                    DIM -= RHO_L[k] * W2N_L[k] * USTxUSTV * MathF.Sqrt(Pow2(AREAZX_L[k]) + Pow2(AREAZY_L[k]));
                                }

                                //Recurrence formula
                                help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                PIM[kn] = BIM_L[kn] * help;
                                QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                m = 2;
                            }
                        }

                        //Obtain new W-components
                        m = 1;
                        for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                        {
                            int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                            if (m == 2)
                            {
                                W1N_L[k] += (relaxv * (PIM[kn] * W2N_L[k] + QIM[kn] - W1N_L[k]));
                                m--;
                            }
                            else
                            {
                                W2N_L[k] += (relaxv * (PIM[kn] * W1N_L[k + 1] + QIM[kn] - W2N_L[k]));
                                m = 2;
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
                            float[] AN2_L = Program.AN2[i][j];
                            float[] AIM_L = Program.AIM[i][j];
                            ReadOnlySpan<float> AREAZX_L = Program.AREAZXImm[i][j].AsSpan();
                            ReadOnlySpan<float> AREAZY_L = Program.AREAZYImm[i][j].AsSpan();
                            float[] AS1_L = Program.AS1[i][j];
                            float[] AW1_L = Program.AW1[i][j];
                            float[] BIM_L = Program.BIM[i][j];
                            float[] CIM_L = Program.CIM[i][j];
                            float[] F1W_L = Program.F1W[i][j];
                            float[] F2W_L = Program.F2W[i][j];
                            ReadOnlySpan<float> RHO_L = Program.RHO[i][j];
                            double[] W1N_L = Program.W1N[i][j]; double[] W1Ni_L = Program.W1N[i + 1][j]; double[] W1NJ_P_L = Program.W1N[i][j + 1];
                            double[] W2N_L = Program.W2N[i][j]; double[] W2Ni_L = Program.W2N[i - 1][j]; double[] W2NJ_P_L = Program.W2N[i][j - 1];
                            float USTxUSTV = (float)(Program.UST[i][j] * Program.USTV[i][j]);

                            int m = 2;
                            for (int kn = 1; kn <= 2 * (NK_P - 1); kn++)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                                      //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = AW1_L[k] * W2Ni_L[k] + AS1_L[k] * W2NJ_P_L[k] + F1W_L[kn];
                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * W1N_L[k] * USTxUSTV * MathF.Sqrt(Pow2(AREAZX_L[k]) + Pow2(AREAZY_L[k]));
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
                                    DIM = AE2_L[k] * W1Ni_L[k] + AN2_L[k] * W1NJ_P_L[k] + F2W_L[kn];
                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * W2N_L[k] * USTxUSTV * MathF.Sqrt(Pow2(AREAZX_L[k]) + Pow2(AREAZY_L[k]));
                                    }

                                    //Recurrence formula
                                    help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                    PIM[kn] = BIM_L[kn] * help;
                                    QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    m = 2;
                                }
                            }

                            //Obtain new W-components
                            m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    W1N_L[k] += (relaxv * (PIM[kn] * W2N_L[k] + QIM[kn] - W1N_L[k]));
                                    m--;
                                }
                                else
                                {
                                    W2N_L[k] += (relaxv * (PIM[kn] * W1N_L[k + 1] + QIM[kn] - W2N_L[k]));
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
                            float[] AN2_L = Program.AN2[i][j];
                            float[] AIM_L = Program.AIM[i][j];
                            ReadOnlySpan<float> AREAZX_L = Program.AREAZXImm[i][j].AsSpan();
                            ReadOnlySpan<float> AREAZY_L = Program.AREAZYImm[i][j].AsSpan();
                            float[] AS1_L = Program.AS1[i][j];
                            float[] AW1_L = Program.AW1[i][j];
                            float[] BIM_L = Program.BIM[i][j];
                            float[] CIM_L = Program.CIM[i][j];
                            float[] F1W_L = Program.F1W[i][j];
                            float[] F2W_L = Program.F2W[i][j];
                            ReadOnlySpan<float> RHO_L = Program.RHO[i][j];
                            double[] W1N_L = Program.W1N[i][j]; double[] W1Ni_L = Program.W1N[i + 1][j]; double[] W1NJ_P_L = Program.W1N[i][j + 1];
                            double[] W2N_L = Program.W2N[i][j]; double[] W2Ni_L = Program.W2N[i - 1][j]; double[] W2NJ_P_L = Program.W2N[i][j - 1];
                            float USTxUSTV = (float)(Program.UST[i][j] * Program.USTV[i][j]);

                            int m = 2;
                            for (int kn = 1; kn <= 2 * (NK_P - 1); kn++)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                                      //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = AW1_L[k] * W2Ni_L[k] + AS1_L[k] * W2NJ_P_L[k] + F1W_L[kn];
                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * W1N_L[k] * USTxUSTV * MathF.Sqrt(Pow2(AREAZX_L[k]) + Pow2(AREAZY_L[k]));
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
                                    DIM = AE2_L[k] * W1Ni_L[k] + AN2_L[k] * W1NJ_P_L[k] + F2W_L[kn];
                                    if (k == 1)
                                    {
                                        DIM -= RHO_L[k] * W2N_L[k] * USTxUSTV * MathF.Sqrt(Pow2(AREAZX_L[k]) + Pow2(AREAZY_L[k]));
                                    }

                                    //Recurrence formula
                                    help = 1 / (AIM_L[kn] - CIM_L[kn] * PIM[kn - 1]);
                                    PIM[kn] = BIM_L[kn] * help;
                                    QIM[kn] = (DIM + CIM_L[kn] * QIM[kn - 1]) * help;
                                    m = 2;
                                }
                            }

                            //Obtain new W-components
                            m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    W1N_L[k] += (relaxv * (PIM[kn] * W2N_L[k] + QIM[kn] - W1N_L[k]));
                                    m--;
                                }
                                else
                                {
                                    W2N_L[k] += (relaxv * (PIM[kn] * W1N_L[k + 1] + QIM[kn] - W2N_L[k]));
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
