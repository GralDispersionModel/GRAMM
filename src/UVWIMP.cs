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
using System.Threading;

namespace GRAMM_2001
{
    partial class Program
    {
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        public static void UVWIMPcalculate(int NI, int NJ, int NK)
        {
            if (Program.ICU)
            {
                //computation of Coriolis-term
                Parallel.For(2, NI, Program.pOptions, i =>
                {
                    int NK_P = NK; int NJ_P = NJ;
                    for (int j = 2; j <= NJ_P - 1; ++j)
                    {
                        float[] AP0_L = Program.AP0[i][j];
                        double[] U1_L = Program.U1[i][j];
                        double[] U2_L = Program.U2[i][j];
                        double[] V1_L = Program.V1[i][j];
                        double[] V2_L = Program.V2[i][j];
                        double[] W1_L = Program.W1[i][j];
                        double[] W2_L = Program.W2[i][j];
                        double[] DIMU_L = Program.DIMU[i][j];
                        double[] DIMV_L = Program.DIMV[i][j];
                        double[] DIMW_L = Program.DIMW[i][j];
                        
                        float[] F1U_L = Program.F1U[i][j];
                        float[] F2U_L = Program.F2U[i][j];
                        float[] F1V_L = Program.F1V[i][j];
                        float[] F2V_L = Program.F2V[i][j];
                        float[] F1W_L = Program.F1W[i][j];
                        float[] F2W_L = Program.F2W[i][j];
                        float[] RHO_L = Program.RHO[i][j];
                        float[] RHOBZ_L = Program.RHOBZ[i][j];
                        double[] QUN_L = Program.QUN[i][j];
                        double[] QBZ_L = Program.QBZ[i][j];
                        
                        double[] V1N_LR = Program.V1N[i][j];
                        double[] U1N_LR = Program.U1N[i][j];
                        float[] VOL_L = Program.VOL[i][j];
                        double[] UG_L = Program.UG[i][j];
                        double[] VG_L = Program.VG[i][j];
                        double[] W1N_LR = Program.W1N[i][j];
                        double[] TPDX_L = Program.TPDX[i][j];
                        double[] TPDX_M_L = Program.TPDX[i - 1][j];
                        double[] TPDY_L = Program.TPDY[i][j];
                        double[] TPDY_M_L = Program.TPDY[i][j - 1];
                        double[] TBZ_L = Program.TBZ[i][j];
                        double[] TN_L = Program.TN[i][j];

                        double f1; double f2 = 0.0;
                        for (int kn = 1; kn <= 2 * (NK_P - 1); ++kn)
                        {
                            int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                            f1 = RHO_L[k] * (Program.FN * (V1N_LR[k] - VG_L[k]) - Program.FH * W1N_LR[k]) * 0.5F * VOL_L[k];
                            f1 += DIMU_L[kn];
                            if (ICPSI)
                            {
                                f1 += TPDX_M_L[k];
                                f2 = TPDX_L[k];
                            }
                            F1U_L[kn] = (float)(AP0_L[k] * U1_L[k] + f1); // Temp array
                            F2U_L[kn] = (float)(AP0_L[k] * U2_L[k] + f1 + f2); // Temp array
                        }

                        for (int kn = 1; kn <= 2 * (NK_P - 1); ++kn)
                        {
                            int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                            f1 = -RHO_L[k] * Program.FN * (U1N_LR[k] - UG_L[k]) * 0.5F * VOL_L[k];
                            f1 += DIMV_L[kn];
                            if (ICPSI)
                            {
                                f1 += TPDY_M_L[k];
                                f2 = TPDY_L[k];
                            }
                            F1V_L[kn] = (float)(AP0_L[k] * V1_L[k] + f1); // Temp array
                            F2V_L[kn] = (float)(AP0_L[k] * V2_L[k] + f1 + f2); // Temp array
                        }
                        
                        for (int kn = 1; kn <= 2 * (NK_P - 1); ++kn)
                        {
                            int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)

                            f1 = RHO_L[k] * Program.FH * U1N_LR[k] * 0.5F * VOL_L[k];
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


                int range_parallel = (int)(NI / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
                range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
                range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
                //Iterative solution using an implicit scheme and the TDMA or Thomas-Algorithm
                //Parallel.For(2, NI, Program.pOptions, i =>
                Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
                {
                    int NK_P = NK; int NJ_P = NJ;
                    double[] U1N_LR = Program.GrammArrayPool.Rent(Program.U1N[1][1].Length);
                    double[] U2N_LR = Program.GrammArrayPool.Rent(Program.U2N[1][1].Length);
                    double[] V1N_LR = Program.GrammArrayPool.Rent(Program.V1N[1][1].Length);
                    double[] V2N_LR = Program.GrammArrayPool.Rent(Program.V2N[1][1].Length);
                    double[] W1N_LR = Program.GrammArrayPool.Rent(Program.W1N[1][1].Length);
                    double[] W2N_LR = Program.GrammArrayPool.Rent(Program.W2N[1][1].Length);

                    for (int i = range.Item1; i < range.Item2; ++i)
                    {
                        int border = Math.Min(i - range.Item1, range.Item2 - 1 - i);
                        for (int j = 2; j <= NJ_P - 1; ++j)
                        {
                            UIMPKernel(i, j, NK_P, border, U1N_LR, V1N_LR, W1N_LR, U2N_LR, V2N_LR, W2N_LR);
                        }
                    }
                    Program.GrammArrayPool.Return(U1N_LR);
                    Program.GrammArrayPool.Return(V1N_LR);
                    Program.GrammArrayPool.Return(W1N_LR);
                    Program.GrammArrayPool.Return(U2N_LR);
                    Program.GrammArrayPool.Return(V2N_LR);
                    Program.GrammArrayPool.Return(W2N_LR);
                });

                range_parallel = (int)(NI / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
                range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
                range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
                //Iterative solution using an implicit scheme and the TDMA or Thomas-Algorithm
                //Parallel.For(2, NI, Program.pOptions, ih =>
                Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
                {
                    int NK_P = NK; int NJ_P = NJ;
                    double[] U1N_LR = Program.GrammArrayPool.Rent(Program.U1N[1][1].Length);
                    double[] U2N_LR = Program.GrammArrayPool.Rent(Program.U2N[1][1].Length);
                    double[] V1N_LR = Program.GrammArrayPool.Rent(Program.V1N[1][1].Length);
                    double[] V2N_LR = Program.GrammArrayPool.Rent(Program.V2N[1][1].Length);
                    double[] W1N_LR = Program.GrammArrayPool.Rent(Program.W1N[1][1].Length);
                    double[] W2N_LR = Program.GrammArrayPool.Rent(Program.W2N[1][1].Length);

                    for (int i = range.Item2 - 1; i >= range.Item1; --i)
                    {
                        int border = Math.Min(i - range.Item1, range.Item2 - 1 - i);
                        for (int j = NJ_P - 1; j >= 2; --j)
                        {
                            UIMPKernel(i, j, NK_P, border, U1N_LR, V1N_LR, W1N_LR, U2N_LR, V2N_LR, W2N_LR);
                        }
                    }
                    Program.GrammArrayPool.Return(U1N_LR);
                    Program.GrammArrayPool.Return(V1N_LR);
                    Program.GrammArrayPool.Return(W1N_LR);
                    Program.GrammArrayPool.Return(U2N_LR);
                    Program.GrammArrayPool.Return(V2N_LR);
                    Program.GrammArrayPool.Return(W2N_LR);
                });

                range_parallel = (int)(NJ / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
                range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
                range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
                //Iterative solution using an implicit scheme and the TDMA or Thomas-Algorithm
                //Parallel.For(2, NJ, Program.pOptions, jh =>
                Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
                {
                    int NK_P = NK; int NI_P = NI;
                    double[] U1N_LR = Program.GrammArrayPool.Rent(Program.U1N[1][1].Length);
                    double[] U2N_LR = Program.GrammArrayPool.Rent(Program.U2N[1][1].Length);
                    double[] V1N_LR = Program.GrammArrayPool.Rent(Program.V1N[1][1].Length);
                    double[] V2N_LR = Program.GrammArrayPool.Rent(Program.V2N[1][1].Length);
                    double[] W1N_LR = Program.GrammArrayPool.Rent(Program.W1N[1][1].Length);
                    double[] W2N_LR = Program.GrammArrayPool.Rent(Program.W2N[1][1].Length);

                    for (int j = range.Item2 - 1; j >= range.Item1; --j)
                    {
                        int border = Math.Min(j - range.Item1, range.Item2 - 1 - j);
                        for (int i = NI_P - 1; i >= 2; --i)
                        {
                            UIMPKernel(i, j, NK_P, border, U1N_LR, V1N_LR, W1N_LR, U2N_LR, V2N_LR, W2N_LR);
                        }
                    }
                    Program.GrammArrayPool.Return(U1N_LR);
                    Program.GrammArrayPool.Return(V1N_LR);
                    Program.GrammArrayPool.Return(W1N_LR);
                    Program.GrammArrayPool.Return(U2N_LR);
                    Program.GrammArrayPool.Return(V2N_LR);
                    Program.GrammArrayPool.Return(W2N_LR);
                });

                range_parallel = (int)(NJ / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
                range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
                range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
                //Iterative solution using an implicit scheme and the TDMA or Thomas-Algorithm
                //Parallel.For(2, NJ, Program.pOptions, j =>
                Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
                {
                    int NK_P = NK; int NI_P = NI;
                    double[] U1N_LR = Program.GrammArrayPool.Rent(Program.U1N[1][1].Length);
                    double[] U2N_LR = Program.GrammArrayPool.Rent(Program.U2N[1][1].Length);
                    double[] V1N_LR = Program.GrammArrayPool.Rent(Program.V1N[1][1].Length);
                    double[] V2N_LR = Program.GrammArrayPool.Rent(Program.V2N[1][1].Length);
                    double[] W1N_LR = Program.GrammArrayPool.Rent(Program.W1N[1][1].Length);
                    double[] W2N_LR = Program.GrammArrayPool.Rent(Program.W2N[1][1].Length);

                    for (int j = range.Item1; j < range.Item2; ++j)
                    {
                        int border = Math.Min(j - range.Item1, range.Item2 - 1 - j);
                        for (int i = 2; i <= NI_P - 1; ++i)
                        {
                            UIMPKernel(i, j, NK_P, border, U1N_LR, V1N_LR, W1N_LR, U2N_LR, V2N_LR, W2N_LR);
                        }
                    }
                    Program.GrammArrayPool.Return(U1N_LR);
                    Program.GrammArrayPool.Return(V1N_LR);
                    Program.GrammArrayPool.Return(W1N_LR);
                    Program.GrammArrayPool.Return(U2N_LR);
                    Program.GrammArrayPool.Return(V2N_LR);
                    Program.GrammArrayPool.Return(W2N_LR);
                });

                range_parallel = (int)(NI / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
                range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
                range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
                //Iterative solution using an implicit scheme and the TDMA or Thomas-Algorithm
                //Parallel.For(2, NI, Program.pOptions, ih =>
                Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
                {
                    int NK_P = NK; int NJ_P = NJ;
                    double[] U1N_LR = Program.GrammArrayPool.Rent(Program.U1N[1][1].Length);
                    double[] U2N_LR = Program.GrammArrayPool.Rent(Program.U2N[1][1].Length);
                    double[] V1N_LR = Program.GrammArrayPool.Rent(Program.V1N[1][1].Length);
                    double[] V2N_LR = Program.GrammArrayPool.Rent(Program.V2N[1][1].Length);
                    double[] W1N_LR = Program.GrammArrayPool.Rent(Program.W1N[1][1].Length);
                    double[] W2N_LR = Program.GrammArrayPool.Rent(Program.W2N[1][1].Length);

                    for (int i = range.Item2 - 1; i >= range.Item1; --i)
                    {
                        int border = Math.Min(i - range.Item1, range.Item2 - 1 - i);
                        for (int j = 2; j <= NJ_P - 1; ++j)
                        {
                            UIMPKernel(i, j, NK_P, border, U1N_LR, V1N_LR, W1N_LR, U2N_LR, V2N_LR, W2N_LR);
                        }
                    }
                    Program.GrammArrayPool.Return(U1N_LR);
                    Program.GrammArrayPool.Return(V1N_LR);
                    Program.GrammArrayPool.Return(W1N_LR);
                    Program.GrammArrayPool.Return(U2N_LR);
                    Program.GrammArrayPool.Return(V2N_LR);
                    Program.GrammArrayPool.Return(W2N_LR);
                });

                range_parallel = (int)(NI / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
                range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
                range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
                //Iterative solution using an implicit scheme and the TDMA or Thomas-Algorithm
                //Parallel.For(2, NI, Program.pOptions, i =>
                Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
                {
                    int NK_P = NK; int NJ_P = NJ;
                    double[] U1N_LR = Program.GrammArrayPool.Rent(Program.U1N[1][1].Length);
                    double[] U2N_LR = Program.GrammArrayPool.Rent(Program.U2N[1][1].Length);
                    double[] V1N_LR = Program.GrammArrayPool.Rent(Program.V1N[1][1].Length);
                    double[] V2N_LR = Program.GrammArrayPool.Rent(Program.V2N[1][1].Length);
                    double[] W1N_LR = Program.GrammArrayPool.Rent(Program.W1N[1][1].Length);
                    double[] W2N_LR = Program.GrammArrayPool.Rent(Program.W2N[1][1].Length);

                    for (int i = range.Item1; i < range.Item2; ++i)
                    {
                        int border = Math.Min(i - range.Item1, range.Item2 - 1 - i);
                        for (int j = NJ_P - 1; j >= 2; --j)
                        {
                            UIMPKernel(i, j, NK_P, border, U1N_LR, V1N_LR, W1N_LR, U2N_LR, V2N_LR, W2N_LR);
                        }
                    }
                    Program.GrammArrayPool.Return(U1N_LR);
                    Program.GrammArrayPool.Return(V1N_LR);
                    Program.GrammArrayPool.Return(W1N_LR);
                    Program.GrammArrayPool.Return(U2N_LR);
                    Program.GrammArrayPool.Return(V2N_LR);
                    Program.GrammArrayPool.Return(W2N_LR);
                });

                range_parallel = (int)(NJ / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
                range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
                range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
                //Iterative solution using an implicit scheme and the TDMA or Thomas-Algorithm
                //Parallel.For(2, NJ, Program.pOptions, jh =>
                Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
                {
                    int NK_P = NK; int NI_P = NI;
                    double[] U1N_LR = Program.GrammArrayPool.Rent(Program.U1N[1][1].Length);
                    double[] U2N_LR = Program.GrammArrayPool.Rent(Program.U2N[1][1].Length);
                    double[] V1N_LR = Program.GrammArrayPool.Rent(Program.V1N[1][1].Length);
                    double[] V2N_LR = Program.GrammArrayPool.Rent(Program.V2N[1][1].Length);
                    double[] W1N_LR = Program.GrammArrayPool.Rent(Program.W1N[1][1].Length);
                    double[] W2N_LR = Program.GrammArrayPool.Rent(Program.W2N[1][1].Length);

                    for (int j = range.Item2 - 1; j >= range.Item1; --j)
                    {
                        int border = Math.Min(j - range.Item1, range.Item2 - 1 - j);
                        for (int i = 2; i <= NI_P - 1; ++i)
                        {
                            UIMPKernel(i, j, NK_P, border, U1N_LR, V1N_LR, W1N_LR, U2N_LR, V2N_LR, W2N_LR);
                        }
                    }
                    Program.GrammArrayPool.Return(U1N_LR);
                    Program.GrammArrayPool.Return(V1N_LR);
                    Program.GrammArrayPool.Return(W1N_LR);
                    Program.GrammArrayPool.Return(U2N_LR);
                    Program.GrammArrayPool.Return(V2N_LR);
                    Program.GrammArrayPool.Return(W2N_LR);
                });

                range_parallel = (int)(NJ / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
                range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
                range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
                //Iterative solution using an implicit scheme and the TDMA or Thomas-Algorithm
                //Parallel.For(2, NJ, Program.pOptions, j =>
                Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
                {
                    int NK_P = NK; int NI_P = NI;
                    double[] U1N_LR = Program.GrammArrayPool.Rent(Program.U1N[1][1].Length);
                    double[] U2N_LR = Program.GrammArrayPool.Rent(Program.U2N[1][1].Length);
                    double[] V1N_LR = Program.GrammArrayPool.Rent(Program.V1N[1][1].Length);
                    double[] V2N_LR = Program.GrammArrayPool.Rent(Program.V2N[1][1].Length);
                    double[] W1N_LR = Program.GrammArrayPool.Rent(Program.W1N[1][1].Length);
                    double[] W2N_LR = Program.GrammArrayPool.Rent(Program.W2N[1][1].Length);

                    for (int j = range.Item1; j < range.Item2; ++j)
                    {
                        int border = Math.Min(j - range.Item1, range.Item2 - 1 - j);
                        for (int i = 2; i <= NI_P - 1; ++i)
                        {
                            UIMPKernel(i, j, NK_P, border, U1N_LR, V1N_LR, W1N_LR, U2N_LR, V2N_LR, W2N_LR);
                        }
                    }
                    Program.GrammArrayPool.Return(U1N_LR);
                    Program.GrammArrayPool.Return(V1N_LR);
                    Program.GrammArrayPool.Return(W1N_LR);
                    Program.GrammArrayPool.Return(U2N_LR);
                    Program.GrammArrayPool.Return(V2N_LR);
                    Program.GrammArrayPool.Return(W2N_LR);
                });
            }
        }
         
        [MethodImpl(MethodImplOptions.AggressiveOptimization | MethodImplOptions.AggressiveInlining)]
        private static void UIMPKernel(int i, int j, int NK_P, int border, double[] U1N_LRR, double[] V1N_LRR, double[] W1N_LRR, double[] U2N_LRR, double[] V2N_LRR, double[] W2N_LRR)
        {
            Unsafe.SkipInit(out double DIMU); 
            Unsafe.SkipInit(out double DIMV); 
            Unsafe.SkipInit(out double DIMW); 
            Unsafe.SkipInit(out double help);
            Unsafe.SkipInit(out double as1);
            Unsafe.SkipInit(out double aw1);
            Unsafe.SkipInit(out double ae2);
            Unsafe.SkipInit(out double an2);
            int knmax = 2 * (NK_P - 1) + 1;
            Span<double> PIMU = stackalloc double[knmax];
            Span<double> QIMU = stackalloc double[knmax];
            Span<double> PIMV = stackalloc double[knmax];
            Span<double> QIMV = stackalloc double[knmax];
            Span<double> PIMW = stackalloc double[knmax];
            Span<double> QIMW = stackalloc double[knmax];
            
            float relaxv = (float)(Program.RELAXV * Program.Relax_Border_factor[i][j]);
            ReadOnlySpan<float> AE2_L = Program.AE2[i][j];
            ReadOnlySpan<float> AN2_L = Program.AN2[i][j];
            ReadOnlySpan<float> AIM_L = Program.AIM[i][j];
            ReadOnlySpan<float> AREA_L = Program.AREA[i][j];
            ReadOnlySpan<float> AS1_L = Program.AS1[i][j];
            ReadOnlySpan<float> AW1_L = Program.AW1[i][j];
            ReadOnlySpan<float> BIM_L = Program.BIM[i][j];
            ReadOnlySpan<float> CIM_L = Program.CIM[i][j];
            ReadOnlySpan<float> F1U_L = Program.F1U[i][j];
            ReadOnlySpan<float> F2U_L = Program.F2U[i][j];
            ReadOnlySpan<float> F1V_L = Program.F1V[i][j];
            ReadOnlySpan<float> F2V_L = Program.F2V[i][j];
            ReadOnlySpan<float> F1W_L = Program.F1W[i][j];
            ReadOnlySpan<float> F2W_L = Program.F2W[i][j];
            ReadOnlySpan<float> RHO_L = Program.RHO[i][j];
            ReadOnlySpan <double> U1Ni_L  = Program.U1N[i + 1][j];
            ReadOnlySpan <double> U1NJ_P_L = Program.U1N[i][j + 1];
            ReadOnlySpan <double> U2Ni_L = Program.U2N[i - 1][j];
            ReadOnlySpan <double> U2NJ_P_L = Program.U2N[i][j - 1];
            ReadOnlySpan <double> V1Ni_L = Program.V1N[i + 1][j];
            ReadOnlySpan <double> V1NJ_P_L = Program.V1N[i][j + 1];
            ReadOnlySpan <double> V2Ni_L = Program.V2N[i - 1][j];
            ReadOnlySpan <double> V2NJ_P_L = Program.V2N[i][j - 1];
            ReadOnlySpan <double> W1Ni_L = Program.W1N[i + 1][j];
            ReadOnlySpan <double> W1NJ_P_L = Program.W1N[i][j + 1];
            ReadOnlySpan <double> W2Ni_L = Program.W2N[i - 1][j];
            ReadOnlySpan <double> W2NJ_P_L = Program.W2N[i][j - 1];

            double[] U1N_LR;
            double[] V1N_LR;
            double[] W1N_LR;
            double[] U2N_LR;
            double[] V2N_LR;
            double[] W2N_LR;

            //Avoid race conditions at the border cells of the sequential calculated stripes
            if (border < 2)
            {
                U1N_LR = U1N_LRR;
                V1N_LR = V1N_LRR;
                W1N_LR = W1N_LRR;
                U2N_LR = U2N_LRR;
                V2N_LR = V2N_LRR;
                W2N_LR = W2N_LRR;
                Program.CopyArraySourceLen(Program.U1N[i][j], U1N_LR);
                Program.CopyArraySourceLen(Program.V1N[i][j], V1N_LR);
                Program.CopyArraySourceLen(Program.W1N[i][j], W1N_LR);
                Program.CopyArraySourceLen(Program.U2N[i][j], U2N_LR);
                Program.CopyArraySourceLen(Program.V2N[i][j], V2N_LR);
                Program.CopyArraySourceLen(Program.W2N[i][j], W2N_LR);
            }
            else
            {
                U1N_LR = Program.U1N[i][j];
                V1N_LR = Program.V1N[i][j];
                W1N_LR = Program.W1N[i][j];
                U2N_LR = Program.U2N[i][j];
                V2N_LR = Program.V2N[i][j];
                W2N_LR = Program.W2N[i][j];
            }

            float USTxUSTV = (float)(Program.UST[i][j] * Program.USTV[i][j]) * RHO_L[1] * AREA_L[1];
            Unsafe.SkipInit(out double aim);
            Unsafe.SkipInit(out double bim);
            Unsafe.SkipInit(out double cim);

            int m = 2;
            for (int kn = 1; kn < PIMU.Length; ++kn)
            {
                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                //Coefficients for the lower half-cell
                if (m == 2)
                {
                    aw1 = AW1_L[k];
                    as1 = AS1_L[k];

                    // U - Component
                    DIMU = aw1 * U2Ni_L[k] + as1 * U2NJ_P_L[k] + F1U_L[kn];
                    DIMV = aw1 * V2Ni_L[k] + as1 * V2NJ_P_L[k] + F1V_L[kn];
                    DIMW = aw1 * W2Ni_L[k] + as1 * W2NJ_P_L[k] + F1W_L[kn];

                    aim = AIM_L[kn];
                    bim = BIM_L[kn];
                    //Recurrence formula
                    if (k == 1)
                    {
                        help = 1 / aim;
                        DIMU -= U1N_LR[k] * USTxUSTV;
                        PIMU[kn] = bim * help;
                        QIMU[kn] = DIMU * help;

                        DIMV -= V1N_LR[k] * USTxUSTV;
                        PIMV[kn] = bim * help;
                        QIMV[kn] = DIMV * help;

                        DIMW -= W1N_LR[k] * USTxUSTV;
                        PIMW[kn] = bim * help;
                        QIMW[kn] = DIMW * help;
                    }
                    else
                    {
                        cim = CIM_L[kn];
                        help = 1 / (aim - cim * PIMU[kn - 1]);
                        PIMU[kn] = bim * help;
                        QIMU[kn] = (DIMU + cim * QIMU[kn - 1]) * help;

                        help = 1 / (aim - cim * PIMV[kn - 1]);
                        PIMV[kn] = bim * help;
                        QIMV[kn] = (DIMV + cim * QIMV[kn - 1]) * help;

                        help = 1 / (aim - cim * PIMW[kn - 1]);
                        PIMW[kn] = bim * help;
                        QIMW[kn] = (DIMW + cim * QIMW[kn - 1]) * help;
                    }
                    m--;
                }
                else
                {
                    ae2 = AE2_L[k];
                    an2 = AN2_L[k];

                    //Coefficients for the upper half-cell
                    // U - Component
                    DIMU = ae2 * U1Ni_L[k] + an2 * U1NJ_P_L[k] + F2U_L[kn];
                    DIMV = ae2 * V1Ni_L[k] + an2 * V1NJ_P_L[k] + F2V_L[kn];
                    DIMW = ae2 * W1Ni_L[k] + an2 * W1NJ_P_L[k] + F2W_L[kn];

                    aim = AIM_L[kn];
                    bim = BIM_L[kn];
                    if (k == 1)
                    {
                        DIMU -= U2N_LR[k] * USTxUSTV;
                        DIMV -= V2N_LR[k] * USTxUSTV;
                        DIMW -= W2N_LR[k] * USTxUSTV;
                    }
                    //Recurrence formula
                    cim = CIM_L[kn];
                    help = 1 / (aim - cim * PIMU[kn - 1]);
                    PIMU[kn] = bim * help;
                    QIMU[kn] = (DIMU + cim * QIMU[kn - 1]) * help;

                    help = 1 / (aim - cim * PIMV[kn - 1]);
                    PIMV[kn] = bim * help;
                    QIMV[kn] = (DIMV + cim * QIMV[kn - 1]) * help;

                    help = 1 / (aim - cim * PIMW[kn - 1]);
                    PIMW[kn] = bim * help;
                    QIMW[kn] = (DIMW + cim * QIMW[kn - 1]) * help;
                    m = 2;
                }
            }

            //Obtain new UVW-components
            for (int kn = PIMU.Length - 1; kn > 1; --kn)
            {
                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                U2N_LR[k] += (relaxv * (PIMU[kn] * U1N_LR[k + 1] + QIMU[kn] - U2N_LR[k]));
                V2N_LR[k] += (relaxv * (PIMV[kn] * V1N_LR[k + 1] + QIMV[kn] - V2N_LR[k]));
                W2N_LR[k] += (relaxv * (PIMW[kn] * W1N_LR[k + 1] + QIMW[kn] - W2N_LR[k]));

                --kn;
                k = 1 + kn >> 1;
                U1N_LR[k] += (relaxv * (PIMU[kn] * U2N_LR[k] + QIMU[kn] - U1N_LR[k]));
                V1N_LR[k] += (relaxv * (PIMV[kn] * V2N_LR[k] + QIMV[kn] - V1N_LR[k]));
                W1N_LR[k] += (relaxv * (PIMW[kn] * W2N_LR[k] + QIMW[kn] - W1N_LR[k]));
            }
            //Avoid race conditions at the border cells of the sequential calculated stripes
            if (border < 2)
            {
                Program.CopyArrayLockDest(U1N_LR, Program.U1N[i][j]);
                Program.CopyArrayLockDest(V1N_LR, Program.V1N[i][j]);
                Program.CopyArrayLockDest(W1N_LR, Program.W1N[i][j]);
                Program.CopyArrayLockDest(U2N_LR, Program.U2N[i][j]);
                Program.CopyArrayLockDest(V2N_LR, Program.V2N[i][j]);
                Program.CopyArrayLockDest(W2N_LR, Program.W2N[i][j]);
            }
        }  
    }
}
