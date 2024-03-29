﻿#region Copyright
///<remarks>
/// <GRAMM Mesoscale Model>
/// Copyright (C) [2021]  [Markus Kuntner]
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
        /// Calculate U, V and W wind components
        /// </summary>
        /// <param name="NI"></param>
        /// <param name="NJ"></param>
        /// <param name="NK"></param>
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
                        ReadOnlySpan<float> RHO_L = Program.RHO[i][j];
                        float[] RHOBZ_L = Program.RHOBZ[i][j];
                        double[] QUN_L = Program.QUN[i][j];
                        double[] QBZ_L = Program.QBZ[i][j];

                        double[] V1N_LR = Program.V1N[i][j];
                        double[] U1N_LR = Program.U1N[i][j];
                        ReadOnlySpan<float> VOL_L = Program.VOLImm[i][j].AsSpan();
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
                            float rho = RHO_L[k];
                            float vol = VOL_L[k];
                            float ap0 = AP0_L[k];

                            f1 = rho * (Program.FN * (V1N_LR[k] - VG_L[k]) - Program.FH * W1N_LR[k]) * 0.5F * vol;
                            f1 += DIMU_L[kn];
                            f2 = 0;
                            if (ICPSI)
                            {
                                f1 += TPDX_M_L[k]; 
                                f2 = TPDX_L[k];
                            }
                            F1U_L[kn] = (float)(ap0 * U1_L[k] + f1); // Temp array
                            F2U_L[kn] = (float)(ap0 * U2_L[k] + f1 + f2); // Temp array
                        
                            f1 = -rho * Program.FN * (U1N_LR[k] - UG_L[k]) * 0.5F * vol;
                            f1 += DIMV_L[kn];
                            f2 = 0;
                            if (ICPSI)
                            {
                                f1 += TPDY_M_L[k];
                                f2 = TPDY_L[k];
                            }
                            F1V_L[kn] = (float)(ap0 * V1_L[k] + f1); // Temp array
                            F2V_L[kn] = (float)(ap0 * V2_L[k] + f1 + f2); // Temp array
                        
                            f1 = rho * Program.FH * U1N_LR[k] * 0.5F * vol;
                            f2 = 0;
                            //Bouyancy
                            if (Program.ICPSI == false)
                            {
                                f2 = 0.5F * vol * ((TN_L[k] + Program.TBZ1) * (1 + 0.00061 * QUN_L[k]) - (TBZ_L[k] * (1 + 0.00061 * QBZ_L[k]))) /
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
                            F1W_L[kn] = (float)(ap0 * W1_L[k] + f1);
                            F2W_L[kn] = (float)(ap0 * W2_L[k] + f1);
                        }
                    }
                });


                //Iterative solution using an implicit scheme and the TDMA or Thomas-Algorithm
                //Parallel.For(2, NI, Program.pOptions, i =>
                Parallel.ForEach(PartitionerI[Program.StripeCounter++ % PartitionerI.Count], range =>
                {
                    int NK_P = NK; int NJ_P = NJ;
                    UIMPKernel(range, NJ_P, NK_P, false, false, true);
                });
                //Parallel.For(2, NI, Program.pOptions, ih =>
                Parallel.ForEach(PartitionerI[Program.StripeCounter++ % PartitionerI.Count], range =>
                {
                    int NK_P = NK; int NJ_P = NJ;
                    UIMPKernel(range, NJ_P, NK_P, true, true, true);
                });
                //Parallel.For(2, NJ, Program.pOptions, jh =>
                Parallel.ForEach(PartitionerJ[Program.StripeCounter++ % PartitionerJ.Count], range =>
                {
                    int NK_P = NK; int NI_P = NI;
                    UIMPKernel(range, NI_P, NK_P, true, true, false);
                });
                //Parallel.For(2, NJ, Program.pOptions, j =>
                Parallel.ForEach(PartitionerJ[Program.StripeCounter++ % PartitionerJ.Count], range =>
                {
                    int NK_P = NK; int NI_P = NI;
                    UIMPKernel(range, NI_P, NK_P, false, false, false);
                });
                //Parallel.For(2, NI, Program.pOptions, ih =>
                Parallel.ForEach(PartitionerI[Program.StripeCounter++ % PartitionerI.Count], range =>
                {
                    int NK_P = NK; int NJ_P = NJ;
                    UIMPKernel(range, NJ_P, NK_P, true, false, true);
                });
                //Parallel.For(2, NI, Program.pOptions, i =>
                Parallel.ForEach(PartitionerI[Program.StripeCounter++ % PartitionerI.Count], range =>
                {
                    int NK_P = NK; int NJ_P = NJ;
                    UIMPKernel(range, NJ_P, NK_P, false, true, true);
                });
                //Parallel.For(2, NJ, Program.pOptions, jh =>
                Parallel.ForEach(PartitionerJ[Program.StripeCounter++ % PartitionerJ.Count], range =>
                {
                    int NK_P = NK; int NI_P = NI;
                    UIMPKernel(range, NI_P, NK_P, true, false, false);
                });
                //Parallel.For(2, NJ, Program.pOptions, j =>
                Parallel.ForEach(PartitionerJ[Program.StripeCounter++ % PartitionerJ.Count], range =>
                {
                    int NK_P = NK; int NI_P = NI;
                    UIMPKernel(range, NI_P, NK_P, false, false, false);
                });
            }
        }

        /// <summary>
        /// Iterative solution using an implicit scheme and the TDMA or Thomas-Algorithm
        /// </summary>
        /// <param name="range">Range for the outer loop</param>
        /// <param name="InnerLoopNmax">Maximum for the inner loop</param>
        /// <param name="outerLoopDir">true=--, false=++</param>
        /// <param name="innerLoopDir">true=--, false=++</param>
        /// <param name="isOuterLoopI">true: outer loop = i, false: outer loop = j</param>
        [MethodImpl(MethodImplOptions.AggressiveOptimization | MethodImplOptions.AggressiveInlining)]
        [SkipLocalsInit]
        private static void UIMPKernel(System.Tuple<int, int> range, int InnerLoopNmax, int NK_P, bool outerLoopDir, bool innerLoopDir, bool isOuterLoopI)
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
            double[] U1N_LRR = Program.GrammArrayPool.Rent(Program.U1N[1][1].Length);
            double[] U2N_LRR = Program.GrammArrayPool.Rent(Program.U2N[1][1].Length);
            double[] V1N_LRR = Program.GrammArrayPool.Rent(Program.V1N[1][1].Length);
            double[] V2N_LRR = Program.GrammArrayPool.Rent(Program.V2N[1][1].Length);
            double[] W1N_LRR = Program.GrammArrayPool.Rent(Program.W1N[1][1].Length);
            double[] W2N_LRR = Program.GrammArrayPool.Rent(Program.W2N[1][1].Length);
            ReadOnlySpan<float> AE2_L ; 
            ReadOnlySpan<float> AN2_L ; 
            ReadOnlySpan<float> AIM_L ; 
            ReadOnlySpan<float> AREA_L;
            ReadOnlySpan<float> AREAZX_L;
            ReadOnlySpan<float> AREAZY_L;
            ReadOnlySpan<float> AS1_L;
            ReadOnlySpan<float> AW1_L;
            ReadOnlySpan<float> BIM_L;
            ReadOnlySpan<float> CIM_L;
            ReadOnlySpan<float> F1U_L;
            ReadOnlySpan<float> F2U_L;
            ReadOnlySpan<float> F1V_L;
            ReadOnlySpan<float> F2V_L;
            ReadOnlySpan<float> F1W_L;
            ReadOnlySpan<float> F2W_L;
            ReadOnlySpan<float> RHO_L;
            double[] U1Ni_L;
            double[] U1NJ_P_L;
            double[] U2Ni_L;
            double[] U2NJ_P_L;
            double[] V1Ni_L;
            double[] V1NJ_P_L;
            double[] V2Ni_L;
            double[] V2NJ_P_L;
            double[] W1Ni_L;
            double[] W1NJ_P_L;
            double[] W2Ni_L;
            double[] W2NJ_P_L;
            
            for (int outerLoop = range.Item1; outerLoop < range.Item2; ++outerLoop)
            {
                int outerL = outerLoop;
                int border = Math.Min(outerL - range.Item1, range.Item2 - 1 - outerL);
                if (outerLoopDir)
                {
                    outerL = range.Item2 - (outerLoop - range.Item1) - 1;
                }
                
                for (int innerLoop = 2; innerLoop < InnerLoopNmax; ++innerLoop)
                {
                    int innerL = innerLoop;
                    if (innerLoopDir)
                    {
                        innerL = InnerLoopNmax - innerLoop + 1;
                    }
                    int i = 0;
                    int j = 0;
                    if (isOuterLoopI)
                    {
                        i = outerL;
                        j = innerL;
                    }
                    else
                    {
                        i = innerL;
                        j = outerL;
                    }

                    float relaxv = (float)(Program.RELAXV * Program.Relax_Border_factor[i][j]);
                    AE2_L = Program.AE2[i][j];
                    AN2_L = Program.AN2[i][j];
                    AIM_L = Program.AIM[i][j];
                    AREA_L = Program.AREAImm[i][j].AsSpan();
                    AREAZX_L = Program.AREAZXImm[i][j].AsSpan();
                    AREAZY_L = Program.AREAZYImm[i][j].AsSpan();
                    AS1_L = Program.AS1[i][j];
                    AW1_L = Program.AW1[i][j];
                    BIM_L = Program.BIM[i][j];
                    CIM_L = Program.CIM[i][j];
                    F1U_L = Program.F1U[i][j];
                    F2U_L = Program.F2U[i][j];
                    F1V_L = Program.F1V[i][j];
                    F2V_L = Program.F2V[i][j];
                    F1W_L = Program.F1W[i][j];
                    F2W_L = Program.F2W[i][j];
                    RHO_L = Program.RHO[i][j];

                    U1Ni_L   = Program.U1N[i + 1][j];
                    U1NJ_P_L = Program.U1N[i][j + 1];
                    U2Ni_L   = Program.U2N[i - 1][j];
                    U2NJ_P_L = Program.U2N[i][j - 1];
                    V1Ni_L   = Program.V1N[i + 1][j];
                    V1NJ_P_L = Program.V1N[i][j + 1];
                    V2Ni_L   = Program.V2N[i - 1][j];
                    V2NJ_P_L = Program.V2N[i][j - 1];
                    W1Ni_L   = Program.W1N[i + 1][j];
                    W1NJ_P_L = Program.W1N[i][j + 1];
                    W2Ni_L   = Program.W2N[i - 1][j];
                    W2NJ_P_L = Program.W2N[i][j - 1];

                    double[] U1N_LR;
                    double[] V1N_LR;
                    double[] W1N_LR;
                    double[] U2N_LR;
                    double[] V2N_LR;
                    double[] W2N_LR;

                    //Avoid race conditions at the border cells of the sequential calculated stripes
                    if (border < 1)
                    {
                        U1N_LR = U1N_LRR;
                        V1N_LR = V1N_LRR;
                        W1N_LR = W1N_LRR;
                        U2N_LR = U2N_LRR;
                        V2N_LR = V2N_LRR;
                        W2N_LR = W2N_LRR;
                        //No lock here, because other threads are reading
                        FastCopy.CopyArraySourceLen(Program.U1N[i][j], U1N_LR);
                        FastCopy.CopyArraySourceLen(Program.V1N[i][j], V1N_LR);
                        FastCopy.CopyArraySourceLen(Program.W1N[i][j], W1N_LR);
                        FastCopy.CopyArraySourceLen(Program.U2N[i][j], U2N_LR);
                        FastCopy.CopyArraySourceLen(Program.V2N[i][j], V2N_LR);
                        FastCopy.CopyArraySourceLen(Program.W2N[i][j], W2N_LR);
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

                    double USTxUSTVW = Program.UST[i][j] * Program.USTV[i][j] * RHO_L[1];
                    double USTxUSTV  = USTxUSTVW * AREA_L[1];                 
                    Unsafe.SkipInit(out double aim);
                    Unsafe.SkipInit(out double bim);
                    Unsafe.SkipInit(out double cim);

                    int m = 2;
                    for (int kn = 1; kn < PIMU.Length; ++kn)
                    {
                        int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                        int knM1 = kn - 1;
                        //Coefficients for the lower half-cell
                        if (m == 2)
                        {
                            aw1 = AW1_L[k];
                            as1 = AS1_L[k];

                            // UVW - Component
                            if (border < 1)
                            {
                                lock (U2Ni_L.SyncRoot)
                                {
                                    lock (U2NJ_P_L.SyncRoot)
                                    {
                                        DIMU = aw1 * U2Ni_L[k] + as1 * U2NJ_P_L[k] + F1U_L[kn];
                                    }
                                }
                                lock (V2Ni_L.SyncRoot)
                                {
                                    lock (V2NJ_P_L.SyncRoot)
                                    {
                                        DIMV = aw1 * V2Ni_L[k] + as1 * V2NJ_P_L[k] + F1V_L[kn];
                                    }
                                }
                                lock (W2Ni_L.SyncRoot)
                                {
                                    lock (W2NJ_P_L.SyncRoot)
                                    {
                                        DIMW = aw1 * W2Ni_L[k] + as1 * W2NJ_P_L[k] + F1W_L[kn];
                                    }
                                }
                            }
                            else
                            {
                                DIMU = aw1 * U2Ni_L[k] + as1 * U2NJ_P_L[k] + F1U_L[kn];
                                DIMV = aw1 * V2Ni_L[k] + as1 * V2NJ_P_L[k] + F1V_L[kn];
                                DIMW = aw1 * W2Ni_L[k] + as1 * W2NJ_P_L[k] + F1W_L[kn];
                            }
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

                                DIMW -= W1N_LR[k] * USTxUSTVW * MathF.Sqrt(Pow2(AREAZX_L[k]) + Pow2(AREAZY_L[k]));
                                PIMW[kn] = bim * help;
                                QIMW[kn] = DIMW * help;
                            }
                            else
                            {
                                cim = CIM_L[kn];
                                help = 1 / (aim - cim * PIMU[knM1]);
                                PIMU[kn] = bim * help;
                                QIMU[kn] = (DIMU + cim * QIMU[knM1]) * help;

                                help = 1 / (aim - cim * PIMV[knM1]);
                                PIMV[kn] = bim * help;
                                QIMV[kn] = (DIMV + cim * QIMV[knM1]) * help;

                                help = 1 / (aim - cim * PIMW[knM1]);
                                PIMW[kn] = bim * help;
                                QIMW[kn] = (DIMW + cim * QIMW[knM1]) * help;
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
                                DIMW -= W2N_LR[k] * USTxUSTVW * MathF.Sqrt(Pow2(AREAZX_L[k]) + Pow2(AREAZY_L[k]));
                                
                            }
                            //Recurrence formula
                            cim = CIM_L[kn];
                            help = 1 / (aim - cim * PIMU[knM1]);
                            PIMU[kn] = bim * help;
                            QIMU[kn] = (DIMU + cim * QIMU[knM1]) * help;

                            help = 1 / (aim - cim * PIMV[knM1]);
                            PIMV[kn] = bim * help;
                            QIMV[kn] = (DIMV + cim * QIMV[knM1]) * help;

                            help = 1 / (aim - cim * PIMW[knM1]);
                            PIMW[kn] = bim * help;
                            QIMW[kn] = (DIMW + cim * QIMW[knM1]) * help;
                            m = 2;
                        }
                    }

                    //Obtain new UVW-components
                    for (int kn = PIMU.Length - 1; kn > 1; --kn)
                    {
                        int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)
                        int kP1 = k + 1;
                        U2N_LR[k] += (relaxv * (PIMU[kn] * U1N_LR[kP1] + QIMU[kn] - U2N_LR[k]));
                        V2N_LR[k] += (relaxv * (PIMV[kn] * V1N_LR[kP1] + QIMV[kn] - V2N_LR[k]));
                        W2N_LR[k] += (relaxv * (PIMW[kn] * W1N_LR[kP1] + QIMW[kn] - W2N_LR[k]));

                        --kn;
                        k = 1 + kn >> 1;
                        U1N_LR[k] += (relaxv * (PIMU[kn] * U2N_LR[k] + QIMU[kn] - U1N_LR[k]));
                        V1N_LR[k] += (relaxv * (PIMV[kn] * V2N_LR[k] + QIMV[kn] - V1N_LR[k]));
                        W1N_LR[k] += (relaxv * (PIMW[kn] * W2N_LR[k] + QIMW[kn] - W1N_LR[k]));
                    }
                    //Avoid race conditions at the border cells of the sequential calculated stripes
                    if (border < 1)
                    {
                        FastCopy.CopyArrayLockDest(U1N_LR, Program.U1N[i][j]);
                        FastCopy.CopyArrayLockDest(V1N_LR, Program.V1N[i][j]);
                        FastCopy.CopyArrayLockDest(W1N_LR, Program.W1N[i][j]);
                        FastCopy.CopyArrayLockDest(U2N_LR, Program.U2N[i][j]);
                        FastCopy.CopyArrayLockDest(V2N_LR, Program.V2N[i][j]);
                        FastCopy.CopyArrayLockDest(W2N_LR, Program.W2N[i][j]);
                    }
                }
            }
            Program.GrammArrayPool.Return(U1N_LRR);
            Program.GrammArrayPool.Return(V1N_LRR);
            Program.GrammArrayPool.Return(W1N_LRR);
            Program.GrammArrayPool.Return(U2N_LRR);
            Program.GrammArrayPool.Return(V2N_LRR);
            Program.GrammArrayPool.Return(W2N_LRR);
        }
    }
}
