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
using System.Threading.Tasks;
using System.Threading;
using System.Runtime.CompilerServices;

namespace GRAMM_2001
{
    partial class Program
    {
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        public static void Timp_calculate(int NI, int NJ, int NK)
        {
            //computation of radiation terms
            Parallel.For(2, NI, Program.pOptions, i =>
            {
                for (int j = 2; j <= NJ - 1; ++j)
                {
                    ReadOnlySpan<double> DT_SOL_L = Program.DT_SOL[i][j];
                    ReadOnlySpan<double> DT_TERR_L = Program.DT_TERR[i][j];
                    ReadOnlySpan<float> FAC_L = Program.FAC[i][j];
                    ReadOnlySpan<double> NBZKP_L = Program.NBZKP[i][j];
                    ReadOnlySpan<double> PNBZKP_L = Program.PNBZKP[i][j];
                    double[] RADIATION_L = Program.RADIATION[i][j];
                    ReadOnlySpan<float> RHO_L = Program.RHOImm[i][j];
                    ReadOnlySpan<float> VOL_J = Program.VOLImm[i][j];

                    for (int k = 1; k <= NK - 1; ++k)
                    {
                        int KKAP = (int)(Math.Floor(NBZKP_L[k]));
                        int KKAM = KKAP - 1;
                        if (KKAP == 1)
                        {
                            KKAM = 1;
                        }

                        if (Program.ICSTR == true)
                        {
                            RADIATION_L[k] = ((DT_SOL_L[KKAP] + DT_TERR_L[KKAP]) * PNBZKP_L[k] +
                            (DT_SOL_L[KKAM] + DT_TERR_L[KKAM]) * (1 - PNBZKP_L[k])) *
                            FAC_L[k] * VOL_J[k] * RHO_L[k];
                        }
                        else
                        {
                            RADIATION_L[k] = 0.0;
                        }
                    }
                }
            });

            int range_parallel = NI / Program.pOptions.MaxDegreeOfParallelism - (StripeCounter % 6);
            range_parallel = Math.Max(Program.StripeWidth - (StripeCounter % 6), range_parallel); // min. Program.StripeWidth cells per processor
            StripeCounter++;
            range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
            //computation of the new temperature
            //Parallel.For(2, NI, Program.pOptions, i =>
            Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
            {
                double DIM;
                Span<double> PIM = stackalloc double[NK + 1];
                Span<double> QIM = stackalloc double[NK + 1];
                double help;
                double[] TN_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] TN_L;
                double[] TNiM_L;
                double[] TNiP_L;

                for (int i = range.Item1; i < range.Item2; ++i)
                {
                    int border = Math.Min(i - range.Item1, range.Item2 - 1 - i);
                    for (int j = 2; j <= NJ - 1; ++j)
                    {
                        ReadOnlySpan<double> A_PS_L = Program.A_PS[i][j];
                        ReadOnlySpan<float> AEAST_PS_L = Program.AEAST_PS[i][j];
                        ReadOnlySpan<float> ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        ReadOnlySpan<float> ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        ReadOnlySpan<float> AWEST_PS_L = Program.AWEST_PS[i][j];
                        ReadOnlySpan<float> AP0_PS_L = Program.AP0_PS[i][j];
                        ReadOnlySpan<double> B_PS_L = Program.B_PS[i][j];
                        ReadOnlySpan<double> C_PS_L = Program.C_PS[i][j];
                        ReadOnlySpan<float> FAC_L = Program.FAC[i][j];
                        ReadOnlySpan<float> RHO_L = Program.RHOImm[i][j];
                        ReadOnlySpan<double> T_L = Program.T[i][j];
                        ReadOnlySpan<double> RADIATION_L = RADIATION[i][j];
                        double WQU_AWQ = Program.WQU[i][j] - Program.AWQ[i][j];
                        ReadOnlySpan<double> TNjM_L = Program.TN[i][j - 1]; ReadOnlySpan<double> TNjP_L = Program.TN[i][j + 1];
                        //Avoid race conditions at the border cells of the sequential calculated stripes
                        if (border < 2)
                        {
                            TN_L = TN_LR;
                            Program.CopyArrayLockSource(Program.TN[i][j], TN_L);
                            TNiM_L = Program.GrammArrayPool.Rent(Program.NZ1);
                            TNiP_L = Program.GrammArrayPool.Rent(Program.NZ1);
                            Program.CopyArrayLockSource(Program.TN[i - 1][j], TNiM_L);
                            Program.CopyArrayLockSource(Program.TN[i + 1][j], TNiP_L);
                        }
                        else
                        {
                            TN_L = Program.TN[i][j];
                            TNiM_L = Program.TN[i - 1][j];
                            TNiP_L = Program.TN[i + 1][j];
                        }

                        for (int k = 1; k <= NK - 1; ++k)
                        {
                            DIM = AWEST_PS_L[k] * TNiM_L[k] + ASOUTH_PS_L[k] * TNjM_L[k] +
                                  AEAST_PS_L[k] * TNiP_L[k] + ANORTH_PS_L[k] * TNjP_L[k] +
                                  AP0_PS_L[k] * T_L[k] + RADIATION_L[k];
                            if ((k == 1) && (Program.ICTB == true))
                            {
                                DIM -= WQU_AWQ / Program.CPLUFT * FAC_L[k] * RHO_L[k];
                            }

                            //Recurrence formula
                            if (k > 1)
                            {
                                help = 1 / (A_PS_L[k] - C_PS_L[k] * PIM[k - 1]);
                                PIM[k] = B_PS_L[k] * help;
                                QIM[k] = (DIM + C_PS_L[k] * QIM[k - 1]) * help;
                            }
                            else
                            {
                                help = 1 / A_PS_L[k];
                                PIM[k] = B_PS_L[k] * help;
                                QIM[k] = DIM * help;
                            }
                        }
                        //Obtain new T-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            TN_L[k] += (Program.RELAXT * (PIM[k] * TN_L[k + 1] + QIM[k] - TN_L[k]));
                        }
                        if (border < 2)
                        {
                            Program.CopyArrayLockDest(TN_L, Program.TN[i][j]);
                            Program.GrammArrayPool.Return(TNiM_L);
                            Program.GrammArrayPool.Return(TNiP_L);
                        }

                    }
                }
                Program.GrammArrayPool.Return(TN_LR);
            });


            range_parallel = (int)(NI / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2 + 4);
            range_parallel = Math.Max(36 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
            range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
            //Parallel.For(2, NI, Program.pOptions, i1 =>
            Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
            {
                double DIM;
                Span<double> PIM = stackalloc double[NK + 1];
                Span<double> QIM = stackalloc double[NK + 1];
                double help;
                double[] TN_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] TN_L;
                double[] TNiM_L;
                double[] TNiP_L;

                for (int i = range.Item2 - 1; i >= range.Item1; --i)
                {
                    int border = Math.Min(i - range.Item1, range.Item2 - 1 - i);
                    for (int j = NJ - 1; j >= 2; --j)
                    {
                        ReadOnlySpan<double> A_PS_L = Program.A_PS[i][j];
                        ReadOnlySpan<float> AEAST_PS_L = Program.AEAST_PS[i][j];
                        ReadOnlySpan<float> ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        ReadOnlySpan<float> ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        ReadOnlySpan<float> AWEST_PS_L = Program.AWEST_PS[i][j];
                        ReadOnlySpan<float> AP0_PS_L = Program.AP0_PS[i][j];
                        ReadOnlySpan<double> B_PS_L = Program.B_PS[i][j];
                        ReadOnlySpan<double> C_PS_L = Program.C_PS[i][j];
                        ReadOnlySpan<float> FAC_L = Program.FAC[i][j];
                        ReadOnlySpan<float> RHO_L = Program.RHOImm[i][j];
                        ReadOnlySpan<double> T_L = Program.T[i][j];
                        ReadOnlySpan<double> RADIATION_L = RADIATION[i][j];
                        double WQU_AWQ = Program.WQU[i][j] - Program.AWQ[i][j];
                        ReadOnlySpan<double> TNjM_L = Program.TN[i][j - 1]; ReadOnlySpan<double> TNjP_L = Program.TN[i][j + 1];

                        if (border < 2)
                        {
                            TN_L = TN_LR;
                            Program.CopyArrayLockSource(Program.TN[i][j], TN_L);
                            TNiM_L = Program.GrammArrayPool.Rent(Program.NZ1);
                            TNiP_L = Program.GrammArrayPool.Rent(Program.NZ1);
                            Program.CopyArrayLockSource(Program.TN[i - 1][j], TNiM_L);
                            Program.CopyArrayLockSource(Program.TN[i + 1][j], TNiP_L);
                        }
                        else
                        {
                            TN_L = Program.TN[i][j];
                            TNiM_L = Program.TN[i - 1][j];
                            TNiP_L = Program.TN[i + 1][j];
                        }

                        for (int k = 1; k <= NK - 1; ++k)
                        {
                            DIM = AWEST_PS_L[k] * TNiM_L[k] + ASOUTH_PS_L[k] * TNjM_L[k] +
                                  AEAST_PS_L[k] * TNiP_L[k] + ANORTH_PS_L[k] * TNjP_L[k] +
                                  AP0_PS_L[k] * T_L[k] + RADIATION_L[k];
                            if ((k == 1) && (Program.ICTB == true))
                            {
                                DIM -= WQU_AWQ / Program.CPLUFT * FAC_L[k] * RHO_L[k];
                            }

                            //Recurrence formula
                            if (k > 1)
                            {
                                help = 1 / (A_PS_L[k] - C_PS_L[k] * PIM[k - 1]);
                                PIM[k] = B_PS_L[k] * help;
                                QIM[k] = (DIM + C_PS_L[k] * QIM[k - 1]) * help;
                            }
                            else
                            {
                                help = 1 / A_PS_L[k];
                                PIM[k] = B_PS_L[k] * help;
                                QIM[k] = DIM * help;
                            }
                        }

                        //Obtain new T-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            TN_L[k] += (Program.RELAXT * (PIM[k] * TN_L[k + 1] + QIM[k] - TN_L[k]));
                        }
                        if (border < 2)
                        {
                            Program.CopyArrayLockDest(TN_L, Program.TN[i][j]);
                            Program.GrammArrayPool.Return(TNiM_L);
                            Program.GrammArrayPool.Return(TNiP_L);
                        }

                    }
                }
                Program.GrammArrayPool.Return(TN_LR);
            });


            range_parallel = (int)(NI / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2 - 4);
            range_parallel = Math.Max(33 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
            range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
            //Parallel.For(2, NI, Program.pOptions, i1 =>
            Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
            {
                double DIM;
                Span<double> PIM = stackalloc double[NK + 1];
                Span<double> QIM = stackalloc double[NK + 1];
                double help;
                double[] TN_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] TN_L;
                double[] TNiM_L;
                double[] TNiP_L;

                for (int i = range.Item2 - 1; i >= range.Item1; --i)
                {
                    int border = Math.Min(i - range.Item1, range.Item2 - 1 - i);
                    for (int j = 2; j <= NJ - 1; ++j)
                    {
                        ReadOnlySpan<double> A_PS_L = Program.A_PS[i][j];
                        ReadOnlySpan<float> AEAST_PS_L = Program.AEAST_PS[i][j];
                        ReadOnlySpan<float> ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        ReadOnlySpan<float> ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        ReadOnlySpan<float> AWEST_PS_L = Program.AWEST_PS[i][j];
                        ReadOnlySpan<float> AP0_PS_L = Program.AP0_PS[i][j];
                        ReadOnlySpan<double> B_PS_L = Program.B_PS[i][j];
                        ReadOnlySpan<double> C_PS_L = Program.C_PS[i][j];
                        ReadOnlySpan<float> FAC_L = Program.FAC[i][j];
                        ReadOnlySpan<float> RHO_L = Program.RHOImm[i][j];
                        ReadOnlySpan<double> T_L = Program.T[i][j];
                        ReadOnlySpan<double> RADIATION_L = RADIATION[i][j];
                        double WQU_AWQ = Program.WQU[i][j] - Program.AWQ[i][j];
                        ReadOnlySpan<double> TNjM_L = Program.TN[i][j - 1]; ReadOnlySpan<double> TNjP_L = Program.TN[i][j + 1];

                        if (border < 2)
                        {
                            TN_L = TN_LR;
                            Program.CopyArrayLockSource(Program.TN[i][j], TN_L);
                            TNiM_L = Program.GrammArrayPool.Rent(Program.NZ1);
                            TNiP_L = Program.GrammArrayPool.Rent(Program.NZ1);
                            Program.CopyArrayLockSource(Program.TN[i - 1][j], TNiM_L);
                            Program.CopyArrayLockSource(Program.TN[i + 1][j], TNiP_L);
                        }
                        else
                        {
                            TN_L = Program.TN[i][j];
                            TNiM_L = Program.TN[i - 1][j];
                            TNiP_L = Program.TN[i + 1][j];
                        }

                        for (int k = 1; k <= NK - 1; ++k)
                        {
                            DIM = AWEST_PS_L[k] * TNiM_L[k] + ASOUTH_PS_L[k] * TNjM_L[k] +
                                  AEAST_PS_L[k] * TNiP_L[k] + ANORTH_PS_L[k] * TNjP_L[k] +
                                  AP0_PS_L[k] * T_L[k] + RADIATION_L[k];
                            if ((k == 1) && (Program.ICTB == true))
                            {
                                DIM -= WQU_AWQ / Program.CPLUFT * FAC_L[k] * RHO_L[k];
                            }

                            //Recurrence formula
                            if (k > 1)
                            {
                                help = 1 / (A_PS_L[k] - C_PS_L[k] * PIM[k - 1]);
                                PIM[k] = B_PS_L[k] * help;
                                QIM[k] = (DIM + C_PS_L[k] * QIM[k - 1]) * help;
                            }
                            else
                            {
                                help = 1 / A_PS_L[k];
                                PIM[k] = B_PS_L[k] * help;
                                QIM[k] = DIM * help;
                            }
                        }

                        //Obtain new T-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            TN_L[k] += (Program.RELAXT * (PIM[k] * TN_L[k + 1] + QIM[k] - TN_L[k]));
                        }
                        if (border < 2)
                        {
                            Program.CopyArrayLockDest(TN_L, Program.TN[i][j]);
                            Program.GrammArrayPool.Return(TNiM_L);
                            Program.GrammArrayPool.Return(TNiP_L);
                        }

                    }
                }
                Program.GrammArrayPool.Return(TN_LR);
            });

            range_parallel = NI / Program.pOptions.MaxDegreeOfParallelism - (StripeCounter % 6);
            range_parallel = Math.Max(Program.StripeWidth - (StripeCounter % 6), range_parallel); // min. Program.StripeWidth cells per processor
            StripeCounter++;
            range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
                                                           //Parallel.For(2, NI, Program.pOptions, i =>
            Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
            {
                double DIM;
                Span<double> PIM = stackalloc double[NK + 1];
                Span<double> QIM = stackalloc double[NK + 1];
                double help;
                double[] TN_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] TN_L;
                double[] TNiM_L;
                double[] TNiP_L;

                for (int i = range.Item1; i < range.Item2; ++i)
                {
                    int border = Math.Min(i - range.Item1, range.Item2 - 1 - i);
                    for (int j = NJ - 1; j >= 2; --j)
                    {
                        ReadOnlySpan<double> A_PS_L = Program.A_PS[i][j];
                        ReadOnlySpan<float> AEAST_PS_L = Program.AEAST_PS[i][j];
                        ReadOnlySpan<float> ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        ReadOnlySpan<float> ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        ReadOnlySpan<float> AWEST_PS_L = Program.AWEST_PS[i][j];
                        ReadOnlySpan<float> AP0_PS_L = Program.AP0_PS[i][j];
                        ReadOnlySpan<double> B_PS_L = Program.B_PS[i][j];
                        ReadOnlySpan<double> C_PS_L = Program.C_PS[i][j];
                        ReadOnlySpan<float> FAC_L = Program.FAC[i][j];
                        ReadOnlySpan<float> RHO_L = Program.RHOImm[i][j];
                        ReadOnlySpan<double> T_L = Program.T[i][j];
                        ReadOnlySpan<double> RADIATION_L = RADIATION[i][j];
                        double WQU_AWQ = Program.WQU[i][j] - Program.AWQ[i][j];
                        ReadOnlySpan<double> TNjM_L = Program.TN[i][j - 1]; ReadOnlySpan<double> TNjP_L = Program.TN[i][j + 1];

                        if (border < 2)
                        {
                            TN_L = TN_LR;
                            Program.CopyArrayLockSource(Program.TN[i][j], TN_L);
                            TNiM_L = Program.GrammArrayPool.Rent(Program.NZ1);
                            TNiP_L = Program.GrammArrayPool.Rent(Program.NZ1);
                            Program.CopyArrayLockSource(Program.TN[i - 1][j], TNiM_L);
                            Program.CopyArrayLockSource(Program.TN[i + 1][j], TNiP_L);
                        }
                        else
                        {
                            TN_L = Program.TN[i][j];
                            TNiM_L = Program.TN[i - 1][j];
                            TNiP_L = Program.TN[i + 1][j];
                        }

                        for (int k = 1; k <= NK - 1; ++k)
                        {
                            DIM = AWEST_PS_L[k] * TNiM_L[k] + ASOUTH_PS_L[k] * TNjM_L[k] +
                                  AEAST_PS_L[k] * TNiP_L[k] + ANORTH_PS_L[k] * TNjP_L[k] +
                                  AP0_PS_L[k] * T_L[k] + RADIATION_L[k];
                            if ((k == 1) && (Program.ICTB == true))
                            {
                                DIM -= WQU_AWQ / Program.CPLUFT * FAC_L[k] * RHO_L[k];
                            }

                            //Recurrence formula
                            if (k > 1)
                            {
                                help = 1 / (A_PS_L[k] - C_PS_L[k] * PIM[k - 1]);
                                PIM[k] = B_PS_L[k] * help;
                                QIM[k] = (DIM + C_PS_L[k] * QIM[k - 1]) * help;
                            }
                            else
                            {
                                help = 1 / A_PS_L[k];
                                PIM[k] = B_PS_L[k] * help;
                                QIM[k] = DIM * help;
                            }
                        }

                        //Obtain new T-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            TN_L[k] += (Program.RELAXT * (PIM[k] * TN_L[k + 1] + QIM[k] - TN_L[k]));
                        }
                        if (border < 2)
                        {
                            Program.CopyArrayLockDest(TN_L, Program.TN[i][j]);
                            Program.GrammArrayPool.Return(TNiM_L);
                            Program.GrammArrayPool.Return(TNiP_L);
                        }

                    }
                }
                Program.GrammArrayPool.Return(TN_LR);
            });

            range_parallel = NJ / Program.pOptions.MaxDegreeOfParallelism - (StripeCounter % 6);
            range_parallel = Math.Max(Program.StripeWidth - (StripeCounter % 6), range_parallel); // min. Program.StripeWidth cells per processor
            StripeCounter++;
            range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
            //Parallel.For(2, NJ, Program.pOptions, j =>
            Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
            {
                double DIM;
                Span<double> PIM = stackalloc double[NK + 1];
                Span<double> QIM = stackalloc double[NK + 1];
                double help;
                double[] TN_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] TN_L;
                double[] TNjM_L;
                double[] TNjP_L;

                for (int j = range.Item1; j < range.Item2; ++j)
                {
                    int border = Math.Min(j - range.Item1, range.Item2 - 1 - j);
                    for (int i = 2; i <= NI - 1; ++i)
                    {
                        ReadOnlySpan<double> A_PS_L = Program.A_PS[i][j];
                        ReadOnlySpan<float> AEAST_PS_L = Program.AEAST_PS[i][j];
                        ReadOnlySpan<float> ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        ReadOnlySpan<float> ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        ReadOnlySpan<float> AWEST_PS_L = Program.AWEST_PS[i][j];
                        ReadOnlySpan<float> AP0_PS_L = Program.AP0_PS[i][j];
                        ReadOnlySpan<double> B_PS_L = Program.B_PS[i][j];
                        ReadOnlySpan<double> C_PS_L = Program.C_PS[i][j];
                        ReadOnlySpan<float> FAC_L = Program.FAC[i][j];
                        ReadOnlySpan<float> RHO_L = Program.RHOImm[i][j];
                        ReadOnlySpan<double> T_L = Program.T[i][j];
                        ReadOnlySpan<double> RADIATION_L = RADIATION[i][j];
                        double WQU_AWQ = Program.WQU[i][j] - Program.AWQ[i][j];
                        ReadOnlySpan<double> TNiM_L = Program.TN[i - 1][j]; ReadOnlySpan<double> TNiP_L = Program.TN[i + 1][j];

                        if (border < 2)
                        {
                            TN_L = TN_LR;
                            Program.CopyArrayLockSource(Program.TN[i][j], TN_L);
                            TNjM_L = Program.GrammArrayPool.Rent(Program.NZ1);
                            TNjP_L = Program.GrammArrayPool.Rent(Program.NZ1);
                            Program.CopyArrayLockSource(Program.TN[i][j - 1], TNjM_L);
                            Program.CopyArrayLockSource(Program.TN[i][j + 1], TNjP_L);
                        }
                        else
                        {
                            TN_L = Program.TN[i][j];
                            TNjM_L = Program.TN[i][j - 1];
                            TNjP_L = Program.TN[i][j + 1];
                        }

                        for (int k = 1; k <= NK - 1; ++k)
                        {
                            DIM = AWEST_PS_L[k] * TNiM_L[k] + ASOUTH_PS_L[k] * TNjM_L[k] +
                                  AEAST_PS_L[k] * TNiP_L[k] + ANORTH_PS_L[k] * TNjP_L[k] +
                                  AP0_PS_L[k] * T_L[k] + RADIATION_L[k];
                            if ((k == 1) && (Program.ICTB == true))
                            {
                                DIM -= WQU_AWQ / Program.CPLUFT * FAC_L[k] * RHO_L[k];
                            }

                            //Recurrence formula
                            if (k > 1)
                            {
                                help = 1 / (A_PS_L[k] - C_PS_L[k] * PIM[k - 1]);
                                PIM[k] = B_PS_L[k] * help;
                                QIM[k] = (DIM + C_PS_L[k] * QIM[k - 1]) * help;
                            }
                            else
                            {
                                help = 1 / A_PS_L[k];
                                PIM[k] = B_PS_L[k] * help;
                                QIM[k] = DIM * help;
                            }
                        }

                        //Obtain new T-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            TN_L[k] += (Program.RELAXT * (PIM[k] * TN_L[k + 1] + QIM[k] - TN_L[k]));
                        }
                        if (border < 2)
                        {
                            Program.CopyArrayLockDest(TN_L, Program.TN[i][j]);
                            Program.GrammArrayPool.Return(TNjM_L);
                            Program.GrammArrayPool.Return(TNjP_L);
                        }

                    }
                }
                Program.GrammArrayPool.Return(TN_LR);
            });

            range_parallel = (int)(NJ / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2 - 4);
            range_parallel = Math.Max(36 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
            range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
            //Parallel.For(2, NJ, Program.pOptions, j1 =>
            Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
            {
                double DIM;
                Span<double> PIM = stackalloc double[NK + 1];
                Span<double> QIM = stackalloc double[NK + 1];
                double help;
                double[] TN_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] TN_L;
                double[] TNjM_L;
                double[] TNjP_L;

                for (int j = range.Item2 - 1; j >= range.Item1; --j)
                {
                    int border = Math.Min(j - range.Item1, range.Item2 - 1 - j);
                    for (int i = NI - 1; i >= 2; --i)
                    {
                        ReadOnlySpan<double> A_PS_L = Program.A_PS[i][j];
                        ReadOnlySpan<float> AEAST_PS_L = Program.AEAST_PS[i][j];
                        ReadOnlySpan<float> ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        ReadOnlySpan<float> ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        ReadOnlySpan<float> AWEST_PS_L = Program.AWEST_PS[i][j];
                        ReadOnlySpan<float> AP0_PS_L = Program.AP0_PS[i][j];
                        ReadOnlySpan<double> B_PS_L = Program.B_PS[i][j];
                        ReadOnlySpan<double> C_PS_L = Program.C_PS[i][j];
                        ReadOnlySpan<float> FAC_L = Program.FAC[i][j];
                        ReadOnlySpan<float> RHO_L = Program.RHOImm[i][j];
                        ReadOnlySpan<double> T_L = Program.T[i][j];
                        ReadOnlySpan<double> RADIATION_L = RADIATION[i][j];
                        double WQU_AWQ = Program.WQU[i][j] - Program.AWQ[i][j];
                        ReadOnlySpan<double> TNiM_L = Program.TN[i - 1][j]; ReadOnlySpan<double> TNiP_L = Program.TN[i + 1][j];

                        if (border < 2)
                        {
                            TN_L = TN_LR;
                            Program.CopyArrayLockSource(Program.TN[i][j], TN_L);
                            TNjM_L = Program.GrammArrayPool.Rent(Program.NZ1);
                            TNjP_L = Program.GrammArrayPool.Rent(Program.NZ1);
                            Program.CopyArrayLockSource(Program.TN[i][j - 1], TNjM_L);
                            Program.CopyArrayLockSource(Program.TN[i][j + 1], TNjP_L);
                        }
                        else
                        {
                            TN_L = Program.TN[i][j];
                            TNjM_L = Program.TN[i][j - 1];
                            TNjP_L = Program.TN[i][j + 1];
                        }

                        for (int k = 1; k <= NK - 1; ++k)
                        {
                            DIM = AWEST_PS_L[k] * TNiM_L[k] + ASOUTH_PS_L[k] * TNjM_L[k] +
                                  AEAST_PS_L[k] * TNiP_L[k] + ANORTH_PS_L[k] * TNjP_L[k] +
                                  AP0_PS_L[k] * T_L[k] + RADIATION_L[k];
                            if ((k == 1) && (Program.ICTB == true))
                            {
                                DIM -= WQU_AWQ / Program.CPLUFT * FAC_L[k] * RHO_L[k];
                            }

                            //Recurrence formula
                            if (k > 1)
                            {
                                help = 1 / (A_PS_L[k] - C_PS_L[k] * PIM[k - 1]);
                                PIM[k] = B_PS_L[k] * help;
                                QIM[k] = (DIM + C_PS_L[k] * QIM[k - 1]) * help;
                            }
                            else
                            {
                                help = 1 / A_PS_L[k];
                                PIM[k] = B_PS_L[k] * help;
                                QIM[k] = DIM * help;
                            }
                        }

                        //Obtain new T-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            TN_L[k] += (Program.RELAXT * (PIM[k] * TN_L[k + 1] + QIM[k] - TN_L[k]));
                        }
                        if (border < 2)
                        {
                            Program.CopyArrayLockDest(TN_L, Program.TN[i][j]);
                            Program.GrammArrayPool.Return(TNjM_L);
                            Program.GrammArrayPool.Return(TNjP_L);
                        }

                    }
                }
                Program.GrammArrayPool.Return(TN_LR);
            });

            range_parallel = (int)(NJ / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2 + 4);
            range_parallel = Math.Max(33 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
            range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
            //Parallel.For(2, NJ, Program.pOptions, j =>
            Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
            {
                double DIM;
                Span<double> PIM = stackalloc double[NK + 1];
                Span<double> QIM = stackalloc double[NK + 1];
                double help;
                double[] TN_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] TN_L;
                double[] TNjM_L;
                double[] TNjP_L;

                for (int j = range.Item1; j < range.Item2; ++j)
                {
                    int border = Math.Min(j - range.Item1, range.Item2 - 1 - j);
                    for (int i = NI - 1; i >= 2; --i)
                    {
                        ReadOnlySpan<double> A_PS_L = Program.A_PS[i][j];
                        ReadOnlySpan<float> AEAST_PS_L = Program.AEAST_PS[i][j];
                        ReadOnlySpan<float> ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        ReadOnlySpan<float> ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        ReadOnlySpan<float> AWEST_PS_L = Program.AWEST_PS[i][j];
                        ReadOnlySpan<float> AP0_PS_L = Program.AP0_PS[i][j];
                        ReadOnlySpan<double> B_PS_L = Program.B_PS[i][j];
                        ReadOnlySpan<double> C_PS_L = Program.C_PS[i][j];
                        ReadOnlySpan<float> FAC_L = Program.FAC[i][j];
                        ReadOnlySpan<float> RHO_L = Program.RHOImm[i][j];
                        ReadOnlySpan<double> T_L = Program.T[i][j];
                        ReadOnlySpan<double> RADIATION_L = RADIATION[i][j];
                        double WQU_AWQ = Program.WQU[i][j] - Program.AWQ[i][j];
                        ReadOnlySpan<double> TNiM_L = Program.TN[i - 1][j]; ReadOnlySpan<double> TNiP_L = Program.TN[i + 1][j];

                        if (border < 2)
                        {
                            TN_L = TN_LR;
                            Program.CopyArrayLockSource(Program.TN[i][j], TN_L);
                            TNjM_L = Program.GrammArrayPool.Rent(Program.NZ1);
                            TNjP_L = Program.GrammArrayPool.Rent(Program.NZ1);
                            Program.CopyArrayLockSource(Program.TN[i][j - 1], TNjM_L);
                            Program.CopyArrayLockSource(Program.TN[i][j + 1], TNjP_L);
                        }
                        else
                        {
                            TN_L = Program.TN[i][j];
                            TNjM_L = Program.TN[i][j - 1];
                            TNjP_L = Program.TN[i][j + 1];
                        }

                        for (int k = 1; k <= NK - 1; ++k)
                        {
                            DIM = AWEST_PS_L[k] * TNiM_L[k] + ASOUTH_PS_L[k] * TNjM_L[k] +
                                   AEAST_PS_L[k] * TNiP_L[k] + ANORTH_PS_L[k] * TNjP_L[k] +
                                   AP0_PS_L[k] * T_L[k] + RADIATION_L[k];
                            if ((k == 1) && (Program.ICTB == true))
                            {
                                DIM -= WQU_AWQ / Program.CPLUFT * FAC_L[k] * RHO_L[k];
                            }

                            //Recurrence formula
                            if (k > 1)
                            {
                                help = 1 / (A_PS_L[k] - C_PS_L[k] * PIM[k - 1]);
                                PIM[k] = B_PS_L[k] * help;
                                QIM[k] = (DIM + C_PS_L[k] * QIM[k - 1]) * help;
                            }
                            else
                            {
                                help = 1 / A_PS_L[k];
                                PIM[k] = B_PS_L[k] * help;
                                QIM[k] = DIM * help;
                            }
                        }

                        //Obtain new T-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            TN_L[k] += (Program.RELAXT * (PIM[k] * TN_L[k + 1] + QIM[k] - TN_L[k]));
                        }
                        if (border < 2)
                        {
                            Program.CopyArrayLockDest(TN_L, Program.TN[i][j]);
                            Program.GrammArrayPool.Return(TNjM_L);
                            Program.GrammArrayPool.Return(TNjP_L);
                        }

                    }
                }
                Program.GrammArrayPool.Return(TN_LR);
            });

            range_parallel = NJ / Program.pOptions.MaxDegreeOfParallelism - (StripeCounter % 6);
            range_parallel = Math.Max(Program.StripeWidth - (StripeCounter % 6), range_parallel); // min. Program.StripeWidth cells per processor
            StripeCounter++;
            range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
            //Parallel.For(2, NJ, Program.pOptions, j1 =>
            Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
            {
                double DIM;
                Span<double> PIM = stackalloc double[NK + 1];
                Span<double> QIM = stackalloc double[NK + 1];
                double help;
                double[] TN_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] TN_L;
                double[] TNjM_L;
                double[] TNjP_L;

                for (int j = range.Item2 - 1; j >= range.Item1; --j)
                {
                    int border = Math.Min(j - range.Item1, range.Item2 - 1 - j);
                    for (int i = 2; i <= NI - 1; ++i)
                    {
                        ReadOnlySpan<double> A_PS_L = Program.A_PS[i][j];
                        ReadOnlySpan<float> AEAST_PS_L = Program.AEAST_PS[i][j];
                        ReadOnlySpan<float> ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        ReadOnlySpan<float> ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        ReadOnlySpan<float> AWEST_PS_L = Program.AWEST_PS[i][j];
                        ReadOnlySpan<float> AP0_PS_L = Program.AP0_PS[i][j];
                        ReadOnlySpan<double> B_PS_L = Program.B_PS[i][j];
                        ReadOnlySpan<double> C_PS_L = Program.C_PS[i][j];
                        ReadOnlySpan<float> FAC_L = Program.FAC[i][j];
                        ReadOnlySpan<float> RHO_L = Program.RHOImm[i][j];
                        ReadOnlySpan<double> T_L = Program.T[i][j];
                        ReadOnlySpan<double> RADIATION_L = RADIATION[i][j];
                        double WQU_AWQ = Program.WQU[i][j] - Program.AWQ[i][j];
                        ReadOnlySpan<double> TNiM_L = Program.TN[i - 1][j]; ReadOnlySpan<double> TNiP_L = Program.TN[i + 1][j];

                        if (border < 2)
                        {
                            TN_L = TN_LR;
                            Program.CopyArrayLockSource(Program.TN[i][j], TN_L);
                            TNjM_L = Program.GrammArrayPool.Rent(Program.NZ1);
                            TNjP_L = Program.GrammArrayPool.Rent(Program.NZ1);
                            Program.CopyArrayLockSource(Program.TN[i][j - 1], TNjM_L);
                            Program.CopyArrayLockSource(Program.TN[i][j + 1], TNjP_L);
                        }
                        else
                        {
                            TN_L = Program.TN[i][j];
                            TNjM_L = Program.TN[i][j - 1];
                            TNjP_L = Program.TN[i][j + 1];
                        }

                        for (int k = 1; k <= NK - 1; ++k)
                        {
                            DIM = AWEST_PS_L[k] * TNiM_L[k] + ASOUTH_PS_L[k] * TNjM_L[k] +
                                  AEAST_PS_L[k] * TNiP_L[k] + ANORTH_PS_L[k] * TNjP_L[k] +
                                  AP0_PS_L[k] * T_L[k] + RADIATION_L[k];
                            if ((k == 1) && (Program.ICTB == true))
                            {
                                DIM -= WQU_AWQ / Program.CPLUFT * FAC_L[k] * RHO_L[k];
                            }

                            //Recurrence formula
                            if (k > 1)
                            {
                                help = 1 / (A_PS_L[k] - C_PS_L[k] * PIM[k - 1]);
                                PIM[k] = B_PS_L[k] * help;
                                QIM[k] = (DIM + C_PS_L[k] * QIM[k - 1]) * help;
                            }
                            else
                            {
                                help = 1 / A_PS_L[k];
                                PIM[k] = B_PS_L[k] * help;
                                QIM[k] = DIM * help;
                            }
                        }

                        //Obtain new T-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            TN_L[k] += (Program.RELAXT * (PIM[k] * TN_L[k + 1] + QIM[k] - TN_L[k]));
                        }
                        if (border < 2)
                        {
                            Program.CopyArrayLockDest(TN_L, Program.TN[i][j]);
                            Program.GrammArrayPool.Return(TNjM_L);
                            Program.GrammArrayPool.Return(TNjP_L);
                        }

                    }
                }
                Program.GrammArrayPool.Return(TN_LR);
            });

        }
    }
}
