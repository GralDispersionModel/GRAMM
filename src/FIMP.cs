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
        /// <summary>
        /// Calculate the new specific humidity 
        /// </summary>
        /// <param name="NI"></param>
        /// <param name="NJ"></param>
        /// <param name="NK"></param>
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        public static void Fimp_calculate(int NI, int NJ, int NK)
        {
            //computation of the new specific humidity
            //Parallel.For(2, NI, Program.pOptions, i =>
            Parallel.ForEach(PartitionerI[Program.StripeCounter++ % PartitionerI.Count], range =>
            {
                double DIM;
                Span<double> PIM = stackalloc double[NK];
                Span<double> QIM = stackalloc double[NK];
                double help;
                double[] QUN_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] QUNiM_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] QUNiP_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] QUN_L;
                double[] QUNiM_L;
                double[] QUNiP_L;
                ReadOnlySpan<float> AWEST_PS_L;
                ReadOnlySpan<float> ASOUTH_PS_L;
                ReadOnlySpan<float> AEAST_PS_L;
                ReadOnlySpan<float> ANORTH_PS_L;
                ReadOnlySpan<float> AP0_PS_L;
                ReadOnlySpan<double> QUNjM_L;
                ReadOnlySpan<double> QUNjP_L;
                ReadOnlySpan<double> A_PS_L;
                ReadOnlySpan<double> B_PS_L;
                ReadOnlySpan<double> C_PS_L;
                ReadOnlySpan<float> RHO_L;
                ReadOnlySpan<float> AREA_L;

                for (int i = range.Item1; i < range.Item2; ++i)
                {
                    int border = Math.Min(i - range.Item1, range.Item2 - 1 - i);
                    for (int j = 2; j <= NJ - 1; ++j)
                    {
                        AWEST_PS_L  = Program.AWEST_PS[i][j];
                        ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        AEAST_PS_L  = Program.AEAST_PS[i][j];
                        ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        AP0_PS_L    = Program.AP0_PS[i][j];
                        ReadOnlySpan<double> QU_L = Program.QU[i][j];
                        A_PS_L = Program.A_PS[i][j];
                        B_PS_L = Program.B_PS[i][j];
                        C_PS_L = Program.C_PS[i][j];
                        RHO_L = Program.RHO[i][j];
                        //Avoid race conditions at the border cells of the sequential calculated stripes
                        if (border < 1)
                        {
                            QUN_L = QUN_LR;
                            QUNiM_L = QUNiM_LR;
                            QUNiP_L = QUNiP_LR;
                            FastCopy.CopyArrayLockSource(Program.QUN[i][j], QUN_L);
                            FastCopy.CopyArrayLockSource(Program.QUN[i - 1][j], QUNiM_L);
                            FastCopy.CopyArrayLockSource(Program.QUN[i + 1][j], QUNiP_L);
                        }
                        else
                        {
                            QUN_L = Program.QUN[i][j];
                            QUNiM_L = Program.QUN[i - 1][j];
                            QUNiP_L = Program.QUN[i + 1][j];
                        }
                        QUNjM_L = Program.QUN[i][j - 1];
                        QUNjP_L = Program.QUN[i][j + 1];
                        AREA_L = Program.AREAImm[i][j].AsSpan();
                        double UST_L = Program.UST[i][j];
                        float XWQ_L = Program.XWQ[i][j];
                        double QUG_L = Program.QUG[i][j];

                        for (int k = 1; k < QIM.Length; ++k)
                        {
                            DIM = AWEST_PS_L[k] * QUNiM_L[k] + ASOUTH_PS_L[k] * QUNjM_L[k] +
                                      AEAST_PS_L[k] * QUNiP_L[k] + ANORTH_PS_L[k] * QUNjP_L[k] +
                                      AP0_PS_L[k] * QU_L[k];
                            if ((k == 1) && (Program.ICTB == true))
                            {
                                DIM -= XWQ_L * UST_L * (QU_L[k] - QUG_L) * AREA_L[k] / 1000 * RHO_L[k];
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

                        //Obtain new QU-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            help = QUN_L[k];
                            help += Program.RELAXT * (PIM[k] * QUN_L[k + 1] + QIM[k] - help);
                            QUN_L[k] = help;
                        }
                        if (border < 1)
                        {
                            FastCopy.CopyArrayLockDest(QUN_L, Program.QUN[i][j]);
                        }
                    }
                }
                Program.GrammArrayPool.Return(QUN_LR);
                Program.GrammArrayPool.Return(QUNiM_LR);
                Program.GrammArrayPool.Return(QUNiP_LR);
            });

            //Parallel.For(2, NI, Program.pOptions, i1 =>
            Parallel.ForEach(PartitionerI[Program.StripeCounter++ % PartitionerI.Count], range =>
            {
                double DIM;
                Span<double> PIM = stackalloc double[NK];
                Span<double> QIM = stackalloc double[NK];
                double help;
                double[] QUN_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] QUNiM_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] QUNiP_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] QUN_L;
                double[] QUNiM_L;
                double[] QUNiP_L;
                ReadOnlySpan<float> AWEST_PS_L;
                ReadOnlySpan<float> ASOUTH_PS_L;
                ReadOnlySpan<float> AEAST_PS_L;
                ReadOnlySpan<float> ANORTH_PS_L;
                ReadOnlySpan<float> AP0_PS_L;
                ReadOnlySpan<double> QUNjM_L;
                ReadOnlySpan<double> QUNjP_L;
                ReadOnlySpan<double> A_PS_L;
                ReadOnlySpan<double> B_PS_L;
                ReadOnlySpan<double> C_PS_L;
                ReadOnlySpan<float> RHO_L;
                ReadOnlySpan<float> AREA_L;

                for (int i = range.Item2 - 1; i >= range.Item1; --i)
                {
                    int border = Math.Min(i - range.Item1, range.Item2 - 1 - i);
                    for (int j = NJ - 1; j >= 2; --j)
                    {
                        AWEST_PS_L = Program.AWEST_PS[i][j];
                        ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        AEAST_PS_L = Program.AEAST_PS[i][j];
                        ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        AP0_PS_L = Program.AP0_PS[i][j];
                        if (border < 1)
                        {
                            QUN_L = QUN_LR;
                            QUNiM_L = QUNiM_LR;
                            QUNiP_L = QUNiP_LR;
                            FastCopy.CopyArrayLockSource(Program.QUN[i][j], QUN_L);
                            FastCopy.CopyArrayLockSource(Program.QUN[i - 1][j], QUNiM_L);
                            FastCopy.CopyArrayLockSource(Program.QUN[i + 1][j], QUNiP_L);
                        }
                        else
                        {
                            QUN_L = Program.QUN[i][j];
                            QUNiM_L = Program.QUN[i - 1][j];
                            QUNiP_L = Program.QUN[i + 1][j];
                        }

                        QUNjM_L = Program.QUN[i][j - 1];
                        QUNjP_L = Program.QUN[i][j + 1];
                        double[] QU_L = Program.QU[i][j];
                        A_PS_L = Program.A_PS[i][j];
                        B_PS_L = Program.B_PS[i][j];
                        C_PS_L = Program.C_PS[i][j];
                        RHO_L = Program.RHO[i][j];
                        AREA_L = Program.AREAImm[i][j].AsSpan();
                        double UST_L = Program.UST[i][j];
                        float XWQ_L = Program.XWQ[i][j];
                        double QUG_L = Program.QUG[i][j];

                        for (int k = 1; k < QIM.Length; ++k)
                        {
                            DIM = AWEST_PS_L[k] * QUNiM_L[k] + ASOUTH_PS_L[k] * QUNjM_L[k] +
                                      AEAST_PS_L[k] * QUNiP_L[k] + ANORTH_PS_L[k] * QUNjP_L[k] +
                                      AP0_PS_L[k] * QU_L[k];
                            if ((k == 1) && (Program.ICTB == true))
                            {
                                DIM -= XWQ_L * UST_L * (QU_L[k] - QUG_L) * AREA_L[k] / 1000 * RHO_L[k];
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

                        //Obtain new QU-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            help = QUN_L[k];
                            help += Program.RELAXT * (PIM[k] * QUN_L[k + 1] + QIM[k] - help);
                            QUN_L[k] = help;
                        }
                        if (border < 1)
                        {
                            FastCopy.CopyArrayLockDest(QUN_L, Program.QUN[i][j]);
                        }
                    }
                }
                Program.GrammArrayPool.Return(QUN_LR);
                Program.GrammArrayPool.Return(QUNiM_LR);
                Program.GrammArrayPool.Return(QUNiP_LR);
            });

            //Parallel.For(2, NJ, Program.pOptions, j =>
            Parallel.ForEach(PartitionerJ[Program.StripeCounter++ % PartitionerJ.Count], range =>
            {
                double DIM;
                Span<double> PIM = stackalloc double[NK];
                Span<double> QIM = stackalloc double[NK];
                double help;
                double[] QUN_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] QUNjM_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] QUNjP_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] QUN_L;
                double[] QUNjM_L;
                double[] QUNjP_L;
                ReadOnlySpan<float> AWEST_PS_L;
                ReadOnlySpan<float> ASOUTH_PS_L;
                ReadOnlySpan<float> AEAST_PS_L;
                ReadOnlySpan<float> ANORTH_PS_L;
                ReadOnlySpan<float> AP0_PS_L;
                ReadOnlySpan<double> QUNiM_L;
                ReadOnlySpan<double> QUNiP_L;
                ReadOnlySpan<double> A_PS_L;
                ReadOnlySpan<double> B_PS_L;
                ReadOnlySpan<double> C_PS_L;
                ReadOnlySpan<float> RHO_L;
                ReadOnlySpan<float> AREA_L;

                for (int j = range.Item1; j < range.Item2; ++j)
                {
                    int border = Math.Min(j - range.Item1, range.Item2 - 1 - j);
                    for (int i = 2; i <= NI - 1; ++i)
                    {
                        AWEST_PS_L = Program.AWEST_PS[i][j];
                        ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        AEAST_PS_L = Program.AEAST_PS[i][j];
                        ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        AP0_PS_L = Program.AP0_PS[i][j];
                        ReadOnlySpan<double> QU_L = Program.QU[i][j];
                        if (border < 1)
                        {
                            QUN_L = QUN_LR;
                            QUNjM_L = QUNjM_LR;
                            QUNjP_L = QUNjP_LR;
                            FastCopy.CopyArrayLockSource(Program.QUN[i][j], QUN_L);
                            FastCopy.CopyArrayLockSource(Program.QUN[i][j - 1], QUNjM_L);
                            FastCopy.CopyArrayLockSource(Program.QUN[i][j + 1], QUNjP_L);
                        }
                        else
                        {
                            QUN_L = Program.QUN[i][j];
                            QUNjM_L = Program.QUN[i][j - 1];
                            QUNjP_L = Program.QUN[i][j + 1];
                        }
                        QUNiM_L = Program.QUN[i - 1][j];
                        QUNiP_L = Program.QUN[i + 1][j];
                        A_PS_L = Program.A_PS[i][j];
                        B_PS_L = Program.B_PS[i][j];
                        C_PS_L = Program.C_PS[i][j];
                        RHO_L = Program.RHO[i][j];
                        AREA_L = Program.AREAImm[i][j].AsSpan();
                        double UST_L = Program.UST[i][j];
                        float XWQ_L = Program.XWQ[i][j];
                        double QUG_L = Program.QUG[i][j];

                        for (int k = 1; k < QIM.Length; ++k)
                        {
                            DIM = AWEST_PS_L[k] * QUNiM_L[k] + ASOUTH_PS_L[k] * QUNjM_L[k] +
                                      AEAST_PS_L[k] * QUNiP_L[k] + ANORTH_PS_L[k] * QUNjP_L[k] +
                                      AP0_PS_L[k] * QU_L[k];
                            if ((k == 1) && (Program.ICTB == true))
                            {
                                DIM -= XWQ_L * UST_L * (QU_L[k] - QUG_L) * AREA_L[k] / 1000 * RHO_L[k];
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

                        //Obtain new QU-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            help = QUN_L[k];
                            help += Program.RELAXT * (PIM[k] * QUN_L[k + 1] + QIM[k] - help);
                            QUN_L[k] = help;
                        }
                        if (border < 1)
                        {
                            FastCopy.CopyArrayLockDest(QUN_L, Program.QUN[i][j]);
                        }
                    }
                }
                Program.GrammArrayPool.Return(QUN_LR);
                Program.GrammArrayPool.Return(QUNjM_LR);
                Program.GrammArrayPool.Return(QUNjP_LR);
            });

            //Parallel.For(2, NJ, Program.pOptions, j1 =>
            Parallel.ForEach(PartitionerJ[Program.StripeCounter++ % PartitionerJ.Count], range =>
            {
                double DIM;
                Span<double> PIM = stackalloc double[NK];
                Span<double> QIM = stackalloc double[NK];
                double help;
                double[] QUN_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] QUNjM_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] QUNjP_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] QUN_L;
                double[] QUNjM_L;
                double[] QUNjP_L;
                ReadOnlySpan<float> AWEST_PS_L;
                ReadOnlySpan<float> ASOUTH_PS_L;
                ReadOnlySpan<float> AEAST_PS_L;
                ReadOnlySpan<float> ANORTH_PS_L;
                ReadOnlySpan<float> AP0_PS_L;
                ReadOnlySpan<double> QUNiM_L;
                ReadOnlySpan<double> QUNiP_L;
                ReadOnlySpan<double> A_PS_L;
                ReadOnlySpan<double> B_PS_L;
                ReadOnlySpan<double> C_PS_L;
                ReadOnlySpan<float> RHO_L;
                ReadOnlySpan<float> AREA_L;

                for (int j = range.Item2 - 1; j >= range.Item1; --j)
                {
                    int border = Math.Min(j - range.Item1, range.Item2 - 1 - j);
                    for (int i = NI - 1; i >= 2; --i)
                    {
                        AWEST_PS_L = Program.AWEST_PS[i][j];
                        ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        AEAST_PS_L = Program.AEAST_PS[i][j];
                        ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        AP0_PS_L = Program.AP0_PS[i][j];
                        ReadOnlySpan<double> QU_L = Program.QU[i][j];
                        if (border < 1)
                        {
                            QUN_L = QUN_LR;
                            QUNjM_L = QUNjM_LR;
                            QUNjP_L = QUNjP_LR;
                            FastCopy.CopyArrayLockSource(Program.QUN[i][j], QUN_L);
                            FastCopy.CopyArrayLockSource(Program.QUN[i][j - 1], QUNjM_L);
                            FastCopy.CopyArrayLockSource(Program.QUN[i][j + 1], QUNjP_L);
                        }
                        else
                        {
                            QUN_L = Program.QUN[i][j];
                            QUNjM_L = Program.QUN[i][j - 1];
                            QUNjP_L = Program.QUN[i][j + 1];
                        }

                        QUNiM_L = Program.QUN[i - 1][j];
                        QUNiP_L = Program.QUN[i + 1][j];
                        A_PS_L = Program.A_PS[i][j];
                        B_PS_L = Program.B_PS[i][j];
                        C_PS_L = Program.C_PS[i][j];
                        RHO_L = Program.RHO[i][j];
                        AREA_L = Program.AREAImm[i][j].AsSpan();
                        double UST_L = Program.UST[i][j];
                        float XWQ_L = Program.XWQ[i][j];
                        double QUG_L = Program.QUG[i][j];

                        for (int k = 1; k < QIM.Length; ++k)
                        {
                            DIM = AWEST_PS_L[k] * QUNiM_L[k] + ASOUTH_PS_L[k] * QUNjM_L[k] +
                                       AEAST_PS_L[k] * QUNiP_L[k] + ANORTH_PS_L[k] * QUNjP_L[k] +
                                       AP0_PS_L[k] * QU_L[k];
                            if ((k == 1) && (Program.ICTB == true))
                            {
                                DIM -= XWQ_L * UST_L * (QU_L[k] - QUG_L) * AREA_L[k] / 1000 * RHO_L[k];
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

                        //Obtain new QU-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            help = QUN_L[k];
                            help += Program.RELAXT * (PIM[k] * QUN_L[k + 1] + QIM[k] - help);
                            QUN_L[k] = help;
                        }
                        if (border < 1)
                        {
                            FastCopy.CopyArrayLockDest(QUN_L, Program.QUN[i][j]);
                        }
                    }
                }
                Program.GrammArrayPool.Return(QUN_LR);
                Program.GrammArrayPool.Return(QUNjM_LR);
                Program.GrammArrayPool.Return(QUNjP_LR);
            });

            //Parallel.For(2, NI, Program.pOptions, i1 =>
            Parallel.ForEach(PartitionerI[Program.StripeCounter++ % PartitionerI.Count], range =>
            {
                double DIM;
                Span<double> PIM = stackalloc double[NK];
                Span<double> QIM = stackalloc double[NK];
                double help;
                double[] QUN_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] QUNiM_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] QUNiP_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] QUN_L;
                double[] QUNiM_L;
                double[] QUNiP_L;
                ReadOnlySpan<float> AWEST_PS_L;
                ReadOnlySpan<float> ASOUTH_PS_L;
                ReadOnlySpan<float> AEAST_PS_L;
                ReadOnlySpan<float> ANORTH_PS_L;
                ReadOnlySpan<float> AP0_PS_L;
                ReadOnlySpan<double> QUNjM_L;
                ReadOnlySpan<double> QUNjP_L;
                ReadOnlySpan<double> A_PS_L;
                ReadOnlySpan<double> B_PS_L;
                ReadOnlySpan<double> C_PS_L;
                ReadOnlySpan<float> RHO_L;
                ReadOnlySpan<float> AREA_L;

                for (int i = range.Item2 - 1; i >= range.Item1; --i)
                {
                    int border = Math.Min(i - range.Item1, range.Item2 - 1 - i);
                    for (int j = 2; j <= NJ - 1; ++j)
                    {
                        AWEST_PS_L = Program.AWEST_PS[i][j];
                        ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        AEAST_PS_L = Program.AEAST_PS[i][j];
                        ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        AP0_PS_L = Program.AP0_PS[i][j];
                        ReadOnlySpan<double> QU_L = Program.QU[i][j];
                        if (border < 1)
                        {
                            QUN_L = QUN_LR;
                            QUNiM_L = QUNiM_LR;
                            QUNiP_L = QUNiP_LR;
                            FastCopy.CopyArrayLockSource(Program.QUN[i][j], QUN_L);
                            FastCopy.CopyArrayLockSource(Program.QUN[i - 1][j], QUNiM_L);
                            FastCopy.CopyArrayLockSource(Program.QUN[i + 1][j], QUNiP_L);
                        }
                        else
                        {
                            QUN_L = Program.QUN[i][j];
                            QUNiM_L = Program.QUN[i - 1][j];
                            QUNiP_L = Program.QUN[i + 1][j];
                        }

                        QUNjM_L = Program.QUN[i][j - 1];
                        QUNjP_L = Program.QUN[i][j + 1];
                        A_PS_L = Program.A_PS[i][j];
                        B_PS_L = Program.B_PS[i][j];
                        C_PS_L = Program.C_PS[i][j];
                        RHO_L = Program.RHO[i][j];
                        AREA_L = Program.AREAImm[i][j].AsSpan();
                        double UST_L = Program.UST[i][j];
                        float XWQ_L = Program.XWQ[i][j];
                        double QUG_L = Program.QUG[i][j];

                        for (int k = 1; k < QIM.Length; ++k)
                        {
                            DIM = AWEST_PS_L[k] * QUNiM_L[k] + ASOUTH_PS_L[k] * QUNjM_L[k] +
                                  AEAST_PS_L[k] * QUNiP_L[k] + ANORTH_PS_L[k] * QUNjP_L[k] +
                                  AP0_PS_L[k] * QU_L[k];
                            if ((k == 1) && (Program.ICTB == true))
                            {
                                DIM -= XWQ_L * UST_L * (QU_L[k] - QUG_L) * AREA_L[k] / 1000 * RHO_L[k];
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

                        //Obtain new QU-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            help = QUN_L[k];
                            help += Program.RELAXT * (PIM[k] * QUN_L[k + 1] + QIM[k] - help);
                            QUN_L[k] = help;
                        }
                        if (border < 1)
                        {
                            FastCopy.CopyArrayLockDest(QUN_L, Program.QUN[i][j]);
                        }
                    }
                }
                Program.GrammArrayPool.Return(QUN_LR);
                Program.GrammArrayPool.Return(QUNiM_LR);
                Program.GrammArrayPool.Return(QUNiP_LR);
            });

            //Parallel.For(2, NI, Program.pOptions, i =>
            Parallel.ForEach(PartitionerI[Program.StripeCounter++ % PartitionerI.Count], range =>
            {
                double DIM;
                Span<double> PIM = stackalloc double[NK];
                Span<double> QIM = stackalloc double[NK];
                double help;
                double[] QUN_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] QUNiM_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] QUNiP_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] QUN_L;
                double[] QUNiM_L;
                double[] QUNiP_L;
                ReadOnlySpan<float> AWEST_PS_L;
                ReadOnlySpan<float> ASOUTH_PS_L;
                ReadOnlySpan<float> AEAST_PS_L;
                ReadOnlySpan<float> ANORTH_PS_L;
                ReadOnlySpan<float> AP0_PS_L;
                ReadOnlySpan<double> QUNjM_L;
                ReadOnlySpan<double> QUNjP_L;
                ReadOnlySpan<double> A_PS_L;
                ReadOnlySpan<double> B_PS_L;
                ReadOnlySpan<double> C_PS_L;
                ReadOnlySpan<float> RHO_L;
                ReadOnlySpan<float> AREA_L;

                for (int i = range.Item1; i < range.Item2; ++i)
                {
                    int border = Math.Min(i - range.Item1, range.Item2 - 1 - i);
                    for (int j = NJ - 1; j >= 2; --j)
                    {
                        AWEST_PS_L = Program.AWEST_PS[i][j];
                        ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        AEAST_PS_L = Program.AEAST_PS[i][j];
                        ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        AP0_PS_L = Program.AP0_PS[i][j];
                        ReadOnlySpan<double> QU_L = Program.QU[i][j];
                        if (border < 1)
                        {
                            QUN_L = QUN_LR;
                            QUNiM_L = QUNiM_LR;
                            QUNiP_L = QUNiP_LR;
                            FastCopy.CopyArrayLockSource(Program.QUN[i][j], QUN_L);
                            FastCopy.CopyArrayLockSource(Program.QUN[i - 1][j], QUNiM_L);
                            FastCopy.CopyArrayLockSource(Program.QUN[i + 1][j], QUNiP_L);
                        }
                        else
                        {
                            QUN_L = Program.QUN[i][j];
                            QUNiM_L = Program.QUN[i - 1][j];
                            QUNiP_L = Program.QUN[i + 1][j];
                        }

                        QUNjM_L = Program.QUN[i][j - 1];
                        QUNjP_L = Program.QUN[i][j + 1];
                        A_PS_L = Program.A_PS[i][j];
                        B_PS_L = Program.B_PS[i][j];
                        C_PS_L = Program.C_PS[i][j];
                        RHO_L = Program.RHO[i][j];
                        AREA_L = Program.AREAImm[i][j].AsSpan();
                        double UST_L = Program.UST[i][j];
                        float XWQ_L = Program.XWQ[i][j];
                        double QUG_L = Program.QUG[i][j];

                        for (int k = 1; k < QIM.Length; ++k)
                        {
                            DIM = AWEST_PS_L[k] * QUNiM_L[k] + ASOUTH_PS_L[k] * QUNjM_L[k] +
                                      AEAST_PS_L[k] * QUNiP_L[k] + ANORTH_PS_L[k] * QUNjP_L[k] +
                                      AP0_PS_L[k] * QU_L[k];
                            if ((k == 1) && (Program.ICTB == true))
                            {
                                DIM -= XWQ_L * UST_L * (QU_L[k] - QUG_L) * AREA_L[k] / 1000 * RHO_L[k];
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

                        //Obtain new QU-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            help = QUN_L[k];
                            help += Program.RELAXT * (PIM[k] * QUN_L[k + 1] + QIM[k] - help);
                            QUN_L[k] = help;
                        }
                        if (border < 1)
                        {
                            FastCopy.CopyArrayLockDest(QUN_L, Program.QUN[i][j]);
                        }
                    }
                }
                Program.GrammArrayPool.Return(QUN_LR);
                Program.GrammArrayPool.Return(QUNiM_LR);
                Program.GrammArrayPool.Return(QUNiP_LR);
            });

            //Parallel.For(2, NJ, Program.pOptions, j =>
            Parallel.ForEach(PartitionerJ[Program.StripeCounter++ % PartitionerJ.Count], range =>
            {
                double DIM;
                Span<double> PIM = stackalloc double[NK];
                Span<double> QIM = stackalloc double[NK];
                double help;
                double[] QUN_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] QUNjM_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] QUNjP_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] QUN_L;
                double[] QUNjM_L;
                double[] QUNjP_L;
                ReadOnlySpan<float> AWEST_PS_L;
                ReadOnlySpan<float> ASOUTH_PS_L;
                ReadOnlySpan<float> AEAST_PS_L;
                ReadOnlySpan<float> ANORTH_PS_L;
                ReadOnlySpan<float> AP0_PS_L;

                for (int j = range.Item1; j < range.Item2; ++j)
                {
                    int border = Math.Min(j - range.Item1, range.Item2 - 1 - j);
                    for (int i = NI - 1; i >= 2; --i)
                    {
                        AWEST_PS_L = Program.AWEST_PS[i][j];
                        ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        AEAST_PS_L = Program.AEAST_PS[i][j];
                        ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        AP0_PS_L = Program.AP0_PS[i][j];
                        ReadOnlySpan<double> QU_L = Program.QU[i][j];
                        if (border < 1)
                        {
                            QUN_L = QUN_LR;
                            QUNjM_L = QUNjM_LR;
                            QUNjP_L = QUNjP_LR;
                            FastCopy.CopyArrayLockSource(Program.QUN[i][j], QUN_L);
                            FastCopy.CopyArrayLockSource(Program.QUN[i][j - 1], QUNjM_L);
                            FastCopy.CopyArrayLockSource(Program.QUN[i][j + 1], QUNjP_L);
                        }
                        else
                        {
                            QUN_L = Program.QUN[i][j];
                            QUNjM_L = Program.QUN[i][j - 1];
                            QUNjP_L = Program.QUN[i][j + 1];
                        }
                        ReadOnlySpan<double> QUNiM_L = Program.QUN[i - 1][j];
                        ReadOnlySpan<double> QUNiP_L = Program.QUN[i + 1][j];
                        ReadOnlySpan<double> A_PS_L = Program.A_PS[i][j];
                        ReadOnlySpan<double> B_PS_L = Program.B_PS[i][j];
                        ReadOnlySpan<double> C_PS_L = Program.C_PS[i][j];
                        ReadOnlySpan<float> RHO_L = Program.RHO[i][j];
                        ReadOnlySpan<float> AREA_L = Program.AREAImm[i][j].AsSpan();
                        double UST_L = Program.UST[i][j];
                        float XWQ_L = Program.XWQ[i][j];
                        double QUG_L = Program.QUG[i][j];

                        for (int k = 1; k < QIM.Length; ++k)
                        {
                            DIM = AWEST_PS_L[k] * QUNiM_L[k] + ASOUTH_PS_L[k] * QUNjM_L[k] +
                                      AEAST_PS_L[k] * QUNiP_L[k] + ANORTH_PS_L[k] * QUNjP_L[k] +
                                      AP0_PS_L[k] * QU_L[k];
                            if ((k == 1) && (Program.ICTB == true))
                            {
                                DIM -= XWQ_L * UST_L * (QU_L[k] - QUG_L) * AREA_L[k] / 1000 * RHO_L[k];
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

                        //Obtain new QU-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            help = QUN_L[k];
                            help += Program.RELAXT * (PIM[k] * QUN_L[k + 1] + QIM[k] - help);
                            QUN_L[k] = help;
                        }
                        if (border < 1)
                        {
                            FastCopy.CopyArrayLockDest(QUN_L, Program.QUN[i][j]);
                        }
                    }
                }
                Program.GrammArrayPool.Return(QUN_LR);
                Program.GrammArrayPool.Return(QUNjM_LR);
                Program.GrammArrayPool.Return(QUNjP_LR);
            });

            //Parallel.For(2, NJ, Program.pOptions, j1 =>
            Parallel.ForEach(PartitionerJ[Program.StripeCounter++ % PartitionerJ.Count], range =>
            {
                double DIM;
                Span<double> PIM = stackalloc double[NK];
                Span<double> QIM = stackalloc double[NK];
                double help;
                double[] QUN_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] QUNjM_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] QUNjP_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                double[] QUN_L;
                double[] QUNjM_L;
                double[] QUNjP_L;
                ReadOnlySpan<float> AWEST_PS_L;
                ReadOnlySpan<float> ASOUTH_PS_L;
                ReadOnlySpan<float> AEAST_PS_L;
                ReadOnlySpan<float> ANORTH_PS_L;
                ReadOnlySpan<float> AP0_PS_L;

                for (int j = range.Item2 - 1; j >= range.Item1; --j)
                {
                    int border = Math.Min(j - range.Item1, range.Item2 - 1 - j);
                    for (int i = 2; i <= NI - 1; ++i)
                    {
                        AWEST_PS_L = Program.AWEST_PS[i][j];
                        ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        AEAST_PS_L = Program.AEAST_PS[i][j];
                        ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        AP0_PS_L = Program.AP0_PS[i][j];
                        ReadOnlySpan<double> QU_L = Program.QU[i][j];
                        if (border < 1)
                        {
                            QUN_L = QUN_LR;
                            QUNjM_L = QUNjM_LR;
                            QUNjP_L = QUNjP_LR;
                            FastCopy.CopyArrayLockSource(Program.QUN[i][j], QUN_L);
                            FastCopy.CopyArrayLockSource(Program.QUN[i][j - 1], QUNjM_L);
                            FastCopy.CopyArrayLockSource(Program.QUN[i][j + 1], QUNjP_L);
                        }
                        else
                        {
                            QUN_L = Program.QUN[i][j];
                            QUNjM_L = Program.QUN[i][j - 1];
                            QUNjP_L = Program.QUN[i][j + 1];
                        }
                        ReadOnlySpan<double> QUNiM_L = Program.QUN[i - 1][j];
                        ReadOnlySpan<double> QUNiP_L = Program.QUN[i + 1][j];
                        ReadOnlySpan<double> A_PS_L = Program.A_PS[i][j];
                        ReadOnlySpan<double> B_PS_L = Program.B_PS[i][j];
                        ReadOnlySpan<double> C_PS_L = Program.C_PS[i][j];
                        ReadOnlySpan<float> RHO_L = Program.RHO[i][j];
                        ReadOnlySpan<float> AREA_L = Program.AREAImm[i][j].AsSpan();
                        double UST_L = Program.UST[i][j];
                        float XWQ_L = Program.XWQ[i][j];
                        double QUG_L = Program.QUG[i][j];

                        for (int k = 1; k < QIM.Length; ++k)
                        {
                            DIM = AWEST_PS_L[k] * QUNiM_L[k] + ASOUTH_PS_L[k] * QUNjM_L[k] +
                                      AEAST_PS_L[k] * QUNiP_L[k] + ANORTH_PS_L[k] * QUNjP_L[k] +
                                      AP0_PS_L[k] * QU_L[k];
                            if ((k == 1) && (Program.ICTB == true))
                            {
                                DIM -= XWQ_L * UST_L * (QU_L[k] - QUG_L) * AREA_L[k] / 1000 * RHO_L[k];
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

                        //Obtain new QU-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            help = QUN_L[k];
                            help += Program.RELAXT * (PIM[k] * QUN_L[k + 1] + QIM[k] - help);
                            QUN_L[k] = help;
                        }
                        if (border < 1)
                        {
                            FastCopy.CopyArrayLockDest(QUN_L, Program.QUN[i][j]);
                        }
                    }
                }
                Program.GrammArrayPool.Return(QUN_LR);
                Program.GrammArrayPool.Return(QUNjM_LR);
                Program.GrammArrayPool.Return(QUNjP_LR);
            });
        }
    }
}
