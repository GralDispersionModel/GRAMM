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

namespace GRAMM_2001
{
    partial class Program
    {
        public static void Fimp_calculate(int NI, int NJ, int NK)
        {

            int range_parallel = (int)(NI / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
            range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
            range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
                                                           //computation of the new specific humidity
                                                           //Parallel.For(2, NI, Program.pOptions, i =>
            Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
            {
                double DIM;
                double[] PIM = new double[NK + 1];
                double[] QIM = new double[NK + 1];
                double help;

                for (int i = range.Item1; i < range.Item2; i++)
                {
                    for (int j = 2; j <= NJ - 1; j++)
                    {
                        float[] AWEST_PS_L = Program.AWEST_PS[i][j];
                        float[] ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        float[] AEAST_PS_L = Program.AEAST_PS[i][j];
                        float[] ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        float[] AP0_PS_L = Program.AP0_PS[i][j];
                        double[] QU_L = Program.QU[i][j];
                        double[] A_PS_L = Program.A_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        float[] RHO_L = Program.RHO[i][j];
                        double[] QUN_L = Program.QUN[i][j];
                        float[] AREA_L = Program.AREA[i][j];
                        double UST_L = Program.UST[i][j];
                        float XWQ_L = Program.XWQ[i][j];
                        double QUG_L = Program.QUG[i][j];

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            DIM = AWEST_PS_L[k] * Program.QUN[i - 1][j][k] + ASOUTH_PS_L[k] * Program.QUN[i][j - 1][k] +
                                AEAST_PS_L[k] * Program.QUN[i + 1][j][k] + ANORTH_PS_L[k] * Program.QUN[i][j + 1][k] +
                                AP0_PS_L[k] * QU_L[k];
                            if ((k == 1) && (Program.ICTB == true)) DIM -= XWQ_L * UST_L * (QU_L[k] - QUG_L) * AREA_L[k] / 1000 * RHO_L[k];

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
                            help += (Program.RELAXT * (PIM[k] * QUN_L[k + 1] + QIM[k] - help));
                            QUN_L[k] = help;

                        }
                    }
                }
            });

            range_parallel = (int)(NI / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
            range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
            range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
                                                           //Parallel.For(2, NI, Program.pOptions, i1 =>
            Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
            {
                double DIM;
                double[] PIM = new double[NK + 1];
                double[] QIM = new double[NK + 1];
                double help;

                for (int i1 = range.Item1; i1 < range.Item2; i1++)
                {
                    int i = NI - i1 + 1;
                    for (int j = NJ - 1; j >= 2; j--)
                    {
                        float[] AWEST_PS_L = Program.AWEST_PS[i][j];
                        float[] ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        float[] AEAST_PS_L = Program.AEAST_PS[i][j];
                        float[] ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        float[] AP0_PS_L = Program.AP0_PS[i][j];
                        double[] QU_L = Program.QU[i][j];
                        double[] A_PS_L = Program.A_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        float[] RHO_L = Program.RHO[i][j];
                        double[] QUN_L = Program.QUN[i][j];
                        float[] AREA_L = Program.AREA[i][j];
                        double UST_L = Program.UST[i][j];
                        float XWQ_L = Program.XWQ[i][j];
                        double QUG_L = Program.QUG[i][j];

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            DIM = AWEST_PS_L[k] * Program.QUN[i - 1][j][k] + ASOUTH_PS_L[k] * Program.QUN[i][j - 1][k] +
                                AEAST_PS_L[k] * Program.QUN[i + 1][j][k] + ANORTH_PS_L[k] * Program.QUN[i][j + 1][k] +
                                AP0_PS_L[k] * QU_L[k];
                            if ((k == 1) && (Program.ICTB == true)) DIM -= XWQ_L * UST_L * (QU_L[k] - QUG_L) * AREA_L[k] / 1000 * RHO_L[k];

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
                            help += (Program.RELAXT * (PIM[k] * QUN_L[k + 1] + QIM[k] - help));
                            QUN_L[k] = help;

                        }
                    }
                }
            });

            range_parallel = (int)(NI / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2 + 6);
            range_parallel = Math.Max(36 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
            range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
                                                           //Parallel.For(2, NI, Program.pOptions, i1 =>
            Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
            {
                double DIM;
                double[] PIM = new double[NK + 1];
                double[] QIM = new double[NK + 1];
                double help;

                for (int i1 = range.Item1; i1 < range.Item2; i1++)
                {
                    int i = NI - i1 + 1;
                    for (int j = 2; j <= NJ - 1; j++)
                    {
                        float[] AWEST_PS_L = Program.AWEST_PS[i][j];
                        float[] ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        float[] AEAST_PS_L = Program.AEAST_PS[i][j];
                        float[] ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        float[] AP0_PS_L = Program.AP0_PS[i][j];
                        double[] QU_L = Program.QU[i][j];
                        double[] A_PS_L = Program.A_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        float[] RHO_L = Program.RHO[i][j];
                        double[] QUN_L = Program.QUN[i][j];
                        float[] AREA_L = Program.AREA[i][j];
                        double UST_L = Program.UST[i][j];
                        float XWQ_L = Program.XWQ[i][j];
                        double QUG_L = Program.QUG[i][j];

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            DIM = AWEST_PS_L[k] * Program.QUN[i - 1][j][k] + ASOUTH_PS_L[k] * Program.QUN[i][j - 1][k] +
                                AEAST_PS_L[k] * Program.QUN[i + 1][j][k] + ANORTH_PS_L[k] * Program.QUN[i][j + 1][k] +
                                AP0_PS_L[k] * QU_L[k];
                            if ((k == 1) && (Program.ICTB == true)) DIM -= XWQ_L * UST_L * (QU_L[k] - QUG_L) * AREA_L[k] / 1000 * RHO_L[k];

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
                            help += (Program.RELAXT * (PIM[k] * QUN_L[k + 1] + QIM[k] - help));
                            QUN_L[k] = help;

                        }
                    }
                }
            });

            range_parallel = (int)(NI / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2 + 3);
            range_parallel = Math.Max(33 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
            range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
                                                           //Parallel.For(2, NI, Program.pOptions, i =>
            Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
            {
                double DIM;
                double[] PIM = new double[NK + 1];
                double[] QIM = new double[NK + 1];
                double help;

                for (int i = range.Item1; i < range.Item2; i++)
                {
                    for (int j = NJ - 1; j >= 2; j--)
                    {
                        float[] AWEST_PS_L = Program.AWEST_PS[i][j];
                        float[] ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        float[] AEAST_PS_L = Program.AEAST_PS[i][j];
                        float[] ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        float[] AP0_PS_L = Program.AP0_PS[i][j];
                        double[] QU_L = Program.QU[i][j];
                        double[] A_PS_L = Program.A_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        float[] RHO_L = Program.RHO[i][j];
                        double[] QUN_L = Program.QUN[i][j];
                        float[] AREA_L = Program.AREA[i][j];
                        double UST_L = Program.UST[i][j];
                        float XWQ_L = Program.XWQ[i][j];
                        double QUG_L = Program.QUG[i][j];

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            DIM = AWEST_PS_L[k] * Program.QUN[i - 1][j][k] + ASOUTH_PS_L[k] * Program.QUN[i][j - 1][k] +
                                AEAST_PS_L[k] * Program.QUN[i + 1][j][k] + ANORTH_PS_L[k] * Program.QUN[i][j + 1][k] +
                                AP0_PS_L[k] * QU_L[k];
                            if ((k == 1) && (Program.ICTB == true)) DIM -= XWQ_L * UST_L * (QU_L[k] - QUG_L) * AREA_L[k] / 1000 * RHO_L[k];

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
                            help += (Program.RELAXT * (PIM[k] * QUN_L[k + 1] + QIM[k] - help));
                            QUN_L[k] = help;

                        }
                    }
                }
            });

            range_parallel = (int)(NJ / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
            range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
            range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
                                                           //Parallel.For(2, NJ, Program.pOptions, j =>
            Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
            {
                double DIM;
                double[] PIM = new double[NK + 1];
                double[] QIM = new double[NK + 1];
                double help;

                for (int j = range.Item1; j < range.Item2; j++)
                {
                    for (int i = 2; i <= NI - 1; i++)
                    {
                        float[] AWEST_PS_L = Program.AWEST_PS[i][j];
                        float[] ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        float[] AEAST_PS_L = Program.AEAST_PS[i][j];
                        float[] ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        float[] AP0_PS_L = Program.AP0_PS[i][j];
                        double[] QU_L = Program.QU[i][j];
                        double[] A_PS_L = Program.A_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        float[] RHO_L = Program.RHO[i][j];
                        double[] QUN_L = Program.QUN[i][j];
                        float[] AREA_L = Program.AREA[i][j];
                        double UST_L = Program.UST[i][j];
                        float XWQ_L = Program.XWQ[i][j];
                        double QUG_L = Program.QUG[i][j];

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            DIM = AWEST_PS_L[k] * Program.QUN[i - 1][j][k] + ASOUTH_PS_L[k] * Program.QUN[i][j - 1][k] +
                                AEAST_PS_L[k] * Program.QUN[i + 1][j][k] + ANORTH_PS_L[k] * Program.QUN[i][j + 1][k] +
                                AP0_PS_L[k] * QU_L[k];
                            if ((k == 1) && (Program.ICTB == true)) DIM -= XWQ_L * UST_L * (QU_L[k] - QUG_L) * AREA_L[k] / 1000 * RHO_L[k];

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
                            help += (Program.RELAXT * (PIM[k] * QUN_L[k + 1] + QIM[k] - help));
                            QUN_L[k] = help;

                        }
                    }
                }
            });

            range_parallel = (int)(NJ / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
            range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
            range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
                                                           //Parallel.For(2, NJ, Program.pOptions, j1 =>
            Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
            {
                double DIM;
                double[] PIM = new double[NK + 1];
                double[] QIM = new double[NK + 1];
                double help;

                for (int j1 = range.Item1; j1 < range.Item2; j1++)
                {
                    int j = NJ - j1 + 1;

                    for (int i = NI - 1; i >= 2; i--)
                    {
                        float[] AWEST_PS_L = Program.AWEST_PS[i][j];
                        float[] ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        float[] AEAST_PS_L = Program.AEAST_PS[i][j];
                        float[] ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        float[] AP0_PS_L = Program.AP0_PS[i][j];
                        double[] QU_L = Program.QU[i][j];
                        double[] A_PS_L = Program.A_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        float[] RHO_L = Program.RHO[i][j];
                        double[] QUN_L = Program.QUN[i][j];
                        float[] AREA_L = Program.AREA[i][j];
                        double UST_L = Program.UST[i][j];
                        float XWQ_L = Program.XWQ[i][j];
                        double QUG_L = Program.QUG[i][j];

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            DIM = AWEST_PS_L[k] * Program.QUN[i - 1][j][k] + ASOUTH_PS_L[k] * Program.QUN[i][j - 1][k] +
                                AEAST_PS_L[k] * Program.QUN[i + 1][j][k] + ANORTH_PS_L[k] * Program.QUN[i][j + 1][k] +
                                AP0_PS_L[k] * QU_L[k];
                            if ((k == 1) && (Program.ICTB == true)) DIM -= XWQ_L * UST_L * (QU_L[k] - QUG_L) * AREA_L[k] / 1000 * RHO_L[k];

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
                            help += (Program.RELAXT * (PIM[k] * QUN_L[k + 1] + QIM[k] - help));
                            QUN_L[k] = help;

                        }
                    }
                }
            });

            range_parallel = (int)(NJ / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
            range_parallel = Math.Max(36 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
            range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
                                                           //Parallel.For(2, NJ, Program.pOptions, j =>
            Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
            {
                double DIM;
                double[] PIM = new double[NK + 1];
                double[] QIM = new double[NK + 1];
                double help;

                for (int j = range.Item1; j < range.Item2; j++)
                {
                    for (int i = NI - 1; i >= 2; i--)
                    {
                        float[] AWEST_PS_L = Program.AWEST_PS[i][j];
                        float[] ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        float[] AEAST_PS_L = Program.AEAST_PS[i][j];
                        float[] ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        float[] AP0_PS_L = Program.AP0_PS[i][j];
                        double[] QU_L = Program.QU[i][j];
                        double[] A_PS_L = Program.A_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        float[] RHO_L = Program.RHO[i][j];
                        double[] QUN_L = Program.QUN[i][j];
                        float[] AREA_L = Program.AREA[i][j];
                        double UST_L = Program.UST[i][j];
                        float XWQ_L = Program.XWQ[i][j];
                        double QUG_L = Program.QUG[i][j];

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            DIM = AWEST_PS_L[k] * Program.QUN[i - 1][j][k] + ASOUTH_PS_L[k] * Program.QUN[i][j - 1][k] +
                                AEAST_PS_L[k] * Program.QUN[i + 1][j][k] + ANORTH_PS_L[k] * Program.QUN[i][j + 1][k] +
                                AP0_PS_L[k] * QU_L[k];
                            if ((k == 1) && (Program.ICTB == true)) DIM -= XWQ_L * UST_L * (QU_L[k] - QUG_L) * AREA_L[k] / 1000 * RHO_L[k];

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
                            help += (Program.RELAXT * (PIM[k] * QUN_L[k + 1] + QIM[k] - help));
                            QUN_L[k] = help;

                        }
                    }
                }
            });

            range_parallel = (int)(NJ / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2 + 6);
            range_parallel = Math.Max(36 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
            range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
                                                           //Parallel.For(2, NJ, Program.pOptions, j1 =>
            Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
            {
                double DIM;
                double[] PIM = new double[NK + 1];
                double[] QIM = new double[NK + 1];
                double help;

                for (int j1 = range.Item1; j1 < range.Item2; j1++)
                {
                    int j = NJ - j1 + 1;
                    for (int i = 2; i <= NI - 1; i++)
                    {
                        float[] AWEST_PS_L = Program.AWEST_PS[i][j];
                        float[] ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        float[] AEAST_PS_L = Program.AEAST_PS[i][j];
                        float[] ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        float[] AP0_PS_L = Program.AP0_PS[i][j];
                        double[] QU_L = Program.QU[i][j];
                        double[] A_PS_L = Program.A_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        float[] RHO_L = Program.RHO[i][j];
                        double[] QUN_L = Program.QUN[i][j];
                        float[] AREA_L = Program.AREA[i][j];
                        double UST_L = Program.UST[i][j];
                        float XWQ_L = Program.XWQ[i][j];
                        double QUG_L = Program.QUG[i][j];

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            DIM = AWEST_PS_L[k] * Program.QUN[i - 1][j][k] + ASOUTH_PS_L[k] * Program.QUN[i][j - 1][k] +
                                AEAST_PS_L[k] * Program.QUN[i + 1][j][k] + ANORTH_PS_L[k] * Program.QUN[i][j + 1][k] +
                                AP0_PS_L[k] * QU_L[k];
                            if ((k == 1) && (Program.ICTB == true)) DIM -= XWQ_L * UST_L * (QU_L[k] - QUG_L) * AREA_L[k] / 1000 * RHO_L[k];

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
                            help += (Program.RELAXT * (PIM[k] * QUN_L[k + 1] + QIM[k] - help));
                            QUN_L[k] = help;

                        }
                    }
                }
            });
        }
    }
}
