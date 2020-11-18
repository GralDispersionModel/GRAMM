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
using System.Runtime.CompilerServices;

namespace GRAMM_2001
{
    partial class Program
    {
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        public static void WatVap_calculate(int NI, int NJ, int NK)
        {
            //computation of temperature change due to latent heat
            int range_parallel = (int)(NI / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
            range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
            range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
                                                           //computation of the new specific humidity
                                                           //Parallel.For(2, NI, Program.pOptions, i =>
            Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
            {
                for (int i = range.Item1; i < range.Item2; i++)
                {
                    for (int j = 2; j <= NJ - 1; j++)
                    {
                        double[] WAT_VAPN_L = Program.WAT_VAPN[i][j];
                        double[] WAT_VAP_L = Program.WAT_VAP[i][j];
                        double[] TABS_L = Program.TABS[i][j];
                        double[] QUN_L = Program.QUN[i][j];
                        double[] TN_L = Program.TN[i][j];
                        float[] FACTOR_L = Program.FACTOR[i][j];
                        double WAT_VAP_integral = 0.0;
                        for (int k = 1; k < NK; k++)
                        {
                            //saturation vapour pressure - Magnus formulae
                            double EGSAT = 611.2 * Math.Exp(17.269 * (TABS_L[k] - 273.16) / (TABS_L[k] - 30.04));
                            //maximum water vapour content possible -> assuming condensation above 80% relative humidity
                            double QGSAT = EGSAT / 461.5 / TABS_L[k];
                            //change in temperature
                            if (QUN_L[k] > QGSAT * 1000.0)
                            {
                                //condensation
                                double WATVAP_COND = (QUN_L[k] - QGSAT * 1000.0);
                                TN_L[k] += WATVAP_COND * 2.5 * FACTOR_L[k];
                                WAT_VAP_L[k] += WATVAP_COND;
                                QUN_L[k] -= WATVAP_COND;
                            }
                            else
                            {
                                //evaporation
                                double WATVAP_EVA = Math.Min(WAT_VAP_L[k], QGSAT * 1000.0 - QUN_L[k]);
                                TN_L[k] -= WATVAP_EVA * 2.5 * FACTOR_L[k];
                                WAT_VAP_L[k] -= WATVAP_EVA;
                                WAT_VAP_L[k] = Math.Max(WAT_VAP_L[k], 0.0);
                                QUN_L[k] += WATVAP_EVA;
                            }
                            WAT_VAP_integral += WAT_VAP_L[k];
                        }
                        //Program.CLOUDS[i][j] = 0;
                        if (WAT_VAP_integral > 1.0)
                        {
                            //Program.CLOUDS[i][j] = 1;
                        }
                    }
                }
            });

            range_parallel = (int)(NI / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
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
                        double[] A_PS_L = Program.A_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        double[] WAT_VAPN_L = Program.WAT_VAPN[i][j];
                        double[] WAT_VAP_L = Program.WAT_VAP[i][j];

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            DIM = AWEST_PS_L[k] * Program.WAT_VAPN[i - 1][j][k] + ASOUTH_PS_L[k] * Program.WAT_VAPN[i][j - 1][k] +
                                AEAST_PS_L[k] * Program.WAT_VAPN[i + 1][j][k] + ANORTH_PS_L[k] * Program.WAT_VAPN[i][j + 1][k] +
                                AP0_PS_L[k] * WAT_VAP_L[k];

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

                        //Obtain new Water vapour-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            help = WAT_VAPN_L[k];
                            help += (Program.RELAXT * (PIM[k] * WAT_VAPN_L[k + 1] + QIM[k] - help));
                            WAT_VAPN_L[k] = help;

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
                        double[] A_PS_L = Program.A_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        double[] WAT_VAPN_L = Program.WAT_VAPN[i][j];
                        double[] WAT_VAP_L = Program.WAT_VAP[i][j];

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            DIM = AWEST_PS_L[k] * Program.WAT_VAPN[i - 1][j][k] + ASOUTH_PS_L[k] * Program.WAT_VAPN[i][j - 1][k] +
                                AEAST_PS_L[k] * Program.WAT_VAPN[i + 1][j][k] + ANORTH_PS_L[k] * Program.WAT_VAPN[i][j + 1][k] +
                                AP0_PS_L[k] * WAT_VAP_L[k];

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

                        //Obtain new Water vapour-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            help = WAT_VAPN_L[k];
                            help += (Program.RELAXT * (PIM[k] * WAT_VAPN_L[k + 1] + QIM[k] - help));
                            WAT_VAPN_L[k] = help;

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
                        double[] A_PS_L = Program.A_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        double[] WAT_VAPN_L = Program.WAT_VAPN[i][j];
                        double[] WAT_VAP_L = Program.WAT_VAP[i][j];

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            DIM = AWEST_PS_L[k] * Program.WAT_VAPN[i - 1][j][k] + ASOUTH_PS_L[k] * Program.WAT_VAPN[i][j - 1][k] +
                                AEAST_PS_L[k] * Program.WAT_VAPN[i + 1][j][k] + ANORTH_PS_L[k] * Program.WAT_VAPN[i][j + 1][k] +
                                AP0_PS_L[k] * WAT_VAP_L[k];

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

                        //Obtain new Water vapour-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            help = WAT_VAPN_L[k];
                            help += (Program.RELAXT * (PIM[k] * WAT_VAPN_L[k + 1] + QIM[k] - help));
                            WAT_VAPN_L[k] = help;

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
                        double[] A_PS_L = Program.A_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        double[] WAT_VAPN_L = Program.WAT_VAPN[i][j];
                        double[] WAT_VAP_L = Program.WAT_VAP[i][j];

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            DIM = AWEST_PS_L[k] * Program.WAT_VAPN[i - 1][j][k] + ASOUTH_PS_L[k] * Program.WAT_VAPN[i][j - 1][k] +
                                AEAST_PS_L[k] * Program.WAT_VAPN[i + 1][j][k] + ANORTH_PS_L[k] * Program.WAT_VAPN[i][j + 1][k] +
                                AP0_PS_L[k] * WAT_VAP_L[k];

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

                        //Obtain new Water vapour-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            help = WAT_VAPN_L[k];
                            help += (Program.RELAXT * (PIM[k] * WAT_VAPN_L[k + 1] + QIM[k] - help));
                            WAT_VAPN_L[k] = help;

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
                        double[] A_PS_L = Program.A_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        double[] WAT_VAPN_L = Program.WAT_VAPN[i][j];
                        double[] WAT_VAP_L = Program.WAT_VAP[i][j];

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            DIM = AWEST_PS_L[k] * Program.WAT_VAPN[i - 1][j][k] + ASOUTH_PS_L[k] * Program.WAT_VAPN[i][j - 1][k] +
                                AEAST_PS_L[k] * Program.WAT_VAPN[i + 1][j][k] + ANORTH_PS_L[k] * Program.WAT_VAPN[i][j + 1][k] +
                                AP0_PS_L[k] * WAT_VAP_L[k];

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

                        //Obtain new Water vapour-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            help = WAT_VAPN_L[k];
                            help += (Program.RELAXT * (PIM[k] * WAT_VAPN_L[k + 1] + QIM[k] - help));
                            WAT_VAPN_L[k] = help;

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
                        double[] A_PS_L = Program.A_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        double[] WAT_VAPN_L = Program.WAT_VAPN[i][j];
                        double[] WAT_VAP_L = Program.WAT_VAP[i][j];

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            DIM = AWEST_PS_L[k] * Program.WAT_VAPN[i - 1][j][k] + ASOUTH_PS_L[k] * Program.WAT_VAPN[i][j - 1][k] +
                                AEAST_PS_L[k] * Program.WAT_VAPN[i + 1][j][k] + ANORTH_PS_L[k] * Program.WAT_VAPN[i][j + 1][k] +
                                AP0_PS_L[k] * WAT_VAP_L[k];

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

                        //Obtain new Water vapour-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            help = WAT_VAPN_L[k];
                            help += (Program.RELAXT * (PIM[k] * WAT_VAPN_L[k + 1] + QIM[k] - help));
                            WAT_VAPN_L[k] = help;

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
                        double[] A_PS_L = Program.A_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        double[] WAT_VAPN_L = Program.WAT_VAPN[i][j];
                        double[] WAT_VAP_L = Program.WAT_VAP[i][j];

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            DIM = AWEST_PS_L[k] * Program.WAT_VAPN[i - 1][j][k] + ASOUTH_PS_L[k] * Program.WAT_VAPN[i][j - 1][k] +
                                AEAST_PS_L[k] * Program.WAT_VAPN[i + 1][j][k] + ANORTH_PS_L[k] * Program.WAT_VAPN[i][j + 1][k] +
                                AP0_PS_L[k] * WAT_VAP_L[k];

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

                        //Obtain new Water vapour-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            help = WAT_VAPN_L[k];
                            help += (Program.RELAXT * (PIM[k] * WAT_VAPN_L[k + 1] + QIM[k] - help));
                            WAT_VAPN_L[k] = help;

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
                        double[] A_PS_L = Program.A_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        double[] WAT_VAPN_L = Program.WAT_VAPN[i][j];
                        double[] WAT_VAP_L = Program.WAT_VAP[i][j];

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            DIM = AWEST_PS_L[k] * Program.WAT_VAPN[i - 1][j][k] + ASOUTH_PS_L[k] * Program.WAT_VAPN[i][j - 1][k] +
                                AEAST_PS_L[k] * Program.WAT_VAPN[i + 1][j][k] + ANORTH_PS_L[k] * Program.WAT_VAPN[i][j + 1][k] +
                                AP0_PS_L[k] * WAT_VAP_L[k];

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

                        //Obtain new Water vapour-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            help = WAT_VAPN_L[k];
                            help += (Program.RELAXT * (PIM[k] * WAT_VAPN_L[k + 1] + QIM[k] - help));
                            WAT_VAPN_L[k] = help;

                        }
                    }
                }
            });
        }
    }
}
