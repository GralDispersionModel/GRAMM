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
        public static void Timp_calculate(int NI, int NJ, int NK)
        {
            //computation of radiation terms
            Parallel.For(2, NI, Program.pOptions, i =>
            {
                for (int j = 2; j <= NJ - 1; j++)
                {
                    double[] DT_SOL_L = Program.DT_SOL[i][j];
                    double[] DT_TERR_L = Program.DT_TERR[i][j];
                    float[] FAC_L = Program.FAC[i][j];
                    double[] NBZKP_L = Program.NBZKP[i][j];
                    double[] PNBZKP_L = Program.PNBZKP[i][j];
                    double[] RADIATION_L = Program.RADIATION[i][j];
                    float[] RHO_L = Program.RHO[i][j];
                    float[] VOL_J = Program.VOL[i][j];

                    for (int k = 1; k <= NK - 1; k++)
                    {
                        int KKAP = (int)(Math.Floor(NBZKP_L[k]));
                        int KKAM = KKAP - 1;
                        if (KKAP == 1) KKAM = 1;
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

            int range_parallel = (int)(NI / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
            range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
            range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
            //computation of the new temperature
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
                        double[] A_PS_L = Program.A_PS[i][j];
                        float[] AEAST_PS_L = Program.AEAST_PS[i][j];
                        float[] ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        float[] ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        float[] AWEST_PS_L = Program.AWEST_PS[i][j];
                        float[] AP0_PS_L = Program.AP0_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        float[] FAC_L = Program.FAC[i][j];
                        float[] RHO_L = Program.RHO[i][j];
                        double[] T_L = Program.T[i][j];
                        double[] TN_L = Program.TN[i][j];
                        double[] RADIATION_L = RADIATION[i][j];
                        double WQU_AWQ = Program.WQU[i][j] - Program.AWQ[i][j];
                        double[] TNmi_L = Program.TN[i - 1][j]; double[] TNi_L = Program.TN[i + 1][j];
                        double[] TNmj_L = Program.TN[i][j - 1]; double[] TNj_L = Program.TN[i][j + 1];

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            DIM = AWEST_PS_L[k] * TNmi_L[k] + ASOUTH_PS_L[k] * TNmj_L[k] +
                                AEAST_PS_L[k] * TNi_L[k] + ANORTH_PS_L[k] * TNj_L[k] +
                                AP0_PS_L[k] * T_L[k] + RADIATION_L[k];
                            if ((k == 1) && (Program.ICTB == true)) DIM -= WQU_AWQ / Program.CPLUFT * FAC_L[k] * RHO_L[k];

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
                    }
                }
            });


            range_parallel = (int)(NI / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2 + 4);
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
                    for (int j = NJ - 1; j >= 2; j--)
                    {
                        double[] A_PS_L = Program.A_PS[i][j];
                        float[] AEAST_PS_L = Program.AEAST_PS[i][j];
                        float[] ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        float[] ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        float[] AWEST_PS_L = Program.AWEST_PS[i][j];
                        float[] AP0_PS_L = Program.AP0_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        float[] FAC_L = Program.FAC[i][j];
                        float[] RHO_L = Program.RHO[i][j];
                        double[] T_L = Program.T[i][j];
                        double[] TN_L = Program.TN[i][j];
                        double[] RADIATION_L = RADIATION[i][j];
                        double WQU_AWQ = Program.WQU[i][j] - Program.AWQ[i][j];
                        double[] TNmi_L = Program.TN[i - 1][j]; double[] TNi_L = Program.TN[i + 1][j];
                        double[] TNmj_L = Program.TN[i][j - 1]; double[] TNj_L = Program.TN[i][j + 1];

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            DIM = AWEST_PS_L[k] * TNmi_L[k] + ASOUTH_PS_L[k] * TNmj_L[k] +
                                AEAST_PS_L[k] * TNi_L[k] + ANORTH_PS_L[k] * TNj_L[k] +
                                AP0_PS_L[k] * T_L[k] + RADIATION_L[k];
                            if ((k == 1) && (Program.ICTB == true)) DIM -= WQU_AWQ / Program.CPLUFT * FAC_L[k] * RHO_L[k];

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
                    }
                }
            });


            range_parallel = (int)(NI / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2 - 4);
            range_parallel = Math.Max(33 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
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
                        double[] A_PS_L = Program.A_PS[i][j];
                        float[] AEAST_PS_L = Program.AEAST_PS[i][j];
                        float[] ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        float[] ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        float[] AWEST_PS_L = Program.AWEST_PS[i][j];
                        float[] AP0_PS_L = Program.AP0_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        float[] FAC_L = Program.FAC[i][j];
                        float[] RHO_L = Program.RHO[i][j];
                        double[] T_L = Program.T[i][j];
                        double[] TN_L = Program.TN[i][j];
                        double[] RADIATION_L = RADIATION[i][j];
                        double WQU_AWQ = Program.WQU[i][j] - Program.AWQ[i][j];
                        double[] TNmi_L = Program.TN[i - 1][j]; double[] TNi_L = Program.TN[i + 1][j];
                        double[] TNmj_L = Program.TN[i][j - 1]; double[] TNj_L = Program.TN[i][j + 1];

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            DIM = AWEST_PS_L[k] * TNmi_L[k] + ASOUTH_PS_L[k] * TNmj_L[k] +
                                AEAST_PS_L[k] * TNi_L[k] + ANORTH_PS_L[k] * TNj_L[k] +
                                AP0_PS_L[k] * T_L[k] + RADIATION_L[k];
                            if ((k == 1) && (Program.ICTB == true)) DIM -= WQU_AWQ / Program.CPLUFT * FAC_L[k] * RHO_L[k];

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
                    }
                }
            });

            range_parallel = (int)(NI / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
            range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
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
                       double[] A_PS_L = Program.A_PS[i][j];
                       float[] AEAST_PS_L = Program.AEAST_PS[i][j];
                       float[] ANORTH_PS_L = Program.ANORTH_PS[i][j];
                       float[] ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                       float[] AWEST_PS_L = Program.AWEST_PS[i][j];
                       float[] AP0_PS_L = Program.AP0_PS[i][j];
                       double[] B_PS_L = Program.B_PS[i][j];
                       double[] C_PS_L = Program.C_PS[i][j];
                       float[] FAC_L = Program.FAC[i][j];
                       float[] RHO_L = Program.RHO[i][j];
                       double[] T_L = Program.T[i][j];
                       double[] TN_L = Program.TN[i][j];
                       double[] RADIATION_L = RADIATION[i][j];
                       double WQU_AWQ = Program.WQU[i][j] - Program.AWQ[i][j];
                       double[] TNmi_L = Program.TN[i - 1][j]; double[] TNi_L = Program.TN[i + 1][j];
                       double[] TNmj_L = Program.TN[i][j - 1]; double[] TNj_L = Program.TN[i][j + 1];

                       for (int k = 1; k <= NK - 1; k++)
                       {
                           DIM = AWEST_PS_L[k] * TNmi_L[k] + ASOUTH_PS_L[k] * TNmj_L[k] +
                               AEAST_PS_L[k] * TNi_L[k] + ANORTH_PS_L[k] * TNj_L[k] +
                               AP0_PS_L[k] * T_L[k] + RADIATION_L[k];
                           if ((k == 1) && (Program.ICTB == true)) DIM -= WQU_AWQ / Program.CPLUFT * FAC_L[k] * RHO_L[k];

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
                        double[] A_PS_L = Program.A_PS[i][j];
                        float[] AEAST_PS_L = Program.AEAST_PS[i][j];
                        float[] ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        float[] ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        float[] AWEST_PS_L = Program.AWEST_PS[i][j];
                        float[] AP0_PS_L = Program.AP0_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        float[] FAC_L = Program.FAC[i][j];
                        float[] RHO_L = Program.RHO[i][j];
                        double[] T_L = Program.T[i][j];
                        double[] TN_L = Program.TN[i][j];
                        double[] RADIATION_L = RADIATION[i][j];
                        double WQU_AWQ = Program.WQU[i][j] - Program.AWQ[i][j];
                        double[] TNmi_L = Program.TN[i - 1][j]; double[] TNi_L = Program.TN[i + 1][j];
                        double[] TNmj_L = Program.TN[i][j - 1]; double[] TNj_L = Program.TN[i][j + 1];

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            DIM = AWEST_PS_L[k] * TNmi_L[k] + ASOUTH_PS_L[k] * TNmj_L[k] +
                                AEAST_PS_L[k] * TNi_L[k] + ANORTH_PS_L[k] * TNj_L[k] +
                                AP0_PS_L[k] * T_L[k] + RADIATION_L[k];
                            if ((k == 1) && (Program.ICTB == true)) DIM -= WQU_AWQ / Program.CPLUFT * FAC_L[k] * RHO_L[k];

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
                    }
                }
            });

            range_parallel = (int)(NJ / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2 - 4);
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
                    for (int i = NI - 1; i >= 2; i--)
                    {
                        double[] A_PS_L = Program.A_PS[i][j];
                        float[] AEAST_PS_L = Program.AEAST_PS[i][j];
                        float[] ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        float[] ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        float[] AWEST_PS_L = Program.AWEST_PS[i][j];
                        float[] AP0_PS_L = Program.AP0_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        float[] FAC_L = Program.FAC[i][j];
                        float[] RHO_L = Program.RHO[i][j];
                        double[] T_L = Program.T[i][j];
                        double[] TN_L = Program.TN[i][j];
                        double[] RADIATION_L = RADIATION[i][j];
                        double WQU_AWQ = Program.WQU[i][j] - Program.AWQ[i][j];
                        double[] TNmi_L = Program.TN[i - 1][j]; double[] TNi_L = Program.TN[i + 1][j];
                        double[] TNmj_L = Program.TN[i][j - 1]; double[] TNj_L = Program.TN[i][j + 1];

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            DIM = AWEST_PS_L[k] * TNmi_L[k] + ASOUTH_PS_L[k] * TNmj_L[k] +
                                AEAST_PS_L[k] * TNi_L[k] + ANORTH_PS_L[k] * TNj_L[k] +
                                AP0_PS_L[k] * T_L[k] + RADIATION_L[k];
                            if ((k == 1) && (Program.ICTB == true)) DIM -= WQU_AWQ / Program.CPLUFT * FAC_L[k] * RHO_L[k];

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
                    }
                }
            });

            range_parallel = (int)(NJ / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2 + 4);
            range_parallel = Math.Max(33 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
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
                        double[] A_PS_L = Program.A_PS[i][j];
                        float[] AEAST_PS_L = Program.AEAST_PS[i][j];
                        float[] ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        float[] ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        float[] AWEST_PS_L = Program.AWEST_PS[i][j];
                        float[] AP0_PS_L = Program.AP0_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        float[] FAC_L = Program.FAC[i][j];
                        float[] RHO_L = Program.RHO[i][j];
                        double[] T_L = Program.T[i][j];
                        double[] TN_L = Program.TN[i][j];
                        double[] RADIATION_L = RADIATION[i][j];
                        double WQU_AWQ = Program.WQU[i][j] - Program.AWQ[i][j];
                        double[] TNmi_L = Program.TN[i - 1][j]; double[] TNi_L = Program.TN[i + 1][j];
                        double[] TNmj_L = Program.TN[i][j - 1]; double[] TNj_L = Program.TN[i][j + 1];

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            DIM = AWEST_PS_L[k] * TNmi_L[k] + ASOUTH_PS_L[k] * TNmj_L[k] +
                                AEAST_PS_L[k] * TNi_L[k] + ANORTH_PS_L[k] * TNj_L[k] +
                                AP0_PS_L[k] * T_L[k] + RADIATION_L[k];
                            if ((k == 1) && (Program.ICTB == true)) DIM -= WQU_AWQ / Program.CPLUFT * FAC_L[k] * RHO_L[k];

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
                    for (int i = 2; i <= NI - 1; i++)
                    {
                        double[] A_PS_L = Program.A_PS[i][j];
                        float[] AEAST_PS_L = Program.AEAST_PS[i][j];
                        float[] ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        float[] ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        float[] AWEST_PS_L = Program.AWEST_PS[i][j];
                        float[] AP0_PS_L = Program.AP0_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        float[] FAC_L = Program.FAC[i][j];
                        float[] RHO_L = Program.RHO[i][j];
                        double[] T_L = Program.T[i][j];
                        double[] TN_L = Program.TN[i][j];
                        double[] RADIATION_L = RADIATION[i][j];
                        double WQU_AWQ = Program.WQU[i][j] - Program.AWQ[i][j];
                        double[] TNmi_L = Program.TN[i - 1][j]; double[] TNi_L = Program.TN[i + 1][j];
                        double[] TNmj_L = Program.TN[i][j - 1]; double[] TNj_L = Program.TN[i][j + 1];

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            DIM = AWEST_PS_L[k] * TNmi_L[k] + ASOUTH_PS_L[k] * TNmj_L[k] +
                                AEAST_PS_L[k] * TNi_L[k] + ANORTH_PS_L[k] * TNj_L[k] +
                                AP0_PS_L[k] * T_L[k] + RADIATION_L[k];
                            if ((k == 1) && (Program.ICTB == true)) DIM -= WQU_AWQ / Program.CPLUFT * FAC_L[k] * RHO_L[k];

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
                    }
                }
            });

        }
    }
}
