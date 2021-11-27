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
        public static void Tkeimp_calculate(int NI, int NJ, int NK)
        {

            int range_parallel = NI / Program.pOptions.MaxDegreeOfParallelism - (StripeCounter % 6);
            range_parallel = Math.Max(Program.StripeWidth - (StripeCounter % 6), range_parallel); // min. Program.StripeWidth cells per processor
            StripeCounter++;
            range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
                                                           //computation of the new turbulent kinetic energy
                                                           //Parallel.For(2, NI, Program.pOptions, i =>
            Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
            {
                double[] PIM = new double[NK + 1];
                double[] QIM = new double[NK + 1];
                double help;
                for (int i = range.Item1; i < range.Item2; i++)
                {
                    for (int j = 2; j <= NJ - 1; j++)
                    {
                        float[] AEAST_PS_L = Program.AEAST_PS[i][j];
                        float[] ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        float[] ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        float[] AWEST_PS_L = Program.AWEST_PS[i][j];
                        float[] AP0_PS_L = Program.AP0_PS[i][j];
                        double[] A_PS_L = Program.A_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        double[] DISS_L = Program.DISS[i][j];
                        float[] RHO_L = Program.RHOImm[i][j];
                        double[] T_L = Program.T[i][j];
                        double[] TE_L = Program.TE[i][j];
                        double[] TEN_L = Program.TEN[i][j];
                        double[] U_L = Program.U[i][j];
                        double[] V_L = Program.V[i][j];
                        double[] VISV_L = Program.VISV[i][j];
                        float[] VOL_L = Program.VOLImm[i][j];
                        double[] W_L = Program.W[i][j];
                        float[] ZSP_L = Program.ZSPImm[i][j];
                        double ZAX_i = 1 / (Program.ZAX[i] + Program.ZAX[i - 1]);
                        double ZAX_j = 1 / (Program.ZAY[j] + Program.ZAY[j - 1]);
                        double DIM;

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            //Production-terms of turbulent kinetic energy
                            double DUDX = (Program.U[i + 1][j][k] - Program.U[i - 1][j][k]) * ZAX_i;
                            double DVDX = (Program.V[i + 1][j][k] - Program.V[i - 1][j][k]) * ZAX_i;
                            double DWDX = (Program.W[i + 1][j][k] - Program.W[i - 1][j][k]) * ZAX_i;
                            double DUDY = (Program.U[i][j + 1][k] - Program.U[i][j - 1][k]) * ZAX_j;
                            double DVDY = (Program.V[i][j + 1][k] - Program.V[i][j - 1][k]) * ZAX_j;
                            double DWDY = (Program.W[i][j + 1][k] - Program.W[i][j - 1][k]) * ZAX_j;

                            double DUDZ = 0;
                            double DVDZ = 0;
                            double DWDZ = 0;
                            double DTDZ = 0;

                            if (k > 1)
                            {
                                help = 1 / (ZSP_L[k + 1] + ZSP_L[k - 1]);
                                DUDZ = (U_L[k + 1] - U_L[k - 1]) * help;
                                DVDZ = (V_L[k + 1] - V_L[k - 1]) * help;
                                DWDZ = (W_L[k + 1] - W_L[k - 1]) * help;
                                DTDZ = (T_L[k + 1] - T_L[k - 1]) * help - 0.00065;
                            }
                            else
                            {
                                help = 1 / (ZSP_L[k + 1] + ZSP_L[k]);
                                DUDZ = (U_L[k + 1] - U_L[k]) * help;
                                DVDZ = (V_L[k + 1] - V_L[k]) * help;
                                DWDZ = (W_L[k + 1] - W_L[k]) * help;
                                DTDZ = (T_L[k + 1] - T_L[k]) * help - 0.00065;
                            }

                            double PROTE = VISV_L[k] / Program.PRTE * VOL_L[k] *
                                (0.5F * (DUDX * DUDX + DVDY * DVDY + DWDZ * DWDZ) +
                                 (DUDY + DVDX) * (DUDY + DVDX) +
                                 (DUDZ + DWDX) * (DUDZ + DWDX) +
                                 (DVDZ + DWDY) * (DVDZ + DWDY) -
                                 DTDZ * Program.GERD / (T_L[k] + Program.TBZ1) / Program.PRTE) - VOL_L[k] * DISS_L[k] * RHO_L[k];

                            DIM = AWEST_PS_L[k] * Program.TEN[i - 1][j][k] + ASOUTH_PS_L[k] * Program.TEN[i][j - 1][k] +
                                AEAST_PS_L[k] * Program.TEN[i + 1][j][k] + ANORTH_PS_L[k] * Program.TEN[i][j + 1][k] +
                                AP0_PS_L[k] * TE_L[k] + PROTE;

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

                        //Obtain new TKE-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            TEN_L[k] += (Program.RELAXT * (PIM[k] * TEN_L[k + 1] + QIM[k] - TEN_L[k]));
                        }
                    }
                }
            });

            range_parallel = NI / Program.pOptions.MaxDegreeOfParallelism - (StripeCounter % 6);
            range_parallel = Math.Max(Program.StripeWidth - (StripeCounter % 6), range_parallel); // min. Program.StripeWidth cells per processor
            StripeCounter++;
            range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
                                                           //Parallel.For(2, NI, Program.pOptions, i1 =>
            Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
            {
                double[] PIM = new double[NK + 1];
                double[] QIM = new double[NK + 1];
                double help;

                for (int i1 = range.Item1; i1 < range.Item2; i1++)
                {
                    int i = NI - i1 + 1;
                    for (int j = NJ - 1; j >= 2; j--)
                    {
                        float[] AEAST_PS_L = Program.AEAST_PS[i][j];
                        float[] ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        float[] ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        float[] AWEST_PS_L = Program.AWEST_PS[i][j];
                        float[] AP0_PS_L = Program.AP0_PS[i][j];
                        double[] A_PS_L = Program.A_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        double[] DISS_L = Program.DISS[i][j];
                        float[] RHO_L = Program.RHOImm[i][j];
                        double[] T_L = Program.T[i][j];
                        double[] TE_L = Program.TE[i][j];
                        double[] TEN_L = Program.TEN[i][j];
                        double[] U_L = Program.U[i][j];
                        double[] V_L = Program.V[i][j];
                        double[] VISV_L = Program.VISV[i][j];
                        float[] VOL_L = Program.VOLImm[i][j];
                        double[] W_L = Program.W[i][j];
                        float[] ZSP_L = Program.ZSPImm[i][j];
                        double ZAX_i = 1 / (Program.ZAX[i] + Program.ZAX[i - 1]);
                        double ZAX_j = 1 / (Program.ZAY[j] + Program.ZAY[j - 1]);
                        double DIM;

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            //Production-terms of turbulent kinetic energy
                            double DUDX = (Program.U[i + 1][j][k] - Program.U[i - 1][j][k]) * ZAX_i;
                            double DVDX = (Program.V[i + 1][j][k] - Program.V[i - 1][j][k]) * ZAX_i;
                            double DWDX = (Program.W[i + 1][j][k] - Program.W[i - 1][j][k]) * ZAX_i;
                            double DUDY = (Program.U[i][j + 1][k] - Program.U[i][j - 1][k]) * ZAX_j;
                            double DVDY = (Program.V[i][j + 1][k] - Program.V[i][j - 1][k]) * ZAX_j;
                            double DWDY = (Program.W[i][j + 1][k] - Program.W[i][j - 1][k]) * ZAX_j;

                            double DUDZ = 0;
                            double DVDZ = 0;
                            double DWDZ = 0;
                            double DTDZ = 0;

                            if (k > 1)
                            {
                                help = 1 / (ZSP_L[k + 1] + ZSP_L[k - 1]);
                                DUDZ = (U_L[k + 1] - U_L[k - 1]) * help;
                                DVDZ = (V_L[k + 1] - V_L[k - 1]) * help;
                                DWDZ = (W_L[k + 1] - W_L[k - 1]) * help;
                                DTDZ = (T_L[k + 1] - T_L[k - 1]) * help - 0.00065;
                            }
                            else
                            {
                                help = 1 / (ZSP_L[k + 1] + ZSP_L[k]);
                                DUDZ = (U_L[k + 1] - U_L[k]) * help;
                                DVDZ = (V_L[k + 1] - V_L[k]) * help;
                                DWDZ = (W_L[k + 1] - W_L[k]) * help;
                                DTDZ = (T_L[k + 1] - T_L[k]) * help - 0.00065;
                            }

                            double PROTE = VISV_L[k] / Program.PRTE * VOL_L[k] *
                                (0.5F * (DUDX * DUDX + DVDY * DVDY + DWDZ * DWDZ) +
                                 (DUDY + DVDX) * (DUDY + DVDX) +
                                 (DUDZ + DWDX) * (DUDZ + DWDX) +
                                 (DVDZ + DWDY) * (DVDZ + DWDY) -
                                 DTDZ * Program.GERD / (T_L[k] + Program.TBZ1) / Program.PRTE) - VOL_L[k] * DISS_L[k] * RHO_L[k];

                            DIM = AWEST_PS_L[k] * Program.TEN[i - 1][j][k] + ASOUTH_PS_L[k] * Program.TEN[i][j - 1][k] +
                                AEAST_PS_L[k] * Program.TEN[i + 1][j][k] + ANORTH_PS_L[k] * Program.TEN[i][j + 1][k] +
                                AP0_PS_L[k] * TE_L[k] + PROTE;

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

                        //Obtain new TKE-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            TEN_L[k] += (Program.RELAXT * (PIM[k] * TEN_L[k + 1] + QIM[k] - TEN_L[k]));
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
                double[] PIM = new double[NK + 1];
                double[] QIM = new double[NK + 1];
                double help;

                for (int i1 = range.Item1; i1 < range.Item2; i1++)
                {
                    int i = NI - i1 + 1;
                    for (int j = 2; j <= NJ - 1; j++)
                    {
                        float[] AEAST_PS_L = Program.AEAST_PS[i][j];
                        float[] ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        float[] ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        float[] AWEST_PS_L = Program.AWEST_PS[i][j];
                        float[] AP0_PS_L = Program.AP0_PS[i][j];
                        double[] A_PS_L = Program.A_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        double[] DISS_L = Program.DISS[i][j];
                        float[] RHO_L = Program.RHOImm[i][j];
                        double[] T_L = Program.T[i][j];
                        double[] TE_L = Program.TE[i][j];
                        double[] TEN_L = Program.TEN[i][j];
                        double[] U_L = Program.U[i][j];
                        double[] V_L = Program.V[i][j];
                        double[] VISV_L = Program.VISV[i][j];
                        float[] VOL_L = Program.VOLImm[i][j];
                        double[] W_L = Program.W[i][j];
                        float[] ZSP_L = Program.ZSPImm[i][j];
                        double ZAX_i = 1 / (Program.ZAX[i] + Program.ZAX[i - 1]);
                        double ZAX_j = 1 / (Program.ZAY[j] + Program.ZAY[j - 1]);
                        double DIM;

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            //Production-terms of turbulent kinetic energy
                            double DUDX = (Program.U[i + 1][j][k] - Program.U[i - 1][j][k]) * ZAX_i;
                            double DVDX = (Program.V[i + 1][j][k] - Program.V[i - 1][j][k]) * ZAX_i;
                            double DWDX = (Program.W[i + 1][j][k] - Program.W[i - 1][j][k]) * ZAX_i;
                            double DUDY = (Program.U[i][j + 1][k] - Program.U[i][j - 1][k]) * ZAX_j;
                            double DVDY = (Program.V[i][j + 1][k] - Program.V[i][j - 1][k]) * ZAX_j;
                            double DWDY = (Program.W[i][j + 1][k] - Program.W[i][j - 1][k]) * ZAX_j;

                            double DUDZ = 0;
                            double DVDZ = 0;
                            double DWDZ = 0;
                            double DTDZ = 0;

                            if (k > 1)
                            {
                                help = 1 / (ZSP_L[k + 1] + ZSP_L[k - 1]);
                                DUDZ = (U_L[k + 1] - U_L[k - 1]) * help;
                                DVDZ = (V_L[k + 1] - V_L[k - 1]) * help;
                                DWDZ = (W_L[k + 1] - W_L[k - 1]) * help;
                                DTDZ = (T_L[k + 1] - T_L[k - 1]) * help - 0.00065;
                            }
                            else
                            {
                                help = 1 / (ZSP_L[k + 1] + ZSP_L[k]);
                                DUDZ = (U_L[k + 1] - U_L[k]) * help;
                                DVDZ = (V_L[k + 1] - V_L[k]) * help;
                                DWDZ = (W_L[k + 1] - W_L[k]) * help;
                                DTDZ = (T_L[k + 1] - T_L[k]) * help - 0.00065;
                            }

                            double PROTE = VISV_L[k] / Program.PRTE * VOL_L[k] *
                                (0.5F * (DUDX * DUDX + DVDY * DVDY + DWDZ * DWDZ) +
                                 (DUDY + DVDX) * (DUDY + DVDX) +
                                 (DUDZ + DWDX) * (DUDZ + DWDX) +
                                 (DVDZ + DWDY) * (DVDZ + DWDY) -
                                 DTDZ * Program.GERD / (T_L[k] + Program.TBZ1) / Program.PRTE) - VOL_L[k] * DISS_L[k] * RHO_L[k];

                            DIM = AWEST_PS_L[k] * Program.TEN[i - 1][j][k] + ASOUTH_PS_L[k] * Program.TEN[i][j - 1][k] +
                                AEAST_PS_L[k] * Program.TEN[i + 1][j][k] + ANORTH_PS_L[k] * Program.TEN[i][j + 1][k] +
                                AP0_PS_L[k] * TE_L[k] + PROTE;

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

                        //Obtain new TKE-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            TEN_L[k] += (Program.RELAXT * (PIM[k] * TEN_L[k + 1] + QIM[k] - TEN_L[k]));
                        }
                    }
                }
            });

            range_parallel = (int)(NI / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2 + 3);
            range_parallel = Math.Max(33 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
            range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
                                                           // Parallel.For(2, NI, Program.pOptions, i =>
            Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
            {
                double[] PIM = new double[NK + 1];
                double[] QIM = new double[NK + 1];
                double help;

                for (int i = range.Item1; i < range.Item2; i++)
                {
                    for (int j = NJ - 1; j >= 2; j--)
                    {
                        float[] AEAST_PS_L = Program.AEAST_PS[i][j];
                        float[] ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        float[] ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        float[] AWEST_PS_L = Program.AWEST_PS[i][j];
                        float[] AP0_PS_L = Program.AP0_PS[i][j];
                        double[] A_PS_L = Program.A_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        double[] DISS_L = Program.DISS[i][j];
                        float[] RHO_L = Program.RHOImm[i][j];
                        double[] T_L = Program.T[i][j];
                        double[] TE_L = Program.TE[i][j];
                        double[] TEN_L = Program.TEN[i][j];
                        double[] U_L = Program.U[i][j];
                        double[] V_L = Program.V[i][j];
                        double[] VISV_L = Program.VISV[i][j];
                        float[] VOL_L = Program.VOLImm[i][j];
                        double[] W_L = Program.W[i][j];
                        float[] ZSP_L = Program.ZSPImm[i][j];
                        double ZAX_i = 1 / (Program.ZAX[i] + Program.ZAX[i - 1]);
                        double ZAX_j = 1 / (Program.ZAY[j] + Program.ZAY[j - 1]);
                        double DIM;

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            //Production-terms of turbulent kinetic energy
                            double DUDX = (Program.U[i + 1][j][k] - Program.U[i - 1][j][k]) * ZAX_i;
                            double DVDX = (Program.V[i + 1][j][k] - Program.V[i - 1][j][k]) * ZAX_i;
                            double DWDX = (Program.W[i + 1][j][k] - Program.W[i - 1][j][k]) * ZAX_i;
                            double DUDY = (Program.U[i][j + 1][k] - Program.U[i][j - 1][k]) * ZAX_j;
                            double DVDY = (Program.V[i][j + 1][k] - Program.V[i][j - 1][k]) * ZAX_j;
                            double DWDY = (Program.W[i][j + 1][k] - Program.W[i][j - 1][k]) * ZAX_j;

                            double DUDZ = 0;
                            double DVDZ = 0;
                            double DWDZ = 0;
                            double DTDZ = 0;

                            if (k > 1)
                            {
                                help = 1 / (ZSP_L[k + 1] + ZSP_L[k - 1]);
                                DUDZ = (U_L[k + 1] - U_L[k - 1]) * help;
                                DVDZ = (V_L[k + 1] - V_L[k - 1]) * help;
                                DWDZ = (W_L[k + 1] - W_L[k - 1]) * help;
                                DTDZ = (T_L[k + 1] - T_L[k - 1]) * help - 0.00065;
                            }
                            else
                            {
                                help = 1 / (ZSP_L[k + 1] + ZSP_L[k]);
                                DUDZ = (U_L[k + 1] - U_L[k]) * help;
                                DVDZ = (V_L[k + 1] - V_L[k]) * help;
                                DWDZ = (W_L[k + 1] - W_L[k]) * help;
                                DTDZ = (T_L[k + 1] - T_L[k]) * help - 0.00065;
                            }

                            double PROTE = VISV_L[k] / Program.PRTE * VOL_L[k] *
                                (0.5F * (DUDX * DUDX + DVDY * DVDY + DWDZ * DWDZ) +
                                 (DUDY + DVDX) * (DUDY + DVDX) +
                                 (DUDZ + DWDX) * (DUDZ + DWDX) +
                                 (DVDZ + DWDY) * (DVDZ + DWDY) -
                                 DTDZ * Program.GERD / (T_L[k] + Program.TBZ1) / Program.PRTE) - VOL_L[k] * DISS_L[k] * RHO_L[k];

                            DIM = AWEST_PS_L[k] * Program.TEN[i - 1][j][k] + ASOUTH_PS_L[k] * Program.TEN[i][j - 1][k] +
                                AEAST_PS_L[k] * Program.TEN[i + 1][j][k] + ANORTH_PS_L[k] * Program.TEN[i][j + 1][k] +
                                AP0_PS_L[k] * TE_L[k] + PROTE;

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

                        //Obtain new TKE-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            TEN_L[k] += (Program.RELAXT * (PIM[k] * TEN_L[k + 1] + QIM[k] - TEN_L[k]));
                        }
                    }
                }
            });

            range_parallel = NJ / Program.pOptions.MaxDegreeOfParallelism - (StripeCounter % 6);
            range_parallel = Math.Max(Program.StripeWidth - (StripeCounter % 6), range_parallel); // min. Program.StripeWidth cells per processor
            StripeCounter++;
            range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
                                                           //Parallel.For(2, NJ, Program.pOptions, j =>
            Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
            {
                double[] PIM = new double[NK + 1];
                double[] QIM = new double[NK + 1];
                double help;

                for (int j = range.Item1; j < range.Item2; j++)
                {
                    for (int i = 2; i <= NI - 1; i++)
                    {
                        float[] AEAST_PS_L = Program.AEAST_PS[i][j];
                        float[] ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        float[] ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        float[] AWEST_PS_L = Program.AWEST_PS[i][j];
                        float[] AP0_PS_L = Program.AP0_PS[i][j];
                        double[] A_PS_L = Program.A_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        double[] DISS_L = Program.DISS[i][j];
                        float[] RHO_L = Program.RHOImm[i][j];
                        double[] T_L = Program.T[i][j];
                        double[] TE_L = Program.TE[i][j];
                        double[] TEN_L = Program.TEN[i][j];
                        double[] U_L = Program.U[i][j];
                        double[] V_L = Program.V[i][j];
                        double[] VISV_L = Program.VISV[i][j];
                        float[] VOL_L = Program.VOLImm[i][j];
                        double[] W_L = Program.W[i][j];
                        float[] ZSP_L = Program.ZSPImm[i][j];
                        double ZAX_i = 1 / (Program.ZAX[i] + Program.ZAX[i - 1]);
                        double ZAX_j = 1 / (Program.ZAY[j] + Program.ZAY[j - 1]);
                        double DIM;

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            //Production-terms of turbulent kinetic energy
                            double DUDX = (Program.U[i + 1][j][k] - Program.U[i - 1][j][k]) * ZAX_i;
                            double DVDX = (Program.V[i + 1][j][k] - Program.V[i - 1][j][k]) * ZAX_i;
                            double DWDX = (Program.W[i + 1][j][k] - Program.W[i - 1][j][k]) * ZAX_i;
                            double DUDY = (Program.U[i][j + 1][k] - Program.U[i][j - 1][k]) * ZAX_j;
                            double DVDY = (Program.V[i][j + 1][k] - Program.V[i][j - 1][k]) * ZAX_j;
                            double DWDY = (Program.W[i][j + 1][k] - Program.W[i][j - 1][k]) * ZAX_j;

                            double DUDZ = 0;
                            double DVDZ = 0;
                            double DWDZ = 0;
                            double DTDZ = 0;

                            if (k > 1)
                            {
                                help = 1 / (ZSP_L[k + 1] + ZSP_L[k - 1]);
                                DUDZ = (U_L[k + 1] - U_L[k - 1]) * help;
                                DVDZ = (V_L[k + 1] - V_L[k - 1]) * help;
                                DWDZ = (W_L[k + 1] - W_L[k - 1]) * help;
                                DTDZ = (T_L[k + 1] - T_L[k - 1]) * help - 0.00065;
                            }
                            else
                            {
                                help = 1 / (ZSP_L[k + 1] + ZSP_L[k]);
                                DUDZ = (U_L[k + 1] - U_L[k]) * help;
                                DVDZ = (V_L[k + 1] - V_L[k]) * help;
                                DWDZ = (W_L[k + 1] - W_L[k]) * help;
                                DTDZ = (T_L[k + 1] - T_L[k]) * help - 0.00065;
                            }

                            double PROTE = VISV_L[k] / Program.PRTE * VOL_L[k] *
                                (0.5F * (DUDX * DUDX + DVDY * DVDY + DWDZ * DWDZ) +
                                 (DUDY + DVDX) * (DUDY + DVDX) +
                                 (DUDZ + DWDX) * (DUDZ + DWDX) +
                                 (DVDZ + DWDY) * (DVDZ + DWDY) -
                                 DTDZ * Program.GERD / (T_L[k] + Program.TBZ1) / Program.PRTE) - VOL_L[k] * DISS_L[k] * RHO_L[k];

                            DIM = AWEST_PS_L[k] * Program.TEN[i - 1][j][k] + ASOUTH_PS_L[k] * Program.TEN[i][j - 1][k] +
                                AEAST_PS_L[k] * Program.TEN[i + 1][j][k] + ANORTH_PS_L[k] * Program.TEN[i][j + 1][k] +
                                AP0_PS_L[k] * TE_L[k] + PROTE;

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

                        //Obtain new TKE-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            TEN_L[k] += (Program.RELAXT * (PIM[k] * TEN_L[k + 1] + QIM[k] - TEN_L[k]));
                        }
                    }
                }
            });

            range_parallel = NJ / Program.pOptions.MaxDegreeOfParallelism - (StripeCounter % 6);
            range_parallel = Math.Max(Program.StripeWidth - (StripeCounter % 6), range_parallel); // min. Program.StripeWidth cells per processor
            StripeCounter++;
            range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
                                                           //Parallel.For(2, NJ, Program.pOptions, j1 =>
            Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
            {
                double[] PIM = new double[NK + 1];
                double[] QIM = new double[NK + 1];
                double help;

                for (int j1 = range.Item1; j1 < range.Item2; j1++)
                {
                    int j = NJ - j1 + 1;
                    for (int i = NI - 1; i >= 2; i--)
                    {
                        float[] AEAST_PS_L = Program.AEAST_PS[i][j];
                        float[] ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        float[] ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        float[] AWEST_PS_L = Program.AWEST_PS[i][j];
                        float[] AP0_PS_L = Program.AP0_PS[i][j];
                        double[] A_PS_L = Program.A_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        double[] DISS_L = Program.DISS[i][j];
                        float[] RHO_L = Program.RHOImm[i][j];
                        double[] T_L = Program.T[i][j];
                        double[] TE_L = Program.TE[i][j];
                        double[] TEN_L = Program.TEN[i][j];
                        double[] U_L = Program.U[i][j];
                        double[] V_L = Program.V[i][j];
                        double[] VISV_L = Program.VISV[i][j];
                        float[] VOL_L = Program.VOLImm[i][j];
                        double[] W_L = Program.W[i][j];
                        float[] ZSP_L = Program.ZSPImm[i][j];
                        double ZAX_i = 1 / (Program.ZAX[i] + Program.ZAX[i - 1]);
                        double ZAX_j = 1 / (Program.ZAY[j] + Program.ZAY[j - 1]);
                        double DIM;

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            //Production-terms of turbulent kinetic energy
                            double DUDX = (Program.U[i + 1][j][k] - Program.U[i - 1][j][k]) * ZAX_i;
                            double DVDX = (Program.V[i + 1][j][k] - Program.V[i - 1][j][k]) * ZAX_i;
                            double DWDX = (Program.W[i + 1][j][k] - Program.W[i - 1][j][k]) * ZAX_i;
                            double DUDY = (Program.U[i][j + 1][k] - Program.U[i][j - 1][k]) * ZAX_j;
                            double DVDY = (Program.V[i][j + 1][k] - Program.V[i][j - 1][k]) * ZAX_j;
                            double DWDY = (Program.W[i][j + 1][k] - Program.W[i][j - 1][k]) * ZAX_j;

                            double DUDZ = 0;
                            double DVDZ = 0;
                            double DWDZ = 0;
                            double DTDZ = 0;

                            if (k > 1)
                            {
                                help = 1 / (ZSP_L[k + 1] + ZSP_L[k - 1]);
                                DUDZ = (U_L[k + 1] - U_L[k - 1]) * help;
                                DVDZ = (V_L[k + 1] - V_L[k - 1]) * help;
                                DWDZ = (W_L[k + 1] - W_L[k - 1]) * help;
                                DTDZ = (T_L[k + 1] - T_L[k - 1]) * help - 0.00065;
                            }
                            else
                            {
                                help = 1 / (ZSP_L[k + 1] + ZSP_L[k]);
                                DUDZ = (U_L[k + 1] - U_L[k]) * help;
                                DVDZ = (V_L[k + 1] - V_L[k]) * help;
                                DWDZ = (W_L[k + 1] - W_L[k]) * help;
                                DTDZ = (T_L[k + 1] - T_L[k]) * help - 0.00065;
                            }

                            double PROTE = VISV_L[k] / Program.PRTE * VOL_L[k] *
                                (0.5F * (DUDX * DUDX + DVDY * DVDY + DWDZ * DWDZ) +
                                 (DUDY + DVDX) * (DUDY + DVDX) +
                                 (DUDZ + DWDX) * (DUDZ + DWDX) +
                                 (DVDZ + DWDY) * (DVDZ + DWDY) -
                                 DTDZ * Program.GERD / (T_L[k] + Program.TBZ1) / Program.PRTE) - VOL_L[k] * DISS_L[k] * RHO_L[k];

                            DIM = AWEST_PS_L[k] * Program.TEN[i - 1][j][k] + ASOUTH_PS_L[k] * Program.TEN[i][j - 1][k] +
                                AEAST_PS_L[k] * Program.TEN[i + 1][j][k] + ANORTH_PS_L[k] * Program.TEN[i][j + 1][k] +
                                AP0_PS_L[k] * TE_L[k] + PROTE;

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

                        //Obtain new TKE-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            TEN_L[k] += (Program.RELAXT * (PIM[k] * TEN_L[k + 1] + QIM[k] - TEN_L[k]));
                        }
                    }
                }
            });

            range_parallel = (int)(NJ / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2 + 6);
            range_parallel = Math.Max(36 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
            range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
                                                           //Parallel.For(2, NJ, Program.pOptions, j =>
            Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
            {
                double[] PIM = new double[NK + 1];
                double[] QIM = new double[NK + 1];
                double help;

                for (int j = range.Item1; j < range.Item2; j++)
                {
                    for (int i = NI - 1; i >= 2; i--)
                    {
                        float[] AEAST_PS_L = Program.AEAST_PS[i][j];
                        float[] ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        float[] ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        float[] AWEST_PS_L = Program.AWEST_PS[i][j];
                        float[] AP0_PS_L = Program.AP0_PS[i][j];
                        double[] A_PS_L = Program.A_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        double[] DISS_L = Program.DISS[i][j];
                        float[] RHO_L = Program.RHOImm[i][j];
                        double[] T_L = Program.T[i][j];
                        double[] TE_L = Program.TE[i][j];
                        double[] TEN_L = Program.TEN[i][j];
                        double[] U_L = Program.U[i][j];
                        double[] V_L = Program.V[i][j];
                        double[] VISV_L = Program.VISV[i][j];
                        float[] VOL_L = Program.VOLImm[i][j];
                        double[] W_L = Program.W[i][j];
                        float[] ZSP_L = Program.ZSPImm[i][j];
                        double ZAX_i = 1 / (Program.ZAX[i] + Program.ZAX[i - 1]);
                        double ZAX_j = 1 / (Program.ZAY[j] + Program.ZAY[j - 1]);
                        double DIM;

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            //Production-terms of turbulent kinetic energy
                            double DUDX = (Program.U[i + 1][j][k] - Program.U[i - 1][j][k]) * ZAX_i;
                            double DVDX = (Program.V[i + 1][j][k] - Program.V[i - 1][j][k]) * ZAX_i;
                            double DWDX = (Program.W[i + 1][j][k] - Program.W[i - 1][j][k]) * ZAX_i;
                            double DUDY = (Program.U[i][j + 1][k] - Program.U[i][j - 1][k]) * ZAX_j;
                            double DVDY = (Program.V[i][j + 1][k] - Program.V[i][j - 1][k]) * ZAX_j;
                            double DWDY = (Program.W[i][j + 1][k] - Program.W[i][j - 1][k]) * ZAX_j;

                            double DUDZ = 0;
                            double DVDZ = 0;
                            double DWDZ = 0;
                            double DTDZ = 0;

                            if (k > 1)
                            {
                                help = 1 / (ZSP_L[k + 1] + ZSP_L[k - 1]);
                                DUDZ = (U_L[k + 1] - U_L[k - 1]) * help;
                                DVDZ = (V_L[k + 1] - V_L[k - 1]) * help;
                                DWDZ = (W_L[k + 1] - W_L[k - 1]) * help;
                                DTDZ = (T_L[k + 1] - T_L[k - 1]) * help - 0.00065;
                            }
                            else
                            {
                                help = 1 / (ZSP_L[k + 1] + ZSP_L[k]);
                                DUDZ = (U_L[k + 1] - U_L[k]) * help;
                                DVDZ = (V_L[k + 1] - V_L[k]) * help;
                                DWDZ = (W_L[k + 1] - W_L[k]) * help;
                                DTDZ = (T_L[k + 1] - T_L[k]) * help - 0.00065;
                            }

                            double PROTE = VISV_L[k] / Program.PRTE * VOL_L[k] *
                                (0.5F * (DUDX * DUDX + DVDY * DVDY + DWDZ * DWDZ) +
                                 (DUDY + DVDX) * (DUDY + DVDX) +
                                 (DUDZ + DWDX) * (DUDZ + DWDX) +
                                 (DVDZ + DWDY) * (DVDZ + DWDY) -
                                 DTDZ * Program.GERD / (T_L[k] + Program.TBZ1) / Program.PRTE) - VOL_L[k] * DISS_L[k] * RHO_L[k];

                            DIM = AWEST_PS_L[k] * Program.TEN[i - 1][j][k] + ASOUTH_PS_L[k] * Program.TEN[i][j - 1][k] +
                                AEAST_PS_L[k] * Program.TEN[i + 1][j][k] + ANORTH_PS_L[k] * Program.TEN[i][j + 1][k] +
                                AP0_PS_L[k] * TE_L[k] + PROTE;

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

                        //Obtain new TKE-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            TEN_L[k] += (Program.RELAXT * (PIM[k] * TEN_L[k + 1] + QIM[k] - TEN_L[k]));
                        }
                    }
                }
            });

            range_parallel = (int)(NJ / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2 + 3);
            range_parallel = Math.Max(33 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
            range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
                                                           //Parallel.For(2, NJ, Program.pOptions, j1 =>
            Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
            {
                double[] PIM = new double[NK + 1];
                double[] QIM = new double[NK + 1];
                double help;

                for (int j1 = range.Item1; j1 < range.Item2; j1++)
                {
                    int j = NJ - j1 + 1;
                    for (int i = 2; i <= NI - 1; i++)
                    {
                        float[] AEAST_PS_L = Program.AEAST_PS[i][j];
                        float[] ANORTH_PS_L = Program.ANORTH_PS[i][j];
                        float[] ASOUTH_PS_L = Program.ASOUTH_PS[i][j];
                        float[] AWEST_PS_L = Program.AWEST_PS[i][j];
                        float[] AP0_PS_L = Program.AP0_PS[i][j];
                        double[] A_PS_L = Program.A_PS[i][j];
                        double[] B_PS_L = Program.B_PS[i][j];
                        double[] C_PS_L = Program.C_PS[i][j];
                        double[] DISS_L = Program.DISS[i][j];
                        float[] RHO_L = Program.RHOImm[i][j];
                        double[] T_L = Program.T[i][j];
                        double[] TE_L = Program.TE[i][j];
                        double[] TEN_L = Program.TEN[i][j];
                        double[] U_L = Program.U[i][j];
                        double[] V_L = Program.V[i][j];
                        double[] VISV_L = Program.VISV[i][j];
                        float[] VOL_L = Program.VOLImm[i][j];
                        double[] W_L = Program.W[i][j];
                        float[] ZSP_L = Program.ZSPImm[i][j];
                        double ZAX_i = 1 / (Program.ZAX[i] + Program.ZAX[i - 1]);
                        double ZAX_j = 1 / (Program.ZAY[j] + Program.ZAY[j - 1]);
                        double DIM;

                        for (int k = 1; k <= NK - 1; k++)
                        {
                            //Production-terms of turbulent kinetic energy
                            double DUDX = (Program.U[i + 1][j][k] - Program.U[i - 1][j][k]) * ZAX_i;
                            double DVDX = (Program.V[i + 1][j][k] - Program.V[i - 1][j][k]) * ZAX_i;
                            double DWDX = (Program.W[i + 1][j][k] - Program.W[i - 1][j][k]) * ZAX_i;
                            double DUDY = (Program.U[i][j + 1][k] - Program.U[i][j - 1][k]) * ZAX_j;
                            double DVDY = (Program.V[i][j + 1][k] - Program.V[i][j - 1][k]) * ZAX_j;
                            double DWDY = (Program.W[i][j + 1][k] - Program.W[i][j - 1][k]) * ZAX_j;

                            double DUDZ = 0;
                            double DVDZ = 0;
                            double DWDZ = 0;
                            double DTDZ = 0;

                            if (k > 1)
                            {
                                help = 1 / (ZSP_L[k + 1] + ZSP_L[k - 1]);
                                DUDZ = (U_L[k + 1] - U_L[k - 1]) * help;
                                DVDZ = (V_L[k + 1] - V_L[k - 1]) * help;
                                DWDZ = (W_L[k + 1] - W_L[k - 1]) * help;
                                DTDZ = (T_L[k + 1] - T_L[k - 1]) * help - 0.00065;
                            }
                            else
                            {
                                help = 1 / (ZSP_L[k + 1] + ZSP_L[k]);
                                DUDZ = (U_L[k + 1] - U_L[k]) * help;
                                DVDZ = (V_L[k + 1] - V_L[k]) * help;
                                DWDZ = (W_L[k + 1] - W_L[k]) * help;
                                DTDZ = (T_L[k + 1] - T_L[k]) * help - 0.00065;
                            }

                            double PROTE = VISV_L[k] / Program.PRTE * VOL_L[k] *
                                (0.5F * (DUDX * DUDX + DVDY * DVDY + DWDZ * DWDZ) +
                                 (DUDY + DVDX) * (DUDY + DVDX) +
                                 (DUDZ + DWDX) * (DUDZ + DWDX) +
                                 (DVDZ + DWDY) * (DVDZ + DWDY) -
                                 DTDZ * Program.GERD / (T_L[k] + Program.TBZ1) / Program.PRTE) - VOL_L[k] * DISS_L[k] * RHO_L[k];

                            DIM = AWEST_PS_L[k] * Program.TEN[i - 1][j][k] + ASOUTH_PS_L[k] * Program.TEN[i][j - 1][k] +
                                AEAST_PS_L[k] * Program.TEN[i + 1][j][k] + ANORTH_PS_L[k] * Program.TEN[i][j + 1][k] +
                                AP0_PS_L[k] * TE_L[k] + PROTE;

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

                        //Obtain new TKE-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            TEN_L[k] += (Program.RELAXT * (PIM[k] * TEN_L[k + 1] + QIM[k] - TEN_L[k]));
                        }
                    }
                }
            });
        }
    }
}
