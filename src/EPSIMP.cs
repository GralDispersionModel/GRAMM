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
        public static void Epsimp_calculate(int NI, int NJ, int NK)
        {
            int range_parallel = NI / Program.pOptions.MaxDegreeOfParallelism - (StripeCounter % 6);
            range_parallel = Math.Max(Program.StripeWidth - (StripeCounter % 6), range_parallel); // min. Program.StripeWidth cells per processor
            StripeCounter++;
            range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
                                                           //computation of the new turbulent kinetic energy
                                                           // Parallel.For(2, NI, Program.pOptions, i =>
            Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
            {
                Span<double> PIM = stackalloc double[NK + 1];
                Span<double> QIM = stackalloc double[NK + 1];
                double help;

                for (int i = range.Item1; i < range.Item2; i++)
                {
                    for (int j = 2; j <= NJ - 1; j++)
                    {
                        for (int k = 1; k <= NK - 1; k++)
                        {
                            //Production-terms of turbulent kinetic energy
                            double DUDX = (Program.U[i + 1][j][k] - Program.U[i - 1][j][k]) / (Program.ZAX[i] + Program.ZAX[i - 1]);
                            double DVDX = (Program.V[i + 1][j][k] - Program.V[i - 1][j][k]) / (Program.ZAX[i] + Program.ZAX[i - 1]);
                            double DWDX = (Program.W[i + 1][j][k] - Program.W[i - 1][j][k]) / (Program.ZAX[i] + Program.ZAX[i - 1]);
                            double DUDY = (Program.U[i][j + 1][k] - Program.U[i][j - 1][k]) / (Program.ZAY[j] + Program.ZAY[j - 1]);
                            double DVDY = (Program.V[i][j + 1][k] - Program.V[i][j - 1][k]) / (Program.ZAY[j] + Program.ZAY[j - 1]);
                            double DWDY = (Program.W[i][j + 1][k] - Program.W[i][j - 1][k]) / (Program.ZAY[j] + Program.ZAY[j - 1]);

                            double DUDZ = 0;
                            double DVDZ = 0;
                            double DWDZ = 0;
                            double DTDZ = 0;

                            if (k > 1)
                            {
                                DUDZ = (Program.U[i][j][k + 1] - Program.U[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]);
                                DVDZ = (Program.V[i][j][k + 1] - Program.V[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]);
                                DWDZ = (Program.W[i][j][k + 1] - Program.W[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]);
                                DTDZ = (Program.T[i][j][k + 1] - Program.T[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]) - 0.00065;
                            }
                            else
                            {
                                DUDZ = (Program.U[i][j][k + 1] - Program.U[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]);
                                DVDZ = (Program.V[i][j][k + 1] - Program.V[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]);
                                DWDZ = (Program.W[i][j][k + 1] - Program.W[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]);
                                DTDZ = (Program.T[i][j][k + 1] - Program.T[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]) - 0.00065;
                            }

                            //production-terms for the dissipation rate
                            double PSTRESS = Program.VISV[i][j][k] *
                                (0.5F * (DUDX * DUDX + DVDY * DVDY + DWDZ * DWDZ) +
                                 (DUDY + DVDX) * (DUDY + DVDX) +
                                 (DUDZ + DWDX) * (DUDZ + DWDX) +
                                 (DVDZ + DWDY) * (DVDZ + DWDY)) * Program.RHO[i][j][k];
                            double PBUOY = -Program.VISV[i][j][k] * DTDZ * Program.GERD / (Program.T[i][j][k] + Program.TBZ1) / Program.PRTE * Program.RHO[i][j][k];

                            double DIM = Program.AWEST_PS[i][j][k] * Program.DISSN[i - 1][j][k] + Program.ASOUTH_PS[i][j][k] * Program.DISSN[i][j - 1][k] +
                                Program.AEAST_PS[i][j][k] * Program.DISSN[i + 1][j][k] + Program.ANORTH_PS[i][j][k] * Program.DISSN[i][j + 1][k] +
                                Program.AP0_PS[i][j][k] * Program.DISS[i][j][k] +
                                (Program.DISS[i][j][k] / Math.Max(Math.Abs(Program.TE[i][j][k]), 0.01) * (1.44 * PSTRESS + 0.4 * PBUOY - 1.92 * Program.DISS[i][j][k] * Program.RHO[i][j][k])) * Program.VOLImm[i][j][k];

                            //Recurrence formula
                            if (k > 1)
                            {
                                help = 1 / (Program.A_PS[i][j][k] - Program.C_PS[i][j][k] * PIM[k - 1]);
                                PIM[k] = Program.B_PS[i][j][k] * help;
                                QIM[k] = (DIM + Program.C_PS[i][j][k] * QIM[k - 1]) * help;
                            }
                            else
                            {
                                help = 1 / Program.A_PS[i][j][k];
                                PIM[k] = Program.B_PS[i][j][k] * help;
                                QIM[k] = DIM * help;
                            }
                        }

                        //Obtain new TKE-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            help = Program.DISSN[i][j][k];
                            help += (Program.RELAXT * (PIM[k] * Program.DISSN[i][j][k + 1] + QIM[k] - help));
                            Program.DISSN[i][j][k] = help;

                        }
                    }
                }
            });

            range_parallel = NI / Program.pOptions.MaxDegreeOfParallelism - (StripeCounter % 6);
            range_parallel = Math.Max(Program.StripeWidth - (StripeCounter % 6), range_parallel); // min. Program.StripeWidth cells per processor
            StripeCounter++;
            range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
                                                           // Parallel.For(2, NI, Program.pOptions, i1 =>
            Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
            {
                Span<double> PIM = stackalloc double[NK + 1];
                Span<double> QIM = stackalloc double[NK + 1];
                double help;

                for (int i1 = range.Item1; i1 < range.Item2; i1++)
                {
                    int i = NI - i1 + 1;
                    for (int j = NJ - 1; j >= 2; j--)
                    {
                        for (int k = 1; k <= NK - 1; k++)
                        {
                            //Production-terms of turbulent kinetic energy
                            double DUDX = (Program.U[i + 1][j][k] - Program.U[i - 1][j][k]) / (Program.ZAX[i] + Program.ZAX[i - 1]);
                            double DVDX = (Program.V[i + 1][j][k] - Program.V[i - 1][j][k]) / (Program.ZAX[i] + Program.ZAX[i - 1]);
                            double DWDX = (Program.W[i + 1][j][k] - Program.W[i - 1][j][k]) / (Program.ZAX[i] + Program.ZAX[i - 1]);
                            double DUDY = (Program.U[i][j + 1][k] - Program.U[i][j - 1][k]) / (Program.ZAY[j] + Program.ZAY[j - 1]);
                            double DVDY = (Program.V[i][j + 1][k] - Program.V[i][j - 1][k]) / (Program.ZAY[j] + Program.ZAY[j - 1]);
                            double DWDY = (Program.W[i][j + 1][k] - Program.W[i][j - 1][k]) / (Program.ZAY[j] + Program.ZAY[j - 1]);

                            double DUDZ = 0;
                            double DVDZ = 0;
                            double DWDZ = 0;
                            double DTDZ = 0;

                            if (k > 1)
                            {
                                DUDZ = (Program.U[i][j][k + 1] - Program.U[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]);
                                DVDZ = (Program.V[i][j][k + 1] - Program.V[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]);
                                DWDZ = (Program.W[i][j][k + 1] - Program.W[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]);
                                DTDZ = (Program.T[i][j][k + 1] - Program.T[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]) - 0.00065;
                            }
                            else
                            {
                                DUDZ = (Program.U[i][j][k + 1] - Program.U[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]);
                                DVDZ = (Program.V[i][j][k + 1] - Program.V[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]);
                                DWDZ = (Program.W[i][j][k + 1] - Program.W[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]);
                                DTDZ = (Program.T[i][j][k + 1] - Program.T[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]) - 0.00065;
                            }

                            //production-terms for the dissipation rate
                            double PSTRESS = Program.VISV[i][j][k] *
                                (0.5F * (DUDX * DUDX + DVDY * DVDY + DWDZ * DWDZ) +
                                 (DUDY + DVDX) * (DUDY + DVDX) +
                                 (DUDZ + DWDX) * (DUDZ + DWDX) +
                                 (DVDZ + DWDY) * (DVDZ + DWDY)) * Program.RHO[i][j][k];
                            double PBUOY = -Program.VISV[i][j][k] * DTDZ * Program.GERD / (Program.T[i][j][k] + Program.TBZ1) / Program.PRTE * Program.RHO[i][j][k];

                            double DIM = Program.AWEST_PS[i][j][k] * Program.DISSN[i - 1][j][k] + Program.ASOUTH_PS[i][j][k] * Program.DISSN[i][j - 1][k] +
                                Program.AEAST_PS[i][j][k] * Program.DISSN[i + 1][j][k] + Program.ANORTH_PS[i][j][k] * Program.DISSN[i][j + 1][k] +
                                Program.AP0_PS[i][j][k] * Program.DISS[i][j][k] +
                                (Program.DISS[i][j][k] / Math.Max(Math.Abs(Program.TE[i][j][k]), 0.01) * (1.44 * PSTRESS + 0.4 * PBUOY - 1.92 * Program.DISS[i][j][k] * Program.RHO[i][j][k])) * Program.VOLImm[i][j][k];

                            //Recurrence formula
                            if (k > 1)
                            {
                                help = 1 / (Program.A_PS[i][j][k] - Program.C_PS[i][j][k] * PIM[k - 1]);
                                PIM[k] = Program.B_PS[i][j][k] * help;
                                QIM[k] = (DIM + Program.C_PS[i][j][k] * QIM[k - 1]) * help;
                            }
                            else
                            {
                                help = 1 / Program.A_PS[i][j][k];
                                PIM[k] = Program.B_PS[i][j][k] * help;
                                QIM[k] = DIM * help;
                            }
                        }

                        //Obtain new TKE-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            help = Program.DISSN[i][j][k];
                            help += (Program.RELAXT * (PIM[k] * Program.DISSN[i][j][k + 1] + QIM[k] - help));
                            Program.DISSN[i][j][k] = help;

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
                Span<double> PIM = stackalloc double[NK + 1];
                Span<double> QIM = stackalloc double[NK + 1];
                double help;

                for (int i1 = range.Item1; i1 < range.Item2; i1++)
                {
                    int i = NI - i1 + 1;
                    for (int j = 2; j <= NJ - 1; j++)
                    {
                        for (int k = 1; k <= NK - 1; k++)
                        {
                            //Production-terms of turbulent kinetic energy
                            double DUDX = (Program.U[i + 1][j][k] - Program.U[i - 1][j][k]) / (Program.ZAX[i] + Program.ZAX[i - 1]);
                            double DVDX = (Program.V[i + 1][j][k] - Program.V[i - 1][j][k]) / (Program.ZAX[i] + Program.ZAX[i - 1]);
                            double DWDX = (Program.W[i + 1][j][k] - Program.W[i - 1][j][k]) / (Program.ZAX[i] + Program.ZAX[i - 1]);
                            double DUDY = (Program.U[i][j + 1][k] - Program.U[i][j - 1][k]) / (Program.ZAY[j] + Program.ZAY[j - 1]);
                            double DVDY = (Program.V[i][j + 1][k] - Program.V[i][j - 1][k]) / (Program.ZAY[j] + Program.ZAY[j - 1]);
                            double DWDY = (Program.W[i][j + 1][k] - Program.W[i][j - 1][k]) / (Program.ZAY[j] + Program.ZAY[j - 1]);

                            double DUDZ = 0;
                            double DVDZ = 0;
                            double DWDZ = 0;
                            double DTDZ = 0;

                            if (k > 1)
                            {
                                DUDZ = (Program.U[i][j][k + 1] - Program.U[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]);
                                DVDZ = (Program.V[i][j][k + 1] - Program.V[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]);
                                DWDZ = (Program.W[i][j][k + 1] - Program.W[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]);
                                DTDZ = (Program.T[i][j][k + 1] - Program.T[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]) - 0.00065;
                            }
                            else
                            {
                                DUDZ = (Program.U[i][j][k + 1] - Program.U[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]);
                                DVDZ = (Program.V[i][j][k + 1] - Program.V[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]);
                                DWDZ = (Program.W[i][j][k + 1] - Program.W[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]);
                                DTDZ = (Program.T[i][j][k + 1] - Program.T[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]) - 0.00065;
                            }

                            //production-terms for the dissipation rate
                            double PSTRESS = Program.VISV[i][j][k] *
                                (0.5F * (DUDX * DUDX + DVDY * DVDY + DWDZ * DWDZ) +
                                 (DUDY + DVDX) * (DUDY + DVDX) +
                                 (DUDZ + DWDX) * (DUDZ + DWDX) +
                                 (DVDZ + DWDY) * (DVDZ + DWDY)) * Program.RHO[i][j][k];
                            double PBUOY = -Program.VISV[i][j][k] * DTDZ * Program.GERD / (Program.T[i][j][k] + Program.TBZ1) / Program.PRTE * Program.RHO[i][j][k];

                            double DIM = Program.AWEST_PS[i][j][k] * Program.DISSN[i - 1][j][k] + Program.ASOUTH_PS[i][j][k] * Program.DISSN[i][j - 1][k] +
                                Program.AEAST_PS[i][j][k] * Program.DISSN[i + 1][j][k] + Program.ANORTH_PS[i][j][k] * Program.DISSN[i][j + 1][k] +
                                Program.AP0_PS[i][j][k] * Program.DISS[i][j][k] +
                                (Program.DISS[i][j][k] / Math.Max(Math.Abs(Program.TE[i][j][k]), 0.01) * (1.44 * PSTRESS + 0.4 * PBUOY - 1.92 * Program.DISS[i][j][k] * Program.RHO[i][j][k])) * Program.VOLImm[i][j][k];

                            //Recurrence formula
                            if (k > 1)
                            {
                                help = 1 / (Program.A_PS[i][j][k] - Program.C_PS[i][j][k] * PIM[k - 1]);
                                PIM[k] = Program.B_PS[i][j][k] * help;
                                QIM[k] = (DIM + Program.C_PS[i][j][k] * QIM[k - 1]) * help;
                            }
                            else
                            {
                                help = 1 / Program.A_PS[i][j][k];
                                PIM[k] = Program.B_PS[i][j][k] * help;
                                QIM[k] = DIM * help;
                            }
                        }

                        //Obtain new TKE-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            help = Program.DISSN[i][j][k];
                            help += (Program.RELAXT * (PIM[k] * Program.DISSN[i][j][k + 1] + QIM[k] - help));
                            Program.DISSN[i][j][k] = help;

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
                Span<double> PIM = stackalloc double[NK + 1];
                Span<double> QIM = stackalloc double[NK + 1];
                double help;

                for (int i = range.Item1; i < range.Item2; i++)
                {
                    for (int j = NJ - 1; j >= 2; j--)
                    {
                        for (int k = 1; k <= NK - 1; k++)
                        {
                            //Production-terms of turbulent kinetic energy
                            double DUDX = (Program.U[i + 1][j][k] - Program.U[i - 1][j][k]) / (Program.ZAX[i] + Program.ZAX[i - 1]);
                            double DVDX = (Program.V[i + 1][j][k] - Program.V[i - 1][j][k]) / (Program.ZAX[i] + Program.ZAX[i - 1]);
                            double DWDX = (Program.W[i + 1][j][k] - Program.W[i - 1][j][k]) / (Program.ZAX[i] + Program.ZAX[i - 1]);
                            double DUDY = (Program.U[i][j + 1][k] - Program.U[i][j - 1][k]) / (Program.ZAY[j] + Program.ZAY[j - 1]);
                            double DVDY = (Program.V[i][j + 1][k] - Program.V[i][j - 1][k]) / (Program.ZAY[j] + Program.ZAY[j - 1]);
                            double DWDY = (Program.W[i][j + 1][k] - Program.W[i][j - 1][k]) / (Program.ZAY[j] + Program.ZAY[j - 1]);

                            double DUDZ = 0;
                            double DVDZ = 0;
                            double DWDZ = 0;
                            double DTDZ = 0;

                            if (k > 1)
                            {
                                DUDZ = (Program.U[i][j][k + 1] - Program.U[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]);
                                DVDZ = (Program.V[i][j][k + 1] - Program.V[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]);
                                DWDZ = (Program.W[i][j][k + 1] - Program.W[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]);
                                DTDZ = (Program.T[i][j][k + 1] - Program.T[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]) - 0.00065;
                            }
                            else
                            {
                                DUDZ = (Program.U[i][j][k + 1] - Program.U[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]);
                                DVDZ = (Program.V[i][j][k + 1] - Program.V[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]);
                                DWDZ = (Program.W[i][j][k + 1] - Program.W[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]);
                                DTDZ = (Program.T[i][j][k + 1] - Program.T[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]) - 0.00065;
                            }

                            //production-terms for the dissipation rate
                            double PSTRESS = Program.VISV[i][j][k] *
                                (0.5F * (DUDX * DUDX + DVDY * DVDY + DWDZ * DWDZ) +
                                 (DUDY + DVDX) * (DUDY + DVDX) +
                                 (DUDZ + DWDX) * (DUDZ + DWDX) +
                                 (DVDZ + DWDY) * (DVDZ + DWDY)) * Program.RHO[i][j][k];
                            double PBUOY = -Program.VISV[i][j][k] * DTDZ * Program.GERD / (Program.T[i][j][k] + Program.TBZ1) / Program.PRTE * Program.RHO[i][j][k];

                            double DIM = Program.AWEST_PS[i][j][k] * Program.DISSN[i - 1][j][k] + Program.ASOUTH_PS[i][j][k] * Program.DISSN[i][j - 1][k] +
                                Program.AEAST_PS[i][j][k] * Program.DISSN[i + 1][j][k] + Program.ANORTH_PS[i][j][k] * Program.DISSN[i][j + 1][k] +
                                Program.AP0_PS[i][j][k] * Program.DISS[i][j][k] +
                                (Program.DISS[i][j][k] / Math.Max(Math.Abs(Program.TE[i][j][k]), 0.01) * (1.44 * PSTRESS + 0.4 * PBUOY - 1.92 * Program.DISS[i][j][k] * Program.RHO[i][j][k])) * Program.VOLImm[i][j][k];

                            //Recurrence formula
                            if (k > 1)
                            {
                                help = 1 / (Program.A_PS[i][j][k] - Program.C_PS[i][j][k] * PIM[k - 1]);
                                PIM[k] = Program.B_PS[i][j][k] * help;
                                QIM[k] = (DIM + Program.C_PS[i][j][k] * QIM[k - 1]) * help;
                            }
                            else
                            {
                                help = 1 / Program.A_PS[i][j][k];
                                PIM[k] = Program.B_PS[i][j][k] * help;
                                QIM[k] = DIM * help;
                            }
                        }

                        //Obtain new TKE-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            help = Program.DISSN[i][j][k];
                            help += (Program.RELAXT * (PIM[k] * Program.DISSN[i][j][k + 1] + QIM[k] - help));
                            Program.DISSN[i][j][k] = help;

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
                Span<double> PIM = stackalloc double[NK + 1];
                Span<double> QIM = stackalloc double[NK + 1];
                double help;

                for (int j = range.Item1; j < range.Item2; j++)
                {
                    for (int i = 2; i <= NI - 1; i++)
                    {
                        for (int k = 1; k <= NK - 1; k++)
                        {
                            //Production-terms of turbulent kinetic energy
                            double DUDX = (Program.U[i + 1][j][k] - Program.U[i - 1][j][k]) / (Program.ZAX[i] + Program.ZAX[i - 1]);
                            double DVDX = (Program.V[i + 1][j][k] - Program.V[i - 1][j][k]) / (Program.ZAX[i] + Program.ZAX[i - 1]);
                            double DWDX = (Program.W[i + 1][j][k] - Program.W[i - 1][j][k]) / (Program.ZAX[i] + Program.ZAX[i - 1]);
                            double DUDY = (Program.U[i][j + 1][k] - Program.U[i][j - 1][k]) / (Program.ZAY[j] + Program.ZAY[j - 1]);
                            double DVDY = (Program.V[i][j + 1][k] - Program.V[i][j - 1][k]) / (Program.ZAY[j] + Program.ZAY[j - 1]);
                            double DWDY = (Program.W[i][j + 1][k] - Program.W[i][j - 1][k]) / (Program.ZAY[j] + Program.ZAY[j - 1]);

                            double DUDZ = 0;
                            double DVDZ = 0;
                            double DWDZ = 0;
                            double DTDZ = 0;

                            if (k > 1)
                            {
                                DUDZ = (Program.U[i][j][k + 1] - Program.U[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]);
                                DVDZ = (Program.V[i][j][k + 1] - Program.V[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]);
                                DWDZ = (Program.W[i][j][k + 1] - Program.W[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]);
                                DTDZ = (Program.T[i][j][k + 1] - Program.T[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]) - 0.00065;
                            }
                            else
                            {
                                DUDZ = (Program.U[i][j][k + 1] - Program.U[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]);
                                DVDZ = (Program.V[i][j][k + 1] - Program.V[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]);
                                DWDZ = (Program.W[i][j][k + 1] - Program.W[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]);
                                DTDZ = (Program.T[i][j][k + 1] - Program.T[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]) - 0.00065;
                            }

                            //production-terms for the dissipation rate
                            double PSTRESS = Program.VISV[i][j][k] *
                                (0.5F * (DUDX * DUDX + DVDY * DVDY + DWDZ * DWDZ) +
                                 (DUDY + DVDX) * (DUDY + DVDX) +
                                 (DUDZ + DWDX) * (DUDZ + DWDX) +
                                 (DVDZ + DWDY) * (DVDZ + DWDY)) * Program.RHO[i][j][k];
                            double PBUOY = -Program.VISV[i][j][k] * DTDZ * Program.GERD / (Program.T[i][j][k] + Program.TBZ1) / Program.PRTE * Program.RHO[i][j][k];

                            double DIM = Program.AWEST_PS[i][j][k] * Program.DISSN[i - 1][j][k] + Program.ASOUTH_PS[i][j][k] * Program.DISSN[i][j - 1][k] +
                                Program.AEAST_PS[i][j][k] * Program.DISSN[i + 1][j][k] + Program.ANORTH_PS[i][j][k] * Program.DISSN[i][j + 1][k] +
                                Program.AP0_PS[i][j][k] * Program.DISS[i][j][k] +
                                (Program.DISS[i][j][k] / Math.Max(Math.Abs(Program.TE[i][j][k]), 0.01) * (1.44 * PSTRESS + 0.4 * PBUOY - 1.92 * Program.DISS[i][j][k] * Program.RHO[i][j][k])) * Program.VOLImm[i][j][k];

                            //Recurrence formula
                            if (k > 1)
                            {
                                help = 1 / (Program.A_PS[i][j][k] - Program.C_PS[i][j][k] * PIM[k - 1]);
                                PIM[k] = Program.B_PS[i][j][k] * help;
                                QIM[k] = (DIM + Program.C_PS[i][j][k] * QIM[k - 1]) * help;
                            }
                            else
                            {
                                help = 1 / Program.A_PS[i][j][k];
                                PIM[k] = Program.B_PS[i][j][k] * help;
                                QIM[k] = DIM * help;
                            }
                        }

                        //Obtain new TKE-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            help = Program.DISSN[i][j][k];
                            help += (Program.RELAXT * (PIM[k] * Program.DISSN[i][j][k + 1] + QIM[k] - help));
                            Program.DISSN[i][j][k] = help;

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
                Span<double> PIM = stackalloc double[NK + 1];
                Span<double> QIM = stackalloc double[NK + 1];
                double help;

                for (int j1 = range.Item1; j1 < range.Item2; j1++)
                {
                    int j = NJ - j1 + 1;
                    for (int i = NI - 1; i >= 2; i--)
                    {
                        for (int k = 1; k <= NK - 1; k++)
                        {
                            //Production-terms of turbulent kinetic energy
                            double DUDX = (Program.U[i + 1][j][k] - Program.U[i - 1][j][k]) / (Program.ZAX[i] + Program.ZAX[i - 1]);
                            double DVDX = (Program.V[i + 1][j][k] - Program.V[i - 1][j][k]) / (Program.ZAX[i] + Program.ZAX[i - 1]);
                            double DWDX = (Program.W[i + 1][j][k] - Program.W[i - 1][j][k]) / (Program.ZAX[i] + Program.ZAX[i - 1]);
                            double DUDY = (Program.U[i][j + 1][k] - Program.U[i][j - 1][k]) / (Program.ZAY[j] + Program.ZAY[j - 1]);
                            double DVDY = (Program.V[i][j + 1][k] - Program.V[i][j - 1][k]) / (Program.ZAY[j] + Program.ZAY[j - 1]);
                            double DWDY = (Program.W[i][j + 1][k] - Program.W[i][j - 1][k]) / (Program.ZAY[j] + Program.ZAY[j - 1]);

                            double DUDZ = 0;
                            double DVDZ = 0;
                            double DWDZ = 0;
                            double DTDZ = 0;

                            if (k > 1)
                            {
                                DUDZ = (Program.U[i][j][k + 1] - Program.U[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]);
                                DVDZ = (Program.V[i][j][k + 1] - Program.V[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]);
                                DWDZ = (Program.W[i][j][k + 1] - Program.W[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]);
                                DTDZ = (Program.T[i][j][k + 1] - Program.T[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]) - 0.00065;
                            }
                            else
                            {
                                DUDZ = (Program.U[i][j][k + 1] - Program.U[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]);
                                DVDZ = (Program.V[i][j][k + 1] - Program.V[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]);
                                DWDZ = (Program.W[i][j][k + 1] - Program.W[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]);
                                DTDZ = (Program.T[i][j][k + 1] - Program.T[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]) - 0.00065;
                            }

                            //production-terms for the dissipation rate
                            double PSTRESS = Program.VISV[i][j][k] *
                                (0.5F * (DUDX * DUDX + DVDY * DVDY + DWDZ * DWDZ) +
                                 (DUDY + DVDX) * (DUDY + DVDX) +
                                 (DUDZ + DWDX) * (DUDZ + DWDX) +
                                 (DVDZ + DWDY) * (DVDZ + DWDY)) * Program.RHO[i][j][k];
                            double PBUOY = -Program.VISV[i][j][k] * DTDZ * Program.GERD / (Program.T[i][j][k] + Program.TBZ1) / Program.PRTE * Program.RHO[i][j][k];

                            double DIM = Program.AWEST_PS[i][j][k] * Program.DISSN[i - 1][j][k] + Program.ASOUTH_PS[i][j][k] * Program.DISSN[i][j - 1][k] +
                                Program.AEAST_PS[i][j][k] * Program.DISSN[i + 1][j][k] + Program.ANORTH_PS[i][j][k] * Program.DISSN[i][j + 1][k] +
                                Program.AP0_PS[i][j][k] * Program.DISS[i][j][k] +
                                (Program.DISS[i][j][k] / Math.Max(Math.Abs(Program.TE[i][j][k]), 0.01) * (1.44 * PSTRESS + 0.4 * PBUOY - 1.92 * Program.DISS[i][j][k] * Program.RHO[i][j][k])) * Program.VOLImm[i][j][k];

                            //Recurrence formula
                            if (k > 1)
                            {
                                help = 1 / (Program.A_PS[i][j][k] - Program.C_PS[i][j][k] * PIM[k - 1]);
                                PIM[k] = Program.B_PS[i][j][k] * help;
                                QIM[k] = (DIM + Program.C_PS[i][j][k] * QIM[k - 1]) * help;
                            }
                            else
                            {
                                help = 1 / Program.A_PS[i][j][k];
                                PIM[k] = Program.B_PS[i][j][k] * help;
                                QIM[k] = DIM * help;
                            }
                        }

                        //Obtain new TKE-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            help = Program.DISSN[i][j][k];
                            help += (Program.RELAXT * (PIM[k] * Program.DISSN[i][j][k + 1] + QIM[k] - help));
                            Program.DISSN[i][j][k] = help;

                        }
                    }
                }
            });

            range_parallel = (int)(NJ / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2 + 6);
            range_parallel = Math.Max(36 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
            range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
                                                           // Parallel.For(2, NJ, Program.pOptions, j =>
            Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
            {
                Span<double> PIM = stackalloc double[NK + 1];
                Span<double> QIM = stackalloc double[NK + 1];
                double help;

                for (int j = range.Item1; j < range.Item2; j++)
                {
                    for (int i = NI - 1; i >= 2; i--)
                    {
                        for (int k = 1; k <= NK - 1; k++)
                        {
                            //Production-terms of turbulent kinetic energy
                            double DUDX = (Program.U[i + 1][j][k] - Program.U[i - 1][j][k]) / (Program.ZAX[i] + Program.ZAX[i - 1]);
                            double DVDX = (Program.V[i + 1][j][k] - Program.V[i - 1][j][k]) / (Program.ZAX[i] + Program.ZAX[i - 1]);
                            double DWDX = (Program.W[i + 1][j][k] - Program.W[i - 1][j][k]) / (Program.ZAX[i] + Program.ZAX[i - 1]);
                            double DUDY = (Program.U[i][j + 1][k] - Program.U[i][j - 1][k]) / (Program.ZAY[j] + Program.ZAY[j - 1]);
                            double DVDY = (Program.V[i][j + 1][k] - Program.V[i][j - 1][k]) / (Program.ZAY[j] + Program.ZAY[j - 1]);
                            double DWDY = (Program.W[i][j + 1][k] - Program.W[i][j - 1][k]) / (Program.ZAY[j] + Program.ZAY[j - 1]);

                            double DUDZ = 0;
                            double DVDZ = 0;
                            double DWDZ = 0;
                            double DTDZ = 0;

                            if (k > 1)
                            {
                                DUDZ = (Program.U[i][j][k + 1] - Program.U[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]);
                                DVDZ = (Program.V[i][j][k + 1] - Program.V[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]);
                                DWDZ = (Program.W[i][j][k + 1] - Program.W[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]);
                                DTDZ = (Program.T[i][j][k + 1] - Program.T[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]) - 0.00065;
                            }
                            else
                            {
                                DUDZ = (Program.U[i][j][k + 1] - Program.U[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]);
                                DVDZ = (Program.V[i][j][k + 1] - Program.V[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]);
                                DWDZ = (Program.W[i][j][k + 1] - Program.W[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]);
                                DTDZ = (Program.T[i][j][k + 1] - Program.T[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]) - 0.00065;
                            }

                            //production-terms for the dissipation rate
                            double PSTRESS = Program.VISV[i][j][k] *
                                (0.5F * (DUDX * DUDX + DVDY * DVDY + DWDZ * DWDZ) +
                                 (DUDY + DVDX) * (DUDY + DVDX) +
                                 (DUDZ + DWDX) * (DUDZ + DWDX) +
                                 (DVDZ + DWDY) * (DVDZ + DWDY)) * Program.RHO[i][j][k];
                            double PBUOY = -Program.VISV[i][j][k] * DTDZ * Program.GERD / (Program.T[i][j][k] + Program.TBZ1) / Program.PRTE * Program.RHO[i][j][k];

                            double DIM = Program.AWEST_PS[i][j][k] * Program.DISSN[i - 1][j][k] + Program.ASOUTH_PS[i][j][k] * Program.DISSN[i][j - 1][k] +
                                Program.AEAST_PS[i][j][k] * Program.DISSN[i + 1][j][k] + Program.ANORTH_PS[i][j][k] * Program.DISSN[i][j + 1][k] +
                                Program.AP0_PS[i][j][k] * Program.DISS[i][j][k] +
                                (Program.DISS[i][j][k] / Math.Max(Math.Abs(Program.TE[i][j][k]), 0.01) * (1.44 * PSTRESS + 0.4 * PBUOY - 1.92 * Program.DISS[i][j][k] * Program.RHO[i][j][k])) * Program.VOLImm[i][j][k];

                            //Recurrence formula
                            if (k > 1)
                            {
                                help = 1 / (Program.A_PS[i][j][k] - Program.C_PS[i][j][k] * PIM[k - 1]);
                                PIM[k] = Program.B_PS[i][j][k] * help;
                                QIM[k] = (DIM + Program.C_PS[i][j][k] * QIM[k - 1]) * help;
                            }
                            else
                            {
                                help = 1 / Program.A_PS[i][j][k];
                                PIM[k] = Program.B_PS[i][j][k] * help;
                                QIM[k] = DIM * help;
                            }
                        }

                        //Obtain new TKE-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            help = Program.DISSN[i][j][k];
                            help += (Program.RELAXT * (PIM[k] * Program.DISSN[i][j][k + 1] + QIM[k] - help));
                            Program.DISSN[i][j][k] = help;

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
                Span<double> PIM = stackalloc double[NK + 1];
                Span<double> QIM = stackalloc double[NK + 1];
                double help;

                for (int j1 = range.Item1; j1 < range.Item2; j1++)
                {
                    int j = NJ - j1 + 1;

                    for (int i = 2; i <= NI - 1; i++)
                    {
                        for (int k = 1; k <= NK - 1; k++)
                        {
                            //Production-terms of turbulent kinetic energy
                            double DUDX = (Program.U[i + 1][j][k] - Program.U[i - 1][j][k]) / (Program.ZAX[i] + Program.ZAX[i - 1]);
                            double DVDX = (Program.V[i + 1][j][k] - Program.V[i - 1][j][k]) / (Program.ZAX[i] + Program.ZAX[i - 1]);
                            double DWDX = (Program.W[i + 1][j][k] - Program.W[i - 1][j][k]) / (Program.ZAX[i] + Program.ZAX[i - 1]);
                            double DUDY = (Program.U[i][j + 1][k] - Program.U[i][j - 1][k]) / (Program.ZAY[j] + Program.ZAY[j - 1]);
                            double DVDY = (Program.V[i][j + 1][k] - Program.V[i][j - 1][k]) / (Program.ZAY[j] + Program.ZAY[j - 1]);
                            double DWDY = (Program.W[i][j + 1][k] - Program.W[i][j - 1][k]) / (Program.ZAY[j] + Program.ZAY[j - 1]);

                            double DUDZ = 0;
                            double DVDZ = 0;
                            double DWDZ = 0;
                            double DTDZ = 0;

                            if (k > 1)
                            {
                                DUDZ = (Program.U[i][j][k + 1] - Program.U[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]);
                                DVDZ = (Program.V[i][j][k + 1] - Program.V[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]);
                                DWDZ = (Program.W[i][j][k + 1] - Program.W[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]);
                                DTDZ = (Program.T[i][j][k + 1] - Program.T[i][j][k - 1]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k - 1]) - 0.00065;
                            }
                            else
                            {
                                DUDZ = (Program.U[i][j][k + 1] - Program.U[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]);
                                DVDZ = (Program.V[i][j][k + 1] - Program.V[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]);
                                DWDZ = (Program.W[i][j][k + 1] - Program.W[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]);
                                DTDZ = (Program.T[i][j][k + 1] - Program.T[i][j][k]) / (Program.ZSPImm[i][j][k + 1] + Program.ZSPImm[i][j][k]) - 0.00065;
                            }

                            //production-terms for the dissipation rate
                            double PSTRESS = Program.VISV[i][j][k] *
                                (0.5F * (DUDX * DUDX + DVDY * DVDY + DWDZ * DWDZ) +
                                 (DUDY + DVDX) * (DUDY + DVDX) +
                                 (DUDZ + DWDX) * (DUDZ + DWDX) +
                                 (DVDZ + DWDY) * (DVDZ + DWDY)) * Program.RHO[i][j][k];
                            double PBUOY = -Program.VISV[i][j][k] * DTDZ * Program.GERD / (Program.T[i][j][k] + Program.TBZ1) / Program.PRTE * Program.RHO[i][j][k];

                            double DIM = Program.AWEST_PS[i][j][k] * Program.DISSN[i - 1][j][k] + Program.ASOUTH_PS[i][j][k] * Program.DISSN[i][j - 1][k] +
                                Program.AEAST_PS[i][j][k] * Program.DISSN[i + 1][j][k] + Program.ANORTH_PS[i][j][k] * Program.DISSN[i][j + 1][k] +
                                Program.AP0_PS[i][j][k] * Program.DISS[i][j][k] +
                                (Program.DISS[i][j][k] / Math.Max(Math.Abs(Program.TE[i][j][k]), 0.01) * (1.44 * PSTRESS + 0.4 * PBUOY - 1.92 * Program.DISS[i][j][k] * Program.RHO[i][j][k])) * Program.VOLImm[i][j][k];

                            //Recurrence formula
                            if (k > 1)
                            {
                                help = 1 / (Program.A_PS[i][j][k] - Program.C_PS[i][j][k] * PIM[k - 1]);
                                PIM[k] = Program.B_PS[i][j][k] * help;
                                QIM[k] = (DIM + Program.C_PS[i][j][k] * QIM[k - 1]) * help;
                            }
                            else
                            {
                                help = 1 / Program.A_PS[i][j][k];
                                PIM[k] = Program.B_PS[i][j][k] * help;
                                QIM[k] = DIM * help;
                            }
                        }

                        //Obtain new TKE-components
                        for (int k = NK - 1; k >= 1; k--)
                        {
                            help = Program.DISSN[i][j][k];
                            help += (Program.RELAXT * (PIM[k] * Program.DISSN[i][j][k + 1] + QIM[k] - help));
                            Program.DISSN[i][j][k] = help;

                        }
                    }
                }
            });
        }
    }
}
