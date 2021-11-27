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

        /// <summary>
        /// Solve the non-hydrostaic pressure
        /// </summary>
        /// <param name="NI"></param>
        /// <param name="NJ"></param>
        /// <param name="NK"></param>
        public static void Primp_calculate(int NI, int NJ, int NK)
        {
            if (Program.INUMS == 1)
            {
                int range_parallel = NI / Program.pOptions.MaxDegreeOfParallelism - (StripeCounter % 6);
                range_parallel = Math.Max(Program.StripeWidth - (StripeCounter % 6), range_parallel); // min. Program.StripeWidth cells per processor
                StripeCounter++;
                range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel

                //Parallel.For(2, NI, Program.pOptions, i =>
                Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
                {
                    int NK_P = NK; int NJ_P = NJ;
                    Span<double> PIM = stackalloc double[2 * NK + 2];
                    Span<double> QIM = stackalloc double[2 * NK + 2];
                    double[] DP_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double[] DPZ_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double[] DPXi_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double[] DPYj_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double TERMP = 0;
                    double DIM, RHOAIM_LL, RHOAJM_LL, RHOiAIM_LL;

                    for (int i = range.Item1; i < range.Item2; ++i)
                    {
                        int border = Math.Min(i - range.Item1, range.Item2 - 1 - i);
                        for (int j = 2; j <= NJ_P - 1; ++j)
                        {
                            int iP1 = i + 1;
                            int jP1 = j + 1;
                            ReadOnlySpan<double> AB_L = Program.AB[i][j];
                            ReadOnlySpan<double> AE_L = Program.AE[i][j];
                            ReadOnlySpan<double> AN_L = Program.AN[i][j];
                            ReadOnlySpan<double> AP_L = Program.AP[i][j];
                            ReadOnlySpan<double> AS_L = Program.AS[i][j];
                            ReadOnlySpan<double> AT_L = Program.AT[i][j];
                            ReadOnlySpan<double> AW_L = Program.AW[i][j];
                            /*
                        double[] AIM_L    = Program.AIM[i][j];
                        double[] AIMi_L   = Program.AIM[iP1][j];
                        double[] AIMj_L   = Program.AIM[i][jP1];
							 */
                            ReadOnlySpan<float> AIM_L = Program.AP0[i][j];
                            ReadOnlySpan<float> AIMi_L = Program.AP0[iP1][j];
                            ReadOnlySpan<float> AIMj_L = Program.AP0[i][jP1];
                            double[] DP_L;
                            double[] DPZ_L;
                            double[] DPXi_L;
                            double[] DPYj_L;
                            //Avoid race conditions at the border cells of the sequential calculated stripes
                            if (border < 2)
                            {
                                DP_L = DP_LR;
                                DPZ_L = DPZ_LR;
                                DPXi_L = DPXi_LR;
                                DPYj_L = DPYj_LR;
                                Program.CopyArrayLockSource(Program.DP[i][j], DP_L);
                                Program.CopyArrayLockSource(Program.DPZ[i][j], DPZ_L);
                                Program.CopyArrayLockSource(Program.DPX[iP1][j], DPXi_L);
                                Program.CopyArrayLockSource(Program.DPY[i][jP1], DPYj_L);
                            }
                            else
                            {
                                DP_L = Program.DP[i][j];
                                DPZ_L = Program.DPZ[i][j];
                                DPXi_L = Program.DPX[iP1][j];
                                DPYj_L = Program.DPY[i][jP1];
                            }
                            ReadOnlySpan<double> DPi_L = Program.DP[iP1][j];
                            ReadOnlySpan<double> DPj_L = Program.DP[i][jP1];
                            ReadOnlySpan<double> DPX_L = Program.DPX[i][j];
                            //ReadOnlySpan<double> DPXi_L = Program.DPX[iP1][j];
                            ReadOnlySpan<double> DPY_L = Program.DPY[i][j];
                            //ReadOnlySpan<double> DPYj_L = Program.DPY[i][jP1];
                            //ReadOnlySpan<double> DPZ_L = Program.DPZ[i][j];
                            ReadOnlySpan<double> DPZi_L = Program.DPZ[iP1][j];
                            ReadOnlySpan<double> DPZj_L = Program.DPZ[i][jP1];
                            ReadOnlySpan<float> AREAX_L = Program.AREAXImm[i][j];
                            ReadOnlySpan<float> AREAXi_L = Program.AREAXImm[iP1][j];
                            ReadOnlySpan<float> AREAY_L = Program.AREAYImm[i][j];
                            ReadOnlySpan<float> AREAYj_L = Program.AREAYImm[i][jP1];
                            ReadOnlySpan<float> AREAZ_L = Program.AREAZImm[i][j];
                            ReadOnlySpan<float> AREAXYZ_L = Program.AREAXYZImm[i][j];
                            ReadOnlySpan<float> AREAZX_L = Program.AREAZXImm[i][j];
                            ReadOnlySpan<float> AREAZXi_L = Program.AREAZXImm[iP1][j];
                            ReadOnlySpan<float> AREAZY_L = Program.AREAZYImm[i][j];
                            ReadOnlySpan<float> AREAZYj_L = Program.AREAZYImm[i][jP1];
                            ReadOnlySpan<float> RHO_L = Program.RHOImm[i][j];
                            ReadOnlySpan<float> RHOi_L = Program.RHOImm[iP1][j];
                            ReadOnlySpan<float> RHOj_L = Program.RHOImm[i][jP1];
                            ReadOnlySpan<float> SUXi_L = Program.SUX[iP1][j];
                            ReadOnlySpan<float> SUXYZ_L = Program.SUXYZ[i][j];
                            ReadOnlySpan<float> SUZ_L = Program.SUZ[i][j];
                            ReadOnlySpan<float> SUY_L = Program.SUY[i][jP1];

                            int m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; --kn)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)

                                //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = SUZ_L[k] * AREAZ_L[k] - AE_L[kn] * DPXi_L[k - 1] -
                                        AN_L[kn] * DPYj_L[k - 1] - AW_L[kn] * DPX_L[k] -
                                        AS_L[kn] * DPY_L[k];

                                    //Recurrence formula
                                    TERMP = 1 / (AP_L[kn] - AT_L[kn] * PIM[kn + 1]);
                                    PIM[kn] = AB_L[kn] * TERMP;
                                    QIM[kn] = (DIM + AT_L[kn] * QIM[kn + 1]) * TERMP;
                                    m--;
                                }
                                //Coefficients for the upper half-cell
                                else
                                {
                                    if (kn < 2 * NK_P)
                                    {
                                        DIM = SUXYZ_L[k] * AREAXYZ_L[k] + AW_L[kn] * DPX_L[k] +
                                            AS_L[kn] * DPY_L[k] + AE_L[kn] * DPXi_L[k] +
                                            AN_L[kn] * DPYj_L[k];

                                        //Recurrence formula
                                        TERMP = 1 / (AP_L[kn] - AT_L[kn] * PIM[kn + 1]);
                                        PIM[kn] = AB_L[kn] * TERMP;
                                        QIM[kn] = (DIM + AT_L[kn] * QIM[kn + 1]) * TERMP;
                                    }

                                    m = 2;
                                }
                            }

                            //Obtain new P-components
                            m = 2;
                            DP_L[0] = 0;
                            DP_L[NK_P] = 0;
                            DPZ_L[NK_P] = 0;
                            DPXi_L[NK_P] = 0;
                            DPYj_L[NK_P] = 0;

                            for (int kn1 = 1; kn1 <= 2 * (NK_P - 1); ++kn1)
                            {
                                int k1 = 1 + kn1 >> 1;  // (int) (kn1 * 0.5F + 0.5F)

                                if (m == 2)
                                {
                                    DPZ_L[k1] = PIM[kn1] * DP_L[k1 - 1] + QIM[kn1];
                                    m--;
                                }
                                else
                                {
                                    DP_L[k1] = PIM[kn1] * DPZ_L[k1] + QIM[kn1];
                                    m = 2;
                                }
                            }

                            //compute neighbours
                            for (int k2 = 1; k2 <= NK_P - 1; ++k2)
                            {
                                int kn2 = 2 * k2 - 1;
                                int k2P = k2 + 1;
                                RHOAIM_LL = RHO_L[k2] / AIM_L[k2];
                                RHOAJM_LL = RHOj_L[k2] / AIMj_L[k2];
                                RHOiAIM_LL = RHOi_L[k2] / AIMi_L[k2];

                                DPXi_L[k2] = 1 / AREAXi_L[k2] / (RHOAIM_LL + RHOiAIM_LL) *
                                    (SUXi_L[k2] + RHOAIM_LL * (DP_L[k2] * (AREAX_L[k2] +
                                                                            AREAZX_L[k2]) - DPZ_L[k2P] * AREAZX_L[k2P]) - RHOiAIM_LL
                                     * (DPZi_L[k2] * AREAZXi_L[k2] - DPi_L[k2] * (AREAXi_L[k2] + AREAZXi_L[k2])));


                                DPYj_L[k2] = 1 / AREAYj_L[k2] / (RHOAIM_LL + RHOAJM_LL) *
                                    (SUY_L[k2] + RHOAIM_LL * (DP_L[k2] * (AREAY_L[k2] +
                                                                          AREAZY_L[k2]) - DPZ_L[k2P] * AREAZY_L[k2P]) - RHOAJM_LL
                                     * (DPZj_L[k2] * AREAZYj_L[k2] - DPj_L[k2] * (AREAYj_L[k2] + AREAZYj_L[k2])));
                            }

                            if (border < 2)
                            {
                                Program.CopyArrayLockDest(DP_L, Program.DP[i][j]);
                                Program.CopyArrayLockDest(DPZ_L, Program.DPZ[i][j]);
                                Program.CopyArrayLockDest(DPXi_L, Program.DPX[iP1][j]);
                                Program.CopyArrayLockDest(DPYj_L, Program.DPY[i][jP1]);
                            }
                        }
                    }
                    Program.GrammArrayPool.Return(DP_LR);
                    Program.GrammArrayPool.Return(DPZ_LR);
                    Program.GrammArrayPool.Return(DPXi_LR);
                    Program.GrammArrayPool.Return(DPYj_LR);
                });
            }

            if (Program.INUMS == 2)
            {
                int range_parallel = NI / Program.pOptions.MaxDegreeOfParallelism - (StripeCounter % 6);
                range_parallel = Math.Max(Program.StripeWidth - (StripeCounter % 6), range_parallel); // min. Program.StripeWidth cells per processor
                StripeCounter++;
                range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel

                //Parallel.For(2, NI, Program.pOptions, i1 =>
                Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
                {
                    int NK_P = NK; int NJ_P = NJ;
                    Span<double> PIM = stackalloc double[2 * NK + 2];
                    Span<double> QIM = stackalloc double[2 * NK + 2];
                    double[] DP_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double[] DPZ_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double[] DPX_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double[] DPY_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double TERMP = 0;
                    double DIM, RHOAIM_LL, RHOAJM_LL, RHOiAIM_LL;

                    for (int i = range.Item2 - 1; i >= range.Item1; --i)
                    {
                        int border = Math.Min(i - range.Item1, range.Item2 - 1 - i);
                        for (int j = NJ_P - 1; j >= 2; --j)
                        {
                            int iP1 = i + 1;
                            int jP1 = j + 1;
                            ReadOnlySpan<double> AB_L = Program.AB[i][j];
                            ReadOnlySpan<double> AE_L = Program.AE[i][j];
                            ReadOnlySpan<double> AN_L = Program.AN[i][j];
                            ReadOnlySpan<double> AP_L = Program.AP[i][j];
                            ReadOnlySpan<double> AS_L = Program.AS[i][j];
                            ReadOnlySpan<double> AT_L = Program.AT[i][j];
                            ReadOnlySpan<double> AW_L = Program.AW[i][j];
                            /*
                        double[] AIM_L    = Program.AIM[i][j];
                        double[] AIMmi_L   = Program.AIM[i - 1][j];
                        double[] AIMj_L   = Program.AIM[i][jP1];
                        double[] AIMmj_L  = Program.AIM[i][j - 1];
            				 */
                            ReadOnlySpan<float> AIM_L = Program.AP0[i][j];
                            ReadOnlySpan<float> AIMmi_L = Program.AP0[i - 1][j];
                            ReadOnlySpan<float> AIMj_L = Program.AP0[i][jP1];
                            ReadOnlySpan<float> AIMmj_L = Program.AP0[i][j - 1];

                            double[] DP_L;
                            double[] DPZ_L;
                            double[] DPX_L;
                            double[] DPY_L;

                            if (border < 2)
                            {
                                DP_L = DP_LR;
                                DPZ_L = DPZ_LR;
                                DPX_L = DPX_LR;
                                DPY_L = DPY_LR;
                                Program.CopyArrayLockSource(Program.DP[i][j], DP_L);
                                Program.CopyArrayLockSource(Program.DPZ[i][j], DPZ_L);
                                Program.CopyArrayLockSource(Program.DPX[i][j], DPX_L);
                                Program.CopyArrayLockSource(Program.DPY[i][j], DPY_L);
                            }
                            else
                            {
                                DP_L = Program.DP[i][j];
                                DPZ_L = Program.DPZ[i][j];
                                DPX_L = Program.DPX[i][j];
                                DPY_L = Program.DPY[i][j];
                            }
                            //double[] DP_L = Program.DP[i][j];
                            ReadOnlySpan<double> DPmi_L = Program.DP[i - 1][j];
                            ReadOnlySpan<double> DPmj_L = Program.DP[i][j - 1];
                            //ReadOnlySpan<double> DPX_L = Program.DPX[i][j];
                            ReadOnlySpan<double> DPXi_L = Program.DPX[iP1][j];
                            //ReadOnlySpan<double> DPY_L = Program.DPY[i][j];
                            ReadOnlySpan<double> DPYj_L = Program.DPY[i][jP1];
                            //ReadOnlySpan<double> DPZ_L = Program.DPZ[i][j];
                            ReadOnlySpan<double> DPZmj_L = Program.DPZ[i][j - 1];
                            ReadOnlySpan<double> DPZmi_L = Program.DPZ[i - 1][j];
                            ReadOnlySpan<float> AREAX_L = Program.AREAXImm[i][j];
                            ReadOnlySpan<float> AREAXmi_L = Program.AREAXImm[i - 1][j];
                            ReadOnlySpan<float> AREAY_L = Program.AREAYImm[i][j];
                            ReadOnlySpan<float> AREAYmj_L = Program.AREAYImm[i][j - 1];
                            ReadOnlySpan<float> AREAZ_L = Program.AREAZImm[i][j];
                            ReadOnlySpan<float> AREAXYZ_L = Program.AREAXYZImm[i][j];
                            ReadOnlySpan<float> AREAZX_L = Program.AREAZXImm[i][j];
                            ReadOnlySpan<float> AREAZXi_L = Program.AREAZXImm[i - 1][j];
                            ReadOnlySpan<float> AREAZY_L = Program.AREAZYImm[i][j];
                            ReadOnlySpan<float> AREAZYmj_L = Program.AREAZYImm[i][j - 1];
                            ReadOnlySpan<float> RHO_L = Program.RHOImm[i][j];
                            ReadOnlySpan<float> RHOmi_L = Program.RHOImm[i - 1][j];
                            ReadOnlySpan<float> RHOmj_L = Program.RHOImm[i][j - 1];
                            ReadOnlySpan<float> SUXYZ_L = Program.SUXYZ[i][j];
                            ReadOnlySpan<float> SUX_L = Program.SUX[i][j];
                            ReadOnlySpan<float> SUY_L = Program.SUY[i][j];
                            ReadOnlySpan<float> SUZ_L = Program.SUZ[i][j];

                            int m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; --kn)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)

                                //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = SUZ_L[k] * AREAZ_L[k] - AE_L[kn] * DPXi_L[k - 1] -
                                        AN_L[kn] * DPYj_L[k - 1] - AW_L[kn] * DPX_L[k] -
                                        AS_L[kn] * DPY_L[k];
                                    //Recurrence formula
                                    TERMP = 1 / (AP_L[kn] - AT_L[kn] * PIM[kn + 1]);
                                    PIM[kn] = AB_L[kn] * TERMP;
                                    QIM[kn] = (DIM + AT_L[kn] * QIM[kn + 1]) * TERMP;
                                    m--;
                                }
                                //Coefficients for the upper half-cell
                                else
                                {
                                    if (kn < 2 * NK_P)
                                    {
                                        DIM = SUXYZ_L[k] * AREAXYZ_L[k] + AW_L[kn] * DPX_L[k] +
                                            AS_L[kn] * DPY_L[k] + AE_L[kn] * DPXi_L[k] +
                                            AN_L[kn] * DPYj_L[k];

                                        //Recurrence formula
                                        TERMP = 1 / (AP_L[kn] - AT_L[kn] * PIM[kn + 1]);
                                        PIM[kn] = AB_L[kn] * TERMP;
                                        QIM[kn] = (DIM + AT_L[kn] * QIM[kn + 1]) * TERMP;
                                    }
                                    m = 2;
                                }
                            }

                            //Obtain new P-components
                            m = 2;
                            DP_L[0] = 0;
                            DP_L[NK_P] = 0;
                            DPZ_L[NK_P] = 0;
                            DPX_L[NK_P] = 0;
                            DPY_L[NK_P] = 0;

                            for (int kn1 = 1; kn1 <= 2 * (NK_P - 1); ++kn1)
                            {
                                int k1 = 1 + kn1 >> 1;  // (int) (kn1 * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    DPZ_L[k1] = PIM[kn1] * DP_L[k1 - 1] + QIM[kn1];
                                    m--;
                                }
                                else
                                {
                                    DP_L[k1] = PIM[kn1] * DPZ_L[k1] + QIM[kn1];
                                    m = 2;
                                }
                            }

                            //compute neighbours
                            for (int k2 = 1; k2 <= NK_P - 1; ++k2)
                            {
                                int kn2 = 2 * k2 - 1;
                                int k2P = k2 + 1;
                                RHOAIM_LL = RHO_L[k2] / AIM_L[k2];
                                RHOAJM_LL = RHOmj_L[k2] / AIMmj_L[k2];
                                RHOiAIM_LL = RHOmi_L[k2] / AIMmi_L[k2];

                                DPX_L[k2] = 1 / AREAX_L[k2] / (RHOiAIM_LL + RHOAIM_LL) *
                                    (SUX_L[k2] + RHOiAIM_LL * (DPmi_L[k2] * (AREAXmi_L[k2] +
                                                                             AREAZXi_L[k2]) - DPZmi_L[k2P] * AREAZXi_L[k2P]) - RHOAIM_LL
                                     * (DPZ_L[k2] * AREAZX_L[k2] - DP_L[k2] * (AREAX_L[k2] + AREAZX_L[k2])));


                                DPY_L[k2] = 1 / AREAY_L[k2] / (RHOAJM_LL + RHOAIM_LL) *
                                    (SUY_L[k2] + RHOAJM_LL * (DPmj_L[k2] * (AREAYmj_L[k2] +
                                                                            AREAZYmj_L[k2]) - DPZmj_L[k2P] * AREAZYmj_L[k2P]) - RHOAIM_LL
                                     * (DPZ_L[k2] * AREAZY_L[k2] - DP_L[k2] * (AREAY_L[k2] + AREAZY_L[k2])));
                            }

                            if (border < 2)
                            {
                                Program.CopyArrayLockDest(DP_L, Program.DP[i][j]);
                                Program.CopyArrayLockDest(DPZ_L, Program.DPZ[i][j]);
                                Program.CopyArrayLockDest(DPX_L, Program.DPX[i][j]);
                                Program.CopyArrayLockDest(DPY_L, Program.DPY[i][j]);
                            }
                        }
                    }
                    Program.GrammArrayPool.Return(DP_LR);
                    Program.GrammArrayPool.Return(DPZ_LR);
                    Program.GrammArrayPool.Return(DPX_LR);
                    Program.GrammArrayPool.Return(DPY_LR);
                });
            }

            if (Program.INUMS == 3)
            {

                int range_parallel = NJ / Program.pOptions.MaxDegreeOfParallelism - (StripeCounter % 6);
                range_parallel = Math.Max(Program.StripeWidth - (StripeCounter % 6), range_parallel); // min. Program.StripeWidth cells per processor
                StripeCounter++;
                range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
                                                               //Parallel.For(2, NJ, Program.pOptions, j1 =>
                Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
                {
                    int NK_P = NK; int NI_P = NI;
                    Span<double> PIM = stackalloc double[2 * NK_P + 2];
                    Span<double> QIM = stackalloc double[2 * NK_P + 2];
                    double[] DP_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double[] DPZ_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double[] DPX_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double[] DPY_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double TERMP = 0;
                    double DIM, RHOAIM_LL, RHOAJM_LL, RHOiAIM_LL;

                    for (int j = range.Item2 - 1; j >= range.Item1; --j)
                    {
                        int border = Math.Min(j - range.Item1, range.Item2 - 1 - j);
                        for (int i = NI_P - 1; i >= 2; --i)
                        {
                            int iP1 = i + 1;
                            int jP1 = j + 1;
                            ReadOnlySpan<double> AB_L = Program.AB[i][j];
                            ReadOnlySpan<double> AE_L = Program.AE[i][j];
                            ReadOnlySpan<double> AN_L = Program.AN[i][j];
                            ReadOnlySpan<double> AP_L = Program.AP[i][j];
                            ReadOnlySpan<double> AS_L = Program.AS[i][j];
                            ReadOnlySpan<double> AT_L = Program.AT[i][j];
                            ReadOnlySpan<double> AW_L = Program.AW[i][j];
                            /*
                        double[] AIM_L    = Program.AIM[i][j];
                        double[] AIMmi_L   = Program.AIM[i - 1][j];
                        double[] AIMj_L   = Program.AIM[i][jP1];
                        double[] AIMmj_L  = Program.AIM[i][j - 1];
            				 */
                            ReadOnlySpan<float> AIM_L = Program.AP0[i][j];
                            ReadOnlySpan<float> AIMmi_L = Program.AP0[i - 1][j];
                            ReadOnlySpan<float> AIMj_L = Program.AP0[i][jP1];
                            ReadOnlySpan<float> AIMmj_L = Program.AP0[i][j - 1];

                            double[] DP_L;
                            double[] DPZ_L;
                            double[] DPX_L;
                            double[] DPY_L;

                            if (border < 2)
                            {
                                DP_L = DP_LR;
                                DPZ_L = DPZ_LR;
                                DPX_L = DPX_LR;
                                DPY_L = DPY_LR;
                                Program.CopyArrayLockSource(Program.DP[i][j], DP_L);
                                Program.CopyArrayLockSource(Program.DPZ[i][j], DPZ_L);
                                Program.CopyArrayLockSource(Program.DPX[i][j], DPX_L);
                                Program.CopyArrayLockSource(Program.DPY[i][j], DPY_L);
                            }
                            else
                            {
                                DP_L = Program.DP[i][j];
                                DPZ_L = Program.DPZ[i][j];
                                DPX_L = Program.DPX[i][j];
                                DPY_L = Program.DPY[i][j];
                            }

                            //ReadOnlySpan<double> DP_L = Program.DP[i][j];
                            ReadOnlySpan<double> DPmi_L = Program.DP[i - 1][j];
                            ReadOnlySpan<double> DPmj_L = Program.DP[i][j - 1];
                            //ReadOnlySpan<double> DPX_L = Program.DPX[i][j];
                            ReadOnlySpan<double> DPXi_L = Program.DPX[iP1][j];
                            //ReadOnlySpan<double> DPY_L = Program.DPY[i][j];
                            ReadOnlySpan<double> DPYj_L = Program.DPY[i][jP1];
                            //ReadOnlySpan<double> DPZ_L = Program.DPZ[i][j];
                            ReadOnlySpan<double> DPZmi_L = Program.DPZ[i - 1][j];
                            ReadOnlySpan<double> DPZmj_L = Program.DPZ[i][j - 1];
                            ReadOnlySpan<float> AREAX_L = Program.AREAXImm[i][j];
                            ReadOnlySpan<float> AREAXmi_L = Program.AREAXImm[i - 1][j];
                            ReadOnlySpan<float> AREAY_L = Program.AREAYImm[i][j];
                            ReadOnlySpan<float> AREAYmj_L = Program.AREAYImm[i][j - 1];
                            ReadOnlySpan<float> AREAZ_L = Program.AREAZImm[i][j];
                            ReadOnlySpan<float> AREAXYZ_L = Program.AREAXYZImm[i][j];
                            ReadOnlySpan<float> AREAZX_L = Program.AREAZXImm[i][j];
                            ReadOnlySpan<float> AREAZXi_L = Program.AREAZXImm[i - 1][j];
                            ReadOnlySpan<float> AREAZY_L = Program.AREAZYImm[i][j];
                            ReadOnlySpan<float> AREAZYmj_L = Program.AREAZYImm[i][j - 1];
                            ReadOnlySpan<float> RHO_L = Program.RHOImm[i][j];
                            ReadOnlySpan<float> RHOmi_L = Program.RHOImm[i - 1][j];
                            ReadOnlySpan<float> RHOmj_L = Program.RHOImm[i][j - 1];
                            ReadOnlySpan<float> SUXYZ_L = Program.SUXYZ[i][j];
                            ReadOnlySpan<float> SUX_L = Program.SUX[i][j];
                            ReadOnlySpan<float> SUY_L = Program.SUY[i][j];
                            ReadOnlySpan<float> SUZ_L = Program.SUZ[i][j];

                            int m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; --kn)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)

                                //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = SUZ_L[k] * AREAZ_L[k] - AE_L[kn] * DPXi_L[k - 1] -
                                        AN_L[kn] * DPYj_L[k - 1] - AW_L[kn] * DPX_L[k] -
                                        AS_L[kn] * DPY_L[k];

                                    //Recurrence formula
                                    TERMP = 1 / (AP_L[kn] - AT_L[kn] * PIM[kn + 1]);
                                    PIM[kn] = AB_L[kn] * TERMP;
                                    QIM[kn] = (DIM + AT_L[kn] * QIM[kn + 1]) * TERMP;
                                    m--;
                                }
                                //Coefficients for the upper half-cell
                                else
                                {
                                    if (kn < 2 * NK_P)
                                    {
                                        DIM = SUXYZ_L[k] * AREAXYZ_L[k] + AW_L[kn] * DPX_L[k] +
                                            AS_L[kn] * DPY_L[k] + AE_L[kn] * DPXi_L[k] +
                                            AN_L[kn] * DPYj_L[k];
                                        //Recurrence formula
                                        TERMP = 1 / (AP_L[kn] - AT_L[kn] * PIM[kn + 1]);
                                        PIM[kn] = AB_L[kn] * TERMP;
                                        QIM[kn] = (DIM + AT_L[kn] * QIM[kn + 1]) * TERMP;
                                    }
                                    m = 2;
                                }
                            }

                            //Obtain new P-components
                            m = 2;
                            DP_L[0] = 0;
                            DP_L[NK_P] = 0;
                            DPZ_L[NK_P] = 0;
                            DPX_L[NK_P] = 0;
                            DPY_L[NK_P] = 0;

                            for (int kn1 = 1; kn1 <= 2 * (NK_P - 1); ++kn1)
                            {
                                int k1 = 1 + kn1 >> 1;  // (int) (kn1 * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    DPZ_L[k1] = PIM[kn1] * DP_L[k1 - 1] + QIM[kn1];
                                    m--;
                                }
                                else
                                {
                                    DP_L[k1] = PIM[kn1] * DPZ_L[k1] + QIM[kn1];
                                    m = 2;
                                }
                            }

                            //compute neighbours
                            for (int k2 = 1; k2 <= NK_P - 1; ++k2)
                            {
                                int kn2 = 2 * k2 - 1;
                                int k2P = k2 + 1;
                                RHOAIM_LL = RHO_L[k2] / AIM_L[k2];
                                RHOAJM_LL = RHOmj_L[k2] / AIMmj_L[k2];
                                RHOiAIM_LL = RHOmi_L[k2] / AIMmi_L[k2];

                                DPX_L[k2] = 1 / AREAX_L[k2] / (RHOiAIM_LL + RHOAIM_LL) *
                                    (SUX_L[k2] + RHOiAIM_LL * (DPmi_L[k2] * (AREAXmi_L[k2] +
                                                                             AREAZXi_L[k2]) - DPZmi_L[k2P] * AREAZXi_L[k2P]) - RHOAIM_LL
                                     * (DPZ_L[k2] * AREAZX_L[k2] - DP_L[k2] * (AREAX_L[k2] + AREAZX_L[k2])));

                                DPY_L[k2] = 1 / AREAY_L[k2] / (RHOAJM_LL + RHOAIM_LL) *
                                    (SUY_L[k2] + RHOAJM_LL * (DPmj_L[k2] * (AREAYmj_L[k2] +
                                                                            AREAZYmj_L[k2]) - DPZmj_L[k2P] * AREAZYmj_L[k2P]) - RHOAIM_LL
                                     * (DPZ_L[k2] * AREAZY_L[k2] - DP_L[k2] * (AREAY_L[k2] + AREAZY_L[k2])));
                            }

                            if (border < 2)
                            {
                                Program.CopyArrayLockDest(DP_L, Program.DP[i][j]);
                                Program.CopyArrayLockDest(DPZ_L, Program.DPZ[i][j]);
                                Program.CopyArrayLockDest(DPX_L, Program.DPX[i][j]);
                                Program.CopyArrayLockDest(DPY_L, Program.DPY[i][j]);
                            }
                        }
                    }
                    Program.GrammArrayPool.Return(DP_LR);
                    Program.GrammArrayPool.Return(DPZ_LR);
                    Program.GrammArrayPool.Return(DPX_LR);
                    Program.GrammArrayPool.Return(DPY_LR);
                });
            }

            if (Program.INUMS == 4)
            {
                int range_parallel = NJ / Program.pOptions.MaxDegreeOfParallelism - (StripeCounter % 6);
                range_parallel = Math.Max(Program.StripeWidth - (StripeCounter % 6), range_parallel); // min. Program.StripeWidth cells per processor
                StripeCounter++;
                range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel

                //Parallel.For(2, NJ, Program.pOptions, j =>
                Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
                {
                    int NK_P = NK; int NI_P = NI;
                    Span<double> PIM = stackalloc double[2 * NK_P + 2];
                    Span<double> QIM = stackalloc double[2 * NK_P + 2];
                    double[] DP_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double[] DPZ_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double[] DPXi_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double[] DPYj_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double TERMP = 0;
                    double DIM, RHOAIM_LL, RHOAJM_LL, RHOiAIM_LL;

                    for (int j = range.Item1; j < range.Item2; ++j)
                    {
                        int border = Math.Min(j - range.Item1, range.Item2 - 1 - j);
                        for (int i = 2; i <= NI_P - 1; ++i)
                        {
                            int iP1 = i + 1;
                            int jP1 = j + 1;
                            ReadOnlySpan<double> AB_L = Program.AB[i][j];
                            ReadOnlySpan<double> AE_L = Program.AE[i][j];
                            ReadOnlySpan<double> AN_L = Program.AN[i][j];
                            ReadOnlySpan<double> AP_L = Program.AP[i][j];
                            ReadOnlySpan<double> AS_L = Program.AS[i][j];
                            ReadOnlySpan<double> AT_L = Program.AT[i][j];
                            ReadOnlySpan<double> AW_L = Program.AW[i][j];
                            /*
                        double[] AIM_L    = Program.AIM[i][j];
                        double[] AIMi_L   = Program.AIM[iP1][j];
                        double[] AIMj_L   = Program.AIM[i][jP1];
            				 */
                            ReadOnlySpan<float> AIM_L = Program.AP0[i][j];
                            ReadOnlySpan<float> AIMi_L = Program.AP0[iP1][j];
                            ReadOnlySpan<float> AIMj_L = Program.AP0[i][jP1];

                            double[] DP_L;
                            double[] DPZ_L;
                            double[] DPXi_L;
                            double[] DPYj_L;

                            if (border < 2)
                            {
                                DP_L = DP_LR;
                                DPZ_L = DPZ_LR;
                                DPXi_L = DPXi_LR;
                                DPYj_L = DPYj_LR;
                                Program.CopyArrayLockSource(Program.DP[i][j], DP_L);
                                Program.CopyArrayLockSource(Program.DPZ[i][j], DPZ_L);
                                Program.CopyArrayLockSource(Program.DPX[iP1][j], DPXi_L);
                                Program.CopyArrayLockSource(Program.DPY[i][jP1], DPYj_L);
                            }
                            else
                            {
                                DP_L = Program.DP[i][j];
                                DPZ_L = Program.DPZ[i][j];
                                DPXi_L = Program.DPX[iP1][j];
                                DPYj_L = Program.DPY[i][jP1];
                            }
                            //ReadOnlySpan<double> DP_L = Program.DP[i][j];
                            ReadOnlySpan<double> DPi_L = Program.DP[iP1][j];
                            ReadOnlySpan<double> DPj_L = Program.DP[i][jP1];
                            ReadOnlySpan<double> DPX_L = Program.DPX[i][j];
                            //ReadOnlySpan<double> DPXi_L = Program.DPX[iP1][j];
                            ReadOnlySpan<double> DPY_L = Program.DPY[i][j];
                            //ReadOnlySpan<double> DPYj_L = Program.DPY[i][jP1];
                            //ReadOnlySpan<double> DPZ_L = Program.DPZ[i][j];
                            ReadOnlySpan<double> DPZi_L = Program.DPZ[iP1][j];
                            ReadOnlySpan<double> DPZj_L = Program.DPZ[i][jP1];
                            ReadOnlySpan<float> AREAX_L = Program.AREAXImm[i][j];
                            ReadOnlySpan<float> AREAXi_L = Program.AREAXImm[iP1][j];
                            ReadOnlySpan<float> AREAY_L = Program.AREAYImm[i][j];
                            ReadOnlySpan<float> AREAYj_L = Program.AREAYImm[i][jP1];
                            ReadOnlySpan<float> AREAZ_L = Program.AREAZImm[i][j];
                            ReadOnlySpan<float> AREAXYZ_L = Program.AREAXYZImm[i][j];
                            ReadOnlySpan<float> AREAZX_L = Program.AREAZXImm[i][j];
                            ReadOnlySpan<float> AREAZXi_L = Program.AREAZXImm[iP1][j];
                            ReadOnlySpan<float> AREAZY_L = Program.AREAZYImm[i][j];
                            ReadOnlySpan<float> AREAZYj_L = Program.AREAZYImm[i][jP1];
                            ReadOnlySpan<float> RHO_L = Program.RHOImm[i][j];
                            ReadOnlySpan<float> RHOi_L = Program.RHOImm[iP1][j];
                            ReadOnlySpan<float> RHOj_L = Program.RHOImm[i][jP1];
                            ReadOnlySpan<float> SUXi_L = Program.SUX[iP1][j];
                            ReadOnlySpan<float> SUXYZ_L = Program.SUXYZ[i][j];
                            ReadOnlySpan<float> SUZ_L = Program.SUZ[i][j];
                            ReadOnlySpan<float> SUY_L = Program.SUY[i][jP1];

                            int m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; --kn)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)

                                //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = SUZ_L[k] * AREAZ_L[k] - AE_L[kn] * DPXi_L[k - 1] -
                                        AN_L[kn] * DPYj_L[k - 1] - AW_L[kn] * DPX_L[k] -
                                        AS_L[kn] * DPY_L[k];

                                    //Recurrence formula
                                    TERMP = 1 / (AP_L[kn] - AT_L[kn] * PIM[kn + 1]);
                                    PIM[kn] = AB_L[kn] * TERMP;
                                    QIM[kn] = (DIM + AT_L[kn] * QIM[kn + 1]) * TERMP;
                                    m--;
                                }
                                //Coefficients for the upper half-cell
                                else
                                {
                                    if (kn < 2 * NK_P)
                                    {
                                        DIM = SUXYZ_L[k] * AREAXYZ_L[k] + AW_L[kn] * DPX_L[k] +
                                            AS_L[kn] * DPY_L[k] + AE_L[kn] * DPXi_L[k] +
                                            AN_L[kn] * DPYj_L[k];

                                        //Recurrence formula
                                        TERMP = 1 / (AP_L[kn] - AT_L[kn] * PIM[kn + 1]);
                                        PIM[kn] = AB_L[kn] * TERMP;
                                        QIM[kn] = (DIM + AT_L[kn] * QIM[kn + 1]) * TERMP;
                                    }
                                    m = 2;
                                }
                            }
                            //Obtain new P-components
                            m = 2;
                            DP_L[0] = 0;
                            DP_L[NK_P] = 0;
                            DPZ_L[NK_P] = 0;
                            DPXi_L[NK_P] = 0;
                            DPYj_L[NK_P] = 0;

                            for (int kn1 = 1; kn1 <= 2 * (NK_P - 1); ++kn1)
                            {
                                int k1 = 1 + kn1 >> 1;  // (int) (kn1 * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    DPZ_L[k1] = PIM[kn1] * DP_L[k1 - 1] + QIM[kn1];
                                    m--;
                                }
                                else
                                {
                                    DP_L[k1] = PIM[kn1] * DPZ_L[k1] + QIM[kn1];
                                    m = 2;
                                }
                            }

                            //compute neighbours
                            for (int k2 = 1; k2 <= NK_P - 1; ++k2)
                            {
                                int kn2 = 2 * k2 - 1;
                                int k2P = k2 + 1;
                                RHOAIM_LL = RHO_L[k2] / AIM_L[k2];
                                RHOAJM_LL = RHOj_L[k2] / AIMj_L[k2];
                                RHOiAIM_LL = RHOi_L[k2] / AIMi_L[k2];

                                DPXi_L[k2] = 1 / AREAXi_L[k2] / (RHOAIM_LL + RHOiAIM_LL) *
                                    (SUXi_L[k2] + RHOAIM_LL * (DP_L[k2] * (AREAX_L[k2] +
                                                                           AREAZX_L[k2]) - DPZ_L[k2P] * AREAZX_L[k2P]) - RHOiAIM_LL
                                     * (DPZi_L[k2] * AREAZXi_L[k2] - DPi_L[k2] * (AREAXi_L[k2] + AREAZXi_L[k2])));

                                DPYj_L[k2] = 1 / AREAYj_L[k2] / (RHOAIM_LL + RHOAJM_LL) *
                                    (SUY_L[k2] + RHOAIM_LL * (DP_L[k2] * (AREAY_L[k2] +
                                                                          AREAZY_L[k2]) - DPZ_L[k2P] * AREAZY_L[k2P]) - RHOAJM_LL
                                     * (DPZj_L[k2] * AREAZYj_L[k2] - DPj_L[k2] * (AREAYj_L[k2] + AREAZYj_L[k2])));
                            }
                            if (border < 2)
                            {
                                Program.CopyArrayLockDest(DP_L, Program.DP[i][j]);
                                Program.CopyArrayLockDest(DPZ_L, Program.DPZ[i][j]);
                                Program.CopyArrayLockDest(DPXi_L, Program.DPX[iP1][j]);
                                Program.CopyArrayLockDest(DPYj_L, Program.DPY[i][jP1]);
                            }
                        }
                    }
                    Program.GrammArrayPool.Return(DP_LR);
                    Program.GrammArrayPool.Return(DPZ_LR);
                    Program.GrammArrayPool.Return(DPXi_LR);
                    Program.GrammArrayPool.Return(DPYj_LR);
                });
            }

            if (Program.INUMS == 5)
            {
                int range_parallel = NI / Program.pOptions.MaxDegreeOfParallelism - (StripeCounter % 6);
                range_parallel = Math.Max(Program.StripeWidth - (StripeCounter % 6), range_parallel); // min. Program.StripeWidth cells per processor
                StripeCounter++;
                range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
                                                               //  Parallel.For(2, NI, Program.pOptions, i1 =>
                Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
                {
                    int NK_P = NK; int NJ_P = NJ;
                    Span<double> PIM = stackalloc double[2 * NK_P + 2];
                    Span<double> QIM = stackalloc double[2 * NK_P + 2];
                    double[] DP_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double[] DPZ_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double[] DPX_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double[] DPYj_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double TERMP = 0;
                    double DIM, RHOAIM_LL, RHOAJM_LL, RHOiAIM_LL;

                    for (int i = range.Item2 - 1; i >= range.Item1; --i)
                    {
                        int border = Math.Min(i - range.Item1, range.Item2 - 1 - i);
                        for (int j = 2; j <= NJ_P - 1; ++j)
                        {
                            int iP1 = i + 1;
                            int jP1 = j + 1;
                            ReadOnlySpan<double> AB_L = Program.AB[i][j];
                            ReadOnlySpan<double> AE_L = Program.AE[i][j];
                            ReadOnlySpan<double> AN_L = Program.AN[i][j];
                            ReadOnlySpan<double> AP_L = Program.AP[i][j];
                            ReadOnlySpan<double> AS_L = Program.AS[i][j];
                            ReadOnlySpan<double> AT_L = Program.AT[i][j];
                            ReadOnlySpan<double> AW_L = Program.AW[i][j];
                            /*
                        double[] AIM_L    = Program.AIM[i][j];
                        double[] AIMmi_L   = Program.AIM[i - 1][j];
                        double[] AIMj_L   = Program.AIM[i][jP1];
                             */
                            ReadOnlySpan<float> AIM_L = Program.AP0[i][j];
                            ReadOnlySpan<float> AIMmi_L = Program.AP0[i - 1][j];
                            ReadOnlySpan<float> AIMj_L = Program.AP0[i][jP1];

                            double[] DP_L;
                            double[] DPZ_L;
                            double[] DPX_L;
                            double[] DPYj_L;

                            if (border < 2)
                            {
                                DP_L = DP_LR;
                                DPZ_L = DPZ_LR;
                                DPX_L = DPX_LR;
                                DPYj_L = DPYj_LR;
                                Program.CopyArrayLockSource(Program.DP[i][j], DP_L);
                                Program.CopyArrayLockSource(Program.DPZ[i][j], DPZ_L);
                                Program.CopyArrayLockSource(Program.DPX[i][j], DPX_L);
                                Program.CopyArrayLockSource(Program.DPY[i][jP1], DPYj_L);
                            }
                            else
                            {
                                DP_L = Program.DP[i][j];
                                DPZ_L = Program.DPZ[i][j];
                                DPX_L = Program.DPX[i][j];
                                DPYj_L = Program.DPY[i][jP1];
                            }

                            //ReadOnlySpan<double> DP_L = Program.DP[i][j];
                            ReadOnlySpan<double> DPj_L = Program.DP[i][jP1];
                            ReadOnlySpan<double> DPmi_L = Program.DP[i - 1][j];
                            //ReadOnlySpan<double> DPX_L = Program.DPX[i][j];
                            ReadOnlySpan<double> DPXi_L = Program.DPX[iP1][j];
                            ReadOnlySpan<double> DPY_L = Program.DPY[i][j];
                            // ReadOnlySpan<double> DPYj_L = Program.DPY[i][jP1];
                            //ReadOnlySpan<double> DPZ_L = Program.DPZ[i][j];
                            ReadOnlySpan<double> DPZmi_L = Program.DPZ[i - 1][j];
                            ReadOnlySpan<double> DPZj_L = Program.DPZ[i][jP1];
                            ReadOnlySpan<float> AREAX_L = Program.AREAXImm[i][j];
                            ReadOnlySpan<float> AREAXmi_L = Program.AREAXImm[i - 1][j];
                            ReadOnlySpan<float> AREAY_L = Program.AREAYImm[i][j];
                            ReadOnlySpan<float> AREAYj_L = Program.AREAYImm[i][jP1];
                            ReadOnlySpan<float> AREAZ_L = Program.AREAZImm[i][j];
                            ReadOnlySpan<float> AREAXYZ_L = Program.AREAXYZImm[i][j];
                            ReadOnlySpan<float> AREAZX_L = Program.AREAZXImm[i][j];
                            ReadOnlySpan<float> AREAZXi_L = Program.AREAZXImm[i - 1][j];
                            ReadOnlySpan<float> AREAZY_L = Program.AREAZYImm[i][j];
                            ReadOnlySpan<float> AREAZYj_L = Program.AREAZYImm[i][jP1];
                            ReadOnlySpan<float> RHO_L = Program.RHOImm[i][j];
                            ReadOnlySpan<float> RHOmi_L = Program.RHOImm[i - 1][j];
                            ReadOnlySpan<float> RHOj_L = Program.RHOImm[i][jP1];
                            ReadOnlySpan<float> SUXYZ_L = Program.SUXYZ[i][j];
                            ReadOnlySpan<float> SUX_L = Program.SUX[i][j];
                            ReadOnlySpan<float> SUY_L = Program.SUY[i][jP1];
                            ReadOnlySpan<float> SUZ_L = Program.SUZ[i][j];

                            int m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; --kn)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)

                                //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = SUZ_L[k] * AREAZ_L[k] - AE_L[kn] * DPXi_L[k - 1] -
                                        AN_L[kn] * DPYj_L[k - 1] - AW_L[kn] * DPX_L[k] -
                                        AS_L[kn] * DPY_L[k];

                                    //Recurrence formula
                                    TERMP = 1 / (AP_L[kn] - AT_L[kn] * PIM[kn + 1]);
                                    PIM[kn] = AB_L[kn] * TERMP;
                                    QIM[kn] = (DIM + AT_L[kn] * QIM[kn + 1]) * TERMP;
                                    m--;
                                }
                                //Coefficients for the upper half-cell
                                else
                                {
                                    if (kn < 2 * NK_P)
                                    {
                                        DIM = SUXYZ_L[k] * AREAXYZ_L[k] + AW_L[kn] * DPX_L[k] +
                                            AS_L[kn] * DPY_L[k] + AE_L[kn] * DPXi_L[k] +
                                            AN_L[kn] * DPYj_L[k];

                                        //Recurrence formula
                                        TERMP = 1 / (AP_L[kn] - AT_L[kn] * PIM[kn + 1]);
                                        PIM[kn] = AB_L[kn] * TERMP;
                                        QIM[kn] = (DIM + AT_L[kn] * QIM[kn + 1]) * TERMP;
                                    }
                                    m = 2;
                                }
                            }
                            //Obtain new P-components
                            m = 2;
                            DP_L[0] = 0;
                            DP_L[NK_P] = 0;
                            DPZ_L[NK_P] = 0;
                            DPX_L[NK_P] = 0;
                            DPYj_L[NK_P] = 0;

                            for (int kn1 = 1; kn1 <= 2 * (NK_P - 1); ++kn1)
                            {
                                int k1 = 1 + kn1 >> 1;  // (int) (kn1 * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    DPZ_L[k1] = PIM[kn1] * DP_L[k1 - 1] + QIM[kn1];
                                    m--;
                                }
                                else
                                {
                                    DP_L[k1] = PIM[kn1] * DPZ_L[k1] + QIM[kn1];
                                    m = 2;
                                }
                            }

                            //compute neighbours
                            for (int k2 = 1; k2 <= NK_P - 1; ++k2)
                            {
                                int kn2 = 2 * k2 - 1;
                                int k2P = k2 + 1;
                                RHOAIM_LL = RHO_L[k2] / AIM_L[k2];
                                RHOAJM_LL = RHOj_L[k2] / AIMj_L[k2];
                                RHOiAIM_LL = RHOmi_L[k2] / AIMmi_L[k2];

                                DPX_L[k2] = 1 / AREAX_L[k2] / (RHOiAIM_LL + RHOAIM_LL) *
                                    (SUX_L[k2] + RHOiAIM_LL * (DPmi_L[k2] * (AREAXmi_L[k2] +
                                                                             AREAZXi_L[k2]) - DPZmi_L[k2P] * AREAZXi_L[k2P]) - RHOAIM_LL
                                     * (DPZ_L[k2] * AREAZX_L[k2] - DP_L[k2] * (AREAX_L[k2] + AREAZX_L[k2])));

                                DPYj_L[k2] = 1 / AREAYj_L[k2] / (RHOAIM_LL + RHOAJM_LL) *
                                    (SUY_L[k2] + RHOAIM_LL * (DP_L[k2] * (AREAY_L[k2] +
                                                                          AREAZY_L[k2]) - DPZ_L[k2P] * AREAZY_L[k2P]) - RHOAJM_LL
                                     * (DPZj_L[k2] * AREAZYj_L[k2] - DPj_L[k2] * (AREAYj_L[k2] + AREAZYj_L[k2])));
                            }
                            if (border < 2)
                            {
                                Program.CopyArrayLockDest(DP_L, Program.DP[i][j]);
                                Program.CopyArrayLockDest(DPZ_L, Program.DPZ[i][j]);
                                Program.CopyArrayLockDest(DPX_L, Program.DPX[i][j]);
                                Program.CopyArrayLockDest(DPYj_L, Program.DPY[i][jP1]);
                            }
                        }
                    }
                    Program.GrammArrayPool.Return(DP_LR);
                    Program.GrammArrayPool.Return(DPZ_LR);
                    Program.GrammArrayPool.Return(DPX_LR);
                    Program.GrammArrayPool.Return(DPYj_LR);
                });
            }

            if (Program.INUMS == 6)
            {
                int range_parallel = NI / Program.pOptions.MaxDegreeOfParallelism - (StripeCounter % 6);
                range_parallel = Math.Max(Program.StripeWidth - (StripeCounter % 6), range_parallel); // min. Program.StripeWidth cells per processor
                StripeCounter++;
                range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
                                                               //Parallel.For(2, NI, Program.pOptions, i =>
                Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
                {
                    int NK_P = NK; int NJ_P = NJ;
                    Span<double> PIM = stackalloc double[2 * NK_P + 2];
                    Span<double> QIM = stackalloc double[2 * NK_P + 2];
                    double[] DP_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double[] DPZ_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double[] DPXi_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double[] DPY_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double TERMP = 0;
                    double DIM, RHOAIM_LL, RHOAJM_LL, RHOiAIM_LL;

                    for (int i = range.Item1; i < range.Item2; ++i)
                    {
                        int border = Math.Min(i - range.Item1, range.Item2 - 1 - i);
                        for (int j = NJ_P - 1; j >= 2; --j)
                        {
                            int iP1 = i + 1;
                            int jP1 = j + 1;
                            ReadOnlySpan<double> AB_L = Program.AB[i][j];
                            ReadOnlySpan<double> AE_L = Program.AE[i][j];
                            ReadOnlySpan<double> AN_L = Program.AN[i][j];
                            ReadOnlySpan<double> AP_L = Program.AP[i][j];
                            ReadOnlySpan<double> AS_L = Program.AS[i][j];
                            ReadOnlySpan<double> AT_L = Program.AT[i][j];
                            ReadOnlySpan<double> AW_L = Program.AW[i][j];
                            /*
                        double[] AIM_L    = Program.AIM[i][j];
                        double[] AIMi_L   = Program.AIM[iP1][j];
                        double[] AIMmj_L  = Program.AIM[i][j - 1];
            				 */
                            ReadOnlySpan<float> AIM_L = Program.AP0[i][j];
                            ReadOnlySpan<float> AIMi_L = Program.AP0[iP1][j];
                            ReadOnlySpan<float> AIMmj_L = Program.AP0[i][j - 1];

                            double[] DP_L;
                            double[] DPZ_L;
                            double[] DPXi_L;
                            double[] DPY_L;

                            if (border < 2)
                            {
                                DP_L = DP_LR;
                                DPZ_L = DPZ_LR;
                                DPXi_L = DPXi_LR;
                                DPY_L = DPY_LR;
                                Program.CopyArrayLockSource(Program.DP[i][j], DP_L);
                                Program.CopyArrayLockSource(Program.DPZ[i][j], DPZ_L);
                                Program.CopyArrayLockSource(Program.DPX[iP1][j], DPXi_L);
                                Program.CopyArrayLockSource(Program.DPY[i][j], DPY_L);
                            }
                            else
                            {
                                DP_L = Program.DP[i][j];
                                DPZ_L = Program.DPZ[i][j];
                                DPXi_L = Program.DPX[iP1][j];
                                DPY_L = Program.DPY[i][j];
                            }

                            //ReadOnlySpan<double> DP_L = Program.DP[i][j];
                            ReadOnlySpan<double> DPmj_L = Program.DP[i][j - 1];
                            ReadOnlySpan<double> DPi_L = Program.DP[iP1][j];
                            ReadOnlySpan<double> DPX_L = Program.DPX[i][j];
                            //ReadOnlySpan<double> DPXi_L = Program.DPX[iP1][j];
                            //ReadOnlySpan<double> DPY_L = Program.DPY[i][j];
                            ReadOnlySpan<double> DPYj_L = Program.DPY[i][jP1];
                            //ReadOnlySpan<double> DPZ_L = Program.DPZ[i][j];
                            ReadOnlySpan<double> DPZmj_L = Program.DPZ[i][j - 1];
                            ReadOnlySpan<double> DPZi_L = Program.DPZ[iP1][j];
                            ReadOnlySpan<float> AREAX_L = Program.AREAXImm[i][j];
                            ReadOnlySpan<float> AREAXi_L = Program.AREAXImm[iP1][j];
                            ReadOnlySpan<float> AREAY_L = Program.AREAYImm[i][j];
                            ReadOnlySpan<float> AREAYmj_L = Program.AREAYImm[i][j - 1];
                            ReadOnlySpan<float> AREAZ_L = Program.AREAZImm[i][j];
                            ReadOnlySpan<float> AREAXYZ_L = Program.AREAXYZImm[i][j];
                            ReadOnlySpan<float> AREAZX_L = Program.AREAZXImm[i][j];
                            ReadOnlySpan<float> AREAZXi_L = Program.AREAZXImm[iP1][j];
                            ReadOnlySpan<float> AREAZY_L = Program.AREAZYImm[i][j];
                            ReadOnlySpan<float> AREAZYmj_L = Program.AREAZYImm[i][j - 1];
                            ReadOnlySpan<float> RHO_L = Program.RHOImm[i][j];
                            ReadOnlySpan<float> RHOi_L = Program.RHOImm[iP1][j];
                            ReadOnlySpan<float> RHOmj_L = Program.RHOImm[i][j - 1];
                            ReadOnlySpan<float> SUXYZ_L = Program.SUXYZ[i][j];
                            ReadOnlySpan<float> SUX_L = Program.SUX[i][j];
                            ReadOnlySpan<float> SUXi_L = Program.SUX[iP1][j];
                            ReadOnlySpan<float> SUY_L = Program.SUY[i][j];
                            ReadOnlySpan<float> SUZ_L = Program.SUZ[i][j];

                            int m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; --kn)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)

                                //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = SUZ_L[k] * AREAZ_L[k] - AE_L[kn] * DPXi_L[k - 1] -
                                        AN_L[kn] * DPYj_L[k - 1] - AW_L[kn] * DPX_L[k] -
                                        AS_L[kn] * DPY_L[k];
                                    //Recurrence formula
                                    TERMP = 1 / (AP_L[kn] - AT_L[kn] * PIM[kn + 1]);
                                    PIM[kn] = AB_L[kn] * TERMP;
                                    QIM[kn] = (DIM + AT_L[kn] * QIM[kn + 1]) * TERMP;
                                    m--;
                                }
                                //Coefficients for the upper half-cell
                                else
                                {
                                    if (kn < 2 * NK_P)
                                    {
                                        DIM = SUXYZ_L[k] * AREAXYZ_L[k] + AW_L[kn] * DPX_L[k] +
                                            AS_L[kn] * DPY_L[k] + AE_L[kn] * DPXi_L[k] +
                                            AN_L[kn] * DPYj_L[k];

                                        //Recurrence formula
                                        TERMP = 1 / (AP_L[kn] - AT_L[kn] * PIM[kn + 1]);
                                        PIM[kn] = AB_L[kn] * TERMP;
                                        QIM[kn] = (DIM + AT_L[kn] * QIM[kn + 1]) * TERMP;
                                    }
                                    m = 2;
                                }
                            }
                            //Obtain new P-components
                            m = 2;
                            DP_L[0] = 0;
                            DP_L[NK_P] = 0;
                            DPZ_L[NK_P] = 0;
                            DPXi_L[NK_P] = 0;
                            DPY_L[NK_P] = 0;

                            for (int kn1 = 1; kn1 <= 2 * (NK_P - 1); ++kn1)
                            {
                                int k1 = 1 + kn1 >> 1;  // (int) (kn1 * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    DPZ_L[k1] = PIM[kn1] * DP_L[k1 - 1] + QIM[kn1];
                                    m--;
                                }
                                else
                                {
                                    DP_L[k1] = PIM[kn1] * DPZ_L[k1] + QIM[kn1];
                                    m = 2;
                                }
                            }

                            //compute neighbours
                            for (int k2 = 1; k2 <= NK_P - 1; ++k2)
                            {
                                int kn2 = 2 * k2 - 1;
                                int k2P = k2 + 1;
                                RHOAIM_LL = RHO_L[k2] / AIM_L[k2];
                                RHOAJM_LL = RHOmj_L[k2] / AIMmj_L[k2];
                                RHOiAIM_LL = RHOi_L[k2] / AIMi_L[k2];

                                DPXi_L[k2] = 1 / AREAXi_L[k2] / (RHOAIM_LL + RHOiAIM_LL) *
                                    (SUXi_L[k2] + RHOAIM_LL * (DP_L[k2] * (AREAX_L[k2] +
                                                                           AREAZX_L[k2]) - DPZ_L[k2P] * AREAZX_L[k2P]) - RHOiAIM_LL
                                     * (DPZi_L[k2] * AREAZXi_L[k2] - DPi_L[k2] * (AREAXi_L[k2] + AREAZXi_L[k2])));

                                DPY_L[k2] = 1 / AREAY_L[k2] / (RHOAJM_LL + RHOAIM_LL) *
                                    (SUY_L[k2] + RHOAJM_LL * (DPmj_L[k2] * (AREAYmj_L[k2] +
                                                                            AREAZYmj_L[k2]) - DPZmj_L[k2P] * AREAZYmj_L[k2P]) - RHOAIM_LL
                                     * (DPZ_L[k2] * AREAZY_L[k2] - DP_L[k2] * (AREAY_L[k2] + AREAZY_L[k2])));
                            }
                            if (border < 2)
                            {
                                Program.CopyArrayLockDest(DP_L, Program.DP[i][j]);
                                Program.CopyArrayLockDest(DPZ_L, Program.DPZ[i][j]);
                                Program.CopyArrayLockDest(DPXi_L, Program.DPX[iP1][j]);
                                Program.CopyArrayLockDest(DPY_L, Program.DPY[i][j]);
                            }
                        }
                    }
                    Program.GrammArrayPool.Return(DP_LR);
                    Program.GrammArrayPool.Return(DPZ_LR);
                    Program.GrammArrayPool.Return(DPXi_LR);
                    Program.GrammArrayPool.Return(DPY_LR);
                });
            }

            if (Program.INUMS == 7)
            {
                int range_parallel = NJ / Program.pOptions.MaxDegreeOfParallelism - (StripeCounter % 6);
                range_parallel = Math.Max(Program.StripeWidth - (StripeCounter % 6), range_parallel); // min. Program.StripeWidth cells per processor
                StripeCounter++;
                range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
                                                               //Parallel.For(2, NJ, Program.pOptions, j1 =>
                Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
                {
                    int NK_P = NK; int NI_P = NI;
                    Span<double> PIM = stackalloc double[2 * NK_P + 2];
                    Span<double> QIM = stackalloc double[2 * NK_P + 2];
                    double[] DP_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double[] DPZ_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double[] DPXi_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double[] DPY_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double TERMP = 0;
                    double DIM, RHOAIM_LL, RHOAJM_LL, RHOiAIM_LL;

                    for (int j = range.Item2 - 1; j >= range.Item1; --j)
                    {
                        int border = Math.Min(j - range.Item1, range.Item2 - 1 - j);
                        for (int i = 2; i <= NI_P - 1; ++i)
                        {
                            int iP1 = i + 1;
                            int jP1 = j + 1;
                            ReadOnlySpan<double> AB_L = Program.AB[i][j];
                            ReadOnlySpan<double> AE_L = Program.AE[i][j];
                            ReadOnlySpan<double> AN_L = Program.AN[i][j];
                            ReadOnlySpan<double> AP_L = Program.AP[i][j];
                            ReadOnlySpan<double> AS_L = Program.AS[i][j];
                            ReadOnlySpan<double> AT_L = Program.AT[i][j];
                            ReadOnlySpan<double> AW_L = Program.AW[i][j];
                            /*
                            double[] AIM_L    = Program.AIM[i][j];
                            double[] AIMi_L   = Program.AIM[iP1][j];
                            double[] AIMmj_L  = Program.AIM[i][j - 1];
                             */
                            ReadOnlySpan<float> AIM_L = Program.AP0[i][j];
                            ReadOnlySpan<float> AIMi_L = Program.AP0[iP1][j];
                            ReadOnlySpan<float> AIMmj_L = Program.AP0[i][j - 1];

                            double[] DP_L;
                            double[] DPZ_L;
                            double[] DPXi_L;
                            double[] DPY_L;

                            if (border < 2)
                            {
                                DP_L = DP_LR;
                                DPZ_L = DPZ_LR;
                                DPXi_L = DPXi_LR;
                                DPY_L = DPY_LR;
                                Program.CopyArrayLockSource(Program.DP[i][j], DP_L);
                                Program.CopyArrayLockSource(Program.DPZ[i][j], DPZ_L);
                                Program.CopyArrayLockSource(Program.DPX[iP1][j], DPXi_L);
                                Program.CopyArrayLockSource(Program.DPY[i][j], DPY_L);
                            }
                            else
                            {
                                DP_L = Program.DP[i][j];
                                DPZ_L = Program.DPZ[i][j];
                                DPXi_L = Program.DPX[iP1][j];
                                DPY_L = Program.DPY[i][j];
                            }

                            //ReadOnlySpan<double> DP_L = Program.DP[i][j];
                            ReadOnlySpan<double> DPmj_L = Program.DP[i][j - 1];
                            ReadOnlySpan<double> DPi_L = Program.DP[iP1][j];
                            ReadOnlySpan<double> DPX_L = Program.DPX[i][j];
                            //ReadOnlySpan<double> DPXi_L = Program.DPX[iP1][j];
                            //ReadOnlySpan<double> DPY_L = Program.DPY[i][j];
                            ReadOnlySpan<double> DPYj_L = Program.DPY[i][jP1];
                            //ReadOnlySpan<double> DPZ_L = Program.DPZ[i][j];
                            ReadOnlySpan<double> DPZi_L = Program.DPZ[iP1][j];
                            ReadOnlySpan<double> DPZmj_L = Program.DPZ[i][j - 1];
                            ReadOnlySpan<float> AREAX_L = Program.AREAXImm[i][j];
                            ReadOnlySpan<float> AREAXi_L = Program.AREAXImm[iP1][j];
                            ReadOnlySpan<float> AREAY_L = Program.AREAYImm[i][j];
                            ReadOnlySpan<float> AREAYmj_L = Program.AREAYImm[i][j - 1];
                            ReadOnlySpan<float> AREAZ_L = Program.AREAZImm[i][j];
                            ReadOnlySpan<float> AREAXYZ_L = Program.AREAXYZImm[i][j];
                            ReadOnlySpan<float> AREAZX_L = Program.AREAZXImm[i][j];
                            ReadOnlySpan<float> AREAZXi_L = Program.AREAZXImm[iP1][j];
                            ReadOnlySpan<float> AREAZY_L = Program.AREAZYImm[i][j];
                            ReadOnlySpan<float> AREAZYmj_L = Program.AREAZYImm[i][j - 1];
                            ReadOnlySpan<float> RHO_L = Program.RHOImm[i][j];
                            ReadOnlySpan<float> RHOi_L = Program.RHOImm[iP1][j];
                            ReadOnlySpan<float> RHOmj_L = Program.RHOImm[i][j - 1];
                            ReadOnlySpan<float> SUXYZ_L = Program.SUXYZ[i][j];
                            ReadOnlySpan<float> SUX_L = Program.SUX[i][j];
                            ReadOnlySpan<float> SUXi_L = Program.SUX[iP1][j];
                            ReadOnlySpan<float> SUY_L = Program.SUY[i][j];
                            ReadOnlySpan<float> SUZ_L = Program.SUZ[i][j];

                            int m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; --kn)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)

                                //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = SUZ_L[k] * AREAZ_L[k] - AE_L[kn] * DPXi_L[k - 1] -
                                        AN_L[kn] * DPYj_L[k - 1] - AW_L[kn] * DPX_L[k] -
                                        AS_L[kn] * DPY_L[k];
                                    //Recurrence formula
                                    TERMP = 1 / (AP_L[kn] - AT_L[kn] * PIM[kn + 1]);
                                    PIM[kn] = AB_L[kn] * TERMP;
                                    QIM[kn] = (DIM + AT_L[kn] * QIM[kn + 1]) * TERMP;
                                    m--;
                                }
                                //Coefficients for the upper half-cell
                                else
                                {
                                    if (kn < 2 * NK_P)
                                    {
                                        DIM = SUXYZ_L[k] * AREAXYZ_L[k] + AW_L[kn] * DPX_L[k] +
                                            AS_L[kn] * DPY_L[k] + AE_L[kn] * DPXi_L[k] +
                                            AN_L[kn] * DPYj_L[k];
                                        //Recurrence formula
                                        TERMP = 1 / (AP_L[kn] - AT_L[kn] * PIM[kn + 1]);
                                        PIM[kn] = AB_L[kn] * TERMP;
                                        QIM[kn] = (DIM + AT_L[kn] * QIM[kn + 1]) * TERMP;
                                    }
                                    m = 2;
                                }
                            }
                            //Obtain new P-components
                            m = 2;
                            DP_L[0] = 0;
                            DP_L[NK_P] = 0;
                            DPZ_L[NK_P] = 0;
                            DPXi_L[NK_P] = 0;
                            DPY_L[NK_P] = 0;

                            for (int kn1 = 1; kn1 <= 2 * (NK_P - 1); ++kn1)
                            {
                                int k1 = 1 + kn1 >> 1;  // (int) (kn1 * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    DPZ_L[k1] = PIM[kn1] * DP_L[k1 - 1] + QIM[kn1];
                                    m--;
                                }
                                else
                                {
                                    DP_L[k1] = PIM[kn1] * DPZ_L[k1] + QIM[kn1];
                                    m = 2;
                                }
                            }

                            //compute neighbours
                            for (int k2 = 1; k2 <= NK_P - 1; ++k2)
                            {
                                int kn2 = 2 * k2 - 1;
                                int k2P = k2 + 1;
                                RHOAIM_LL = RHO_L[k2] / AIM_L[k2];
                                RHOAJM_LL = RHOmj_L[k2] / AIMmj_L[k2];
                                RHOiAIM_LL = RHOi_L[k2] / AIMi_L[k2];

                                DPXi_L[k2] = 1 / AREAXi_L[k2] / (RHOAIM_LL + RHOiAIM_LL) *
                                    (SUXi_L[k2] + RHOAIM_LL * (DP_L[k2] * (AREAX_L[k2] +
                                                                           AREAZX_L[k2]) - DPZ_L[k2P] * AREAZX_L[k2P]) - RHOiAIM_LL
                                     * (DPZi_L[k2] * AREAZXi_L[k2] - DPi_L[k2] * (AREAXi_L[k2] + AREAZXi_L[k2])));


                                DPY_L[k2] = 1 / AREAY_L[k2] / (RHOAJM_LL + RHOAIM_LL) *
                                    (SUY_L[k2] + RHOAJM_LL * (DPmj_L[k2] * (AREAYmj_L[k2] +
                                                                            AREAZYmj_L[k2]) - DPZmj_L[k2P] * AREAZYmj_L[k2P]) - RHOAIM_LL
                                     * (DPZ_L[k2] * AREAZY_L[k2] - DP_L[k2] * (AREAY_L[k2] + AREAZY_L[k2])));
                            }
                            if (border < 2)
                            {
                                Program.CopyArrayLockDest(DP_L, Program.DP[i][j]);
                                Program.CopyArrayLockDest(DPZ_L, Program.DPZ[i][j]);
                                Program.CopyArrayLockDest(DPXi_L, Program.DPX[iP1][j]);
                                Program.CopyArrayLockDest(DPY_L, Program.DPY[i][j]);
                            }
                        }
                    }
                    Program.GrammArrayPool.Return(DP_LR);
                    Program.GrammArrayPool.Return(DPZ_LR);
                    Program.GrammArrayPool.Return(DPXi_LR);
                    Program.GrammArrayPool.Return(DPY_LR);
                });
            }

            if (Program.INUMS == 8)
            {
                int range_parallel = NJ / Program.pOptions.MaxDegreeOfParallelism - (StripeCounter % 6);
                range_parallel = Math.Max(Program.StripeWidth - (StripeCounter % 6), range_parallel); // min. Program.StripeWidth cells per processor
                StripeCounter++;
                range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
                //Parallel.For(2, NJ, Program.pOptions, j =>
                Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
                {
                    int NK_P = NK; int NI_P = NI;
                    Span<double> PIM = stackalloc double[2 * NK_P + 2];
                    Span<double> QIM = stackalloc double[2 * NK_P + 2];
                    double[] DP_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double[] DPZ_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double[] DPX_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double[] DPYj_LR = Program.GrammArrayPool.Rent(Program.NZ1);
                    double TERMP = 0;
                    double DIM, RHOAIM_LL, RHOAJM_LL, RHOiAIM_LL;

                    for (int j = range.Item1; j < range.Item2; ++j)
                    {
                        int border = Math.Min(j - range.Item1, range.Item2 - 1 - j);
                        for (int i = NI_P - 1; i >= 2; --i)
                        {
                            int iP1 = i + 1;
                            int jP1 = j + 1;
                            ReadOnlySpan<double> AB_L = Program.AB[i][j];
                            ReadOnlySpan<double> AE_L = Program.AE[i][j];
                            ReadOnlySpan<double> AN_L = Program.AN[i][j];
                            ReadOnlySpan<double> AP_L = Program.AP[i][j];
                            ReadOnlySpan<double> AS_L = Program.AS[i][j];
                            ReadOnlySpan<double> AT_L = Program.AT[i][j];
                            ReadOnlySpan<double> AW_L = Program.AW[i][j];
                            /*
                        double[] AIM_L    = Program.AIM[i][j];
                        double[] AIMmi_L   = Program.AIM[i - 1][j];
                        double[] AIMj_L   = Program.AIM[i][jP1];
                			 */
                            ReadOnlySpan<float> AIM_L = Program.AP0[i][j];
                            ReadOnlySpan<float> AIMmi_L = Program.AP0[i - 1][j];
                            ReadOnlySpan<float> AIMj_L = Program.AP0[i][jP1];

                            double[] DP_L;
                            double[] DPZ_L;
                            double[] DPX_L;
                            double[] DPYj_L;

                            if (border < 2)
                            {
                                DP_L = DP_LR;
                                DPZ_L = DPZ_LR;
                                DPX_L = DPX_LR;
                                DPYj_L = DPYj_LR;
                                Program.CopyArrayLockSource(Program.DP[i][j], DP_L);
                                Program.CopyArrayLockSource(Program.DPZ[i][j], DPZ_L);
                                Program.CopyArrayLockSource(Program.DPX[i][j], DPX_L);
                                Program.CopyArrayLockSource(Program.DPY[i][jP1], DPYj_L);
                            }
                            else
                            {
                                DP_L = Program.DP[i][j];
                                DPZ_L = Program.DPZ[i][j];
                                DPX_L = Program.DPX[i][j];
                                DPYj_L = Program.DPY[i][jP1];
                            }

                            //ReadOnlySpan<double> DP_L = Program.DP[i][j];
                            ReadOnlySpan<double> DPj_L = Program.DP[i][jP1];
                            ReadOnlySpan<double> DPmi_L = Program.DP[i - 1][j];
                            //ReadOnlySpan<double> DPX_L = Program.DPX[i][j];
                            ReadOnlySpan<double> DPXi_L = Program.DPX[iP1][j];
                            ReadOnlySpan<double> DPY_L = Program.DPY[i][j];
                            //ReadOnlySpan<double> DPYj_L = Program.DPY[i][jP1];
                            //ReadOnlySpan<double> DPZ_L = Program.DPZ[i][j];
                            ReadOnlySpan<double> DPZmi_L = Program.DPZ[i - 1][j];
                            ReadOnlySpan<double> DPZj_L = Program.DPZ[i][jP1];
                            ReadOnlySpan<float> AREAX_L = Program.AREAXImm[i][j];
                            ReadOnlySpan<float> AREAXmi_L = Program.AREAXImm[i - 1][j];
                            ReadOnlySpan<float> AREAY_L = Program.AREAYImm[i][j];
                            ReadOnlySpan<float> AREAYj_L = Program.AREAYImm[i][jP1];
                            ReadOnlySpan<float> AREAZ_L = Program.AREAZImm[i][j];
                            ReadOnlySpan<float> AREAXYZ_L = Program.AREAXYZImm[i][j];
                            ReadOnlySpan<float> AREAZX_L = Program.AREAZXImm[i][j];
                            ReadOnlySpan<float> AREAZXi_L = Program.AREAZXImm[i - 1][j];
                            ReadOnlySpan<float> AREAZY_L = Program.AREAZYImm[i][j];
                            ReadOnlySpan<float> AREAZYj_L = Program.AREAZYImm[i][jP1];
                            ReadOnlySpan<float> RHO_L = Program.RHOImm[i][j];
                            ReadOnlySpan<float> RHOmi_L = Program.RHOImm[i - 1][j];
                            ReadOnlySpan<float> RHOj_L = Program.RHOImm[i][jP1];
                            ReadOnlySpan<float> SUXYZ_L = Program.SUXYZ[i][j];
                            ReadOnlySpan<float> SUX_L = Program.SUX[i][j];
                            ReadOnlySpan<float> SUY_L = Program.SUY[i][jP1];
                            ReadOnlySpan<float> SUZ_L = Program.SUZ[i][j];

                            int m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; --kn)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)

                                //Coefficients for the lower half-cell
                                if (m == 2)
                                {
                                    DIM = SUZ_L[k] * AREAZ_L[k] - AE_L[kn] * DPXi_L[k - 1] -
                                        AN_L[kn] * DPYj_L[k - 1] - AW_L[kn] * DPX_L[k] -
                                        AS_L[kn] * DPY_L[k];

                                    //Recurrence formula
                                    TERMP = 1 / (AP_L[kn] - AT_L[kn] * PIM[kn + 1]);
                                    PIM[kn] = AB_L[kn] * TERMP;
                                    QIM[kn] = (DIM + AT_L[kn] * QIM[kn + 1]) * TERMP;
                                    m--;
                                }
                                //Coefficients for the upper half-cell
                                else
                                {
                                    if (kn < 2 * NK_P)
                                    {
                                        DIM = SUXYZ_L[k] * AREAXYZ_L[k] + AW_L[kn] * DPX_L[k] +
                                            AS_L[kn] * DPY_L[k] + AE_L[kn] * DPXi_L[k] +
                                            AN_L[kn] * DPYj_L[k];

                                        //Recurrence formula
                                        TERMP = 1 / (AP_L[kn] - AT_L[kn] * PIM[kn + 1]);
                                        PIM[kn] = AB_L[kn] * TERMP;
                                        QIM[kn] = (DIM + AT_L[kn] * QIM[kn + 1]) * TERMP;
                                    }
                                    m = 2;
                                }
                            }
                            //Obtain new P-components
                            m = 2;
                            DP_L[0] = 0;
                            DP_L[NK_P] = 0;
                            DPZ_L[NK_P] = 0;
                            DPX_L[NK_P] = 0;
                            DPYj_L[NK_P] = 0;

                            for (int kn1 = 1; kn1 <= 2 * (NK_P - 1); ++kn1)
                            {
                                int k1 = 1 + kn1 >> 1;  // (int) (kn1 * 0.5F + 0.5F)
                                if (m == 2)
                                {
                                    DPZ_L[k1] = PIM[kn1] * DP_L[k1 - 1] + QIM[kn1];
                                    m--;
                                }
                                else
                                {
                                    DP_L[k1] = PIM[kn1] * DPZ_L[k1] + QIM[kn1];
                                    m = 2;
                                }
                            }

                            //compute neighbours
                            for (int k2 = 1; k2 <= NK_P - 1; ++k2)
                            {
                                int kn2 = 2 * k2 - 1;
                                int k2P = k2 + 1;
                                RHOAIM_LL = RHO_L[k2] / AIM_L[k2];
                                RHOAJM_LL = RHOj_L[k2] / AIMj_L[k2];
                                RHOiAIM_LL = RHOmi_L[k2] / AIMmi_L[k2];

                                DPX_L[k2] = 1 / AREAX_L[k2] / (RHOiAIM_LL + RHOAIM_LL) *
                                    (SUX_L[k2] + RHOiAIM_LL * (DPmi_L[k2] * (AREAXmi_L[k2] +
                                                                             AREAZXi_L[k2]) - DPZmi_L[k2P] * AREAZXi_L[k2P]) - RHOAIM_LL
                                     * (DPZ_L[k2] * AREAZX_L[k2] - DP_L[k2] * (AREAX_L[k2] + AREAZX_L[k2])));

                                DPYj_L[k2] = 1 / AREAYj_L[k2] / (RHOAIM_LL + RHOAJM_LL) *
                                    (SUY_L[k2] + RHOAIM_LL * (DP_L[k2] * (AREAY_L[k2] +
                                                                          AREAZY_L[k2]) - DPZ_L[k2P] * AREAZY_L[k2P]) - RHOAJM_LL
                                     * (DPZj_L[k2] * AREAZYj_L[k2] - DPj_L[k2] * (AREAYj_L[k2] + AREAZYj_L[k2])));
                            }
                            if (border < 2)
                            {
                                Program.CopyArrayLockDest(DP_L, Program.DP[i][j]);
                                Program.CopyArrayLockDest(DPZ_L, Program.DPZ[i][j]);
                                Program.CopyArrayLockDest(DPX_L, Program.DPX[i][j]);
                                Program.CopyArrayLockDest(DPYj_L, Program.DPY[i][jP1]);
                            }
                        }
                    }
                    Program.GrammArrayPool.Return(DP_LR);
                    Program.GrammArrayPool.Return(DPZ_LR);
                    Program.GrammArrayPool.Return(DPX_LR);
                    Program.GrammArrayPool.Return(DPYj_LR);
                });
            }
        }
    }
}
