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
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        public static int Primp_calculate(int NI, int NJ, int NK)
        {
            // define the prsssure solver iteration number
            int iteration = 0;
            
            for (; iteration < 9; iteration++)
            {
                if (iteration == 1)
                {
                    //Parallel.For(2, NI, Program.pOptions, i =>
                    Parallel.ForEach(PartitionerI[Program.StripeCounter++ % PartitionerI.Count], range =>
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
                        ReadOnlySpan<double> DPi_L;
                        ReadOnlySpan<double> DPj_L;
                        ReadOnlySpan<double> DPX_L;
                        ReadOnlySpan<double> DPY_L;
                        ReadOnlySpan<double> DPZi_L;
                        ReadOnlySpan<double> DPZj_L;
                        ReadOnlySpan<float> AREAX_L;
                        ReadOnlySpan<float> AREAXi_L;
                        ReadOnlySpan<float> AREAY_L;
                        ReadOnlySpan<float> AREAYj_L;
                        ReadOnlySpan<float> AREAZ_L;
                        ReadOnlySpan<float> AREAXYZ_L;
                        ReadOnlySpan<float> AREAZX_L;
                        ReadOnlySpan<float> AREAZXi_L;
                        ReadOnlySpan<float> AREAZY_L;
                        ReadOnlySpan<float> AREAZYj_L;
                        ReadOnlySpan<float> RHO_L;
                        ReadOnlySpan<float> RHOi_L;
                        ReadOnlySpan<float> RHOj_L;
                        ReadOnlySpan<float> SUXi_L;
                        ReadOnlySpan<float> SUXYZ_L;
                        ReadOnlySpan<float> SUZ_L;
                        ReadOnlySpan<float> SUY_L;
                        ReadOnlySpan<double> AB_L;
                        ReadOnlySpan<double> AE_L;
                        ReadOnlySpan<double> AN_L;
                        ReadOnlySpan<double> AP_L;
                        ReadOnlySpan<double> AS_L;
                        ReadOnlySpan<double> AT_L;
                        ReadOnlySpan<double> AW_L;
                        ReadOnlySpan<float> AIM_L;
                        ReadOnlySpan<float> AIMi_L;
                        ReadOnlySpan<float> AIMj_L;

                        for (int i = range.Item1; i < range.Item2; ++i)
                        {
                            int border = Math.Min(i - range.Item1, range.Item2 - 1 - i);
                            double[][] AB_Li = Program.AB[i];
                            double[][] AE_Li = Program.AE[i];
                            double[][] AN_Li = Program.AN[i];
                            double[][] AP_Li = Program.AP[i];
                            double[][] AS_Li = Program.AS[i];
                            double[][] AT_Li = Program.AT[i];
                            double[][] AW_Li = Program.AW[i];
                            float[][] AIM_Li = Program.AP0[i];
                            float[][] AIMi_Li = Program.AP0[i + 1];
                            
                            for (int j = 2; j <= NJ_P - 1; ++j)
                            {
                                int iP1 = i + 1;
                                int jP1 = j + 1;
                                AB_L = AB_Li[j];
                                AE_L = AE_Li[j];
                                AN_L = AN_Li[j];
                                AP_L = AP_Li[j];
                                AS_L = AS_Li[j];
                                AT_L = AT_Li[j];
                                AW_L = AW_Li[j];
                                /*
                                double[] AIM_L    = Program.AIM[i][j];
                                double[] AIMi_L   = Program.AIM[iP1][j];
                                double[] AIMj_L   = Program.AIM[i][jP1];
                                 */
                                AIM_L =  AIM_Li[j];
                                AIMi_L = AIMi_Li[j];
                                AIMj_L = AIM_Li[jP1];
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
                                    FastCopy.CopyArrayLockSource(Program.DP[i][j], DP_L);
                                    FastCopy.CopyArrayLockSource(Program.DPZ[i][j], DPZ_L);
                                    FastCopy.CopyArrayLockSource(Program.DPX[iP1][j], DPXi_L);
                                    FastCopy.CopyArrayLockSource(Program.DPY[i][jP1], DPYj_L);
                                }
                                else
                                {
                                    DP_L = Program.DP[i][j];
                                    DPZ_L = Program.DPZ[i][j];
                                    DPXi_L = Program.DPX[iP1][j];
                                    DPYj_L = Program.DPY[i][jP1];
                                }
                                DPi_L       = Program.DP[iP1][j];
                                DPj_L       = Program.DP[i][jP1];
                                DPX_L       = Program.DPX[i][j];
                                DPY_L       = Program.DPY[i][j];
                                DPZi_L      = Program.DPZ[iP1][j];
                                DPZj_L      = Program.DPZ[i][jP1];
                                AREAX_L     = Program.AREAXImm[i][j].AsSpan();
                                AREAXi_L    = Program.AREAXImm[iP1][j].AsSpan();
                                AREAY_L     = Program.AREAYImm[i][j].AsSpan();
                                AREAYj_L    = Program.AREAYImm[i][jP1].AsSpan();
                                AREAZ_L     = Program.AREAZImm[i][j].AsSpan();
                                AREAXYZ_L   = Program.AREAXYZImm[i][j].AsSpan();
                                AREAZX_L    = Program.AREAZXImm[i][j].AsSpan();
                                AREAZXi_L   = Program.AREAZXImm[iP1][j].AsSpan();
                                AREAZY_L    = Program.AREAZYImm[i][j].AsSpan();
                                AREAZYj_L   = Program.AREAZYImm[i][jP1].AsSpan();
                                RHO_L       = Program.RHO[i][j];
                                RHOi_L      = Program.RHO[iP1][j].AsSpan();
                                RHOj_L      = Program.RHO[i][jP1].AsSpan();
                                SUXi_L      = Program.SUX[iP1][j];
                                SUXYZ_L     = Program.SUXYZ[i][j];
                                SUZ_L       = Program.SUZ[i][j];
                                SUY_L       = Program.SUY[i][jP1];

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
                                    FastCopy.CopyArrayLockDest(DP_L, Program.DP[i][j]);
                                    FastCopy.CopyArrayLockDest(DPZ_L, Program.DPZ[i][j]);
                                    FastCopy.CopyArrayLockDest(DPXi_L, Program.DPX[iP1][j]);
                                    FastCopy.CopyArrayLockDest(DPYj_L, Program.DPY[i][jP1]);
                                }
                            }
                        }
                        Program.GrammArrayPool.Return(DP_LR);
                        Program.GrammArrayPool.Return(DPZ_LR);
                        Program.GrammArrayPool.Return(DPXi_LR);
                        Program.GrammArrayPool.Return(DPYj_LR);
                    });
                }
                else if (iteration == 2)
                {
                    //Parallel.For(2, NI, Program.pOptions, i1 =>
                    Parallel.ForEach(PartitionerI[Program.StripeCounter++ % PartitionerI.Count], range =>
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
                        ReadOnlySpan<double> DPmi_L;
                        ReadOnlySpan<double> DPmj_L;
                        ReadOnlySpan<double> DPXi_L;
                        ReadOnlySpan<double> DPYj_L;
                        ReadOnlySpan<double> DPZmj_L;
                        ReadOnlySpan<double> DPZmi_L;
                        ReadOnlySpan<float> AREAX_L;
                        ReadOnlySpan<float> AREAXmi_L;
                        ReadOnlySpan<float> AREAY_L;
                        ReadOnlySpan<float> AREAYmj_L;
                        ReadOnlySpan<float> AREAZ_L;
                        ReadOnlySpan<float> AREAXYZ_L;
                        ReadOnlySpan<float> AREAZX_L;
                        ReadOnlySpan<float> AREAZXi_L;
                        ReadOnlySpan<float> AREAZY_L;
                        ReadOnlySpan<float> AREAZYmj_L;
                        ReadOnlySpan<float> RHO_L;
                        ReadOnlySpan<float> RHOmi_L;
                        ReadOnlySpan<float> RHOmj_L;
                        ReadOnlySpan<float> SUXYZ_L;
                        ReadOnlySpan<float> SUX_L;
                        ReadOnlySpan<float> SUY_L;
                        ReadOnlySpan<float> SUZ_L;
                        ReadOnlySpan<double> AB_L;
                        ReadOnlySpan<double> AE_L;
                        ReadOnlySpan<double> AN_L;
                        ReadOnlySpan<double> AP_L;
                        ReadOnlySpan<double> AS_L;
                        ReadOnlySpan<double> AT_L;
                        ReadOnlySpan<double> AW_L;
                        ReadOnlySpan<float> AIM_L;
                        ReadOnlySpan<float> AIMmi_L;
                        ReadOnlySpan<float> AIMj_L;
                        ReadOnlySpan<float> AIMmj_L;

                        for (int i = range.Item2 - 1; i >= range.Item1; --i)
                        {
                            int border = Math.Min(i - range.Item1, range.Item2 - 1 - i);
                            double[][] AB_Li = Program.AB[i];
                            double[][] AE_Li = Program.AE[i];
                            double[][] AN_Li = Program.AN[i];
                            double[][] AP_Li = Program.AP[i];
                            double[][] AS_Li = Program.AS[i];
                            double[][] AT_Li = Program.AT[i];
                            double[][] AW_Li = Program.AW[i];
                            float[][] AIM_Li = Program.AP0[i];
                            float[][] AIMmi_Li = Program.AP0[i - 1];
                           
                            for (int j = NJ_P - 1; j >= 2; --j)
                            {
                                int iP1 = i + 1;
                                int jP1 = j + 1;
                                AB_L = AB_Li[j];
                                AE_L = AE_Li[j];
                                AN_L = AN_Li[j];
                                AP_L = AP_Li[j];
                                AS_L = AS_Li[j];
                                AT_L = AT_Li[j];
                                AW_L = AW_Li[j];
                                /*
                                double[] AIM_L    = Program.AIM[i][j];
                                double[] AIMi_L   = Program.AIM[iP1][j];
                                double[] AIMj_L   = Program.AIM[i][jP1];
                                 */
                                AIM_L =  AIM_Li[j];
                                AIMmi_L = AIMmi_Li[j];
                                AIMj_L = AIM_Li[jP1];
                                AIMmj_L = AIM_Li[j - 1];

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
                                    FastCopy.CopyArrayLockSource(Program.DP[i][j], DP_L);
                                    FastCopy.CopyArrayLockSource(Program.DPZ[i][j], DPZ_L);
                                    FastCopy.CopyArrayLockSource(Program.DPX[i][j], DPX_L);
                                    FastCopy.CopyArrayLockSource(Program.DPY[i][j], DPY_L);
                                }
                                else
                                {
                                    DP_L = Program.DP[i][j];
                                    DPZ_L = Program.DPZ[i][j];
                                    DPX_L = Program.DPX[i][j];
                                    DPY_L = Program.DPY[i][j];
                                }
                                //double[] DP_L = Program.DP[i][j];
                                DPmi_L     = Program.DP[i - 1][j];
                                DPmj_L     = Program.DP[i][j - 1];
                                DPXi_L     = Program.DPX[iP1][j];
                                DPYj_L     = Program.DPY[i][jP1];
                                DPZmj_L    = Program.DPZ[i][j - 1];
                                DPZmi_L    = Program.DPZ[i - 1][j];
                                AREAX_L    = Program.AREAXImm[i][j].AsSpan();
                                AREAXmi_L  = Program.AREAXImm[i - 1][j].AsSpan();
                                AREAY_L    = Program.AREAYImm[i][j].AsSpan();
                                AREAYmj_L  = Program.AREAYImm[i][j - 1].AsSpan();
                                AREAZ_L    = Program.AREAZImm[i][j].AsSpan();
                                AREAXYZ_L  = Program.AREAXYZImm[i][j].AsSpan();
                                AREAZX_L   = Program.AREAZXImm[i][j].AsSpan();
                                AREAZXi_L  = Program.AREAZXImm[i - 1][j].AsSpan();
                                AREAZY_L   = Program.AREAZYImm[i][j].AsSpan();
                                AREAZYmj_L = Program.AREAZYImm[i][j - 1].AsSpan();
                                RHO_L      = Program.RHO[i][j];
                                RHOmi_L    = Program.RHO[i - 1][j].AsSpan();
                                RHOmj_L    = Program.RHO[i][j - 1].AsSpan();
                                SUXYZ_L    = Program.SUXYZ[i][j];
                                SUX_L      = Program.SUX[i][j];
                                SUY_L      = Program.SUY[i][j];
                                SUZ_L      = Program.SUZ[i][j];

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
                                    FastCopy.CopyArrayLockDest(DP_L, Program.DP[i][j]);
                                    FastCopy.CopyArrayLockDest(DPZ_L, Program.DPZ[i][j]);
                                    FastCopy.CopyArrayLockDest(DPX_L, Program.DPX[i][j]);
                                    FastCopy.CopyArrayLockDest(DPY_L, Program.DPY[i][j]);
                                }
                            }
                        }
                        Program.GrammArrayPool.Return(DP_LR);
                        Program.GrammArrayPool.Return(DPZ_LR);
                        Program.GrammArrayPool.Return(DPX_LR);
                        Program.GrammArrayPool.Return(DPY_LR);
                    });
                }
                else if (iteration == 3)
                {
                    Parallel.ForEach(PartitionerJ[Program.StripeCounter++ % PartitionerJ.Count], range =>
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
                        ReadOnlySpan<double> DPmi_L;
                        ReadOnlySpan<double> DPmj_L;
                        ReadOnlySpan<double> DPXi_L;
                        ReadOnlySpan<double> DPYj_L;
                        ReadOnlySpan<double> DPZmi_L;
                        ReadOnlySpan<double> DPZmj_L;
                        ReadOnlySpan<float> AREAX_L;
                        ReadOnlySpan<float> AREAXmi_L;
                        ReadOnlySpan<float> AREAY_L;
                        ReadOnlySpan<float> AREAYmj_L;
                        ReadOnlySpan<float> AREAZ_L;
                        ReadOnlySpan<float> AREAXYZ_L;
                        ReadOnlySpan<float> AREAZX_L;
                        ReadOnlySpan<float> AREAZXi_L;
                        ReadOnlySpan<float> AREAZY_L;
                        ReadOnlySpan<float> AREAZYmj_L;
                        ReadOnlySpan<float> RHO_L;
                        ReadOnlySpan<float> RHOmi_L;
                        ReadOnlySpan<float> RHOmj_L;
                        ReadOnlySpan<float> SUXYZ_L;
                        ReadOnlySpan<float> SUX_L;
                        ReadOnlySpan<float> SUY_L;
                        ReadOnlySpan<float> SUZ_L;
                        ReadOnlySpan<double> AB_L;
                        ReadOnlySpan<double> AE_L;
                        ReadOnlySpan<double> AN_L;
                        ReadOnlySpan<double> AP_L;
                        ReadOnlySpan<double> AS_L;
                        ReadOnlySpan<double> AT_L;
                        ReadOnlySpan<double> AW_L;
                        ReadOnlySpan<float> AIM_L;
                        ReadOnlySpan<float> AIMmi_L;
                        ReadOnlySpan<float> AIMj_L;
                        ReadOnlySpan<float> AIMmj_L;

                        for (int j = range.Item2 - 1; j >= range.Item1; --j)
                        {
                            int border = Math.Min(j - range.Item1, range.Item2 - 1 - j);
                            for (int i = NI_P - 1; i >= 2; --i)
                            {
                                int iP1 = i + 1;
                                int jP1 = j + 1;
                                AB_L = Program.AB[i][j];
                                AE_L = Program.AE[i][j];
                                AN_L = Program.AN[i][j];
                                AP_L = Program.AP[i][j];
                                AS_L = Program.AS[i][j];
                                AT_L = Program.AT[i][j];
                                AW_L = Program.AW[i][j];
                                /*
                                double[] AIM_L    = Program.AIM[i][j];
                                double[] AIMmi_L   = Program.AIM[i - 1][j];
                                double[] AIMj_L   = Program.AIM[i][jP1];
                                double[] AIMmj_L  = Program.AIM[i][j - 1];
                                 */
                                AIM_L = Program.AP0[i][j];
                                AIMmi_L = Program.AP0[i - 1][j];
                                AIMj_L = Program.AP0[i][jP1];
                                AIMmj_L = Program.AP0[i][j - 1];

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
                                    FastCopy.CopyArrayLockSource(Program.DP[i][j], DP_L);
                                    FastCopy.CopyArrayLockSource(Program.DPZ[i][j], DPZ_L);
                                    FastCopy.CopyArrayLockSource(Program.DPX[i][j], DPX_L);
                                    FastCopy.CopyArrayLockSource(Program.DPY[i][j], DPY_L);
                                }
                                else
                                {
                                    DP_L = Program.DP[i][j];
                                    DPZ_L = Program.DPZ[i][j];
                                    DPX_L = Program.DPX[i][j];
                                    DPY_L = Program.DPY[i][j];
                                }

                                DPmi_L    = Program.DP[i - 1][j];
                                DPmj_L    = Program.DP[i][j - 1];
                                DPXi_L    = Program.DPX[iP1][j];
                                DPYj_L    = Program.DPY[i][jP1];
                                DPZmi_L   = Program.DPZ[i - 1][j];
                                DPZmj_L   = Program.DPZ[i][j - 1];
                                AREAX_L    = Program.AREAXImm[i][j].AsSpan();
                                AREAXmi_L  = Program.AREAXImm[i - 1][j].AsSpan();
                                AREAY_L    = Program.AREAYImm[i][j].AsSpan();
                                AREAYmj_L  = Program.AREAYImm[i][j - 1].AsSpan();
                                AREAZ_L    = Program.AREAZImm[i][j].AsSpan();
                                AREAXYZ_L  = Program.AREAXYZImm[i][j].AsSpan();
                                AREAZX_L   = Program.AREAZXImm[i][j].AsSpan();
                                AREAZXi_L  = Program.AREAZXImm[i - 1][j].AsSpan();
                                AREAZY_L   = Program.AREAZYImm[i][j].AsSpan();
                                AREAZYmj_L = Program.AREAZYImm[i][j - 1].AsSpan();
                                RHO_L      = Program.RHO[i][j];
                                RHOmi_L    = Program.RHO[i - 1][j].AsSpan();
                                RHOmj_L    = Program.RHO[i][j - 1].AsSpan();
                                SUXYZ_L    = Program.SUXYZ[i][j];
                                SUX_L      = Program.SUX[i][j];
                                SUY_L      = Program.SUY[i][j];
                                SUZ_L      = Program.SUZ[i][j];

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
                                    FastCopy.CopyArrayLockDest(DP_L, Program.DP[i][j]);
                                    FastCopy.CopyArrayLockDest(DPZ_L, Program.DPZ[i][j]);
                                    FastCopy.CopyArrayLockDest(DPX_L, Program.DPX[i][j]);
                                    FastCopy.CopyArrayLockDest(DPY_L, Program.DPY[i][j]);
                                }
                            }
                        }
                        Program.GrammArrayPool.Return(DP_LR);
                        Program.GrammArrayPool.Return(DPZ_LR);
                        Program.GrammArrayPool.Return(DPX_LR);
                        Program.GrammArrayPool.Return(DPY_LR);
                    });
                }
                else if (iteration == 4)
                {
                    Parallel.ForEach(PartitionerJ[Program.StripeCounter++ % PartitionerJ.Count], range =>
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
                        ReadOnlySpan<double> DPi_L;
                        ReadOnlySpan<double> DPj_L;
                        ReadOnlySpan<double> DPX_L;
                        ReadOnlySpan<double> DPY_L;
                        ReadOnlySpan<double> DPZi_L;
                        ReadOnlySpan<double> DPZj_L;
                        ReadOnlySpan<float> AREAX_L;
                        ReadOnlySpan<float> AREAXi_L;
                        ReadOnlySpan<float> AREAY_L;
                        ReadOnlySpan<float> AREAYj_L;
                        ReadOnlySpan<float> AREAZ_L;
                        ReadOnlySpan<float> AREAXYZ_L;
                        ReadOnlySpan<float> AREAZX_L;
                        ReadOnlySpan<float> AREAZXi_L;
                        ReadOnlySpan<float> AREAZY_L;
                        ReadOnlySpan<float> AREAZYj_L;
                        ReadOnlySpan<float> RHO_L;
                        ReadOnlySpan<float> RHOi_L;
                        ReadOnlySpan<float> RHOj_L;
                        ReadOnlySpan<float> SUXi_L;
                        ReadOnlySpan<float> SUXYZ_L;
                        ReadOnlySpan<float> SUZ_L;
                        ReadOnlySpan<float> SUY_L;
                        ReadOnlySpan<double> AB_L;
                        ReadOnlySpan<double> AE_L;
                        ReadOnlySpan<double> AN_L;
                        ReadOnlySpan<double> AP_L;
                        ReadOnlySpan<double> AS_L;
                        ReadOnlySpan<double> AT_L;
                        ReadOnlySpan<double> AW_L;
                        ReadOnlySpan<float> AIM_L;
                        ReadOnlySpan<float> AIMi_L;
                        ReadOnlySpan<float> AIMj_L;
                        
                        for (int j = range.Item1; j < range.Item2; ++j)
                        {
                            int border = Math.Min(j - range.Item1, range.Item2 - 1 - j);
                            for (int i = 2; i <= NI_P - 1; ++i)
                            {
                                int iP1 = i + 1;
                                int jP1 = j + 1;
                                AB_L = Program.AB[i][j];
                                AE_L = Program.AE[i][j];
                                AN_L = Program.AN[i][j];
                                AP_L = Program.AP[i][j];
                                AS_L = Program.AS[i][j];
                                AT_L = Program.AT[i][j];
                                AW_L = Program.AW[i][j];
                                /*
                            double[] AIM_L    = Program.AIM[i][j];
                            double[] AIMi_L   = Program.AIM[iP1][j];
                            double[] AIMj_L   = Program.AIM[i][jP1];
                                 */
                                AIM_L = Program.AP0[i][j];
                                AIMi_L = Program.AP0[iP1][j];
                                AIMj_L = Program.AP0[i][jP1];

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
                                    FastCopy.CopyArrayLockSource(Program.DP[i][j], DP_L);
                                    FastCopy.CopyArrayLockSource(Program.DPZ[i][j], DPZ_L);
                                    FastCopy.CopyArrayLockSource(Program.DPX[iP1][j], DPXi_L);
                                    FastCopy.CopyArrayLockSource(Program.DPY[i][jP1], DPYj_L);
                                }
                                else
                                {
                                    DP_L = Program.DP[i][j];
                                    DPZ_L = Program.DPZ[i][j];
                                    DPXi_L = Program.DPX[iP1][j];
                                    DPYj_L = Program.DPY[i][jP1];
                                }
                                
                                DPi_L     = Program.DP[iP1][j];
                                DPj_L     = Program.DP[i][jP1];
                                DPX_L     = Program.DPX[i][j];
                                DPY_L     = Program.DPY[i][j];
                                DPZi_L    = Program.DPZ[iP1][j];
                                DPZj_L    = Program.DPZ[i][jP1];
                                AREAX_L   = Program.AREAXImm[i][j].AsSpan();
                                AREAXi_L  = Program.AREAXImm[iP1][j].AsSpan();
                                AREAY_L   = Program.AREAYImm[i][j].AsSpan();
                                AREAYj_L  = Program.AREAYImm[i][jP1].AsSpan();
                                AREAZ_L   = Program.AREAZImm[i][j].AsSpan();
                                AREAXYZ_L = Program.AREAXYZImm[i][j].AsSpan();
                                AREAZX_L  = Program.AREAZXImm[i][j].AsSpan();
                                AREAZXi_L = Program.AREAZXImm[iP1][j].AsSpan();
                                AREAZY_L  = Program.AREAZYImm[i][j].AsSpan();
                                AREAZYj_L = Program.AREAZYImm[i][jP1].AsSpan();
                                RHO_L     = Program.RHO[i][j];
                                RHOi_L    = Program.RHO[iP1][j].AsSpan();
                                RHOj_L    = Program.RHO[i][jP1].AsSpan();
                                SUXi_L    = Program.SUX[iP1][j];
                                SUXYZ_L   = Program.SUXYZ[i][j];
                                SUZ_L     = Program.SUZ[i][j];
                                SUY_L     = Program.SUY[i][jP1];

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
                                    FastCopy.CopyArrayLockDest(DP_L, Program.DP[i][j]);
                                    FastCopy.CopyArrayLockDest(DPZ_L, Program.DPZ[i][j]);
                                    FastCopy.CopyArrayLockDest(DPXi_L, Program.DPX[iP1][j]);
                                    FastCopy.CopyArrayLockDest(DPYj_L, Program.DPY[i][jP1]);
                                }
                            }
                        }
                        Program.GrammArrayPool.Return(DP_LR);
                        Program.GrammArrayPool.Return(DPZ_LR);
                        Program.GrammArrayPool.Return(DPXi_LR);
                        Program.GrammArrayPool.Return(DPYj_LR);
                    });
                }
                else if (iteration == 5)
                {
                    Parallel.ForEach(PartitionerI[Program.StripeCounter++ % PartitionerI.Count], range =>
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
                        ReadOnlySpan<double> DPj_L;
                        ReadOnlySpan<double> DPmi_L;
                        ReadOnlySpan<double> DPXi_L;
                        ReadOnlySpan<double> DPY_L;
                        ReadOnlySpan<double> DPZmi_L;
                        ReadOnlySpan<double> DPZj_L;
                        ReadOnlySpan<float> AREAX_L;
                        ReadOnlySpan<float> AREAXmi_L;
                        ReadOnlySpan<float> AREAY_L;
                        ReadOnlySpan<float> AREAYj_L;
                        ReadOnlySpan<float> AREAZ_L;
                        ReadOnlySpan<float> AREAXYZ_L;
                        ReadOnlySpan<float> AREAZX_L;
                        ReadOnlySpan<float> AREAZXi_L;
                        ReadOnlySpan<float> AREAZY_L;
                        ReadOnlySpan<float> AREAZYj_L;
                        ReadOnlySpan<float> RHO_L;
                        ReadOnlySpan<float> RHOmi_L;
                        ReadOnlySpan<float> RHOj_L;
                        ReadOnlySpan<float> SUXYZ_L;
                        ReadOnlySpan<float> SUX_L;
                        ReadOnlySpan<float> SUY_L;
                        ReadOnlySpan<float> SUZ_L;
                        ReadOnlySpan<double> AB_L;
                        ReadOnlySpan<double> AE_L;
                        ReadOnlySpan<double> AN_L;
                        ReadOnlySpan<double> AP_L;
                        ReadOnlySpan<double> AS_L;
                        ReadOnlySpan<double> AT_L;
                        ReadOnlySpan<double> AW_L;
                        ReadOnlySpan<float> AIM_L;
                        ReadOnlySpan<float> AIMmi_L;
                        ReadOnlySpan<float> AIMj_L;
                        
                        for (int i = range.Item2 - 1; i >= range.Item1; --i)
                        {
                            int border = Math.Min(i - range.Item1, range.Item2 - 1 - i);
                            double[][] AB_Li = Program.AB[i];
                            double[][] AE_Li = Program.AE[i];
                            double[][] AN_Li = Program.AN[i];
                            double[][] AP_Li = Program.AP[i];
                            double[][] AS_Li = Program.AS[i];
                            double[][] AT_Li = Program.AT[i];
                            double[][] AW_Li = Program.AW[i];
                            float[][] AIM_Li = Program.AP0[i];
                            float[][] AIMmi_Li = Program.AP0[i - 1];
                            float[][] AIMj_Li = Program.AP0[i];
                            for (int j = 2; j <= NJ_P - 1; ++j)
                            {
                                int iP1 = i + 1;
                                int jP1 = j + 1;
                                AB_L = AB_Li[j];
                                AE_L = AE_Li[j];
                                AN_L = AN_Li[j];
                                AP_L = AP_Li[j];
                                AS_L = AS_Li[j];
                                AT_L = AT_Li[j];
                                AW_L = AW_Li[j];
                                /*
                                double[] AIM_L    = Program.AIM[i][j];
                                double[] AIMi_L   = Program.AIM[iP1][j];
                                double[] AIMj_L   = Program.AIM[i][j - 1];
                                 */
                                AIM_L =  AIM_Li[j];
                                AIMmi_L = AIMmi_Li[j];
                                AIMj_L = AIM_Li[jP1];

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
                                    FastCopy.CopyArrayLockSource(Program.DP[i][j], DP_L);
                                    FastCopy.CopyArrayLockSource(Program.DPZ[i][j], DPZ_L);
                                    FastCopy.CopyArrayLockSource(Program.DPX[i][j], DPX_L);
                                    FastCopy.CopyArrayLockSource(Program.DPY[i][jP1], DPYj_L);
                                }
                                else
                                {
                                    DP_L = Program.DP[i][j];
                                    DPZ_L = Program.DPZ[i][j];
                                    DPX_L = Program.DPX[i][j];
                                    DPYj_L = Program.DPY[i][jP1];
                                }

                                DPj_L     = Program.DP[i][jP1];
                                DPmi_L    = Program.DP[i - 1][j];
                                DPXi_L    = Program.DPX[iP1][j];
                                DPY_L     = Program.DPY[i][j];
                                DPZmi_L   = Program.DPZ[i - 1][j];
                                DPZj_L    = Program.DPZ[i][jP1];
                                AREAX_L   = Program.AREAXImm[i][j].AsSpan();
                                AREAXmi_L = Program.AREAXImm[i - 1][j].AsSpan();
                                AREAY_L   = Program.AREAYImm[i][j].AsSpan();
                                AREAYj_L  = Program.AREAYImm[i][jP1].AsSpan();
                                AREAZ_L   = Program.AREAZImm[i][j].AsSpan();
                                AREAXYZ_L = Program.AREAXYZImm[i][j].AsSpan();
                                AREAZX_L  = Program.AREAZXImm[i][j].AsSpan();
                                AREAZXi_L = Program.AREAZXImm[i - 1][j].AsSpan();
                                AREAZY_L  = Program.AREAZYImm[i][j].AsSpan();
                                AREAZYj_L = Program.AREAZYImm[i][jP1].AsSpan();
                                RHO_L     = Program.RHO[i][j];
                                RHOmi_L   = Program.RHO[i - 1][j].AsSpan();
                                RHOj_L    = Program.RHO[i][jP1].AsSpan();
                                SUXYZ_L   = Program.SUXYZ[i][j];
                                SUX_L     = Program.SUX[i][j];
                                SUY_L     = Program.SUY[i][jP1];
                                SUZ_L     = Program.SUZ[i][j];

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
                                    FastCopy.CopyArrayLockDest(DP_L, Program.DP[i][j]);
                                    FastCopy.CopyArrayLockDest(DPZ_L, Program.DPZ[i][j]);
                                    FastCopy.CopyArrayLockDest(DPX_L, Program.DPX[i][j]);
                                    FastCopy.CopyArrayLockDest(DPYj_L, Program.DPY[i][jP1]);
                                }
                            }
                        }
                        Program.GrammArrayPool.Return(DP_LR);
                        Program.GrammArrayPool.Return(DPZ_LR);
                        Program.GrammArrayPool.Return(DPX_LR);
                        Program.GrammArrayPool.Return(DPYj_LR);
                    });
                }
                else if (iteration == 6)
                {
                    //Parallel.For(2, NI, Program.pOptions, i =>
                    Parallel.ForEach(PartitionerI[Program.StripeCounter++ % PartitionerI.Count], range =>
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
                        ReadOnlySpan<double> DPmj_L;
                        ReadOnlySpan<double> DPi_L;
                        ReadOnlySpan<double> DPX_L;
                        ReadOnlySpan<double> DPYj_L;
                        ReadOnlySpan<double> DPZmj_L;
                        ReadOnlySpan<double> DPZi_L;
                        ReadOnlySpan<float> AREAX_L;
                        ReadOnlySpan<float> AREAXi_L;
                        ReadOnlySpan<float> AREAY_L;
                        ReadOnlySpan<float> AREAYmj_L;
                        ReadOnlySpan<float> AREAZ_L;
                        ReadOnlySpan<float> AREAXYZ_L;
                        ReadOnlySpan<float> AREAZX_L;
                        ReadOnlySpan<float> AREAZXi_L;
                        ReadOnlySpan<float> AREAZY_L;
                        ReadOnlySpan<float> AREAZYmj_L;
                        ReadOnlySpan<float> RHO_L;
                        ReadOnlySpan<float> RHOi_L;
                        ReadOnlySpan<float> RHOmj_L;
                        ReadOnlySpan<float> SUXYZ_L;
                        ReadOnlySpan<float> SUX_L;
                        ReadOnlySpan<float> SUXi_L;
                        ReadOnlySpan<float> SUY_L;
                        ReadOnlySpan<float> SUZ_L;
                        ReadOnlySpan<double> AB_L;
                        ReadOnlySpan<double> AE_L;
                        ReadOnlySpan<double> AN_L;
                        ReadOnlySpan<double> AP_L;
                        ReadOnlySpan<double> AS_L;
                        ReadOnlySpan<double> AT_L;
                        ReadOnlySpan<double> AW_L;
                        ReadOnlySpan<float> AIM_L;
                        ReadOnlySpan<float> AIMi_L;
                        ReadOnlySpan<float> AIMmj_L;

                        for (int i = range.Item1; i < range.Item2; ++i)
                        {
                            int border = Math.Min(i - range.Item1, range.Item2 - 1 - i);
                            double[][] AB_Li = Program.AB[i];
                            double[][] AE_Li = Program.AE[i];
                            double[][] AN_Li = Program.AN[i];
                            double[][] AP_Li = Program.AP[i];
                            double[][] AS_Li = Program.AS[i];
                            double[][] AT_Li = Program.AT[i];
                            double[][] AW_Li = Program.AW[i];
                            float[][] AIM_Li = Program.AP0[i];
                            float[][] AIMi_Li = Program.AP0[i + 1];
                            float[][] AIMj_Li = Program.AP0[i];
                            for (int j = NJ_P - 1; j >= 2; --j)
                            {
                                int iP1 = i + 1;
                                int jP1 = j + 1;
                                AB_L = AB_Li[j];
                                AE_L = AE_Li[j];
                                AN_L = AN_Li[j];
                                AP_L = AP_Li[j];
                                AS_L = AS_Li[j];
                                AT_L = AT_Li[j];
                                AW_L = AW_Li[j];
                                /*
                                double[] AIM_L    = Program.AIM[i][j];
                                double[] AIMi_L   = Program.AIM[iP1][j];
                                double[] AIMmj_L  = Program.AIM[i][j - 1];
                                 */
                                AIM_L = AIM_Li[j];
                                AIMi_L = AIMi_Li[j];
                                AIMmj_L = AIM_Li[j - 1];

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
                                    FastCopy.CopyArrayLockSource(Program.DP[i][j], DP_L);
                                    FastCopy.CopyArrayLockSource(Program.DPZ[i][j], DPZ_L);
                                    FastCopy.CopyArrayLockSource(Program.DPX[iP1][j], DPXi_L);
                                    FastCopy.CopyArrayLockSource(Program.DPY[i][j], DPY_L);
                                }
                                else
                                {
                                    DP_L = Program.DP[i][j];
                                    DPZ_L = Program.DPZ[i][j];
                                    DPXi_L = Program.DPX[iP1][j];
                                    DPY_L = Program.DPY[i][j];
                                }

                                DPmj_L    = Program.DP[i][j - 1];
                                DPi_L     = Program.DP[iP1][j];
                                DPX_L     = Program.DPX[i][j];
                                DPYj_L    = Program.DPY[i][jP1];
                                DPZmj_L   = Program.DPZ[i][j - 1];
                                DPZi_L    = Program.DPZ[iP1][j];
                                AREAX_L   = Program.AREAXImm[i][j].AsSpan();
                                AREAXi_L  = Program.AREAXImm[iP1][j].AsSpan();
                                AREAY_L   = Program.AREAYImm[i][j].AsSpan();
                                AREAYmj_L = Program.AREAYImm[i][j - 1].AsSpan();
                                AREAZ_L   = Program.AREAZImm[i][j].AsSpan();
                                AREAXYZ_L = Program.AREAXYZImm[i][j].AsSpan();
                                AREAZX_L  = Program.AREAZXImm[i][j].AsSpan();
                                AREAZXi_L = Program.AREAZXImm[iP1][j].AsSpan();
                                AREAZY_L  = Program.AREAZYImm[i][j].AsSpan();
                                AREAZYmj_L= Program.AREAZYImm[i][j - 1].AsSpan();
                                RHO_L     = Program.RHO[i][j];
                                RHOi_L    = Program.RHO[iP1][j].AsSpan();
                                RHOmj_L   = Program.RHO[i][j - 1].AsSpan();
                                SUXYZ_L   = Program.SUXYZ[i][j];
                                SUX_L     = Program.SUX[i][j];
                                SUXi_L    = Program.SUX[iP1][j];
                                SUY_L     = Program.SUY[i][j];
                                SUZ_L     = Program.SUZ[i][j];

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
                                    FastCopy.CopyArrayLockDest(DP_L, Program.DP[i][j]);
                                    FastCopy.CopyArrayLockDest(DPZ_L, Program.DPZ[i][j]);
                                    FastCopy.CopyArrayLockDest(DPXi_L, Program.DPX[iP1][j]);
                                    FastCopy.CopyArrayLockDest(DPY_L, Program.DPY[i][j]);
                                }
                            }
                        }
                        Program.GrammArrayPool.Return(DP_LR);
                        Program.GrammArrayPool.Return(DPZ_LR);
                        Program.GrammArrayPool.Return(DPXi_LR);
                        Program.GrammArrayPool.Return(DPY_LR);
                    });
                }
                else if (iteration == 7)
                {
                    //Parallel.For(2, NJ, Program.pOptions, j1 =>
                    Parallel.ForEach(PartitionerJ[Program.StripeCounter++ % PartitionerJ.Count], range =>
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
                        ReadOnlySpan<double> DPmj_L;
                        ReadOnlySpan<double> DPi_L;
                        ReadOnlySpan<double> DPX_L;
                        ReadOnlySpan<double> DPYj_L;
                        ReadOnlySpan<double> DPZi_L;
                        ReadOnlySpan<double> DPZmj_L;
                        ReadOnlySpan<float> AREAX_L;
                        ReadOnlySpan<float> AREAXi_L;
                        ReadOnlySpan<float> AREAY_L;
                        ReadOnlySpan<float> AREAYmj_L;
                        ReadOnlySpan<float> AREAZ_L;
                        ReadOnlySpan<float> AREAXYZ_L;
                        ReadOnlySpan<float> AREAZX_L;
                        ReadOnlySpan<float> AREAZXi_L;
                        ReadOnlySpan<float> AREAZY_L;
                        ReadOnlySpan<float> AREAZYmj_L;
                        ReadOnlySpan<float> RHO_L;
                        ReadOnlySpan<float> RHOi_L;
                        ReadOnlySpan<float> RHOmj_L;
                        ReadOnlySpan<float> SUXYZ_L;
                        ReadOnlySpan<float> SUX_L;
                        ReadOnlySpan<float> SUXi_L;
                        ReadOnlySpan<float> SUY_L;
                        ReadOnlySpan<float> SUZ_L;
                        ReadOnlySpan<double> AB_L;
                        ReadOnlySpan<double> AE_L;
                        ReadOnlySpan<double> AN_L;
                        ReadOnlySpan<double> AP_L;
                        ReadOnlySpan<double> AS_L;
                        ReadOnlySpan<double> AT_L;
                        ReadOnlySpan<double> AW_L;
                        ReadOnlySpan<float> AIM_L;
                        ReadOnlySpan<float> AIMi_L;
                        ReadOnlySpan<float> AIMmj_L;

                        for (int j = range.Item2 - 1; j >= range.Item1; --j)
                        {
                            int border = Math.Min(j - range.Item1, range.Item2 - 1 - j);
                            for (int i = 2; i <= NI_P - 1; ++i)
                            {
                                int iP1 = i + 1;
                                int jP1 = j + 1;
                                AB_L = Program.AB[i][j];
                                AE_L = Program.AE[i][j];
                                AN_L = Program.AN[i][j];
                                AP_L = Program.AP[i][j];
                                AS_L = Program.AS[i][j];
                                AT_L = Program.AT[i][j];
                                AW_L = Program.AW[i][j];
                                /*
                                double[] AIM_L    = Program.AIM[i][j];
                                double[] AIMi_L   = Program.AIM[iP1][j];
                                double[] AIMmj_L  = Program.AIM[i][j - 1];
                                 */
                                AIM_L = Program.AP0[i][j];
                                AIMi_L = Program.AP0[iP1][j];
                                AIMmj_L = Program.AP0[i][j - 1];

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
                                    FastCopy.CopyArrayLockSource(Program.DP[i][j], DP_L);
                                    FastCopy.CopyArrayLockSource(Program.DPZ[i][j], DPZ_L);
                                    FastCopy.CopyArrayLockSource(Program.DPX[iP1][j], DPXi_L);
                                    FastCopy.CopyArrayLockSource(Program.DPY[i][j], DPY_L);
                                }
                                else
                                {
                                    DP_L = Program.DP[i][j];
                                    DPZ_L = Program.DPZ[i][j];
                                    DPXi_L = Program.DPX[iP1][j];
                                    DPY_L = Program.DPY[i][j];
                                }

                                DPmj_L     = Program.DP[i][j - 1];
                                DPi_L      = Program.DP[iP1][j];
                                DPX_L      = Program.DPX[i][j];
                                DPYj_L     = Program.DPY[i][jP1];
                                DPZi_L     = Program.DPZ[iP1][j];
                                DPZmj_L    = Program.DPZ[i][j - 1];
                                AREAX_L    = Program.AREAXImm[i][j].AsSpan();
                                AREAXi_L   = Program.AREAXImm[iP1][j].AsSpan();
                                AREAY_L    = Program.AREAYImm[i][j].AsSpan();
                                AREAYmj_L  = Program.AREAYImm[i][j - 1].AsSpan();
                                AREAZ_L    = Program.AREAZImm[i][j].AsSpan();
                                AREAXYZ_L  = Program.AREAXYZImm[i][j].AsSpan();
                                AREAZX_L   = Program.AREAZXImm[i][j].AsSpan();
                                AREAZXi_L  = Program.AREAZXImm[iP1][j].AsSpan();
                                AREAZY_L   = Program.AREAZYImm[i][j].AsSpan();
                                AREAZYmj_L = Program.AREAZYImm[i][j - 1].AsSpan();
                                RHO_L      = Program.RHO[i][j];
                                RHOi_L     = Program.RHO[iP1][j].AsSpan();
                                RHOmj_L    = Program.RHO[i][j - 1].AsSpan();
                                SUXYZ_L    = Program.SUXYZ[i][j];
                                SUX_L      = Program.SUX[i][j];
                                SUXi_L     = Program.SUX[iP1][j];
                                SUY_L      = Program.SUY[i][j];
                                SUZ_L      = Program.SUZ[i][j];

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
                                    FastCopy.CopyArrayLockDest(DP_L, Program.DP[i][j]);
                                    FastCopy.CopyArrayLockDest(DPZ_L, Program.DPZ[i][j]);
                                    FastCopy.CopyArrayLockDest(DPXi_L, Program.DPX[iP1][j]);
                                    FastCopy.CopyArrayLockDest(DPY_L, Program.DPY[i][j]);
                                }
                            }
                        }
                        Program.GrammArrayPool.Return(DP_LR);
                        Program.GrammArrayPool.Return(DPZ_LR);
                        Program.GrammArrayPool.Return(DPXi_LR);
                        Program.GrammArrayPool.Return(DPY_LR);
                    });
                }
                else if (iteration == 8)
                {
                    //Parallel.For(2, NJ, Program.pOptions, j =>
                    Parallel.ForEach(PartitionerJ[Program.StripeCounter++ % PartitionerJ.Count], range =>
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
                        ReadOnlySpan<double> DPj_L;
                        ReadOnlySpan<double> DPmi_L;
                        ReadOnlySpan<double> DPXi_L;
                        ReadOnlySpan<double> DPY_L;
                        ReadOnlySpan<double> DPZmi_L;
                        ReadOnlySpan<double> DPZj_L;
                        ReadOnlySpan<float> AREAX_L;
                        ReadOnlySpan<float> AREAXmi_L;
                        ReadOnlySpan<float> AREAY_L;
                        ReadOnlySpan<float> AREAYj_L;
                        ReadOnlySpan<float> AREAZ_L;
                        ReadOnlySpan<float> AREAXYZ_L;
                        ReadOnlySpan<float> AREAZX_L;
                        ReadOnlySpan<float> AREAZXi_L;
                        ReadOnlySpan<float> AREAZY_L;
                        ReadOnlySpan<float> AREAZYj_L;
                        ReadOnlySpan<float> RHO_L;
                        ReadOnlySpan<float> RHOmi_L;
                        ReadOnlySpan<float> RHOj_L;
                        ReadOnlySpan<float> SUXYZ_L;
                        ReadOnlySpan<float> SUX_L;
                        ReadOnlySpan<float> SUY_L;
                        ReadOnlySpan<float> SUZ_L;
                        ReadOnlySpan<double> AB_L;
                        ReadOnlySpan<double> AE_L;
                        ReadOnlySpan<double> AN_L;
                        ReadOnlySpan<double> AP_L;
                        ReadOnlySpan<double> AS_L;
                        ReadOnlySpan<double> AT_L;
                        ReadOnlySpan<double> AW_L;
                        ReadOnlySpan<float> AIM_L;
                        ReadOnlySpan<float> AIMmi_L;
                        ReadOnlySpan<float> AIMj_L;

                        for (int j = range.Item1; j < range.Item2; ++j)
                        {
                            int border = Math.Min(j - range.Item1, range.Item2 - 1 - j);
                            for (int i = NI_P - 1; i >= 2; --i)
                            {
                                int iP1 = i + 1;
                                int jP1 = j + 1;
                                AB_L = Program.AB[i][j];
                                AE_L = Program.AE[i][j];
                                AN_L = Program.AN[i][j];
                                AP_L = Program.AP[i][j];
                                AS_L = Program.AS[i][j];
                                AT_L = Program.AT[i][j];
                                AW_L = Program.AW[i][j];
                                /*
                            double[] AIM_L    = Program.AIM[i][j];
                            double[] AIMmi_L   = Program.AIM[i - 1][j];
                            double[] AIMj_L   = Program.AIM[i][jP1];
                                 */
                                AIM_L = Program.AP0[i][j];
                                AIMmi_L = Program.AP0[i - 1][j];
                                AIMj_L = Program.AP0[i][jP1];

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
                                    FastCopy.CopyArrayLockSource(Program.DP[i][j], DP_L);
                                    FastCopy.CopyArrayLockSource(Program.DPZ[i][j], DPZ_L);
                                    FastCopy.CopyArrayLockSource(Program.DPX[i][j], DPX_L);
                                    FastCopy.CopyArrayLockSource(Program.DPY[i][jP1], DPYj_L);
                                }
                                else
                                {
                                    DP_L = Program.DP[i][j];
                                    DPZ_L = Program.DPZ[i][j];
                                    DPX_L = Program.DPX[i][j];
                                    DPYj_L = Program.DPY[i][jP1];
                                }

                                DPj_L     = Program.DP[i][jP1];
                                DPmi_L    = Program.DP[i - 1][j];
                                DPXi_L    = Program.DPX[iP1][j];
                                DPY_L     = Program.DPY[i][j];
                                DPZmi_L   = Program.DPZ[i - 1][j];
                                DPZj_L    = Program.DPZ[i][jP1];
                                AREAX_L   = Program.AREAXImm[i][j].AsSpan();
                                AREAXmi_L = Program.AREAXImm[i - 1][j].AsSpan();
                                AREAY_L   = Program.AREAYImm[i][j].AsSpan();
                                AREAYj_L  = Program.AREAYImm[i][jP1].AsSpan();
                                AREAZ_L   = Program.AREAZImm[i][j].AsSpan();
                                AREAXYZ_L = Program.AREAXYZImm[i][j].AsSpan();
                                AREAZX_L  = Program.AREAZXImm[i][j].AsSpan();
                                AREAZXi_L = Program.AREAZXImm[i - 1][j].AsSpan();
                                AREAZY_L  = Program.AREAZYImm[i][j].AsSpan();
                                AREAZYj_L = Program.AREAZYImm[i][jP1].AsSpan();
                                RHO_L     = Program.RHO[i][j];
                                RHOmi_L   = Program.RHO[i - 1][j].AsSpan();
                                RHOj_L    = Program.RHO[i][jP1].AsSpan();
                                SUXYZ_L   = Program.SUXYZ[i][j];
                                SUX_L      = Program.SUX[i][j];
                                SUY_L     = Program.SUY[i][jP1];
                                SUZ_L     = Program.SUZ[i][j];

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
                                    FastCopy.CopyArrayLockDest(DP_L, Program.DP[i][j]);
                                    FastCopy.CopyArrayLockDest(DPZ_L, Program.DPZ[i][j]);
                                    FastCopy.CopyArrayLockDest(DPX_L, Program.DPX[i][j]);
                                    FastCopy.CopyArrayLockDest(DPYj_L, Program.DPY[i][jP1]);
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
            return iteration;
        }
    }
}
