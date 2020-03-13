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
        public static void Primp_calculate(int NI, int NJ, int NK)
        {
            if (Program.INUMS == 1)
            {
                int range_parallel = (int)(NI / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
                range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
                range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel

                //Parallel.For(2, NI, Program.pOptions, i =>
                Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
                {
                    int NK_P = NK; int NJ_P = NJ;
                    double[] PIM = new double[2 * NK + 2];
                    double[] QIM = new double[2 * NK + 2];
                    double TERMP = 0;
                    double DIM, RHOAIM_LL, RHOAJM_LL, RHOiAIM_LL;
                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        for (int j = 2; j <= NJ_P - 1; j++)
                        {
                            double[] AB_L = Program.AB[i][j];
                            double[] AE_L = Program.AE[i][j];
                            double[] AN_L = Program.AN[i][j];
                            double[] AP_L = Program.AP[i][j];
                            double[] AS_L = Program.AS[i][j];
                            double[] AT_L = Program.AT[i][j];
                            double[] AW_L = Program.AW[i][j];
                            /*
                        double[] AIM_L    = Program.AIM[i][j];
                        double[] AIMi_L   = Program.AIM[i + 1][j];
                        double[] AIMj_L   = Program.AIM[i][j + 1];
							 */
                            float[] AIM_L = Program.AP0[i][j];
                            float[] AIMi_L = Program.AP0[i + 1][j];
                            float[] AIMj_L = Program.AP0[i][j + 1];

                            double[] DP_L = Program.DP[i][j];
                            double[] DPi_L = Program.DP[i + 1][j];
                            double[] DPj_L = Program.DP[i][j + 1];
                            double[] DPX_L = Program.DPX[i][j];
                            double[] DPXi_L = Program.DPX[i + 1][j];
                            double[] DPY_L = Program.DPY[i][j];
                            double[] DPYj_L = Program.DPY[i][j + 1];
                            double[] DPZ_L = Program.DPZ[i][j];
                            double[] DPZi_L = Program.DPZ[i + 1][j];
                            double[] DPZj_L = Program.DPZ[i][j + 1];
                            float[] AREAX_L = Program.AREAX[i][j];
                            float[] AREAXi_L = Program.AREAX[i + 1][j];
                            float[] AREAY_L = Program.AREAY[i][j];
                            float[] AREAYj_L = Program.AREAY[i][j + 1];
                            float[] AREAZ_L = Program.AREAZ[i][j];
                            float[] AREAXYZ_L = Program.AREAXYZ[i][j];
                            float[] AREAZX_L = Program.AREAZX[i][j];
                            float[] AREAZXi_L = Program.AREAZX[i + 1][j];
                            float[] AREAZY_L = Program.AREAZY[i][j];
                            float[] AREAZYj_L = Program.AREAZY[i][j + 1];
                            float[] RHO_L = Program.RHO[i][j];
                            float[] RHOi_L = Program.RHO[i + 1][j];
                            float[] RHOj_L = Program.RHO[i][j + 1];
                            float[] SUXi_L = Program.SUX[i + 1][j];
                            float[] SUXYZ_L = Program.SUXYZ[i][j];
                            float[] SUZ_L = Program.SUZ[i][j];
                            float[] SUY_L = Program.SUY[i][j + 1];

                            int m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)

                                //Coefficients for the lower half-cell
                                if ((m - 1) == 1)
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

                            for (int kn1 = 1; kn1 <= 2 * (NK_P - 1); kn1++)
                            {
                                int k1 = 1 + kn1 >> 1;  // (int) (kn1 * 0.5F + 0.5F)

                                if ((m - 1) == 1)
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
                            for (int k2 = 1; k2 <= NK_P - 1; k2++)
                            {
                                int kn2 = 2 * k2 - 1;
                                RHOAIM_LL = RHO_L[k2] / AIM_L[k2];
                                RHOAJM_LL = RHOj_L[k2] / AIMj_L[k2];
                                RHOiAIM_LL = RHOi_L[k2] / AIMi_L[k2];

                                DPXi_L[k2] = 1 / AREAXi_L[k2] / (RHOAIM_LL + RHOiAIM_LL) *
                                    (SUXi_L[k2] + RHOAIM_LL * (DP_L[k2] * (AREAX_L[k2] +
                                                                            AREAZX_L[k2]) - DPZ_L[k2 + 1] * AREAZX_L[k2 + 1]) - RHOiAIM_LL
                                     * (DPZi_L[k2] * AREAZXi_L[k2] - DPi_L[k2] * (AREAXi_L[k2] + AREAZXi_L[k2])));


                                DPYj_L[k2] = 1 / AREAYj_L[k2] / (RHOAIM_LL + RHOAJM_LL) *
                                    (SUY_L[k2] + RHOAIM_LL * (DP_L[k2] * (AREAY_L[k2] +
                                                                          AREAZY_L[k2]) - DPZ_L[k2 + 1] * AREAZY_L[k2 + 1]) - RHOAJM_LL
                                     * (DPZj_L[k2] * AREAZYj_L[k2] - DPj_L[k2] * (AREAYj_L[k2] + AREAZYj_L[k2])));
                            }
                        }
                    }
                });
            }

            if (Program.INUMS == 2)
            {
                int range_parallel = (int)(NI / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
                range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
                range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel

                //Parallel.For(2, NI, Program.pOptions, i1 =>
                Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
                {
                    int NK_P = NK; int NJ_P = NJ;
                    double[] PIM = new double[2 * NK + 2];
                    double[] QIM = new double[2 * NK + 2];
                    double TERMP = 0;
                    double DIM, RHOAIM_LL, RHOAJM_LL, RHOiAIM_LL;
                    for (int i1 = range.Item1; i1 < range.Item2; i1++)
                    {
                        int i = NI - i1 + 1;

                        for (int j = NJ_P - 1; j >= 2; j--)
                        {
                            double[] AB_L = Program.AB[i][j];
                            double[] AE_L = Program.AE[i][j];
                            double[] AN_L = Program.AN[i][j];
                            double[] AP_L = Program.AP[i][j];
                            double[] AS_L = Program.AS[i][j];
                            double[] AT_L = Program.AT[i][j];
                            double[] AW_L = Program.AW[i][j];
                            /*
                        double[] AIM_L    = Program.AIM[i][j];
                        double[] AIMmi_L   = Program.AIM[i - 1][j];
                        double[] AIMj_L   = Program.AIM[i][j + 1];
                        double[] AIMmj_L  = Program.AIM[i][j - 1];
            				 */
                            float[] AIM_L = Program.AP0[i][j];
                            float[] AIMmi_L = Program.AP0[i - 1][j];
                            float[] AIMj_L = Program.AP0[i][j + 1];
                            float[] AIMmj_L = Program.AP0[i][j - 1];

                            double[] DP_L = Program.DP[i][j];
                            double[] DPmi_L = Program.DP[i - 1][j];
                            double[] DPmj_L = Program.DP[i][j - 1];
                            double[] DPX_L = Program.DPX[i][j];
                            double[] DPXi_L = Program.DPX[i + 1][j];
                            double[] DPY_L = Program.DPY[i][j];
                            double[] DPYj_L = Program.DPY[i][j + 1];
                            double[] DPZ_L = Program.DPZ[i][j];
                            double[] DPZmj_L = Program.DPZ[i][j - 1];
                            double[] DPZmi_L = Program.DPZ[i - 1][j];
                            float[] AREAX_L = Program.AREAX[i][j];
                            float[] AREAXmi_L = Program.AREAX[i - 1][j];
                            float[] AREAY_L = Program.AREAY[i][j];
                            float[] AREAYmj_L = Program.AREAY[i][j - 1];
                            float[] AREAZ_L = Program.AREAZ[i][j];
                            float[] AREAXYZ_L = Program.AREAXYZ[i][j];
                            float[] AREAZX_L = Program.AREAZX[i][j];
                            float[] AREAZXi_L = Program.AREAZX[i - 1][j];
                            float[] AREAZY_L = Program.AREAZY[i][j];
                            float[] AREAZYmj_L = Program.AREAZY[i][j - 1];
                            float[] RHO_L = Program.RHO[i][j];
                            float[] RHOmi_L = Program.RHO[i - 1][j];
                            float[] RHOmj_L = Program.RHO[i][j - 1];
                            float[] SUXYZ_L = Program.SUXYZ[i][j];
                            float[] SUX_L = Program.SUX[i][j];
                            float[] SUY_L = Program.SUY[i][j];
                            float[] SUZ_L = Program.SUZ[i][j];

                            int m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)

                                //Coefficients for the lower half-cell
                                if ((m - 1) == 1)
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

                            for (int kn1 = 1; kn1 <= 2 * (NK_P - 1); kn1++)
                            {
                                int k1 = 1 + kn1 >> 1;  // (int) (kn1 * 0.5F + 0.5F)

                                if ((m - 1) == 1)
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
                            for (int k2 = 1; k2 <= NK_P - 1; k2++)
                            {
                                int kn2 = 2 * k2 - 1;
                                RHOAIM_LL = RHO_L[k2] / AIM_L[k2];
                                RHOAJM_LL = RHOmj_L[k2] / AIMmj_L[k2];
                                RHOiAIM_LL = RHOmi_L[k2] / AIMmi_L[k2];

                                DPX_L[k2] = 1 / AREAX_L[k2] / (RHOiAIM_LL + RHOAIM_LL) *
                                    (SUX_L[k2] + RHOiAIM_LL * (DPmi_L[k2] * (AREAXmi_L[k2] +
                                                                             AREAZXi_L[k2]) - DPZmi_L[k2 + 1] * AREAZXi_L[k2 + 1]) - RHOAIM_LL
                                     * (DPZ_L[k2] * AREAZX_L[k2] - DP_L[k2] * (AREAX_L[k2] + AREAZX_L[k2])));


                                DPY_L[k2] = 1 / AREAY_L[k2] / (RHOAJM_LL + RHOAIM_LL) *
                                    (SUY_L[k2] + RHOAJM_LL * (DPmj_L[k2] * (AREAYmj_L[k2] +
                                                                            AREAZYmj_L[k2]) - DPZmj_L[k2 + 1] * AREAZYmj_L[k2 + 1]) - RHOAIM_LL
                                     * (DPZ_L[k2] * AREAZY_L[k2] - DP_L[k2] * (AREAY_L[k2] + AREAZY_L[k2])));
                            }
                        }
                    }
                });
            }

            if (Program.INUMS == 3)
            {

                int range_parallel = (int)(NJ / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
                range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
                range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
                                                               //Parallel.For(2, NJ, Program.pOptions, j1 =>
                Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
                {
                    int NK_P = NK; int NI_P = NI;
                    double[] PIM = new double[2 * NK_P + 2];
                    double[] QIM = new double[2 * NK_P + 2];
                    double TERMP = 0;
                    double DIM, RHOAIM_LL, RHOAJM_LL, RHOiAIM_LL;
                    for (int j1 = range.Item1; j1 < range.Item2; j1++)
                    {
                        int j = NJ - j1 + 1;
                        for (int i = NI_P - 1; i >= 2; i--)
                        {
                            double[] AB_L = Program.AB[i][j];
                            double[] AE_L = Program.AE[i][j];
                            double[] AN_L = Program.AN[i][j];
                            double[] AP_L = Program.AP[i][j];
                            double[] AS_L = Program.AS[i][j];
                            double[] AT_L = Program.AT[i][j];
                            double[] AW_L = Program.AW[i][j];
                            /*
                        double[] AIM_L    = Program.AIM[i][j];
                        double[] AIMmi_L   = Program.AIM[i - 1][j];
                        double[] AIMj_L   = Program.AIM[i][j + 1];
                        double[] AIMmj_L  = Program.AIM[i][j - 1];
            				 */
                            float[] AIM_L = Program.AP0[i][j];
                            float[] AIMmi_L = Program.AP0[i - 1][j];
                            float[] AIMj_L = Program.AP0[i][j + 1];
                            float[] AIMmj_L = Program.AP0[i][j - 1];

                            double[] DP_L = Program.DP[i][j];
                            double[] DPmi_L = Program.DP[i - 1][j];
                            double[] DPmj_L = Program.DP[i][j - 1];
                            double[] DPX_L = Program.DPX[i][j];
                            double[] DPXi_L = Program.DPX[i + 1][j];
                            double[] DPY_L = Program.DPY[i][j];
                            double[] DPYj_L = Program.DPY[i][j + 1];
                            double[] DPZ_L = Program.DPZ[i][j];
                            double[] DPZmi_L = Program.DPZ[i - 1][j];
                            double[] DPZmj_L = Program.DPZ[i][j - 1];
                            float[] AREAX_L = Program.AREAX[i][j];
                            float[] AREAXmi_L = Program.AREAX[i - 1][j];
                            float[] AREAY_L = Program.AREAY[i][j];
                            float[] AREAYmj_L = Program.AREAY[i][j - 1];
                            float[] AREAZ_L = Program.AREAZ[i][j];
                            float[] AREAXYZ_L = Program.AREAXYZ[i][j];
                            float[] AREAZX_L = Program.AREAZX[i][j];
                            float[] AREAZXi_L = Program.AREAZX[i - 1][j];
                            float[] AREAZY_L = Program.AREAZY[i][j];
                            float[] AREAZYmj_L = Program.AREAZY[i][j - 1];
                            float[] RHO_L = Program.RHO[i][j];
                            float[] RHOmi_L = Program.RHO[i - 1][j];
                            float[] RHOmj_L = Program.RHO[i][j - 1];
                            float[] SUXYZ_L = Program.SUXYZ[i][j];
                            float[] SUX_L = Program.SUX[i][j];
                            float[] SUY_L = Program.SUY[i][j];
                            float[] SUZ_L = Program.SUZ[i][j];

                            int m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)

                                //Coefficients for the lower half-cell
                                if ((m - 1) == 1)
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

                            for (int kn1 = 1; kn1 <= 2 * (NK_P - 1); kn1++)
                            {
                                int k1 = 1 + kn1 >> 1;  // (int) (kn1 * 0.5F + 0.5F)

                                if ((m - 1) == 1)
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
                            for (int k2 = 1; k2 <= NK_P - 1; k2++)
                            {
                                int kn2 = 2 * k2 - 1;
                                RHOAIM_LL = RHO_L[k2] / AIM_L[k2];
                                RHOAJM_LL = RHOmj_L[k2] / AIMmj_L[k2];
                                RHOiAIM_LL = RHOmi_L[k2] / AIMmi_L[k2];

                                DPX_L[k2] = 1 / AREAX_L[k2] / (RHOiAIM_LL + RHOAIM_LL) *
                                    (SUX_L[k2] + RHOiAIM_LL * (DPmi_L[k2] * (AREAXmi_L[k2] +
                                                                             AREAZXi_L[k2]) - DPZmi_L[k2 + 1] * AREAZXi_L[k2 + 1]) - RHOAIM_LL
                                     * (DPZ_L[k2] * AREAZX_L[k2] - DP_L[k2] * (AREAX_L[k2] + AREAZX_L[k2])));


                                DPY_L[k2] = 1 / AREAY_L[k2] / (RHOAJM_LL + RHOAIM_LL) *
                                    (SUY_L[k2] + RHOAJM_LL * (DPmj_L[k2] * (AREAYmj_L[k2] +
                                                                            AREAZYmj_L[k2]) - DPZmj_L[k2 + 1] * AREAZYmj_L[k2 + 1]) - RHOAIM_LL
                                     * (DPZ_L[k2] * AREAZY_L[k2] - DP_L[k2] * (AREAY_L[k2] + AREAZY_L[k2])));
                            }
                        }
                    }
                });
            }

            if (Program.INUMS == 4)
            {
                int range_parallel = (int)(NJ / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
                range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
                range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel

                //Parallel.For(2, NJ, Program.pOptions, j =>
                Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
                {
                    int NK_P = NK; int NI_P = NI;
                    double[] PIM = new double[2 * NK_P + 2];
                    double[] QIM = new double[2 * NK_P + 2];
                    double TERMP = 0;
                    double DIM, RHOAIM_LL, RHOAJM_LL, RHOiAIM_LL;

                    for (int j = range.Item1; j < range.Item2; j++)
                    {
                        for (int i = 2; i <= NI_P - 1; i++)
                        {
                            double[] AB_L = Program.AB[i][j];
                            double[] AE_L = Program.AE[i][j];
                            double[] AN_L = Program.AN[i][j];
                            double[] AP_L = Program.AP[i][j];
                            double[] AS_L = Program.AS[i][j];
                            double[] AT_L = Program.AT[i][j];
                            double[] AW_L = Program.AW[i][j];
                            /*
                        double[] AIM_L    = Program.AIM[i][j];
                        double[] AIMi_L   = Program.AIM[i + 1][j];
                        double[] AIMj_L   = Program.AIM[i][j + 1];
            				 */
                            float[] AIM_L = Program.AP0[i][j];
                            float[] AIMi_L = Program.AP0[i + 1][j];
                            float[] AIMj_L = Program.AP0[i][j + 1];

                            double[] DP_L = Program.DP[i][j];
                            double[] DPi_L = Program.DP[i + 1][j];
                            double[] DPj_L = Program.DP[i][j + 1];
                            double[] DPX_L = Program.DPX[i][j];
                            double[] DPXi_L = Program.DPX[i + 1][j];
                            double[] DPY_L = Program.DPY[i][j];
                            double[] DPYj_L = Program.DPY[i][j + 1];
                            double[] DPZ_L = Program.DPZ[i][j];
                            double[] DPZi_L = Program.DPZ[i + 1][j];
                            double[] DPZj_L = Program.DPZ[i][j + 1];
                            float[] AREAX_L = Program.AREAX[i][j];
                            float[] AREAXi_L = Program.AREAX[i + 1][j];
                            float[] AREAY_L = Program.AREAY[i][j];
                            float[] AREAYj_L = Program.AREAY[i][j + 1];
                            float[] AREAZ_L = Program.AREAZ[i][j];
                            float[] AREAXYZ_L = Program.AREAXYZ[i][j];
                            float[] AREAZX_L = Program.AREAZX[i][j];
                            float[] AREAZXi_L = Program.AREAZX[i + 1][j];
                            float[] AREAZY_L = Program.AREAZY[i][j];
                            float[] AREAZYj_L = Program.AREAZY[i][j + 1];
                            float[] RHO_L = Program.RHO[i][j];
                            float[] RHOi_L = Program.RHO[i + 1][j];
                            float[] RHOj_L = Program.RHO[i][j + 1];
                            float[] SUXi_L = Program.SUX[i + 1][j];
                            float[] SUXYZ_L = Program.SUXYZ[i][j];
                            float[] SUZ_L = Program.SUZ[i][j];
                            float[] SUY_L = Program.SUY[i][j + 1];

                            int m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)

                                //Coefficients for the lower half-cell
                                if ((m - 1) == 1)
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

                            for (int kn1 = 1; kn1 <= 2 * (NK_P - 1); kn1++)
                            {
                                int k1 = 1 + kn1 >> 1;  // (int) (kn1 * 0.5F + 0.5F)

                                if ((m - 1) == 1)
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
                            for (int k2 = 1; k2 <= NK_P - 1; k2++)
                            {
                                int kn2 = 2 * k2 - 1;
                                RHOAIM_LL = RHO_L[k2] / AIM_L[k2];
                                RHOAJM_LL = RHOj_L[k2] / AIMj_L[k2];
                                RHOiAIM_LL = RHOi_L[k2] / AIMi_L[k2];

                                DPXi_L[k2] = 1 / AREAXi_L[k2] / (RHOAIM_LL + RHOiAIM_LL) *
                                    (SUXi_L[k2] + RHOAIM_LL * (DP_L[k2] * (AREAX_L[k2] +
                                                                           AREAZX_L[k2]) - DPZ_L[k2 + 1] * AREAZX_L[k2 + 1]) - RHOiAIM_LL
                                     * (DPZi_L[k2] * AREAZXi_L[k2] - DPi_L[k2] * (AREAXi_L[k2] + AREAZXi_L[k2])));


                                DPYj_L[k2] = 1 / AREAYj_L[k2] / (RHOAIM_LL + RHOAJM_LL) *
                                    (SUY_L[k2] + RHOAIM_LL * (DP_L[k2] * (AREAY_L[k2] +
                                                                          AREAZY_L[k2]) - DPZ_L[k2 + 1] * AREAZY_L[k2 + 1]) - RHOAJM_LL
                                     * (DPZj_L[k2] * AREAZYj_L[k2] - DPj_L[k2] * (AREAYj_L[k2] + AREAZYj_L[k2])));
                            }
                        }
                    }
                });
            }

            if (Program.INUMS == 5)
            {
                int range_parallel = (int)(NI / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
                range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
                range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
                                                               //  Parallel.For(2, NI, Program.pOptions, i1 =>
                Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
                {
                    int NK_P = NK; int NJ_P = NJ;
                    double[] PIM = new double[2 * NK_P + 2];
                    double[] QIM = new double[2 * NK_P + 2];
                    double TERMP = 0;
                    double DIM, RHOAIM_LL, RHOAJM_LL, RHOiAIM_LL;

                    for (int i1 = range.Item1; i1 < range.Item2; i1++)
                    {
                        int i = NI - i1 + 1;

                        for (int j = 2; j <= NJ_P - 1; j++)
                        {
                            double[] AB_L = Program.AB[i][j];
                            double[] AE_L = Program.AE[i][j];
                            double[] AN_L = Program.AN[i][j];
                            double[] AP_L = Program.AP[i][j];
                            double[] AS_L = Program.AS[i][j];
                            double[] AT_L = Program.AT[i][j];
                            double[] AW_L = Program.AW[i][j];
                             /*
                         double[] AIM_L    = Program.AIM[i][j];
                         double[] AIMmi_L   = Program.AIM[i - 1][j];
                         double[] AIMj_L   = Program.AIM[i][j + 1];
                              */
                            float[] AIM_L = Program.AP0[i][j];
                            float[] AIMmi_L = Program.AP0[i - 1][j];
                            float[] AIMj_L = Program.AP0[i][j + 1];

                            double[] DP_L = Program.DP[i][j];
                            double[] DPj_L = Program.DP[i][j + 1];
                            double[] DPmi_L = Program.DP[i - 1][j];
                            double[] DPX_L = Program.DPX[i][j];
                            double[] DPXi_L = Program.DPX[i + 1][j];
                            double[] DPY_L = Program.DPY[i][j];
                            double[] DPYj_L = Program.DPY[i][j + 1];
                            double[] DPZ_L = Program.DPZ[i][j];
                            double[] DPZmi_L = Program.DPZ[i - 1][j];
                            double[] DPZj_L = Program.DPZ[i][j + 1];
                            float[] AREAX_L = Program.AREAX[i][j];
                            float[] AREAXmi_L = Program.AREAX[i - 1][j];
                            float[] AREAY_L = Program.AREAY[i][j];
                            float[] AREAYj_L = Program.AREAY[i][j + 1];
                            float[] AREAZ_L = Program.AREAZ[i][j];
                            float[] AREAXYZ_L = Program.AREAXYZ[i][j];
                            float[] AREAZX_L = Program.AREAZX[i][j];
                            float[] AREAZXi_L = Program.AREAZX[i - 1][j];
                            float[] AREAZY_L = Program.AREAZY[i][j];
                            float[] AREAZYj_L = Program.AREAZY[i][j + 1];
                            float[] RHO_L = Program.RHO[i][j];
                            float[] RHOmi_L = Program.RHO[i - 1][j];
                            float[] RHOj_L = Program.RHO[i][j + 1];
                            float[] SUXYZ_L = Program.SUXYZ[i][j];
                            float[] SUX_L = Program.SUX[i][j];
                            float[] SUY_L = Program.SUY[i][j + 1];
                            float[] SUZ_L = Program.SUZ[i][j];

                            int m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)

                                 //Coefficients for the lower half-cell
                                 if ((m - 1) == 1)
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

                            for (int kn1 = 1; kn1 <= 2 * (NK_P - 1); kn1++)
                            {
                                int k1 = 1 + kn1 >> 1;  // (int) (kn1 * 0.5F + 0.5F)

                                 if ((m - 1) == 1)
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
                             for (int k2 = 1; k2 <= NK_P - 1; k2++)
                            {
                                int kn2 = 2 * k2 - 1;
                                RHOAIM_LL = RHO_L[k2] / AIM_L[k2];
                                RHOAJM_LL = RHOj_L[k2] / AIMj_L[k2];
                                RHOiAIM_LL = RHOmi_L[k2] / AIMmi_L[k2];

                                DPX_L[k2] = 1 / AREAX_L[k2] / (RHOiAIM_LL + RHOAIM_LL) *
                                    (SUX_L[k2] + RHOiAIM_LL * (DPmi_L[k2] * (AREAXmi_L[k2] +
                                                                             AREAZXi_L[k2]) - DPZmi_L[k2 + 1] * AREAZXi_L[k2 + 1]) - RHOAIM_LL
                                     * (DPZ_L[k2] * AREAZX_L[k2] - DP_L[k2] * (AREAX_L[k2] + AREAZX_L[k2])));


                                DPYj_L[k2] = 1 / AREAYj_L[k2] / (RHOAIM_LL + RHOAJM_LL) *
                                    (SUY_L[k2] + RHOAIM_LL * (DP_L[k2] * (AREAY_L[k2] +
                                                                          AREAZY_L[k2]) - DPZ_L[k2 + 1] * AREAZY_L[k2 + 1]) - RHOAJM_LL
                                     * (DPZj_L[k2] * AREAZYj_L[k2] - DPj_L[k2] * (AREAYj_L[k2] + AREAZYj_L[k2])));
                            }
                        }
                    }
                });
            }

            if (Program.INUMS == 6)
            {
                int range_parallel = (int)(NI / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
                range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
                range_parallel = Math.Min(NI, range_parallel); // if NI < range_parallel
                                                               //Parallel.For(2, NI, Program.pOptions, i =>
                Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
                {
                    int NK_P = NK; int NJ_P = NJ;
                    double[] PIM = new double[2 * NK_P + 2];
                    double[] QIM = new double[2 * NK_P + 2];
                    double TERMP = 0;
                    double DIM, RHOAIM_LL, RHOAJM_LL, RHOiAIM_LL;
                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        for (int j = NJ_P - 1; j >= 2; j--)
                        {
                            double[] AB_L = Program.AB[i][j];
                            double[] AE_L = Program.AE[i][j];
                            double[] AN_L = Program.AN[i][j];
                            double[] AP_L = Program.AP[i][j];
                            double[] AS_L = Program.AS[i][j];
                            double[] AT_L = Program.AT[i][j];
                            double[] AW_L = Program.AW[i][j];
                            /*
                        double[] AIM_L    = Program.AIM[i][j];
                        double[] AIMi_L   = Program.AIM[i + 1][j];
                        double[] AIMmj_L  = Program.AIM[i][j - 1];
            				 */
                            float[] AIM_L = Program.AP0[i][j];
                            float[] AIMi_L = Program.AP0[i + 1][j];
                            float[] AIMmj_L = Program.AP0[i][j - 1];

                            double[] DP_L = Program.DP[i][j];
                            double[] DPmj_L = Program.DP[i][j - 1];
                            double[] DPi_L = Program.DP[i + 1][j];
                            double[] DPX_L = Program.DPX[i][j];
                            double[] DPXi_L = Program.DPX[i + 1][j];
                            double[] DPY_L = Program.DPY[i][j];
                            double[] DPYj_L = Program.DPY[i][j + 1];
                            double[] DPZ_L = Program.DPZ[i][j];
                            double[] DPZmj_L = Program.DPZ[i][j - 1];
                            double[] DPZi_L = Program.DPZ[i + 1][j];
                            float[] AREAX_L = Program.AREAX[i][j];
                            float[] AREAXi_L = Program.AREAX[i + 1][j];
                            float[] AREAY_L = Program.AREAY[i][j];
                            float[] AREAYmj_L = Program.AREAY[i][j - 1];
                            float[] AREAZ_L = Program.AREAZ[i][j];
                            float[] AREAXYZ_L = Program.AREAXYZ[i][j];
                            float[] AREAZX_L = Program.AREAZX[i][j];
                            float[] AREAZXi_L = Program.AREAZX[i + 1][j];
                            float[] AREAZY_L = Program.AREAZY[i][j];
                            float[] AREAZYmj_L = Program.AREAZY[i][j - 1];
                            float[] RHO_L = Program.RHO[i][j];
                            float[] RHOi_L = Program.RHO[i + 1][j];
                            float[] RHOmj_L = Program.RHO[i][j - 1];
                            float[] SUXYZ_L = Program.SUXYZ[i][j];
                            float[] SUX_L = Program.SUX[i][j];
                            float[] SUXi_L = Program.SUX[i + 1][j];
                            float[] SUY_L = Program.SUY[i][j];
                            float[] SUZ_L = Program.SUZ[i][j];

                            int m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)

                                //Coefficients for the lower half-cell
                                if ((m - 1) == 1)
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

                            for (int kn1 = 1; kn1 <= 2 * (NK_P - 1); kn1++)
                            {
                                int k1 = 1 + kn1 >> 1;  // (int) (kn1 * 0.5F + 0.5F)

                                if ((m - 1) == 1)
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
                            for (int k2 = 1; k2 <= NK_P - 1; k2++)
                            {
                                int kn2 = 2 * k2 - 1;
                                RHOAIM_LL = RHO_L[k2] / AIM_L[k2];
                                RHOAJM_LL = RHOmj_L[k2] / AIMmj_L[k2];
                                RHOiAIM_LL = RHOi_L[k2] / AIMi_L[k2];

                                DPXi_L[k2] = 1 / AREAXi_L[k2] / (RHOAIM_LL + RHOiAIM_LL) *
                                    (SUXi_L[k2] + RHOAIM_LL * (DP_L[k2] * (AREAX_L[k2] +
                                                                           AREAZX_L[k2]) - DPZ_L[k2 + 1] * AREAZX_L[k2 + 1]) - RHOiAIM_LL
                                     * (DPZi_L[k2] * AREAZXi_L[k2] - DPi_L[k2] * (AREAXi_L[k2] + AREAZXi_L[k2])));


                                DPY_L[k2] = 1 / AREAY_L[k2] / (RHOAJM_LL + RHOAIM_LL) *
                                    (SUY_L[k2] + RHOAJM_LL * (DPmj_L[k2] * (AREAYmj_L[k2] +
                                                                            AREAZYmj_L[k2]) - DPZmj_L[k2 + 1] * AREAZYmj_L[k2 + 1]) - RHOAIM_LL
                                     * (DPZ_L[k2] * AREAZY_L[k2] - DP_L[k2] * (AREAY_L[k2] + AREAZY_L[k2])));
                            }
                        }
                    }
                });
            }

            if (Program.INUMS == 7)
            {
                int range_parallel = (int)(NJ / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
                range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
                range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
                                                               //Parallel.For(2, NJ, Program.pOptions, j1 =>
                Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
               {
                   int NK_P = NK; int NI_P = NI;
                   double[] PIM = new double[2 * NK_P + 2];
                   double[] QIM = new double[2 * NK_P + 2];
                   double TERMP = 0;
                   double DIM, RHOAIM_LL, RHOAJM_LL, RHOiAIM_LL;

                   for (int j1 = range.Item1; j1 < range.Item2; j1++)
                   {
                       int j = NJ - j1 + 1;
                       for (int i = 2; i <= NI_P - 1; i++)
                       {
                           double[] AB_L = Program.AB[i][j];
                           double[] AE_L = Program.AE[i][j];
                           double[] AN_L = Program.AN[i][j];
                           double[] AP_L = Program.AP[i][j];
                           double[] AS_L = Program.AS[i][j];
                           double[] AT_L = Program.AT[i][j];
                           double[] AW_L = Program.AW[i][j];
                        /*
                        double[] AIM_L    = Program.AIM[i][j];
                        double[] AIMi_L   = Program.AIM[i + 1][j];
                        double[] AIMmj_L  = Program.AIM[i][j - 1];
            			 */
                           float[] AIM_L = Program.AP0[i][j];
                           float[] AIMi_L = Program.AP0[i + 1][j];
                           float[] AIMmj_L = Program.AP0[i][j - 1];

                           double[] DP_L = Program.DP[i][j];
                           double[] DPmj_L = Program.DP[i][j - 1];
                           double[] DPi_L = Program.DP[i + 1][j];
                           double[] DPX_L = Program.DPX[i][j];
                           double[] DPXi_L = Program.DPX[i + 1][j];
                           double[] DPY_L = Program.DPY[i][j];
                           double[] DPYj_L = Program.DPY[i][j + 1];
                           double[] DPZ_L = Program.DPZ[i][j];
                           double[] DPZi_L = Program.DPZ[i + 1][j];
                           double[] DPZmj_L = Program.DPZ[i][j - 1];
                           float[] AREAX_L = Program.AREAX[i][j];
                           float[] AREAXi_L = Program.AREAX[i + 1][j];
                           float[] AREAY_L = Program.AREAY[i][j];
                           float[] AREAYmj_L = Program.AREAY[i][j - 1];
                           float[] AREAZ_L = Program.AREAZ[i][j];
                           float[] AREAXYZ_L = Program.AREAXYZ[i][j];
                           float[] AREAZX_L = Program.AREAZX[i][j];
                           float[] AREAZXi_L = Program.AREAZX[i + 1][j];
                           float[] AREAZY_L = Program.AREAZY[i][j];
                           float[] AREAZYmj_L = Program.AREAZY[i][j - 1];
                           float[] RHO_L = Program.RHO[i][j];
                           float[] RHOi_L = Program.RHO[i + 1][j];
                           float[] RHOmj_L = Program.RHO[i][j - 1];
                           float[] SUXYZ_L = Program.SUXYZ[i][j];
                           float[] SUX_L = Program.SUX[i][j];
                           float[] SUXi_L = Program.SUX[i + 1][j];
                           float[] SUY_L = Program.SUY[i][j];
                           float[] SUZ_L = Program.SUZ[i][j];

                           int m = 1;
                           for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                           {
                               int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)

                            //Coefficients for the lower half-cell
                            if ((m - 1) == 1)
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

                           for (int kn1 = 1; kn1 <= 2 * (NK_P - 1); kn1++)
                           {
                               int k1 = 1 + kn1 >> 1;  // (int) (kn1 * 0.5F + 0.5F)

                            if ((m - 1) == 1)
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
                        for (int k2 = 1; k2 <= NK_P - 1; k2++)
                           {
                               int kn2 = 2 * k2 - 1;
                               RHOAIM_LL = RHO_L[k2] / AIM_L[k2];
                               RHOAJM_LL = RHOmj_L[k2] / AIMmj_L[k2];
                               RHOiAIM_LL = RHOi_L[k2] / AIMi_L[k2];

                               DPXi_L[k2] = 1 / AREAXi_L[k2] / (RHOAIM_LL + RHOiAIM_LL) *
                                   (SUXi_L[k2] + RHOAIM_LL * (DP_L[k2] * (AREAX_L[k2] +
                                                                          AREAZX_L[k2]) - DPZ_L[k2 + 1] * AREAZX_L[k2 + 1]) - RHOiAIM_LL
                                    * (DPZi_L[k2] * AREAZXi_L[k2] - DPi_L[k2] * (AREAXi_L[k2] + AREAZXi_L[k2])));


                               DPY_L[k2] = 1 / AREAY_L[k2] / (RHOAJM_LL + RHOAIM_LL) *
                                   (SUY_L[k2] + RHOAJM_LL * (DPmj_L[k2] * (AREAYmj_L[k2] +
                                                                           AREAZYmj_L[k2]) - DPZmj_L[k2 + 1] * AREAZYmj_L[k2 + 1]) - RHOAIM_LL
                                    * (DPZ_L[k2] * AREAZY_L[k2] - DP_L[k2] * (AREAY_L[k2] + AREAZY_L[k2])));
                           }
                       }
                   }
               });
            }

            if (Program.INUMS == 8)
            {
                int range_parallel = (int)(NJ / Program.pOptions.MaxDegreeOfParallelism - (ITIME % 3) * 2);
                range_parallel = Math.Max(30 - (ITIME % 3) * 2, range_parallel); // min. 30 steps per processor
                range_parallel = Math.Min(NJ, range_parallel); // if NI < range_parallel
                //Parallel.For(2, NJ, Program.pOptions, j =>
                Parallel.ForEach(Partitioner.Create(2, NJ, range_parallel), range =>
                {
                    int NK_P = NK; int NI_P = NI;
                    double[] PIM = new double[2 * NK_P + 2];
                    double[] QIM = new double[2 * NK_P + 2];
                    double TERMP = 0;
                    double DIM, RHOAIM_LL, RHOAJM_LL, RHOiAIM_LL;
                    for (int j = range.Item1; j < range.Item2; j++)
                    {
                        for (int i = NI_P - 1; i >= 2; i--)
                        {
                            double[] AB_L = Program.AB[i][j];
                            double[] AE_L = Program.AE[i][j];
                            double[] AN_L = Program.AN[i][j];
                            double[] AP_L = Program.AP[i][j];
                            double[] AS_L = Program.AS[i][j];
                            double[] AT_L = Program.AT[i][j];
                            double[] AW_L = Program.AW[i][j];
                            /*
                        double[] AIM_L    = Program.AIM[i][j];
                        double[] AIMmi_L   = Program.AIM[i - 1][j];
                        double[] AIMj_L   = Program.AIM[i][j + 1];
                			 */
                            float[] AIM_L = Program.AP0[i][j];
                            float[] AIMmi_L = Program.AP0[i - 1][j];
                            float[] AIMj_L = Program.AP0[i][j + 1];

                            double[] DP_L = Program.DP[i][j];
                            double[] DPj_L = Program.DP[i][j + 1];
                            double[] DPmi_L = Program.DP[i - 1][j];
                            double[] DPX_L = Program.DPX[i][j];
                            double[] DPXi_L = Program.DPX[i + 1][j];
                            double[] DPY_L = Program.DPY[i][j];
                            double[] DPYj_L = Program.DPY[i][j + 1];
                            double[] DPZ_L = Program.DPZ[i][j];
                            double[] DPZmi_L = Program.DPZ[i - 1][j];
                            double[] DPZj_L = Program.DPZ[i][j + 1];
                            float[] AREAX_L = Program.AREAX[i][j];
                            float[] AREAXmi_L = Program.AREAX[i - 1][j];
                            float[] AREAY_L = Program.AREAY[i][j];
                            float[] AREAYj_L = Program.AREAY[i][j + 1];
                            float[] AREAZ_L = Program.AREAZ[i][j];
                            float[] AREAXYZ_L = Program.AREAXYZ[i][j];
                            float[] AREAZX_L = Program.AREAZX[i][j];
                            float[] AREAZXi_L = Program.AREAZX[i - 1][j];
                            float[] AREAZY_L = Program.AREAZY[i][j];
                            float[] AREAZYj_L = Program.AREAZY[i][j + 1];
                            float[] RHO_L = Program.RHO[i][j];
                            float[] RHOmi_L = Program.RHO[i - 1][j];
                            float[] RHOj_L = Program.RHO[i][j + 1];
                            float[] SUXYZ_L = Program.SUXYZ[i][j];
                            float[] SUX_L = Program.SUX[i][j];
                            float[] SUY_L = Program.SUY[i][j + 1];
                            float[] SUZ_L = Program.SUZ[i][j];

                            int m = 1;
                            for (int kn = 2 * (NK_P - 1); kn >= 1; kn--)
                            {
                                int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)

                                //Coefficients for the lower half-cell
                                if ((m - 1) == 1)
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

                            for (int kn1 = 1; kn1 <= 2 * (NK_P - 1); kn1++)
                            {
                                int k1 = 1 + kn1 >> 1;  // (int) (kn1 * 0.5F + 0.5F)

                                if ((m - 1) == 1)
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
                            for (int k2 = 1; k2 <= NK_P - 1; k2++)
                            {
                                int kn2 = 2 * k2 - 1;
                                RHOAIM_LL = RHO_L[k2] / AIM_L[k2];
                                RHOAJM_LL = RHOj_L[k2] / AIMj_L[k2];
                                RHOiAIM_LL = RHOmi_L[k2] / AIMmi_L[k2];

                                DPX_L[k2] = 1 / AREAX_L[k2] / (RHOiAIM_LL + RHOAIM_LL) *
                                    (SUX_L[k2] + RHOiAIM_LL * (DPmi_L[k2] * (AREAXmi_L[k2] +
                                                                             AREAZXi_L[k2]) - DPZmi_L[k2 + 1] * AREAZXi_L[k2 + 1]) - RHOAIM_LL
                                     * (DPZ_L[k2] * AREAZX_L[k2] - DP_L[k2] * (AREAX_L[k2] + AREAZX_L[k2])));


                                DPYj_L[k2] = 1 / AREAYj_L[k2] / (RHOAIM_LL + RHOAJM_LL) *
                                    (SUY_L[k2] + RHOAIM_LL * (DP_L[k2] * (AREAY_L[k2] +
                                                                          AREAZY_L[k2]) - DPZ_L[k2 + 1] * AREAZY_L[k2 + 1]) - RHOAJM_LL
                                     * (DPZj_L[k2] * AREAZYj_L[k2] - DPj_L[k2] * (AREAYj_L[k2] + AREAZYj_L[k2])));
                            }
                        }
                    }
                });
            }
        }
    }
}
