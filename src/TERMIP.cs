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

using System.Threading.Tasks;
using System.Runtime.CompilerServices;

namespace GRAMM_2001
{
    partial class Program
    {
        /// <summary>
        /// Procedure calculating the terms for the non-hydrostatic pressure equation
        /// </summary>
        /// <param name="NI"></param>
        /// <param name="NJ"></param>
        /// <param name="NK"></param>
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        public static void TERMIPterms(int NI, int NJ, int NK)
        {
            Parallel.For(2, NI, Program.pOptions, i =>
            {
                double AIM_Rez = 0, RHO = 0, RHO_AIM = 0, AIM_RezM1 = 0, RHOM1 = 0, RHOM1_AIMM1 = 0;
                double AREAX_LL = 0, AREAY_LL = 0, AREAZX_LL = 0, AREAZY_LL = 0, AREASUMX_LL = 0, AREASUMY_LL = 0;
                for (int j = 2; j <= NJ - 1; j++)
                {
                    double[] AB_L = Program.AB[i][j];
                    double[] AE_L = Program.AE[i][j];
                    double[] AN_L = Program.AN[i][j];
                    double[] AP_L = Program.AP[i][j];
                    float[] AP0_L = Program.AP0[i][j];
                    double[] AS_L = Program.AS[i][j];
                    double[] AT_L = Program.AT[i][j];
                    double[] AW_L = Program.AW[i][j];
                    float[] AREA_L = Program.AREA[i][j];
                    float[] AREAX_L = Program.AREAXImm[i][j];
                    float[] AREAY_L = Program.AREAY[i][j];
                    float[] AREAZX_L = Program.AREAZX[i][j];
                    float[] AREAZY_L = Program.AREAZY[i][j];
                    float[] RHO_L = Program.RHO[i][j];

                    int m = 2;

                    for (int kn = 1; kn <= 2 * NK; kn++)
                    {
                        int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)

                        //As array access is a slow process, certain often-used array-values are used as constant
                        if (kn % 2 != 0) // ODD numbers 1,3,5,..- new k value - set local constants new
                        {
                            //double AIM_Rez = Program.AIM_Rez[i][j][kn];
                            AIM_Rez = 1 / AP0_L[k];
                            RHO = RHO_L[k];
                            RHO_AIM = RHO * AIM_Rez;
                            //double AIM_RezM1 = Program.AIM_Rez[i][j][kn-1];
                            AIM_RezM1 = 1 / AP0_L[k - 1];
                            RHOM1 = Program.RHO[i][j][k - 1];
                            RHOM1_AIMM1 = RHOM1 * AIM_RezM1;
                            AREAX_LL = AREAX_L[k];
                            AREAY_LL = AREAY_L[k];
                            AREAZX_LL = AREAZX_L[k];
                            AREAZY_LL = AREAZY_L[k];
                            AREASUMX_LL = AREAX_LL + AREAZX_LL;
                            AREASUMY_LL = AREAY_LL + AREAZY_LL;
                        }

                        //Assemble TDMA coefficients
                        //Half-cell below
                        if ((m - 1) == 1)
                        {
                            m = m - 1;

                            //coefficients of pressure equation
                            if (kn > 1)
                            {
                                AP_L[kn] = Pow2(AREA_L[k]) * (RHOM1_AIMM1 + RHO_AIM)
                                    + (Pow2(AREAZX_LL) + Pow2(AREAZY_LL)) * (RHOM1_AIMM1 + RHO_AIM);
                                AB_L[kn] = Pow2(AREA_L[k]) * (RHOM1_AIMM1)
                                    + (AREAZX_LL * (AREAX_L[k - 1] + AREAZX_L[k - 1]) +
                                       AREAZY_LL * (AREAY_L[k - 1] + AREAZY_L[k - 1]))
                                    * (RHOM1_AIMM1);
                                AE_L[kn] = Program.AREAXImm[i + 1][j][k - 1] * AREAZX_LL * RHOM1_AIMM1;
                                AN_L[kn] = Program.AREAY[i][j + 1][k - 1] * AREAZY_LL * RHOM1_AIMM1;
                            }
                            else
                            {
                                AP_L[kn] = (Pow2(AREA_L[k]) + Pow2(AREAZX_LL) + Pow2(AREAZY_LL)) * RHO_AIM;
                                AB_L[kn] = 0;
                                AE_L[kn] = 0;
                                AN_L[kn] = 0;
                            }
                            AT_L[kn] = Pow2(AREA_L[k]) * RHO_AIM
                                + (AREAZX_LL * AREASUMX_LL +
                                   AREAZY_LL * AREASUMY_LL)
                                * RHO_AIM;
                            AW_L[kn] = AREAX_LL * AREAZX_LL * RHO_AIM;
                            AS_L[kn] = AREAY_LL * AREAZY_LL * RHO_AIM;
                        }
                        else
                        {
                            m = 2;

                            AP_L[kn] = Pow2(AREA_L[k]) * (RHO_AIM + RHO_AIM)
                                + (Pow2(AREASUMX_LL) + Pow2(AREASUMY_LL)) *
                                (RHO_AIM + RHO_AIM);
                            AB_L[kn] = Pow2(AREA_L[k]) * (RHO_AIM)
                                + (AREAZX_LL * AREASUMX_LL +
                                   AREAZY_LL * AREASUMY_LL)
                                * RHO_AIM;
                            AT_L[kn] = Pow2(AREA_L[k]) * RHO_AIM
                                + (AREAZX_L[k + 1] * AREASUMX_LL +
                                   AREAZY_L[k + 1] * AREASUMY_LL)
                                * RHO_AIM;
                            AW_L[kn] = AREAX_LL * (AREASUMX_LL) * (RHO_AIM);
                            AS_L[kn] = AREAY_LL * (AREASUMY_LL) * (RHO_AIM);
                            AE_L[kn] = Program.AREAXImm[i + 1][j][k] * AREASUMX_LL * RHO_AIM;
                            AN_L[kn] = Program.AREAY[i][j + 1][k] * AREASUMY_LL * RHO_AIM;
                        }
                    }
                }
            });
        }
    }
}
