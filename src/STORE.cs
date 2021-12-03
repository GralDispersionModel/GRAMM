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
using System.Linq;
using System.Runtime.CompilerServices;
using System.Threading.Tasks;

namespace GRAMM_2001
{
    partial class Program
    {
        /// <summary>
        /// Update border cells and boundary values
        /// </summary>
        /// <param name="NI"></param>
        /// <param name="NJ"></param>
        /// <param name="NK"></param>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void STOREcalculate(int NI, int NJ, int NK)
        {

            //Boundary values for turbulent kinetic energy and dissipation
            Parallel.For(2, NI, Program.pOptions, i =>
            {
                for (int j = 2; j <= NJ - 1; j++)
                {
                    for (int k = 1; k <= NK - 1; k++)
                    {
                        if (j == 1)
                        {
                            Program.TEN[i][1][k] = Program.TEN[i][2][k];
                            Program.TEN[i][NJ][k] = Program.TEN[i][NJ - 1][k];
                            Program.DISSN[i][1][k] = Program.DISSN[i][2][k];
                            Program.DISSN[i][NJ][k] = Program.DISSN[i][NJ - 1][k];
                        }

                        if (i == 1)
                        {
                            Program.TEN[1][j][k] = Program.TEN[2][j][k];
                            Program.TEN[NI][j][k] = Program.TEN[NI - 1][j][k];
                            Program.DISSN[1][j][k] = Program.DISSN[2][j][k];
                            Program.DISSN[NI][j][k] = Program.DISSN[NI - 1][j][k];
                        }
                    }
                }
            });

            /*
            //Smoothing the velocity fields with a 5-point filter
            if (Program.ISTAT >= 2)
            {
                for (int i = 2; i <= NI - 1; i++)
                {
                    for (int j = 2; j <= NJ - 1; j++)
                    {
                        for (int k = 1; k <= NK; k++)
                        {
                            Program.U1N[i][j][k] = 0.2 * (Program.U1N[i + 1][j][k] + Program.U1N[i - 1][j][k] + Program.U1N[i][j - 1][k] + Program.U1N[i][j + 1][k] + Program.U1N[i][j][k]);
                            Program.V1N[i][j][k] = 0.2 * (Program.V1N[i + 1][j][k] + Program.V1N[i - 1][j][k] + Program.V1N[i][j - 1][k] + Program.V1N[i][j + 1][k] + Program.V1N[i][j][k]);
                            Program.U2N[i][j][k] = 0.2 * (Program.U2N[i + 1][j][k] + Program.U2N[i - 1][j][k] + Program.U2N[i][j - 1][k] + Program.U2N[i][j + 1][k] + Program.U2N[i][j][k]);
                            Program.V2N[i][j][k] = 0.2 * (Program.V2N[i + 1][j][k] + Program.V2N[i - 1][j][k] + Program.V2N[i][j - 1][k] + Program.V2N[i][j + 1][k] + Program.V2N[i][j][k]);
                        }
                    }
                }
            }
            */

            //Update old and actual values
            //Parallel.For(1, NI + 1, Program.pOptions, i =>
            Parallel.ForEach(Partitioner.Create(1, NI + 1, (int)((NI + 1) / Program.pOptions.MaxDegreeOfParallelism)), range =>
            {
                for (int i = range.Item1; i < range.Item2; i++)
                {
                    for (int j = 1; j <= NJ; j++)
                    {
                        double[] U1N_L = Program.U1N[i][j];
                        double[] U2N_L = Program.U2N[i][j];
                        double[] V1N_L = Program.V1N[i][j];
                        double[] V2N_L = Program.V2N[i][j];
                        double[] W1N_L = Program.W1N[i][j];
                        double[] W2N_L = Program.W2N[i][j];
                        double[] UN_L = Program.UN[i][j];
                        double[] VN_L = Program.VN[i][j];
                        double[] WN_L = Program.WN[i][j];
                        double[] U1_L = Program.U1[i][j];
                        double[] V1_L = Program.V1[i][j];
                        double[] W1_L = Program.W1[i][j];
                        double[] U2_L = Program.U2[i][j];
                        double[] V2_L = Program.V2[i][j];
                        double[] W2_L = Program.W2[i][j];
                        double[] U_L = Program.U[i][j];
                        double[] V_L = Program.V[i][j];
                        double[] W_L = Program.W[i][j];
                        double[] T_L = Program.T[i][j];
                        double[] TE_L = Program.TE[i][j];
                        double[] TN_L = Program.TN[i][j];
                        double[] TBZ_L = Program.TBZ[i][j];
                        double[] TEN_L = Program.TEN[i][j];
                        double[] QU_L = Program.QU[i][j];
                        double[] QUN_L = Program.QUN[i][j];
                        double[] QBZ_L = Program.QBZ[i][j];
                        double[] DISS_L = Program.DISS[i][j];
                        double[] DISSN_L = Program.DISSN[i][j];
                        double[] TABS_L = Program.TABS[i][j];
                        float[] FAC_L = Program.FAC[i][j];
                        double[] WAT_VAPN_L = Program.WAT_VAPN[i][j];
                        double[] WAT_VAP_L = Program.WAT_VAP[i][j];
                        double[] UG_L = Program.UG[i][j];
                        double[] VG_L = Program.VG[i][j];
                        float[] UG1_L = Program.UG1[i][j];
                        float[] VG1_L = Program.VG1[i][j];
                        float[] UG2_L = Program.UG2[i][j];
                        float[] VG2_L = Program.VG2[i][j];
                        ReadOnlySpan<float> ZSP_L = Program.ZSPImm[i][j].AsSpan();

                        for (int k = 1; k <= NK; k++)
                        {
                            if (i == 1) // just one times necessary
                            {
                                Program.U1N[1][1][k] = 0.5F * (Program.UN[1][2][k] + Program.UN[2][1][k]);
                                Program.U2N[1][1][k] = 0.5F * (Program.UN[1][2][k] + Program.UN[2][1][k]);
                                Program.V1N[1][1][k] = 0.5F * (Program.VN[1][2][k] + Program.VN[2][1][k]);
                                Program.V2N[1][1][k] = 0.5F * (Program.VN[1][2][k] + Program.VN[2][1][k]);
                                Program.W1N[1][1][k] = 0.5F * (Program.WN[1][2][k] + Program.WN[2][1][k]);
                                Program.W2N[1][1][k] = 0.5F * (Program.WN[1][2][k] + Program.WN[2][1][k]);

                                Program.U1N[1][NJ][k] = 0.5F * (Program.UN[1][NJ - 1][k] + Program.UN[2][NJ][k]);
                                Program.U2N[1][NJ][k] = 0.5F * (Program.UN[1][NJ - 1][k] + Program.UN[2][NJ][k]);
                                Program.V1N[1][NJ][k] = 0.5F * (Program.VN[1][NJ - 1][k] + Program.VN[2][NJ][k]);
                                Program.V2N[1][NJ][k] = 0.5F * (Program.VN[1][NJ - 1][k] + Program.VN[2][NJ][k]);
                                Program.W1N[1][NJ][k] = 0.5F * (Program.WN[1][NJ - 1][k] + Program.WN[2][NJ][k]);
                                Program.W2N[1][NJ][k] = 0.5F * (Program.WN[1][NJ - 1][k] + Program.WN[2][NJ][k]);

                                Program.U1N[NI][NJ][k] = 0.5F * (Program.UN[NI][NJ - 1][k] + Program.UN[NI - 1][NJ][k]);
                                Program.U2N[NI][NJ][k] = 0.5F * (Program.UN[NI][NJ - 1][k] + Program.UN[NI - 1][NJ][k]);
                                Program.V1N[NI][NJ][k] = 0.5F * (Program.VN[NI][NJ - 1][k] + Program.VN[NI - 1][NJ][k]);
                                Program.V2N[NI][NJ][k] = 0.5F * (Program.VN[NI][NJ - 1][k] + Program.VN[NI - 1][NJ][k]);
                                Program.W1N[NI][NJ][k] = 0.5F * (Program.WN[NI][NJ - 1][k] + Program.WN[NI - 1][NJ][k]);
                                Program.W2N[NI][NJ][k] = 0.5F * (Program.WN[NI][NJ - 1][k] + Program.WN[NI - 1][NJ][k]);

                                Program.U1N[NI][1][k] = 0.5F * (Program.UN[NI][2][k] + Program.UN[NI - 1][1][k]);
                                Program.U2N[NI][1][k] = 0.5F * (Program.UN[NI][2][k] + Program.UN[NI - 1][1][k]);
                                Program.V1N[NI][1][k] = 0.5F * (Program.VN[NI][2][k] + Program.VN[NI - 1][1][k]);
                                Program.V2N[NI][1][k] = 0.5F * (Program.VN[NI][2][k] + Program.VN[NI - 1][1][k]);
                                Program.W1N[NI][1][k] = 0.5F * (Program.WN[NI][2][k] + Program.WN[NI - 1][1][k]);
                                Program.W2N[NI][1][k] = 0.5F * (Program.WN[NI][2][k] + Program.WN[NI - 1][1][k]);
                            }

                            U1_L[k] = U1N_L[k];
                            V1_L[k] = V1N_L[k];
                            W1_L[k] = W1N_L[k];
                            U2_L[k] = U2N_L[k];
                            V2_L[k] = V2N_L[k];
                            W2_L[k] = W2N_L[k];


                            UN_L[k] = 0.5F * (U1N_L[k] + U2N_L[k]);
                            VN_L[k] = 0.5F * (V1N_L[k] + V2N_L[k]);
                            WN_L[k] = 0.5F * (W1N_L[k] + W2N_L[k]);
                            U_L[k] = UN_L[k];
                            V_L[k] = VN_L[k];
                            W_L[k] = WN_L[k];

                            T_L[k] = TN_L[k];
                            QU_L[k] = QUN_L[k];
                            TE_L[k] = TEN_L[k];
                            DISS_L[k] = DISSN_L[k];

                            TABS_L[k] = (T_L[k] + Program.TBZ1) / FAC_L[k]; //updated values for absolute temperature 

                            WAT_VAP_L[k] = WAT_VAPN_L[k];
                        }
                    }
                }
            });

            //Computation of values at receptor points
            if (Program.recexist == true)
            {
                for (int ianz = 0; ianz < Program.Urec.Count(); ianz++)
                {
                    double wind = 0;
                    double z1 = 0;
                    double u1 = 0;
                    double v1 = 0;
                    if (knrec[ianz] == 1)
                    {
                        z1 = 0;
                        u1 = 0;
                        v1 = 0;
                    }
                    else
                    {
                        z1 = Program.ZSPImm[Program.inrec[ianz]][Program.jnrec[ianz]][Program.knrec[ianz] - 1] - Program.AHImm[Program.inrec[ianz]][Program.jnrec[ianz]];
                        u1 = Program.UN[Program.inrec[ianz]][Program.jnrec[ianz]][Program.knrec[ianz] - 1];
                        v1 = Program.VN[Program.inrec[ianz]][Program.jnrec[ianz]][Program.knrec[ianz] - 1];
                    }

                    Linear_interpolation(z1, Program.ZSPImm[Program.inrec[ianz]][Program.jnrec[ianz]][Program.knrec[ianz]] - Program.AHImm[Program.inrec[ianz]][Program.jnrec[ianz]], Program.Zrec[ianz], u1, Program.UN[Program.inrec[ianz]][Program.jnrec[ianz]][Program.knrec[ianz]], ref wind);
                    Program.Urec[ianz] += wind * Program.DT / 3600;
                    Linear_interpolation(z1, Program.ZSPImm[Program.inrec[ianz]][Program.jnrec[ianz]][Program.knrec[ianz]] - Program.AHImm[Program.inrec[ianz]][Program.jnrec[ianz]], Program.Zrec[ianz], v1, Program.VN[Program.inrec[ianz]][Program.jnrec[ianz]][Program.knrec[ianz]], ref wind);
                    Program.Vrec[ianz] += wind * Program.DT / 3600;

                    double T2 = 0;
                    Linear_interpolation(0.0, Program.ZSPImm[Program.inrec[ianz]][Program.jnrec[ianz]][1] - Program.AHImm[Program.inrec[ianz]][Program.jnrec[ianz]], 2.0, Program.TB[Program.inrec[ianz]][Program.jnrec[ianz]][2] - 273.0, Program.TABS[Program.inrec[ianz]][Program.jnrec[ianz]][1] - 273.0, ref T2);
                    Program.Trec[ianz] += T2 * Program.DT / 3600;
                    Program.Globradrec[ianz] += Program.GLOBRAD[Program.inrec[ianz]][Program.jnrec[ianz]] * Program.DT / 3600;
                    Program.Soilheatfluxrec[ianz] += (Program.ALAMBDA[Program.inrec[ianz]][Program.jnrec[ianz]] * (Program.TB[Program.inrec[ianz]][Program.jnrec[ianz]][3] - Program.TB[Program.inrec[ianz]][Program.jnrec[ianz]][2]) / Program.DWB[3] -
                                     Program.ALAMBDA[Program.inrec[ianz]][Program.jnrec[ianz]] * (Program.TB[Program.inrec[ianz]][Program.jnrec[ianz]][2] - Program.TB[Program.inrec[ianz]][Program.jnrec[ianz]][1]) / Program.DWB[2]) * Program.DT / 3600;
                    //influence of snow cover, except in urban areas
                    if (Program.SNOW[Program.inrec[ianz]][Program.jnrec[ianz]] > 0.15 && Program.SNOW[Program.inrec[ianz]][Program.jnrec[ianz]] > Program.Z0[Program.inrec[ianz]][Program.jnrec[ianz]] && Program.ALAMBDA[Program.inrec[ianz]][Program.jnrec[ianz]] != 4.0)
                    {
                        Program.Soilheatfluxrec[ianz] = 0.0;
                    }
                    Program.Longradrec[ianz] += (Program.RL[Program.inrec[ianz]][Program.jnrec[ianz]] - Program.EPSG[Program.inrec[ianz]][Program.jnrec[ianz]] * Program.SIGMA * Pow4(Program.TB[Program.inrec[ianz]][Program.jnrec[ianz]][2])) * Program.DT / 3600;
                    Program.Sensheatfluxrec[ianz] += Program.WQU[Program.inrec[ianz]][Program.jnrec[ianz]] / Program.DDXImm[1] / Program.DDYImm[1] * Program.DT / 3600;
                    double ALW1 = Program.ALW - 2300 * (Program.TB[Program.inrec[ianz]][Program.jnrec[ianz]][2] - 273.15);
                    Program.Latheatfluxrec[ianz] += Program.XWQ[Program.inrec[ianz]][Program.jnrec[ianz]] * Program.RHO[Program.inrec[ianz]][Program.jnrec[ianz]][1] * Program.UST[Program.inrec[ianz]][Program.jnrec[ianz]] * ALW1 * (Program.QU[Program.inrec[ianz]][Program.jnrec[ianz]][1] - Program.QUG[Program.inrec[ianz]][Program.jnrec[ianz]]) * 0.001 * Program.DT / 3600;
                }
            }

            //updated values for absolute temperature - see next loop
            //            Parallel.For(1, NI + 1, Program.pOptions, i =>
            //            {
            //                for (int j = 1; j <= NJ; j++)
            //                {
            //                	double[] TABS_L = Program.TABS[i][j];
            //                	double[] T_L    = Program.T[i][j];
            //                	double[] FAC_L  = Program.FAC[i][j];
            //                	
            //                    for (int k = 1; k <= NK - 1; k++)
            //                    {
            //                    	TABS_L[k] = (T_L[k] + Program.TBZ1) / FAC_L[k];
            //                    }
            //                }
            //            });

            //updated values for soil and absolute temperature
            Parallel.For(1, NI + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= NJ; j++)
                {
                    double[] TBA_L = Program.TBA[i][j];
                    double[] TB_L = Program.TB[i][j];
                    double[] TBN_L = Program.TBN[i][j];

                    for (int k = 1; k <= Program.NZB - 1; k++)
                    {
                        TBA_L[k] = TB_L[k];
                        TB_L[k] = TBN_L[k];
                    }
                    Program.TG[i][j] = Math.Round(TB_L[2], 2);
                }
            });

            //Turbulent eddy viscosity model
            if (Program.ICTE == true)
            {
                //Parallel.For(2, NI, Program.pOptions, i =>
                Parallel.ForEach(Partitioner.Create(2, NI, (int)(NI / Program.pOptions.MaxDegreeOfParallelism)), range =>
                  {
                      for (int i = range.Item1; i < range.Item2; i++)
                      {
                          for (int j = 2; j <= NJ - 1; j++)
                          {
                              double[] VISV_L = Program.VISV[i][j];
                              double[] DISSN_L = Program.DISSN[i][j];
                              double[] TE_L = Program.TE[i][j];
                              double[] TEN_L = Program.TEN[i][j];

                              for (int k = 1; k <= NK - 1; k++)
                              {
                                  if (TE_L[k] <= 0)
                                  {
                                      TE_L[k] = 0;
                                  }

                                  VISV_L[k] += Program.RELAXT * (0.09 * Pow2(TEN_L[k]) / DISSN_L[k] - VISV_L[k]);
                                  VISV_L[k] = Math.Max(Program.VISEL, VISV_L[k]);
                                  VISV_L[k] = Math.Min(50, VISV_L[k]);
                              }
                          }
                      }
                  });
            }
            else
            {
                //computation of turbulent viscosity using an algebraic mixing length model (Flassak, 1990)
                //Parallel.For(2, NI, Program.pOptions, i =>
                Parallel.ForEach(Partitioner.Create(2, NI, (int)(NI / Program.pOptions.MaxDegreeOfParallelism)), range =>
                {
                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        for (int j = 2; j <= NJ - 1; j++)
                        {
                            ReadOnlySpan<float> ZSP_L = Program.ZSPImm[i][j].AsSpan();
                            double[] U_L = Program.U[i][j];
                            double[] V_L = Program.V[i][j];
                            double[] VISV_L = Program.VISV[i][j];
                            double[] RITSCH_L = Program.RITSCH[i][j];

                            for (int k = 1; k <= NK - 1; k++)
                            {
                                Program.ZI[i][j] = 0;
                                double DZZ = (ZSP_L[k + 1] - ZSP_L[k]);
                                double ALUN = Math.Abs(Program.FN) / (0.007 * Program.UST[i][j]);
                                double ALBL = Program.CK * (ZSP_L[k] - Program.AHImm[i][j]) /
                                   (1 + Program.CK * (ZSP_L[k] - ZSP_L[1]) * ALUN);
                                double DVDZ = 0.5 * Math.Sqrt(Pow2((U_L[k + 1] - U_L[k]) / DZZ) + Pow2((V_L[k + 1] - V_L[k]) / DZZ));
                                double FRITSCH = 0;
                                if ((RITSCH_L[k] > 0) && (RITSCH_L[k] <= 0.33))
                                {
                                    FRITSCH = Pow2(1 - 3 * RITSCH_L[k]);
                                }
                                else if ((RITSCH_L[k] <= 0) && (RITSCH_L[k] > -0.048))
                                {
                                    FRITSCH = Math.Pow(1 + 3 * RITSCH_L[k], -2);
                                }
                                else if (RITSCH_L[k] <= -0.048)
                                {
                                    FRITSCH = 2.08 * Math.Pow(-RITSCH_L[k], 0.33);
                                }

                                VISV_L[k] = ALBL * ALBL * DVDZ * FRITSCH;
                                VISV_L[k] = Math.Max(Program.VISEL, VISV_L[k]);

                                if (k == 1)
                                {
                                    VISV_L[k] = 0.4 * Program.UST[i][j] * (ZSP_L[k] - Program.AHImm[i][j]);
                                }
                            }
                        }
                    }
                });
            }

            //homogenous Von Neumann boundary-conditions for TE, VIS, W
            Parallel.For(1, NI + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= NJ; j++)
                {
                    double[] TE_L = Program.TE[i][j];
                    double[] DISS_L = Program.DISS[i][j];
                    double[] VISV_L = Program.VISV[i][j];

                    for (int k = 1; k <= NK; k++)
                    {
                        TE_L[NK] = TE_L[NK - 1];
                        DISS_L[NK] = DISS_L[NK - 1];
                        VISV_L[NK] = VISV_L[NK - 1];

                        if (j == 1) // just one times
                        {
                            Program.TE[i][NJ][k] = Program.TE[i][NJ - 1][k];
                            Program.TE[i][1][k] = Program.TE[i][2][k];

                            Program.DISS[i][NJ][k] = Program.DISS[i][NJ - 1][k];
                            Program.DISS[i][1][k] = Program.DISS[i][2][k];

                            Program.VISV[i][NJ][k] = Program.VISV[i][NJ - 1][k];
                            Program.VISV[i][1][k] = Program.VISV[i][2][k];

                            Program.ZI[i][1] = Program.ZI[i][2];
                            Program.ZI[i][NJ] = Program.ZI[i][NJ - 1];

                            Program.RSOLG[i][1] = Program.RSOLG[i][2];
                            Program.RSOLG[i][NJ] = Program.RSOLG[i][NJ - 1];
                        }

                        if (i == 1)
                        {
                            Program.TE[NI][j][k] = Program.TE[NI - 1][j][k];
                            Program.TE[1][j][k] = Program.TE[2][j][k];

                            Program.DISS[NI][j][k] = Program.DISS[NI - 1][j][k];
                            Program.DISS[1][j][k] = Program.DISS[2][j][k];

                            Program.VISV[NI][j][k] = Program.VISV[NI - 1][j][k];
                            Program.VISV[1][j][k] = Program.VISV[2][j][k];

                            Program.VISV[1][1][k] = 0.5F * (Program.VISV[2][1][k] + Program.VISV[1][2][k]);
                            Program.VISV[NI][1][k] = 0.5F * (Program.VISV[NI - 1][1][k] + Program.VISV[NI][2][k]);
                            Program.VISV[1][NJ][k] = 0.5F * (Program.VISV[2][NJ][k] + Program.VISV[1][NJ - 1][k]);
                            Program.VISV[NI][NJ][k] = 0.5F * (Program.VISV[NI - 1][NJ][k] + Program.VISV[NI][NJ - 1][k]);

                            Program.ZI[1][1] = Program.ZI[2][2];
                            Program.ZI[1][NJ] = Program.ZI[2][NJ - 1];
                            Program.ZI[NI][NJ] = Program.ZI[NI - 1][NJ - 1];
                            Program.ZI[NI][1] = Program.ZI[NI - 1][2];
                            Program.ZI[1][j] = Program.ZI[2][j];
                            Program.ZI[NI][j] = Program.ZI[NI - 1][j];

                            Program.RSOLG[1][1] = Program.RSOLG[2][2];
                            Program.RSOLG[1][NJ] = Program.RSOLG[2][NJ - 1];
                            Program.RSOLG[NI][NJ] = Program.RSOLG[NI - 1][NJ - 1];
                            Program.RSOLG[NI][1] = Program.RSOLG[NI - 1][2];
                            Program.RSOLG[1][j] = Program.RSOLG[2][j];
                            Program.RSOLG[NI][j] = Program.RSOLG[NI - 1][j];
                        }

                    }
                }
            });

            //nudging at the top boundary
            Parallel.For(1, NI + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= NJ; j++)
                {
                    double[] UN_L = Program.UN[i][j];
                    double[] VN_L = Program.VN[i][j];
                    double[] WN_L = Program.WN[i][j];
                    double[] U1N_L = Program.U1N[i][j];
                    double[] U2N_L = Program.U2N[i][j];
                    double[] V1N_L = Program.V1N[i][j];
                    double[] V2N_L = Program.V2N[i][j];
                    double[] W1N_L = Program.W1N[i][j];
                    double[] W2N_L = Program.W2N[i][j];

                    for (int k = NK - 3; k <= NK; k++)
                    {
                        double ATTENU = 1.5;
                        double AALPHA = 1 - Math.Tanh(ATTENU * (NK - k));
                        UN_L[k] -= AALPHA * (UN_L[k] - Program.UG[i][j][k]);
                        VN_L[k] -= AALPHA * (VN_L[k] - Program.VG[i][j][k]);
                        WN_L[k] -= AALPHA * (WN_L[k]);

                        U1N_L[k] = UN_L[k];
                        V1N_L[k] = VN_L[k];
                        U2N_L[k] = UN_L[k];
                        V2N_L[k] = VN_L[k];
                        W1N_L[k] = WN_L[k];
                        W2N_L[k] = WN_L[k];
                    }
                }
            });

            //mass-fluxes at the cell faces
            Parallel.For(1, NI + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= NJ; j++)
                {
                    for (int k = 1; k <= NK; k++)
                    {
                        if (j == 1)
                        {
                            float RHO3 = Program.RHO[i][1][k];
                            float RHO4 = Program.RHO[i][NJ][k];

                            Program.U2NRHO[i][1][k] = (float)(Program.U[i][1][k] * RHO3);
                            Program.U1NRHO[i][NJ][k] = (float)(Program.U[i][NJ][k] * RHO4);

                            Program.V2NRHO[i][1][k] = (float)(Program.V[i][1][k] * RHO3);
                            Program.V1NRHO[i][NJ][k] = (float)(Program.V[i][NJ][k] * RHO4);

                            Program.W2NRHO[i][1][k] = (float)(Program.W[i][1][k] * RHO3);
                            Program.W1NRHO[i][NJ][k] = (float)(Program.W[i][NJ][k] * RHO4);
                        }

                        if (i == 1)
                        {
                            float RHO1 = Program.RHO[1][j][k];
                            float RHO2 = Program.RHO[NI][j][k];

                            Program.U2NRHO[1][j][k] = (float)(Program.U[1][j][k] * RHO1);
                            Program.U1NRHO[NI][j][k] = (float)(Program.U[NI][j][k] * RHO2);

                            Program.V2NRHO[1][j][k] = (float)(Program.V[1][j][k] * RHO1);
                            Program.V1NRHO[NI][j][k] = (float)(Program.V[NI][j][k] * RHO2);

                            Program.W2NRHO[1][j][k] = (float)(Program.W[1][j][k] * RHO1);
                            Program.W1NRHO[NI][j][k] = (float)(Program.W[NI][j][k] * RHO2);
                        }
                    }
                }
            });
        }

        //linear interpolation between two points
        public static void Linear_interpolation(double x1, double x2, double x3, double val1, double val2, ref double val3)
        {
            if (x2 - x1 != 0)
            {
                val3 = val1 + (val2 - val1) / (x2 - x1) * (x3 - x1);
            }
            else
            {
                val3 = val1;
            }
        }

    }
}
