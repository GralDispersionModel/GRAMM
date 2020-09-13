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
        public static void CALCPR_calculate(int NI, int NJ, int NK)
        {
            //compute mass-fluxes at the cell faces
            Parallel.For(1, NI + 1, Program.pOptions, i =>
            {
                int NK_P = NK; int NJ_P = NJ;

                for (int j = 1; j <= NJ_P; j++)
                {
                    float[] RHO_L = Program.RHO[i][j];
                    float[] U1NRHO_L = Program.U1NRHO[i][j];
                    float[] U2NRHO_L = Program.U2NRHO[i][j];
                    float[] V1NRHO_L = Program.V1NRHO[i][j];
                    float[] V2NRHO_L = Program.V2NRHO[i][j];
                    float[] W1NRHO_L = Program.W1NRHO[i][j];
                    float[] W2NRHO_L = Program.W2NRHO[i][j];
                    double[] U1N_L = Program.U1N[i][j];
                    double[] V1N_L = Program.V1N[i][j];
                    double[] W1N_L = Program.W1N[i][j];
                    double[] U2N_L = Program.U2N[i][j];
                    double[] V2N_L = Program.V2N[i][j];
                    double[] W2N_L = Program.W2N[i][j];
                    float RHO;

                    U1NRHO_L[0] = 0;
                    U2NRHO_L[0] = 0;
                    V1NRHO_L[0] = 0;
                    V2NRHO_L[0] = 0;
                    W1NRHO_L[0] = 0;
                    W2NRHO_L[0] = 0;

                    Program.DPX[i][j][0] = 0;
                    Program.DPY[i][j][0] = 0;
                    Program.DPZ[i][j][0] = 0;
                    Program.DP[i][j][0] = 0;

                    for (int k = 1; k <= NK_P; k++)
                    {
                        RHO = RHO_L[k];
                        U1NRHO_L[k] = (float)(U1N_L[k] * RHO);
                        U2NRHO_L[k] = (float)(U2N_L[k] * RHO);
                        V1NRHO_L[k] = (float)(V1N_L[k] * RHO);
                        V2NRHO_L[k] = (float)(V2N_L[k] * RHO);
                        W1NRHO_L[k] = (float)(W1N_L[k] * RHO);
                        W2NRHO_L[k] = (float)(W2N_L[k] * RHO);

                        /*
                        Program.DPX[i][j][k] += Program.TPX[i][j][k];
                        Program.DPY[i][j][k] += Program.TPY[i][j][k];
                        Program.DPZ[i][j][k] += Program.TPZ[i][j][k];
                        Program.DP[i][j][k] += Program.TP[i][j][k];
                        */
                    }
                }

            });

            //compute mass-divergences   
            Parallel.For(2, NI + 1, Program.pOptions, i =>      
            {
                int NK_P = NK; int NJ_P = NJ;
                //Console.WriteLine(" Partitoner: " + Convert.ToString(range.Item1) +"/" + Convert.ToString(range.Item2));                             	
                for (int j = 2; j <= NJ_P; j++)
                {
                    float[] U1NRHO_L = Program.U1NRHO[i][j];
                    float[] U2NRHO_L = Program.U2NRHO[i][j];
                    float[] V1NRHO_L = Program.V1NRHO[i][j];
                    float[] V2NRHO_L = Program.V2NRHO[i][j];
                    float[] W1NRHO_L = Program.W1NRHO[i][j];
                    float[] W2NRHO_L = Program.W2NRHO[i][j];
                    float[] SUX_L = Program.SUX[i][j];
                    float[] SUY_L = Program.SUY[i][j];
                    float[] SUZ_L = Program.SUZ[i][j];
                    float[] SUXYZ_L = Program.SUXYZ[i][j];
                    float[] AREA_L = Program.AREA[i][j];
                    float[] AREAZ_L = Program.AREAZ[i][j];
                    float[] AREAY_L = Program.AREAY[i][j];
                    float[] AREAX_L = Program.AREAX[i][j];
                    float[] AREAZX_L = Program.AREAZX[i][j];
                    float[] AREAZY_L = Program.AREAZY[i][j];
                    float[] AREAXYZ_L = Program.AREAXYZ[i][j];

                    for (int k = 1; k <= NK_P; k++)
                    {
                        //mass-divergence in east-west direction
                        if ((j < NJ_P) && (k < NK_P)) SUX_L[k] = (Program.U2NRHO[i - 1][j][k] - U1NRHO_L[k]);

                        //mass-divergence in south-north direction
                        if ((i < NI) && (k < NK_P)) SUY_L[k] = (Program.V2NRHO[i][j - 1][k] - V1NRHO_L[k]);

                        //mass-divergence in the z-direction
                        if ((i < NI) && (j < NJ_P)) SUZ_L[k] = ((W2NRHO_L[k - 1] - W1NRHO_L[k]) * AREA_L[k] +
                                                                            (U2NRHO_L[k - 1] - U1NRHO_L[k]) * AREAZX_L[k] +
                                                                            (V2NRHO_L[k - 1] - V1NRHO_L[k]) * AREAZY_L[k]) /
                                                                            AREAZ_L[k];

                        //mass-divergence between the two half-cells
                        if ((i < NI) && (j < NJ_P) && (k < NK_P))
                            SUXYZ_L[k] = ((U1NRHO_L[k] - U2NRHO_L[k]) * AREAX_L[k] +
                                                    (V1NRHO_L[k] - V2NRHO_L[k]) * AREAY_L[k] +
                                                    (W1NRHO_L[k] - W2NRHO_L[k]) * AREA_L[k] +
                                                    (U1NRHO_L[k] - U2NRHO_L[k]) * AREAZX_L[k] +
                                                    (V1NRHO_L[k] - V2NRHO_L[k]) * AREAZY_L[k]) /
                                                    AREAXYZ_L[k];

                        //round-off errors cause the pressure equation to produce meaningless gradients
                        //cutting off the last digits solves this problem largely (Oettl, Sept 2015)
                        if (Math.Abs(SUX_L[k]) < 0.00001) SUX_L[k] = 0;
                        if (Math.Abs(SUY_L[k]) < 0.00001) SUY_L[k] = 0;
                        if (Math.Abs(SUZ_L[k]) < 0.00001) SUZ_L[k] = 0;
                        if (Math.Abs(SUXYZ_L[k]) < 0.00001) SUXYZ_L[k] = 0;
                    }
                }
            });

            //            //compute total mass-divergence
            //            double SUM = 0;
            //            for (int i = 2; i <= NI; i++)
            //            {
            //                for (int j = 2; j <= NJ_P; j++)
            //                {
            //                	double[] SUX_L    = Program.SUX[i][j];
            //                	double[] SUY_L    = Program.SUY[i][j];
            //                	double[] SUZ_L    = Program.SUZ[i][j];
            //                	double[] SUXYZ_L  = Program.SUXYZ[i][j];
            //                    for (int k = 1; k <= NK_P; k++)
            //                    {
            //                        SUM += Math.Abs(SUX_L[k] + SUY_L[k] + SUZ_L[k] + SUXYZ_L[k]);
            //                    }
            //                }
            //            }

            //solve the non-hydrostatic pressure equation iteratively using the TDMA or Thomas-algorithm
            Program.INUMS = 0;
            while (Program.INUMS < 9)
            {
                Program.INUMS++;
                if (Program.ICPN == true) Primp_calculate(NI, NJ, NK);
            }

            //compute pressure gradients to correct wind speeds
            Parallel.For(2, NI, Program.pOptions, i =>
            {
                int NK_P = NK; int NJ_P = NJ;
                for (int j = 2; j <= NJ_P - 1; j++)
                {
                    float[] AREA_L = Program.AREA[i][j];
                    float[] AREAX_L = Program.AREAX[i][j];
                    float[] AREAY_L = Program.AREAY[i][j];
                    float[] AREAZX_L = Program.AREAZX[i][j];
                    float[] AREAZY_L = Program.AREAZY[i][j];

                    double[] DP_L = Program.DP[i][j];
                    double[] DPZ_L = Program.DPZ[i][j];
                    float[] AP0_L = Program.AP0[i][j];
                    double[] U1N_L = Program.U1N[i][j];
                    double[] V1N_L = Program.V1N[i][j];
                    double[] W1N_L = Program.W1N[i][j];
                    double[] U2N_L = Program.U2N[i][j];
                    double[] V2N_L = Program.V2N[i][j];
                    double[] W2N_L = Program.W2N[i][j];
                    double[] DDP1DX_L = Program.DDP1DX[i][j];
                    double[] DDP2DX_L = Program.DDP2DX[i][j];
                    double[] DDP1DY_L = Program.DDP1DY[i][j];
                    double[] DDP2DY_L = Program.DDP2DY[i][j];
                    double[] DDP1DZ_L = Program.DDP1DZ[i][j];
                    double[] DDP2DZ_L = Program.DDP2DZ[i][j];
                    double[] DPX_L = Program.DPX[i][j];
                    double[] DPY_L = Program.DPY[i][j];
                    double DP_LL = 0, DPZ_LL = 0, DPZp_LL = 0, f1 = 0, f2 = 0;
                    int m = 2;

                    for (int kn = 1; kn <= 2 * (NK_P - 1); kn++)
                    {
                        int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)

                        if (kn % 2 != 0) // ODD numbers 1,3,5,..- new k value - set local constants new
                        {
                            DP_LL = DP_L[k];
                            DPZ_LL = DPZ_L[k];
                            DPZp_LL = DPZ_L[k + 1];
                            f1 = (AREAX_L[k] + AREAZX_L[k]) * DP_LL;
                            f2 = (AREAY_L[k] + AREAZY_L[k]) * DP_LL;

                            //Pressure gradients

                            DDP1DX_L[k] = (AREAX_L[k] * DPX_L[k] - f1 + AREAZX_L[k] * DPZ_LL);

                            DDP2DX_L[k] = (-Program.AREAX[i + 1][j][k] * Program.DPX[i + 1][j][k] + f1 - AREAZX_L[k + 1] * DPZp_LL);

                            DDP1DY_L[k] = (AREAY_L[k] * DPY_L[k] - f2 + AREAZY_L[k] * DPZ_LL);

                            DDP2DY_L[k] = (-Program.AREAY[i][j + 1][k] * Program.DPY[i][j + 1][k] + f2 - AREAZY_L[k + 1] * DPZp_LL);

                            DDP1DZ_L[k] = AREA_L[k] * (DPZ_LL - DP_LL);

                            DDP2DZ_L[k] = AREA_L[k + 1] * (DP_LL - DPZp_LL);
                        }

                        //Velocity corrections
                        if ((m - 1) == 1)
                        {
                            m--;
                            float temp = 1 / AP0_L[k];
                            U1N_L[k] += DDP1DX_L[k] * temp;
                            V1N_L[k] += DDP1DY_L[k] * temp;
                            W1N_L[k] += DDP1DZ_L[k] * temp;
                        }
                        else
                        {
                            m = 2;
                            float temp = 1 / AP0_L[k];
                            U2N_L[k] += DDP2DX_L[k] * temp;
                            V2N_L[k] += DDP2DY_L[k] * temp;
                            W2N_L[k] += DDP2DZ_L[k] * temp;

                            //ACHTUNG
                            /*
                            if (i == 1) U1N_L[k] = U2N_L[k];
                            if (i == NI - 1) U2N_L[k] = U1N_L[k];
                            if (j == 1) V1N_L[k] = V2N_L[k];
                            if (j == NJ - 1) V2N_L[k] = V1N_L[k];
                            */

                            Program.UN[i][j][k] = 0.5F * (U1N_L[k] + U2N_L[k]);
                            Program.VN[i][j][k] = 0.5F * (V1N_L[k] + V2N_L[k]);
                            Program.WN[i][j][k] = 0.5F * (W1N_L[k] + W2N_L[k]);
                        }
                    }
                }
            });

            Parallel.For(1, NI + 1, Program.pOptions, i =>
            {
                int NK_P = NK; int NJ_P = NJ;
                for (int j = 1; j <= NJ_P; j++)
                {
                    float[] RHO_L = Program.RHO[i][j];
                    float[] U1NRHO_L = Program.U1NRHO[i][j];
                    float[] U2NRHO_L = Program.U2NRHO[i][j];
                    float[] V1NRHO_L = Program.V1NRHO[i][j];
                    float[] V2NRHO_L = Program.V2NRHO[i][j];
                    float[] W1NRHO_L = Program.W1NRHO[i][j];
                    float[] W2NRHO_L = Program.W2NRHO[i][j];
                    double[] U1N_L = Program.U1N[i][j];
                    double[] V1N_L = Program.V1N[i][j];
                    double[] W1N_L = Program.W1N[i][j];
                    double[] U2N_L = Program.U2N[i][j];
                    double[] V2N_L = Program.V2N[i][j];
                    double[] W2N_L = Program.W2N[i][j];

                    for (int k = 1; k <= NK_P; k++)
                    {
                        float RHO_LL = RHO_L[k];
                        U1NRHO_L[k] = (float)(U1N_L[k] * RHO_LL);
                        U2NRHO_L[k] = (float)(U2N_L[k] * RHO_LL);
                        V1NRHO_L[k] = (float)(V1N_L[k] * RHO_LL);
                        V2NRHO_L[k] = (float)(V2N_L[k] * RHO_LL);
                        W1NRHO_L[k] = (float)(W1N_L[k] * RHO_LL);
                        W2NRHO_L[k] = (float)(W2N_L[k] * RHO_LL);
                    }

                    U1NRHO_L[0] = 0;
                    U2NRHO_L[0] = 0;
                    V2NRHO_L[0] = 0;
                    V2NRHO_L[0] = 0;
                    W1NRHO_L[0] = 0;
                    W2NRHO_L[0] = 0;
                }
            });

            //Final computation of total mass-divergence after velocity correction
            Parallel.For(2, NI + 1, Program.pOptions, i =>
            {
                int NK_P = NK; int NJ_P = NJ;
                //Console.WriteLine(" Partitoner: " + Convert.ToString(range.Item1) +"/" + Convert.ToString(range.Item2));
                for (int j = 2; j <= NJ_P; j++)
                    {
                        float[] SUX_L = Program.SUX[i][j];
                        float[] SUY_L = Program.SUY[i][j];
                        float[] SUZ_L = Program.SUZ[i][j];
                        float[] SUXYZ_L = Program.SUXYZ[i][j];

                        float[] AREA_L = Program.AREA[i][j];
                        float[] AREAZ_L = Program.AREAZ[i][j];
                        float[] AREAX_L = Program.AREAX[i][j];
                        float[] AREAY_L = Program.AREAY[i][j];
                        float[] AREAZY_L = Program.AREAZY[i][j];
                        float[] AREAZX_L = Program.AREAZX[i][j];
                        float[] AREAXYZ_L = Program.AREAXYZ[i][j];

                        float[] U1NRHO_L = Program.U1NRHO[i][j];
                        float[] U2NRHO_L = Program.U2NRHO[i][j];
                        float[] V1NRHO_L = Program.V1NRHO[i][j];
                        float[] V2NRHO_L = Program.V2NRHO[i][j];
                        float[] W1NRHO_L = Program.W1NRHO[i][j];
                        float[] W2NRHO_L = Program.W2NRHO[i][j];


                        for (int k = 1; k <= NK_P - 1; k++)
                        {
                            //mass-divergence in east-west direction
                            if ((j < NJ_P) && (k < NK_P)) SUX_L[k] = (Program.U2NRHO[i - 1][j][k] - U1NRHO_L[k]);

                            //mass-divergence in south-north direction
                            if ((i < NI) && (k < NK_P)) SUY_L[k] = (Program.V2NRHO[i][j - 1][k] - V1NRHO_L[k]);

                            //mass-divergence in the z-direction
                            if ((i < NI) && (j < NJ_P)) SUZ_L[k] = ((W2NRHO_L[k - 1] - W1NRHO_L[k]) * AREA_L[k] +
                                                                            (U2NRHO_L[k - 1] - U1NRHO_L[k]) * AREAZX_L[k] +
                                                                            (V2NRHO_L[k - 1] - V1NRHO_L[k]) * AREAZY_L[k]) /
                                                                            AREAZ_L[k];

                            //mass-divergence between the two half-cells
                            if ((i < NI) && (j < NJ_P) && (k < NK_P))
                                SUXYZ_L[k] = ((U1NRHO_L[k] - U2NRHO_L[k]) * AREAX_L[k] +
                                                    (V1NRHO_L[k] - V2NRHO_L[k]) * AREAY_L[k] +
                                                    (W1NRHO_L[k] - W2NRHO_L[k]) * AREA_L[k] +
                                                    (U1NRHO_L[k] - U2NRHO_L[k]) * AREAZX_L[k] +
                                                    (V1NRHO_L[k] - V2NRHO_L[k]) * AREAZY_L[k]) /
                                                    AREAXYZ_L[k];
                        }
                    }
            });

            //compute total mass-divergence
            Program.SUMG = 0;
            float sum = 0;
            object obj = new object(); // Kuntner 14052018: use parallel.foreach()
            int range_parallel = (int)((NI - 2) / Program.pOptions.MaxDegreeOfParallelism);
            range_parallel = Math.Min(NI - 2, range_parallel); // if NI < range_parallel
            Parallel.ForEach(Partitioner.Create(2, NI, range_parallel), range =>
            {
                float sum_i = 0;
                for (int i = range.Item1; i < range.Item2; ++i)
                //for (int i = 2; i < NI - 1; i++) // Kuntner15022017: change <= NI to < NI-1, same at NJ
                {
                    for (int j = 2; j < NJ - 1; ++j)
                    {
                        float[] SUX_L = Program.SUX[i][j];
                        float[] SUY_L = Program.SUY[i][j];
                        float[] SUZ_L = Program.SUZ[i][j];
                        float[] SUXYZ_L = Program.SUXYZ[i][j];
                        for (int k = 1; k <= NK; ++k)
                        {
                            sum_i += Math.Abs(SUX_L[k] + SUY_L[k] + SUZ_L[k] + SUXYZ_L[k]);
                        }
                    }
                }
                lock (obj)
                {
                    sum += sum_i;
                }
            });
            Program.SUMG = sum; // divergence complete
            obj = null;

            //Temporal trend of the total mass-divergence -> used for the dynamic time-step calculation
            Program.IDIV++;
            if (Program.IDIV == 11) Program.IDIV = 1;
            Program.MASSOURCE[Program.IDIV] = Math.Round(Program.SUMG * 0.1, 0) * 10;

            if (Program.IDIV == 10)
            {
                Int32 ROFFSETANFANG = Convert.ToInt32(Math.Floor(Program.MASSOURCE[1] - 50));
                double QUADRMIN = 10000000;
                for (int ISTEIGUNG = -100; ISTEIGUNG <= 100; ISTEIGUNG++)
                {
                    for (int ROFFSET = ROFFSETANFANG; ROFFSET <= ROFFSETANFANG + 200; ROFFSET++)
                    {
                        double QUADRABW = 0;
                        double RSTEIGUNG = 0;
                        double temp;
                        for (int NIJ = 1; NIJ <= 10; NIJ++)
                        {
                            RSTEIGUNG = (float)ISTEIGUNG * 0.1;
                            temp = RSTEIGUNG * (float)NIJ + ROFFSET - Program.MASSOURCE[NIJ];
                            QUADRABW += temp * temp;
                        }
                        if (QUADRABW <= QUADRMIN)
                        {
                            QUADRMIN = Math.Min(QUADRMIN, QUADRABW);
                            Program.STEIGUNG = RSTEIGUNG;
                        }
                    }
                }
            }

            // 5.4.2017 Ku  
            if (Program.IDIV_Up == 0) // first time
            {
                Program.MASSOURCE_Old = Program.SUMG;
                Program.MASSOURCE_Act = Program.SUMG;
            }
            if (Program.IDIV_Up < 330) Program.IDIV_Up++;
            else
            {
                Program.IDIV_Up = 30;
                Program.IDIV_LockUp = 0;
            }
            double a = Math.Exp(-1D / 20.9);
            Program.MASSOURCE_Act = Program.MASSOURCE_Act * a + Program.SUMG * (1 - a);
            Program.MASSOURCE_Queue.Enqueue(Program.SUMG);
            if (Program.IDIV_Up > 25)
            {
                Program.MASSOURCE_Old = Program.MASSOURCE_Old * a + Program.MASSOURCE_Queue.Dequeue() * (1 - a);
            }
            //Console.WriteLine(MASSOURCE_Act.ToString()+"/"+MASSOURCE_Old.ToString());
            //5.4.2017 Ku 
        }
    }
}
