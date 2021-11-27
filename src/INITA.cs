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
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using System.Collections.Immutable;

namespace GRAMM_2001
{
    partial class Program
    {
        /// <summary>
        /// Init basic values
        /// </summary>
        /// <param name="NI"></param>
        /// <param name="NJ"></param>
        /// <param name="NK"></param>
        public static void INITA(int NI, int NJ, int NK)
        {
            double B0 = 38.00699;
            double B1 = 6546.307;
            double B2 = -3.257572;
            double B3 = -0.003163528;

            //Stretching-factor in the vertical
            Program.STRETCH = (Program.ZSPImm[1][1][3] - Program.ZSPImm[1][1][2]) / (Program.ZSPImm[1][1][2] - Program.ZSPImm[1][1][1]);

            //Coriolisfrequencies
            Program.FN = Math.PI / 21600 * Math.Sin(Program.BGRAD * Math.PI / 180);
            Program.FH = Math.PI / 21600 * Math.Cos(Program.BGRAD * Math.PI / 180) * 0;

            //Inital values for the saturation pressure
            Program.PSAT = new double[212]{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,102.9,113.3,
            124.7,136.9,150.4,165.1,180.9,198.1,216.9,237.3,259.4,283.3,309.0,338.0,368.0,
            401.0,437.0,475.0,517.0,562.0,611.0,657.0,706.0,758.0,813.0,872.0,
            935.0,1002,1073,1148,1228,1312,1402,1497,1598,1705,
            1817,1937,2063,2196,2337,2486,2643,2808,2982,3170,
            3360,3560,3780,4000,4240,4490,4750,5030,5320,5620,
            5940,6270,6620,6990,7370,7777,8198,8639,9100,9582,
            10086,10612,11162,11736,12335,12961,13613,14293,15002,15741,
            16511,17313,18147,19016,19920,20860,21840,22860,23910,25010,
            26150,27330,28560,29840,31160,32530,33960,35430,36960,38550,
            40190,41890,43650,45470,47360,49310,51330,53420,55570,57800,
            60110,62490,64950,67490,70110};

            //Saturation temperature
            for (int n = 1; n <= 161; n++)
            {
                double TSAT = 149.95 + (float)n;
                double PSATN = 100000 * Math.Exp(B0 - B1 / TSAT + B2 * Math.Log(TSAT) + B3 * TSAT);
                if (Program.PSAT[n] == 0)
                {
                    Program.PSAT[n] = PSATN;
                }
            }

            //Gravitational acceleration, general gas constant, heat capacity of air by constant pressure, evaporation heat of water
            Program.GERD = 9.81;
            Program.GASCON = 287;
            Program.CPLUFT = 1000;
            Program.ALW = 2500000;

            //Initialization procedure
            //Program.ISOL = 1; //Switch for radiation model: 1=thin clouds 2=thick clouds
            if (Program.ISTAT < 2)
            {
                Temp_INIT(NI, NJ, NK);
            }

            if (Program.ISTAT < 2)
            {
                Parallel.For(1, NI + 1, Program.pOptions, i =>
                {
                    for (int j = 1; j <= NJ; j++)
                    {
                        float[] RHO_L = Program.RHO[i][j];
                        float[] RHOBZ_L = Program.RHOBZ[i][j];
                        double[] W1_L = Program.W1[i][j];
                        double[] W2_L = Program.W2[i][j];
                        double[] W1N_L = Program.W1N[i][j];
                        double[] W2N_L = Program.W2N[i][j];
                        double[] PN_L = Program.PN[i][j];
                        double[] DP_L = Program.DP[i][j];
                        double[] TE_L = Program.TE[i][j];
                        double[] TEN_L = Program.TEN[i][j];
                        double[] DT_SOL_L = Program.DT_SOL[i][j];
                        double[] DT_TERR_L = Program.DT_TERR[i][j];

                        for (int k = 1; k <= NK; k++)
                        {
                            //densitiy
                            RHO_L[k] = (float)(Math.Round(RHOBZ_L[k], 3));

                            //vertical velocities
                            W1_L[k] = 0;
                            W2_L[k] = 0;
                            W1N_L[k] = 0;
                            W2N_L[k] = 0;

                            //non-hydrostatic pressure
                            PN_L[k] = 0;
                            DP_L[k] = 0;

                            //turbulent-kinetic energy
                            TE_L[k] = 0;
                            TEN_L[k] = 0;

                            //solar and terrestrial radiation heating rates
                            DT_SOL_L[k] = 0;
                            DT_TERR_L[k] = 0;
                        }
                    }
                });
                RHOImm = ImmutableArray.Create(Program.RHO);

                //Alpha: parameter used for nudging variables towards large-scale values at lateral boundaries
                for (int k = 1; k <= NK; k++)
                {
                    Program.ALPHA[k] = 0.05;
                }

                //van Karman constant
                Program.CK = 0.35;

                //minimum turbulent viscosity
                Program.VISEL = 0.05;

                //heat capacity soil
                Program.CPBOD = 900;

                //Stefan-Bolzmann constant
                Program.SIGMA = 5.6697e-8;
            }
        }

        // Topographical preprocessor: find U valleys and bassins - June 2018 M. Kuntner
        private static void Analyze_Topography()
        {
            Console.Write("Analyzing topography....  ");

            List<float> TPI_300 = new List<float>();
            List<float> TPI_2000 = new List<float>();

            float Radius_inner = 1.5F * Program.DDXImm[1];
            float Radius_outer_1 = 800;
            float Radius_outer_2 = 4000;

            Read_TPI_Settings(ref Radius_inner, ref Radius_outer_1, ref Radius_outer_2);

            if (Radius_inner < Program.DDXImm[1] || Radius_outer_1 < Program.DDXImm[1] || Radius_outer_2 < 4 * Program.DDXImm[1]) // if settings not valid
            {
                Set_Bassins_to_AHmin();
                Console.WriteLine("failed");
                return;
            }

            int Radius_300 = (int)Math.Max(1, Radius_inner / Program.DDXImm[1] + 0.5F);
            int Radius_2000 = (int)Math.Max(1, Radius_outer_2 / Program.DDXImm[1] + 0.5F);

            float[][] TPI300 = CreateArray<float[]>(NX1, () => new float[NY1]);
            float[][] TPI2000 = CreateArray<float[]>(NX1, () => new float[NY1]);
            float[][] Slope300_Max = CreateArray<float[]>(NX1, () => new float[NY1]);
            float[][] Slope300_Min = CreateArray<float[]>(NX1, () => new float[NY1]);

            int x_start = Math.Max(nr_cell_smooth, Radius_300);
            int x_end = Math.Min(NX - nr_cell_smooth, NX - Radius_300) + 1;
            int y_end = Math.Min(NY - nr_cell_smooth, NY - Radius_300) + 1;

            object obj = new object();
            Parallel.For(x_start, x_end, Program.pOptions, i =>
            {
                for (int j = x_start; j < y_end; j++)
                {
                    float slopemax = 0;
                    float slopemin = 0;
                    TPI300[i][j] = TPI_Calc(0F, Radius_inner, i, j, ref slopemax, ref slopemin);
                    Slope300_Max[i][j] = (float)(Math.Round(slopemax));
                    Slope300_Min[i][j] = (float)(Math.Round(slopemin));
                    lock (obj)
                    {
                        TPI_300.Add(TPI300[i][j]);
                    }
                }
            });

            x_start = Math.Max(nr_cell_smooth, Radius_2000);
            x_end = Math.Min(NX - nr_cell_smooth, NX - Radius_2000) + 1;
            y_end = Math.Min(NY - nr_cell_smooth, NY - Radius_2000) + 1;

            Parallel.For(x_start, x_end, Program.pOptions, i =>
            {
                for (int j = x_start; j < y_end; j++)
                {
                    float slopemax = 0;
                    float slopemin = 0;
                    TPI2000[i][j] = TPI_Calc(Radius_outer_1, Radius_outer_2, i, j, ref slopemax, ref slopemin);
                    lock (obj)
                    {
                        TPI_2000.Add(TPI2000[i][j]);
                    }
                    //Slope300[i][j] = slope;
                }
            });

            obj = null;

            double TPI300_Mean = 0;
            double STDev300 = StdDev(TPI_300, ref TPI300_Mean);

            double TPI2000_Mean = 0;
            double STDev2000 = StdDev(TPI_2000, ref TPI2000_Mean);

            Console.WriteLine("mean near " + Math.Round(TPI300_Mean, 2).ToString() + " mean large " + Math.Round(TPI2000_Mean, 2).ToString());

            float[][] TPI = CreateArray<float[]>(NX1, () => new float[NY1]);

            // Write_ASRII("TPI_SlopeMax.txt", Slope300_Max);

            Parallel.For(1, NX, Program.pOptions, i =>
            {
                for (int j = 1; j < NY; j++)
                {
                    TPI300[i][j] = (int)((TPI300[i][j] - TPI300_Mean) / STDev300 * 100F + 0.5F); // normalize TPI300
                    TPI2000[i][j] = (int)((TPI2000[i][j] - TPI2000_Mean) / STDev2000 * 100F + 0.5F); // normalize TPI2000

                    int TPI_300_Stdi = (int)(TPI300[i][j]);
                    int TPI_2000_Stdi = (int)(TPI2000[i][j]);

                    float slope = Slope300_Max[i][j];

                    // classify TPI values
                    if (TPI_300_Stdi > -100 && TPI_300_Stdi < 100 && TPI_2000_Stdi > -100 && TPI_2000_Stdi < 100 && slope <= 5F)
                    {
                        TPI[i][j] = 5; // Broad flat areas
                    }
                    else if (TPI_300_Stdi > -100 && TPI_300_Stdi < 100 && TPI_2000_Stdi > -100 && TPI_2000_Stdi < 100 && slope > 5F)
                    {
                        TPI[i][j] = 6; // Broad open slopes
                    }
                    else if (TPI_300_Stdi > -100 && TPI_300_Stdi < 100 && TPI_2000_Stdi >= 100)
                    {
                        TPI[i][j] = 7; // flat ridge tops
                    }
                    else if (TPI_300_Stdi > -100 && TPI_300_Stdi < 100 && TPI_2000_Stdi <= -100 && slope <= 5F)
                    {
                        TPI[i][j] = 4; // U Shape valley
                    }
                    else if (TPI_300_Stdi <= -100 && TPI_2000_Stdi > -100 && TPI_2000_Stdi < 100)
                    {
                        TPI[i][j] = 2; // local valley in planes
                    }
                    else if (TPI_300_Stdi >= 100 && TPI_2000_Stdi > -100 && TPI_2000_Stdi < 100)
                    {
                        TPI[i][j] = 9; // local ridges in planes
                    }
                    else if (TPI_300_Stdi <= -100 && TPI_2000_Stdi >= 100)
                    {
                        TPI[i][j] = 3;  // Upland incised drainages
                    }
                    else if (TPI_300_Stdi <= -100 && TPI_2000_Stdi <= -100)
                    {
                        TPI[i][j] = 1; // V Shape valley
                    }
                    else if (TPI_300_Stdi >= 100 && TPI_2000_Stdi >= 100)
                    {
                        TPI[i][j] = 10;  // Mountain Top
                    }
                    else if (TPI_300_Stdi >= 100 && TPI_2000_Stdi <= -100)
                    {
                        TPI[i][j] = 8;  // local hilltop within broad valleys
                    }
                    else
                    {
                        TPI[i][j] = 2; // not found -> set to local valley in plain
                    }

                    if (Program.FW[i][j] > 0.99) // water bodies
                    {
                        TPI[i][j] = 4; // U-valley at large water bodies
                    }
                }
            });

            // Low pass filter to slopemin
            Parallel.For(0, NX, Program.pOptions, i =>
            {
                for (int j = 0; j < NY; j++)
                {
                    TPI2000[i][j] = Slope300_Min[i][j];
                }
            });

            Parallel.For(2, NX - 2, Program.pOptions, i =>
            {
                for (int j = 2; j < NY - 2; j++)
                {
                    double temp = 0.06321 * (0.354 * TPI2000[i - 2][j + 2] +
                    0.354 * TPI2000[i - 2][j - 2] + 0.354 * TPI2000[i + 2][j + 2] +
                    0.354 * TPI2000[i + 2][j - 2] + 0.447 * TPI2000[i - 2][j + 1] +
                    0.447 * TPI2000[i - 2][j - 1] + 0.447 * TPI2000[i + 2][j + 1] +
                    0.447 * TPI2000[i + 2][j - 1] + 0.447 * TPI2000[i - 1][j + 2] +
                    0.447 * TPI2000[i - 1][j - 2] + 0.447 * TPI2000[i + 1][j + 2] +
                    0.447 * TPI2000[i + 1][j - 2] + 0.5 * TPI2000[i][j + 2] +
                    0.5 * TPI2000[i][j - 2] + 0.5 * TPI2000[i + 2][j] +
                    0.5 * TPI2000[i - 2][j] + TPI2000[i - 1][j] +
                    TPI2000[i + 1][j] + TPI2000[i][j + 1] +
                    TPI2000[i][j - 1] +
                    0.707 * TPI2000[i - 1][j + 1] + 0.707 * TPI2000[i + 1][j + 1] + 0.707 * TPI2000[i - 1][j - 1] + 0.707 * TPI2000[i + 1][j - 1] +
                    2 * TPI2000[i][j]);
                    Slope300_Min[i][j] = (float)Math.Round(temp, 1);
                }
            });


            Write_ASRII("TPI_SlopeMin.txt", Slope300_Min);


            // remove little islands of class 1, 4 and 5 
            // 1st copy TPI
            Parallel.For(0, NX, Program.pOptions, i =>
            {
                for (int j = 0; j < NY; j++)
                {
                    TPI2000[i][j] = TPI[i][j];
                }
            });

            // 2nd remove islands
            for (int i = 0; i < NX; i++)
            {
                for (int j = 0; j < NY; j++)
                {
                    if (i < 4 || j < 4 || i > NX - 4 || j > NY - 4) // no classes 1, 4 and 5 at the border
                    {
                        if (TPI[i][j] == 1 || TPI[i][j] == 4 || TPI[i][j] == 5)
                        {
                            TPI[i][j] = 6;
                        }
                    }
                    else
                    {
                        if (TPI[i][j] == 1)
                        {
                            if (Math.Round(Mean_TPI_Value(i, j, TPI2000, 1)) > 2F)
                            {
                                TPI[i][j] = 2;
                            }
                        }
                        if (TPI[i][j] == 4)
                        {
                            if (Math.Round(Mean_TPI_Value(i, j, TPI2000, 4)) != 4)
                            {
                                TPI[i][j] = 6;
                            }
                        }
                        if (TPI[i][j] == 5)
                        {
                            if (Math.Round(Mean_TPI_Value(i, j, TPI2000, 5)) != 5)
                            {
                                TPI[i][j] = 6;
                            }
                        }
                    }
                    Program.TPI[i][j] = TPI[i][j];
                }
            }

            Write_ASRII("TPI_STDI.txt", TPI);

            if (Math.Abs(TPI2000_Mean) < 20 && Math.Abs(TPI300_Mean) < 10) // the standardization to STDEV is not valid at large mean values
            {
                float frac = 0.01F;
                for (int i = 1; i < NX; i++)
                {
                    for (int j = 1; j < NY; j++)
                    {
                        float _TPI = TPI[i][j];
                        if ((_TPI == 1 || _TPI == 4 || _TPI == 5) && Slope300_Min[i][j] < 10) // flat V valleys, U-valleys, broad flat aereas
                        //if ((_TPI == 4 || _TPI == 5) && Slope300_Min[i][j] < 5) // U-valleys, broad flat aereas
                        {
                            AH_Bassins[i][j] = Program.AHImm[i][j] + frac; // take the original height for that cells
                            frac += 0.01F;
                        }
                        else
                        {
                            AH_Bassins[i][j] = -1000F;
                        }
                    }
                }

                //Write_ASRII("TPI_Valleys1.txt", AH_Bassins);

                // now enlarge the cells 
                bool finish = false; // flag, that all cells are modified
                float AH_upper = (float)AHMIN; // upper limit to enlarge the arrays -> increase slowly
                bool[][] enlarged = CreateArray<bool[]>(NX1, () => new bool[NY1]);

                while (finish == false) // as long, as some cells are not modified
                {
                    finish = true;
                    AH_upper += Program.DDXImm[1] / 6F;

                    Parallel.For(0, NX, Program.pOptions, i => // set enlarge[][] to true
                    {
                        for (int j = 0; j < NY; j++)
                        {
                            enlarged[i][j] = true;
                        }
                    });

                    for (int i = 1; i < NX - 1; i++)
                    {
                        for (int j = 1; j < NY - 1; j++)
                        {
                            if (AH_Bassins[i][j] < AH_upper && AH_Bassins[i][j] > -900 && enlarged[i][j] == true) // enlarge this cell if possible for 1 cell to all directions
                            {
                                for (int ii = i - 1; ii < i + 2; ii++)
                                {
                                    for (int jj = j - 1; jj < j + 2; jj++)
                                    {
                                        if (AH_Bassins[ii][jj] < -900F && enlarged[ii][jj] == true
                                        /*&& (Program.AHImm[ii][jj] < (AH_Bassins[i][j] + 250))  // influence up to 250 m
                                        && (Program.AHImm[ii][jj] > (AH_Bassins[i][j] - 50))   // limit the influence of an higher bassin to a lower valley
                                        && (TPI[ii][jj] < 9)*/)                       // do not enlarge at hilltops
                                        {
                                            if (AH_Bassins[ii][jj] < -900F)
                                            {
                                                enlarged[ii][jj] = false; // do not enlarge this cell once more in this loop 
                                            }
                                            AH_Bassins[ii][jj] = AH_Bassins[i][j];
                                        }
                                    }
                                }
                            }
                            if (AH_Bassins[i][j] < -900 /*&& AH_upper < (AHMAX - 200)*/) // still not filled cells available
                            {
                                finish = false;
                            }
                        }
                    }
                }

                // Low pass filter 
                Parallel.For(2, NX - 2, Program.pOptions, i =>
                {
                    for (int j = 2; j < NY - 2; j++)
                    {
                        double temp = 0.06321 * (0.354 * AH_Bassins[i - 2][j + 2] +
                        0.354 * AH_Bassins[i - 2][j - 2] + 0.354 * AH_Bassins[i + 2][j + 2] +
                        0.354 * AH_Bassins[i + 2][j - 2] + 0.447 * AH_Bassins[i - 2][j + 1] +
                        0.447 * AH_Bassins[i - 2][j - 1] + 0.447 * AH_Bassins[i + 2][j + 1] +
                        0.447 * AH_Bassins[i + 2][j - 1] + 0.447 * AH_Bassins[i - 1][j + 2] +
                        0.447 * AH_Bassins[i - 1][j - 2] + 0.447 * AH_Bassins[i + 1][j + 2] +
                        0.447 * AH_Bassins[i + 1][j - 2] + 0.5 * AH_Bassins[i][j + 2] +
                        0.5 * AH_Bassins[i][j - 2] + 0.5 * AH_Bassins[i + 2][j] +
                        0.5 * AH_Bassins[i - 2][j] + AH_Bassins[i - 1][j] +
                        AH_Bassins[i + 1][j] + AH_Bassins[i][j + 1] +
                        AH_Bassins[i][j - 1] +
                        0.707 * AH_Bassins[i - 1][j + 1] + 0.707 * AH_Bassins[i + 1][j + 1] + 0.707 * AH_Bassins[i - 1][j - 1] + 0.707 * AH_Bassins[i + 1][j - 1] +
                        2 * AH_Bassins[i][j]);

                        if (temp < AH_Bassins[i][j])
                        {
                            AH_Bassins[i][j] = (float)(temp);
                        }
                    }
                });


                Parallel.For(1, NX + 1, Program.pOptions, i =>
                {
                    for (int j = 1; j <= NY; j++)
                    {
                        if (AH_Bassins[i][j] < AHMIN)
                        {
                            AH_Bassins[i][j] = (float)AHMIN;
                        }
                        else if (AH_Bassins[i][j] > Program.AHImm[i][j]) // AH_Bassin must not be higher than Program.AHImm[][]
                        {
                            Parallel.For(0, NX, Program.pOptions, ii => // set complete Bassin to lowest cell height inside the bassin
                            {
                                for (int jj = 0; jj < NY; jj++)
                                {
                                    if (Math.Abs(AH_Bassins[ii][jj] - AH_Bassins[i][j]) < 2 * float.Epsilon)
                                    {
                                        AH_Bassins[i][j] = Program.AHImm[i][j];
                                    }
                                }
                            });
                        }
                    }
                });

            }
            else // the topography preprocessing failed
            {
                Set_Bassins_to_AHmin();
                Console.WriteLine("The mean TPI value exceeds the permissible value range - try an appropriate setting by using TPI_Settings.txt");
            }

            //Write_ASRII("TPI_Base.txt", AH_Bassins);

            TPI_300.Clear();
            TPI_2000.Clear();
            TPI300 = null;
            TPI2000 = null;
            Slope300_Max = null;
        }

        static private float TPI_Calc(float R_in, float R_out, int xm, int ym, ref float slopemax, ref float slopemin)
        {
            float Mean_inner = 0; int Count_inner = 0;
            float Mean_outer = 0; int Count_outer = 0;
            float AH_max = -100000; int ind_xmax = 0; int ind_ymax = 0;
            float AH_min = 100000; int ind_xmin = 0; int ind_ymin = 0;

            int r_out_ind = Math.Max(1, (int)(R_out / Program.DDXImm[1] + 1.5F)); // outer donut circle

            for (int xi = xm - r_out_ind; xi <= xm + r_out_ind; xi++)
            {
                for (int yi = ym - r_out_ind; yi <= ym + r_out_ind; yi++)
                {
                    if (xi == xm && yi == ym) // inner point
                    {
                        if (xi > 0 && yi > 0 && xi < NX1 && yi < NY1)
                        {
                            Mean_inner += Program.AHImm[xi][yi];
                            Count_inner++;
                            if (Program.AHImm[xi][yi] > AH_max)
                            {
                                AH_max = Program.AHImm[xi][yi];
                                ind_xmax = xi;
                                ind_ymax = yi;
                            }
                            if (Program.AHImm[xi][yi] < AH_min)
                            {
                                AH_min = Program.AHImm[xi][yi];
                                ind_xmin = xi;
                                ind_ymin = yi;
                            }
                        }
                    }
                    else // outer circle
                    {

                        float r = Distance_between_cells(xm, ym, xi, yi);

                        if (r > R_in && r < R_out) // inside circle?
                        {
                            if (xi > 0 && yi > 0 && xi < NX1 && yi < NY1)
                            {
                                Mean_outer += Program.AHImm[xi][yi];
                                Count_outer++;
                                if (Program.AHImm[xi][yi] > AH_max)
                                {
                                    AH_max = Program.AHImm[xi][yi];
                                    ind_xmax = xi;
                                    ind_ymax = yi;
                                }
                                if (Program.AHImm[xi][yi] < AH_min)
                                {
                                    AH_min = Program.AHImm[xi][yi];
                                    ind_xmin = xi;
                                    ind_ymin = yi;
                                }
                            }
                        }

                    }
                }
            }

            float TPI = 0;
            //Console.Write(Count_outer.ToString() + " ");
            if (Count_outer > 0)
            {
                Mean_outer = Mean_outer / Count_outer;
                TPI = Mean_inner - Mean_outer;
                //Console.Write((AH_max - AH_min).ToString() + " ");
                //slope = (float) (Math.Atan(Math.Abs(AH_max - AH_min) / (R_out * 2)) * 180F / Math.PI); // slope is calculated with a diameter of 1.5 * Program.DDXImm[1]
                // find the maximum slope from cell xm/ym to cells xmin/ymin or xmax/ymax
                float rmin = Distance_between_cells(xm, ym, ind_xmin, ind_ymin);
                float sl1 = 0;
                if (rmin > 1)
                {
                    sl1 = (float)(Math.Atan(Math.Abs(Program.AHImm[xm][ym] - AH_min) / rmin) * 180F / Math.PI);
                }
                float rmax = Distance_between_cells(xm, ym, ind_xmax, ind_ymax);
                float sl2 = 0;
                if (rmax > 1)
                {
                    sl2 = (float)(Math.Atan(Math.Abs(Program.AHImm[xm][ym] - AH_max) / rmax) * 180F / Math.PI);
                }

                slopemax = Math.Max(sl1, sl2);
                slopemin = Math.Min(sl1, sl2);
            }

            return TPI;
        }

        static private float Distance_between_cells(int x1, int y1, int x2, int y2)
        {
            float dx = (x1 - x2) * Program.DDXImm[1];
            float dy = (y1 - y2) * Program.DDXImm[1];
            return (float)(Math.Sqrt(Pow2(dx) + Pow2(dy)));
        }

        static private float Mean_TPI_Value(int x, int y, float[][] TPI, int TPI_Check)
        {
            float TPI_mean = 0;
            int count = 0;
            int dx = (int)(Math.Max(1, 300 / Program.DDXImm[1] + 0.5));

            for (int i = x - dx; i < x + dx; i++)
            {
                for (int j = y - dx; j < y + dx; j++)
                {
                    if (i > 0 && j > 0 && i < NX && j < NY)
                    {
                        if (TPI[i][j] - TPI_Check < 0.5)
                        {
                            TPI_mean += TPI[i][j];
                            count++;
                        }
                        else
                        {
                            TPI_mean += TPI_Check + 2;
                            count++;
                        }
                    }
                }
            }
            float mean = TPI[x][y];
            if (count > 0)
            {
                mean = TPI_mean / count;
            }
            return (float)(mean);
        }
        static private double StdDev(List<float> Height, ref double Avg)
        {
            double sd = 0;
            if (Height.Count > 0)
            {
                double avg = Height.Average();
                //double sumOfSquaresOfDifferences = test.Select(val => (val - avg) * (val - avg)).Sum();
                double sumOfSquaresOfDifferences = Height.Sum(val => (val - avg) * (val - avg));
                sd = Math.Sqrt(sumOfSquaresOfDifferences / Height.Count); // standard-deviation
                Avg = avg;
            }
            return sd;
        }

        //optional: write ASRII txt file
        static private void Write_ASRII(string filename, float[][] values)
        {
            try
            {
                CultureInfo ic = CultureInfo.InvariantCulture;
                using (StreamWriter wt = new StreamWriter(filename))
                {
                    wt.WriteLine("ncols         " + Program.NX.ToString(CultureInfo.InvariantCulture));
                    wt.WriteLine("nrows         " + Program.NY.ToString(CultureInfo.InvariantCulture));
                    wt.WriteLine("xllcorner     " + Program.GRAMM_West.ToString(CultureInfo.InvariantCulture));
                    wt.WriteLine("yllcorner     " + Program.GRAMM_South.ToString(CultureInfo.InvariantCulture));
                    wt.WriteLine("cellsize      " + Program.DDXImm[1].ToString(CultureInfo.InvariantCulture));
                    wt.WriteLine("NODATA_value  " + "-9999");

                    for (int jj = Program.NY; jj >= 1; jj--)
                    {
                        for (int o = 1; o <= Program.NX; o++)
                        {
                            wt.Write((values[o][jj]).ToString(ic) + " ");
                        }
                        wt.WriteLine();
                    }
                }
            }
            catch { }
        }//optional: write ASRII txt file

        static private bool Read_TPI_Settings(ref float r1, ref float r2, ref float r3)
        {
            try
            {
                if (File.Exists("TPI_Settings.txt"))
                {
                    using (FileStream fs = new FileStream("TPI_Settings.txt", FileMode.Open, FileAccess.Read, FileShare.Read))
                    {
                        using (StreamReader myreader = new StreamReader(fs))
                        {
                            string[] text = new string[10];
                            text = myreader.ReadLine().Split(new char[] { ' ', ';', '!', '\t' });
                            text[0] = text[0].Replace(".", decsep);
                            r1 = Convert.ToSingle(text[0]);

                            text = myreader.ReadLine().Split(new char[] { ' ', ';', '!', '\t' });
                            text[0] = text[0].Replace(".", decsep);
                            r2 = Convert.ToSingle(text[0]);

                            text = myreader.ReadLine().Split(new char[] { ' ', ';', '!', '\t' });
                            text[0] = text[0].Replace(".", decsep);
                            r3 = Convert.ToSingle(text[0]);
                        }
                    }
                    return true;
                }
                else
                {
                    return false;
                }
            }
            catch
            {
                return false;
            }
        }

        static private void Set_Bassins_to_AHmin()
        {
            Parallel.For(1, NX, Program.pOptions, i =>
            {
                for (int j = 1; j < NY; j++)
                {
                    AH_Bassins[i][j] = (float)AHMIN;
                }
            });
        }

    }
}
