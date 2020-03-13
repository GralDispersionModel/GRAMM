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

/*
 * This routine computes the initial wind- and temperature fields based on either
 * the file meteopgt.all representing a single point measurement and stability class or
 * on detailled profile and point measurements using a specific input format -> the file name is free in this case
*/
using System;
using System.IO;
using System.Threading.Tasks;

namespace GRAMM_2001
{
    partial class Program
    {
        public static void Temp_INIT(int NI, int NJ, int NK)
        {
            //local variables declaration block
            Random rnd = new Random();
            string VARI;
            string VAIR;
            string ZEILE;
            string FNAME;
            string[] STATION = new string[51];
            string PRONAM;
            string INTIM;
            string INDAT;
            Int32 O = 0;
            //    Int32 P = 0;
            Int16 INPUT;
            double[] DISTX = new double[52];
            double[] DISTY = new double[52];
            double[][] DISTZ = Program.CreateArray<double[]>(52, () => new double[52]);
            double[][] MARK = Program.CreateArray<double[]>(52, () => new double[52]);
            double[][] TEMPI = Program.CreateArray<double[]>(52, () => new double[52]);
            double[][] WINDU = Program.CreateArray<double[]>(52, () => new double[52]);
            double[][] WINDV = Program.CreateArray<double[]>(52, () => new double[52]);
            double[] HUG = new double[52];
            double blh0;
            double blh;
            double USTinit;
            Boolean LOGWIND;
            Boolean meteopgtexist = false;
            double TIMESERIES = 0;
            double SECTORWIDTH = 10;
            double WINDDIR = 0;
            double WINDGE = 0;
            double WU = 0;
            double WV = 0;
            double WIND = 0;
            double WINDI = 0;
            double GRADIENT = 0;
            double moist_adiabatic = 1.0;  //correction factor to determine the correct adiabatic temperature gradient for moist air
            double TMAX = -1000;
            double HMAX = -1;
            double windexpon = 1;
            double HEIGHT = 0;
            Int32 INDI;
            Int32 INDJ;
            Int32 INDK;
            double ZABST;
            Int32 II;
            Int32 JJ;
            Int32 KK;
            double MINH;
            double SUMGEW = 0;
            Int32 IDO = 0;
            Int32 IDU = 0;
            Int32 IDOM = 0;
            Int32 IDON = 0;
            Int32 IDUM = 0;
            Int32 IDUN = 0;
            double UNTO;
            double UNTU;
            double DIFF;
            double ZSPMAX;
            Int32 N2 = 0;
            Int32 M2 = 0;
            double GEW = 0;
            double TFIX = 0;
            double ZSPFIX = 0;
            double UNT = 0;
            double DIFFST;
            double USTR;
            double VSTR;
            double USTR1;
            double VSTR1;
            double USTR2;
            double VSTR2;
            double PUNTEN;
            //double POBEN;
            double PMEER;

            //increment actual computed flow field situation
            Program.IWETTER++;

            //read file GRAMMin.dat for some basic information about the way of initialization (use meteopgt.all or not)
            if (Program.IWETTER == 1) // initialize data
            {
                Program.GRAMMin_File_Read(Program.IWETTER);
                Program.Z0[Program.inrec[0]][Program.jnrec[0]] = (float)(Program.Rauigkeit);
            }
            else // refresh (switching steady state output on/off 
                Program.GRAMMin_File_Read(Program.IWETTER);

            //check whether original number of weather situations in meteopgt.all is exceeded
            if ((Program.IWETTER > Program.meteopgt_nr) && (Program.meteopgt_nr != 0))
            {
                Console.WriteLine();
                Console.WriteLine("GRAMM simulations finished. Press any key to continue...");
                if (Program.IOUTPUT <= 0) 			// if not a SOUNDPLAN Project
                    Console.ReadKey(true); 	// wait for a key input
                Environment.Exit(0); 		// Exit console
            }

            // 11.4.17 Ku use arguments
            if ((Program.IWETTER > IWetter_Console_Last))
            {
                Console.WriteLine();
                Console.Write("GRAMM simulations for weather situations ");
                Console.Write(IWetter_Console_First.ToString() + " to " + IWetter_Console_Last.ToString());
                Console.WriteLine(" finished. Press any key to continue...");
                if (Program.IOUTPUT <= 0) 			// if not a SOUNDPLAN Project
                    Console.ReadKey(true); 	// wait for a key input
                Environment.Exit(0); 		// Exit console
            }

            //Write actually computed flow situation to file DispNrGramm.txt -> used in the GUI 
            Program.Counter++;
            if ((Program.Counter > 0) && (IWetter_Console_First <= 1)) // Write status for one instance only
            {
                Program.Counter = 0;
                using (StreamWriter mywriter = new StreamWriter("DispNrGramm.txt"))
                {
                    mywriter.WriteLine(Convert.ToString(Program.IWETTER));
                }
            }

            //Set all temperatures to zero
            Parallel.For(1, NI + 1, Program.pOptions, i =>
                {
                    for (int j = 1; j <= NJ; j++)
                    {
                        for (int k = 1; k <= NK; k++)
                        {
                            Program.T[i][j][k] = 0;
                            Program.TBZ[i][j][k] = 0;
                            Program.U[i][j][k] = 0;
                            Program.V[i][j][k] = 0;
                        }
                    }
                });

            //Set all temporary fields used for establishing the initial wind field to zero
            Parallel.For(0, 51, Program.pOptions, i =>
            {
                for (int j = 0; j <= 50; j++)
                {
                    TEMPI[i][j] = 0;
                    WINDU[i][j] = 0;
                    WINDV[i][j] = 0;
                    DISTZ[i][j] = 0;
                    DISTX[i] = 0;
                    DISTY[i] = 0;
                }
            });

            //Get number of vertical profiles
            Int32 L = 0;
            Int32 M = 0;

            //GRADIENT defines the neutral moist-adiabatic temperature gradient in case of transient simulations
            GRADIENT = Program.TGRAD;
            moist_adiabatic = -GRADIENT * 105.0; //the factor 105 leads to a slightly stable atmosphere


            if ((Program.METEO == "Y") || (Program.METEO == "y"))
            {
                /*
                 * HEREAFTER A SIMPLE DIAGNOSTIC MODEL IS USED TO OBTAIN THE INITIAL STATES OF
                 * TEMPERATURE AND WIND BASED ON A SINGLE POINT OBSERVATION AND A STABILITY CLASS
                 * AS STORED IN THE FILE METEOPGT.ALL
                */

                //get the specific meteorological conditions
                Console.WriteLine("Met.-Situation Nr.: " + Convert.ToString(Program.IWETTER));
                meteopgtexist = File.Exists("meteopgt.all");
                if (meteopgtexist == true)
                {
                    //in case of dynamic sun rise, intermediate flow field files are appended to meteopgt.all
                    if (Program.meteopgt_nr > 0)
                    {
                        Meteopgtall.meteopgtall_generate(Program.meteopgt_nr, Program.TLIMIT2, Program.IOUTPUT);
                    }

                    Program.IPGT = 1;
                    try
                    {
                        using (FileStream fs = new FileStream("meteopgt.all", FileMode.Open, FileAccess.Read, FileShare.Read))
                        {
                            using (StreamReader myreader = new StreamReader(fs))
                            {
                                string[] text = new string[10];
                                text = myreader.ReadLine().Split(new char[] { ' ', ';', ',', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                                Program.ANEMO = Convert.ToDouble(text[0].Replace(".", Program.decsep));
                                TIMESERIES = Convert.ToDouble(text[1].Replace(".", Program.decsep));
                                SECTORWIDTH = Convert.ToDouble(text[2].Replace(".", Program.decsep));
                                for (int inid = 1; inid <= Program.IWETTER; inid++)
                                {
                                    text = myreader.ReadLine().Split(new char[] { ' ', ';', ',', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                                }
                                text = myreader.ReadLine().Split(new char[] { ' ', ';', ',', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                                WINDDIR = Convert.ToDouble(text[0].Replace(".", Program.decsep));
                                WINDGE = Math.Max(Convert.ToDouble(text[1].Replace(".", Program.decsep)), 0.001);
                                Program.AKLA = Convert.ToInt16(text[2]);
                                Program.Windspeed_meteopgt = WINDGE;

                                // Set DTMAX in dependece of the SC and the initial wind speed 
                                DTMAX = max_time_step_original; //13.4.2017 Ku
                                                                //in case of low wind speed initialization and convective conditions, GRAMM easily becomes unstable
                                if ((Program.AKLA < 3) && (WINDGE < 1))
                                {
                                    Program.DTMAX = Math.Max(1.5, Program.max_time_step_original * 0.25); //4.4.2017 Ku time step original as basis, Math.Max 1,5 s
                                }
                                //in case of higher wind speed initialization, the divergence increases at large time steps
                                if (WINDGE > 0.2) //3.4.2017 Ku
                                {   // max time step = Rastersize/4/Windge
                                    Program.DTMAX = Math.Min(Program.DTMAX, Math.Max(1.5, Program.DDX[1] * 0.25 / WINDGE)); //3.4.2017 Ku
                                }
                                // in case of a retry
                                if (Program.computation_retry > 0)
                                    DTMAX = Math.Max(1.0, DTMAX / (2 * Program.computation_retry)); // try with a reduced max. time step 3.4.2017 Ku - Max(1.0)
                                if (Program.computation_retry > 1) //13.4.2017 Ku try a reduction of relax factors 
                                {
                                    RELAXT *= 0.85;
                                    RELAXV *= 0.85;
                                }
                            }
                        }
                    }
                    catch
                    {
                        Console.WriteLine("Error when reading file meteopgt.all - Execution stopped");
                        Environment.Exit(0);
                    }
                }
                else
                {
                    Console.WriteLine("File meteopgt.all is missing - Execution stopped");
                    Environment.Exit(0);
                }

                //RANDOM wind direction within 10° Sector width to avoid finger like structures with point sources
                if (TIMESERIES == 0)
                {
                    WINDDIR = WINDDIR - SECTORWIDTH / 20 + SECTORWIDTH / 10 * rnd.NextDouble();
                }
                WINDDIR *= 10;
                Console.WriteLine("Wind direction: " + Convert.ToString(Math.Round(WINDDIR, 0)));
                WINDDIR = (270 - WINDDIR) * Math.PI / 180;
                WU = WINDGE * Math.Cos(WINDDIR);
                WV = WINDGE * Math.Sin(WINDDIR);
                Console.WriteLine("Wind speed : " + Convert.ToString(Math.Round(WINDGE, 2)) + "m/s");
                Console.WriteLine("U-component: " + Convert.ToString(Math.Round(WU, 2)) + "m/s");
                Console.WriteLine("V-component: " + Convert.ToString(Math.Round(WV, 2)) + "m/s");
                Console.WriteLine("Stability class: " + Convert.ToString(Program.AKLA));

                //The simulation time is modified according to the stability class and wind speed
                if (Program.TLIMIT2 <= 1)
                {
                    Program.TLIMIT = Math.Sqrt(Math.Pow(Program.DDX[3] * NI, 2) + Math.Pow(Program.DDY[3] * NJ, 2)) / WINDGE * Program.TLIMIT2;
                    Program.DTI = Program.TLIMIT;
                }
                else if ((Program.ISTAT == 0) && (Program.IPGT == 1))
                {
                    //Initial temperature gradient
                    GRADIENT = -0.0065;
                    moist_adiabatic = -GRADIENT * 105.0; //the factor 105 leads to a slightly stable atmosphere
                    moist_adiabatic = 0.9;  //leads to more stable stratification

                    //integration time dependent on stability class and wind speed
                    Meteopgtall.Integrationtime(Program.AKLA, WINDGE, Program.TLIMIT2, Program.meteopgt_nr, ref Program.DTI);
                    Program.TLIMIT = Program.DTI;

                }
                else if ((Program.ISTAT == 0) && (Program.IPGT == 0))
                {
                    if ((WINDGE <= 2) && (Program.AKLA != 2))
                    {
                        Program.DTI = Program.TLIMIT2 * 2;
                        Program.TLIMIT = Program.TLIMIT2 * 2;
                    }
                    else
                    {
                        Program.DTI = Math.Min(Program.TLIMIT2, 400);
                        Program.TLIMIT = Math.Min(Program.TLIMIT2, 400);
                    }
                }
                else if (Program.ISTAT != 0)
                {
                    Program.TLIMIT = Program.TLIMIT2;
                }
                else
                {
                    Program.DTI = Math.Min(Program.TLIMIT2, 12000);
                    Program.TLIMIT = Math.Min(Program.TLIMIT2, 12000);
                }

                //Determination of initial values of surface and soil temperatures
                if (Program.IPGT == 1)
                {
                    if (Program.ICATAFORC != 1)
                    {
                        /*The soil temperature in 1m depth should be around 5-10°C at a latitude around 47N and a sea level of around 300m
                         * Around 2000m above sea level permafrost exists in the Alps (=0°C)
                         * Further, permafrost is found at latitudes larger/smaller than 70N/70S -> dependency on the latitude is considered
                        */
                        Program.TBINIT1 = 293 - 0.005 * Math.Pow(Program.BGRAD, 2) + 0.006 * Math.Abs(Program.BGRAD);
                        if (ISTAT == 0)
                        {
                            Program.TINIT = Program.TBINIT1;
                            Program.TBINIT = Program.TBINIT1;
                        }
                        if (Program.AKLA == 1)
                        {
                            Program.Obini = Math.Min(1 / (-0.37 * Math.Pow(Program.Z0[Program.inrec[0]][Program.jnrec[0]] * 100, -0.55)), -4);
                        }
                        else if (Program.AKLA == 2)
                        {
                            Program.Obini = Math.Min(1 / (-0.12 * Math.Pow(Program.Z0[Program.inrec[0]][Program.jnrec[0]] * 100, -0.50)), -4);
                        }
                        else if (Program.AKLA == 3)
                        {
                            Program.Obini = Math.Min(1 / (-0.067 * Math.Pow(Program.Z0[Program.inrec[0]][Program.jnrec[0]] * 100, -0.56)), -4);
                        }
                        else if (Program.AKLA == 4)
                        {
                            Program.Obini = 9999;
                        }
                        else if (Program.AKLA == 5)
                        {
                            Program.Obini = Math.Max(1 / (0.02 * Math.Pow(Program.Z0[Program.inrec[0]][Program.jnrec[0]] * 100, -0.30)), 4);
                        }
                        else if (Program.AKLA == 6)
                        {
                            Program.Obini = Math.Max(1 / (0.05 * Math.Pow(Program.Z0[Program.inrec[0]][Program.jnrec[0]] * 100, -0.50)), 4);
                        }
                        else if (Program.AKLA == 7)
                        {
                            Program.Obini = Math.Max(1 / (0.2 * Math.Pow(Program.Z0[Program.inrec[0]][Program.jnrec[0]] * 100, -0.55)), 4);
                        }
                    }
                    else
                    {
                        if (ISTAT == 0)
                        {
                            if (Program.AKLA == 1)
                            {
                                GRADIENT = -0.00985;
                                moist_adiabatic = -GRADIENT * 105.0; //the factor 105 leads to a slightly stable atmosphere
                                Program.TINIT = 300;
                                Program.TBINIT1 = 300;
                                Program.TBINIT = 300;
                                Program.Obini = Math.Min(1 / (-0.37 * Math.Pow(Program.Z0[Program.inrec[0]][Program.jnrec[0]] * 100, -0.55)), -4);
                            }
                            else if (Program.AKLA == 2)
                            {
                                GRADIENT = -0.00985;
                                moist_adiabatic = -GRADIENT * 105.0; //the factor 105 leads to a slightly stable atmosphere
                                Program.TINIT = 295;
                                Program.TBINIT1 = 295;
                                Program.TBINIT = 295;
                                Program.Obini = Math.Min(1 / (-0.12 * Math.Pow(Program.Z0[Program.inrec[0]][Program.jnrec[0]] * 100, -0.50)), -4);
                            }
                            else if (Program.AKLA == 3)
                            {
                                GRADIENT = -0.00985;
                                moist_adiabatic = -GRADIENT * 105.0; //the factor 105 leads to a slightly stable atmosphere
                                Program.TINIT = 280;
                                Program.TBINIT1 = 280;
                                Program.TBINIT = 280;
                                Program.Obini = Math.Min(1 / (-0.067 * Math.Pow(Program.Z0[Program.inrec[0]][Program.jnrec[0]] * 100, -0.56)), -4);
                            }
                            else if (Program.AKLA == 4)
                            {
                                GRADIENT = -0.00985;
                                moist_adiabatic = -GRADIENT * 105.0; //the factor 105 leads to a slightly stable atmosphere
                                Program.TINIT = 280;
                                Program.TBINIT1 = 280;
                                Program.TBINIT = 280;
                                Program.Obini = 9999;
                            }
                            else if (Program.AKLA == 5)
                            {
                                GRADIENT = -0.00985;
                                moist_adiabatic = -GRADIENT * 105.0; //the factor 105 leads to a slightly stable atmosphere
                                Program.TINIT = 280;
                                Program.TBINIT1 = 280;
                                Program.TBINIT = 280;
                                Program.Obini = Math.Max(1 / (0.02 * Math.Pow(Program.Z0[Program.inrec[0]][Program.jnrec[0]] * 100, -0.30)), 4);
                            }
                            else if (Program.AKLA == 6)
                            {
                                GRADIENT = -0.00985;
                                moist_adiabatic = -GRADIENT * 105.0; //the factor 105 leads to a slightly stable atmosphere
                                Program.TINIT = 275;
                                Program.TBINIT1 = 275;
                                Program.TBINIT = 275;
                                Program.Obini = Math.Max(1 / (0.05 * Math.Pow(Program.Z0[Program.inrec[0]][Program.jnrec[0]] * 100, -0.50)), 4);
                            }
                            else if (Program.AKLA == 7)
                            {
                                GRADIENT = -0.00985;
                                moist_adiabatic = -GRADIENT * 105.0; //the factor 105 leads to a slightly stable atmosphere
                                Program.TINIT = 270;
                                Program.TBINIT1 = 270;
                                Program.TBINIT = 270;
                                Program.Obini = Math.Max(1 / (0.2 * Math.Pow(Program.Z0[Program.inrec[0]][Program.jnrec[0]] * 100, -0.55)), 4);
                            }
                        }
                    }
                    Console.WriteLine("Initial Obukhov length: " + Convert.ToString(Math.Round(Program.Obini, 3)) + "m");
                    Console.WriteLine("Roughness length: " + Convert.ToString(Math.Round(Program.Z0[Program.inrec[0]][Program.jnrec[0]], 3)) + "m");
                }

                //set initial values for relative humidity in %
                if (((Program.METEO == "Y") || (Program.METEO == "y")) && (Program.ISTAT == 0))
                {
                    if (Program.AKLA == 1) Program.QUINIT = 0.3;
                    if (Program.AKLA == 2) Program.QUINIT = 0.4;
                    if (Program.AKLA == 3) Program.QUINIT = 0.5;
                    if (Program.AKLA == 4) Program.QUINIT = 0.5;
                    if (Program.AKLA == 5) Program.QUINIT = 0.7;
                    if (Program.AKLA == 6) Program.QUINIT = 0.7;
                    if (Program.AKLA == 7) Program.QUINIT = 0.8;
                }

                //vertical temperature profile->not sure if needed here??
                for (int k = 1; k <= NK; k++)
                {
                    Program.T[Program.AHMINI][Program.AHMINJ][k] = Program.TINIT + GRADIENT * (Program.ZSP[Program.AHMINI][Program.AHMINJ][k] - 0);
                    if (Program.ZSP[Program.AHMINI][Program.AHMINJ][k] - Program.AHMIN > Program.ZNEUT)
                    {
                        Program.T[Program.AHMINI][Program.AHMINJ][k] = TMAX;
                    }
                    if (Program.T[Program.AHMINI][Program.AHMINJ][k] > TMAX)
                    {
                        TMAX = Math.Max(TMAX, Program.T[Program.AHMINI][Program.AHMINJ][k]);
                        HMAX = Math.Max(HMAX, Program.ZSP[Program.AHMINI][Program.AHMINJ][k]);
                    }
                }

                //get exponent for power-law wind profile
                if (Program.Obini < 0)
                {
                    windexpon = Math.Max(0.35 - 0.4 * Math.Pow(Math.Abs(Program.Obini), -0.15), 0.05);
                }
                else
                {
                    windexpon = 0.56 * Math.Pow(Program.Obini, -0.15);
                }

                //computation of boundary-layer height (Hanna 1982)
                USTinit = (WINDGE + 0.15) * 0.35 / Math.Log((Program.ZSP[3][3][1] - Program.AH[3][3]) / Program.Rauigkeit);
                if ((Program.Obini >= 0) && (Program.Obini < 100))
                    blh0 = Math.Min(0.4 * Math.Sqrt(USTinit * Program.Obini / Program.FN), 2000);
                else if ((Program.Obini < 0) && (Program.Obini > -100))
                    blh0 = 800;
                else
                    blh0 = Math.Min(0.2 * USTinit / Program.FN, 2000);
                blh = blh0;
                Console.WriteLine("Initial Boundary-Layer height: " + Convert.ToString(Math.Round(blh0, 0)) + "m");
                Console.WriteLine("Fricition velocity: " + Convert.ToString(Math.Round(USTinit, 2)) + "m/s");

                //vertical temperature profile if (Program.AKLA==7)
                Inversion_Height = 400;
                Wind_Velocity = (float)WINDGE;
                {
                    float factor_inv = 1f;
                    if (Wind_Velocity < 0.5F)
                        factor_inv = 1;
                    else if (Wind_Velocity >= 0.5F && Wind_Velocity < 1)
                        factor_inv = 0.6f;
                    else if (Wind_Velocity >= 1 && Wind_Velocity < 1.5)
                        factor_inv = 0.4f;
                    else if (Wind_Velocity >= 1.5)
                        factor_inv = 0.2f;

                    if (Program.Wind_Velocity >= 0.35F)
                    {
                        Inversion_Height *= factor_inv;
                    }
                    else
                    {
                        Inversion_Height = (float)(0.33 * (Program.AHMAX - Program.AHMIN));
                        Inversion_Height = (float)(Math.Max(400, Inversion_Height));
                    }
                    //INV_HEIGHT = (float)(Math.Max(400, INV_HEIGHT));
                    if (Program.AKLA == 7)
                        Console.WriteLine("Inversion height: " + Convert.ToString(Math.Round(Inversion_Height, 0)) + "m");
                }

                Parallel.For(1, NI + 1, Program.pOptions, i =>
                {
                    for (int j = 1; j <= NJ; j++)
                    {
                        for (int k = 1; k <= NK; k++)
                        {
                            Program.T[i][j][k] = Program.TINIT + GRADIENT * Program.ZSP[i][j][k];

                            if (Program.AKLA == 7)
                            {
                                if (Program.Wind_Velocity >= 0.35F)
                                {
                                    // Kuntner 28.6.2018: initialization like SC6, but lower temperature
                                    Program.T[i][j][k] = Program.TINIT + 2 + GRADIENT * Program.ZSP[i][j][k];
                                }
                                else // low wind velocities
                                {
                                    //in case of stability class 7 (strongly stable) ans low wind velocities it is assumed that strong ground inversions have already developed
                                    //these suppress the development of cold air drainage flows in this zone
                                    //the initial air temperature at ground is set 5K below the surface temperature -> probably not correct because the surface temperature is set equal to the air temperature in all cases in INITB.cs
                                    if ((Program.ZSP[i][j][k] - Program.AHMIN) <= Inversion_Height)
                                    {
                                        Program.T[i][j][k] = Program.TINIT - 5 - 0.005 * Program.AHMIN + 0.01 * (Program.ZSP[i][j][k] - Program.AHMIN);
                                    }
                                    else if ((Program.ZSP[i][j][k] - Program.AHMIN) > Inversion_Height)
                                    {
                                        Program.T[i][j][k] = Program.TINIT - 5 - 0.005 * Program.AHMIN + 0.01 * Inversion_Height + GRADIENT * (Math.Max(0, Program.ZSP[i][j][k] - Program.AHMIN - Inversion_Height));
                                    }
                                    else if ((Program.ZSP[i][j][k] - Program.AHMIN) > Program.ZNEUT)
                                    {
                                        Program.T[i][j][k] = Program.T[i][j][k - 1] - 0.003 * (Program.ZSP[i][j][k] - Program.ZSP[i][j][k - 1]);
                                    }
                                }
                            }

                            //in case of stability class 6 (moderatly stable) cold air drainage flows are assumed to be strongest.
                            else if (Program.AKLA == 6)
                            {
                                Program.T[i][j][k] = Program.TINIT + 10 + GRADIENT * Program.ZSP[i][j][k];

                                if ((Program.ZSP[i][j][k] - Program.AHMIN) > Program.ZNEUT)
                                {
                                    Program.T[i][j][k] = Program.T[i][j][k - 1] - 0.003 * (Program.ZSP[i][j][k] - Program.ZSP[i][j][k - 1]);
                                }
                            }
                            //in case of stability class 1 (strongly convective) temperature need to be high enough, otherwise thunderstorms develop
                            else if (Program.AKLA == 1)
                            {
                                Program.T[i][j][k] = Program.TINIT + 15 + GRADIENT * Program.ZSP[i][j][k];

                                if ((Program.ZSP[i][j][k] - Program.AHMIN) > Program.ZNEUT)
                                {
                                    Program.T[i][j][k] = Program.T[i][j][k - 1] - 0.003 * (Program.ZSP[i][j][k] - Program.ZSP[i][j][k - 1]);
                                }
                            }
                            else if (Program.AKLA == 2 || Program.AKLA == 3)
                            {
                                Program.T[i][j][k] = Program.TINIT + 5 + GRADIENT * Program.ZSP[i][j][k];

                                if ((Program.ZSP[i][j][k] - Program.AHMIN) > Program.ZNEUT)
                                {
                                    Program.T[i][j][k] = Program.T[i][j][k - 1] - 0.003 * (Program.ZSP[i][j][k] - Program.ZSP[i][j][k - 1]);
                                }
                            }

                            if (Program.ISTAT == 1)
                            {
                                if ((Program.ZSP[i][j][k] - Program.AHMIN) > Program.ZNEUT)
                                {
                                    Program.T[i][j][k] = Program.TINIT - 0.09 * Math.Max(Program.ZNEUT - Program.AHMIN, 0) + GRADIENT * (Program.ZSP[i][j][k] - Program.ZNEUT - Program.AHMIN);
                                }
                                else
                                {
                                    Program.T[i][j][k] = Program.TINIT - 0.09 * (Program.ZSP[i][j][k] - Program.AHMIN);
                                }
                            }

                            Program.U[i][j][k] = WU * Math.Pow((Program.ZSP[i][j][k] - Program.AH[i][j]) / Program.ANEMO, windexpon);
                            Program.V[i][j][k] = WV * Math.Pow((Program.ZSP[i][j][k] - Program.AH[i][j]) / Program.ANEMO, windexpon);
                        }
                    }
                });
            }
            /*
             * FROM HERE ON THE DIAGNOSTIC WIND FIELD MODEL DEVELOPED BY OETTL, 2000 IS IMPLEMENTED
             * IT GENERATES THE INITIAL STATES FOR TEMPERATURE AND WIND BASED ON PROFILE AND POINT OBSERVATIONS
             * NOTE THAT IT REQUIRES A SPECIFIC FILE FORMAT, WHICH IS GENERATED BY THE MODEL ITSELF AS ALL
             * OBSERVATIONS HAVE TO BE PROVIDED VIA THE KEYBOARD. THE FILE IS NAMED BY THE USER AND THEN STORED FOR 
             * FURTHER USAGE
             * TO INVOKE THIS MODEL THE FIRST LINE IN THE FILE GRAMMin.dat NEEDS TO BE "n" INSTEAD OF "y"
            */
            else
            {
                Console.WriteLine();
                Console.WriteLine("********************************************");
                Console.WriteLine("********************************************");
                Console.WriteLine("**                                        **");
                Console.WriteLine("**                                        **");
                Console.WriteLine("**       START OF TEMPERATURE INPUT       **");
                Console.WriteLine("**                                        **");
                Console.WriteLine("**                                        **");
                Console.WriteLine("********************************************");
                Console.WriteLine("********************************************");
                Console.WriteLine(" ");
                Console.WriteLine(" ");
                Console.WriteLine("         AT LEAST ONE VERTICAL PROFILE      ");
                Console.WriteLine("         IS NECESSARY FOR INTERPOLATION     ");
                Console.WriteLine(" ");
                Console.WriteLine(" ");
                Console.WriteLine("Should a stable stratification be provided");
                Console.Write("starting from a certain height (Y)es (N)o:  ");
                VARI = Console.ReadLine();
                Console.WriteLine(" ");

                double GRADIENT1 = 0;
                HEIGHT = 0;
                if ((VARI == "Y") || (VARI == "y"))
                {
                    Console.WriteLine(" ");
                    Console.Write("starting at which height:  ");
                    HEIGHT = Convert.ToDouble(Console.ReadLine().Replace(".", Program.decsep));
                    Console.WriteLine(" ");
                    Console.Write("Gradient [K/m]:  ");
                    GRADIENT1 = Convert.ToDouble(Console.ReadLine().Replace(".", Program.decsep));
                    Console.WriteLine(" ");
                }

                Console.WriteLine(" ");
                Console.Write("READ FROM EXISTING FILE (Y)ES (N)O : ");
                VARI = Console.ReadLine();
                Console.WriteLine(" ");

                //HEREAFTER OBSERVATIONS CAN BE PROVIDED BY THE USER VIA THE KEYBOARD
                if ((VARI != "Y") && (VARI != "y"))
                {
                    Console.WriteLine(" ");
                    Console.Write("   PROJECT NAME                :  ");
                    PRONAM = Console.ReadLine();
                    Console.WriteLine(" ");
                    Console.Write("   INPUT DATE [DAY.MONTH.YEAR] :  ");
                    INDAT = Console.ReadLine();
                    Console.WriteLine(" ");
                    Console.Write("   INPUT TIME [HOUR:MINUTE]    :  ");
                    INTIM = Console.ReadLine();
                GOTO98:
                    Console.WriteLine(" ");
                    Console.WriteLine(" ");
                    Console.WriteLine("   INPUT OF POINT MEASUREMENTS = 1");
                    Console.WriteLine("   INPUT OF VERTICAL PROFILES  = 2");
                    Console.Write("   SELECT NUMBER....");
                    INPUT = Convert.ToInt16(Console.ReadLine());
                    Console.WriteLine(" ");

                    //INPUT A POINT OBSERVATION
                    if (INPUT == 1)
                    {
                        Console.WriteLine("********************************************");
                        Console.WriteLine("********************************************");
                        Console.WriteLine("**                                        **");
                        Console.WriteLine("**                                        **");
                        Console.WriteLine("**     INPUT OF ONE POINT MEASUREMENT     **");
                        Console.WriteLine("**                                        **");
                        Console.WriteLine("**                                        **");
                        Console.WriteLine("********************************************");
                        Console.WriteLine("********************************************");
                        Console.WriteLine(" ");
                        Console.WriteLine(" ");
                        L++;
                        Console.Write("       MONITORING STATION               : ");
                        STATION[L] = Console.ReadLine();
                        Console.WriteLine(" ");
                    WESTERNBOUNDARY:
                        Console.Write("       DISTANCE FROM WESTERN BOUNDARY  [M]: ");
                        DISTX[L] = Convert.ToDouble(Console.ReadLine().Replace(".", Program.decsep));
                        Console.WriteLine(" ");
                        if (DISTX[L] > (Program.X[NI] + Program.DDX[NI])) goto WESTERNBOUNDARY;
                        SOUTHERNBOUNDARY:
                        Console.Write("       DISTANCE FROM SOUTHERN BOUNDARY  [M]: ");
                        DISTY[L] = Convert.ToDouble(Console.ReadLine().Replace(".", Program.decsep));
                        Console.WriteLine(" ");
                        if (DISTY[L] > (Program.Y[NJ] + Program.DDY[NJ])) goto SOUTHERNBOUNDARY;
                        SEALEVEL:
                        Console.Write("       HEIGHT ABOVE SEA LEVEL        [M]: ");
                        DISTZ[L][0] = Convert.ToDouble(Console.ReadLine().Replace(".", Program.decsep));
                        Console.WriteLine(" ");
                        if ((DISTZ[L][0] > Program.Z[NK + 1]) || (DISTZ[L][0] < Program.Z[1])) goto SEALEVEL;
                        Console.Write("       HEIGHT ABOVE GROUND LEVEL     [M]: ");
                        HUG[L] = Convert.ToDouble(Console.ReadLine().Replace(".", Program.decsep));
                        Console.WriteLine(" ");
                        Console.Write("       TEMPERATURE             [CELSIUS]: ");
                        TEMPI[L][0] = Convert.ToDouble(Console.ReadLine().Replace(".", Program.decsep));
                        Console.WriteLine(" ");
                        Console.Write("       WIND SPEED                  [M/S]: ");
                        WIND = Convert.ToDouble(Console.ReadLine().Replace(".", Program.decsep));
                        Console.WriteLine(" ");
                        Console.Write("       WIND DIRECTION (N=0,S=180)  [DEG]: ");
                        WINDI = Convert.ToDouble(Console.ReadLine().Replace(".", Program.decsep));
                        Console.WriteLine(" ");
                        Console.Write("       ADD TO VERTICAL PROFILES ? (Y)ES (N)O : ");
                        VAIR = Console.ReadLine();
                        Console.WriteLine(" ");
                        TEMPI[L][0] += 273.15;
                        WINDU[L][0] = -WIND * Math.Sin(WINDI * Math.PI / 180);
                        WINDV[L][0] = -WIND * Math.Cos(WINDI * Math.PI / 180);

                        //CALCULATION OF CELL INDICES
                        INDI = 0;
                        INDJ = 0;
                        INDK = 0;
                        ZABST = 10000;
                        for (int i = 1; i <= NI; i++)
                        {
                            for (int j = 1; j <= NJ; j++)
                            {
                                for (int k = 1; k <= NK; k++)
                                {
                                    if ((DISTX[L] == (Program.X[i] + Program.DDX[i] * 0.5F)) &&
                                       (DISTY[L] == (Program.Y[j] + Program.DDY[j] * 0.5F)) &&
                                        DISTZ[L][0] == Program.ZSP[i][j][k])
                                    {
                                        INDI = i;
                                        INDJ = j;
                                        INDK = k;
                                    }
                                    if ((DISTX[L] >= Program.X[i]) && (DISTX[L] < Program.X[i + 1])) II = i;
                                    if ((DISTY[L] >= Program.Y[j]) && (DISTY[L] < Program.Y[j + 1]))
                                    {
                                        JJ = j;
                                        MINH = Math.Abs(Program.ZSP[i][j][k] - DISTZ[L][0]);
                                        if (MINH <= ZABST)
                                        {
                                            ZABST = MINH;
                                            KK = k;
                                        }
                                    }
                                }
                            }
                        }
                        //INPUT TO BE CONTROLLED BY THE USER
                        Console.WriteLine(" ");
                        Console.WriteLine("   DISTANCE X            = " + Convert.ToString(DISTX[L]));
                        Console.WriteLine("   DISTANCE Y            = " + Convert.ToString(DISTY[L]));
                        Console.WriteLine("   ABSOLUTE HEIGHT       = " + Convert.ToString(DISTZ[L][0]));
                        Console.WriteLine("   RELATIVE HEIGHT       = " + Convert.ToString(HUG[L]));
                        Console.WriteLine("   TEMPERATURE IN KELVIN = " + Convert.ToString(TEMPI[L][0]));
                        Console.WriteLine("   WIND IN WEST/EAST     = " + Convert.ToString(WINDU[L][0]));
                        Console.WriteLine("   WIND IN SOUTH/NORTH   = " + Convert.ToString(WINDV[L][0]));
                        Console.WriteLine(" ");
                        Console.Write("   INPUT CORRECT ??     (C)ONTINUE ?    (A)BORT ? ");
                        VARI = Console.ReadLine();
                        if ((VARI == "A") || (VARI == "a"))
                        {
                            L--;
                            goto GOTO98;
                        }
                        if ((INDI != 0) && (INDJ != 0) && (INDK != 0))
                        {
                            Program.T[INDI][INDJ][INDK] = TEMPI[L][0];
                            Program.U[INDI][INDJ][INDK] = WINDU[L][0];
                            Program.V[INDI][INDJ][INDK] = WINDV[L][0];
                        }
                        if ((VAIR == "Y") || (VAIR == "y"))
                        {
                            DISTZ[L][1] = DISTZ[L][0];
                            TEMPI[L][1] = TEMPI[L][0];
                            TEMPI[L][0] = 0;
                        }
                        Console.WriteLine(" ");
                        Console.Write("   MORE INPUT ??     (Y)ES ?    (N)O? ");
                        VARI = Console.ReadLine();
                        if ((VARI == "Y") || (VARI == "y"))
                            goto GOTO98;
                    }
                    //INPUT A VERTICAL PROFILE
                    else if (INPUT == 2)
                    {
                        Console.WriteLine(" ");
                        Console.WriteLine(" ");
                        Console.WriteLine("********************************************");
                        Console.WriteLine("********************************************");
                        Console.WriteLine("**                                        **");
                        Console.WriteLine("**                                        **");
                        Console.WriteLine("**     INPUT OF ONE VERTICAL  PROFILE     **");
                        Console.WriteLine("**                                        **");
                        Console.WriteLine("**                                        **");
                        Console.WriteLine("********************************************");
                        Console.WriteLine("********************************************");
                        Console.WriteLine(" ");
                        Console.WriteLine(" ");
                        L++;
                        M = 0;
                        Console.Write("       MONITORING STATION               : ");
                        STATION[L] = Console.ReadLine();
                        Console.WriteLine(" ");
                    WESTERNBOUNDARY_V:
                        Console.Write("       DISTANCE FROM WESTERN BOUNDARY  [M]: ");
                        DISTX[L] = Convert.ToDouble(Console.ReadLine().Replace(".", Program.decsep));
                        Console.WriteLine(" ");
                        if (DISTX[L] > (Program.X[NI] + Program.DDX[NI])) goto WESTERNBOUNDARY_V;
                        SOUTHERNBOUNDARY_V:
                        Console.Write("       DISTANCE FROM SOUTHERN BOUNDARY  [M]: ");
                        DISTY[L] = Convert.ToDouble(Console.ReadLine().Replace(".", Program.decsep));
                        Console.WriteLine(" ");
                        if (DISTY[L] > (Program.Y[NJ] + Program.DDY[NJ])) goto SOUTHERNBOUNDARY_V;
                        Console.WriteLine("      HEIGHT ABOVE GROUND LEVEL");
                        Console.Write("       FOR CLOSEST POINT TO GROUND   [M]: ");
                        HUG[L] = Convert.ToDouble(Console.ReadLine().Replace(".", Program.decsep));
                    GOTO96:
                        Console.WriteLine(" ");
                        M++;
                    SEALEVEL_V:
                        Console.Write("       HEIGHT ABOVE SEA LEVEL        [M]: ");
                        DISTZ[L][M] = Convert.ToDouble(Console.ReadLine().Replace(".", Program.decsep));
                        Console.WriteLine(" ");
                        if ((DISTZ[L][M] > Program.Z[NK + 1])) goto SEALEVEL_V;
                        Console.Write("       TEMPERATURE             [CELSIUS]: ");
                        TEMPI[L][M] = Convert.ToDouble(Console.ReadLine().Replace(".", Program.decsep));
                        Console.WriteLine(" ");
                        Console.Write("       WIND SPEED                  [M/S]: ");
                        WIND = Convert.ToDouble(Console.ReadLine().Replace(".", Program.decsep));
                        Console.WriteLine(" ");
                        Console.Write("       WIND DIRECTION (N=0,S=180)  [DEG]: ");
                        WINDI = Convert.ToDouble(Console.ReadLine().Replace(".", Program.decsep));
                        Console.WriteLine(" ");
                        TEMPI[L][M] += 273.15;
                        WINDU[L][M] = -WIND * Math.Sin(WINDI * Math.PI / 180);
                        WINDV[L][M] = -WIND * Math.Cos(WINDI * Math.PI / 180);

                        //CALCULATION OF CELL INDICES
                        INDI = 0;
                        INDJ = 0;
                        INDK = 0;
                        ZABST = 10000;
                        for (int i = 1; i <= NI; i++)
                        {
                            for (int j = 1; j <= NJ; j++)
                            {
                                for (int k = 1; k <= NK; k++)
                                {
                                    if ((DISTX[L] == (Program.X[i] + Program.DDX[i] * 0.5F)) &&
                                       (DISTY[L] == (Program.Y[j] + Program.DDY[j] * 0.5F)) &&
                                        DISTZ[L][M] == Program.ZSP[i][j][k])
                                    {
                                        INDI = i;
                                        INDJ = j;
                                        INDK = k;
                                    }
                                    if ((DISTX[L] >= Program.X[i]) && (DISTX[L] < Program.X[i + 1])) II = i;
                                    if ((DISTY[L] >= Program.Y[j]) && (DISTY[L] < Program.Y[j + 1]))
                                    {
                                        JJ = j;
                                        MINH = Math.Abs(Program.ZSP[i][j][k] - DISTZ[L][M]);
                                        if (MINH <= ZABST)
                                        {
                                            ZABST = MINH;
                                            KK = k;
                                        }
                                    }
                                }
                            }
                        }

                        //INPUT TO BE CONTROLLED BY THE USER
                        Console.WriteLine(" ");
                        Console.WriteLine("   DISTANCE X            = " + Convert.ToString(DISTX[L]));
                        Console.WriteLine("   DISTANCE Y            = " + Convert.ToString(DISTY[L]));
                        Console.WriteLine("   ABSOLUTE HEIGHT       = " + Convert.ToString(DISTZ[L][M]));
                        Console.WriteLine("   RELATIVE HEIGHT       = " + Convert.ToString(HUG[L]));
                        Console.WriteLine("   TEMPERATURE IN KELVIN = " + Convert.ToString(TEMPI[L][M]));
                        Console.WriteLine("   WIND IN WEST/EAST     = " + Convert.ToString(WINDU[L][M]));
                        Console.WriteLine("   WIND IN SOUTH/NORTH   = " + Convert.ToString(WINDV[L][M]));
                        Console.WriteLine(" ");
                        Console.Write("   INPUT CORRECT ??  (E)ND ?   (C)ONTINUE ?  (A)BORT ? ");
                        VARI = Console.ReadLine();
                        if ((VARI == "A") || (VARI == "a"))
                        {
                            for (int j = 1; j <= M; j++)
                            {
                                TEMPI[L][j] = 0;
                                WINDU[L][j] = 0;
                                WINDV[L][j] = 0;
                            }
                            L--;
                            goto GOTO98;
                        }
                        else if ((VARI == "C") || (VARI == "c"))
                        {
                            if ((INDI != 0) && (INDJ != 0) && (INDK != 0))
                            {
                                Program.T[INDI][INDJ][INDK] = TEMPI[L][M];
                                Program.U[INDI][INDJ][INDK] = WINDU[L][M];
                                Program.V[INDI][INDJ][INDK] = WINDV[L][M];
                            }
                            goto GOTO96;
                        }
                        Console.WriteLine(" ");
                        Console.Write("   MORE INPUT ??     (Y)ES ?    (N)O? ");
                        VARI = Console.ReadLine();
                        if ((VARI == "Y") || (VARI == "y"))
                            goto GOTO98;
                    }
                    else
                        goto GOTO98;

                    //SAVE INPUT FILE
                    Console.WriteLine(" ");
                    Console.Write("   SAVE INPUT TO FILE ??   (Y)ES ?  (N)O ? ");
                    VARI = Console.ReadLine();
                    if ((VARI == "Y") || (VARI == "y"))
                    {
                        Console.WriteLine(" ");
                        Console.Write("   FILENAME ? ");
                        FNAME = Console.ReadLine();
                        using (StreamWriter mywriter = new StreamWriter(FNAME))
                        {
                            mywriter.WriteLine(PRONAM);
                            mywriter.WriteLine(INDAT);
                            mywriter.WriteLine(INTIM);
                            mywriter.WriteLine(" NUMBER OF MONITORING STATIONS: " + Convert.ToString(L).PadLeft(2, '0'));
                            for (int i = 1; i <= L; i++)
                            {
                                mywriter.WriteLine("*");
                                if (STATION[i] != " ")
                                {
                                    if (HUG[i] != 0)
                                    {
                                        mywriter.WriteLine(" MONITORING STATION  : " + Convert.ToString(STATION[i].PadRight(20)) + "(" + i.ToString("00") + ")");
                                        mywriter.WriteLine(" DISTANCE FROM W     : " + Convert.ToString(DISTX[i].ToString("0.00")).PadLeft(8).Replace(Program.decsep, ".") + "(" + i.ToString("00") + ")");
                                        mywriter.WriteLine(" DISTANCE FROM S     : " + Convert.ToString(DISTY[i].ToString("0.00")).PadLeft(8).Replace(Program.decsep, ".") + "(" + i.ToString("00") + ")");
                                        mywriter.WriteLine(" HEIGHT ABOVE GROUND : " + Convert.ToString(HUG[i].ToString("0.00")).PadLeft(8).Replace(Program.decsep, ".") + "(" + i.ToString("00") + ")");
                                    }
                                    for (int j = 0; j <= 51; j++)
                                    {
                                        if (DISTZ[i][j] != 0)
                                        {
                                            mywriter.WriteLine(" HEIGHT ABOVE SEA LEV: " + Convert.ToString(DISTZ[i][j].ToString("0.00")).PadLeft(8).Replace(Program.decsep, ".") + "(" + i.ToString("00") + "," + j.ToString("00") + ")");
                                            mywriter.WriteLine(" TEMPERATURE      [C]: " + Convert.ToString(TEMPI[i][j].ToString("0.00")).PadLeft(8).Replace(Program.decsep, ".") + "(" + i.ToString("00") + "," + j.ToString("00") + ")");
                                            mywriter.WriteLine(" U-COMPONENT OF WIND : " + Convert.ToString(WINDU[i][j].ToString("0.00")).PadLeft(8).Replace(Program.decsep, ".") + "(" + i.ToString("00") + "," + j.ToString("00") + ")");
                                            mywriter.WriteLine(" V-COMPONENT OF WIND : " + Convert.ToString(WINDV[i][j].ToString("0.00")).PadLeft(8).Replace(Program.decsep, ".") + "(" + i.ToString("00") + "," + j.ToString("00") + ")");
                                        }
                                    }
                                }
                            }
                            mywriter.WriteLine("**");

                            for (int i = 1; i <= NI; i++)
                            {
                                for (int j = 1; j <= NJ; j++)
                                {
                                    for (int k = 1; k <= NK; k++)
                                    {
                                        if ((Program.T[i][j][k] != 0) || (Program.U[i][j][k] != 0) || (Program.V[i][j][k] != 0))
                                        {
                                            mywriter.WriteLine("T(I,J,K)   : " + Convert.ToString(Program.T[i][j][k].ToString("0.00")).PadLeft(8).Replace(Program.decsep, ".") + "," + i.ToString() + "," + j.ToString() + "," + k.ToString("00"));
                                            mywriter.WriteLine("U-COMPONENT: " + Convert.ToString(Program.U[i][j][k].ToString("0.00")).PadLeft(8).Replace(Program.decsep, ".") + "," + i.ToString() + "," + j.ToString() + "," + k.ToString("00"));
                                            mywriter.WriteLine("V-COMPONENT: " + Convert.ToString(Program.V[i][j][k].ToString("0.00")).PadLeft(8).Replace(Program.decsep, ".") + "," + i.ToString() + "," + j.ToString() + "," + k.ToString("00"));
                                        }
                                    }
                                }
                            }
                        } // using (mywriter) 

                        Console.WriteLine(" ");
                        Console.Write("   WROTE INPUT TO FILE " + FNAME);
                    }
                }
                //HEREAFTER AN ALREADY EXISTING FILE IS OPENED
                else
                {
                    Console.WriteLine(" ");
                    Console.WriteLine("************************************************");
                    Console.WriteLine("**                                            **");
                    Console.WriteLine("**         READ FROM EXISTING FILE            **");
                    Console.WriteLine("**                                            **");
                    Console.WriteLine("************************************************");
                    Console.WriteLine(" ");
                    Console.Write("   FILENAME ?? ");
                    FNAME = Console.ReadLine();
                    StreamReader myreader = new StreamReader(FNAME);
                    try
                    {
                        PRONAM = myreader.ReadLine();
                        INDAT = myreader.ReadLine();
                        INTIM = myreader.ReadLine();
                        string[] text = new string[10];
                        text = myreader.ReadLine().Split(new char[] { ',', ':', '(', ')' });
                        L = Convert.ToInt32(text[1]);
                        ZEILE = myreader.ReadLine();
                        for (int i = 1; i <= L; i++)
                        {
                            text = myreader.ReadLine().Split(new char[] { ',', ':', '(', ')' });
                            STATION[i] = text[1];
                            text = myreader.ReadLine().Split(new char[] { ',', ':', '(', ')' });
                            DISTX[i] = Convert.ToDouble(text[1].Replace(".", Program.decsep));
                            text = myreader.ReadLine().Split(new char[] { ',', ':', '(', ')' });
                            DISTY[i] = Convert.ToDouble(text[1].Replace(".", Program.decsep));
                            text = myreader.ReadLine().Split(new char[] { ',', ':', '(', ')' });
                            HUG[i] = Convert.ToDouble(text[1].Replace(".", Program.decsep));
                            for (int j = 0; j <= 50; j++)
                            {
                                text = myreader.ReadLine().Split(new char[] { ',', ':', '(', ')' });
                                if (text[0] == "*")
                                    break;
                                if (text[0] == "**")
                                    goto FINISH;
                                II = Convert.ToInt32(text[2]);
                                JJ = Convert.ToInt32(text[3]);
                                DISTZ[II][JJ] = Convert.ToDouble(text[1].Replace(".", Program.decsep));

                                //CALCULATION OF CELL INDICES
                                INDI = 0;
                                INDJ = 0;
                                INDK = 0;
                                ZABST = 10000;
                                for (int iii = 1; iii <= NI; iii++)
                                {
                                    for (int jjj = 1; jjj <= NJ; jjj++)
                                    {
                                        for (int kkk = 1; kkk <= NK; kkk++)
                                        {
                                            if ((DISTX[i] >= Program.X[iii]) && (DISTX[i] < Program.X[iii + 1]))
                                            {
                                                int IIII = iii;
                                                if ((DISTY[i] >= Program.Y[jjj]) && (DISTY[i] < Program.X[jjj + 1]))
                                                {
                                                    int JJJJ = jjj;
                                                    MINH = Math.Abs(Program.ZSP[iii][jjj][kkk] - DISTZ[II][JJ]);
                                                    if (MINH <= ZABST)
                                                    {
                                                        ZABST = MINH;
                                                        int KKKK = kkk;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                text = myreader.ReadLine().Split(new char[] { ',', ':', '(', ')' });
                                TEMPI[II][JJ] = Convert.ToDouble(text[1].Replace(".", Program.decsep));
                                text = myreader.ReadLine().Split(new char[] { ',', ':', '(', ')' });
                                WINDU[II][JJ] = Convert.ToDouble(text[1].Replace(".", Program.decsep));
                                text = myreader.ReadLine().Split(new char[] { ',', ':', '(', ')' });
                                WINDV[II][JJ] = Convert.ToDouble(text[1].Replace(".", Program.decsep));
                            }
                        }
                    FINISH:
                        try
                        {
                            for (int i = 1; i <= 100000; i++)
                            {
                                text = myreader.ReadLine().Split(new char[] { ',', ':', '(', ')' });
                                II = Convert.ToInt32(text[2]);
                                int J = Convert.ToInt32(text[3]);
                                int K = Convert.ToInt32(text[4]);
                                Program.T[II][J][K] = Convert.ToDouble(text[1].Replace(".", Program.decsep));
                                if ((II != 1) || (J != 1) || (K != 1)) Program.T[II][J][K] = 0;
                                text = myreader.ReadLine().Split(new char[] { ',', ':', '(', ')' });
                                II = Convert.ToInt32(text[2]);
                                J = Convert.ToInt32(text[3]);
                                K = Convert.ToInt32(text[4]);
                                Program.U[II][J][K] = Convert.ToDouble(text[1].Replace(".", Program.decsep));
                                if ((II != 1) || (J != 1) || (K != 1)) Program.U[II][J][K] = 0;
                                text = myreader.ReadLine().Split(new char[] { ',', ':', '(', ')' });
                                II = Convert.ToInt32(text[2]);
                                J = Convert.ToInt32(text[3]);
                                K = Convert.ToInt32(text[4]);
                                Program.V[II][J][K] = Convert.ToDouble(text[1].Replace(".", Program.decsep));
                                if ((II != 1) || (J != 1) || (K != 1)) Program.V[II][J][K] = 0;
                            }
                        }
                        catch
                        { }
                        myreader.Close();
                        myreader.Dispose();
                    }
                    catch
                    {
                        Console.WriteLine("Unable to open file: " + FNAME);
                        Environment.Exit(0);
                    }

                    Console.WriteLine(" ");
                    Console.Write("   ADD MORE MONITORING STATIONS ? (Y)ES (N)O : -> THIS FUNCTIONALITY IS NOT SUPPORTED WITH THE CURRENT VERSION  ");
                    VARI = Console.ReadLine();
                    if ((VARI == "Y") || (VARI == "y"))
                    //here should be a jump to the goto GOTO98 statement, which is impossible in C#
                    { }
                }
                //INTERPOLATION PROCEDURE

                //Check, whether observation is identical with a cell-centre
                Parallel.For(1, NI + 1, Program.pOptions, i =>
                {
                    for (int j = 1; j <= NJ; j++)
                    {
                        for (int k = 1; k <= NK; k++)
                        {
                            if (Program.T[i][j][k] != 0) Console.WriteLine("T(I,J,K)   : " + Convert.ToString(Program.T[i][j][k].ToString("0.00")).PadLeft(8).Replace(Program.decsep, ".") + "," + i.ToString() + "," + j.ToString() + "," + k.ToString());
                            if (Program.U[i][j][k] != 0) Console.WriteLine("U(I,J,K)   : " + Convert.ToString(Program.U[i][j][k].ToString("0.00")).PadLeft(8).Replace(Program.decsep, ".") + "," + i.ToString() + "," + j.ToString() + "," + k.ToString());
                            if (Program.V[i][j][k] != 0) Console.WriteLine("V(I,J,K)   : " + Convert.ToString(Program.V[i][j][k].ToString("0.00")).PadLeft(8).Replace(Program.decsep, ".") + "," + i.ToString() + "," + j.ToString() + "," + k.ToString());
                        }
                    }
                });

                //Select terrain-following or height dependent wind interpolation
                LOGWIND = false;
                Console.WriteLine(" ");
                Console.WriteLine(" ");
                Console.Write(" WIND INTERPOLATION TO SAME HEIGHT (Y)ES (N)O:  ");
                VARI = Console.ReadLine();
                if ((VARI == "Y") || (VARI == "y"))
                {
                    LOGWIND = true;
                }
                Console.WriteLine(" ");
                Console.WriteLine(" ");
                Console.WriteLine("*****  INTERPOLATION  OF TEMPERATURE  *****");
                for (int i = 1; i <= NI; i++)
                {
                    for (int j = 1; j <= NJ; j++)
                    {
                        for (int k = 1; k <= NK; k++)
                        {
                            if (Program.T[i][j][k] == 0)
                            {
                                //Mark all cells next above and below every cell to be interpolated
                                SUMGEW = 0;
                                IDO = 0;
                                IDU = 0;
                                IDOM = 0;
                                IDON = 0;
                                IDUM = 0;
                                IDUN = 0;
                                for (int n = 1; n <= L + 1; n++)
                                {
                                    MARK[IDON][IDOM] = 1;
                                    MARK[IDUN][IDUM] = 1;
                                    IDOM = 0;
                                    IDON = 0;
                                    IDUM = 0;
                                    IDUN = 0;
                                    UNTO = -100000;
                                    UNTU = 100000;
                                    for (int m = 1; m <= 50; m++)
                                    {
                                        MARK[n][m] = 0;
                                        if (TEMPI[n][m] == 0)
                                            break;
                                        DIFF = Program.ZSP[i][j][k] - DISTZ[n][m];
                                        if ((DIFF <= 0) && (DIFF > UNTO))
                                        {
                                            UNTO = DIFF;
                                            IDON = n;
                                            IDOM = m;
                                            IDO = 1;
                                        }
                                        if ((DIFF > 0) && (DIFF < UNTU))
                                        {
                                            UNTU = DIFF;
                                            IDUN = n;
                                            IDUM = m;
                                            IDU = 1;
                                        }
                                    }
                                }
                                //IF NO VALUE ABOVE EXISTS
                                if (IDO == 0)
                                {
                                    ZSPMAX = -1000;
                                    for (int n = 1; n <= L; n++)
                                    {
                                        for (int m = 1; m <= 50; m++)
                                        {
                                            if ((MARK[n][m] == 1) && (m > 1) && (TEMPI[n][m - 1] != 0) && (DISTZ[n][m] > ZSPMAX))
                                            {
                                                N2 = n;
                                                M2 = m - 1;
                                                MARK[n][m] = 0;
                                                ZSPMAX = DISTZ[n][m];
                                            }
                                        }
                                        MARK[N2][M2] = 1;
                                        MARK[N2][M2 + 1] = 1;
                                    }
                                }
                                //IF NO VALUE BELOW EXISTS
                                if (IDU == 0)
                                {
                                    for (int n = 1; n <= L; n++)
                                    {
                                        IDU = 0;
                                        for (int m = 1; m <= 50; m++)
                                        {
                                            if ((IDU == 0) && (MARK[n][m] == 1) && (TEMPI[n][m + 1] != 0))
                                            {
                                                MARK[n][m + 1] = 1;
                                                IDU = 1;
                                            }
                                        }
                                    }
                                }
                                //Computation of temperature                            
                                for (int n = 1; n <= L; n++)
                                {
                                    for (int m = 1; m <= 50; m++)
                                    {
                                        if (MARK[n][m] == 1)
                                        {
                                            double TEMP1 = TEMPI[n][m];
                                            double GEW1 = Math.Pow(DISTX[n] - Program.X[i] - Program.DDX[i] * 0.5, 2) + Math.Pow(DISTY[n] - Program.Y[j] - Program.DDY[j] * 0.5, 2);
                                            Int32 MAUT = m;
                                            for (int o = n; o <= L; o++)
                                            {
                                                for (int p = (MAUT + 1); p <= 50; p++)
                                                {
                                                    if (MARK[o][p] == 1)
                                                    {
                                                        if ((DISTZ[n][m] - DISTZ[o][p] == 0))
                                                            break;
                                                        double TEMP2 = TEMPI[o][p];
                                                        double GEW2 = Math.Pow(DISTX[o] - Program.X[i] - Program.DDX[i] * 0.5, 2) + Math.Pow(DISTY[o] - Program.Y[j] - Program.DDY[j] * 0.5, 2);
                                                        if ((GEW1 == 0) && (GEW2 == 0))
                                                            GEW = 0.00000000001;
                                                        else if ((n == 0) && (GEW1 != 0) && (GEW2 != 0))
                                                            GEW = Math.Sqrt(GEW1 + GEW2);
                                                        else
                                                            GEW = GEW1 + GEW2;
                                                        SUMGEW += 1 / GEW;
                                                        Program.T[i][j][k] += ((Program.ZSP[i][j][k] - DISTZ[n][m]) * (TEMP1 - TEMP2) / (DISTZ[n][m] - DISTZ[o][p]) + TEMP1) / GEW;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                //Weighting factors
                                if (SUMGEW != 0) Program.T[i][j][k] /= SUMGEW;
                                SUMGEW = 1;

                                //Influence of ground measurements
                                for (int n = 1; n <= L; n++)
                                {
                                    if ((TEMPI[n][0] != 0) && (Math.Abs(DISTZ[n][0] - Program.ZSP[i][j][k]) < 20))
                                    {
                                        double DUMMY = Math.Sqrt(Math.Pow(DISTX[n] - Program.X[i] - Program.DDX[i] * 0.5, 2) + Math.Pow(DISTY[n] - Program.Y[j] - Program.DDY[j] * 0.5, 2));
                                        if (DUMMY == 0)
                                            GEW = 1 / 10000;
                                        else
                                            GEW = DUMMY / 500;
                                        SUMGEW += 1 / GEW;
                                        Program.T[i][j][k] += TEMPI[n][0] / GEW;
                                    }
                                }
                                //Weighting factors
                                if (SUMGEW != 0) Program.T[i][j][k] /= SUMGEW;
                            }
                        }
                    }
                }
                for (int k = 1; k <= NK; k++)
                {
                    if ((Program.ZSP[Program.AHMINI][Program.AHMINJ][k] >= HEIGHT) && (HEIGHT != 0))
                    {
                        TFIX = Program.T[Program.AHMINI][Program.AHMINJ][k - 1];
                        ZSPFIX = Program.ZSP[Program.AHMINI][Program.AHMINJ][k - 1];
                        break;
                    }
                }
                Parallel.For(1, NI + 1, Program.pOptions, i =>
                {
                    for (int j = 1; j <= NJ; j++)
                    {
                        for (int k = 1; k <= NK; k++)
                        {
                            if ((Program.ZSP[i][j][k] >= HEIGHT) && (HEIGHT != 0))
                                Program.T[i][j][k] = TFIX + GRADIENT1 * (Program.ZSP[i][j][k] - ZSPFIX);
                        }
                    }
                });
                //check if each cell is assigned a temperature
                M = 0;
                O = 0;
                for (int i = 1; i <= NI; i++)
                {
                    for (int j = 1; j <= NJ; j++)
                    {
                        for (int k = 1; k <= NK; k++)
                        {
                            if (Program.T[i][j][k] == 0) M++;
                            O++;
                        }
                    }
                }
                Console.WriteLine("NUMBER OF UNDEFINED TEMPERATURES           : " + M.ToString());
                Console.WriteLine("TOTAL NUMBER OF CELLS: " + O.ToString());
                if (M != 0)
                {
                    Console.WriteLine("INTERPOLATION FAILED!!");
                    Environment.Exit(0);
                }

                //Interpolation of windcomponents
                Console.WriteLine();
                Console.WriteLine("***** INTERPOLATION OF WINDCOMPONENTS *****");
                Console.WriteLine();
                //interpolate wind fields terrain-following
                if (LOGWIND == false)
                {
                    for (int i = 1; i <= NI; i++)
                    {
                        for (int j = 1; j <= NJ; j++)
                        {
                            for (int k = 1; k <= NK; k++)
                            {
                                Program.W[i][j][k] = 0;
                                SUMGEW = 0;
                                if ((Program.U[i][j][k] != 0) && (Program.V[i][j][k] != 0))
                                    continue;
                                UNT = Program.ZSP[i][j][k] - Program.AH[i][j];
                                for (int n = 1; n <= L; n++)
                                {
                                    IDOM = 0;
                                    IDUM = 0;
                                    IDON = 0;
                                    IDUN = 0;
                                    IDU = 0;
                                    IDO = 0;
                                    UNTO = -100000;
                                    UNTU = 100000;
                                    DIFFST = 0;
                                    for (int m = 0; m <= 50; m++)
                                    {
                                        if ((WINDU[n][m] != 0) || (WINDV[n][m] != 0))
                                        {
                                            DIFFST = Math.Max(DISTZ[n][0], DISTZ[n][1]) - HUG[n];
                                            DIFF = (Program.ZSP[i][j][k] - Program.AH[i][j]) - (DISTZ[n][m] - DIFFST);
                                            if ((DIFF <= 0) && (DIFF > UNTO))
                                            {
                                                IDON = n;
                                                IDOM = m;
                                                IDO = 1;
                                                UNTO = DIFF;
                                            }
                                            if ((DIFF > 0) && (DIFF < UNTU))
                                            {
                                                IDUN = n;
                                                IDUM = m;
                                                IDU = 1;
                                                UNTU = DIFF;
                                            }
                                        }
                                    }
                                    //interpolation within the Prandtl layer
                                    double GEW1 = 0;
                                    if ((IDO == 0) && (IDU == 1) && (UNT <= 70))
                                    {
                                        USTR = WINDU[IDUN][IDUM] * Math.Pow((Program.ZSP[i][j][k] - Program.AH[i][j]) / (DISTZ[IDUN][IDUM] - DIFFST), 0.25);
                                        GEW1 = Math.Pow((DISTX[IDUN] - Program.X[i] - Program.DDX[i] * 0.5F), 2) + Math.Pow((DISTY[IDUN] - Program.Y[j] - Program.DDY[j] * 0.5F), 2);
                                        if (GEW1 == 0) GEW1 = 0.000000000001;
                                        Program.U[i][j][k] += USTR / GEW1;
                                        SUMGEW += 1 / GEW1;
                                        VSTR = WINDV[IDUN][IDUM] * Math.Pow((Program.ZSP[i][j][k] - Program.AH[i][j]) / (DISTZ[IDUN][IDUM] - DIFFST), 0.25);
                                        Program.V[i][j][k] += VSTR / GEW1;
                                    }
                                    else if ((IDO == 1) && (IDU == 0))
                                    {
                                        USTR = WINDU[IDON][IDOM] * Math.Pow((Program.ZSP[i][j][k] - Program.AH[i][j]) / (DISTZ[IDON][IDOM] - DIFFST), 0.25);
                                        GEW1 = Math.Pow((DISTX[IDON] - Program.X[i] - Program.DDX[i] * 0.5F), 2) + Math.Pow((DISTY[IDON] - Program.Y[j] - Program.DDY[j] * 0.5F), 2);
                                        if (GEW1 == 0) GEW1 = 0.000000000001;
                                        Program.U[i][j][k] += USTR / GEW1;
                                        SUMGEW += 1 / GEW1;
                                        VSTR = WINDV[IDON][IDOM] * Math.Pow((Program.ZSP[i][j][k] - Program.AH[i][j]) / (DISTZ[IDON][IDOM] - DIFFST), 0.25);
                                        Program.V[i][j][k] += VSTR / GEW1;
                                    }
                                    else if ((IDO == 1) && (IDU == 1))
                                    {
                                        USTR = WINDU[IDUN][IDUM] + (WINDU[IDON][IDOM] - WINDU[IDUN][IDUM]) /
                                            (DISTZ[IDON][IDOM] - DISTZ[IDUN][IDUM]) * (Program.ZSP[i][j][k] - Program.AH[i][j] - DISTZ[IDUN][IDUM] + DIFFST);
                                        GEW1 = Math.Pow((DISTX[IDON] - Program.X[i] - Program.DDX[i] * 0.5F), 2) + Math.Pow((DISTY[IDON] - Program.Y[j] - Program.DDY[j] * 0.5F), 2);
                                        if (GEW1 == 0) GEW1 = 0.000000000001;
                                        Program.U[i][j][k] += USTR / GEW1;
                                        SUMGEW += 1 / GEW1;
                                        VSTR = WINDV[IDUN][IDUM] + (WINDV[IDON][IDOM] - WINDV[IDUN][IDUM]) /
                                            (DISTZ[IDON][IDOM] - DISTZ[IDUN][IDUM]) * (Program.ZSP[i][j][k] - Program.AH[i][j] - DISTZ[IDUN][IDUM] + DIFFST);
                                        Program.V[i][j][k] += VSTR / GEW1;
                                    }
                                    //interpolation above Prandtl-layer
                                    else if ((IDO == 0) && (IDU == 1) && (UNT > 70) && (IDUM != 0))
                                    {
                                        GEW1 = Math.Pow((DISTX[IDON] - Program.X[i] - Program.DDX[i] * 0.5F), 2) + Math.Pow((DISTY[IDON] - Program.Y[j] - Program.DDY[j] * 0.5F), 2);
                                        if (GEW1 == 0) GEW1 = 0.000000000001;
                                        Program.U[i][j][k] += WINDU[IDUN][IDUM] / GEW1;
                                        Program.V[i][j][k] += WINDV[IDUN][IDUM] / GEW1;
                                        SUMGEW += 1 / GEW1;
                                    }
                                }
                                //weighting factors
                                if (SUMGEW != 0) Program.U[i][j][k] /= SUMGEW;
                                if (SUMGEW != 0) Program.V[i][j][k] /= SUMGEW;
                            }
                        }
                    }
                }
                //interpolate wind fields height dependent
                else
                {
                    for (int i = 1; i <= NI; i++)
                    {
                        for (int j = 1; j <= NJ; j++)
                        {
                            for (int k = 1; k <= NK; k++)
                            {
                                if ((Program.U[i][j][k] == 0) && (Program.V[i][j][k] == 0))
                                {
                                    //mark all cells closest above and below every cell to be interpolated
                                    IDOM = 0;
                                    IDUM = 0;
                                    IDON = 0;
                                    IDUN = 0;
                                    IDU = 0;
                                    IDO = 0;
                                    SUMGEW = 0;
                                    for (int n = 1; n <= L + 1; n++)
                                    {
                                        MARK[IDON][IDOM] = 1;
                                        MARK[IDUN][IDUM] = 1;
                                        IDOM = 0;
                                        IDUM = 0;
                                        IDON = 0;
                                        IDUN = 0;
                                        UNTO = -100000;
                                        UNTU = 100000;
                                        for (int m = 0; m <= 50; m++)
                                        {
                                            MARK[n][m] = 0;
                                            if ((WINDU[n][m] != 0) || (WINDV[n][m] != 0))
                                            {
                                                DIFF = Program.ZSP[i][j][k] - DISTZ[n][m];
                                                if ((DIFF <= 0) && (DIFF > UNTO))
                                                {
                                                    IDON = n;
                                                    IDOM = m;
                                                    IDO = 1;
                                                    UNTO = DIFF;
                                                }
                                                if ((DIFF > 0) && (DIFF < UNTU))
                                                {
                                                    IDUN = n;
                                                    IDUM = m;
                                                    IDU = 1;
                                                    UNTU = DIFF;
                                                }
                                            }
                                        }
                                    }
                                    //computation of windcomponents in between two measurements
                                    double GEW1 = 0;
                                    double GEW2 = 0;
                                    double DUMMY = 0;
                                    Int32 MAUT = 0;
                                    for (int n = 1; n <= L; n++)
                                    {
                                        for (int m = 0; m <= 50; m++)
                                        {
                                            if (MARK[n][m] == 1)
                                            {
                                                USTR1 = WINDU[n][m];
                                                VSTR1 = WINDV[n][m];
                                                GEW1 = Math.Pow((DISTX[n] - Program.X[i] - Program.DDX[i] * 0.5F), 2) + Math.Pow((DISTY[n] - Program.Y[j] - Program.DDY[j] * 0.5F), 2);
                                                DUMMY = 0;
                                                MAUT = m;
                                                for (int o = n; o <= L; o++)
                                                {
                                                    if (DUMMY == 1) MAUT = 0;
                                                    DUMMY = 1;
                                                    for (int p = MAUT + 1; p <= 50; p++)
                                                    {
                                                        if (MARK[o][p] == 1)
                                                        {
                                                            if ((DISTZ[n][m] - DISTZ[o][p]) == 0) break;
                                                            USTR2 = WINDU[o][p];
                                                            VSTR2 = WINDV[o][p];
                                                            GEW2 = Math.Pow((DISTX[o] - Program.X[i] - Program.DDX[i] * 0.5F), 2) + Math.Pow((DISTY[o] - Program.Y[j] - Program.DDY[j] * 0.5F), 2);
                                                            if ((GEW1 == 0) && (GEW2 == 0)) GEW = 0.0000000001;
                                                            else if ((n == 0) && (GEW1 != 0) && (GEW2 != 0)) GEW = Math.Sqrt(GEW1 + GEW2);
                                                            else
                                                                GEW = GEW1 + GEW2;
                                                            SUMGEW += 1 / GEW;
                                                            Program.U[i][j][k] += ((Program.ZSP[i][j][k] - DISTZ[n][m]) * (USTR1 - USTR2) / (DISTZ[n][m] - DISTZ[o][p]) + USTR1) / GEW;
                                                            Program.V[i][j][k] += ((Program.ZSP[i][j][k] - DISTZ[n][m]) * (VSTR1 - VSTR2) / (DISTZ[n][m] - DISTZ[o][p]) + USTR1) / GEW;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    //weighting factors
                                    if (SUMGEW != 0) Program.U[i][j][k] /= SUMGEW;
                                    if (SUMGEW != 0) Program.V[i][j][k] /= SUMGEW;

                                    //if no value above exists - use no gradient for interpolation
                                    if (IDO == 0)
                                    {
                                        double DISTZMAX = -100000;
                                        for (int n = 1; n <= L; n++)
                                        {
                                            //only 1 station is used
                                            for (int m = 0; m <= 50; m++)
                                            {
                                                if ((IDO == 0) && (MARK[n][m] == 1) && (DISTZ[n][m] > DISTZMAX))
                                                {
                                                    DISTZMAX = DISTZ[n][m];
                                                    Program.U[i][j][k] = WINDU[n][m];
                                                    Program.V[i][j][k] = WINDV[n][m];
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                //Smoothing the velocity fields with a 5-point filter
                Console.WriteLine();
                Console.WriteLine("***** SMOOTHING OF VELOCITYS *****");
                Console.WriteLine();
                for (int i = 2; i <= NI - 1; i++)
                {
                    for (int j = 2; j <= NJ - 1; j++)
                    {
                        for (int k = 1; k <= NK; k++)
                        {
                            //velocity smoothing in corners
                            Program.U[1][1][k] = 0.333 * (Program.U[1][2][k] + Program.U[2][1][k] + Program.U[1][1][k]);
                            Program.V[1][1][k] = 0.333 * (Program.V[1][2][k] + Program.V[2][1][k] + Program.V[1][1][k]);
                            Program.U[1][NJ][k] = 0.333 * (Program.U[1][NJ - 1][k] + Program.U[2][NJ][k] + Program.U[1][NJ][k]);
                            Program.V[1][NJ][k] = 0.333 * (Program.V[1][NJ - 1][k] + Program.V[2][NJ][k] + Program.V[1][NJ][k]);
                            Program.U[NI][NJ][k] = 0.333 * (Program.U[NI - 1][NJ][k] + Program.U[NI][NJ - 1][k] + Program.U[NI][NJ][k]);
                            Program.V[NI][NJ][k] = 0.333 * (Program.V[NI - 1][NJ][k] + Program.V[NI][NJ - 1][k] + Program.V[NI][NJ][k]);
                            Program.U[NI][1][k] = 0.333 * (Program.U[NI - 1][1][k] + Program.U[NI][2][k] + Program.U[NI][1][k]);
                            Program.V[NI][1][k] = 0.333 * (Program.V[NI - 1][1][k] + Program.V[NI][2][k] + Program.V[NI][1][k]);

                            //velocity smoothing along border lines
                            Program.U[1][j][k] = 0.25 * (Program.U[1][j - 1][k] + Program.U[1][j + 1][k] + Program.U[2][j][k] + Program.U[1][j][k]);
                            Program.V[1][j][k] = 0.25 * (Program.V[1][j - 1][k] + Program.V[1][j + 1][k] + Program.V[2][j][k] + Program.V[1][j][k]);
                            Program.U[i][NJ][k] = 0.25 * (Program.U[i - 1][NJ][k] + Program.U[i + 1][NJ][k] + Program.U[i][NJ - 1][k] + Program.U[i][NJ][k]);
                            Program.V[i][NJ][k] = 0.25 * (Program.V[i - 1][NJ][k] + Program.V[i + 1][NJ][k] + Program.V[i][NJ - 1][k] + Program.V[i][NJ][k]);
                            Program.U[NI][j][k] = 0.25 * (Program.U[NI][j + 1][k] + Program.U[NI][j - 1][k] + Program.U[NI - 1][j][k] + Program.U[NI][j][k]);
                            Program.V[NI][j][k] = 0.25 * (Program.V[NI][j + 1][k] + Program.V[NI][j - 1][k] + Program.V[NI - 1][j][k] + Program.V[NI][j][k]);
                            Program.U[i][1][k] = 0.25 * (Program.U[i + 1][1][k] + Program.U[i - 1][1][k] + Program.U[i][2][k] + Program.U[i][1][k]);
                            Program.V[i][1][k] = 0.25 * (Program.V[i + 1][1][k] + Program.V[i - 1][1][k] + Program.V[i][2][k] + Program.V[i][1][k]);

                            //velocity smoothing inside
                            Program.U[i][j][k] = 0.2 * (Program.U[i + 1][j][k] + Program.U[i - 1][j][k] + Program.U[i][j - 1][k] + Program.U[i][j + 1][k] + Program.U[i][j][k]);
                            Program.V[i][j][k] = 0.2 * (Program.V[i + 1][j][k] + Program.V[i - 1][j][k] + Program.V[i][j - 1][k] + Program.V[i][j + 1][k] + Program.V[i][j][k]);
                        }
                    }
                }

                //round values
                Parallel.For(1, NI + 1, Program.pOptions, i =>
                {
                    for (int j = 2; j <= NJ - 1; j++)
                    {
                        for (int k = 1; k <= NK; k++)
                        {
                            Program.U[i][j][k] = Math.Round(Program.U[i][j][k], 3);
                            Program.V[i][j][k] = Math.Round(Program.V[i][j][k], 3);
                        }
                    }
                });

                //check if each cell is assigned a wind speed
                M = 0;
                O = 0;
                for (int i = 1; i <= NI; i++)
                {
                    for (int j = 1; j <= NJ; j++)
                    {
                        for (int k = 1; k <= NK; k++)
                        {
                            if ((Program.U[i][j][k] == 0) && (Program.V[i][j][k] == 0)) M++;
                            O++;
                        }
                    }
                }
                Console.WriteLine("NUMBER OF UNDEFINED WIND SPEEDS           : " + M.ToString());
                Console.WriteLine("TOTAL NUMBER OF CELLS: " + O.ToString());
                if (M != 0)
                {
                    Console.WriteLine("INTERPOLATION FAILED!!");
                    Environment.Exit(0);
                }
            }

            //pressure for the model top
            if ((Program.METEO == "Y") || (Program.METEO == "y"))
            {
                PUNTEN = 99000;
                PMEER = 101300;
            }
            else
            {
                Console.WriteLine();
                Console.Write("	SURFACE PRESSURE OF LOWEST POINT IN MODEL =  ");
                PUNTEN = Convert.ToDouble(Console.ReadLine().Replace(".", Program.decsep));
                Console.WriteLine(" ");
                Console.Write("   SURFACE PRESSURE AT SEA LEVEL            =  ");
                PMEER = Convert.ToDouble(Console.ReadLine().Replace(".", Program.decsep));
            }
            Program.POBEN = PMEER * Math.Exp(-Program.ZSP[2][2][NK] / 8000);

        GOTO1234:

            Parallel.For(1, NI + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= NJ; j++)
                {
                    for (int k = 1; k <= NK; k++)
                    {
                        //harmonic mean temperature
                        double THARM = 0;
                        for (int m = k; m <= NK - 1; m++)
                        {
                            THARM += (Program.ZSP[i][j][m + 1] - Program.ZSP[i][j][m]) / (Program.T[i][j][m] + Program.T[i][j][m + 1]) * 2;
                        }
                        //pressure profile
                        Program.PBZ[i][j][k] = (float)(Program.POBEN * Math.Exp(Program.GERD / Program.GASCON * THARM));
                    }
                }
            });

            if (Math.Abs((Program.PBZ[Program.AHMINI][Program.AHMINJ][1] - PUNTEN) / PUNTEN * 100) >= 0.01)
            {
                if (Program.PBZ[Program.AHMINI][Program.AHMINJ][1] > PUNTEN)
                {
                    Program.POBEN -= (Program.PBZ[Program.AHMINI][Program.AHMINJ][1] - PUNTEN) * 0.1;
                    goto GOTO1234;
                }
                else if (Program.PBZ[Program.AHMINI][Program.AHMINJ][1] < PUNTEN)
                {
                    Program.POBEN -= (Program.PBZ[Program.AHMINI][Program.AHMINJ][1] - PUNTEN) * 0.1;
                    goto GOTO1234;
                }
            }
            TMAX = 0;
            double TMIN = 374;
            for (int i = 1; i <= NI; i++)
            {
                for (int j = 1; j <= NJ; j++)
                {
                    for (int k = 1; k <= NK; k++)
                    {
                        //harmonic mean temperature
                        double THARM = 0;
                        for (int m = k; m <= NK - 1; m++)
                        {
                            THARM += (Program.ZSP[i][j][m + 1] - Program.ZSP[i][j][m]) / (Program.T[i][j][m] + Program.T[i][j][m + 1]) * 2;
                        }

                        //pressure profile
                        Program.PBZ[i][j][k] = (float)(Program.POBEN * Math.Exp(Program.GERD / Program.GASCON * THARM));
                        if (Program.ISTAT == 1)
                        {
                            Program.FACTOR[i][j][k] = (float)(Math.Pow(PMEER / Program.PBZ[i][j][k], 0.287));
                        }
                        else
                        {
                            Program.FACTOR[i][j][k] = (float)(Math.Pow(PMEER / Program.PBZ[i][j][k], 0.287 * moist_adiabatic));
                        }
                        Program.PBZ[i][j][k] = (float)(Math.Round(Program.PBZ[i][j][k], 1));

                        //absolute temperature
                        Program.T[i][j][k] = Math.Round(Program.T[i][j][k], 2);
                        Program.TABS[i][j][k] = (float)(Program.T[i][j][k]);
                        if (Program.TABS[i][j][k] < TMIN) TMIN = Program.TABS[i][j][k];
                        if (Program.TABS[i][j][k] > TMAX) TMAX = Program.TABS[i][j][k];

                        //potential temperature
                        Program.T[i][j][k] *= Program.FACTOR[i][j][k];

                        //convective initialisation needs to be avoided due to numerical reasons
                        Program.T[i][j][k] = Math.Max(Program.T[i][j][k], Program.T[i][j][k - 1]);

                        //basic state of potential temperature
                        Program.TBZ[i][j][k] = Program.T[i][j][k];

                        //pressure profile for radiation model
                        Program.PBZZ[k] = Program.PBZ[2][2][k];

                        //densitiy profile for radiation model
                        Program.RHOBZ[i][j][k] = (float)(Program.PBZ[i][j][k] / Program.GASCON / Program.TBZ[i][j][k] / Program.FACTOR[i][j][k]);
                        Program.RHOBZZ[k] = Program.RHOBZ[2][2][k];
                    }
                }
            }

            Int32 I = Program.AHMINI;
            Int32 J1 = Program.AHMINJ;
            Console.WriteLine();
            for (int k = NK; k >= 1; k--)
            {
                Console.WriteLine("  HEIGHT : " + Convert.ToString(Math.Round(Program.ZSP[I][J1][k], 1).ToString("0.0")).PadLeft(6) +
                    "   Tpot = " + Convert.ToString(Math.Round(Program.T[I][J1][k], 2).ToString("0.00")).PadLeft(6) +
                    "   T = " + Convert.ToString(Math.Round(Program.TABS[I][J1][k], 2).ToString("0.00")).PadLeft(6) +
                    "   P = " + Convert.ToString(Math.Round(Program.PBZ[I][J1][k], 2).ToString("0")).PadLeft(9) +
                    "   RHO = " + Convert.ToString(Math.Round(Program.RHOBZ[I][J1][k], 4).ToString("0.000")).PadLeft(6));
            }

            Console.WriteLine();
            Console.WriteLine("   MAXIMUM TEMPERATURE : " + Convert.ToString(Math.Round(TMAX, 2)).PadLeft(6));
            Console.WriteLine("   MINIMUM TEMPERATURE : " + Convert.ToString(Math.Round(TMIN, 2)).PadLeft(6));
            Console.WriteLine();
            if (TMIN < 120)
            {
                Console.WriteLine(" TEMPERATURE BELOW 120K !!");
                Environment.Exit(0);
            }

            //Cloud scheme for radiation model 1=thin / 2=thick clouds
            //Program.ISOL = 1;

            //Profiles of the radiation model
            Int32 NPROF = 30;
            Program.ZPROF[1] = Math.Max(Program.AHMIN, 0);
            for (int n = 1; n <= NPROF; n++)
            {
                Int32 INDO = 0;
                Int32 INDU = 0;
                if (n > 1) Program.ZPROF[n] = Program.ZPROF[n - 1] + 500;
                if (Program.ZPROF[n] > Program.ZSP[Program.AHMINI][Program.AHMINJ][NK])
                {
                    //pressure profile above GRAMM domain
                    Program.TPROF[n] = Program.TBZ[Program.AHMINI][Program.AHMINJ][NK];
                    //temperature profile above GRAMM domain
                    Program.PPROF[n] = PMEER * Math.Exp(-Program.ZPROF[n] / 8000);
                }
                else if (Program.ZPROF[n] <= Program.ZSP[Program.AHMINI][Program.AHMINJ][NK])
                {
                    //temperature and pressure profile within GRAMM domain
                    UNTU = -100000;
                    UNTO = 100000;
                    INDO = 0;
                    INDU = 0;
                    double GRAD = 0;
                    double GRADP = 0;
                    for (int m = 1; m <= NK; m++)
                    {
                        DIFF = Program.ZSP[Program.AHMINI][Program.AHMINJ][m] - Program.ZPROF[n];
                        if ((DIFF >= 0) && (DIFF < UNTO))
                        {
                            UNTO = DIFF;
                            INDO = m;
                        }
                        if ((DIFF < 0) && (DIFF > UNTU))
                        {
                            UNTU = DIFF;
                            INDU = m;
                        }
                    }
                    if (INDO == 0)
                    {
                        GRAD = (Program.TBZ[Program.AHMINI][Program.AHMINJ][INDU] - Program.TBZ[Program.AHMINI][Program.AHMINJ][INDU - 1]) /
                            (Program.ZSP[Program.AHMINI][Program.AHMINJ][INDU] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDU - 1]);
                        GRADP = (Program.PBZ[Program.AHMINI][Program.AHMINJ][INDU] - Program.PBZ[Program.AHMINI][Program.AHMINJ][INDU - 1]) /
                            (Program.ZSP[Program.AHMINI][Program.AHMINJ][INDU] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDU - 1]);
                        Program.TPROF[n] = Program.TBZ[Program.AHMINI][Program.AHMINJ][INDU] + (Program.ZPROF[n] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDU]) * GRAD;
                        Program.PPROF[n] = Program.PBZ[Program.AHMINI][Program.AHMINJ][INDU] + (Program.ZPROF[n] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDU]) * GRADP;
                    }
                    else if (INDU == 0)
                    {
                        GRAD = (Program.TBZ[Program.AHMINI][Program.AHMINJ][INDO] - Program.TBZ[Program.AHMINI][Program.AHMINJ][INDO + 1]) /
                            (Program.ZSP[Program.AHMINI][Program.AHMINJ][INDO] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDO + 1]);
                        GRADP = (Program.PBZ[Program.AHMINI][Program.AHMINJ][INDO] - Program.PBZ[Program.AHMINI][Program.AHMINJ][INDO + 1]) /
                            (Program.ZSP[Program.AHMINI][Program.AHMINJ][INDO] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDO + 1]);
                        Program.TPROF[n] = Program.TBZ[Program.AHMINI][Program.AHMINJ][INDO] + (Program.ZPROF[n] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDO]) * GRAD;
                        Program.PPROF[n] = Program.PBZ[Program.AHMINI][Program.AHMINJ][INDO] + (Program.ZPROF[n] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDO]) * GRADP;
                    }
                    else
                    {
                        GRAD = (Program.TBZ[Program.AHMINI][Program.AHMINJ][INDO] - Program.TBZ[Program.AHMINI][Program.AHMINJ][INDU]) /
                            (Program.ZSP[Program.AHMINI][Program.AHMINJ][INDO] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDU]);
                        GRADP = (Program.PBZ[Program.AHMINI][Program.AHMINJ][INDO] - Program.PBZ[Program.AHMINI][Program.AHMINJ][INDU]) /
                            (Program.ZSP[Program.AHMINI][Program.AHMINJ][INDO] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDU]);
                        Program.TPROF[n] = Program.TBZ[Program.AHMINI][Program.AHMINJ][INDO] + (Program.ZPROF[n] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDO]) * GRAD;
                        Program.PPROF[n] = Program.PBZ[Program.AHMINI][Program.AHMINJ][INDO] + (Program.ZPROF[n] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDO]) * GRADP;
                    }
                }
                //compute absolute temperature based on absolute temperature
                Program.TPROF[n] /= Math.Pow(PMEER / Program.PPROF[n], 0.287);
                if (Program.TPROF[n] < 153.15) Program.TPROF[n] = 153.15;
                Program.VNORM[n] = 36.5419617 + 4.8939118 * (Program.ZPROF[n] * 0.001) +
                    4.1091542 * Math.Pow(Program.ZPROF[n] * 0.001, 2) - 0.1456879 * Math.Pow(Program.ZPROF[n] * 0.001, 3) + 0.0149291 * Math.Pow(Program.ZPROF[n] * 0.001, 4);
                Program.VNORM[n] *= 1000;

                //compute absolute humidity based on specific humidity
                double TBZN = Program.TPROF[n] - 153.15;
                TBZN = Math.Max(0, Math.Min(209.0, TBZN));
                int TBZNINT = Convert.ToInt32(Math.Floor(TBZN));
                double PDST = Program.PSAT[TBZNINT + 1] + (Program.PSAT[TBZNINT + 2] - Program.PSAT[TBZNINT + 1]) * (TBZN - (float)TBZNINT);

                //water vapour in the atmosphere for the radiation model
                Program.QVAP[n] = 18.02 / 28.96 * PDST / (Program.PPROF[n] / Program.QUINIT - PDST);

                //water content in clouds
                Program.QCLD[n] = 0;

                //water content in rain
                Program.QRAIN[n] = 0;

                //water content in ice
                Program.QICE[n] = 0;

                //water content in snow
                Program.QSNOW[n] = 0;
            }
        }
    }
}
