#region Copyright
///<remarks>
/// <GRAMM Mesoscale Model>
/// Copyright (C) [2019] [Dietmar Oettl, Markus Kuntner]
/// This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
/// the Free Software Foundation version 3 of the License
/// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
/// You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
///</remarks>
#endregion

using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Threading.Tasks;

namespace GRAMM_2001
{
    /* GRAZ MESOSCALE MODELL GRAMM
    COMPREHENSIVE DESCRIPTION CAN BE FOUND IN OETTL, 2016
    THE GRAMM MODEL HAS BEEN ORIGINALLY DEVELOPED BY RAIMUND ALMBAUER AROUND 1990.
    THE RADIATION MODEL IS BASED ON THE MODEL SUGGESTED BY SOMIESKI 1988 AND HAS BEEN IMPLEMENTED AROUND 1990 BY ROBERT ALMER
    THE RE-INITIALIZATION PROCEDURE HAS BEEN IMPLEMENTED BY ULRICH UHRNER IN 2013
    THE FOLLOWING MODEL PARTS HAVE BEEN IMPLEMENTED BY DIETMAR OETTL SINCE 1997:
    -DIAGNOSTIC MODEL FOR INITIALIZATION OF WIND AND TEMPERATURE FIELDS USING POINT AND PROFILE OBSERVATIONS
    -FULL IMPLICIT SOLVERS OF ALL CONSERVATION EQUATIONS USING TDMA ALGORITHM TO IMPROVE NUMERICAL STABILITIES IN HIGHLY COMPLEX TERRAIN
    -STANDARD K-EPS MODEL
    -PARALLELISATION
    -PORTING THE CODE FROM FORTRAN TO C#
    */

    partial class Program
    {
        /*INPUT FILES :
        BASIC DOMAIN INFORMATION GRAMM.geb
        MAIN CONTROL PARAM FILEs iin.dat or IIN.dat
        GEOMETRY DATA ggeom.asc
        LANDUSE DATA landuse.asc
        METEOROLOGICAL DATA meteopgt.all
        MAX. NUMBER OF CPUs Max_Proc.txt
        */

        static void Main(string[] args)
        {
            GdalConfiguration.ConfigureOgr();
            GdalConfiguration.ConfigureGdal();

            // Arguments: "First Situation" "Final Situation" "Max. Time Step" "RelaxV" "RelaxT"
            // if "First Situation" = "?" -> print Info & Stop

            int p = (int)Environment.OSVersion.Platform;

            if ((p == 4) || (p == 6) || (p == 128))
            {
                //Console.WriteLine ("Running on Unix");
                unix = true;
            }
            else
            {
                //Console.WriteLine ("NOT running on Unix");
            }

            //WRITE GRAMM VERSION INFORMATION TO SCREEN
            Console.WriteLine("");
            Console.WriteLine("+------------------------------------------------------+");
            Console.WriteLine("|                                                      |");
            string Info =     "+         > > G R A M M VERSION: 21.01Alpha < <        +";
            Console.WriteLine(Info);
            if (unix)
            {
                Console.WriteLine("|                      L I N U X                     |");
            }
#if NETCOREAPP2_1 || NETCOREAPP2_0 || NETCOREAPP3_0
Console.WriteLine("| .Net Core Version |");
#endif
            Console.WriteLine("+------------------------------------------------------+");
            Console.WriteLine(" ");

            //show licence terms
            ShowCopyright(args);

            // 11.04.17 Ku use arguments
            Console_Arguments(args);

            //User defined decimal seperator
            decsep = NumberFormatInfo.CurrentInfo.NumberDecimalSeparator;

            //read number of grid cells stored in the file "GRAMM.geb"
            Read_Gramm_Geb();

            //check if chemistry is envoked
            if (File.Exists("chemistry.txt"))
            {
                Read_Chemistry();
                //@Johannes: at this place call the chemical mechanism, that feeds back the number of spezies to be computed
            }

            // Write to "Logfile_GRAMMCore"
            try
            {
                ProgramWriters.LogfileGrammCoreWrite(new String('-', 80));
                ProgramWriters.LogfileGrammCoreWrite(Info);
                ProgramWriters.LogfileGrammCoreWrite("Computation started at: " + DateTime.Now.ToString());
                ProgramWriters.LogfileGrammCoreWrite("Computation folder: " + Directory.GetCurrentDirectory());
                Info = "Application hash code: " + GetAppHashCode();
                ProgramWriters.LogfileGrammCoreWrite(Info);
            }
            catch { }

            //Allocate Memory -> Define Arrays
            Define_Arrays();
            List<int>[] Month_List = new List<int>[3]; // 3 types of month lists
            Month_List[0] = new List<int>() { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
            Month_List[1] = new List<int>() { 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 2 };
            Month_List[2] = new List<int>() { 6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5 };

            if (NX * NY < 50 * 50)
            {
                TerminalThreshold = 12;
            }
            else if (NX * NY < 200 * 200)
            {
                TerminalThreshold = 6;
            }
            else
            {
                TerminalThreshold = 1;
            }

            Console.WriteLine(".");
            //END OF VARIABLE DECLARATION BLOCK

            //SET NUMBER OF FLOW FIELD SITUATION TO ZERO
            IWETTER = 0;

            //OPEN THE MAIN CONTROL FILE "IIN.dat"
            Read_IIN_Dat();

            //When using ERA5 data GRAMM automatically starts with a number for the weather situation based on the existing *wnd files
            if (ISTAT >= 2)
            {
                string wdir = Directory.GetCurrentDirectory();
                DirectoryInfo d = new DirectoryInfo(wdir);
                FileInfo[] Files = d.GetFiles("*.wnd");
                IWETTER = Files.Length;
            }

            //Total simulation time
            DTI = TLIMIT2;
            DTMAX = DT;
            max_time_step_original = DTMAX;

            if (MaxTimeStep_Console > 0.99)
                max_time_step_original = MaxTimeStep_Console;
            if (Relaxv_Console > 0.0099)
                RELAXV = Relaxv_Console;
            if (Relaxt_Console > 0.0099)
                RELAXT = Relaxt_Console;

            // Read number of max. used processors
            IPROC = 1;
            Max_Proc_File_Read();

            //Convert relative humidity in %
            QUINIT *= 0.01;

            //Get year, month, day, hour, and minute of the simulation start
            if (NDATUM.ToString().Length == 8)
            {
                IJAHR4digits = Convert.ToInt32(Math.Floor(NDATUM / 10000d));
                IJAHR = IJAHR4digits - Convert.ToInt32(Math.Floor(IJAHR4digits / 100d)) * 100;
                IMON = Convert.ToInt32(Math.Floor(NDATUM / 100d)) - IJAHR4digits * 100;
                ITAG = NDATUM - IJAHR4digits * 10000 - IMON * 100;
            }
            else
            {
                IJAHR = Convert.ToInt32(Math.Floor(NDATUM / 10000d));
                IJAHR4digits = IJAHR;
                IMON = Convert.ToInt32(Math.Floor(NDATUM / 100d)) - IJAHR * 100;
                ITAG = NDATUM - IJAHR * 10000 - IMON * 100;
            }
            ISTU = Convert.ToInt32(Math.Floor(ISTUNDE / 100d));
            IMIN = ISTUNDE - 100 * ISTU;

            //Set flags for various compuation options
            ICU = false;
            ICV = false;
            ICW = false;
            ICPN = false;
            ICT = false;
            ICPH = false;
            IFOU = false;
            ICQU = false;
            ICPSI = false;
            ICTE = false;
            ICSTR = false;
            ICPR = false;
            ICBR = false;
            ICTB = false;
            ICGW = false;

            Int32 INUMM = 0;
            if (IFLAGS1 / 1000000 == 1) ICU = true;
            if (ICU) INUMM += 1000000;
            if ((IFLAGS1 - INUMM + 1) / 100000 == 1) ICV = true;
            if (ICV) INUMM += 100000;
            if ((IFLAGS1 - INUMM + 1) / 10000 == 1) ICW = true;
            if (ICW) INUMM += 10000;
            if ((IFLAGS1 - INUMM + 1) / 1000 == 1) ICPN = true;
            if (ICPN) INUMM += 1000;
            if ((IFLAGS1 - INUMM + 1) / 100 == 1) ICT = true;
            if (ICT) INUMM += 100;
            if ((IFLAGS1 - INUMM + 1) / 10 == 1) ICGW = true;
            if (ICGW) INUMM += 10;
            if ((IFLAGS1 - INUMM) == 1) IFOU = true;

            INUMM = 0;
            if (IFLAGS2 / 1000000 == 1) ICBR = true;
            if (ICBR) INUMM += 1000000;
            if ((IFLAGS2 - INUMM + 1) / 100000 == 1) ICPR = true;
            if (ICPR) INUMM += 100000;
            if ((IFLAGS2 - INUMM + 1) / 10000 == 1) ICQU = true;
            if (ICQU) INUMM += 10000;
            if ((IFLAGS2 - INUMM + 1) / 1000 == 1) ICPSI = true;
            if (ICPSI) INUMM += 1000;
            if ((IFLAGS2 - INUMM + 1) / 100 == 1) ICTE = true;
            if (ICTE) INUMM += 100;
            if ((IFLAGS2 - INUMM + 1) / 10 == 1) ICTB = true;
            if (ICTB) INUMM += 10;
            if ((IFLAGS2 - INUMM) == 1) ICSTR = true;

            //COMPUTE MODEL GRID
            Console.WriteLine(" *** GENERATING MODEL GRID *** ");
            GEOM();

            //INQUIRE IF RECEPTOR POINTS ARE SET
            Read_Receptor_Dat();

            //Analyze_Topography(); // find U valleys and bassins

            Relaxv_ori = RELAXV;
            Relaxt_ori = RELAXT;
            Relax_Border_factor[0][0] = -4; // flag, that factor must be computed

            ProgramWriters.LogfileGrammCoreInfo();

            // init radiation model 
            RadiationModel = new RadiationCalculation(); 

        //Loop_______________________________________________________________________________________________
        NEXTWEATHERSITUATION:

            clear_arrays();

            //INITIALIZE FIELDS
            Console.WriteLine();
            Console.WriteLine();
            Console.WriteLine(" *** INITIALIZATION *** ");
            INITA(NX, NY, NZ);

            Int32 ISTUD = Convert.ToInt32(Math.Floor(TLIMIT / 3600d));
            double AMIND = (TLIMIT / 3600 - (float)ISTUD) * 60;
            Int32 IMIND = Convert.ToInt32(Math.Floor(AMIND));
            double ASECD = (AMIND - (float)IMIND) * 60;
            Int32 ISECD = Convert.ToInt32(Math.Floor(ASECD));

            if (ISTAT >= 2)
            {
                //setting UTC time for reading ERA5 data
                dateUTC = new DateTime(Program.IJAHR4digits, Program.IMON, Program.ITAG, Program.ISTU, Program.IMIN, 0);
                //GRAMM uses sun time -> UTC has to be transferred
                if (Longitude > 180.0)
                    Longitude -= 180;
                if (Longitude < -180.0)
                    Longitude += 180;
                ISTU += Convert.ToInt32(Longitude * 12 / 180.0);
                if (ISTU >= 24)
                {
                    ISTU -= 24;
                    DateTime dateref = new DateTime(Program.IJAHR4digits, Program.IMON, Program.ITAG, Program.ISTU, Program.IMIN, 0);
                    dateref = dateref.AddDays(1);
                    ITAG = dateref.Day;
                    IMON = dateref.Month;
                }
            }

            if (ISTAT == 0)
            {
                Intermed_Threshold = IOUTPUT / 2; // 1st intermediate output after IOUTPUT / 2 (default 1800 s) for SC 5, 6 and SC 7
                if (AKLA < 5)
                {
                    Intermed_Threshold = IOUTPUT; // 1st intermediate output after IOUTPUT s
                }
            }
            else
            {
                Intermed_Threshold = IOUTPUT;
            }
            int Radiation_Threshold = 1800; // compute new radiation each 1800 s

            //set further initial values
            if (ISTAT < 2)
                INITB.INIT(NX, NY, NZ);

            Relax_Border(); // set reduction factors for the border cells

            //initialize the radiation model
            if (((METEO != "Y") && (METEO != "y")) || (Program.ISTAT != 0)) TLIMIT = TLIMIT2;
            Boolean LGEOM = false;
            Boolean LGEOMW = false;
            Boolean LGEOMR = false;
            INIT_Radiation(ref LGEOM, ref LGEOMW, ref LGEOMR, ref ISTUD, ref AMIND, ref IMIND, ref ASECD, ref ISECD, Month_List);

            Console.WriteLine(ITAG.ToString() + "." + IMON.ToString() + " - " + ISTU.ToString() + ":" + IMIND.ToString("D2"));

            REALTIME = 0;
            ITIME = 0;
            INUMS = 0;
            DIVSUM = 0;
            IDIV = 0;
            DT = 1.5;
            for (int i = 0; i < 11; i++)
            {
                MASSOURCE[i] = 0;
            }

            //START OF THE INTEGRATION LOOP
            while (REALTIME < TLIMIT)
            {
                //number of iteration
                ITIME++;

                if (ITIME % 20 == 0) Max_Proc_File_Read(); // read MaxProc at 20th time step

                //total simulation expressed in seconds
                TJETZT = (float)ISTU * 3600 + (float)IMIN * 60 + REALTIME;

                //write normalised actual time used in the GUI for visualisation of the modelling progress
                if ((ITIME % 10) == 0 && (IWetter_Console_First <= 1)) // 11.4.17 Ku use arguments
                {
                    try
                    {
                        using (StreamWriter w = new StreamWriter("PercentGramm.txt"))
                        {
                            double tnorm = REALTIME * TLIMIT2 / TLIMIT;
                            w.WriteLine(tnorm.ToString());
                        }
                    }
                    catch { }
                }


                //Online output of fields for GUI
                if (GRAMM_Online_flag || (ITIME % 10d) == 0)
                {
                    GRAMM_Online_flag = false; // reset flag
                    GrammOnline(NX, NY, NZ); // flag is at GRAMMOnline set to true, if only output is necessary
                }

                //Implicit solution algorith for conservation equations (SIMPLE - Patankar, 1980)
                if (REALTIME <= DTI)
                {
                    if (SOLUTION(NX, NY, NZ) == false && ISTAT == 0) // numerical problems using overall massdivergence
                    {
                        computation_retry++;
                        if (computation_retry < 3) // try 3 times -> otherwise let the app crash
                        {
                            IWETTER--; //try same situation
                            TLIMIT += TLIMIT2;
                            DTI += TLIMIT2;
                            clear_arrays();
                            GEOM();
                            goto NEXTWEATHERSITUATION;
                        }

                    }
                }

                //new radiation data
                if (ITIME == 1)
                {
                    LGEOM = true;
                    LGEOMW = true;
                    LGEOMR = false;
                }
                else
                {
                    LGEOM = false;
                    LGEOMW = false;
                    LGEOMR = false;
                }

                if (ICSTR == true)
                {
                    if ((ITIME == 1) || ((ITIME % IRAD) == 0))
                    {
                        double TJETZT1 = TJETZT;
                        int IMIN_RAD = IMIN;
                        int ITAG_RAD = ITAG;
                        int IMON_RAD = IMON;
                        if (ISTAT == 0)
                        {
                            ISTUD = Convert.ToInt32(Math.Floor((TJETZT + 10) / 3600d));
                        }
                        else if (ISTAT == 1)
                        {
                            ISTUD = ISTU;
                            TJETZT1 = (float)ISTU * 3600 + (float)IMIN * 60;
                        }
                        else if (ISTAT >= 2)
                        {
                            //ISTUD = Convert.ToInt32(Math.Floor((TJETZT + 10) / 3600d));
                            //GRAMM uses sun time -> UTC has to be transferred
                            DateTime dateref = new DateTime(Program.IJAHR4digits, Program.IMON, Program.ITAG, Program.ISTU, Program.IMIN, 0);
                            dateref = dateref.AddSeconds(REALTIME);
                            IMON_RAD = dateref.Month;
                            ITAG_RAD = dateref.Day;
                            ISTUD = dateref.Hour;
                            IMIND = dateref.Minute;
                            TJETZT1 = ISTUD * 3600;
                        }
                        AMIND = (TJETZT / 3600 - (float)ISTUD) * 60;
                        IMIND = Convert.ToInt32(Math.Floor(AMIND));
                        ASECD = (AMIND - (float)IMIND) * 60;
                        ISECD = Convert.ToInt32(Math.Floor(ASECD));

                        for (int II = 1; II <= 1000; II++)
                        {
                            if (TJETZT1 >= 86400) TJETZT1 -= 86400;
                        }
                        RadiationModel.RADIATRAD(LGEOM, LGEOMW, LGEOMR, ITAG_RAD, IMON_RAD, IJAHR, TJETZT1, NX, NY, NZ);
                        Console.WriteLine(ITAG_RAD.ToString() + "." + IMON_RAD.ToString() + " - " + ISTUD.ToString() + ":" + IMIND.ToString("D2"));
                    }

                    //dynamic sun (global radiation) every 1800 seconds
                    if (((METEO == "Y") || (METEO == "y")) && (Program.ISTAT == 0))
                    {
                        if (Program.AKLA < 4)
                        {
                            if (REALTIME > Radiation_Threshold)
                            {
                                Radiation_Threshold += 1800; // compute new radiation each 1800 s
                                double TJETZT1 = TJETZT;
                                ISTUD = Convert.ToInt32(Math.Floor((TJETZT + 10) / 3600d));
                                AMIND = (TJETZT / 3600 - (float)ISTUD) * 60;
                                IMIND = Convert.ToInt32(Math.Floor(AMIND));
                                ASECD = (AMIND - (float)IMIND) * 60;
                                ISECD = Convert.ToInt32(Math.Floor(ASECD));

                                for (int II = 1; II <= 1000; II++)
                                {
                                    if (TJETZT1 >= 86400) TJETZT1 -= 86400;
                                }
                                RadiationModel.RADIATRAD(LGEOM, LGEOMW, LGEOMR, ITAG, IMON, IJAHR, TJETZT1, NX, NY, NZ);
                                Console.WriteLine(ITAG.ToString() + "." + IMON.ToString() + " - " + ISTUD.ToString() + ":" + IMIND.ToString("D2"));
                            }
                        }
                    }
                }

                //store last time step
                STOREcalculate(NX, NY, NZ);

                //intermediate output
                if (((METEO == "Y") || (METEO == "y")) && (Program.ISTAT == 0))
                {
                    //only in case of steady-state simulations using meteopgt.all no intermediate output is provided
                    if (Program.meteopgt_nr == 0)
                    {
                        IOUTPUT = 100000000;
                        Intermed_Threshold = 100000000;
                    }
                }

                Int16 IHOUR = Convert.ToInt16(Math.Floor(TJETZT / IOUTPUT));
                if (REALTIME > Intermed_Threshold && DTI > 3602) // if threshold is exceeded && simulation time > 1 h -> create intermed. outputs
                {
                    Console.Write(" INTERMEDIATE OUTPUT ");

                    // set new Threshold value
                    if (Intermed_Threshold < IOUTPUT)
                    {
                        Intermed_Threshold = IOUTPUT * 2; // 2nd intermediate output after 2 hours
                    }
                    else
                    {
                        Intermed_Threshold += IOUTPUT;
                    }

                    OUTPUT(NX, NY, NZ, true); // intermediate output

                    //generating meteopgt.all and mettimeseries.dat for ERA5/GUI postprocessing
                    if (Program.ISTAT == 2 || Program.ISTAT == 4)
                    {
                        Meteopgtall.GRALGenerateMeteoFilesForERA5();
                    }

                    Console.WriteLine();
                    Console.WriteLine();

                    if (((METEO != "Y") && (METEO != "y")) || (Program.ISTAT != 0))
                    {
                        IWETTER++;
                    }
                }

                //increase simulation time
                REALTIME += DT;
                IHOURO = IHOUR;

                //save intermediate surface flow fields after 75% of the total simulatin time (VDI 3783-7)
                if ((REALTIME >= DTI * 0.75) && (REALTIME <= (DTI * 0.75 + DT * 2)))
                {
                    Parallel.For(1, NX + 1, Program.pOptions, i =>
                    {
                        for (int j = 1; j <= NY; j++)
                        {
                            U_TEMP[i][j] = (float)(U[i][j][1]);
                            V_TEMP[i][j] = (float)(V[i][j][1]);
                            W_TEMP[i][j] = (float)(W[i][j][1]);
                        }
                    });
                }
            }

            computation_retry = 0; // reset compuation retry counter
                                   //Ultimate output at the end of each situation
            Console.WriteLine("");
            Console.Write(" MMAIN : OUT ");
            if ((METEO != "Y") && (METEO != "y"))
            {
                OUTPUT(NX, NY, NZ, false); // final output
                Console.WriteLine();
            }
            else
            {
                OUTPUT(NX, NY, NZ, false); // final output
                Console.WriteLine();
                goto NEXTWEATHERSITUATION;
            }
        }

        //module to initialze a jagged array
        public static T[] CreateArray<T>(int cnt, Func<T> itemCreator)
        {
            T[] result = new T[cnt];
            for (int i = 0; i < result.Length; i++)
            {
                result[i] = itemCreator();
            }
            return result;
        }

        //licence GPL 3 terms
        static void ShowCopyright(string[] args)
        {
            Console.WriteLine("[GRAMM] Copyright (C) <2019> <Dietmar Oettl, Markus Kuntner>");
            Console.WriteLine("This program comes with ABSOLUTELY NO WARRANTY; for details start GRAMM with a startup parameter ‘show_w’");
            Console.WriteLine("This is free software, and you are welcome to redistribute it under certain conditions; start GRAMM with a startup parameter ‘show_c’ for details. )");

            if (args.Length > 0 && args[0].Contains("show_w"))
            {
                Console.WriteLine("This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of" +
                " MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.");
            }
            else if (args.Length > 0 && args[0].Contains("show_c"))
            {
                Console.WriteLine("This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by" +
                "the Free Software Foundation, either version 3 of the License.");
                Console.WriteLine("You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.");
            }
            Console.WriteLine();
        }

        private static bool Console_Arguments(string[] args)
        {
            int _off = 0;
            if (args.Length > 0)
            {
                if (args[0].Contains("?") == true || args[0].ToUpper().Contains("HELP") == true) // Info & stop
                {
                    Console.WriteLine("GRAMM console arguments: 'Working Directory' 'First Situation' 'Final Situation' 'Max. Time Step' 'RelaxV' 'RelaxT'");
                    Console.WriteLine("or");
                    Console.WriteLine("GRAMM console arguments: 'First Situation' 'Final Situation' 'Max. Time Step' 'RelaxV' 'RelaxT'");
                    Environment.Exit(0);
                }

                int temp;
                if (Int32.TryParse(args[0], out temp) == false) // not a valid number
                {
                    if (Directory.Exists(args[0]) == true) // arg[0] = Directory!
                    {
                        Directory.SetCurrentDirectory(args[0]);
                        _off = 1;
                    }
                }
            }
            if (args.Length > (1 + _off)) // + 2 arguments -> first and last weather situation
            {
                if (Int32.TryParse(args[0 + _off], out IWetter_Console_First))
                    Int32.TryParse(args[1 + _off], out IWetter_Console_Last);
                if (IWetter_Console_Last < IWetter_Console_First || IWetter_Console_First < 1)
                {
                    IWetter_Console_First = 0;
                    IWetter_Console_Last = 9999999;
                }
                // Max. Time Step
                if (args.Length > (2 + _off))
                    if (Double.TryParse(args[2 + _off], NumberStyles.AllowDecimalPoint, CultureInfo.InvariantCulture, out MaxTimeStep_Console) == false)
                        MaxTimeStep_Console = 0;
                // RelaxV
                if (args.Length > (3 + _off))
                    if (Double.TryParse(args[3 + _off], NumberStyles.AllowDecimalPoint, CultureInfo.InvariantCulture, out Relaxv_Console) == false)
                        Relaxv_Console = 0;
                // RelaxT
                if (args.Length > (4 + _off))
                    if (Double.TryParse(args[4 + _off], NumberStyles.AllowDecimalPoint, CultureInfo.InvariantCulture, out Relaxt_Console) == false)
                        Relaxt_Console = 0;
            }

            if (IWetter_Console_First > 0) // 11.04.17 Ku use arguments
            {
                Console.WriteLine("");
                Console.WriteLine("Starting with situation " + IWetter_Console_First.ToString() +
                " ...final situation " + IWetter_Console_Last.ToString() + " ");
            }
            if (IWetter_Console_First <= 1) // 11.04.17 Ku use arguments
            {
                try
                {
                    if (File.Exists("albeq.dat"))
                        File.Delete("albeq.dat");
                }
                catch { }
            }
            else
            {
                Console.Write("waiting for albeq.dat");
                while (File.Exists("albeq.dat") == false)
                {
                    System.Threading.Thread.Sleep(2000);
                    Console.Write(".");
                }
                System.Threading.Thread.Sleep(2000);
                Console.WriteLine("");
            }
            return true;
        }


        //computes wind direction
        public static double winkel(double u, double v)
        {
            double winkel = 0;

            if (v == 0)
                winkel = 90;
            else
                winkel = Math.Abs(Math.Atan(u / v)) * 180 / Math.PI;

            if ((v > 0) && (u <= 0)) winkel = 180 - winkel;
            if ((v >= 0) && (u > 0)) winkel = 180 + winkel;
            if ((v < 0) && (u >= 0)) winkel = 360 - winkel;

            return winkel;
        }

        private static int set_month_type(double windspeed, float AKLA, double Brgrad) // set start month for radiation search
        {
            int month_setting = 0; // Jan to Dec
            if (Brgrad > 0) // northern hemisphere
            {
                if (Windspeed_meteopgt > 2 || AKLA < 5) // higher wind vel. or labile or neutral SC -> start radiation search at Mar
                    month_setting = 1;
                else if (Windspeed_meteopgt > 4) // very high wind vel. -> start radiation search at Jun
                    month_setting = 2;
            }
            else // southern hemisphere
            {
                month_setting = 2; // June to Mai
                if (Windspeed_meteopgt > 2 || AKLA < 5) // higher wind vel. or labile or neutral SC -> start radiation search at Mar
                    month_setting = 1;
                else if (Windspeed_meteopgt > 4) // very high wind vel. -> start radiation search at Jan
                    month_setting = 0;
            }
            return month_setting;
        }

        /// <summary>
        /// Get a Hash code of the running app
        /// </summary>
        public static string GetAppHashCode()
        {
            string filename = System.Diagnostics.Process.GetCurrentProcess().MainModule.FileName;
            string hashstring = string.Empty;

            if (File.Exists(filename))
            {
                try
                {
                    using (System.Security.Cryptography.SHA1 sha = System.Security.Cryptography.SHA1.Create())
                    {
                        byte[] hash;

                        using (FileStream stream = File.OpenRead(filename))
                        {
                            hash = sha.ComputeHash(stream);
                        }

                        foreach (byte item in hash)
                        {
                            hashstring += item.ToString("x2");
                        }
                    }
                }
                catch { }
            }
            return hashstring;
        }

    }
}
