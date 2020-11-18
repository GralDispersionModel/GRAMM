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
using System.IO;

namespace GRAMM_2001
{
    class Meteopgtall
    {
        //computes the correct filenumber for intermediate GRAMM flow field files
        public static void meteopgtall_calculate(int Ori_Weather_nr, int IWETTER, double DTI, double TIME, int OUTPUT, double TLIMIT2, ref int Filenumber)
        {
            //compute index of intermediate file
            int IHOUR = Convert.ToInt32(Math.Floor(TIME / OUTPUT));

            //compute filenumber
            Filenumber = Ori_Weather_nr;
            List<int> counter = new List<int>();
            //reading meteopgt.all

            try
            {
                using (FileStream fs = new FileStream("meteopgt.all", FileMode.Open, FileAccess.Read, FileShare.Read))
                {
                    using (StreamReader myreader = new StreamReader(fs))
                    {
                        string[] text = new string[10];
                        text = myreader.ReadLine().Split(new char[] { ' ', ';', ',', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                        text = myreader.ReadLine().Split(new char[] { ' ', ';', ',', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                        for (int inid = 1; inid < Program.IWETTER; inid++)
                        {
                            text = myreader.ReadLine().Split(new char[] { ' ', ';', ',', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                            double WINDGE = Math.Max(Convert.ToDouble(text[1].Replace(".", Program.decsep)), 0.001);
                            int AKLA = Convert.ToInt16(text[2]);
                            double ITIME = TLIMIT2;
                            Integrationtime(AKLA, WINDGE, TLIMIT2, Ori_Weather_nr, ref ITIME);

                            //number of additional weather situations
                            int addnumb = Convert.ToInt32(Math.Ceiling(ITIME / OUTPUT - 1));
                            for (int i = 1; i <= addnumb; i++)
                            {
                                counter.Add(1);
                            }
                        }
                    }
                }
            }
            catch
            {
                Console.WriteLine("Error when reading file meteopgt.all during intermediate output of GRAMM flow fields - Execution stopped");
                Environment.Exit(0);
            }

            Filenumber += counter.Count + IHOUR;
        }

        //computes a new meteopgt.all with additional weather situations at the end of the original file, where intermediate GRAMM flow fields are stored
        public static void meteopgtall_generate(int Ori_Weather_nr, double TLIMIT2, int OUTPUT)
        {
            List<string> texttoadd = new List<string>();
            bool extended_meteofile = false;
            //reading meteopgt.all

            try
            {
                using (FileStream fs = new FileStream("meteopgt.all", FileMode.Open, FileAccess.Read, FileShare.Read))
                {
                    using (StreamReader myreader = new StreamReader(fs))
                    {

                        string[] text = new string[10];
                        text = myreader.ReadLine().Split(new char[] { ' ', ';', ',', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                        text = myreader.ReadLine().Split(new char[] { ' ', ';', ',', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                        for (int inid = 1; inid <= Ori_Weather_nr; inid++)
                        {
                            text = myreader.ReadLine().Split(new char[] { ' ', ';', ',', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                            double WINDGE = Math.Max(Convert.ToDouble(text[1].Replace(".", Program.decsep)), 0.001);
                            double WINDDIR = Convert.ToDouble(text[0].Replace(".", Program.decsep)) * 10;
                            int AKLA = Convert.ToInt16(text[2]);
                            double ITIME = TLIMIT2;
                            Integrationtime(AKLA, WINDGE, TLIMIT2, Ori_Weather_nr, ref ITIME);
                            //number of additional weather situations
                            int addnumb = Convert.ToInt32(Math.Ceiling(ITIME / OUTPUT - 1));
                            //Console.WriteLine(inid.ToString()+"/" + addnumb.ToString() +"/" + ITIME.ToString() + "/" + OUTPUT.ToString());
                            text[1] = text[1].Replace(".", string.Empty);
                            for (int i = 1; i <= addnumb; i++)
                            {
                                texttoadd.Add(Convert.ToString(Convert.ToInt32(WINDDIR)) + "." + text[1] + "," + inid.ToString() + "." + i.ToString() + "," + text[2] + ",0");
                            }
                        }
                        if (myreader.EndOfStream)
                            extended_meteofile = true;
                    }
                }
            }
            catch
            {
                Console.WriteLine("Error when reading file meteopgt.all during intermediate output of GRAMM flow fields - Execution stopped");
                Environment.Exit(0);
            }

            //writing extended meteopgt.all, only if it has not already been done before
            if (extended_meteofile == true)
            {
                try
                {
                    using (StreamWriter mywriter = File.AppendText("meteopgt.all"))
                    {
                        for (int inid = 0; inid < texttoadd.Count; inid++)
                        {
                            mywriter.WriteLine(texttoadd[inid]);
                        }
                    }
                }
                catch
                {
                    Console.WriteLine("Error when writing extended file meteopgt.all for intermediate output of GRAMM flow fields - Execution stopped");
                    Environment.Exit(0);
                }
            }
        }

        //computes the GRAMM integration time in dependence on wind speed and stability class
        public static void Integrationtime(int stability, double windspeed, double TLIMIT2, int dynamic, ref double DTI)
        {
            //New time scheme OETTL, 10 Feb 2017 to increase wind field variability
            DTI = TLIMIT2;

            if (stability < 4)
            {
                if (windspeed < 3)
                {
                    DTI = TLIMIT2 * 1.0;
                }
                else if (windspeed < 5)
                {
                    DTI = TLIMIT2 * 0.5;
                }
                else if (windspeed < 7)
                {
                    DTI = TLIMIT2 * 0.3;
                }
                else
                {
                    DTI = TLIMIT2 * 0.17;
                }

                //in case of dynamic sun rise simulation
                if (dynamic > 0)
                {
                    DTI = 21600;  //should be exactly 6 hours as the simulation starts at 6 o'clock in the morning
                }
            }
            if (stability == 6)
            {
                if (windspeed <= 2.0)
                {
                    //in case of dynamic sun rise simulation
                    if (dynamic > 0)
                    {
                        DTI = TLIMIT2 * 3.0;
                    }
                    else
                    {
                        DTI = Math.Min(TLIMIT2 * 1.5, 3600);
                    }
                }
                /*
                else if (windspeed < 1.0)
                {
                    DTI = TLIMIT2 * 1.5;
                }
                */
                else
                {
                    DTI = TLIMIT2 * 1.0;
                }
            }
            if (stability == 7)
            {
                if (windspeed < 0.5F)
                {
                    //in case of dynamic sun rise simulation
                    if (dynamic > 0)
                    {
                        DTI = Math.Min(TLIMIT2 * 2.0, 7200);
                    }
                    else
                    {
                        DTI = Math.Min(TLIMIT2 * 2.0, 3600);
                    }
                }
                else if (windspeed < 1.0)
                {
                    //in case of dynamic sun rise simulation
                    if (dynamic > 0)
                    {
                        DTI = Math.Min(TLIMIT2 * 1.5, 7200);
                    }
                    else
                    {
                        DTI = Math.Min(TLIMIT2 * 1.5, 3600);
                    }
                }
                else
                {
                    DTI = Math.Min(TLIMIT2 * 1.0, 3600);
                }
            }

            if (stability == 4)
            {
                DTI = Math.Min(TLIMIT2 * 0.17, 600);
            }
        }
        /// <summary>
        /// Create meteopgt.all and mettimeseries.dat for usage within ERA5 applications (i.e. postprocessing tools of the GUI)
        /// </summary>
        public static void GRALGenerateMeteoFilesForERA5()
        {
            System.Globalization.CultureInfo ic = System.Globalization.CultureInfo.InvariantCulture;
            //generate the gral input files mettimeseries.dat and meteopgt.all
            try
            {
                DateTime dateref = Program.dateUTC;
                dateref = dateref.AddSeconds(Program.REALTIME);

                //at the beginning of the ERA5 simulations the default files generated by the GUI need to be deleted, such that new files with correct time information can be established
                //the meteo information is arbitrary and due to the generated wnd and scl files not used by any subsequent GRAL simulation
                if (Program.IWETTER == 1)
                {
                    File.Delete("mettimeseries.dat");
                    File.Delete("meteopgt.all");
                    using (StreamWriter myWriterAll = new StreamWriter("meteopgt.all", true))
                    {
                        myWriterAll.WriteLine("10,0,0,");
                        myWriterAll.WriteLine("Wind direction sector,Wind speed class,stability class, frequency");
                    }
                }
                using (StreamWriter myWriterMet = new StreamWriter("mettimeseries.dat", true))
                {
                    using (StreamWriter myWriterAll = new StreamWriter("meteopgt.all", true))
                    {
                        string vel = Math.Round(((int)(Program.IWETTER / 360) * 0.1) + 0.1, 2).ToString(ic);
                        int dirDeg = (Program.IWETTER % 360);
                        string dirDecDeg = Math.Round(dirDeg * 0.1, 1).ToString(ic);
                        string _dateTime = dateref.Day.ToString(ic) + "." + dateref.Month.ToString(ic) + "," + dateref.Hour.ToString(ic) + ",";
                        myWriterMet.WriteLine(_dateTime + vel + "," + dirDecDeg + ",4");
                        myWriterAll.WriteLine(dirDecDeg + "," + vel + ",4" + ",0.05556");
                    }
                }
            }
            catch
            {

            }
        }
    }
}
