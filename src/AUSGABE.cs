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
using System.IO;
using System.IO.Compression;
using System.Linq;

namespace GRAMM_2001
{
    partial class Program
    {

        /// <summary>
        /// Write result files to the calculation folder
        /// </summary>
        /// <param name="NI"></param>
        /// <param name="NJ"></param>
        /// <param name="NK"></param>
        /// <param name="intermediate"></param>
        public static void OUTPUT(int NI, int NJ, int NK, bool intermediate)
        {
            //Check for quasi-steady-state requirement (VDI 3783-7) -> the lateral region where topography is smoothed is not considered
            if (Program.REALTIME > 0.8 * Program.DTI)
            {
                StreamWriter writesteadystate = null; // create streamwriter-container


                if (Program.WriteSteadyState) // write data for steady state criterion to a file
                {
                    string filename = (Convert.ToString(Program.IWETTER).PadLeft(5, '0') + "_steady_state.txt");
                    try
                    {
                        writesteadystate = new StreamWriter(filename);
                        writesteadystate.WriteLine("ncols         " + Convert.ToString(NI - 2 * Program.nr_cell_smooth));
                        writesteadystate.WriteLine("nrows         " + Convert.ToString(NJ - 2 * Program.nr_cell_smooth));
                        writesteadystate.WriteLine("xllcorner     " + Convert.ToString(Program.IKOOA + Program.nr_cell_smooth * Program.DDXImm[1]));
                        writesteadystate.WriteLine("yllcorner     " + Convert.ToString(Program.JKOOA + Program.nr_cell_smooth * Program.DDXImm[1]));
                        writesteadystate.WriteLine("cellsize      " + Convert.ToString(Program.DDXImm[1]));
                        writesteadystate.WriteLine("NODATA_value  " + "-9999");
                    }
                    catch { }
                }

                double UTREFFER = 0;
                double VTREFFER = 0;
                double WTREFFER = 0;
                for (int j = NJ - Program.nr_cell_smooth; j >= 1 + Program.nr_cell_smooth; j--)
                {
                    string steadystateline = "";
                    for (int i = 1 + Program.nr_cell_smooth; i <= NI - Program.nr_cell_smooth; i++)
                    {
                        Byte steadystatevalue = 0;
                        if ((Math.Abs(Program.U_TEMP[i][j] - Program.U[i][j][1]) <= 0.35) || (Math.Abs((Program.U_TEMP[i][j] - Program.U[i][j][1]) / (Math.Max(Math.Abs(Program.U[i][j][1]), 0.001))) <= 0.1))
                        {
                            UTREFFER++;
                            steadystatevalue = 1;
                        }
                        if ((Math.Abs(Program.V_TEMP[i][j] - Program.V[i][j][1]) <= 0.35) || (Math.Abs((Program.V_TEMP[i][j] - Program.V[i][j][1]) / (Math.Max(Math.Abs(Program.V[i][j][1]), 0.001))) <= 0.1))
                        {
                            VTREFFER++;
                            steadystatevalue += 2;
                        }
                        if ((Math.Abs(Program.W_TEMP[i][j] - Program.W[i][j][1]) <= 0.35) || (Math.Abs((Program.W_TEMP[i][j] - Program.W[i][j][1]) / (Math.Max(Math.Abs(Program.W[i][j][1]), 0.001))) <= 0.1))
                        {
                            WTREFFER++;
                            steadystatevalue += 4;
                        }

                        if (writesteadystate != null) // write file!
                        {
                            if (j == 1 + Program.nr_cell_smooth && i == 1 + Program.nr_cell_smooth) // Avoid blank raster data set for GUI
                            {
                                steadystateline += Convert.ToString((double)steadystatevalue + 0.01).Replace(Program.decsep, ".") + " ";
                            }
                            else
                            {
                                steadystateline += Convert.ToString(steadystatevalue) + " ";
                            }
                        }
                    }

                    if (writesteadystate != null)
                    {
                        writesteadystate.WriteLine(steadystateline);
                    }

                }

                if (writesteadystate != null)
                {
                    try
                    {
                        writesteadystate.Close();
                        writesteadystate.Dispose();
                    }
                    catch { }
                }



                UTREFFER /= ((NI - 2 * Program.nr_cell_smooth) * (NJ - 2 * Program.nr_cell_smooth));
                VTREFFER /= ((NI - 2 * Program.nr_cell_smooth) * (NJ - 2 * Program.nr_cell_smooth));
                WTREFFER /= ((NI - 2 * Program.nr_cell_smooth) * (NJ - 2 * Program.nr_cell_smooth));
            }

            //Check for possible numerical instabilities
            if (Program.MASSOURCE_FIRST < Program.SUMG)
            {
                //compute maximum wind speeds
                float UMAX = 0;
                float VMAX = 0;
                float WMAX = 0;
                for (int i = 2; i <= NI - 1; i++)
                {
                    for (int j = 2; j <= NJ - 1; j++)
                    {
                        for (int k = 1; k <= NK - 1; k++)
                        {
                            UMAX = (float)Math.Max(Math.Abs(Program.U[i][j][k]), UMAX);
                            VMAX = (float)Math.Max(Math.Abs(Program.V[i][j][k]), VMAX);
                            WMAX = (float)Math.Max(Math.Abs(Program.W[i][j][k]), WMAX);
                        }
                    }
                }

                try
                {
                    if ((UMAX > 50) || (VMAX > 50) || (WMAX > 20))
                    {
                        try
                        {
                            string err = "Final mass divergence is larger than initial mass divergence." + "  Possible numerical instabilities detected for flow field: " +
                                             Program.IWETTER.ToString();
                            ProgramWriters.LogfileProblemreportWrite(err);
                            err = "Maximum wind speeds in east/west/vertical direction: " + UMAX.ToString("0000.0") + ";" + VMAX.ToString("0000.0") + ";" + WMAX.ToString("0000.0");
                            ProgramWriters.LogfileProblemreportWrite(err);
                        }
                        catch { }
                    }
                }
                catch { }
            }

            //differentiating between intermediate output and final output in case of meteopgt.all input
            string wndfilename;
            if (intermediate == false)
            {
                wndfilename = (Convert.ToString(Program.IWETTER).PadLeft(5, '0') + ".wnd");
            }
            else
            {
                wndfilename = (Convert.ToString(Program.IWETTER).PadLeft(5, '0') + ".wnd");
                //intermediate output in case of meteopgt.all
                if (((METEO == "Y") || (METEO == "y")) && (Program.ISTAT == 0))
                {
                    //compute fictious weather number
                    int FILENUMBER = 0;
                    double time_real = Program.REALTIME;
                    if (time_real < Program.IOUTPUT) // 1st situation should be stored after 1800 s
                    {
                        time_real = Program.IOUTPUT + 2;
                    }
                    Meteopgtall.meteopgtall_calculate(Program.meteopgt_nr, Program.IWETTER, Program.DTI, time_real, Program.IOUTPUT, Program.TLIMIT2, ref FILENUMBER);
                    wndfilename = (Convert.ToString(FILENUMBER).PadLeft(5, '0') + ".wnd");
                }
            }

            Console.Write(wndfilename + "  "); // write windfile name to the console

            BinaryWriter writer = new BinaryWriter(File.Open(wndfilename, FileMode.Create));
            int header = -1;
            Int16 dummy;
            float GRAMMhorgridsize = (float)Program.DDXImm[1];

            //there are two different formats: IOUTPUT = 0 (standard output for GRAL-GUI users) and IOUTPUT = 3 for SOUNDPLAN USERS
            if (Program.IOUT == 0)
            {
                writer.Write(header);
                writer.Write(NI);
                writer.Write(NJ);
                writer.Write(NK);
                writer.Write(GRAMMhorgridsize);
                for (int i = 1; i <= NI; i++)
                {
                    for (int j = 1; j <= NJ; j++)
                    {
                        for (int k = 1; k <= NK; k++)
                        {
                            try
                            {
                                dummy = Convert.ToInt16(Math.Max(Int16.MinValue, Math.Min(Int16.MaxValue, Program.U[i][j][k] * 100)));
                            }
                            catch
                            {
                                //Console.WriteLine("U " + e.Message.ToString() + Program.U[i][j][k].ToString());
                                dummy = Int16.MaxValue;
                            }
                            writer.Write(dummy);

                            try
                            {
                                dummy = Convert.ToInt16(Math.Max(Int16.MinValue, Math.Min(Int16.MaxValue, Program.V[i][j][k] * 100)));
                            }
                            catch
                            {
                                //Console.WriteLine("V " + e.Message.ToString() + Program.V[i][j][k].ToString());
                                dummy = Int16.MaxValue;
                            }
                            writer.Write(dummy);
                            try
                            {
                                dummy = Convert.ToInt16(Math.Max(Int16.MinValue, Math.Min(Int16.MaxValue,Program.W[i][j][k] * 100)));
                            }
                            catch
                            {
                                //Console.WriteLine("W " + e.Message.ToString() + Program.W[i][j][k].ToString());
                                dummy = Int16.MaxValue;
                            }
                            writer.Write(dummy);
                        }
                    }
                }
            }
            if (Program.IOUT == 3)
            {
                for (int i = 1; i <= NI; i++)
                {
                    for (int j = 1; j <= NJ; j++)
                    {
                        for (int k = 1; k <= NK; k++)
                        {
                            writer.Write((float)Program.U[i][j][k]);
                            writer.Write((float)Program.V[i][j][k]);
                            writer.Write((float)Program.W[i][j][k]);
                        }
                    }
                }
            }

            writer.Close();
            writer.Dispose();

            //output for friction velocity and Obukhov length
            // write a Zip file
            string stabclassfilename = (Convert.ToString(Program.IWETTER).PadLeft(5, '0') + ".scl");

            if (intermediate == false)
            {
                stabclassfilename = (Convert.ToString(Program.IWETTER).PadLeft(5, '0') + ".scl");
            }
            else
            {
                stabclassfilename = (Convert.ToString(Program.IWETTER).PadLeft(5, '0') + ".scl");
                //intermediate output in case of meteopgt.all
                if (((METEO == "Y") || (METEO == "y")) && (Program.ISTAT == 0))
                {
                    //compute fictious weather number
                    int FILENUMBER = 0;
                    double time_real = Program.REALTIME;
                    if (time_real < Program.IOUTPUT) // 1st situation should be stored after 1800 s
                    {
                        time_real = Program.IOUTPUT + 2;
                    }
                    Meteopgtall.meteopgtall_calculate(Program.meteopgt_nr, Program.IWETTER, Program.DTI, time_real, Program.IOUTPUT, Program.TLIMIT2, ref FILENUMBER);
                    stabclassfilename = (Convert.ToString(FILENUMBER).PadLeft(5, '0') + ".scl");
                }
            }

            try
            {
                Console.Write(stabclassfilename + "  "); // write windfile name to the console

                using (FileStream zipToOpen = new FileStream(stabclassfilename, FileMode.Create))
                {
                    using (ZipArchive archive = new ZipArchive(zipToOpen, ZipArchiveMode.Update))
                    {
                        string ustarfilename = (Convert.ToString(Program.IWETTER).PadLeft(5, '0') + ".ust");
                        ZipArchiveEntry write_entry1 = archive.CreateEntry(ustarfilename);
                        using (writer = new BinaryWriter(write_entry1.Open()))
                        {
                            writer.Write(header);
                            writer.Write(NI);
                            writer.Write(NJ);
                            writer.Write(NK);
                            writer.Write(GRAMMhorgridsize);
                            for (int i = 1; i <= NI; i++)
                            {
                                for (int j = 1; j <= NJ; j++)
                                {
                                    writer.Write(Convert.ToInt16(Math.Max(Int16.MinValue, Math.Min(Int16.MaxValue, Program.UST[i][j] * 1000))));
                                }
                            }
                        }

                        string obukhovfilename = (Convert.ToString(Program.IWETTER).PadLeft(5, '0') + ".obl");
                        ZipArchiveEntry write_entry2 = archive.CreateEntry(obukhovfilename);
                        using (writer = new BinaryWriter(write_entry2.Open()))
                        {
                            writer.Write(header);
                            writer.Write(NI);
                            writer.Write(NJ);
                            writer.Write(NK);
                            writer.Write(GRAMMhorgridsize);
                            for (int i = 1; i <= NI; i++)
                            {
                                for (int j = 1; j <= NJ; j++)
                                {
                                    writer.Write(Convert.ToInt16(Math.Max(Int16.MinValue, Math.Min(Int16.MaxValue, Program.OL[i][j] * 10))));
                                }
                            }
                        }

                        //computation and ouput of stability classes
                        string stabilityfile = (Convert.ToString(Program.IWETTER).PadLeft(5, '0') + ".scl");
                        ZipArchiveEntry write_entry3 = archive.CreateEntry(stabilityfile);
                        using (writer = new BinaryWriter(write_entry3.Open()))
                        {
                            writer.Write(header);
                            writer.Write(NI);
                            writer.Write(NJ);
                            writer.Write(NK);
                            writer.Write(GRAMMhorgridsize);
                            for (int i = 1; i <= NI; i++)
                            {
                                for (int j = 1; j <= NJ; j++)
                                {
                                    writer.Write((Int16) Math.Max(Int16.MinValue, Math.Min(Int16.MaxValue, Program.stabilityclass[i][j])));
                                }
                            }
                        }
                    } // archive
                } // Zip File
            } // catch
            catch { }

            Program.Max_Proc_File_Read(); // read number of max. Processors

            if (recexist == true && IWetter_Console_First == 0) // write receptor wind fields; not if multi-instances are used
            {
                try
                {
                    using (StreamWriter wr = new StreamWriter("GRAMM.dat", true))
                    {
                        for (int ianz = 0; ianz < Urec.Count(); ianz++)
                        {
                            double angle = winkel(Urec[ianz], Vrec[ianz]);
                            wr.Write(Math.Sqrt(Math.Pow(Urec[ianz], 2) + Math.Pow(Vrec[ianz], 2)).ToString("0.00").Replace(decsep, ".") +
                                     "," + angle.ToString("0").Replace(decsep, ".") + "," + Trec[ianz].ToString("0.0").Replace(decsep, ".") + "," + Globradrec[ianz].ToString("0").Replace(decsep, ".")
                                     + "," + Longradrec[ianz].ToString("0").Replace(decsep, ".") + "," + Soilheatfluxrec[ianz].ToString("0").Replace(decsep, ".")
                                     + "," + Sensheatfluxrec[ianz].ToString("0").Replace(decsep, ".") + "," + Latheatfluxrec[ianz].ToString("0").Replace(decsep, ".") + ",");
                            Urec[ianz] = 0;
                            Vrec[ianz] = 0;
                            Trec[ianz] = 0;
                            Globradrec[ianz] = 0;
                            Longradrec[ianz] = 0;
                            Soilheatfluxrec[ianz] = 0;
                            Sensheatfluxrec[ianz] = 0;
                            Latheatfluxrec[ianz] = 0;
                        }
                        wr.WriteLine();
                    }
                }
                catch { }
            }

        }

    }
}
