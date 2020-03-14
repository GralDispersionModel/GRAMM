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
using System.Threading.Tasks;

namespace GRAMM_2001
{
    partial class Program
    {
        public static void Bords_calculate(int NI, int NJ, int NK)
        {
            //damping factor
            double ATTENU = 2;
            double ATTENU2 = 0.1; //stronger damping factor
            double ATTENU1 = 5.0; //weaker damping factor
            double WINDDIR1 = 0;
            double WINDDIR2 = 0;
            double WINDGE1 = 0;
            double WINDGE2 = 0;

            //transient boundary conditions using ECMWF ERA 5 GRIB2 DATA
            if (REALTIME > REINITIALIZATION_Threshold && (Program.ISTAT == 2 || Program.ISTAT == 4))
            {
                REINITIALIZATION_Threshold += REINITIALIZATION;
                ERA5BOUND_Threshold += ERA5BOUND;

                if (Program.ISTAT == 2 || Program.ISTAT == 4)
                {
                    DateTime dateref = Program.dateUTC;
                    dateref = dateref.AddSeconds(REALTIME);
                    List<float> XERA5 = new List<float>();
                    List<float> YERA5 = new List<float>();
                    List<float> ZERA5 = new List<float>();
                    List<float> OROERA5 = new List<float>();
                    List<float> UERA5 = new List<float>();
                    List<float> VERA5 = new List<float>();
                    List<float> WERA5 = new List<float>();
                    List<float> TERA5 = new List<float>();
                    List<float> QERA5 = new List<float>();
                    List<float> WATERCLOUDERA5 = new List<float>();
                    List<float> SDERA5 = new List<float>();
                    List<float> CCERA5 = new List<float>();
                    List<float> LCCERA5 = new List<float>();
                    List<float> MCCERA5 = new List<float>();
                    List<float> HCCERA5 = new List<float>();
                    List<float> PERA5 = new List<float>();
                    List<float> STERA5 = new List<float>();
                    List<float> SWERA5 = new List<float>();
                    List<float> MSLPERA5 = new List<float>();
                    List<float> SEATEMPERA5 = new List<float>();

                    //read ERA5 grib2 files
                    ERA5_Read(dateref, ref ERA5_date1, ref ERA5_date2, ref XERA5, ref YERA5, ref ZERA5, ref OROERA5, ref UERA5, ref VERA5, ref WERA5, ref TERA5,
                        ref QERA5, ref SDERA5, ref CCERA5, ref PERA5, ref STERA5, ref SWERA5, ref MSLPERA5, ref LCCERA5, ref MCCERA5, ref HCCERA5, ref WATERCLOUDERA5, ref SEATEMPERA5);

                    //re-initialization using ERA5 data
                    ERA5_InitializeOntoGRAMMgrid(XERA5, YERA5, ZERA5, UERA5, VERA5, WERA5, TERA5, QERA5, PERA5, SDERA5, CCERA5, STERA5, SWERA5, MSLPERA5, LCCERA5, MCCERA5, HCCERA5, WATERCLOUDERA5, SEATEMPERA5);

                    if (ICSTR)
                    {
                        dateref = new DateTime(Program.IJAHR4digits, Program.IMON, Program.ITAG, Program.ISTU, Program.IMIN, 0);
                        dateref = dateref.AddSeconds(REALTIME);
                        int IMON_RAD = dateref.Month;
                        int ITAG_RAD = dateref.Day;
                        int ISTUD = dateref.Hour;
                        int IMIND = dateref.Minute;
                        double TJETZT1 = ISTUD * 3600;
                        double AMIND = (TJETZT / 3600 - (float)ISTUD) * 60;
                        IMIND = Convert.ToInt32(Math.Floor(AMIND));
                        double ASECD = (AMIND - (float)IMIND) * 60;
                        int ISECD = Convert.ToInt32(Math.Floor(ASECD));
                        RadiationModel.RADIATRAD(false, false, false, ITAG_RAD, IMON_RAD, IJAHR, TJETZT1, NX, NY, NZ);
                        Console.WriteLine(ITAG_RAD.ToString() + "." + IMON_RAD.ToString() + "  -  " + ISTUD.ToString() + ":" + IMIND.ToString("D2"));
                    }
                }
            }
            else if (REALTIME > ERA5BOUND_Threshold && (Program.ISTAT == 2 || Program.ISTAT == 4))
            {
                ERA5BOUND_Threshold += ERA5BOUND;

                if (Program.ISTAT == 2 || Program.ISTAT == 4)
                {
                    DateTime dateref = Program.dateUTC;
                    dateref = dateref.AddSeconds(REALTIME);
                    List<float> XERA5 = new List<float>();
                    List<float> YERA5 = new List<float>();
                    List<float> ZERA5 = new List<float>();
                    List<float> OROERA5 = new List<float>();
                    List<float> UERA5 = new List<float>();
                    List<float> VERA5 = new List<float>();
                    List<float> WERA5 = new List<float>();
                    List<float> TERA5 = new List<float>();
                    List<float> QERA5 = new List<float>();
                    List<float> WATERCLOUDERA5 = new List<float>();
                    List<float> SDERA5 = new List<float>();
                    List<float> CCERA5 = new List<float>();
                    List<float> LCCERA5 = new List<float>();
                    List<float> MCCERA5 = new List<float>();
                    List<float> HCCERA5 = new List<float>();
                    List<float> PERA5 = new List<float>();
                    List<float> STERA5 = new List<float>();
                    List<float> SWERA5 = new List<float>();
                    List<float> MSLPERA5 = new List<float>();
                    List<float> SEATEMPERA5 = new List<float>();

                    //read ERA5 grib2 files
                    ERA5_Read(dateref, ref ERA5_date1, ref ERA5_date2, ref XERA5, ref YERA5, ref ZERA5, ref OROERA5, ref UERA5, ref VERA5, ref WERA5, ref TERA5,
                        ref QERA5, ref SDERA5, ref CCERA5, ref PERA5, ref STERA5, ref SWERA5, ref MSLPERA5, ref LCCERA5, ref MCCERA5, ref HCCERA5, ref WATERCLOUDERA5, ref SEATEMPERA5);

                    //calculate boundary conditions using ERA5 data
                    ERA5_BoundaryConditions2(XERA5, YERA5, ZERA5, UERA5, VERA5, WERA5, TERA5, QERA5, PERA5, SDERA5, CCERA5, STERA5, SWERA5, MSLPERA5, LCCERA5, MCCERA5, HCCERA5, WATERCLOUDERA5, SEATEMPERA5);
                }
            }

            //in case of non-steady-state option: read file meteopgt.all to get geostrophic wind for forcing
            Int16 IHOUR = Convert.ToInt16(Math.Floor(TJETZT / IOUTPUT));
            if (((REALTIME == 0) || (REALTIME > 0) && ((ITIME % IRAD) == 0)) && (Program.ISTAT == 1))
            {
                bool meteopgtexist = File.Exists("meteopgt.all");
                if (meteopgtexist == true)
                {
                    int inid = 1;
                    try
                    {
                        using (FileStream fs = new FileStream("mettimeseries.dat", FileMode.Open, FileAccess.Read, FileShare.Read))
                        {
                            using (StreamReader myreader = new StreamReader(fs))
                            {
                                string[] text = new string[10];
                                string[] datum = new string[2];
                                for (; inid < Program.IWETTER; inid++)
                                {
                                    text = myreader.ReadLine().Split(new char[] { ' ', ';', ',', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                                }
                                text = myreader.ReadLine().Split(new char[] { ' ', ';', ',', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                                WINDDIR1 = Convert.ToDouble(text[3].Replace(".", Program.decsep));
                                WINDGE1 = Math.Max(Convert.ToDouble(text[2].Replace(".", Program.decsep)), 0.001);

                                datum = text[0].Split(new char[] { '.', ':' }, StringSplitOptions.RemoveEmptyEntries);
                                ITAG = Convert.ToInt32(datum[0]);
                                IMON = Convert.ToInt32(datum[1]);
                                ISTU = Convert.ToInt32(text[1]);
                                //Console.WriteLine(ITAG.ToString() + "," + IMON.ToString() + "," + ISTU.ToString());

                                text = myreader.ReadLine().Split(new char[] { ' ', ';', ',', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                                WINDDIR2 = Convert.ToDouble(text[3].Replace(".", Program.decsep));
                                WINDGE2 = Math.Max(Convert.ToDouble(text[2].Replace(".", Program.decsep)), 0.001);

                                WINDDIR1 *= 10;
                                WINDDIR2 *= 10;
                                WINDDIR1 = (270 - WINDDIR1) * Math.PI / 180;
                                WINDDIR2 = (270 - WINDDIR2) * Math.PI / 180;

                                WU1 = WINDGE1 * Math.Cos(WINDDIR1);
                                WV1 = WINDGE1 * Math.Sin(WINDDIR1);
                                WU2 = WINDGE2 * Math.Cos(WINDDIR1);
                                WV2 = WINDGE2 * Math.Sin(WINDDIR1);
                            }
                        }
                    }
                    catch
                    {
                        Console.WriteLine("Error when reading file mettimeseries.dat in line " + inid.ToString() + "- Execution stopped");
                        Environment.Exit(0);
                    }
                }
            }

            //transient forcing using mettimeseries.dat
            if (Program.ISTAT == 1)
            {
                double ANEMO_I = 1.0 / 10000.0;
                double TLIMIT_L = 1.0 / Program.TLIMIT2 * (Program.REALTIME - Program.TLIMIT + Program.TLIMIT2);
                Parallel.For(1, NI + 1, Program.pOptions, i =>
                  {
                      for (int j = 1; j <= NJ; j++)
                      {
                          float[] ZSP_L = Program.ZSP[i][j];
                          double[] UG_L = Program.UG[i][j];
                          double[] VG_L = Program.VG[i][j];
                          double AH = Program.AH[i][j];
                          double Delta_L = 0;

                          for (int k = 1; k <= NK; k++)
                          {
                              Delta_L = Math.Pow(Math.Min(ZSP_L[k] - AH, 10000) * ANEMO_I, 0.25);

                              UG_L[k] = WU1 - (WU1 - WU2) * TLIMIT_L;
                              VG_L[k] = WV1 - (WV1 - WV2) * TLIMIT_L;


                              Program.USWN[j][k] = UG_L[k] * Delta_L;
                              Program.VSWN[j][k] = VG_L[k] * Delta_L;
                              Program.USEN[j][k] = Program.USWN[j][k];
                              Program.VSEN[j][k] = Program.VSWN[j][k];
                              Program.USSN[i][k] = UG_L[k] * Delta_L;
                              Program.VSSN[i][k] = VG_L[k] * Delta_L;
                              Program.USNN[i][k] = Program.USSN[i][k];
                              Program.VSNN[i][k] = Program.VSSN[i][k];

                          }
                      }
                  });
            }

            //compute boundary values at the western border
            Parallel.For(1, NJ + 1, Program.pOptions, j =>
            {
                for (int k = 1; k <= NK; k++)
                {
                    Program.U1N[1][j][k] = Program.U1N[2][j][k];
                    Program.V1N[1][j][k] = Program.V1N[2][j][k];
                    Program.UN[1][j][k] = Program.U1N[2][j][k];
                    Program.VN[1][j][k] = Program.V1N[2][j][k];

                    Program.WN[1][j][k] = 0;
                    Program.TBN[1][j][2] = Program.TBN[2][j][2];
                    Program.QUN[1][j][k] = Program.QUN[2][j][k];

                    if (Program.IBOUND == 1 && Program.ISTAT <= 1)
                    {
                        Program.UN[1][j][k] = Program.U1N[2][j][k];
                        Program.TN[1][j][k] = Program.TN[2][j][k];
                    }
                    else if (Program.IBOUND == 6 && Program.ISTAT <= 1)
                    {
                        for (int i = 1; i <= 6; i++)
                        {
                            int i1 = Math.Min(i, NI);
                            double AALPHA = Math.Pow(2.73, -ATTENU * (float)i1);
                            Program.UN[i1][j][k] = Program.U1N[i1][j][k] - AALPHA * (Program.U1N[i1][j][k] - Program.USWN[j][k]);
                            Program.VN[i1][j][k] = Program.V1N[i1][j][k] - AALPHA * (Program.V1N[i1][j][k] - Program.VSWN[j][k]);
                        }
                        Program.WN[1][j][k] = 0;
                        Program.TN[1][j][k] = Program.TN[2][j][k];
                    }
                    else if (Program.IBOUND == 6 && Program.ISTAT >= 2)
                    {

                        for (int i = 1; i <= 1; i++)
                        {
                            int i1 = Math.Min(i, NI);
                            Linear_interpolation(1, NZ, k, 0.5, 2.0, ref ATTENU);
                            double AALPHA = Math.Pow(2.73, -ATTENU * (float)i1);
                            Program.UN[i1][j][k] = Program.UN[i1 + 1][j][k] - AALPHA * (Program.UN[i1 + 1][j][k] - Program.USWN[j][k]);
                            Program.VN[i1][j][k] = Program.VN[i1 + 1][j][k] - AALPHA * (Program.VN[i1 + 1][j][k] - Program.VSWN[j][k]);
                            Program.TN[i1][j][k] = Program.TN[i1 + 1][j][k] - AALPHA * (Program.TN[i1 + 1][j][k] - Program.TSWN[j][k]);
                        }

                        Program.WN[1][j][k] = 0;
                        //Program.VN[1][j][k] = Program.VSWN[j][k];
                        //Program.TN[1][j][k] = Program.TSWN[j][k];
                        Program.QUN[1][j][k] = Program.QUSWN[j][k];
                    }

                    Program.U1N[1][j][k] = Program.UN[1][j][k];
                    Program.U2N[1][j][k] = Program.UN[1][j][k];
                    Program.V1N[1][j][k] = Program.VN[1][j][k];
                    Program.V2N[1][j][k] = Program.VN[1][j][k];
                    Program.W1N[1][j][k] = Program.WN[1][j][k];
                    Program.W2N[1][j][k] = Program.WN[1][j][k];
                }
            });

            //compute boundary values at the eastern border
            Parallel.For(1, NJ + 1, Program.pOptions, j =>
            {
                for (int k = 1; k <= NK; k++)
                {

                    Program.U2N[NI][j][k] = Program.U2N[NI - 1][j][k];
                    Program.V2N[NI][j][k] = Program.V2N[NI - 1][j][k];
                    Program.UN[NI][j][k] = Program.U2N[NI - 1][j][k];
                    Program.VN[NI][j][k] = Program.V2N[NI - 1][j][k];


                    Program.WN[NI][j][k] = 0;
                    Program.TBN[NI][j][2] = Program.TBN[NI - 1][j][2];
                    Program.QUN[NI][j][k] = Program.QUN[NI - 1][j][k];

                    if (Program.IBOUND == 1 && Program.ISTAT <= 1)
                    {
                        Program.UN[NI][j][k] = Program.U2N[NI - 1][j][k];
                        Program.TN[NI][j][k] = Program.TN[NI - 1][j][k];
                    }
                    else if (Program.IBOUND == 6 && Program.ISTAT <= 1)
                    {
                        for (int i = NI - 5; i <= NI; i++)
                        {
                            int i1 = Math.Max(i, 1);
                            double AALPHA = Math.Pow(2.73, -ATTENU * (float)(NI - i1 + 1));
                            Program.UN[i1][j][k] = Program.U2N[i1][j][k] - AALPHA * (Program.U2N[i1][j][k] - Program.USEN[j][k]);
                            Program.VN[i1][j][k] = Program.V2N[i1][j][k] - AALPHA * (Program.V2N[i1][j][k] - Program.VSEN[j][k]);
                        }
                        Program.WN[NI][j][k] = 0;
                        Program.TN[NI][j][k] = Program.TN[NI - 1][j][k];
                    }
                    else if (Program.IBOUND == 6 && Program.ISTAT >= 2)
                    {

                        for (int i = NI; i <= NI; i++)
                        {
                            int i1 = Math.Max(i, 1);
                            Linear_interpolation(1, NZ, k, 0.5, 2.0, ref ATTENU);
                            double AALPHA = Math.Pow(2.73, -ATTENU * (float)(NI - i1 + 1));
                            Program.UN[i1][j][k] = Program.UN[i1 - 1][j][k] - AALPHA * (Program.UN[i1 - 1][j][k] - Program.USEN[j][k]);
                            Program.VN[i1][j][k] = Program.VN[i1 - 1][j][k] - AALPHA * (Program.VN[i1 - 1][j][k] - Program.VSEN[j][k]);
                            Program.TN[i1][j][k] = Program.TN[i1 - 1][j][k] - AALPHA * (Program.TN[i1 - 1][j][k] - Program.TSEN[j][k]);
                        }

                        Program.WN[NI][j][k] = 0;
                        //Program.VN[NI][j][k] = Program.VSEN[j][k];
                        //Program.TN[NI][j][k] = Program.TSEN[j][k];
                        Program.QUN[NI][j][k] = Program.QUSEN[j][k];
                    }

                    Program.U1N[NI][j][k] = Program.UN[NI][j][k];
                    Program.U2N[NI][j][k] = Program.UN[NI][j][k];
                    Program.V1N[NI][j][k] = Program.VN[NI][j][k];
                    Program.V2N[NI][j][k] = Program.VN[NI][j][k];
                    Program.W1N[NI][j][k] = Program.WN[NI][j][k];
                    Program.W2N[NI][j][k] = Program.WN[NI][j][k];


                }
            });

            //compute boundary values at the southern border
            Parallel.For(1, NI + 1, Program.pOptions, i =>
            {
                for (int k = 1; k <= NK; k++)
                {
                    Program.U1N[i][1][k] = Program.U1N[i][2][k];
                    Program.V1N[i][1][k] = Program.V1N[i][2][k];
                    Program.UN[i][1][k] = Program.U1N[i][2][k];
                    Program.VN[i][1][k] = Program.V1N[i][2][k];
                    Program.WN[i][1][k] = 0;
                    Program.TBN[i][1][2] = Program.TBN[i][2][2];
                    Program.QUN[i][1][k] = Program.QUN[i][2][k];

                    if (Program.IBOUND == 1 && Program.ISTAT <= 1)
                    {
                        Program.VN[i][1][k] = Program.V1N[i][2][k];
                        Program.TN[i][1][k] = Program.TN[i][2][k];
                    }
                    else if (Program.IBOUND == 6 && Program.ISTAT <= 1)
                    {
                        for (int j = 1; j <= 6; j++)
                        {
                            int j1 = Math.Min(j, NJ);
                            double AALPHA = Math.Pow(2.73, -ATTENU * (float)j1);
                            Program.VN[i][j1][k] = Program.V1N[i][j1][k] - AALPHA * (Program.V1N[i][j1][k] - Program.VSSN[i][k]);
                            Program.UN[i][j1][k] = Program.U1N[i][j1][k] - AALPHA * (Program.U1N[i][j1][k] - Program.USSN[i][k]);
                        }
                        Program.WN[i][1][k] = 0;
                        Program.TN[i][1][k] = Program.TN[i][2][k];
                    }
                    else if (Program.IBOUND == 6 && Program.ISTAT >= 2)
                    {

                        for (int j = 1; j <= 1; j++)
                        {
                            int j1 = Math.Min(j, NJ);
                            Linear_interpolation(1, NZ, k, 0.5, 2.0, ref ATTENU);
                            double AALPHA = Math.Pow(2.73, -ATTENU * (float)j1);
                            Program.VN[i][j1][k] = Program.VN[i][j1 + 1][k] - AALPHA * (Program.VN[i][j1 + 1][k] - Program.VSSN[i][k]);
                            Program.UN[i][j1][k] = Program.UN[i][j1 + 1][k] - AALPHA * (Program.UN[i][j1 + 1][k] - Program.USSN[i][k]);
                            Program.TN[i][j1][k] = Program.TN[i][j1 + 1][k] - AALPHA * (Program.TN[i][j1 + 1][k] - Program.TSSN[i][k]);
                        }

                        Program.WN[i][1][k] = 0;
                        //Program.UN[i][1][k] = Program.USSN[i][k];
                        //Program.TN[i][1][k] = Program.TSSN[i][k];
                        Program.QUN[i][1][k] = Program.QUSSN[i][k];
                    }

                    Program.U1N[i][1][k] = Program.UN[i][1][k];
                    Program.U2N[i][1][k] = Program.UN[i][1][k];
                    Program.V1N[i][1][k] = Program.VN[i][1][k];
                    Program.V2N[i][1][k] = Program.VN[i][1][k];
                    Program.W1N[i][1][k] = Program.WN[i][1][k];
                    Program.W2N[i][1][k] = Program.WN[i][1][k];

                }
            });

            //compute boundary values at the northern border
            Parallel.For(1, NI + 1, Program.pOptions, i =>
            {
                for (int k = 1; k <= NK; k++)
                {
                    Program.U2N[i][NJ][k] = Program.U2N[i][NJ - 1][k];
                    Program.V2N[i][NJ][k] = Program.V2N[i][NJ - 1][k];
                    Program.UN[i][NJ][k] = Program.U2N[i][NJ - 1][k];
                    Program.VN[i][NJ][k] = Program.V2N[i][NJ - 1][k];
                    Program.WN[i][NJ][k] = 0;
                    Program.TBN[i][NJ][2] = Program.TBN[i][NJ - 1][2];
                    Program.QUN[i][NJ][k] = Program.QUN[i][NJ - 1][k];

                    if (Program.IBOUND == 1 && Program.ISTAT <= 1)
                    {
                        Program.VN[i][NJ][k] = Program.V2N[i][NJ - 1][k];
                        Program.TN[i][NJ][k] = Program.TN[i][NJ - 1][k];
                    }
                    else if (Program.IBOUND == 6 && Program.ISTAT <= 1)
                    {
                        for (int j = NJ - 5; j <= NJ; j++)
                        {
                            int j1 = Math.Max(j, 1);
                            double AALPHA = Math.Pow(2.73, -ATTENU * (float)(NJ - j1 + 1));
                            Program.VN[i][j1][k] = Program.V2N[i][j1][k] - AALPHA * (Program.V2N[i][j1][k] - Program.VSNN[i][k]);
                            Program.UN[i][j1][k] = Program.U2N[i][j1][k] - AALPHA * (Program.U2N[i][j1][k] - Program.USNN[i][k]);
                        }
                        Program.WN[i][NJ][k] = 0;
                        Program.TN[i][NJ][k] = Program.TN[i][NJ - 1][k];
                    }
                    else if (Program.IBOUND == 6 && Program.ISTAT >= 2)
                    {

                        for (int j = NJ; j <= NJ; j++)
                        {
                            int j1 = Math.Max(j, 1);
                            Linear_interpolation(1, NZ, k, 0.5, 2.0, ref ATTENU);
                            double AALPHA = Math.Pow(2.73, -ATTENU * (float)(NJ - j1 + 1));
                            Program.VN[i][j1][k] = Program.VN[i][j1 - 1][k] - AALPHA * (Program.VN[i][j1 - 1][k] - Program.VSNN[i][k]);
                            Program.UN[i][j1][k] = Program.UN[i][j1 - 1][k] - AALPHA * (Program.UN[i][j1 - 1][k] - Program.USNN[i][k]);
                            Program.TN[i][j1][k] = Program.TN[i][j1 - 1][k] - AALPHA * (Program.TN[i][j1 - 1][k] - Program.TSNN[i][k]);
                        }

                        Program.WN[i][NJ][k] = 0;
                        //Program.UN[i][NJ][k] = Program.USNN[i][k];
                        //Program.TN[i][NJ][k] = Program.TSNN[i][k];
                        Program.QUN[i][NJ][k] = Program.QUSNN[i][k];
                    }

                    Program.U1N[i][NJ][k] = Program.UN[i][NJ][k];
                    Program.U2N[i][NJ][k] = Program.UN[i][NJ][k];
                    Program.V1N[i][NJ][k] = Program.VN[i][NJ][k];
                    Program.V2N[i][NJ][k] = Program.VN[i][NJ][k];
                    Program.W1N[i][NJ][k] = Program.WN[i][NJ][k];
                    Program.W2N[i][NJ][k] = Program.WN[i][NJ][k];

                }
            });

            Parallel.For(1, NJ + 1, Program.pOptions, j =>
            {
                for (int k = 1; k <= NK; k++)
                {
                    Program.T[1][j][k] = Program.TN[1][j][k];
                    Program.QU[1][j][k] = Program.QUN[1][j][k];

                    Program.T[NI][j][k] = Program.TN[NI][j][k];
                    Program.QU[NI][j][k] = Program.QUN[NI][j][k];
                }
            });


            Parallel.For(1, NI + 1, Program.pOptions, i =>
            {
                for (int k = 1; k <= NK; k++)
                {
                    Program.T[i][1][k] = Program.TN[i][1][k];
                    Program.QU[i][1][k] = Program.QUN[i][1][k];

                    Program.T[i][NJ][k] = Program.TN[i][NJ][k];
                    Program.QU[i][NJ][k] = Program.QUN[i][NJ][k];
                }
            });

            //Corner values
            for (int k = 1; k <= NK; k++)
            {
                //Corner points 1,1
                Program.UN[1][1][k] = Program.UN[2][2][k];
                Program.VN[1][1][k] = Program.VN[2][2][k];
                Program.WN[1][1][k] = 0;
                Program.TN[1][1][k] = Program.TN[2][2][k];
                Program.QUN[1][1][k] = Program.QUN[2][2][k];

                //Corner points 1,NJ
                Program.UN[1][NJ][k] = Program.UN[2][NJ - 1][k];
                Program.VN[1][NJ][k] = Program.VN[2][NJ - 1][k];
                Program.WN[1][NJ][k] = 0;
                Program.TN[1][NJ][k] = Program.TN[2][NJ - 1][k];
                Program.QUN[1][NJ][k] = Program.QUN[2][NJ - 1][k];

                //Corner points NI,1
                Program.UN[NI][1][k] = Program.UN[NI - 1][2][k];
                Program.VN[NI][1][k] = Program.VN[NI - 1][2][k];
                Program.WN[NI][1][k] = 0;
                Program.TN[NI][1][k] = Program.TN[NI - 1][2][k];
                Program.QUN[NI][1][k] = Program.QUN[NI - 1][2][k];

                //Corner points NI,NJ
                Program.UN[NI][NJ][k] = Program.UN[NI][NJ - 1][k];
                Program.UN[NI][NJ][k] = Program.UN[NI - 1][NJ][k];
                Program.VN[NI][NJ][k] = Program.VN[NI][NJ - 1][k];
                Program.VN[NI][NJ][k] = Program.VN[NI - 1][NJ][k];
                Program.WN[NI][NJ][k] = 0;
                Program.TN[NI][NJ][k] = Program.TN[NI - 1][NJ - 1][k];
                Program.QUN[NI][NJ][k] = Program.QUN[NI - 1][NJ - 1][k];
            }
        }

    }
}
