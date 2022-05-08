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
using System.Collections.Immutable;
using System.Runtime.CompilerServices;
using System.Collections.Concurrent;

namespace GRAMM_2001
{
    partial class Program
    {
        /// <summary>
        /// Set the boundary condítions
        /// </summary>
        /// <param name="NI"></param>
        /// <param name="NJ"></param>
        /// <param name="NK"></param>
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        public static void Bords_calculate(int NI, int NJ, int NK)
        {
            //damping factor
            //double ATTENU2 = 0.1; //stronger damping factor
            //double ATTENU1 = 5.0; //weaker damping factor
            double WINDDIR1 = 0;
            double WINDDIR2 = 0;
            double WINDGE1 = 0;
            double WINDGE2 = 0;

            //in case of non-steady-state option: read file meteopgt.all to get geostrophic wind for forcing
            int IHOUR = Convert.ToInt32(Math.Floor(TJETZT / IOUTPUT));
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
                          ImmutableArray<float> ZSP_L = Program.ZSPImm[i][j];
                          double[] UG_L = Program.UG[i][j];
                          double[] VG_L = Program.VG[i][j];
                          double AH = Program.AHImm[i][j];
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
            //Parallel.For(1, NJ + 1, Program.pOptions, j =>
            Parallel.ForEach(Partitioner.Create(1, NJ + 1, NJ / Program.pOptions.MaxDegreeOfParallelism), range =>
            {
                double ATTENU = 2;
                for (int j = range.Item1; j < range.Item2; ++j)
                {
                    Array.Clear(Program.WN[1][j]);
                    FastCopy.CopyArraySourceLen(Program.U1N[2][j], Program.U1N[1][j]);
                    FastCopy.CopyArraySourceLen(Program.V1N[2][j], Program.V1N[1][j]);
                    FastCopy.CopyArraySourceLen(Program.U1N[2][j],  Program.UN[1][j]);
                    FastCopy.CopyArraySourceLen(Program.V1N[2][j],  Program.VN[1][j]);
                    FastCopy.CopyArraySourceLen(Program.QUN[2][j], Program.QUN[1][j]);
                    Program.TBN[1][j][2] = Program.TBN[2][j][2];

                    for (int k = 1; k <= NK; k++)
                    {
                        if (Program.IBOUND == 1 && Program.ISTAT <= 1)
                        {
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
                                ATTENU = Linear_interpolation(1, NZ, k, 0.5, 2.0);
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
                    }
                    FastCopy.CopyArraySourceLen(Program.UN[1][j], Program.U1N[1][j]);
                    FastCopy.CopyArraySourceLen(Program.UN[1][j], Program.U2N[1][j]);
                    FastCopy.CopyArraySourceLen(Program.VN[1][j], Program.V1N[1][j]);
                    FastCopy.CopyArraySourceLen(Program.VN[1][j], Program.V2N[1][j]);
                    FastCopy.CopyArraySourceLen(Program.WN[1][j], Program.W1N[1][j]);
                    FastCopy.CopyArraySourceLen(Program.WN[1][j], Program.W2N[1][j]);
                }
            });

            //compute boundary values at the eastern border
            //Parallel.For(1, NJ + 1, Program.pOptions, j =>
            Parallel.ForEach(Partitioner.Create(1, NJ + 1, NJ / Program.pOptions.MaxDegreeOfParallelism), range =>
            {
                double ATTENU = 2;
                for (int j = range.Item1; j < range.Item2; ++j)
                {
                    Array.Clear(Program.WN[NI][j]);
                    FastCopy.CopyArraySourceLen(Program.U2N[NI - 1][j], Program.U2N[NI][j]);
                    FastCopy.CopyArraySourceLen(Program.V2N[NI - 1][j], Program.V2N[NI][j]);
                    FastCopy.CopyArraySourceLen(Program.U2N[NI - 1][j],  Program.UN[NI][j]);
                    FastCopy.CopyArraySourceLen(Program.V2N[NI - 1][j],  Program.VN[NI][j]);
                    FastCopy.CopyArraySourceLen(Program.QUN[NI - 1][j], Program.QUN[NI][j]);
                    Program.TBN[NI][j][2] = Program.TBN[NI - 1][j][2];

                    for (int k = 1; k <= NK; k++)
                    {
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
                                ATTENU = Linear_interpolation(1, NZ, k, 0.5, 2.0);
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
                    }
                    FastCopy.CopyArraySourceLen(Program.UN[NI][j], Program.U1N[NI][j]);
                    FastCopy.CopyArraySourceLen(Program.UN[NI][j], Program.U2N[NI][j]);
                    FastCopy.CopyArraySourceLen(Program.VN[NI][j], Program.V1N[NI][j]);
                    FastCopy.CopyArraySourceLen(Program.VN[NI][j], Program.V2N[NI][j]);
                    FastCopy.CopyArraySourceLen(Program.WN[NI][j], Program.W1N[NI][j]);
                    FastCopy.CopyArraySourceLen(Program.WN[NI][j], Program.W2N[NI][j]);
                }
            });

            //compute boundary values at the southern border
            //Parallel.For(1, NI + 1, Program.pOptions, i =>
            Parallel.ForEach(Partitioner.Create(1, NI + 1, NI / Program.pOptions.MaxDegreeOfParallelism), range =>
            {
                double ATTENU = 2;
                for (int i = range.Item1; i < range.Item2; ++i)
                {
                    Array.Clear(Program.WN[i][1]);
                    FastCopy.CopyArraySourceLen(Program.U1N[i][2], Program.U1N[i][1]);
                    FastCopy.CopyArraySourceLen(Program.V1N[i][2], Program.V1N[i][1]);
                    FastCopy.CopyArraySourceLen(Program.U1N[i][2],  Program.UN[i][1]);
                    FastCopy.CopyArraySourceLen(Program.V1N[i][2],  Program.VN[i][1]);
                    FastCopy.CopyArraySourceLen(Program.QUN[i][2], Program.QUN[i][1]);
                    Program.TBN[i][1][2] = Program.TBN[i][2][2];

                    for (int k = 1; k <= NK; k++)
                    {
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
                                ATTENU = Linear_interpolation(1, NZ, k, 0.5, 2.0);
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
                    }
                    FastCopy.CopyArraySourceLen(Program.UN[i][1], Program.U1N[i][1]);
                    FastCopy.CopyArraySourceLen(Program.UN[i][1], Program.U2N[i][1]);
                    FastCopy.CopyArraySourceLen(Program.VN[i][1], Program.V1N[i][1]);
                    FastCopy.CopyArraySourceLen(Program.VN[i][1], Program.V2N[i][1]);
                    FastCopy.CopyArraySourceLen(Program.WN[i][1], Program.W1N[i][1]);
                    FastCopy.CopyArraySourceLen(Program.WN[i][1], Program.W2N[i][1]);
                }
            });

            //compute boundary values at the northern border
            //Parallel.For(1, NI + 1, Program.pOptions, i =>
            Parallel.ForEach(Partitioner.Create(1, NI + 1, NI / Program.pOptions.MaxDegreeOfParallelism), range =>
            {
                double ATTENU = 2;
                for (int i = range.Item1; i < range.Item2; ++i)
                {
                    Array.Clear(Program.WN[i][NJ]);
                    FastCopy.CopyArraySourceLen(Program.U2N[i][NJ - 1], Program.U2N[i][NJ]);
                    FastCopy.CopyArraySourceLen(Program.V2N[i][NJ - 1], Program.V2N[i][NJ]);
                    FastCopy.CopyArraySourceLen(Program.U2N[i][NJ - 1],  Program.UN[i][NJ]);
                    FastCopy.CopyArraySourceLen(Program.V2N[i][NJ - 1],  Program.VN[i][NJ]);
                    FastCopy.CopyArraySourceLen(Program.QUN[i][NJ - 1], Program.QUN[i][NJ]);
                    Program.TBN[i][NJ][2] = Program.TBN[i][NJ - 1][2];

                    for (int k = 1; k <= NK; k++)
                    {
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
                                ATTENU = Linear_interpolation(1, NZ, k, 0.5, 2.0);
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
                    }
                    FastCopy.CopyArraySourceLen(Program.UN[i][NJ], Program.U1N[i][NJ]);
                    FastCopy.CopyArraySourceLen(Program.UN[i][NJ], Program.U2N[i][NJ]);
                    FastCopy.CopyArraySourceLen(Program.VN[i][NJ], Program.V1N[i][NJ]);
                    FastCopy.CopyArraySourceLen(Program.VN[i][NJ], Program.V2N[i][NJ]);
                    FastCopy.CopyArraySourceLen(Program.WN[i][NJ], Program.W1N[i][NJ]);
                    FastCopy.CopyArraySourceLen(Program.WN[i][NJ], Program.W2N[i][NJ]);
                }
            });

            //Parallel.For(1, NJ + 1, Program.pOptions, j =>
            Parallel.ForEach(Partitioner.Create(1, NJ + 1, NJ / Program.pOptions.MaxDegreeOfParallelism), range =>
            {
                for (int j = range.Item1; j < range.Item2; ++j)
                {
                    FastCopy.CopyArraySourceLen(Program.TN[1][j],     Program.T[1][j]);
                    FastCopy.CopyArraySourceLen(Program.QUN[1][j],   Program.QU[1][j]);
                    FastCopy.CopyArraySourceLen(Program.TN[NI][j],   Program.T[NI][j]);
                    FastCopy.CopyArraySourceLen(Program.QUN[NI][j], Program.QU[NI][j]);
                }
            });


            //Parallel.For(1, NI + 1, Program.pOptions, i =>
            Parallel.ForEach(Partitioner.Create(1, NI + 1, NI / Program.pOptions.MaxDegreeOfParallelism), range =>
            {
                for (int i = range.Item1; i < range.Item2; ++i)
                {
                    FastCopy.CopyArraySourceLen(Program.TN[i][1],     Program.T[i][1]);
                    FastCopy.CopyArraySourceLen(Program.QUN[i][1],   Program.QU[i][1]);
                    FastCopy.CopyArraySourceLen(Program.TN[i][NJ],   Program.T[i][NJ]);
                    FastCopy.CopyArraySourceLen(Program.QUN[i][NJ], Program.QU[i][NJ]);
                }
            });

            //Copy corner values
            //Corner points 1,1
            FastCopy.CopyArraySourceLen( Program.UN[2][2],  Program.UN[1][1]);
            FastCopy.CopyArraySourceLen( Program.VN[2][2],  Program.VN[1][1]);
            FastCopy.CopyArraySourceLen( Program.TN[2][2],  Program.TN[1][1]);
            FastCopy.CopyArraySourceLen(Program.QUN[2][2], Program.QUN[1][1]);
            Array.Clear(Program.WN[1][1]);

            //Corner points 1,NJ
            FastCopy.CopyArraySourceLen( Program.UN[2][NJ - 1],  Program.UN[1][NJ]);
            FastCopy.CopyArraySourceLen( Program.VN[2][NJ - 1],  Program.VN[1][NJ]);
            FastCopy.CopyArraySourceLen( Program.TN[2][NJ - 1],  Program.TN[1][NJ]);
            FastCopy.CopyArraySourceLen(Program.QUN[2][NJ - 1], Program.QUN[1][NJ]);
            Array.Clear(Program.WN[1][NJ]);

            //Corner points NI,1
            FastCopy.CopyArraySourceLen( Program.UN[NI - 1][2],  Program.UN[NI][1]);
            FastCopy.CopyArraySourceLen( Program.VN[NI - 1][2],  Program.VN[NI][1]);
            FastCopy.CopyArraySourceLen( Program.TN[NI - 1][2],  Program.TN[NI][1]);
            FastCopy.CopyArraySourceLen(Program.QUN[NI - 1][2], Program.QUN[NI][1]);
            Array.Clear(Program.WN[NI][1]);
            
            //Corner points NI,NJ
            FastCopy.CopyArraySourceLen( Program.UN[NI - 1][NJ],  Program.UN[NI][NJ]);
            FastCopy.CopyArraySourceLen( Program.VN[NI - 1][NJ],  Program.VN[NI][NJ]);
            FastCopy.CopyArraySourceLen( Program.TN[NI - 1][NJ - 1],  Program.TN[NI][NJ]);
            FastCopy.CopyArraySourceLen(Program.QUN[NI - 1][NJ - 1], Program.QUN[NI][NJ]);
            Array.Clear(Program.WN[NI][NJ]);
        }

    }
}
