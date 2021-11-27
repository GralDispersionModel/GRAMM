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
using System.Threading.Tasks;
using System.Collections.Immutable;

namespace GRAMM_2001
{
    partial class Program
    {
        /// <summary>
        /// Read the file ggeom.asc
        /// </summary>
        public static void GEOM()
        {
            //read values from file ggeom.asc
            Console.Write("Reading ggeom.asc....");
            string[] Is_Binary = new string[1];

            using (StreamReader isbin = new StreamReader("ggeom.asc"))
            {
                Is_Binary = isbin.ReadLine().Split(new char[] { ' ', ',', '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
            }

            Console.Write(".");

            int NX = 0;
            int NY = 0;
            int NZ = 0;
            int n = 0;

            if (Convert.ToDouble(Is_Binary[0]) < 0) // read binary format of ggeom.asc
            {
                using (BinaryReader readbin = new BinaryReader(File.Open("ggeom.asc", FileMode.Open, FileAccess.Read, FileShare.Read)))
                {
                    // read 1st line inclusive carriage return and line feed
                    byte[] header;
                    header = readbin.ReadBytes(6);

                    //obtain array size in x,y,z direction
                    NX = readbin.ReadInt32();
                    NY = readbin.ReadInt32();
                    NZ = readbin.ReadInt32();

                    // read AH[] array
                    for (int j = 1; j < NY + 1; j++)
                    {
                        for (int i = 1; i < NX + 1; i++)
                        {
                            Program.AH[i][j] = readbin.ReadSingle();
                        }
                    }
                    // read ZSP[] array
                    for (int k = 1; k < NZ + 1; k++)
                    {
                        for (int j = 1; j < NY + 1; j++)
                        {
                            for (int i = 1; i < NX + 1; i++)
                            {
                                Program.ZSP[i][j][k] = readbin.ReadSingle();
                            }
                        }
                    }

                    //obtain X, Y, and Z
                    for (int i = 1; i < NX + 2; i++)
                    {
                        Program.X[i] = readbin.ReadSingle();
                    }
                    for (int i = 1; i < NY + 2; i++)
                    {
                        Program.Y[i] = readbin.ReadSingle();
                    }
                    for (int i = 1; i < NZ + 2; i++)
                    {
                        Program.Z[i] = readbin.ReadSingle();
                    }

                    //obtain grid volumes
                    for (int k = 1; k < NZ + 1; k++)
                    {
                        for (int j = 1; j < NY + 1; j++)
                        {
                            for (int i = 1; i < NX + 1; i++)
                            {
                                Program.VOL[i][j][k] = readbin.ReadSingle();
                            }
                        }
                    }
                    Console.Write(".");

                    //obtain areas in x-direction
                    for (int k = 1; k < NZ + 1; k++)
                    {
                        for (int j = 1; j < NY + 1; j++)
                        {
                            for (int i = 1; i < NX + 2; i++)
                            {
                                Program.AREAX[i][j][k] = readbin.ReadSingle();
                            }
                        }
                    }
                    Console.Write(".");

                    //obtain areas in y-direction
                    for (int k = 1; k < NZ + 1; k++)
                    {
                        for (int j = 1; j < NY + 2; j++)
                        {
                            for (int i = 1; i < NX + 1; i++)
                            {
                                Program.AREAY[i][j][k] = readbin.ReadSingle();
                            }
                        }
                    }
                    Console.Write(".");

                    //obtain projection of z-area in x-direction
                    for (int k = 1; k < NZ + 2; k++)
                    {
                        for (int j = 1; j < NY + 1; j++)
                        {
                            for (int i = 1; i < NX + 1; i++)
                            {
                                Program.AREAZX[i][j][k] = readbin.ReadSingle();
                            }
                        }
                    }
                    Console.Write(".");

                    //obtain projection of z-area in y-direction
                    for (int k = 1; k < NZ + 2; k++)
                    {
                        for (int j = 1; j < NY + 1; j++)
                        {
                            for (int i = 1; i < NX + 1; i++)
                            {
                                Program.AREAZY[i][j][k] = readbin.ReadSingle();
                            }
                        }
                    }
                    Console.Write(".");

                    //obtain area in z-direction
                    for (int k = 1; k < NZ + 2; k++)
                    {
                        for (int j = 1; j < NY + 1; j++)
                        {
                            for (int i = 1; i < NX + 1; i++)
                            {
                                Program.AREAZ[i][j][k] = readbin.ReadSingle();
                            }
                        }
                    }
                    //obtain grid cells sizes in x-direction
                    for (int i = 1; i < NX + 1; i++)
                    {
                        Program.DDX[i] = readbin.ReadSingle();
                    }
                    //obtain grid cells sizes in y-direction
                    for (int i = 1; i < NY + 1; i++)
                    {
                        Program.DDY[i] = readbin.ReadSingle();
                    }
                    //obtain distances of neighbouring grid cells in x-direction
                    for (int i = 1; i < NX + 1; i++)
                    {
                        Program.ZAX[i] = readbin.ReadSingle();
                    }
                    //obtain distances of neighbouring grid cells in y-direction
                    for (int i = 1; i < NY + 1; i++)
                    {
                        Program.ZAY[i] = readbin.ReadSingle();
                    }

                    //obtain western and southern borders of the model domain and the angle (not used anymore) between the main model domain axis and north
                    IKOOA = readbin.ReadInt32();
                    JKOOA = readbin.ReadInt32();
                    double winkel = readbin.ReadDouble(); // angle not used
                } // Using

            }
            else // read ascii format of ggeom.asc
            {
                using (FileStream fs = new FileStream("ggeom.asc", FileMode.Open, FileAccess.Read, FileShare.Read))
                {
                    using (StreamReader reader = new StreamReader(fs))
                    {
                        int count = 0;

                        Is_Binary = reader.ReadLine().Split(new char[] { ' ', ',', '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
                        count = count + Is_Binary.Length;

                        //obtain array sizes in x,y,z direction
                        NX = Convert.ToInt32(Is_Binary[0]);
                        NY = Convert.ToInt32(Is_Binary[1]);
                        NZ = Convert.ToInt32(Is_Binary[2]);

                        //obtain X, Y, and Z
                        for (int i = 3; i < 3 + NX + 1; i++)
                        {
                            Program.X[i - 2] = Convert.ToSingle(Is_Binary[i].Replace(".", Program.decsep));
                        }
                        n = 0;
                        for (int i = 3 + NX + 1; i < 3 + NX + 1 + NY + 1; i++)
                        {
                            n++;
                            Program.Y[n] = Convert.ToSingle(Is_Binary[i].Replace(".", Program.decsep));
                        }
                        n = 0;
                        for (int i = 3 + NX + 1 + NY + 1; i < Is_Binary.Length; i++)
                        {
                            n++;
                            Program.Z[n] = Convert.ToSingle(Is_Binary[i].Replace(".", Program.decsep));
                        }

                        int ix = 0;
                        int iy = 1;
                        //obtain surface heights
                        n = 0;
                        Is_Binary = reader.ReadLine().Split(new char[] { ' ', ',', '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
                        Console.Write(".");
                        for (int j = 1; j < NY + 1; j++)
                        {
                            for (int i = 1; i < NX + 1; i++)
                            {
                                Program.AH[i][j] = Convert.ToSingle(Is_Binary[n].Replace(".", Program.decsep));
                                n++;
                            }
                        }
                        //obtain grid volumes
                        n = 0;
                        Is_Binary = reader.ReadLine().Split(new char[] { ' ', ',', '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
                        Console.Write(".");
                        for (int k = 1; k < NZ + 1; k++)
                        {
                            for (int j = 1; j < NY + 1; j++)
                            {
                                for (int i = 1; i < NX + 1; i++)
                                {
                                    Program.VOL[i][j][k] = Convert.ToSingle(Is_Binary[n].Replace(".", Program.decsep));
                                    n++;
                                }
                            }
                        }
                        //obtain areas in x-direction
                        n = 0;
                        Is_Binary = reader.ReadLine().Split(new char[] { ' ', ',', '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
                        Console.Write(".");
                        for (int k = 1; k < NZ + 1; k++)
                        {
                            for (int j = 1; j < NY + 1; j++)
                            {
                                for (int i = 1; i < NX + 2; i++)
                                {
                                    Program.AREAX[i][j][k] = Convert.ToSingle(Is_Binary[n].Replace(".", Program.decsep));
                                    n++;
                                }
                            }
                        }
                        //obtain areas in y-direction
                        n = 0;
                        Is_Binary = reader.ReadLine().Split(new char[] { ' ', ',', '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
                        Console.Write(".");
                        for (int k = 1; k < NZ + 1; k++)
                        {
                            for (int j = 1; j < NY + 2; j++)
                            {
                                for (int i = 1; i < NX + 1; i++)
                                {
                                    Program.AREAY[i][j][k] = Convert.ToSingle(Is_Binary[n].Replace(".", Program.decsep));
                                    n++;
                                }
                            }
                        }
                        //obtain projection of z-area in x-direction
                        n = 0;
                        Is_Binary = reader.ReadLine().Split(new char[] { ' ', ',', '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
                        Console.Write(".");
                        for (int k = 1; k < NZ + 2; k++)
                        {
                            for (int j = 1; j < NY + 1; j++)
                            {
                                for (int i = 1; i < NX + 1; i++)
                                {
                                    Program.AREAZX[i][j][k] = Convert.ToSingle(Is_Binary[n].Replace(".", Program.decsep));
                                    n++;
                                }
                            }
                        }
                        //obtain projection of z-area in y-direction
                        n = 0;
                        Is_Binary = reader.ReadLine().Split(new char[] { ' ', ',', '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
                        Console.Write(".");
                        for (int k = 1; k < NZ + 2; k++)
                        {
                            for (int j = 1; j < NY + 1; j++)
                            {
                                for (int i = 1; i < NX + 1; i++)
                                {
                                    Program.AREAZY[i][j][k] = Convert.ToSingle(Is_Binary[n].Replace(".", Program.decsep));
                                    n++;
                                }
                            }
                        }
                        //obtain area in z-direction
                        n = 0;
                        Is_Binary = reader.ReadLine().Split(new char[] { ' ', ',', '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
                        Console.Write(".");
                        for (int k = 1; k < NZ + 2; k++)
                        {
                            for (int j = 1; j < NY + 1; j++)
                            {
                                for (int i = 1; i < NX + 1; i++)
                                {
                                    Program.AREAZ[i][j][k] = Convert.ToSingle(Is_Binary[n].Replace(".", Program.decsep));
                                    n++;
                                }
                            }
                        }
                        //obtain cell-heights
                        n = 0;
                        Is_Binary = reader.ReadLine().Split(new char[] { ' ', ',', '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
                        Console.Write(".");
                        for (int k = 1; k < NZ + 1; k++)
                        {
                            for (int j = 1; j < NY + 1; j++)
                            {
                                for (int i = 1; i < NX + 1; i++)
                                {
                                    Program.ZSP[i][j][k] = Convert.ToSingle(Is_Binary[n].Replace(".", Program.decsep));
                                    n++;
                                }
                            }
                        }
                        //obtain grid cells sizes in x-direction
                        n = 0;
                        Is_Binary = reader.ReadLine().Split(new char[] { ' ', ',', '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
                        Console.Write(".");
                        for (int i = 1; i < NX + 1; i++)
                        {
                            Program.DDX[i] = Convert.ToSingle(Is_Binary[n].Replace(".", Program.decsep));
                            n++;
                        }
                        //obtain grid cells sizes in y-direction
                        n = 0;
                        Is_Binary = reader.ReadLine().Split(new char[] { ' ', ',', '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
                        Console.Write(".");
                        for (int i = 1; i < NY + 1; i++)
                        {
                            Program.DDY[i] = Convert.ToSingle(Is_Binary[n].Replace(".", Program.decsep));
                            n++;
                        }
                        //obtain distances of neighbouring grid cells in x-direction
                        n = 0;
                        Is_Binary = reader.ReadLine().Split(new char[] { ' ', ',', '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
                        Console.Write(".");
                        for (int i = 1; i < NX + 1; i++)
                        {
                            Program.ZAX[i] = Convert.ToSingle(Is_Binary[n].Replace(".", Program.decsep));
                            n++;
                        }
                        //obtain distances of neighbouring grid cells in y-direction
                        n = 0;
                        Is_Binary = reader.ReadLine().Split(new char[] { ' ', ',', '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
                        Console.Write(".");
                        for (int i = 1; i < NY + 1; i++)
                        {
                            Program.ZAY[i] = Convert.ToSingle(Is_Binary[n].Replace(".", Program.decsep));
                            n++;
                        }
                        //obtain western and southern borders of the model domain and the angle (not used anymore) between the main model domain axis and north
                        n = 0;
                        Is_Binary = reader.ReadLine().Split(new char[] { ' ', ',', '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);

                        Program.IKOOA = Convert.ToInt32(Is_Binary[0]);
                        Program.JKOOA = Convert.ToInt32(Is_Binary[1]);
                        double winkel = Convert.ToDouble(Is_Binary[2].Replace(".", Program.decsep));
                    }
                }
            }

            Console.Write("I");
            Program.AHImm = ImmutableArray.Create(Program.AH);
            Program.ZSPImm = ImmutableArray.Create(Program.ZSP);
            Program.VOLImm = ImmutableArray.Create(Program.VOL);
            Program.DDXImm = ImmutableArray.Create(Program.DDX);
            Program.DDYImm = ImmutableArray.Create(Program.DDY);
            Program.AREAImm = ImmutableArray.Create(Program.AREA);
            Program.AREAXImm = ImmutableArray.Create(Program.AREAX);
            Program.AREAYImm = ImmutableArray.Create(Program.AREAY);
            Program.AREAZXImm = ImmutableArray.Create(Program.AREAZX);
            Program.AREAZYImm = ImmutableArray.Create(Program.AREAZY);
            Program.AREAZImm = ImmutableArray.Create(Program.AREAZ);
            Program.AREAXYZImm = ImmutableArray.Create(Program.AREAXYZ);

            Console.WriteLine("I");

            //check whether array dimensions given in the input file GRAMM.geb meets the one of the input file ggeom.asc
            if (NX != Program.NX)
            {
                Console.WriteLine("Field dimensions in x-direction of 'GRAMM.geb' and 'ggeom.asc' are different");
                return;
            }
            if (NY != Program.NY)
            {
                Console.WriteLine("Field dimensions in y-direction of 'GRAMM.geb' and 'ggeom.asc' are different");
                return;
            }
            if (NZ != Program.NZ)
            {
                Console.WriteLine("Field dimensions in z-direction of 'GRAMM.geb' and 'ggeom.asc' are different");
                return;
            }

            //definition of vertical cells of the soil model
            double DZB0 = 0.02 / 2.3;
            Program.DWB[1] = 0;
            for (n = 2; n <= Program.NZB; n++)
            {
                Program.DWB[n] = (float)(DZB0 * Math.Pow(2.3, n - 2));
            }

            for (n = 2; n <= Program.NZB - 1; n++)
            {
                Program.DZB[n] = 0.5F * (Program.DWB[n] + Program.DWB[n + 1]);
            }

            Program.DZB[1] = Program.DZB[2];
            Program.DZB[Program.NZB] = Program.DZB[Program.NZB - 1];

            //projection of z-surface in z-direction
            Parallel.For(1, Program.NX + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= Program.NY; j++)
                {
                    for (int k = 1; k <= Program.NZ; k++)
                    {
                        Program.AREA[i][j][k] = Program.DDX[i] * Program.DDYImm[j];
                    }
                }
            });

            //computation of minimum and maximum surface elevation and of the area separating the two half-cells
            Program.AHMIN = 100000;
            Program.AHMAX = -10000;
            double AHMIN_D = 100000;
            Parallel.For(1, Program.NX + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= Program.NY; j++)
                {
                    if (Program.AH[i][j] < Program.AHMIN)
                    {
                        Program.AHMIN = Program.AH[i][j];
                        Program.AHMINI = i;
                        Program.AHMINJ = j;
                    }
                    if (i > 2 && j > 2 && i < Program.NX && j < Program.NY)
                    {
                        if (Program.AH[i][j] < AHMIN_D)
                        {
                            AHMIN_D = Program.AH[i][j];
                            Program.AHMINI_D = i;
                            Program.AHMINJ_D = j;
                        }
                    }
                    if (Program.AH[i][j] > Program.AHMAX)
                    {
                        Program.AHMAX = Program.AH[i][j];
                    }
                    for (int k = 1; k <= Program.NZ; k++)
                    {
                        Program.AREAXYZ[i][j][k] = (float)Math.Sqrt(Pow2(Program.AREAX[i][j][k] + Program.AREAZX[i][j][k]) + Pow2(Program.AREAY[i][j][k] + Program.AREAZY[i][j][k]) + Pow2(Program.AREA[i][j][k]));
                    }
                }
            });
            Console.WriteLine("Minimum and maximum surface elevations: " + Convert.ToString(Math.Round(Program.AHMIN, 0)) + "m  " + Convert.ToString(Math.Round(Program.AHMAX, 0)) + "m ");

            //computation of geometry data for the radiation model
            Parallel.For(1, Program.NX + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= Program.NY; j++)
                {
                    Program.KST[i][j] = 2;
                    for (int k = 1; k <= Program.NZ; k++)
                    {
                        double DELZP = 0.5F * (Program.Z[k + 1] - Program.Z[k]);
                        double DELZM;
                        if (k > 1)
                        {
                            DELZM = 0.5F * (Program.Z[k] - Program.Z[k - 1]);
                        }
                        else
                        {
                            DELZM = 0;
                        }
                        Program.ZZ[1] = Program.Z[1];
                        if (k < Program.NZ)
                        {
                            Program.ZZ[k + 1] = 0.5F * (Program.Z[k + 1] + Program.Z[k]);
                        }

                        if ((Program.AH[i][j] >= (Program.Z[k] - DELZM)) && (Program.AH[i][j] < (Program.Z[k] - DELZP)))
                        {
                            Program.KST[i][j] = Math.Max(k + 1, 1);
                        }

                        for (int kk = 1; kk <= Program.NZ; kk++)
                        {
                            if (k > 1)
                            {
                                if ((Program.ZSP[i][j][kk] >= Program.ZZ[k - 1]) && (Program.ZSP[i][j][kk] < Program.ZZ[k]))
                                {
                                    Program.NBZKP[i][j][kk] = k;
                                    Program.PNBZKP[i][j][kk] = (Program.ZSP[i][j][kk] - Program.ZZ[k - 1]) / (Program.ZZ[k] - Program.ZZ[k - 1]);
                                }
                                else if (Program.ZSP[i][j][kk] >= Program.ZZ[Program.NZ])
                                {
                                    Program.NBZKP[i][j][kk] = Program.NZ;
                                    Program.PNBZKP[i][j][kk] = 1;
                                }
                            }
                            else
                            {
                                if (Program.ZSP[i][j][kk] < Program.ZZ[k])
                                {
                                    Program.NBZKP[i][j][kk] = k;
                                    Program.PNBZKP[i][j][kk] = 1;
                                }
                            }
                        }
                    }
                }
            });

        }


    }
}
