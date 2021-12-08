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
using System.Linq;
using System.Collections.Immutable;

namespace GRAMM_2001
{
    partial class Program
    {
        /// <summary>
        /// Read the GRAMM.geb file
        /// </summary>
        private static void Read_Gramm_Geb()
        {
            try
            {
                using (FileStream fs = new FileStream("GRAMM.geb", FileMode.Open, FileAccess.Read, FileShare.Read))
                {
                    using (StreamReader myreader = new StreamReader(fs))
                    {
                        string[] text = new string[10];
                        text = myreader.ReadLine().Split(new char[] { '!' });
                        text[0] = text[0].Trim();
                        text[0] = text[0].Replace(".", decsep);
                        Program.NX = Convert.ToInt32(text[0]);
                        text = myreader.ReadLine().Split(new char[] { '!' });
                        text[0] = text[0].Trim();
                        text[0] = text[0].Replace(".", decsep);
                        Program.NY = Convert.ToInt32(text[0]);
                        text = myreader.ReadLine().Split(new char[] { '!' });
                        text[0] = text[0].Trim();
                        text[0] = text[0].Replace(".", decsep);
                        Program.NZ = Convert.ToInt32(text[0]);

                        if (myreader.EndOfStream == false)
                        {
                            text = myreader.ReadLine().Split(new char[] { '!' });
                            text[0] = text[0].Trim();
                            text[0] = text[0].Replace(".", decsep);
                            Program.GRAMM_West = Convert.ToSingle(text[0]);
                        }
                        if (myreader.EndOfStream == false)
                        {
                            text = myreader.ReadLine().Split(new char[] { '!' });
                        }
                        if (myreader.EndOfStream == false)
                        {
                            text = myreader.ReadLine().Split(new char[] { '!' });
                            text[0] = text[0].Trim();
                            text[0] = text[0].Replace(".", decsep);
                            Program.GRAMM_South = Convert.ToSingle(text[0]);
                        }
                    }
                }
            }
            catch
            {
                Console.WriteLine("Error when reading domain data from file 'GRAMM.geb' - Execution stopped - press ESC to continue");
                while (!(Console.KeyAvailable && Console.ReadKey(true).Key == ConsoleKey.Escape))
                {
                    ;
                }

                Environment.Exit(0);
            }
        }

        private static void Read_Chemistry()
        {
            try
            {
                using (FileStream fs = new FileStream("chemistry.txt", FileMode.Open, FileAccess.Read, FileShare.Read))
                {
                    using (StreamReader myreader = new StreamReader(fs))
                    {
                        string[] text = new string[10];
                        text = myreader.ReadLine().Split(new char[] { '!' });
                        text[0] = text[0].Trim();
                        Program.chemistry_mechanism = text[0];

                        text = myreader.ReadLine().Split(new char[] { '!' });
                        text[0] = text[0].Trim();
                        Program.Update_Chemistry = Convert.ToSingle(text[0]);
                        Program.Update_Chemistry_Threshold = Program.Update_Chemistry;

                        Program.chemistry = true;
                    }
                }
            }
            catch
            {
                Console.WriteLine("Error when reading file 'chemistry.txt' - Execution stopped - press ESC to continue");
                while (!(Console.KeyAvailable && Console.ReadKey(true).Key == ConsoleKey.Escape))
                {
                    ;
                }

                Environment.Exit(0);
            }
        }

        /// <summary>
        /// Define the size of all used arrays
        /// </summary>
        private static void Define_Arrays()
        {
            NX1 = NX + 1;
            NY1 = NY + 1;
            NZ1 = NZ + 1;
            NX2 = NX + 2;
            NY2 = NY + 2;
            NZ2 = NZ + 2;
            NZB = 7;
            IMETSTR = 0;
            IHOURO = 0;
            NPROFMAX = 51;

            Console.Write("Allocate memory.");

            AH = CreateArray<float[]>(NX1, () => new float[NY1]);
            AH_Bassins = CreateArray<float[]>(NX1, () => new float[NY1]);
            TPI = CreateArray<float[]>(NX1, () => new float[NY1]);
            VOL = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            AREAX = CreateArray<float[][]>(NX2, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            AREAY = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY2, () => new float[NZ1]));
            AREAZ = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ2]));
            AREA = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            AREAXYZ = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            AREAZX = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ2]));
            AREAZY = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ2]));
            AHE = CreateArray<float[][]>(NX2, () => CreateArray<float[]>(NY2, () => new float[NZ2]));
            ZSP = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            
            AHImm   = new ImmutableArray<float>[NX1];
            ZSPImm  = CreateArray<ImmutableArray<float>[]>(NX1, () => new ImmutableArray<float>[NY1]);
            VOLImm  = CreateArray<ImmutableArray<float>[]>(NX1, () => new ImmutableArray<float>[NY1]);
            AREAImm = CreateArray<ImmutableArray<float>[]>(NX1, () => new ImmutableArray<float>[NY1]);
            AREAXImm = CreateArray<ImmutableArray<float>[]>(NX2, () => new ImmutableArray<float>[NY1]);
            AREAYImm = CreateArray<ImmutableArray<float>[]>(NX1, () => new ImmutableArray<float>[NY2]);
            AREAZXImm = CreateArray<ImmutableArray<float>[]>(NX1, () => new ImmutableArray<float>[NY1]);
            AREAZYImm = CreateArray<ImmutableArray<float>[]>(NX1, () => new ImmutableArray<float>[NY1]);
            AREAZImm = CreateArray<ImmutableArray<float>[]>(NX1, () => new ImmutableArray<float>[NY1]);
            AREAXYZImm = CreateArray<ImmutableArray<float>[]>(NX1, () => new ImmutableArray<float>[NY1]);
            
            DDX = new float[NX1];
            DDY = new float[NY1];
            ZAX = new float[NX1];
            ZAY = new float[NY1];
            X = new float[NX2];
            Y = new float[NY2];
            Z = new float[NZ2];
            //LAND = CreateArray<double[]>(NX2, () => new double[NY2]);
            PBZZ = new double[NZ1];
            RHOBZZ = new double[NZ1];
            TABS = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            RSOLG = CreateArray<double[]>(NX1, () => new double[NY1]);
            GLOBRAD = CreateArray<double[]>(NX1, () => new double[NY1]);
            DT_SOL = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            DT_TERR = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            RL = CreateArray<double[]>(NX1, () => new double[NY1]);
            TG = CreateArray<double[]>(NX1, () => new double[NY1]);
            EPSG = CreateArray<double[]>(NX1, () => new double[NY1]);
            KST = CreateArray<int[]>(NX1, () => new int[NY1]);
            ZZ = new float[NZ1];
            Console.Write(".");
            ALBEDO = CreateArray<double[]>(NX1, () => new double[NY1]);
            CLOUDS = CreateArray<double[]>(NX1, () => new double[NY1]);
            SNOW = CreateArray<double[]>(NX1, () => new double[NY1]);
            ZPROF = new double[NPROFMAX];
            PPROF = new double[NPROFMAX];
            TPROF = new double[NPROFMAX];
            VNORM = new double[NPROFMAX];
            QVAP = new double[NPROFMAX];
            QCLD = new double[NPROFMAX];
            QRAIN = new double[NPROFMAX];
            QICE = new double[NPROFMAX];
            QSNOW = new double[NPROFMAX];
            U = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            V = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            W = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            RHO = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            UN = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            VN = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            WN = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            U1 = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            V1 = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            W1 = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            U2 = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            V2 = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            W2 = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            U1N = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            V1N = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            W1N = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            U2N = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            V2N = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            W2N = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            PN = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            DP = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            Console.Write(".");
            DDP1DX = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            DDP1DY = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            DDP1DZ = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            DDP2DX = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            DDP2DY = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            DDP2DZ = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            DPX = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            DPY = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            DPZ = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            TPDX = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            TPDY = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            TP = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            TPX = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            TPY = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            TPZ = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            SUX = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            SUY = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            SUZ = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            SUXYZ = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            U1NRHO = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            V1NRHO = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            W1NRHO = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            U2NRHO = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            V2NRHO = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            W2NRHO = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            PBZ = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            RHOBZ = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            PNBZKP = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            Console.Write(".");
            NBZKP = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            UG1 = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            UG2 = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            VG1 = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            VG2 = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            UG = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            VG = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            VISH = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            VISV = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            QU = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            QUN = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            WAT_VAP = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            WAT_VAPN = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            QBZ = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            QUG = CreateArray<double[]>(NX1, () => new double[NY1]);
            FW = CreateArray<double[]>(NX1, () => new double[NY1]);
            T = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            TN = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            
            TBZ = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            FACTOR = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            DISS = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            DISSN = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            TE = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            TEN = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            Console.Write(".");
            ZI = CreateArray<float[]>(NX1, () => new float[NY1]);
            UST = CreateArray<double[]>(NX1, () => new double[NY1]);
            USTV = CreateArray<double[]>(NX1, () => new double[NY1]);
            TST = CreateArray<double[]>(NX1, () => new double[NY1]);
            OL = CreateArray<float[]>(NX1, () => new float[NY1]);
            XWQ = CreateArray<float[]>(NX1, () => new float[NY1]);
            Z0 = CreateArray<float[]>(NX1, () => new float[NY1]);
            FAC = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            RITSCH = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ]));
            WQU = CreateArray<double[]>(NX1, () => new double[NY1]);
            AWQ = CreateArray<float[]>(NX1, () => new float[NY1]);
            TB = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZB + 1]));
            TBN = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZB + 1]));
            TBA = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZB + 1]));
            DWB = new float[NZB + 1];
            DZB = new float[NZB + 1];
            RHOB = CreateArray<float[]>(NX1, () => new float[NY1]);
            ALAMBDA = CreateArray<float[]>(NX1, () => new float[NY1]);
            DRY_AREA = CreateArray<bool[]>(NX1, () => new bool[NY1]);
            ALPHA = new double[NZ1];
            Console.Write(".");
            A_PS = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            B_PS = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            C_PS = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            ASOUTH_PS = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            AWEST_PS = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            ANORTH_PS = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            AEAST_PS = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            AP0_PS = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            //SU_PS = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ1]));
            AIM = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[2 * NZ1]));
            BIM = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[2 * NZ1]));
            CIM = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[2 * NZ1]));
            AW1 = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            AS1 = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            AE2 = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            AN2 = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            AP0 = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[NZ1]));
            Console.Write(".");
            //SU = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[2 * NZ1]));
            DIMU = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[2 * NZ1]));
            F1U = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[2 * NZ1]));
            F2U = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[2 * NZ1]));
            F1V = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[2 * NZ1]));
            F2V = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[2 * NZ1]));
            F1W = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[2 * NZ1]));
            F2W = CreateArray<float[][]>(NX1, () => CreateArray<float[]>(NY1, () => new float[2 * NZ1]));
            DIMV = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[2 * NZ1]));
            DIMW = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[2 * NZ1]));
            AP = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[2 * NZ1]));
            AB = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[2 * NZ1]));
            AT = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[2 * NZ1]));
            AS = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[2 * NZ1]));
            AN = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[2 * NZ1]));
            AW = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[2 * NZ1]));
            AE = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[2 * NZ1]));
            RADIATION = CreateArray<double[][]>(NX1, () => CreateArray<double[]>(NY1, () => new double[NZ]));
            Relax_Border_factor = CreateArray<float[]>(NX1, () => new float[NY1]);

            VDATA1 = CreateArray<float[]>(NX1, () => new float[NY1]);
            VDATA2 = CreateArray<float[]>(NX1, () => new float[NY1]);
            VDATA3 = CreateArray<float[]>(NX1, () => new float[NY1]);
            VDATA4 = CreateArray<float[]>(NX1, () => new float[NY1]);
            VDATA8 = CreateArray<float[]>(NX1, () => new float[NY1]);
            VDATA5 = new float[NZ1];
            VDATA6 = new float[NZ1];
            VDATA7 = new float[NZ1];
            VDATA9 = new float[NZ1];
            MASSOURCE = new double[11];

            Console.Write(".");
            U_TEMP = CreateArray<float[]>(NX1, () => new float[NY1]);
            V_TEMP = CreateArray<float[]>(NX1, () => new float[NY1]);
            W_TEMP = CreateArray<float[]>(NX1, () => new float[NY1]);
            USE1 = CreateArray<double[]>(NY1, () => new double[NZ1]);
            VSE = CreateArray<double[]>(NY1, () => new double[NZ1]);
            WSE = CreateArray<double[]>(NY1, () => new double[NZ1]);
            TSE = CreateArray<double[]>(NY1, () => new double[NZ1]);
            QUSE = CreateArray<double[]>(NY1, () => new double[NZ1]);
            USEN = CreateArray<double[]>(NY1, () => new double[NZ1]);
            VSEN = CreateArray<double[]>(NY1, () => new double[NZ1]);
            WSEN = CreateArray<double[]>(NY1, () => new double[NZ1]);
            TSEN = CreateArray<double[]>(NY1, () => new double[NZ1]);
            QUSEN = CreateArray<double[]>(NY1, () => new double[NZ1]);
            USW = CreateArray<double[]>(NY1, () => new double[NZ1]);
            VSW = CreateArray<double[]>(NY1, () => new double[NZ1]);
            WSW = CreateArray<double[]>(NY1, () => new double[NZ1]);
            TSW = CreateArray<double[]>(NY1, () => new double[NZ1]);
            QUSW = CreateArray<double[]>(NY1, () => new double[NZ1]);
            USWN = CreateArray<double[]>(NY1, () => new double[NZ1]);
            VSWN = CreateArray<double[]>(NY1, () => new double[NZ1]);
            WSWN = CreateArray<double[]>(NY1, () => new double[NZ1]);
            TSWN = CreateArray<double[]>(NY1, () => new double[NZ1]);
            QUSWN = CreateArray<double[]>(NY1, () => new double[NZ1]);
            Console.Write(".");
            USN = CreateArray<double[]>(NX1, () => new double[NZ1]);
            VSN = CreateArray<double[]>(NX1, () => new double[NZ1]);
            WSN = CreateArray<double[]>(NX1, () => new double[NZ1]);
            TSN = CreateArray<double[]>(NX1, () => new double[NZ1]);
            QUSN = CreateArray<double[]>(NX1, () => new double[NZ1]);
            USNN = CreateArray<double[]>(NX1, () => new double[NZ1]);
            VSNN = CreateArray<double[]>(NX1, () => new double[NZ1]);
            WSNN = CreateArray<double[]>(NX1, () => new double[NZ1]);
            TSNN = CreateArray<double[]>(NX1, () => new double[NZ1]);
            QUSNN = CreateArray<double[]>(NX1, () => new double[NZ1]);
            USS = CreateArray<double[]>(NX1, () => new double[NZ1]);
            VSS = CreateArray<double[]>(NX1, () => new double[NZ1]);
            WSS = CreateArray<double[]>(NX1, () => new double[NZ1]);
            TSS = CreateArray<double[]>(NX1, () => new double[NZ1]);
            QUSS = CreateArray<double[]>(NX1, () => new double[NZ1]);
            USSN = CreateArray<double[]>(NX1, () => new double[NZ1]);
            VSSN = CreateArray<double[]>(NX1, () => new double[NZ1]);
            WSSN = CreateArray<double[]>(NX1, () => new double[NZ1]);
            TSSN = CreateArray<double[]>(NX1, () => new double[NZ1]);
            QUSSN = CreateArray<double[]>(NX1, () => new double[NZ1]);
            stabilityclass = CreateArray<Int16[]>(NX1, () => new Int16[NY1]);
        }

        /// <summary>
        /// Read the file IIN.dat
        /// </summary>
        private static void Read_IIN_Dat()
        {
            try
            {
                using (FileStream fs = new FileStream("IIN.dat", FileMode.Open, FileAccess.Read, FileShare.Read))
                {
                    using (StreamReader myreader = new StreamReader(fs))
                    {
                        string[] text = new string[10];
                        text = myreader.ReadLine().Split(new char[] { '!', ':' });
                        text[1] = text[1].Trim();
                        text[1] = text[1].Replace(".", decsep);
                        NDATUM = Convert.ToInt32(text[1]);
                        text = myreader.ReadLine().Split(new char[] { '!', ':' });
                        text[1] = text[1].Trim();
                        text[1] = text[1].Replace(".", decsep);
                        ISTUNDE = Convert.ToInt32(text[1]);
                        text = myreader.ReadLine().Split(new char[] { '!', ':' });
                        text[1] = text[1].Trim();
                        text[1] = text[1].Replace(".", decsep);
                        DT = Convert.ToDouble(text[1]);
                        text = myreader.ReadLine().Split(new char[] { '!', ':' });
                        text[1] = text[1].Trim();
                        text[1] = text[1].Replace(".", decsep);
                        TLIMIT2 = Convert.ToDouble(text[1]);
                        text = myreader.ReadLine().Split(new char[] { '!', ':' });
                        text[1] = text[1].Trim();
                        text[1] = text[1].Replace(".", decsep);
                        IOUTPUT = Convert.ToInt32(text[1]);
                        text = myreader.ReadLine().Split(new char[] { '!', ':' });
                        text[1] = text[1].Trim();
                        text[1] = text[1].Replace(".", decsep);
                        DELW = Convert.ToDouble(text[1]);
                        text = myreader.ReadLine().Split(new char[] { '!', ':' });
                        text[1] = text[1].Trim();
                        text[1] = text[1].Replace(".", decsep);
                        QUINIT = Convert.ToDouble(text[1]);
                        text = myreader.ReadLine().Split(new char[] { '!', ':' });
                        text[1] = text[1].Trim();
                        text[1] = text[1].Replace(".", decsep);
                        ZSEEH = Convert.ToDouble(text[1]);
                        text = myreader.ReadLine().Split(new char[] { '!', ':' });
                        text[1] = text[1].Trim();
                        text[1] = text[1].Replace(".", decsep);
                        TINIT = Convert.ToDouble(text[1]);
                        text = myreader.ReadLine().Split(new char[] { '!', ':' });
                        text[1] = text[1].Trim();
                        text[1] = text[1].Replace(".", decsep);
                        TGRAD = Convert.ToDouble(text[1]);
                        text = myreader.ReadLine().Split(new char[] { '!', ':' });
                        text[1] = text[1].Trim();
                        text[1] = text[1].Replace(".", decsep);
                        ZNEUT = Convert.ToDouble(text[1]);
                        text = myreader.ReadLine().Split(new char[] { '!', ':' });
                        text[1] = text[1].Trim();
                        text[1] = text[1].Replace(".", decsep);
                        TBINIT = Convert.ToDouble(text[1]);
                        text = myreader.ReadLine().Split(new char[] { '!', ':' });
                        text[1] = text[1].Trim();
                        text[1] = text[1].Replace(".", decsep);
                        TBINIT1 = Convert.ToDouble(text[1]);
                        text = myreader.ReadLine().Split(new char[] { '!', ':' });
                        text[1] = text[1].Trim();
                        text[1] = text[1].Replace(".", decsep);
                        BGRAD = Convert.ToDouble(text[1]);
                        text = myreader.ReadLine().Split(new char[] { '!', ':' });
                        text[1] = text[1].Trim();
                        text[1] = text[1].Replace(".", decsep);
                        IRAD = Convert.ToInt32(text[1]);
                        text = myreader.ReadLine().Split(new char[] { '!', ':' });
                        text[1] = text[1].Trim();
                        text[1] = text[1].Replace(".", decsep);
                        IDEBUG = Convert.ToInt32(text[1]);
                        text = myreader.ReadLine().Split(new char[] { '!', ':' });
                        text[1] = text[1].Trim();
                        text[1] = text[1].Replace(".", decsep);
                        IFLAGS1 = Convert.ToInt32(text[1]);
                        text = myreader.ReadLine().Split(new char[] { '!', ':' });
                        text[1] = text[1].Trim();
                        text[1] = text[1].Replace(".", decsep);
                        IFLAGS2 = Convert.ToInt32(text[1]);
                        text = myreader.ReadLine().Split(new char[] { '!', ':' });
                        text[1] = text[1].Trim();
                        text[1] = text[1].Replace(".", decsep);
                        DIA = Convert.ToDouble(text[1]);
                        text = myreader.ReadLine().Split(new char[] { '!', ':' });
                        text[1] = text[1].Trim();
                        text[1] = text[1].Replace(".", decsep);
                        DIAI = Convert.ToDouble(text[1]);
                        text = myreader.ReadLine().Split(new char[] { '!', ':' });
                        text[1] = text[1].Trim();
                        text[1] = text[1].Replace(".", decsep);
                        RELAXV = Convert.ToDouble(text[1]);
                        text = myreader.ReadLine().Split(new char[] { '!', ':' });
                        text[1] = text[1].Trim();
                        text[1] = text[1].Replace(".", decsep);
                        RELAXT = Convert.ToDouble(text[1]);
                        text = myreader.ReadLine().Split(new char[] { '!', ':' });
                        text[1] = text[1].Trim();
                        text[1] = text[1].Replace(".", decsep);
                        ICATAFORC = Convert.ToInt32(text[1]);
                        text = myreader.ReadLine().Split(new char[] { '!', ':' });
                        text[1] = text[1].Trim();
                        text[1] = text[1].Replace(".", decsep);
                        ERA5BOUND = Convert.ToSingle(text[1]) * 3600;
                        ERA5BOUND_Threshold = ERA5BOUND;
                        text = myreader.ReadLine().Split(new char[] { '!', ':' });
                        text[1] = text[1].Trim();
                        text[1] = text[1].Replace(".", decsep);
                        IBOUND = Convert.ToInt32(text[1]);
                        text = myreader.ReadLine().Split(new char[] { '!', ':' });
                        text[1] = text[1].Trim();
                        text[1] = text[1].Replace(".", decsep);
                        ISTAT = Convert.ToInt32(text[1]);
                        text = myreader.ReadLine().Split(new char[] { '!', ':' });
                        text[1] = text[1].Trim();
                        text[1] = text[1].Replace(".", decsep);
                        REINITIALIZATION = Convert.ToSingle(text[1]) * 3600;
                        REINITIALIZATION_Threshold = REINITIALIZATION;
                        text = myreader.ReadLine().Split(new char[] { '!', ':' });
                        text[1] = text[1].Trim();
                        text[1] = text[1].Replace(".", decsep);
                        NESTGRAMM = Convert.ToInt32(text[1]);
                        text = myreader.ReadLine().Split(new char[] { '!', ':' });
                        text[1] = text[1].Trim();
                        text[1] = text[1].Replace(".", decsep);
                        IOUT = Convert.ToInt32(text[1]);
                    }
                }
            }
            catch
            {
                Console.WriteLine("File 'IIN.dat' not found - Execution stopped  - press ESC to continue");
                while (!(Console.KeyAvailable && Console.ReadKey(true).Key == ConsoleKey.Escape))
                {
                    ;
                }

                Environment.Exit(0);
            }
        }

        /// <summary>
        /// Read the file Receptor.dat
        /// </summary>
        private static void Read_Receptor_Dat()
        {
            try
            {
                recexist = File.Exists("Receptor_GRAMM.dat");
                if (recexist == true)
                {
                    FileStream fs = new FileStream("Receptor_GRAMM.dat", FileMode.Open, FileAccess.Read, FileShare.Read);
                    StreamReader myreader = new StreamReader(fs);
                    try
                    {
                        string[] text = new string[10];
                        text = myreader.ReadLine().Split(new char[] { ',', ' ' });
                        text[0] = text[0].Trim();
                        text[0] = text[0].Replace(".", decsep);
                        int irecall = Convert.ToInt32(text[0]);
                        int n = 0;

                        while (myreader.EndOfStream == false)
                        {
                            text = myreader.ReadLine().Split(new char[] { ',', ' ' });

                            if (text.Length > 3)
                            {
                                text[0] = text[0].Trim();
                                text[0] = text[0].Replace(".", decsep);
                                double _xrec = Convert.ToDouble(text[1].Replace(".", decsep));
                                double _yrec = Convert.ToDouble(text[2].Replace(".", decsep));
                                double _zrec = Convert.ToDouble(text[3].Replace(".", decsep));
                                int _inrec = Convert.ToInt32(Math.Floor((_xrec - IKOOA) / DDX[1]) + 1);
                                int _jnrec = Convert.ToInt32(Math.Floor((_yrec - JKOOA) / DDY[1]) + 1);
                                n++;

                                if ((_inrec >= NX) || (_inrec <= 0))
                                {
                                    Console.WriteLine("  Receptor nr. " + n.ToString() + " out of Domain in X - direction!");
                                }
                                else if ((_jnrec >= NY) || (_jnrec <= 0))
                                {
                                    Console.WriteLine("  Receptor nr. " + n.ToString() + " out of Domain in Y - direction!");
                                }
                                else // Receptor points inside model domain
                                {
                                    double zdummy = 0;
                                    int _knrec = 0;
                                    for (int k = 1; k <= NZ; k++)
                                    {
                                        zdummy = ZSP[_inrec][_jnrec][k] - AH[_inrec][_jnrec];
                                        if (zdummy > _zrec)
                                        {
                                            _knrec = k;
                                            break;
                                        }
                                    }
                                    Xrec.Add(_xrec);
                                    Yrec.Add(_yrec);
                                    Zrec.Add(_zrec);
                                    Urec.Add(0);
                                    Vrec.Add(0);
                                    Trec.Add(0);
                                    Globradrec.Add(0);
                                    Longradrec.Add(0);
                                    Soilheatfluxrec.Add(0);
                                    Sensheatfluxrec.Add(0);
                                    Latheatfluxrec.Add(0);
                                    inrec.Add(_inrec);
                                    jnrec.Add(_jnrec);
                                    knrec.Add(_knrec);
                                }
                            }
                        }
                    }
                    catch
                    {
                        Console.WriteLine("Error when reading receptor data from file 'Receptor_GRAMM.dat' - Execution stopped - press ESC to continue");
                        while (!(Console.KeyAvailable && Console.ReadKey(true).Key == ConsoleKey.Escape))
                        {
                            ;
                        }

                        Environment.Exit(0);
                    }
                    myreader.Close();
                    myreader.Dispose();
                    fs.Dispose();

                    // Create File GRAMM.dat 
                    if (recexist == true && IWetter_Console_First == 0 && File.Exists("Gramm.dat") == false) //  not if multi-instances are used
                    {
                        try
                        {
                            using (StreamWriter wr = new StreamWriter("GRAMM.dat", false))
                            {
                                for (int ianz = 0; ianz < Xrec.Count(); ianz++)
                                {
                                    wr.Write(" Receptor" + (ianz + 1).ToString() + ", x:" + Xrec[ianz] + ", y:" + Yrec[ianz] + ", z:" + Zrec[ianz] + ", , , , ,");
                                }
                                wr.WriteLine();
                                for (int ianz = 0; ianz < Xrec.Count(); ianz++)
                                {
                                    wr.Write(" V, Dir, Temp, Solrad, Terrrad, Soilflux, Sensflux, Latentflux, ");
                                }
                                wr.WriteLine();
                            }
                        }
                        catch { }
                    }
                    if (Urec.Count == 0)
                    {
                        Urec.Add(0);
                        Vrec.Add(0);
                        Trec.Add(0);
                        inrec.Add(1);
                        jnrec.Add(1);
                        knrec.Add(1);
                    }
                }
                else
                {
                    Urec.Add(0);
                    Vrec.Add(0);
                    Trec.Add(0);
                    inrec.Add(1);
                    jnrec.Add(1);
                    knrec.Add(1);
                }
            }
            catch
            {
            }
        }

        /// <summary>
        /// Reset array values to 0
        /// </summary>
        private static void clear_arrays()
        {
            Read_IIN_Dat();
            IDIV = 0;
            IDIV_Up = 0;        //5.4.2017 Ku Reset Counter for MASSOURCE Queue
            IDIV_LockDown2 = 30;
            IDIV_LockDown = 0;
            IDIV_LockUp = 0;
            IDIV_PingPong = 0;
            Program.STEIGUNG = 0;
            IDIV_LockRelax = 0; //5.4.2017 Ku Reset Counter for Relax decrease
            RELAXV = Relaxv_ori; //5.4.2017 Ku Reset relax
            RELAXT = Relaxt_ori; //5.4.2017 Ku Reset relax
            Program.SUMG = 0;
            Program.DIVSUM = 0;
            Program.MASSOURCE_Old = 0;
            Program.MASSOURCE_Act = 0;
            Program.MASSOURCE_FIRST = 0;
            Program.MASSOURCE_Queue.Clear();
            Divergence_Min = 10e9;
            Program.TerminalOut = 0;
            Program.POBEN = 0;
            Program.TBZ1 = 0;

            
            ClearJaggedArray(MASSOURCE);
            ClearJaggedArray(PBZZ);
            ClearJaggedArray(RHOBZZ);
            ClearJaggedArray(RSOLG);
            ClearJaggedArray(GLOBRAD);
            ClearJaggedArray(DT_SOL);
            ClearJaggedArray(DT_TERR);
            ClearJaggedArray(RL);
            ClearJaggedArray(TG);
            ClearJaggedArray(ALBEDO);
            ClearJaggedArray(CLOUDS);
            ClearJaggedArray(SNOW);
            ClearJaggedArray(ZPROF);
            ClearJaggedArray(PPROF);
            ClearJaggedArray(TPROF);
            ClearJaggedArray(VNORM);
            ClearJaggedArray(QVAP);
            ClearJaggedArray(QCLD);
            ClearJaggedArray(QRAIN);
            ClearJaggedArray(QICE);
            ClearJaggedArray(QSNOW);
            ClearJaggedArray(EPSG);
            ClearJaggedArray(UST);
            ClearJaggedArray(USTV);
            ClearJaggedArray(TST);
            ClearJaggedArray(XWQ);
            ClearJaggedArray(WQU);
            ClearJaggedArray(AWQ);
            ClearJaggedArray(RHOB);
            ClearJaggedArray(RHO);
            ClearJaggedArray(Program.RHOBZ);
            ClearJaggedArray(ALAMBDA);

            ClearJaggedArray(TBN);
            ClearJaggedArray(TBA);
            ClearJaggedArray(TN);
            ClearJaggedArray(TE);
            ClearJaggedArray(TBZ);

            ClearJaggedArray(U);
            ClearJaggedArray(V);
            ClearJaggedArray(W);
            ClearJaggedArray(U1);
            ClearJaggedArray(V1);
            ClearJaggedArray(W1);
            ClearJaggedArray(U2);
            ClearJaggedArray(V2);
            ClearJaggedArray(W2);
            ClearJaggedArray(U1N);
            ClearJaggedArray(V1N);
            ClearJaggedArray(W1N);
            ClearJaggedArray(U2N);
            ClearJaggedArray(V2N);
            ClearJaggedArray(W2N);
            ClearJaggedArray(DIMU);
            ClearJaggedArray(DIMV);
            ClearJaggedArray(DIMW);
            ClearJaggedArray(A_PS);
            ClearJaggedArray(B_PS);
            ClearJaggedArray(C_PS);
            ClearJaggedArray(ASOUTH_PS);
            ClearJaggedArray(AEAST_PS);
            ClearJaggedArray(ANORTH_PS);
            ClearJaggedArray(AWEST_PS);
            ClearJaggedArray(AP0_PS);
            ClearJaggedArray(AIM);
            ClearJaggedArray(BIM);
            ClearJaggedArray(CIM);
            ClearJaggedArray(U1NRHO);
            ClearJaggedArray(V1NRHO);
            ClearJaggedArray(W1NRHO);
            ClearJaggedArray(U2NRHO);
            ClearJaggedArray(V2NRHO);
            ClearJaggedArray(W2NRHO);
            ClearJaggedArray(SUX);
            ClearJaggedArray(SUY);
            ClearJaggedArray(SUZ);
            ClearJaggedArray(SUXYZ);
            ClearJaggedArray(DPX);
            ClearJaggedArray(DPY);
            ClearJaggedArray(DPZ);
            ClearJaggedArray(DP);
            ClearJaggedArray(PN);
            ClearJaggedArray(PBZ);
            ClearJaggedArray(PNBZKP);
            ClearJaggedArray(NBZKP);
            ClearJaggedArray(FACTOR);
            ClearJaggedArray(DISS);
            ClearJaggedArray(DISSN);
            ClearJaggedArray(DDP1DX);
            ClearJaggedArray(DDP1DY);
            ClearJaggedArray(DDP1DZ);
            ClearJaggedArray(DDP2DX);
            ClearJaggedArray(DDP2DY);
            ClearJaggedArray(DDP2DZ);
            ClearJaggedArray(T);
            ClearJaggedArray(TB);
            ClearJaggedArray(TEN);
            ClearJaggedArray(TABS);
            ClearJaggedArray(TP);
            ClearJaggedArray(TPDY);
            ClearJaggedArray(TPDX);
            ClearJaggedArray(TPX);
            ClearJaggedArray(TPY);
            ClearJaggedArray(TPZ);
            ClearJaggedArray(VISV);
            ClearJaggedArray(VISH);
            ClearJaggedArray(QU);
            ClearJaggedArray(QUN);
            ClearJaggedArray(FAC);
            ClearJaggedArray(RITSCH);
            ClearJaggedArray(USEN);
            ClearJaggedArray(VSEN);
            ClearJaggedArray(WSEN);
            ClearJaggedArray(USW);
            ClearJaggedArray(VSW);
            ClearJaggedArray(WSW);
            ClearJaggedArray(TSEN);
            ClearJaggedArray(WSWN);
            ClearJaggedArray(VSWN);
            ClearJaggedArray(USWN);
            ClearJaggedArray(QUSW);
            ClearJaggedArray(QUSWN);
            ClearJaggedArray(QUSE);
            ClearJaggedArray(TSWN);
            ClearJaggedArray(USSN);
            ClearJaggedArray(VSSN);
            ClearJaggedArray(WSSN);
            ClearJaggedArray(TSSN);
            ClearJaggedArray(TSS);
            ClearJaggedArray(USS);
            ClearJaggedArray(VSS);
            ClearJaggedArray(WSS);
            ClearJaggedArray(QUSNN);
            ClearJaggedArray(QUSSN);
            ClearJaggedArray(QUSS);
            ClearJaggedArray(USNN);
            ClearJaggedArray(VSNN);
            ClearJaggedArray(WSNN);
            ClearJaggedArray(TSNN);
            ClearJaggedArray(USN);
            ClearJaggedArray(VSN);
            ClearJaggedArray(WSN);
            ClearJaggedArray(TSN);
            ClearJaggedArray(AB);
            ClearJaggedArray(AE);
            ClearJaggedArray(AN);
            ClearJaggedArray(AP);
            ClearJaggedArray(AS);
            ClearJaggedArray(AT);
            ClearJaggedArray(AW);
            ClearJaggedArray(F1U);
            ClearJaggedArray(F2U);
            ClearJaggedArray(F1V);
            ClearJaggedArray(F2V);
            ClearJaggedArray(F1W);
            ClearJaggedArray(F2W);
            ClearJaggedArray(QBZ);
            ClearJaggedArray(UG);
            ClearJaggedArray(VG);
            ClearJaggedArray(AW1);
            ClearJaggedArray(AS1);
            ClearJaggedArray(AE2);
            ClearJaggedArray(AN2);
            ClearJaggedArray(Program.RADIATION);
            ClearJaggedArray(OL);
            ClearJaggedArray(stabilityclass);
            ClearJaggedArray(Program.QUG);
            ClearJaggedArray(AP0);
            ClearJaggedArray(WAT_VAP);
            ClearJaggedArray(WAT_VAPN);
            ClearJaggedArray(FW);
            ClearJaggedArray(ZI);
            ClearJaggedArray(Z0);
        }
        /// <summary>
        /// Clear jagged array
        /// </summary>
        private static void ClearJaggedArray(double[][][] array)
        {
            for (int i = 0; i < array.Length; i++)
            {
                for (int j = 0; j < array[i].Length; j++)
                {
                    Array.Clear(array[i][j], 0, array[i][j].Length);
                }
            }
        }
        /// <summary>
        /// Clear jagged array
        /// </summary>
        private static void ClearJaggedArray(float[][][] array)
        {
            for (int i = 0; i < array.Length; i++)
            {
                for (int j = 0; j < array[i].Length; j++)
                {
                    Array.Clear(array[i][j], 0, array[i][j].Length);
                }
            }
        }
        /// <summary>
        /// Clear jagged array
        /// </summary>
        private static void ClearJaggedArray(double[][] array)
        {
            for (int i = 0; i < array.Length; i++)
            {
                Array.Clear(array[i], 0, array[i].Length);
            }
        }
        /// <summary>
        /// Clear jagged array
        /// </summary>
        private static void ClearJaggedArray(float[][] array)
        {
            for (int i = 0; i < array.Length; i++)
            {
                Array.Clear(array[i], 0, array[i].Length);
            }
        }
        /// <summary>
        /// Clear jagged array
        /// </summary>
        private static void ClearJaggedArray(short[][] array)
        {
            for (int i = 0; i < array.Length; i++)
            {
                Array.Clear(array[i], 0, array[i].Length);
            }
        }
        /// <summary>
        /// Clear jagged array
        /// </summary>
        private static void ClearJaggedArray(int[][] array)
        {
            for (int i = 0; i < array.Length; i++)
            {
                Array.Clear(array[i], 0, array[i].Length);
            }
        }
        /// <summary>
        /// Clear jagged array
        /// </summary>
        private static void ClearJaggedArray(double[] array)
        {
            Array.Clear(array, 0, array.Length);
        }
        /// <summary>
        /// Clear jagged array
        /// </summary>
        private static void ClearJaggedArray(float[] array)
        {
            Array.Clear(array, 0, array.Length);
        }

        /// <summary>
        /// Read the number of max. used processor cores
        /// </summary>
        public static void Max_Proc_File_Read()
        {
            //Set the maximum number of threads to be used in each parallelized region
            try
            {
                using (FileStream fs = new FileStream("Max_Proc.txt", FileMode.Open, FileAccess.Read, FileShare.Read))
                {
                    using (StreamReader myreader = new StreamReader(fs))
                    {
                        string text = myreader.ReadLine();
                        int IPROC = Convert.ToInt32(text);
                        IPROC = Math.Min(IPROC, Math.Min(NX, NY) / 2); // limit IPROC to 1/2 of lowest cell number count
                        IPROC = Math.Min(Environment.ProcessorCount, IPROC); // limit IPROC to Environment.ProcessorCount
                        pOptions.MaxDegreeOfParallelism = IPROC;
                    }
                }
            }
            catch
            { }
        }

        /// <summary>
        /// Read the file GRAMMin.dat, if iWetter == 1 then intialize computation
        /// </summary>
        /// <param name="iWetter"></param>
        public static void GRAMMin_File_Read(int iWetter)
        {
            try
            {
                using (FileStream fs = new FileStream("GRAMMin.dat", FileMode.Open, FileAccess.Read, FileShare.Read))
                {
                    using (StreamReader myreader = new StreamReader(fs))
                    {
                        string text;
                        text = myreader.ReadLine();
                        if (text.Contains("Version") == true) // new version of GRAMMin.dat
                        {
                            text = myreader.ReadLine();
                            if (iWetter == 1)
                            {
                                Program.METEO = text;
                            }

                            text = myreader.ReadLine();
                            if (iWetter == 1)
                            {
                                Program.Rauigkeit = Convert.ToDouble(text.Replace(".", Program.decsep));
                            }

                            string[] text5 = new string[10];
                            text5 = myreader.ReadLine().Split(new char[] { ' ', ',', '\t', ';', '!' }, StringSplitOptions.RemoveEmptyEntries);
                            if (iWetter == 1)
                            {
                                Program.IWETTER = Convert.ToInt16(text5[0]);
                                Program.nr_cell_smooth = Convert.ToInt16(text5[1]);
                                if (IWetter_Console_First > 0) // 11.4.17 Ku use arguments
                                {
                                    Program.IWETTER = IWetter_Console_First;
                                }
                            }

                            if (!myreader.EndOfStream)
                            {
                                text = myreader.ReadLine();
                                if (text.Contains("yes"))
                                {
                                    Program.WriteSteadyState = true;
                                    if (iWetter == 1)
                                    {
                                        Console.WriteLine("Write _steady_state.txt: yes");
                                    }
                                }
                                else
                                {
                                    Program.WriteSteadyState = false;
                                    if (iWetter == 1)
                                    {
                                        Console.WriteLine("Write _steady_state.txt: no");
                                    }
                                }
                            }
                            else
                            {
                                Program.WriteSteadyState = false;
                            }

                            //read original number of weather situations stored in meteopgt.all
                            Program.meteopgt_nr = 0;
                            if (!myreader.EndOfStream)
                            {
                                text5 = myreader.ReadLine().Split(new char[] { ' ', ',', '\t', ';', '!' }, StringSplitOptions.RemoveEmptyEntries);
                                Program.meteopgt_nr = Convert.ToInt16(text5[0]);
                            }
                        }
                        else if (iWetter == 1) // just the first time
                        {
                            text = text.Replace(".", Program.decsep);
                            Program.METEO = text;
                            text = myreader.ReadLine();
                            Program.Rauigkeit = Convert.ToDouble(text.Replace(".", Program.decsep));

                            text = myreader.ReadLine();
                            Program.IWETTER = Convert.ToInt16(text);

                            if (IWetter_Console_First > 0) // 11.4.17 Ku use arguments
                            {
                                Program.IWETTER = IWetter_Console_First;
                            }
                        }
                    }
                }
            }
            catch
            {
                Console.WriteLine("Error when reading file 'GRAMMin.dat' - Execution stopped");
                Environment.Exit(0);
            }
        }

    }
}
