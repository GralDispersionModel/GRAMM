#region Copyright
///<remarks>
/// <GRAMM Mesoscale Model>
/// Copyright (C) [2020]  [Dietmar Oettl, Markus Kuntner]
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

namespace GRAMM_2001
{
    ///<summary>
    ///This class computes shortwave incoming solar radiation and longwave terrestrial radiation
    ///</summary>
    class RadiationCalculation
    {
        ///<summary>
        ///Number of days since 1 Jan
        ///</summary>
        private int Na = 0;
        ///<summary>
        ///Time after midnight in seconds
        ///</summary>
        private double TimeR = 0;
        ///<summary>
        ///Vertical index of the highest surface cell
        ///</summary>
        private int KStMax = 0;
        ///<summary>
        ///General lower threshold
        ///</summary>
        private const double Eps = 0.00005;
        ///<summary>
        ///Number of sectors used to compute the horizons of the surroundings (should be divisible by 4)
        ///</summary>
        private const int NumberOmega = 40;
        ///<summary>
        ///Variation of the solar flux over the year
        ///</summary>
        private float GPhi = 0;
        ///<summary>
        ///Sinus of sun heigth (He)
        ///</summary>
        private float Mye = 0;
        ///<summary>
        ///Sinus of sun heigth (H)
        ///</summary>
        private float My = 0;
        ///<summary>
        ///Cosinus He(Mye)
        ///</summary>
        private float CosHe = 0;
        ///<summary>
        ///Cosinus He(My)
        ///</summary>
        private float CosH = 0;
        ///<summary>
        ///Angle of sun height [rad]
        ///</summary>
        private float He = 0;
        ///<summary>
        ///Azimuth of the sun [rad]
        ///</summary>
        private float Eta = 0;
        ///<summary>
        ///Conversion factor PI/180 in rad
        ///</summary>
        private const double GrRa = 0.0174532925199;
        ///<summary>
        ///Conversion factor 180/PI in deg
        ///</summary>
        private const double RaGr = 57.2957795131;
        ///<summary>
        ///reference pressure
        ///</summary>
        private const double P0 = 1.01325E5;
        ///<summary>
        ///Eulerian constant
        ///</summary>
        private const double Eulerian = 2.7182818285;
        ///<summary>
        ///extraterrestrial solar radiation density [W/m²]
        ///</summary>
        private const double IG0 = 1367;
        ///<summary>
        ///heat capacity of air [Joule/kg Kelvin]
        ///</summary>
        private const double Cp = 1004;
        ///<summary>
        ///mass-concentration of CO2
        ///</summary>
        private const double QCo2 = 0.00052;
        ///<summary>
        ///W-integral of black cloud [kg/m²]
        ///</summary>
        private const double Wcrit = 0.03;
        ///<summary>
        ///Epsilon at the top of the atmosphere
        ///</summary>
        private const double EpsTop = 0.2;
        ///<summary>
        ///Temperature at the top of the atmosphere
        ///</summary>
        private const double TTop = 216.65;
        ///<summary>
        ///Default surface albedo
        ///</summary>
        private const double AgNorm = 0.6;
        ///<summary>
        ///Intermediate value for computing radiation density
        ///</summary>
        private double A0 = 0;
        ///<summary>
        ///Inclination angle of the surface fall-line (=line of greates slope)
        ///</summary>
        private readonly float[][] Alfa;
        ///<summary>
        ///Azimuth angle of fall-line (=line of greates slope)
        ///</summary>
        private readonly float[][] Beta;
        // North = 0
        ///<summary>
        ///1-Q = Sky-view factor (Horizontüberhöhung) in the radiation model
        ///</summary>
        private readonly double[][] Q;
        ///<summary>
        /// Integrated visibility in the radiation model
        ///</summary>
        private readonly double[][] YG;
        ///<summary>
        ///Integrated water vapour in the radiation model
        ///</summary>
        private readonly double[][] R;
        ///<summary>
        ///influence of clouds on the solar incoming radiation
        ///</summary>
        private readonly double[][] W1rad;
        //private  double[][][][] W1rad = Program.CreateArray<double[][][]>(1, () => Program.CreateArray<double[][]>(1, () => Program.CreateArray<double[]>(1, () => new double[1])));
        ///<summary>
        /// influence of rain on the solar incoming radiation
        ///</summary>
        private readonly double[][] W2rad;
        ///<summary>
        ///? something in the radiation model
        ///</summary>
        private readonly double[][] Wrad;
        ///<summary>
        ///Lower cell height of the radiation model
        ///</summary>
        private readonly double[] Dz1 = new double[1];
        ///<summary>
        ///Upper cell height of the radiation model
        ///</summary>
        private readonly double[] Dz2 = new double[1];
        ///<summary>
        /// Surface albedo in the radiation model
        ///</summary>
        private readonly double[][] Ag;
        ///<summary>
        ///Linke turbidity coefficient (~1 - 8)
        ///</summary>
        private readonly double[] Tstern = new double[1];
        ///<summary>
        ///Cloud parameter
        ///</summary>
        private readonly float[] CloudES = new float[1];
        //private  double[][][] eS = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        /// Cloud parameter
        ///</summary>
        private readonly float[] CloudESeS = new float[1];
        //private  double[][][] eSeS = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Shortwave solar radiation
        ///</summary>
        private readonly double[] IG = new double[1];
        ///<summary>
        ///Total solar radiation density
        ///</summary>
        private readonly double[][][] EG;
        ///<summary>
        ///Direct solar radiation density
        ///</summary>
        private readonly double[][][] SG;
        ///<summary>
        ///Diffusive solar radiation density
        ///</summary>
        private readonly double[][][] DG;
        ///<summary>
        ///Transmission function
        ///</summary>
        private readonly double[] Ns = new double[1];
        ///<summary>
        ///Transmission function
        ///</summary>
        private readonly double[] Nd = new double[1];
        ///<summary>
        ///Temporary array for cloud
        ///</summary>
        private readonly double[] Tau = new double[1];
        ///<summary>
        /// Temporary array in the radiation model
        ///</summary>
        private readonly double[][][] Arad;
        ///<summary>
        ///Temporary array in the radiation model
        ///</summary>
        private readonly float[] Myx;
        ///<summary>
        /// Integrated water vapour
        ///</summary>
        private readonly double[][] Rp;
        ///<summary>
        ///Integrated CO2
        ///</summary>
        private readonly double[][] CGp;
        ///<summary>
        /// Downward terrestrial radiation
        ///</summary>
        private readonly double[] LStrich = new double[1];
        ///<summary>
        ///Total incoming terrestrial radiation
        ///</summary>
        private readonly double[][] RTerrG;
        ///<summary>
        ///Temporary array in the radiation model
        ///</summary>
        private readonly double[][][] EpsAp;
        ///<summary>
        ///Temporary array in the radiation model
        ///</summary>
        private readonly double[][][] EpsAm;
        ///<summary>
        ///Temperature of the black cloud (at top of the atmosphere)
        ///</summary>
        private readonly double[][][] Tab;
        ///<summary>
        ///Emissivity? of the black cloud (at top of the atmosphere)
        ///</summary>
        private readonly double[][][] Eab;

        ///<summary>
        ///Constructor: array declarations for the radiation model
        ///</summary>
        public RadiationCalculation()
        {
            Alfa = Program.CreateArray<float[]>(Program.NX1, () => new float[Program.NY1]);
            Beta = Program.CreateArray<float[]>(Program.NX1, () => new float[Program.NY1]);
            Q = Program.CreateArray<double[]>(Program.NX1, () => new double[Program.NY1]);
            YG = Program.CreateArray<double[]>(Program.NZ1, () => new double[4]);
            R = Program.CreateArray<double[]>(Program.NZ1, () => new double[4]);
            W1rad = Program.CreateArray<double[]>(Program.NZ1, () => new double[4]);
            Wrad = Program.CreateArray<double[]>(Program.NZ1, () => new double[4]);
            W2rad = Program.CreateArray<double[]>(Program.NZ1, () => new double[4]);
            Dz1 = new double[Program.NZ1];
            Dz2 = new double[Program.NZ1];
            Ag = Program.CreateArray<double[]>(Program.NX1, () => new double[Program.NY1]);
            Tstern = new double[Program.NZ1];
            CloudES = new float[Program.NZ1];
            CloudESeS = new float[Program.NZ1];
            //eS = Program.CreateArray<double[][]>(Program.NX1, () => Program.CreateArray<double[]>(Program.NY1, () => new double[Program.NZ1]));
            //eSeS = Program.CreateArray<double[][]>(Program.NX1, () => Program.CreateArray<double[]>(Program.NY1, () => new double[Program.NZ1]));
            IG = new double[Program.NZ1];
            EG = Program.CreateArray<double[][]>(Program.NX1, () => Program.CreateArray<double[]>(Program.NY1, () => new double[Program.NZ1]));
            SG = Program.CreateArray<double[][]>(Program.NX1, () => Program.CreateArray<double[]>(Program.NY1, () => new double[Program.NZ1]));
            DG = Program.CreateArray<double[][]>(Program.NX1, () => Program.CreateArray<double[]>(Program.NY1, () => new double[Program.NZ1]));
            Ns = new double[Program.NZ1];
            Nd = new double[Program.NZ1];
            Tau = new double[Program.NZ1];
            //nS = Program.CreateArray<double[][]>(Program.NX1, () => Program.CreateArray<double[]>(Program.NY1, () => new double[Program.NZ1]));
            //nD = Program.CreateArray<double[][]>(Program.NX1, () => Program.CreateArray<double[]>(Program.NY1, () => new double[Program.NZ1]));
            //Tau = Program.CreateArray<double[][]>(Program.NX1, () => Program.CreateArray<double[]>(Program.NY1, () => new double[Program.NZ1]));
            Arad = Program.CreateArray<double[][]>(Program.NX1, () => Program.CreateArray<double[]>(Program.NY1, () => new double[Program.NZ1]));
            Myx = new float[Program.NZ1];
            Rp = Program.CreateArray<double[]>(Program.NZ1, () => new double[4]);
            CGp = Program.CreateArray<double[]>(Program.NZ1, () => new double[4]);
            LStrich = new double[Program.NZ1];
            //L_Strich = Program.CreateArray<double[][]>(Program.NX1, () => Program.CreateArray<double[]>(Program.NY1, () => new double[Program.NZ1]));
            RTerrG = Program.CreateArray<double[]>(Program.NX1, () => new double[Program.NY1]);
            EpsAp = Program.CreateArray<double[][]>(Program.NX1, () => Program.CreateArray<double[]>(Program.NY1, () => new double[Program.NZ1]));
            EpsAm = Program.CreateArray<double[][]>(Program.NX1, () => Program.CreateArray<double[]>(Program.NY1, () => new double[Program.NZ1]));
            Tab = Program.CreateArray<double[][]>(Program.NX1, () => Program.CreateArray<double[]>(Program.NY1, () => new double[Program.NZ1]));
            Eab = Program.CreateArray<double[][]>(Program.NX1, () => Program.CreateArray<double[]>(Program.NY1, () => new double[Program.NZ1]));
        }

        ///<summary>
        ///This routine computes shortwave incoming solar radiation and longwave terrestrial radiation
        ///</summary>
        ///<param name="lGeom">Flag for deciding whether geometry data should be computed or not></param>
        ///<param name="lGeomWr">Flag for deciding whether geometry data should be stored in the file "albeq.dat" or not></param>
        ///<param name="lGeomRe">Flag for deciding whether geometry data should be read from file "albeq.dat" or not></param>
        ///<param name="iTag">Day for calculation></param>
        ///<param name="iMon">Month for calculation></param>
        ///<param name="iJahr">Year of calculation></param>
        ///<param name="Zeit">Time after midnight in seconds></param>
        ///<param name="NI">Number of cells in x direction></param>
        ///<param name="NJ">Number of cells in y direction></param>
        ///<param name="NK">Number of cells in z direction></param>
        public void RADIATRAD(bool lGeom, bool lGeomWr, bool lGeomRe, int iTag, int iMon, int iJahr, double Zeit, int NI, int NJ, int NK)
        {
            /*
             * This routine computes shortwave incoming solar radiation and longwave terrestrial radiation
             * as well as atmospheric heating rates
             * lGeom:       Flag for deciding whether geometry data should be computed or not
             * LGeomWr:     Flag for deciding whether geometry data should be stored in the file "albeq.dat" or not
             * lGeomRe:     Flag for deciding whether geometry data should be read from file "albeq.dat" or not
             * iTag, iMon,
             * iJahr:       Date
             * Zeit:        Time in seconds after midnight
             * KAlb:        Classification of albedo
             * iTheta:      Sector containing the suns azimuth
             * XP, XW, XE,..Coordinates of the surface cell and its neighbours
             * DOmega:      Angle in [rad] of one sector
             * Omega:       Average angle of actual sector in [rad]
             * Theta:       Excess-angle of cell [rad]
             * ThetaB:      Excess-angle of surface cell [rad]
             * ThetaOld:    Theta of previous loop
            */

            int iTheta; // int Kalb;

            if (NK != 3) Console.Write("RADIATION started");
            if ((Program.ISOL != 1) && (Program.ISOL != 2))
            {
                Console.WriteLine("===========================================");
                Console.WriteLine("Selected radiation scheme must be '1' or '2'");
                Console.WriteLine("===========================================");
                Environment.Exit(0);
            }

            string Date1 = "0101" + iJahr.ToString("0000");
            string Date2 = iTag.ToString("00") + iMon.ToString("00") + iJahr.ToString("0000");
            Na = 1 + DateNum(Date2) - DateNum(Date1);

            TimeR = Zeit;

            // Highest ground cell
            KStMax = 1;
            for (int i = 1; i <= NI; i++)
            {
                for (int j = 1; j <= NJ; j++)
                {
                    KStMax = Math.Max(KStMax, Program.KST[i][j]);
                }
            }

            //Computing geometry (inclination, exposition and projection of surface cells)
            if (lGeom == true && File.Exists("albeq.dat") == false) // if possible, read albeq.dat
            {

                Parallel.For(1, NI + 1, Program.pOptions, i =>
                {
                    for (int j = 1; j <= NJ; j++)
                    {
                        Q[i][j] = 0;
                        Alfa[i][j] = 0;
                        Beta[i][j] = 0;
                    }
                });

                Parallel.For(2, NI, Program.pOptions, i =>
                {
                    for (int j = 2; j <= NJ - 1; j++)
                    {
                        double XP, XW, XE, YP, YS, YN, ZE, ZN, ZW, ZS;
                        XP = Program.X[i];
                        XW = Program.X[i - 1];
                        XE = Program.X[i + 1];
                        YP = Program.Y[j];
                        YS = Program.Y[j - 1];
                        YN = Program.Y[j + 1];
                        ZW = Program.AH[i - 1][j];
                        ZE = Program.AH[i + 1][j];
                        ZS = Program.AH[i][j - 1];
                        ZN = Program.AH[i][j + 1];

                        //Inclination and exposition
                        float HX = (float)((ZE - ZW) / (XE - XW)); // inclination gradients of the ground
                        float HY = (float)((ZN - ZS) / (YN - YS));

                        Alfa[i][j] = MathF.Atan(MathF.Sqrt(HX * HX + HY * HY)); //inclination of the ground cell

                        if (Math.Abs(Alfa[i][j]) < Eps)
                            Beta[i][j] = 0; // direction of the ground cell 0 = north
                        else
                            Beta[i][j] = MathF.Atan2(HX, HY);

                        Q[i][j] = CLQ(i, j, Alfa[i][j], Beta[i][j]);  //sky view factor  = 0 in case of flat terrain (1 - Q = horizontal elevation)
                    }
                });

                if (NK != 3) Console.Write(" / geometry computed");
                //save geometry data
                if (lGeomWr == true)
                {
                    try
                    {
                        using (BinaryWriter wr = new BinaryWriter(File.Open("albeq.dat", FileMode.Create)))
                        {
                            for (int i = 1; i <= Program.NX; i++)
                                for (int j = 1; j <= Program.NY; j++)
                                {
                                    {
                                        wr.Write((double)Alfa[i][j]); // write as double to keep compatibility to version 20.01 and earlier
                                        wr.Write((double)Beta[i][j]);
                                        wr.Write(Q[i][j]);
                                    }
                                }
                            wr.Write(KStMax);
                        }
                        if (NK != 3) Console.Write(" / geometry saved");
                    }
                    catch
                    {
                        Console.WriteLine("Unable to write geometry data of the radiation model to file 'albeq.dat'");
                        File.Delete("albeq.dat");
                    }
                }
            }


            //read geometry data
            if (lGeomRe == true || File.Exists("albeq.dat"))
            {
                try
                {
                    using (BinaryReader br = new BinaryReader(File.Open("albeq.dat", FileMode.Open, FileAccess.Read, FileShare.Read)))
                    {
                        for (int i = 1; i <= Program.NX; i++)
                            for (int j = 1; j <= Program.NY; j++)
                            {
                                {
                                    Alfa[i][j] = (float)br.ReadDouble();
                                    Beta[i][j] = (float)br.ReadDouble();
                                    Q[i][j] = br.ReadDouble();
                                }
                            }
                        KStMax = br.ReadInt32();
                    }
                    if (NK != 3) Console.Write(" / geometry loaded");
                }
                catch
                {
                    Console.WriteLine("Unable to read geometry data of the radiation model from file 'albeq.dat'");
                }
            }

            //sun azimuth, angle, and height
            ClMyEt();  // calculates eta and He, a0, CosH, CosHe, My, Mye, GPhi

            //surface albedo
            ClAg();

            //vertical integral values
            ClYrWc();

            //general height dependent values
            CLTStern(); // uses a0

            double DOmega, Omega;
            DOmega = 2 * Math.PI / NumberOmega;
            iTheta = (int)((Eta + Math.PI) / (2 * Math.PI / NumberOmega)) + 1;     //sector where eta lies
            if (iTheta > NumberOmega) iTheta = 1;
            Omega = DOmega * (iTheta - 1) + DOmega * 0.5;                               //actual angle

            //Console.WriteLine("Itime in hours" + (timeR/3600).ToString() + " Omega" + Omega.ToString() + " Eta "+ Eta.ToString());

            //main loop for computing solar radiation
            Parallel.For(2, NI, Program.pOptions, i =>
            {
                double SinIe = 0;                            //sun angle above surface (Mye) [rad]
                double SinI = 0;
                double ThetaB = 0;
                double ThetaOld = 0;
                double Theta = 0;
                //            for (int i = 2; i <= Program.NX - 1; i++)
                //            {
                for (int j = 2; j <= NJ - 1; j++)
                {
                    // sun angle over the surface
                    SinIe = Mye * MathF.Cos(Alfa[i][j]) + CosHe * MathF.Sin(Alfa[i][j]) * MathF.Cos(Eta - Beta[i][j]);
                    SinI = My * MathF.Cos(Alfa[i][j]) + CosH * MathF.Sin(Alfa[i][j]) * MathF.Cos(Eta - Beta[i][j]);

                    ThetaB = CLThet(i, j, Program.AH[i][j], Omega, DOmega, ThetaB);
                    ThetaOld = -90 * GrRa;

                    for (int k = Program.KST[i][j] - 1; k <= Program.NZ; k++)
                    {
                        ClAlbedo(i, j, k, Program.KST[i][j] - 1, ThetaB, SinI, SinIe);
                        if (ThetaOld >= He)
                            Theta = CLThet(i, j, Program.ZZ[k], Omega, DOmega, Theta);
                        else
                            Theta = -90 * GrRa;

                        ThetaOld = Theta;

                        if (Program.ISOL == 1) // thin clouds
                        {
                            ClSol_1(i, j, k, Theta);
                        }
                        else // thick clouds
                        {
                            ClSol_2(i, j, k, Theta);
                        }
                    }
                    RSolGround(i, j, Program.KST[i][j] - 1, ThetaB, SinIe);
                    dTdt_Sol(i, j);
                }
            });

            if (NK != 3) Console.WriteLine(" / radiation computed");
            //terrestrial downward longwave radiation
            Cl_LStrich();

            //Epsilon Above
            for (int k = 2; k <= NK; k++)
            {
                for (int i = 1; i <= NI; i++)
                {
                    for (int j = 1; j <= NJ; j++)
                    {
                        Cl_EpsA(i, j, k, ref EpsAp[i][j][k], ref EpsAm[i][j][k], ref Eab[i][j][k], ref Tab[i][j][k]);
                    }
                }
            }

            //Terrestrial radiation
            for (int i = 1; i <= NI; i++)
            {
                for (int j = 1; j <= NJ; j++)
                {
                    Cl_RterrG(i, j);
                    Cl_dTdtterr(i, j);
                }
            }

            //border values          
            for (int k = 1; k <= NK; k++)
            {
                Parallel.For(2, NI, Program.pOptions, i =>
                {
                    //                for (int i = 2; i <= Program.NX - 1; i++)
                    //                {
                    Arad[i][1][k] = Arad[i][2][k];
                    SG[i][1][k] = SG[i][2][k];
                    DG[i][1][k] = DG[i][2][k];
                    EG[i][1][k] = EG[i][2][k];
                    Program.DT_SOL[i][1][k] = Program.DT_SOL[i][2][k];
                    Program.DT_TERR[i][1][k] = Program.DT_TERR[i][2][k];

                    Arad[i][Program.NY][k] = Arad[i][Program.NY - 1][k];
                    SG[i][Program.NY][k] = SG[i][Program.NY - 1][k];
                    DG[i][Program.NY][k] = DG[i][Program.NY - 1][k];
                    EG[i][Program.NY][k] = EG[i][Program.NY - 1][k];
                    Program.DT_SOL[i][Program.NY][k] = Program.DT_SOL[i][Program.NY - 1][k];
                    Program.DT_TERR[i][Program.NY][k] = Program.DT_TERR[i][Program.NY - 1][k];
                });
                Parallel.For(2, NJ, Program.pOptions, j =>
                {
                    //                for (int j = 2; j <= Program.NY - 1; j++)
                    //                {
                    Arad[1][j][k] = Arad[2][j][k];
                    SG[1][j][k] = SG[2][j][k];
                    DG[1][j][k] = DG[2][j][k];
                    EG[1][j][k] = EG[2][j][k];
                    Program.DT_SOL[1][j][k] = Program.DT_SOL[2][j][k];
                    Program.DT_TERR[1][j][k] = Program.DT_TERR[2][j][k];

                    Arad[Program.NX][j][k] = Arad[Program.NX - 1][j][k];
                    SG[Program.NX][j][k] = SG[Program.NX - 1][j][k];
                    DG[Program.NX][j][k] = DG[Program.NX - 1][j][k];
                    EG[Program.NX][j][k] = EG[Program.NX - 1][j][k];
                    Program.DT_SOL[Program.NX][j][k] = Program.DT_SOL[Program.NX - 1][j][k];
                    Program.DT_TERR[Program.NX][j][k] = Program.DT_TERR[Program.NX - 1][j][k];
                });
            }
        }

        ///<summary>
        ///Number of days since 1.1.1600
        ///</summary>
        private int DateNum(string Date)
        {
            int DateNum = 0;
            int Day, Month, Year;
            Day = Convert.ToInt32(Date.Substring(0, 2));
            Month = Convert.ToInt32(Date.Substring(2, 2));
            Year = Convert.ToInt32(Date.Substring(4, 4));
            //            if (Month > 2)
            //            {
            //                Month -= 3;
            //            }
            //            else
            //            { 
            //                Month += 9;
            //                Year--;
            //            }
            //            Year -= 1600;
            //            DateNum = Convert.ToInt32((Convert.ToInt32((Year / 100)) * 146097) / 4)
            //                    + Convert.ToInt32(((Year % 100) * 1461) / 4)
            //                    + Convert.ToInt32(((Month * 153) + 2) / 5) + Day + 59;
            //
            DateTime firstdate = new DateTime(1600, 1, 1);
            DateTime seconddate = new DateTime(Year, Month, Day);
            TimeSpan ts = seconddate - firstdate;
            DateNum = ts.Days; // calculates the number of days between two dates
            return DateNum;
        }

        ///<summary>
        ///Calculate the projection to the ground cell
        ///</summary>
        private double CLQ(int i, int j, double Alfai, double Betai)
        {
            //Computation of the projection to the ground cell
            // Alfai = inclination of the ground cell
            // Betai azimuth of the ground cell
            float Omega, DOmega;
            //Projection of surface cell
            DOmega = 2F * MathF.PI / (float)NumberOmega;
            double QI = 0;
            double Theta = 0;
            for (int N = 1; N <= NumberOmega; N++)
            {
                Omega = DOmega * (N - 1) + DOmega * 0.5F;                                    //compute sectors
                float AlfEff = -Alfa[i][j] * MathF.Cos(Beta[i][j] - Omega);                 //actual angle
                                                                                            //Compute horizons of surroundings
                Theta = CLThet(i, j, Program.AH[i][j], Omega, DOmega, Theta);
                float TheS = (float)Math.Max(Theta, AlfEff);
                double DQi = MathF.Cos(AlfEff) * 0.5F * (MathF.Pow(MathF.Sin(TheS), 2) - MathF.Pow(MathF.Sin(AlfEff), 2))
                    - MathF.Sin(AlfEff) * 0.5F * ((MathF.Sin(2 * TheS) - MathF.Sin(2 * AlfEff)) * 0.5
                    + TheS - AlfEff);
                QI += DQi;                                                                               //sum up projection
            }

            QI *= DOmega / Math.PI;

            return QI;
        }

        ///<summary>
        ///search horizon for radiation or shadow estimation
        ///</summary>
        private double CLThet(int I, int J, float Hoehe, double Omega, double DOmega, double Theta)
        {
            // search horizon for radiation or shadow estimation
            int IAkt, JAkt, IZ, JZ, Nenner;
            double Dx1, Dx2, Dy1, Dy2,ThetaI;
            float  DXrad, DYrad, XSum, YSum, ZSum;
            bool Ende;
            float OM1 = (float)(Omega - DOmega * 0.5);
            float OM2 = (float)(Omega + DOmega * 0.5);
            Theta = 0;

            if ((Omega > Math.PI * 1.75) || (Omega <= Math.PI * 0.25)) // northern sector
            {
                IAkt = I;             //angle in the north sector
                JZ = J + 1;
                Ende = false;

                while ((JZ <= Program.NY) && (Ende != true))
                {
                    DYrad = Program.Y[JZ] - Program.Y[J];   //distance JZ-J
                    Dx1 = MathF.Tan(OM1) * DYrad;            //distance in x-directions of OM1
                    Dx2 = MathF.Tan(OM2) * DYrad;            //distance in x-directions of OM2
                    if (((Program.X[I] + Dx2) >= Program.X[Program.NX]) || (Program.X[I] + Dx1) <= Program.X[1])
                    {
                        Ende = true; 						// sector outside the computation area
                        Theta = Math.Max(Theta, 0);
                    }
                    else
                    {
                        if (Program.X[IAkt] < (Program.X[I] + Dx1)) // counter less OM1
                        {
                            while (Program.X[IAkt] < (Program.X[I] + Dx1))
                                IAkt++;
                            IZ = IAkt;
                        }
                        else // counter equal or larger OM1
                        {
                            while (Program.X[IAkt] >= (Program.X[I] + Dx1))
                                IAkt--;
                            IZ = IAkt + 1;
                            IAkt++;
                        }

                        if (Program.X[IZ] > (Program.X[I] + Dx2))     //no point between OM1 and OM2
                        {
                            XSum = (Program.X[IZ] + Program.X[IZ - 1]) * 0.5F - Program.X[I];
                            ZSum = (Program.AH[IZ][JZ] + Program.AH[IZ - 1][JZ]) * 0.5F - Hoehe;
                        }
                        else                                          //numbers up to OM2
                        {
                            XSum = 0;
                            ZSum = 0;
                            Nenner = 0;
                            while (Program.X[IZ] <= (Program.X[I] + Dx2))
                            {
                                Nenner++;
                                XSum += Program.X[IZ] - Program.X[I];
                                ZSum += Program.AH[IZ][JZ] - Hoehe;
                                IZ++;
                            }
                            XSum /= Nenner;
                            ZSum /= Nenner;
                        }
                        ThetaI = MathF.Atan2(ZSum, MathF.Sqrt(XSum * XSum + DYrad * DYrad));
                        Theta = Math.Max(ThetaI, Theta);
                    }
                    JZ++;
                }
            }
            else if ((Omega > Math.PI * 0.25) && (Omega <= Math.PI * 0.75)) // eastern sector
            {
                JAkt = J;
                IZ = I + 1;
                Ende = false;
                while ((IZ <= Program.NX) && (Ende == false))
                {
                    DXrad = Program.X[IZ] - Program.X[I];
                    Dy1 = MathF.Tan(MathF.PI * 0.5F - OM1) * DXrad;
                    Dy2 = MathF.Tan(MathF.PI * 0.5F - OM2) * DXrad;
                    if (((Program.Y[J] + Dy2) <= Program.Y[1]) || (Program.Y[J] + Dy1) >= Program.Y[Program.NY])
                    {
                        Ende = true;
                        Theta = Math.Max(Theta, 0);
                    }
                    else
                    {
                        if (Program.Y[JAkt] > (Program.Y[J] + Dy1))
                        {
                            while (Program.Y[JAkt] > (Program.Y[J] + Dy1))
                                JAkt--;
                            JZ = JAkt;
                        }
                        else
                        {
                            while (Program.Y[JAkt] <= (Program.Y[J] + Dy1))
                                JAkt++;
                            JZ = JAkt - 1;
                            JAkt--;
                        }
                        if (Program.Y[JZ] < (Program.Y[J] + Dy2))
                        {
                            YSum = (Program.Y[JZ] + Program.Y[JZ + 1]) * 0.5F - Program.Y[J];
                            ZSum = (Program.AH[IZ][JZ] + Program.AH[IZ][JZ + 1]) * 0.5F - Hoehe;
                        }
                        else
                        {
                            YSum = 0;
                            ZSum = 0;
                            Nenner = 0;
                            while (Program.Y[JZ] >= (Program.Y[J] + Dy2))
                            {
                                Nenner++;
                                YSum += Program.Y[JZ] - Program.Y[J];
                                ZSum += Program.AH[IZ][JZ] - Hoehe;
                                JZ--;
                            }
                            YSum /= Nenner;
                            ZSum /= Nenner;
                        }
                        ThetaI = MathF.Atan2(ZSum, MathF.Sqrt(YSum * YSum + DXrad * DXrad));
                        Theta = Math.Max(ThetaI, Theta);
                    }
                    IZ++;
                }
            }
            else if ((Omega > Math.PI * 0.75) && (Omega <= Math.PI * 1.25))
            {
                IAkt = I;                                           //angle in the southern sector
                JZ = J - 1;
                Ende = false;
                while ((JZ >= 1) && (Ende == false))
                {
                    DYrad = Program.Y[J] - Program.Y[JZ];
                    Dx1 = MathF.Tan(MathF.PI - OM1) * DYrad;
                    Dx2 = MathF.Tan(MathF.PI - OM2) * DYrad;
                    if (((Program.X[I] + Dx2) <= Program.X[1]) || (Program.X[I] + Dx1) >= Program.X[Program.NX])
                    {
                        Ende = true;
                        Theta = Math.Max(Theta, 0);
                    }
                    else
                    {
                        if (Program.X[IAkt] > (Program.X[I] + Dx1))
                        {
                            while (Program.X[IAkt] > (Program.X[I] + Dx1))
                                IAkt--;
                            IZ = IAkt;
                        }
                        else
                        {
                            while (Program.X[IAkt] <= (Program.X[I] + Dx1))
                                IAkt++;
                            IZ = IAkt - 1;
                            IAkt--;
                        }
                        if (Program.X[IZ] < (Program.X[I] + Dx2))
                        {
                            XSum = (Program.X[IZ] + Program.X[IZ + 1]) * 0.5F - Program.X[I];
                            ZSum = (Program.AH[IZ][JZ] + Program.AH[IZ + 1][JZ]) * 0.5F - Hoehe;
                        }
                        else
                        {
                            XSum = 0;
                            ZSum = 0;
                            Nenner = 0;
                            while (Program.X[IZ] >= (Program.X[I] + Dx2))
                            {
                                Nenner++;
                                XSum += Program.X[IZ] - Program.X[I];
                                ZSum += Program.AH[IZ][JZ] - Hoehe;
                                IZ--;
                            }
                            XSum /= Nenner;
                            ZSum /= Nenner;
                        }
                        ThetaI = MathF.Atan2(ZSum, MathF.Sqrt(XSum * XSum + DYrad * DYrad));
                        Theta = Math.Max(ThetaI, Theta);
                    }
                    JZ--;
                }
            }
            else if ((Omega > Math.PI * 1.25) && (Omega <= Math.PI * 1.75))
            {
                JAkt = J;                                            //angle in the western sector
                IZ = I - 1;
                Ende = false;
                while ((IZ >= 1) && (Ende == false))
                {
                    DXrad = Program.X[I] - Program.X[IZ];
                    Dy1 = -MathF.Tan(MathF.PI * 1.5F - OM1) * DXrad;
                    Dy2 = -MathF.Tan(MathF.PI * 1.5F - OM2) * DXrad;
                    if (((Program.Y[J] + Dy2) >= Program.Y[Program.NY]) || (Program.Y[J] + Dy1) <= Program.Y[1])
                    {
                        Ende = true;
                        Theta = Math.Max(Theta, 0);
                    }
                    else
                    {
                        if (Program.Y[JAkt] < (Program.Y[J] + Dy1))
                        {
                            while (Program.Y[JAkt] < (Program.Y[J] + Dy1))
                                JAkt++;
                            JZ = JAkt;
                        }
                        else
                        {
                            while (Program.Y[JAkt] >= (Program.Y[J] + Dy1))
                                JAkt--;
                            JZ = JAkt + 1;
                            JAkt++;
                        }
                        if (Program.Y[JZ] > (Program.Y[J] + Dy2))
                        {
                            YSum = (Program.Y[JZ] + Program.Y[JZ - 1]) * 0.5F - Program.Y[J];
                            ZSum = (Program.AH[IZ][JZ] + Program.AH[IZ][JZ - 1]) * 0.5F - Hoehe;
                        }
                        else
                        {
                            YSum = 0;
                            ZSum = 0;
                            Nenner = 0;
                            while (Program.Y[JZ] <= (Program.Y[J] + Dy2))
                            {
                                Nenner++;
                                YSum += Program.Y[JZ] - Program.Y[J];
                                ZSum += Program.AH[IZ][JZ] - Hoehe;
                                JZ++;
                            }
                            YSum /= Nenner;
                            ZSum /= Nenner;
                        }
                        ThetaI = MathF.Atan2(ZSum, MathF.Sqrt(YSum * YSum + DXrad * DXrad));
                        Theta = Math.Max(ThetaI, Theta);
                    }
                    IZ--;
                }
            }
            else
            {
                Console.WriteLine("Error when computing Theta in the radiation model");
                Environment.Exit(0);
            }

            return Theta;
        }

        ///<summary>
        ///Sun height, angle, and azimuth
        ///</summary>
        private double ClMyEt()
        {
            float PhiThe;
            float Delta;       //sun declination [rad]
            float Psirad;      //sun-hour angle [rad]
            double CosEta;      //cosinus of sun azimuth [rad]

            PhiThe = 2 * MathF.PI * ((float)Na - 2.84F) / 365;
            GPhi = 1.0006F + 0.03343F * MathF.Cos(PhiThe) + 0.0011F * MathF.Sin(PhiThe);   //radiation flux variation
            Delta = MathF.Asin(MathF.Sin(23.5F * MathF.PI / 180) * MathF.Sin(2 * MathF.PI * (Na - 80) / 365));   //sun declination
            Psirad = (float)(TimeR / 43200 * Math.PI - Math.PI);   // horizontal sun angle in hours 0:00 = -Pi, 12:00 = 0, 24:00 = Pi
            
            float latitudeRadiant = (float) (Program.BGRAD * GrRa);
            Mye = MathF.Sin(latitudeRadiant) * MathF.Sin(Delta) + MathF.Cos(latitudeRadiant) * MathF.Cos(Delta) * MathF.Cos(Psirad);    //sinus of sun angle
            My = MathF.Max(0.02F, Mye);					  // needed to avoid errors if sun under the horizon
            CosHe = MathF.Sqrt(1 - Mye * Mye);            // cosinus of sun angle
            CosH = MathF.Sqrt(1 - My * My);
            He = MathF.Asin(Mye);                         // vertical sun angle
            A0 = 0.04 + 0.065 * MathF.Exp(-0.18F / My);

            if ((MathF.Abs(CosHe) < Eps) || (latitudeRadiant == MathF.PI * 0.5F) || (MathF.Abs(Psirad) < Eps))
            {
                Eta = 0; // at the pole, at noon
            }
            else
            {
                CosEta = (Mye * MathF.Sin(latitudeRadiant) - MathF.Sin(Delta)) / (CosHe * MathF.Cos(latitudeRadiant));
                if (Math.Abs(CosEta) <= 1 + Eps)
                {
                    if (CosEta > 1)
                        CosEta = 1;
                    if (CosEta < -1)
                        CosEta = -1;
                    Eta = (float)Math.Acos(CosEta);
                    if (Psirad < 0)
                        Eta = -Eta;
                }
                else
                {
                    Console.WriteLine("Error when computing Eta in the radiation model (-1<coseta<1 not fulfilled)");
                    Environment.Exit(0);
                }
            }

            if (Program.BGRAD < 0) // on the southern hemisphere - set the sun position to the north
            {
                //Console.WriteLine(Eta.ToString());
                if (Eta < 0 && Eta > -Math.PI / 2)
                {
                    Eta = -MathF.PI - Eta;
                }
                if (Eta > 0 && Eta < Math.PI / 2)
                {
                    Eta = MathF.PI - Eta;
                }

                if (Eta > 2 * Math.PI)
                    Eta -= 2 * MathF.PI;

                if (Eta < Eps)
                    Eta = MathF.PI;
            }
            //Console.WriteLine("ETA " + Eta.ToString());
            return He;
        }

        ///<summary>
        ///Surface albedo
        ///</summary>
        private double ClAg()
        {
            Parallel.For(1, Program.NX + 1, Program.pOptions, I =>
            {
                //            for (int I = 1; I <= Program.NX; I++)
                //            {
                for (int J = 1; J <= Program.NY; J++)
                {
                    Ag[I][J] = Program.ALBEDO[I][J];
                    if (Ag[I][J] == 0.08)
                        Ag[I][J] = MathF.Min(1, MathF.Max(0.03F, -0.0139F + MathF.Tan(MathF.PI * 0.5F - He)));    //Flassak p 133
                }
            });
            return He;
        }

        ///<summary>
        ///Vertical integrals of haze, etc.
        ///</summary>
        private double ClYrWc()
        {
            //N,               ! Zaehlervariable
            //K                ! Hoehenindex der Zelle
            double[] Rhop = new double[Program.NPROFMAX];      //Density [kg/m³]
            double[] QRp = new double[Program.NPROFMAX];       //Intermediate value for calculating Rp F50a
            double[] QCp = new double[Program.NPROFMAX];       //Intermediate value for calculating Cp F50b
            double[] Vnz = new double[Program.NPROFMAX];       //Intermediate value
            double[] QVapz = new double[Program.NPROFMAX];     //Intermediate value
            double[] QCldz = new double[Program.NPROFMAX];     //Intermediate value
            //double[][][] QCldz = CreateArray<double[][]>(Program.Program.NX1, () => CreateArray<double[]>(Program.NY1, () => new double[Program.NPROFMAX]));     //Intermediate value
            double[] QIcez = new double[Program.NPROFMAX];     //Intermediate value


            for (int N = 1; N <= 30; N++)
            {
                Rhop[N] = Program.PPROF[N] / 287 / Program.TPROF[N];             //general gas equation
                QRp[N] = Rhop[N] * Program.QVAP[N] * Program.PPROF[N] / P0;      //F50a
                QCp[N] = Rhop[N] * QCo2 * Program.PPROF[N] / P0;                 //F50b
                Vnz[N] = 570 * (1 / MathF.Pow((float)Program.VNORM[N], 1.5F));  //F13a
                QVapz[N] = Rhop[N] * Program.QVAP[N];                            //F13b
                QCldz[N] = Rhop[N] * (Program.QCLD[N] + Program.QRAIN[N]);       //F13c
                QIcez[N] = Rhop[N] * (Program.QICE[N] + Program.QSNOW[N]);       //F13d

                /*
                Program.ZPROF[1] = Math.Max(Program.AHMIN, 0);
                double UNTO;
                double UNTU;
                double DIFF;
                Int32 INDO = 0;
                Int32 INDU = 0;
                for (int i = 1; i <= Program.NX; i++)
                {
                    for (int j = 1; j <= Program.NY; j++)
                    {
                        if (N > 1) Program.ZPROF[N] = Program.ZPROF[N - 1] + 500;
                        if (Program.ZPROF[N] > Program.ZSP[i][j][Program.NZ])
                        {
                            QCLD[N] = 0.0;
                        }
                        else if (Program.ZPROF[N] <= Program.ZSP[i][j][Program.NZ])
                        {
                            //temperature and pressure profile within GRAMM domain
                            UNTU = -100000;
                            UNTO = 100000;
                            INDO = 0;
                            INDU = 0;
                            double GRAD = 0;
                            for (int m = 1; m <= Program.NZ; m++)
                            {
                                DIFF = Program.ZSP[i][j][m] - Program.ZPROF[N];
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
                                GRAD = (Program.WAT_VAP[i][j][INDU] - Program.WAT_VAP[i][j][INDU - 1]) /
                                    (Program.ZSP[i][j][INDU] - Program.ZSP[i][j][INDU - 1]);
                                QCLD[N] = Program.WAT_VAP[i][j][INDU] + (Program.ZPROF[N] - Program.ZSP[i][j][INDU]) * GRAD;
                            }
                            else if (INDU == 0)
                            {
                                GRAD = (Program.WAT_VAP[i][j][INDO] - Program.WAT_VAP[i][j][INDO + 1]) /
                                    (Program.ZSP[i][j][INDO] - Program.ZSP[i][j][INDO + 1]);
                                QCLD[N] = Program.WAT_VAP[i][j][INDO] + (Program.ZPROF[N] - Program.ZSP[i][j][INDO]) * GRAD;
                            }
                            else
                            {
                                GRAD = (Program.WAT_VAP[i][j][INDO] - Program.WAT_VAP[i][j][INDU]) /
                                    (Program.ZSP[i][j][INDO] - Program.ZSP[i][j][INDU]);
                                QCLD[N] = Program.WAT_VAP[i][j][INDO] + (Program.ZPROF[N] - Program.ZSP[i][j][INDO]) * GRAD;
                            }
                        }
                        QCldz[i][j][N] = 0 * Rhop[N] * (Math.Max(QCLD[N], 0) * 0.001 + QRAIN[N]);     //F13c
                    }
                }
                */
            }

            // it seems, that there is a dependency, so parallelisation is not possible!
            for (int K = 1; K <= Program.NZ; K++) // start of integration
            {
                if (K == 1) // Border of the 1st cell
                {
                    Dz1[K] = 0;                                          //borders of first cell
                    Dz2[K] = (Program.ZZ[K] + Program.ZZ[K + 1]) * 0.5;
                }
                else if (K == Program.NZ) // border of the last cell k
                {
                    Dz1[K] = (Program.ZZ[K] + Program.ZZ[K - 1]) * 0.5;    //borders of last cell
                    Dz2[K] = Program.ZZ[K] + (Program.ZZ[K] - Program.ZZ[K - 1]) * 0.5;
                }
                else
                {
                    Dz1[K] = (Program.ZZ[K] + Program.ZZ[K - 1]) * 0.5;    //borders of the k-th cell
                    Dz2[K] = (Program.ZZ[K] + Program.ZZ[K + 1]) * 0.5;
                }
                Program.TABS[1][1][K] = rInteg(Program.ZPROF, Program.TPROF, Program.NPROFMAX, Dz1[K], Dz2[K]) / (Dz2[K] - Dz1[K]);   //temperature in k-th cell
                Program.PBZZ[K] = rInteg(Program.ZPROF, Program.PPROF, Program.NPROFMAX, Dz1[K], Dz2[K]) / (Dz2[K] - Dz1[K]);   //pressure in k-th cell

                Program.RHOBZZ[K] = Program.PBZZ[K] / 287 / Program.TABS[1][1][K]; // gas equation
                YG[K][3] = rInteg(Program.ZPROF, Vnz, Program.NPROFMAX, Program.ZZ[K], Program.ZPROF[30]);      // valus direction up 
                R[K][3] = rInteg(Program.ZPROF, QVapz, Program.NPROFMAX, Program.ZZ[K], Program.ZPROF[30]);
                W2rad[K][3] = rInteg(Program.ZPROF, QIcez, Program.NPROFMAX, Program.ZZ[K], Program.ZPROF[30]);
                Rp[K][3] = rInteg(Program.ZPROF, QRp, Program.NPROFMAX, Program.ZZ[K], Program.ZPROF[30]);
                CGp[K][3] = rInteg(Program.ZPROF, QCp, Program.NPROFMAX, Program.ZZ[K], Program.ZPROF[30]);

                YG[K][1] = YG[1][3] - YG[K][3];        // values directon down
                R[K][1] = R[1][3] - R[K][3];
                W2rad[K][1] = W2rad[1][3] - W2rad[K][3];
                Rp[K][1] = Rp[1][3] - Rp[K][3];
                CGp[K][1] = CGp[1][3] - CGp[K][3];

                YG[K][2] = rInteg(Program.ZPROF, Vnz, Program.NPROFMAX, Dz1[K], Dz2[K]); // value at the cell
                R[K][2] = rInteg(Program.ZPROF, QVapz, Program.NPROFMAX, Dz1[K], Dz2[K]);
                W2rad[K][2] = rInteg(Program.ZPROF, QIcez, Program.NPROFMAX, Dz1[K], Dz2[K]);
                Rp[K][2] = rInteg(Program.ZPROF, QRp, Program.NPROFMAX, Dz1[K], Dz2[K]);
                CGp[K][2] = rInteg(Program.ZPROF, QCp, Program.NPROFMAX, Dz1[K], Dz2[K]);

                W1rad[K][3] = rInteg(Program.ZPROF, QCldz, Program.NPROFMAX, Program.ZZ[K], Program.ZPROF[30]);
                W1rad[K][1] = W1rad[1][3] - W1rad[K][3];
                W1rad[K][2] = rInteg(Program.ZPROF, QCldz, Program.NPROFMAX, Dz1[K], Dz2[K]);
                Wrad[K][3] = W1rad[K][3] + 0.5F * W2rad[K][3]; // F13e
                Wrad[K][1] = W1rad[K][1] + 0.5F * W2rad[K][1];
                Wrad[K][2] = W1rad[K][2] + 0.5F * W2rad[K][2];

                /*
                for (int i = 1; i <= Program.NX; i++)
                {
                    for (int j = 1; j <= Program.NY; j++)
                    {
                        W1rad[i][j][K][3] = rInteg(Program.ZPROF, QCldz, Program.NPROFMAX, Program.ZZ[K], Program.ZPROF[30]);
                        W1rad[i][j][K][1] = W1rad[i][j][1][3] - W1rad[i][j][K][3];
                        W1rad[i][j][K][2] = rInteg(Program.ZPROF, QCldz, Program.NPROFMAX, Dz1[K], Dz2[K]);

                        Wrad[i][j][K][3] = W1rad[i][j][K][3] + 0.5F * W2rad[K][3]; // F13e
                        Wrad[i][j][K][1] = W1rad[i][j][K][1] + 0.5F * W2rad[K][1];
                        Wrad[i][j][K][2] = W1rad[i][j][K][2] + 0.5F * W2rad[K][2];
                    }
                }
                */
            }

            return 0d;
        }

        ///<summary>
        ///integration over a field of points
        ///</summary>
        private double rInteg(double[] Xvalue, double[] Yvalue, int nMax, double X1, double X2)
        {
            double rInteg = 0;
            double Y = 0;
            int N = 2;
            while (X1 > Xvalue[N])
            {
                N++;
            }
            if (X2 < Xvalue[N])
            {
                Y = Yvalue[N - 1] + ((X1 + X2) * 0.5F - Xvalue[N - 1]) * ((Yvalue[N] - Yvalue[N - 1]) / (Xvalue[N] - Xvalue[N - 1]));
                rInteg = Y * (X2 - X1);
            }
            else
            {
                Y = Yvalue[N - 1] + ((X1 + Xvalue[N]) * 0.5F - Xvalue[N - 1]) * ((Yvalue[N] - Yvalue[N - 1]) / (Xvalue[N] - Xvalue[N - 1]));
                rInteg = Y * (Xvalue[N] - X1);
                N++;
                while (Xvalue[N] < X2)
                {
                    Y = Yvalue[N - 1] + ((Xvalue[N] + Xvalue[N - 1]) * 0.5F - Xvalue[N - 1]) * ((Yvalue[N] - Yvalue[N - 1]) / (Xvalue[N] - Xvalue[N - 1]));
                    rInteg += Y * (Xvalue[N] - Xvalue[N - 1]);
                    N++;
                }
                Y = Yvalue[N - 1] + ((X2 + Xvalue[N - 1]) * 0.5F - Xvalue[N - 1]) * ((Yvalue[N] - Yvalue[N - 1]) / (Xvalue[N] - Xvalue[N - 1]));
                rInteg += Y * (X2 - Xvalue[N - 1]);
            }

            return rInteg;
        }

        ///<summary>
        ///General height dependent values
        ///</summary>
        private double CLTStern()
        {
            double Tau1, Tau2, TauS;

            for (int k = 1; k <= Program.NZ; k++)
            {
                Tstern[k] = 0.3 + (MathF.Pow(0.1F * (float)R[k][3], 0.2F)) * (0.75 + My) +
                    (8.0 * P0 / Program.PBZZ[k]) * (YG[k][3] + 0.1F * Program.PBZZ[k] / P0 + 0.02);
                IG[k] = GPhi * IG0 * Program.Solar_reduction_factor * MathF.Exp((float)(-1 * ((Tstern[k] / My) * A0)));

                float tmp = 1000F * (float)W1rad[k][3] + 1;
                if (MathF.Log(tmp) == 0)
                {
                    Tau1 = 0;
                }
                else
                {
                    Tau1 = 0.07 * MathF.Pow(MathF.Log(tmp), 3.95F);
                }
                Tau2 = 56.6 * W2rad[k][3];
                Tau[k] = Tau1 + Tau2;
                TauS = Math.Min(32, Tau[k]);
                CloudES[k] = MathF.Pow((float)(1 - TauS / 32), 4);
                CloudESeS[k] = Math.Max(0.1f, CloudES[k]);

                Ns[k] = 0.346 * (CloudESeS[k] - 0.1) + 0.85 * MathF.Pow(CloudESeS[k] - 0.1F, 2);
                Nd[k] = 0.3 * MathF.Sin(MathF.PI * (CloudES[k] + 0.25F)) + 0.36F * CloudES[k] - 0.15F;
                /*
                for (int i = 1; i <= Program.NX; i++)
                {
                    for (int j = 1; j <= Program.NY; j++)
                    {
                        Tau2 = 56.6 * W2rad[k][3];
                        Tau[i][j][k] = Tau1 + Tau2;
                        TauS = Math.Min(32, Tau[i][j][k]);
                        Program.CloudeS[i][j][k] = Math.Pow(1 - TauS / 32, 4);
                        Program.CloudeSeS[i][j][k] = Math.Max(0.1, Program.CloudeS[i][j][k]);

                        nS[i][j][k] = 0.346 * (Program.CloudeSeS[i][j][k] - 0.1) + 0.85 * Math.Pow(Program.CloudeSeS[i][j][k] - 0.1, 2);
                        nD[i][j][k] = 0.3 * Math.Sin(Math.PI * (Program.CloudeS[i][j][k] + 0.25)) + 0.36 * Program.CloudeS[i][j][k] - 0.15;
                    }
                }
                */
            }
            return 0d;
        }

        ///<summary>
        ///Computing albedo in each cell
        ///</summary>
        private double ClAlbedo(int i, int j, int k, int KSti, double Theta, double SinI, double SinIe)
        {
            /*
             * i,j,k		  : cell number 
             * KSti          E: cell-number after which no surface exists
               WUnt           : integrated water downward
               TClearUP,..    : transmission-coefficients
               Ac             : cloud-albedo
               AgMod          : modified surface albedo
               Theta          : excess angle of surface cell
               AcZw           : intermediate value for computing Ac
            */
            // negative values are possible, if Mye is used in the following formulasand the sun is under the horizon (t=20:00 Phi = 47°)
            // a square root of a negative number willcause an error
            // therefore My is used and limited by 0.02

            double TClearUP, TClearDown, Ac, AgMod;
            float AcZw, WUnt, TCUp, TCDown;

            Myx[k] = My * CloudES[k] + 0.5F * (1 - CloudES[k]);
            WUnt = (float) (Wrad[KSti][3] - Wrad[k][3]);
            TClearUP = Math.Pow(Program.PBZZ[k] / Program.PBZZ[KSti], 0.25);
            TClearDown = TClearUP;
            TCDown = 1.2F * MathF.Pow(Myx[k] + 0.2F, 0.5F) * (1.3F * MathF.Pow(1000F * WUnt + 0.1F, -0.33F) - 0.06F);
            TCDown = MathF.Min(1, TCDown);
            TCDown = MathF.Max(0, TCDown);
            TCUp = 1.3F * MathF.Pow(1000F * WUnt + 0.1F, -0.33F) - 0.06F;
            TCUp = MathF.Min(1, TCUp);
            TCUp = MathF.Max(0, TCUp);
            AcZw = 2600 * WUnt * (1.2F - Myx[k]);

            if (AcZw <= Eulerian)
                Ac = 0;
            else
                Ac = 0.42 * MathF.Log(MathF.Log(AcZw));

            if ((SinIe < 0) || (He < Theta))
                AgMod = Ag[i][j] * (1 - CloudES[KSti]);
            else
                AgMod = Ag[i][j] * (1 - CloudES[KSti] + CloudES[KSti] * SinI / My);

            Arad[i][j][k] = Ac + (TClearUP * TClearDown * TCUp * TCDown * AgMod) / (1 - Ac * AgMod);

            return Ac;
        }

        ///<summary>
        ///Solar scheme 1 (thin clouds)
        ///</summary>
        private double ClSol_1(int i, int j, int k, double Theta)
        {
            /*
             * i,j,k		  : cell number 
             * Tauj           : cloud parameter F32
               fA             : correction factor F31
               jzw            : intermediate value for computing j32 F32
               j32            : F32
               Theta          : excess angle of surface cell
            */
            double fA, jzw, j32;
            float Tauj;
            jzw = 0;

            fA = 1 + (1 - 0.5F * CloudES[k]) * (Arad[i][j][k] - 0.2);

            if (Tau[k] < 2)
            {
                jzw = ((0.07 * Tstern[k] - 0.115) / My) * (Program.PBZZ[k] / P0);
                j32 = Math.Max(0, jzw);
            }
            else
            {
                Tauj = MathF.Max(32, (float)Tau[k]);
                j32 = 45.25 * MathF.Pow(Tauj, -1.5F);
            }

            if (He < Theta)
            {
                SG[i][j][k] = 0;
            }
            else
            {
                SG[i][j][k] = IG[k] * Mye * Ns[k] * (Program.PBZZ[k] / P0) * (1 - 0.233 * Program.CLOUDS[i][j] - 0.415 * Program.Pow2(Program.CLOUDS[i][j]));
            }

            DG[i][j][k] = IG[k] * My * (Nd[k] + j32) * fA * (Program.PBZZ[k] / P0) * (1 - 0.233 * Program.CLOUDS[i][j] - 0.415 * Program.Pow2(Program.CLOUDS[i][j]));
            EG[i][j][k] = SG[i][j][k] + DG[i][j][k];

            return jzw;
        }

        ///<summary>
        ///Solar scheme 2 (thick clouds)
        ///</summary>
        private double ClSol_2(int i, int j, int k, double Theta)
        {
            /*
             * fAE            : correction factor F39
               jzw            : intermediate value for computing j38 F38
               j38            : F38
               Tc             : cloud transmission factor F40
               SStrich        : intermediate value for computing direct solar radiation F36
               DStrich        : intermediate value for computing diffusive solar radiation F37
               Theta          : excess angle of surface cell
            */
            double fAE, jzw, j38, SStrich, DStrich;
            float Tc;
            jzw = 0;

            jzw = ((0.07 * Tstern[k] - 0.115) / My) * (Program.PBZZ[k] / P0);
            j38 = Math.Max(0, jzw);
            fAE = 1 + (1 - 0.9 * CloudES[k]) * (Arad[i][j][k] - 0.2);
            Tc = 1.2F * MathF.Pow(Mye + 0.2F, 0.5F) * (1.3F * MathF.Pow((float)(1000F * Wrad[k][3] + 0.1F), -0.33F) - 0.06F);
            Tc = MathF.Min(1, Tc);
            Tc = MathF.Max(0, Tc);

            if (He < Theta)
            {
                SStrich = 0;
            }
            else
            {
                SStrich = IG[k] * Mye;
            }

            DStrich = IG[k] * My * j38;
            EG[i][j][k] = (SStrich + DStrich) * Tc * fAE;

            if (He < Theta)
            {
                SG[i][j][k] = 0;
            }
            else
            {
                SG[i][j][k] = 0.8 * EG[i][j][k] * (CloudESeS[k] - 0.1) / 0.9;
            }

            DG[i][j][k] = EG[i][j][k] - SG[i][j][k];

            return jzw;
        }

        ///<summary>
        ///Solar radiation at the surface
        ///</summary>
        private double RSolGround(int i, int j, int k, double Theta, double SinIe)
        {
            /*
               Theta          : excess angle of surface cell
               XH             : intermediate value
               SGG            : direct solar radiation at the surface
               DGSG           : direct share of diffusive solar radiation at the surface
               DGDG           : isotropic share of diffusive solar radiation at the surface
               BGG            : reflective share of solar radiation at the surface
            */
            double XH, SGG, DGSG, DGDG, BGG;

            XH = 1.5 * My * CloudES[k] * (0.65 + 0.04 * Tstern[k] * Program.PBZZ[k] / P0 - 0.6 * CloudES[k]);

            if ((SinIe < 0) || (He < Theta))
            {
                SGG = 0;
                DGSG = 0;
            }
            else
            {
                SGG = SG[i][j][k] * SinIe / Mye;
                DGSG = XH * SinIe * DG[i][j][k] / My;
            }

            DGDG = (1 - XH) * (1 - Q[i][j]) * DG[i][j][k];
            BGG = Ag[i][j] * Q[i][j] * (SG[i][j][k] + DG[i][j][k]);

            Program.GLOBRAD[i][j] = (SGG + DGSG + DGDG + BGG);
            Program.RSOLG[i][j] = (1 - Ag[i][j]) * Program.GLOBRAD[i][j];

            //influence of snow cover, except in urban areas and water bodies
            if (Program.SNOW[i][j] > 0.15 && Program.SNOW[i][j] > Program.Z0[i][j] && Program.ALAMBDA[i][j] != 4.0 && Program.FW[i][j] < 0.99)
            {
                Program.RSOLG[i][j] = (1 - 0.6) * Program.GLOBRAD[i][j];
            }

            return XH;
        }

        ///<summary>
        ///Heating rate caused by solar radiation
        ///</summary>
        private double dTdt_Sol(int i, int j)
        {
            /*
               k                : height index
               IceDens          : densitiy of ice in air [kg/m³]
               WatDens          : densitiy of water in air [kg/m³]
               Dz               : height of scalar cell [m]
               EzPlus           : total solar radiation at position z+z/2
               EzMin            : total solar radiation at position z-z/2
               AzMin            : albedo at position z-z/2
               AcDown           : absorption function of clouds for downward radiation
               AcUp             : absorption function of clouds for upward radiation
            */
            double IceDens, WatDens, Dz, EzMin, EzPlus, AzMin, AcDown;
            double AcUp = 0;

            for (int k = Program.KST[i][j]; k <= Program.NZ - 1; k++)
            {
                Dz = (Dz2[k] - Dz1[k]);
                IceDens = W2rad[k][2] / Dz;
                WatDens = W1rad[k][2] / Dz;
                AzMin = Arad[i][j][k];

                if (SG[i][j][k] != 0)
                {
                    EzMin = EG[i][j][k] + (Dz1[k] - Program.ZZ[k]) * ((EG[i][j][k + 1] - EG[i][j][k]) / (Program.ZZ[k + 1] - Program.ZZ[k]));
                    EzPlus = EG[i][j][k] + (Dz2[k] - Program.ZZ[k]) * ((EG[i][j][k + 1] - EG[i][j][k]) / (Program.ZZ[k + 1] - Program.ZZ[k]));
                }
                else
                {
                    EzMin = EG[i][j][k - 1] + (Dz1[k] - Program.ZZ[k - 1]) * ((EG[i][j][k] - EG[i][j][k - 1]) / (Program.ZZ[k] - Program.ZZ[k - 1]));
                    EzPlus = EG[i][j][k - 1] + (Dz2[k] - Program.ZZ[k - 1]) * ((EG[i][j][k] - EG[i][j][k - 1]) / (Program.ZZ[k] - Program.ZZ[k - 1]));
                }

                if ((IceDens < 0.000001) && (WatDens < 0.00001))
                {
                    Program.DT_SOL[i][j][k] = 1 / (Program.RHOBZZ[k] * Cp * Dz) * Math.Abs(EzPlus - EzMin) * 0.7 * (1 + AzMin);
                }
                else
                {
                    float tmp = (float)(1 + 1000 * Wrad[k][2]);
                    AcDown = 0.025 * MathF.Pow((Myx[k] + Myx[k + 1]) * 0.5F + 0.1F, 0.5F) * MathF.Log(tmp);
                    AcUp = 0.0194 * MathF.Log(tmp);
                    Program.DT_SOL[i][j][k] = 1 / (Program.RHOBZZ[k] * Cp * Dz) * Math.Abs(EzPlus * AcDown + EzMin * AzMin * AcUp);
                }
            }

            Program.DT_SOL[i][j][Program.NZ] = Program.DT_SOL[i][j][Program.NZ - 1];
            Program.DT_SOL[i][j][Program.KST[i][j] - 1] = Program.DT_SOL[i][j][Program.KST[i][j]];

            return AcUp;
        }

        ///<summary>
        ///Longwave downward terrestrial radiation at specific heights
        ///</summary>
        private double Cl_LStrich()
        {
            /*
               k                : height index
               kk               : summation height index
               WMin             : W-value below kk
               WPlus            : W-value above kk
               CGpMin, CGpPlus, RpMin, RpPlus             : analogu to WPlus
               WPlusG           : W-value at the border kk and kk+1
               WMinG            : W-value at the border kk and kk-1
               CGpPlusG, CGpMinG, RpMinG, RpPlusG         : analogu to WPlus
               EpsPlusG         : Epsilon-value at the upper border
               EpsMin G         : Epsilon-value at the lower border
               W_Last           : Epsilon-value of previous cell
            */
            int kk;
            double WPlus, CGpMin, CGpPlus, RpMin, RpPlus, W_Last;
            double WMin = 0;
            double EpsPlus = 0;
            double EpsMin = 0;

            for (int k = 1; k <= KStMax; k++)
            {
                kk = k;
                LStrich[k] = 0;
                W_Last = 0;
                while ((kk <= Program.NZ) && (W_Last < Wcrit))
                {
                    if (kk == k)
                    {
                        WMin = Wrad[kk][1] + (Dz1[kk] - Program.ZZ[kk + 1]) * ((Wrad[kk + 1][1] - Wrad[kk][1]) / (Program.ZZ[kk] - Program.ZZ[kk + 1]));
                        CGpMin = CGp[kk][1] + (Dz1[kk] - Program.ZZ[kk + 1]) * ((CGp[kk + 1][1] - CGp[kk][1]) / (Program.ZZ[kk] - Program.ZZ[kk + 1]));
                        RpMin = Rp[kk][1] + (Dz1[kk] - Program.ZZ[kk + 1]) * ((Rp[kk + 1][1] - Rp[kk][1]) / (Program.ZZ[kk] - Program.ZZ[kk + 1]));
                    }
                    else
                    {
                        WMin = Wrad[kk - 1][1] + (Dz1[kk] - Program.ZZ[kk - 1]) * ((Wrad[kk][1] - Wrad[kk - 1][1]) / (Program.ZZ[kk] - Program.ZZ[kk - 1]));
                        CGpMin = CGp[kk - 1][1] + (Dz1[kk] - Program.ZZ[kk - 1]) * ((CGp[kk][1] - CGp[kk - 1][1]) / (Program.ZZ[kk] - Program.ZZ[kk - 1]));
                        RpMin = Rp[kk - 1][1] + (Dz1[kk] - Program.ZZ[kk - 1]) * ((Rp[kk][1] - Rp[kk - 1][1]) / (Program.ZZ[kk] - Program.ZZ[kk - 1]));
                    }

                    if (kk == Program.NZ)
                    {
                        WPlus = Wrad[kk - 1][1] + (Dz2[kk] - Program.ZZ[kk - 1]) * ((Wrad[kk][1] - Wrad[kk - 1][1]) / (Program.ZZ[kk] - Program.ZZ[kk - 1]));
                        CGpPlus = CGp[kk - 1][1] + (Dz2[kk] - Program.ZZ[kk - 1]) * ((CGp[kk][1] - CGp[kk - 1][1]) / (Program.ZZ[kk] - Program.ZZ[kk - 1]));
                        RpPlus = Rp[kk - 1][1] + (Dz2[kk] - Program.ZZ[kk - 1]) * ((Rp[kk][1] - Rp[kk - 1][1]) / (Program.ZZ[kk] - Program.ZZ[kk - 1]));
                    }
                    else
                    {
                        WPlus = Wrad[kk][1] + (Dz2[kk] - Program.ZZ[kk]) * ((Wrad[kk + 1][1] - Wrad[kk][1]) / (Program.ZZ[kk + 1] - Program.ZZ[kk]));
                        CGpPlus = CGp[kk][1] + (Dz2[kk] - Program.ZZ[kk]) * ((CGp[kk + 1][1] - CGp[kk][1]) / (Program.ZZ[kk + 1] - Program.ZZ[kk]));
                        RpPlus = Rp[kk][1] + (Dz2[kk] - Program.ZZ[kk]) * ((Rp[kk + 1][1] - Rp[kk][1]) / (Program.ZZ[kk + 1] - Program.ZZ[kk]));
                    }

                    WPlus -= Wrad[k][1];
                    CGpPlus -= CGp[k][1];
                    RpPlus -= Rp[k][1];
                    WMin -= Wrad[k][1];
                    CGpMin -= CGp[k][1];
                    RpMin -= Rp[k][1];

                    EpsPlus = Cl_Eps(EpsPlus, RpPlus, CGpPlus, WPlus);
                    EpsMin = Cl_Eps(EpsMin, RpMin, CGpMin, WMin);

                    LStrich[k] += Program.SIGMA * Math.Pow(Program.TABS[1][1][kk], 4) * (EpsPlus - EpsMin);
                    W_Last = Wrad[kk][2];
                    kk++;
                }
            }

            return WMin;
        }

        ///<summary>
        ///Epsilon computation by Somielski
        ///</summary>
        private double Cl_Eps(double Epsi, double Rpi, double CGpi, double Wi)
        {
            /*
             * Epsilon computation by Somielski 
             * Epsi        A: Epsilon F55
             * Rpi         E: Integral water vapour
             * CGpi        E: Integral CO2
             * Wi          E: Integral water
             * Eps_r        : Epsilon water vapour
             * Eps_CO2      : Epsilon water carbon dioxid
             * Eps_c        : Epsilon clouds
             * Eps_clear    : Intermediate value
            */
            double Eps_r, Eps_CO2, Eps_c, Eps_clear;

            Rpi = (Rpi * (1 / 1000)) * 100;
            float rpi = (float)Rpi;

            //CGpi must not be negative
            CGpi = Math.Max((CGpi * (1 / 1.963)) * 100, 0);

            float tmp = MathF.Log(rpi);
            if (Rpi < 0.0183)
                Eps_r = 0.113 * MathF.Log(1 + 12.63F * rpi);
            else if (tmp < -3)
                Eps_r = 0.104 * MathF.Log(rpi) + 0.44;
            else if (tmp < -1.5)
                Eps_r = 0.121 * MathF.Log(rpi) + 0.491;
            else if (tmp < -1.0)
                Eps_r = 0.146 * MathF.Log(rpi) + 0.527;
            else if (tmp < 0)
                Eps_r = 0.161 * MathF.Log(rpi) + 0.542;
            else
                Eps_r = 0.136 * MathF.Log(rpi) + 0.542;

            Eps_CO2 = 0.185 * (1 - MathF.Exp(-0.39F * MathF.Pow((float)CGpi, 0.4F)));
            Eps_clear = 1.03 * (Eps_r + Eps_CO2);
            Eps_c = 1 - MathF.Exp((float)(-130 * Wi));

            Epsi = 1 - (1 - Eps_clear) * (1 - Eps_c);

            return Epsi;
        }

        ///<summary>
        ///Epsilon upward
        ///</summary>
        private double Cl_EpsA(int i, int j, int k, ref double EpsApi, ref double EpsAmi, ref double Eabi, ref double Tabi)
        {
            /*
             * k         : height index
             * EpsApi    : Epsilon from the cell-top upwards
             * EpsAmi    : Epsilon from the cell-bootom upwards
             * Eabi      : Epsilon above
             * Tabi      : Temperature above
             * WSearch   : Integral of W until a black cloud
             * CGpSearch : Integral of CGp until a black cloud
             * RpSearch  : Integral of Rp until a black cloud
             * WOutp     : W between cell-top until a black cloud
             * WOutm     : W between cell-bottom until a black cloud
             * CGpOutp   : CGp between cell-top until a black cloud
             * CGpOutm   : CGp between cell-bottom until a black cloud
             * RpOutp    : Rp between cell-top until a black cloud
             * RpOutm    : Rp between cell-bottom until a black cloud
            */
            double WSearch, CGpSearch, RpSearch, WOutp, WOutm, CGpOutp, CGpOutm, RpOutp, RpOutm;
            int kk = k + 1;
            WSearch = 0;
            CGpSearch = 0;
            RpSearch = 0;
            CGpOutm = 0;
            CGpOutp = 0;
            RpOutm = 0;
            RpOutp = 0;
            WOutp = 0;
            WOutm = 0;

            while ((kk <= Program.NZ) && (WSearch < Wcrit))
            {
                WSearch += Wrad[kk][2];
                CGpSearch += CGp[kk][2];
                RpSearch += Rp[kk][2];
                kk++;
            }

            if (WSearch >= Wcrit)
            {
                if (Wrad[kk - 1][2] < Wcrit)
                {
                    WOutp = 0;
                }
                else
                {
                    WOutp = WSearch - Wrad[kk - 1][2];
                    CGpOutp = CGpSearch - CGp[kk - 1][2];
                    RpOutp = RpSearch - Rp[kk - 1][2];
                }
            }
            else
            {
                WOutp = Wrad[k][3] - Wrad[k][2] / (Dz2[k] - Dz1[k]) * (Dz2[k] - Program.ZZ[k]);
                CGpOutp = CGp[k][3] - CGp[k][2] / (Dz2[k] - Dz1[k]) * (Dz2[k] - Program.ZZ[k]);
                RpOutp = Rp[k][3] - Rp[k][2] / (Dz2[k] - Dz1[k]) * (Dz2[k] - Program.ZZ[k]);
            }

            WOutm = WOutp + Wrad[k][2];
            CGpOutm = CGpOutp + CGp[k][2];
            RpOutm = RpOutp + Rp[k][2];

            EpsApi = Cl_Eps(EpsApi, RpOutp, CGpOutp, WOutp);
            EpsAmi = Cl_Eps(EpsAmi, RpOutm, CGpOutm, WOutm);

            if (WSearch >= Wcrit)
            {
                Eabi = 1;
                Tabi = Program.TABS[1][1][kk - 1];
            }
            else
            {
                Eabi = EpsTop;
                Tabi = TTop;
            }

            return Tabi;
        }

        ///<summary>
        ///Terrestrial surface radiation
        ///</summary>
        private void Cl_RterrG(int i, int j)
        {
            //emissivity of the atmosphere in dependence on the cloudyness -> Master theses Manzl, Uni Innsbruck, 2010
            double EpsAtm = 0.7 * (1 - Math.Pow(Program.CLOUDS[i][j], 2.0)) + 0.90 * Math.Pow(Program.CLOUDS[i][j], 2.0);

            Program.RL[i][j] = LStrich[Program.KST[i][j] - 1] * (Q[i][j] + 0.03 * MathF.Sin(Alfa[i][j])) +               ////ACHTUNG Q statt 1-Q; Öttl, Dez 18
                EpsAtm * Program.SIGMA * Math.Pow(0.5 * (Program.TABS[i][j][1] + Program.TABS[i][j][10]), 4) * (1 - Q[i][j]);                                                      ////ACHTUNG 1-Q statt Q; Öttl, Dez 18
            RTerrG[i][j] = Program.EPSG[i][j] * (Program.RL[i][j] - Program.SIGMA * Math.Pow(Program.TG[i][j], 4));
        }

        ///<summary>
        ///Atmospheric heating rates due to terrestrial radiation
        ///</summary>
        private void Cl_dTdtterr(int i, int j)
        {
            /*
             * Dz      :Cell-height [m]
             * EpsAp   :Epsilon from cell-top upwards
             * EpsAm   :Epsilon from cell-bottom upwards
             * Eab     :Epsilon above
             * Tab     :Temperature above
             * EpsBp   :Epsilon from cell-top downwards
             * EpsBm   :Epsilon from cell-bottom downwards
             * Ebel    :Epsilon below
             * Tbel    :Temperature below
             */
            double Dz, EpsBp, EpsBm, Ebel, Tbel;
            EpsBm = 0;
            EpsBp = 0;
            Ebel = 0;
            Tbel = 0;

            for (int k = Program.KST[i][j]; k <= Program.NZ; k++)
            {
                Dz = Dz2[k] - Dz1[k];
                Cl_EpsB(i, j, k, ref EpsBp, ref EpsBm, ref Ebel, ref Tbel);

                Program.DT_TERR[i][j][k] = (1 / (Program.RHOBZZ[k] * Cp * Dz)) *
                    ((Eab[i][j][k] * Program.SIGMA * MathF.Pow((float)Tab[i][j][k], 4) - Program.SIGMA * MathF.Pow((float)Program.TABS[1][1][k], 4)) * (EpsAm[i][j][k] - EpsAp[i][j][k]) +
                    (Ebel * Program.SIGMA * MathF.Pow((float)Tbel, 4) + (1 - Ebel) * Program.RL[i][j] - Program.SIGMA * MathF.Pow((float)Program.TABS[1][1][k], 4)) * (EpsBp - EpsBm));
            }

            Program.DT_TERR[i][j][Program.KST[i][j] - 1] = Program.DT_TERR[i][j][Program.KST[i][j]];
        }

        //Epsilon downward
        private double Cl_EpsB(int i, int j, int k, ref double EpsBp, ref double EpsBm, ref double Ebel, ref double Tbel)
        {
            /*
             * k         : height index
             * EpsBp     : Epsilon from the cell-top downwards
             * EpsBm     : Epsilon from the cell-bootom downwards
             * Ebel      : Epsilon below
             * Tbel      : Temperature below
             * WSearch   : Integral of W until a black cloud
             * CGpSearch : Integral of CGp until a black cloud
             * RpSearch  : Integral of Rp until a black cloud
             * WOutp     : W between cell-top until a black cloud
             * WOutm     : W between cell-bottom until a black cloud
             * CGpOutp   : CGp between cell-top until a black cloud
             * CGpOutm   : CGp between cell-bottom until a black cloud
             * RpOutp    : Rp between cell-top until a black cloud
             * RpOutm    : Rp between cell-bottom until a black cloud
            */
            double WSearch, CGpSearch, RpSearch, WOutp, WOutm, CGpOutp, CGpOutm, RpOutp, RpOutm;

            int kk = k - 1;
            WSearch = 0;
            CGpSearch = 0;
            RpSearch = 0;
            CGpOutm = 0;
            CGpOutp = 0;
            RpOutm = 0;
            RpOutp = 0;
            WOutp = 0;
            WOutm = 0;
            while ((kk >= (Program.KST[i][j] - 1)) && (WSearch < Wcrit))
            {
                WSearch += Wrad[kk][2];
                CGpSearch += CGp[kk][2];
                RpSearch += Rp[kk][2];
                kk--;
            }
            if (WSearch >= Wcrit)
            {
                if (Wrad[kk + 1][2] < Wcrit)
                {
                    WOutp = 0;
                }
                else
                {
                    WOutp = WSearch - Wrad[kk + 1][2];
                    CGpOutp = CGpSearch - CGp[kk + 1][2];
                    RpOutp = RpSearch - Rp[kk + 1][2];
                }
            }
            else
            {
                WOutp = WSearch;
                CGpOutp = CGpSearch;
                RpOutp = RpSearch;
            }
            WOutm = WOutp + Wrad[k][2];
            CGpOutm = CGpOutp + CGp[k][2];
            RpOutm = RpOutp + Rp[k][2];

            EpsBp = Cl_Eps(EpsBp, RpOutp, CGpOutp, WOutp);
            EpsBm = Cl_Eps(EpsBm, RpOutm, CGpOutm, WOutm);

            if (WSearch >= Wcrit)
            {
                Ebel = 1;
                Tbel = Program.TABS[1][1][kk + 1];
            }
            else
            {
                Ebel = Program.EPSG[i][j];
                Tbel = Program.TG[i][j];
            }

            return Tbel;
        }
    }
}
