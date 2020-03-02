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
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace GRAMM_CSharp_Test
{
    partial class Program
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

        public static int Na = 0;                   //Number of days since 1 Jan
        public static double timeR = 0;             //Time after midnight in seconds
        public static int KStMax = 0;               //Vertical index of the highest surface cell
        public const double Eps = 0.00005;          //General lower threshold
        public const int nOmega = 40;               //Number of sectors used to compute the horizons of the surroundings (should be divisible by 4)
        public static float GPhi = 0;              //Variation of the solar flux over the year
        public static float Mye = 0;               //Sinus of sun heigth (He)
        public static float My = 0;                //Sinus of sun heigth (H)
        public static float CosHe = 0;             //Cosinus He(Mye)
        public static float CosH = 0;              //Cosinus He(My)
        public static float He = 0;                //Angle of sun height [rad]
        public static float Eta = 0;               //Azimuth of the sun [rad]
        public const double GrRa = 0.0174532925199;               //Conversion factor PI/180 in rad
        public const double RaGr = 57.2957795131;                 //Conversion factor 180/PI in deg
        public const double p0 = 1.01325E5;                       //reference pressure
        public const double eHoch1 = 2.7182818285;                //Eulerian constant
        public const double IG0 = 1367;                           //extraterrestrial solar radiation density [W/m²]
        public const double cp = 1004;                            //heat capacity of air [Joule/kg Kelvin]
        public const double QCo2 = 0.00052;                       //mass-concentration of CO2
        public const double Wcrit = 0.03;                         //W-integral of black cloud [kg/m²]
        public const double EpsTop = 0.2;                         //Epsilon at the top of the atmosphere
        public const double TTop = 216.65;                        //Temperature at the top of the atmosphere
        public const double AgNorm = 0.6;                         //Default surface albedo
        public static double a0 = 0;                              //Intermediate value for computing radiation density
//        public static double SinIe = 0;                         //sun angle above surface (Mye) [rad]
//        public static double SinI = 0;                          //sun angle above surface (My) [rad]

        public static void RADIATRAD(bool lGeom,bool lGeomWr,bool lGeomRe,int iTag,int iMon,int iJahr,double Zeit, int NI, int NJ, int NK)
        {
            int iTheta; // int Kalb;

            if (NK != 3) Console.Write("RADIATION started");
            if((Program.ISOL != 1) && (Program.ISOL != 2))
            {
                Console.WriteLine("===========================================");
                Console.WriteLine("Selected radiation scheme must equal 1 or 2");
                Console.WriteLine("===========================================");
                Environment.Exit(0);
            }

            string Date1 = "0101" + iJahr.ToString("0000");
            string Date2 = iTag.ToString("00") + iMon.ToString("00") + iJahr.ToString("0000");
            Na = 1 + DateNum(Date2) - DateNum(Date1);

            timeR = Zeit;

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
            if(lGeom == true  && File.Exists("albeq.dat") == false) // if possible, read albeq.dat
            {

            	Parallel.For(1, NI+1, Program.pOptions, i =>
                {
                   for (int j = 1; j <= NJ; j++)
                    {
                        Program.Q[i][j] = 0;
                        Program.Alfa[i][j] = 0;
                        Program.Beta[i][j] = 0;
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
                        float HX = (float) ((ZE - ZW) / (XE - XW)); // inclination gradients of the ground
                        float HY = (float) ((ZN - ZS) / (YN - YS));
                        
                        Program.Alfa[i][j] = MathF.Atan(MathF.Sqrt(HX * HX + HY * HY)); //inclination of the ground cell 
                        
                        if (Math.Abs(Program.Alfa[i][j]) < Eps)
                            Program.Beta[i][j] = 0; // direction of the ground cell 0 = north
                        else
                            Program.Beta[i][j] = MathF.Atan2(HX, HY);
                        
                        Program.Q[i][j] = CLQ(i, j, Program.Alfa[i][j], Program.Beta[i][j]);  //sky view factor  = 0 in case of flat terrain (1 - Q = Horizontüberhöhung)
                    }
            	  });

                if (NK != 3) Console.Write(" / geometry computed");
                //save geometry data
                if(lGeomWr==true)
                {
                	try
                	{
                		using (BinaryWriter wr = new BinaryWriter(File.Open("albeq.dat", FileMode.Create)))
                		{
                			for (int i = 1; i <= Program.NX; i++)
                				for (int j = 1; j <= Program.NY; j++)
                			{
                				{
                					wr.Write((double) Program.Alfa[i][j]);
                					wr.Write((double) Program.Beta[i][j]);
                					wr.Write(Program.Q[i][j]);
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
            					Program.Alfa[i][j] = (float) br.ReadDouble();
            					Program.Beta[i][j] = (float) br.ReadDouble();
            					Program.Q[i][j] = br.ReadDouble();
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
            DOmega = 2 * Math.PI / nOmega;
            iTheta = (int)((Eta + Math.PI) / (2 * Math.PI / nOmega)) + 1;     //sector where eta lies
            if (iTheta > nOmega) iTheta = 1;
            Omega = DOmega * (iTheta - 1) + DOmega * 0.5;                               //actual angle
 			
			//Console.WriteLine("Itime in hours" + (timeR/3600).ToString() + " Omega" + Omega.ToString() + " Eta "+ Eta.ToString());

            //main loop for computing solar radiation
			Parallel.For(2, NI, Program.pOptions, i =>
            { 
      		   double SinIe = 0;                            //sun angle above surface (Mye) [rad]
       		   double SinI=0;
      		   double ThetaB=0;
     		   double ThetaOld = 0;
     		   double Theta = 0;
//            for (int i = 2; i <= Program.NX - 1; i++)
//            {
                for (int j = 2; j <= NJ - 1; j++)
                {
                	// sun angle over the surface
                    SinIe = Mye * MathF.Cos(Program.Alfa[i][j]) + CosHe * MathF.Sin(Program.Alfa[i][j]) * MathF.Cos(Eta - Program.Beta[i][j]);
                    SinI = My * MathF.Cos(Program.Alfa[i][j]) + CosH * MathF.Sin(Program.Alfa[i][j]) * MathF.Cos(Eta - Program.Beta[i][j]);
                    
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
            
            if (NK != 3)  Console.WriteLine(" / radiation computed");
            //terrestrial downward longwave radiation
            Cl_LStrich();

            //Epsilon Above
            for (int k = 2; k <= NK; k++)
            {
                for (int i = 1; i <= NI; i++)
                {
                    for (int j = 1; j <= NJ; j++)
                    {
                        Cl_EpsA(i, j, k, ref Program.EpsAp[i][j][k], ref Program.EpsAm[i][j][k], ref Program.Eab[i][j][k], ref Program.Tab[i][j][k]);
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
                    Program.Arad[i][1][k] = Program.Arad[i][2][k];
                    Program.SG[i][1][k] = Program.SG[i][2][k];
                    Program.DG[i][1][k] = Program.DG[i][2][k];
                    Program.EG[i][1][k] = Program.EG[i][2][k];
                    Program.DT_SOL[i][1][k] = Program.DT_SOL[i][2][k];
                    Program.DT_TERR[i][1][k] = Program.DT_TERR[i][2][k];

                    Program.Arad[i][Program.NY][k] = Program.Arad[i][Program.NY - 1][k];
                    Program.SG[i][Program.NY][k] = Program.SG[i][Program.NY - 1][k];
                    Program.DG[i][Program.NY][k] = Program.DG[i][Program.NY - 1][k];
                    Program.EG[i][Program.NY][k] = Program.EG[i][Program.NY - 1][k];
                    Program.DT_SOL[i][Program.NY][k] = Program.DT_SOL[i][Program.NY - 1][k];
                    Program.DT_TERR[i][Program.NY][k] = Program.DT_TERR[i][Program.NY - 1][k];
            	});
            	Parallel.For(2, NJ, Program.pOptions, j =>
                {
//                for (int j = 2; j <= Program.NY - 1; j++)
//                {
                    Program.Arad[1][j][k] = Program.Arad[2][j][k];
                    Program.SG[1][j][k] = Program.SG[2][j][k];
                    Program.DG[1][j][k] = Program.DG[2][j][k];
                    Program.EG[1][j][k] = Program.EG[2][j][k];
                    Program.DT_SOL[1][j][k] = Program.DT_SOL[2][j][k];
                    Program.DT_TERR[1][j][k] = Program.DT_TERR[2][j][k];

                    Program.Arad[Program.NX][j][k] = Program.Arad[Program.NX - 1][j][k];
                    Program.SG[Program.NX][j][k] = Program.SG[Program.NX - 1][j][k];
                    Program.DG[Program.NX][j][k] = Program.DG[Program.NX - 1][j][k];
                    Program.EG[Program.NX][j][k] = Program.EG[Program.NX - 1][j][k];
                    Program.DT_SOL[Program.NX][j][k] = Program.DT_SOL[Program.NX - 1][j][k];
                    Program.DT_TERR[Program.NX][j][k] = Program.DT_TERR[Program.NX - 1][j][k];
            	 });
            }
        }        

        //Number of days since 1.1.1600
        public static int DateNum(string Date)
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
            DateTime firstdate  = new DateTime(1600,1,1);
			DateTime seconddate = new DateTime(Year,Month,Day);
            TimeSpan ts = seconddate - firstdate;
            DateNum = ts.Days; // calculates the number of days between two dates
            return DateNum;
        }

        
        private static double CLQ(int i, int j, double Alfai, double Betai)
        {
        	//Computation of the projection to the ground cell
        	// Alfai = inclination of the ground cell
        	// Betai azimuth of the ground cell
        	double DOmega, Omega;
        	//Projection of surface cell
        	DOmega = 2 * Math.PI / (float)nOmega;
        	double QI = 0;
        	double Theta = 0;
        	for (int N = 1; N <= nOmega; N++)
        	{
        		Omega = DOmega * (N - 1) + DOmega * 0.5;                                            //compute sectors
        		float AlfEff = (float) (-Program.Alfa[i][j] * Math.Cos(Program.Beta[i][j] - Omega));         //actual angle
        		//Compute horizons of surroundings
        		Theta = CLThet(i, j, Program.AH[i][j], Omega, DOmega, Theta);
        		float TheS = (float) Math.Max(Theta, AlfEff);
        		double DQi = Math.Cos(AlfEff) * 0.5F * (MathF.Pow(MathF.Sin(TheS), 2) - MathF.Pow(MathF.Sin(AlfEff),2))
        			- MathF.Sin(AlfEff) * 0.5F * ((MathF.Sin(2 * TheS) - MathF.Sin(2 * AlfEff)) * 0.5
        			+TheS - AlfEff);
        		QI += DQi;                                                                               //sum up projection
        	}

        	QI *= DOmega / Math.PI;
			
        	return QI;
		}
        
        //Horizons of surroundings
        public static double CLThet(int I, int J, float Hoehe, double Omega, double DOmega, double Theta)
        {
        	// search horizon for radiation or shadow estimation
            int IAkt, JAkt, IZ, JZ, Nenner;
            double DXrad, DYrad, Dx1, Dx2, Dy1, Dy2, XSum, YSum, ZSum, ThetaI;
            bool Ende;
            float OM1 = (float) (Omega - DOmega * 0.5);
            float OM2 = (float) (Omega + DOmega * 0.5);
            Theta = 0;

            if((Omega > Math.PI * 1.75) || (Omega <= Math.PI * 0.25)) // northern sector
            {
                IAkt = I;             //angle in the north sector
                JZ = J + 1;
                Ende = false;
                
                while((JZ <= Program.NY) && (Ende != true))
                {
                    DYrad = Program.Y[JZ] - Program.Y[J];   //distance JZ-J
                    Dx1 = MathF.Tan(OM1) * DYrad;            //distance in x-directions of OM1
                    Dx2 = MathF.Tan(OM2) * DYrad;            //distance in x-directions of OM2
                    if(((Program.X[I]+Dx2)>=Program.X[Program.NX])||(Program.X[I]+Dx1)<=Program.X[1])
                    {
                        Ende = true; 						// sector outside the computation area
                        Theta = Math.Max(Theta, 0);
                    }
                    else
                    {
                        if(Program.X[IAkt] < (Program.X[I] + Dx1)) // counter less OM1
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
                            while(Program.X[IZ] <= (Program.X[I] + Dx2))
                            {
                                Nenner++;
                                XSum += Program.X[IZ] - Program.X[I];
                                ZSum += Program.AH[IZ][JZ] - Hoehe;
                                IZ++;
                            }
                            XSum /= Nenner;
                            ZSum /= Nenner;
                        }
                        ThetaI = Math.Atan2(ZSum, Math.Sqrt(XSum * XSum + DYrad * DYrad));
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
                while((IZ<=Program.NX)&&(Ende==false))
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
                        ThetaI = Math.Atan2(ZSum, Math.Sqrt(YSum * YSum + DXrad * DXrad));
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
                        ThetaI = Math.Atan2(ZSum, Math.Sqrt(XSum * XSum + DYrad * DYrad));
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
                        ThetaI = Math.Atan2(ZSum, Math.Sqrt(YSum * YSum + DXrad * DXrad));
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

        //Sun height, angle, and azimuth
        public static double ClMyEt()
        {
            float PhiThe;
            float Delta;       //sun declination [rad]
            float Psirad;      //sun-hour angle [rad]
            double CosEta;      //cosinus of sun azimuth [rad]

            PhiThe = 2 * MathF.PI * ((float)Na - 2.84F) / 365;
            GPhi = 1.0006F + 0.03343F * MathF.Cos(PhiThe) + 0.0011F * MathF.Sin(PhiThe);   //radiation flux variation
            Delta = MathF.Asin(MathF.Sin(23.5F * MathF.PI / 180) * MathF.Sin(2 * MathF.PI * (Na - 80) / 365));   //sun declination
            Psirad = (float) (timeR / 43200 * Math.PI - Math.PI);   // horizontal sun angle in hours 0:00 = -Pi, 12:00 = 0, 24:00 = Pi
            
            Mye = (float) (Math.Sin(Program.BGRAD * GrRa) * MathF.Sin(Delta) + Math.Cos(Program.BGRAD * GrRa) * MathF.Cos(Delta) * MathF.Cos(Psirad));    //sinus of sun angle
            My = MathF.Max(0.02F, Mye);					  // needed to avoid errors if sun under the horizon
            CosHe = MathF.Sqrt(1 - Mye * Mye);             // cosinus of sun angle
            CosH = MathF.Sqrt(1 - My * My);
            He = MathF.Asin(Mye);                          // vertical sun angle
            a0 = 0.04 + 0.065 * MathF.Exp(-0.18F / My);

            if ((MathF.Abs(CosHe) < Eps) || (Program.BGRAD * GrRa == MathF.PI * 0.5F) || (MathF.Abs(Psirad) < Eps))
            {
            	Eta = 0; // at the pole, at noon
            }
            else
            {
                CosEta = (Mye * Math.Sin(Program.BGRAD * GrRa) - MathF.Sin(Delta)) / (CosHe * Math.Cos(Program.BGRAD * GrRa));
                if (Math.Abs(CosEta) <= 1 + Eps)
                {
                    if(CosEta>1)
                        CosEta = 1;
                    if (CosEta < -1)
                        CosEta = -1;
                    Eta = (float) Math.Acos(CosEta);
                    if (Psirad < 0)
                        Eta = -Eta;
                }
                else
                {
                    Console.WriteLine("Error when computing Eta in the radiation model (-1<coseta<1 not fulfilled)");
                    Environment.Exit(0);
                }
            }
            
            if (Program.BGRAD < 0) // Kuntner 22.3.2017 :on the southern hemisphere - set the sun position to the north
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
            		Eta -= 2* MathF.PI;
            	
            	if (Eta < Eps)
            		Eta = MathF.PI;
            }
            //Console.WriteLine("ETA " + Eta.ToString());
            return He;
        }

        //Surface albedo
        public static double ClAg()
        {
        	Parallel.For(1, Program.NX + 1, Program.pOptions, I =>
            {
//            for (int I = 1; I <= Program.NX; I++)
//            {
                for (int J = 1; J <= Program.NY; J++)
                {
                    Program.Ag[I][J] = Program.ALBEDO[I][J];
                    if (Program.Ag[I][J] == 0.08)
                        Program.Ag[I][J] = MathF.Min(1, MathF.Max(0.03F, -0.0139F + MathF.Tan(MathF.PI * 0.5F - He)));    //Flassak p 133
                }
        	});
            return He;   
        }

        //Vertical integrals of haze, etc.
        public static double ClYrWc()
        {
        	//N,               ! Zaehlervariable
     	    //K                ! Hoehenindex der Zelle
        	double[] Rhop = new double[Program.NPROFMAX];      //Density [kg/m³]
            double[] QRp = new double[Program.NPROFMAX];       //Intermediate value for calculating Rp F50a
            double[] QCp = new double[Program.NPROFMAX];       //Intermediate value for calculating Cp F50b
            double[] Vnz = new double[Program.NPROFMAX];       //Intermediate value
            double[] QVapz = new double[Program.NPROFMAX];     //Intermediate value
            double[] QCldz = new double[Program.NPROFMAX];     //Intermediate value
            //double[][][] QCldz = CreateArray<double[][]>(Program.NX1, () => CreateArray<double[]>(Program.NY1, () => new double[Program.NPROFMAX]));     //Intermediate value
            double[] QIcez = new double[Program.NPROFMAX];     //Intermediate value

            
            for (int N = 1; N <= 30; N++)
            {
                Rhop[N] = Program.PPROF[N] / 287 / Program.TPROF[N];             //general gas equation
                QRp[N] = Rhop[N] * Program.QVAP[N] * Program.PPROF[N] / p0;      //F50a
                QCp[N] = Rhop[N] * QCo2 * Program.PPROF[N] / p0;                 //F50b
                Vnz[N] = 570 * (1 / MathF.Pow((float) Program.VNORM[N], 1.5F));  //F13a
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
                            Program.QCLD[N] = 0.0;
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
                                Program.QCLD[N] = Program.WAT_VAP[i][j][INDU] + (Program.ZPROF[N] - Program.ZSP[i][j][INDU]) * GRAD;
                            }
                            else if (INDU == 0)
                            {
                                GRAD = (Program.WAT_VAP[i][j][INDO] - Program.WAT_VAP[i][j][INDO + 1]) /
                                    (Program.ZSP[i][j][INDO] - Program.ZSP[i][j][INDO + 1]);
                                Program.QCLD[N] = Program.WAT_VAP[i][j][INDO] + (Program.ZPROF[N] - Program.ZSP[i][j][INDO]) * GRAD;
                            }
                            else
                            {
                                GRAD = (Program.WAT_VAP[i][j][INDO] - Program.WAT_VAP[i][j][INDU]) /
                                    (Program.ZSP[i][j][INDO] - Program.ZSP[i][j][INDU]);
                                Program.QCLD[N] = Program.WAT_VAP[i][j][INDO] + (Program.ZPROF[N] - Program.ZSP[i][j][INDO]) * GRAD;
                            }
                        }
                        QCldz[i][j][N] = 0 * Rhop[N] * (Math.Max(Program.QCLD[N], 0) * 0.001 + Program.QRAIN[N]);     //F13c
                    }
                }
                */
            }
 			
            // it seems, that there is a dependency, so parallelisation is not possible!
            for (int K = 1; K <= Program.NZ; K++) // start of integration
            {
                if(K==1) // Border of the 1st cell
                {
                    Program.Dz1[K] = 0;                                          //borders of first cell
                    Program.Dz2[K] = (Program.ZZ[K] + Program.ZZ[K + 1]) * 0.5;
                }
                else if(K==Program.NZ) // border of the last cell k
                {
                    Program.Dz1[K] = (Program.ZZ[K] + Program.ZZ[K - 1]) * 0.5;    //borders of last cell
                    Program.Dz2[K] = Program.ZZ[K] + (Program.ZZ[K] - Program.ZZ[K - 1]) * 0.5;
                }
                else
                {
                    Program.Dz1[K] = (Program.ZZ[K] + Program.ZZ[K - 1]) * 0.5;    //borders of the k-th cell
                    Program.Dz2[K] = (Program.ZZ[K] + Program.ZZ[K + 1]) * 0.5;
                }
                Program.TABS[1][1][K] = rInteg(Program.ZPROF, Program.TPROF, Program.NPROFMAX, Program.Dz1[K], Program.Dz2[K]) / (Program.Dz2[K] - Program.Dz1[K]);   //temperature in k-th cell
                Program.PBZZ[K] = rInteg(Program.ZPROF, Program.PPROF, Program.NPROFMAX, Program.Dz1[K], Program.Dz2[K]) / (Program.Dz2[K] - Program.Dz1[K]);   //pressure in k-th cell

                Program.RHOBZZ[K] = Program.PBZZ[K] / 287 / Program.TABS[1][1][K]; // gas equation
                Program.YG[K][3] = rInteg(Program.ZPROF, Vnz, Program.NPROFMAX, Program.ZZ[K], Program.ZPROF[30]);      // valus direction up 
                Program.R[K][3] = rInteg(Program.ZPROF, QVapz, Program.NPROFMAX, Program.ZZ[K], Program.ZPROF[30]);
                Program.W2rad[K][3] = rInteg(Program.ZPROF, QIcez, Program.NPROFMAX, Program.ZZ[K], Program.ZPROF[30]);
                Program.Rp[K][3] = rInteg(Program.ZPROF, QRp, Program.NPROFMAX, Program.ZZ[K], Program.ZPROF[30]);
                Program.CGp[K][3] = rInteg(Program.ZPROF, QCp, Program.NPROFMAX, Program.ZZ[K], Program.ZPROF[30]);

                Program.YG[K][1] = Program.YG[1][3] - Program.YG[K][3];        // values directon down
                Program.R[K][1] = Program.R[1][3] - Program.R[K][3];
                Program.W2rad[K][1] = Program.W2rad[1][3] - Program.W2rad[K][3];
                Program.Rp[K][1] = Program.Rp[1][3] - Program.Rp[K][3];
                Program.CGp[K][1] = Program.CGp[1][3] - Program.CGp[K][3];

                Program.YG[K][2] = rInteg(Program.ZPROF, Vnz, Program.NPROFMAX, Program.Dz1[K], Program.Dz2[K]); // value at the cell
                Program.R[K][2] = rInteg(Program.ZPROF, QVapz, Program.NPROFMAX, Program.Dz1[K], Program.Dz2[K]);
                Program.W2rad[K][2] = rInteg(Program.ZPROF, QIcez, Program.NPROFMAX, Program.Dz1[K], Program.Dz2[K]);
                Program.Rp[K][2] = rInteg(Program.ZPROF, QRp, Program.NPROFMAX, Program.Dz1[K], Program.Dz2[K]);
                Program.CGp[K][2] = rInteg(Program.ZPROF, QCp, Program.NPROFMAX, Program.Dz1[K], Program.Dz2[K]);

                Program.W1rad[K][3] = rInteg(Program.ZPROF, QCldz, Program.NPROFMAX, Program.ZZ[K], Program.ZPROF[30]);
                Program.W1rad[K][1] = Program.W1rad[1][3] - Program.W1rad[K][3];
                Program.W1rad[K][2] = rInteg(Program.ZPROF, QCldz, Program.NPROFMAX, Program.Dz1[K], Program.Dz2[K]);
                Program.Wrad[K][3] = Program.W1rad[K][3] + 0.5F * Program.W2rad[K][3]; // F13e
                Program.Wrad[K][1] = Program.W1rad[K][1] + 0.5F * Program.W2rad[K][1];
                Program.Wrad[K][2] = Program.W1rad[K][2] + 0.5F * Program.W2rad[K][2];

                /*
                for (int i = 1; i <= Program.NX; i++)
                {
                    for (int j = 1; j <= Program.NY; j++)
                    {
                        Program.W1rad[i][j][K][3] = rInteg(Program.ZPROF, QCldz, Program.NPROFMAX, Program.ZZ[K], Program.ZPROF[30]);
                        Program.W1rad[i][j][K][1] = Program.W1rad[i][j][1][3] - Program.W1rad[i][j][K][3];
                        Program.W1rad[i][j][K][2] = rInteg(Program.ZPROF, QCldz, Program.NPROFMAX, Program.Dz1[K], Program.Dz2[K]);

                        Program.Wrad[i][j][K][3] = Program.W1rad[i][j][K][3] + 0.5F * Program.W2rad[K][3]; // F13e
                        Program.Wrad[i][j][K][1] = Program.W1rad[i][j][K][1] + 0.5F * Program.W2rad[K][1];
                        Program.Wrad[i][j][K][2] = Program.W1rad[i][j][K][2] + 0.5F * Program.W2rad[K][2];
                    }
                }
                */
             } 
            
            return 0d;   
        }

        //integration over a field of points
        public static double rInteg(double[] Xvalue,double[] Yvalue,int nMax,double X1,double X2)
        {
            double rInteg = 0;
            double Y = 0;
            int N = 2;
            while(X1>Xvalue[N])
            {
                N++;
            }
            if(X2<Xvalue[N])
            {
                Y = Yvalue[N - 1] + ((X1 + X2) * 0.5F - Xvalue[N - 1]) * ((Yvalue[N] - Yvalue[N - 1]) / (Xvalue[N] - Xvalue[N - 1]));
                rInteg = Y * (X2 - X1);
            }
            else
            {
                Y = Yvalue[N - 1] + ((X1 + Xvalue[N]) * 0.5F - Xvalue[N - 1]) * ((Yvalue[N] - Yvalue[N - 1]) / (Xvalue[N] - Xvalue[N - 1]));
                rInteg = Y * (Xvalue[N] - X1);
                N++;
                while(Xvalue[N]<X2)
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

        //General height dependent values
        public static double CLTStern()
        {
         	double Tau1, Tau2, TauS;
            
         	for (int k = 1; k <= Program.NZ;k++ )
            {
                Program.Tstern[k] = 0.3 + (MathF.Pow(0.1F * (float) Program.R[k][3], 0.2F)) * (0.75 + My) +
                    (8.0 * p0 / Program.PBZZ[k]) * (Program.YG[k][3] + 0.1F * Program.PBZZ[k] / p0 + 0.02);
                Program.IG[k] = GPhi * IG0 * Solar_reduction_factor * MathF.Exp((float) (-1 * ((Program.Tstern[k] / My) * a0)));

                float tmp = 1000F * (float) Program.W1rad[k][3] + 1;
                if (MathF.Log(tmp) == 0)
                {
                    Tau1 = 0;
                }
                else
                {
                    Tau1 = 0.07 * MathF.Pow(MathF.Log(tmp), 3.95F);
                }
                Tau2 = 56.6 * Program.W2rad[k][3];
                Program.Tau[k] = Tau1 + Tau2;
                TauS = Math.Min(32, Program.Tau[k]);
                Program.CloudES[k] = MathF.Pow((float) (1 - TauS / 32), 4);
                Program.CloudESeS[k] = Math.Max(0.1f, Program.CloudES[k]);

                Program.nS[k] = 0.346 * (Program.CloudESeS[k] - 0.1) + 0.85 * MathF.Pow((float) Program.CloudESeS[k] - 0.1F, 2);
                Program.nD[k] = 0.3 * MathF.Sin(MathF.PI * (Program.CloudES[k] + 0.25F)) + 0.36F * Program.CloudES[k] - 0.15F;

                /*
                for (int i = 1; i <= Program.NX; i++)
                {
                    for (int j = 1; j <= Program.NY; j++)
                    {
                        Tau2 = 56.6 * Program.W2rad[k][3];
                        Program.Tau[i][j][k] = Tau1 + Tau2;
                        TauS = Math.Min(32, Program.Tau[i][j][k]);
                        Program.CloudeS[i][j][k] = Math.Pow(1 - TauS / 32, 4);
                        Program.CloudeSeS[i][j][k] = Math.Max(0.1, Program.CloudeS[i][j][k]);

                        Program.nS[i][j][k] = 0.346 * (Program.CloudeSeS[i][j][k] - 0.1) + 0.85 * Math.Pow(Program.CloudeSeS[i][j][k] - 0.1, 2);
                        Program.nD[i][j][k] = 0.3 * Math.Sin(Math.PI * (Program.CloudeS[i][j][k] + 0.25)) + 0.36 * Program.CloudeS[i][j][k] - 0.15;
                    }
                }
                */
            }

            return 0d;
        }

        //Computing albedo in each cell
        public static double ClAlbedo(int i, int j, int k, int KSti, double Theta, double SinI, double SinIe)
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
		   
		   double WUnt, TClearUP, TClearDown, TCUp, TCDown, Ac, AgMod;
           float AcZw;

		   Program.Myx[k] = My * Program.CloudES[k] + 0.5F * (1 - Program.CloudES[k]);
		   WUnt = Program.Wrad[KSti][3] - Program.Wrad[k][3];
		   TClearUP = Math.Pow(Program.PBZZ[k] / Program.PBZZ[KSti], 0.25);
		   TClearDown = TClearUP;
		   TCDown = 1.2 * MathF.Pow((float) Program.Myx[k] + 0.2F, 0.5F) * (1.3 * MathF.Pow(1000F * (float) WUnt + 0.1F, -0.33F) - 0.06);
		   TCDown = Math.Min(1, TCDown);
		   TCDown = Math.Max(0, TCDown);
		   TCUp = 1.3 * MathF.Pow(1000F * (float) WUnt + 0.1F, -0.33F) - 0.06;
		   TCUp = Math.Min(1, TCUp);
		   TCUp = Math.Max(0, TCUp);
		   AcZw = (float) (2600 * WUnt * (1.2 - Program.Myx[k]));
		  
		   if (AcZw <= eHoch1)
		   	Ac = 0;
		   else
		   	Ac = 0.42 * MathF.Log(MathF.Log(AcZw));
		   
		   if ((SinIe < 0) || (He < Theta))
		   	AgMod = Program.Ag[i][j] * (1 - Program.CloudES[KSti]);
		   else
		   	AgMod = Program.Ag[i][j] * (1 - Program.CloudES[KSti] + Program.CloudES[KSti] * SinI / My);
		   
		   Program.Arad[i][j][k] = Ac + (TClearUP * TClearDown * TCUp * TCDown * AgMod) / (1 - Ac * AgMod);

		   return Ac;
        }

        //solar scheme 1 (thin clouds)
        public static double ClSol_1(int i, int j, int k, double Theta)
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

            fA = 1 + (1 - 0.5F * Program.CloudES[k]) * (Program.Arad[i][j][k] - 0.2);
            
            if(Program.Tau[k]<2)
            {
                jzw = ((0.07 * Program.Tstern[k] - 0.115) / My) * (Program.PBZZ[k] / p0);
                j32 = Math.Max(0, jzw);
            }
            else
            {
                Tauj = MathF.Max(32, (float) Program.Tau[k]);
                j32 = 45.25 * MathF.Pow(Tauj, -1.5F);
            }
            
            if(He<Theta)
            {
                Program.SG[i][j][k] = 0;
            }
            else
            {
                Program.SG[i][j][k] = Program.IG[k] * Mye * Program.nS[k] * (Program.PBZZ[k] / p0) * (1 - 0.233 * Program.CLOUDS[i][j] - 0.415 * Pow2(Program.CLOUDS[i][j]));
            }

            Program.DG[i][j][k] = Program.IG[k] * My * (Program.nD[k] + j32) * fA * (Program.PBZZ[k] / p0) * (1 - 0.233 * Program.CLOUDS[i][j] - 0.415 * Pow2(Program.CLOUDS[i][j]));
            Program.EG[i][j][k] = Program.SG[i][j][k] + Program.DG[i][j][k];

            return jzw;
        }

        //solar scheme 2 (thick clouds)
        public static double ClSol_2(int i, int j, int k, double Theta)
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
            double fAE, jzw, j38, Tc, SStrich, DStrich;
            jzw = 0;

            jzw = ((0.07 * Program.Tstern[k] - 0.115) / My) * (Program.PBZZ[k] / p0);
            j38 = Math.Max(0, jzw);
            fAE = 1 + (1 - 0.9 * Program.CloudES[k]) * (Program.Arad[i][j][k] - 0.2);
            Tc = 1.2 * MathF.Pow((float) Mye + 0.2F, 0.5F) * (1.3 * MathF.Pow((float) (1000F * Program.Wrad[k][3] + 0.1F), -0.33F) - 0.06);
            Tc = Math.Min(1, Tc);
            Tc = Math.Max(0, Tc);

            if (He < Theta)
            {
                SStrich = 0;
            }
            else
            {
                SStrich = Program.IG[k] * Mye;
            }
            
            DStrich = Program.IG[k] * My * j38;
            Program.EG[i][j][k] = (SStrich + DStrich) * Tc * fAE;
            
            if (He < Theta)
            {
                Program.SG[i][j][k] = 0;
            }
            else
            {
                Program.SG[i][j][k] = 0.8 * Program.EG[i][j][k] * (Program.CloudESeS[k] - 0.1) / 0.9;
            }

            Program.DG[i][j][k] = Program.EG[i][j][k] - Program.SG[i][j][k];
            

            return jzw;
        }

        //Solar radiation at the surface
        public static double RSolGround(int i, int j, int k, double Theta, double SinIe)
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
            
            XH = 1.5 * My * Program.CloudES[k] * (0.65 + 0.04 * Program.Tstern[k] * Program.PBZZ[k] / p0 - 0.6 * Program.CloudES[k]);
            
            if((SinIe<0)||(He<Theta))
            {
                SGG = 0;
                DGSG = 0;
            }
            else
            {
                SGG = Program.SG[i][j][k] * SinIe / Mye;
                DGSG = XH * SinIe * Program.DG[i][j][k] / My;
            }
            
            DGDG = (1 - XH) * (1 - Program.Q[i][j]) * Program.DG[i][j][k];
            BGG = Program.Ag[i][j] * Program.Q[i][j] * (Program.SG[i][j][k] + Program.DG[i][j][k]);

            Program.GLOBRAD[i][j] =(SGG + DGSG + DGDG + BGG);
            Program.RSOLG[i][j] = (1 - Program.Ag[i][j]) * Program.GLOBRAD[i][j];

            //influence of snow cover, except in urban areas and water bodies
            if (Program.SNOW[i][j] > 0.15 && Program.SNOW[i][j] > Program.Z0[i][j] && Program.ALAMBDA[i][j] != 4.0 && Program.FW[i][j] < 0.99)
            {
                Program.RSOLG[i][j] = (1 - 0.6) * Program.GLOBRAD[i][j];
            }

            return XH;
        }

        //Heating rate caused by solar radiation
        public static double dTdt_Sol(int i, int j)
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
                Dz = (Program.Dz2[k] - Program.Dz1[k]);
                IceDens = Program.W2rad[k][2] / Dz;
                WatDens = Program.W1rad[k][2] / Dz;
                AzMin = Program.Arad[i][j][k];
                
                if(Program.SG[i][j][k]!=0)
                {
                    EzMin = Program.EG[i][j][k] + (Program.Dz1[k] - Program.ZZ[k]) * ((Program.EG[i][j][k + 1] - Program.EG[i][j][k]) / (Program.ZZ[k + 1] - Program.ZZ[k]));
                    EzPlus = Program.EG[i][j][k] + (Program.Dz2[k] - Program.ZZ[k]) * ((Program.EG[i][j][k + 1] - Program.EG[i][j][k]) / (Program.ZZ[k + 1] - Program.ZZ[k]));
                }
                else
                {
                    EzMin = Program.EG[i][j][k - 1] + (Program.Dz1[k] - Program.ZZ[k - 1]) * ((Program.EG[i][j][k] - Program.EG[i][j][k - 1]) / (Program.ZZ[k] - Program.ZZ[k - 1]));
                    EzPlus = Program.EG[i][j][k - 1] + (Program.Dz2[k] - Program.ZZ[k - 1]) * ((Program.EG[i][j][k] - Program.EG[i][j][k - 1]) / (Program.ZZ[k] - Program.ZZ[k - 1]));
                }
                
                if((IceDens < 0.000001) && (WatDens < 0.00001))
                {
                    Program.DT_SOL[i][j][k] = 1 / (Program.RHOBZZ[k] * cp * Dz) * Math.Abs(EzPlus - EzMin) * 0.7 * (1 + AzMin);
                }
                else
                {
                    float tmp = (float) (1 + 1000 * Program.Wrad[k][2]);
                    AcDown = 0.025 * Math.Pow((Program.Myx[k] + Program.Myx[k + 1]) * 0.5F + 0.1, 0.5F) * MathF.Log(tmp);
                    AcUp = 0.0194 * MathF.Log(tmp);
                    Program.DT_SOL[i][j][k] = 1 / (Program.RHOBZZ[k] * cp * Dz) * Math.Abs(EzPlus * AcDown + EzMin * AzMin * AcUp);
                }
            }
            
            Program.DT_SOL[i][j][Program.NZ] = Program.DT_SOL[i][j][Program.NZ - 1];
            Program.DT_SOL[i][j][Program.KST[i][j] - 1] = Program.DT_SOL[i][j][Program.KST[i][j]];
            
            return AcUp;
        }

        //longwave downward terrestrial radiation at specific heights
        public static double Cl_LStrich()
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
                Program.L_Strich[k] = 0;
                W_Last = 0;
                while ((kk <= Program.NZ) && (W_Last < Wcrit))
                {
                    if (kk == k)
                    {
                        WMin = Program.Wrad[kk][1] + (Program.Dz1[kk] - Program.ZZ[kk + 1]) * ((Program.Wrad[kk + 1][1] - Program.Wrad[kk][1]) / (Program.ZZ[kk] - Program.ZZ[kk + 1]));
                        CGpMin = Program.CGp[kk][1] + (Program.Dz1[kk] - Program.ZZ[kk + 1]) * ((Program.CGp[kk + 1][1] - Program.CGp[kk][1]) / (Program.ZZ[kk] - Program.ZZ[kk + 1]));
                        RpMin = Program.Rp[kk][1] + (Program.Dz1[kk] - Program.ZZ[kk + 1]) * ((Program.Rp[kk + 1][1] - Program.Rp[kk][1]) / (Program.ZZ[kk] - Program.ZZ[kk + 1]));
                    }
                    else
                    {
                        WMin = Program.Wrad[kk - 1][1] + (Program.Dz1[kk] - Program.ZZ[kk - 1]) * ((Program.Wrad[kk][1] - Program.Wrad[kk - 1][1]) / (Program.ZZ[kk] - Program.ZZ[kk - 1]));
                        CGpMin = Program.CGp[kk - 1][1] + (Program.Dz1[kk] - Program.ZZ[kk - 1]) * ((Program.CGp[kk][1] - Program.CGp[kk - 1][1]) / (Program.ZZ[kk] - Program.ZZ[kk - 1]));
                        RpMin = Program.Rp[kk - 1][1] + (Program.Dz1[kk] - Program.ZZ[kk - 1]) * ((Program.Rp[kk][1] - Program.Rp[kk - 1][1]) / (Program.ZZ[kk] - Program.ZZ[kk - 1]));
                    }

                    if (kk == Program.NZ)
                    {
                        WPlus = Program.Wrad[kk - 1][1] + (Program.Dz2[kk] - Program.ZZ[kk - 1]) * ((Program.Wrad[kk][1] - Program.Wrad[kk - 1][1]) / (Program.ZZ[kk] - Program.ZZ[kk - 1]));
                        CGpPlus = Program.CGp[kk - 1][1] + (Program.Dz2[kk] - Program.ZZ[kk - 1]) * ((Program.CGp[kk][1] - Program.CGp[kk - 1][1]) / (Program.ZZ[kk] - Program.ZZ[kk - 1]));
                        RpPlus = Program.Rp[kk - 1][1] + (Program.Dz2[kk] - Program.ZZ[kk - 1]) * ((Program.Rp[kk][1] - Program.Rp[kk - 1][1]) / (Program.ZZ[kk] - Program.ZZ[kk - 1]));
                    }
                    else
                    {
                        WPlus = Program.Wrad[kk][1] + (Program.Dz2[kk] - Program.ZZ[kk]) * ((Program.Wrad[kk + 1][1] - Program.Wrad[kk][1]) / (Program.ZZ[kk + 1] - Program.ZZ[kk]));
                        CGpPlus = Program.CGp[kk][1] + (Program.Dz2[kk] - Program.ZZ[kk]) * ((Program.CGp[kk + 1][1] - Program.CGp[kk][1]) / (Program.ZZ[kk + 1] - Program.ZZ[kk]));
                        RpPlus = Program.Rp[kk][1] + (Program.Dz2[kk] - Program.ZZ[kk]) * ((Program.Rp[kk + 1][1] - Program.Rp[kk][1]) / (Program.ZZ[kk + 1] - Program.ZZ[kk]));
                    }

                    WPlus -= Program.Wrad[k][1];
                    CGpPlus -= Program.CGp[k][1];
                    RpPlus -= Program.Rp[k][1];
                    WMin -= Program.Wrad[k][1];
                    CGpMin -= Program.CGp[k][1];
                    RpMin -= Program.Rp[k][1];

                    EpsPlus = Cl_Eps(EpsPlus, RpPlus, CGpPlus, WPlus);
                    EpsMin = Cl_Eps(EpsMin, RpMin, CGpMin, WMin);

                    Program.L_Strich[k] += Program.SIGMA * Math.Pow(Program.TABS[1][1][kk], 4) * (EpsPlus - EpsMin);
                    W_Last = Program.Wrad[kk][2];
                    kk++;
                }
            }

            return WMin;
        }

        //Epsilon
        public static double Cl_Eps(double Epsi,double Rpi,double CGpi,double Wi)
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
            float rpi = (float) Rpi;

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

            Eps_CO2 = 0.185 * (1 - MathF.Exp(-0.39F * MathF.Pow((float) CGpi, 0.4F)));
            Eps_clear = 1.03 * (Eps_r + Eps_CO2);
            Eps_c = 1 - MathF.Exp((float) (-130 * Wi));

            Epsi = 1 - (1 - Eps_clear) * (1 - Eps_c);

            return Epsi;
        }

        //Epsilon upward
        public static double Cl_EpsA(int i, int j, int k, ref double EpsApi, ref double EpsAmi, ref double Eabi, ref double Tabi)
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
            
            while((kk <= Program.NZ) && (WSearch < Wcrit))
            {
                WSearch += Program.Wrad[kk][2];
                CGpSearch += Program.CGp[kk][2];
                RpSearch += Program.Rp[kk][2];
                kk++;
            }
            
            if(WSearch >= Wcrit)
            {
                if(Program.Wrad[kk-1][2] < Wcrit)
                {
                    WOutp = 0;
                }
                else
                {
                    WOutp = WSearch - Program.Wrad[kk - 1][2];
                    CGpOutp = CGpSearch - Program.CGp[kk - 1][2];
                    RpOutp = RpSearch - Program.Rp[kk - 1][2];
                }
            }
            else
            {
                WOutp = Program.Wrad[k][3] - Program.Wrad[k][2] / (Program.Dz2[k] - Program.Dz1[k]) * (Program.Dz2[k] - Program.ZZ[k]);
                CGpOutp = Program.CGp[k][3] - Program.CGp[k][2] / (Program.Dz2[k] - Program.Dz1[k]) * (Program.Dz2[k] - Program.ZZ[k]);
                RpOutp = Program.Rp[k][3] - Program.Rp[k][2] / (Program.Dz2[k] - Program.Dz1[k]) * (Program.Dz2[k] - Program.ZZ[k]);
            }
            
            WOutm = WOutp + Program.Wrad[k][2];
            CGpOutm = CGpOutp + Program.CGp[k][2];
            RpOutm = RpOutp + Program.Rp[k][2];
            
            EpsApi = Cl_Eps(EpsApi, RpOutp, CGpOutp, WOutp);
            EpsAmi = Cl_Eps(EpsAmi, RpOutm, CGpOutm, WOutm);
            
            if(WSearch>=Wcrit)
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

        //Terrestrial surface radiation
        public static double Cl_RterrG(int i, int j)
        {
            double OK = 0;

            //emissivity of the atmosphere in dependence on the cloudyness -> Master theses Manzl, Uni Innsbruck, 2010
            double EpsAtm = 0.7 * (1 - Math.Pow(Program.CLOUDS[i][j], 2.0)) + 0.90 * Math.Pow(Program.CLOUDS[i][j], 2.0);

            Program.RL[i][j] = Program.L_Strich[Program.KST[i][j] - 1] * (Program.Q[i][j] + 0.03 * MathF.Sin(Program.Alfa[i][j])) +               ////ACHTUNG Q statt 1-Q; Öttl, Dez 18
                EpsAtm * Program.SIGMA * Math.Pow(0.5 * (Program.TABS[i][j][1] + Program.TABS[i][j][10]), 4) * (1 - Program.Q[i][j]);                                                      ////ACHTUNG 1-Q statt Q; Öttl, Dez 18
            Program.RTerrG[i][j] = Program.EPSG[i][j] * (Program.RL[i][j] - Program.SIGMA * Math.Pow(Program.TG[i][j], 4));

            return OK;
        }

        //Atmospheric heating rates due to terrestrial radiation
        public static double Cl_dTdtterr(int i, int j)
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
            double OK = 0;
            double Dz, EpsBp, EpsBm, Ebel, Tbel;
            EpsBm = 0;
            EpsBp = 0;
            Ebel = 0;
            Tbel = 0;

            for (int k = Program.KST[i][j]; k <= Program.NZ;k++ )
            {
                Dz = Program.Dz2[k] - Program.Dz1[k];
                Cl_EpsB(i, j, k, ref EpsBp, ref EpsBm, ref Ebel, ref Tbel);
                
                Program.DT_TERR[i][j][k] = (1 / (Program.RHOBZZ[k] * cp * Dz)) *
                    ((Program.Eab[i][j][k] * Program.SIGMA * MathF.Pow( (float) Program.Tab[i][j][k], 4) - Program.SIGMA * MathF.Pow( (float)Program.TABS[1][1][k], 4)) * (Program.EpsAm[i][j][k] - Program.EpsAp[i][j][k]) +
                    (Ebel * Program.SIGMA * MathF.Pow((float) Tbel, 4) + (1 - Ebel) * Program.RL[i][j] - Program.SIGMA * MathF.Pow( (float)Program.TABS[1][1][k], 4)) * (EpsBp - EpsBm));
            }
            
            Program.DT_TERR[i][j][Program.KST[i][j] - 1] = Program.DT_TERR[i][j][Program.KST[i][j]];

            return OK;
        }

        //Epsilon downward
        public static double Cl_EpsB(int i, int j, int k, ref double EpsBp, ref double EpsBm, ref double Ebel, ref double Tbel)
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
            while ((kk >= (Program.KST[i][j]-1)) && (WSearch < Wcrit))
            {
                WSearch += Program.Wrad[kk][2];
                CGpSearch += Program.CGp[kk][2];
                RpSearch += Program.Rp[kk][2];
                kk--;
            }
            if (WSearch >= Wcrit)
            {
                if (Program.Wrad[kk + 1][2] < Wcrit)
                {
                    WOutp = 0;
                }
                else
                {
                    WOutp = WSearch - Program.Wrad[kk + 1][2];
                    CGpOutp = CGpSearch - Program.CGp[kk + 1][2];
                    RpOutp = RpSearch - Program.Rp[kk + 1][2];
                }
            }
            else
            {
                WOutp = WSearch;
                CGpOutp = CGpSearch;
                RpOutp = RpSearch;
            }
            WOutm = WOutp + Program.Wrad[k][2];
            CGpOutm = CGpOutp + Program.CGp[k][2];
            RpOutm = RpOutp + Program.Rp[k][2];

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
