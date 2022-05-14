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
using System.Globalization;

namespace GRAMM_2001
{
    public class INITB
    {
        /// <summary>
        /// Init specific values
        /// </summary>
        /// <param name="NI"></param>
        /// <param name="NJ"></param>
        /// <param name="NK"></param>
        public static void INIT(int NI, int NJ, int NK)
        {
            double UINIT = Program.U[1][1][NK];
            double VINIT = Program.V[1][1][NK];

            //factor used to initialize turbulent kinetic energy
            double TURBIN = 0.01;
            double TEINIT = 0;

            //initial turbulent kinetic energy
            if ((UINIT == 0) && (VINIT == 0))
            {
                TEINIT = 0.01;
            }
            else
            {
                TEINIT = TURBIN * (Math.Pow(UINIT, 2) + Math.Pow(VINIT, 2));
            }

            //turbulent eddy viscosity
            double VISINIT = 0;
            if (TEINIT > 0)
            {
                VISINIT = Program.CK * Math.Sqrt(TEINIT) * 10;
            }
            else if (TEINIT == 0)
            {
                VISINIT = Program.VISEL;
            }

            //computing specific humidity from relative humidity

            Parallel.For(1, NI + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= NJ; j++)
                {
                    float[] PBZ_L = Program.PBZ[i][j];
                    float[] FAC_L = Program.FAC[i][j];
                    float[] FACT_L = Program.FACTOR[i][j];
                    double[] QBZ_L = Program.QBZ[i][j];

                    for (int k = 1; k <= NK; k++)
                    {
                        double TBZN = Math.Round(Program.TABS[i][j][k] - 153.15, 3);
                        //                        if (Convert.ToInt32(TBZN) > 209) TBZN = 209;
                        //                        Int32 TBZNINT = Convert.ToInt32(TBZN);
                        TBZN = Math.Max(0, Math.Min(209.999, TBZN));
                        int TBZNINT = Convert.ToInt32(Math.Floor(TBZN));
                        double PDST = Program.PSAT[TBZNINT + 1] + (Program.PSAT[TBZNINT + 2] - Program.PSAT[TBZNINT + 1]) * (TBZN - (float)TBZNINT);

                        // initial specific humidity
                        QBZ_L[k] = 18.02F / 28.96F * PDST / (PBZ_L[k] / Program.QUINIT - PDST) * 1000;

                        if (Program.AKLA == 7 && Program.METEO.ToUpper() == "Y" && Program.ISTAT == 0)
                        {
                            //if (Program.AHImm[i][j] - 30 < Program.AH_Bassins[i][j])
                            if (Program.TPI[i][j] == 1 || Program.TPI[i][j] == 2 || Program.TPI[i][j] == 4 || Program.TPI[i][j] == 5)
                            {
                                if (Program.Wind_Velocity >= 0.35F) // lower humidity at the ground and above the inversion height
                                {
                                    QBZ_L[k] *= 0.6;
                                }
                                else // at low wind velocities: higher humidity at the ground 
                                {
                                    QBZ_L[k] = 18.02F / 28.96F * PDST / (PBZ_L[k] / 0.95 - PDST) * 1000;
                                }
                            }
                            else if ((Program.ZSPImm[i][j][k] - Program.AH_Bassins[i][j]) > Program.Inversion_Height)
                            {
                                QBZ_L[k] *= 0.4;
                            }
                        }

                        FAC_L[k] = FACT_L[k];
                    }
                }
            });

            //temperature value used to improve numerical accuracy of the solution algorithm for temperature
            Program.TBZ1 = Program.TBZ[2][2][1];

            Parallel.For(1, NI + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= NJ; j++)
                {
                    double[] T_L = Program.T[i][j];
                    double[] TE_L = Program.TE[i][j];
                    double[] TN_L = Program.TN[i][j];
                    double[] TEN_L = Program.TEN[i][j];
                    double[] DISS_L = Program.DISS[i][j];
                    double[] DISSN_L = Program.DISSN[i][j];
                    double[] VISH_L = Program.VISH[i][j];
                    double[] VISV_L = Program.VISV[i][j];

                    for (int k = 1; k <= NK; k++)
                    {
                        //specific humidity - actual time step
                        Program.QU[i][j][k] = Program.QBZ[i][j][k];

                        //specific humidity - next time step
                        Program.QUN[i][j][k] = Program.QBZ[i][j][k];

                        //reference temperature is subtracted to improve numerical accuracy
                        T_L[k] -= Program.TBZ1;
                        TN_L[k] = T_L[k];

                        //turbulent kinetic energy
                        TE_L[k] = 0.01;
                        TEN_L[k] = 0.01;

                        //dissipation
                        DISS_L[k] = 0.0001F;
                        DISSN_L[k] = 0.0001F;

                        //turbulent eddy viscosity
                        VISH_L[k] = VISINIT;
                        VISV_L[k] = VISINIT;
                    }
                }
            });

            //Values at south and north borders
            Parallel.For(1, NI + 1, Program.pOptions, i =>
            {
                for (int k = 1; k <= NK; k++)
                {
                    int j = 1;
                    //large-scale values at the southern border for the next time step
                    Program.TSSN[i][k] = Program.T[i][j][k];
                    //large-scale values at the southern border for the actual time step
                    Program.TSS[i][k] = Program.T[i][j][k];

                    Program.QUSSN[i][k] = Program.QBZ[i][j][k];
                    Program.QUSS[i][k] = Program.QBZ[i][j][k];

                    Program.WSSN[i][k] = 0;
                    Program.WSS[i][k] = 0;

                    j = NJ;
                    Program.TSNN[i][k] = Program.T[i][j][k];
                    Program.TSN[i][k] = Program.T[i][j][k];

                    Program.QUSNN[i][k] = Program.QBZ[i][j][k];
                    Program.QUSN[i][k] = Program.QBZ[i][j][k];

                    Program.WSNN[i][k] = 0;
                    Program.WSN[i][k] = 0;
                }
            });

            //Values at west and east borders
            Parallel.For(1, NJ + 1, Program.pOptions, j =>
            {
                for (int k = 1; k <= NK; k++)
                {
                    int i = 1;
                    Program.TSWN[j][k] = Program.T[i][j][k];
                    Program.TSW[j][k] = Program.T[i][j][k];

                    Program.QUSWN[j][k] = Program.QBZ[i][j][k];
                    Program.QUSW[j][k] = Program.QBZ[i][j][k];

                    Program.WSWN[j][k] = 0;
                    Program.WSW[j][k] = 0;

                    i = NI;
                    Program.TSEN[j][k] = Program.T[i][j][k];
                    Program.TSE[j][k] = Program.T[i][j][k];

                    Program.QUSEN[j][k] = Program.QBZ[i][j][k];
                    Program.QUSE[j][k] = Program.QBZ[i][j][k];

                    Program.WSEN[j][k] = 0;
                    Program.WSE[j][k] = 0;
                }
            });

            CustomAirSoilInit cInit = new CustomAirSoilInit(Program.TINIT, Program.TBINIT, Program.TBINIT1, Program.QUINIT);
            if (cInit.SnowHeightThreshold < 100000)
            {
                Console.WriteLine(" Snow altitude threshold: {0} m", cInit.SnowHeightThreshold);
            }
            if (Math.Abs(cInit.DeltaTSurfaceWater) > 0)
            {
                Console.WriteLine(" Temperature difference to ground temperature: {0} °C", cInit.DeltaTSurfaceWater);
            }

            //boundary values for surface parameters
            Parallel.For(1, NI + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= NJ; j++)
                {
                    //friction velocity
                    Program.UST[i][j] = 0;
                    if ((Program.METEO == "Y") || (Program.METEO == "y"))
                    {
                        Program.UST[i][j] = 0.15;
                    }

                    Program.USTV[i][j] = 0;
                    if ((Program.METEO == "Y") || (Program.METEO == "y"))
                    {
                        Program.USTV[i][j] = 0.15 * Program.CK / Math.Log((Program.ZSPImm[3][3][1] - Program.AHImm[3][3]) / Program.Rauigkeit);
                    }

                    //initial values for stability classes
                    if ((Program.METEO == "Y") || (Program.METEO == "y"))
                    {
                        Program.stabilityclass[i][j] = Program.AKLA;
                    }

                    //characteristic potential temperature
                    Program.TST[i][j] = 0;

                    //Obukhov length
                    Program.OL[i][j] = 10;

                    Program.XWQ[i][j] = 0;
                    Program.WQU[i][j] = 0;
                    Program.RSOLG[i][j] = 0;
                    Program.GLOBRAD[i][j] = 0;
                    Program.RL[i][j] = 0;
                    Program.AWQ[i][j] = 0;
                    Program.QUG[i][j] = Program.QU[i][j][1];

                    //soil parameters
                    double DUMMY = 1;
                    for (int kb = 1; kb <= Program.NZB; kb++)
                    {
                        //surface temperature
                        /*
                         * The surface temperature is initialized in dependence on the height above sea level and the latitude (ÖTTL, 22 JULY 2015)
                        */
                        DUMMY -= Program.DWB[kb];
                        double tbInit = Program.TBINIT;
                        double tbInitDepth = Program.TBINIT1;
                        
                        //in case of snow covered terrain
                        if (Program.AHImm[i][j] > cInit.SnowHeightThreshold)
                        {
                            tbInit = Math.Min(270, Program.TBINIT);
                            tbInitDepth = Math.Min(272, Program.TBINIT1);
                        }
                        else
                        {
                            tbInit = Program.TBINIT;
                            tbInitDepth = Program.TBINIT1;
                        }
                        
                        double DELTA_TB = cInit.TAir2m - tbInit; //this is the desired temperature difference between the surface and the air temperature at the ground
                        //Dependency on height above sea level -> note that the dependency on latitude is taken into account in the calculation of TBINIT1 (tempint.cs)
                        Program.TB[i][j][kb] = (tbInitDepth + cInit.SoilTempGradient * Program.AHImm[i][j]) + (Program.TABS[i][j][1] - DELTA_TB - (tbInitDepth + cInit.SoilTempGradient * Program.AHImm[i][j])) * DUMMY;

                        /*
                         * !Oettl, 7.Juni 2013 - Versuch die oberste Schicht der Erdbodentemperatur in Abhängigkeit von der AKLA abweichend zur Lufttemperatur darüber einzustellen
                         * !damit der sensible Wärmefluss zu Beginn der Berechnung sofort hoch genug ist
                         * !da bei berechnungen für Voitsberg und Mürzzuschlag zu starke kaltluftabflüsse berechnet wurden, am 3.9.13 wieder zurückgenommen
                        */
                        if ((Program.METEO == "Y") || (Program.METEO == "y"))
                        {
                            if (kb < Program.NZB)
                            {
                                if (Program.AKLA == 7)
                                {
                                    Program.TB[i][j][kb] = (tbInitDepth + cInit.SoilTempGradient * Program.AHImm[i][j]) + (Program.TABS[Program.AHMINI][Program.AHMINJ][1] - (tbInitDepth)) * DUMMY - 2;
                                }
                                if (Program.AKLA == 6)
                                {
                                    Program.TB[i][j][kb] = (tbInitDepth + cInit.SoilTempGradient * Program.AHImm[i][j]) + (Program.TABS[i][j][1] - (tbInitDepth + cInit.SoilTempGradient * Program.AHImm[i][j])) * DUMMY;
                                }
                                if (Program.AKLA == 5)
                                {
                                    Program.TB[i][j][kb] = (tbInitDepth + cInit.SoilTempGradient * Program.AHImm[i][j]) + (Program.TABS[i][j][1] - (tbInitDepth + cInit.SoilTempGradient * Program.AHImm[i][j])) * DUMMY;
                                }

                                /*
                                if (Program.AKLA > 2)
                                {
                                    //CORINE CLASS 111 CONTINUOUS URBAN FABRIC
                                    if (Program.Z0[i][j] >= 1.2)
                                    {
                                        Program.TB[i][j][kb] += 1 * (Program.NZB - kb) / ((float) Program.NZB);
                                    }
                                    //CORINE CLASS 112 DISCONTINUOUS URBAN FABRIC
                                    if ((Program.Z0[i][j] == 0.5F) && (Program.ALAMBDA[i][j] == 4F))
                                    {
                                        Program.TB[i][j][kb] += 0.5F * (Program.NZB - kb) / ((float) Program.NZB);
                                    }
                                    //CORINE CLASS 121 INUDSTRIAL UNITS
                                    if ((Program.Z0[i][j] == 0.5F) && (Program.ALAMBDA[i][j] == 3.5F))
                                    {
                                        Program.TB[i][j][kb] += 0.2F * (Program.NZB - kb) / ((float) Program.NZB);
                                    }
                                }
                                */

                                //surface temperature should be less or equal the air temperature at ground at all heights
                                // if (Program.AKLA == 7)
                                // {
                                //     Program.TB[i][j][kb] = Math.Min(Program.TB[i][j][kb], Program.TABS[i][j][1] - 0.5F);

                                //     // avoid high wind velocioties at layer 1 and 2
                                //     if (Program.Wind_Velocity < 0.35F)
                                //     {
                                //         Program.TB[i][j][kb] = Math.Max(Program.TB[i][j][kb], Program.TABS[i][j][1] - 5);
                                //     }
                                // }
                                // else
                                // {
                                //     Program.TB[i][j][kb] = Math.Min(Program.TB[i][j][kb], Program.TABS[i][j][1]);
                                // }
                            }
                            else
                            {
                                Program.TB[i][j][kb] = tbInitDepth - 0.005 * Program.AHImm[i][j];
                            }
                        }

                        Program.TG[i][j] = Program.TB[i][j][2];

                        Program.TBN[i][j][kb] = Program.TB[i][j][kb] - 0.1F;
                        Program.TBA[i][j][kb] = Program.TB[i][j][kb] + 0.1F;
                    }
                }
            });

            //initial velocity fields
            Parallel.For(1, NI + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= NJ; j++)
                {
                    double[] U_L = Program.U[i][j];
                    double[] V_L = Program.V[i][j];
                    double[] UG_L = Program.UG[i][j];

                    for (int k = NK; k >= 1; k--)
                    {
                        Program.U1NRHO[i][j][k] = (float)(U_L[k] * Program.RHO[i][j][k]);
                        Program.U2NRHO[i][j][k] = (float)(U_L[k] * Program.RHO[i][j][k]);
                        Program.U1[i][j][k] = U_L[k];
                        Program.U2[i][j][k] = U_L[k];
                        Program.U1N[i][j][k] = U_L[k];
                        Program.U2N[i][j][k] = U_L[k];
                        Program.UN[i][j][k] = U_L[k];
                        if (Program.UN[i][j][k] > 50)
                        {
                            Console.WriteLine("U-component: " + Convert.ToString(Math.Round(Program.UN[i][j][k], 1)) + "m/s");
                        }

                        //Geostrophic wind estimation
                        if (k < NK)
                        {
                            if (Math.Sqrt(Math.Pow(U_L[k], 2) + Math.Pow(V_L[k], 2)) > Math.Sqrt(Math.Pow(UG_L[k + 1], 2) + Math.Pow(Program.VG[i][j][k + 1], 2)))
                            {
                                UG_L[k] = U_L[k];
                                Program.UG1[i][j][k] = (float)(U_L[k]);
                                Program.UG2[i][j][k] = (float)(U_L[k]);
                                Program.VG[i][j][k] = (float)(V_L[k]);
                                Program.VG1[i][j][k] = (float)(V_L[k]);
                                Program.VG2[i][j][k] = (float)(V_L[k]);
                            }
                            else
                            {
                                UG_L[k] = UG_L[k + 1];
                                Program.UG1[i][j][k] = (float)(UG_L[k + 1]);
                                Program.UG2[i][j][k] = (float)(UG_L[k + 1]);
                                Program.VG[i][j][k] = (float)(Program.VG[i][j][k + 1]);
                                Program.VG1[i][j][k] = (float)(Program.VG[i][j][k + 1]);
                                Program.VG2[i][j][k] = (float)(Program.VG[i][j][k + 1]);
                            }
                        }
                        else
                        {
                            UG_L[k] = U_L[k];
                            Program.UG1[i][j][k] = (float)(U_L[k]);
                            Program.UG2[i][j][k] = (float)(U_L[k]);
                            Program.VG[i][j][k] = (float)(V_L[k]);
                            Program.VG1[i][j][k] = (float)(V_L[k]);
                            Program.VG2[i][j][k] = (float)(V_L[k]);
                        }

                        /*
                        if ((Program.METEO == "Y") || (Program.METEO == "y"))
                        {
                            UG_L[k] = U_L[NK];
                            Program.UG1[i][j][k] = U_L[NK];
                            Program.UG2[i][j][k] = U_L[NK];
                        }
                        else
                        {
                            UG_L[k] = U_L[k];
                            Program.UG1[i][j][k] = U_L[k];
                            Program.UG2[i][j][k] = U_L[k];
                        }
                        */

                        //large-scale values at the boundaries
                        Program.USN[i][k] = Program.U[i][NJ][k];
                        Program.USNN[i][k] = Program.U[i][NJ][k];

                        Program.USS[i][k] = Program.U[i][1][k];
                        Program.USSN[i][k] = Program.U[i][1][k];

                        if (i == 1) // avoid deadlock
                        {
                            Program.USE1[j][k] = Program.U[NI][j][k];
                            Program.USEN[j][k] = Program.U[NI][j][k];
                            Program.USW[j][k] = Program.U[1][j][k];
                            Program.USWN[j][k] = Program.U[1][j][k];
                        }


                        //set wind speeds initially to zero
                        /*
                        if (((Program.METEO == "Y") || (Program.METEO == "y")) && Program.ISTAT == 1)
                        {
                            Program.U1NRHO[i][j][k] = 0;
                            Program.U2NRHO[i][j][k] = 0;
                            Program.U1[i][j][k] = 0;
                            Program.U2[i][j][k] = 0;
                            Program.U1N[i][j][k] = 0;
                            Program.U2N[i][j][k] = 0;
                            Program.UN[i][j][k] = 0;
                        }
                        */

                    }
                }
            });

            Parallel.For(1, NI + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= NJ; j++)
                {

                    double[] V_L = Program.V[i][j];


                    for (int k = NK; k >= 1; k--)
                    {
                        Program.V1NRHO[i][j][k] = (float)(V_L[k] * Program.RHO[i][j][k]);
                        Program.V2NRHO[i][j][k] = (float)(V_L[k] * Program.RHO[i][j][k]);
                        Program.V1[i][j][k] = V_L[k];
                        Program.V2[i][j][k] = V_L[k];
                        Program.V1N[i][j][k] = V_L[k];
                        Program.V2N[i][j][k] = V_L[k];
                        Program.VN[i][j][k] = V_L[k];

                        //Geostrophic balance is assumed when verical profiles are used as met-input
                        /*
                        if ((Program.METEO == "Y") || (Program.METEO == "y"))
                        {
                            Program.VG[i][j][k] = V_L[NK];
                            Program.VG1[i][j][k] = V_L[NK];
                            Program.VG2[i][j][k] = V_L[NK];
                        }
                        else
                        {
                            Program.VG[i][j][k] = V_L[k];
                            Program.VG1[i][j][k] = V_L[k];
                            Program.VG2[i][j][k] = V_L[k];
                        }     
                        */

                        //large-scale values at the boundaries
                        Program.VSN[i][k] = Program.V[i][NJ][k];
                        Program.VSNN[i][k] = Program.V[i][NJ][k];

                        Program.VSS[i][k] = Program.V[i][1][k];
                        Program.VSSN[i][k] = Program.V[i][1][k];

                        if (i == 1) // avoid deadlock
                        {
                            Program.VSE[j][k] = Program.V[NI][j][k];
                            Program.VSEN[j][k] = Program.V[NI][j][k];

                            Program.VSW[j][k] = Program.V[1][j][k];
                            Program.VSWN[j][k] = Program.V[1][j][k];
                        }


                        //set wind speeds initially to zero
                        /*
                        if (((Program.METEO == "Y") || (Program.METEO == "y")) && Program.ISTAT == 1)
                        {
                            Program.V1NRHO[i][j][k] = 0;
                            Program.V2NRHO[i][j][k] = 0;
                            Program.V1[i][j][k] = 0;
                            Program.V2[i][j][k] = 0;
                            Program.V1N[i][j][k] = 0;
                            Program.V2N[i][j][k] = 0;
                            Program.VN[i][j][k] = 0;
                        }
                        */

                    }
                }
            });

            //reading land-use data
            double DUMMY1 = 0;
            double DUMMY2 = 0;
            double DUMMY3 = 0;
            double DUMMY4 = 0;
            double DUMMY5 = 0;
            double DUMMY6 = 0;
            if ((Program.METEO == "Y") || (Program.METEO == "y"))
            {
                DUMMY1 = 0.1;
                DUMMY3 = 0.15;
                DUMMY4 = 0.000001;
                DUMMY5 = 2;
                DUMMY6 = 0.9;
            }
            
            if (File.Exists("landuse.asc"))
            {
                string[] text = new string[(NI + 2) * (NJ + 2)];
                try
                {
                    using (FileStream fs = new FileStream("landuse.asc", FileMode.Open, FileAccess.Read, FileShare.Read))
                    {
                        using (StreamReader r = new StreamReader(fs))
                        {
                            text = Convert.ToString(r.ReadLine()).Split(new char[] { ' ', ',', ';' }, StringSplitOptions.RemoveEmptyEntries);
                            int n = 0;
                            for (int j = 1; j < NJ + 1; j++)
                            {
                                for (int i = 1; i < NI + 1; i++)
                                {
                                    Program.RHOB[i][j] = Convert.ToSingle(text[n].Replace(".", Program.decsep));
                                    n++;
                                    // in case of snow, set soil density to snow value
                                    if (Program.AHImm[i][j] > cInit.SnowHeightThreshold)
                                    {
                                       Program.RHOB[i][j] = Math.Max(1 / 0.0000005F / 900, Program.RHOB[i][j]);
                                    }
                                }
                            }
                            text = Convert.ToString(r.ReadLine()).Split(new char[] { ' ', ',', ';' }, StringSplitOptions.RemoveEmptyEntries);
                            n = 0;
                            for (int j = 1; j < NJ + 1; j++)
                            {
                                for (int i = 1; i < NI + 1; i++)
                                {
                                    Program.ALAMBDA[i][j] = Convert.ToSingle(text[n].Replace(".", Program.decsep));
                                    n++;
                                    // in case of snow, set low surface albedo to snow value
                                    if (Program.AHImm[i][j] > cInit.SnowHeightThreshold && Program.ALAMBDA[i][j] < 1)
                                    {
                                        Program.ALAMBDA[i][j] = 1F;
                                    }
                                }
                            }
                            text = Convert.ToString(r.ReadLine()).Split(new char[] { ' ', ',', ';' }, StringSplitOptions.RemoveEmptyEntries);
                            n = 0;
                            for (int j = 1; j < NJ + 1; j++)
                            {
                                for (int i = 1; i < NI + 1; i++)
                                {
                                    Program.Z0[i][j] = Convert.ToSingle(text[n].Replace(".", Program.decsep));
                                    n++;
                                    // in case of snow cover, set green areas to a lower roughness lenght
                                    if (Program.AHImm[i][j] > cInit.SnowHeightThreshold && Program.Z0[i][j] < 0.4F)
                                    {
                                       Program.Z0[i][j] = Math.Min(Program.Z0[i][j], 0.01F);
                                    }
                                }
                            }
                            text = Convert.ToString(r.ReadLine()).Split(new char[] { ' ', ',', ';' }, StringSplitOptions.RemoveEmptyEntries);
                            n = 0;
                            for (int j = 1; j < NJ + 1; j++)
                            {
                                for (int i = 1; i < NI + 1; i++)
                                {
                                    Program.FW[i][j] = Convert.ToSingle(text[n].Replace(".", Program.decsep));
                                    n++;
                                    // in case of snow cover, set soil moisture to snow
                                    if (Program.AHImm[i][j] > cInit.SnowHeightThreshold && Program.FW[i][j] < 0.1F)
                                    {
                                        Program.FW[i][j] = 0.1F;
                                    }
                                }
                            }
                            text = Convert.ToString(r.ReadLine()).Split(new char[] { ' ', ',', ';' }, StringSplitOptions.RemoveEmptyEntries);
                            n = 0;
                            for (int j = 1; j < NJ + 1; j++)
                            {
                                for (int i = 1; i < NI + 1; i++)
                                {
                                    Program.EPSG[i][j] = Convert.ToDouble(text[n].Replace(".", Program.decsep));
                                    n++;
                                    // in case of snow cover, set surface emissivity to snow
                                    if (Program.AHImm[i][j] > cInit.SnowHeightThreshold)
                                    {
                                        Program.EPSG[i][j] = 0.95F;
                                    }
                                }
                            }
                            text = Convert.ToString(r.ReadLine()).Split(new char[] { ' ', ',', ';' }, StringSplitOptions.RemoveEmptyEntries);
                            n = 0;
                            for (int j = 1; j < NJ + 1; j++)
                            {
                                for (int i = 1; i < NI + 1; i++)
                                {
                                    Program.ALBEDO[i][j] = Convert.ToDouble(text[n].Replace(".", Program.decsep));
                                    n++;
                                    // in case of snow cover, set surface albedo to snow
                                    if (Program.AHImm[i][j] > cInit.SnowHeightThreshold)
                                    {
                                        Program.ALBEDO[i][j] = 0.6F;
                                    }
                                }
                            }

                            //clouds are set to zero in the current version
                            for (int j = 1; j < NJ + 1; j++)
                            {
                                for (int i = 1; i < NI + 1; i++)
                                {
                                    Program.CLOUDS[i][j] = 0;
                                    //Program.SNOW[i][j] = 0;
                                }
                            }

                            //diagnostic evaluation of anthropogenic heat fluxes
                            /*
                            for (int j = 1; j < NJ + 1; j++)
                            {
                                for (int i = 1; i < NI + 1; i++)
                                {
                                    //CORINE CLASS 111 CONTINUOUS URBAN FABRIC
                                    if (Program.Z0[i][j] >= 1.2)
                                        Program.AWQ[i][j] = 40 * Program.DDXImm[i] * Program.DDYImm[j];
                                    //CORINE CLASS 112 DISCONTINUOUS URBAN FABRIC
                                    if ((Program.Z0[i][j] == 0.5F) && (Program.ALAMBDA[i][j] == 4))
                                        Program.AWQ[i][j] = 20 * Program.DDXImm[i] * Program.DDYImm[j];
                                    //CORINE CLASS 121 INUDSTRIAL UNITS
                                    if ((Program.Z0[i][j] == 0.5F) && (Program.ALAMBDA[i][j] == 3.5))
                                        Program.AWQ[i][j] = 50 * Program.DDXImm[i] * Program.DDYImm[j];
                                }
                            }
                            */
                        }
                    }
                }
                catch
                {
                    Console.WriteLine("Error when reading file 'landuse.asc' - Execution stopped");
                    Environment.Exit(0);
                }
            }
            else
            {
                if ((Program.METEO == "Y") || (Program.METEO == "y"))
                { }
                else
                {
                    Console.WriteLine();
                    Console.WriteLine(" FILE 'landuse.asc' DOES NOT EXIST!");
                    Console.WriteLine(" PLEASE PROVIDE SOIL DATA FOR YOUR ");
                    Console.WriteLine(" MODEL AREA ...");
                    Console.WriteLine();
                    Console.WriteLine();
                    Console.WriteLine(" ***********************************");
                    Console.WriteLine(" *             ALBEDO              *");
                    Console.WriteLine(" ***********************************");
                    Console.WriteLine();
                    Console.WriteLine(" EXAMPLE VALUES: ");
                    Console.WriteLine();
                    Console.WriteLine("      CITY               0.20");
                    Console.WriteLine("      SUBURB             0.23");
                    Console.WriteLine("      FOREST             0.10");
                    Console.WriteLine("      FARMLAND           0.22");
                    Console.WriteLine("      SPARSE VEGETATION  0.22");
                    Console.WriteLine("      WATER              0.08");
                    Console.WriteLine();
                    Console.Write(" ALBEDO = ");
                    DUMMY1 = Convert.ToDouble(Console.ReadLine().Replace(".", Program.decsep));

                    Console.WriteLine();
                    Console.WriteLine();
                    Console.WriteLine(" ***********************************");
                    Console.WriteLine(" *      AERODYNAMIC ROUGHNESS      *");
                    Console.WriteLine(" ***********************************");
                    Console.WriteLine();
                    Console.WriteLine(" EXAMPLE VALUES: ");
                    Console.WriteLine();
                    Console.WriteLine("      CITY               0.80");
                    Console.WriteLine("      SUBURB             0.50");
                    Console.WriteLine("      FOREST             0.50");
                    Console.WriteLine("      FARMLAND           0.05");
                    Console.WriteLine("      SPARSE VEGETATION  0.01");
                    Console.WriteLine("      WATER              0.0001F");
                    Console.WriteLine();
                    Console.Write(" AERODYNAMIC ROUGHNESS = ");
                    DUMMY2 = Convert.ToDouble(Console.ReadLine().Replace(".", Program.decsep));

                    Console.WriteLine();
                    Console.WriteLine();
                    Console.WriteLine(" ***********************************");
                    Console.WriteLine(" *        MOISTURE PARAMETER       *");
                    Console.WriteLine(" ***********************************");
                    Console.WriteLine();
                    Console.WriteLine(" EXAMPLE VALUES: ");
                    Console.WriteLine();
                    Console.WriteLine("      CITY               0.05");
                    Console.WriteLine("      SUBURB             0.10");
                    Console.WriteLine("      FOREST             0.20");
                    Console.WriteLine("      FARMLAND           0.15");
                    Console.WriteLine("      SPARSE VEGETATION  0.05");
                    Console.WriteLine("      WATER              1.00");
                    Console.WriteLine();
                    Console.Write(" MOISTURE PARAMETER = ");
                    DUMMY3 = Convert.ToDouble(Console.ReadLine().Replace(".", Program.decsep));

                    Console.WriteLine();
                    Console.WriteLine();
                    Console.WriteLine(" ***************************************");
                    Console.WriteLine(" * CONDUCTIBILITY OF TEMPERATURE M^2/S *");
                    Console.WriteLine(" ***************************************");
                    Console.WriteLine();
                    Console.WriteLine(" EXAMPLE VALUES: ");
                    Console.WriteLine();
                    Console.WriteLine("      CITY               2.0E-6");
                    Console.WriteLine("      SUBURB             1.3E-6");
                    Console.WriteLine("      FOREST             0.8E-6");
                    Console.WriteLine("      FARMLAND           0.7E-6");
                    Console.WriteLine("      SPARSE VEGETATION  1.0E-6");
                    Console.WriteLine("      WATER              1.0E-6");
                    Console.WriteLine();
                    Console.Write(" CONDUCTIBILITY OF TEMPERATURE = ");
                    DUMMY4 = Convert.ToDouble(Console.ReadLine().Replace(".", Program.decsep));

                    Console.WriteLine();
                    Console.WriteLine();
                    Console.WriteLine(" *************************************");
                    Console.WriteLine(" *    THERMAL CONDUCTIVITY W/(MK)    *");
                    Console.WriteLine(" *************************************");
                    Console.WriteLine();
                    Console.WriteLine(" EXAMPLE VALUES: ");
                    Console.WriteLine();
                    Console.WriteLine("      CITY               4.68");
                    Console.WriteLine("      SUBURB             2.86");
                    Console.WriteLine("      FOREST             0.94");
                    Console.WriteLine("      FARMLAND           2.00");
                    Console.WriteLine("      SPARSE VEGETATION  2.68");
                    Console.WriteLine("      WATER            100.00");
                    Console.WriteLine();
                    Console.Write(" THERMAL CONDUCTIVITY = ");
                    DUMMY5 = Convert.ToDouble(Console.ReadLine().Replace(".", Program.decsep));

                    Console.WriteLine();
                    Console.WriteLine();
                    Console.WriteLine(" ***********************************");
                    Console.WriteLine(" *       EMISSIVITY OF SOIL        *");
                    Console.WriteLine(" ***********************************");
                    Console.WriteLine();
                    Console.WriteLine(" EXAMPLE VALUES: ");
                    Console.WriteLine();
                    Console.WriteLine("      CITY               0.88");
                    Console.WriteLine("      SUBURB             0.88");
                    Console.WriteLine("      FOREST             0.92");
                    Console.WriteLine("      FARMLAND           0.90");
                    Console.WriteLine("      SPARSE VEGETATION  0.90");
                    Console.WriteLine("      WATER              0.98");
                    Console.WriteLine();
                    Console.Write(" EMISSIVITY = ");
                    DUMMY6 = Convert.ToDouble(Console.ReadLine().Replace(".", Program.decsep));
                }
                for (int j = 1; j < NJ + 1; j++)
                {
                    for (int i = 1; i < NI + 1; i++)
                    {
                        Program.ALBEDO[i][j] = DUMMY1;
                        Program.FW[i][j] = (float)DUMMY3;
                        Program.Z0[i][j] = (float)(DUMMY2);
                        if ((Program.METEO == "Y") || (Program.METEO == "y"))
                        {
                            Program.Z0[i][j] = (float)(Program.Rauigkeit);
                        }

                        Program.EPSG[i][j] = DUMMY6;
                        Program.ALAMBDA[i][j] = (float)(DUMMY5);
                        Program.RHOB[i][j] = (float)(Program.ALAMBDA[i][j] / DUMMY4 / Program.CPBOD);
                    }
                }
            }

            //set water temperatures for water bodies
            if ((Program.METEO == "Y") || (Program.METEO == "y"))
            {
                Parallel.For(1, NI + 1, Program.pOptions, i =>
                {
                    for (int j = 1; j <= NJ; j++)
                    {
                        if (Program.FW[i][j] > 0.99)
                        {
                            int k = 1;
                            double DUMMY = 1;
                            for (int kb = 1; kb <= Program.NZB; kb++)
                            {
                                // user defined deltaT for water bodies
                                if (Math.Abs(cInit.DeltaTSurfaceWater) > 0.2)
                                {
                                    Program.TB[i][j][kb] = Math.Max(270, (Program.TBINIT1 + cInit.SoilTempGradient * Program.AHImm[i][j]) + cInit.DeltaTSurfaceWater);
                                }
                                else
                                {   // default deltaT for water bodies
                                    DUMMY -= Program.DWB[kb];
                                    if (kb < Program.NZB)
                                    {
                                        if (Program.AKLA == 7)
                                        {
                                            Program.TB[i][j][kb] = (Program.TBINIT1 + cInit.SoilTempGradient * Program.AHImm[i][j]) + (Program.TABS[i][j][1] + 10 - (Program.TBINIT1 + cInit.SoilTempGradient * Program.AHImm[i][j])) * DUMMY;
                                        }

                                        if (Program.AKLA == 6)
                                        {
                                            Program.TB[i][j][kb] = (Program.TBINIT1 + cInit.SoilTempGradient * Program.AHImm[i][j]) + (Program.TABS[i][j][1] + 5 - (Program.TBINIT1 + cInit.SoilTempGradient * Program.AHImm[i][j])) * DUMMY;
                                        }

                                        if (Program.AKLA == 5)
                                        {
                                            Program.TB[i][j][kb] = (Program.TBINIT1 + cInit.SoilTempGradient * Program.AHImm[i][j]) + (Program.TABS[i][j][1] + 2 - (Program.TBINIT1 + cInit.SoilTempGradient * Program.AHImm[i][j])) * DUMMY;
                                        }

                                        if (Program.AKLA == 1)
                                        {
                                            Program.TB[i][j][kb] = (Program.TBINIT1 + cInit.SoilTempGradient * Program.AHImm[i][j]) + (Program.TABS[i][j][1] - 10 - (Program.TBINIT1 + cInit.SoilTempGradient * Program.AHImm[i][j])) * DUMMY;
                                        }

                                        if (Program.AKLA == 2)
                                        {
                                            Program.TB[i][j][kb] = (Program.TBINIT1 + cInit.SoilTempGradient * Program.AHImm[i][j]) + (Program.TABS[i][j][1] - 5 - (Program.TBINIT1 + cInit.SoilTempGradient * Program.AHImm[i][j])) * DUMMY;
                                        }

                                        if (Program.AKLA == 3)
                                        {
                                            Program.TB[i][j][kb] = (Program.TBINIT1 + cInit.SoilTempGradient * Program.AHImm[i][j]) + (Program.TABS[i][j][1] - 2 - (Program.TBINIT1 + cInit.SoilTempGradient * Program.AHImm[i][j])) * DUMMY;
                                        }
                                    }
                                    else
                                    {
                                        Program.TB[i][j][kb] = Program.TBINIT1 + cInit.SoilTempGradient * Program.AHImm[i][j];
                                    }
                                }
                                Program.TG[i][j] = Program.TB[i][j][2];
                                Program.TBN[i][j][kb] = Program.TB[i][j][kb];
                                Program.TBA[i][j][kb] = Program.TB[i][j][kb];
                            }
                        }
                    }
                });
            }

            Int32 IZELL = 0;
            double ARSUM = 0;
            double VOLSUM = 0;
            for (int i = 2; i < NI; i++)
            {
                for (int j = 2; j < NJ; j++)
                {
                    ARSUM += Program.DDXImm[i] * Program.DDYImm[j] * 0.000001;
                    for (int k = 1; k < NK; k++)
                    {
                        IZELL++;
                        VOLSUM += Program.VOLImm[i][j][k] * 0.000000001;
                    }
                }
            }
            WriteZ0ToFile();
            Console.WriteLine(" IINIT : NUMBER OF CALCULATED CELLS    : " + IZELL.ToString());
            Console.WriteLine(" IINIT : SIZE OF CALCULATED AREA       : " + ARSUM.ToString("0") + " km**2");
            Console.WriteLine(" IINIT : SIZE OF CALCULATED VOLUME     : " + VOLSUM.ToString("0") + " km**3");
        }

        /// <summary>
        /// Write used Z0 values to Roughness.txt
        /// </summary>
        private static void WriteZ0ToFile()
        {
            try
            {
                using (StreamWriter wt = new StreamWriter("RoughnessUsed.txt"))
                {
                    wt.WriteLine("ncols         " + Program.NX.ToString(CultureInfo.InvariantCulture));
                    wt.WriteLine("nrows         " + Program.NY.ToString(CultureInfo.InvariantCulture));
                    wt.WriteLine("xllcorner     " + Program.IKOOA.ToString(CultureInfo.InvariantCulture));
                    wt.WriteLine("yllcorner     " + Program.JKOOA.ToString(CultureInfo.InvariantCulture));
                    wt.WriteLine("cellsize      " + Program.DDXImm[1].ToString(CultureInfo.InvariantCulture));
                    wt.WriteLine("NODATA_value  " + "-9999");

                    for (int y = Program.NY; y >= 1; y--)
                    {
                        for (int x = 1; x <= Program.NX; x++)
                        {
                            wt.Write(Math.Round(Program.Z0[x][y], 5).ToString(CultureInfo.InvariantCulture) + " ");
                        }
                        wt.WriteLine();
                    }
                }
            }
            catch
            { }
        }

        /// <summary>
        /// Read the snow file and set snow to 1 if the height is exceeded<br></br>
        /// <b>Returns user defined</b> (snowHeightThreshold)(airTemp at 2 m)(surfaceTemp)(surfaceTemp 1m depth)(relative humidity)
        /// </summary>
        public static (double, double, double, double, double, double, double) ReadCustomInit()
        {
            double snowHeight = 1000000;
            double tAir2m = Program.TINIT;
            double tSurface = Program.TBINIT;
            double tSurface1m = Program.TBINIT1;
            double relHumidity = Program.QUINIT;
            double inversionHeight = 10000;
            double userdefinedAirTemp = 0;
            
            if (File.Exists("CustomInit.txt"))
            {
                try
                {
                    using (FileStream fs = new FileStream("CustomInit.txt", FileMode.Open, FileAccess.Read, FileShare.Read))
                    {
                        string[] text;
                        using (StreamReader r = new StreamReader(fs))
                        {
                            r.ReadLine();
                            r.ReadLine();
                            for (int inid = 1; inid < Program.IWETTER; inid++)
                            {
                                r.ReadLine();
                            }
                            text = Convert.ToString(r.ReadLine()).Split(new char[] { ' ', ',', ';' }, StringSplitOptions.RemoveEmptyEntries);
                            if (text.Length > 4)
                            {
                                double.TryParse(text[4], NumberStyles.AllowDecimalPoint, CultureInfo.InvariantCulture, out snowHeight);
                            }
                            if (text.Length > 5)
                            {
                                double.TryParse(text[5], NumberStyles.AllowDecimalPoint, CultureInfo.InvariantCulture, out double temp);
                                if (temp > 0)
                                {
                                    tAir2m = temp;
                                    userdefinedAirTemp = 1;
                                }
                            }
                            if (text.Length > 6)
                            {
                                double.TryParse(text[6], NumberStyles.AllowDecimalPoint, CultureInfo.InvariantCulture, out double temp);
                                if (temp > 0)
                                {
                                    tSurface = temp;
                                }
                            }
                            if (text.Length > 7)
                            {
                                double.TryParse(text[7], NumberStyles.AllowDecimalPoint, CultureInfo.InvariantCulture, out double temp);
                                if (temp > 0)
                                {
                                    tSurface1m = temp;
                                }
                            }
                            if (text.Length > 8)
                            {
                                double.TryParse(text[8], NumberStyles.AllowDecimalPoint, CultureInfo.InvariantCulture, out double temp);
                                if (temp > 0)
                                {
                                    relHumidity = temp;
                                }
                            }
                            if (text.Length > 9)
                            {
                                double.TryParse(text[9], NumberStyles.AllowDecimalPoint, CultureInfo.InvariantCulture, out double temp);
                                if (temp > 0)
                                {
                                    inversionHeight = temp;
                                }
                            }
                        }
                    }
                }
                catch{}
            }
            //the array snow is not used, except for receptors (Ku 01/2022)
            for (int x = 0; x < Program.NX1; x++)
            {
                for (int y = 0; y < Program.NY1; y++)
                {
                    if (Program.AH[x][y] < snowHeight)
                    {
                        // no snow
                        Program.SNOW[x][y] = 0;
                    }
                    else
                    {
                        //set snow to 1
                        Program.SNOW[x][y] = 1;
                    }
                }
            }
            return (snowHeight, tAir2m, userdefinedAirTemp, tSurface, tSurface1m, relHumidity, inversionHeight);
        }
    }
}
