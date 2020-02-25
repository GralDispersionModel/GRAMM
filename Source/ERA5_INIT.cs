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
using Grib.Api;
using OSGeo.OSR;
using System.IO;
using System.Collections.Concurrent;

namespace GRAMM_CSharp_Test
{
    partial class Program
    {
        public static void ERA5_Read(DateTime dateref, ref DateTime date1, ref DateTime date2, ref List<float> X, ref List<float> Y, ref List<float> Z, ref List<float> ORO, ref List<float> U, ref List<float> V,
            ref List<float> W, ref List<float> T, ref List<float> Q, ref List<float> SD, ref List<float> CC, ref List<float> P, ref List<float> ST, ref List<float> SW, ref List<float> MSLP
            , ref List<float> LCC, ref List<float> MCC, ref List<float> HCC, ref List<float> CLOUDWATER, ref List<float> SEATEMP)
        {
            Console.WriteLine(" READING ERA5 GRIB DATA");
            //dateref = date for which ERA5 data should be read
            //date1 and date2 = dates before and after dateref used interpolation

            /*
            Parameter ID    Parameter                       Unit
            3               potential temperature           K
            34              sea surface temperature         K
            39              soil wetness 0-7cm              m³/m³
            42              soil wetness 72-189cm           m³/m³
            53              Montgomery potential            m²/s²   ->   g*z+cp*T
            54              pressure                        Pa
            129             geopotential                    m²/s²
            130             temperature                     K
            131             u component of wind             m/s
            132             v component of wind             m/s
            133             specific humidity               kg/kg
            135             vertical velocity               m/s
            139             soil temperature 0-7cm          K
            141             snow depth                      m of water equivalent
            151             mean sea level pressure         Pa
            152             logarithm of surface pressure   Pa
            164             total cloud cover               (0-1)
            186             lower cloud cover (<2km)        (0-1)
            187             medium cloud cover (2-6km)      (0-1)
            188             high cloud cover (>6km)         (0-1)
            236             soil temperature 72-189cm       K
            246             cloud water content             kg/kg
            247             cloud ice content               kg/kg
            228002          orography                       m
            */

            //search for grib file
            DirectoryInfo d = new DirectoryInfo(@"..\Computation");
            FileInfo[] Files = d.GetFiles("*.grib");
            string str = Convert.ToString(Files[0]);

            //transforming GRAMM boundaries into ERA5 coordinate system
            string GRAMM_projection = "";
            try
            {
                using (FileStream fs = new FileStream("GRAMM-coordinatesystem.txt", FileMode.Open, FileAccess.Read, FileShare.Read))
                {
                    using (StreamReader myreader = new StreamReader(fs))
                    {
                        GRAMM_projection = myreader.ReadLine();
                    }
                }
            }
            catch { }

            
            SpatialReference ERA5pro = new SpatialReference("");
            ERA5pro.SetWellKnownGeogCS("WGS84");
            SpatialReference GRAMMpro = new SpatialReference("");
            StringBuilder EPSG = new StringBuilder();
            EPSG.Append(GRAMM_projection);
            int GRAMMEPGS = Convert.ToInt16(GRAMM_projection);
            EPSG.Insert(0, "EPSG:");
            GRAMM_projection = EPSG.ToString();
            GRAMMpro.ImportFromEPSG(GRAMMEPGS);
            GRAMMpro.SetWellKnownGeogCS(GRAMM_projection);
            CoordinateTransformation ct = new CoordinateTransformation(ERA5pro, GRAMMpro);
            CoordinateTransformation ctrev = new CoordinateTransformation(GRAMMpro, ERA5pro);

            double[] p = new double[3];
            //get longitude of GRAMM domain
            double longitude = 0.5 * (Program.GRAMM_West + Program.GRAMM_West + Program.DDX[1] * Program.NX);
            p[1] = (float)Program.GRAMM_South; p[0] = (float)longitude; p[2] = 0;
            ctrev.TransformPoint(p);
            Program.Longitude = (float)p[0];

            float timeint = 0;
            float timepoint = 0;

            //determining the dates before and after the reference date
            string fname = Files[0].Name;
            bool found = false;
            TimeSpan span = new TimeSpan();
            for (int i = 0; i < Files.Length; i++)
            {
                fname = Files[i].Name;
                if (found == false)
                {
                    using (GribFile file = new GribFile(fname))
                    {
                        foreach (GribMessage msg in file)
                        {
                            DateTime date = msg.ReferenceTime;
                            span = date.Subtract(dateref);
                            timeint = Math.Abs((float)span.TotalMinutes);
                            if (DateTime.Compare(dateref, date) < 0 && timeint < 400)
                            {
                                date2 = date;
                                found = true;
                                break;
                            }
                            date1 = date;
                        }
                        file.Rewind();
                        TimeSpan span = date2.Subtract(date1);
                        timeint = (float)span.TotalSeconds;
                        span = dateref.Subtract(date1);
                        timepoint = (float)span.TotalSeconds;
                    }
                }
            }

            //reading all grib files
            int[] count = new int[20];
            for (int i = 0; i < Files.Length; i++)
            {
                fname = Files[i].Name;
                using (GribFile file = new GribFile(fname))
                {                    
                    double dummy = 0.0;

                    foreach (GribMessage msg in file)
                    {
                        int level = msg.Level;
                        DateTime date = msg.ReferenceTime;

                        if (DateTime.Compare(date, date1) == 0 || DateTime.Compare(date, date2) == 0)
                        {
                            foreach (GeoSpatialValue var in msg.GeoSpatialValues)
                            {
                                p[0] = (float)var.Longitude; p[1] = (float)var.Latitude; p[2] = 0;
                                try
                                {
                                    ct.TransformPoint(p);
                                    if (p[0] >= Program.GRAMM_West - 30000 && p[0] <= Program.GRAMM_West + 30000 + Program.DDX[1] * Program.NX &&
                                        p[1] >= Program.GRAMM_South - 30000 && p[1] <= Program.GRAMM_South + 30000 + Program.DDY[1] * Program.NY)
                                    {
                                        if (DateTime.Compare(date, date1) == 0)
                                        {
                                            //total cloud cover
                                            if (msg.ShortName.Equals("tcc"))
                                            {
                                                CC.Add((float)var.Value);
                                            }
                                            //low cloud cover
                                            if (msg.ShortName.Equals("lcc"))
                                            {
                                                LCC.Add((float)var.Value);
                                            }
                                            //medium cloud cover
                                            if (msg.ShortName.Equals("mcc"))
                                            {
                                                MCC.Add((float)var.Value);
                                            }
                                            //high cloud cover
                                            if (msg.ShortName.Equals("hcc"))
                                            {
                                                HCC.Add((float)var.Value);
                                            }
                                            //geopotential
                                            if (msg.ShortName.Equals("z"))
                                            {
                                                if (level == 0)
                                                {
                                                    //orography
                                                    ORO.Add((float)(var.Value / 9.8));
                                                }
                                                else
                                                {
                                                    Z.Add((float)(var.Value / 9.8));
                                                }
                                            }
                                            //U velocity and pressure level -> note that the U velocity as others must be provided at pressure levels
                                            if (msg.ShortName.Equals("u"))
                                            {
                                                X.Add((float)(p[0] - Program.GRAMM_West));
                                                Y.Add((float)(p[1] - Program.GRAMM_South));
                                                U.Add((float)var.Value);
                                                P.Add(level);
                                            }
                                            //V velocity
                                            if (msg.ShortName.Equals("v"))
                                            {
                                                V.Add((float)var.Value);
                                            }
                                            //W velocity
                                            if (msg.ShortName.Equals("w"))
                                            {
                                                W.Add((float)var.Value);
                                            }
                                            //Temperature
                                            if (msg.ShortName.Equals("t"))
                                            {
                                                T.Add((float)var.Value);
                                            }
                                            //specific humidity
                                            if (msg.ShortName.Equals("q"))
                                            {
                                                Q.Add((float)var.Value);
                                            }
                                            //specific cloud liquid water content
                                            if (msg.ShortName.Equals("clwc"))
                                            {
                                                CLOUDWATER.Add((float)var.Value);
                                            }
                                            //snow depth
                                            if (msg.ShortName.Equals("sd"))
                                            {
                                                SD.Add((float)var.Value);
                                            }
                                            //soil temperature
                                            if (msg.ShortName.Equals("stl1"))
                                            {
                                                ST.Add((float)var.Value);
                                            }
                                            //soil wetness
                                            if (msg.ShortName.Equals("swvl1"))
                                            {
                                                SW.Add((float)var.Value);
                                            }
                                            //mean sea level pressure
                                            if (msg.ShortName.Equals("msl"))
                                            {
                                                MSLP.Add((float)var.Value);
                                            }
                                            //mean sea surface temperature
                                            if (msg.ShortName.Equals("sst"))
                                            {
                                                SEATEMP.Add((float)var.Value);
                                            }
                                        }
                                        else
                                        {
                                            //carring out temporal interpolation
                                            //total cloud cover
                                            if (msg.ShortName.Equals("tcc"))
                                            {
                                                Linear_interpolation(0.0, timeint, timepoint, CC[count[0]], var.Value, ref dummy);
                                                CC[count[0]] = (float)dummy;
                                                count[0]++;
                                            }
                                            //low cloud cover
                                            if (msg.ShortName.Equals("lcc"))
                                            {
                                                Linear_interpolation(0.0, timeint, timepoint, LCC[count[11]], var.Value, ref dummy);
                                                LCC[count[11]] = (float)dummy;
                                                count[11]++;
                                            }
                                            //medium cloud cover
                                            if (msg.ShortName.Equals("mcc"))
                                            {
                                                Linear_interpolation(0.0, timeint, timepoint, MCC[count[12]], var.Value, ref dummy);
                                                MCC[count[12]] = (float)dummy;
                                                count[12]++;
                                            }
                                            //high cloud cover
                                            if (msg.ShortName.Equals("hcc"))
                                            {
                                                Linear_interpolation(0.0, timeint, timepoint, HCC[count[13]], var.Value, ref dummy);
                                                HCC[count[13]] = (float)dummy;
                                                count[13]++;
                                            }
                                            //geopotential -> note that the z coordinate varies as all quantities are delivered on pressure levels
                                            if (msg.ShortName.Equals("z"))
                                            {
                                                if (level > 0)
                                                {
                                                    Linear_interpolation(0.0, timeint, timepoint, Z[count[1]], var.Value / 9.8, ref dummy);
                                                    Z[count[1]] = (float)dummy;
                                                    count[1]++;
                                                }
                                            }
                                            //U velocity
                                            if (msg.ShortName.Equals("u"))
                                            {
                                                Linear_interpolation(0.0, timeint, timepoint, U[count[2]], var.Value, ref dummy);
                                                U[count[2]] = (float)dummy;
                                                count[2]++;
                                            }
                                            //V velocity
                                            if (msg.ShortName.Equals("v"))
                                            {
                                                Linear_interpolation(0.0, timeint, timepoint, V[count[3]], var.Value, ref dummy);
                                                V[count[3]] = (float)dummy;
                                                count[3]++;
                                            }
                                            // W velocity
                                            if (msg.ShortName.Equals("w"))
                                            {
                                                Linear_interpolation(0.0, timeint, timepoint, W[count[9]], var.Value, ref dummy);
                                                W[count[9]] = (float)dummy;
                                                count[9]++;
                                            }
                                            //Temperature
                                            if (msg.ShortName.Equals("t"))
                                            {
                                                Linear_interpolation(0.0, timeint, timepoint, T[count[4]], var.Value, ref dummy);
                                                T[count[4]] = (float)dummy;
                                                count[4]++;
                                            }
                                            //specific humidity
                                            if (msg.ShortName.Equals("q"))
                                            {
                                                Linear_interpolation(0.0, timeint, timepoint, Q[count[5]], var.Value, ref dummy);
                                                Q[count[5]] = (float)dummy;
                                                count[5]++;
                                            }
                                            //specific cloud liquid water content
                                            if (msg.ShortName.Equals("clwc"))
                                            {
                                                Linear_interpolation(0.0, timeint, timepoint, CLOUDWATER[count[14]], var.Value, ref dummy);
                                                CLOUDWATER[count[14]] = (float)dummy;
                                                count[14]++;
                                            }
                                            //snow depth
                                            if (msg.ShortName.Equals("sd"))
                                            {
                                                Linear_interpolation(0.0, timeint, timepoint, SD[count[6]], var.Value, ref dummy);
                                                SD[count[6]] = (float)dummy;
                                                count[6]++;
                                            }
                                            //soil temperature
                                            if (msg.ShortName.Equals("stl1"))
                                            {
                                                Linear_interpolation(0.0, timeint, timepoint, ST[count[7]], var.Value, ref dummy);
                                                ST[count[7]] = (float)dummy;
                                                count[7]++;
                                            }
                                            //soil wetness
                                            if (msg.ShortName.Equals("swvl1"))
                                            {
                                                Linear_interpolation(0.0, timeint, timepoint, SW[count[8]], var.Value, ref dummy);
                                                SW[count[8]] = (float)dummy;
                                                count[8]++;
                                            }
                                            //soil wetness
                                            if (msg.ShortName.Equals("msl"))
                                            {
                                                Linear_interpolation(0.0, timeint, timepoint, MSLP[count[10]], var.Value, ref dummy);
                                                MSLP[count[10]] = (float)dummy;
                                                count[10]++;
                                            }
                                            //soil wetness
                                            if (msg.ShortName.Equals("sst"))
                                            {
                                                Linear_interpolation(0.0, timeint, timepoint, SEATEMP[count[14]], var.Value, ref dummy);
                                                SEATEMP[count[14]] = (float)dummy;
                                                count[14]++;
                                            }
                                        }
                                    }
                                }
                                catch { }
                            }
                        }
                    }
                }
            }
        }

        //interpolate any 3D meteorological data onto the GRAMM grid
        public static void ERA5_InitializeOntoGRAMMgrid(List<float> X, List<float> Y, List<float> Z, List<float> U, List<float> V, List<float> W, List<float> T, List<float> Q, List<float> P, List<float> SD, List<float> CC, List<float> ST, List<float> SW, List<float> MSLP, List<float> LCC, List<float> MCC, List<float> HCC, List<float> CLOUDWATER, List<float> SEATEMP)
        {
            Console.WriteLine(" INTERPOLATING ERA5 DATA");

            double TMAX = 0;
            double TMIN = 374;
            double[][] MSLPinterim = CreateArray<double[]>(NX1, () => new double[NY1]);
            object obj = new object();

            //resetting variables
            Parallel.For(1, NX + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= NY; j++)
                {
                    double[] U_L = Program.U[i][j];
                    double[] V_L = Program.V[i][j];
                    double[] TABS_L = Program.TABS[i][j];
                    double[] QU_L = Program.QU[i][j];
                    double[] WATVAP_L = Program.WAT_VAP[i][j];
                    float[] PBZ_L = Program.PBZ[i][j];
                    Program.FW[i][j] = 0;
                    Program.CLOUDS[i][j] = 0;
                    Program.SNOW[i][j] = 0;
                    Program.TB[i][j][2] = 0;
                    for (int k = 1; k <= NZ; k++)
                    {
                        TABS_L[k] = 0;
                        U_L[k] = 0;
                        V_L[k] = 0;
                        QU_L[k] = 0;
                        WATVAP_L[k] = 0;
                        PBZ_L[k] = 0;
                    }
                }
            });
            Parallel.For(1, Program.NX + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= Program.NY; j++)
                {
                    double[] U_L = Program.U[i][j];
                    double[] V_L = Program.V[i][j];
                    double[] T_L = Program.T[i][j];
                    double[] TABS_L = Program.TABS[i][j];
                    double[] TBZ_L = Program.TBZ[i][j];
                    double[] QU_L = Program.QU[i][j];
                    double[] WATVAP_L = Program.WAT_VAP[i][j];
                    float[] PBZ_L = Program.PBZ[i][j];
                    double[] QBZ_L = Program.QBZ[i][j];
                    double[] TB_L = Program.TB[i][j];
                    float[] ZSP_L = Program.ZSP[i][j];
                    float[] FACTOR_L = Program.FACTOR[i][j];
                    float[] FAC_L = Program.FAC[i][j];
                    float AH_L = Program.AH[i][j];
                    for (int k = 1; k <= Program.NZ; k++)
                    {
                        //3D variables
                        double DUMMY1 = 0;
                        double DUMMY2 = 0;
                        double GEW = 0;
                        double SUMGEW = 0.0;
                        double SUMGEWSP = 0.0;
                        for (int n = 0; n < T.Count; n++)  //loop over all ERA5 cells within GRAMM domain
                        {
                            DUMMY1 = Math.Sqrt(Math.Pow(X[n] - Program.X[i] - Program.DDX[i] * 0.5, 2) + Math.Pow(Y[n] - Program.Y[j] - Program.DDY[j] * 0.5, 2));
                            if (DUMMY1 < 30000.0)
                            {
                                DUMMY1 = Math.Max(Math.Pow(DUMMY1, 2), 1.0);
                                DUMMY2 = Math.Max(Math.Pow(Z[n] + AH_L - ZSP_L[k], 4) * 0.001, 1.0);
                                GEW = DUMMY1 * DUMMY2 / 1000000.0;
                                SUMGEW += 1 / GEW;
                                TABS_L[k] += T[n] / GEW;
                                U_L[k] += U[n] / GEW;
                                V_L[k] += V[n] / GEW;
                                QU_L[k] += Q[n] * 1000.0 / GEW;
                                WATVAP_L[k] += CLOUDWATER[n] * 1000.0 / GEW;

                                if (ZSP_L[k] >= 0.0)
                                {
                                    SUMGEWSP += 1 / GEW;
                                    PBZ_L[k] += (float)(P[n] * 100.0 / GEW);
                                }
                            }
                        }
                        for (int n = 0; n < MSLP.Count; n++)  //additional loop for pressure to include surface pressure
                        {
                            DUMMY1 = Math.Sqrt(Math.Pow(X[n] - Program.X[i] - Program.DDX[i] * 0.5, 2) + Math.Pow(Y[n] - Program.Y[j] - Program.DDY[j] * 0.5, 2));
                            if (DUMMY1 < 30000.0)
                            {
                                DUMMY1 = Math.Max(Math.Pow(DUMMY1, 2), 1.0);
                                DUMMY2 = Math.Max(Math.Pow(0.0 - ZSP_L[k], 4) * 0.001, 1.0);
                                GEW = DUMMY1 * DUMMY2 / 1000000.0;
                                SUMGEWSP += 1 / GEW;

                                //Pressure calculation, if terrain is below sea level
                                if (ZSP_L[k] < 0.0)
                                {
                                    DUMMY2 = 9.81 * 0.029 / (8.31 * T[n]);
                                    DUMMY1 = MSLP[n] * Math.Exp(-ZSP_L[k] * DUMMY2);
                                    PBZ_L[k] += (float)(DUMMY1 / GEW);
                                }
                                else
                                {
                                    PBZ_L[k] += (float)(MSLP[n] / GEW);
                                }
                            }
                        }
                        //Weighting factors
                        if (SUMGEW != 0)
                        {
                            TABS_L[k] /= SUMGEW;
                            U_L[k] /= SUMGEW;
                            V_L[k] /= SUMGEW;
                            QU_L[k] /= SUMGEW;
                            QBZ_L[k] = QU_L[k];
                            WATVAP_L[k] /= SUMGEW;

                            //ERA 5 data only for initialization of temperature, but no large scale forcing via the geostrophic wind
                            if (Program.ISTAT >= 3)
                            {
                                U_L[k] = 0;
                                V_L[k] = 0;
                            }

                            //PBZ_L[k] /= (float)SUMGEWSP;
                            //FACTOR_L[k] = (float)(Math.Pow(MSLP[0] / PBZ_L[k], 0.287));
                            //FAC_L[k] = FACTOR_L[k];
                            //potential temperature
                            //T_L[k] = TABS_L[k] * FACTOR_L[k];
                            //hydrostatic balance
                            //TBZ_L[k] = T_L[k];

                            lock (obj)
                            {
                                if (TABS_L[k] < TMIN) TMIN = TABS_L[k];
                                if (TABS_L[k] > TMAX) TMAX = TABS_L[k];
                            }
                        }

                        //2D variables
                        if (k == 1)
                        {
                            SUMGEW = 0.0;
                            DUMMY1 = 0;
                            for (int n = 0; n < CC.Count; n++)  //loop over all ERA5 cells within GRAMM domain
                            {
                                DUMMY1 = Math.Sqrt(Math.Pow(X[n] - Program.X[i] - Program.DDX[i] * 0.5, 2) + Math.Pow(Y[n] - Program.Y[j] - Program.DDY[j] * 0.5, 2));
                                if (DUMMY1 < 30000.0)
                                {
                                    DUMMY1 = Math.Max(Math.Pow(DUMMY1, 2), 1.0);
                                    GEW = DUMMY1 / 1000000;
                                    SUMGEW += 1 / GEW;
                                    TB_L[2] += ST[n] / GEW;
                                    Program.FW[i][j] += SW[n] / GEW;
                                    //take height of clouds and topography into account
                                    if (AH_L > 6000.0)
                                    {
                                        Program.CLOUDS[i][j] += HCC[n] / GEW;
                                    }
                                    else if (AH_L <= 6000.0 && AH_L > 2000.0)
                                    {
                                        Program.CLOUDS[i][j] += Math.Max(HCC[n], MCC[n]) / GEW;
                                    }
                                    else
                                    {
                                        Program.CLOUDS[i][j] += Math.Max(HCC[n], Math.Max(MCC[n], LCC[n])) / GEW;
                                        ////Program.CLOUDS[i][j] += CC[n] / GEW;
                                    }

                                    Program.SNOW[i][j] += SD[n] * 10 / GEW;
                                    MSLPinterim[i][j] += MSLP[n] / GEW;
                                }
                            }
                            //Weighting factors
                            if (SUMGEW != 0)
                            {
                                TB_L[2] /= SUMGEW;
                                Program.FW[i][j] /= SUMGEW;
                                Program.CLOUDS[i][j] /= SUMGEW;
                                Program.SNOW[i][j] /= SUMGEW;
                                MSLPinterim[i][j] /= SUMGEW;

                                //no clouds considered in this case
                                if (Program.ISTAT == 3)
                                {
                                    Program.CLOUDS[i][j] = 0;
                                }
                            }
                        }
                    }
                }
            });

            /*
            //computation of how much water vapour exists as cloud liquid water
            Parallel.For(1, Program.NX + 1, Program.pOptions, i =>
            {
                for (int j = 2; j <= Program.NY - 1; j++)
                {
                    double[] WAT_VAPN_L = Program.WAT_VAPN[i][j];
                    double[] WAT_VAP_L = Program.WAT_VAP[i][j];
                    double[] TABS_L = Program.TABS[i][j];
                    double[] QUN_L = Program.QUN[i][j];
                    double[] QU_L = Program.QU[i][j];
                    for (int k = 1; k < Program.NZ; k++)
                    {
                        //saturation vapour pressure - Magnus formulae
                        double EGSAT = 611.2 * Math.Exp(17.269 * (TABS_L[k] - 273.16) / (TABS_L[k] - 30.04));
                        //maximum water vapour content possible -> assuming condensation above 80% relative humidity
                        double QGSAT = EGSAT / 461.5 / TABS_L[k];
                        //change in temperature
                        if (QUN_L[k] > QGSAT * 1000.0)
                        {
                            //condensation
                            double WATVAP_COND = (QUN_L[k] - QGSAT * 1000.0);
                            WAT_VAPN_L[k] += WATVAP_COND;
                            QUN_L[k] -= WATVAP_COND;
                            QU_L[k] = QUN_L[k];
                            WAT_VAP_L[k] = WAT_VAPN_L[k];
                        }
                        else
                        {
                            //evaporation
                            double WATVAP_EVA = Math.Min(WAT_VAP_L[k], QGSAT * 1000.0 - QUN_L[k]);
                            WAT_VAPN_L[k] -= WATVAP_EVA;
                            WAT_VAPN_L[k] = Math.Max(WAT_VAPN_L[k], 0.0);
                            QUN_L[k] += WATVAP_EVA;
                            QU_L[k] = QUN_L[k];
                            WAT_VAP_L[k] = WAT_VAPN_L[k];
                        }
                    }
                }
            });
            */

            //calculating the horizontally homegenuous basic state of p, T and rho according to COSMO model 
            //seeking for a close approximation of an idealized atmosphere to the actual observed one
            //for simplification the profile above the lowest topography is used
            double betamin = 42.0;
            double Tsurfacemin = 288.0;
            double fehlermin = 10000000.0;
            for (int beta = 42; beta <= 42; beta++)
            {
                for (int Tsurface = 210; Tsurface <= 350; Tsurface++)
                {
                    double fehler = 0.0;
                    for (int k = 1; k <= Program.NZ; k++)
                    {
                        double term = Math.Sqrt(1.0 - 2.0 * (float)beta * 9.81 * Program.ZSP[AHMINI][AHMINJ][k] / Program.GASCON / Math.Pow((float)Tsurface, 2));
                        fehler += Math.Pow(Program.TABS[AHMINI][AHMINJ][k] - (float)Tsurface * term, 2);
                    }
                    if (fehler < fehlermin)
                    {
                        betamin = beta;
                        Tsurfacemin = Tsurface;
                        fehlermin = fehler;
                    }
                }
            }
            Parallel.For(1, Program.NX + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= Program.NY; j++)
                {
                    double[] T_L = Program.T[i][j];
                    double[] TABS_L = Program.TABS[i][j];
                    double[] TBZ_L = Program.TBZ[i][j];
                    float[] ZSP_L = Program.ZSP[i][j];
                    float[] PBZ_L = Program.PBZ[i][j];
                    float[] RHOBZ_L = Program.RHOBZ[i][j];
                    float[] RHO_L = Program.RHO[i][j];
                    float[] FACTOR_L = Program.FACTOR[i][j];
                    float[] FAC_L = Program.FAC[i][j];

                    for (int k = 1; k <= Program.NZ; k++)
                    {
                        double term = Math.Sqrt(1.0 - 2.0 * betamin * 9.81 * ZSP_L[k] / Program.GASCON / Math.Pow(Tsurfacemin, 2));
                        TBZ_L[k] = Tsurfacemin * term;

                        TABS_L[k] = TBZ_L[k];

                        PBZ_L[k] = (float)(100000.0 * Math.Exp(-Tsurfacemin / betamin * (1.0 - term)));
                        RHOBZ_L[k] = (float)(PBZ_L[k] / TBZ_L[k] / Program.GASCON);
                        RHO_L[k] = RHOBZ_L[k];
                        FACTOR_L[k] = (float)(Math.Pow(100000.0 / PBZ_L[k], 0.287 * 1.0));
                        FAC_L[k] = FACTOR_L[k];
                        TBZ_L[k] *= FACTOR_L[k];

                        T_L[k] = TABS_L[k] * FACTOR_L[k];

                        //densitiy and pressure profile for radiation model
                        Program.RHOBZZ[k] = Program.RHOBZ[2][2][k];
                        Program.PBZZ[k] = Program.PBZ[2][2][k];
                    }
                }
            });

            //calculating the basic and initial state of Q, for simplification the profile above the lowest topography is used
            Parallel.For(1, Program.NX + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= Program.NY; j++)
                {
                    double[] QU_L = Program.QU[i][j];
                    double[] QBZ_L = Program.QBZ[i][j];
                    float[] ZSP_L = Program.ZSP[i][j];
                    float[] ZSPmin_L = Program.ZSP[AHMINI][AHMINJ];
                    double[] QUmin_L = Program.QU[AHMINI][AHMINJ];

                    for (int k = 1; k <= Program.NZ; k++)
                    {
                        int n1 = -1;
                        double QUref = 0;
                        for (int n = 1; n <= Program.NZ; n++)
                        {
                            if (ZSPmin_L[n] >= ZSP_L[k])
                            {
                                n1 = n;
                                break;
                            }
                        }
                        if (n1 < 0)
                        {
                            n1 = Program.NZ;
                        }
                        n1 = Math.Max(n1, 2);
                        Linear_interpolation(ZSPmin_L[n1 - 1] - AHMIN, ZSPmin_L[n1] - AHMIN, ZSP_L[k] - AHMIN, QUmin_L[n1 - 1], QUmin_L[n1], ref QUref);
                        QU_L[k] = QUref;
                        QBZ_L[k] = QU_L[k];
                    }
                }
            });

            //The soil temperature in 1m depth should be around 5-10°C at a latitude around 47N and a sea level of around 300m
            //Around 2000m above sea level permafrost exists in the Alps (=0°C)
            //Further, permafrost is found at latitudes larger/smaller than 70N/70S -> dependency on the latitude is considered
            //

            Program.TBINIT1 = 293 - 0.005 * Math.Pow(Program.BGRAD, 2) + 0.006 * Math.Abs(Program.BGRAD);

            //boundary values for surface parameters
            Parallel.For(1, Program.NX + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= Program.NY; j++)
                {
                    //friction velocity
                    Program.UST[i][j] = 0;
                    Program.USTV[i][j] = 0;
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
                    Program.TB[i][j][1] = Program.TB[i][j][2];
                    for (int kb = 3; kb <= Program.NZB; kb++)
                    {
                        //The surface temperature is initialized in dependence on the height above sea level and the latitude (ÖTTL, 22 JULY 2015)

                        DUMMY -= Program.DWB[kb];
                        //Dependency on height above sea level -> note that the dependency on latitude is taken into account in the calculation of TBINIT1 (tempint.cs)
                        double DELTA_TB = Program.TABS[i][j][1] - Program.TB[i][j][1]; //this is the desired temperature difference between the surface and the air temperature at the ground
                        Program.TB[i][j][kb] = (Program.TBINIT1 - 0.005 * Program.AH[i][j]) + (Program.TABS[i][j][1] - DELTA_TB - (Program.TBINIT1 - 0.005 * Program.AH[i][j])) * DUMMY;
                        Program.TG[i][j] = Program.TB[i][j][2];
                        Program.TBN[i][j][kb] = Program.TB[i][j][kb] - 0.1F;
                        Program.TBA[i][j][kb] = Program.TB[i][j][kb] + 0.1F;
                    }
                }
            });

            //Computing mass conservative ERA5 wind field at GRAMM grid
            //compute mass-divergences         
            Parallel.ForEach(Partitioner.Create(2, NX + 1, (int)((NX + 1) / Program.pOptions.MaxDegreeOfParallelism)), range =>
            {
                int NK_P = NZ; int NJ_P = NY;
                //Console.WriteLine(" Partitoner: " + Convert.ToString(range.Item1) +"/" + Convert.ToString(range.Item2));                             	
                for (int i = range.Item1; i < range.Item2; i++)
                {
                    for (int j = 2; j <= NJ_P; j++)
                    {
                        float[] SUX_L = Program.SUX[i][j];
                        float[] SUY_L = Program.SUY[i][j];
                        float[] SUZ_L = Program.SUZ[i][j];
                        float[] SUXYZ_L = Program.SUXYZ[i][j];
                        float[] AREA_L = Program.AREA[i][j];
                        float[] AREAZ_L = Program.AREAZ[i][j];
                        float[] AREAY_L = Program.AREAY[i][j];
                        float[] AREAX_L = Program.AREAX[i][j];
                        float[] AREAZX_L = Program.AREAZX[i][j];
                        float[] AREAZY_L = Program.AREAZY[i][j];
                        float[] AREAXYZ_L = Program.AREAXYZ[i][j];

                        for (int k = 1; k <= NK_P; k++)
                        {
                            //mass-divergence in east-west direction
                            if ((j < NJ_P) && (k < NK_P)) SUX_L[k] = (float)(Program.U[i - 1][j][k] * Program.RHOBZ[i - 1][j][k] - Program.U[i][j][k] * Program.RHOBZ[i][j][k]);

                            //mass-divergence in south-north direction
                            if ((i < NX) && (k < NK_P)) SUY_L[k] = (float)(Program.V[i][j - 1][k] * Program.RHOBZ[i][j - 1][k] - Program.V[i][j][k] * Program.RHOBZ[i][j][k]);

                            //mass-divergence in the z-direction
                            if ((i < NX) && (j < NJ_P)) SUZ_L[k] = (float)((Program.U[i][j][k - 1] * Program.RHOBZ[i][j][k - 1] - Program.U[i][j][k] * Program.RHOBZ[i][j][k]) * AREAZX_L[k] + (Program.V[i][j][k - 1] * Program.RHOBZ[i][j][k - 1] - Program.V[i][j][k] * Program.RHOBZ[i][j][k]) * AREAZY_L[k]) / AREAZ_L[k];

                            //mass-divergence between the two half-cells
                            if ((i < NX) && (j < NJ_P) && (k < NK_P))
                                SUXYZ_L[k] = 0;

                            //round-off errors cause the pressure equation to produce meaningless gradients
                            //cutting off the last digits solves this problem largely (Oettl, Sept 2015)
                            if (Math.Abs(SUX_L[k]) < 0.00001) SUX_L[k] = 0;
                            if (Math.Abs(SUY_L[k]) < 0.00001) SUY_L[k] = 0;
                            if (Math.Abs(SUZ_L[k]) < 0.00001) SUZ_L[k] = 0;
                            if (Math.Abs(SUXYZ_L[k]) < 0.00001) SUXYZ_L[k] = 0;
                        }
                    }
                }
            });

            TERMIVterms(NX, NY, NZ);
            TERMIPSterms(NX, NY, NZ);
            TERMIPterms(NX, NY, NZ);

            Program.INUMS = 0;
            while (Program.INUMS < 20)
            {
                Program.INUMS++;
                Primp_calculate(NX, NY, NZ);
            }

            //compute pressure gradients to correct wind speeds
            Parallel.For(2, NX, Program.pOptions, i =>
            {
                int NK_P = NZ; int NJ_P = NY;
                for (int j = 2; j <= NJ_P - 1; j++)
                {
                    float[] AREA_L = Program.AREA[i][j];
                    float[] AREAX_L = Program.AREAX[i][j];
                    float[] AREAY_L = Program.AREAY[i][j];
                    float[] AREAZX_L = Program.AREAZX[i][j];
                    float[] AREAZY_L = Program.AREAZY[i][j];

                    double[] DP_L = Program.DP[i][j];
                    double[] DPZ_L = Program.DPZ[i][j];
                    float[] AP0_L = Program.AP0[i][j];
                    double[] U1N_L = Program.U1N[i][j];
                    double[] V1N_L = Program.V1N[i][j];
                    double[] W1N_L = Program.W1N[i][j];
                    double[] U2N_L = Program.U2N[i][j];
                    double[] V2N_L = Program.V2N[i][j];
                    double[] W2N_L = Program.W2N[i][j];
                    double[] DDP1DX_L = Program.DDP1DX[i][j];
                    double[] DDP2DX_L = Program.DDP2DX[i][j];
                    double[] DDP1DY_L = Program.DDP1DY[i][j];
                    double[] DDP2DY_L = Program.DDP2DY[i][j];
                    double[] DDP1DZ_L = Program.DDP1DZ[i][j];
                    double[] DDP2DZ_L = Program.DDP2DZ[i][j];
                    double[] DPX_L = Program.DPX[i][j];
                    double[] DPY_L = Program.DPY[i][j];
                    double DP_LL = 0, DPZ_LL = 0, DPZp_LL = 0, f1 = 0, f2 = 0;
                    int m = 2;

                    for (int kn = 1; kn <= 2 * (NK_P - 1); kn++)
                    {
                        int k = 1 + kn >> 1;  // (int) (kn * 0.5F + 0.5F)

                        if (kn % 2 != 0) // ODD numbers 1,3,5,..- new k value - set local constants new
                        {
                            DP_LL = DP_L[k];
                            DPZ_LL = DPZ_L[k];
                            DPZp_LL = DPZ_L[k + 1];
                            f1 = (AREAX_L[k] + AREAZX_L[k]) * DP_LL;
                            f2 = (AREAY_L[k] + AREAZY_L[k]) * DP_LL;

                            //Pressure gradients

                            DDP1DX_L[k] = (AREAX_L[k] * DPX_L[k] - f1 + AREAZX_L[k] * DPZ_LL);

                            DDP2DX_L[k] = (-Program.AREAX[i + 1][j][k] * Program.DPX[i + 1][j][k] + f1 - AREAZX_L[k + 1] * DPZp_LL);

                            DDP1DY_L[k] = (AREAY_L[k] * DPY_L[k] - f2 + AREAZY_L[k] * DPZ_LL);

                            DDP2DY_L[k] = (-Program.AREAY[i][j + 1][k] * Program.DPY[i][j + 1][k] + f2 - AREAZY_L[k + 1] * DPZp_LL);

                            DDP1DZ_L[k] = AREA_L[k] * (DPZ_LL - DP_LL);

                            DDP2DZ_L[k] = AREA_L[k + 1] * (DP_LL - DPZp_LL);
                        }

                        //Velocity corrections
                        if ((m - 1) == 1)
                        {
                            m--;
                            float temp = 1 / AP0_L[k];
                            U1N_L[k] = Program.U[i][j][k] + DDP1DX_L[k] * temp;
                            V1N_L[k] = Program.V[i][j][k] + DDP1DY_L[k] * temp;
                            W1N_L[k] += DDP1DZ_L[k] * temp;
                        }
                        else
                        {
                            m = 2;
                            float temp = 1 / AP0_L[k];
                            U2N_L[k] = Program.U[i][j][k] + DDP2DX_L[k] * temp;
                            V2N_L[k] = Program.V[i][j][k] + DDP2DY_L[k] * temp;
                            W2N_L[k] += DDP2DZ_L[k] * temp;

                            Program.U[i][j][k] = 0.5F * (U1N_L[k] + U2N_L[k]);
                            Program.V[i][j][k] = 0.5F * (V1N_L[k] + V2N_L[k]);
                            Program.W[i][j][k] = 0.5F * (W1N_L[k] + W2N_L[k]);
                        }
                    }
                }
            });

            //Geostrophic wind estimation
            //Part I: Geostrophic wind at sea level
            double xgrad = 0;
            double ygrad = 0;
            int numbx = 0;
            int numby = 0;
            Parallel.For(0, MSLP.Count, Program.pOptions, i =>
            {
                for (int j = 0; j < MSLP.Count; j++)
                {
                    if (i != j && Math.Abs(X[i] - X[j]) > 10000)
                    {
                        xgrad += (MSLP[i] - MSLP[j]) / (X[i] - X[j]);
                        numbx++;
                    }
                    if (i != j && Math.Abs(Y[i] - Y[j]) > 10000)
                    {
                        ygrad += (MSLP[i] - MSLP[j]) / (Y[i] - Y[j]);
                        numby++;
                    }
                }
            });
            xgrad /= numbx;
            ygrad /= numby;
            double UGsea = -ygrad / 1.2 / Program.FN;
            double VGsea = xgrad / 1.2 / Program.FN;

            //PartII: Geostrophic wind above mixing layer
            int numb = 0;
            double UGtemp = 0;
            double VGtemp = 0;
            Parallel.For(1, Program.NX + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= Program.NY; j++)
                {
                    double[] U_L = Program.U[i][j];
                    double[] V_L = Program.V[i][j];
                    float[] ZSP_L = Program.ZSP[i][j];
                    for (int k = Program.NZ; k >= 1; k--)
                    {
                        if (ZSP_L[k] - Program.AHMAX > 2000.0 && ZSP_L[k] < 10000.0)
                        {
                            UGtemp += U_L[k];
                            VGtemp += V_L[k];
                            numb++;
                        }
                    }
                }
            });
            double UGfree = UGtemp / numb;
            double VGfree = VGtemp / numb;

            Parallel.For(1, Program.NX + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= Program.NY; j++)
                {
                    double[] UG_L = Program.UG[i][j];
                    double[] VG_L = Program.VG[i][j];
                    float[] UG1_L = Program.UG1[i][j];
                    float[] VG1_L = Program.VG1[i][j];
                    float[] UG2_L = Program.UG2[i][j];
                    float[] VG2_L = Program.VG2[i][j];
                    double[] U_L = Program.U[i][j];
                    double[] V_L = Program.V[i][j];
                    double[] W_L = Program.W[i][j];
                    double[] UN_L = Program.UN[i][j];
                    double[] VN_L = Program.VN[i][j];
                    double[] U1_L = Program.U1[i][j];
                    double[] V1_L = Program.V1[i][j];
                    double[] U2_L = Program.U2[i][j];
                    double[] V2_L = Program.V2[i][j];
                    double[] U1N_L = Program.U1N[i][j];
                    double[] V1N_L = Program.V1N[i][j];
                    double[] U2N_L = Program.U2N[i][j];
                    double[] V2N_L = Program.V2N[i][j];
                    float[] ZSP_L = Program.ZSP[i][j];
                    for (int k = Program.NZ; k >= 1; k--)
                    {
                        /*
                        if (ZSP_L[k] > AHMAX + 2000)
                        {
                            UG_L[k] = UGfree;
                            VG_L[k] = VGfree;
                        }
                        else if(ZSP_L[k] <= AHMAX + 2000 && ZSP_L[k] >= AHMAX)
                        {
                            Linear_interpolation(0, AHMAX + 2000, ZSP_L[k], 0, UGfree, ref UG_L[k]);
                            Linear_interpolation(0, AHMAX + 2000, ZSP_L[k], 0, VGfree, ref VG_L[k]);
                        }
                        else
                        {
                            UG_L[k] = 0;
                            VG_L[k] = 0;
                        }

                        UG1_L[k] = (float)UG_L[k];
                        UG2_L[k] = (float)UG_L[k];
                        VG1_L[k] = (float)VG_L[k];
                        VG2_L[k] = (float)VG_L[k];
                        U_L[k] = UG_L[k];
                        V_L[k] = VG_L[k];
                        Program.U1NRHO[i][j][k] = (float)(U_L[k] * Program.RHO[i][j][k]);
                        Program.U2NRHO[i][j][k] = (float)(U_L[k] * Program.RHO[i][j][k]);
                        U1_L[k] = U_L[k];
                        U2_L[k] = U_L[k];
                        U1N_L[k] = U_L[k];
                        U2N_L[k] = U_L[k];
                        UN_L[k] = U_L[k];
                        Program.V1NRHO[i][j][k] = (float)(V_L[k] * Program.RHO[i][j][k]);
                        Program.V2NRHO[i][j][k] = (float)(V_L[k] * Program.RHO[i][j][k]);
                        V1_L[k] = V_L[k];
                        V2_L[k] = V_L[k];
                        V1N_L[k] = V_L[k];
                        V2N_L[k] = V_L[k];
                        VN_L[k] = V_L[k];
                        */

                        //ERA 5 data only for initialization of temperature, but no large scale forcing via the geostrophic wind
                        if (Program.ISTAT >= 2)
                        {
                            UG_L[k] = 0;
                            VG_L[k] = 0;
                        }

                        //TEST
                        /*
                        UG_L[k] = (float)UGfree;
                        VG_L[k] = (float)VGfree;
                        UG1_L[k] = (float)UGfree;
                        UG2_L[k] = (float)UGfree;
                        VG1_L[k] = (float)VGfree;
                        VG2_L[k] = (float)VGfree;
                        */
                        Program.UN[i][j][k] = U_L[k];
                        Program.U1N[i][j][k] = U_L[k];
                        Program.U1[i][j][k] = U_L[k];
                        Program.U2N[i][j][k] = U_L[k];
                        Program.U2[i][j][k] = U_L[k];
                        Program.VN[i][j][k] = V_L[k];
                        Program.V1N[i][j][k] = V_L[k];
                        Program.V1[i][j][k] = V_L[k];
                        Program.V2N[i][j][k] = V_L[k];
                        Program.V2[i][j][k] = V_L[k];
                        Program.WN[i][j][k] = W_L[k];
                        Program.W1N[i][j][k] = W_L[k];
                        Program.W1[i][j][k] = W_L[k];
                        Program.W2N[i][j][k] = W_L[k];
                        Program.W2[i][j][k] = W_L[k];
                    }
                }
            });

            //reading land-use data
            double DUMMY11 = 0;
            double DUMMY21 = 0;
            double DUMMY4 = 0;
            double DUMMY5 = 0;
            double DUMMY6 = 0;

            if (File.Exists("landuse.asc"))
            {
                string[] text = new string[(Program.NX + 2) * (Program.NY + 2)];
                try
                {
                    using (FileStream fs = new FileStream("landuse.asc", FileMode.Open, FileAccess.Read, FileShare.Read))
                    {
                        using (StreamReader r = new StreamReader(fs))
                        {
                            text = Convert.ToString(r.ReadLine()).Split(new char[] { ' ', ',', ';' }, StringSplitOptions.RemoveEmptyEntries);
                            int n = 0;
                            for (int j = 1; j < Program.NY + 1; j++)
                            {
                                for (int i = 1; i < Program.NX + 1; i++)
                                {
                                    Program.RHOB[i][j] = Convert.ToSingle(text[n].Replace(".", Program.decsep));
                                    n++;
                                }
                            }
                            text = Convert.ToString(r.ReadLine()).Split(new char[] { ' ', ',', ';' }, StringSplitOptions.RemoveEmptyEntries);
                            n = 0;
                            for (int j = 1; j < Program.NY + 1; j++)
                            {
                                for (int i = 1; i < Program.NX + 1; i++)
                                {
                                    Program.ALAMBDA[i][j] = Convert.ToSingle(text[n].Replace(".", Program.decsep));
                                    n++;
                                }
                            }
                            text = Convert.ToString(r.ReadLine()).Split(new char[] { ' ', ',', ';' }, StringSplitOptions.RemoveEmptyEntries);
                            n = 0;
                            for (int j = 1; j < Program.NY + 1; j++)
                            {
                                for (int i = 1; i < Program.NX + 1; i++)
                                {
                                    Program.Z0[i][j] = Convert.ToSingle(text[n].Replace(".", Program.decsep));
                                    n++;
                                }
                            }
                            //the default soil wetness based on CORINE is read and it is determined, whether this is a dry area or not
                            text = Convert.ToString(r.ReadLine()).Split(new char[] { ' ', ',', ';' }, StringSplitOptions.RemoveEmptyEntries);
                            n = 0;
                            for (int j = 1; j < Program.NY + 1; j++)
                            {
                                for (int i = 1; i < Program.NX + 1; i++)
                                {
                                    Program.DRY_AREA[i][j] = false;
                                    if (Convert.ToSingle(text[n].Replace(".", Program.decsep)) < 0.05)
                                    {
                                        Program.DRY_AREA[i][j] = true;
                                    }
                                    n++;
                                }
                            }
                            text = Convert.ToString(r.ReadLine()).Split(new char[] { ' ', ',', ';' }, StringSplitOptions.RemoveEmptyEntries);
                            n = 0;
                            for (int j = 1; j < Program.NY + 1; j++)
                            {
                                for (int i = 1; i < Program.NX + 1; i++)
                                {
                                    Program.EPSG[i][j] = Convert.ToDouble(text[n].Replace(".", Program.decsep));
                                    n++;
                                }
                            }
                            text = Convert.ToString(r.ReadLine()).Split(new char[] { ' ', ',', ';' }, StringSplitOptions.RemoveEmptyEntries);
                            n = 0;
                            for (int j = 1; j < Program.NY + 1; j++)
                            {
                                for (int i = 1; i < Program.NX + 1; i++)
                                {
                                    Program.ALBEDO[i][j] = Convert.ToDouble(text[n].Replace(".", Program.decsep));
                                    n++;
                                }
                            }
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
                DUMMY11 = Convert.ToDouble(Console.ReadLine().Replace(".", Program.decsep));

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
                DUMMY21 = Convert.ToDouble(Console.ReadLine().Replace(".", Program.decsep));

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
                for (int j = 1; j < Program.NY + 1; j++)
                {
                    for (int i = 1; i < Program.NX + 1; i++)
                    {
                        Program.ALBEDO[i][j] = DUMMY11;
                        Program.Z0[i][j] = (float)(DUMMY21);
                        Program.EPSG[i][j] = DUMMY6;
                        Program.ALAMBDA[i][j] = (float)(DUMMY5);
                        Program.RHOB[i][j] = (float)(Program.ALAMBDA[i][j] / DUMMY4 / Program.CPBOD);
                    }
                }
            }

            //set wetness for water bodies and dry areas (due to the coarse ERA5 data resolution, this is necessary to better account for these areas)
            Parallel.For(1, Program.NX + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= Program.NY; j++)
                {
                    if (Program.ALAMBDA[i][j] > 99.0)
                    {
                        for (int kb = 1; kb <= Program.NZB; kb++)
                        {
                            Program.FW[i][j] = 1.0;

                            if (Program.ALBEDO[i][j] == 0.081 && kb >= 2)
                            {
                                //sea surface temperature homogenous within 1 m
                                Program.TB[i][j][kb] = Program.TB[i][j][kb - 1];
                                Program.TG[i][j] = Program.TB[i][j][2];
                                Program.TBN[i][j][kb] = Program.TB[i][j][kb] - 0.1F;
                                Program.TBA[i][j][kb] = Program.TB[i][j][kb] + 0.1F;
                            }
                        }
                    }
                    if (Program.DRY_AREA[i][j] == true)
                    {
                        Program.FW[i][j] = 0.001;
                    }
                }
            });

            Int32 I = Program.AHMINI;
            Int32 J1 = Program.AHMINJ;
            Console.WriteLine();
            for (int k = Program.NZ; k >= 1; k--)
            {
                Console.WriteLine("  HEIGHT : " + Convert.ToString(Math.Round(Program.ZSP[I][J1][k], 1).ToString("0.0")).PadLeft(6) +
                    "   Tpot = " + Convert.ToString(Math.Round(Program.T[I][J1][k], 2).ToString("0.00")).PadLeft(6) +
                    "   T = " + Convert.ToString(Math.Round(Program.TABS[I][J1][k], 2).ToString("0.00")).PadLeft(6) +
                    "   P = " + Convert.ToString(Math.Round(Program.PBZ[I][J1][k], 2).ToString("0")).PadLeft(9) +
                    "   RHO = " + Convert.ToString(Math.Round(Program.RHOBZ[I][J1][k], 4).ToString("0.000")).PadLeft(6));
            }

            Console.WriteLine();
            Console.WriteLine("   MAXIMUM TEMPERATURE : " + Convert.ToString(Math.Round(TMAX, 2)).PadLeft(6));
            Console.WriteLine("   MINIMUM TEMPERATURE : " + Convert.ToString(Math.Round(TMIN, 2)).PadLeft(6));
            Console.WriteLine();
            if (TMIN < 120)
            {
                Console.WriteLine(" TEMPERATURE BELOW 120K !!");
                Environment.Exit(0);
            }

            //Profiles of the radiation model
            Int32 NPROF = 30;
            Program.ZPROF[1] = Math.Max(Program.AHMIN, 0);
            double UNTO = 0;
            double UNTU = 0;
            double DIFF = 0;
            for (int n = 1; n <= NPROF; n++)
            {
                Int32 INDO = 0;
                Int32 INDU = 0;
                if (n > 1) Program.ZPROF[n] = Program.ZPROF[n - 1] + 500;
                if (Program.ZPROF[n] > Program.ZSP[Program.AHMINI][Program.AHMINJ][Program.NZ])
                {
                    //pressure profile above GRAMM domain
                    Program.TPROF[n] = Program.TBZ[Program.AHMINI][Program.AHMINJ][Program.NZ];
                    //temperature profile above GRAMM domain
                    Program.PPROF[n] = MSLP[0] * Math.Exp(-Program.ZPROF[n] / 8000);
                    //cloud water content above GRAMM domain
                    Program.QCLD[n] = 0.0;
                    //water vapour content above GRAMM domain
                    Program.QVAP[n] = 0.0;
                }
                else if (Program.ZPROF[n] <= Program.ZSP[Program.AHMINI][Program.AHMINJ][Program.NZ])
                {
                    //temperature and pressure profile within GRAMM domain
                    UNTU = -100000;
                    UNTO = 100000;
                    INDO = 0;
                    INDU = 0;
                    double GRAD = 0;
                    double GRADP = 0;
                    for (int m = 1; m <= Program.NZ; m++)
                    {
                        DIFF = Program.ZSP[Program.AHMINI][Program.AHMINJ][m] - Program.ZPROF[n];
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
                        GRAD = (Program.TBZ[Program.AHMINI][Program.AHMINJ][INDU] - Program.TBZ[Program.AHMINI][Program.AHMINJ][INDU - 1]) /
                            (Program.ZSP[Program.AHMINI][Program.AHMINJ][INDU] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDU - 1]);
                        Program.TPROF[n] = Program.TBZ[Program.AHMINI][Program.AHMINJ][INDU] + (Program.ZPROF[n] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDU]) * GRAD;

                        GRADP = (Program.PBZ[Program.AHMINI][Program.AHMINJ][INDU] - Program.PBZ[Program.AHMINI][Program.AHMINJ][INDU - 1]) /
                            (Program.ZSP[Program.AHMINI][Program.AHMINJ][INDU] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDU - 1]);
                        Program.PPROF[n] = Program.PBZ[Program.AHMINI][Program.AHMINJ][INDU] + (Program.ZPROF[n] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDU]) * GRADP;


                        GRAD = (Program.QU[Program.AHMINI][Program.AHMINJ][INDU] - Program.QU[Program.AHMINI][Program.AHMINJ][INDU - 1]) /
                            (Program.ZSP[Program.AHMINI][Program.AHMINJ][INDU] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDU - 1]);
                        Program.QVAP[n] = (Program.QU[Program.AHMINI][Program.AHMINJ][INDU] + (Program.ZPROF[n] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDU]) * GRAD) * 0.001;
                    }
                    else if (INDU == 0)
                    {
                        GRAD = (Program.TBZ[Program.AHMINI][Program.AHMINJ][INDO] - Program.TBZ[Program.AHMINI][Program.AHMINJ][INDO + 1]) /
                            (Program.ZSP[Program.AHMINI][Program.AHMINJ][INDO] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDO + 1]);
                        Program.TPROF[n] = Program.TBZ[Program.AHMINI][Program.AHMINJ][INDO] + (Program.ZPROF[n] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDO]) * GRAD;

                        GRADP = (Program.PBZ[Program.AHMINI][Program.AHMINJ][INDO] - Program.PBZ[Program.AHMINI][Program.AHMINJ][INDO + 1]) /
                            (Program.ZSP[Program.AHMINI][Program.AHMINJ][INDO] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDO + 1]);
                        Program.PPROF[n] = Program.PBZ[Program.AHMINI][Program.AHMINJ][INDO] + (Program.ZPROF[n] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDO]) * GRADP;

                        GRAD = (Program.QU[Program.AHMINI][Program.AHMINJ][INDO] - Program.QU[Program.AHMINI][Program.AHMINJ][INDO + 1]) /
                            (Program.ZSP[Program.AHMINI][Program.AHMINJ][INDO] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDO + 1]);
                        Program.QVAP[n] = (Program.QU[Program.AHMINI][Program.AHMINJ][INDO] + (Program.ZPROF[n] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDO]) * GRAD) * 0.001;

                    }
                    else
                    {
                        GRAD = (Program.TBZ[Program.AHMINI][Program.AHMINJ][INDO] - Program.TBZ[Program.AHMINI][Program.AHMINJ][INDU]) /
                            (Program.ZSP[Program.AHMINI][Program.AHMINJ][INDO] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDU]);
                        Program.TPROF[n] = Program.TBZ[Program.AHMINI][Program.AHMINJ][INDO] + (Program.ZPROF[n] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDO]) * GRAD;

                        GRADP = (Program.PBZ[Program.AHMINI][Program.AHMINJ][INDO] - Program.PBZ[Program.AHMINI][Program.AHMINJ][INDU]) /
                            (Program.ZSP[Program.AHMINI][Program.AHMINJ][INDO] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDU]);
                        Program.PPROF[n] = Program.PBZ[Program.AHMINI][Program.AHMINJ][INDO] + (Program.ZPROF[n] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDO]) * GRADP;

                        GRAD = (Program.QU[Program.AHMINI][Program.AHMINJ][INDO] - Program.QU[Program.AHMINI][Program.AHMINJ][INDU]) /
                            (Program.ZSP[Program.AHMINI][Program.AHMINJ][INDO] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDU]);
                        Program.QVAP[n] = (Program.QU[Program.AHMINI][Program.AHMINJ][INDO] + (Program.ZPROF[n] - Program.ZSP[Program.AHMINI][Program.AHMINJ][INDO]) * GRAD) * 0.001;

                    }
                }
                //compute absolute temperature based on absolute temperature
                Program.TPROF[n] /= Math.Pow(MSLP[0] / Program.PPROF[n], 0.287);
                if (Program.TPROF[n] < 153.15) Program.TPROF[n] = 153.15;
                Program.VNORM[n] = 36.5419617 + 4.8939118 * (Program.ZPROF[n] * 0.001) +
                    4.1091542 * Math.Pow(Program.ZPROF[n] * 0.001, 2) - 0.1456879 * Math.Pow(Program.ZPROF[n] * 0.001, 3) + 0.0149291 * Math.Pow(Program.ZPROF[n] * 0.001, 4);
                Program.VNORM[n] *= 1000;

                //water content in clouds for the radiation model
                Program.QCLD[n] = 0.0;

                //water content in rain for the radiation model
                Program.QRAIN[n] = 0;

                //water content in ice for the radiation model
                Program.QICE[n] = 0;

                //water content in snow for the radiation model
                Program.QSNOW[n] = 0;
            }

            Parallel.For(1, Program.NX + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= Program.NY; j++)
                {
                    float[] RHO_L = Program.RHO[i][j];
                    float[] RHOBZ_L = Program.RHOBZ[i][j];
                    double[] W1_L = Program.W1[i][j];
                    double[] W2_L = Program.W2[i][j];
                    double[] W1N_L = Program.W1N[i][j];
                    double[] W2N_L = Program.W2N[i][j];
                    double[] PN_L = Program.PN[i][j];
                    double[] DP_L = Program.DP[i][j];
                    double[] TE_L = Program.TE[i][j];
                    double[] TEN_L = Program.TEN[i][j];
                    double[] DT_SOL_L = Program.DT_SOL[i][j];
                    double[] DT_TERR_L = Program.DT_TERR[i][j];

                    for (int k = 1; k <= Program.NZ; k++)
                    {
                        //densitiy
                        RHO_L[k] = (float)(Math.Round(RHOBZ_L[k], 3));

                        //vertical velocities
                        W1_L[k] = 0;
                        W2_L[k] = 0;
                        W1N_L[k] = 0;
                        W2N_L[k] = 0;

                        //non-hydrostatic pressure
                        PN_L[k] = 0;
                        DP_L[k] = 0;

                        //turbulent-kinetic energy
                        TE_L[k] = 0;
                        TEN_L[k] = 0;

                        //solar and terrestrial radiation heating rates
                        DT_SOL_L[k] = 0;
                        DT_TERR_L[k] = 0;
                    }
                }
            });

            //Alpha: parameter used for nudging variables towards large-scale values at lateral boundaries
            for (int k = 1; k <= Program.NZ; k++)
            {
                Program.ALPHA[k] = 0.05;
            }

            //van Karman constant
            Program.CK = 0.35;

            //minimum turbulent viscosity
            Program.VISEL = 0.05;

            //heat capacity soil
            Program.CPBOD = 900;

            //Stefan-Bolzmann constant
            Program.SIGMA = 5.6697e-8;

            double UINIT = Program.U[1][1][Program.NZ];
            double VINIT = Program.V[1][1][Program.NZ];

            //factor used to initialize turbulent kinetic energy
            double TURBIN = 0.01;
            double TEINIT = 0;

            //turbulent Prandtl-number
            Program.PRTE = 0.9;

            //initial turbulent kinetic energy
            if ((UINIT == 0) && (VINIT == 0))
                TEINIT = 0.01;
            else
                TEINIT = TURBIN * (Math.Pow(UINIT, 2) + Math.Pow(VINIT, 2));

            //turbulent eddy viscosity
            double VISINIT = 0;
            if (TEINIT > 0)
                VISINIT = Program.CK * Math.Sqrt(TEINIT) * 10;
            else if (TEINIT == 0)
                VISINIT = Program.VISEL;
            

            //computing specific humidity from relative humidity
            Parallel.For(1, Program.NX + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= Program.NY; j++)
                {
                    float[] PBZ_L = Program.PBZ[i][j];
                    float[] FAC_L = Program.FAC[i][j];
                    float[] FACT_L = Program.FACTOR[i][j];
                    double[] QBZ_L = Program.QBZ[i][j];

                    for (int k = 1; k <= Program.NZ; k++)
                    {
                        double TBZN2 = Math.Round(Program.TABS[i][j][k] - 153.15, 3);
                        //                        if (Convert.ToInt32(TBZN) > 209) TBZN = 209;
                        //                        Int32 TBZNINT = Convert.ToInt32(TBZN);
                        TBZN2 = Math.Max(0, Math.Min(209.999, TBZN2));
                        int TBZNINT2 = Convert.ToInt32(Math.Floor(TBZN2));
                        double PDST2 = Program.PSAT[TBZNINT2 + 1] + (Program.PSAT[TBZNINT2 + 2] - Program.PSAT[TBZNINT2 + 1]) * (TBZN2 - (float)TBZNINT2);

                        // initial specific humidity
                        //QBZ_L[k] = 18.02F / 28.96F * PDST2 / (PBZ_L[k] / Program.QUINIT - PDST2) * 1000;

                        FAC_L[k] = FACT_L[k];
                    }
                }
            });

            //temperature value used to improve numerical accuracy of the solution algorithm for temperature
            Program.TBZ1 = Program.TBZ[2][2][1];

            Parallel.For(1, Program.NX + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= Program.NY; j++)
                {
                    double[] T_L = Program.T[i][j];
                    double[] TE_L = Program.TE[i][j];
                    double[] TN_L = Program.TN[i][j];
                    double[] TEN_L = Program.TEN[i][j];
                    double[] DISS_L = Program.DISS[i][j];
                    double[] DISSN_L = Program.DISSN[i][j];
                    double[] VISH_L = Program.VISH[i][j];
                    double[] VISV_L = Program.VISV[i][j];

                    for (int k = 1; k <= Program.NZ; k++)
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
            Parallel.For(1, Program.NX + 1, Program.pOptions, i =>
            {
                for (int k = 1; k <= Program.NZ; k++)
                {
                    int j = 1;
                    //large-scale values at the southern border for the next time step
                    Program.TSSN[i][k] = Program.T[i][j][k];

                    Program.QUSSN[i][k] = Program.QBZ[i][j][k];

                    Program.USSN[i][k] = Program.U[i][j][k];

                    Program.VSSN[i][k] = Program.V[i][j][k];

                    Program.WSSN[i][k] = 0;

                    j = Program.NY;
                    Program.TSNN[i][k] = Program.T[i][j][k];

                    Program.QUSNN[i][k] = Program.QBZ[i][j][k];

                    Program.USNN[i][k] = Program.U[i][j][k];

                    Program.VSNN[i][k] = Program.V[i][j][k];

                    Program.WSNN[i][k] = 0;
                }
            });

            //Values at west and east borders
            Parallel.For(1, Program.NY + 1, Program.pOptions, j =>
            {
                for (int k = 1; k <= Program.NZ; k++)
                {
                    int i = 1;
                    Program.TSWN[j][k] = Program.T[i][j][k];

                    Program.QUSWN[j][k] = Program.QBZ[i][j][k];

                    Program.USWN[j][k] = Program.U[i][j][k];

                    Program.VSWN[j][k] = Program.V[i][j][k];

                    Program.WSWN[j][k] = 0;

                    i = Program.NX;
                    Program.TSEN[j][k] = Program.T[i][j][k];

                    Program.QUSEN[j][k] = Program.QBZ[i][j][k];

                    Program.USEN[j][k] = Program.U[i][j][k];

                    Program.VSEN[j][k] = Program.V[i][j][k];

                    Program.WSEN[j][k] = 0;
                }
            });

            //boundary values for surface parameters
            Parallel.For(1, Program.NX + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= Program.NY; j++)
                {
                    //friction velocity
                    Program.UST[i][j] = 0.15;
                    Program.USTV[i][j] = 0.15 * Program.CK / Math.Log((Program.ZSP[3][3][1] - Program.AH[3][3]) / Program.Rauigkeit);

                    //initial values for stability classes
                    Program.stabilityclass[i][j] = Program.AKLA;

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
                }
            });

            Int32 IZELL = 0;
            double ARSUM = 0;
            double VOLSUM = 0;
            for (int i = 2; i < Program.NX; i++)
            {
                for (int j = 2; j < Program.NY; j++)
                {
                    ARSUM += Program.DDX[i] * Program.DDY[j] * 0.000001;
                    for (int k = 1; k < Program.NZ; k++)
                    {
                        IZELL++;
                        VOLSUM += Program.VOL[i][j][k] * 0.000000001;
                    }
                }
            }
            Console.WriteLine(" IINIT : NUMBER OF CALCULATED CELLS    : " + IZELL.ToString());
            Console.WriteLine(" IINIT : SIZE OF CALCULATED AREA       : " + ARSUM.ToString("0") + " km**2");
            Console.WriteLine(" IINIT : SIZE OF CALCULATED VOLUME     : " + VOLSUM.ToString("0") + " km**3");
        }

        //calculate boundary conditions
        public static void ERA5_BoundaryConditions(List<float> XERA5, List<float> YERA5, List<float> ZERA5, List<float> UERA5, List<float> VERA5, List<float> WERA5, List<float> TERA5, List<float> QERA5, List<float> PERA5, List<float> SDERA5, List<float> CCERA5, List<float> STERA5, List<float> SWERA5, List<float> MSLPERA5, List<float> LCCERA5, List<float> MCCERA5, List<float> HCCERA5, List<float> CLOUDWATERERA5, List<float> SEATEMPERA5)
        {
            Console.WriteLine(" RE-INITIALIZATION USING ERA5 DATA");
            double[][][] Uinterim = CreateArray<double[][]>(Program.NX1, () => CreateArray<double[]>(Program.NY1, () => new double[Program.NZ1]));
            double[][][] Vinterim = CreateArray<double[][]>(Program.NX1, () => CreateArray<double[]>(Program.NY1, () => new double[Program.NZ1]));
            double[][][] Tinterim = CreateArray<double[][]>(Program.NX1, () => CreateArray<double[]>(Program.NY1, () => new double[Program.NZ1]));
            double[][][] QUinterim = CreateArray<double[][]>(Program.NX1, () => CreateArray<double[]>(Program.NY1, () => new double[Program.NZ1]));
            double[][][] CWinterim = CreateArray<double[][]>(Program.NX1, () => CreateArray<double[]>(Program.NY1, () => new double[Program.NZ1]));
            double[][] MSLPinterim = CreateArray<double[]>(NX1, () => new double[NY1]);
            object obj = new object();
            Parallel.For(1, Program.NX + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= Program.NY; j++)
                {
                    double[] U_L = Uinterim[i][j];
                    double[] V_L = Vinterim[i][j];
                    double[] T_L = Tinterim[i][j];
                    double[] QU_L = QUinterim[i][j];
                    double[] CW_L = CWinterim[i][j];
                    float[] ZSP_L = Program.ZSP[i][j];
                    float[] FACTOR_L = Program.FACTOR[i][j];
                    float AH_L = Program.AH[i][j];
                    //set all to zero
                    Program.FW[i][j] = 0.0;
                    Program.CLOUDS[i][j] = 0.0;
                    Program.SNOW[i][j] = 0.0;
                    for (int k = 1; k <= Program.NZ; k++)
                    {
                        //set to zero
                        //3D variables
                        double DUMMY1 = 0;
                        double DUMMY2 = 0;
                        double GEW = 0;
                        double SUMGEW = 0.0;
                        for (int n = 0; n < TERA5.Count; n++)  //loop over all ERA5 cells within GRAMM domain
                        {
                            DUMMY1 = Math.Sqrt(Math.Pow(XERA5[n] - Program.X[i] - Program.DDX[i] * 0.5, 2) + Math.Pow(YERA5[n] - Program.Y[j] - Program.DDY[j] * 0.5, 2));
                            if (DUMMY1 < 30000.0)
                            {
                                DUMMY1 = Math.Max(Math.Pow(DUMMY1, 2), 1.0);
                                DUMMY2 = Math.Max(Math.Pow(ZERA5[n] + AH_L - ZSP_L[k], 4) * 0.001, 1.0);
                                GEW = DUMMY1 * DUMMY2 / 1000000.0;
                                SUMGEW += 1 / GEW;
                                T_L[k] += TERA5[n] / GEW;
                                U_L[k] += UERA5[n] / GEW;
                                V_L[k] += VERA5[n] / GEW;
                                QU_L[k] += QERA5[n] * 1000.0 / GEW;
                                CW_L[k] += CLOUDWATERERA5[n] * 1000.0 / GEW;
                            }
                        }
                        //Weighting factors
                        if (SUMGEW != 0)
                        {
                            T_L[k] /= SUMGEW;
                            U_L[k] /= SUMGEW;
                            V_L[k] /= SUMGEW;
                            QU_L[k] /= SUMGEW;
                            CW_L[k] /= SUMGEW;
                        }

                        //2D variables
                        if (k == 1)
                        {
                            SUMGEW = 0.0;
                            DUMMY1 = 0;
                            for (int n = 0; n < CCERA5.Count; n++)  //loop over all ERA5 cells within GRAMM domain
                            {
                                DUMMY1 = Math.Sqrt(Math.Pow(XERA5[n] - Program.X[i] - Program.DDX[i] * 0.5, 2) + Math.Pow(YERA5[n] - Program.Y[j] - Program.DDY[j] * 0.5, 2));
                                if (DUMMY1 < 30000.0)
                                {
                                    DUMMY1 = Math.Max(Math.Pow(DUMMY1, 2), 1.0);
                                    GEW = DUMMY1 / 1000000;
                                    SUMGEW += 1 / GEW;
                                    Program.FW[i][j] += SWERA5[n] / GEW;
                                    //take height of clouds and topography into account
                                    if (AH_L > 6000.0)
                                    {
                                        Program.CLOUDS[i][j] += HCCERA5[n] / GEW;
                                    }
                                    else if (AH_L <= 6000.0 && AH_L > 2000.0)
                                    {
                                        Program.CLOUDS[i][j] += Math.Max(HCCERA5[n], MCCERA5[n]) / GEW;
                                    }
                                    else
                                    {
                                        Program.CLOUDS[i][j] += Math.Max(HCCERA5[n], Math.Max(MCCERA5[n], LCCERA5[n])) / GEW;
                                    }
                                    Program.SNOW[i][j] += SDERA5[n] * 10.0 / GEW;
                                    MSLPinterim[i][j] += MSLPERA5[n] / GEW;
                                }
                            }
                            //Weighting factors
                            if (SUMGEW != 0)
                            {
                                Program.FW[i][j] /= SUMGEW;
                                Program.CLOUDS[i][j] /= SUMGEW;
                                Program.SNOW[i][j] /= SUMGEW;
                                MSLPinterim[i][j] /= SUMGEW;
                                if (Program.ISTAT == 3)
                                {
                                    Program.CLOUDS[i][j] = 0;
                                }
                            }
                        }
                    }
                }
            });

            //calculating the basic state of p, T and rho according to COSMO model, horizontally homegenuous, steady-state
            //seeking for a close approximation of an idealized atmosphere to the actual observed one
            double betamin = 42.0;
            double Tsurfacemin = 288.0;
            double fehlermin = 10000000.0;
            for (int beta = 42; beta <= 42; beta++)
            {
                for (int Tsurface = 210; Tsurface <= 350; Tsurface++)
                {
                    double fehler = 0.0;
                    for (int k = 1; k <= Program.NZ; k++)
                    {
                        double term = Math.Sqrt(1.0 - 2.0 * (float)beta * 9.81 * Program.ZSP[AHMINI][AHMINJ][k] / Program.GASCON / Math.Pow((float)Tsurface, 2));
                        fehler += Math.Pow(Tinterim[AHMINI][AHMINJ][k] - (float)Tsurface * term, 2);
                    }
                    if (fehler < fehlermin)
                    {
                        betamin = beta;
                        Tsurfacemin = Tsurface;
                        fehlermin = fehler;
                    }
                }
            }
            Parallel.For(1, Program.NX + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= Program.NY; j++)
                {
                    double[] Tinter_L = Tinterim[i][j];
                    float[] FACTOR_L = Program.FACTOR[i][j];
                    float[] ZSP_L = Program.ZSP[i][j];
                    double[] TBZ_L = Program.TBZ[i][j];

                    for (int k = 1; k <= Program.NZ; k++)
                    {
                        double term = Math.Sqrt(1.0 - 2.0 * betamin * 9.81 * ZSP_L[k] / Program.GASCON / Math.Pow(Tsurfacemin, 2));
                        double Thelp= Tsurfacemin * term * FACTOR_L[k] - Program.TBZ1;
                        TBZ_L[k] = Thelp + Program.TBZ1;
                        Tinter_L[k] = Tinter_L[k] * FACTOR_L[k] - Program.TBZ1;
                    }
                }
            });

            //calculating the transient state of Q, for simplification the profile above the lowest topography is used
            Parallel.For(1, Program.NX + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= Program.NY; j++)
                {
                    double[] QU_L = QUinterim[i][j];
                    float[] ZSP_L = Program.ZSP[i][j];
                    float[] ZSPmin_L = Program.ZSP[AHMINI][AHMINJ];
                    double[] QUmin_L = QUinterim[AHMINI][AHMINJ];

                    for (int k = 1; k <= Program.NZ; k++)
                    {
                        int n1 = -1;
                        double QUref = 0;
                        for (int n = 1; n <= Program.NZ; n++)
                        {
                            if (ZSPmin_L[n] >= ZSP_L[k])
                            {
                                n1 = n;
                                break;
                            }
                        }
                        if (n1 < 0)
                        {
                            n1 = Program.NZ;
                        }
                        n1 = Math.Max(n1, 2);
                        Linear_interpolation(ZSPmin_L[n1 - 1] - AHMIN, ZSPmin_L[n1] - AHMIN, ZSP_L[k] - AHMIN, QUmin_L[n1 - 1], QUmin_L[n1], ref QUref);
                        QU_L[k] = QUref;
                    }
                }
            });

            /*
            //Geostrophic wind estimation
            //Part I: Geostrophic wind at sea level
            double xgrad = 0;
            double ygrad = 0;
            int numbx = 0;
            int numby = 0;
            Parallel.For(0, MSLPERA5.Count, Program.pOptions, i =>
            {
                for (int j = 0; j < MSLPERA5.Count; j++)
                {
                    if (i != j && Math.Abs(XERA5[i] - XERA5[j]) > 10000)
                    {
                        xgrad += (MSLPERA5[i] - MSLPERA5[j]) / (XERA5[i] - XERA5[j]);
                        numbx++;
                    }
                    if (i != j && Math.Abs(YERA5[i] - YERA5[j]) > 10000)
                    {
                        ygrad += (MSLPERA5[i] - MSLPERA5[j]) / (YERA5[i] - YERA5[j]);
                        numby++;
                    }
                }
            });
            xgrad /= numbx;
            ygrad /= numby;
            double UGsea = -ygrad / 1.2 / Program.FN;
            double VGsea = xgrad / 1.2 / Program.FN;

            //PartII: Geostrophic wind above mixing layer
            int numb = 0;
            double UGtemp = 0;
            double VGtemp = 0;
            Parallel.For(1, Program.NX + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= Program.NY; j++)
                {
                    double[] U_L = Uinterim[i][j];
                    double[] V_L = Vinterim[i][j];
                    float[] ZSP_L = Program.ZSP[i][j];
                    for (int k = Program.NZ; k >= 1; k--)
                    {
                        if (ZSP_L[k] - Program.AHMAX > 2000.0 && ZSP_L[k] < 10000.0)
                        {
                            UGtemp += U_L[k];
                            VGtemp += V_L[k];
                            numb++;
                        }
                    }
                }
            });
            double UGfree = UGtemp / numb;
            double VGfree = VGtemp / numb;
            */

            Parallel.For(1, Program.NX + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= Program.NY; j++)
                {
                    double[] U_L = Uinterim[i][j];
                    double[] V_L = Vinterim[i][j];
                    double[] W_L = Program.W[i][j];
                    double[] UG_L = Program.UG[i][j];
                    double[] VG_L = Program.VG[i][j];
                    float[] UG1_L = Program.UG1[i][j];
                    float[] VG1_L = Program.VG1[i][j];
                    float[] UG2_L = Program.UG2[i][j];
                    float[] VG2_L = Program.VG2[i][j];
                    float[] ZSP_L = Program.ZSP[i][j];
                    for (int k = Program.NZ; k >= 1; k--)
                    {
                        /*
                        if (ZSP_L[k] > AHMAX + 2000)
                        {
                            U_L[k] = UGfree;
                            V_L[k] = VGfree;
                        }
                        else if (ZSP_L[k] <= AHMAX + 2000 && ZSP_L[k] >= AHMAX)
                        {
                            Linear_interpolation(0, AHMAX + 2000, ZSP_L[k], 0, UGfree, ref U_L[k]);
                            Linear_interpolation(0, AHMAX + 2000, ZSP_L[k], 0, VGfree, ref V_L[k]);
                        }
                        else
                        {
                            U_L[k] = 0;
                            V_L[k] = 0;
                        }
                        UG1_L[k] = (float)UG_L[k];
                        UG2_L[k] = (float)UG_L[k];
                        VG1_L[k] = (float)VG_L[k];
                        VG2_L[k] = (float)VG_L[k];
                        */

                        //TEST
                        /*
                        UG_L[k] = (float)UGfree;
                        VG_L[k] = (float)VGfree;
                        UG1_L[k] = (float)UGfree;
                        UG2_L[k] = (float)UGfree;
                        VG1_L[k] = (float)VGfree;
                        VG2_L[k] = (float)VGfree;
                        */
                        Program.UN[i][j][k] = U_L[k];
                        Program.U1N[i][j][k] = U_L[k];
                        Program.U1[i][j][k] = U_L[k];
                        Program.U2N[i][j][k] = U_L[k];
                        Program.U2[i][j][k] = U_L[k];
                        Program.U[i][j][k] = U_L[k];
                        Program.V[i][j][k] = V_L[k];
                        Program.VN[i][j][k] = V_L[k];
                        Program.V1N[i][j][k] = V_L[k];
                        Program.V1[i][j][k] = V_L[k];
                        Program.V2N[i][j][k] = V_L[k];
                        Program.V2[i][j][k] = V_L[k];
                        W_L[k] = 0;
                        Program.WN[i][j][k] = W_L[k];
                        Program.W1N[i][j][k] = W_L[k];
                        Program.W1[i][j][k] = W_L[k];
                        Program.W2N[i][j][k] = W_L[k];
                        Program.W2[i][j][k] = W_L[k];
                        

                        //Program.TN[i][j][k] = TBZ[i][j][k] - Program.TBZ1;
                        //Program.T[i][j][k] = TBZ[i][j][k] - Program.TBZ1;
                    }
                }
            });

            //Values at south and north borders
            Parallel.For(1, Program.NX + 1, Program.pOptions, i =>
            {
                for (int k = 1; k <= Program.NZ; k++)
                {
                    int j = 1;
                    Program.TSSN[i][k] = Tinterim[i][j][k];                    
                    Program.QUSSN[i][k] = QUinterim[i][j][k];

                    Program.WSSN[i][k] = 0;

                    Program.VSSN[i][k] = Vinterim[i][j][k];

                    Program.USSN[i][k] = Uinterim[i][j][k];

                    /*
                    Program.WAT_VAPN[i][j][k] = CWinterim[i][j][k];
                    */

                    j = Program.NY;
                    Program.TSNN[i][k] = Tinterim[i][j][k];
                    Program.QUSNN[i][k] = QUinterim[i][j][k];

                    Program.WSNN[i][k] = 0;

                    Program.USNN[i][k] = Uinterim[i][j][k];

                    Program.VSNN[i][k] = Vinterim[i][j][k];

                    /*
                    Program.WAT_VAPN[i][j][k] = CWinterim[i][j][k];
                    */
                }
            });

            //Values at west and east borders
            Parallel.For(1, Program.NY + 1, Program.pOptions, j =>
            {
                for (int k = 1; k <= Program.NZ; k++)
                {
                    int i = 1;
                    Program.TSWN[j][k] = Tinterim[i][j][k];
                    Program.QUSWN[j][k] = QUinterim[i][j][k];

                    Program.WSWN[j][k] = 0;

                    Program.VSWN[j][k] = Vinterim[i][j][k];

                    Program.USWN[j][k] = Uinterim[i][j][k];

                    /*
                    Program.WAT_VAPN[i][j][k] = CWinterim[i][j][k];
                    */

                    i = Program.NX;
                    Program.TSEN[j][k] = Tinterim[i][j][k];
                    Program.QUSEN[j][k] = QUinterim[i][j][k];

                    Program.WSEN[j][k] = 0;

                    Program.USEN[j][k] = Uinterim[i][j][k];

                    Program.VSEN[j][k] = Vinterim[i][j][k];

                    /*
                    Program.WAT_VAPN[i][j][k] = CWinterim[i][j][k];
                    */
                }
            });            

            //set wetness for water bodies and dry areas (due to the coarse ERA5 data resolution, this is necessary to better account for water bodies and dry areas)
            Parallel.For(1, Program.NX + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= Program.NY; j++)
                {
                    if (Program.ALAMBDA[i][j] > 99.0)
                    {
                        for (int kb = 1; kb <= Program.NZB; kb++)
                        {
                            Program.FW[i][j] = 1.0;
                        }
                    }
                    if (Program.DRY_AREA[i][j] == true)
                    {
                        Program.FW[i][j] = 0.001;
                    }
                }
            });
        }

        //calculate boundary conditions
        public static void ERA5_BoundaryConditions2(List<float> XERA5, List<float> YERA5, List<float> ZERA5, List<float> UERA5, List<float> VERA5, List<float> WERA5, List<float> TERA5, List<float> QERA5, List<float> PERA5, List<float> SDERA5, List<float> CCERA5, List<float> STERA5, List<float> SWERA5, List<float> MSLPERA5, List<float> LCCERA5, List<float> MCCERA5, List<float> HCCERA5, List<float> CLOUDWATERERA5, List<float> SEATEMPERA5)
        {
            Console.WriteLine(" UPDATING BOUNDARY VALUES USING ERA5 DATA");
            double[][][] Uinterim = CreateArray<double[][]>(Program.NX1, () => CreateArray<double[]>(Program.NY1, () => new double[Program.NZ1]));
            double[][][] Vinterim = CreateArray<double[][]>(Program.NX1, () => CreateArray<double[]>(Program.NY1, () => new double[Program.NZ1]));
            double[][][] Tinterim = CreateArray<double[][]>(Program.NX1, () => CreateArray<double[]>(Program.NY1, () => new double[Program.NZ1]));
            double[][][] QUinterim = CreateArray<double[][]>(Program.NX1, () => CreateArray<double[]>(Program.NY1, () => new double[Program.NZ1]));
            double[][][] CWinterim = CreateArray<double[][]>(Program.NX1, () => CreateArray<double[]>(Program.NY1, () => new double[Program.NZ1]));
            double[][] MSLPinterim = CreateArray<double[]>(NX1, () => new double[NY1]);
            object obj = new object();
            Parallel.For(1, Program.NX + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= Program.NY; j++)
                {
                    double[] U_L = Uinterim[i][j];
                    double[] V_L = Vinterim[i][j];
                    double[] T_L = Tinterim[i][j];
                    double[] QU_L = QUinterim[i][j];
                    double[] CW_L = CWinterim[i][j];
                    float[] ZSP_L = Program.ZSP[i][j];
                    float[] FACTOR_L = Program.FACTOR[i][j];
                    float AH_L = Program.AH[i][j];
                    //set all to zero
                    Program.FW[i][j] = 0.0;
                    Program.CLOUDS[i][j] = 0.0;
                    Program.SNOW[i][j] = 0.0;
                    for (int k = 1; k <= Program.NZ; k++)
                    {
                        //set to zero
                        //3D variables
                        double DUMMY1 = 0;
                        double DUMMY2 = 0;
                        double GEW = 0;
                        double SUMGEW = 0.0;
                        for (int n = 0; n < TERA5.Count; n++)  //loop over all ERA5 cells within GRAMM domain
                        {
                            DUMMY1 = Math.Sqrt(Math.Pow(XERA5[n] - Program.X[i] - Program.DDX[i] * 0.5, 2) + Math.Pow(YERA5[n] - Program.Y[j] - Program.DDY[j] * 0.5, 2));
                            if (DUMMY1 < 30000.0)
                            {
                                DUMMY1 = Math.Max(Math.Pow(DUMMY1, 2), 1.0);
                                DUMMY2 = Math.Max(Math.Pow(ZERA5[n] + AH_L - ZSP_L[k], 4) * 0.001, 1.0);
                                GEW = DUMMY1 * DUMMY2 / 1000000.0;
                                SUMGEW += 1 / GEW;
                                T_L[k] += TERA5[n] / GEW;
                                U_L[k] += UERA5[n] / GEW;
                                V_L[k] += VERA5[n] / GEW;
                                QU_L[k] += QERA5[n] * 1000.0 / GEW;
                                CW_L[k] += CLOUDWATERERA5[n] * 1000.0 / GEW;
                            }
                        }
                        //Weighting factors
                        if (SUMGEW != 0)
                        {
                            T_L[k] /= SUMGEW;
                            U_L[k] /= SUMGEW;
                            V_L[k] /= SUMGEW;
                            QU_L[k] /= SUMGEW;
                            CW_L[k] /= SUMGEW;
                        }

                        //2D variables
                        if (k == 1)
                        {
                            SUMGEW = 0.0;
                            DUMMY1 = 0;
                            for (int n = 0; n < CCERA5.Count; n++)  //loop over all ERA5 cells within GRAMM domain
                            {
                                DUMMY1 = Math.Sqrt(Math.Pow(XERA5[n] - Program.X[i] - Program.DDX[i] * 0.5, 2) + Math.Pow(YERA5[n] - Program.Y[j] - Program.DDY[j] * 0.5, 2));
                                if (DUMMY1 < 30000.0)
                                {
                                    DUMMY1 = Math.Max(Math.Pow(DUMMY1, 2), 1.0);
                                    GEW = DUMMY1 / 1000000;
                                    SUMGEW += 1 / GEW;
                                    Program.FW[i][j] += SWERA5[n] / GEW;
                                    //take height of clouds and topography into account
                                    if (AH_L > 6000.0)
                                    {
                                        Program.CLOUDS[i][j] += HCCERA5[n] / GEW;
                                    }
                                    else if (AH_L <= 6000.0 && AH_L > 2000.0)
                                    {
                                        Program.CLOUDS[i][j] += Math.Max(HCCERA5[n], MCCERA5[n]) / GEW;
                                    }
                                    else
                                    {
                                        Program.CLOUDS[i][j] += Math.Max(HCCERA5[n], Math.Max(MCCERA5[n], LCCERA5[n])) / GEW;
                                    }
                                    Program.SNOW[i][j] += SDERA5[n] * 10.0 / GEW;
                                    MSLPinterim[i][j] += MSLPERA5[n] / GEW;
                                }
                            }
                            //Weighting factors
                            if (SUMGEW != 0)
                            {
                                Program.FW[i][j] /= SUMGEW;
                                Program.CLOUDS[i][j] /= SUMGEW;
                                Program.SNOW[i][j] /= SUMGEW;
                                MSLPinterim[i][j] /= SUMGEW;
                                if (Program.ISTAT == 3)
                                {
                                    Program.CLOUDS[i][j] = 0;
                                }
                            }
                        }
                    }
                }
            });

            //calculating the basic state of p, T and rho according to COSMO model, horizontally homegenuous, steady-state
            //seeking for a close approximation of an idealized atmosphere to the actual observed one
            double betamin = 42.0;
            double Tsurfacemin = 288.0;
            double fehlermin = 10000000.0;
            for (int beta = 42; beta <= 42; beta++)
            {
                for (int Tsurface = 210; Tsurface <= 350; Tsurface++)
                {
                    double fehler = 0.0;
                    for (int k = 1; k <= Program.NZ; k++)
                    {
                        double term = Math.Sqrt(1.0 - 2.0 * (float)beta * 9.81 * Program.ZSP[AHMINI][AHMINJ][k] / Program.GASCON / Math.Pow((float)Tsurface, 2));
                        fehler += Math.Pow(Tinterim[AHMINI][AHMINJ][k] - (float)Tsurface * term, 2);
                    }
                    if (fehler < fehlermin)
                    {
                        betamin = beta;
                        Tsurfacemin = Tsurface;
                        fehlermin = fehler;
                    }
                }
            }
            Parallel.For(1, Program.NX + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= Program.NY; j++)
                {
                    double[] Tinter_L = Tinterim[i][j];
                    float[] FACTOR_L = Program.FACTOR[i][j];
                    float[] ZSP_L = Program.ZSP[i][j];
                    double[] TBZ_L = Program.TBZ[i][j];

                    for (int k = 1; k <= Program.NZ; k++)
                    {
                        double term = Math.Sqrt(1.0 - 2.0 * betamin * 9.81 * ZSP_L[k] / Program.GASCON / Math.Pow(Tsurfacemin, 2));
                        double Thelp = Tsurfacemin * term * FACTOR_L[k] - Program.TBZ1;
                        TBZ_L[k] = Thelp + Program.TBZ1;
                        Tinter_L[k] = Tinter_L[k] * FACTOR_L[k] - Program.TBZ1;
                    }
                }
            });

            //calculating the transient state of Q, for simplification the profile above the lowest topography is used
            Parallel.For(1, Program.NX + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= Program.NY; j++)
                {
                    double[] QU_L = QUinterim[i][j];
                    float[] ZSP_L = Program.ZSP[i][j];
                    float[] ZSPmin_L = Program.ZSP[AHMINI][AHMINJ];
                    double[] QUmin_L = QUinterim[AHMINI][AHMINJ];

                    for (int k = 1; k <= Program.NZ; k++)
                    {
                        int n1 = -1;
                        double QUref = 0;
                        for (int n = 1; n <= Program.NZ; n++)
                        {
                            if (ZSPmin_L[n] >= ZSP_L[k])
                            {
                                n1 = n;
                                break;
                            }
                        }
                        if (n1 < 0)
                        {
                            n1 = Program.NZ;
                        }
                        n1 = Math.Max(n1, 2);
                        Linear_interpolation(ZSPmin_L[n1 - 1] - AHMIN, ZSPmin_L[n1] - AHMIN, ZSP_L[k] - AHMIN, QUmin_L[n1 - 1], QUmin_L[n1], ref QUref);
                        QU_L[k] = QUref;
                    }
                }
            });

            //Values at south and north borders
            Parallel.For(1, Program.NX + 1, Program.pOptions, i =>
            {
                for (int k = 1; k <= Program.NZ; k++)
                {
                    int j = 1;
                    Program.TSSN[i][k] = Tinterim[i][j][k];
                    Program.QUSSN[i][k] = QUinterim[i][j][k];

                    Program.WSSN[i][k] = 0;

                    Program.VSSN[i][k] = Vinterim[i][j][k];

                    Program.USSN[i][k] = Uinterim[i][j][k];

                    /*
                    Program.WAT_VAPN[i][j][k] = CWinterim[i][j][k];
                    */

                    j = Program.NY;
                    Program.TSNN[i][k] = Tinterim[i][j][k];
                    Program.QUSNN[i][k] = QUinterim[i][j][k];

                    Program.WSNN[i][k] = 0;

                    Program.USNN[i][k] = Uinterim[i][j][k];

                    Program.VSNN[i][k] = Vinterim[i][j][k];

                    /*
                    Program.WAT_VAPN[i][j][k] = CWinterim[i][j][k];
                    */
                }
            });

            //Values at west and east borders
            Parallel.For(1, Program.NY + 1, Program.pOptions, j =>
            {
                for (int k = 1; k <= Program.NZ; k++)
                {
                    int i = 1;
                    Program.TSWN[j][k] = Tinterim[i][j][k];
                    Program.QUSWN[j][k] = QUinterim[i][j][k];

                    Program.WSWN[j][k] = 0;

                    Program.VSWN[j][k] = Vinterim[i][j][k];

                    Program.USWN[j][k] = Uinterim[i][j][k];

                    /*
                    Program.WAT_VAPN[i][j][k] = CWinterim[i][j][k];
                    */

                    i = Program.NX;
                    Program.TSEN[j][k] = Tinterim[i][j][k];
                    Program.QUSEN[j][k] = QUinterim[i][j][k];

                    Program.WSEN[j][k] = 0;

                    Program.USEN[j][k] = Uinterim[i][j][k];

                    Program.VSEN[j][k] = Vinterim[i][j][k];
                }
            });            

            //set wetness for water bodies and dry areas (due to the coarse ERA5 data resolution, this is necessary to better account for water bodies and dry areas )
            Parallel.For(1, Program.NX + 1, Program.pOptions, i =>
            {
                for (int j = 1; j <= Program.NY; j++)
                {
                    if (Program.ALAMBDA[i][j] > 99.0)
                    {
                        for (int kb = 1; kb <= Program.NZB; kb++)
                        {
                            Program.FW[i][j] = 1.0;
                        }
                    }
                    if (Program.DRY_AREA[i][j] == true)
                    {
                        Program.FW[i][j] = 0.001;
                    }
                }
            });
        }
    }
}
