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
using System.Globalization;
using System.IO;

namespace GRAMM_2001
{
    partial class Program
    {
        /// <summary>
        /// Write GRAMM online data 
        /// </summary>
        /// <param name="NI"></param>
        /// <param name="NJ"></param>
        /// <param name="NK"></param>
        public static void GrammOnline(int NI, int NJ, int NK)
        { // Schreibe Berechnungsdaten für GRAMM Online

            string config1 = "";
            string config2 = "";
            string filename = "";

            if (IWetter_Console_First > 1) // 11.4.17 Ku - do not write GRAMM Online for multiple instances
            {
                return;
            }

            for (int Ausgabe = 1; Ausgabe <= 29; Ausgabe++) // insgesamt 29 Fälle
            {
                switch (Ausgabe) // für jeden Fall den Filenamen festlegen
                {
                    case 1: //ausgabe des u,v-windfeldes fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        filename = "uv";
                        break;
                    case 2: // ausgabe der u-windfeldkomponente fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        filename = "u";
                        break;
                    case 3: // ausgabe der horizontalen windgeschwindigkeit fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        filename = "speed";
                        break;
                    case 4: // ausgabe der horizontalen windgeschwindigkeit v fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        filename = "v";
                        break;
                    case 5: // ausgabe der vertikalen windgeschwindigkeit fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        filename = "w";
                        break;
                    case 6: // ausgabe der absoluten Temperatur fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        filename = "tabs";
                        break;
                    case 7: // ausgabe der potentiellen Temperatur fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        filename = "tpot";
                        break;
                    case 8: // ausgabe der spez. Feuchte fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        filename = "hum";
                        break;
                    case 9: // ausgabe des nicht hydrostat. Drucks fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        filename = "nhp";
                        break;
                    case 10: // ausgabe der turbulenten kinetischen Energie fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        filename = "tke";
                        break;
                    case 11: // ausgabe der Dissipation fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        filename = "dis";
                        break;
                    case 12: //ausgabe der Globalstrahlung fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        filename = "glob";
                        break;
                    case 13: //ausgabe der terrestrischen Ausstrahlung fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        filename = "terr";
                        break;
                    case 14: // ausgabe des sensiblen Waermeflusses fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        filename = "sensheat";
                        break;
                    case 15: // ausgabe des latenten Waermeflusses fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        filename = "latheat";
                        break;
                    case 16: // ausgabe der Schubspannungsgeschwindigkeit fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        filename = "fricvel";
                        break;
                    case 17: // ausgabe der inversen Monin-Obukhov Laenge fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        filename = "inverseMO";
                        break;
                    case 18: // ausgabe der Bodentemperatur fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        filename = "surfTemp";
                        break;
                    case 19: // ausgabe der Stabilitätsklasse fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        filename = "stabilityclass";
                        break;
                    case 20: // ausgabe der horizontalen Windgeschwindigkeit fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        filename = "vp_Wind-speed";
                        break;
                    case 21: // ausgabe der West/Ost Windkomponente fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        filename = "vp_u-wind-component";
                        break;
                    case 22: // ausgabe der Nord/sued Windkomponente fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        filename = "vp_v-wind-component";
                        break;
                    case 23: // ausgabe der vertikalen Windkomponente fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        filename = "vp_w-wind-component";
                        break;
                    case 24: // ausgabe der Temperatur fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        filename = "vp_Temperature";
                        break;
                    case 25: //  ausgabe der potentiellen Temperatur fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        filename = "vp_potential-Temperature";
                        break;
                    case 26: // ausgabe der Feuchtigkeit fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        filename = "vp_Humidity";
                        break;
                    case 27: // ausgabe des nicht-hydrostatischen Drucks fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        filename = "vp_non-hydr-pressure";
                        break;
                    case 28: // ausgabe der turbulenten kinetischen Energie fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        filename = "vp_Turbulent-kinetic-energy";
                        break;
                    case 29: // ausgabe der Dissipation fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        filename = "vp_Dissipation";
                        break;
                }

                if (File.Exists(filename + ".txt")) // muss file ausgeschrieben werden?
                {
                    GRAMM_Online_flag = true;
                    try
                    {
                        using (StreamReader rd = new StreamReader(filename + ".txt")) // dann einlesen des gewuenschten Ausgabelevels ueber Grund
                        {
                            config1 = rd.ReadLine();   // versuche string zu lesen - bei Fehler ans Ende des using-Statements
                            config2 = rd.ReadLine();   // auch zweite Zeile vorhanden?
                        }

                        int LX = 1; // Defaultwert
                        Int32.TryParse(config1, out LX); // Layer oder LX
                        int LY = 1; // default
                        Int32.TryParse(config2, out LY); // immer LY

                        string writefile;
                        if (Ausgabe < 20) // man könnte auf den Textteil "vp_" prüfen, ist unabhängiger, aber nicht so schnell
                        {
                            writefile = filename + "_GRAMM.txt";
                        }
                        else
                        {
                            writefile = "GRAMM-" + filename + ".txt";
                        }

                        using (StreamWriter wt = new StreamWriter(writefile)) // ausgabe des windfeldes - using = try-catch
                        {
                            if (Ausgabe < 20)
                            {
                                wt.WriteLine("ncols         " + NI.ToString(CultureInfo.InvariantCulture));
                                wt.WriteLine("nrows         " + NJ.ToString(CultureInfo.InvariantCulture));
                                wt.WriteLine("xllcorner     " + Program.IKOOA.ToString(CultureInfo.InvariantCulture));
                                wt.WriteLine("yllcorner     " + Program.JKOOA.ToString(CultureInfo.InvariantCulture));
                                wt.WriteLine("cellsize      " + Program.DDXImm[1].ToString(CultureInfo.InvariantCulture));
                                wt.WriteLine("NODATA_value  " + "-9999");

                                for (int jj = NJ; jj >= 1; jj--)
                                {
                                    for (int o = 1; o <= NI; o++)
                                    {
                                        if (Ausgabe == 1) // Sonderfall Ausgabe 1
                                        {
                                            //wt.Write(Convert.ToString(Math.Round(Program.U[o][jj][LX], 3)).Replace(Program.decsep, ".") + " " + Convert.ToString(Math.Round(Program.V[o][jj][LX], 3)).Replace(Program.decsep, ".") + " ");
                                            wt.Write(Math.Round(Program.U[o][jj][LX], 3).ToString(CultureInfo.InvariantCulture) + " " + Math.Round(Program.V[o][jj][LX], 3).ToString(System.Globalization.CultureInfo.InvariantCulture) + " ");
                                        }
                                        else
                                        {
                                            //wt.Write(Convert.ToString(value(Ausgabe,o,jj,LX,LY)).Replace(Program.decsep, ".") + " ");
                                            wt.Write(value(Ausgabe, o, jj, LX, LY).ToString(CultureInfo.InvariantCulture) + " ");
                                        }
                                    }
                                    wt.WriteLine();
                                }
                            }

                            else // Ausgabe ab Variante 20
                            {
                                for (int k = 1; k <= Program.NZ; k++)
                                {
                                    //wt.Write(Convert.ToString(value(Ausgabe,0,k,LX,LY)).Replace(Program.decsep, ".") + " " + Convert.ToString(value(Ausgabe,1,k,LX,LY)).Replace(Program.decsep, ".") + " ");
                                    wt.WriteLine(Math.Round(Program.ZSPImm[LX][LY][k], 3).ToString(CultureInfo.InvariantCulture) + " " + value(Ausgabe, 1, k, LX, LY).ToString(CultureInfo.InvariantCulture) + " ");
                                }
                            }

                        } // StreamWriter Ende
                    }
                    catch
                    {
                    }
                } // FileExist() Ende

            } // Ausgabenschleife Ende

        } // GRAMMONLINE Ende

        private static double value(int Ausgabe, int o, int jj, int LX, int LY)
        { // Ermittle Ausgabewerte und Zahlen-Rundung
            try
            {
                double val;
                switch (Ausgabe)
                {
                    case 1: // Ausgabe des u,v-windfeldes fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        return 0;

                    case 2: // Ausgabe der u-windfeldkomponente fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        return Math.Round(Program.WAT_VAPN[o][jj][LX], 3);

                    case 3: // Ausgabe der horizontalen windgeschwindigkeit fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        return Math.Round(Math.Sqrt(Program.U[o][jj][LX] * Program.U[o][jj][LX] + Program.V[o][jj][LX] * Program.V[o][jj][LX]), 3);

                    case 4: //Ausgabe der horizontalen windgeschwindigkeit v fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        return Math.Round(Program.V[o][jj][LX], 3);

                    case 5: //Ausgabe der vertikalen windgeschwindigkeit fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        return Math.Round(Program.W[o][jj][LX], 3);

                    case 6: //Ausgabe der absoluten Temperatur fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        return Math.Round(Program.TABS[o][jj][LX] - 273.15, 3);

                    case 7: // Ausgabe der potentiellen Temperatur fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        return Math.Round(Program.TN[o][jj][LX] + Program.TBZ1, 3);

                    case 8: //ausgabe der spez. Feuchte fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        return Math.Round(Program.QUN[o][jj][LX], 3);

                    case 9: //ausgabe des nicht hydrostat. Drucks fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        if (Program.ICPSI)
                        {
                            return Math.Round(Program.TP[o][jj][LX], 3);
                        }
                        else
                        {
                            return Math.Round(Program.DP[o][jj][LX], 3);
                        }

                    case 10: //ausgabe der turbulenten kinetischen Energie fuer GRAMM Online Analyse mit Benutzeroberflaeche
                             //return Math.Round(Program.TEN[o][jj][LX], 3);
                        return Math.Round(Program.CLOUDS[o][jj], 3);

                    case 11: //ausgabe der Dissipation fuer GRAMM Online Analyse mit Benutzeroberflaeche
                             //return Math.Round(Program.DISSN[o][jj][LX], 3);
                        return Math.Round(Program.SNOW[o][jj], 3);

                    case 12: //ausgabe der Globalstrahlung fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        return Math.Round(Program.GLOBRAD[o][jj], 3);

                    case 13: // ausgabe der terrestrischen Ausstrahlung fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        return Math.Round(Program.EPSG[o][jj] * (Program.RL[o][jj] - Program.SIGMA * Math.Pow(Program.TB[o][jj][2], 4)), 3);

                    case 14: //  ausgabe des sensiblen Waermeflusses fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        val = Program.WQU[o][jj] / Program.DDXImm[o] / Program.DDYImm[jj];
                        if (double.IsInfinity(val) || double.IsNaN(val))
                        {
                            val = -9999; // Error
                        }

                        return Math.Round(val, 3);

                    case 15: // ausgabe des latenten Waermeflusses fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        return Math.Round(Program.XWQ[o][jj] * Program.RHO[o][jj][1] * Program.UST[o][jj] * Program.ALW * (Program.QU[o][jj][1] - Program.QUG[o][jj]) / 1000, 3);

                    case 16: //  ausgabe der Schubspannungsgeschwindigkeit fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        return Math.Round(Program.UST[o][jj], 3);

                    case 17: // ausgabe der inversen Monin-Obukhov Laenge fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        val = 1 / Program.OL[o][jj];
                        if (double.IsInfinity(val) || double.IsNaN(val))
                        {
                            val = -9999; // Error
                        }

                        return Math.Round(val, 3);

                    case 18: // ausgabe der Bodentemperatur fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        return Math.Round(Program.TB[o][jj][2], 3);

                    case 19: // ausgabe der Stabilitätsklasse fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        return Program.stabilityclass[o][jj];

                    case 20: // ausgabe der horizontalen Windgeschwindigkeit fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        return Math.Round(Math.Sqrt(Math.Pow(Program.U[LX][LY][jj], 2) + Math.Pow(Program.V[LX][LY][jj], 2)), 3);

                    case 21: // ausgabe der West/Ost Windkomponente fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        return Math.Round(Program.U[LX][LY][jj], 3);

                    case 22: // ausgabe der Nord/sued Windkomponente fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        return Math.Round(Program.V[LX][LY][jj], 3);

                    case 23: // ausgabe der vertikalen Windkomponente fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        return Math.Round(Program.W[LX][LY][jj], 3);

                    case 24: // ausgabe der Temperatur fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        return Math.Round(Program.TABS[LX][LY][jj] - 273.15, 3);

                    case 25: // ausgabe der potentiellen Temperatur fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        return Math.Round(Program.TN[LX][LY][jj] + Program.TBZ1, 3);

                    case 26: // ausgabe der Feuchtigkeit fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        return Math.Round(Program.QUN[LX][LY][jj], 3);

                    case 27: // ausgabe des nicht-hydrostatischen Drucks fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        if (Program.ICPSI)
                        {
                            return Math.Round(Program.TP[LX][LY][jj], 3);
                        }
                        else
                        {
                            return Math.Round(Program.DP[LX][LY][jj], 3);
                        }

                    case 28: // ausgabe der turbulenten kinetischen Energie fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        return Math.Round(Program.TEN[LX][LY][jj], 3);

                    case 29: // ausgabe der Dissipation fuer GRAMM Online Analyse mit Benutzeroberflaeche
                        return Math.Round(Program.DISSN[LX][LY][jj], 3);

                    default:
                        return -9999; // ERROR
                }
            }
            catch
            {
                return -9999; // ERROR
            }
        } // value() Ende

    }
}