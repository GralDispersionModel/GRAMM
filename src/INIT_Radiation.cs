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


namespace GRAMM_2001
{
    partial class Program
    {
        public static void INIT_Radiation(ref bool LGEOM, ref bool LGEOMW, ref bool LGEOMR, ref int ISTUD, ref double AMIND,
        ref int IMIND, ref double ASECD, ref int ISECD, List<int>[] Month_List)
        {
            int NX_Radiat_Guess = Math.Max(NX / 4, 3 * nr_cell_smooth);
            int NY_Radiat_Guess = Math.Max(NY / 4, 3 * nr_cell_smooth);
            Program.ISOL = 1; // thin clouds -> default value
            Solar_reduction_factor = 1.0F; // full solar radiation = default value

            if (((METEO == "Y") || (METEO == "y")) && (Program.ISTAT == 0))
            {
                if (IMETSTR == 1)
                {
                    LGEOM = false;
                    LGEOMW = false;
                    LGEOMR = false;
                    IMETSTR = 1;
                }
                else
                {
                    LGEOM = true;
                    LGEOMW = true;
                    LGEOMR = false;
                }
                if (ISTAT == 0)
                {
                    //option: dynamic sun rise for convective classes (1-3)
                    if (Program.meteopgt_nr > 0)
                    {
                        if (AKLA == 1)
                        {
                            ITAG = 21;
                            if (BGRAD > 0)
                                IMON = 6;
                            else
                                IMON = 12;
                            IJAHR = 2006;
                            ISTU = 6;
                            IMIN = 0;

                            int month_setting = set_month_type(Windspeed_meteopgt, AKLA, BGRAD);
                            bool exit = false;
                            float rsolg_max = 0;

                            //find the month and time of the year, where the radiation equals xx W/m? as required by SC 1
                            foreach (int imonth_var in Month_List[month_setting])
                            {
                                for (int iday_var = 6; iday_var < 26; iday_var += 2)
                                {
                                    float rsolg_daymax = 0;
                                    for (int istu_var = 4; istu_var < 14; istu_var += 2)
                                    {
                                        TJETZT = (float)istu_var * 3600 + (float)IMIN * 60;
                                        ISTUD = Convert.ToInt32(Math.Floor((TJETZT + 10) / 3600d));
                                        AMIND = (TJETZT / 3600 - (float)istu_var) * 60;
                                        IMIND = Convert.ToInt32(Math.Floor(AMIND));
                                        ASECD = (AMIND - (float)IMIND) * 60;
                                        ISECD = Convert.ToInt32(Math.Floor(ASECD));

                                        RadiationModel.RADIATRAD(LGEOM, false, false, iday_var, imonth_var, IJAHR, TJETZT, NX_Radiat_Guess, NY_Radiat_Guess, 3);
                                        Console.Write(iday_var.ToString() + "." + imonth_var.ToString() + "  -  " + istu_var.ToString() + ":00        ");
                                        Console.SetCursorPosition(0, Console.CursorTop);

                                        float rsolg_mean = Mean_Array_Value(GLOBRAD, NX_Radiat_Guess, NY_Radiat_Guess);
                                        rsolg_daymax = Math.Max(rsolg_daymax, rsolg_mean);

                                        if (Windspeed_meteopgt <= 2)
                                        {
                                            if (rsolg_mean > 800)
                                            {
                                                ISTU = istu_var;
                                                ITAG = iday_var;
                                                IMON = imonth_var;
                                                iday_var = 31;
                                                exit = true;
                                                break;
                                            }
                                        }
                                        else if (Windspeed_meteopgt > 2)
                                        {
                                            if (rsolg_mean > 950)
                                            {
                                                ISTU = istu_var;
                                                ITAG = iday_var;
                                                IMON = imonth_var;
                                                iday_var = 31;
                                                exit = true;
                                                break;
                                            }
                                        }
                                    }

                                    // find the day with the max. radiation
                                    rsolg_max = Math.Max(rsolg_max, rsolg_daymax);
                                    if ((rsolg_daymax + 1) < rsolg_max)
                                    {
                                        ISTU = 12;
                                        ITAG = iday_var;
                                        IMON = imonth_var;
                                        iday_var = 31;
                                        exit = true;
                                    }
                                }
                                if (exit)
                                    break;
                            }
                            //reset hour to 6 o'clock as the sun rises dynamically
                            ISTU = 6;
                        }
                        else if (AKLA == 2)
                        {
                            ITAG = 21;
                            IMON = 3;
                            IJAHR = 2006;
                            ISTU = 6;
                            IMIN = 0;

                            int month_setting = set_month_type(Windspeed_meteopgt, AKLA, BGRAD);
                            bool exit = false;

                            //find the month and time of the year, where the radiation equals xx W/m? as required by SC 2
                            foreach (int imonth_var in Month_List[month_setting])
                            {
                                for (int iday_var = 1; iday_var < 29; iday_var += 2)
                                {
                                    for (int istu_var = 4; istu_var < 14; istu_var += 2)
                                    {
                                        TJETZT = (float)istu_var * 3600 + (float)IMIN * 60;
                                        ISTUD = Convert.ToInt32(Math.Floor((TJETZT + 10) / 3600d));
                                        AMIND = (TJETZT / 3600 - (float)istu_var) * 60;
                                        IMIND = Convert.ToInt32(Math.Floor(AMIND));
                                        ASECD = (AMIND - (float)IMIND) * 60;
                                        ISECD = Convert.ToInt32(Math.Floor(ASECD));

                                        RadiationModel.RADIATRAD(LGEOM, false, false, iday_var, imonth_var, IJAHR, TJETZT, NX_Radiat_Guess, NY_Radiat_Guess, 3);
                                        Console.Write(iday_var.ToString() + "." + imonth_var.ToString() + "  -  " + istu_var.ToString() + ":00        ");
                                        Console.SetCursorPosition(0, Console.CursorTop);

                                        float rsolg_mean = Mean_Array_Value(GLOBRAD, NX_Radiat_Guess, NY_Radiat_Guess);

                                        if (Windspeed_meteopgt <= 2)
                                        {
                                            if ((rsolg_mean > 200) && (rsolg_mean <= 600))
                                            {
                                                ISTU = istu_var;
                                                ITAG = iday_var;
                                                IMON = imonth_var;
                                                iday_var = 31;
                                                exit = true;
                                                break;
                                            }
                                        }
                                        else if ((Windspeed_meteopgt > 2) && (Windspeed_meteopgt <= 3))
                                        {
                                            if ((rsolg_mean > 690) && (rsolg_mean <= 900))
                                            {
                                                ISTU = istu_var;
                                                ITAG = iday_var;
                                                IMON = imonth_var;
                                                iday_var = 31;
                                                exit = true;
                                                break;
                                            }
                                        }
                                        else if (Windspeed_meteopgt > 3)
                                        {
                                            if (rsolg_mean > 690)
                                            {
                                                ISTU = istu_var;
                                                ITAG = iday_var;
                                                IMON = imonth_var;
                                                iday_var = 31;
                                                exit = true;
                                                break;
                                            }
                                        }
                                    }
                                }
                                if (exit)
                                    break;
                            }

                            //reset hour to 12 o'clock as the sun goes down dynamically
                            ISTU = 12;
                        }
                        else if (AKLA == 3)
                        {
                            ITAG = 21;
                            IMON = 3;
                            IJAHR = 2006;
                            ISTU = 6;
                            IMIN = 0;

                            int month_setting = set_month_type(Windspeed_meteopgt, AKLA, BGRAD);
                            bool exit = false;

                            //find the month and time of the year, where the radiation equals xx W/m? as required by SC 3
                            foreach (int imonth_var in Month_List[month_setting])
                            {
                                for (int iday_var = 1; iday_var < 29; iday_var += 2)
                                {
                                    for (int istu_var = 4; istu_var < 13; istu_var += 2)
                                    {
                                        TJETZT = (float)istu_var * 3600 + (float)IMIN * 60;
                                        ISTUD = Convert.ToInt32(Math.Floor((TJETZT + 10) / 3600d));
                                        AMIND = (TJETZT / 3600 - (float)istu_var) * 60;
                                        IMIND = Convert.ToInt32(Math.Floor(AMIND));
                                        ASECD = (AMIND - (float)IMIND) * 60;
                                        ISECD = Convert.ToInt32(Math.Floor(ASECD));
                                        RadiationModel.RADIATRAD(LGEOM, false, false, iday_var, imonth_var, IJAHR, TJETZT, 3, 3, 3);
                                        Console.Write(iday_var.ToString() + "." + imonth_var.ToString() + "  -  " + istu_var.ToString() + ":00        ");
                                        Console.SetCursorPosition(0, Console.CursorTop);

                                        if (Windspeed_meteopgt <= 2)
                                        {
                                            if ((GLOBRAD[2][2] > 100) && (GLOBRAD[2][2] <= 200))
                                            {
                                                ISTU = istu_var;
                                                ITAG = iday_var;
                                                IMON = imonth_var;
                                                iday_var = 31;
                                                exit = true;
                                                break;
                                            }
                                        }
                                        else if ((Windspeed_meteopgt > 2) && (Windspeed_meteopgt <= 3))
                                        {
                                            if ((GLOBRAD[2][2] > 150) && (GLOBRAD[2][2] <= 650))
                                            {
                                                ISTU = istu_var;
                                                ITAG = iday_var;
                                                IMON = imonth_var;
                                                iday_var = 31;
                                                exit = true;
                                                break;
                                            }
                                        }
                                        else if ((Windspeed_meteopgt > 3) && (Windspeed_meteopgt <= 5))
                                        {
                                            if ((GLOBRAD[2][2] > 180) && (GLOBRAD[2][2] <= 650))
                                            {
                                                ISTU = istu_var;
                                                ITAG = iday_var;
                                                IMON = imonth_var;
                                                iday_var = 31;
                                                exit = true;
                                                break;
                                            }
                                        }
                                        else if (Windspeed_meteopgt > 5)
                                        {
                                            if (GLOBRAD[2][2] > 900)
                                            {
                                                ISTU = istu_var;
                                                ITAG = iday_var;
                                                IMON = imonth_var;
                                                iday_var = 31;
                                                exit = true;
                                                break;
                                            }
                                        }
                                    }
                                }
                                if (exit)
                                    break;
                            }

                            //reset hour to 9 o'clock as the sun proceedes dynamically
                            ISTU = 9;
                        }
                    } // option dynamic sunrise end
                    else
                    {
                        if (AKLA == 1)
                        {
                            ITAG = 21;
                            if (BGRAD > 0)
                                IMON = 6;
                            else
                                IMON = 12;
                            IJAHR = 2006;
                            ISTU = 12;
                            IMIN = 0;

                            int month_setting = set_month_type(Windspeed_meteopgt, AKLA, BGRAD);
                            bool exit = false;
                            float rsolg_max = 0;

                            //find the month and time of the year, where the radiation equals xx W/m? as required by SC 1
                            foreach (int imonth_var in Month_List[month_setting])
                            {
                                for (int iday_var = 6; iday_var < 26; iday_var += 2)
                                {
                                    float rsolg_daymax = 0;
                                    for (int istu_var = 4; istu_var < 14; istu_var += 2)
                                    {
                                        TJETZT = (float)istu_var * 3600 + (float)IMIN * 60;
                                        ISTUD = Convert.ToInt32(Math.Floor((TJETZT + 10) / 3600d));
                                        AMIND = (TJETZT / 3600 - (float)istu_var) * 60;
                                        IMIND = Convert.ToInt32(Math.Floor(AMIND));
                                        ASECD = (AMIND - (float)IMIND) * 60;
                                        ISECD = Convert.ToInt32(Math.Floor(ASECD));
                                        RadiationModel.RADIATRAD(LGEOM, false, false, iday_var, imonth_var, IJAHR, TJETZT, NX_Radiat_Guess, NY_Radiat_Guess, 3);
                                        Console.Write(iday_var.ToString() + "." + imonth_var.ToString() + "  -  " + istu_var.ToString() + ":00        ");
                                        Console.SetCursorPosition(0, Console.CursorTop);

                                        float rsolg_mean = Mean_Array_Value(GLOBRAD, NX_Radiat_Guess, NY_Radiat_Guess);
                                        rsolg_daymax = Math.Max(rsolg_daymax, rsolg_mean);

                                        if (Windspeed_meteopgt <= 2)
                                        {
                                            if (rsolg_mean > 800)
                                            {
                                                ISTU = istu_var;
                                                ITAG = iday_var;
                                                IMON = imonth_var;
                                                iday_var = 31;
                                                exit = true;
                                                break;
                                            }
                                        }
                                        else if (Windspeed_meteopgt > 2)
                                        {
                                            if (rsolg_mean > 950)
                                            {
                                                ISTU = istu_var;
                                                ITAG = iday_var;
                                                IMON = imonth_var;
                                                iday_var = 31;
                                                exit = true;
                                                break;
                                            }
                                        }
                                    }

                                    // find the day with the max. radiation
                                    rsolg_max = Math.Max(rsolg_max, rsolg_daymax);
                                    if ((rsolg_daymax + 1) < rsolg_max)
                                    {
                                        ISTU = 12;
                                        ITAG = iday_var;
                                        IMON = imonth_var;
                                        iday_var = 31;
                                        exit = true;
                                    }
                                }
                                if (exit)
                                    break;
                            }
                        }
                        else if (AKLA == 2)
                        {
                            ITAG = 21;
                            IMON = 3;
                            IJAHR = 2006;
                            ISTU = 12;
                            IMIN = 0;

                            int month_setting = set_month_type(Windspeed_meteopgt, AKLA, BGRAD);
                            bool exit = false;

                            //find the month and time of the year, where the radiation equals xx W/m? as required by SC 2
                            foreach (int imonth_var in Month_List[month_setting])
                            {
                                for (int iday_var = 1; iday_var < 29; iday_var += 2)
                                {
                                    for (int istu_var = 4; istu_var < 14; istu_var += 2)
                                    {
                                        TJETZT = (float)istu_var * 3600 + (float)IMIN * 60;
                                        ISTUD = Convert.ToInt32(Math.Floor((TJETZT + 10) / 3600d));
                                        AMIND = (TJETZT / 3600 - (float)istu_var) * 60;
                                        IMIND = Convert.ToInt32(Math.Floor(AMIND));
                                        ASECD = (AMIND - (float)IMIND) * 60;
                                        ISECD = Convert.ToInt32(Math.Floor(ASECD));

                                        RadiationModel.RADIATRAD(LGEOM, false, false, iday_var, imonth_var, IJAHR, TJETZT, NX_Radiat_Guess, NY_Radiat_Guess, 3);
                                        Console.Write(iday_var.ToString() + "." + imonth_var.ToString() + "  -  " + istu_var.ToString() + ":00        ");
                                        Console.SetCursorPosition(0, Console.CursorTop);

                                        float rsolg_mean = Mean_Array_Value(GLOBRAD, NX_Radiat_Guess, NY_Radiat_Guess);

                                        if (Windspeed_meteopgt <= 2)
                                        {
                                            if ((rsolg_mean > 200) && (rsolg_mean <= 600))
                                            {
                                                ISTU = istu_var;
                                                ITAG = iday_var;
                                                IMON = imonth_var;
                                                iday_var = 31;
                                                exit = true;
                                                break;
                                            }
                                        }
                                        else if ((Windspeed_meteopgt > 2) && (Windspeed_meteopgt <= 3))
                                        {
                                            if ((rsolg_mean > 690) && (rsolg_mean <= 900))
                                            {
                                                ISTU = istu_var;
                                                ITAG = iday_var;
                                                IMON = imonth_var;
                                                iday_var = 31;
                                                exit = true;
                                                break;
                                            }
                                        }
                                        else if (Windspeed_meteopgt > 3)
                                        {
                                            if (rsolg_mean > 690)
                                            {
                                                ISTU = istu_var;
                                                ITAG = iday_var;
                                                IMON = imonth_var;
                                                iday_var = 31;
                                                exit = true;
                                                break;
                                            }
                                        }
                                    }
                                }
                                if (exit)
                                    break;
                            }
                        }
                        else if (AKLA == 3)
                        {
                            ITAG = 21;
                            IMON = 3;
                            IJAHR = 2006;
                            ISTU = 10;
                            IMIN = 0;

                            int month_setting = set_month_type(Windspeed_meteopgt, AKLA, BGRAD);
                            bool exit = false;

                            //find the month and time of the year, where the radiation equals xx W/m? as required by SC 3
                            foreach (int imonth_var in Month_List[month_setting])
                            {
                                for (int iday_var = 1; iday_var < 29; iday_var += 2)
                                {
                                    for (int istu_var = 4; istu_var < 13; istu_var += 2)
                                    {
                                        TJETZT = (float)istu_var * 3600 + (float)IMIN * 60;
                                        ISTUD = Convert.ToInt32(Math.Floor((TJETZT + 10) / 3600d));
                                        AMIND = (TJETZT / 3600 - (float)istu_var) * 60;
                                        IMIND = Convert.ToInt32(Math.Floor(AMIND));
                                        ASECD = (AMIND - (float)IMIND) * 60;
                                        ISECD = Convert.ToInt32(Math.Floor(ASECD));
                                        RadiationModel.RADIATRAD(LGEOM, false, false, iday_var, imonth_var, IJAHR, TJETZT, 3, 3, 3);
                                        Console.Write(iday_var.ToString() + "." + imonth_var.ToString() + "  -  " + istu_var.ToString() + ":00        ");
                                        Console.SetCursorPosition(0, Console.CursorTop);

                                        if (Windspeed_meteopgt <= 2)
                                        {
                                            if ((GLOBRAD[2][2] > 100) && (GLOBRAD[2][2] <= 200))
                                            {
                                                ISTU = istu_var;
                                                ITAG = iday_var;
                                                IMON = imonth_var;
                                                iday_var = 31;
                                                exit = true;
                                                break;
                                            }
                                        }
                                        else if ((Windspeed_meteopgt > 2) && (Windspeed_meteopgt <= 3))
                                        {
                                            if ((GLOBRAD[2][2] > 150) && (GLOBRAD[2][2] <= 650))
                                            {
                                                ISTU = istu_var;
                                                ITAG = iday_var;
                                                IMON = imonth_var;
                                                iday_var = 31;
                                                exit = true;
                                                break;
                                            }
                                        }
                                        else if ((Windspeed_meteopgt > 3) && (Windspeed_meteopgt <= 5))
                                        {
                                            if ((GLOBRAD[2][2] > 180) && (GLOBRAD[2][2] <= 650))
                                            {
                                                ISTU = istu_var;
                                                ITAG = iday_var;
                                                IMON = imonth_var;
                                                iday_var = 31;
                                                exit = true;
                                                break;
                                            }
                                        }
                                        else if (Windspeed_meteopgt > 5)
                                        {
                                            if (GLOBRAD[2][2] > 900)
                                            {
                                                ISTU = istu_var;
                                                ITAG = iday_var;
                                                IMON = imonth_var;
                                                iday_var = 31;
                                                exit = true;
                                                break;
                                            }
                                        }
                                    }
                                }
                                if (exit)
                                    break;
                            }
                        }
                    } // option dynamic sunrise

                    if (AKLA == 4)
                    {
                        ITAG = 21;
                        IMON = 3;
                        IJAHR = 2006;
                        ISTU = 12;
                        IMIN = 0;
                        Program.ISOL = 2; // 04/07/2018 Kuntner: thick clouds + lower solar radiation
                        Solar_reduction_factor = 0.1F; // try with reduced solar radiation

                        int month_setting = set_month_type(Windspeed_meteopgt, AKLA, BGRAD);
                        bool exit = false;

                        //find the month and time of the year, where the radiation equals 100 W/m? as required by SC 4
                        foreach (int imonth_var in Month_List[month_setting])
                        {
                            for (int iday_var = 1; iday_var < 29; iday_var += 2)
                            {
                                for (int istu_var = 1; istu_var < 16; istu_var += 2)
                                {
                                    TJETZT = (float)istu_var * 3600 + (float)IMIN * 60;
                                    ISTUD = Convert.ToInt32(Math.Floor((TJETZT + 10) / 3600d));
                                    AMIND = (TJETZT / 3600 - (float)istu_var) * 60;
                                    IMIND = Convert.ToInt32(Math.Floor(AMIND));
                                    ASECD = (AMIND - (float)IMIND) * 60;
                                    ISECD = Convert.ToInt32(Math.Floor(ASECD));
                                    RadiationModel.RADIATRAD(LGEOM, false, false, iday_var, imonth_var, IJAHR, TJETZT, NX_Radiat_Guess, NY_Radiat_Guess, 3);
                                    Console.Write(iday_var.ToString() + "." + imonth_var.ToString() + "  -  " + istu_var.ToString() + ":00        ");
                                    Console.SetCursorPosition(0, Console.CursorTop);

                                    float rsolg_mean = Mean_Array_Value(GLOBRAD, NX_Radiat_Guess, NY_Radiat_Guess);

                                    if ((rsolg_mean > 45) && (rsolg_mean <= 150))
                                    {
                                        ISTU = istu_var;
                                        ITAG = iday_var;
                                        IMON = imonth_var;
                                        iday_var = 31;
                                        istu_var = 16;
                                        exit = true;
                                        break;
                                    }
                                }
                                if (Solar_reduction_factor < 0.95F) // found no situation -> increase solar_reduction_factor
                                {
                                    Solar_reduction_factor += 0.05F;
                                }
                            }
                            if (exit)
                                break;
                        }
                    }
                    else if (AKLA == 5)
                    {
                        ITAG = 21;
                        IMON = 12;
                        IJAHR = 2006;
                        ISTU = 18;
                        IMIN = 0;
                    }
                    else if (AKLA == 6)
                    {
                        ITAG = 21;
                        IMON = 12;
                        IJAHR = 2006;
                        ISTU = 18;
                        IMIN = 0;
                    }
                    else if (AKLA == 7)
                    {
                        ITAG = 21;
                        IMON = 12;
                        IJAHR = 2006;
                        ISTU = 18;
                        IMIN = 0;
                    }
                }
                TJETZT = (float)ISTU * 3600 + (float)IMIN * 60;
                Console.WriteLine("  ** STR **                                ");
                ISTUD = Convert.ToInt32(Math.Floor((TJETZT + 10) / 3600d));
                AMIND = (TJETZT / 3600 - (float)ISTUD) * 60;
                IMIND = Convert.ToInt32(Math.Floor(AMIND));
                ASECD = (AMIND - (float)IMIND) * 60;
                ISECD = Convert.ToInt32(Math.Floor(ASECD));
                RadiationModel.RADIATRAD(LGEOM, LGEOMW, LGEOMR, ITAG, IMON, IJAHR, TJETZT, NX, NY, NZ);
            }
        }

        public static float Mean_Array_Value(double[][] Array, int NX, int NY)
        {
            float mean = 0;
            int count = 0;
            int startx = nr_cell_smooth;
            if (startx > (NX - 5))
            {
                startx = Math.Max(1, (NX - 5));
            }
            int starty = nr_cell_smooth;
            if (starty > (NY - 5))
            {
                starty = Math.Max(1, (NY - 5));
            }

            for (int i = startx; i < NX; i++)
            {
                for (int j = starty; j < NY; j++)
                {
                    count++;
                    mean += (float)Array[i][j];
                }
            }
            return (float)(mean / count);
        }

    }
}