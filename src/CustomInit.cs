#region Copyright
///<remarks>
/// <GRAMM Mesoscale Model>
/// Copyright (C) [2022]  [Markus Kuntner]
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
    /// <summary>
    /// All user defined custom air temperature and soil temperature values
    /// </summary>
    public struct CustomAirSoilInit
    {
        /// <summary>
        /// Soil is set to snow cover above this absolute height 
        /// </summary>
        public readonly double SnowHeightThreshold;
        /// <summary>
        /// Initial air temperature at sea-level in 2m height above ground
        /// </summary>
        public readonly double TAir2m;
        /// <summary>
        /// Initial soil temperature at sea-level
        /// </summary>
        public readonly double TSurface;
        /// <summary>
        /// Initial soil temperature at sea-level in 1 m depth
        /// </summary>
        public readonly double TSurface1m;
        /// <summary>
        /// Delta T to surface for water bodies
        /// </summary>
        public readonly double DeltaTSurfaceWater;
        /// <summary>
        /// Relative air humidity 
        /// </summary>
        public readonly double RelHumidity;
        /// <summary>
        /// Inversion height above lowest terrain cell within the GRAMM domain 
        /// </summary>
        public readonly double InversionHeight;
        /// <summary>
        /// Flag, if the air temperature has been set by the user
        /// </summary>
        public readonly bool UserdefinedAirTemp;
        /// <summary>
        /// Gradient of the air temperature - default -0.0065 K
        /// </summary>
        public readonly double AirTempGradient;
        /// <summary>
        /// Gradient of the air temperature below the inversion height - default +0.01 K
        /// </summary>
        public readonly double AirTempGradientBelowInversionHeight;
        /// <summary>
        /// Gradient of the soil temperature - default -0.005 K
        /// </summary>
        public readonly double SoilTempGradient;

        /// <summary>
        /// Init Custom Init Struct
        /// </summary>
        public CustomAirSoilInit(double Tinit, double TBinit, double TBinit1, double Qinit)
        {
            //set default parameters
            SnowHeightThreshold = 1000000;
            TAir2m = Tinit;
            TSurface = TBinit;
            TSurface1m = TBinit1;
            RelHumidity = Qinit;
            InversionHeight = 10000;
            UserdefinedAirTemp = false;
            AirTempGradient = -0.0065;
            AirTempGradientBelowInversionHeight = 0.01;
            SoilTempGradient = -0.005;
            DeltaTSurfaceWater = 0;

            // Read the custom init file and set snow to 1 if the height is exceeded<br></br>
            // <b>Set user defined custom init values
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
                                double.TryParse(text[4], NumberStyles.AllowDecimalPoint, CultureInfo.InvariantCulture, out double temp);
                                if (temp > 0)
                                {
                                    SnowHeightThreshold = temp;
                                }
                            }
                            if (text.Length > 5)
                            {
                                double.TryParse(text[5], NumberStyles.AllowDecimalPoint, CultureInfo.InvariantCulture, out double temp);
                                if (temp > 0)
                                {
                                    TAir2m = temp;
                                    UserdefinedAirTemp = true;
                                }
                            }
                            if (text.Length > 6)
                            {
                                double.TryParse(text[6], NumberStyles.AllowDecimalPoint, CultureInfo.InvariantCulture, out double temp);
                                if (temp > 0)
                                {
                                    TSurface = temp;
                                }
                            }
                            if (text.Length > 7)
                            {
                                double.TryParse(text[7], NumberStyles.AllowDecimalPoint, CultureInfo.InvariantCulture, out double temp);
                                if (temp > 0)
                                {
                                    TSurface1m = temp;
                                }
                            }
                            if (text.Length > 8)
                            {
                                double.TryParse(text[8], NumberStyles.AllowDecimalPoint | NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out double temp);
                                if (temp != 0)
                                {
                                    DeltaTSurfaceWater = temp;
                                }
                            }
                            if (text.Length > 9)
                            {
                                double.TryParse(text[9], NumberStyles.AllowDecimalPoint, CultureInfo.InvariantCulture, out double temp);
                                if (temp > 0)
                                {
                                    RelHumidity = temp;
                                }
                            }
                            if (text.Length > 10)
                            {
                                double.TryParse(text[10], NumberStyles.AllowDecimalPoint, CultureInfo.InvariantCulture, out double temp);
                                if (temp > 0)
                                {
                                    InversionHeight = temp;
                                }
                            }
                            if (text.Length > 11)
                            {
                                double.TryParse(text[11], NumberStyles.AllowDecimalPoint | NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out double temp);
                                if (temp != 0)
                                {
                                    AirTempGradient = temp;
                                }
                            }
                            if (text.Length > 12)
                            {
                                double.TryParse(text[12], NumberStyles.AllowDecimalPoint | NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out double temp);
                                if (temp != 0)
                                {
                                    AirTempGradientBelowInversionHeight = temp;
                                }
                            }
                            if (text.Length > 13)
                            {
                                double.TryParse(text[13], NumberStyles.AllowDecimalPoint | NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out double temp);
                                if (temp != 0)
                                {
                                    SoilTempGradient = temp;
                                }
                            }
                        }
                    }
                }
                catch { }
            }
        }

    }
}