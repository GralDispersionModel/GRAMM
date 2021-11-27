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

namespace GRAMM_2001
{
    public class ProgramWriters
    {
        private static int month_old = 0;
        private static int day_old = 0;
        private static int hour_old = 0;
        private static int minute_old = 0;


        /// <summary>
        ///Output of GRAMM logging file
        /// </summary>
        public static void LogfileGrammCoreInfo()
        {
            // Write additional Data to the Log-File
            LogfileGrammCoreWrite("");

            string err = "GRAMM domain area";
            LogfileGrammCoreWrite(err);
            err = "  West  (abs): " + Program.GRAMM_West.ToString();
            LogfileGrammCoreWrite(err);
            err = "  East  (abs): " + Convert.ToString((Int32)(Program.GRAMM_West + Program.NX * Program.DDXImm[1]));
            LogfileGrammCoreWrite(err);
            err = "  North (abs): " + Convert.ToString((Int32)(Program.GRAMM_South + Program.NY * Program.DDYImm[1]));
            LogfileGrammCoreWrite(err);
            err = "  South (abs): " + Program.GRAMM_South.ToString();
            LogfileGrammCoreWrite(err);
            err = "  Latitude [°]: " + Program.BGRAD.ToString();
            LogfileGrammCoreWrite(err);

        } // Write Logfile		

        /// <summary>
        /// Output to GRAMM log file
        /// </summary>
        /// <param name="a"></param>
        public static void LogfileGrammCoreWrite(string a)
        {
            try
            {
                using (StreamWriter sw = new StreamWriter("Logfile_GRAMMCore.txt", true))
                {
                    sw.WriteLine(a);
                    sw.Flush();
                }
            }
            catch { }
        }

        /// <summary>
        ///Output of GRAMM Problem report
        /// </summary>
        public static void LogfileProblemreportWrite(string a)
        {
            a = "GRAMM Error: " + a;
            try
            {
                using (StreamWriter sw = new StreamWriter("Problemreport_GRAMM.txt", true))
                {
                    sw.WriteLine(a);
                    sw.Flush();
                }
            }
            catch { }
            LogfileGrammCoreWrite(a); // Write error to LogfileCore
        }

    }
}
