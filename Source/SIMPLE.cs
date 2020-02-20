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
using System.Collections.Concurrent;
//using System.Diagnostics;

namespace GRAMM_CSharp_Test
{
    partial class Program
    {
        public static bool SOLUTION(int NI, int NJ, int NK)
        {
            //Check for numerical problems using overall massdivergence
            if (Program.MASSOURCE[Program.IDIV] > 1000) // 13.4.2017 Ku Check if Divergence increased to a very high value 
            	Program.Divergence_Min = Math.Min(Program.Divergence_Min, Program.MASSOURCE[Program.IDIV]);
           
            if((Program.MASSOURCE[Program.IDIV]*0.001 >= 500000) || (Program.MASSOURCE[Program.IDIV] > Program.Divergence_Min * 50))
            {
            	if (Program.computation_retry >= 3) // try 4 times
            	{
            		try
            		{
                        string err = "Numerical instabilities detected for flow field: " + Program.IWETTER.ToString();
                        ProgramWriters.LogfileProblemreportWrite(err);
            		}
            		catch { }
            		
            		Program.REALTIME = Program.TLIMIT; // exit this situation
            	}
                for (int i = 0; i < 11;i++ )
                    Program.MASSOURCE[i] = 0;
                return false; // error
            }
  
            //Calculation of several terms for the implicit scheme according to Patankar, 1980
            //Note that the wind speeds and the non-hydrostatic pressure are calculated using half-cells,
            //while pot. temperature, humidity, turbulent-kinetic energy and dissipation are modelled based on the full-cell
            if (Program.REALTIME < Program.DTI)
            {            	
  				TERMIVterms(NI, NJ, NK);
            	
                TERMIPSterms(NI, NJ, NK);
            					
        		TERMIPterms(NI, NJ, NK); 

                //Prandtl-layer quantities
                if(Program.ICPR) Prandtl_calculate(NI, NJ, NK);

                //Soil temperature (as the last 5 cells near the borders are not calculated, the routine just makes sense if model domains have more than 10 grid cells in the horizontal directions
                if ((NI > 10) && (NJ > 10) && (Program.ICTB)) Calctb_calculate(NI, NJ, NK);

                //Solve momentum equations
				if (Program.pOptions.MaxDegreeOfParallelism > 5) // more cores - divide the work -> try to avoid false sharing             	
                {
//                	int cores = Program.pOptions.MaxDegreeOfParallelism;
//                	Program.pOptions.MaxDegreeOfParallelism = Convert.ToInt32(Math.Ceiling(cores / 3.0f));
                	
                	Parallel.Invoke(Program.pOptions,
                    () => UIMPcalculate(NI, NJ, NK),
                    () => VIMPcalculate(NI, NJ, NK),
                    () => WVELcalculate(NI, NJ, NK));

//                	Program.pOptions.MaxDegreeOfParallelism = cores; // reset core-count
                }
                else
                {
                	UIMPcalculate(NI, NJ, NK);
                	VIMPcalculate(NI, NJ, NK);
                	WVELcalculate(NI, NJ, NK);
                }
                 
                //Boundary conditions
                Bords_calculate(NI, NJ, NK);

                //Solve the non-hydrostatic pressure equation
                if (Program.ICPN == true) CALCPR_calculate(NI, NJ, NK);

                //save initial mass-divergence to check for numerical instabilities at the end of the simulation
                if (Program.REALTIME <= (1.5 * Program.DT))
                {
                    Program.MASSOURCE_FIRST = Program.SUMG;
                }

                
                if (Program.ICT == true && Program.ICQU == true)
                {
                	Parallel.Invoke(Program.pOptions,
                    () => Timp_calculate(NI, NJ, NK),
                    () => Fimp_calculate(NI, NJ, NK));
                    //computation of water vapour content
                    if (Program.ICQU == true && Program.ISTAT != 0) WatVap_calculate(NI, NJ, NK);
                }
                else
                {
                	//computation of pot. temperature
                	if (Program.ICT == true) Timp_calculate(NI, NJ, NK);
					//computation of humidity
                	if (Program.ICQU == true) Fimp_calculate(NI, NJ, NK);
                    //computation of water vapour content
                    if (Program.ICQU == true && Program.ISTAT != 0) WatVap_calculate(NI, NJ, NK);
                }

                //computation of thermal pressure
                //if (Program.ICPSI == true) THERMPRESScalculate(NI, NJ, NK);
                
                //computation of turbulent kinetic energy and dissipation rate
                if (Program.ICTE == true)
                {
                	Parallel.Invoke(Program.pOptions,
                    () => Tkeimp_calculate(NI, NJ, NK),
                    () => Epsimp_calculate(NI, NJ, NK));
                }

                //advection/diffusion of passive scalars for the chmemical reaction mechanism + call chemistry subroutine
                if (Program.ISTAT == 2 && chemistry == true && REALTIME > Update_Chemistry_Threshold)
                {
                    Update_Chemistry_Threshold += Update_Chemistry;
                    Parallel.For(0, Program.NSPEZ, Program.pOptions, n =>
                    {
                        PSimp_calculate(NI, NJ, NK, n);
                    });
                    //@Johannes: insert here the call for the chemical reaction mechanism 
                }


                //time-step adjustion
                if (Program.ISTAT == 0)
                {
                    Program.IDIV_LockDown2--;
                    if (Program.IDIV_LockDown > 0) // 5.4.2017 Ku avoid increase of time step after a decrease
                    {
                        Program.IDIV_LockDown--;
                    }
                    else
                    {
                        Program.STEIGUNG = Math.Round(Program.STEIGUNG, 2);
                        if ((Program.IDIV == 10) && (Program.STEIGUNG <= 0))
                            if (REALTIME <= DTI * 0.75) // do not adjust time step if steady state criterion is measured
                            {
                                Program.DT = Math.Min(Program.DT + 0.5, Program.DTMAX); //3.4.2017 Ku
                                Program.IDIV_LockUp = 0;
                                Program.IDIV_LockDown2 = 40;
                            }
                    }

                    //5.4.2017 Ku
                    if (((Program.MASSOURCE_Old * 1.01 < Program.MASSOURCE_Act) || ((Program.MASSOURCE_Old + 500) < Program.MASSOURCE_Act)) && IDIV_LockUp <= 3 && IDIV_LockDown == 0 && IDIV_LockDown2 < 0)
                    {
                        Program.DT = Math.Max(Program.DT * 0.80, 1.5); //5.4.2017 Ku
                        RELAXV = Math.Max(RELAXV * 0.80, 0.05);
                        RELAXT = Math.Max(RELAXT * 0.80, 0.05);
                        Program.IDIV_LockDown = 11;
                        Program.IDIV_LockDown2 = 20;
                        Program.IDIV_LockUp++;
                        Program.IDIV_PingPong++;

                        if (Program.IDIV_PingPong > 6 * Math.Pow(Program.Relaxv_ori / Program.RELAXV, 2) && Program.IDIV_LockRelax < 6)
                        {   // Set relax factors down and increase lockDown value if downsteps occur periodically
                            Program.DT = Math.Max(Program.DT * 0.80, 1.5); //5.4.2017 Ku
                            RELAXV = Math.Max(RELAXV * 0.80, 0.05);
                            RELAXT = Math.Max(RELAXT * 0.80, 0.05);
                            Program.IDIV_LockRelax++;
                            Program.IDIV_PingPong = 0;
                            Program.IDIV_LockDown = 20;
                        }
                        else if (IDIV_LockUp > 2 && IDIV_LockRelax < 6)
                        {
                            RELAXV = Math.Max(RELAXV * 0.80, 0.05);
                            RELAXT = Math.Max(RELAXT * 0.80, 0.05);
                            Program.IDIV_LockRelax++;
                        }
                    }
                }
                //5.4.2017 Ku
                

                //2.9.2019 Öt
                //dynamic time step and relaxation factor adjustments for GRAMM full transient simulations
                if (Program.ISTAT != 0)
                {
                    //Console.WriteLine(Program.IDIV.ToString() + "," + Program.MASSOURCE[Program.IDIV].ToString() + "," + Program.MASSOURCE[1].ToString() + "," + Program.MASSOURCE_FIRST.ToString());
                    /*
                    if ((Program.IDIV == 5) && (Program.MASSOURCE[Program.IDIV] > Program.MASSOURCE[1]) && (Program.MASSOURCE_FIRST < Program.MASSOURCE[Program.IDIV]))
                    {
                        Program.DT = Math.Max(Program.DT * 0.70, 1.5); //5.4.2017 Ku
                        RELAXV = Math.Max(RELAXV * 0.70, 0.05);
                        RELAXT = Math.Max(RELAXT * 0.70, 0.05);
                        Program.IDIV = 0;
                    }
                    else if ((Program.IDIV == 5) && (Program.MASSOURCE[Program.IDIV] <= Program.MASSOURCE[1]))
                    {

                        Program.DT = Math.Min(Program.DT + 0.5, Program.DTMAX);
                        RELAXV = Math.Min(RELAXV * 1.01, 0.15);
                        RELAXT = Math.Min(RELAXT * 1.01, 0.15);
                        Program.IDIV = 0;
                    }
                    */
                    Program.DT = Program.DTMAX;
                }
                //2.9.2019 Öt

                Program.TerminalOut++;
                if (Program.TerminalOut >= TerminalThreshold) // if Counter is even
                {
                	//output to screen
                	string LOGT = "-";
                	string LOGH = "-";
                	string LOGU = "-";
                	string LOGV = "-";
                	string LOGW = "-";
                	if (Program.ICT) LOGT = "+";
                	if (Program.ICQU) LOGH = "+";
                	if (Program.ICU) LOGU = "+";
                	if (Program.ICV) LOGV = "+";
                	if (Program.ICW) LOGW = "+";
                	Console.WriteLine("-------------------------------------------------------------------------------------");
                	Console.WriteLine("WEATHER-SIT. TIME[s]  TIMESTEP[s]  ENDTIME[s]  PRESS-ITERATIONS  DIVERGENCE U V W T H");
                	Console.WriteLine(Program.IWETTER.ToString().PadLeft(10) + "/" + (1+computation_retry).ToString().PadLeft(1) + " " + Program.REALTIME.ToString("0.0").PadLeft(7) + "  " +
                	                  Program.DT.ToString("0.0").PadLeft(11) + "  " + Program.DTI.ToString("0").PadLeft(10) + "               " +
                	                  (Program.INUMS - 1).ToString("0").PadLeft(3) + "  " + (Program.SUMG * 0.001).ToString("0.000").PadLeft(10) + " " +
                	                  LOGU + " " + LOGV + " " + LOGW + " " + LOGT + " " + LOGH);
                }                
            }
            
			return true; // computation OK
        }

		public static void Relax_Border()
        {

            if (Program.Relax_Border_factor[0][0] < 0)
            {
                int NX = Program.NX;
                int NY = Program.NY;

                int BorderE = NX - Program.nr_cell_smooth;
                int BorderS = NY - Program.nr_cell_smooth;
                int Border_Expand = (int)(Program.nr_cell_smooth * 1.8);

                for (int i = 0; i <= NX; ++i)
                {
                    for (int j = 0; j <= NY; ++j)
                    {
                        double factor = 1;
                        if (i >= Border_Expand && j >= Border_Expand && i <= (NX - Border_Expand) && j <= (NY - Border_Expand))
                        {
                            factor = 1;
                        }
                        else if (i < Program.nr_cell_smooth)
                        {
                            factor = 1 - Math.Abs((AH[Program.nr_cell_smooth][j] - AH[1][j]) / (2 * Program.DDX[1] * Program.nr_cell_smooth));
                        }
                        else if (j < Program.nr_cell_smooth)
                        {
                            factor = 1 - Math.Abs((AH[i][Program.nr_cell_smooth] - AH[i][1]) / (2 * Program.DDY[1] * Program.nr_cell_smooth));
                        }
                        else if (i > BorderE)
                        {
                            factor = 1 - Math.Abs((AH[NX - Program.nr_cell_smooth][j] - AH[NX][j]) / (2 * Program.DDX[1] * Program.nr_cell_smooth));
                        }
                        else if (j > BorderS)
                        {
                            factor = 1 - Math.Abs((AH[i][NY - Program.nr_cell_smooth] - AH[i][NY]) / (2 * Program.DDY[1] * Program.nr_cell_smooth));
                        }

                        // check transition area
                        if (i < Border_Expand)
                        {
                            factor *= Math.Max(0.5, (float)i / Border_Expand);
                        }
                        else if (j < Border_Expand)
                        {
                            factor *= Math.Max(0.5, (float)j / Border_Expand);
                        }
                        else if (i > NX - Border_Expand)
                        {
                            factor *= Math.Max(0.5, (float)(NX - i) / Border_Expand);
                        }
                        else if (j > NY - Border_Expand)
                        {
                            factor *= Math.Max(0.5, (float)(NY - j) / Border_Expand);
                        }

                        factor = Math.Min(1F, Math.Max(factor, 0.05F));
                        Program.Relax_Border_factor[i][j] = (float)factor;
                    }
                }
            }
        }

    }
}
