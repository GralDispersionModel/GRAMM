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
using System.Data;
using System.Globalization;
#if __ECMWF__
using OSGeo.OSR;
#endif

namespace GRAMM_CSharp_Test
{
    /*  GRAZ MESOSCALE MODELL GRAMM
        COMPREHENSIVE DESCRIPTION CAN BE FOUND IN OETTL, 2016
        THE GRAMM MODEL HAS BEEN ORIGINALLY DEVELOPED BY RAIMUND ALMBAUER AROUND 1990.
        THE RADIATION MODEL IS BASED ON THE MODEL SUGGESTED BY SOMIESKI 1988 AND HAS BEEN IMPLEMENTED AROUND 1990 BY ROBERT ALMER
        THE RE-INITIALIZATION PROCEDURE HAS BEEN IMPLEMENTED BY ULRICH UHRNER IN 2013
        THE FOLLOWING MODEL PARTS HAVE BEEN IMPLEMENTED BY DIETMAR OETTL SINCE 1997:
        -DIAGNOSTIC MODEL FOR INITIALIZATION OF WIND AND TEMPERATURE FIELDS USING POINT AND PROFILE OBSERVATIONS
        -FULL IMPLICIT SOLVERS OF ALL CONSERVATION EQUATIONS USING TDMA ALGORITHM TO IMPROVE NUMERICAL STABILITIES IN HIGHLY COMPLEX TERRAIN
        -STANDARD K-EPS MODEL
        -PARALLELISATION
        -PORTING THE CODE FROM FORTRAN TO C#
	 */

    partial class Program
    {
        /*INPUT FILES :
          BASIC DOMAIN INFORMATION        GRAMM.geb
          MAIN CONTROL PARAM FILEs        iin.dat  or IIN.dat
          GEOMETRY DATA                   ggeom.asc
          LANDUSE DATA                    landuse.asc
          METEOROLOGICAL DATA             meteopgt.all
          MAX. NUMBER OF CPUs             Max_Proc.txt
		 */

        //global variables declaration block
        //note that arrays are defined as jagged-arrays as they have a much better computational performance than matrices

        public static ParallelOptions pOptions = new ParallelOptions();                //sets the maximum number of cores to be used in the parallelisation
        public static string decsep;                                                   //global decimal separator of the system
        public static Int32 NX;                                                        //number of cells in x-direction
        public static Int32 NY;                                                        //number of cells in y-direction
        public static Int32 NZ;                                                        //number of cells in z-direction
        public static Int32 NX1;                                                       //number of cells in x-direction +1
        public static Int32 NY1;                                                       //number of cells in y-direction +1
        public static Int32 NZ1;                                                       //number of cells in z-direction +1
        public static Int32 NX2;                                                       //number of cells in x-direction +2
        public static Int32 NY2;                                                       //number of cells in y-direction +2
        public static Int32 NZ2;                                                       //number of cells in z-direction +2
        public static Int32 NZB;                                                       //NUMBER OF VERTICAL GRID CELLS OF THE SOIL MODEL
        public static Int32 NPROFMAX;                                                  //NUMBER OF VERTICAL GRID CELLS OF THE RADIATION MODEL
        public static Int16 IMETSTR;                                                   //Flag needed that geometry data for the radiation model is read just once
        public static float[][] AH = CreateArray<float[]>(1, () => new float[1]);                                               //Height of the surface
        public static float[][] AH_Bassins = CreateArray<float[]>(1, () => new float[1]);                                       //Height of the bassins and valleys 
        public static float[][] TPI = CreateArray<float[]>(1, () => new float[1]);                                       //Height of the bassins and valleys 
        public static Int16 IHOURO;                                                    //Flag needed for the intermediate output of GRAMM flowfields
        public static float[][][] VOL = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));           //Volume of grid cells
        public static float[][][] AREAX = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));         //Area of the grid cell in x-direction
        public static float[][][] AREAY = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));         //Area of the grid cell in y-direction
        public static float[][][] AREAZ = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));         //Bottom area of the grid cell
        public static float[][][] AREA = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));          //Projection of the ground area of the grid cell in z-direction
        public static float[][][] AREAXYZ = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));       //Area between the two halfs of the grid cell
        public static float[][][] AREAZX = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));        //Projection of the ground area of the grid cell in x-direction
        public static float[][][] AREAZY = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));        //Projection of the ground area of the grid cell in y-direction
        public static float[][][] AHE = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));           //Heights of the corner points of each grid cell
        public static float[][][] ZSP = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));           //Height of the centre point of each grid cell
        public static float[] DDX = new float[1];                                                                               //Horizontal grid size in x-direction
        public static float[] DDY = new float[1];                                                                               //Horizontal grid size in y-direction
        public static float[] ZAX = new float[1];                                                                               //Distance between neighbouring grid cells in x-direction
        public static float[] ZAY = new float[1];                                                                               //Distance between neighbouring grid cells in y-direction
        public static float[] X = new float[1];                                                                                 //Distance of grid cell centre from western domain border
        public static float[] Y = new float[1];                                                                                 //Distance of grid cell centre from southern domain border
        public static float[] Z = new float[1];                                                                                 //Temporary array for generating the terrain-following grid
                                                                                                                                //public static double[][] LAND = CreateArray<double[]>(1, () => new double[1]);                                            //CORINE land-use categories
        public static double[] PBZZ = new double[1];                                                                              //Base state for pressure used in the radiation model
        public static double[] RHOBZZ = new double[1];                                                                            //Base state for density used in the radiation model
        public static double[][][] TABS = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));        //Absolute temperature in K
        public static double[][] RSOLG = CreateArray<double[]>(1, () => new double[1]);                                           //Solar incoming radiation without reflected radiation
        public static double[][] GLOBRAD = CreateArray<double[]>(1, () => new double[1]);                                         //Globar incoming radiation
        public static double[][][] DT_SOL = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));      //Solar heating rate
        public static double[][] RL = CreateArray<double[]>(1, () => new double[1]);                                              //Longwave outgoing radiation
        public static double[][][] DT_TERR = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));     //Terrestrial heating rate
        public static double[][] TG = CreateArray<double[]>(1, () => new double[1]);                                              //Absolute soil temperature at the second level beneath the surface
        public static double[][] EPSG = CreateArray<double[]>(1, () => new double[1]);                                            //Surface emissivity
        public static int[][] KST = CreateArray<int[]>(1, () => new int[1]);                                                   //Defines the cell in vertical direction from which onward the radiation model is computed
        public static float[] ZZ = new float[1];                                                                                  //Average value of two neighbouring cell heights
        public static double[][] ALBEDO = CreateArray<double[]>(1, () => new double[1]);                                          //Surface albedo
        public static double[][] CLOUDS = CreateArray<double[]>(1, () => new double[1]);                                          //Cloudyness 0 < x < 1
        public static double[][] SNOW = CreateArray<double[]>(1, () => new double[1]);                                            //Snow cover 0/1
        public static double[] ZPROF = new double[1];                                                                             //Cell heights for the radiation model
        public static double[] TPROF = new double[1];                                                                             //Vertical temperature profile in the radiation model
        public static double[] PPROF = new double[1];                                                                             //Vertical pressure profile in the radiation model
        public static double[] VNORM = new double[1];                                                                             //standard optical length
        public static double[] QVAP = new double[1];                                                                              //water content in the atmosphere used in radiation model
        public static double[] QCLD = new double[1];                                                                              //water content in clouds used in radiation model
        public static double[] QRAIN = new double[1];                                                                             //water content in rain used in radiation model
        public static double[] QICE = new double[1];                                                                              //water content in ice used in radiation model
        public static double[] QSNOW = new double[1];                                                                             //water content in snow used in radiation model
        public static double[][][] U = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));           //Average wind speed in x-direction over the two half-cells at time t
        public static double[][][] V = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));           //Average wind speed in y-direction over the two half-cells at time t
        public static double[][][] W = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));           //Average wind speed in z-direction over the two half-cells at time t
        public static float[][][] RHO = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));             //Density of the air
        public static double[][][] UN = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Average wind speed in x-direction over the two half-cells at time t+1
        public static double[][][] VN = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Average wind speed in y-direction over the two half-cells at time t+1
        public static double[][][] WN = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Average wind speed in z-direction over the two half-cells at time t+1
        public static double[][][] U1 = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Average wind speed in x-direction in the lower half-cells at time t
        public static double[][][] V1 = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Average wind speed in y-direction in the lower half-cells at time t
        public static double[][][] W1 = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Average wind speed in z-direction in the lower half-cells at time t
        public static double[][][] U2 = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Average wind speed in x-direction in the upper half-cells at time t
        public static double[][][] V2 = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Average wind speed in y-direction in the upper half-cells at time t
        public static double[][][] W2 = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Average wind speed in z-direction in the upper half-cells at time t
        public static double[][][] U1N = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));         //Average wind speed in x-direction in the lower half-cells at time t+1
        public static double[][][] V1N = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));         //Average wind speed in y-direction in the lower half-cells at time t+1
        public static double[][][] W1N = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));         //Average wind speed in z-direction in the lower half-cells at time t+1
        public static double[][][] U2N = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));         //Average wind speed in x-direction in the upper half-cells at time t+1
        public static double[][][] V2N = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));         //Average wind speed in y-direction in the upper half-cells at time t+1
        public static double[][][] W2N = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));         //Average wind speed in z-direction in the upper half-cells at time t+1
        public static double[][][] PN = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Non-hydrostatic pressure in the cell-centre
        public static float[][] Relax_Border_factor = CreateArray<float[]>(1, () => new float[1]);                                //Relaxation reduction factor at the border cells

        public static double[][][] DP = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Non-hydrostatic pressure in the cell-centre
        public static double[][][] DPX = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));         //Non-hydrostatic pressure in x-direction
        public static double[][][] DPY = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));         //Non-hydrostatic pressure in y-direction
        public static double[][][] DPZ = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));         //Non-hydrostatic pressure in z-direction
        public static double[][][] DDP1DX = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));      //Pressure gradient used for the velocity correction in x-direction
        public static double[][][] DDP1DY = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));      //Pressure gradient used for the velocity correction in y-direction
        public static double[][][] DDP1DZ = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));      //Pressure gradient used for the velocity correction in z-direction
        public static double[][][] DDP2DX = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));      //Pressure gradient used for the velocity correction in x-direction
        public static double[][][] DDP2DY = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));      //Pressure gradient used for the velocity correction in y-direction
        public static double[][][] DDP2DZ = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));      //Pressure gradient used for the velocity correction in z-direction          //Non-hydrostatic pressure in the cell-centre
        public static double[][][] TPDX = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));        //thermic pressure gradient in x-direction
        public static double[][][] TPDY = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));        //thermic pressure gradient in y-direction
        public static double[][][] TP = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //thermal pressure
        public static double[][][] TPX = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));         //thermal pressure in x-direction
        public static double[][][] TPY = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));         //thermal pressure in y-direction
        public static double[][][] TPZ = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));         //thermal pressure in z-direction
        public static double POBEN;                                                                                               //pressure at the top of the model domain
        public static float[][][] SUX = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));             //Mass divergence in x-direction
        public static float[][][] SUY = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));             //Mass divergence in y-direction
        public static float[][][] SUZ = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));             //Mass divergence in z-direction
        public static float[][][] SUXYZ = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));           //Mass divergence in between the two half-cells
        public static float[][][] U1NRHO = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));          //Mass-flux in x-direction
        public static float[][][] U2NRHO = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));          //Mass-flux in y-direction
        public static float[][][] V1NRHO = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));          //Mass-flux in z-direction
        public static float[][][] V2NRHO = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));          //Mass-flux in x-direction
        public static float[][][] W1NRHO = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));          //Mass-flux in y-direction
        public static float[][][] W2NRHO = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));          //Mass-flux in z-direction
        public static float[][][] PBZ = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));         //Hydrostatic pressure
        public static float[][][] RHOBZ = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));           //Air densitiy at the simulation start
        public static double[][][] PNBZKP = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));      //Temporary field in the calculation of the radiation fluxes
        public static double[][][] NBZKP = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));       //Temporary field in the calculation of the radiation fluxes
        public static float[][][] UG1 = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));             //Geostrophic wind in x-direction at time t-1
        public static float[][][] VG1 = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));             //Geostrophic wind in y-direction at time t-1
        public static float[][][] UG2 = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));             //Geostrophic wind in x-direction at time t+1
        public static float[][][] VG2 = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));             //Geostrophic wind in y-direction at time t+1
        public static double[][][] UG = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Geostrophic wind in x-direction at time t
        public static double[][][] VG = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Geostrophic wind in y-direction at time t
        public static double[][][] VISV = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));        //Turbulent viscosity in z-direction
        public static double[][][] VISH = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));        //Turbulent viscosity in horizontal-directions
        public static double[][][] QU = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Specific humidity at time t
        public static double[][][] QUN = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));         //Specific humidity at time t+1
        public static double[][][] WAT_VAP = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));     //Content of cloud water in the atmosphere in g/kg at time t
        public static double[][][] WAT_VAPN = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));    //Content of cloud water in the atmosphere in g/kg at time t+1
        public static double[][][] QBZ = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));         //Basic state of specific humidity at the simulation start
        public static double[][] QUG = CreateArray<double[]>(1, () => new double[1]);                                             //Soil moisture content
        public static double[][] FW = CreateArray<double[]>(1, () => new double[1]);                                              //Specific soil moisture parameter (e.g. water = 1)
        public static double[][][] T = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));           //Potential temperature of air at time t
        public static double[][][] TN = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Potential temperature of air at time t+1
        public static float[][][] FACTOR = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));          //Conversion factor between potential and absolute temperature of air
        public static double[][][] DISS = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));        //Turbulent dissipation rate at time t
        public static double[][][] DISSN = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));       //Turbulent dissipation rate at time t+1
        public static double[][][] TE = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Turbulent kinetic energy at time t
        public static double[][][] TEN = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));         //Turbulent kinetic energy at time t+1
        public static float[][] ZI = CreateArray<float[]>(1, () => new float[1]);                                                 //Mixing height
        public static double[][] UST = CreateArray<double[]>(1, () => new double[1]);                                             //Fricition velocity
        public static double[][] USTV = CreateArray<double[]>(1, () => new double[1]);                                            //Friction velocity normalized by the wind speed
        public static double[][] TST = CreateArray<double[]>(1, () => new double[1]);                                             //Potential temperature scale
        public static float[][] OL = CreateArray<float[]>(1, () => new float[1]);                                                 //Obukhov-length
        public static float[][] XWQ = CreateArray<float[]>(1, () => new float[1]);                                                //Temporary field used for calculating the sensible and latent heat fluxes
        public static float[][] Z0 = CreateArray<float[]>(1, () => new float[1]);                                                 //Aerodynamic surface roughness
        public static float[][][] FAC = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));             //Conversion factor between potential and absolute temperature of air
        public static double[][][] RITSCH = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));      //Gradient Richardson number
        public static double[][] WQU = CreateArray<double[]>(1, () => new double[1]);                                             //Sensible heat flux
        public static float[][] AWQ = CreateArray<float[]>(1, () => new float[1]);                                                //Anthropogenic heat flux
        public static double[][][] TB = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Absolute soil temperature at time t
        public static double[][][] TBN = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));         //Absolute soil temperature at time t+1
        public static double[][][] TBA = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));         //Absolute soil temperature at time t-1
        public static double[][][] TBZ = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));         //Basic state of potential temperature at the simulation start
        public static float[] DWB = new float[1];                                                                                 //Vertical cell heights of the soil model
        public static float[] DZB = new float[1];                                                                                 //Vertical distance between neighbouring cells of the soil model
        public static float[][] RHOB = CreateArray<float[]>(1, () => new float[1]);                                               //Soil density
        public static float[][] ALAMBDA = CreateArray<float[]>(1, () => new float[1]);                                            //Surface albedo
        public static double[] ALPHA = new double[1];                                                                             //Damping parameter used at the model boundaries
        public static double[][][] A_PS = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));        //Temporary field used in the TDMA solver for passive scalars
        public static double[][][] B_PS = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));        //Temporary field used in the TDMA solver for passive scalars
        public static double[][][] C_PS = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));        //Temporary field used in the TDMA solver for passive scalars
        public static float[][][] AWEST_PS = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));    //Temporary field used in the TDMA solver for passive scalars
        public static float[][][] ASOUTH_PS = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));   //Temporary field used in the TDMA solver for passive scalars
        public static float[][][] AEAST_PS = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));    //Temporary field used in the TDMA solver for passive scalars
        public static float[][][] ANORTH_PS = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));   //Temporary field used in the TDMA solver for passive scalars
        public static float[][][] AP0_PS = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));      //Temporary field used in the TDMA solver for passive scalars
                                                                                                                              //public static double[][][] SU_PS = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));       //Temporary field used in the TDMA solver for passive scalars
        public static float[][][] AIM = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));         //Temporary field used in the TDMA solver for velocities
        public static float[][][] BIM = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));         //Temporary field used in the TDMA solver for velocities
        public static float[][][] CIM = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));         //Temporary field used in the TDMA solver for velocities
        public static float[][][] AW1 = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));         //Temporary field used in the TDMA solver for velocities
        public static float[][][] AS1 = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));         //Temporary field used in the TDMA solver for velocities
        public static float[][][] AE2 = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));         //Temporary field used in the TDMA solver for velocities
        public static float[][][] AN2 = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));         //Temporary field used in the TDMA solver for velocities
                                                                                                                              //public static double[][][] SU = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));        //Temporary field used in the TDMA solver for velocities
        public static float[][][] AP0 = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));             //Temporary field used in the TDMA solver for velocities
        public static float[][][] F1U = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));             //Temporary field for UIMP.  VIMP, WIMP
        public static float[][][] F2U = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));             //Temporary field for UIMP.  VIMP, WIMP
        public static float[][][] F1V = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));             //Temporary field for UIMP.  VIMP, WIMP
        public static float[][][] F2V = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));             //Temporary field for UIMP.  VIMP, WIMP
        public static float[][][] F1W = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));             //Temporary field for UIMP.  VIMP, WIMP
        public static float[][][] F2W = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));             //Temporary field for UIMP.  VIMP, WIMP
        public static double[][][] DIMU = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));        //Temporary field used in the TDMA solver for velocities
        public static double[][][] DIMV = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));        //Temporary field used in the TDMA solver for velocities
        public static double[][][] DIMW = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));        //Temporary field used in the TDMA solver for velocities
        public static double[][][] AP = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Temporary field used in the TDMA solver for non-hydrostatic pressure
        public static double[][][] AB = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Temporary field used in the TDMA solver for non-hydrostatic pressure
        public static double[][][] AT = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Temporary field used in the TDMA solver for non-hydrostatic pressure
        public static double[][][] AW = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Temporary field used in the TDMA solver for non-hydrostatic pressure
        public static double[][][] AE = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Temporary field used in the TDMA solver for non-hydrostatic pressure
        public static double[][][] AS = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Temporary field used in the TDMA solver for non-hydrostatic pressure
        public static double[][][] AN = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Temporary field used in the TDMA solver for non-hydrostatic pressure
        public static float[][] VDATA1 = CreateArray<float[]>(1, () => new float[1]);                                            //Temporary field used for online output of prognostic variables
        public static float[][] VDATA2 = CreateArray<float[]>(1, () => new float[1]);                                            //Temporary field used for online output of prognostic variables
        public static float[][] VDATA3 = CreateArray<float[]>(1, () => new float[1]);                                            //Temporary field used for online output of prognostic variables
        public static float[][] VDATA4 = CreateArray<float[]>(1, () => new float[1]);                                            //Temporary field used for online output of prognostic variables
        public static float[][] VDATA8 = CreateArray<float[]>(1, () => new float[1]);                                            //Temporary field used for online output of prognostic variables
        public static float[] VDATA5 = new float[1];                                                                              //Temporary field used for online output of prognostic variables
        public static float[] VDATA6 = new float[1];                                                                              //Temporary field used for online output of prognostic variables
        public static float[] VDATA7 = new float[1];                                                                              //Temporary field used for online output of prognostic variables
        public static float[] VDATA9 = new float[1];                                                                              //Temporary field used for online output of prognostic variables
        public static double[] MASSOURCE = new double[1];                                                                         //Control variable for the temporal development of the entire mass divergence
        public static float[][] U_TEMP = CreateArray<float[]>(1, () => new float[1]);                                             //Control variable for checking flow convergence in the lowest layer
        public static float[][] V_TEMP = CreateArray<float[]>(1, () => new float[1]);                                             //Control variable for checking flow convergence in the lowest layer
        public static float[][] W_TEMP = CreateArray<float[]>(1, () => new float[1]);                                             //Control variable for checking flow convergence in the lowest layer
        public static double[][] USE1 = CreateArray<double[]>(1, () => new double[1]);                                            //Large scale value of u-wind at eastern boundary at time t
        public static double[][] VSE = CreateArray<double[]>(1, () => new double[1]);                                             //Large scale value of v-wind at eastern boundary at time t
        public static double[][] WSE = CreateArray<double[]>(1, () => new double[1]);                                             //Large scale value of w-wind at eastern boundary at time t
        public static double[][] TSE = CreateArray<double[]>(1, () => new double[1]);                                             //Large scale value of pot. temp. at eastern boundary at time t
        public static double[][] QUSE = CreateArray<double[]>(1, () => new double[1]);                                            //Large scale value of spec. humidity at eastern boundary at time t
        public static double[][] USEN = CreateArray<double[]>(1, () => new double[1]);                                            //Large scale value of u-wind at eastern boundary at time t+1
        public static double[][] VSEN = CreateArray<double[]>(1, () => new double[1]);                                            //Large scale value of v-wind at eastern boundary at time t+1
        public static double[][] WSEN = CreateArray<double[]>(1, () => new double[1]);                                            //Large scale value of w-wind at eastern boundary at time t+1
        public static double[][] TSEN = CreateArray<double[]>(1, () => new double[1]);                                            //Large scale value of pot. temp. at eastern boundary at time t+1
        public static double[][] QUSEN = CreateArray<double[]>(1, () => new double[1]);                                           //Large scale value of spec. humidity at eastern boundary at time t+1
        public static double[][] USW = CreateArray<double[]>(1, () => new double[1]);                                             //Large scale value of u-wind at western boundary at time t
        public static double[][] VSW = CreateArray<double[]>(1, () => new double[1]);                                             //Large scale value of v-wind at western boundary at time t
        public static double[][] WSW = CreateArray<double[]>(1, () => new double[1]);                                             //Large scale value of w-wind at western boundary at time t
        public static double[][] TSW = CreateArray<double[]>(1, () => new double[1]);                                             //Large scale value of pot. temp. at western boundary at time t
        public static double[][] QUSW = CreateArray<double[]>(1, () => new double[1]);                                            //Large scale value of spec. humidity at western boundary at time t
        public static double[][] USWN = CreateArray<double[]>(1, () => new double[1]);                                            //Large scale value of u-wind at western boundary at time t+1
        public static double[][] VSWN = CreateArray<double[]>(1, () => new double[1]);                                            //Large scale value of v-wind at western boundary at time t+1
        public static double[][] WSWN = CreateArray<double[]>(1, () => new double[1]);                                            //Large scale value of w-wind at western boundary at time t+1
        public static double[][] TSWN = CreateArray<double[]>(1, () => new double[1]);                                            //Large scale value of pot. temp. at western boundary at time t+1
        public static double[][] QUSWN = CreateArray<double[]>(1, () => new double[1]);                                           //Large scale value of spec. humidity at western boundary at time t+1
        public static double[][] USN = CreateArray<double[]>(1, () => new double[1]);                                             //Large scale value of u-wind at northern boundary at time t
        public static double[][] VSN = CreateArray<double[]>(1, () => new double[1]);                                             //Large scale value of v-wind at northern boundary at time t
        public static double[][] WSN = CreateArray<double[]>(1, () => new double[1]);                                             //Large scale value of w-wind at northern boundary at time t
        public static double[][] TSN = CreateArray<double[]>(1, () => new double[1]);                                             //Large scale value of pot. temp. at northern boundary at time t
        public static double[][] QUSN = CreateArray<double[]>(1, () => new double[1]);                                            //Large scale value of spec. humidity at northern boundary at time t
        public static double[][] USNN = CreateArray<double[]>(1, () => new double[1]);                                            //Large scale value of u-wind at northern boundary at time t+1
        public static double[][] VSNN = CreateArray<double[]>(1, () => new double[1]);                                            //Large scale value of v-wind at northern boundary at time t+1
        public static double[][] WSNN = CreateArray<double[]>(1, () => new double[1]);                                            //Large scale value of w-wind at northern boundary at time t+1
        public static double[][] TSNN = CreateArray<double[]>(1, () => new double[1]);                                            //Large scale value of pot. temp. at northern boundary at time t+1
        public static double[][] QUSNN = CreateArray<double[]>(1, () => new double[1]);                                           //Large scale value of spec. humidity at northern boundary at time t+1
        public static double[][] USS = CreateArray<double[]>(1, () => new double[1]);                                             //Large scale value of u-wind at southern boundary at time t
        public static double[][] VSS = CreateArray<double[]>(1, () => new double[1]);                                             //Large scale value of v-wind at southern boundary at time t
        public static double[][] WSS = CreateArray<double[]>(1, () => new double[1]);                                             //Large scale value of w-wind at southern boundary at time t
        public static double[][] TSS = CreateArray<double[]>(1, () => new double[1]);                                             //Large scale value of pot. temp. at southern boundary at time t
        public static double[][] QUSS = CreateArray<double[]>(1, () => new double[1]);                                            //Large scale value of spec. humidity at southern boundary at time t
        public static double[][] USSN = CreateArray<double[]>(1, () => new double[1]);                                            //Large scale value of u-wind at southern boundary at time t+1
        public static double[][] VSSN = CreateArray<double[]>(1, () => new double[1]);                                            //Large scale value of v-wind at southern boundary at time t+1
        public static double[][] WSSN = CreateArray<double[]>(1, () => new double[1]);                                            //Large scale value of w-wind at southern boundary at time t+1
        public static double[][] TSSN = CreateArray<double[]>(1, () => new double[1]);                                            //Large scale value of pot. temp. at southern boundary at time t+1
        public static double[][] QUSSN = CreateArray<double[]>(1, () => new double[1]);                                           //Large scale value of spec. humidity at southern boundary at time t+1
        public static double[][] Alfa = CreateArray<double[]>(1, () => new double[1]);                                            //Inclination angle of the surface fall-line (=line of greates slope)
        public static double[][] Beta = CreateArray<double[]>(1, () => new double[1]);                                            //Azimuth angle of fall-line (=line of greates slope); North = 0
        public static double[][] Q = CreateArray<double[]>(1, () => new double[1]);                                               //1-Q = Sky-view factor (Horizontüberhöhung) in the radiation model
        public static double[][] YG = CreateArray<double[]>(1, () => new double[1]);                                              //Integrated visibility in the radiation model
        public static double[][] R = CreateArray<double[]>(1, () => new double[1]);                                               //Integrated water vapour in the radiation model
        public static double[][][][] W1rad = CreateArray<double[][][]>(1, () => CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1])));  //influence of clouds on the solar incoming radiation
        public static double[][] W2rad = CreateArray<double[]>(1, () => new double[1]);                                           //influence of rain on the solar incoming radiation
        public static double[][][][] Wrad = CreateArray<double[][][]>(1, () => CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1])));   //? something in the radiation model
        public static double[][][] RADIATION = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));   //Temporary Radiation in TIMP
        public static double[] Dz1 = new double[1];                                                                               //Lower cell height of the radiation model
        public static double[] Dz2 = new double[1];                                                                               //Upper cell height of the radiation model
        public static double[][] Ag = CreateArray<double[]>(1, () => new double[1]);                                              //Surface albedo in the radiation model
        public static double[] Tstern = new double[1];                                                                            //Linke turbidity coefficient (~1 - 8)
        public static double[][][] eS = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Cloud parameter
        public static double[][][] eSeS = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));        //Cloud parameter
        public static double[] IG = new double[1];                                                                                //Shortwave solar radiation
        public static double[][][] EG = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Total solar radiation density
        public static double[][][] SG = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Direct solar radiation density
        public static double[][][] DG = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Diffusive solar radiation density
        public static double[][][] nS = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Transmission function
        public static double[][][] nD = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));          //Transmission function
        public static double[][][] Tau = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));         //Temporary field for cloud
        public static double[][][] Arad = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));        //Temporary field in the radiatio model
        public static double[] Myx = new double[1];                                                                               //Temporary field in the radiatio model
        public static double[][] Rp = CreateArray<double[]>(1, () => new double[1]);                                              //Integrated water vapour
        public static double[][] CGp = CreateArray<double[]>(1, () => new double[1]);                                             //Integrated CO2
        public static double[][][] L_Strich = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));    //Downward terrestrial radiation
        public static double[][] RTerrG = CreateArray<double[]>(1, () => new double[1]);                                          //Total incoming terrestrial radiation
        public static double[][][] EpsAp = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));       //Temporary field in the radiatio model
        public static double[][][] EpsAm = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));       //Temporary field in the radiatio model
        public static double[][][] Tab = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));         //Temperature of the black cloud (at top of the atmosphere)
        public static double[][][] Eab = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));         //Emissivity? of the black cloud (at top of the atmosphere)
        public static Int16[][] stabilityclass = CreateArray<Int16[]>(1, () => new Int16[1]);                                     //Stability class (1-7)
        public static float Solar_reduction_factor = 1.0F;                                                                        //Factor to reduce the solar radiation

        public static Int32 NDATUM;                   //Starting date of computation
        public static Int32 ISTUNDE;                  //Starting hour of computation
        public static double DT;                      //Maximum allowed time step
        public static double TLIMIT2;                 //Total modelling time
        public static Int32 IOUTPUT;                  //Write intermediate output files after IOUTPUT seconds
        public static double DELW;                    //Historic value not used currently
        public static double QUINIT;                  //Initial relative humidity
        public static double ZSEEH;                   //Height of the lowest compuation height - not used anymore
        public static double TINIT;                   //Temperature at 2m height
        public static double TGRAD;                   //Gradient of potential temperature - not used anymore
        public static double ZNEUT;                   //Height of the neutral boundary layer - not used anymore
        public static double TBINIT;                  //Surface temperature
        public static double TBINIT1;                 //Soil temperature in 1m depth -> Note that TBINIT1, TBINIT and TINIT are not used by their absolute values but by their differences to define the initial state of the atmospheric and soil temperatures
        public static double BGRAD;                   //Latitude
        public static Int32 IRAD;                     //Number of timesteps after which the radiation is updated
        public static Int32 IDEBUG;                   //Debug level 0 none, 3 highest
        public static Int32 IFLAGS1;                  //Compute u v w pn t gw fo
        public static Int32 IFLAGS2;                  //Compute br pr qu psi te tb str
        public static double DIA;                     //Flag to determine whether initale state of the model is computed using the diagnostic model
        public static double DIAI;                    //Not used anymore - former flag to switch between explicit and implicit solutaion algorithm
        public static double RELAXV;                  //Underrelaxation factor for momentum equations
        public static double RELAXT;                  //Underrelaxation factor for passive scalar equations
        public static Int32 ICATAFORC;                //Forcing of catabatic flows with -40/-25/-15 W/m² AKLA 7/6/5
        public static float ERA5BOUND;                //Defines the interval after which the boundary conditions are updated with certain ERA5 data
        public static float ERA5BOUND_Threshold;      //Counter for new boundary conditions with certain ERA5 data
        public static Int32 IBOUND;                   //BOUNDARY CONDITION (1,6,10,11,12)
        public static Int32 ISTAT;                    //NON-STEADYSTATE/FORCING W. "METEOPGT.ALL" YES=1
        public static Int32 MM5IFU;                   //NON-STEADYSTATE/FORCING W. IFU-MM5 DATA YES=1 -> not yet implemented in this version
        public static Int32 NESTGRAMM;                //GRAMM IN GRAMM NESTING YES=1 -> not yet implemented in this version
        public static Int32 IOUT;                     //flow out as binary-stream, little_endian, with header 0, or 3 without header
        public static double DTI;                     //Total simulation time
        public static double DTMAX;                   //Maximum allowed time step
        public static Int32 IPROC;                    //Maximum number of processors/threads used in each parallelized procedure
        public static Int32 IJAHR;                    //Year of the simulation (only for transient simulations)
        public static Int32 IJAHR4digits;             //Year of the simulation (only for transient simulations)
        public static Int32 IMON;                     //Month of the simulation start (only for transient simulations)
        public static Int32 ITAG;                     //Day of the simulation start (only for transient simulations)
        public static Int32 ISTU;                     //Hour of the simulation start (only for transient simulations)
        public static Int32 IMIN;                     //Minute of the simulation start (only for transient simulations)
        public static Int32 IKOOA;                    //Western border of model domain
        public static Int32 JKOOA;                    //Southern border of model domain
        public static float AHMIN;                   //Minimum surface elevation
        public static float AHMAX;                   //Maximum surface elevation
        public static float Inversion_Height;        //Inversion height
        public static float Wind_Velocity;            //Velocity
        public static Int32 AHMINI;                   //Cell indices of minimum surface elevation
        public static Int32 AHMINJ;                   //Cell indices of minimum surface elevation
        public static Int32 AHMINI_D;                   //Cell indices of minimum surface elevation within calculation domain
        public static Int32 AHMINJ_D;                   //Cell indices of minimum surface elevation within calculation domain
        public static List<int> inrec = new List<int>();   //Index in x-direction of receptor points
        public static List<int> jnrec = new List<int>();   //Index in y-direction of receptor points
        public static List<int> knrec = new List<int>();   //Index in z-direction of receptor points
        public static List<double> Xrec = new List<double>();  //x-coordinate of receptor points
        public static List<double> Yrec = new List<double>();  //y-coordinate of receptor points
        public static List<double> Zrec = new List<double>();  //z-coordinate of receptor points
        public static List<double> Urec = new List<double>();  //u-component at receptor points
        public static List<double> Vrec = new List<double>();  //v-component at receptor points
        public static List<double> Trec = new List<double>();  //temperature at receptor points
        public static List<double> Globradrec = new List<double>();  //global radiation at receptor points
        public static List<double> Longradrec = new List<double>();  //longwave outgoing radiation at receptor points
        public static List<double> Soilheatfluxrec = new List<double>();  //soil heat flux at receptor points
        public static List<double> Sensheatfluxrec = new List<double>();  //sensible heat flux at receptor points
        public static List<double> Latheatfluxrec = new List<double>();  //latent heat flux at receptor points
        public static double STRETCH;                 //Stretching-factor in the vertical
        public static double FN;                      //Coriolis frequency normal to the earths surface
        public static double FH;                      //Coriolis frequency parallel to the earths surface
        public static double[] PSAT = new double[1];  //saturation pressure
        public static double GERD;                    //Gravitational acceleration
        public static double GASCON;                  //General gas constant
        public static double CPLUFT;                  //Heat capacity of air for constant pressure
        public static double ALW;                     //Evaporation heat of water
        public static Int16 ISOL;                     //Switch for radiation model: 1=thin clouds 2=thick clouds
        public static int IWETTER;                    //Number of the actual calculated flow field situation
        public static string METEO;                   //Character determining whether GRAMM is driven by the file meteopgt.all or not
        public static double ANEMO;                   //Anemometer heat when using meteopgt.all as input file
        public static Int16 AKLA;                     //Stability class when using meteopgt.all as input file
        public static double Windspeed_meteopgt;      //Wind speed as stored in meteopgt.all used for determining the correct global radiation
        public static double TLIMIT;                  //Modified total simulation time in dependence on the stability class
        public static double Obini;                   //Initial Obukhov length
        public static double Rauigkeit;               //Initial surface-roughness length
        public static double CK;                      //van Karman constant
        public static double VISEL;                   //minimum turbulent viscosity
        public static double CPBOD;                   //heat capacity of soil
        public static double SIGMA;                   //Stefan Bolzmann constant
        public static double TBZ1;                    //temperature value used to improve numerical accuracy of the solution algorithm for temperature
        public static Int16 IPGT;                     //Flag used to check whether the file meteopgt.all is used as input file
        public static double TJETZT;                  //specific modelling time format used for the radiation model
        public static double REALTIME;                //integrated modelling time
        public static Int16 ITIME;                    //total time steps
        public static Int16 INUMS;                    //number of pressure-iterations
        public static double DIVSUM, SUMG;            //total mass divergence
        public static Int16 IDIV;                     //counter for computing the trend of the mass-divergence used to adjust the time-step
        public static double STEIGUNG;                //temporal trend of the total mass-divergence used to compute the actual time-step
        public static double MASSOURCE_FIRST;         //initial total mass-divergence at the beginning of each computation
        public static double PRTE;                    //turbulent Prandtl-number

        public static Boolean ICU;                    //Flag switching the compuation of the u-component on/off
        public static Boolean ICV;                    //Flag switching the compuation of the v-component on/off
        public static Boolean ICW;                    //Flag switching the compuation of the w-component on/off
        public static Boolean ICPN;                   //Flag switching the compuation of the non-hydrostatic pressure on/off
        public static Boolean ICT;                    //Flag switching the compuation of the pot. temperature on/off
        public static Boolean ICPH;                   //Flag switching the compuation of the hydrostatic pressure on/off (not used)
        public static Boolean IFOU;                   //Flag switching the compuation of the Fourier-Transformed upper boundary conditions on/off (not used)
        public static Boolean ICQU;                   //Flag switching the compuation of the spec. humidity on/off
        public static Boolean ICPSI;                  //Flag switching the compuation of a passive scalar on/off (not used)
        public static Boolean ICTE;                   //Flag switching the compuation of the turbulent kinetic energy on/off
        public static Boolean ICSTR;                  //Flag switching the compuation of the radiation model on/off
        public static Boolean ICPR;                   //Flag switching the compuation of the surface layer model on/off
        public static Boolean ICBR;                   //Flag switching the compuation of the boundary conditions on/off
        public static Boolean ICTB;                   //Flag switching the compuation of the surface temperature on/off
        public static Boolean ICGW;                   //Flag switching the compuation of the geostrophic winds from the large scale model when run in nesting mode on/off (not used)
        public static Boolean recexist;               //Flag set when receptor points are set (file receptor_GRAMM.dat)
        public static Boolean WriteSteadyState = false;  // Flag if Steady-State criterion should be written to a file

        public static int Counter = 0;                // Counter for writing "DispGRAMM.txt"
        public static int TerminalOut = 0;            // Counter for Terminal Output
        public static int TerminalThreshold = 1;      // Threshold for terminal output
        public static int nr_cell_smooth = 0;         // Number of cells used at the lateral boundaries for smoothing the orography

        public static bool GRAMM_Online_flag = true;
        public static bool unix = false;
        public static int computation_retry = 0;              // flag, if a compuatation is repeated caused by numerical instabilities
        public static int meteopgt_nr = 0;                    // original number of weather situations in the file meteopgt.all
        public static double max_time_step_original = 0; // the original maximum time step

        public static Queue<double> MASSOURCE_Queue = new Queue<double>(); // 5.4.2017 Ku Massource queue
        public static int IDIV_Up;                                         // 5.4.2017 Ku IDIV up Counter
        public static int IDIV_LockDown;                                   // 5.4.2017 Ku IDIV Lock to increase time step
        public static int IDIV_LockUp;                                     // 5.4.2017 Ku IDIV Lock to increase time step
        public static int IDIV_LockDown2;                                  // 5.4.2017 Ku IDIV Lock to increase time step
        public static int IDIV_LockRelax;                                  // 5.4.2017 Ku IDIV Lock to avoit too low relax factors
        public static int IDIV_PingPong;                                   // 6.4.2017 Ku IDIV Lock to avoit PingPong effects
        public static double MASSOURCE_Act;                                // 5.4.2017 Ku MASSOURCE actual
        public static double MASSOURCE_Old;                                // 5.4.2017 Ku MASSOURCE old
        public static double MASSOURCE_minusone;                           // 5.4.2017 ÖT MASSOURCE previous time step
        private static int IWetter_Console_First = 0;                      // 11.4.2017 Ku first weather situation from console
        private static int IWetter_Console_Last = 9999999;                 // 11.4.2017 Ku last weather situation from console
        private static double Relaxv_Console = 0;                          // 25.4.2017 Ku Relaxv from console
        private static double Relaxt_Console = 0;                          // 25.4.2017 Ku Relaxt from console
        private static double MaxTimeStep_Console = 0;                     // 25.4.2017 Ku MaxTimeStep from console

        public static double WU1 = 0;       //large-scale forcing using meteopgt.all
        public static double WU2 = 0;       //large-scale forcing using meteopgt.all
        public static double WV1 = 0;       //large-scale forcing using meteopgt.all
        public static double WV2 = 0;       //large-scale forcing using meteopgt.all

        public static double Relaxv_ori;
        public static double Relaxt_ori;
        public static double Divergence_Min = 10e9;
        public static float GRAMM_West = 0;
        public static float GRAMM_South = 0;

        public static DateTime ERA5_date1 = new DateTime();  //first date used for interpolation ERA5 data to force/initialize GRAMM
        public static DateTime ERA5_date2 = new DateTime();  //second date used for interpolation ERA5 data to force/initialize GRAMM
        public static float Longitude = 47;                  //avearge longitude of GRAMM domain

        public static int Intermed_Threshold;                //time interval for intermediate GRAMM flow field output
        public static float REINITIALIZATION;                //Defines the time interval in seconds after which GRAMM is completely re-initialized with ERA5 data
        public static float REINITIALIZATION_Threshold;      //Counter for re-initialization

        static void Main(string[] args)
        {
            #if __ECMWF__
            GdalConfiguration.ConfigureOgr();
            GdalConfiguration.ConfigureGdal();
            #endif
            // Arguments: "First Situation" "Final Situation" "Max. Time Step" "RelaxV" "RelaxT" 
            // if "First Situation" = "?" -> print Info & Stop

            int p = (int)Environment.OSVersion.Platform;

            if ((p == 4) || (p == 6) || (p == 128))
            {
                //Console.WriteLine ("Running on Unix");
                unix = true;
            }
            else
            {
                //Console.WriteLine ("NOT running on Unix");
            }

            //WRITE GRAMM VERSION INFORMATION TO SCREEN
            Console.WriteLine("");
            Console.WriteLine("+------------------------------------------------------+");
            Console.WriteLine("|                                                      |");
            string Info = "+  > >         G R A M M VERSION: 20.01          < <   +";
            Console.WriteLine(Info);
            if (unix)
            {
                Console.WriteLine("|                     L I N U X                        |");
            }
#if NETCOREAPP2_1 || NETCOREAPP2_0 || NETCOREAPP3_0
			    Console.WriteLine("|                 .Net Core Version                    |");
#endif
            Console.WriteLine("+------------------------------------------------------+");
            Console.WriteLine("");

            //show licence terms
            ShowCopyright(args);

            // 11.04.17 Ku use arguments
            Console_Arguments(args);

            //User defined decimal seperator
            decsep = NumberFormatInfo.CurrentInfo.NumberDecimalSeparator;

            //read number of grid cells stored in the file "GRAMM.geb"
            Read_Gramm_Geb();

            // Write to "Logfile_GRAMMCore"
            try
            {
                ProgramWriters.LogfileGrammCoreWrite(new String('-', 80));
                ProgramWriters.LogfileGrammCoreWrite(Info);
                ProgramWriters.LogfileGrammCoreWrite("Computation started at: " + DateTime.Now.ToString());
                ProgramWriters.LogfileGrammCoreWrite("Computation folder:     " + Directory.GetCurrentDirectory());
                Info = "Application hash code:  " + GetAppHashCode();
                ProgramWriters.LogfileGrammCoreWrite(Info);
            }
            catch { }

            //Allocate Memory -> Define Arrays
            Define_Arrays();
            List<int>[] Month_List = new List<int>[3]; // 3 types of month lists
            Month_List[0] = new List<int>() { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
            Month_List[1] = new List<int>() { 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 2 };
            Month_List[2] = new List<int>() { 6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5 };

            if (NX * NY < 50 * 50)
            {
                TerminalThreshold = 12;
            }
            else if (NX * NY < 200 * 200)
            {
                TerminalThreshold = 6;
            }
            else
            {
                TerminalThreshold = 1;
            }

            Console.WriteLine(".");
            //END OF VARIABLE DECLARATION BLOCK

            //SET NUMBER OF FLOW FIELD SITUATION TO ZERO
            IWETTER = 0;

            //OPEN THE MAIN CONTROL FILE "IIN.dat"
            Read_IIN_Dat();

            //When using ERA5 data GRAMM automatically starts with a number for the weather situation based on the existing *wnd files
            if (ISTAT >= 2)
            {
                DirectoryInfo d = new DirectoryInfo(@"..\Computation");
                FileInfo[] Files = d.GetFiles("*.wnd");
                IWETTER = Files.Length;
            }

            //Total simulation time
            DTI = TLIMIT2;
            DTMAX = DT;
            max_time_step_original = DTMAX;

            if (MaxTimeStep_Console > 0.99)
                max_time_step_original = MaxTimeStep_Console;
            if (Relaxv_Console > 0.0099)
                RELAXV = Relaxv_Console;
            if (Relaxt_Console > 0.0099)
                RELAXT = Relaxt_Console;

            // Read number of max. used processors
            IPROC = 1;
            Max_Proc_File_Read();

            //Convert relative humidity in %
            QUINIT *= 0.01;

            //Get year, month, day, hour, and minute of the simulation start
            if (NDATUM.ToString().Length == 8)
            {
                IJAHR4digits= Convert.ToInt32(Math.Floor(NDATUM / 10000d));
                IJAHR = IJAHR4digits - Convert.ToInt32(Math.Floor(IJAHR4digits / 100d)) * 100;
                IMON = Convert.ToInt32(Math.Floor(NDATUM / 100d)) - IJAHR4digits * 100;
                ITAG = NDATUM - IJAHR4digits * 10000 - IMON * 100;
            }
            else
            {
                IJAHR = Convert.ToInt32(Math.Floor(NDATUM / 10000d));
                IJAHR4digits = IJAHR;
                IMON = Convert.ToInt32(Math.Floor(NDATUM / 100d)) - IJAHR * 100;
                ITAG = NDATUM - IJAHR * 10000 - IMON * 100;
            }
            ISTU = Convert.ToInt32(Math.Floor(ISTUNDE / 100d));
            IMIN = ISTUNDE - 100 * ISTU;

            //Set flags for various compuation options
            ICU = false;
            ICV = false;
            ICW = false;
            ICPN = false;
            ICT = false;
            ICPH = false;
            IFOU = false;
            ICQU = false;
            ICPSI = false;
            ICTE = false;
            ICSTR = false;
            ICPR = false;
            ICBR = false;
            ICTB = false;
            ICGW = false;

            Int32 INUMM = 0;
            if (IFLAGS1 / 1000000 == 1) ICU = true;
            if (ICU) INUMM += 1000000;
            if ((IFLAGS1 - INUMM + 1) / 100000 == 1) ICV = true;
            if (ICV) INUMM += 100000;
            if ((IFLAGS1 - INUMM + 1) / 10000 == 1) ICW = true;
            if (ICW) INUMM += 10000;
            if ((IFLAGS1 - INUMM + 1) / 1000 == 1) ICPN = true;
            if (ICPN) INUMM += 1000;
            if ((IFLAGS1 - INUMM + 1) / 100 == 1) ICT = true;
            if (ICT) INUMM += 100;
            if ((IFLAGS1 - INUMM + 1) / 10 == 1) ICGW = true;
            if (ICGW) INUMM += 10;
            if ((IFLAGS1 - INUMM) == 1) IFOU = true;

            INUMM = 0;
            if (IFLAGS2 / 1000000 == 1) ICBR = true;
            if (ICBR) INUMM += 1000000;
            if ((IFLAGS2 - INUMM + 1) / 100000 == 1) ICPR = true;
            if (ICPR) INUMM += 100000;
            if ((IFLAGS2 - INUMM + 1) / 10000 == 1) ICQU = true;
            if (ICQU) INUMM += 10000;
            if ((IFLAGS2 - INUMM + 1) / 1000 == 1) ICPSI = true;
            if (ICPSI) INUMM += 1000;
            if ((IFLAGS2 - INUMM + 1) / 100 == 1) ICTE = true;
            if (ICTE) INUMM += 100;
            if ((IFLAGS2 - INUMM + 1) / 10 == 1) ICTB = true;
            if (ICTB) INUMM += 10;
            if ((IFLAGS2 - INUMM) == 1) ICSTR = true;

            //COMPUTE MODEL GRID
            Console.WriteLine("  *** GENERATING MODEL GRID *** ");
            GEOM();

            //INQUIRE IF RECEPTOR POINTS ARE SET
            Read_Receptor_Dat();

            Analyze_Topography(); // find U valleys and bassins

            Relaxv_ori = RELAXV;
            Relaxt_ori = RELAXT;
            Relax_Border_factor[0][0] = -4; // flag, that factor must be computed

            ProgramWriters.LogfileGrammCoreInfo();

            //Loop_______________________________________________________________________________________________
            NEXTWEATHERSITUATION:

            clear_arrays();

            //INITIALIZE FIELDS
            Console.WriteLine();
            Console.WriteLine();
            Console.WriteLine("  *** INITIALIZATION *** ");
            INITA(NX, NY, NZ);

            Int32 ISTUD = Convert.ToInt32(Math.Floor(TLIMIT / 3600d));
            double AMIND = (TLIMIT / 3600 - (float)ISTUD) * 60;
            Int32 IMIND = Convert.ToInt32(Math.Floor(AMIND));
            double ASECD = (AMIND - (float)IMIND) * 60;
            Int32 ISECD = Convert.ToInt32(Math.Floor(ASECD));

            if (ISTAT >= 2)
            {
                //GRAMM uses sun time -> UTC has to be transferred
                if (Longitude > 180.0)
                    Longitude -= 180;
                if (Longitude < -180.0)
                    Longitude += 180;
                ISTU += Convert.ToInt32(Longitude * 12 / 180.0);
                if (ISTU >= 24)
                {
                    ISTU -= 24;
                    DateTime dateref = new DateTime(Program.IJAHR4digits, Program.IMON, Program.ITAG, Program.ISTU, Program.IMIN, 0);
                    dateref = dateref.AddDays(1);
                    ITAG = dateref.Day;
                    IMON = dateref.Month;
                }
            }

            if (ISTAT == 0)
            {
                Intermed_Threshold = IOUTPUT / 2;  // 1st intermediate output after IOUTPUT / 2 (default 1800 s) for SC 5, 6 and SC 7
                if (AKLA < 5)
                {
                    Intermed_Threshold = IOUTPUT;  // 1st intermediate output after IOUTPUT s
                }
            }
            else
            {
                Intermed_Threshold = IOUTPUT;
            }
            int Radiation_Threshold = 1800; // compute new radiation each 1800 s

            //set further initial values
            if (ISTAT < 2)
                INITB.INIT(NX, NY, NZ);

            Relax_Border(); // set reduction factors for the border cells

            //initialize the radiation model
            if (((METEO != "Y") && (METEO != "y")) || (Program.ISTAT != 0)) TLIMIT = TLIMIT2;
            Boolean LGEOM = false;
            Boolean LGEOMW = false;
            Boolean LGEOMR = false;
            INIT_Radiation(ref LGEOM, ref LGEOMW, ref LGEOMR, ref ISTUD, ref AMIND, ref IMIND, ref ASECD, ref ISECD, Month_List);
            
            Console.WriteLine(ITAG.ToString() + "." + IMON.ToString() + "  -  " + ISTU.ToString() + ":" + IMIND.ToString("D2"));

            REALTIME = 0;
            ITIME = 0;
            INUMS = 0;
            DIVSUM = 0;
            IDIV = 0;
            DT = 1.5;
            for (int i = 0; i < 11; i++)
            {
                MASSOURCE[i] = 0;
            }

            //START OF THE INTEGRATION LOOP
            while (REALTIME < TLIMIT)
            {
                //number of iteration
                ITIME++;

                if (ITIME % 20 == 0) Max_Proc_File_Read(); // read MaxProc at 20th time step

                //total simulation expressed in seconds
                TJETZT = (float)ISTU * 3600 + (float)IMIN * 60 + REALTIME;

                //write normalised actual time used in the GUI for visualisation of the modelling progress
                if ((ITIME % 10) == 0 && (IWetter_Console_First <= 1)) // 11.4.17 Ku use arguments
                {
                    try
                    {
                        using (StreamWriter w = new StreamWriter("PercentGramm.txt"))
                        {
                            double tnorm = REALTIME * TLIMIT2 / TLIMIT;
                            w.WriteLine(tnorm.ToString());
                        }
                    }
                    catch { }
                }


                //Online output of fields for GUI
                if (GRAMM_Online_flag || (ITIME % 10d) == 0)
                {
                    GRAMM_Online_flag = false; // reset flag
                    GrammOnline(NX, NY, NZ);   // flag is at GRAMMOnline set to true, if only output is necessary
                }

                //Implicit solution algorith for conservation equations (SIMPLE - Patankar, 1980)
                if (REALTIME <= DTI)
                {
                    if (SOLUTION(NX, NY, NZ) == false && ISTAT == 0) //  numerical problems using overall massdivergence
                    {
                        computation_retry++;
                        if (computation_retry < 3) // try 3 times -> otherwise let the app crash
                        {
                            IWETTER--; //try same situation
                            TLIMIT += TLIMIT2;
                            DTI += TLIMIT2;
                            clear_arrays();
                            GEOM();
                            goto NEXTWEATHERSITUATION;
                        }

                    }
                }

                //new radiation data
                if (ITIME == 1)
                {
                    LGEOM = true;
                    LGEOMW = true;
                    LGEOMR = false;
                }
                else
                {
                    LGEOM = false;
                    LGEOMW = false;
                    LGEOMR = false;
                }

                if (ICSTR == true)
                {
                    if ((ITIME == 1) || ((ITIME % IRAD) == 0))
                    {
                        double TJETZT1 = TJETZT;
                        int IMIN_RAD = IMIN;
                        int ITAG_RAD = ITAG;
                        int IMON_RAD = IMON;
                        if (ISTAT == 0)
                        {
                            ISTUD = Convert.ToInt32(Math.Floor((TJETZT + 10) / 3600d));
                        }
                        else if(ISTAT == 1)
                        {
                            ISTUD = ISTU;
                            TJETZT1 = (float)ISTU * 3600 + (float)IMIN * 60;
                        }
                        else if (ISTAT >= 2)
                        {
                            //ISTUD = Convert.ToInt32(Math.Floor((TJETZT + 10) / 3600d));
                            //GRAMM uses sun time -> UTC has to be transferred
                            DateTime dateref = new DateTime(Program.IJAHR4digits, Program.IMON, Program.ITAG, Program.ISTU, Program.IMIN, 0);
                            dateref = dateref.AddSeconds(REALTIME);
                            IMON_RAD = dateref.Month;
                            ITAG_RAD = dateref.Day;
                            ISTUD = dateref.Hour;
                            IMIND = dateref.Minute;
                            TJETZT1 = ISTUD * 3600;
                        }
                        AMIND = (TJETZT / 3600 - (float)ISTUD) * 60;
                        IMIND = Convert.ToInt32(Math.Floor(AMIND));
                        ASECD = (AMIND - (float)IMIND) * 60;
                        ISECD = Convert.ToInt32(Math.Floor(ASECD));

                        for (int II = 1; II <= 1000; II++)
                        {
                            if (TJETZT1 >= 86400) TJETZT1 -= 86400;
                        }
                        RADIATRAD(LGEOM, LGEOMW, LGEOMR, ITAG_RAD, IMON_RAD, IJAHR, TJETZT1, NX, NY, NZ);
                        Console.WriteLine(ITAG_RAD.ToString() + "." + IMON_RAD.ToString() + "  -  " + ISTUD.ToString() + ":" + IMIND.ToString("D2"));
                    }

                    //dynamic sun (global radiation) every 1800 seconds
                    if (((METEO == "Y") || (METEO == "y")) && (Program.ISTAT == 0))
                    {
                        if (Program.AKLA < 4)
                        {
                            if (REALTIME > Radiation_Threshold)
                            {
                                Radiation_Threshold += 1800; // compute new radiation each 1800 s
                                double TJETZT1 = TJETZT;
                                ISTUD = Convert.ToInt32(Math.Floor((TJETZT + 10) / 3600d));
                                AMIND = (TJETZT / 3600 - (float)ISTUD) * 60;
                                IMIND = Convert.ToInt32(Math.Floor(AMIND));
                                ASECD = (AMIND - (float)IMIND) * 60;
                                ISECD = Convert.ToInt32(Math.Floor(ASECD));

                                for (int II = 1; II <= 1000; II++)
                                {
                                    if (TJETZT1 >= 86400) TJETZT1 -= 86400;
                                }
                                RADIATRAD(LGEOM, LGEOMW, LGEOMR, ITAG, IMON, IJAHR, TJETZT1, NX, NY, NZ);
                                Console.WriteLine(ITAG.ToString() + "." + IMON.ToString() + "  -  " + ISTUD.ToString() + ":" + IMIND.ToString("D2"));
                            }
                        }
                    }
                }

                //store last time step
                STOREcalculate(NX, NY, NZ);

                //intermediate output
                if (((METEO == "Y") || (METEO == "y")) && (Program.ISTAT == 0))
                {
                    //only in case of steady-state simulations using meteopgt.all no intermediate output is provided
                    if (Program.meteopgt_nr == 0)
                    {
                        IOUTPUT = 100000000;
                        Intermed_Threshold = 100000000;
                    }
                }

                Int16 IHOUR = Convert.ToInt16(Math.Floor(TJETZT / IOUTPUT));
                if (REALTIME > Intermed_Threshold && DTI > 3602) // if threshold is exceeded && simulation time > 1 h -> create intermed. outputs
                {
                    Console.Write(" INTERMEDIATE OUTPUT ");

                    // set new Threshold value
                    if (Intermed_Threshold < IOUTPUT)
                    {
                        Intermed_Threshold = IOUTPUT * 2; // 2nd intermediate output after 2 hours
                    }
                    else
                    {
                        Intermed_Threshold += IOUTPUT;
                    }

                    OUTPUT(NX, NY, NZ, true); // intermediate output

                    Console.WriteLine();
                    Console.WriteLine();

                    if (((METEO != "Y") && (METEO != "y")) || (Program.ISTAT != 0))
                    {
                        IWETTER++;
                    }
                }

                //increase simulation time
                REALTIME += DT;
                IHOURO = IHOUR;

                //save intermediate surface flow fields after 75% of the total simulatin time (VDI 3783-7)
                if ((REALTIME >= DTI * 0.75) && (REALTIME <= (DTI * 0.75 + DT * 2)))
                {
                    Parallel.For(1, NX + 1, Program.pOptions, i =>
                    {
                        for (int j = 1; j <= NY; j++)
                        {
                            U_TEMP[i][j] = (float)(U[i][j][1]);
                            V_TEMP[i][j] = (float)(V[i][j][1]);
                            W_TEMP[i][j] = (float)(W[i][j][1]);
                        }
                    });
                }
            }

            computation_retry = 0; // reset compuation retry counter
                                   //Ultimate output at the end of each situation
            Console.WriteLine("");
            Console.Write(" MMAIN : OUT ");
            if ((METEO != "Y") && (METEO != "y"))
            {
                OUTPUT(NX, NY, NZ, false); // final output
                Console.WriteLine();
            }
            else
            {
                OUTPUT(NX, NY, NZ, false); // final output
                Console.WriteLine();
                goto NEXTWEATHERSITUATION;
            }
        }

        //module to initialze a jagged array
        public static T[] CreateArray<T>(int cnt, Func<T> itemCreator)
        {
            T[] result = new T[cnt];
            for (int i = 0; i < result.Length; i++)
            {
                result[i] = itemCreator();
            }
            return result;
        }

        //licence GPL 3 terms
        static void ShowCopyright(string[] args)
        {
            Console.WriteLine("[GRAMM]  Copyright (C) <2019>  <Dietmar Oettl, Markus Kuntner>");
            Console.WriteLine("This program comes with ABSOLUTELY NO WARRANTY; for details start GRAMM with a startup parameter ‘show_w’");
            Console.WriteLine("This is free software, and you are welcome to redistribute it under certain conditions; start GRAMM with a startup parameter ‘show_c’ for details. )");

            if (args.Length > 0 && args[0].Contains("show_w"))
            {
                Console.WriteLine("This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of" +
                                  " MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.");
            }
            else if (args.Length > 0 && args[0].Contains("show_c"))
            {
                Console.WriteLine("This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by" +
                                  "the Free Software Foundation, either version 3 of the License.");
                Console.WriteLine("You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.");
            }
            Console.WriteLine();
        }

        private static bool Console_Arguments(string[] args)
        {
            int _off = 0;
            if (args.Length > 0)
            {
                if (args[0].Contains("?") == true || args[0].ToUpper().Contains("HELP") == true) // Info & stop
                {
                    Console.WriteLine("GRAMM console arguments: 'Working Directory' 'First Situation' 'Final Situation' 'Max. Time Step' 'RelaxV' 'RelaxT'");
                    Console.WriteLine("or");
                    Console.WriteLine("GRAMM console arguments: 'First Situation' 'Final Situation' 'Max. Time Step' 'RelaxV' 'RelaxT'");
                    Environment.Exit(0);
                }

                int temp;
                if (Int32.TryParse(args[0], out temp) == false) // not a valid number
                {
                    if (Directory.Exists(args[0]) == true) // arg[0] = Directory!
                    {
                        Directory.SetCurrentDirectory(args[0]);
                        _off = 1;
                    }
                }
            }
            if (args.Length > (1 + _off)) //  + 2 arguments -> first and last weather situation
            {
                if (Int32.TryParse(args[0 + _off], out IWetter_Console_First))
                    Int32.TryParse(args[1 + _off], out IWetter_Console_Last);
                if (IWetter_Console_Last < IWetter_Console_First || IWetter_Console_First < 1)
                {
                    IWetter_Console_First = 0;
                    IWetter_Console_Last = 9999999;
                }
                // Max. Time Step
                if (args.Length > (2 + _off))
                    if (Double.TryParse(args[2 + _off], NumberStyles.AllowDecimalPoint, CultureInfo.InvariantCulture, out MaxTimeStep_Console) == false)
                        MaxTimeStep_Console = 0;
                // RelaxV
                if (args.Length > (3 + _off))
                    if (Double.TryParse(args[3 + _off], NumberStyles.AllowDecimalPoint, CultureInfo.InvariantCulture, out Relaxv_Console) == false)
                        Relaxv_Console = 0;
                // RelaxT
                if (args.Length > (4 + _off))
                    if (Double.TryParse(args[4 + _off], NumberStyles.AllowDecimalPoint, CultureInfo.InvariantCulture, out Relaxt_Console) == false)
                        Relaxt_Console = 0;
            }

            if (IWetter_Console_First > 0) // 11.04.17 Ku use arguments
            {
                Console.WriteLine("");
                Console.WriteLine("Starting with situation " + IWetter_Console_First.ToString() +
                                  " ...final situation " + IWetter_Console_Last.ToString() + "    ");
            }
            if (IWetter_Console_First <= 1) // 11.04.17 Ku use arguments
            {
                try
                {
                    if (File.Exists("albeq.dat"))
                        File.Delete("albeq.dat");
                }
                catch { }
            }
            else
            {
                Console.Write("waiting for albeq.dat");
                while (File.Exists("albeq.dat") == false)
                {
                    System.Threading.Thread.Sleep(2000);
                    Console.Write(".");
                }
                System.Threading.Thread.Sleep(2000);
                Console.WriteLine("");
            }
            return true;
        }
        

        //computes wind direction
        public static double winkel(double u, double v)
        {
            double winkel = 0;

            if (v == 0)
                winkel = 90;
            else
                winkel = Math.Abs(Math.Atan(u / v)) * 180 / Math.PI;

            if ((v > 0) && (u <= 0)) winkel = 180 - winkel;
            if ((v >= 0) && (u > 0)) winkel = 180 + winkel;
            if ((v < 0) && (u >= 0)) winkel = 360 - winkel;

            return winkel;
        }

        private static int set_month_type(double windspeed, float AKLA, double Brgrad) // set start month for radiation search
        {
            int month_setting = 0; // Jan to Dec
            if (Brgrad > 0) // northern hemisphere
            {
                if (Windspeed_meteopgt > 2 || AKLA < 5) // higher wind vel. or labile or neutral SC -> start radiation search at Mar
                    month_setting = 1;
                else if (Windspeed_meteopgt > 4) // very high wind vel. -> start radiation search at Jun
                    month_setting = 2;
            }
            else // southern hemisphere
            {
                month_setting = 2;  // June to Mai
                if (Windspeed_meteopgt > 2 || AKLA < 5) // higher wind vel.  or labile or neutral SC -> start radiation search at Mar
                    month_setting = 1;
                else if (Windspeed_meteopgt > 4) // very high wind vel. -> start radiation search at Jan
                    month_setting = 0;
            }
            return month_setting;
        }

        /// <summary>
        /// Get a Hash code of the running app
        /// </summary>
        public static string GetAppHashCode()
        {
            string filename = System.Diagnostics.Process.GetCurrentProcess().MainModule.FileName;
            string hashstring = string.Empty;

            if (File.Exists(filename))
            {
                try
                {
                    using (System.Security.Cryptography.SHA1 sha = System.Security.Cryptography.SHA1.Create())
                    {
                        byte[] hash;

                        using (FileStream stream = File.OpenRead(filename))
                        {
                            hash = sha.ComputeHash(stream);
                        }

                        foreach (byte item in hash)
                        {
                            hashstring += item.ToString("x2");
                        }
                    }
                }
                catch { }
            }
            return hashstring;
        }

    }
}
