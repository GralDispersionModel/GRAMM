#region Copyright
///<remarks>
/// <GRAMM Mesoscale Model>
/// Copyright (C) [2019] [Dietmar Oettl, Markus Kuntner]
/// This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
/// the Free Software Foundation version 3 of the License
/// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
/// You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
///</remarks>
#endregion

using System;
using System.Collections.Generic;
using System.Threading.Tasks;

namespace GRAMM_2001
{
    partial class Program
    {
        //global variables declaration block
        //note that arrays are defined as jagged-arrays as they have a much better computational performance than matrices

        ///<summary>
        ///sets the maximum number of cores to be used in the parallelisation
        ///</summary>
        public static ParallelOptions pOptions = new ParallelOptions();
        ///<summary>
        ///global decimal separator of the system
        ///</summary>
        public static string decsep;
        ///<summary>
        ///number of cells in x-direction
        ///</summary>
        public static Int32 NX;
        ///<summary>
        ///number of cells in y-direction
        ///</summary>
        public static Int32 NY;
        ///<summary>
        ///number of cells in z-direction
        ///</summary>
        public static Int32 NZ;
        ///<summary>
        ///number of cells in x-direction +1
        ///</summary>
        public static Int32 NX1;
        ///<summary>
        ///number of cells in y-direction +1
        ///</summary>
        public static Int32 NY1;
        ///<summary>
        ///number of cells in z-direction +1
        ///</summary>
        public static Int32 NZ1;
        ///<summary>
        ///number of cells in x-direction +2
        ///</summary>
        public static Int32 NX2;
        ///<summary>
        ///number of cells in y-direction +2
        ///</summary>
        public static Int32 NY2;
        ///<summary>
        ///number of cells in z-direction +2
        ///</summary>
        public static Int32 NZ2;
        ///<summary>
        ///NUMBER OF VERTICAL GRID CELLS OF THE SOIL MODEL
        ///</summary>
        public static Int32 NZB;
        ///<summary>
        /// NUMBER OF VERTICAL GRID CELLS OF THE RADIATION MODEL
        ///</summary>
        public static Int32 NPROFMAX;
        ///<summary>
        ///Flag needed that geometry data for the radiation model is read just once
        ///</summary>
        public static Int16 IMETSTR;
        ///<summary>
        ///Height of the surface
        ///</summary>
        public static float[][] AH = CreateArray<float[]>(1, () => new float[1]);
        ///<summary>
        /// Height of the bassins and valleys
        ///</summary>
        public static float[][] AH_Bassins = CreateArray<float[]>(1, () => new float[1]);
        ///<summary>
        /// Height of the bassins and valleys
        ///</summary>
        public static float[][] TPI = CreateArray<float[]>(1, () => new float[1]);
        ///<summary>
        ///Flag needed for the intermediate output of GRAMM flowfields
        ///</summary>
        public static Int16 IHOURO;
        ///<summary>
        /// Volume of grid cells
        ///</summary>
        public static float[][][] VOL = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Area of the grid cell in x-direction
        ///</summary>
        public static float[][][] AREAX = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Area of the grid cell in y-direction
        ///</summary>
        public static float[][][] AREAY = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Bottom area of the grid cell
        ///</summary>
        public static float[][][] AREAZ = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Projection of the ground area of the grid cell in z-direction
        ///</summary>
        public static float[][][] AREA = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Area between the two halfs of the grid cell
        ///</summary>
        public static float[][][] AREAXYZ = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        /// Projection of the ground area of the grid cell in x-direction
        ///</summary>
        public static float[][][] AREAZX = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        /// Projection of the ground area of the grid cell in y-direction
        ///</summary>
        public static float[][][] AREAZY = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        /// Heights of the corner points of each grid cell
        ///</summary>
        public static float[][][] AHE = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        /// Height of the centre point of each grid cell
        ///</summary>
        public static float[][][] ZSP = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Horizontal grid size in x-direction
        ///</summary>
        public static float[] DDX = new float[1];
        ///<summary>
        ///Horizontal grid size in y-direction
        ///</summary>
        public static float[] DDY = new float[1];
        ///<summary>
        ///Distance between neighbouring grid cells in x-direction
        ///</summary>
        public static float[] ZAX = new float[1];
        ///<summary>
        ///Distance between neighbouring grid cells in y-direction
        ///</summary>
        public static float[] ZAY = new float[1];
        ///<summary>
        /// Distance of grid cell centre from western domain border
        ///</summary>
        public static float[] X = new float[1];
        ///<summary>
        /// Distance of grid cell centre from southern domain border
        ///</summary>
        public static float[] Y = new float[1];
        ///<summary>
        /// Temporary array for generating the terrain-following grid
        ///</summary>
        public static float[] Z = new float[1];
        ///<summary>
        ///CORINE land-use categories
        ///</summary>
        public static double[][] LAND = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        /// Base state for pressure used in the radiation model
        ///</summary>
        public static double[] PBZZ = new double[1];
        ///<summary>
        ///Base state for density used in the radiation model
        ///</summary>
        public static double[] RHOBZZ = new double[1];
        ///<summary>
        /// Absolute temperature in K
        ///</summary>
        public static double[][][] TABS = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        /// Solar incoming radiation without reflected radiation
        ///</summary>
        public static double[][] RSOLG = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Globar incoming radiation
        ///</summary>
        public static double[][] GLOBRAD = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Solar heating rate
        ///</summary>
        public static double[][][] DT_SOL = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        /// Longwave outgoing radiation
        ///</summary>
        public static double[][] RL = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Terrestrial heating rate
        ///</summary>
        public static double[][][] DT_TERR = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        /// Absolute soil temperature at the second level beneath the surface
        ///</summary>
        public static double[][] TG = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Surface emissivity
        ///</summary>
        public static double[][] EPSG = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Defines the cell in vertical direction from which onward the radiation model is computed
        ///</summary>
        public static int[][] KST = CreateArray<int[]>(1, () => new int[1]);
        ///<summary>
        ///Average value of two neighbouring cell heights
        ///</summary>
        public static float[] ZZ = new float[1];
        ///<summary>
        ///Surface albedo
        ///</summary>
        public static double[][] ALBEDO = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        /// Flag characterizing dry areas
        ///</summary>
        public static bool[][] DRY_AREA = CreateArray<bool[]>(1, () => new bool[1]);
        ///<summary>
        ///Cloudyness 0 - x - 1
        ///</summary>
        public static double[][] CLOUDS = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Snow cover 0/1
        ///</summary>
        public static double[][] SNOW = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Cell heights for the radiation model
        ///</summary>
        public static double[] ZPROF = new double[1];
        ///<summary>
        ///Vertical temperature profile in the radiation model
        ///</summary>
        public static double[] TPROF = new double[1];
        ///<summary>
        ///Vertical pressure profile in the radiation model
        ///</summary>
        public static double[] PPROF = new double[1];
        ///<summary>
        ///standard optical length
        ///</summary>
        public static double[] VNORM = new double[1];
        ///<summary>
        /// water content in the atmosphere used in radiation model
        ///</summary>
        public static double[] QVAP = new double[1];
        ///<summary>
        /// water content in clouds used in radiation model
        ///</summary>
        public static double[] QCLD = new double[1];
        ///<summary>
        ///water content in rain used in radiation model
        ///</summary>
        public static double[] QRAIN = new double[1];
        ///<summary>
        /// water content in ice used in radiation model
        ///</summary>
        public static double[] QICE = new double[1];
        ///<summary>
        ///water content in snow used in radiation model
        ///</summary>
        public static double[] QSNOW = new double[1];
        ///<summary>
        /// Average wind speed in x-direction over the two half-cells at time t
        ///</summary>
        public static double[][][] U = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        /// Average wind speed in y-direction over the two half-cells at time t
        ///</summary>
        public static double[][][] V = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        /// Average wind speed in z-direction over the two half-cells at time t
        ///</summary>
        public static double[][][] W = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Density of the air
        ///</summary>
        public static float[][][] RHO = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Average wind speed in x-direction over the two half-cells at time t+1
        ///</summary>
        public static double[][][] UN = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Average wind speed in y-direction over the two half-cells at time t+1
        ///</summary>
        public static double[][][] VN = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Average wind speed in z-direction over the two half-cells at time t+1
        ///</summary>
        public static double[][][] WN = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Average wind speed in x-direction in the lower half-cells at time t
        ///</summary>
        public static double[][][] U1 = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Average wind speed in y-direction in the lower half-cells at time t
        ///</summary>
        public static double[][][] V1 = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Average wind speed in z-direction in the lower half-cells at time t
        ///</summary>
        public static double[][][] W1 = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Average wind speed in x-direction in the upper half-cells at time t
        ///</summary>
        public static double[][][] U2 = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Average wind speed in y-direction in the upper half-cells at time t
        ///</summary>
        public static double[][][] V2 = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Average wind speed in z-direction in the upper half-cells at time t
        ///</summary>
        public static double[][][] W2 = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Average wind speed in x-direction in the lower half-cells at time t+1
        ///</summary>
        public static double[][][] U1N = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Average wind speed in y-direction in the lower half-cells at time t+1
        ///</summary>
        public static double[][][] V1N = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Average wind speed in z-direction in the lower half-cells at time t+1
        ///</summary>
        public static double[][][] W1N = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Average wind speed in x-direction in the upper half-cells at time t+1
        ///</summary>
        public static double[][][] U2N = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Average wind speed in y-direction in the upper half-cells at time t+1
        ///</summary>
        public static double[][][] V2N = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Average wind speed in z-direction in the upper half-cells at time t+1
        ///</summary>
        public static double[][][] W2N = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Non-hydrostatic pressure in the cell-centre
        ///</summary>
        public static double[][][] PN = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        /// Relaxation reduction factor at the border cells
        ///</summary>
        public static float[][] Relax_Border_factor = CreateArray<float[]>(1, () => new float[1]);

        ///<summary>
        ///Non-hydrostatic pressure in the cell-centre
        ///</summary>
        public static double[][][] DP = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Non-hydrostatic pressure in x-direction
        ///</summary>
        public static double[][][] DPX = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Non-hydrostatic pressure in y-direction
        ///</summary>
        public static double[][][] DPY = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Non-hydrostatic pressure in z-direction
        ///</summary>
        public static double[][][] DPZ = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Pressure gradient used for the velocity correction in x-direction
        ///</summary>
        public static double[][][] DDP1DX = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Pressure gradient used for the velocity correction in y-direction
        ///</summary>
        public static double[][][] DDP1DY = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Pressure gradient used for the velocity correction in z-direction
        ///</summary>
        public static double[][][] DDP1DZ = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Pressure gradient used for the velocity correction in x-direction
        ///</summary>
        public static double[][][] DDP2DX = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Pressure gradient used for the velocity correction in y-direction
        ///</summary>
        public static double[][][] DDP2DY = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Pressure gradient used for the velocity correction in z-directionNon-hydrostatic pressure in the cell-centre
        ///</summary>
        public static double[][][] DDP2DZ = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        /// thermic pressure gradient in x-direction
        ///</summary>
        public static double[][][] TPDX = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        /// thermic pressure gradient in y-direction
        ///</summary>
        public static double[][][] TPDY = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///thermal pressure
        ///</summary>
        public static double[][][] TP = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///thermal pressure in x-direction
        ///</summary>
        public static double[][][] TPX = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///thermal pressure in y-direction
        ///</summary>
        public static double[][][] TPY = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///thermal pressure in z-direction
        ///</summary>
        public static double[][][] TPZ = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        /// pressure at the top of the model domain
        ///</summary>
        public static double POBEN;
        ///<summary>
        ///Mass divergence in x-direction
        ///</summary>
        public static float[][][] SUX = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Mass divergence in y-direction
        ///</summary>
        public static float[][][] SUY = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Mass divergence in z-direction
        ///</summary>
        public static float[][][] SUZ = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        /// Mass divergence in between the two half-cells
        ///</summary>
        public static float[][][] SUXYZ = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Mass-flux in x-direction
        ///</summary>
        public static float[][][] U1NRHO = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Mass-flux in y-direction
        ///</summary>
        public static float[][][] U2NRHO = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Mass-flux in z-direction
        ///</summary>
        public static float[][][] V1NRHO = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Mass-flux in x-direction
        ///</summary>
        public static float[][][] V2NRHO = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Mass-flux in y-direction
        ///</summary>
        public static float[][][] W1NRHO = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Mass-flux in z-direction
        ///</summary>
        public static float[][][] W2NRHO = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Hydrostatic pressure
        ///</summary>
        public static float[][][] PBZ = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        /// Air densitiy at the simulation start
        ///</summary>
        public static float[][][] RHOBZ = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Temporary array in the calculation of the radiation fluxes
        ///</summary>
        public static double[][][] PNBZKP = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Temporary array in the calculation of the radiation fluxes
        ///</summary>
        public static double[][][] NBZKP = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Geostrophic wind in x-direction at time t-1
        ///</summary>
        public static float[][][] UG1 = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Geostrophic wind in y-direction at time t-1
        ///</summary>
        public static float[][][] VG1 = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Geostrophic wind in x-direction at time t+1
        ///</summary>
        public static float[][][] UG2 = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Geostrophic wind in y-direction at time t+1
        ///</summary>
        public static float[][][] VG2 = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Geostrophic wind in x-direction at time t
        ///</summary>
        public static double[][][] UG = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Geostrophic wind in y-direction at time t
        ///</summary>
        public static double[][][] VG = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        /// Turbulent viscosity in z-direction
        ///</summary>
        public static double[][][] VISV = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        /// Turbulent viscosity in horizontal-directions
        ///</summary>
        public static double[][][] VISH = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Specific humidity at time t
        ///</summary>
        public static double[][][] QU = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Specific humidity at time t+1
        ///</summary>
        public static double[][][] QUN = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Content of cloud water in the atmosphere in g/kg at time t
        ///</summary>
        public static double[][][] WAT_VAP = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        /// Content of cloud water in the atmosphere in g/kg at time t+1
        ///</summary>
        public static double[][][] WAT_VAPN = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Basic state of specific humidity at the simulation start
        ///</summary>
        public static double[][][] QBZ = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Soil moisture content
        ///</summary>
        public static double[][] QUG = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        /// Specific soil moisture parameter (e.g. water = 1)
        ///</summary>
        public static double[][] FW = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        /// Potential temperature of air at time t
        ///</summary>
        public static double[][][] T = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Potential temperature of air at time t+1
        ///</summary>
        public static double[][][] TN = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Passive scalar at time t
        ///</summary>
        public static double[][][][] PS = CreateArray<double[][][]>(1, () => CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1])));
        ///<summary>
        ///Passive scalar at time t+1
        ///</summary>
        public static double[][][][] PSN = CreateArray<double[][][]>(1, () => CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1])));
        ///<summary>
        ///temporary field for a passive scalar
        ///</summary>
        public static double[][][] PStemp = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///temporary field for a passive scalar
        ///</summary>
        public static double[][][] PSNtemp = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));

        ///<summary>
        ///Conversion factor between potential and absolute temperature of air
        ///</summary>
        public static float[][][] FACTOR = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        /// Turbulent dissipation rate at time t
        ///</summary>
        public static double[][][] DISS = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Turbulent dissipation rate at time t+1
        ///</summary>
        public static double[][][] DISSN = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Turbulent kinetic energy at time t
        ///</summary>
        public static double[][][] TE = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Turbulent kinetic energy at time t+1
        ///</summary>
        public static double[][][] TEN = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Mixing height
        ///</summary>
        public static float[][] ZI = CreateArray<float[]>(1, () => new float[1]);
        ///<summary>
        ///Fricition velocity
        ///</summary>
        public static double[][] UST = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Friction velocity normalized by the wind speed
        ///</summary>
        public static double[][] USTV = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Potential temperature scale
        ///</summary>
        public static double[][] TST = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Obukhov-length
        ///</summary>
        public static float[][] OL = CreateArray<float[]>(1, () => new float[1]);
        ///<summary>
        ///Temporary array used for calculating the sensible and latent heat fluxes
        ///</summary>
        public static float[][] XWQ = CreateArray<float[]>(1, () => new float[1]);
        ///<summary>
        ///Aerodynamic surface roughness
        ///</summary>
        public static float[][] Z0 = CreateArray<float[]>(1, () => new float[1]);
        ///<summary>
        ///Conversion factor between potential and absolute temperature of air
        ///</summary>
        public static float[][][] FAC = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Gradient Richardson number
        ///</summary>
        public static double[][][] RITSCH = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Sensible heat flux
        ///</summary>
        public static double[][] WQU = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Anthropogenic heat flux
        ///</summary>
        public static float[][] AWQ = CreateArray<float[]>(1, () => new float[1]);
        ///<summary>
        ///Absolute soil temperature at time t
        ///</summary>
        public static double[][][] TB = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Absolute soil temperature at time t+1
        ///</summary>
        public static double[][][] TBN = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Absolute soil temperature at time t-1
        ///</summary>
        public static double[][][] TBA = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Basic state of potential temperature at the simulation start
        ///</summary>
        public static double[][][] TBZ = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        /// Vertical cell heights of the soil model
        ///</summary>
        public static float[] DWB = new float[1];
        ///<summary>
        /// Vertical distance between neighbouring cells of the soil model
        ///</summary>
        public static float[] DZB = new float[1];
        ///<summary>
        ///Soil density
        ///</summary>
        public static float[][] RHOB = CreateArray<float[]>(1, () => new float[1]);
        ///<summary>
        ///Surface albedo
        ///</summary>
        public static float[][] ALAMBDA = CreateArray<float[]>(1, () => new float[1]);
        ///<summary>
        ///Damping parameter used at the model boundaries
        ///</summary>
        public static double[] ALPHA = new double[1];
        ///<summary>
        /// Temporary array used in the TDMA solver for passive scalars
        ///</summary>
        public static double[][][] A_PS = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        /// Temporary array used in the TDMA solver for passive scalars
        ///</summary>
        public static double[][][] B_PS = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        /// Temporary array used in the TDMA solver for passive scalars
        ///</summary>
        public static double[][][] C_PS = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        /// Temporary array used in the TDMA solver for passive scalars
        ///</summary>
        public static float[][][] AWEST_PS = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Temporary array used in the TDMA solver for passive scalars
        ///</summary>
        public static float[][][] ASOUTH_PS = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        /// Temporary array used in the TDMA solver for passive scalars
        ///</summary>
        public static float[][][] AEAST_PS = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Temporary array used in the TDMA solver for passive scalars
        ///</summary>
        public static float[][][] ANORTH_PS = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Temporary array used in the TDMA solver for passive scalars
        ///</summary>
        public static float[][][] AP0_PS = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Temporary array used in the TDMA solver for passive scalars
        ///</summary>
        public static double[][][] SU_PS = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Temporary array used in the TDMA solver for velocities
        ///</summary>
        public static float[][][] AIM = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Temporary array used in the TDMA solver for velocities
        ///</summary>
        public static float[][][] BIM = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Temporary array used in the TDMA solver for velocities
        ///</summary>
        public static float[][][] CIM = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Temporary array used in the TDMA solver for velocities
        ///</summary>
        public static float[][][] AW1 = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Temporary array used in the TDMA solver for velocities
        ///</summary>
        public static float[][][] AS1 = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Temporary array used in the TDMA solver for velocities
        ///</summary>
        public static float[][][] AE2 = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Temporary array used in the TDMA solver for velocities
        ///</summary>
        public static float[][][] AN2 = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        /// Temporary array used in the TDMA solver for velocities
        ///</summary>
        public static double[][][] SU = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Temporary array used in the TDMA solver for velocities
        ///</summary>
        public static float[][][] AP0 = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Temporary array for UIMP.VIMP, WIMP
        ///</summary>
        public static float[][][] F1U = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Temporary array for UIMP.VIMP, WIMP
        ///</summary>
        public static float[][][] F2U = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Temporary array for UIMP.VIMP, WIMP
        ///</summary>
        public static float[][][] F1V = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Temporary array for UIMP.VIMP, WIMP
        ///</summary>
        public static float[][][] F2V = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Temporary array for UIMP.VIMP, WIMP
        ///</summary>
        public static float[][][] F1W = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        ///Temporary array for UIMP.VIMP, WIMP
        ///</summary>
        public static float[][][] F2W = CreateArray<float[][]>(1, () => CreateArray<float[]>(1, () => new float[1]));
        ///<summary>
        /// Temporary array used in the TDMA solver for velocities
        ///</summary>
        public static double[][][] DIMU = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        /// Temporary array used in the TDMA solver for velocities
        ///</summary>
        public static double[][][] DIMV = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        /// Temporary array used in the TDMA solver for velocities
        ///</summary>
        public static double[][][] DIMW = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Temporary array used in the TDMA solver for non-hydrostatic pressure
        ///</summary>
        public static double[][][] AP = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Temporary array used in the TDMA solver for non-hydrostatic pressure
        ///</summary>
        public static double[][][] AB = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Temporary array used in the TDMA solver for non-hydrostatic pressure
        ///</summary>
        public static double[][][] AT = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Temporary array used in the TDMA solver for non-hydrostatic pressure
        ///</summary>
        public static double[][][] AW = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Temporary array used in the TDMA solver for non-hydrostatic pressure
        ///</summary>
        public static double[][][] AE = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Temporary array used in the TDMA solver for non-hydrostatic pressure
        ///</summary>
        public static double[][][] AS = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Temporary array used in the TDMA solver for non-hydrostatic pressure
        ///</summary>
        public static double[][][] AN = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Temporary array used for online output of prognostic variables
        ///</summary>
        public static float[][] VDATA1 = CreateArray<float[]>(1, () => new float[1]);
        ///<summary>
        ///Temporary array used for online output of prognostic variables
        ///</summary>
        public static float[][] VDATA2 = CreateArray<float[]>(1, () => new float[1]);
        ///<summary>
        ///Temporary array used for online output of prognostic variables
        ///</summary>
        public static float[][] VDATA3 = CreateArray<float[]>(1, () => new float[1]);
        ///<summary>
        ///Temporary array used for online output of prognostic variables
        ///</summary>
        public static float[][] VDATA4 = CreateArray<float[]>(1, () => new float[1]);
        ///<summary>
        ///Temporary array used for online output of prognostic variables
        ///</summary>
        public static float[][] VDATA8 = CreateArray<float[]>(1, () => new float[1]);
        ///<summary>
        /// Temporary array used for online output of prognostic variables
        ///</summary>
        public static float[] VDATA5 = new float[1];
        ///<summary>
        /// Temporary array used for online output of prognostic variables
        ///</summary>
        public static float[] VDATA6 = new float[1];
        ///<summary>
        /// Temporary array used for online output of prognostic variables
        ///</summary>
        public static float[] VDATA7 = new float[1];
        ///<summary>
        /// Temporary array used for online output of prognostic variables
        ///</summary>
        public static float[] VDATA9 = new float[1];
        ///<summary>
        ///Control variable for the temporal development of the entire mass divergence
        ///</summary>
        public static double[] MASSOURCE = new double[1];
        ///<summary>
        ///Control variable for checking flow convergence in the lowest layer
        ///</summary>
        public static float[][] U_TEMP = CreateArray<float[]>(1, () => new float[1]);
        ///<summary>
        ///Control variable for checking flow convergence in the lowest layer
        ///</summary>
        public static float[][] V_TEMP = CreateArray<float[]>(1, () => new float[1]);
        ///<summary>
        ///Control variable for checking flow convergence in the lowest layer
        ///</summary>
        public static float[][] W_TEMP = CreateArray<float[]>(1, () => new float[1]);
        ///<summary>
        ///Large scale value of u-wind at eastern boundary at time t
        ///</summary>
        public static double[][] USE1 = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of v-wind at eastern boundary at time t
        ///</summary>
        public static double[][] VSE = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of w-wind at eastern boundary at time t
        ///</summary>
        public static double[][] WSE = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of pot. temp. at eastern boundary at time t
        ///</summary>
        public static double[][] TSE = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of spec. humidity at eastern boundary at time t
        ///</summary>
        public static double[][] QUSE = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of u-wind at eastern boundary at time t+1
        ///</summary>
        public static double[][] USEN = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of v-wind at eastern boundary at time t+1
        ///</summary>
        public static double[][] VSEN = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of w-wind at eastern boundary at time t+1
        ///</summary>
        public static double[][] WSEN = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of pot. temp. at eastern boundary at time t+1
        ///</summary>
        public static double[][] TSEN = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        /// Large scale value of spec. humidity at eastern boundary at time t+1
        ///</summary>
        public static double[][] QUSEN = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of u-wind at western boundary at time t
        ///</summary>
        public static double[][] USW = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of v-wind at western boundary at time t
        ///</summary>
        public static double[][] VSW = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of w-wind at western boundary at time t
        ///</summary>
        public static double[][] WSW = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of pot. temp. at western boundary at time t
        ///</summary>
        public static double[][] TSW = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of spec. humidity at western boundary at time t
        ///</summary>
        public static double[][] QUSW = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of u-wind at western boundary at time t+1
        ///</summary>
        public static double[][] USWN = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of v-wind at western boundary at time t+1
        ///</summary>
        public static double[][] VSWN = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of w-wind at western boundary at time t+1
        ///</summary>
        public static double[][] WSWN = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of pot. temp. at western boundary at time t+1
        ///</summary>
        public static double[][] TSWN = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        /// Large scale value of spec. humidity at western boundary at time t+1
        ///</summary>
        public static double[][] QUSWN = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of u-wind at northern boundary at time t
        ///</summary>
        public static double[][] USN = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of v-wind at northern boundary at time t
        ///</summary>
        public static double[][] VSN = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of w-wind at northern boundary at time t
        ///</summary>
        public static double[][] WSN = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of pot. temp. at northern boundary at time t
        ///</summary>
        public static double[][] TSN = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of spec. humidity at northern boundary at time t
        ///</summary>
        public static double[][] QUSN = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of u-wind at northern boundary at time t+1
        ///</summary>
        public static double[][] USNN = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of v-wind at northern boundary at time t+1
        ///</summary>
        public static double[][] VSNN = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of w-wind at northern boundary at time t+1
        ///</summary>
        public static double[][] WSNN = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of pot. temp. at northern boundary at time t+1
        ///</summary>
        public static double[][] TSNN = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        /// Large scale value of spec. humidity at northern boundary at time t+1
        ///</summary>
        public static double[][] QUSNN = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of u-wind at southern boundary at time t
        ///</summary>
        public static double[][] USS = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of v-wind at southern boundary at time t
        ///</summary>
        public static double[][] VSS = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of w-wind at southern boundary at time t
        ///</summary>
        public static double[][] WSS = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of pot. temp. at southern boundary at time t
        ///</summary>
        public static double[][] TSS = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of spec. humidity at southern boundary at time t
        ///</summary>
        public static double[][] QUSS = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of u-wind at southern boundary at time t+1
        ///</summary>
        public static double[][] USSN = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of v-wind at southern boundary at time t+1
        ///</summary>
        public static double[][] VSSN = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of w-wind at southern boundary at time t+1
        ///</summary>
        public static double[][] WSSN = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        ///Large scale value of pot. temp. at southern boundary at time t+1
        ///</summary>
        public static double[][] TSSN = CreateArray<double[]>(1, () => new double[1]);
        ///<summary>
        /// Large scale value of spec. humidity at southern boundary at time t+1
        ///</summary>
        public static double[][] QUSSN = CreateArray<double[]>(1, () => new double[1]);
        //public static double[][][][] Wrad = CreateArray<double[][][]>(1, () => CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1])));
        ///<summary>
        ///Temporary Radiation in TIMP
        ///</summary>
        public static double[][][] RADIATION = CreateArray<double[][]>(1, () => CreateArray<double[]>(1, () => new double[1]));
        ///<summary>
        ///Stability class (1-7)
        ///</summary>
        public static Int16[][] stabilityclass = CreateArray<Int16[]>(1, () => new Int16[1]);
        ///<summary>
        ///Factor to reduce the solar radiation
        ///</summary>
        public static float Solar_reduction_factor = 1.0F;

        ///<summary>
        ///Starting date of computation
        ///</summary>
        public static Int32 NDATUM;
        ///<summary>
        /// Starting hour of computation
        ///</summary>
        public static Int32 ISTUNDE;
        ///<summary>
        /// Maximum allowed time step
        ///</summary>
        public static double DT;
        ///<summary>
        ///Total modelling time
        ///</summary>
        public static double TLIMIT2;
        ///<summary>
        /// Write intermediate output files after IOUTPUT seconds
        ///</summary>
        public static Int32 IOUTPUT;
        ///<summary>
        ///Historic value not used currently
        ///</summary>
        public static double DELW;
        ///<summary>
        /// Initial relative humidity
        ///</summary>
        public static double QUINIT;
        ///<summary>
        ///Height of the lowest compuation height - not used anymore
        ///</summary>
        public static double ZSEEH;
        ///<summary>
        ///Temperature at 2m height
        ///</summary>
        public static double TINIT;
        ///<summary>
        ///Gradient of potential temperature - not used anymore
        ///</summary>
        public static double TGRAD;
        ///<summary>
        ///Height of the neutral boundary layer - not used anymore
        ///</summary>
        public static double ZNEUT;
        ///<summary>
        /// Surface temperature
        ///</summary>
        public static double TBINIT;
        ///<summary>
        ///Soil temperature in 1m depth -> Note that TBINIT1, TBINIT and TINIT are not used by their absolute values but by their differences to define the initial state of the atmospheric and soil temperatures
        ///</summary>
        public static double TBINIT1;
        ///<summary>
        ///Latitude
        ///</summary>
        public static double BGRAD;
        ///<summary>
        ///Number of timesteps after which the radiation is updated
        ///</summary>
        public static Int32 IRAD;
        ///<summary>
        ///Debug level 0 none, 3 highest
        ///</summary>
        public static Int32 IDEBUG;
        ///<summary>
        /// Compute u v w pn t gw fo
        ///</summary>
        public static Int32 IFLAGS1;
        ///<summary>
        /// Compute br pr qu psi te tb str
        ///</summary>
        public static Int32 IFLAGS2;
        ///<summary>
        ///Flag to determine whether initale state of the model is computed using the diagnostic model
        ///</summary>
        public static double DIA;
        ///<summary>
        ///Not used anymore - former flag to switch between explicit and implicit solutaion algorithm
        ///</summary>
        public static double DIAI;
        ///<summary>
        /// Underrelaxation factor for momentum equations
        ///</summary>
        public static double RELAXV;
        ///<summary>
        /// Underrelaxation factor for passive scalar equations
        ///</summary>
        public static double RELAXT;
        ///<summary>
        ///Forcing of catabatic flows with -40/-25/-15 W/m² AKLA 7/6/5
        ///</summary>
        public static Int32 ICATAFORC;
        ///<summary>
        ///Defines the interval after which the boundary conditions are updated with certain ERA5 data
        ///</summary>
        public static float ERA5BOUND;
        ///<summary>
        ///Counter for new boundary conditions with certain ERA5 data
        ///</summary>
        public static float ERA5BOUND_Threshold;
        ///<summary>
        ///BOUNDARY CONDITION (1,6,10,11,12)
        ///</summary>
        public static Int32 IBOUND;
        ///<summary>
        ///NON-STEADYSTATE/FORCING: Default = 0, "METEOPGT.ALL" = 1, GRAMM transient = 2
        ///</summary>
        public static Int32 ISTAT;
        ///<summary>
        ///NON-STEADYSTATE/FORCING W. IFU-MM5 DATA YES=1 -> not yet implemented in this version
        ///</summary>
        public static Int32 MM5IFU;
        ///<summary>
        ///GRAMM IN GRAMM NESTING YES=1 -> not yet implemented in this version
        ///</summary>
        public static Int32 NESTGRAMM;
        ///<summary>
        ///flow out as binary-stream, little_endian, with header 0, or 3 without header
        ///</summary>
        public static Int32 IOUT;
        ///<summary>
        ///Total simulation time
        ///</summary>
        public static double DTI;
        ///<summary>
        ///Maximum allowed time step
        ///</summary>
        public static double DTMAX;
        ///<summary>
        ///Maximum number of processors/threads used in each parallelized procedure
        ///</summary>
        public static Int32 IPROC;
        ///<summary>
        ///Year of the simulation (only for transient simulations)
        ///</summary>
        public static Int32 IJAHR;
        ///<summary>
        ///Year of the simulation (only for transient simulations)
        ///</summary>
        public static Int32 IJAHR4digits;
        ///<summary>
        ///Month of the simulation start (only for transient simulations)
        ///</summary>
        public static Int32 IMON;
        ///<summary>
        ///Day of the simulation start (only for transient simulations)
        ///</summary>
        public static Int32 ITAG;
        ///<summary>
        ///Hour of the simulation start (only for transient simulations)
        ///</summary>
        public static Int32 ISTU;
        ///<summary>
        ///Minute of the simulation start (only for transient simulations)
        ///</summary>
        public static Int32 IMIN;
        ///<summary>
        ///Western border of model domain
        ///</summary>
        public static Int32 IKOOA;
        ///<summary>
        ///Southern border of model domain
        ///</summary>
        public static Int32 JKOOA;
        ///<summary>
        ///Minimum surface elevation
        ///</summary>
        public static float AHMIN;
        ///<summary>
        ///Maximum surface elevation
        ///</summary>
        public static float AHMAX;
        ///<summary>
        /// Inversion height
        ///</summary>
        public static float Inversion_Height;
        ///<summary>
        ///Velocity
        ///</summary>
        public static float Wind_Velocity;
        ///<summary>
        ///Cell indices of minimum surface elevation
        ///</summary>
        public static Int32 AHMINI;
        ///<summary>
        ///Cell indices of minimum surface elevation
        ///</summary>
        public static Int32 AHMINJ;
        ///<summary>
        ///Cell indices of minimum surface elevation within calculation domain
        ///</summary>
        public static Int32 AHMINI_D;
        ///<summary>
        ///Cell indices of minimum surface elevation within calculation domain
        ///</summary>
        public static Int32 AHMINJ_D;
        ///<summary>
        ///Index in x-direction of receptor points
        ///</summary>
        public static List<int> inrec = new List<int>();
        ///<summary>
        ///Index in y-direction of receptor points
        ///</summary>
        public static List<int> jnrec = new List<int>();
        ///<summary>
        ///Index in z-direction of receptor points
        ///</summary>
        public static List<int> knrec = new List<int>();
        ///<summary>
        ///x-coordinate of receptor points
        ///</summary>
        public static List<double> Xrec = new List<double>();
        ///<summary>
        ///y-coordinate of receptor points
        ///</summary>
        public static List<double> Yrec = new List<double>();
        ///<summary>
        ///z-coordinate of receptor points
        ///</summary>
        public static List<double> Zrec = new List<double>();
        ///<summary>
        ///u-component at receptor points
        ///</summary>
        public static List<double> Urec = new List<double>();
        ///<summary>
        ///v-component at receptor points
        ///</summary>
        public static List<double> Vrec = new List<double>();
        ///<summary>
        ///temperature at receptor points
        ///</summary>
        public static List<double> Trec = new List<double>();
        ///<summary>
        ///global radiation at receptor points
        ///</summary>
        public static List<double> Globradrec = new List<double>();
        ///<summary>
        ///longwave outgoing radiation at receptor points
        ///</summary>
        public static List<double> Longradrec = new List<double>();
        ///<summary>
        ///soil heat flux at receptor points
        ///</summary>
        public static List<double> Soilheatfluxrec = new List<double>();
        ///<summary>
        ///sensible heat flux at receptor points
        ///</summary>
        public static List<double> Sensheatfluxrec = new List<double>();
        ///<summary>
        ///latent heat flux at receptor points
        ///</summary>
        public static List<double> Latheatfluxrec = new List<double>();
        ///<summary>
        ///Stretching-factor in the vertical
        ///</summary>
        public static double STRETCH;
        ///<summary>
        /// Coriolis frequency normal to the earths surface
        ///</summary>
        public static double FN;
        ///<summary>
        /// Coriolis frequency parallel to the earths surface
        ///</summary>
        public static double FH;
        ///<summary>
        ///saturation pressure
        ///</summary>
        public static double[] PSAT = new double[1];
        ///<summary>
        ///Gravitational acceleration
        ///</summary>
        public static double GERD;
        ///<summary>
        /// General gas constant
        ///</summary>
        public static double GASCON;
        ///<summary>
        /// Heat capacity of air for constant pressure
        ///</summary>
        public static double CPLUFT;
        ///<summary>
        ///Evaporation heat of water
        ///</summary>
        public static double ALW;
        ///<summary>
        ///Switch for radiation model: 1=thin clouds 2=thick clouds
        ///</summary>
        public static Int16 ISOL;
        ///<summary>
        ///Number of the actual calculated flow field situation
        ///</summary>
        public static int IWETTER;
        ///<summary>
        ///Character determining whether GRAMM is driven by the file meteopgt.all or not
        ///</summary>
        public static string METEO;
        ///<summary>
        ///Anemometer heat when using meteopgt.all as input file
        ///</summary>
        public static double ANEMO;
        ///<summary>
        ///Stability class when using meteopgt.all as input file
        ///</summary>
        public static Int16 AKLA;
        ///<summary>
        ///Wind speed as stored in meteopgt.all used for determining the correct global radiation
        ///</summary>
        public static double Windspeed_meteopgt;
        ///<summary>
        /// Modified total simulation time in dependence on the stability class
        ///</summary>
        public static double TLIMIT;
        ///<summary>
        ///Initial Obukhov length
        ///</summary>
        public static double Obini;
        ///<summary>
        /// Initial surface-roughness length
        ///</summary>
        public static double Rauigkeit;
        ///<summary>
        /// van Karman constant
        ///</summary>
        public static double CK;
        ///<summary>
        ///minimum turbulent viscosity
        ///</summary>
        public static double VISEL;
        ///<summary>
        ///heat capacity of soil
        ///</summary>
        public static double CPBOD;
        ///<summary>
        ///Stefan Bolzmann constant
        ///</summary>
        public static double SIGMA;
        ///<summary>
        ///temperature value used to improve numerical accuracy of the solution algorithm for temperature
        ///</summary>
        public static double TBZ1;
        ///<summary>
        ///Flag used to check whether the file meteopgt.all is used as input file
        ///</summary>
        public static Int16 IPGT;
        ///<summary>
        /// specific modelling time format used for the radiation model
        ///</summary>
        public static double TJETZT;
        ///<summary>
        ///integrated modelling time
        ///</summary>
        public static double REALTIME;
        ///<summary>
        ///total time steps
        ///</summary>
        public static Int16 ITIME;
        ///<summary>
        ///number of pressure-iterations
        ///</summary>
        public static Int16 INUMS;
        ///<summary>
        ///total mass divergence
        ///</summary>
        public static double DIVSUM, SUMG;
        ///<summary>
        ///counter for computing the trend of the mass-divergence used to adjust the time-step
        ///</summary>
        public static Int16 IDIV;
        ///<summary>
        ///temporal trend of the total mass-divergence used to compute the actual time-step
        ///</summary>
        public static double STEIGUNG;
        ///<summary>
        ///initial total mass-divergence at the beginning of each computation
        ///</summary>
        public static double MASSOURCE_FIRST;
        ///<summary>
        ///turbulent Prandtl-number
        ///</summary>
        public static double PRTE;

        ///<summary>
        ///Flag switching the computation of the u-component on/off
        ///</summary>
        public static Boolean ICU;
        ///<summary>
        ///Flag switching the computation of the v-component on/off
        ///</summary>
        public static Boolean ICV;
        ///<summary>
        ///Flag switching the computation of the w-component on/off
        ///</summary>
        public static Boolean ICW;
        ///<summary>
        ///Flag switching the computation of the non-hydrostatic pressure on/off
        ///</summary>
        public static Boolean ICPN;
        ///<summary>
        ///Flag switching the computation of the pot. temperature on/off
        ///</summary>
        public static Boolean ICT;
        ///<summary>
        ///Flag switching the computation of the hydrostatic pressure on/off (not used)
        ///</summary>
        public static Boolean ICPH;
        ///<summary>
        ///Flag switching the computation of the Fourier-Transformed upper boundary conditions on/off (not used)
        ///</summary>
        public static Boolean IFOU;
        ///<summary>
        ///Flag switching the computation of the spec. humidity on/off
        ///</summary>
        public static Boolean ICQU;
        ///<summary>
        /// Flag switching the computation of a passive scalar on/off (not used)
        ///</summary>
        public static Boolean ICPSI;
        ///<summary>
        ///Flag switching the computation of the turbulent kinetic energy on/off
        ///</summary>
        public static Boolean ICTE;
        ///<summary>
        /// Flag switching the computation of the radiation model on/off
        ///</summary>
        public static Boolean ICSTR;
        ///<summary>
        ///Flag switching the computation of the surface layer model on/off
        ///</summary>
        public static Boolean ICPR;
        ///<summary>
        ///Flag switching the computation of the boundary conditions on/off
        ///</summary>
        public static Boolean ICBR;
        ///<summary>
        ///Flag switching the computation of the surface temperature on/off
        ///</summary>
        public static Boolean ICTB;
        ///<summary>
        ///Flag switching the computation of the geostrophic winds from the large scale model when run in nesting mode on/off (not used)
        ///</summary>
        public static Boolean ICGW;
        ///<summary>
        /// Flag set when receptor points are set (file receptor_GRAMM.dat)
        ///</summary>
        public static Boolean recexist;
        ///<summary>
        ///Flag if Steady-State criterion should be written to a file
        ///</summary>
        public static Boolean WriteSteadyState = false;

        ///<summary>
        ///Counter for writing "DispGRAMM.txt"
        ///</summary>
        public static int Counter = 0;
        ///<summary>
        ///Counter for Terminal Output
        ///</summary>
        public static int TerminalOut = 0;
        ///<summary>
        ///Threshold for terminal output
        ///</summary>
        public static int TerminalThreshold = 1;
        ///<summary>
        ///Number of cells used at the lateral boundaries for smoothing the orography
        ///</summary>
        public static int nr_cell_smooth = 0;

        ///<summary>
        /// Write online data
        ///</summary>

        public static bool GRAMM_Online_flag = true;
        ///<summary>
        /// Running in linux?
        ///</summary>
        public static bool unix = false;
        ///<summary>
        /// flag, if a computation is repeated caused by numerical instabilities
        ///</summary>
        public static int computation_retry = 0;
        ///<summary>
        ///original number of weather situations in the file meteopgt.all
        ///</summary>
        public static int meteopgt_nr = 0;
        ///<summary>
        ///the original maximum time step
        ///</summary>
        public static double max_time_step_original = 0;

        ///<summary>
        ///5.4.2017 Ku Massource queue
        ///</summary>
        public static Queue<double> MASSOURCE_Queue = new Queue<double>();
        ///<summary>
        ///5.4.2017 Ku IDIV up Counter
        ///</summary>
        public static int IDIV_Up;
        ///<summary>
        /// 5.4.2017 Ku IDIV Lock to increase time step
        ///</summary>
        public static int IDIV_LockDown;
        ///<summary>
        ///5.4.2017 Ku IDIV Lock to increase time step
        ///</summary>
        public static int IDIV_LockUp;
        ///<summary>
        ///5.4.2017 Ku IDIV Lock to increase time step
        ///</summary>
        public static int IDIV_LockDown2;
        ///<summary>
        ///5.4.2017 Ku IDIV Lock to avoit too low relax factors
        ///</summary>
        public static int IDIV_LockRelax;
        ///<summary>
        /// 6.4.2017 Ku IDIV Lock to avoit PingPong effects
        ///</summary>
        public static int IDIV_PingPong;
        ///<summary>
        ///5.4.2017 Ku MASSOURCE actual
        ///</summary>
        public static double MASSOURCE_Act;
        ///<summary>
        ///5.4.2017 Ku MASSOURCE old
        ///</summary>
        public static double MASSOURCE_Old;
        ///<summary>
        ///5.4.2017 ÖT MASSOURCE previous time step
        ///</summary>
        public static double MASSOURCE_minusone;
        ///<summary>
        ///11.4.2017 Ku first weather situation from console
        ///</summary>
        private static int IWetter_Console_First = 0;
        ///<summary>
        /// 11.4.2017 Ku last weather situation from console
        ///</summary>
        private static int IWetter_Console_Last = 9999999;
        ///<summary>
        ///25.4.2017 Ku Relaxv from console
        ///</summary>
        private static double Relaxv_Console = 0;
        ///<summary>
        ///25.4.2017 Ku Relaxt from console
        ///</summary>
        private static double Relaxt_Console = 0;
        ///<summary>
        /// 25.4.2017 Ku MaxTimeStep from console
        ///</summary>
        private static double MaxTimeStep_Console = 0;


        ///<summary>
        ///large-scale forcing using meteopgt.all
        ///</summary>
        public static double WU1 = 0;

        ///<summary>
        ///large-scale forcing using meteopgt.all
        ///</summary>
        public static double WU2 = 0;

        ///<summary>
        ///large-scale forcing using meteopgt.all
        ///</summary>
        public static double WV1 = 0;

        ///<summary>
        ///large-scale forcing using meteopgt.all
        ///</summary>
        public static double WV2 = 0;

        ///<summary>
        /// Original RelaxV value
        ///</summary>
        public static double Relaxv_ori;
        ///<summary>
        /// Original RelaxT value
        ///</summary>
        public static double Relaxt_ori;
        ///<summary>
        /// Min divergence
        ///</summary>
        public static double Divergence_Min = 10e9;
        ///<summary>
        /// Western border of GRAMM domain
        ///</summary>
        public static float GRAMM_West = 0;
        ///<summary>
        /// Southern border of GRAMM domain
        ///</summary>
        public static float GRAMM_South = 0;

        ///<summary>
        ///UTC date used for interpolation ERA5 data to force/initialize GRAMM
        ///</summary>
        public static DateTime dateUTC = new DateTime();
        ///<summary>
        ///first date used for interpolation ERA5 data to force/initialize GRAMM
        ///</summary>
        public static DateTime ERA5_date1 = new DateTime();
        ///<summary>
        ///second date used for interpolation ERA5 data to force/initialize GRAMM
        ///</summary>
        public static DateTime ERA5_date2 = new DateTime();
        ///<summary>
        /// avearge longitude of GRAMM domain
        ///</summary>
        public static float Longitude = 47;

        ///<summary>
        ///time interval for intermediate GRAMM flow field output
        ///</summary>
        public static int Intermed_Threshold;
        ///<summary>
        ///Defines the time interval in seconds after which GRAMM is completely re-initialized with ERA5 data
        ///</summary>
        public static float REINITIALIZATION;
        ///<summary>
        ///Counter for re-initialization
        ///</summary>
        public static float REINITIALIZATION_Threshold;

        ///<summary>
        ///flag determining whether chemistry is computed or not
        ///</summary>
        public static bool chemistry = false;
        ///<summary>
        ///chemical mechanism
        ///</summary>
        public static string chemistry_mechanism;
        ///<summary>
        ///number of chemical species for which advection and diffusion has to be computed
        ///</summary>
        public static int NSPEZ = 20;
        ///<summary>
        ///time interval, after which the chemical solver is called
        ///</summary>
        public static float Update_Chemistry = 60;
        ///<summary>
        ///counter for calling chemical solver
        ///</summary>
        public static float Update_Chemistry_Threshold = 60;
        ///<summary>
        ///Calculation of solar radiation
        ///</summary>
        public static RadiationCalculation RadiationModel;
    }
}
