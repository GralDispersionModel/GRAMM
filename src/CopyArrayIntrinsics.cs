#region Copyright
///<remarks>
/// <GRAMM Mesoscale Model>
/// Copyright (C) [2021]  [Markus Kuntner]
/// This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
/// the Free Software Foundation version 3 of the License
/// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
/// You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
///</remarks>
#endregion

using System;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;
using System.Runtime.CompilerServices;
using System.Threading;

namespace GRAMM_2001
{
    partial class Program
    {
        /// <summary>
        /// Copy an array using AVX intrinsics and lock the source array
        /// </summary>
        /// <param name="copyFrom">Source</param>
        /// <param name="copyTo">Destination</param>
        /// <param name="length">Number of elements</param>
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        public static void CopyArrayLockSource(double[] copyFrom, double[] copyTo)
        {
            lock (copyFrom.SyncRoot)
            {
                int i = 0;
                // if (Avx.IsSupported)
                // {
                //     fixed (double* source = copyFrom, dest = copyTo)
                //     {
                //         int len = copyFrom.Length - Vector256<double>.Count;
                //         for (; i < len; i += Vector256<double>.Count)
                //         {
                //             Avx.Store(dest + i, Avx.LoadVector256(source + i));
                //         }
                //     }
                // }
                for (; i < copyFrom.Length; i++)
                {
                    copyTo[i] = copyFrom[i];
                }
            }
        }

        /// <summary>
        /// Copy an array using AVX intrinsics 
        /// </summary>
        /// <param name="copyFrom">Source</param>
        /// <param name="copyTo">Destination</param>
        /// <param name="length">Number of elements</param>
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        public static void CopyArray(double[] copyFrom, double[] copyTo, int length)
        {
            int i = 0;
            // if (Avx.IsSupported)
            // {
            //     fixed (double* source = copyFrom, dest = copyTo)
            //     {
            //         int len = length - Vector256<double>.Count;
            //         for (; i < len; i += Vector256<double>.Count)
            //         {
            //             Avx.Store(dest + i, Avx.LoadVector256(source + i));
            //         }
            //     }
            // }
            for (; i < length; i++)
            {
                copyTo[i] = copyFrom[i];
            }
        }
        /// <summary>
        /// Copy an array using the source length
        /// </summary>
        /// <param name="copyFrom">Source</param>
        /// <param name="copyTo">Destination</param>
        /// <param name="length">Number of elements</param>
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        public static void CopyArraySourceLen(double[] copyFrom, double[] copyTo)
        {
            int i = 0;
            for (; i < copyFrom.Length; i++)
            {
                copyTo[i] = copyFrom[i];
            }
        }
        /// <summary>
        /// Copy an array using the destination length
        /// </summary>
        /// <param name="copyFrom">Source</param>
        /// <param name="copyTo">Destination</param>
        /// <param name="length">Number of elements</param>
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        public static void CopyArrayDestLen(double[] copyFrom, double[] copyTo)
        {
            int i = 0;
            for (; i < copyTo.Length; i++)
            {
                copyTo[i] = copyFrom[i];
            }
        }

        /// <summary>
        /// Copy an array using AVX intrinsics and the lenght of the destination array
        /// </summary>
        /// <param name="copyFrom">Source</param>
        /// <param name="copyTo">Destination</param>
        /// <param name="length">Number of elements</param>
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        public static void CopyArrayLockDest(double[] copyFrom, double[] copyTo)
        {
            lock (copyTo.SyncRoot)
            {
                int i = 0;
                // if (Avx.IsSupported)
                // {
                //     fixed (double* source = copyFrom, dest = copyTo)
                //     {
                //         int len = length - Vector256<double>.Count;
                //         for (; i < len; i += Vector256<double>.Count)
                //         {
                //             Avx.Store(dest + i, Avx.LoadVector256(source + i));
                //         }
                //     }
                // }
                for (; i < copyTo.Length; i++)
                {
                    copyTo[i] = copyFrom[i];
                }
            }
        }
    }
}
