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

using System.Runtime.CompilerServices;
using System;

namespace GRAMM_2001
{
    /// <summary>
    /// Fast array copying using MemoryCopy
    /// </summary>
    public class FastCopy
    {
        /// <summary>
        /// Copy an array using MemoryCopy and lock the source array
        /// </summary>
        /// <param name="source">Source</param>
        /// <param name="destination">Destination</param>
        /// <param name="length">Number of elements</param>
        [MethodImpl(MethodImplOptions.AggressiveOptimization | MethodImplOptions.AggressiveInlining)]
        public static unsafe void CopyArrayLockSource(double[] source, double[] destination)
        {
            lock (source.SyncRoot)
            {
                fixed (double* srcPtr = source, destPtr = destination)
                {
                    long length = source.LongLength * sizeof(double);
                    Buffer.MemoryCopy(srcPtr, destPtr, length, length);
                }
            }
        }

        /// <summary>
        /// Copy an array using MemoryCopy
        /// </summary>
        /// <param name="source">Source</param>
        /// <param name="destination">Destination</param>
        /// <param name="length">Number of elements</param>
        [MethodImpl(MethodImplOptions.AggressiveOptimization | MethodImplOptions.AggressiveInlining)]
        public static void CopyArray(double[] source, double[] destination, int length)
        {
            for (int i = 0; i < length; i++)
            {
                destination[i] = source[i];
            }
        }
        /// <summary>
        /// Copy an array using the source length
        /// </summary>
        /// <param name="source">Source</param>
        /// <param name="destination">Destination</param>
        /// <param name="length">Number of elements</param>
        [MethodImpl(MethodImplOptions.AggressiveOptimization | MethodImplOptions.AggressiveInlining)]
        public static unsafe void CopyArraySourceLen(double[] source, double[] destination)
        {
            fixed (double* srcPtr = source, destPtr = destination)
            {
                long length = source.LongLength * sizeof(double);
                Buffer.MemoryCopy(srcPtr, destPtr, length, length);
            }
        }
        /// <summary>
        /// Copy an array using the destination length
        /// </summary>
        /// <param name="source">Source</param>
        /// <param name="destination">Destination</param>
        /// <param name="length">Number of elements</param>
        [MethodImpl(MethodImplOptions.AggressiveOptimization | MethodImplOptions.AggressiveInlining)]
        public static void CopyArrayDestLen(double[] source, double[] destination)
        {
            int i = 0;
            for (; i < destination.Length; i++)
            {
                destination[i] = source[i];
            }
        }

        /// <summary>
        /// Copy an array using MemoryCopy and the lenght of the destination array
        /// </summary>
        /// <param name="source">Source</param>
        /// <param name="destination">Destination</param>
        /// <param name="length">Number of elements</param>
        [MethodImpl(MethodImplOptions.AggressiveOptimization | MethodImplOptions.AggressiveInlining)]
        public static unsafe void CopyArrayLockDest(double[] source, double[] destination)
        {
            lock (destination.SyncRoot)
            {
                fixed (double* srcPtr = source, destPtr = destination)
                {
                    long length = destination.LongLength * sizeof(double);
                    Buffer.MemoryCopy(srcPtr, destPtr, length, length);
                }        
            }
        }
    }
}
