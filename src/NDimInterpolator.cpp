/*******************************************************************************
 * This file is part of pyNFrac
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * pyNFrac is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * pyNFrac is distributed in the hope that it will be useful,
 * but WITOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with pyNFrac. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

/**
 * @file NDimInterpolator.cpp
 *
 * @brief NDimInterpolator implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "NDimInterpolator.hpp"

/**
 * @brief N-dimensional linear interpolator.
 *
 * @param dataTable N-dimensional double precision floating point array on which
 * to interpolate.
 * @param paramArrays N 1-dimensional data arrays corresponding to the parameter
 * values.
 * @param tableSize N lenghts of the parameter data arrays.
 * @param values N data values marking the parameter values for which we want
 * the interpolated value.
 * @param ndim Number of dimensions N.
 * @param index Array containing an initial guess for the lower indices of the
 * parameter interval the value resides in (if -1, the interpolator will
 * determine these values). When the function returns, this array always
 * contains these indices.
 * @param boxValues Position of the interpolation point within the parameter
 * intervals (if -1., the interpolator will determine these values). When the
 * function returns, this array always contains these positions.
 * @param extrapolate If true, this function will extrapolate for parameter
 * values outside the given parameter ranges. If false, the values for the
 * outermost parameter value will be used.
 * @return Value at the requested point in the N-dimensional parameter space.
 */
double NDimInterpolator::interpolate(double *dataTable, double **paramArrays,
                                     unsigned int *tableSize, double *values,
                                     unsigned int ndim, int *index,
                                     double *boxValues, bool extrapolate) {

  /* hunt for location in dataTable (aka find lower indexes) and convert values
   * to range [0..1] within box around location
   * will not be done if index and boxValues are already supplied */
  int idim = ndim - 1;
  unsigned int outOfRange = 0;
  if (index[0] == -1) {
    while (idim + 1) {

      hunt(paramArrays[idim], tableSize[idim], values[idim], &index[idim]);

      /* if value out of range, do necessary stuff for limiting or extrapolating
       */
      /* value lower than lower limit */
      if (index[idim] == -1) {
        /* make sure databox is within table */
        index[idim] = 0;
        /* if extrapolating not desired, set to lower limit - otherwise flag
         * outOfRange */
        if (!extrapolate)
          values[idim] = paramArrays[idim][0];
        else
          outOfRange = 1;
      }
      /* value higher than upper limit */
      else if (index[idim] == tableSize[idim] - 1) {
        /* make sure databox is within table */
        index[idim] = tableSize[idim] - 2;
        /* if extrapolating not desired, set to upper limit - otherwise flag
         * outOfRange */
        if (!extrapolate)
          values[idim] = paramArrays[idim][tableSize[idim] - 1];
        else
          outOfRange = 1;
      }
      boxValues[idim] =
          (values[idim] - paramArrays[idim][index[idim]]) /
          (paramArrays[idim][index[idim] + 1] - paramArrays[idim][index[idim]]);
      idim--;
    }
  }
  /* extract "box" of data from table, surrounding the location */
  unsigned int i, j, position, factor;
  unsigned int boxLength = 1 << ndim; // 1<<ndim = pow(2,ndim)
  unsigned int *bits = new unsigned int[ndim];
  double *dataBox = new double[boxLength];
  for (i = 0; i < boxLength; i++) {
    convertToBits(i, bits, ndim);
    position = 0;
    for (idim = 0; idim < ndim; idim++) {
      factor = 1;
      for (j = idim + 1; j < ndim; j++)
        factor *= tableSize[j];
      position += (index[idim] + bits[ndim - 1 - idim]) * factor;
    }
    dataBox[i] = dataTable[position];
  }

  /* do interpolation = "collapse" box dimension per dimension onto first
   * element */
  int separation = 1, ipos;
  idim = ndim;
  while (idim) {
    for (ipos = 0; ipos < boxLength; ipos += 2 * separation) {
      /* extrapolate or interpolate */
      if (outOfRange)
        dataBox[ipos] =
            dataBox[ipos] +
            (dataBox[ipos + separation] - dataBox[ipos]) * boxValues[idim - 1];
      else
        dataBox[ipos] = dataBox[ipos] * (1 - boxValues[idim - 1]) +
                        dataBox[ipos + separation] * boxValues[idim - 1];
    }
    separation *= 2;
    idim--;
  }

  double result = dataBox[0];

  delete[] dataBox;
  delete[] bits;

  return result;
}

/**
 * @brief Find the lower index of the interval containing the given value in the
 * given array with given length.
 *
 * @param xx Parameter array.
 * @param n Length of the parameter array.
 * @param x Value to locate.
 * @param jlo Index of the lower bound of the interval containing the value. If
 * this value is initialized to a value in the range [0, n], it is used as an
 * initial guess.
 */
void NDimInterpolator::hunt(double *xx, int n, double x, int *jlo) {
  int jm, jhi, inc, nm1;
  int ascnd;
  nm1 = n - 1;
  ascnd =
      (xx[nm1] >= xx[0]); // True if ascending order of table, false otherwise.

  if (*jlo <= -1 || *jlo > nm1) {
    // Input guess not useful. Go immediately to bisection
    *jlo = -1;
    jhi = n;
  } else {
    // Set the hunting increment.
    inc = 1;

    // Hunt up:
    if ((x >= xx[*jlo]) == ascnd) {

      if (*jlo == nm1)
        return;

      jhi = (*jlo) + 1;

      while ((x >= xx[jhi]) == ascnd) { // Not done hunting,
        *jlo = jhi;
        inc += inc; // so double the increment
        jhi = (*jlo) + inc;
        if (jhi > nm1) { // Done hunting, since off end of table.
          jhi = n;
          break;
        } // Try again.
      }   // Done hunting, value bracketed.
    }

    // Hunt down:
    else {

      if (*jlo == 0) {
        *jlo = -1;
        return;
      }

      jhi = (*jlo)--;

      while ((x < xx[*jlo]) == ascnd) { // Not done hunting,
        jhi = (*jlo);
        inc <<= 1;        // so double the increment
        if (inc >= jhi) { // Done hunting, since off end of table.
          *jlo = -1;
          break;
        } else
          *jlo = jhi - inc;
      } // and try again.
    }   // Done hunting, value bracketed.
  }

  // Hunt is done, so begin the final bisection phase:
  while ((jhi - (*jlo)) != 1) {
    jm = (jhi + (*jlo)) >> 1;
    if ((x >= xx[jm]) == ascnd)
      *jlo = jm;
    else
      jhi = jm;
  }

  if (x == xx[nm1])
    *jlo = n - 2;
  if (x == xx[0])
    *jlo = 0;
}

/**
 * @brief Convert an unsigned integer to its binary representation.
 *
 * @param number Number to convert.
 * @param buffer Buffer to store the resulting bits in.
 * @param bufSize Size of the buffer.
 */
void NDimInterpolator::convertToBits(unsigned int number, unsigned int *buffer,
                                     unsigned int bufSize) {
  unsigned int i;
  for (i = 0; i < bufSize; number >>= 1)
    buffer[i++] = (number & 1) ? 1 : 0;
}
