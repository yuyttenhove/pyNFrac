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
 * @file NDimInterpolator.hpp
 *
 * @brief N-dimensional interpolator.
 *
 * This interpolator was based on the ND interpolator written by Joeri Schroyen
 * for Hyplot: https://sourceforge.net/projects/hyplot/.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef NDIMINTERPOLATOR_HPP
#define NDIMINTERPOLATOR_HPP

/**
 * @brief N-dimensional interpolator.
 */
class NDimInterpolator {
public:
  static double interpolate(double *dataTable, double **paramArrays,
                            unsigned int *tableSize, double *values,
                            unsigned int ndim, int *index, double *boxValues,
                            bool extrapolate);

  static void hunt(double *xx, int n, double x, int *jlo);
  static void convertToBits(unsigned int number, unsigned int *buffer,
                            unsigned int bufSize);
};

#endif // NDIMINTERPOLATOR_HPP
