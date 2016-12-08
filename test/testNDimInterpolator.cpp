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
 * @file testNDimInterpolator.cpp
 *
 * @brief Unit test for the NDimInterpolator class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "../src/NDimInterpolator.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

/*! @brief Number of values for the first test parameter. */
#define N1 10
/*! @brief Number of values for the second test parameter. */
#define N2 10
/*! @brief Number of values for the third test parameter. */
#define N3 10
/*! @brief Number of values for the fourth test parameter. */
#define N4 10
/*! @brief Number of values for the fifth test parameter. */
#define N5 10

/**
 * @brief 2D test interpolation function.
 *
 * @param x First parameter.
 * @param y Second parameter.
 * @return Test function value.
 */
double f2(double x, double y) { return 3. * x + 2. * y - 7.; }

/**
 * @brief 5D test interpolation function.
 *
 * @param x1 First parameter.
 * @param x2 Second parameter.
 * @param x3 Third parameter.
 * @param x4 Fourth parameter.
 * @param x5 Fifth parameter.
 * @return Test function value.
 */
double f5(double x1, double x2, double x3, double x4, double x5) {
  return 25. * x1 + 4. * x2 - 17. * x3 + x4 - 0.5 * x5 + 124.;
}

/**
 * @brief Get a random uniform double precision floating point value.
 *
 * @return Random uniform double precision floating point value in the range
 * [0., 1.].
 */
double rand_double() { return ((double)std::rand()) / ((double)RAND_MAX); }

/**
 * @brief Unit test for the NDimInterpolator class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  double x1[N1];
  double x2[N2];
  double x3[N3];
  double x4[N4];
  double x5[N5];
  double d1 = 1. / N1;
  for (unsigned int i = 0; i < N1; ++i) {
    x1[i] = i * d1;
  }
  double d2 = 1. / N2;
  for (unsigned int i = 0; i < N2; ++i) {
    x2[i] = i * d2;
  }
  double d3 = 1. / N3;
  for (unsigned int i = 0; i < N3; ++i) {
    x3[i] = i * d3;
  }
  double d4 = 1. / N4;
  for (unsigned int i = 0; i < N4; ++i) {
    x4[i] = i * d4;
  }
  double d5 = 1. / N5;
  for (unsigned int i = 0; i < N5; ++i) {
    x5[i] = i * d5;
  }

  double f2D[N1 * N2];
  double f5D[N1 * N2 * N3 * N4 * N5];
  for (unsigned int i1 = 0; i1 < N1; ++i1) {
    for (unsigned int i2 = 0; i2 < N2; ++i2) {
      f2D[i1 * N2 + i2] = f2(x1[i1], x2[i2]);
      for (unsigned int i3 = 0; i3 < N3; ++i3) {
        for (unsigned int i4 = 0; i4 < N4; ++i4) {
          for (unsigned int i5 = 0; i5 < N5; ++i5) {
            f5D[i1 * N2 * N3 * N4 * N5 + i2 * N3 * N4 * N5 + i3 * N4 * N5 +
                i4 * N5 + i5] = f5(x1[i1], x2[i2], x3[i3], x4[i4], x5[i5]);
          }
        }
      }
    }
  }

  double *params2[2] = {x1, x2};
  double *params5[5] = {x1, x2, x3, x4, x5};
  unsigned int sizes2[2] = {N1, N2};
  unsigned int sizes5[5] = {N1, N2, N3, N4, N5};
  for (unsigned int i = 0; i < 100; ++i) {
    double values2[2] = {rand_double(), rand_double()};
    int index2[2] = {-1, -1};
    double boxValues2[2] = {-1., -1.};
    double xtest = NDimInterpolator::interpolate(f2D, params2, sizes2, values2,
                                                 2, index2, boxValues2, true);

    double tval2 = f2(values2[0], values2[1]);
    if (std::abs(xtest - tval2) > 1.e-15 * std::abs(xtest + tval2)) {
      std::cerr << "Wrong interpolated value: " << xtest << " (" << tval2
                << ", relative difference: "
                << std::abs(xtest - tval2) / std::abs(xtest + tval2) << ")!"
                << std::endl;
      abort();
    }

    double values5[5] = {rand_double(), rand_double(), rand_double(),
                         rand_double(), rand_double()};
    int index5[5] = {-1, -1, -1, -1, -1};
    double boxValues5[5] = {-1., -1., -1., -1., -1.};
    xtest = NDimInterpolator::interpolate(f5D, params5, sizes5, values5, 5,
                                          index5, boxValues5, true);

    double tval5 =
        f5(values5[0], values5[1], values5[2], values5[3], values5[4]);
    if (std::abs(xtest - tval5) > 1.e-15 * std::abs(xtest + tval5)) {
      std::cerr << "Wrong interpolated value: " << xtest << " (" << tval5
                << ", relative difference: "
                << std::abs(xtest - tval5) / std::abs(xtest + tval5) << ")!"
                << std::endl;
      abort();
    }
  }

  return 0;
}
