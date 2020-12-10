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
 * @file testNeutralFracTable.cpp
 *
 * @brief Unit test for the NeutralFracTable class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "../src/NeutralFracTable.hpp"
#include "NeutralFracTableDataLocation.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

/**
 * @brief Unit test for the NeutralFracTable class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
    // check if we can reconstruct a single file name
    std::string name = NeutralFracTable::get_filename(-99., 0., 2., 1.e-3);

    std::stringstream refname;
    refname << NEUTRALFRACTABLEDATALOCATION
          << "/RadLoss_-99.00_0.00_2.00_1.00e-03.neutral";
    if (name != refname.str()) {
    std::cerr << "Wrong name: " << name << std::endl;
    abort();
    }

    // now check if we can open all files
    double FeH[7] = {-99., -4., -2., -1., -0.5, 0., 0.5};
    double MgFe[3] = {-0.55, 0., 0.47};
    double z[11] = {0., 0.2, 0.5, 1., 2., 4., 10., 10.25, 10.5, 10.65, 11.};
    double rho[8] = {1.e-9, 1.e-6, 1.e-4, 1.e-3, 1.e-2, 1.e-1, 1., 1.e2};
    for (unsigned int iFeH = 0; iFeH < 7; ++iFeH) {
    for (unsigned int iMgFe = 0; iMgFe < 3; ++iMgFe) {
        for (unsigned int iz = 0; iz < 11; ++iz) {
            for (unsigned int irho = 0; irho < 8; ++irho) {
                // skip files with [Fe/H] == -99 and MgFe != 0., since these do not
                // exist
                if (FeH[iFeH] != -99. || MgFe[iMgFe] == 0.) {
                    name = NeutralFracTable::get_filename(FeH[iFeH], MgFe[iMgFe], z[iz],
                                                          rho[irho]);
                    std::ifstream file(name);
                    if (!file) {
                        std::cerr << "Unable to open file: " << name << "!" << std::endl;
                        abort();
                    }
                }
            }
        }
    }
    }

  // now check some values
  NeutralFracTable table;

  std::ofstream ofile("testcurve.txt");
  double Tarr[351], narr[351];
  table.get_curve(3, 1, 7, 4, Tarr, narr, 351);
  double gasdens = 2.10197728002e-26;
  double fac = rho[4] / gasdens;
  bool wrong = false;
  double FeHval = FeH[3];
  double MgFeval = MgFe[1];
  double rhoval = rho[4] / fac;
  double zval = z[7];
  for (unsigned int i = 0; i < 351; ++i) {
    double nfracH =
        table.get_neutral_fraction(FeHval, MgFeval, rhoval, zval, Tarr[i]);
    ofile << FeHval << "\t" << MgFeval << "\t" << rhoval << "\t" << zval << "\t"
          << Tarr[i] << "\t" << narr[i] << "\t" << nfracH << "\n";
    if (std::abs(narr[i] - nfracH) > 1.e-15 * std::abs(narr[i] + nfracH)) {
      std::cerr << "Error: wrong neutral fraction: " << nfracH << " ("
                << narr[i] << ", relative difference "
                << std::abs(narr[i] - nfracH) / std::abs(narr[i] + nfracH)
                << ")!" << std::endl;
      wrong = true;
    }
  }
  ofile.close();

  if (wrong) {
    return 1;
  }

  return 0;
}
