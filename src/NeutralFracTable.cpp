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
 * @file NeutralFracTable.cpp
 *
 * @brief NeutralFracTable implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "NeutralFracTable.hpp"
#include "NeutralFracTableDataLocation.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

/**
 * @brief Constructor.
 */
NeutralFracTable::NeutralFracTable()
    : _feh_array{-99., -4., -2., -1., -0.5, 0., 0.5},
      _mgfe_array{-0.55, 0., 0.47},
      _z_array{0., 0.2, 0.5, 1., 2., 4., 10., 10.25, 10.5, 10.65, 11.},
      _dens_array{1.e-9, 1.e-6, 1.e-4, 1.e-3, 1.e-2, 1.e-1, 1., 1.e2} {
  for (unsigned int iFeH = 0; iFeH < 7; ++iFeH) {
    for (unsigned int iMgFe = 0; iMgFe < 3; ++iMgFe) {
      for (unsigned int iz = 0; iz < 11; ++iz) {
        for (unsigned int irho = 0; irho < 8; ++irho) {
          // we don't have data files for these values
          if (_feh_array[iFeH] != -99. || _mgfe_array[iMgFe] == 0.) {
            std::string filename =
                get_filename(_feh_array[iFeH], _mgfe_array[iMgFe], _z_array[iz],
                             _dens_array[irho]);
            std::ifstream file(filename);
            std::string line;
            getline(file, line);
            // extract density conversion factor
            std::istringstream lstream(line);
            lstream >> _densityConversionFactors[iFeH][iMgFe];
            for (unsigned int itemp = 0; itemp < NEUTRALFRACTABLE_NTEMP;
                 ++itemp) {
              getline(file, line);
              std::istringstream lstream(line);
              lstream >> _temp_array[itemp] >>
                  _neutralFracTable[itemp][iFeH][iMgFe][iz][irho];
            }
          } else {
            for (unsigned int itemp = 0; itemp < NEUTRALFRACTABLE_NTEMP;
                 ++itemp) {
              _neutralFracTable[itemp][iFeH][iMgFe][iz][irho] = 0.;
            }
          }
        }
      }
    }
  }
}

/**
 * @brief Get the name of the file containing the data for the given parameter
 * values.
 *
 * @param FeH [Fe/H] value.
 * @param MgFe [Mg/Fe] value.
 * @param z Redshift value.
 * @param rho Density value (in cm^-3).
 * @return Name of the file containing the data for those parameter values.
 */
std::string NeutralFracTable::get_filename(double FeH, double MgFe, double z,
                                           double rho) {
  std::stringstream stream;
  stream << NEUTRALFRACTABLEDATALOCATION << "/RadLoss_";
  stream << std::setprecision(2) << std::fixed << FeH << "_";
  stream << std::setprecision(2) << std::fixed << MgFe << "_";
  stream << std::setprecision(2) << std::fixed << z << "_";
  stream << std::setprecision(2) << std::scientific << rho << ".neutral";
  return stream.str();
}

/**
 * @brief Get the neutral fraction curve for the given parameters indices.
 *
 * @param iFeH [Fe/H] value index.
 * @param iMgFe [Mg/Fe] value index.
 * @param iz Redshift value index.
 * @param irho Number density value index.
 * @param Tarr Empty temperature array of size size that will be filled.
 * @param narr Empty neutral fraction array of size size that will be filled.
 * @param size Size of the empty arrays.
 * @return Number of elements written to the arrays by this routine.
 */
unsigned int NeutralFracTable::get_curve(unsigned int iFeH, unsigned int iMgFe,
                                         unsigned int iz, unsigned int irho,
                                         double *Tarr, double *narr,
                                         unsigned int size) {
  unsigned int readsize = NEUTRALFRACTABLE_NTEMP;
  if (readsize > size) {
    readsize = size;
  }
  for (unsigned int itemp = 0; itemp < readsize; ++itemp) {
    Tarr[itemp] = _temp_array[itemp];
    narr[itemp] = _neutralFracTable[itemp][iFeH][iMgFe][iz][irho];
  }

  return readsize;
}
