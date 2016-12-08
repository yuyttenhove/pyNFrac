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
 * @file NeutralFracTable.hpp
 *
 * @brief Neutral fraction table.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef NEUTRALFRACTABLE_HPP
#define NEUTRALFRACTABLE_HPP

/*! @brief Number of temperature values. */
#define NEUTRALFRACTABLE_NTEMP 351
/*! @brief Number of [Fe/H] values. */
#define NEUTRALFRACTABLE_NFEH 7
/*! @brief Number of [Mg/Fe] values. */
#define NEUTRALFRACTABLE_NMGFE 3
/*! @brief Number of redshift values. */
#define NEUTRALFRACTABLE_NZ 11
/*! @brief Number of density values. */
#define NEUTRALFRACTABLE_NDENS 8

#include <string>

/**
 * @brief Neutral fraction table.
 */
class NeutralFracTable {
private:
  /*! @brief Temperature array. */
  double _temp_array[NEUTRALFRACTABLE_NTEMP];

  /*! @brief [Fe/H] array. */
  double _feh_array[NEUTRALFRACTABLE_NFEH];

  /*! @brief [Mg/Fe] array. */
  double _mgfe_array[NEUTRALFRACTABLE_NMGFE];

  /*! @brief Redshift array. */
  double _z_array[NEUTRALFRACTABLE_NZ];

  /*! @brief Density array. */
  double _dens_array[NEUTRALFRACTABLE_NDENS];

  /*! @brief Neutral fraction table. */
  double _neutralFracTable[NEUTRALFRACTABLE_NTEMP * NEUTRALFRACTABLE_NFEH *
                           NEUTRALFRACTABLE_NMGFE * NEUTRALFRACTABLE_NZ *
                           NEUTRALFRACTABLE_NDENS];

  /*! @brief Density conversion factors. */
  double
      _densityConversionFactors[NEUTRALFRACTABLE_NFEH * NEUTRALFRACTABLE_NMGFE];

public:
  NeutralFracTable();

  static std::string get_filename(double FeH, double MgFe, double z,
                                  double rho);

  unsigned int get_curve(unsigned int iFeH, unsigned int iMgFe, unsigned int iz,
                         unsigned int irho, double *Tarr, double *narr,
                         unsigned int size);

  double get_neutral_fraction(double FeH, double MgFe, double rho, double z,
                              double T);
};

#endif // NEUTRALFRACTABLE_HPP
