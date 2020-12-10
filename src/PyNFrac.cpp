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
 * @file PyNFrac.cpp
 *
 * @brief Python module functions.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "NeutralFracTable.hpp"

#include <cstdlib>
#include <iostream>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace p = boost::python;
namespace np = boost::python::numpy;

/**
 * @brief Get the length of a 1D numpy.ndarray.
 *
 * If the array is not 1D, an error is thrown.
 *
 * @param a numpy.ndarray.
 * @return Size of the 1D array.
 */
static unsigned int get_size_1D_array(np::ndarray &a) {
    p::tuple ashape = p::extract<p::tuple>(a.attr("shape"));
    const unsigned int adim = p::len(ashape);
    if (adim != 1) {
        std::cerr << "Wrong shape for 1D array (dimension = " << adim << ")!"
                  << std::endl;
        abort();
    }
    return p::extract<unsigned int>(ashape[0]);
}

/**
 * @brief Get the neutral fractions corresponding to the given parameter values.
 *
 * @param FeHarr numpy.ndarray with [Fe/H] values.
 * @param MgFearr numpy.ndarray with [Mg/Fe] values.
 * @param rhoarr numpy.ndarray with density values (in g cm^-3).
 * @param Tarr numpy.ndarray with temperature values (in K).
 * @param z Redshift value.
 * @return numpy.ndarray containing the neutral fractions of hydrogen.
 */
static np::ndarray
get_neutral_fractions(np::ndarray &FeHarr,
                      np::ndarray &MgFearr,
                      np::ndarray &rhoarr,
                      np::ndarray &Tarr, double z) {
    const unsigned int nFeH = get_size_1D_array(FeHarr);
    const unsigned int nMgFe = get_size_1D_array(MgFearr);
    const unsigned int nrho = get_size_1D_array(rhoarr);
    const unsigned int nT = get_size_1D_array(Tarr);

    if (nFeH != nMgFe || nFeH != nrho || nFeH != nT) {
        std::cerr << "Not all arrays have the same size: " << nFeH << ", " << nMgFe
                  << ", " << nrho << ", " << nT << "!" << std::endl;
        abort();
    }

    // we create the ndarray that will store the results
    p::tuple shape = p::make_tuple(nFeH);
    np::ndarray result = np::zeros(shape, np::dtype::get_builtin<double>());

    // create the NeutralFracTable
    NeutralFracTable table;

    // now fill the array
    for (unsigned int i = 0; i < nFeH; ++i) {
        double FeH = p::extract<double>(FeHarr[i]);
        double MgFe = p::extract<double>(MgFearr[i]);
        double rho = p::extract<double>(rhoarr[i]);
        double T = p::extract<double>(Tarr[i]);

            result[i] = table.get_neutral_fraction(FeH, MgFe, rho, z, T);
    }

    return result;
}

/**
 * @brief Python module exposure.
 */
BOOST_PYTHON_MODULE (pyNFrac) {
    Py_Initialize();
    np::initialize();
    def("get_neutral_fractions", &get_neutral_fractions);
}
