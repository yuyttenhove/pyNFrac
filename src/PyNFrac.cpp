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

#include <string>

#include <boost/python/def.hpp>
#include <boost/python/module.hpp>

static std::string say_hello() { return "Hello!"; }

BOOST_PYTHON_MODULE(pyNFrac) { boost::python::def("say_hello", &say_hello); }
