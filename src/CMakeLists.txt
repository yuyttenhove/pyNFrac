################################################################################
# This file is part of pyNFrac
# Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#
# pyNFrac is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyNFrac is distributed in the hope that it will be useful,
# but WITOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with pyNFrac. If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Configure data files
set(NFRACTARFILE ${PROJECT_SOURCE_DIR}/data/neutral_frac_tables.tar.gz)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/data)
execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf ${NFRACTARFILE}
                WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/data/)
set(NEUTRALFRACTABLEDATALOCATION ${PROJECT_BINARY_DIR}/data)
configure_file(${PROJECT_SOURCE_DIR}/src/NeutralFracTableDataLocation.hpp.in
               ${PROJECT_BINARY_DIR}/src/NeutralFracTableDataLocation.hpp
               @ONLY)

# Add program source files
set(PYNFRAC_SOURCES
    # Python module
    PyNFrac.cpp

    # NeutralFracTable
    NeutralFracTable.cpp
    NeutralFracTable.hpp
)

add_library(pyNFrac SHARED ${PYNFRAC_SOURCES})

# Remove the 'lib' prefix from the library name.
set_target_properties(pyNFrac PROPERTIES PREFIX "")
