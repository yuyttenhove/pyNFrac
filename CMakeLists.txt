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

cmake_minimum_required(VERSION 3.5)

project(pyNFrac)

# Set project version
set(PYNFRAC_VERSION_MAJOR 1)
set(PYNFRAC_VERSION_MINOR 1)

# Enable unit testing using CTest
enable_testing()

# Macro to safely set compiler flags
macro(add_compiler_flag FLAG FATAL)
    set(BACKUP_CXX_FLAGS ${CMAKE_CXX_FLAGS})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAG}")
    try_compile(FLAG_WORKS ${PROJECT_BINARY_DIR}
            ${PROJECT_BINARY_DIR}/mintest.cpp)
    if(NOT FLAG_WORKS)
        if(${FATAL} STREQUAL "REQUIRED")
            message(FATAL_ERROR "The compiler does not support the ${FLAG} flag!")
        else(${FATAL} STREQUAL "REQUIRED")
            message(STATUS "Not using unsupported compiler flag ${FLAG}.")
            set(CMAKE_CXX_FLAGS ${BACKUP_CXX_FLAGS})
        endif(${FATAL} STREQUAL "REQUIRED")
    endif(NOT FLAG_WORKS)
endmacro(add_compiler_flag)

# Add minimal cpp file for flag testing
execute_process(COMMAND echo "int main(int argc, char **argv){ return 0; }"
                OUTPUT_FILE ${PROJECT_BINARY_DIR}/mintest.cpp)

# Enable C++11
add_compiler_flag("-std=c++11" REQUIRED)

# Find python - required dependency
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)

# Find Boost with required components
set(BOOST_COMPONENTS
    python3${Python3_VERSION_MINOR}
    numpy3${Python3_VERSION_MINOR}
)
find_package(Boost COMPONENTS ${BOOST_COMPONENTS} REQUIRED)

# Set environment
include_directories(${PROJECT_BINARY_DIR}/src)
include_directories( ${Python3_INCLUDE_DIRS} )
include_directories( ${Boost_INCLUDE_DIR} )
link_libraries(${Boost_LIBRARIES} ${Python3_LIBRARIES})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/rundir)

# Add source files
add_subdirectory(src)

# Add unit tests
add_subdirectory(test)

# Generate documentation
find_package(Doxygen)
if(DOXYGEN_FOUND)
    message(STATUS
            "Doxygen found! You can generate documentation using 'make doc'")
    # Configure the Doxyfile (sets the correct output path)
    configure_file(${PROJECT_SOURCE_DIR}/Doxyfile.in
            ${PROJECT_BINARY_DIR}/Doxyfile @ONLY)
    # Add a command that will generate the doxygen documentation
    add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/doc/html/index.html
            COMMAND ${DOXYGEN_EXECUTABLE}
            ${PROJECT_BINARY_DIR}/Doxyfile
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM)
    # Add a custom target that calls the command defined above
    add_custom_target(doc DEPENDS ${PROJECT_BINARY_DIR}/doc/html/index.html)
else(DOXYGEN_FOUND)
    # Could not find doxygen. Inform the user that documentation will not be
    # available.
    message(WARNING "Doxygen not found. 'make doc' will not work!")
endif(DOXYGEN_FOUND)
