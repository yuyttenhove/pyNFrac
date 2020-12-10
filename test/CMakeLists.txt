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

### Convenient macros to automate unit test generation #########################
### Do not touch the code below unless you know what you're doing! #############

# Add a new unit test
# A new target with the test sources is constructed, and a CTest test with the
# same name is created. The new test is also added to the global list of test
# contained in the check target
macro(add_unit_test)
    set(oneValueArgs NAME)
    set(multiValueArgs SOURCES LIBS)
    cmake_parse_arguments(TEST "${options}" "${oneValueArgs}"
                               "${multiValueArgs}" ${ARGN})
    message(STATUS "generating " ${TEST_NAME})
    add_executable(${TEST_NAME} EXCLUDE_FROM_ALL ${TEST_SOURCES})
    set_target_properties(${TEST_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                          ${PROJECT_BINARY_DIR}/rundir/test)
    target_link_libraries(${TEST_NAME} ${TEST_LIBS})

    add_test(NAME ${TEST_NAME}
             WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/rundir/test
             COMMAND ${TEST_NAME})

    set(TESTNAMES ${TESTNAMES} ${TEST_NAME})
endmacro(add_unit_test)

macro(add_python_unit_test)
    set(oneValueArgs NAME)
    cmake_parse_arguments(TEST "${options}" "${oneValueArgs}"
                               "${multiValueArgs}" ${ARGN})
    message(STATUS "generating python test " ${TEST_NAME})

    # copy Python test script
    configure_file(${PROJECT_SOURCE_DIR}/test/${TEST_NAME}
                   ${PROJECT_BINARY_DIR}/rundir/test/${TEST_NAME} @ONLY)
    add_test(NAME ${TEST_NAME}
             WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/rundir/test
             COMMAND ${Python3_EXECUTABLE} ${TEST_NAME})

    set(TESTNAMES ${TESTNAMES} ${TEST_NAME})
endmacro(add_python_unit_test)

## Configure the load_module script that is used to load the Python library
set(PYTHON_MODULEPATH ${PROJECT_BINARY_DIR}/rundir)
configure_file(${PROJECT_SOURCE_DIR}/test/load_module.py.in
               ${PROJECT_BINARY_DIR}/rundir/test/load_module.py @ONLY)

### Actual unit test generation ################################################
### Add new unit tests below ###################################################

## Unit test for NeutralFracTable
set(TESTNEUTRALFRACTABLE_SOURCES
    testNeutralFracTable.cpp

    ../src/NDimInterpolator.cpp
    ../src/NDimInterpolator.hpp
    ../src/NeutralFracTable.cpp
    ../src/NeutralFracTable.hpp
)
add_unit_test(NAME testNeutralFracTable
              SOURCES ${TESTNEUTRALFRACTABLE_SOURCES})

## Unit test for NDimInterpolator
set(TESTNDIMINTERPOLATOR_SOURCES
    testNDimInterpolator.cpp

    ../src/NDimInterpolator.cpp
    ../src/NDimInterpolator.hpp
)
add_unit_test(NAME testNDimInterpolator
              SOURCES ${TESTNDIMINTERPOLATOR_SOURCES})

## Unit test for pyNFrac
add_python_unit_test(NAME test_pyNFrac.py)

### Done adding unit tests. Create the 'make check' target #####################
### Do not touch these lines unless you know what you're doing! ################
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
                        DEPENDS ${TESTNAMES})
