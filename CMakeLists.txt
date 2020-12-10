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

# Find python and Boost - both are required dependencies
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(Boost COMPONENTS python38 numpy38 REQUIRED)

# Without this, any build libraries automatically have names "lib{x}.so"
set(CMAKE_SHARED_MODULE_PREFIX "")

# Set environment
include_directories(${PROJECT_BINARY_DIR}/src)
include_directories( ${Python_INCLUDE_DIRS} )
include_directories( ${Boost_INCLUDE_DIR} )
link_libraries(${Boost_LIBRARIES} ${Python_LIBRARIES})
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

# Add a shared module - modules are intended to be imported at runtime.
# - This is where you add the source files
#add_library(pyNFrac SHARED src/pyNFrac.cpp src/NDimInterpolator.cpp src/NDimInterpolator.hpp src/NeutralFracTable.cpp src/NeutralFracTable.hpp src/NeutralFracTableDataLocation.hpp)

# Set up the libraries and header search paths for this target
#target_link_libraries(pyNFrac ${Boost_LIBRARIES})
#target_include_directories(pyNFrac PRIVATE ${PYTHON_INCLUDE_DIRS})

#set_target_properties( pyNFrac PROPERTIES PREFIX "" OUTPUT_NAME pyNFrac)