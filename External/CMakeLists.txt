# ******************************************************************************
# MIT License
#
# Copyright (c) 2024 Moritz Laupichler
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# ******************************************************************************


# Add RoutingKit library (needs to be built beforehand!) and include headers from git submodule.
# Use interface target to bundle sources and headers.
find_library(RoutingKit_LIBRARY routingkit PATHS RoutingKit/lib)
if (NOT RoutingKit_LIBRARY)
    message(FATAL_ERROR "RoutingKit library (e.g. libroutingkit.so) not found. Did you forget to build the library (see ReadMe)?")
endif ()
set(RoutingKit_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/RoutingKit/include/")
add_library(routingkit INTERFACE)
target_link_libraries(routingkit INTERFACE ${RoutingKit_LIBRARY})
# include as SYSTEM to avoid warnings for features deprecated in C++20 as RoutingKit is C++17
target_include_directories(routingkit SYSTEM INTERFACE ${RoutingKit_INCLUDE_DIR})

# Add fast-cpp-csv-parser header-only library from git submodule.
add_library(fast_cpp_csv_parser INTERFACE)
target_include_directories(fast_cpp_csv_parser SYSTEM INTERFACE fast-cpp-csv-parser)

# Add kassert header-only library from git submodule.
add_library(kassert INTERFACE)
target_include_directories(kassert SYSTEM INTERFACE kassert/include)
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(KASSERT_ASSERTION_LEVEL 40)
else ()
    set(KASSERT_ASSERTION_LEVEL 20)
endif ()
target_compile_definitions(kassert INTERFACE -DKASSERT_ASSERTION_LEVEL=${KASSERT_ASSERTION_LEVEL})

# Add vectorclass header-only library from git submodule.
add_library(vectorclass INTERFACE)
target_include_directories(vectorclass SYSTEM INTERFACE vectorclass)


# Add bundled rapidxml header only library as interface target
add_library(rapidxml INTERFACE)
target_include_directories(rapidxml SYSTEM INTERFACE rapidxml)

# Add bundled randomc header only library as interface target
add_library(randomc INTERFACE)
target_include_directories(randomc SYSTEM INTERFACE randomc)

# Add bundled stocc header only library as interface target
add_library(stocc INTERFACE)
target_include_directories(stocc SYSTEM INTERFACE stocc)

add_subdirectory(ttlsa)