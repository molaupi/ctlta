# Specify all dependencies fetched at configure time using FetchContent
include(FetchContent)

# Declare kassert dependency
FetchContent_Declare(
        kassert
        GIT_REPOSITORY https://github.com/kamping-site/kassert
        GIT_TAG v0.1.0
)

# Declare fast-cpp-csv-parser dependency
FetchContent_Declare(
        fast_cpp_csv_parser
        GIT_REPOSITORY https://github.com/ben-strasser/fast-cpp-csv-parser
)

# Declare vectorlclass dependency
FetchContent_Declare(
        vectorclass
        GIT_REPOSITORY https://github.com/vectorclass/version2.git
)

# Declare proj dependency
FetchContent_Declare(
        proj
        URL https://download.osgeo.org/proj/proj-9.5.0.tar.gz
        URL_MD5 ac46b4e31562890d012ea6b31e579cf6
)

# Fetch kassert
message("Fetching kassert library...")
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(KASSERT_ASSERTION_LEVEL 40)
else ()
    set(KASSERT_ASSERTION_LEVEL 20)
endif ()
FetchContent_MakeAvailable(kassert)

# Fetch fast-cpp-csv-parser (header only library, create interface target)
FetchContent_MakeAvailable(fast_cpp_csv_parser)
FetchContent_GetProperties(fast_cpp_csv_parser SOURCE_DIR fast_cpp_csv_parser_SOURCE_DIR)
add_library(fast_cpp_csv_parser INTERFACE)
target_include_directories(fast_cpp_csv_parser SYSTEM INTERFACE ${fast_cpp_csv_parser_SOURCE_DIR})

# Fetch vectorclass (header only library, create interface target)
FetchContent_MakeAvailable(vectorclass)
FetchContent_GetProperties(vectorclass SOURCE_DIR vectorclass_SOURCE_DIR)
add_library(vectorclass INTERFACE)
target_include_directories(vectorclass SYSTEM INTERFACE ${vectorclass_SOURCE_DIR})

# Fetch proj
message("Fetching proj library...")
set(BUILD_APPS OFF)
set(BUILD_TESTING OFF)
set(ENABLE_CURL OFF)
set(ENABLE_TIFF OFF)
set(TESTING_USE_NETWORK OFF)
FetchContent_MakeAvailable(proj)