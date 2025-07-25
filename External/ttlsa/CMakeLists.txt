# Add library target for external TTL implementation
add_library(ttlsa STATIC src/road_network.cpp src/util.cpp)
target_include_directories(ttlsa PUBLIC include)

set(TTLSA_MAX_VALLEY_PATH_LENGTH 0 CACHE STRING "Choose the maximum number of edges in one shortcut. Set to 0 to use triangles like regular CCH.")
target_compile_definitions(ttlsa PUBLIC TTLSA_MAX_VALLEY_PATH_LENGTH=${TTLSA_MAX_VALLEY_PATH_LENGTH})

if((NUM_THREADS GREATER 1))
    target_compile_definitions(ttlsa PRIVATE NUM_THREADS=${NUM_THREADS})
else()
    target_compile_definitions(ttlsa PRIVATE NUM_THREADS=1)
endif()


