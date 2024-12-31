# Add library target for external TTL implementation
add_library(ttlsa STATIC src/road_network.cpp src/util.cpp)
target_include_directories(ttlsa PUBLIC include)


if((NUM_THREADS GREATER 1))
    target_compile_definitions(ttlsa PRIVATE NUM_THREADS=${NUM_THREADS})
else()
    target_compile_definitions(ttlsa PRIVATE NUM_THREADS=1)
endif()