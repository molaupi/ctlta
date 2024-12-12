# Add library target for external TTL implementation
add_library(ttlsa src/road_network.cpp src/util.cpp)
target_include_directories(ttlsa PUBLIC include)