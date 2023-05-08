cmake_minimum_required(VERSION 3.1)
project(stdev)

# Find Boost libraries
find_package(Boost REQUIRED)

# Add executable and link Boost libraries
add_executable(stdev stdev.cpp)
target_link_libraries(stdev Boost::boost)
