cmake_minimum_required(VERSION 3.1)
project(Hurst)

# Find Boost libraries
find_package(Boost REQUIRED)

# Add executable and link Boost libraries
add_executable(Hurst Hurst.cpp)
target_link_libraries(Hurst Boost::boost)
