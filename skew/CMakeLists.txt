cmake_minimum_required(VERSION 3.1)
project(skew)

# Find Boost libraries
find_package(Boost REQUIRED)

# Add executable and link Boost libraries
add_executable(skew skew.cpp)
target_link_libraries(skew Boost::boost)
