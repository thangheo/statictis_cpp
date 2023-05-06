cmake_minimum_required(VERSION 3.1)
project(var)

# Find Boost libraries
find_package(Boost REQUIRED)

# Add executable and link Boost libraries
add_executable(var var.cpp)
target_link_libraries(var Boost::boost)
