cmake_minimum_required(VERSION 3.1)
project(kurt)

# Find Boost libraries
find_package(Boost REQUIRED)

# Add executable and link Boost libraries
add_executable(kurt kurt.cpp)
target_link_libraries(kurt Boost::boost)
