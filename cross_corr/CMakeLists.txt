cmake_minimum_required(VERSION 3.1)
project(cross)

# Find Boost libraries
find_package(Boost REQUIRED)

# Add executable and link Boost libraries
add_executable(cross cross.cpp)
target_link_libraries(cross Boost::boost)
