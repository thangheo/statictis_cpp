cmake_minimum_required(VERSION 3.1)
project(hist)

# Find Boost libraries
find_package(Boost REQUIRED)

# Add executable and link Boost libraries
add_executable(hist hist.cpp)
target_link_libraries(hist Boost::boost)
