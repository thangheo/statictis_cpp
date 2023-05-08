cmake_minimum_required(VERSION 3.1)
project(psd)

# Find Boost libraries
find_package(Boost REQUIRED)

# Add executable and link Boost libraries
add_executable(psd psd.cpp)
target_link_libraries(psd Boost::boost)
