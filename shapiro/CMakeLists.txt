cmake_minimum_required(VERSION 3.1)
project(shapiro)

# Find Boost libraries
find_package(Boost REQUIRED)

# Add executable and link Boost libraries
add_executable(shapiro shapiro.cpp)
target_link_libraries(shapiro Boost::boost)
