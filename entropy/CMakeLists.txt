cmake_minimum_required(VERSION 3.1)
project(entropy)

# Find Boost libraries
find_package(Boost REQUIRED)

# Add executable and link Boost libraries
add_executable(entropy entropy.cpp)
target_link_libraries(entropy Boost::boost)
