cmake_minimum_required(VERSION 3.1)
project(corr)

# Find Boost libraries
find_package(Boost REQUIRED)


# Add executable and link Boost libraries
add_executable(corr corr.cpp)
target_link_libraries(corr Boost::boost)
