cmake_minimum_required(VERSION 3.1)
project(covar)

# Find Boost libraries
find_package(Boost REQUIRED)


# Add executable and link Boost libraries
add_executable(covar covar.cpp)
target_link_libraries(covar Boost::boost)
