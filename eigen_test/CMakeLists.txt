cmake_minimum_required(VERSION 3.5)

project(EigenExample)

# Find Eigen
find_package(Eigen3 REQUIRED)

# Add the executable
add_executable(eigen_example eigen.cpp)

# Link against Eigen
target_link_libraries(eigen_example Eigen3::Eigen)
