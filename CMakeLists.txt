# cmake needs this line
cmake_minimum_required(VERSION 3.1)

# Define project name
project(test)

# find_package( Eigen3 REQUIRED )
find_package(ifopt REQUIRED)
include_directories( EIGEN3_INCLUDE_DIR ifopt_INCLUDE_DIRS)
include_directories(/home/damon/softwares/eigen-master)

add_executable( ifopt_test ifopt_test.cpp )
target_link_libraries( ifopt_test ifopt::ifopt_core ifopt::ifopt_ipopt)
