# cmake needs this line
cmake_minimum_required(VERSION 3.1)

# Define project name
project(test)

# find_package( Eigen3 REQUIRED )
find_package(ifopt REQUIRED)
find_package(Armadillo REQUIRED)


include_directories(ARMADILLO_INCLUDE_DIRS)
include_directories( EIGEN3_INCLUDE_DIR ifopt_INCLUDE_DIRS)
include_directories(/home/damon/softwares/eigen-master)

# add_executable( ifopt_test ifopt_test.cpp )
# add_executable( test test.cpp )
# add_executable(findShortestPath findShortestPath.cpp)
# add_executable(mainHD mainHerdDouble.cpp)
# target_link_libraries(mainHD -O2 -larmadillo)
#  add_executable(test1 motionPlanForDefOpenForm.cpp)
#  target_link_libraries(test1 -O2 -larmadillo)
 add_executable(test2 interSecLineLine.cpp)
 target_link_libraries(test2 -O2 -larmadillo)
#target_link_libraries(findShortestPath Armadillo)
# target_link_libraries( ifopt_test ifopt::ifopt_core ifopt::ifopt_ipopt)
