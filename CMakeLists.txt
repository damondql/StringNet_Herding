# cmake needs this line
cmake_minimum_required(VERSION 3.1)

# Define project name
project(test)

set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")

# find_package( Eigen3 REQUIRED )
find_package(ifopt REQUIRED)
find_package(Armadillo REQUIRED)
find_library(GUROBI_LIBRARY
    NAMES gurobi gurobi912
    HINTS ${GUROBI_DIR} $ENV{GUROBI_HOME}
    PATH_SUFFIXES lib)

include_directories(ARMADILLO_INCLUDE_DIRS)
include_directories( EIGEN3_INCLUDE_DIR ifopt_INCLUDE_DIRS)
include_directories(/home/damon/softwares/eigen-master)
include_directories(/opt/gurobi912/linux64/include)


# add_executable( ifopt_test ifopt_test.cpp )
# add_executable( test test.cpp )
# add_executable(findShortestPath findShortestPath.cpp)
# add_executable(mainHD mainHerdDouble.cpp)
# target_link_libraries(mainHD "/opt/gurobi912/linux64/lib/libgurobi_c++.a" "/opt/gurobi912/linux64/lib/libgurobi91.so" -O2 -larmadillo)
#  add_executable(test1 motionPlanForDefOpenForm.cpp)
#  target_link_libraries(test1 -O2 -larmadillo)
#  add_executable(test2 interSecLineLine.cpp)
#  target_link_libraries(test2 -O2 -larmadillo)
#  add_executable(test3 motionPlanForDefOpenForm.cpp)
#  target_link_libraries(test3 "/opt/gurobi912/linux64/lib/libgurobi_c++.a" "/opt/gurobi912/linux64/lib/libgurobi91.so" -O2 -larmadillo)
#  add_executable(test3 defGoalAssignMIQP.cpp)
#  target_link_libraries(test3 "/opt/gurobi912/linux64/lib/libgurobi_c++.a" "/opt/gurobi912/linux64/lib/libgurobi91.so" -O2 -larmadillo)
# add_executable(test4 defInitDesirePos.cpp)
#  target_link_libraries(test4 "/opt/gurobi912/linux64/lib/libgurobi_c++.a" "/opt/gurobi912/linux64/lib/libgurobi91.so" -O2 -larmadillo)
add_executable(test5 controlAttacker4.cpp)
target_link_libraries(test5 "/opt/gurobi912/linux64/lib/libgurobi_c++.a" "/opt/gurobi912/linux64/lib/libgurobi91.so" -O2 -larmadillo)
#target_link_libraries(findShortestPath Armadillo)
# target_link_libraries( ifopt_test ifopt::ifopt_core ifopt::ifopt_ipopt)