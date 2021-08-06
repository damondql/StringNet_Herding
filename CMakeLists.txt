# cmake needs this line
cmake_minimum_required(VERSION 3.1)
 
# Define project name
project(test)
 
set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")
 
find_package(Armadillo REQUIRED)
find_library(GUROBI_LIBRARY
   NAMES gurobi gurobi912
   HINTS ${GUROBI_DIR} $ENV{GUROBI_HOME}
   PATH_SUFFIXES lib)
 
include_directories(ARMADILLO_INCLUDE_DIRS)
include_directories(/home/dasclab/SupportingLib/gurobi912/linux64/include)
 
 
# add_executable(mainHD mainHerdDouble.cpp)
# target_link_libraries(mainHD "/home/dasclab/SupportingLib/gurobi912/linux64/lib/libgurobi_c++.a" "/home/dasclab/SupportingLib/gurobi912/linux64/lib/libgurobi91.so" -O2 -larmadillo)
 
 
# add_executable(mainHDExp mainHerdDoubleExperiment.cpp)
# target_link_libraries(mainHDExp "/home/dasclab/SupportingLib/gurobi912/linux64/lib/libgurobi_c++.a" "/home/dasclab/SupportingLib/gurobi912/linux64/lib/libgurobi91.so" -O2 -larmadillo)

add_executable(test test.cpp)
target_link_libraries(test "/home/dasclab/SupportingLib/gurobi912/linux64/lib/libgurobi_c++.a" "/home/dasclab/SupportingLib/gurobi912/linux64/lib/libgurobi91.so" -O2 -larmadillo)