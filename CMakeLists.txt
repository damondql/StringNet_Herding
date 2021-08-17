# cmake needs this line
cmake_minimum_required(VERSION 3.1)
 
# Define project name
project(test)
 
set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")
 
if(CMAKE_HOST_WIN32)
    set(libname "libconfig")
else()
    set(libname "config")
endif()


find_package(Armadillo REQUIRED)
find_library(GUROBI_LIBRARY
   NAMES gurobi gurobi912
   HINTS ${GUROBI_DIR} $ENV{GUROBI_HOME}
   PATH_SUFFIXES lib)
find_library(LIBCONFIGPP_LOCATION libconfig++)
 
include_directories(ARMADILLO_INCLUDE_DIRS)
include_directories(/home/dasclab/SupportingLib/gurobi912/linux64/include)
# include_directories(/home/dasclab/SupportingLib/libconfig-1.7.3/lib)
 
 
# add_executable(mainHD mainHerdDouble.cpp)
# target_link_libraries(mainHD "/home/dasclab/SupportingLib/gurobi912/linux64/lib/libgurobi_c++.a" "/home/dasclab/SupportingLib/gurobi912/linux64/lib/libgurobi91.so" -O2 -larmadillo)

# add_executable(mainHD mainHerdDoubleNew.cpp)
# target_link_libraries(mainHD "/home/dasclab/SupportingLib/gurobi912/linux64/lib/libgurobi_c++.a" "/home/dasclab/SupportingLib/gurobi912/linux64/lib/libgurobi91.so" -O2 -larmadillo)
 
 
# add_executable(mainHDExp mainHerdDoubleExperiment.cpp)
# target_link_libraries(mainHDExp "/home/dasclab/SupportingLib/gurobi912/linux64/lib/libgurobi_c++.a" "/home/dasclab/SupportingLib/gurobi912/linux64/lib/libgurobi91.so" -O2 -larmadillo)

# add_executable(test test.cpp)
# target_link_libraries(test "/home/dasclab/SupportingLib/gurobi912/linux64/lib/libgurobi_c++.a" "/home/dasclab/SupportingLib/gurobi912/linux64/lib/libgurobi91.so" -O2 -larmadillo)

# add_executable(grbTest mip1_c++.cpp)
# target_link_libraries(grbTest "/home/dasclab/SupportingLib/gurobi912/linux64/lib/libgurobi_c++.a" "/home/dasclab/SupportingLib/gurobi912/linux64/lib/libgurobi91.so" -O2 -larmadillo)
 
add_executable(c++_example4	config_example.cpp )
target_include_directories(c++_example4 PRIVATE ${CMAKE_SOURCE_DIR}/lib)
target_link_libraries(c++_example4 ${libname}++ )

add_executable(test1	readConfig.cpp )
target_include_directories(test1 PRIVATE ${CMAKE_SOURCE_DIR}/lib)
target_link_libraries(test1 ${libname}++ )