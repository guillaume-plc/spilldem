cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(SpillDEM VERSION 0.1)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

configure_file(src/SpillDEMConfig.h.in SpillDEM.h)
find_package(GDAL REQUIRED)

# add executable
add_executable(spilldem src/main.cpp)

# include dirs
target_include_directories(spilldem PUBLIC ${PROJECT_BINARY_DIR} ${GDAL_INCLUDE_DIRS})
target_link_libraries(spilldem ${GDAL_LIBRARIES})