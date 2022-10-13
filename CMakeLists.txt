cmake_minimum_required(VERSION 3.9)
project(FHO_PreTable LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
find_package(OpenMP REQUIRED)

add_executable(FHO_PreTable main.cpp shockwave.cpp shockwave.h FHO.cpp FHO.h vib_const.h species.h random.cpp random.h)
target_link_libraries(FHO_PreTable PRIVATE OpenMP::OpenMP_CXX)