cmake_minimum_required(VERSION 3.5)
set(CMAKE_OSX_ARCHITECTURES "x86_64" CACHE STRING "") #Apple Silicon
project(Panorama)
list(APPEND CMAKE_FRAMEWORK_PATH /Library/Frameworks) #Mac, why not auto?
find_package(Imagine REQUIRED COMPONENTS LinAlg Images)

add_executable(Panorama Panorama.cpp
        Panorama.cpp
        Panorama.cpp
        Panorama.cpp)
target_link_libraries(Panorama Imagine::LinAlg Imagine::Images)
