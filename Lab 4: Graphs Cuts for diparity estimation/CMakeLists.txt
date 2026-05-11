cmake_minimum_required(VERSION 3.12)
project(GraphCuts)
list(APPEND CMAKE_FRAMEWORK_PATH /Library/Frameworks) #Mac, why not auto?
find_package(Imagine REQUIRED COMPONENTS Images)

add_executable(exampleGC exampleGC.cpp maxflow/graph.cpp)
add_executable(GCDisparity GCDisparity.cpp maxflow/graph.cpp)
target_link_libraries(GCDisparity Imagine::Images)
