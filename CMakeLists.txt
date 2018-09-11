cmake_minimum_required(VERSION 2.8.11)
project(EqnPlotter)

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")

## User defined library/include paths
include(PkgConfigPath.cmake)

## mio
set(MIO_INCLUDE_DIR "${CMAKE_CURRENT_LIST_DIR}/modules")
include_directories(${MIO_INCLUDE_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${MIO_INCLUDE_DIR}/mio/cmake/Modules")

## Setup Release and Debug variables
include("${MIO_INCLUDE_DIR}/mio/cmake/DefaultConfigTypes.cmake")

## OpenGL
# OPENGL_FOUND, OPENGL_INCLUDE_DIR, OPENGL_LIBRARIES
find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIR})

## Qt5
find_package(Qt5 5 REQUIRED COMPONENTS Widgets Core OpenGL Svg Concurrent PrintSupport Xml)

## Qwt
find_package(Qwt REQUIRED)
include_directories(${QWT_INCLUDE_DIRS})

## OpenCV
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
set(OPENCV_LIBRARIES ${OpenCV_LIBS})
message(STATUS "Found OpenCV version ${OpenCV_VERSION_MAJOR}.${OpenCV_VERSION_MINOR}.${OpenCV_VERSION_PATCH}")

## VTK
option(WITH_VTK "Build with VTK support" OFF)
if(WITH_VTK)
  find_package(VTK QUIET)
  if(VTK_FOUND)
    message(STATUS "Building with the VTK support")
    include(${VTK_USE_FILE})
    add_definitions(-DHAVE_VTK=${VTK_FOUND})
  else()
    message(STATUS "Could not find VTK. Building without VTK support.")
  endif()
endif()

## libmio_adv_slider_widget
add_subdirectory(${MIO_INCLUDE_DIR}/mio/qt/adv_slider_widget ${CMAKE_CURRENT_BINARY_DIR}/mio/qt/adv_slider_widget)
add_definitions(-DHAVE_ADV_SLIDER_WIDGET)

## libmio_qwt_graph
qt5_wrap_cpp(mio_qwt_graph_HEADERS_MOC ${MIO_INCLUDE_DIR}/mio/qt/qwt_graph.h)
add_library(mio_qwt_graph ${mio_qwt_graph_HEADERS_MOC})
target_link_libraries(mio_qwt_graph Qt5::Widgets Qt5::PrintSupport Qt5::Svg Qt5::OpenGL
                                   ${QWT_LIBRARIES} ${OPENGL_LIBRARIES})

set(EQN_PLOTTER_DIR "${CMAKE_CURRENT_LIST_DIR}/modules/eqn_plotter")
include_directories(${EQN_PLOTTER_DIR})
add_executable(eqn_plotter "${EQN_PLOTTER_DIR}/eqn_plotter.cpp")
target_link_libraries(eqn_plotter mio_adv_slider_widget mio_qwt_graph ${OPENCV_LIBRARIES})

