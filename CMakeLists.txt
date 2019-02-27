CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(topopt99)

INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/LibSL-small
  ${PROJECT_SOURCE_DIR}/LibSL-small/src/
  ${PROJECT_SOURCE_DIR}/LibSL-small/src/LibSL
  ${PROJECT_SOURCE_DIR}/eigen-git-mirror/
)

CONFIGURE_FILE(
  "${CMAKE_CURRENT_SOURCE_DIR}/path.h.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/path.h"
  )

SET(SOURCES 
  topopt.cpp
  LibSL-small/src/LibSL/Math/Math.cpp
  LibSL-small/src/LibSL/Math/Vertex.cpp
  LibSL-small/src/LibSL/Mesh/Mesh.cpp
  LibSL-small/src/LibSL/Mesh/MeshFormat_stl.cpp
  LibSL-small/src/LibSL/Mesh/VertexFormat_dynamic.cpp
  LibSL-small/src/LibSL/Image/Image.cpp
  LibSL-small/src/LibSL/Image/ImageFormat_TGA.cpp
  LibSL-small/src/LibSL/Image/tga.cpp
  LibSL-small/src/LibSL/System/System.cpp
  LibSL-small/src/LibSL/CppHelpers/CppHelpers.cpp
)

ADD_EXECUTABLE(topopt99 ${SOURCES})

if(WIN32)
target_link_libraries(topopt99 shlwapi)
endif(WIN32)
