cmake_minimum_required(VERSION 3.10)

# set the project name
project(BFieldModels)

# Set build type and optimisations (assume unix like platform)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")


# Find all executables
file(GLOB MAINS ${PROJECT_SOURCE_DIR}/mains/*.cpp)

# Add sources
file(GLOB SRC_FILES ${PROJECT_SOURCE_DIR}/src/*)

# Set include path
include_directories(${PROJECT_SOURCE_DIR}/include)

# Setup all programs in mains folder
foreach(MAIN ${MAINS})
   get_filename_component(MAIN_NAME ${MAIN} NAME_WLE )
   add_executable(${MAIN_NAME} ${MAIN} ${SRC_FILES})
   #install(${MAIN_NAME})
endforeach(MAIN ${MAINS})

