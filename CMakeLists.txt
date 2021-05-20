cmake_minimum_required (VERSION 2.9)

project (jpacTriangle)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_FLAGS "-fPIC -O3") 
set(CMAKE_BUILD_TYPE "Release")

# INSTALLATION SETTINGS
set( LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )

# Make sure we can find ROOT
execute_process(COMMAND root-config --prefix 
                COMMAND tr -d \\n
                RESULT_VARIABLE RC_RESULT 
                OUTPUT_VARIABLE ROOT_PREFIX)

list(APPEND CMAKE_MODULE_PATH "${ROOT_PREFIX}/share/root/cmake")

## Connect ROOT
find_package(ROOT REQUIRED MathMore)
if (ROOT_FOUND)
    message("-- ROOT found!")
    include_directories(${ROOT_INCLUDE_DIRS})
    link_directories(${ROOT_LIBRARY_DIRS})
else()
    message(SEND_ERROR "-- ROOT not found!")
endif()

## Connect BOOST
find_package(Boost REQUIRED)
if (Boost_FOUND)
    include_directories(${BOOST_INCLUDE_DIRS})
    link_directories(${BOOST_LIBRARY_DIRS})
else()
    message(SEND_ERROR "-- BOOST not found!")
endif()

# BUILD CUBATURE LIBRARY
include_directories("cubature")
file(GLOB CUB_INC "cubature/cubature.h")
file(GLOB CUB_SRC "cubature/*cubature.c")

# Find the jpacStyle library
find_library(JSTYLELIB NAMES jpacStyle libjpacStyle 
                       HINTS "$ENV{JPACSTYLE}/lib")
if (JSTYLELIB)
    message("-- jpacStyle found! (${JSTYLELIB})")
    include_directories("$ENV{JPACSTYLE}/include")
    link_libraries(${JSTYLELIB})
else()
    message("-- jpacStyle not found! Executables will not be installed.")
endif()

# BUILD LIBRARY FROM LOCAL FiLES
include_directories("include")
include_directories("src")
file(GLOB_RECURSE INC "include/*.hpp")
file(GLOB_RECURSE SRC "src/*.cpp")
add_library( jpacTriangle SHARED ${INC} ${CUB_INC} ${SRC} ${CUB_SRC})

# INSTALLATION SETTINGS
set( LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )

install(TARGETS jpacTriangle
        LIBRARY DESTINATION "${LIBRARY_OUTPUT_DIRECTORY}" )

# GRAB EXECUTABLES
include_directories("executables")
file(GLOB EXE_FILES "executables/*.cpp")
foreach( exefile ${EXE_FILES} )
get_filename_component( exename ${exefile} NAME_WE)
add_executable( ${exename} ${exefile} )
target_link_libraries( ${exename} jpacTriangle)
target_link_libraries( ${exename} ${JSTYLELIB})
target_link_libraries( ${exename} ${ROOT_LIBRARIES})
target_link_libraries( ${exename} ${BOOST_LIBRARIES})
endforeach( exefile ${EXE_FILES} )
