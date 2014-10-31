cmake_minimum_required (VERSION 2.8.11)
project (YarnProject)


# General Configuration: find all the required libraries.
	list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

#===================================================================
# compile option
option(USE_DEBUG "Turn on the debug mode" OFF)
option(USE_64BIT_ARCH "Compile the 64bit executable" OFF)

# Options for building individual tools
option(BUILD_RODSOUND         	"Compile the RodSound Target" ON)
option(BUILD_SIMULATOR     		"Compile the Simulator Target" ON)
option(BUILD_VISUALIZER         "Compile the Visualizer Target" OFF)

#===================================================================
## Compiler
#===================================================================
# set compiler flags for debug/release
if ( USE_DEBUG )
    MESSAGE(STATUS "Using debug mode")
    add_definitions(-DDEBUG)
    set(CMAKE_BUILD_TYPE Debug)
else ( USE_DEBUG )
    MESSAGE(STATUS "Using release mode")
    set(CMAKE_BUILD_TYPE Release)
endif ( USE_DEBUG )

#===================================================================
## Libraries
#===================================================================

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

include_directories(include)

add_subdirectory(src)

IF(APPLE)
	include(CheckCXXCompilerFlag)
	CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
	if(COMPILER_SUPPORTS_CXX11)
	    set(CMAKE_CXX_FLAGS_RELEASE "-std=c++11 ${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "" FORCE)
	    message(STATUS ${CMAKE_CXX_FLAGS_RELEASE})
	else()
	    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
	endif()
	set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++0x")
	set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
else()
	include(CheckCXXCompilerFlag)
	CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
	if(COMPILER_SUPPORTS_CXX11)
	    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
	else()
	        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
	endif()
endif()

# Recurse into the "Hello" and "Demo" subdirectories. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.
# add_subdirectory (src)