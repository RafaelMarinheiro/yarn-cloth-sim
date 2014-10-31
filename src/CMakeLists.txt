#===================================================================
## RodSound
set(RODSOUND_APP_HPP "${CMAKE_CURRENT_SOURCE_DIR}/RodSound/Headers/RodSoundApp.h")
set(RODSOUND_APP_CPP "${CMAKE_CURRENT_SOURCE_DIR}/RodSound/Source/RodSoundApp.cpp")
set(RODSOUND_MAIN_CPP "${CMAKE_CURRENT_SOURCE_DIR}/RodSound/Source/main.cpp")

file(GLOB RODSOUND_HEADERS "RodSound/Headers/*.h")
file(GLOB RODSOUND_SOURCE "RodSound/Source/*.cpp")

set(RODSOUND_HEADERS_DIR "RodSound/Headers/")

list(REMOVE_ITEM RODSOUND_HEADERS ${RODSOUND_APP_HPP})
list(REMOVE_ITEM RODSOUND_SOURCE ${RODSOUND_APP_CPP})
list(REMOVE_ITEM RODSOUND_SOURCE ${RODSOUND_MAIN_CPP})

## Simulator
set(SIMULATOR_APP_HPP "${CMAKE_CURRENT_SOURCE_DIR}/Simulator/Headers/SimulatorApp.h")
set(SIMULATOR_APP_CPP "${CMAKE_CURRENT_SOURCE_DIR}/Simulator/Source/SimulatorApp.cpp")
set(SIMULATOR_MAIN_CPP "${CMAKE_CURRENT_SOURCE_DIR}/Simulator/Source/main.cpp")

file(GLOB SIMULATOR_HEADERS "Simulator/Headers/*.h")
file(GLOB SIMULATOR_SOURCE "Simulator/Source/*.cpp")

list(REMOVE_ITEM SIMULATOR_HEADERS ${SIMULATOR_APP_HPP})
list(REMOVE_ITEM SIMULATOR_SOURCE ${SIMULATOR_APP_CPP})
list(REMOVE_ITEM SIMULATOR_SOURCE ${SIMULATOR_MAIN_CPP})

set(SIMULATOR_HEADERS_DIR "Simulator/Headers/")

## Visualizer
set(VISUALIZER_APP_HPP "Visualizer/Headers/VisualizerApp.h")
set(VISUALIZER_APP_CPP "Visualizer/Source/VisualizerApp.cpp")

file(GLOB VISUALIZER_HEADERS "Visualizer/Headers/*.h")
file(GLOB VISUALIZER_SOURCE "Visualizer/Source/*.cpp")

list(REMOVE_ITEM VISUALIZER_HEADERS ${VISUALIZER_APP_HPP})
list(REMOVE_ITEM VISUALIZER_SOURCE ${VISUALIZER_APP_CPP})

set(VISUALIZER_HEADERS_DIR "Visualizer/Headers/")

## Common 
file(GLOB COMMON_HEADERS "Common/Headers/*.h")
file(GLOB COMMON_SOURCE "Common/Source/*.cpp")
set(COMMON_HEADERS_DIR "Common/Headers/")

#===================================================================

IF (BUILD_RODSOUND)
	FIND_PACKAGE(Qt4 COMPONENTS QtCore QtXml QtOpenGL QtGui REQUIRED)
	find_package(QGLViewer REQUIRED)
	find_package(OpenGL REQUIRED)

	find_package(Boost COMPONENTS system filesystem REQUIRED)

	set(QT_USE_QTXML 1)
	set(QT_USE_QTOPENGL 1)
	set(CMAKE_AUTOMOC ON)
	include(${QT_USE_FILE})

	add_executable(RodSound
					${RODSOUND_MAIN_CPP}
					${RODSOUND_HEADERS} ${RODSOUND_SOURCE}
					${COMMON_HEADERS}	 ${COMMON_SOURCE}
				    ${SIMULATOR_HEADERS} ${SIMULATOR_SOURCE}
				   # ${VISUALIZER_HEADERS} ${VISUALIZER_SOURCE}
				   )

	target_include_directories(RodSound
								PUBLIC ${QT_USE_FILE}
								PUBLIC ${QGLViewer_INCLUDE_DIR}
								PUBLIC ${COMMON_HEADERS_DIR}
								PUBLIC ${SIMULATOR_HEADERS_DIR}
								PUBLIC ${RODSOUND_HEADERS_DIR}
								)

	target_link_libraries(RodSound
							${Boost_LIBRARIES}
							${QT_LIBRARIES}
					        ${QGLViewer_LIBRARY}
					        ${OPENGL_LIBRARIES}
					        )
ENDIF(BUILD_RODSOUND)

IF (BUILD_SIMULATOR)
	FIND_PACKAGE(Qt4 COMPONENTS QtCore QtXml QtOpenGL QtGui REQUIRED)
	find_package(QGLViewer REQUIRED)
	find_package(OpenGL REQUIRED)

	find_package(Boost COMPONENTS system filesystem REQUIRED)

	set(QT_USE_QTXML 1)
	set(QT_USE_QTOPENGL 1)
	set(CMAKE_AUTOMOC ON)
	include(${QT_USE_FILE})

	add_executable(Simulator
					# ${RODSOUND_HEADERS} ${RODSOUND_SOURCE}
					${SIMULATOR_MAIN_CPP}
					${COMMON_HEADERS}	 ${COMMON_SOURCE}
				    ${SIMULATOR_HEADERS} ${SIMULATOR_SOURCE}
				    # ${VISUALIZER_HEADERS} ${VISUALIZER_SOURCE}
				   )

	# add_precompiled_header(Simulator "Simulator/Headers/Simulator_Prefix.pch")
	target_include_directories(Simulator
								PUBLIC ${QT_USE_FILE}
								PUBLIC ${QGLViewer_INCLUDE_DIR}
								PUBLIC ${COMMON_HEADERS_DIR}
								PUBLIC ${SIMULATOR_HEADERS_DIR}
								)

	target_link_libraries(Simulator
							${Boost_LIBRARIES}
							${QT_LIBRARIES}
					        ${QGLViewer_LIBRARY}
					        ${OPENGL_LIBRARIES}
					        )
ENDIF(BUILD_SIMULATOR)

# IF (BUILD_VISUALIZER)
# 	add_executable(Visualizer ${VISUALIZER_APP_HPP} ${VISUALIZER_APP_CPP}
# 				   ${RODSOUND_HEADERS} ${RODSOUND_SOURCE}
# 				   ${SIMULATOR_HEADERS} ${SIMULATOR_SOURCE}
# 				   ${VISUALIZER_HEADERS} ${VISUALIZER_SOURCE}
# 				   )

# 	target_include_directories(Visualizer PUBLIC ${RODSOUND_HEADERS_DIR} ${SIMULATOR_HEADERS_DIR} ${VISUALIZER_HEADERS_DIR}
# 								${CINDER_INCLUDE_DIRS})
# 	target_link_libraries(Visualizer ${CINDER_LIBRARIES})
# 	add_macos_framework(RodSound Cocoa)
# 	add_macos_framework(RodSound OpenGL)
# 	add_macos_framework(RodSound CoreVideo)
# 	add_macos_framework(RodSound QTKit)
# 	add_macos_framework(RodSound Accelerate)
# 	add_macos_framework(RodSound AudioToolbox)
# 	add_macos_framework(RodSound AudioUnit)
# 	add_macos_framework(RodSound CoreAudio)
# ENDIF(BUILD_VISUALIZER)