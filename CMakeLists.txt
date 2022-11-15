cmake_minimum_required (VERSION 3.11)
project (SBNfit)
# The version number.
set (SBNfit_VERSION_MAJOR 3)
set (SBNfit_VERSION_MINOR 0)

if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
          set (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
    endif ()
    else ()
	    set (CMAKE_CXX_STANDARD 17)
endif ()

message("@main CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
unset (CMAKE_CXX_FLAGS)

option(USE_GPU "USE GPU" OFF)
if(USE_GPU)
	message("")
	message("Using GPU")
	message("")
	add_compile_definitions(USE_GPU)
else()
	message("")
	message("Not using GPU")
	message("")
endif()



FIND_PACKAGE(ROOT COMPONENTS MathCore MathMore)# COMPONENTS Core Cint RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread)

if(ROOT_FOUND)
	message("Found ROOT using CMAKE FIND_PACKAGE")
	include(${ROOT_USE_FILE})
ELSE()
	message("ROOT not compiled using CMAKE, use root-config")
	#Get our ROOT libraries and such
	FIND_PROGRAM(ROOT_CONFIG root-config)
	EXEC_PROGRAM(${ROOT_CONFIG} ARGS --cflags --glibs OUTPUT_VARIABLE ROOT_CXX_FLAGS)
	message("ROOT cxx flags set to: ${ROOT_CXX_FLAGS}")
ENDIF()

#Get our GSL libraries
FIND_PROGRAM(GSL_CONFIG gsl-config)
EXEC_PROGRAM(${GSL_CONFIG} ARGS --cflags OUTPUT_VARIABLE GSL_CXX_FLAGS)
EXEC_PROGRAM(${GSL_CONFIG} ARGS --libs OUTPUT_VARIABLE GSL_LIBS)

set(EIGEN_CONFIG "$ENV{EIGEN_INC}")

message("@main CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
unset (CMAKE_CXX_FLAGS)

message("GSL_LIBS: ${GSL_LIBS}")
message("GSL_CXX_FLAGS: ${GSL_CXX_FLAGS}")

if(USE_GPU)
    set(CMAKE_MAGMA_CXX_FLAGS "-$ENV{MAGMA_CFLAGS} $ENV{MAGMA_LIB}")
    message("@main CMAKE_MAGMA_CXX_FLAGS: ${CMAKE_MAGMA_CXX_FLAGS}")
    set(CMAKE_GPU_CXX_FLAGS "-ta=tesla:nollvm -Mcudalib=curand -Minfo=accel -acc -w ${CMAKE_MAGMA_CXX_FLAGS}" )
ENDIF()

#gopt is forcing LLVM compilation
# -noswitcherror only for pgcc
#set (CMAKE_SO_CXX_FLAGS "-noswitcherror ${ROOT_CXX_FLAGS} ${GSL_LIBS} ${CMAKE_CXX_FLAGS} -w -Wno-deprecated-declarations")
set (CMAKE_SO_CXX_FLAGS " ${ROOT_CXX_FLAGS} ${GSL_LIBS} ${GSL_CXX_FLAGS} ${CMAKE_CXX_FLAGS} -w -Wno-deprecated-declarations")
set (CMAKE_ACC_CXX_FLAGS "${CMAKE_GPU_CXX_FLAGS} ${CMAKE_SO_CXX_FLAGS} -g -I ${EIGEN_CONFIG}")
set (CMAKE_PCC_CXX_FLAGS "-ta=multicore -Minfo=accel -acc -w ${CMAKE_SO_CXX_FLAGS}")

set (CMAKE_CXX_FLAGS "${CMAKE_ACC_CXX_FLAGS}")

message("@main CMAKE_SO_CXX_FLAGS: ${CMAKE_SO_CXX_FLAGS}")
message("@main CMAKE_ACC_CXX_FLAGS: ${CMAKE_ACC_CXX_FLAGS}")
message("@main CMAKE_ACC_PXX_FLAGS: ${CMAKE_PCC_CXX_FLAGS}")
message("@main CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

include_directories(${PROJECT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inc/)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/tinyxml/)

add_subdirectory(libio)

unset(CMAKE_CXX_FLAGS)
set(CMAKE_CXX_FLAGS "${CMAKE_ACC_CXX_FLAGS}")
message("@main CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

add_subdirectory(tinyxml)
add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(bin)
add_subdirectory(unit)
add_subdirectory(HPCfunctionality)

