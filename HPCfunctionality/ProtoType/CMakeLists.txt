
option                      (diy_thread     "Enable diy threading"                 OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-march=native -mtune=native -std=gnu++17 -Wall -O2 -DNDEBUG -L${HDF5_LIB_DIRS} -L$ENV{SBNFIT_LIBDIR}")# -fPermissive")

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


include(CheckIPOSupported)
check_ipo_supported(RESULT supported OUTPUT error)

if( supported )
    message(STATUS "IPO / LTO enabled")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    message(STATUS "IPO / LTO not supported: <${error}>")
endif()


## AddressSanitize
set(CMAKE_C_FLAGS_ASAN
    "-fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g -O1"
    CACHE STRING "Flags used by the C compiler during AddressSanitizer builds."
    FORCE)
set(CMAKE_CXX_FLAGS_ASAN
    "-fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g -O1"
    CACHE STRING "Flags used by the C++ compiler during AddressSanitizer builds."
    FORCE)

 #LeakSanitizer
set(CMAKE_C_FLAGS_LSAN
    "-fsanitize=leak -fno-omit-frame-pointer -g -O1"
    CACHE STRING "Flags used by the C compiler during LeakSanitizer builds."
    FORCE)
set(CMAKE_CXX_FLAGS_LSAN
    "-fsanitize=leak -fno-omit-frame-pointer -g -O1"
    CACHE STRING "Flags used by the C++ compiler during LeakSanitizer builds."
    FORCE)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)



set(CMAKE_CXX_STANDARD 17)
 # MPI
#find_package                (MPI REQUIRED)
    set                     (libraries ${libraries} ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})

#find_package(MPI REQUIRED)

set (LOCAL_MPI_DIR "/uboone/app/users/markrl/programs/MPICH/mpich-install")
set (LOCAL_MPI_INC "${LOCAL_MPI_DIR}/include") 
set (LOCAL_MPI_LIB "${LOCAL_MPI_DIR}/lib") 


#find_package                (HDF5)
set                     (libraries ${libraries} ${HDF5_C_LIBRARIES})
#set                     (libraries ${libraries} ${HDF5_C_LIBRARIES} ${NOMAD_INCLUDE_DIRS}/../lib/libnomad.so)

find_package(GSL REQUIRED)
set                     (libraries ${libraries} ${GSL_LIBRARIES})
FIND_PROGRAM(GSL_CONFIG gsl-config)
EXEC_PROGRAM(${GSL_CONFIG} ARGS --cflags OUTPUT_VARIABLE GSL_CXX_FLAGS)
EXEC_PROGRAM(${GSL_CONFIG} ARGS --libs OUTPUT_VARIABLE GSL_LIBS)



add_definitions         (-DDIY_NO_THREADS)
# DIY
find_path                   (DIY_INCLUDE_DIRS diy/types.h)
#find_path                   (NOMAD_INCLUDE_DIRS nomad.hpp)
#find_path                   (SGTE_INCLUDE_DIRS sgtelib.hpp)

set (MY_PUBLIC_HEADERS 
    ${CMAKE_CURRENT_SOURCE_DIR}/include/GenericBlock.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/opts.h
    )

# Include dirs
set                         (CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem")


include_directories         (${CMAKE_CURRENT_SOURCE_DIR}/include
                            ${DIY_INCLUDE_DIRS}
                            ${GSL_INCLUDE_DIRS}
                            ${LOCAL_MPI_DIR}/include
                            SYSTEM ${MPI_INCLUDE_PATH}
                            ${HDF5_INCLUDE_DIRS}
                            ${EIGEN3_INCLUDE_DIR}
                            ${HIGHFIVE_DIR}/include 
                            )
                            #${SGTE_INCLUDE_DIRS}
                            #${NOMAD_INCLUDE_DIRS}
                        
message(STATUS "GARR Lbraries, check HDF5: ${libraries}")
message(STATUS "GARR HDF5 include dirs, check HDF5: ${HDF5_INCLUDE_DIRS}")
message(STATUS "CMAKE_MODULE_PATH : ${CMAKE_MODULE_PATH}")

##set (ROOT_LIBS "-I${ROOTSYS}/include -L${ROOTSYS}/lib -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTDataFrame -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc")
#Remove some niche ones from above. ? not sure useful. can replace with root-config?
set (ROOT_LIBS "-I${ROOTSYS}/include -L${ROOTSYS}/lib -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad  -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc")


# -lstdc++fs is for filesystem
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  ${GSL_LIBS} ${GSL_CXX_FLAGS} ${ROOT_LIBS}  -lmpi -L${LOCAL_MPI_LIB} -I${LOCAL_MPI_INC} -I${HDF5_INCLUDE_DIRS} -L${HDF5_LIB_DIRS}  -I/uboone/app/users/markrl/programs/HighFive/HighFive/include -I/uboone/app/users/markrl/programs/DIY/diy/include -L/nashome/i/isafa/usr/share/eigen3/   -Wunused-variable -Wno-deprecated-declarations  -Wnarrowing -Wreturn-type  "  )

message(STATUS "CXX_FLAGS FINAL PLACE HERE YES LOOK: ${CMAKE_CXX_FLAGS}")

# Figure out correct library suffix
if                          (APPLE)
  set                       (CMAKE_SHARED_LIBRARY_SUFFIX ".dylib")
  message                   ("Shared library suffix set to dylib")
else                        ()
  set                       (CMAKE_SHARED_LIBRARY_SUFFIX ".so")
  message                   ("Shared library suffix set to so")
endif                       ()


# Test if user has supplied -DHIGHFIVE_DIR                        
if (NOT HIGHFIVE_DIR)
    message                 (FATAL_ERROR " Need HighFive in order to build program --- check -DHIGHFIVE_DIR")
endif                       ()

if (NOT ROOTSYS)
    message                 (FATAL_ERROR " Unfortunately need ROOT in order to build program --- check -DROOTSYS")
endif                       ()


find_file(HAVE_H5_H
    NAMES H5File.hpp
    HINTS "${HIGHFIVE_DIR}/include/highfive" 
    )
if(NOT HAVE_H5_H)
    message(FATAL_ERROR "Specified HighFive dir not found --- check -DHIGHFIVE_DIR")
endif()

# If all found, add build target
add_subdirectory        (src)

set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "6")
set(CPACK_PACKAGE_VERSION_PATCH "6")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_PACKAGE_FILE_NAME
  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(CPACK_SOURCE_IGNORE_FILES
    "Docker;benchmark;test;.git;~;build;${CPACK_SOURCE_IGNORE_FILES}")
include(CPack)

add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
