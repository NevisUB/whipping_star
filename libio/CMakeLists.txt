include(${ROOT_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR}/inc)

unset(CMAKE_CXX_FLAGS)
set(CMAKE_CXX_FLAGS " ${CMAKE_SO_CXX_FLAGS}")
message("@libio CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

#--- create eventweight library
ROOT_GENERATE_DICTIONARY(G__EventWeight eventweight.h LINKDEF eventweight_LinkDef.h)

#---Create a shared library with geneated dictionary
add_library(EventWeight SHARED eventweight.cxx G__EventWeight.cxx)
target_link_libraries(EventWeight ${ROOT_LIBRARIES})


#--- create DeNan library ---
ROOT_GENERATE_DICTIONARY(DeNan_Dict denan.h LINKDEF denan_LinkDef.h)
add_library(SBNfit_DeNan SHARED denan.cxx DeNan_Dict.cxx)
target_link_libraries(SBNfit_DeNan ${ROOT_LIBRARIES})
