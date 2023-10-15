
#add_library (SBNfitlib SBNcovariance.cxx SBNgenerate SBNconfig.cxx SBNspec.cxx SBNfit.cxx SBNchi.cxx prob.cxx SBNfit3pN.cxx SBNosc.cxx SBNprob.cxx SBNcls.cxx SBNfeld.cxx )
#add_library (SBNfitlib SBNcovariance.cxx SBNgenerate SBNconfig.cxx SBNspec.cxx SBNfit.cxx SBNchi.cxx prob.cxx SBNfit3pN.cxx SBNosc.cxx SBNprob.cxx SBNcls.cxx SBNfeld.cxx HIVEplot.cxx SBNconditional.cxx bayes.cxx)
add_library (SBNfitlib SBNcovariance.cxx SBNgenerate SBNconfig.cxx SBNspec.cxx SBNfit.cxx SBNchi.cxx prob.cxx SBNfit3pN.cxx SBNosc.cxx SBNprob.cxx SBNcls.cxx SBNfeld.cxx SBNconditional.cxx bayes.cxx)

include_directories ( ${CMAKE_SOURCE_DIR}/inc/)
include_directories ( ${CMAKE_SOURCE_DIR}/src/)
include_directories ( ${CMAKE_SOURCE_DIR}/tinyxml/)
include_directories( ${EIGEN3_INCLUDE_DIR} )

#ROOT_GENERATE_DICTIONARY(LINKDEF mdict.h )
tarGet_link_libraries(SBNfitlib ${ROOT_LIBRARIES})
tarGet_link_libraries(SBNfitlib gsl gslcblas)
tarGet_link_libraries(SBNfitlib ${GSL_LIBRARIES})
tarGet_link_libraries(SBNfitlib tinyxmllib)
