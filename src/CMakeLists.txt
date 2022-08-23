
add_library (SBNfitlib SBNcovariance.cxx SBNgenerate.cxx SBNconfig.cxx SBNspec.cxx SBNfit.cxx SBNchi.cxx prob.cxx SBNfit3pN.cxx SBNosc.cxx SBNprob.cxx SBNcls.cxx SBNfeld.cxx HIVEplot.cxx SBNconditional.cxx SBNsinglephoton.cxx bayes.cxx)

# state that anybody linking to us needs to include the current source dir
# to find header file, while we don't.
target_include_directories(SBNfitlib
          INTERFACE ${PROJECT_SOURCE_DIR}/inc/
          INTERFACE ${PROJECT_SOURCE_DIR}/src/
          )

include_directories ( ${PROJECT_SOURCE_DIR}/inc/)
include_directories ( ${PROJECT_SOURCE_DIR}/src/)
include_directories ( ${PROJECT_SOURCE_DIR}/tinyxml/)

#ROOT_GENERATE_DICTIONARY(LINKDEF mdict.h )
tarGet_link_libraries(SBNfitlib ${ROOT_LIBRARIES})
tarGet_link_libraries(SBNfitlib gsl gslcblas)
tarGet_link_libraries(SBNfitlib ${GSL_LIBRARIES})
tarGet_link_libraries(SBNfitlib sbnfit_tinyxmllib)
