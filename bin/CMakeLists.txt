include_directories(${CMAKE_SOURCE_DIR}/inc)

add_executable (sbnfit_make_covariance sbnfit_make_covariance.cxx)
target_link_libraries (sbnfit_make_covariance SBNfitlib tinyxmllib EventWeight)

add_executable (sbnfit_feldman_cousins sbnfit_feldman_cousins.cxx)
target_link_libraries (sbnfit_feldman_cousins SBNfitlib tinyxmllib EventWeight)

add_executable (sbnfit_make_spec sbnfit_make_spec.cxx)
target_link_libraries (sbnfit_make_spec SBNfitlib tinyxmllib EventWeight)

add_executable (sbnfit_uboone_scaling_fc sbnfit_uboone_scaling_fc.cxx)
target_link_libraries (sbnfit_uboone_scaling_fc SBNfitlib tinyxmllib EventWeight)


add_executable (sbnfit_lee_frequentist_study sbnfit_lee_frequentist_study.cxx)
target_link_libraries (sbnfit_lee_frequentist_study SBNfitlib tinyxmllib EventWeight)

add_executable (sbnfit_scale_spec sbnfit_scale_spec.cxx)
target_link_libraries (sbnfit_scale_spec SBNfitlib tinyxmllib EventWeight)


