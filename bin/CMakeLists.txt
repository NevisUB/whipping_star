include_directories(${CMAKE_SOURCE_DIR}/inc)
include_directories( ${EIGEN3_INCLUDE_DIR} )

#add_executable (sbnfit_make_covariance sbnfit_make_covariance.cxx)
#target_link_libraries (sbnfit_make_covariance SBNfitlib tinyxmllib EventWeight)

#add_executable (sbnfit_plot_covariance sbnfit_plot_covariance.cxx)
#target_link_libraries (sbnfit_plot_covariance SBNfitlib tinyxmllib EventWeight)

#add_executable (sbnfit_conditional_constraint sbnfit_conditional_constraint.cxx)
#target_link_libraries (sbnfit_conditional_constraint SBNfitlib tinyxmllib EventWeight)

#add_executable (sbnfit_hive_plot sbnfit_hive_plot.cxx)
#target_link_libraries (sbnfit_hive_plot SBNfitlib tinyxmllib EventWeight)

#add_executable (sbnfit_fix_fractional sbnfit_fix_fractional.cxx)
#target_link_libraries (sbnfit_fix_fractional SBNfitlib tinyxmllib EventWeight)

#add_executable (sbnfit_flat_fractional sbnfit_flat_fractional.cxx)
#target_link_libraries (sbnfit_flat_fractional SBNfitlib tinyxmllib EventWeight)

#Davio's DL code, not on master
#add_executable (DLtest DLDis.cc)
#target_link_libraries (DLtest SBNfitlib tinyxmllib EventWeight)


#add_executable (sbnfit_merge_fractional sbnfit_merge_fractional.cxx)
#target_link_libraries (sbnfit_merge_fractional SBNfitlib tinyxmllib EventWeight)

#add_executable (sbnfit_make_spec sbnfit_make_spec.cxx)
#target_link_libraries (sbnfit_make_spec SBNfitlib tinyxmllib EventWeight)

add_executable (sbnfit_wilks sbnfit_wilks.cxx)
target_link_libraries (sbnfit_wilks SBNfitlib sbnfit_tinyxmllib EventWeight)

#add_executable (sbnfit_feldman_cousins sbnfit_feldman_cousins.cxx)
#target_link_libraries (sbnfit_feldman_cousins SBNfitlib tinyxmllib EventWeight)

#add_executable (sbnfit_fraction_fit sbnfit_fraction_fit.cxx)
#target_link_libraries (sbnfit_fraction_fit SBNfitlib tinyxmllib EventWeight)

#add_executable (sbnfit_fraction_fit_v1 sbnfit_fraction_fit_v1.cxx)
#target_link_libraries (sbnfit_fraction_fit_v1 SBNfitlib tinyxmllib EventWeight)
#add_executable (sbnfit_fraction_fit_v2 sbnfit_fraction_fit_v2.cxx)
#target_link_libraries (sbnfit_fraction_fit_v2 SBNfitlib tinyxmllib EventWeight)

#add_executable (sbnfit_uboone_scaling_fc sbnfit_uboone_scaling_fc.cxx)
#target_link_libraries (sbnfit_uboone_scaling_fc SBNfitlib tinyxmllib EventWeight)

#add_executable (sbnfit_lee_frequentist_study sbnfit_lee_frequentist_study.cxx)
#target_link_libraries (sbnfit_lee_frequentist_study SBNfitlib tinyxmllib EventWeight)

#add_executable (sbnfit_scale_spec sbnfit_scale_spec.cxx)
#target_link_libraries (sbnfit_scale_spec SBNfitlib tinyxmllib EventWeight)

#add_executable (sbnfit_pi0_tests sbnfit_pi0_tests.cxx)
#target_link_libraries (sbnfit_pi0_tests SBNfitlib tinyxmllib EventWeight)

