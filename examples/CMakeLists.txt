include_directories(${PROJECT_SOURCE_DIR}/inc)

add_executable (example1 example_1.cxx)
target_link_libraries (example1 SBNfitlib sbnfit_tinyxmllib EventWeight)

add_executable (example2 example_2.cxx)
target_link_libraries (example2 SBNfitlib sbnfit_tinyxmllib EventWeight)

add_executable (example3 example_3.cxx)
target_link_libraries (example3 SBNfitlib sbnfit_tinyxmllib EventWeight)

add_executable (exampletest example_test.cxx)
target_link_libraries (exampletest SBNfitlib sbnfit_tinyxmllib EventWeight)


#add_executable (exampleI example_i.cxx)
#target_link_libraries (exampleI SBNfitlib sbnfit_tinyxmllib EventWeight)

#add_executable (exampleI2 example_i2.cxx)
#target_link_libraries (exampleI2 SBNfitlib sbnfit_tinyxmllib EventWeight)

#add_executable (example4 example_4.cxx)
#target_link_libraries (example4 SBNfitlib sbnfit_tinyxmllib EventWeight)

file(COPY ${PROJECT_SOURCE_DIR}/xml/example.xml DESTINATION ${PROJECT_SOURCE_DIR}/build/examples/)




