include_directories(${PROJECT_SOURCE_DIR}/inc)

add_executable (unit1 unit1.cxx)
target_link_libraries (unit1 SBNfitlib sbnfit_tinyxmllib EventWeight)

add_executable (unit2 unit2.cxx)
target_link_libraries (unit2 SBNfitlib sbnfit_tinyxmllib EventWeight)



file(COPY ${PROJECT_SOURCE_DIR}/unit/unit1a.root DESTINATION ${PROJECT_SOURCE_DIR}/build/unit/)
file(COPY ${PROJECT_SOURCE_DIR}/unit/unit1b.root DESTINATION ${PROJECT_SOURCE_DIR}/build/unit/)
file(COPY ${PROJECT_SOURCE_DIR}/unit/unit1c.root DESTINATION ${PROJECT_SOURCE_DIR}/build/unit/)

file(COPY ${PROJECT_SOURCE_DIR}/unit/unit1a_matrix.root DESTINATION ${PROJECT_SOURCE_DIR}/build/unit/)
file(COPY ${PROJECT_SOURCE_DIR}/unit/unit1b_matrix.root DESTINATION ${PROJECT_SOURCE_DIR}/build/unit/)
file(COPY ${PROJECT_SOURCE_DIR}/unit/unit1c_matrix.root DESTINATION ${PROJECT_SOURCE_DIR}/build/unit/)

file(COPY ${PROJECT_SOURCE_DIR}/xml/unit/unit1a.xml DESTINATION ${PROJECT_SOURCE_DIR}/build/unit/)
file(COPY ${PROJECT_SOURCE_DIR}/xml/unit/unit1b.xml DESTINATION ${PROJECT_SOURCE_DIR}/build/unit/)
file(COPY ${PROJECT_SOURCE_DIR}/xml/unit/unit1c.xml DESTINATION ${PROJECT_SOURCE_DIR}/build/unit/)

