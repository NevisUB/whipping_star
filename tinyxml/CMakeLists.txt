# The extension is already found. Any number of sources could be listed here.
add_library (sbnfit_tinyxmllib  tinystr.cpp tinyxml.cpp tinyxmlerror.cpp tinyxmlparser.cpp)

# Make sure the compiler can find include files for our Hello library
# when other libraries or executables link to Hello
#target_include_directories (sbnfit_tinyxmllib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
include_directories ( ${CMAKE_CURRENT_SOURCE_DIR})
