install(FILES ${MY_PUBLIC_HEADERS}
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include
    COMPONENT devel)

add_executable              (pureeigen pureeigen.cpp)

target_link_libraries (pureeigen LINK_PUBLIC ${libraries} ${GSL_LIBS} SBNfitlib Eigen3::Eigen "-lstdc++fs" )

install(TARGETS pureeigen
    DESTINATION ${CMAKE_INSTALL_PREFIX/bin}
    RUNTIME DESTINATION bin
    )
