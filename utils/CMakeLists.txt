add_executable(utils_makeCAFs makeCAFsNDGAr.cxx)
target_include_directories(utils_makeCAFs PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../duneanaobj/)
#target_include_directories(utils_makeCAFs PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(utils_makeCAFs duneanaobj)

install(TARGETS utils_makeCAFs DESTINATION lib)
