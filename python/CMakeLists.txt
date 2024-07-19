add_subdirectory(src)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/python/apps/ DESTINATION ${MACH3_PYTHONPATH}/..)