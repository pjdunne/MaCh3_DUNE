add_subdirectory(interface)

# Install Python Modules somewhere usable
install(DIRECTORY ${CMAKE_SOURCE_DIR}/python/src/modules/ DESTINATION ${MACH3_PYTHONPATH})