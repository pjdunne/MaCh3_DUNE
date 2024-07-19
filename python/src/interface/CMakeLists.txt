pybind11_add_module(_pyMaCh3 SHARED
                                pyMaCh3.cpp
                                MaCh3Instance.cpp)


target_include_directories(_pyMaCh3 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(_pyMaCh3 PRIVATE
                                     MaCh3 SamplePDFDUNE splinesDUNE pybind11::module ROOT::ROOT ${Python3_INCLUDE_DIRS})

install(TARGETS _pyMaCh3 
        LIBRARY DESTINATION ${MACH3_PYTHONPATH})

# Install Python Modules somewhere usable
install(DIRECTORY ${CMAKE_SOURCE_DIR}/python/src/interface/modules/ DESTINATION ${MACH3_PYTHONPATH})