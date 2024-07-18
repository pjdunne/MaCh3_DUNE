pybind11_add_module(_pyMaCh3 SHARED MaCh3Instance.cpp)


target_include_directories(_pyMaCh3 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(_pyMaCh3 PRIVATE
                                     MaCh3 SamplePDFDUNE splinesDUNE pybind11::module ROOT::ROOT ${Python3_INCLUDE_DIRS})

install(TARGETS _pyMaCh3 
        LIBRARY DESTINATION ${MACH3_PYTHONPATH})

