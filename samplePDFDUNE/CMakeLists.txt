set(SamplePDFDUNE_implementation_files
  samplePDFBlob.cpp
	samplePDFDUNEBase.cpp
	samplePDFDUNEBaseND.cpp) 

set(SamplePDFDUNE_header_files 
	samplePDFDUNEBase.h
	samplePDFDUNEBaseND.h)

add_library(SamplePDFDUNE SHARED ${SamplePDFDUNE_implementation_files})
target_include_directories(SamplePDFDUNE PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(SamplePDFDUNE MaCh3 splinesDUNE Eigen3::Eigen TreePlayer)

if(NOT ${CPU_ONLY})
  set_target_properties(SamplePDFDUNE PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  target_link_libraries(SamplePDFDUNE probGpu)
else()
endif()

install(TARGETS SamplePDFDUNE DESTINATION lib)
