set(SamplePDFDUNE_implementation_files 
	samplePDFDUNEBase.cpp) 
        

set(SamplePDFDUNE_header_files 
	samplePDFDUNEBase.h)

add_library(SamplePDFDUNE SHARED ${SamplePDFDUNE_implementation_files})

target_include_directories(SamplePDFDUNE PRIVATE $ENV{MACH3CORE})
target_compile_features(SamplePDFDUNE PUBLIC)

if(NOT ${CPU_ONLY})
  set_target_properties(SamplePDFDUNE PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  add_dependencies(SamplePDFDUNE SamplePDF Splines Covariance probGpu)
else()
  add_dependencies(SamplePDFDUNE SamplePDF Splines Covariance)
endif()

install(TARGETS SamplePDFDUNE DESTINATION lib)
