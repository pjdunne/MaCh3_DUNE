set(HEADERS
  samplePDFDUNEAtmBase.h
  samplePDFDUNEBeamFDBase.h
  samplePDFDUNEBeamNDBase.h
  samplePDFDUNEBeamNDGarBase.h
)

add_library(SamplePDFDUNE SHARED
  samplePDFDUNEAtmBase.cpp
  samplePDFDUNEBeamFDBase.cpp
  samplePDFDUNEBeamNDBase.cpp
  samplePDFDUNEBeamNDGarBase.cpp
)

set_target_properties(SamplePDFDUNE PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME SamplePDFDUNE)

if(NOT CPU_ONLY)
  set_target_properties(SamplePDFDUNE PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()

target_link_libraries(SamplePDFDUNE MaCh3::All MaCh3DUNECompilerOptions duneanaobj::All)

target_include_directories(SamplePDFDUNE PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
  $<INSTALL_INTERFACE:include>)

install(TARGETS SamplePDFDUNE
        EXPORT mach3dune-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/samplePDF)

add_library(MACH3DUNE::SamplePDFDUNE ALIAS SamplePDFDUNE)
