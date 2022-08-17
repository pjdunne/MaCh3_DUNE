set(APPS 
	 SamplePDFDUNE_test)


foreach(a ${APPS})
	  add_executable(${a} ${a}.cpp)

	  #set_property(TARGET ${a} PROPERTY CUDA_SEPARABLE_COMPILATION ON)
	  add_dependencies(${a} throwParms SamplePDF SamplePDFDUNE Covariance MCMC Splines Manager)
	  target_include_directories(${a} PRIVATE ${CMAKE_SOURCE_DIR})
          #set_source_files_properties(${a} PROPERTIES COMPILE_OPTIONS)
          link_directories(${CMAKE_BINARY_DIR}/lib)	  
	  target_link_libraries(${a} -L${CMAKE_BINARY_DIR}/lib SamplePDFDUNE throwParms SamplePDF Covariance Splines Manager ThreeProb_1.00 ${CMAKE_BINARY_DIR}/lib/libconfig++.so.9.1.3 ${ROOT_LD_FLAGS})
	  install(TARGETS ${a} DESTINATION bin)
endforeach()



