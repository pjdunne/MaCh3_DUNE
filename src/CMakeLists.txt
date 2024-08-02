set(APPS EventRatesDUNE SigmaVariationsDUNE LLHScansDUNE jointFitDUNE Selections)

foreach(a ${APPS})
	  add_executable(${a} ${a}.cpp)
	  target_link_libraries(${a} SamplePDFDUNE)
	  install(TARGETS ${a} DESTINATION bin)
endforeach()
