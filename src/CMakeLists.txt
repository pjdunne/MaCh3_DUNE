set(APPS EventRatesDUNE SigmaVariationsDUNE jointFitDUNE oscPlots test Oscillator testAtmDUNE DUNE_LLHscan_osc)

foreach(a ${APPS})
	  add_executable(${a} ${a}.cpp)
	  target_link_libraries(${a} SamplePDFDUNE)
	  install(TARGETS ${a} DESTINATION bin)
endforeach()
