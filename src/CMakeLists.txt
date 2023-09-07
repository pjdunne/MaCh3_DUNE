set(APPS EventRatesDUNE SigmaVariationsDUNE jointFitDUNE oscPlots test Oscillator testAtmDUNE DUNE_LLHscan_osc DUNE_LLHscan_systs DUNE_atm_chain)

foreach(a ${APPS})
	  add_executable(${a} ${a}.cpp)
	  target_link_libraries(${a} SamplePDFDUNE)
	  install(TARGETS ${a} DESTINATION bin)
endforeach()
