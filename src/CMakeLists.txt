set(APPS EventRatesDUNE SigmaVariationsDUNE LLHScansDUNE jointFitDUNE)

foreach(a ${APPS})
	  add_executable(${a} ${a}.cpp)
	  target_link_libraries(${a} SamplePDFDUNE)
	  install(TARGETS ${a} DESTINATION bin)
endforeach()

add_executable(LukeMessAround LukeMessAround.cpp)
    target_link_libraries(LukeMessAround yaml-cpp Eigen3::Eigen SamplePDFDUNE TreePlayer)
    install(TARGETS LukeMessAround DESTINATION bin)