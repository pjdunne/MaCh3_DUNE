set(splinesDUNE_implementation_files 
  splinesDUNE.cpp) 

set(splinesDUNE_header_files 
	splinesDUNE.h)

add_library(splinesDUNE SHARED ${splinesDUNE_implementation_files})
target_include_directories(splinesDUNE PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(splinesDUNE MaCh3)

install(TARGETS splinesDUNE DESTINATION lib)
