find_package(GSL REQUIRED)

set(PROGRAMS
  inversionMethod
  rejectionMethod
  metropolisMethod
  metropolisEquilibrium
  semiclassicalQuantizationLJ
  semiclassicalQuantizationMorse)

foreach(program ${PROGRAMS})
  add_executable(${program} ${program}.cpp)
  target_link_libraries(${program} ${CORELIBS})
endforeach(program)

set(PROGRAMS_GSL
  semiclassicalQuantizationLJ-GSL)
foreach(program ${PROGRAMS_GSL})
  add_executable(${program} ${program}.cpp)
  target_link_libraries(${program} ${CORELIBS})
  target_link_libraries(${program} ${GSL_LIBRARIES})
endforeach(program)
