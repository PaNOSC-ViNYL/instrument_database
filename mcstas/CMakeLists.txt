cmake_minimum_required(VERSION 3.20.0)
project(mcstas_instruments
  VERSION 0.1.0
  LANGUAGES C
  )

include(GNUInstallDirs)
include(FetchContent)

find_package (Python REQUIRED)

string(TOUPPER ${PROJECT_NAME} component_instruments)
set(component_instruments ${component_instruments})

# mcstas code generates plenty of warnings!!! 
set(CMAKE_C_FLAGS "-O3 -Wno-error -Wno-all -Wno-extra -Wno-unused-result -Wno-format-truncation -Wno-format-overflow -Wno-format")

if(NOT DEFINED ENV{MCSTAS})
  message(FATAL_ERROR "MCSTAS environment not defined")
else()
  set(MCSTAS $ENV{MCSTAS} CACHE PATH "Path of the mcstas installation [MANDATORY]")
endif()
set(PROXY "--proxy=$ENV{http_proxy}" CACHE STRING "http proxy")
############################# MCSTAS
message(STATUS "MCSTAS: ${MCSTAS}")
#include_directories(/usr/share/mcstas/2.6/libs/mcpl)
#link_directories(/usr/share/mcstas/2.6/libs/mcpl/)
include_directories(${MCSTAS}/libs/mcpl)
link_directories(${MCSTAS}/libs/mcpl/)

#------------------------------ Python environment with isolated version of Mcstasscript
add_custom_command(OUTPUT python_env.done
  COMMENT "Creating python virtual environment"
  COMMAND ${Python_EXECUTABLE} -m venv ${CMAKE_CURRENT_BINARY_DIR}/python_env
  COMMAND touch python_env.done
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  VERBATIM
  )

add_custom_command(OUTPUT mcstasscript.done
  COMMENT "Install fresh release of McStasScript"
  DEPENDS python_env.done
  DEPENDS requirements.txt
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/python_env/bin/pip ${PROXY} install -r ${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt
  COMMAND touch mcstasscript.done
  VERBATIM
  )


#------------------------------ Function to convert and compile mcstas instruments
set(ILL_INSTRUMENTS ${CMAKE_CURRENT_SOURCE_DIR}/ILL/instruments/)
set(ILL_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/ILL/sources/)
set(ILL_INSTRUMENTS_BUILD ${CMAKE_CURRENT_BINARY_DIR}/ILL/instruments/)
set(TEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/Test_outputs)
file(MAKE_DIRECTORY ${TEST_DIR})

  
# instrument_file should be relative to ILL_INSTRUMENTS
function(run_mcstas instrument_file)
  cmake_path(SET instrument_file_full_path ${ILL_INSTRUMENTS}/${instrument_file})

  cmake_path(REMOVE_FILENAME instrument_file OUTPUT_VARIABLE dir)

  cmake_path(GET instrument_file FILENAME outname)
  cmake_path(REMOVE_EXTENSION outname)

  set(TEST_REFERENCE_DIR ${ILL_INSTRUMENTS}/${dir}/test/)

  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${dir})
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${dir}/${outname}.c
  	COMMENT "creation of ${outname}.c for version: ${version}"
  	DEPENDS mcstasscript.done
	DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/scripts/py_to_instr.py
	DEPENDS ${instrument_file_full_path}
	DEPENDS ${ILL_SOURCES} # any source file, trigger the recompilation of all the instruments...
  	COMMAND ${CMAKE_CURRENT_BINARY_DIR}/python_env/bin/python3 ${CMAKE_CURRENT_SOURCE_DIR}/scripts/py_to_instr.py ${instrument_file_full_path}
  	COMMAND mcstas -t -o ${outname}.c  ${outname}.instr 
  	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${dir}
  	VERBATIM
  	)
  add_executable(${outname}.out ${CMAKE_CURRENT_BINARY_DIR}/${dir}/${outname}.c)
  target_link_libraries(${outname}.out PUBLIC m mcpl)

  file(MAKE_DIRECTORY ${TEST_DIR}/${outname}/)
  add_test(NAME ${outname}
	COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_test_simulation.sh ${CMAKE_CURRENT_BINARY_DIR} ${outname} ${TEST_DIR}
	)
  set_tests_properties(${outname} PROPERTIES LABELS mcstas)
  
  add_test(NAME ${outname}_output
	COMMAND  ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_test.sh ${CMAKE_CURRENT_BINARY_DIR} ${outname} ${TEST_REFERENCE_DIR} ${TEST_DIR}/${outname} 

	)
endfunction()


#------------------------------ List of instruments to be compiled
run_mcstas(D22/HEAD/D22_quick.py)
run_mcstas(ThALES/HEAD/ThALES.py)

#------------------------------ Testing
enable_testing()
