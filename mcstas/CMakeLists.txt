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
set(MCSTAS $ENV{MCSTAS} CACHE PATH "Path of the mcstas installation [MANDATORY]")
set(PROXY "--proxy=$ENV{http_proxy}" CACHE STRING "http proxy")
############################# MCSTAS
message(STATUS "MCSTAS: ${MCSTAS}")
#include_directories(/usr/share/mcstas/2.6/libs/mcpl)
#link_directories(/usr/share/mcstas/2.6/libs/mcpl/)
include_directories(${MCSTAS}/libs/mcpl)
link_directories(${MCSTAS}/libs/mcpl/)

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
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/python_env/bin/pip ${PROXY} install -r ${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt
  COMMAND touch mcstasscript.done
  VERBATIM
  )

#set(ENV{PYTHONPATH} ${CMAKE_CURRENT_SOURCE_DIR}/ILL/instruments/D22/HEAD/)
#message(${CMAKE_CURRENT_BINARY_DIR})

set(ILL_INSTRUMENTS ${CMAKE_CURRENT_SOURCE_DIR}/ILL/instruments/)
# instrument_file should be relative to ILL_INSTRUMENTS
function(run_mcstas instrument_file)
  cmake_path(SET instrument_file_full_path ${ILL_INSTRUMENTS}/${instrument_file})

  cmake_path(REMOVE_FILENAME instrument_file OUTPUT_VARIABLE dir)

  cmake_path(GET instrument_file FILENAME outname)
  cmake_path(REMOVE_EXTENSION outname)

  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${dir})
  add_custom_command(OUTPUT ${outname}.c
  	COMMENT "creation of ${outname}.c for version: ${version}"
  	DEPENDS mcstasscript.done
  	COMMAND ${CMAKE_CURRENT_BINARY_DIR}/python_env/bin/python3 ${CMAKE_CURRENT_SOURCE_DIR}/scripts/py_to_instr.py ${instrument_file_full_path} ${outname}
  	COMMAND mcstas -t -o ${CMAKE_CURRENT_BINARY_DIR}/${outname}.c -I ${CMAKE_CURRENT_SOURCE_DIR}/ ${CMAKE_CURRENT_BINARY_DIR}/${outname}.instr 
  	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${dir}
  	VERBATIM
  	)
  add_executable(${outname}.out ${outname}.c)
  target_link_libraries(${outname}.out PUBLIC m mcpl)
endfunction()


run_mcstas(D22/HEAD/D22_quick.py)

enable_testing()
add_test(D22_quick ${CMAKE_CURRENT_BINARY_DIR}/D22_quick.out  lambda=4.51  -s 554321 -n 200000)
set_tests_properties(D22_quick PROPERTIES LABELS mcstas)
