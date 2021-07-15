cmake_minimum_required(VERSION 3.14.0)
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
############################# MCSTAS
message(STATUS "MCSTAS: ${MCSTAS}")
#include_directories(/usr/share/mcstas/2.6/libs/mcpl)
#link_directories(/usr/share/mcstas/2.6/libs/mcpl/)
include_directories(${MCSTAS}/libs/mcpl)
link_directories(${MCSTAS}/libs/mcpl/)

#set(PROXY "--proxy=http://193.49.43.123:8888")
set(PROXY "--proxy=$ENV{http_proxy}" CACHE STRING "http proxy")


add_custom_target(python_env.done
  COMMENT "Creating python virtual environment"
  COMMAND ${Python_EXECUTABLE} -m venv ${CMAKE_CURRENT_BINARY_DIR}/python_env
  COMMAND touch python_env.done
  VERBATIM
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

add_custom_target(mcstasscript.done
  COMMENT "Install fresh release of McStasScript"
  DEPENDS python_env.done
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/python_env/bin/pip ${PROXY} install -r ${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt
  COMMAND touch mcstasscript.done
  VERBATIM
  )

set(ENV{PYTHONPATH} ${CMAKE_CURRENT_SOURCE_DIR}/ILL/instruments/D22/HEAD/)
message(${CMAKE_CURRENT_BINARY_DIR})

function(run_mcstas outname )
  add_custom_command(OUTPUT ${outname}.c
	COMMENT "creation of ${outname}.c"
	DEPENDS mcstasscript.done
	#DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${outname}.instr #${CMAKE_CURRENT_SOURCE_DIR}/MCPL_input.comp
	COMMAND PYTHONPATH=$ENV{PYTHONPATH} ${CMAKE_CURRENT_BINARY_DIR}/python_env/bin/python3 ${CMAKE_CURRENT_SOURCE_DIR}/scripts/py_to_instr.py ${CMAKE_CURRENT_SOURCE_DIR}/ILL/instruments/D22/HEAD/D22_quick.py D22_quick
	COMMAND mcstas -t -o ${CMAKE_CURRENT_BINARY_DIR}/${outname}.c -I ${CMAKE_CURRENT_SOURCE_DIR}/ ${CMAKE_CURRENT_BINARY_DIR}/${outname}.instr 
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/
	VERBATIM
	)
  add_executable(${outname}.out ${outname}.c)
  target_link_libraries(${outname}.out PUBLIC m mcpl)
endfunction()


run_mcstas(D22_quick)
