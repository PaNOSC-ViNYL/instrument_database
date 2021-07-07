#!/bin/bash
# this script is meant to run a test simulation for a given instrument
# this script is needed in order to clean up the simulation output directory each time the ctest is launched

# directory where to find the mcstas executable
BINARY_DIR=$1
# name of the mcstas executable
outname=$2
# base directory where the test simulation output should go
TEST_DIR=$3

if [ ! -e $TEST_DIR/${outname}/ ]
then
	mkdir ${TEST_DIR}/ -p
else
	rm ${TEST_DIR}/${outname}/ -Rf
fi

cd ${TEST_DIR}
$BINARY_DIR/${outname}.out -s 554321 -n 200000 -d ${outname}/ stage=-1
