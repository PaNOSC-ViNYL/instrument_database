#!/bin/bash
# this script compares the output of the test simulation with a reference committed in the repository

# this script is meant to run a test simulation for a given instrument
BINARY_DIR=$1
outname=$2
TEST_DIR=$3
TEST_REFERENCE_DIR=$4

#cd ${TEST_DIR}
source ${BINARY_DIR}/python_env/bin/activate
echo "PYTHONPATH: $PYTHONPATH"
python3 `dirname $0`/test.py  ${TEST_REFERENCE_DIR} ${TEST_DIR}

