import sys
from mcstasscript.interface.functions import  load_metadata, load_monitor
import numpy as np
import pytest
pytest.main(["--rootdir=/tmp/"])

reference_file_dir = sys.argv[1]
test_file_dir = sys.argv[2]

metadata_list = load_metadata(reference_file_dir)

reference_results = {}
for metadata in metadata_list:
    monitor = load_monitor(metadata,reference_file_dir)
    reference_results[monitor.name] = monitor
metadata_list = load_metadata(test_file_dir)
test_results={}
for metadata in metadata_list:
    monitor = load_monitor(metadata,test_file_dir)
    test_results[monitor.name] = monitor

    #        reference_results[metadata.append(load_monitor(metadata, data_folder_name))

for name in list(test_results):
    if name not in reference_results:
        print("[ERROR] Monitor defined in test but not in reference")
        sys.exit(1)

for name in list(reference_results):
    if name not in test_results:
        print("[ERROR] Monitor defined in reference but not in test")
        sys.exit(1)
    assert np.array_equal(reference_results[name].Intensity, test_results[name].Intensity)
    assert np.array_equal(reference_results[name].Error, test_results[name].Error)
    assert np.array_equal(reference_results[name].Ncount, test_results[name].Ncount)

#sys.exit(1)
