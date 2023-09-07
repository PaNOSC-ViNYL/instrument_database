"""
Reads a D11 NeXus file ...
"""

import sys
import os
import h5py
import numpy as np
import matplotlib.pyplot as plt


file = sys.argv[1]
print("Reading file" + file)


# Read NeXus
f = h5py.File(file, "r")
central_detector = f["entry0"]["D11"]["Detector 1"]["data"]
print(central_detector.shape)
plt.imshow(central_detector, aspect="auto", cmap="seismic", origin="lower")
f.close()

sys.exit(0)
print("NeXus data shapes:")
print(data1.shape, data2.shape, data3.shape, data4.shape, data5.shape)
print()

if plot2D:

    plt.imshow(
        np.transpose(data2[:, :, 0]), aspect="auto", cmap="seismic", origin="lower"
    )
    plt.title("Detector 2 (right)")
    plt.show()

    plt.imshow(
        np.transpose(data3[:, :, 0]), aspect="auto", cmap="seismic", origin="lower"
    )
    plt.title("Detector 3 (left)")
    plt.show()

    plt.imshow(
        np.transpose(data4[:, :, 0]), aspect="auto", cmap="seismic", origin="lower"
    )
    plt.title("Detector 4 (bottom)")
    plt.show()

    plt.imshow(
        np.transpose(data5[:, :, 0]), aspect="auto", cmap="seismic", origin="lower"
    )
    plt.title("Detector 5 (top)")
    plt.show()

print("Center of mass full panels (counting from 1):")
cx1, cy1 = center_of_mass(data1[:, :, 0], 0, 255, 0, 127)
cx2, cy2 = center_of_mass(data2[:, :, 0], 0, 31, 0, 255)
cx3, cy3 = center_of_mass(data3[:, :, 0], 0, 31, 0, 255)
cx4, cy4 = center_of_mass(data4[:, :, 0], 0, 255, 0, 31)
cx5, cy5 = center_of_mass(data5[:, :, 0], 0, 255, 0, 31)
print("Detector 1 = ", cx1, cy1)
print("Detector 2 = ", cx2, cy2)
print("Detector 3 = ", cx3, cy3)
print("Detector 4 = ", cx4, cy4)
print("Detector 5 = ", cx5, cy5)

print()
print("Center of mass using ROI (counting from 1):")
cx1, cy1 = center_of_mass(data1[:, :, 0], 118, 134, 59, 66)
cx2, cy2 = center_of_mass(data2[:, :, 0], 11, 16, 127, 138)
cx3, cy3 = center_of_mass(data3[:, :, 0], 14, 18, 124, 136)
cx4, cy4 = center_of_mass(data4[:, :, 0], 120, 131, 13, 18)
cx5, cy5 = center_of_mass(data5[:, :, 0], 118, 130, 12, 16)
print("Detector 1 = ", cx1, cy1)
print("Detector 2 = ", cx2, cy2)
print("Detector 3 = ", cx3, cy3)
print("Detector 4 = ", cx4, cy4)
print("Detector 5 = ", cx5, cy5)
