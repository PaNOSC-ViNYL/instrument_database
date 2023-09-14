"""
Reads a D11 NeXus file ...
"""


def center_of_mass(data, nx_min, nx_max, ny_min, ny_max):
    cx = 0.0
    cy = 0.0
    nx = 0.0
    ny = 0.0
    for i in range(nx_min, nx_max + 1):
        cx += np.sum(data[i, :]) * i
        nx += np.sum(data[i, :])
    for i in range(ny_min, ny_max + 1):
        cy += np.sum(data[:, i]) * i
        ny += np.sum(data[:, i])
    return cx / nx + 1, cy / ny + 1


import sys
import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
from mcstasscript.interface.functions import load_metadata, load_monitor


file_data = sys.argv[1]
simulation_dir = sys.argv[2]

print("Reading file: " + file_data)


# Read NeXus
f = h5py.File(file_data, "r")
central_detector_data = f["entry0"]["D11"]["Detector 1"]["data"][:, :, 0]
left_detector_data = f["entry0"]["D11"]["Detector 2"]["data"][:, :, 0]
right_detector_data = f["entry0"]["D11"]["Detector 3"]["data"][:, :, 0]

detectors_data = {
    "left": left_detector_data,
    "central": central_detector_data,
    "right": right_detector_data,
}

metadata_list = load_metadata(simulation_dir)
print(left_detector_data.shape, central_detector_data.shape, right_detector_data.shape)
# print(metadata_list)

detectors_simulation = {}
for metadata in metadata_list:
    if metadata.component_name in [
        "detector_left",
        "detector_right",
        "detector_central",
    ]:
        monitor = load_monitor(metadata, simulation_dir)
        if monitor.name == "detector_left":
            detectors_simulation["left"] = monitor
        elif monitor.name == "detector_right":
            detectors_simulation["right"] = monitor
        elif monitor.name == "detector_central":
            detectors_simulation["central"] = monitor
        # print("### MONITOR:")
        # help(monitor)

print(
    detectors_simulation["left"].Intensity.transpose().shape,
    detectors_simulation["central"].Intensity.transpose().shape,
    detectors_simulation["right"].Intensity.transpose().shape,
)
ratio = {}
ratio["left"] = (
    detectors_simulation["left"].Intensity.transpose() + detectors_data["left"]
)
# print(ratio["left"])
# help(detectors_simulation["right"].Intensity)
# print(metadata_list["detector_right"])
# plt.show()


fig, axs = plt.subplots(2, 3)
fig.suptitle("Vertically stacked subplots")
axs[0][0].imshow(
    detectors_data["left"].transpose(), aspect="auto", cmap="seismic", origin="lower"
)
axs[0][1].imshow(
    detectors_data["central"].transpose(), aspect="auto", cmap="seismic", origin="lower"
)
axs[0][2].imshow(
    detectors_data["right"].transpose(), aspect="auto", cmap="seismic", origin="lower"
)

axs[1][0].imshow(
    detectors_simulation["left"].Intensity,
    aspect="auto",
    cmap="seismic",
    origin="lower",
)
axs[1][1].imshow(
    detectors_simulation["central"].Intensity,
    aspect="auto",
    cmap="seismic",
    origin="lower",
)
axs[1][2].imshow(
    detectors_simulation["right"].Intensity,
    aspect="auto",
    cmap="seismic",
    origin="lower",
)
# cbar = fig.colorbar(detectors_simulation["central"])
# cbar.set_label("ZLabel", loc="top")
# cax = plt.axes([0.85, 0.1, 0.075, 0.8])
# plt.colorbar(cax=cax)
plt.show()
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
