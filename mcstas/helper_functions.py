import numpy as np
from scipy.stats import norm
import copy


class PlotHelper:
    """Class helper for comparison plots of two datasets for one 2D detector"""

    class stat:
        """structure containing the stat info"""

        mean = None
        std = None

    @property
    def pixels(self):
        return self.__pixel_size

    def set_pixel_size(self, x, y):
        print("Setting: ", x, y)
        self.__pixel_size[0] = float(x)
        self.__pixel_size[1] = float(y)
        print(self.__pixel_size)

    def set_limits(self, lims):
        if len(lims) == 2:
            self.__xmin = [lims[0]]
            self.__xmax = [lims[1]]
        elif len(lims) == 4:
            self.__xmin = [lims[0], lims[2]]
            self.__xmax = [lims[1], lims[3]]
        else:
            raise RuntimeError("limits should be list of 2 or 4 elements")

        x, wx = self.marginal(axis=1)
        y, wy = self.marginal(axis=0)
        self.__stats[1].mean, self.__stats[1].std = self.calc_stats(x, self.xbins(1))
        self.__stats[0].mean, self.__stats[0].std = self.calc_stats(y, self.xbins(0))

    def __init__(self, data, weights=None):
        """
        data: 2D numpy array
        """
        if data.ndim != 2:
            raise RuntimeError(f"wrong number of dimensions <{data.ndim}>")

        self.__scatter_plot_settings = {"s": 10}
        self.__errorbar_plot_settings = {"fmt": "o"}
        self.__data = None
        self.__weights = None
        self.__stats = None
        self.__pixel_size = [1, 1]
        self.__xmin = None

        self.__data = data
        if weights is None:
            self.__weights = np.ones(data.shape)

        # make 1D statistics
        x, wx = self.marginal(axis=1)
        y, wy = self.marginal(axis=0)

        self.__stats = [self.stat(), self.stat()]  # [x[mean,std], y[mean, std]]
        self.__stats[1].mean, self.__stats[1].std = self.calc_stats(x, self.xbins(1))
        self.__stats[0].mean, self.__stats[0].std = self.calc_stats(y, self.xbins(0))
        self.__norm = np.sum(self.__data)

    @property
    def norm(self):
        return self.__norm

    @property
    def stats(self):
        return self.__stats

    def calc_stats(self, hist, bins=None):
        """Calculate mean and standard deviation for a numpy array representing the entries for an evenly binned histogram"""
        if bins is None:
            bins = np.arange(0, hist.shape[0])
        bin_values = hist

        sumw = np.sum(bin_values)
        if sumw != 0:
            sumwx = np.sum(bins * bin_values)
            sumwx2 = np.sum(np.square(bins) * bin_values)
            m = sumwx / sumw
            std = sumwx2 / sumw - m * m
        else:
            m = 0
            std = 0
        return (m, np.sqrt(std))

    def marginal(self, axis):
        """currently implemented only for the 2D case, not generalized
        axis: axis to marginalize, i.e. axis desappearing
        """
        x = np.sum(self.__data, axis=axis)
        w = np.sqrt(np.sum(np.square(self.__weights), axis=axis))
        return x, w

    def limits(self, nsigma, axis):
        return (
            self.__stats[axis].mean - self.__stats[axis].std * nsigma,
            self.__stats[axis].mean + self.__stats[axis].std * nsigma,
        )

    def xbins(self, axis):
        x = None
        if self.__xmin is not None:
            x = np.linspace(
                self.__xmin[axis], self.__xmax[axis], self.__data.shape[axis]
            )
        else:
            x = np.arange(
                0,
                self.__data.shape[axis] * self.__pixel_size[axis],
                self.__pixel_size[axis],
            )
        return x

    def plot_marginal(self, ax, axis, isData, label=None, plotGauss=False, norm=False):
        if axis not in [0, 1]:
            raise RuntimeError("plot_marginal: axis should be 0 or 1")

        y, w = self.marginal(axis)
        if y.sum() == 0:
            ax.plot(np.NaN, np.NaN, "-", color="none", label=label)
            return
        x = None
        if self.__xmin is not None:
            x = np.linspace(self.__xmin[axis], self.__xmax[axis], y.shape[0])
        else:
            x = np.arange(
                0, y.shape[0] * self.__pixel_size[axis], self.__pixel_size[axis]
            )
        yerr = np.sqrt(y)
        if norm:
            norm = y.sum()
        else:
            norm = 1

        l = None
        if isData is True:
            dx = x[1] - x[0]
            xedges = np.append(x - dx / 2, x[-1] + dx / 2)
            l = ax.stairs(y / norm, xedges, fill=True)
        else:
            l = ax.errorbar(
                x,
                y / norm,
                yerr=yerr / norm,
                label=label,
                **self.__errorbar_plot_settings,
            )
        # ax.set_xlim(np.array(self.limits(3, axis)) * self.__pixel_size[axis])

        if plotGauss:
            x = np.linspace(0, y.shape[0], y.shape[0] * 10)
            g = norm.pdf(x, self.stats[axis].mean, self.stats[axis].std) * y.sum()
            print(x, y)
            print(y.max(), self.stats[axis].mean, self.stats[axis].std)
            ax.plot(x, g)

        return l

    def legend_options_marginal(self):
        return {"bbox_to_anchor": (0, 1.15), "loc": "upper center"}

    def plot_2D(self, ax, axis, label=None, norm=None):
        settings = {
            "cmap": "hsv",
            "origin": "lower",
            "extent": [
                self.__xmin[0],
                self.__xmax[0],
                self.__xmin[1],
                self.__xmax[1],
                # self.__pixel_size[0] * self.__data.shape[0],
                # 0, self.__pixel_size[1] * self.__data.shape[1],
            ],
        }
        print(
            self.__pixel_size,
            self.__data.shape,
        )
        if norm is not None:
            l = ax.imshow(self.__data.T / norm, vmin=0, vmax=1, **settings)
        else:
            l = ax.imshow(self.__data.T, **settings)
        return l


class ComparePlots:
    def __init__(self, ph1, ph2, legend1, legend2):
        """ph1 and ph2: dict("name": PlotHelper)"""
        self._ph1 = ph1
        self._ph2 = ph2
        self._legend1 = legend1
        self._legend2 = legend2

    def plot_2Ds(self, fig2D, plot_positions, axs, colorbar_position: str):
        norm_ph1 = 0
        norm_ph2 = 0

        for d in plot_positions:
            norm_ph1 = norm_ph1 + self._ph1[d].norm
            if self._ph2 is not None and len(self._ph2) > 0:
                norm_ph2 = norm_ph2 + self._ph2[d].norm

        imgs1 = []
        imgs2 = []

        for d in plot_positions:
            ax_sim = axs[plot_positions[d] + "_sim"]
            imgs1.append(self._ph1[d].plot_2D(ax_sim, self._legend1))
            if self._ph2 is not None and len(self._ph2) > 0:
                ax_data = axs[plot_positions[d] + "_data"]
                imgs2.append(self._ph2[d].plot_2D(ax_data, self._legend1))

        from matplotlib import colors

        vmin = 0
        vmax1 = max(image.get_array().max() for image in imgs1)
        norm1 = colors.Normalize(vmin=vmin, vmax=vmax1)
        vmax2 = 0
        if len(imgs2) > 0:
            vmax2 = max(image.get_array().max() for image in imgs2)
            norm2 = colors.Normalize(vmin=vmin, vmax=vmax2)

        for im in imgs1:
            im.set_norm(norm1)
        for im in imgs2:
            im.set_norm(norm2)

        fig2D.colorbar(
            imgs1[0], ax=axs[colorbar_position + "_sim"], orientation="vertical"
        )
        if len(imgs2) > 0:
            fig2D.colorbar(
                imgs2[0], ax=axs[colorbar_position + "_data"], orientation="vertical"
            )
        # print(a.stats)
        # print(a.limits_x(4))
        # axs["x"].set_xlim(a.mean_y - a.std_y * 4, a.mean_y + a.std_y * 4)

    def plot_marginals(self, plot_positions, axs):
        """plot_positions = dict("detactor_name": "position")
        e.g. plot_positions = {"detector_left": "l", "detector_central": "c", "detector_right": "r"}
        axs = fig1D.subplot_mosaic(
        \"\"\"
        lcr
        LCR
        \"\"\",
        # set the height ratios between the rows
        # height_ratios=[1, 1],
        # set the width ratios between the columns
        # width_ratios=[1, 2, 1],
        )
        """

        for d in plot_positions:
            ax_x = axs[plot_positions[d] + "x"]
            ax_y = axs[plot_positions[d] + "y"]
            self._ph1[d].plot_marginal(ax_x, 0, False, label=self._legend1, norm=True)
            self._ph1[d].plot_marginal(ax_y, 1, False, label=self._legend1, norm=True)
            if self._ph2 is not None and len(self._ph2) > 0:
                self._ph2[d].plot_marginal(
                    ax_x, 0, True, label=self._legend2, norm=True
                )
                self._ph2[d].plot_marginal(
                    ax_y, 1, True, label=self._legend2, norm=True
                )

            box = ax_x.get_position()
            # ax_x.set_position([box.x0, box.y0, box.width, box.height * 0.9])
            ax_x.legend(
                [self._legend1, self._legend2], **self._ph1[d].legend_options_marginal()
            )
            ax_y.legend(
                [self._legend1, self._legend2], **self._ph1[d].legend_options_marginal()
            )

            if self._ph2 is not None and len(self._ph2) > 0:
                ax_x.annotate(
                    "Mean x ({}): {:.2f}\nStd. dev.: {:.2f}\nMean x ({}): {:.2f}\nStd. dev.: {:.2f}".format(
                        self._legend1,
                        self._ph1[d].stats[0].mean,
                        self._ph1[d].stats[0].std,
                        self._legend2,
                        self._ph2[d].stats[0].mean,
                        self._ph2[d].stats[0].std,
                    ),
                    (0.5, 1.1),
                    xycoords="axes fraction",
                )
                ax_y.annotate(
                    "Mean y ({}): {:.2f}\nStd. dev.: {:.2f}\nMean y ({}): {:.2f}\nStd. dev.: {:.2f}".format(
                        self._legend1,
                        self._ph1[d].stats[1].mean,
                        self._ph1[d].stats[1].std,
                        self._legend2,
                        self._ph2[d].stats[1].mean,
                        self._ph2[d].stats[1].std,
                    ),
                    (0.5, 1.1),
                    xycoords="axes fraction",
                )
            else:
                ax_x.annotate(
                    "Mean x ({}): {:.2f}\nStd. dev.: {:.2f}".format(
                        self._legend1,
                        self._ph1[d].stats[0].mean,
                        self._ph1[d].stats[0].std,
                    ),
                    (0.5, 1.1),
                    xycoords="axes fraction",
                )
                ax_y.annotate(
                    "Mean y ({}): {:.2f}\nStd. dev.: {:.2f}".format(
                        self._legend1,
                        self._ph1[d].stats[1].mean,
                        self._ph1[d].stats[1].std,
                    ),
                    (0.5, 1.1),
                    xycoords="axes fraction",
                )
