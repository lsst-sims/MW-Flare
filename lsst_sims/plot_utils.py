import matplotlib.pyplot as plt
import numpy as np


__all__ = ["make_distribution_plot", "make_density_plot"]



def make_distribution_plot(xx_in, dd, color='k'):
    xx_dex = np.round(xx_in/dd).astype(int)
    unq, ct = np.unique(xx_dex, return_counts=True)

    sorted_dex = np.argsort(unq)

    return plt.plot(unq[sorted_dex]*dd, ct[sorted_dex], color=color)


def make_density_plot(xx_in, yy_in, dd, cmin=0, cmax=150, dc=25):

    data_x = np.round(xx_in/dd).astype(int)
    data_y = np.round(yy_in/dd).astype(int)

    xmax = data_x.max()
    ymax = data_y.max()
    xmin = data_x.min()
    ymin = data_y.min()
    factor = int(np.power(10.0, np.round(np.log10(ymax-ymin)+1.0)))
    dex_arr = (data_x-xmin)*factor + data_y-ymin

    unq, counts = np.unique(dex_arr, return_counts=True)

    x_unq = xmin + unq//factor
    y_unq = ymin + unq % factor

    grid = {}

    for xx, yy, cc in zip(x_unq, y_unq, counts):
        if xx not in grid:
            grid[xx]= {}

        if yy not in grid[xx]:
            grid[xx][yy] = cc
        else:
            grid[xx][yy] += cc

    xx_arr = []
    yy_arr = []
    ct_arr = []
    for xx in grid:
        for yy in grid[xx]:
            xx_arr.append(xx)
            yy_arr.append(yy)
            ct_arr.append(grid[xx][yy])

    ct_arr = np.array(ct_arr)
    xx_arr = np.array(xx_arr)
    yy_arr = np.array(yy_arr)

    xx_grid, yy_grid = np.mgrid[slice(xx_arr.min()*dd, xx_arr.max()*dd+dd, dd),
                                slice(yy_arr.min()*dd, yy_arr.max()*dd+dd, dd)]

    cc_grid = np.ones(xx_grid.shape)*(-99.0)

    for ix, iy, cc in zip(xx_arr, yy_arr, ct_arr):
        cc_grid[ix-xx_arr.min()][iy-yy_arr.min()] = cc

    cc_masked = np.ma.masked_values(cc_grid, -99.0)
    plt.pcolormesh(xx_grid, yy_grid, cc_masked,
                   cmap=plt.cm.gist_ncar,
                   edgecolor='')

    cb = plt.colorbar()
    plt.clim(cmin,cmax)
    cb.set_ticks(np.arange(cmin,cmax+10,dc))

    return None
