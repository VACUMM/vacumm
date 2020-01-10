"""
Utilities for running tests
"""

import numpy as np
import matplotlib.pyplot as plt


def get_slices_(ndim, axis, index):
    slices = [slice(None)] * ndim
    slices[axis] = index
    return tuple(slices)


def meshgrid(x, y):
    x = np.asarray(x)
    y = np.asarray(y)
    if x.ndim == 1 and y.ndim == 1:
        x, y = np.meshgrid(x, y)
    elif x.ndim == 1:
        x = np.resize(x, y.shape)
    elif y.ndim == 1:
        y = np.repeat(y[:, None], x.shape[1], axis=1)
    return x, y


def meshcells(x, y=None, axis=-1):

    x = np.asarray(x)

    # 1D
    if y is None:
        sh = list(x.shape)
        sh[axis] += 1
        xb = np.empty(sh)
        xb[get_slices_(x.ndim, axis, slice(1, -1))] = 0.5 * (
            x[get_slices_(x.ndim, axis, slice(None, -1))] +
            x[get_slices_(x.ndim, axis, slice(1, None))])
        xb[get_slices_(x.ndim, axis, 0)] = (
            2 * x[get_slices_(x.ndim, axis, 0)] -
            x[get_slices_(x.ndim, axis, 1)])
        xb[get_slices_(x.ndim, axis, -1)] = (
            2 * x[get_slices_(x.ndim, axis, -1)] -
            x[get_slices_(x.ndim, axis, -2)])
        return xb

    # 2D
    x, y = meshgrid(x, y)
    ny, nx = x.shape
    xxb = np.empty((ny+1, nx+1))
    yyb = np.empty((ny+1, nx+1))
    for zzb, z in [(xxb, x), (yyb, y)]:
        zzb[1:-1, 1:-1] = 0.25 * (z[1:, 1:] + z[:-1, 1:] +
                                  z[1:, :-1] + z[:-1, :-1])
        zzb[1:-1, 0] = 2 * zzb[1:-1, 1] - zzb[1:-1, 2]
        zzb[1:-1, -1] = 2 * zzb[1:-1, -2] - zzb[1:-1, -3]
        zzb[0, 1:-1] = 2 * zzb[1, 1:-1] - zzb[2, 1:-1]
        zzb[-1, 1:-1] = 2 * zzb[-2, 1:-1] - zzb[-3, 1:-1]
        zzb[0, 0] = 2 * zzb[1, 1] - zzb[2, 2]
        zzb[-1, -1] = 2 * zzb[-2, -2] - zzb[-3, -3]
        zzb[0, -1] = 2 * zzb[1, -2] - zzb[2, -3]
        zzb[-1, 0] = 2 * zzb[-2, 1] - zzb[-3, 2]
    return xxb, yyb


def rotate_grid(x, y, angle):
    x, y = meshgrid(x, y)
    xc = x.min()
    yc = y.min()
    angle *= np.pi / 180
    ca, sa = np.cos(angle), np.sin(angle)
    xr, yr = np.array([[ca, -sa], [sa, ca]]).dot([x.ravel()-xc, y.ravel()-yc])
    return xr.reshape(x.shape) + xc, yr.reshape(y.shape) + yc


def plot_grid(x, y, color='k', linewidth=.5, **kwargs):
    x, y = meshgrid(x, y)
    plt.plot(x, y, color=color, linewidth=linewidth, **kwargs)
    plt.plot(x.T, y.T, color=color, linewidth=linewidth, **kwargs)
