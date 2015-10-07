# -*- coding: utf8 -*-
"""Utilities imported and slightly adapted from exterbal sources :

    - http://matplotlib.sourceforge.net/examples/pylab_examples/demo_agg_filter.html
"""

import matplotlib.pyplot as plt

import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
from matplotlib.artist import Artist
import matplotlib.transforms as mtransforms

def smooth1d(x, window_len):
    # copied from http://www.scipy.org/Cookbook/SignalSmooth (9,) 9 (24,) (8,)

    s=np.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    w = np.hanning(window_len)
    y=np.convolve(w/w.sum(),s,mode='same')
#    # FIXME: smooth1d with size equal to window_len (9)
#    if x.shape!=y[window_len-1:-window_len+1].shape:
#        print 'x.shape, window_len,  y.shape, y[window_len-1:-window_len+1].shape', x.shape, window_len,  y.shape, y[window_len-1:-window_len+1].shape
#        print 'np.convolve(w/w.sum(),s,mode="valid").shape', np.convolve(w/w.sum(),s,mode='valid').shape
#    else:
#        print 'x.shape, window_len,  y.shape, y[window_len-1:-window_len+1].shape', x.shape, window_len,  y.shape, y[window_len-1:-window_len+1].shape
##    return y[window_len-1:-window_len+1]
    return y[window_len-1:window_len-1+x.shape[0]]

def smooth2d(A, sigma=3):

    window_len = max(int(sigma), 3)*2+1
#    if window_len>=x.shape[0]
    A1 = np.array([smooth1d(x, window_len) for x in np.asarray(A)])
    A2 = np.transpose(A1)
    A3 = np.array([smooth1d(x, window_len) for x in A2])
    A4 = np.transpose(A3)

    return A4




class BaseFilter(object):
    def prepare_image(self, src_image, dpi, pad):
        ny, nx, depth = src_image.shape
        #tgt_image = np.zeros([pad*2+ny, pad*2+nx, depth], dtype="d")
        padded_src = np.zeros([pad*2+ny, pad*2+nx, depth], dtype="d")
        padded_src[pad:-pad, pad:-pad,:] = src_image[:,:,:]

        return padded_src#, tgt_image

    def get_pad(self, dpi):
        return 0

    def __call__(self, im, dpi):
        pad = self.get_pad(dpi)
        padded_src = self.prepare_image(im, dpi, pad)
        tgt_image = self.process_image(padded_src, dpi)
        return tgt_image, -pad, -pad


class OffsetFilter(BaseFilter):
    def __init__(self, offsets=None):
        if offsets is None:
            self.offsets = (0, 0)
        else:
            self.offsets = offsets

    def get_pad(self, dpi):
        return int(max(*self.offsets)/72.*dpi)

    def process_image(self, padded_src, dpi):
        ox, oy = self.offsets
        a1 = np.roll(padded_src, int(ox/72.*dpi), axis=1)
        a2 = np.roll(a1, -int(oy/72.*dpi), axis=0)
        return a2

class GaussianFilter(BaseFilter):
    "simple gauss filter"
    def __init__(self, sigma, alpha=0.5, color=None):
        self.sigma = sigma
        self.alpha = alpha
        if color is None:
            self.color=(0, 0, 0)
        else:
            self.color=color

    def get_pad(self, dpi):
        return int(self.sigma*3/72.*dpi)


    def process_image(self, padded_src, dpi):
        #offsetx, offsety = int(self.offsets[0]), int(self.offsets[1])
        tgt_image = np.zeros_like(padded_src)
        aa = smooth2d(padded_src[:,:,-1]*self.alpha,
                      self.sigma/72.*dpi)
        tgt_image[:,:,-1] = aa
        tgt_image[:,:,:-1] = self.color
        return tgt_image

class DropShadowFilter(BaseFilter):
    def __init__(self, sigma, alpha=0.3, color=None, offsets=None):
        self.gauss_filter = GaussianFilter(sigma, alpha, color)
        self.offset_filter = OffsetFilter(offsets)

    def get_pad(self, dpi):
        return max(self.gauss_filter.get_pad(dpi),
                   self.offset_filter.get_pad(dpi))

    def process_image(self, padded_src, dpi):
        t1 = self.gauss_filter.process_image(padded_src, dpi)
        t2 = self.offset_filter.process_image(t1, dpi)
        return t2

class GrowFilter(BaseFilter):
    "enlarge the area"
    def __init__(self, pixels, color=None, alpha=1.):
        self.pixels = pixels
        if color is None:
            self.color=(1, 1, 1)
        else:
            self.color=color
        self.alpha = alpha

    def __call__(self, im, dpi):
        pad = self.pixels
        ny, nx, depth = im.shape
        new_im = np.empty([pad*2+ny, pad*2+nx, depth], dtype="d")
        alpha = new_im[:,:,3]
        alpha.fill(0)
        alpha[pad:-pad, pad:-pad] = im[:,:,-1]
        alpha2 = np.clip(smooth2d(alpha, self.pixels/72.*dpi) * 5, 0, 1)*self.alpha
        new_im[:,:,-1] = alpha2
        new_im[:,:,:-1] = self.color
        offsetx, offsety = -pad, -pad

        return new_im, offsetx, offsety


class LightFilter(BaseFilter):
    "simple gauss filter"
    def __init__(self, sigma, fraction=0.5, **kwargs):
        self.gauss_filter = GaussianFilter(sigma, alpha=1)
        from matplotlib.colors import LightSource
        self.light_source = LightSource(**kwargs)
        self.fraction = fraction
        #hsv_min_val=0.5,hsv_max_val=0.9,
        #                                hsv_min_sat=0.1,hsv_max_sat=0.1)
    def get_pad(self, dpi):
        return self.gauss_filter.get_pad(dpi)

    def process_image(self, padded_src, dpi):
        t1 = self.gauss_filter.process_image(padded_src, dpi)
        elevation = t1[:,:,3]
        rgb = padded_src[:,:,:3]

        rgb2 = self.light_source.shade_rgb(rgb, elevation,
                                           fraction=self.fraction)
        print 'process'
        tgt = np.empty_like(padded_src)
        tgt[:,:,:3] = rgb2
        tgt[:,:,3] = padded_src[:,:,3]

        return tgt



class FilteredArtistList(Artist):
    """
    A simple container to draw filtered artist.
    """
    def __init__(self, artist_list, filter):
        self._artist_list = artist_list
        self._filter = filter
        Artist.__init__(self)

    def draw(self, renderer):
        renderer.start_rasterizing()
        if hasattr(renderer, 'start_filter'): renderer.start_filter()
        for a in self._artist_list:
            if hasattr(a, 'draw'):
                a.draw(renderer)
        renderer.stop_filter(self._filter)
        renderer.stop_rasterizing()
