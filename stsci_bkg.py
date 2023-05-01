# Данный код заимствован с сайта: https://spacetelescope.github.io/jdat_notebooks/notebooks/background_estimation_imaging/Imaging_Sky_Background_Estimation.html
# И собран мною в одну функцию для учета фона неба.

# Loading nessecary packages
import numpy as np
from scipy import stats, ndimage, interpolate
from astropy import stats as astrostats
from photutils.background import (
    Background2D,  # For estimating the background
    MedianBackground, BiweightLocationBackground, SExtractorBackground,
    BkgIDWInterpolator, BkgZoomInterpolator)  # For interpolating background
from photutils import datasets  # For making simulated data
#from photutils.segmentation import make_source_mask
#from photutils.utils import ShepardIDWInterpolator as idw
from astropy.table import Table
from astropy.convolution import (
    convolve, Box2DKernel, Tophat2DKernel,
    Ring2DKernel, Gaussian2DKernel)
from scipy.ndimage import median_filter
from astropy.modeling import models, fitting
from astropy.nddata.blocks import block_reduce
#from IPython.display import HTML
#from jdaviz import Imviz

import matplotlib.pyplot as plt
import matplotlib as mpl

# Make sources mask

from astropy.convolution import convolve
from photutils.segmentation import make_2dgaussian_kernel
from photutils.segmentation import detect_sources
from photutils.segmentation import detect_threshold
from photutils.utils import circular_footprint
def make_source_mask(data, nsigma, npixels, dilate_size, filter_fwhm): 
    #data = scene - bkg.background
    threshold = detect_threshold(data, nsigma) #nsigma * bkg.background_rms
    #kernel = make_2dgaussian_kernel(filter_fwhm, size=5)  # FWHM = 3.0
    #convolved_data = convolve(data, kernel)
    #segment_map = detect_sources(convolved_data, threshold, npixels=npixels)
    
    #footprint = circular_footprint(dilate_size)
    #mask_2sigma = segment_map.make_source_mask(footprint=footprint)
    

    kernel = Tophat2DKernel(filter_fwhm)
    convolved_data = convolve(data, kernel)
    segment_map = detect_sources(convolved_data, threshold, npixels=npixels)

    #mask_2sigma += segment_map.make_source_mask(footprint=footprint)
    
    return  ' ' , segment_map

# grow mask
def dilate_mask(mask, tophat_size):
    ''' Take a mask and make the masked regions bigger.'''
    area = np.pi*tophat_size**2.
    kernel = Tophat2DKernel(tophat_size)
    dilated_mask = convolve(mask, kernel) >= 1./area
    return dilated_mask

# Iterating mask
def my_background(img, box_size, mask, interp=None, filter_size=1,
                  exclude_percentile=90, c_mask=None):
    #c_mask = (scene==0)
    ''' Run photutils background with SigmaClip and MedianBackground'''
    if interp is None:
        interp = BkgZoomInterpolator()
    return Background2D(img, box_size,
                        sigma_clip=astrostats.SigmaClip(sigma=3.),
                        filter_size=filter_size,
                        bkg_estimator=MedianBackground(),
                        exclude_percentile=exclude_percentile,
                        mask=mask,
                        interpolator=interp,
                        coverage_mask = c_mask,
                       )

class SourceMask:
    def __init__(self, img, nsigma=3., npixels=3):
        ''' Helper for making & dilating a source mask.
             See Photutils docs for make_source_mask.'''
        self.img = img
        self.nsigma = nsigma
        self.npixels = npixels

    def single(self, filter_fwhm=3., tophat_size=5., mask=None):
        '''Mask on a single scale'''
        
        if mask is None:
            plt.show()
            image = self.img
        else:
            image = self.img*(1-mask)
        mask = make_source_mask(image, nsigma=self.nsigma,
                                npixels=self.npixels,
                                dilate_size=1, filter_fwhm=filter_fwhm)
        print(mask.shape)
        return dilate_mask(mask, tophat_size)

    def multiple(self, filter_fwhm=[3.], tophat_size=[3.], mask=None):
        '''Mask repeatedly on different scales'''
        if mask is None:
            self.mask = np.zeros(self.img.shape, dtype=bool)
        for fwhm, tophat in zip(filter_fwhm, tophat_size):
            smask = self.single(filter_fwhm=fwhm, tophat_size=tophat)
            self.mask = self.mask | smask  # Or the masks at each iteration
        return self.mask


def main(scene, add_mask):
    
    c_mask = (scene == 0)

    mask_3sigma = make_source_mask(scene, nsigma=1.5, npixels=1, dilate_size=10, filter_fwhm=3)
    mask = mask_3sigma
    mask = dilate_mask(mask, 11)
    # Создание первичной маски
    #ring = Ring2DKernel(40, 5)
    #filtered = median_filter(scene, footprint=ring.array)
    #difference = scene-filtered
    #smoothed = convolve(difference, Gaussian2DKernel(3))
    #mask = smoothed > 1.*smoothed.std()


    interpolator = BkgIDWInterpolator(n_neighbors=20, power=1, reg=30)
    bkg4 = my_background(scene, box_size=20, filter_size=3, mask=mask | add_mask,
                       interp=interpolator, exclude_percentile=90, c_mask=c_mask)
   
    sm = SourceMask(scene - bkg4.background, nsigma=1.5)
    mask = sm.multiple(filter_fwhm=[ 3, 5, 7 ,10],
                   tophat_size=[ 7, 5, 3, 1])

    interpolator = BkgIDWInterpolator(n_neighbors=20, power=0, reg=30)

    bkg5 = my_background(scene, box_size=6, filter_size=3, mask=mask | add_mask,
                       interp=interpolator, exclude_percentile=90, c_mask=c_mask).background
    
    bkgd = bkg5

    c_mask = (scene == 0)

    mean_masked = bkgd[mask].mean()
    std_masked = bkgd[mask].std()
    stderr_masked = mean_masked/(np.sqrt(len(bkgd[mask]))*std_masked)

    mean_unmasked = bkgd[~(mask | c_mask)].mean()
    std_unmasked = bkgd[~(mask | c_mask)].std()
    stderr_unmasked = mean_unmasked/(np.sqrt(len(bkgd[~(mask| c_mask)]))*std_unmasked)

    diff = mean_masked - mean_unmasked
    significance = diff/np.sqrt(stderr_masked**2 + stderr_unmasked**2)

    print(f"Mean under masked pixels   = {mean_masked:.4f} +- {stderr_masked:.4f}")
    print(f"Mean under unmasked pixels = "
          f"{mean_unmasked:.4f} +- {stderr_unmasked:.4f}")
    print(f"Difference = {diff:.4f} at {significance:.2f} sigma significance")
    return bkg5, mask+c_mask
