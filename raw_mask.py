from stsci_bkg import *
from mask_to_reg import main as m2r
from astropy.io import fits 
import argparse
from pathlib import Path
from scipy.ndimage import maximum_filter as mf
import warnings

def run(args):
    inp = args.input
    out = args.output
    data = fits.getdata(inp) 
    show = args.show == 'True'
    reg = args.region == 'True'
    
    #kernel = make_2dgaussian_kernel(3, size=5)  # FWHM = 3.0
    #data = convolve(data, kernel)
    #sm = SourceMask(data, nsigma=3)
    #mask = sm.multiple(filter_fwhm=[ 3, 5, 7, 10], tophat_size=[ 7, 5, 3, 1])
    from photutils.segmentation import deblend_sources
    #m2r(mask, 'append/NGC5504.reg')


    _, segm = make_source_mask(data, nsigma=3, npixels=10, dilate_size=1, filter_fwhm=3) 
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mask = deblend_sources(data, segm, npixels=10, nlevels=32, contrast=0.1, progress_bar=True)
    

    maskf = mf(mask.data, size=1)
    
    if show:
        plt.figure(figsize=(5,5), dpi=300)
        plt.title('Image mask')
        plt.imshow(maskf, origin='lower', cmap = mask.cmap )
        plt.tight_layout()
        plt.savefig('temp_mask.png')
    fits.writeto(out, 1*maskf, overwrite=True)
    if reg:
        m2r(mask, 'temp_mask.reg')

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    
    parser.add_argument("-input", type=str,  help='Input file')
    parser.add_argument("--output", default='temp_mask.fits',type=str, help='Output file name')
    parser.add_argument("--show", default='True',type=str, help='Show result key. True or False.')
    parser.add_argument("--region", default='False',type=str, help='Create region file result key. True or False.')
    
    args = parser.parse_args()
    run(args)

