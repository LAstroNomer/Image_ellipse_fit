import argparse
import fit_ell
import mask2contures
import subprocess
from pathlib import Path
from astropy.io import fits
import numpy as np
from time import time
import os
#==========================================
q = 0.15
#=========================================
def save_pars_to_data(pars, Ps, file):
    with open(file, 'w') as f:
        print('#Num', 'Xc(pix)', 'Yc(pix)', 'A(pix)', 'e', 'PA(deg)', 'eps/A', 'fit (1=good)', file=f)
        for i, par, p in zip(range(len(Ps)), pars, Ps):
            if p > q:
                key = 0
            else:
                key = 1

            xc = par[0]
            yc = par[1]
            A  = np.max([par[2], par[3]])
            B  = np.min([par[3], par[2]])
            e  = np.sqrt(1 - B**2/A**2)
            PA = np.degrees(par[4]) 
            
            print('%5i %8.2f %8.2f %8.2f %8.3f %8.2f %8.3f %5i' %(i, xc, yc, A, e, PA, p, key), file=f)
    return
    

if __name__ == '__main__':
     parser = argparse.ArgumentParser()
     parser.add_argument("-input", type=str,  help='Input fits')
     parser.add_argument("--outdir", type=str, default='.',  help='Output directory')
     parser.add_argument("--show", default='True',type=str, help='Show result key. True or False.')
     parser.add_argument('--save_temps', default='True', type=str, help='Save temp files key. True or False')

     args = parser.parse_args()
     file = args.input
     if os.path.exists(file):
        print('Input file:', file)
     else:
        print('Err! No such file:', file)
        os._exit(2)

     outdir = args.outdir
     if os.path.exists(outdir):
        print('Out dir:', outdir)
     else:
        key = True
        while key:
            ans = input("No path %s. Do you want to create path? (y/n) \n" %outdir)

            if (ans == "y"):
                subprocess.run("mkdir -p %s" %outdir, shell=True)  
                key = not(key)
            elif (ans == "n"):
                print("Exit...")
                key = not(key)
                os._exit(1)

     show = args.show == 'True'
     print('Show key:', show)
     save = args.save_temps == 'True'
     print('Save temp files:', save)
     
     start = time()
     print('==== START ====')
     print('Creates image mask')     
     subprocess.run(['python3', 'raw_mask.py', '-input', file])
      
     print('Gets contures from mask')
     subprocess.run(['python3', 'mask2contures.py'])
    
     print('Ellipse Fitting')
     image = fits.getdata(file)
     data  = fits.getdata('temp_mask.fits')
     jf = 'temp_contures.json'
     pars, Ps = fit_ell.fit_ell(jf, data)
     save_pars_to_data(pars, Ps, file=Path(outdir,'fit_result.dat'))
     if show:
         fit_ell.plot_fit(data, image, pars, Ps, outdir)
     
     if save:
        subprocess.run(['mv temp_* %s' %outdir], shell=True)
     else:
        subprocess.run(['rm temp_*'], shell=True)
     print('Done! Time: %8.1f s.' %(time() -start))
