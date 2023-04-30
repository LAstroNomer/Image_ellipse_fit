import matplotlib.pyplot as plt
import numpy as np

def main(data, out_mask):

    
    plt.figure()
    cs = plt.contour(data, levels=1) #‘solid’ | ‘dashed’ | ‘dashdot’ | ‘dotted’m = [[15,14,13,12],[14,12,10,8],[13,10,7,4],[12,8,4,0]]
    plt.clf()
    plt.close()

    with open(out_mask, 'w') as m:
        print('# Region file format: DS9 version 4.1', file=m)
        print('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1', file=m)
        print('image',  file=m)
        for cols in cs.collections:  
            p0 = cols.get_paths()
            for p in p0:
                v = p.vertices
                x = v[:,0]
                y = v[:,1]
    
                if len(x) < 50:
                    continue
                print('polygon(', end='', file=m)
                res = []
                for i, xi, yi in zip(range(len(x)), x, y):
                    res.append(str(np.round(xi,2)))
                    res.append(str(np.round(yi,2)))
                res = ','.join(res)
                print(res, file=m, end=') \n')

if __name__ == '__main__':   
    from astropy.io import fits
    data = fits.getdata('append/mask.fits')
    main(data, 'append/test.reg') 
