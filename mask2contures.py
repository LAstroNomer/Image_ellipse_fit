from astropy.io import fits
from matplotlib import pyplot as plt
import numpy as np 
import argparse
import json
import multiprocessing
import tqdm
from itertools import product
import warnings

def plot(data, all_data):
    plt.figure()
    plt.imshow(data, origin='lower')
    for ad in all_data:
        plt.plot(ad[0], ad[1], '-r')
    plt.savefig('temp_deblending_contures.png')

def faster(args):
    i = args[0]
    data = args[1]
    temp_data = np.zeros(data.shape)
    temp_data[np.where(data==i)] = 1
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cs = plt.contour(temp_data, levels=[1])  
    #del data
    #del temp_data
    plt.clf()
    plt.close()
    
    cols = cs.collections[0] 
    #del cs
    p0 = cols.get_paths()
    #print(p0)
    #del cols
    v = p0[0].vertices
    #del p0
    #v = p.vertices
    #del p
    x = v[:,0]
    y = v[:,1]
    return [list(x), list(y)]

def main(data):
    nmax = int(np.max(data))

    #all_data = []
    #for i in range(1, 3):
    
    hs = [[i, data] for i in range(1, nmax+1)]
    with multiprocessing.Pool(1) as pool:
        all_data = list(tqdm.tqdm(pool.imap(faster, hs), total=len(hs)))
        
        return all_data

def cre_dict(p):
    res = dict()
    
    for i in range(len(p)):
        res[i] = p[i]
    return res

if __name__ == '__main__':
    data = fits.getdata('temp_mask.fits')
    p = main(data)

    plot(data, p)        
    dic = cre_dict(p)
    with open('temp_contures.json' ,'w') as t:
        contures = json.dumps(dic, sort_keys=True, indent=4)
        t.write(contures)
