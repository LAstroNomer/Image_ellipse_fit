from astropy.io import fits
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
#from fourier import fourier
from scipy import interpolate
from scipy.stats import median_abs_deviation as mad
import multiprocessing
from pathlib import Path
from astropy.visualization import ZScaleInterval, PercentileInterval
from astropy.visualization import (MinMaxInterval, LogStretch,
                                   ImageNormalize)
import tqdm
import json 

def get_R(alpha):
    R = np.zeros((2,2))
    R[0,0] = R[1,1] = np.cos(alpha)
    R[0, 1] = - np.sin(alpha)
    R[1, 0] =   np.sin(alpha)
    return R

def rotate(x, y, alpha):
    
    R = get_R(alpha)

    res_x = np.array([])
    res_y = np.array([])
    
    for xi, yi in zip(x,y):
        x_, y_ = np.dot(R, [xi, yi])
        res_x = np.append(res_x, x_)
        res_y = np.append(res_y, y_)

    return res_x, res_y

def get_M(a11, a12, a22, a13, a23):
    
    M = np.zeros((2,2))
    M[0,0] = a11
    M[1,0] = M[0,1] = a12
    M[1,1] = a22
    #M[2,2] = 1
    #M[0,2] = M[2,0] = a13    
    #M[1,2] = M[2,1] = a23
    
    return M

def check_inner(data, x, y):
    
    xc = int(np.mean(x))
    yc = int(np.mean(y))

    if data[yc, xc] == 0:
        return False
    return True

def fit(x, y):
    
    X = x
    Y = y

    def cart_model(x, a11, a12, a22, a13, a23):
        return a11*X**2 + a22*Y**2 + 2*a12*X*Y + 2*a13*X + 2*a23*Y - 1

    p, covp = curve_fit(cart_model, xdata=X, ydata=np.zeros(len(X)))
    a11, a12, a22, a13, a23 = p
    #M = get_M(a11, a12, a22, a13, a23)

    theta = 0.5* np.arctan2(2*a12, a11 - a22)
    #R = get_R(theta)

    #M_ = np.dot(R.transpose(), np.dot(M, R))
    
    #for a in M_:
    #    print(a)

    #A = 1/np.sqrt(M[0,0])
    #B = 1/np.sqrt(M[1,1])
    
    D = a11*a22 - a12**2
    x0 = (a12*a23 - a13*a22)/D
    y0 = (a12*a13 - a11*a23)/D

    xc = x - x0
    yc = y - y0

    xr, yr = rotate(xc, yc, -theta)
    def simp_ell(x, a, b):
        return ((xr)/a)**2 + ((yr)/b)**2 - 1
    p, covp = curve_fit(simp_ell, xdata=x, ydata=np.zeros(len(x)))
    #print(p) 
    A = np.max([p[0], p[1]])
    B = np.min([p[1], p[0]])
    #e = np.sqrt(1 - B**2/A**2)
    P = np.sum(simp_ell(x, p[0], p[1])**2)/A

    #P = (np.sqrt(covp[0,0] + covp[1,1])-0.02)*np.hypot(A, B)/1.23 - 1.06 -2*e**2 + e
    #P = (np.sqrt(covp[0,0] + covp[1,1]) - 0.03)/e - 0.02027316
    #e = A*np.sqrt(1 -B**2/A**2) 

    #R = np.hypot(xr, yr)
    #Phi = np.arctan2(yr, xr)
    #Phi[Phi <0] += 2*np.pi
    #t = np.arange(0, 4*np.pi, 0.0001)
    #xir, yir = A*np.cos(t),  B*np.sin(t)
    #R_est = np.hypot(xir, yir)
    #Psi = np.arctan2(yir, xir)
    #Psi[Psi < 0] += 2*np.pi
    #Psi0 = np.array([phi for phi, dr in sorted(zip(Psi, R_est))])
    #R_est = np.array([dr for phi, dr in sorted(zip(Psi, R_est))])
    #print(Phi)
    

    #Phi0 = np.array([phi for phi, dr in sorted(zip(Phi, dR))])
    #dR = np.array([dr for phi, dr in sorted(zip(Phi, dR))])
    #try:    
    #    el = interpolate.interp1d(Psi0, R_est)
    #    R_est = el(Phi)

    #    dR = R - R_est
    #    p, nmax = fourier(dR, Phi) 
    #
    #    Ap = p[ : nmax]
    #    Bp = p[nmax :]
    #except:
    #    return 0, 0, False
    #def model(x, p, nmax):
    #    res_c = np.sum([p[i] * np.cos(i*x) for i in range(nmax)], axis=0)
    #    res_s = np.sum([p[nmax + i] * np.sin(i*x) for i in range(nmax)], axis=0)
    #    return res_c + res_s
    #delta_R = model(Phi, p, nmax) 
    #R_fit = R_est + delta_R
    
    #t = np.linspace(0, 2*np.pi, 1000)
    #de_R = model(t, p, nmax) 
    #plt.figure()
    #plt.plot(xr, yr, 'r--')
    #plt.plot(R_est*np.cos(Phi), R_est*np.sin(Phi), '--b')
    #plt.gca().set_aspect('equal')
    #plt.show()

    t = np.linspace(0, 2*np.pi, 1000)
    if (A > 1000) or (B > 1000):
        return 0 ,0 ,False
    ##print(theta, x0, y0, A, B)
    #xir, yir = rotate(A*np.cos(t+theta),  B*np.sin(t+theta), theta)
    #plt.plot( x0+xir,  y0+yir, '-r')
    #xir, yir =  R_fit*np.cos(Phi),  R_fit*np.sin(Phi)
    #xir, yir = rotate(xir, yir, theta)
    #plt.plot( x0 + xir,  y0 + yir, '-g')
    #print('%5.2f %5.2f %5.2f %5.2f %5.2f' %(x0, y0, A, B, theta-np.pi/2), file=f)
    #print('%5.3f' %(P), file=f)
    return [x0, y0, A, B, theta], P, True


def plot_fit(data, image, pars, Ps, outdir):
    from photutils.utils.colormaps import make_random_cmap
    from matplotlib import colors
    cmap = make_random_cmap(np.max(data) + 1, seed=111)
    cmap.colors[0] = colors.to_rgba('#000000ff')
    def mod(x, a, b, c):
        return a*x**2 + b*x + c
 
    plt.figure(figsize=(10,10), dpi=300)
    m = np.mean(data)
    s = np.std(data)
    plt.imshow(data, cmap=cmap, origin='lower')
    t = np.linspace(0, 2*np.pi, 100)
    ps_mean =  0.15 #np.median(Ps) + mad(Ps)
    #ps_mean =  np.mean(Ps) + np.std(Ps)
    rkey = ykey = True
    for p, ps in zip(pars, Ps):
        x0 = p[0]
        y0 = p[1]
        A  = p[2]
        B  = p[3]
        theta = p[4]
        if (A > 2000) or (B > 2000):
            continue
        xir, yir = rotate( A*np.cos(t - theta),  B*np.sin(t - theta), theta)
        if ps > ps_mean:
            if rkey:
                plt.plot(x0 + xir, y0 + yir, '-r', linewidth=2, label='bad')
                rkey = not(rkey)
            plt.plot(x0 + xir, y0 + yir, '-r', linewidth=2)
        else:
            if ykey:
                plt.plot(x0 + xir, y0 + yir, '-y', linewidth=2, label='good')
                ykey = not(ykey)
            plt.plot(x0 + xir, y0 + yir, '-y', linewidth=2)

    plt.legend(bbox_to_anchor=(1,1), loc="upper left")
    plt.tight_layout()
    plt.savefig(Path(outdir, 'result_ells.pdf'))


    plt.figure(figsize=(5, 5), dpi=300)
    m = np.mean(image)
    s = np.std(image)
    norm = ImageNormalize(image, interval=PercentileInterval(97))
    plt.imshow(image, cmap='Greys_r', norm=norm, origin='lower')
    t = np.linspace(0, 2*np.pi, 100)
    ps_mean =  0.15 #np.median(Ps) + mad(Ps)
    #ps_mean =  np.mean(Ps) + np.std(Ps)
    rkey = ykey = True
    for p, ps in zip(pars, Ps):
        x0 = p[0]
        y0 = p[1]
        A  = p[2]
        B  = p[3]
        theta = p[4]
        if (A > 2000) or (B > 2000):
            continue
        xir, yir = rotate( A*np.cos(t - theta),  B*np.sin(t - theta), theta)
        if ps > ps_mean:
            if rkey:
                plt.plot(x0 + xir, y0 + yir, '-r', linewidth=1, label='bad')
                rkey = not(rkey)
            plt.plot(x0 + xir, y0 + yir, '-r', linewidth=1)
        else:
            if ykey:
                plt.plot(x0 + xir, y0 + yir, '-g', linewidth=1, label='good')
                ykey = not(ykey)
            plt.plot(x0 + xir, y0 + yir, '-g', linewidth=1)

    plt.legend(bbox_to_anchor=(1,1), loc="upper left")
    plt.tight_layout()
    plt.savefig(Path(outdir, 'image_ells.pdf'))


    mp = np.max(Ps)
    width = 0.01
    nbins = int(mp/width)
    plt.figure(figsize=(5,5), dpi=300)
    plt.title('Image fit destribution')
    plt.hist(Ps, bins=nbins)
    plt.xlabel(r'$\frac{\varepsilon}{A}$')
    plt.ylabel('N')
    plt.tight_layout()
    plt.savefig(Path(outdir, 'fit_hist.png'))

    #plt.figure()
    #xs = []
    #ys = []
    #for p, ps in zip(pars, Ps):
    #    A = np.max([p[2], p[3]])
    #    B = np.min([p[2], p[3]])
    #    #x = np.sqrt(1 - B**2/A**2)
    #    x = A
    #    xs.append(x)
    #    ys.append(ps)
    #    plt.plot(x, ps, 'o')
    #xs = np.array(xs)
    #ys = np.array(ys)
    #p, covp = curve_fit(mod, xdata=xs, ydata=ys)
    ##print(p)
    #t = np.linspace(0, 1, 100)
    #plt.plot(t, mod(t, p[0], p[1], p[2]), '-b')
    #plt.xlabel('A')
    #plt.ylabel('p')
    #plt.show()

def faster(conture):
        #conture = contures[name]
        x = np.array(conture[0])
        y = np.array(conture[1])         
        if len(x) > 5:
            if True: #check_inner(data, x, y):
                par, P, key = fit(x, y)
                if key:
                    return [par, P]
        return [[0], [0]]

def fit_ell(jf, data):
    #f = open('fit_test.txt', 'w')
    with open(jf, 'r') as jf0:
        contures = json.load(jf0)
    #pars = []
    #Ps = []
    plt.figure()
    plt.imshow(data)

    cs = [contures[name] for name in contures]
    with multiprocessing.Pool(2) as pool:
        all_data = list(tqdm.tqdm(pool.imap(faster, cs), total=len(cs)))

    pars = [a for a, b in all_data  if len(a) > 1]  
    Ps = [b for a, b in all_data  if len(a) > 1]  
        
    return pars, Ps


def fit_ell0(data, show):
    #f = open('fit_test.txt', 'a')
    plt.figure() 
    cs = plt.contour(data, levels=100)    
    plt.close()
    pars = []
    Ps = []
    for cols in cs.collections: 
        p0 = cols.get_paths()
        
        for p in p0:
            v = p.vertices
            x = v[:,0]
            y = v[:,1]
            if len(x) > 5:
                if check_inner(data, x, y):
                    par, P, key = fit(x, y, show, f)
                    if key:
                        pars.append(par)
                        Ps.append(P)
        break
    return pars, Ps
import argparse

if __name__ == '__main__':
    #parser = argparse.ArgumentParser()
    ###
    #parser.add_argument("-input", type=str,  help='Input fits')
    #parser.add_argument("--show", default='True',type=str, help='Show result')

    #args = parser.parse_args()
    #file = args.input
    #show = args.show == 'True'
    file = 'out_temp_mask.fits'#'diff_deprojected.fits'
    data = fits.getdata(file)
    jf = 'test.json'
    show = True
    pars, Ps = fit_ell(jf, data)
    #pars, Ps = fit_ell0( data, show)
    if show:
        plot_fit(data, pars, Ps)
