#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
from scipy.optimize import brentq

def analyse_isotherm(tab, rho_in, temp):
    if len(rho_in) == 2:
        rho_min, rho_max = rho_in
        rho_vect = np.logspace(np.log10(rho_max),  np.log10(rho_min), 400)
    else:
        rho_vect = rho_in
        rho_min, rho_max = rho_vect.min(), rho_vect.max()
    temp_vect = np.ones(rho_vect.shape)*temp
    Pt_DT_dFx = tab.Pt_DT.dFx(rho_vect, temp_vect)
    Pt_DT_dFxx = tab.Pt_DT.dFxx(rho_vect, temp_vect)
    Pt_DT_dFx_func = lambda x: tab.Pt_DT.dFx(x, temp)
    Pt_DT_dFxx_func = lambda x: tab.Pt_DT.dFxx(x, temp)
    mask_dFx = np.diff((Pt_DT_dFx>0)*1)
    mask_dFxx = np.diff((Pt_DT_dFxx>0)*1)
    out_dFx_p = []
    out_dFx_n = []
    out_dFxx_p = []
    out_dFxx_n = []
    for mmask, func, out in [(mask_dFx<0, Pt_DT_dFx_func, out_dFx_n),
                            (mask_dFx>0, Pt_DT_dFx_func,  out_dFx_p),
                            (mask_dFxx<0, Pt_DT_dFxx_func,  out_dFxx_n),
                            (mask_dFxx>0, Pt_DT_dFxx_func,  out_dFxx_p)]:
        idx = np.nonzero(mask_dFx<0)
        for idx in np.nonzero(mmask)[0]:
            rho_l, rho_r = rho_vect[[idx,idx+1]]
            f_l, f_r =  func(rho_l), func(rho_r)
            if np.isnan(f_l*f_r):
                continue
            elif f_l*f_r>0: # the sign can't be the same!!
                continue

            root = brentq(func, rho_l, rho_r)
            out.append(root)

    return {'dFx': (out_dFx_n, out_dFx_p),
            'dFxx': (out_dFxx_n, out_dFxx_p)}


def get_critical_point(eosmat, x, y, level=0.0, x_min=None, x_max=None, y_min=None, y_max=None):
    from scipy.spatial.distance import cdist
    # get a contours for Pt_DT.dFx == 0 and Pt_DT.dFxx == 0
    Pt_x = get_contour(eosmat.Pt_DT.dFx, x, y, 0.0, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max)
    Pt_xx = get_contour(eosmat.Pt_DT.dFxx, x, y, 0.0, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max)
    # find roughly the intesection point
    distance = cdist(np.log(Pt_x), np.log(Pt_xx), 'euclidean')
    xidx, yidx =  np.unravel_index(np.argmin(distance), distance.shape)
    # localy fit both contours by a polynomial to find the intersection exactly.

    vals0 = Pt_x[xidx+np.arange(-1,2), :]
    vals1 = Pt_xx[yidx+np.arange(-1,2), :]
    poly0 = np.polyfit(vals0[:,0], vals0[:,1], 2)
    poly1 = np.polyfit(vals1[:,0], vals1[:,1], 2)

    # bound for the root intersection
    rho_c_l = vals0[:,0].min()
    rho_c_r = vals0[:,0].max()

    roots = np.roots(poly0 - poly1)
    for root in roots:
        if not np.iscomplex(root) and (rho_c_l<root<rho_c_r):
           rho_c = root
           break
    else:
        raise ValueError("Couldn't find intersection of curves dPt_DT.dFx==0 and dPt_DT.dFxx==0 !")

    temp_c = np.polyval(poly0, rho_c)
    pres_c = float(eosmat.Pt_DT(rho_c, temp_c))
    if pres_c < 0:
        raise ValueError("The critical point can't have a negative pressure!")

    return {'pres_c': pres_c, 'rho_c': rho_c, 'temp_c': temp_c}



def get_contour(func, x, y, level=0.0, x_min=None, x_max=None, y_min=None, y_max=None):
    """
    Calculate contour for a  variable

    Parameters:
    -----------
      - func : for instace EosMaterial.Pt_DT or EosMaterial.Pt_DT.dFxx
      - x [ndarray] (N,) e.g. density 1D array
      - y [ndarray] (N,) e.g. temperature  1D array
    """
    from skimage.measure import find_contours
    x_mask = ~np.isnan(x)
    y_mask = ~np.isnan(y)
    if x_min:
        x_mask *= (x>x_min)
    if x_max:
        x_mask *= (x<x_max)
    if y_min:
        y_mask *= (y>y_min)
    if y_max:
        y_mask *= (y<y_max)

    x = x[x_mask]
    y = y[y_mask]
    X, Y = np.meshgrid(x, y, indices='ij')

    M =  func(X, Y).T
    res = find_contours(M, level=level, fully_connected='low')
    out = []
    for el in res:
        out.append(np.array([
                np.interp(el[:,0], range(len(x)), x),
                np.interp(el[:,1], range(len(y)), y)]).T)

    return sorted(out, key=lambda x: x.shape[0])[-1]



