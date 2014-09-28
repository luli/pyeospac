#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
from scipy.optimize import brentq
from scipy.interpolate import interp1d

def get_critical_point(eosmat, x, y, x_min=None, x_max=None, y_min=None, y_max=None, extra_fields={}):
    from scipy.spatial.distance import cdist
    # get a contours for Pt_DT.dFx == 0 and Pt_DT.dFxx == 0
    Pt = get_contour(eosmat.Pt_DT.dFx, x, y, 0.0, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max)
    Pt_x = get_contour(eosmat.Pt_DT.dFx, x, y, 0.0, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max, extra_fields=[eosmat.Pt_DT])
    Pt_xx = get_contour(eosmat.Pt_DT.dFxx, x, y, 0.0, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max, extra_fields=[eosmat.Pt_DT])

    # out of all the results select the point with the highest temperature
    
    rho_c, temp_c = get_curve_intersection(Pt_x[:,0], Pt_x[:,1], Pt_xx[:,0], Pt_xx[:,1], log=True)
    idx = np.argmax(temp_c) 

    rho_c = rho_c[idx]
    temp_c = temp_c[idx]

    pres_c = float(eosmat.Pt_DT(rho_c, temp_c))
    #if pres_c < 0:
    #    raise ValueError("The critical point can't have a negative pressure!")

    out = {'pres_c': pres_c, 'rho_c': rho_c, 'temp_c': temp_c}
    for key, func in extra_fields.iteritems():
        out[key] = float(func(out['rho_c'], out['temp_c']))
    return out


def saturation_curve(tab, rho, temp, extra_fields={}):
    crit_pt = get_critical_point(tab, rho, temp)
    rho_c = crit_pt['rho_c']
    temp_c = crit_pt['temp_c']
    temp = temp[temp<=temp_c]
    sat_l = []
    sat_v = []
    for idx, temp_val in enumerate(temp):
       res  = maxwell_construction(tab, rho, temp_val, rho_c, full_output=False)
       if res:
           sat_l.append([res['rho_l'], temp_val])
           sat_v.append([res['rho_v'], temp_val])
    sat_l.append([rho_c, temp_c])
    sat_v.append([rho_c, temp_c])


    sat_l = np.array(sat_l)
    sat_v = np.array(sat_v)

    Pneg = get_contour(tab.Pt_DT, rho, temp)
    mask =  (Pneg[:,1]<sat_l[:,1].min())*(Pneg[:,0]>rho_c)

    sat_l = np.concatenate((Pneg[mask,:], sat_l), axis=0)

    sat_l = {'rho': sat_l[:,0], 'temp': sat_l[:,1]}
    sat_v = {'rho': sat_v[:,0], 'temp': sat_v[:,1]}

    out = {'sat_l': sat_l, 'sat_v': sat_v}
    for maj_key in out:
        for key, func in extra_fields.iteritems():
            out[maj_key][key] = func(out[maj_key]['rho'], out[maj_key]['temp'])

    return  out



def maxwell_construction(tab, rho, temp_val, rho_c, full_output=False):
    """
    Compute Maxwell constractions on the isotherm
    Parameters:
    -----------
      - rho [ndarray] should be increasing
      - pres [ndarray]
      - temp [float]
    """
    if (np.diff(rho)<0).any():
        raise ValueError('Density array should be increasing!')
    temp = np.ones(rho.shape)*temp_val
    mask = rho>rho_c
    pres = tab.Pt_DT(rho, temp)
    A = tab.At_DT(rho, temp)
    G = A + pres/rho

    pres_m_eq, G_m_eq = get_curve_intersection(pres[mask], G[mask], pres[~mask], G[~mask])
    if len(pres_m_eq) != 1:
        return None

    dp_func = lambda x:  float(tab.Pt_DT(x, temp_val)) - pres_m_eq
    roots = global_roots(dp_func, rho)
    if len(roots) != 3: 
        # must have exactly 3 roots!
        return None
    rho_l, rho_r = roots[[0,-1]]
    out = {"pres_m_eq": pres_m_eq, 'rho_v': rho_l, 'rho_l': rho_r}
    if full_output:
        pres_m = np.array(pres)
        eint_m = tab.Ut_DT(rho, temp)
        eint_m_l = float(tab.Ut_DT(rho_l, temp_val))
        eint_m_r = float(tab.Ut_DT(rho_r, temp_val))
        mask = (rho>=rho_l)*(rho<=rho_r)
        pres_m[mask] = pres_m_eq
        eint_m[mask] = (eint_m_l + (eint_m_r - eint_m_l)*(1./rho - 1/rho_l)/(1./rho_r - 1./rho_l))[mask]
        out.update(**{'pres_m': pres_m, 'eint_m': eint_m})
    
    return out




def global_roots(func, x, extra_args=None):
    """
    Parameters:
    ----------- 
      - func: function of one argument
      - x [ndarray]: points at which the function will be evaluated
    """
    roots = []
    for idx in range(len(x)-1):
        x_l, x_r = x[idx], x[idx+1]
        if func(x_l)*func(x_r)<0:
                roots.append(brentq(func, x_l, x_r))
    return np.array(roots)




def get_curve_intersection(x1, y1, x2, y2, log=False):
    """Find all intersection of curves 1 and 2 defined by points x, y
    The curve 1 is assumed to be such that for every x value there is only one y value.
    """
    if log:
        x1 = np.log10(x1)
        y1 = np.log10(y1)
        x2 = np.log10(x2)
        y2 = np.log10(y2)

    itp = interp1d(x1, y1, kind='linear',  bounds_error=False)
    x_root = []
    y_root = []
    for idx in range(len(y2)-1):
        x_l, x_r = x2[[idx, idx+1]]
        y2_l, y2_r = y2[[idx, idx+1]]
        a = (y2_r - y2_l)/(x_r-x_l)
        b = y2_l - a*x_l
        err = lambda x: itp(x) - (a*x+b)
        if err(x_l)*err(x_r)<0: # there is an intersection
            x_sol = brentq(err, x_l, x_r)
            y_sol = a*x_sol+b
            x_root.append(x_sol)
            y_root.append(y_sol)
        else:
            continue
    x_root = np.array(x_root)
    y_root = np.array(y_root)
    if log:
        x_root = np.power(10, x_root)
        y_root = np.power(10, y_root)
    return x_root, y_root

def get_contour(func, x, y, level=0.0, x_min=None, x_max=None, y_min=None, y_max=None, extra_fields=[]):
    """
    Calculate contour for a  variable

    Parameters:
    -----------
      - func : for instace EosMaterial.Pt_DT or EosMaterial.Pt_DT.dFxx
      - x [ndarray] (N,) e.g. density 1D array
      - y [ndarray] (N,) e.g. temperature  1D array
    Returns:
    --------
       ndarray (N, 2+len(extra_fields)) the columns correspond to x, y, and extra_fields
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
        x_contour = el[:,0]
        y_contour = el[:,1]
        F = [ np.interp(x_contour, range(len(x)), x),
              np.interp(y_contour, range(len(y)), y)]
        for field in extra_fields:
            F.append(field(x_contour, y_contour))
        out.append(np.array(F).T)

    return sorted(out, key=lambda x: x.shape[0])[-1]



