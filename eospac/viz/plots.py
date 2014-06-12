#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from ..quantities import derived_quantities

K2eV = 1./11640

def eos_plot(mat, name, ax,spec='t', vmin=None, vmax=None,
        linthreshx=3e-5, linthreshy=0.1, nx=300, ny=350, xmax=None, ymax=None):
    """
    Plot a table or a derived quantity
    """

    DT_tables = [el for el in mat.tables if el.endswith('{s}_DT'.format(s=spec))]

    Rmax = np.min([mat.get_table(el)['Rmax'] for el in DT_tables])
    Rmin = np.max([mat.get_table(el)['Rmin'] for el in DT_tables])

    Tmax = np.min([mat.get_table(el)['Tmax'] for el in DT_tables])
    Tmin = np.max([mat.get_table(el)['Tmin'] for el in DT_tables])

    if Rmin == 0:
        Xarr = np.logspace(-5, np.log10(Rmax)-0.1, nx)
        Xarr[0] = 0
    else:
        Xarr = np.logspace(np.log10(Rmin), np.log10(Rmax)-0.1, nx)
    if Tmin == 0:
        Yarr = np.logspace(1, np.log10(Tmax)-0.1, ny)
        #Yarr[0] = 0
    else:
        Yarr = np.logspace(np.log10(Tmin)+0.1, np.log10(Tmax)-0.1, ny)
    print Yarr.min()

    X, Y = np.meshgrid(Xarr, Yarr, indexing='ij')

    if name in mat.tables:
        print 'ok'
        tab = mat.get_table(name, spec)
        F = tab(X,Y)
        if vmax is None:
            vmax = np.percentile(F, 99.5) 
        if vmin is None:
            vmin = np.percentile(F[F>0], 0.5)

        cs = ax.pcolormesh(X, Y*K2eV, F, cmap=plt.cm.jet, norm = LogNorm(),
                vmin=vmin, vmax=vmax)
        levels = np.arange(np.log10(F[F>0].min()), int(np.log10(F.max())))
        logF = np.log10(np.where(F>0, F, F[F>0].min()))
        cl = ax.contour(X, Y/11640, logF, levels, colors='k')
        plt.clabel(cl, fontsize=10, inline=False, fmt='%1.0d', use_clabeltext=True)
        plt.title('Table {0}: {1}'.format(tab['Material_ID'], name.replace('_', '\_')))
        cb = plt.colorbar(cs)
        if F.min()<0:
            min_label = '   (min {0:.0e} GPa)'.format(F.min())
        else:
            min_label = ''
        cb.set_label('{0} [{1}] {2}'.format(tab.label.replace('_', '\_'),
                tab.units, min_label))
        contour_at_1 = False
    else:
        F = mat.derived_quantity[name, spec](X,Y)
        if vmax is None:
            vmax = np.percentile(F, 99) 
        cs = ax.pcolormesh(X, Y*K2eV, F, vmin=0, vmax=vmax, cmap=plt.cm.jet)
        plt.title(ur'Table {0}: total: {1}'.format(
            mat.get_table(DT_tables[0], spec)['Material_ID'],
            derived_quantities[name]['description'].decode('utf8').replace('\n','')))
        plt.colorbar(cs)
        if name in ['gamma_3']:
            cl = ax.contour(X, Y/11640, F, [1, 2,3], colors='k')
        contour_at_1 = True
        if name in ['G3', 'gamma_3']:
            contour_at_1 = False
        if contour_at_1:
            cl = ax.contour(X, Y*K2eV, F,  [1], colors='white', linestyles='solid')


    cl = ax.contourf(X, Y*K2eV, F>0, [0,0.5], colors='white', hatches=['//'])

    ax.set_xscale('symlog', linthreshx=3e-5)
    ax.set_yscale('symlog', linthreshy=0.1)
    if xmax is None:
        ax.set_xlim(0, Xarr.max())
    else:
        ax.set_xlim(0, xmax)
    if ymax is None:
        ax.set_ylim(0, Yarr.max()*K2eV)
    else:
        ax.set_ylim(0, ymax)

    ax.set_xlabel(r'$\rho$ [g.cm$^{-3}$]')
    ax.set_ylabel(r'$T$ [eV]')
    return ax

def plot_eos_table(ax, mat, table_name, spec='t', vmin=None, vmax=None,
         nx=300, ny=350, xmax=None, ymax=None, xmin=None, ymin=None):
    """
    Plot an EoS table
    """

    table_name = table_name.format(s=spec)
    tab = mat.get_table(table_name)

    Rmin, Rmax = tab['Rmin'], tab['Rmax']
    Tmin, Tmax = tab['Tmin'], tab['Tmax']
    if xmin is not None:
        Rmin = xmin
    if ymin is not None:
        Tmin = ymin

    Xarr = np.logspace(np.log10(Rmin), np.log10(Rmax)-0.1, nx)
    Yarr = np.logspace(np.log10(Tmin), np.log10(Tmax)-0.1, ny)

    X, Y = np.meshgrid(Xarr, Yarr, indexing='ij')

    F = tab(X,Y)

    if vmax is None:
        vmax = np.percentile(F, 99.5) 
    if vmin is None:
        vmin = np.percentile(F[F>0], 0.5)

    cs = ax.pcolormesh(X, Y*K2eV, F, cmap=plt.cm.jet, norm = LogNorm(),
            vmin=vmin, vmax=vmax)
    if vmin is not None:
        levels = np.arange(int(np.log10(vmin)), int(np.log10(F.max())))
    else:
        levels = np.arange(np.log10(F[F>0].min()), int(np.log10(F.max())))
    logF = np.log10(np.where(F>0, F, F[F>0].min()))
    cl = ax.contour(X, Y/11640, logF, levels, colors='k')
    plt.clabel(cl, fontsize=10, inline=False, fmt='%1.0d', use_clabeltext=True)
    plt.title('Table {0}: {1}'.format(tab['Material_ID'], table_name.replace('_', '\_')))
    cb = plt.colorbar(cs)
    if F.min()<0:
        min_label = '   (min {0:.0e} GPa)'.format(F.min())
    else:
        min_label = ''
    cb.set_label('{0} [{1}] {2}'.format(tab.label.replace('_', '\_'),
            tab.units, min_label))

    cl = ax.contourf(X, Y*K2eV, F>0, [0,0.5], colors='white', hatches=['//'])

    ax.set_xscale('symlog', linthreshx=3e-5)
    ax.set_yscale('symlog', linthreshy=0.1)
    if xmax is None:
        ax.set_xlim(0, Xarr.max())
    else:
        ax.set_xlim(0, xmax)
    if ymax is None:
        ax.set_ylim(0, Yarr.max()*K2eV)
    else:
        ax.set_ylim(0, ymax)

    ax.set_xlabel(r'$\rho$ [g.cm$^{-3}$]')
    ax.set_ylabel(r'$T$ [eV]')
    return ax

