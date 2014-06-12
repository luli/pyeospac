#!/usr/bin/python
# -*- coding: utf-8 -*-

from evtk.hl import gridToVTK, grid2DToVTK
import numpy as np
import random as rnd

def save2vtk(self, filename, tables=[], quantities=[], z_pos=False,
            grid_pos=True, grid_log=True,
            rho=None, temp=None, spec='t'):
    """
    Save specified tables and quantities of EosMaterial to a vtk format
    for vizualization with VisIt

    Parameters:
    ----------
     - filename: path of the file that will be created
     - tables: list of tables
     - quantities: list of quantities
     - zpos: enforce pozitivity for all fields
     - shift_zero_pt: shift zero point for density and temperature point

    Warning: this probably doensn't work so well at the moment for 3T SESAME, have
    to check what 'spec' variable is doing exactly
    """
    if type(spec) is str:
        spec = [spec]

    if not tables and not quantities:
        raise ValueError('At least table list or quantities list should be provided')
    if rho is None and temp is None:
        grid_provided = False
        for idx, tab_name in enumerate(tables):
            tab = self.get_table(tab_name, spec[0])
            if idx == 0:
                rho = tab['R_Array']
                temp = tab['T_Array']
            else:
                if not np.allclose(rho, tab['R_Array']):
                    raise ValueError("Provided tables {0}, {1} don't use the same density grid!".format(tables[0], tab_name))
                if not np.allclose(temp, tab['T_Array']):
                    raise ValueError("Provided tables {0}, {1} don't use the same temperature grid!".format(tables[0], tab_name))

        # getting a TD grid even if no tables are providied
        if not tables:
            tab_name = filter(lambda x: spec[0]+'_DT' in x, self.tables)[0] 
            tab = self.get_table(tab_name, spec[0])
            rho = tab['R_Array']
            temp = tab['T_Array']
    else:
        grid_provided = True
    if grid_pos or grid_log:
        if rho[0] == 0.0:
            rho[0] = 0.5*rho[1]
        if temp[0] == 0.0:
            temp[0] = 0.5*temp[1]

    # converting 1D grid to 3D
    rho_grid, temp_grid = np.meshgrid(rho, temp, indices='ij')
    rho_grid = rho_grid.reshape(rho_grid.shape + (1,))
    temp_grid = temp_grid.reshape(temp_grid.shape + (1,))
    z_grid = np.zeros(temp_grid.shape) # setting the 3rd coordinate to zero
    data = {}
    for spec_el in spec:
        for tab_name in tables:
            tab_name = tab_name.format(s=spec_el)
            tab = self.get_table(tab_name, spec_el) # redondant but harmless
            if not grid_provided:
                data[tab_name] = tab['F_Array'].reshape(rho_grid.shape).T
            else:
                data[tab_name] = tab(rho_grid, temp_grid)
        for q_name in quantities:
            data[q_name+'_'+spec_el] = self.q[q_name, spec_el](rho_grid, temp_grid)

    data_out = {}
    for key in data:
        vals = data[key]
        if z_pos:
            neg_mask = vals <= 0
            vals[neg_mask] = vals[~neg_mask].min()

        data_out['eospac_'+key] = vals
    if grid_log:
        rho_grid = np.log10(rho_grid)
        temp_grid = np.log10(temp_grid)

    gridToVTK(filename, rho_grid, temp_grid, z_grid,
                        cellData =None, pointData = data_out)


