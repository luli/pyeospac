#!/usr/bin/python
# -*- coding: utf-8 -*-
# Copyright CNRS 2012
# Roman Yurchak (LULI)
# This software is governed by the CeCILL-B license under French law and
# abiding by the rules of distribution of free software.

from .interface import EosMaterial
import numpy as np
from scipy.integrate import ode
from pandas import DataFrame
from copy import deepcopy


# solving the equation
# dUt/drho = Pt/rho^2

def ode_func(rho, eint,  eos_itp):
    #dy/dt = f(y, t)
    args = np.array([rho]), eint
    Pt = eos_itp.Pt_DUt(np.array([rho]), eint)
    return Pt/rho**2

def ode_func_jacobian(eint, rho, eos_itp):
    #dy/dt = f(y, t)
    args = rho, np.array([eint])
    print args
    return eos_itp.Pt_DUt.dFy(*args)/rho**2


class Adiabat:
    def __init__(self, state0, rho, backend='gamma', eos_pars={},
            root_opts={'method': 'hybr'}, check=True, material=None,
            **args):
        """
        Compute the hydrodynamic Rankine-Hugoniot relationships
        with different equations of state

        Parameters:
        -----------
          - state0: initial state
          - rho_lim: the limit for the density.
          - backend: eos backend
          - opt_method: optimization method to pass to scipy.optimize.root
          - eos_pars: EoS backend specific parameters
        """
        state0 = deepcopy(state0)
        res = {'rho': deepcopy(rho)}#, 'eint': np.ones(len(rho))*np.nan}
        # Initialize EoS
        eos_itp = EosMaterial(backend=backend, tables=['[PU]t_DT', 'T_DUt', 'Pt_DUt',
                                                       'Ut_DPt', 'T_DPt', 'St_DT', 'T_DSt'],
                material=material,
                units="eospac", # using eospac units, it allows a better conditioned problem
                options=eos_pars)
        state0 = eos_itp.get_state(**state0)
        S0 = eos_itp.St_DT(state0['rho'], state0['temp'])
        if S0.shape:
            S0 = S0[0]
        S0_vect = S0*np.ones(res['rho'].shape)#*state0['rho']/res['rho']

        res['temp'] = eos_itp.T_DSt(res['rho'], S0_vect)/11640.


        state_res = eos_itp.get_state(**res)

        self.res = DataFrame(state_res)

    @classmethod
    def solve(cls, *pargs, **vargs):
        """Solve Rankine Hugoniot conditions. See RankineHugoniot class for parameters.
        """
        x = cls(*pargs, **vargs)
        return x.res
