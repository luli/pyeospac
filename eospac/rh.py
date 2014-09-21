#!/usr/bin/python
# -*- coding: utf-8 -*-
# Copyright CNRS 2012
# Roman Yurchak (LULI)
# This software is governed by the CeCILL-B license under French law and
# abiding by the rules of distribution of free software.

from .interface import EosMaterial
import numpy as np
from scipy.optimize import root
from pandas import DataFrame
from copy import deepcopy

def _gen_eos_args(state_i, eos_pars):
    """ Generate a dict that we can pass to EosInterface """
    pars = dict( state_i.items() + eos_pars.items())
    for key in ['u_s', 'u_2']:
        if key in pars:
            del pars[key]
    return pars

class RankineHugoniot:
    def __init__(self, state0, state1, backend='gamma', eos_pars={},
            root_opts={'method': 'hybr'}, check=True, material=None,
            **args):
        """
        Compute the hydrodynamic Rankine-Hugoniot relationships
        with different equations of state

        Parameters:
        -----------
          - state0: initial state
          - state1: shocked state
          - backend: eos backend
          - opt_method: optimization method to pass to scipy.optimize.root
          - eos_pars: EoS backend specific parameters
        """
        # Just making sure there are no weird side-effects anywhere...
        state0 = deepcopy(state0)
        state1 = deepcopy(state1)
        # Initialize EoS
        eos_itp = EosMaterial(backend=backend, tables=['[PU]t_DT', 'T_DUt', 'Pt_DUt',
                                                       'Ut_DPt', 'T_DPt'],
                material=material,
                units="eospac", # using eospac units, it allows a better conditioned problem
                options=eos_pars)
        self.eos = eos_itp
        # get all thermodynamic variables for the state 0
        state0 = eos_itp.get_state(**state0)
        if 'rho' in state1:
            state0, state1 = self._solve_ener_rho1(state0, state1, 
                    eos_itp, root_opts)
        elif 'pres' in state1:
            state0, state1 = self._solve_ener_pres1(state0, state1, 
                    eos_itp, root_opts)
        u_s, u_p = self._solve_impulse_mass(state0, state1)

        #self.check_consistency(state0, state1, u_s=u_s, u_p=u_p, 
        out = {'u_s': u_s, 'u_p': u_p}
        for state_idx, state in enumerate([state0, state1]):
            for key in state:
                out[key+str(state_idx)] = state[key]*(state_idx and 1 or np.ones(state1[key].shape))
        self.res = DataFrame(out)
        return

    @classmethod
    def solve(cls, *pargs, **vargs):
        """Solve Rankine Hugoniot conditions. See RankineHugoniot class for parameters.
        """
        x = cls(*pargs, **vargs)
        return x.res

    @staticmethod
    def _solve_ener_rho1(state0, state1, eos_itp, root_opts):
        """
        Implicitly solve the Rankine-Hugoniot energy equation 
              ε₁ - ε₀ = ½(p₁ + p₀) (1/ρ₀ - 1/ρ₁)
        relationship given ρ₁
        The huginiot is bivariate in ρ₁ so this method should be used with care!
        """
        # the normalized version of the energy equation is
        # ε₁/ε₀ - 1 = ½(p₁/ε₀ + p₀/ε₀) (1/ρ₀ - 1/ρ₁)
        rho0, pres0, eint0 = [state0[key] for key in ['rho', 'pres', 'eint']]
        rho1 = state1['rho']

        def objfunc_ener(x, rho0, pres0, eint0, rho1, eos_itp):
            """x = eint1/eint0 in order to condition the problem
            """
            cstate1 = eos_itp.get_state(rho=rho1, eint=x,
                            mode='DU')
            pres1= cstate1['pres']
            err =  (x - eint0) - 0.5*(pres0 + pres1)*(1/rho0 - 1/rho1)
            return err
        # just use 1 with appropriate shape as initial state
        x0 =  eos_itp.Ut_DT(rho1, 300*np.ones(rho1.shape))
        sol = root(objfunc_ener, x0, 
                args=(rho0, pres0, eint0, rho1, eos_itp),
                method=root_opts['method'])

        eint1 = sol.x
        mask_bad = (sol.x<1.0)+(sol.x==x0) # this means we were not able to find a solution
        rho1[mask_bad] = np.nan
        eint1[mask_bad] = np.nan

        state1 = eos_itp.get_state(rho=rho1, eint=eint1, mode='DU')
        # if len 1 array, convert to scalar
        if len(state1['eint']) == 1:
            for key in state1:
                if key != 'rho':
                    state1[key]  = state1[key][0]
        return state0, state1

    @staticmethod
    def _solve_ener_pres1(state0, state1, eos_itp, root_opts):
        """
        Implicitly solve the Rankine-Hugoniot energy equation 
              ε₁ - ε₀ = ½(p₁ + p₀) (1/ρ₀ - 1/ρ₁)
        relationship given P₁
        """
        rho0, pres0, eint0 = [state0[key] for key in ['rho', 'pres', 'eint']]
        pres1 = state1['pres']

        def objfunc_ener(x, rho0, pres0, eint0, pres1, eos_itp):
            """x = pres1 in order to condition the problem
            """
            cstate1 = eos_itp.get_state(rho=x, pres=pres1,
                            mode='DP')
            eint1= cstate1['eint']
            err =  (eint1 - eint0) - 0.5*(pres0 + pres1)*(1/rho0 - 1/x)
            return err
        # just use 1 with appropriate shape as initial state
        #x0 = 0.90*rho0*pres1**0#RankineHugoniot._get_analytical_pres2dens(state0, pres1, game=10)
        x0 = RankineHugoniot._get_analytical_pres2dens(state0, pres1, game=200)
        sol = root(objfunc_ener, x0,
                args=(rho0, pres0, eint0, pres1, eos_itp),
                method=root_opts['method'])

        rho1 = sol.x
        mask_bad = (sol.x<=x0) # this means we were not able to find a solution
        rho1[mask_bad] = np.nan
        pres1[mask_bad] = np.nan

        state1 = eos_itp.get_state(rho=rho1, pres=pres1, mode='DP')
        # if len 1 array, convert to scalar
        if len(state1['rho']) == 1:
            for key in state1:
                if key != 'pres':
                    state1[key]  = state1[key][0]
        return state0, state1

    @staticmethod
    def _get_analytical_pres2dens(state0, pres1, game):
        """
        Compute post-shock density given the post_shock pressure and the initial state.
        Usually used as initialization for a more realistic EoS.
        """
        pres0 = state0['pres']
        rho0 = state0['rho']
        rho1 = rho0*((game+1)*pres1 + (game-1)*pres0)/((game-1)*pres1 + (game+1)*pres0)
        return rho1


    @staticmethod
    def _solve_impulse_mass(state0, state1):
        s0, s1 = state0, state1

        # Now computing shock velocities and particle velocities. 
        u_s_2 = s1['rho']*(s1['pres'] - s0['pres'])/(s0['rho']*(s1['rho'] - s0['rho']))
        u_s = u_s_2**0.5
        u_p = (1 - s0['rho']/s1['rho'])*u_s
        return u_s, u_p

