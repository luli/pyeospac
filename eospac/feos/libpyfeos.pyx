#!/usr/bin/python
# -*- coding: utf-8 -*-
# cython: boundscheck=False
# cython: cdivision=True
# cython: wraparound=False

import numpy as np
cimport numpy as np

from libc.stdlib cimport malloc, free
# to use the Numpy-C-API from Cython
np.import_array()

ctypedef np.int32_t int32
ctypedef np.float64_t float64

cdef extern from "libfeos.h":
    void FEOS_Initialize( int matnum, int printflag)
    void FEOS_Finalize( )
    void FEOS_Get_Calc_Limits( double*, double* )
    void FEOS_Init_Mat( int, int, int, int, double*, int, int*, int*, double* )
    void FEOS_Delete_Mat( int )
    void FEOS_Get_Mat_Par( int, double*, double*, double*, double*, double*, double*,
                        double*, double*, double*, int* )
    void FEOS_Get_SoftSphere_Par( int, double*, double*, double*, double*, double* )
    void FEOS_Get_Energy_Offsets( int, double*, double* )
    void FEOS_Get_Crit_Point( int, double*, double*, double*, double*, double*, double* )
    void FEOS_Get_Binodal( int, double*, double*, double*, double*, double*, double*,
                        double*, double*, double*, double*, double*, double*, double* )
    void FEOS_Get_Spinodal( int, double*, double*, double*, double*, double* )
    void FEOS_Get_EOS_All( int, int, double, double , int*,
                                            double*, double*, double*, double*, double*, double*,
                        double*, double*, double*, double*,
                        double*, double*, double*, double*,
                        double*, double*, double*, double* )
    void FEOS_Get_EOS( int, int, int, double, double, int*,
                                        double*, double*, double*, double*, double*, double* )
    void FEOS_Get_Pmin_Pmax( int, double,
                                            double*, double*, double*, double* )


cpdef _initialize(int matnum, int printflag):
    FEOS_Initialize(matnum, printflag)

cpdef _finalize():
    FEOS_Finalize()

cpdef _get_calc_limits():
    cdef double t_zero
    cdef double rho_zero
    FEOS_Get_Calc_Limits(&t_zero, &rho_zero)
    return t_zero, rho_zero

cpdef _init_mat(int table_handle, int matid, int use_maxwell,
        int use_softsphere,
        np.ndarray[float64, ndim=1, mode='c'] maxwell_temp_arr):
    cdef int Npt_maxwell_temp_arr
    cdef int num_elements
    cdef int num_maxwell_iso
    cdef double maxwell_is_reliable
    Npt_maxwell_temp_arr = len(maxwell_temp_arr)
    FEOS_Init_Mat(table_handle, matid, use_maxwell, use_softsphere,
            <double *> np.PyArray_DATA(maxwell_temp_arr), Npt_maxwell_temp_arr,
            &num_elements, &num_maxwell_iso,
            &maxwell_is_reliable)
    return num_elements, num_maxwell_iso, maxwell_is_reliable


cpdef _get_eos_all(int table_handle, int num_elements, int use_maxwell,
        double [::1] rho, double [::1] temp):
    """
    Vectorized version of Get_EOS_All
    """
    cdef int k, Npt
    cdef double [::1] Pt, Ut, St, At, Zbar,\
                        Pe, Ue, Se, Ae,\
                        Pi, Ui, Si, Ai,\
                        P_tf, U_tf, S_tf, A_tf
    cdef int [::1] bellow_binodal
    cdef dict res
    cdef np.ndarray[float64, ndim=2, mode='c'] Zi

    if len(rho) != len(temp):
        raise ValueError('Rho and temp arrays should be of the same shape!')
    Npt = len(rho)
    Pt, Ut, St, At = np.zeros(Npt), np.zeros(Npt), np.zeros(Npt), np.zeros(Npt)
    Pe, Ue, Se, Ae = np.zeros(Npt), np.zeros(Npt), np.zeros(Npt), np.zeros(Npt)
    Pi, Ui, Si, Ai = np.zeros(Npt), np.zeros(Npt), np.zeros(Npt), np.zeros(Npt)
    P_tf, U_tf, S_tf, A_tf = np.zeros(Npt), np.zeros(Npt), np.zeros(Npt), np.zeros(Npt)
    Zbar, Zi = np.zeros(Npt), np.zeros((Npt, num_elements))
    bellow_binodal = np.zeros(Npt, 'int32')

    E_offsets = _get_energy_offsets(table_handle)

    for k in range(Npt):
        FEOS_Get_EOS_All(table_handle, use_maxwell, rho[k], temp[k], # inputs
                &bellow_binodal[k], &Pt[k], &Ut[k], &St[k], &At[k],\
                <double *> np.PyArray_DATA(Zi[k]), &Zbar[k],\
                &Pe[k], &Ue[k], &Se[k], &Ae[k],\
                &Pi[k], &Ui[k], &Si[k], &Ai[k],\
                &P_tf[k], &U_tf[k], &S_tf[k], &A_tf[k])

    #res =  {'Pt': Pt.base,'Ut': Ut.base, 'St': St.base, 'At': At.base,
    #        'Piz': Pi.base,'Uiz': Ui.base, 'Siz': Si.base, 'Aiz': Ai.base,
    #        'Pec': Pe.base,'Uec': Ue.base, 'Sec': Se.base, 'Aec': Ae.base,
    #        # pure the electronic part is given by the TF without cold
    #        # curve corrections
    #        'Pe': P_tf.base, 'Ue': U_tf.base, 'Se': S_tf.base, 'Ae': A_tf.base,
    #        #'Pe': Pi.base, 'Ue': Ui.base, 'Se': Si.base, 'Ae': Ai.base,
    #        'Pic': Pt.base - P_tf.base ,'Uic': Ut.base - U_tf.base,
    #        'Sic': St.base - S_tf.base, 'Aic': At.base - A_tf.base,
    #        'Zfc': Zbar.base, 'Zi': Zi, 'bellow_binodal': bellow_binodal.base}
    res =  {'Pt': Pt.base,'Ut': Ut.base, 'St': St.base, 'At': At.base,
            'Piz': Pi.base,'Uiz': Ui.base, 'Siz': Si.base, 'Aiz': Ai.base,
            'Pe': Pe.base,'Ue': Ue.base, 'Se': Se.base, 'Ae': Ae.base,
            # pure the electronic part is given by the TF without cold
            # curve corrections
            #'Pe': P_tf.base, 'Ue': U_tf.base, 'Se': S_tf.base, 'Ae': A_tf.base,
            #'Pe': Pi.base, 'Ue': Ui.base, 'Se': Si.base, 'Ae': Ai.base,
            'Pic': Pt.base - P_tf.base ,'Uic': Ut.base - U_tf.base,
            'Sic': St.base - S_tf.base, 'Aic': At.base - A_tf.base,
            'Zfc': Zbar.base, 'Zi': Zi, 'bellow_binodal': bellow_binodal.base}
    #for key in res.keys():
    #    print key, np.min(res[key])
    #for key in E_offsets:
    #    print key, E_offsets[key]
    return res


cpdef _get_eos(int table_handle, int num_elements, int task, int use_maxwell,
        double [::1] rho, double [::1] temp):
    """
    Vectorized version of Get_EOS
    """
    cdef int k, Npt
    cdef double [::1] P, U, S, A, Zbar
    cdef int [::1] bellow_binodal
    cdef dict res
    cdef np.ndarray[float64, ndim=2, mode='c'] Zi

    if len(rho) != len(temp):
        raise ValueError('Rho and temp arrays should be of the same shape!')
    Npt = len(rho)
    P, U, S, A = np.zeros(Npt), np.zeros(Npt), np.zeros(Npt), np.zeros(Npt)
    Zbar, Zi = np.zeros(Npt), np.zeros((Npt, num_elements))
    bellow_binodal = np.zeros(Npt, 'int32')

    for k in range(Npt):
        FEOS_Get_EOS(table_handle, task, use_maxwell, rho[k], temp[k], # inputs
                &bellow_binodal[k], &P[k], &U[k], &S[k], &A[k],\
                <double *> np.PyArray_DATA(Zi[k]), &Zbar[k])

    res =  {'P': P.base,'U': U.base, 'S': S.base, 'A': A.base,
            'Zbar': Zbar.base, 'Zi': Zi, 'bellow_binodal': bellow_binodal.base}
    return res


cpdef _delete_material(int table_handle):
    FEOS_Delete_Mat(table_handle)


cpdef _get_mat_par(int table_handle, int num_elements):
    cdef np.ndarray[float64, ndim=1, mode='c'] A, Z, X
    cdef double Atot, Ztot, Xtot, temp_ref, rho_ref, bulk_mod_ref
    cdef int matid
    A, Z, X = np.zeros(num_elements), np.zeros(num_elements), np.zeros(num_elements)

    FEOS_Get_Mat_Par(table_handle, <double *> np.PyArray_DATA(A),
            <double *> np.PyArray_DATA(Z), <double *> np.PyArray_DATA(X),
            &Atot, &Ztot, &Xtot, &temp_ref, &rho_ref, &bulk_mod_ref, &matid)
    return {'A': A, 'Z': Z, 'X': X,
            'Atot': Atot, 'Ztot': Ztot, 'Xtot': Xtot,
            'temp_ref': temp_ref, 'rho_ref': rho_ref, 'bulk_mod_ref': bulk_mod_ref,
            'matid': matid}


cpdef _get_crit_point(int table_handle):
    cdef double temp, rho, P, H, S, critical_compressibility
    FEOS_Get_Crit_Point(table_handle, &temp, &rho,
            &P, &H, &S, &critical_compressibility)
    return {'critical_temp': temp, 'critical_rho': rho,'critical_P': P,
            'critical_H': H, 'critical_S': S, 
            'critical_compressibility': critical_compressibility}

cpdef _get_soft_shere_par(int table_handle):
    cdef double Ecoh, n, m, A, B
    FEOS_Get_SoftSphere_Par(table_handle, &Ecoh, &n, &m,
            &A, &B)
    return {'softsph_Ecoh': Ecoh, 'softsph_n': n, 'softsph_m': m, 'softsph_A':A, 'softsph_B':B}

cpdef _get_energy_offsets(int table_handle):
    cdef double electron_offs, ion_offs
    FEOS_Get_Energy_Offsets(table_handle, &electron_offs, &ion_offs)
    return {'offs_electron': electron_offs, 'offs_ion': ion_offs}

# to implement Binodal, et Spinodal routines
