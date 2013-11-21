#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
cimport numpy as np

import constants as cst
from libc.stdlib cimport malloc, free
# to use the Numpy-C-API from Cython
np.import_array()


ctypedef np.int32_t int32
ctypedef np.float64_t float64

cdef extern from "ses_defines.h":
    ctypedef int          ses_file_handle
    ctypedef char*        ses_string
    ctypedef char         ses_open_type
    ctypedef int          ses_error_flag
    ctypedef long         ses_material_id
    ctypedef long*        ses_material_id_reference
    ctypedef long         ses_table_id
    ctypedef long*        ses_table_id_reference
    ctypedef long*        ses_table_sizes_reference
    ctypedef char*        ses_label
    ctypedef long         ses_number
    ctypedef long*        ses_number_reference
    ctypedef int          ses_boolean
    ctypedef double       ses_word
    ctypedef double*      ses_word_reference
    ctypedef char         ses_file_type
    ctypedef int          ses_array_order

    #  open, setup, and close, and exit


    cdef ses_error_flag   ses_close(ses_file_handle the_handle)
    cdef ses_boolean      ses_exit()
    cdef ses_file_handle  ses_open(ses_string filename, ses_open_type open_flags)
    cdef ses_error_flag   ses_setup(ses_file_handle the_handle, 
                                        ses_material_id the_mid, ses_table_id the_tid)
    cdef ses_error_flag   ses_write_setup(ses_file_handle the_handle,
            ses_material_id the_mid, ses_table_id the_tid, ses_number nr,
            ses_number nt, ses_number ntab)

    #  write array interface
    cdef ses_error_flag   ses_write_1D(ses_file_handle the_handle, ses_word_reference the_buffer, ses_number dim)
    cdef ses_error_flag   ses_write_2D(ses_file_handle the_handle, ses_word_reference the_buffer,
                                    ses_number dim1, ses_number dim2)
    cdef ses_error_flag   ses_write_3D(ses_file_handle the_handle, ses_word_reference the_buffer,
                                    ses_number dim1, ses_number dim2, ses_number dim3)
    cdef ses_error_flag   ses_write_comments(ses_file_handle the_handle, ses_string the_comments, ses_number dim)
    cdef ses_error_flag   ses_write_next(ses_file_handle the_handle, ses_word_reference the_buffer,
                                    ses_number dim, ses_label the_label)
    cdef ses_error_flag   ses_write_number(ses_file_handle the_handle, ses_number the_buffer)
    cdef ses_error_flag   ses_write_pairs(ses_file_handle the_handle, ses_word_reference buf1,
                                    ses_word_reference buf2, ses_number dim)
    cdef ses_error_flag   ses_write_word(ses_file_handle the_handle, ses_word the_buffer)


   #  setter routines 

    cdef ses_error_flag   ses_set_array_order(ses_file_handle the_handle, ses_array_order the_order)
    cdef ses_error_flag   ses_set_date(ses_file_handle the_handle, long the_date)
    cdef ses_error_flag   ses_set_format(ses_file_handle the_handle, ses_file_type the_file_type)
    cdef ses_error_flag   ses_set_grid(ses_file_handle the_handle, ses_number nr, ses_number nt, ses_number ntab)
    cdef ses_error_flag   ses_set_label(ses_file_handle the_handle, ses_label the_label)
    cdef ses_error_flag   ses_set_material_order(ses_file_handle the_handle)
    cdef ses_error_flag   ses_set_significant_digits(ses_file_handle the_handle, ses_number number_digits)
    cdef ses_error_flag   ses_set_validate(ses_file_handle the_handle)
    cdef ses_error_flag   ses_set_version(ses_file_handle the_handle, long the_version)

    # getter routines
    cdef ses_material_id_reference ses_get_materials(ses_file_handle the_handle, long* size);




cpdef _write_sesbin(ses_string filename, ses_material_id matid, dict prop, dict tab_dict):
    """ Write binary file for a single material

    Parameters
    ----------
      - filename
      - matid: material id
      - prop: a dicttionary with table properties with the following keys
            ['Abar', 'Zbar', 'rho_ref'] 
      - tabs: a list of dicts containing the following information
            'name': Pt, Ut, At, Zt, etc.. all are assumed to be in DT plane
            'shape': 1 or 2
            'comments': str
            'D_Array': density grid 1d array
            'T_Array': temperature grid 1D array
            'P_Array': data for the given variable 2D
            'U_Array': data for the given variable 2D
            'F_Array': data for the given variable 2D
    """
    cdef ses_file_handle f
    cdef ses_table_id tid
    cdef ses_open_type wopen_flag = 'W'
    cdef ses_file_type file_type = 'B'

    f = ses_open(filename, wopen_flag)
    ses_set_format(f, file_type)

    # Write 201 table
    tid = 201
    ses_write_setup(f, matid, tid, 5, 1, -1)
    for key in ['Mean_Atomic_Num', 'Mean_Atomic_Mass',
                'Normal_Density', 'Modulus', 'Exchange_Coeff']:
        if key in prop:
            val = prop[key]
        else:
            val = 0.
        ses_write_word(f, <ses_word> val)

    # Write 301 table
    for tid in sorted(tab_dict):
        if tid in [301, 303, 304, 305]:
            tab = tab_dict[tid]
            nr = len(tab['rho'])
            nt = len(tab['temp'])
            ses_write_setup(f, matid, tid, nr, nt, 3)
            ses_write_number(f, nr)
            ses_write_number(f, nt)
            ses_write_1D(f, <ses_word_reference> np.PyArray_DATA(tab['rho']), nr)
            ses_write_1D(f, <ses_word_reference> np.PyArray_DATA(tab['temp']), nt)
            for var in ['P', 'U', 'A']:
                ses_write_2D(f, <ses_word_reference> np.PyArray_DATA(
                    tab[var].astype('float64', order='C')),
                    nr, nt)
        elif tid in [502, 503, 504, 505, 601, 602, 603, 604, 605]:
            tab = tab_dict[tid]
            rho, temp = tab['rho'].copy(), tab['temp'].copy()
            #if rho[0]  == 0.0: rho[0] = 1e-10
            #if temp[0] == 0.0: temp[0] = 1e-10
            nr = len(rho)
            nt = len(temp)
            ses_write_setup(f, matid, tid, nr, nt, 1)
            ses_write_number(f, nr)
            ses_write_number(f, nt)
            ses_write_1D(f, <ses_word_reference> np.PyArray_DATA((rho)), nr)
            ses_write_1D(f, <ses_word_reference> np.PyArray_DATA((temp)), nt)
            ses_write_2D(f, <ses_word_reference> np.PyArray_DATA(
                (tab['Z']).astype('float64', order='C')),
                nr, nt)
        else:
            raise NotImplemented('Table {0} not implemened!'.format(tid))



    ses_close(f)


#def _get_materials(ses_string filename, float size):
#    cdef ses_file_handle f
#    cdef ses_open_type wopen_flag = 'R'
#    cdef ses_file_type file_type = 'B'
#
#    f = ses_open(filename, wopen_flag)
#    ses_set_format(f, file_type)
#
#    ses_get_materials(f, &size)







