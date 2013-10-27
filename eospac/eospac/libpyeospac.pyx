#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
cimport numpy as np

import constants as cst
from libc.stdlib cimport malloc, free
# to use the Numpy-C-API from Cython
np.import_array()

cdef extern from "eos_universal_types.h":
    ctypedef long EOS_INTEGER
    ctypedef unsigned long long_UNSIGNED
    ctypedef double EOS_REAL
    ctypedef char EOS_CHAR
    cdef enum EOS_BOOLEAN_enum:
        EOS_FALSE=0
        EOS_TRUE=1
        EOS_INVALID=-99999
    ctypedef EOS_BOOLEAN_enum EOS_BOOLEAN
    int EOS_MaxErrMsgLen
    int EOS_OK

ctypedef np.int32_t int32
ctypedef np.float64_t float64

cdef extern from "eos_Interface.h":
    cdef EOS_INTEGER EOS_Cmnt_Len
    cdef EOS_INTEGER EOS_INTERP_EXTRAPOLATED
    void eos_CheckExtrap (EOS_INTEGER *tableHandle,
                                   EOS_INTEGER *nXYPairs,
                                   EOS_REAL *xVals,
                                   EOS_REAL *yVals,
                                   EOS_INTEGER *xyBounds,
                                   EOS_INTEGER *errorCode)

    void eos_CreateTables (EOS_INTEGER *nTables,
                                    EOS_INTEGER tableType[],
                                    EOS_INTEGER matID[],
                                    EOS_INTEGER tableHandles[],
                                    EOS_INTEGER *errorCode)
    void eos_DestroyAll (EOS_INTEGER *errorCode)
    void eos_DestroyTables (EOS_INTEGER *nTables,
                                     EOS_INTEGER tableHandles[],
                                     EOS_INTEGER *errorCode)
    void eos_GetErrorCode (EOS_INTEGER *tableHandle,
                                    EOS_INTEGER *errorCode)
    void eos_GetErrorMessage (EOS_INTEGER *errorCode,
                                       EOS_CHAR errorMsg[])
    void eos_GetMaxDataFileNameLength (EOS_INTEGER *max_length)
    void eos_GetPackedTables (EOS_INTEGER *nTables,
                                       EOS_INTEGER tableHandles[],
                                       EOS_CHAR *packedTables,
                                       EOS_INTEGER *errorCode)
    void eos_GetPackedTablesSize (EOS_INTEGER *nTables,
                                           EOS_INTEGER tableHandles[],
                                           EOS_INTEGER *packedTablesSize,
                                           EOS_INTEGER *errorCode)
    void eos_Interpolate (EOS_INTEGER *tableHandle,
                                   EOS_INTEGER *nXYPairs,
                                   EOS_REAL *xVals,
                                   EOS_REAL *yVals,
                                   EOS_REAL *fVals,
                                   EOS_REAL *dFx,
                                   EOS_REAL *dFy, EOS_INTEGER *errorCode)
    void eos_LoadTables (EOS_INTEGER *nTables,
                                  EOS_INTEGER tableHandles[],
                                  EOS_INTEGER *errorCode)
    void eos_ResetOption (EOS_INTEGER *tableHandle,
                                   const EOS_INTEGER *tableOption,
                                   EOS_INTEGER *errorCode)
    void eos_SetDataFileName_Cwrapper (EOS_INTEGER *tableHandle,
                                                EOS_INTEGER *matID,
                                                EOS_INTEGER *tableType,
                                                EOS_CHAR *fileName,
                                                EOS_INTEGER *errorCode)

    void eos_SetPackedTables (EOS_INTEGER *nTables,
                                       EOS_INTEGER *packedTablesSize,
                                       EOS_CHAR *packedTables,
                                       EOS_INTEGER tableHandles[],
                                       EOS_INTEGER *errorCode)
    void eos_SetOption (EOS_INTEGER *tableHandle,
                                 const EOS_INTEGER *tableOption,
                                 const EOS_REAL *tableOptionVal,
                                 EOS_INTEGER *errorCode)
    void eos_GetTableInfo (EOS_INTEGER *tableHandle,
                                    EOS_INTEGER *numInfoItems,
                                    EOS_INTEGER *infoItems,
                                    EOS_REAL *infoVals,
                                    EOS_INTEGER *errorCode)
    void eos_GetTableCmnts (EOS_INTEGER *tableHandle,
                                     EOS_CHAR *cmntStr,
                                     EOS_INTEGER *errorCode)
    void eos_Mix (EOS_INTEGER *nTables, EOS_INTEGER *tableHandles,
                           EOS_INTEGER *nXYPairs, EOS_REAL *concInMix,
                           EOS_REAL *xVals, EOS_REAL *yVals, EOS_REAL *fVals,
                           EOS_REAL *dFx, EOS_REAL *dFy,
                           EOS_INTEGER *errorCode)
    void eos_GetVersionLength (EOS_INTEGER *length)
    void eos_GetVersion (EOS_CHAR *version)
    void eos_Time (EOS_BOOLEAN *reset, EOS_REAL *wctime,
                            EOS_REAL *cputime, EOS_REAL *cpucycles,
                            EOS_INTEGER *err)

cpdef _create_tables(np.ndarray[int32, ndim=1, mode='c'] table_type,
                    np.ndarray[int32, ndim=1, mode='c'] matid):
    cdef EOS_INTEGER ntables
    cdef EOS_INTEGER error_code
    cdef np.ndarray[int32, ndim=1, mode='c'] table_handles

    ntables = len(table_type)

    if len(table_type) != len(matid):
        raise RuntimeError('Both table_type and matid should have the same size!')

    table_handles = np.ones(ntables, dtype="int32")*-1
    eos_CreateTables(&ntables, <EOS_INTEGER*> np.PyArray_DATA(table_type),
                              <EOS_INTEGER*> np.PyArray_DATA(matid),
                              <EOS_INTEGER*> np.PyArray_DATA(table_handles),
                              &error_code)
    if error_code:
        _print_errors( table_handles )
        raise RuntimeError("Could not create files!")
    else:
        return table_handles

cpdef _load_tables(np.ndarray[int32, ndim=1, mode='c'] table_handles):
    cdef EOS_INTEGER ntables
    cdef EOS_INTEGER error_code

    ntables = len(table_handles)

    eos_LoadTables(&ntables, <EOS_INTEGER*> np.PyArray_DATA(table_handles),
                              &error_code)
    if error_code:
        _print_errors( table_handles )
        raise RuntimeError("Could not load files!")
    else:
        return 0

cpdef _get_error_code(long table_handle):
    cdef EOS_INTEGER error_code
    eos_GetErrorCode(<EOS_INTEGER*> &table_handle, &error_code)
    return error_code

cpdef _get_error_message(long error_code):
    cdef char* error_message
    error_message = <char *> malloc(EOS_MaxErrMsgLen*sizeof(char))
    eos_GetErrorMessage(<EOS_INTEGER*> &error_code,  <EOS_CHAR*> error_message)
    return error_message

cpdef _get_error_list(np.ndarray[int32, ndim=1, mode='c'] table_handles):
    messages_out = {}
    for handle in table_handles:
        error_code = _get_error_code(handle)
        if error_code:
            message = _get_error_message(error_code)
        else:
            message = "ok"
        messages_out[handle] = message
    return messages_out

cpdef _print_errors(np.ndarray[int32, ndim=1, mode='c'] table_handles):
    errors = _get_error_list(table_handles)
    print "EOSPAC traceback"
    print "="*16
    print ' handle    error'
    for handle, msg in errors.iteritems():
        print '   {0}   ->   {1}'.format(handle, msg)
    print "="*16
    return ''

cpdef _interpolate(long table_handle,
                  np.ndarray[float64, ndim=1, mode='c'] xVals,
                  np.ndarray[float64, ndim=1, mode='c'] yVals):
    cdef EOS_INTEGER error_code
    cdef EOS_INTEGER nXYPairs
    cdef np.ndarray[float64, ndim=1, mode='c'] fVals
    cdef np.ndarray[float64, ndim=1, mode='c'] dFx
    cdef np.ndarray[float64, ndim=1, mode='c'] dFy
    cdef np.ndarray[int32, ndim=1, mode='c'] xyBounds
    if  len(xVals) != len(yVals):
        raise ValueError('xVals and yVals arguments should be of the same shape!')
    nXYPairs = len(xVals)
    fVals = np.zeros((nXYPairs,), dtype='float64')
    dFx = np.zeros((nXYPairs,), dtype='float64')
    dFy = np.zeros((nXYPairs,), dtype='float64')

    eos_Interpolate (<EOS_INTEGER*> &table_handle,
                     <EOS_INTEGER*> &nXYPairs,
                     <EOS_REAL*> np.PyArray_DATA(xVals),
                     <EOS_REAL*> np.PyArray_DATA(yVals),
                     <EOS_REAL*> np.PyArray_DATA(fVals),
                     <EOS_REAL*> np.PyArray_DATA(dFx),
                     <EOS_REAL*> np.PyArray_DATA(dFy),
                     &error_code)
    if error_code:
        if error_code == EOS_INTERP_EXTRAPOLATED:
            xyBounds = _check_extrap(table_handle, xVals, yVals)
            for error in np.unique(xyBounds):
                if error:
                    default_error_labels =  cst.errors.keys()
                    default_error_ids =  np.array(cst.errors.values())
                    print 'Error {0} occured {1} times at indices {2}'.format(
                            default_error_labels[np.argmin(np.abs(default_error_ids - error))],
                            np.sum(xyBounds==error),
                            np.nonzero(xyBounds==error))


        raise RuntimeError( _get_error_message(error_code) )
    return fVals, dFx, dFy

cpdef _mix(np.ndarray[int32, ndim=1, mode='c'] table_handles,
                  np.ndarray[float64, ndim=1, mode='c'] xVals,
                  np.ndarray[float64, ndim=1, mode='c'] yVals,
                  np.ndarray[float64, ndim=2, mode='c'] concInMix2D,
                  ):
    cdef EOS_INTEGER error_code
    cdef EOS_INTEGER nXYPairs
    cdef EOS_INTEGER nTables
    cdef np.ndarray[float64, ndim=1, mode='c'] concInMix1D
    cdef np.ndarray[float64, ndim=1, mode='c'] fVals
    cdef np.ndarray[float64, ndim=1, mode='c'] dFx
    cdef np.ndarray[float64, ndim=1, mode='c'] dFy
    if  len(xVals) != len(yVals):
        raise ValueError('xVals and yVals arguments should be of the same shape!')
    nXYPairs = len(xVals)
    nTables = len(table_handles)
    fVals = np.zeros((nXYPairs,), dtype='float64')
    dFx = np.zeros((nXYPairs,), dtype='float64')
    dFy = np.zeros((nXYPairs,), dtype='float64')

    concInMix1D = concInMix2D.flatten()

    eos_Mix (&nTables,
             <EOS_INTEGER*> np.PyArray_DATA(table_handles),
             &nXYPairs,
             <EOS_REAL*> np.PyArray_DATA(concInMix1D),
             <EOS_REAL*> np.PyArray_DATA(xVals),
             <EOS_REAL*> np.PyArray_DATA(yVals),
             <EOS_REAL*> np.PyArray_DATA(fVals),
             <EOS_REAL*> np.PyArray_DATA(dFx),
             <EOS_REAL*> np.PyArray_DATA(dFy),
             &error_code)

    if error_code:
        raise RuntimeError( _get_error_message(error_code))
    return fVals, dFx, dFy

cpdef _get_table_info(long table_handle,
                  np.ndarray[int32, ndim=1, mode='c'] info_items):
    cdef EOS_INTEGER error_code
    cdef EOS_INTEGER num_info_items
    cdef np.ndarray[float64, ndim=1, mode='c'] info_vals
    
    info_vals = np.zeros(len(info_items), dtype='float64')

    num_info_items = len(info_items)
    eos_GetTableInfo (<EOS_INTEGER*> &table_handle,
                      &num_info_items,
                      <EOS_INTEGER*> np.PyArray_DATA(info_items),
                      <EOS_REAL*> np.PyArray_DATA(info_vals),
                      &error_code)
    if error_code:
        raise RuntimeError( _get_error_message(error_code))
    return info_vals

cpdef _set_option(long table_handle,
                      long table_option,
                      double table_option_val):

    cdef EOS_INTEGER error_code
    eos_SetOption (<EOS_INTEGER*> &table_handle,
                   <const EOS_INTEGER*> &table_option,
                   <const EOS_REAL*> &table_option_val,
                   &error_code)
    if error_code:
        raise RuntimeError(_get_error_message(error_code))


cpdef _get_table_cmnts(long table_handle):
    cdef EOS_INTEGER error_code
    cdef char* comment_str
    cdef int comment_len
    cdef np.ndarray[int32, ndim=1, mode='c'] info_list
    info_list = np.array([EOS_Cmnt_Len], dtype='int32')

    comment_len = <int> _get_table_info(table_handle, info_list)[0]
    comment_str = <char *> malloc(comment_len*sizeof(char))
    eos_GetTableCmnts (<EOS_INTEGER*> &table_handle,
                       <EOS_CHAR*> comment_str,
                       &error_code)
    if error_code:
        raise RuntimeError(_get_error_message(error_code))
    return comment_str


cpdef _check_extrap(long table_handle,
                  np.ndarray[float64, ndim=1, mode='c'] xVals,
                  np.ndarray[float64, ndim=1, mode='c'] yVals):
    cdef EOS_INTEGER error_code
    cdef EOS_INTEGER nXYPairs
    cdef np.ndarray[int32, ndim=1, mode='c'] xyBounds
    if  len(xVals) != len(yVals):
        raise ValueError('xVals and yVals arguments should be of the same shape!')
    nXYPairs = len(xVals)
    xyBounds = np.zeros((nXYPairs,), dtype='int32')
    eos_CheckExtrap (<EOS_INTEGER*> &table_handle,
                     &nXYPairs,
                     <EOS_REAL*> np.PyArray_DATA(xVals),
                     <EOS_REAL*> np.PyArray_DATA(yVals),
                     <EOS_INTEGER*> np.PyArray_DATA(xyBounds),
                     &error_code)
    if error_code:
        raise RuntimeError(_get_error_message(error_code))
    return xyBounds

