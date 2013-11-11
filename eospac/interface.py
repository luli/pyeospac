#!/usr/bin/python
# -*- coding: utf-8 -*-

def EosMaterial(material=None, tables=['*'], options={},
            spec=['t'], backend='eospac', units='cgs'):
       """
       Parameters:
       -----------
         - material: int: 4 digit SESAME/FEOS material ID
         - tables: list: ['table1_id', 'table2_id', ..etc]
         - options: dict: {'tableid_regexpr': {'opt1': optval}, etc}
            For example {".t_DT": {'create_tzero': True}} would apply the
            EOS_CREATE_TZERO option to both 'Pt_DT' and 'Ut_DT' tables. The
            tableid_regexpr accepts regular expressions.
            [see  help(re) for more details].
         - backend: the backend that should be used ['gamma', 'EOSPAC', 'FEOS',
                                                    'tabulated']
         - units: physical units: 'cgs', 'feos', or 'eospac'
       """
       if backend=='gamma':
           from . import gamma
           return gamma.GammaMaterial(material=material, tables=tables,
                                options=options, spec=spec, units=units)
       elif backend=='eospac':
           from . import eospac as eospac_c
           return eospac_c.EospacMaterial(material=material, tables=tables,
                                options=options, spec=spec, units=units)
       elif backend=='feos':
           from . import feos
           return feos.FeosMaterial(material=material, tables=tables,
                                options=options, spec=spec, units=units)
