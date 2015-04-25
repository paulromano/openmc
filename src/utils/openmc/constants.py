"""Dictionaries of integer-to-string mappings from openmc/src/constants.F90"""

SURFACE_TYPES = {1: 'x-plane',
                 2: 'y-plane',
                 3: 'z-plane',
                 4: 'plane',
                 5: 'x-cylinder',
                 6: 'y-cylinder',
                 7: 'z-cylinder',
                 8: 'sphere',
                 9: 'x-cone',
                 10: 'y-cone',
                 11: 'z-cone'}

BC_TYPES = {0: 'transmission',
            1: 'vacuum',
            2: 'reflective',
            3: 'periodic'}

FILL_TYPES = {1: 'normal',
              2: 'fill',
              3: 'lattice'}

LATTICE_TYPES = {1: 'rectangular',
                 2: 'hexagonal'}

ESTIMATOR_TYPES = {1: 'analog',
                   2: 'tracklength'}

FILTER_TYPES = {1: 'universe',
                2: 'material',
                3: 'cell',
                4: 'cellborn',
                5: 'surface',
                6: 'mesh',
                7: 'energy',
                8: 'energyout',
                9: 'distribcell'}

SCORE_TYPES = {-1: 'flux',
               -2: 'total',
               -3: 'scatter',
               -4: 'nu-scatter',
               -5: 'scatter-n',
               -6: 'scatter-pn',
               -7: 'nu-scatter-n',
               -8: 'nu-scatter-pn',
               -9: 'transport',
               -10: 'n1n',
               -11: 'n2n',
               -12: 'n3n',
               -13: 'n4n',
               -14: 'absorption',
               -15: 'fission',
               -16: 'nu-fission',
               -17: 'kappa-fission',
               -18: 'current',
               -19: 'flux-yn',
               -20: 'total-yn',
               -21: 'scatter-yn',
               -22: 'nu-scatter-yn',
               -23: 'events',
               1: '(n,total)',
               2: '(n,elastic)',
               4: '(n,level)',
               11: '(n,2nd)',
               16: '(n,2n)',
               17: '(n,3n)',
               18: '(n,fission)',
               19: '(n,f)',
               20: '(n,nf)',
               21: '(n,2nf)',
               22: '(n,na)',
               23: '(n,n3a)',
               24: '(n,2na)',
               25: '(n,3na)',
               28: '(n,np)',
               29: '(n,n2a)',
               30: '(n,2n2a)',
               32: '(n,nd)',
               33: '(n,nt)',
               34: '(n,nHe-3)',
               35: '(n,nd2a)',
               36: '(n,nt2a)',
               37: '(n,4n)',
               38: '(n,3nf)',
               41: '(n,2np)',
               42: '(n,3np)',
               44: '(n,n2p)',
               45: '(n,npa)',
               91: '(n,nc)',
               101: '(n,disappear)',
               102: '(n,gamma)',
               103: '(n,p)',
               104: '(n,d)',
               105: '(n,t)',
               106: '(n,3He)',
               107: '(n,a)',
               108: '(n,2a)',
               109: '(n,3a)',
               111: '(n,2p)',
               112: '(n,pa)',
               113: '(n,t2a)',
               114: '(n,d2a)',
               115: '(n,pd)',
               116: '(n,pt)',
               117: '(n,da)',
               201: '(n,Xn)',
               202: '(n,Xgamma)',
               203: '(n,Xp)',
               204: '(n,Xd)',
               205: '(n,Xt)',
               206: '(n,X3He)',
               207: '(n,Xa)',
               444: '(damage)',
               649: '(n,pc)',
               699: '(n,dc)',
               749: '(n,tc)',
               799: '(n,3Hec)',
               849: '(n,tc)'}
SCORE_TYPES.update({MT: '(n,n' + str(MT-50) + ')' for MT in range(51,91)})
SCORE_TYPES.update({MT: '(n,p' + str(MT-600) + ')' for MT in range(600,649)})
SCORE_TYPES.update({MT: '(n,d' + str(MT-650) + ')' for MT in range(650,699)})
SCORE_TYPES.update({MT: '(n,t' + str(MT-700) + ')' for MT in range(700,749)})
SCORE_TYPES.update({MT: '(n,3He' + str(MT-750) + ')' for MT in range(750,649)})
SCORE_TYPES.update({MT: '(n,a' + str(MT-800) + ')' for MT in range(800,849)})
