
MATERIALS = {
    'Si': {
        'alat': 10.2, 'cation': 'Si', 'anion': 'Si',
        'cation_pp': '../pseudo/Si_ONCV_PBE-1.2.upf', 'anion_pp': '../pseudo/Si_ONCV_PBE-1.2.upf',
        'cation_mass': 28.086, 'anion_mass': 28.086,
    },
    'GaAs': {
        'alat': 10.45, 'cation': 'Ga', 'anion': 'As',
        'cation_pp': '../pseudo/Ga_ONCV_PBE-1.2.upf', 'anion_pp': '../pseudo/As_ONCV_PBE-1.2.upf',
        'cation_mass': 69.723, 'anion_mass': 74.922,
    },
    'AlAs': {
        'alat': 10.59, 'cation': 'Al', 'anion': 'As',
        'cation_pp': '../pseudo/Al_ONCV_PBE-1.2.upf', 'anion_pp': '../pseudo/As_ONCV_PBE-1.2.upf',
        'cation_mass': 26.982, 'anion_mass': 74.922,
    },
    'GaP': {
        'alat': 10.11, 'cation': 'Ga', 'anion': 'P',
        'cation_pp': '../pseudo/Ga_ONCV_PBE-1.2.upf', 'anion_pp': '../pseudo/P_ONCV_PBE-1.2.upf',
        'cation_mass': 69.723, 'anion_mass': 30.974,
    },
    'AlP': {
        'alat': 10.24, 'cation': 'Al', 'anion': 'P',
        'cation_pp': '../pseudo/Al_ONCV_PBE-1.2.upf', 'anion_pp': '../pseudo/P_ONCV_PBE-1.2.upf',
        'cation_mass': 26.982, 'anion_mass': 30.974,
    },
}

PRL89_VALUES = {
    'GaAs': {'alat': 10.45, 'Z_star': 2.00, 'eps_inf': 11.9, 'eps_static': 13.5,
             'chi2': 134, 'gamma14': -0.40, 'gamma14_0': -1.42},
    'AlAs': {'alat': 10.59, 'Z_star': 2.14, 'eps_inf': 9.6, 'eps_static': 11.5,
             'chi2': 64, 'gamma14': -0.10, 'gamma14_0': -1.40},
    'GaP':  {'alat': 10.11, 'Z_star': 2.10, 'eps_inf': 9.4, 'eps_static': 11.2,
             'chi2': 66, 'gamma14': -0.25, 'gamma14_0': -1.35},
    'AlP':  {'alat': 10.24, 'Z_star': 2.24, 'eps_inf': 8.1, 'eps_static': 10.2,
             'chi2': 39, 'gamma14': 0.05, 'gamma14_0': -1.31},
}

EXPT_VALUES = {
    'GaAs': {'alat': 10.68, 'Z_star': 2.07, 'eps_inf': 10.9, 'eps_static': 13.2,
             'chi2': 166, 'gamma14': -0.32, 'gamma14_0': '-'},
    'AlAs': {'alat': 10.69, 'Z_star': 2.18, 'eps_inf': 8.2, 'eps_static': 10.1,
             'chi2': '-', 'gamma14': '-', 'gamma14_0': '-'},
    'GaP':  {'alat': 10.28, 'Z_star': 2.04, 'eps_inf': 9.0, 'eps_static': 11.1,
             'chi2': 74, 'gamma14': -0.18, 'gamma14_0': '-'},
    'AlP':  {'alat': 10.33, 'Z_star': 2.28, 'eps_inf': 7.5, 'eps_static': 9.8,
             'chi2': '-', 'gamma14': '-', 'gamma14_0': '-'},
}
