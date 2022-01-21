study_site_monthly_EIRs = {
    'Dapelogo':     [1, 1, 1, 1, 1, 3, 27, 70, 130, 57, 29, 1],                                                 # EIR = 322
    "Namawala":     [43.8, 68.5, 27.4, 46.6, 49.4, 24.7, 13.7, 11, 11, 2.74, 13.7, 16.5] ,                      # EIR = 329
    "Dielmo":       [10.4, 13, 6, 2.6, 4, 6, 35, 21, 28, 15, 10.4, 8.4],                                        # EIR = 160
    "Sugungum":     [2, 1, 0.1, 0.1, 0.1, 5, 9.625, 10.375, 20, 23.25, 4.125, 4],
    "Matsari":      [4, 1, 0.1, 0.1, 1, 2.5, 7.75, 11.625, 56.25, 27.125, 12, 6],
    'Siaya':        [55*m for m in [0.07, 0.08, 0.04, 0.02, 0.03, 0.04, 0.22, 0.13, 0.18, 0.09, 0.07, 0.05]],   # EIR =  56
    'Laye':         [1, 1, 1, 1, 1, 7, 8, 9, 5, 4, 6, 1],                                                       # EIR =  45
    'Kilifi_South': [3.6, 1.4, 0.0, 2.8, 0.4, 10.3, 10.6, 3.1, 0.4, 0.0, 0.7, 0.8],                             # EIR =  34
    "Ndiop":        [0.39, 0.19, 0.77, 0, 0, 0, 6.4, 2.2, 4.7, 3.9, 0.87, 0.58],                                # EIR =  20
    "Rafin_Marke":  [1, 1, 0.5, 1, 1, 2, 3.875, 7.75, 15.0, 3.875, 1, 1],
    'Kilifi_North': [0.6, 0.25, 0.0, 0.45, 0.05, 1.7, 1.75, 0.5, 0.05, 0.0, 0.1, 0.15],                         # EIR =   5.6
    'Sukuta':       [3*m for m in [0.07, 0.08, 0.04, 0.02, 0.03, 0.04, 0.22, 0.13, 0.18, 0.09, 0.07, 0.05]],    # EIR =   3.1
    'Bakau':        [0.12*m for m in [0.02, 0.01, 0.04, 0.00, 0.00, 0.00, 0.32, 0.11, 0.24, 0.20, 0.04, 0.03]], # EIR =   0.12
    'Nanoro':       [1.2734122280773272, 0.6893925372248239, 0.47196150548615945, 1.3592037219478288,           # EIR =  22
                     2.3805710805016487, 2.6720894266549253, 2.3689851180114783, 1.895147447566689, 1.4730660393929509,
                     1.453070342033863, 2.033218878887547, 2.201760717967919],
    # Nanoro EIR obtained from calibration to Nanoro clinical cases in infants binned every 3 months
    'Navrongo':     [13.86770953, 1.898386789, 2.606851375, 3.874630107, 6.111886693,
                     10.6982627, 27.14209861, 91.23949981, 158.394485, 96.27332713,
                     36.15764605, 20.49506014],  # average of three years (2001-2003) of monthly EIRs reported in Kasasa et al. 2013. Note: fairly similar to Dapelogo EIRs
    '_constant_':   [1.0/12] * 12,
    '_pulse_':      [1] + [0] * 11
    }

sites_with_interventions = ['Dielmo', 'Ndiop', 'Navrongo']

sites_cm_seek = {
    'Dielmo': 0.5,
    'Ndiop': 0.5,
    'Navrongo': 0.2
}


def mAb_vs_EIR(EIR):
    # Rough cut at function from eyeballing a few BinnedReport outputs parsed into antibody fractions
    mAb = 0.9 * (1e-4*EIR*EIR + 0.7*EIR) / ( 0.7*EIR + 2 )
    return min(mAb, 1.0)